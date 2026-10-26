use super::{
    history::Event, math::chol_lower_with_enforce, robust::ROBUST_SENSOR_STAR, Estimator, MatN,
    MeasMat, StarMeas,
};
use crate::util::enforce_pd;
use nalgebra::linalg::Cholesky;

impl Estimator {
    pub fn update_star(&mut self, meas: StarMeas) {
        self.process_star(meas, true, true);
    }

    pub(super) fn process_star(&mut self, meas: StarMeas, record: bool, allow_rewind: bool) {
        let replay_events = if allow_rewind {
            self.rewind_if_needed(meas.t)
        } else {
            None
        };
        let predicted_end_state = if self.seg_has_data {
            // RTS needs the propagation prior at this timestamp after the update
            // mutates `self.state`, so the pre-update state must be owned here.
            Some(self.state.clone())
        } else {
            None
        };
        if std::env::var("TRACE").ok().is_some() {
            eprintln!("[star] meas.t={:.3} state.t={:.3}", meas.t, self.state.t);
        }
        let dq = if allow_rewind {
            meas.q_bi * self.state.q.inverse()
        } else if meas.t < self.state.t {
            if let Some(idx) = self.history.find_snapshot_index(meas.t) {
                meas.q_bi * self.history.snapshots[idx].state.q.inverse()
            } else {
                meas.q_bi * self.state.q.inverse()
            }
        } else {
            meas.q_bi * self.state.q.inverse()
        };
        let angle_axis = dq.scaled_axis();
        let mut h = MeasMat::<3>::zeros();
        for i in 0..3 {
            h[(i, i + 6)] = 1.0;
        }
        let p = self.l_cov * self.l_cov.transpose();
        let base_s = h.clone() * p.clone() * h.transpose() + meas.r_cov;
        let base_s_chol = match Cholesky::new(base_s.clone()) {
            Some(c) => c,
            None => {
                if record {
                    self.history.push_event(Event::Star(meas));
                }
                self.replay_tail(replay_events);
                return;
            }
        };
        let base_s_inv = base_s_chol.inverse();
        let base_d2 = (angle_axis.transpose() * base_s_inv * angle_axis)[(0, 0)];
        self.last_star_base_nis = Some((meas.t, base_d2, 3));
        let meas_r_cov = self.robust_scale(ROBUST_SENSOR_STAR, 3, base_d2, meas.r_cov);
        let s = h.clone() * p.clone() * h.transpose() + meas_r_cov;
        let s_chol = match Cholesky::new(s.clone()) {
            Some(c) => c,
            None => {
                if record {
                    self.history.push_event(Event::Star(meas));
                }
                self.replay_tail(replay_events);
                return;
            }
        };
        let s_inv = s_chol.inverse();
        let d2 = (angle_axis.transpose() * s_inv.clone() * angle_axis)[(0, 0)];
        let gate_sq = self.robust_gate_sq_or(3, self.noise.star_nis_gate_sq);
        if !d2.is_finite() || d2 > gate_sq {
            // An initially contaminated attitude can be overconfident. Three
            // consistent rejected innovations indicate a lost attitude solution,
            // whereas isolated or inconsistent Cauchy outliers do not.
            if angle_axis.iter().all(|v| v.is_finite()) {
                let consistent = self
                    .star_reject_innovation
                    .is_some_and(|previous| (angle_axis - previous).norm() < 0.002);
                self.star_reject_streak = if consistent {
                    self.star_reject_streak + 1
                } else {
                    1
                };
                self.star_reject_innovation = Some(angle_axis);
                if self.star_reject_streak >= 3 {
                    let mut widened = self.cov();
                    for i in 6..9 {
                        widened[(i, i)] += 5.0_f64.to_radians().powi(2);
                    }
                    self.l_cov = chol_lower_with_enforce(widened, "star reacquisition");
                    self.star_reject_streak = 0;
                    self.star_reject_innovation = None;
                }
            }
            if record {
                self.history.push_event(Event::Star(meas));
            }
            self.replay_tail(replay_events);
            return;
        }
        self.star_reject_streak = 0;
        self.star_reject_innovation = None;
        self.update_iae_scale(ROBUST_SENSOR_STAR, 3, d2);

        let k = p.clone() * h.transpose() * s_inv;
        self.apply_state_error(k * angle_axis);
        let k_att = k;
        let i_kh = MatN::identity() - (k_att.clone() * h.clone());
        let mut p_new =
            i_kh * p * i_kh.transpose() + k_att.clone() * meas_r_cov * k_att.transpose();
        enforce_pd(&mut p_new);
        self.l_cov = chol_lower_with_enforce(p_new, "process_star");

        self.maybe_checkpoint_with_prediction(true, predicted_end_state);
        if record {
            self.history.push_event(Event::Star(meas));
        }
        self.replay_tail(replay_events);
    }
}
