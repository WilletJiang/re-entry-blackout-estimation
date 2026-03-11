use super::{
    history::Event, math::chol_lower_with_enforce, robust::ROBUST_SENSOR_STAR, Estimator, MatN,
    MeasMat, StarMeas, VecN,
};
use crate::types::STATE_DIM;
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
                self.replay_tail(replay_events);
                return;
            }
        };
        let base_s_inv = base_s_chol.inverse();
        let base_d2 = (angle_axis.transpose() * base_s_inv * angle_axis)[(0, 0)];
        let meas_r_cov = self.robust_scale(ROBUST_SENSOR_STAR, 3, base_d2, meas.r_cov);
        let s = h.clone() * p.clone() * h.transpose() + meas_r_cov;
        let s_chol = match Cholesky::new(s.clone()) {
            Some(c) => c,
            None => {
                self.replay_tail(replay_events);
                return;
            }
        };
        let s_inv = s_chol.inverse();
        let d2 = (angle_axis.transpose() * s_inv.clone() * angle_axis)[(0, 0)];
        let gate_sq = self.robust_gate_sq_or(3, self.noise.star_nis_gate_sq);
        if !d2.is_finite() || d2 > gate_sq {
            self.replay_tail(replay_events);
            return;
        }
        self.update_iae_scale(ROBUST_SENSOR_STAR, 3, d2);

        let k = p.clone() * h.transpose() * s_inv;
        let dx = k.clone() * angle_axis;
        let mut dx_att_only = VecN::zeros();
        for i in 6..9 {
            dx_att_only[i] = dx[i];
        }
        self.apply_state_error(dx_att_only);

        let mut k_att = k.clone();
        for i in 0..6 {
            for j in 0..3 {
                k_att[(i, j)] = 0.0;
            }
        }
        for i in 9..STATE_DIM {
            for j in 0..3 {
                k_att[(i, j)] = 0.0;
            }
        }
        let i_kh = MatN::identity() - (k_att.clone() * h.clone());
        let mut p_new =
            i_kh * p * i_kh.transpose() + k_att.clone() * meas_r_cov * k_att.transpose();
        enforce_pd(&mut p_new);
        self.l_cov = chol_lower_with_enforce(p_new, "process_star");

        self.maybe_checkpoint(true);
        if record {
            self.history.push_event(Event::Star(meas));
        }
        self.replay_tail(replay_events);
    }
}
