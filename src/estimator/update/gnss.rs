use super::{
    history::Event, math::chol_lower_with_enforce, robust::ROBUST_SENSOR_GNSS_PVT, Estimator,
    GnssMeas, MatN, MeasMat,
};
use crate::types::STATE_DIM;
use crate::util::enforce_pd;
use nalgebra::{linalg::Cholesky, Matrix6, SMatrix, SVector};

type Gain6 = SMatrix<f64, STATE_DIM, 6>;

impl Estimator {
    pub fn update_gnss(&mut self, meas: GnssMeas) {
        self.process_gnss(meas, true, true);
    }

    pub fn force_align_gnss(&mut self, meas: GnssMeas) {
        if let Some(events) = self.rewind_if_needed(meas.t) {
            if std::env::var("TRACE").ok().is_some() {
                eprintln!(
                    "[force_align] rewind to t={:.3}; set r,v to meas; before |r-|",
                    self.state.t
                );
            }
            self.state.r = meas.pos;
            self.state.v = meas.vel;
            let mut p = self.cov();
            for i in 0..3 {
                p[(i, i)] = (300.0_f64).powi(2);
            }
            for i in 0..3 {
                p[(i + 3, i + 3)] = (10.0_f64).powi(2);
            }
            enforce_pd(&mut p);
            self.l_cov = chol_lower_with_enforce(p, "force_align_gnss (rewind)");
            self.history.push_event(Event::Gnss(meas));
            if std::env::var("TRACE").ok().is_some() {
                eprintln!(
                    "[force_align] after set r,v; now replay {} events forward",
                    events.len()
                );
            }
            for e in events {
                self.replay_event(e);
            }
            return;
        }

        self.state.r = meas.pos;
        self.state.v = meas.vel;
        let mut p = self.cov();
        for i in 0..3 {
            p[(i, i)] = (300.0_f64).powi(2);
        }
        for i in 0..3 {
            p[(i + 3, i + 3)] = (10.0_f64).powi(2);
        }
        enforce_pd(&mut p);
        self.l_cov = chol_lower_with_enforce(p, "force_align_gnss (no rewind)");
        self.history.push_event(Event::Gnss(meas));
    }

    pub(super) fn process_gnss(&mut self, meas: GnssMeas, record: bool, allow_rewind: bool) {
        let trace = std::env::var("TRACE").ok().is_some();
        let replay_events = if allow_rewind {
            self.rewind_if_needed(meas.t)
        } else {
            None
        };
        let mut h = MeasMat::<6>::zeros();
        for i in 0..3 {
            h[(i, i)] = 1.0;
            h[(i + 3, i + 3)] = 1.0;
        }
        let z = SVector::<f64, 6>::from_iterator(meas.pos.iter().chain(meas.vel.iter()).cloned());
        let hx = SVector::<f64, 6>::from_iterator(
            self.state.r.iter().chain(self.state.v.iter()).cloned(),
        );
        let innov = z - hx;

        let p_prior = self.l_cov * self.l_cov.transpose();
        self.last_gnss_prior_t = Some(meas.t);
        self.last_gnss_prior_rv_p = Some(p_prior.fixed_view::<6, 6>(0, 0).into_owned());
        self.last_gnss_prior_r = Some(self.state.r);
        self.last_gnss_prior_v = Some(self.state.v);

        let s = h * p_prior * h.transpose() + meas.r_cov;
        let s_chol = match Cholesky::new(s) {
            Some(c) => c,
            None => {
                self.replay_tail(replay_events);
                return;
            }
        };
        let s_inv = s_chol.inverse();
        let base_d2 = (innov.transpose() * s_inv * innov)[(0, 0)];
        let r_cov = self.robust_scale(ROBUST_SENSOR_GNSS_PVT, 6, base_d2, meas.r_cov);

        let s = h * p_prior * h.transpose() + r_cov;
        let s_chol = match Cholesky::new(s) {
            Some(c) => c,
            None => {
                self.replay_tail(replay_events);
                return;
            }
        };
        let s_inv = s_chol.inverse();
        let d2 = (innov.transpose() * s_inv.clone() * innov)[(0, 0)];
        self.last_gnss_nis = Some((meas.t, d2, 6));
        if trace {
            self.trace_gnss_pre(meas.t, &innov, d2);
        }

        let gate_sq = self.robust_gate_sq(6, self.noise.nis_alpha);
        let accepted = d2.is_finite() && d2 <= gate_sq;
        self.last_gnss_accepted = Some((meas.t, accepted));
        if !accepted {
            self.gnss_reject_streak = self.gnss_reject_streak.saturating_add(1);
            if self.gnss_reject_streak >= 3 {
                self.inflate_rv(1_000.0, 100.0);
                self.gnss_reject_streak = 0;
            }
            self.replay_tail(replay_events);
            return;
        }
        self.gnss_reject_streak = 0;
        self.update_iae_scale(ROBUST_SENSOR_GNSS_PVT, 6, d2);

        self.bias_cal.from_pvt = true;
        self.bias_cal_observe_rv(meas.t, meas.pos, meas.vel, &r_cov);

        let mut k = p_prior * h.transpose() * s_inv;
        self.apply_gnss_active_set(meas.t, &innov, &mut k);

        let dx = &k * &innov;
        self.apply_state_error(dx);

        let i_kh = MatN::identity() - k * h;
        let mut p_new = i_kh * p_prior * i_kh.transpose() + k * r_cov * k.transpose();
        enforce_pd(&mut p_new);
        self.l_cov = chol_lower_with_enforce(p_new, "process_gnss");

        if trace {
            self.trace_gnss_post(meas.t);
        }

        let p_post = self.cov();
        self.last_gnss_post_t = Some(meas.t);
        self.last_gnss_post_rv_p = Some(p_post.fixed_view::<6, 6>(0, 0).into_owned());
        self.last_gnss_post_r = Some(self.state.r);
        self.last_gnss_post_v = Some(self.state.v);

        self.update_aero_learning(meas.vel, meas.t);

        self.state.t = meas.t;
        self.maybe_checkpoint(true);
        if record {
            self.history.push_event(Event::Gnss(meas));
        }
        self.replay_tail(replay_events);
    }

    fn apply_gnss_active_set(&self, t: f64, innov: &SVector<f64, 6>, k: &mut Gain6) {
        let (mut bias_updates_enabled, max_ba_update) = self.bias_update_profile(t);
        if self.bias_cal.enabled {
            bias_updates_enabled = false;
        }

        for i in 6..STATE_DIM {
            let active_ba = bias_updates_enabled && (15..18).contains(&i);
            if active_ba {
                continue;
            }
            for j in 0..6 {
                k[(i, j)] = 0.0;
            }
        }

        if !bias_updates_enabled {
            return;
        }

        for i in 15..18 {
            for j in 0..3 {
                k[(i, j)] = 0.0;
            }
        }
        let gain = self.noise.gps_ba_update_gain.clamp(0.0, 1.0);
        if (gain - 1.0).abs() > 0.0 {
            for i in 15..18 {
                for j in 0..6 {
                    k[(i, j)] *= gain;
                }
            }
        }

        let dx_raw = (&*k) * innov;
        for i in 15..18 {
            let abs = dx_raw[i].abs();
            if abs > max_ba_update && abs > 0.0 {
                let s = max_ba_update / abs;
                for j in 0..6 {
                    k[(i, j)] *= s;
                }
            }
        }
    }

    fn trace_gnss_pre(&self, t: f64, innov: &SVector<f64, 6>, d2: f64) {
        let p_rv = self.last_gnss_prior_rv_p.unwrap_or(Matrix6::identity());
        let mut sig = [0.0f64; 6];
        for i in 0..6 {
            sig[i] = p_rv[(i, i)].abs().sqrt();
        }
        let pos_innov = innov.fixed_rows::<3>(0).norm();
        let vel_innov = innov.fixed_rows::<3>(3).norm();
        eprintln!(
            "[GNSS pre] t={:.3} sig_r=({:.2},{:.2},{:.2}) m sig_v=({:.3},{:.3},{:.3}) m/s |dz_r|={:.2} |dz_v|={:.3} NIS={:.2}",
            t,
            sig[0],
            sig[1],
            sig[2],
            sig[3],
            sig[4],
            sig[5],
            pos_innov,
            vel_innov,
            d2
        );
    }

    fn trace_gnss_post(&self, t: f64) {
        let p_post = self.cov();
        let p_rv = p_post.fixed_view::<6, 6>(0, 0).into_owned();
        let mut sig = [0.0f64; 6];
        for i in 0..6 {
            sig[i] = p_rv[(i, i)].abs().sqrt();
        }
        eprintln!(
            "[GNSS post] t={:.3} sig_r=({:.2},{:.2},{:.2}) m sig_v=({:.3},{:.3},{:.3}) m/s",
            t, sig[0], sig[1], sig[2], sig[3], sig[4], sig[5]
        );
    }
}
