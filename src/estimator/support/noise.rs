use super::{Estimator, MatN};

impl Estimator {
    pub(super) fn process_noise(&self, dt: f64) -> MatN {
        let mut q = MatN::zeros();
        let sg = self.noise.gyro_noise;
        let sa = self.noise.accel_noise;
        let sbg = self.noise.gyro_bias_rw;
        let sba = self.noise.accel_bias_rw;
        let s_af = self.noise.accel_residual;
        let (s_ar, s_tr) = {
            let base_f = self.noise.aero_force_rw;
            let base_t = self.noise.aero_torque_rw;
            if let Some((bs, be)) = self.blackout_window {
                let t = self.state.t;
                if t >= bs && t <= be {
                    (
                        self.aero_force_rw_blackout.unwrap_or(base_f),
                        self.aero_torque_rw_blackout.unwrap_or(base_t),
                    )
                } else {
                    (base_f, base_t)
                }
            } else {
                (base_f, base_t)
            }
        };

        let t = self.state.t;
        let (rv_scale, ba_scale) = if let Some((bs, be)) = self.blackout_window {
            if t < bs {
                (self.noise.q_scale_rv_pre, self.noise.q_scale_ba_pre)
            } else if t <= be {
                (self.noise.q_scale_rv_blk, self.noise.q_scale_ba_blk)
            } else {
                (self.noise.q_scale_rv_post, self.noise.q_scale_ba_post)
            }
        } else {
            (1.0, 1.0)
        };
        let rv_scale = rv_scale.max(0.0);
        let ba_scale = ba_scale.max(0.0);

        let in_star_saturation = if let Some((bs, _)) = self.blackout_window {
            t >= bs + 10.0 && t <= bs + 30.0
        } else {
            false
        };

        let sg_effective = if in_star_saturation { sg * 10.0 } else { sg };

        for i in 0..3 {
            q[(i + 6, i + 6)] = sg_effective * sg_effective * dt;
        }

        let s_af_effective = if in_star_saturation { s_af * 5.0 } else { s_af };
        let q_acc = sa * sa + s_af_effective * s_af_effective;
        let dt2 = dt * dt;
        let q_rr = rv_scale * q_acc * (dt2 * dt) / 3.0;
        let q_rv = rv_scale * q_acc * (dt2 * 0.5);
        let q_vv = rv_scale * q_acc * dt;
        for i in 0..3 {
            q[(i, i)] += q_rr;
            q[(i, i + 3)] += q_rv;
            q[(i + 3, i)] += q_rv;
            q[(i + 3, i + 3)] += q_vv;
        }

        for i in 0..3 {
            q[(i + 12, i + 12)] = sbg * sbg * dt;
            q[(i + 15, i + 15)] = ba_scale * sba * sba * dt;
        }

        let phi_f = (-dt / self.dyn_params.tau_aero_force).exp();
        let phi_t = (-dt / self.dyn_params.tau_aero_torque).exp();
        let qd_f = s_ar * s_ar * (1.0 - phi_f * phi_f);
        let qd_t = s_tr * s_tr * (1.0 - phi_t * phi_t);
        for i in 0..3 {
            q[(i + 18, i + 18)] = qd_f;
        }
        for i in 0..3 {
            q[(i + 21, i + 21)] = qd_t;
        }
        q
    }
}
