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

        let (rate_gain, torque_gain) = self.rate_observer_gains(dt);
        let mut gyro_input = super::MeasMat::<3>::zeros().transpose();
        gyro_input
            .fixed_view_mut::<3, 3>(6, 0)
            .copy_from(&(self.state.q.to_rotation_matrix().matrix() * dt));
        gyro_input
            .fixed_view_mut::<3, 3>(9, 0)
            .copy_from(&(nalgebra::Matrix3::identity() * rate_gain));
        gyro_input
            .fixed_view_mut::<3, 3>(21, 0)
            .copy_from(&(self.dyn_params.inertia * torque_gain));
        q += gyro_input * gyro_input.transpose() * (sg * sg / dt);
        let inertia_inv = self
            .dyn_params
            .inertia
            .try_inverse()
            .expect("inertia must be invertible");
        let angular_process =
            inertia_inv * inertia_inv.transpose() * self.noise.torque_noise.powi(2) * dt;
        for i in 0..3 {
            for j in 0..3 {
                q[(i + 9, j + 9)] += angular_process[(i, j)];
            }
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
            q[(i + 21, i + 21)] += qd_t;
        }
        q
    }

    pub(super) fn process_noise_factor(&self, q: &MatN) -> MatN {
        let mut lq = MatN::zeros();

        for axis in 0..3 {
            let pos_idx = axis;
            let vel_idx = axis + 3;
            let q_rr = q[(pos_idx, pos_idx)].max(0.0);
            let q_rv = q[(vel_idx, pos_idx)];
            let q_vv = q[(vel_idx, vel_idx)].max(0.0);

            if q_rr > 0.0 {
                let l11 = q_rr.sqrt();
                let l21 = q_rv / l11;
                let l22_sq = (q_vv - l21 * l21).max(0.0);
                lq[(pos_idx, pos_idx)] = l11;
                lq[(vel_idx, pos_idx)] = l21;
                lq[(vel_idx, vel_idx)] = l22_sq.sqrt();
            } else {
                lq[(vel_idx, vel_idx)] = q_vv.sqrt();
            }
        }

        // Factor only the coupled attitude/rate/torque block. A zero pivot
        // is allowed for exact zero-noise configurations without bias jitter.
        let angular_indices = [6, 7, 8, 9, 10, 11, 21, 22, 23];
        for (row, &i) in angular_indices.iter().enumerate() {
            for (column, &j) in angular_indices.iter().enumerate().take(row + 1) {
                let mut value = q[(i, j)];
                for &k in angular_indices.iter().take(column) {
                    value -= lq[(i, k)] * lq[(j, k)];
                }
                lq[(i, j)] = if row == column {
                    value.max(0.0).sqrt()
                } else if lq[(j, j)] > 0.0 {
                    value / lq[(j, j)]
                } else {
                    0.0
                };
            }
        }
        for diag_idx in 12..15 {
            lq[(diag_idx, diag_idx)] = q[(diag_idx, diag_idx)].max(0.0).sqrt();
        }
        for diag_idx in 15..18 {
            lq[(diag_idx, diag_idx)] = q[(diag_idx, diag_idx)].max(0.0).sqrt();
        }
        for diag_idx in 18..21 {
            lq[(diag_idx, diag_idx)] = q[(diag_idx, diag_idx)].max(0.0).sqrt();
        }

        lq
    }
}
