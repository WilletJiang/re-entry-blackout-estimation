use super::{
    history::{Event, LoggedImu},
    math::{chol_lower_with_enforce, insert_block, qr_r_upper_in_place},
    Estimator, Mat2NByN, MatN,
};
use crate::physics::nominal_body_dynamics_with_residual;
use crate::types::{ImuSample, STATE_DIM};
use crate::util::skew;
use nalgebra::{Matrix3, UnitQuaternion, Vector3};

impl Estimator {
    pub fn propagate(&mut self, imu: ImuSample) {
        self.propagate_internal(imu, true);
    }

    pub(super) fn propagate_internal(&mut self, imu: ImuSample, record: bool) {
        let dt = imu.dt;
        let dt2 = dt * dt;
        let half_dt2 = 0.5 * dt2;
        let bg = self.state.bg;
        let ba = self.state.ba;
        let omega_meas = imu.gyro_meas - bg;
        let accel_meas = imu.accel_meas - ba;

        self.state.omega = omega_meas;

        let dq = UnitQuaternion::from_scaled_axis(omega_meas * dt);
        self.state.q = UnitQuaternion::new_normalize((self.state.q * dq).into_inner());

        let c_bi = self.state.q.to_rotation_matrix();
        let g = self.gravity(&self.state.r);
        let g_grad = self.gravity_gradient(&self.state.r);
        let delta_r = self.state.v * dt;
        let g_avg = g + 0.5 * g_grad * delta_r;
        let linearize_specific_force_b =
            self.nominal_specific_force_b(imu.force_cmd, imu.torque_cmd, accel_meas);

        let aero_correction = if self.aero_learning_enabled {
            self.learned_aero_correction
        } else {
            Vector3::zeros()
        };

        if self.bias_cal.enabled && self.bias_cal.ins_inited {
            let g_ins = self.gravity(&self.bias_cal.ins_r);
            let a_ins = c_bi * accel_meas + g_ins;
            let ins_v_prev = self.bias_cal.ins_v;
            self.bias_cal.ins_r += ins_v_prev * dt + a_ins * half_dt2;
            self.bias_cal.ins_v = ins_v_prev + a_ins * dt;

            let c = c_bi.matrix();
            self.bias_cal.sens_r_ba += self.bias_cal.sens_v_ba * dt - (*c) * half_dt2;
            self.bias_cal.sens_v_ba += -(*c) * dt;
        }

        let specific_force_nav = c_bi * accel_meas + g_avg + aero_correction;
        let v_prev = self.state.v;
        self.state.r += v_prev * dt + specific_force_nav * half_dt2;
        self.state.v = v_prev + specific_force_nav * dt;

        let phi_f = (-dt / self.dyn_params.tau_aero_force).exp();
        let phi_t = (-dt / self.dyn_params.tau_aero_torque).exp();
        self.state.aero_res_force *= phi_f;
        self.state.aero_res_torque *= phi_t;
        self.state.t = imu.t;

        let (phi, qk) = self.linearize(dt, linearize_specific_force_b);
        self.time_update_sqrt(&phi, &qk);
        self.accumulate_segment(&phi, &qk);
        if record {
            self.history
                .push_event(Event::Imu(LoggedImu { sample: imu }));
        }
        self.maybe_checkpoint(false);
    }

    pub(super) fn linearize(&self, dt: f64, specific_force_b: Vector3<f64>) -> (MatN, MatN) {
        let mut phi = MatN::identity();
        let dt2 = dt * dt;
        let half_dt2 = 0.5 * dt2;
        for i in 0..3 {
            phi[(i, i + 3)] = dt;
        }

        let g_grad = self.gravity_gradient(&self.state.r);
        let dv_dr = g_grad * dt;
        insert_block(&mut phi, &dv_dr, 3, 0);

        let c = self.state.q.to_rotation_matrix().matrix().clone();
        let a_nav = c * specific_force_b;
        let dr_datt = -skew(a_nav) * half_dt2;
        insert_block(&mut phi, &dr_datt, 0, 6);
        insert_block(&mut phi, &(-c * half_dt2), 0, 15);
        let dr_daf = (c / self.dyn_params.mass) * half_dt2;
        insert_block(&mut phi, &dr_daf, 0, 18);

        let dv_datt = -skew(a_nav) * dt;
        insert_block(&mut phi, &dv_datt, 3, 6);
        insert_block(&mut phi, &(-c * dt), 3, 15);
        let dv_daf = (c / self.dyn_params.mass) * dt;
        insert_block(&mut phi, &dv_daf, 3, 18);

        for i in 0..3 {
            phi[(i + 6, i + 12)] = -dt;
        }
        let omega_x = -skew(self.state.omega) * dt;
        insert_block(&mut phi, &(Matrix3::identity() + omega_x), 6, 6);

        let phi_f = (-dt / self.dyn_params.tau_aero_force).exp();
        let phi_t = (-dt / self.dyn_params.tau_aero_torque).exp();
        insert_block(&mut phi, &(Matrix3::identity() * phi_f), 18, 18);
        insert_block(&mut phi, &(Matrix3::identity() * phi_t), 21, 21);
        let qk = self.process_noise(dt);
        (phi, qk)
    }

    pub(super) fn time_update_sqrt(&mut self, phi: &MatN, qk: &MatN) {
        let a = phi * self.l_cov;
        let lq = chol_lower_with_enforce(*qk, "time_update_sqrt Qk");

        let mut mt = Mat2NByN::zeros();
        for r in 0..STATE_DIM {
            for c in 0..STATE_DIM {
                mt[(r, c)] = a[(c, r)];
                mt[(r + STATE_DIM, c)] = lq[(c, r)];
            }
        }

        let mut r = qr_r_upper_in_place(&mut mt);
        for i in 0..STATE_DIM {
            if r[(i, i)] < 0.0 {
                for j in i..STATE_DIM {
                    r[(i, j)] = -r[(i, j)];
                }
            }
        }

        self.l_cov = r.transpose();
    }

    #[cfg(feature = "diagnostics")]
    pub fn bench_time_update_sqrt(&mut self, phi: &MatN, qk: &MatN) {
        self.time_update_sqrt(phi, qk);
    }

    fn nominal_specific_force_b(
        &self,
        force_cmd_b: Vector3<f64>,
        torque_cmd_b: Vector3<f64>,
        imu_specific_force_b: Vector3<f64>,
    ) -> Vector3<f64> {
        let Some(table) = self.aero_table.as_ref() else {
            return imu_specific_force_b;
        };
        let nominal = nominal_body_dynamics_with_residual(
            table.as_ref(),
            &self.dyn_params,
            self.state.r,
            self.state.v,
            &self.state.q,
            self.state.omega,
            force_cmd_b,
            torque_cmd_b,
            self.state.aero_res_force,
            self.state.aero_res_torque,
        );
        let specific_force_b = nominal.specific_force_b;
        if specific_force_b.iter().all(|value| value.is_finite()) {
            return specific_force_b;
        }
        imu_specific_force_b
    }
}
