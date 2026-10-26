use super::SimRun;
use crate::physics::NominalBodyDynamics;
use crate::sim::truth::normal_vector;
use crate::types::ImuSample;
use nalgebra::{UnitQuaternion, Vector3};

impl SimRun {
    pub(super) fn update_truth_residuals(&mut self, t: f64) {
        let in_blk = (t >= self.cfg.blackout_start) && (t <= self.cfg.blackout_end);
        let noise = &self.cfg.truth_noise;
        let force_sigma = if in_blk {
            noise.blackout_force_sigma
        } else {
            noise.aero_force_sigma
        };
        let torque_sigma = if in_blk {
            noise.blackout_torque_sigma
        } else {
            noise.aero_torque_sigma
        };
        let n_f = normal_vector(&mut self.truth_rng)
            * (force_sigma * (1.0 - self.phi_f_const * self.phi_f_const).sqrt());
        let n_t = normal_vector(&mut self.truth_rng)
            * (torque_sigma * (1.0 - self.phi_t_const * self.phi_t_const).sqrt());
        self.truth.aero_res_force = self.truth.aero_res_force * self.phi_f_const + n_f;
        self.truth.aero_res_torque = self.truth.aero_res_torque * self.phi_t_const + n_t;
    }

    pub(super) fn propagate_estimator(
        &mut self,
        t: f64,
        tau_b: Vector3<f64>,
        f_b: Vector3<f64>,
        nominal_truth: &NominalBodyDynamics,
    ) {
        let (accel_noise, gyro_noise) = self
            .cfg
            .truth_noise
            .sample_imu_noise(self.cfg.dt, &mut self.imu_rng);
        let acc_meas = nominal_truth.specific_force_b + self.truth.ba + accel_noise;
        let angular_acceleration = self.inertia_inv
            * (nominal_truth.total_torque_b
                - self
                    .truth
                    .omega
                    .cross(&(self.dyn_params.inertia * self.truth.omega)));
        let gyro_meas = self.truth.omega
            + angular_acceleration * (0.5 * self.cfg.dt)
            + self.truth.bg
            + gyro_noise;
        self.est.propagate(ImuSample {
            t: t + self.cfg.dt,
            dt: self.cfg.dt,
            accel_meas: acc_meas,
            gyro_meas,
            torque_cmd: tau_b,
            force_cmd: f_b,
        });
    }

    pub(super) fn propagate_truth(&mut self, nominal_truth: &NominalBodyDynamics) {
        let previous_omega = self.truth.omega;
        self.truth.omega += self.inertia_inv
            * (nominal_truth.total_torque_b
                - self
                    .truth
                    .omega
                    .cross(&(self.dyn_params.inertia * self.truth.omega)))
            * self.cfg.dt;
        self.truth.q = UnitQuaternion::new_normalize(
            (self.truth.q
                * UnitQuaternion::from_scaled_axis(
                    (previous_omega + self.truth.omega) * (0.5 * self.cfg.dt),
                ))
            .into_inner(),
        );
        let force_nav = self.truth.q.to_rotation_matrix() * nominal_truth.specific_force_b
            + self.est.gravity(&self.truth.r);
        self.truth.r += self.truth.v * self.cfg.dt + 0.5 * force_nav * self.cfg.dt * self.cfg.dt;
        self.truth.v += force_nav * self.cfg.dt;
    }

    pub(super) fn record_truth(&mut self, t: f64) {
        self.truth_times.push(t);
        self.truth_hist.push(self.truth.clone());
    }
}
