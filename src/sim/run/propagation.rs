use super::SimRun;
use crate::types::ImuSample;
use nalgebra::{UnitQuaternion, Vector3};
use rand_distr::Distribution;

impl SimRun {
    pub(super) fn update_truth_residuals(&mut self, t: f64) {
        let in_blk = (t >= self.cfg.blackout_start) && (t <= self.cfg.blackout_end);
        let a_rw = if in_blk {
            self.noise
                .aero_force_rw_blackout
                .unwrap_or(self.noise.aero_force_rw)
        } else {
            self.noise.aero_force_rw
        };
        let t_rw = if in_blk {
            self.noise
                .aero_torque_rw_blackout
                .unwrap_or(self.noise.aero_torque_rw)
        } else {
            self.noise.aero_torque_rw
        };
        let s_f = a_rw * (1.0 - self.phi_f_const * self.phi_f_const).sqrt();
        let s_t = t_rw * (1.0 - self.phi_t_const * self.phi_t_const).sqrt();
        let n_f = Vector3::new(
            s_f * self.standard_normal.sample(&mut self.rng),
            s_f * self.standard_normal.sample(&mut self.rng),
            s_f * self.standard_normal.sample(&mut self.rng),
        );
        let n_t = Vector3::new(
            s_t * self.standard_normal.sample(&mut self.rng),
            s_t * self.standard_normal.sample(&mut self.rng),
            s_t * self.standard_normal.sample(&mut self.rng),
        );
        self.truth.aero_res_force = self.truth.aero_res_force * self.phi_f_const + n_f;
        self.truth.aero_res_torque = self.truth.aero_res_torque * self.phi_t_const + n_t;
    }

    pub(super) fn propagate_estimator(
        &mut self,
        t: f64,
        tau_b: Vector3<f64>,
        f_b: Vector3<f64>,
        aero_rt: &crate::physics::AeroResult,
    ) {
        let acc_b = (aero_rt.force_b + f_b + self.truth.aero_res_force) / self.dyn_params.mass;
        let imu_true_acc = acc_b;
        let imu_gyro = self.truth.omega;
        let acc_meas = imu_true_acc
            + self.truth.ba
            + Vector3::new(
                self.noise.accel_noise * self.standard_normal.sample(&mut self.rng),
                self.noise.accel_noise * self.standard_normal.sample(&mut self.rng),
                self.noise.accel_noise * self.standard_normal.sample(&mut self.rng),
            );
        let gyro_meas = imu_gyro
            + self.truth.bg
            + Vector3::new(
                self.noise.gyro_noise * self.standard_normal.sample(&mut self.rng),
                self.noise.gyro_noise * self.standard_normal.sample(&mut self.rng),
                self.noise.gyro_noise * self.standard_normal.sample(&mut self.rng),
            );
        self.est.propagate(ImuSample {
            t: t + self.cfg.dt,
            dt: self.cfg.dt,
            accel_meas: acc_meas,
            gyro_meas,
            torque_cmd: tau_b,
            force_cmd: f_b,
        });
    }

    pub(super) fn propagate_truth(
        &mut self,
        tau_b: Vector3<f64>,
        f_b: Vector3<f64>,
        aero_rt: &crate::physics::AeroResult,
    ) {
        self.truth.omega += self.inertia_inv
            * (tau_b + aero_rt.torque_b + self.truth.aero_res_torque
                - self
                    .truth
                    .omega
                    .cross(&(self.dyn_params.inertia * self.truth.omega)))
            * self.cfg.dt;
        self.truth.q = UnitQuaternion::new_normalize(
            (self.truth.q * UnitQuaternion::from_scaled_axis(self.truth.omega * self.cfg.dt))
                .into_inner(),
        );
        let force_nav = self.truth.q.to_rotation_matrix()
            * (aero_rt.force_b + f_b + self.truth.aero_res_force)
            / self.dyn_params.mass
            + self.est.gravity(&self.truth.r);
        self.truth.r += self.truth.v * self.cfg.dt + 0.5 * force_nav * self.cfg.dt * self.cfg.dt;
        self.truth.v += force_nav * self.cfg.dt;
    }

    pub(super) fn record_truth(&mut self, t: f64) {
        self.truth_times.push(t);
        self.truth_hist.push(self.truth.clone());
    }
}

