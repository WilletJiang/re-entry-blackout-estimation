use super::{
    bias_cal, history::Event, history::HistoryBuffer, history::PropRecord, history::Snapshot,
    Estimator, EstimatorHealth, MatN,
};
use crate::data::AeroTable;
use crate::util::enforce_pd;
use nalgebra::{Matrix3, Vector3};
use std::sync::Arc;

impl Estimator {
    pub fn new(
        noise: crate::types::NoiseParams,
        dyn_params: crate::types::DynamicsParams,
        history_sec: f64,
    ) -> Self {
        let cov0 = bias_cal::initial_cov();
        let l0 = super::math::chol_lower_with_enforce(cov0, "initial_cov");
        let history_snapshot_dt = 0.01;
        let mut est = Self {
            state: crate::types::State::new(),
            l_cov: l0,
            noise,
            dyn_params,
            history: HistoryBuffer::new(history_sec.max(1.0), history_snapshot_dt),
            history_snapshot_dt,
            last_snapshot_t: 0.0,
            seg_phi: MatN::identity(),
            seg_q: MatN::zeros(),
            seg_has_data: false,
            aero_table: None,
            blackout_window: None,
            aero_force_rw_blackout: None,
            aero_torque_rw_blackout: None,
            last_gnss_nis: None,
            last_gnss_prior_t: None,
            last_gnss_prior_rv_p: None,
            last_gnss_prior_r: None,
            last_gnss_prior_v: None,
            last_gnss_post_t: None,
            last_gnss_post_rv_p: None,
            last_gnss_post_r: None,
            last_gnss_post_v: None,
            last_gnss_accepted: None,
            gnss_reject_streak: 0,
            robust_gate_alpha: None,
            robust_r_scale_cap: super::robust::DEFAULT_R_SCALE_CAP,
            robust_student_nu: 0.0,
            iae_gamma: None,
            iae_scales: [1.0; super::robust::ROBUST_SENSOR_COUNT],
            aero_learning_enabled: false,
            last_v_meas: None,
            learned_aero_correction: Vector3::zeros(),
            aero_learning_rate: 0.02,
            allow_bias_precalibration: true,
            precalibration_end_time: 120.0,
            bias_cal: bias_cal::BiasCal::new(true),
        };

        est.history.push_snapshot(Snapshot {
            state: est.state.clone(),
            l_cov: est.l_cov,
            internal: est.snapshot_internal(),
        });
        est.last_snapshot_t = est.state.t;
        est
    }

    pub fn health(&self) -> EstimatorHealth {
        use std::mem::size_of;
        let bytes = self.history.snapshots.capacity() * size_of::<Snapshot>()
            + self.history.props.capacity() * size_of::<Option<PropRecord>>()
            + self.history.events.capacity() * size_of::<Event>();
        EstimatorHealth {
            t: self.state.t,
            last_gnss_nis: self.last_gnss_nis,
            last_gnss_accepted: self.last_gnss_accepted,
            gnss_reject_streak: self.gnss_reject_streak,
            history_snapshots: self.history.snapshots.len(),
            history_events: self.history.events.len(),
            history_props: self.history.props.len(),
            history_capacity_bytes: bytes,
        }
    }

    #[inline]
    pub(super) fn accumulate_segment(&mut self, phi: &MatN, qk: &MatN) {
        let phi_t = phi.transpose();
        self.seg_q = (*phi) * self.seg_q * phi_t + (*qk);
        self.seg_phi = (*phi) * self.seg_phi;
        self.seg_has_data = true;
    }

    #[inline]
    pub(super) fn maybe_checkpoint(&mut self, force: bool) {
        self.maybe_checkpoint_with_prediction(force, None);
    }

    #[inline]
    pub(super) fn maybe_checkpoint_with_prediction(
        &mut self,
        force: bool,
        predicted_end_state: Option<crate::types::State>,
    ) {
        if !force {
            let dt = self.state.t - self.last_snapshot_t;
            if dt < self.history_snapshot_dt - 1e-9 {
                return;
            }
        }
        self.checkpoint(predicted_end_state);
    }

    pub(super) fn checkpoint(&mut self, predicted_end_state: Option<crate::types::State>) {
        let t = self.state.t;
        let internal = self.snapshot_internal();
        // Both the history snapshot and the propagation record need owned copies
        // of the state because the buffers outlive this call.
        let state = self.state.clone();
        let l_cov = self.l_cov;

        if let Some(last) = self.history.snapshots.back_mut() {
            if (t - last.state.t).abs() < 1e-9 {
                last.state = state;
                last.l_cov = l_cov;
                last.internal = internal;
                self.last_snapshot_t = t;
                self.seg_phi = MatN::identity();
                self.seg_q = MatN::zeros();
                self.seg_has_data = false;
                return;
            }
        }

        if self.seg_has_data {
            let predicted_end_state = match predicted_end_state {
                Some(state_before_update) => state_before_update,
                None => {
                    // The same state must be stored both as the snapshot payload
                    // and as the predicted end state for RTS.
                    state.clone()
                }
            };
            self.history.push_prop(PropRecord {
                phi: self.seg_phi,
                q: self.seg_q,
                predicted_end_state,
            });
            self.seg_phi = MatN::identity();
            self.seg_q = MatN::zeros();
            self.seg_has_data = false;
        }

        self.history.push_snapshot(Snapshot {
            state,
            l_cov,
            internal,
        });
        self.last_snapshot_t = t;
    }

    pub fn cov(&self) -> MatN {
        let p = self.l_cov * self.l_cov.transpose();
        #[cfg(feature = "diagnostics")]
        super::diag::maybe_check_covariance(&p, "cov()");
        p
    }

    #[cfg(feature = "diagnostics")]
    pub fn cov_strict(&self) -> MatN {
        let p = self.l_cov * self.l_cov.transpose();
        let l = super::math::chol_lower_with_enforce(p, "cov_strict()");
        l * l.transpose()
    }

    pub fn set_aero_table(&mut self, table: AeroTable) {
        self.aero_table = Some(Arc::new(table));
    }

    pub fn set_aero_table_shared(&mut self, table: Arc<AeroTable>) {
        self.aero_table = Some(table);
    }

    pub fn gravity(&self, r: &Vector3<f64>) -> Vector3<f64> {
        let mu = self.dyn_params.mu_earth;
        let norm = r.norm();
        if norm < 1e-3 {
            return Vector3::zeros();
        }
        -mu * r / norm.powi(3)
    }

    pub(super) fn gravity_gradient(&self, r: &Vector3<f64>) -> Matrix3<f64> {
        let mu = self.dyn_params.mu_earth;
        let r_norm = r.norm();
        if r_norm < 1e-3 {
            return Matrix3::zeros();
        }
        let rr_t = r * r.transpose();
        let r2 = r_norm * r_norm;
        let r5 = r2 * r2 * r_norm;
        (3.0 * rr_t / r5 - Matrix3::identity() / r_norm.powi(3)) * mu
    }

    pub fn set_blackout_residual(
        &mut self,
        start: f64,
        end: f64,
        aero_force_rw: f64,
        aero_torque_rw: Option<f64>,
    ) {
        self.blackout_window = Some((start, end));
        self.aero_force_rw_blackout = Some(aero_force_rw);
        self.aero_torque_rw_blackout = aero_torque_rw;
    }

    pub fn set_robust_tuning(&mut self, gate_p: f64, r_scale_cap: f64, student_nu: f64) {
        self.robust_gate_alpha = Self::sanitize_gate_alpha(gate_p);
        self.robust_r_scale_cap = if r_scale_cap.is_finite() {
            r_scale_cap.max(super::robust::DEFAULT_R_SCALE_CAP.min(1.0))
        } else {
            super::robust::DEFAULT_R_SCALE_CAP
        };
        self.robust_student_nu = if student_nu.is_finite() {
            student_nu.max(0.0)
        } else {
            0.0
        };
        self.reset_iae_scales();
    }

    pub fn enable_iae(&mut self, gamma: f64) {
        let gamma = if gamma.is_finite() {
            gamma.clamp(0.01, 1.0)
        } else {
            0.1
        };
        self.iae_gamma = Some(gamma);
        self.reset_iae_scales();
    }

    pub fn disable_iae(&mut self) {
        self.iae_gamma = None;
        self.reset_iae_scales();
    }

    pub fn inflate_rv(&mut self, pos_scale: f64, vel_scale: f64) {
        let mut p = self.cov();
        for i in 0..3 {
            p[(i, i)] *= pos_scale;
            p[(i + 3, i + 3)] *= vel_scale;
        }
        enforce_pd(&mut p);
        self.l_cov = super::math::chol_lower_with_enforce(p, "inflate_rv");
    }

    pub fn inflate_rv_record(&mut self, pos_scale: f64, vel_scale: f64) {
        let t = self.state.t;
        self.inflate_rv(pos_scale, vel_scale);
        self.maybe_checkpoint(true);
        self.history.push_event(Event::InflateRv {
            t,
            pos_scale,
            vel_scale,
        });
    }

    pub fn enable_aero_learning(&mut self, learning_rate: f64) {
        self.aero_learning_enabled = true;
        self.aero_learning_rate = learning_rate.clamp(0.001, 0.1);
    }

    pub fn disable_aero_learning(&mut self) {
        self.aero_learning_enabled = false;
    }

    pub(super) fn update_aero_learning(&mut self, v_meas: Vector3<f64>, t: f64) {
        if !self.aero_learning_enabled {
            return;
        }

        if let Some((t_prev, v_prev)) = self.last_v_meas {
            let dt = t - t_prev;
            if dt > 0.01 && dt < 2.0 {
                let v_change_meas = v_meas - v_prev;
                let v_change_pred = self.state.v - v_prev;
                let residual_accel = (v_change_meas - v_change_pred) / dt;

                let alpha = self.aero_learning_rate;
                self.learned_aero_correction =
                    (1.0 - alpha) * self.learned_aero_correction + alpha * residual_accel;

                let max_correction = self.noise.aero_learning_max_correction;
                let norm = self.learned_aero_correction.norm();
                if norm > max_correction {
                    self.learned_aero_correction *= max_correction / norm;
                }
            }
        }

        self.last_v_meas = Some((t, v_meas));
    }
}
