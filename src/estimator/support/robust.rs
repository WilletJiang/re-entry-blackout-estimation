use super::{gating::chi_square_gate, Estimator};
use nalgebra::SMatrix;

pub(super) const ROBUST_SENSOR_GNSS_PVT: usize = 0;
pub(super) const ROBUST_SENSOR_GPS_RAW: usize = 1;
pub(super) const ROBUST_SENSOR_STAR: usize = 2;
pub(super) const ROBUST_SENSOR_RADAR: usize = 3;
pub(super) const ROBUST_SENSOR_COUNT: usize = 4;

pub(super) const DEFAULT_R_SCALE_CAP: f64 = 1.0e3;

impl Estimator {
    pub(super) fn robust_gate_sq(&self, dof: usize, fallback_alpha: f64) -> f64 {
        let alpha = self.robust_gate_alpha.unwrap_or(fallback_alpha);
        chi_square_gate(dof, alpha)
    }

    pub(super) fn robust_gate_sq_or(&self, dof: usize, fallback_sq: f64) -> f64 {
        let Some(alpha) = self.robust_gate_alpha else {
            return fallback_sq;
        };
        chi_square_gate(dof, alpha)
    }

    pub(super) fn robust_scale<const M: usize>(
        &self,
        sensor_id: usize,
        dof: usize,
        d2: f64,
        r_cov: SMatrix<f64, M, M>,
    ) -> SMatrix<f64, M, M> {
        let scale = self.robust_covariance_scale(sensor_id, dof, d2);
        r_cov * scale
    }

    pub(super) fn robust_covariance_scale(&self, sensor_id: usize, dof: usize, d2: f64) -> f64 {
        let sensor_scale = self.iae_scales.get(sensor_id).copied().unwrap_or(1.0);
        let student_scale = self.student_t_scale(dof, d2);
        let scale = sensor_scale * student_scale;
        let cap = self.robust_r_scale_cap.max(DEFAULT_R_SCALE_CAP);
        scale.clamp(1.0, cap)
    }

    pub(super) fn update_iae_scale(&mut self, sensor_id: usize, dof: usize, d2: f64) {
        let Some(gamma) = self.iae_gamma else {
            return;
        };
        if dof == 0 || !d2.is_finite() {
            return;
        }

        let cap = self.robust_r_scale_cap.max(DEFAULT_R_SCALE_CAP);
        let innovation_ratio = (d2 / (dof as f64)).max(1.0).min(cap);
        let current_scale = self.iae_scales.get(sensor_id).copied().unwrap_or(1.0);
        let blended_scale = (1.0 - gamma) * current_scale + gamma * innovation_ratio;
        if let Some(scale) = self.iae_scales.get_mut(sensor_id) {
            *scale = blended_scale.clamp(1.0, cap);
        }
    }

    pub(super) fn reset_iae_scales(&mut self) {
        self.iae_scales.fill(1.0);
    }

    pub(super) fn sanitize_gate_alpha(gate_p: f64) -> Option<f64> {
        if !gate_p.is_finite() {
            return None;
        }
        let clamped_gate_p = gate_p.clamp(0.5, 0.999_999);
        Some((1.0 - clamped_gate_p).max(1.0e-6))
    }

    fn student_t_scale(&self, dof: usize, d2: f64) -> f64 {
        if self.robust_student_nu <= 0.0 || dof == 0 || !d2.is_finite() {
            return 1.0;
        }

        let nu = self.robust_student_nu.max(1.0e-6);
        let scale = (nu + d2) / (nu + dof as f64);
        scale.max(1.0)
    }
}
