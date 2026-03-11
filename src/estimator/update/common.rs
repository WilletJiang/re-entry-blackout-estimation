use super::{math::chol_lower_with_enforce, Estimator, MatN, MeasMat, MeasVec};
use crate::util::enforce_pd;
use nalgebra::{linalg::Cholesky, SMatrix};

impl Estimator {
    pub(super) fn measurement_update<const M: usize>(
        &mut self,
        sensor_id: usize,
        innov: MeasVec<M>,
        h: MeasMat<M>,
        r: SMatrix<f64, M, M>,
        dof: usize,
    ) -> bool {
        let p = self.cov();
        let base_s = h.clone() * p.clone() * h.transpose() + r;
        let base_s_chol = match Cholesky::new(base_s.clone()) {
            Some(c) => c,
            None => return false,
        };
        let base_s_inv = base_s_chol.inverse();
        let base_d2 = (innov.transpose() * base_s_inv * innov)[(0, 0)];
        let r = self.robust_scale(sensor_id, dof, base_d2, r);

        let s = h.clone() * p.clone() * h.transpose() + r;
        let s_chol = match Cholesky::new(s.clone()) {
            Some(c) => c,
            None => return false,
        };
        let s_inv = s_chol.inverse();
        let d2 = (innov.transpose() * s_inv.clone() * innov)[(0, 0)];
        let gate_sq = self.robust_gate_sq(dof, self.noise.nis_alpha);
        if d2.is_nan() || d2 > gate_sq {
            return false;
        }
        self.update_iae_scale(sensor_id, dof, d2);
        let k = p * h.transpose() * s_inv;
        let dx = k * innov;

        let mut dx_no_bias = dx.clone();
        for i in 12..18 {
            dx_no_bias[i] = 0.0;
        }
        self.apply_state_error(dx_no_bias);
        let i_kh = MatN::identity() - k * h;
        let mut p_new = i_kh * p * i_kh.transpose() + k * r * k.transpose();
        enforce_pd(&mut p_new);
        self.l_cov = chol_lower_with_enforce(p_new, "measurement_update");
        true
    }
}
