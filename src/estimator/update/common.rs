use super::{gating::chi_square_gate, math::chol_lower_with_enforce, Estimator, MatN, MeasMat, MeasVec};
use crate::util::enforce_pd;
use nalgebra::{linalg::Cholesky, SMatrix};

impl Estimator {
    pub(super) fn measurement_update<const M: usize>(
        &mut self,
        innov: MeasVec<M>,
        h: MeasMat<M>,
        r: SMatrix<f64, M, M>,
        dof: usize,
    ) {
        let p = self.cov();
        let s = h.clone() * p.clone() * h.transpose() + r;
        let s_chol = match Cholesky::new(s.clone()) {
            Some(c) => c,
            None => return,
        };
        let s_inv = s_chol.inverse();
        let d2 = (innov.transpose() * s_inv.clone() * innov)[(0, 0)];
        let gate_sq = chi_square_gate(dof, self.noise.nis_alpha);
        if d2.is_nan() || d2 > gate_sq {
            return;
        }
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
    }
}
