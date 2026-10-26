use super::{Estimator, MatN, VecN};
use crate::util::enforce_pd;
use nalgebra::{linalg::Cholesky, UnitQuaternion};

impl Estimator {
    /// Smooth the estimator history over a fixed lag.
    ///
    /// `lag_sec` is measured backward from the newest snapshot time. The
    /// function returns `None` when the history does not cover the requested lag
    /// or when the backward pass encounters a numerically unrecoverable step.
    pub fn rts_smooth(&self, lag_sec: f64) -> Option<(crate::types::State, MatN)> {
        let mut oldest = None;
        self.visit_smoothed(lag_sec, |state, cov| oldest = Some((state.clone(), *cov)))?;
        oldest
    }

    /// Visit smoothed checkpoints newest first, using only data already received.
    /// No history-sized workspace is allocated. A failed pass returns `None`;
    /// callers must then discard any samples already delivered to the callback.
    pub fn visit_smoothed(
        &self,
        lag_sec: f64,
        mut visit: impl FnMut(&crate::types::State, &MatN),
    ) -> Option<()> {
        if !lag_sec.is_finite() || lag_sec < 0.0 {
            return None;
        }
        let n_snap = self.history.snapshots.len();
        if n_snap == 0 {
            return None;
        }
        let newest = self.history.snapshots.back()?;
        let target_t = newest.state.t - lag_sec;
        let target_idx = self.history.find_snapshot_index(target_t)?;

        let mut smoothed_state = newest.state.clone();
        let mut smoothed_cov = newest.l_cov * newest.l_cov.transpose();

        visit(&smoothed_state, &smoothed_cov);
        let newest_idx = n_snap - 1;
        let eye = MatN::identity();
        let zero = MatN::zeros();

        for k in (target_idx..newest_idx).rev() {
            let snapshot = &self.history.snapshots[k];
            let filtered_cov = snapshot.l_cov * snapshot.l_cov.transpose();
            let prop_k = if k < self.history.props.len() {
                self.history.props[k].as_ref()
            } else {
                None
            };
            let (phi_k, q_k) = if let Some(prop_record) = prop_k {
                (&prop_record.phi, &prop_record.q)
            } else {
                (&eye, &zero)
            };

            let mut p_pred = *phi_k * filtered_cov * phi_k.transpose() + *q_k;
            enforce_pd(&mut p_pred);
            // Solve P_pred A^T = Phi P_filtered; avoid explicitly inverting
            // a covariance whose translation and bias scales differ greatly.
            let factor = Cholesky::new(p_pred)?;
            let a_k = factor.solve(&(*phi_k * filtered_cov)).transpose();

            let predicted_next_state = if let Some(prop_record) = prop_k {
                &prop_record.predicted_end_state
            } else {
                &self.history.snapshots[k + 1].state
            };
            let dx_next = state_error(&smoothed_state, predicted_next_state);
            let dx = a_k * dx_next;
            smoothed_state = vec_to_state(&snapshot.state, dx);

            let mut p_new = filtered_cov + a_k * (smoothed_cov - p_pred) * a_k.transpose();
            enforce_pd(&mut p_new);
            smoothed_cov = p_new;
            visit(&smoothed_state, &smoothed_cov);
        }

        Some(())
    }
}

fn state_error(target: &crate::types::State, base: &crate::types::State) -> VecN {
    let mut v = VecN::zeros();
    v.fixed_rows_mut::<3>(0).copy_from(&(target.r - base.r));
    v.fixed_rows_mut::<3>(3).copy_from(&(target.v - base.v));
    let dq = target.q * base.q.inverse();
    v.fixed_rows_mut::<3>(6).copy_from(&dq.scaled_axis());
    v.fixed_rows_mut::<3>(9)
        .copy_from(&(target.omega - base.omega));
    v.fixed_rows_mut::<3>(12).copy_from(&(target.bg - base.bg));
    v.fixed_rows_mut::<3>(15).copy_from(&(target.ba - base.ba));
    v.fixed_rows_mut::<3>(18)
        .copy_from(&(target.aero_res_force - base.aero_res_force));
    v.fixed_rows_mut::<3>(21)
        .copy_from(&(target.aero_res_torque - base.aero_res_torque));
    v
}

fn vec_to_state(base: &crate::types::State, dx: VecN) -> crate::types::State {
    let mut s = base.clone();
    s.r += dx.fixed_rows::<3>(0).into_owned();
    s.v += dx.fixed_rows::<3>(3).into_owned();
    let dq = UnitQuaternion::from_scaled_axis(dx.fixed_rows::<3>(6).into_owned());
    s.q = UnitQuaternion::new_normalize((dq * s.q).into_inner());
    s.omega += dx.fixed_rows::<3>(9).into_owned();
    s.bg += dx.fixed_rows::<3>(12).into_owned();
    s.ba += dx.fixed_rows::<3>(15).into_owned();
    s.aero_res_force += dx.fixed_rows::<3>(18).into_owned();
    s.aero_res_torque += dx.fixed_rows::<3>(21).into_owned();
    s
}
