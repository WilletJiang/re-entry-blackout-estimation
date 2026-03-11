use super::{math::chol_lower_with_enforce, Estimator, MatN, VecN};
use crate::util::enforce_pd;
use nalgebra::UnitQuaternion;

impl Estimator {
    pub fn rts_smooth(&self, lag_sec: f64) -> Option<(crate::types::State, MatN)> {
        let n_snap = self.history.snapshots.len();
        if n_snap == 0 {
            return None;
        }
        let newest = self.history.snapshots.back().unwrap();
        let target_t = newest.state.t - lag_sec;
        let target_idx = self.history.find_snapshot_index(target_t)?;

        let mut p_f: Vec<MatN> = Vec::with_capacity(n_snap);
        for snap in &self.history.snapshots {
            let p = snap.l_cov * snap.l_cov.transpose();
            p_f.push(p);
        }

        let mut x_s: Vec<crate::types::State> = self
            .history
            .snapshots
            .iter()
            .map(|s| s.state.clone())
            .collect();
        let mut p_s: Vec<MatN> = p_f.clone();

        if n_snap == 1 {
            return Some((x_s[0].clone(), p_s[0]));
        }

        let newest_idx = n_snap - 1;
        let eye = MatN::identity();
        let zero = MatN::zeros();

        for k in (0..newest_idx).rev() {
            let (phi_k, q_k) = if k < self.history.props.len() {
                let pr = &self.history.props[k];
                (&pr.phi, &pr.q)
            } else {
                (&eye, &zero)
            };

            let mut p_pred = *phi_k * p_f[k] * phi_k.transpose() + *q_k;
            enforce_pd(&mut p_pred);
            let p_pred_inv = match p_pred.try_inverse() {
                Some(inv) => inv,
                None => {
                    let l_pred = chol_lower_with_enforce(p_pred, "rts_smooth P_pred");
                    let l_inv = l_pred
                        .try_inverse()
                        .expect("Cholesky factor must be invertible for SPD P_pred");
                    let l_inv_t = l_inv.transpose();
                    l_inv_t * l_inv
                }
            };
            let a_k = p_f[k] * phi_k.transpose() * p_pred_inv;

            let dx_next = state_error(&x_s[k + 1], &self.history.snapshots[k + 1].state);
            let dx = a_k * dx_next;
            x_s[k] = vec_to_state(&self.history.snapshots[k].state, dx);

            let mut p_new = p_f[k] + a_k * (p_s[k + 1] - p_pred) * a_k.transpose();
            enforce_pd(&mut p_new);
            p_s[k] = p_new;
        }

        Some((x_s[target_idx].clone(), p_s[target_idx]))
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
