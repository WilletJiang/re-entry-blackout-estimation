use super::{
    gating::chi_square_gate, history::Event, math::chol_lower_with_enforce,
    robust::ROBUST_SENSOR_GPS_RAW, Estimator, GpsRawMeas, GpsSatMeas, MatN,
};
use crate::types::STATE_DIM;
use crate::util::enforce_pd;
use nalgebra::{linalg::Cholesky, DMatrix, DVector, Matrix3, Matrix6, SymmetricEigen, Vector3};

const GPS_RAW_MAX_ROWS: usize = 64;
const ROW_PR: u8 = 0;
const ROW_RR: u8 = 1;

struct WlsRv {
    pos: Vector3<f64>,
    vel: Vector3<f64>,
    r_cov: Matrix6<f64>,
    nis: f64,
    dof: usize,
}

struct GpsRawModel {
    rows: usize,
    row_kind: [u8; GPS_RAW_MAX_ROWS],
    h: DMatrix<f64>,
    innov: DVector<f64>,
    r_diag: DVector<f64>,
    info_v: Matrix3<f64>,
}

impl GpsRawModel {
    fn build(r_user: Vector3<f64>, v_user: Vector3<f64>, sats: &[GpsSatMeas]) -> Option<Self> {
        let rows = count_rows(sats)?;
        let mut row_kind = [0u8; GPS_RAW_MAX_ROWS];
        let mut h = DMatrix::<f64>::zeros(rows, STATE_DIM);
        let mut innov = DVector::<f64>::zeros(rows);
        let mut r_diag = DVector::<f64>::zeros(rows);
        let mut info_v = Matrix3::<f64>::zeros();

        let mut row = 0usize;
        for sat in sats {
            let dr = sat.pos_ecef - r_user;
            let rho = dr.norm().max(1e-3);
            let los = dr / rho;

            for j in 0..3 {
                h[(row, j)] = -los[j];
            }
            let pr_pred = rho;
            row_kind[row] = ROW_PR;
            innov[row] = sat.pr_m - pr_pred;
            r_diag[row] = sat.pr_sigma * sat.pr_sigma;
            row += 1;

            if let (Some(rr_meas), Some(rr_sig)) = (sat.range_rate_mps, sat.rr_sigma) {
                let rel_v = sat.vel_ecef - v_user;
                let rr_pred = rel_v.dot(&los);
                let h_r = -(rel_v - los * rr_pred) / rho;
                for j in 0..3 {
                    h[(row, j)] = h_r[j];
                }
                for j in 0..3 {
                    h[(row, j + 3)] = -los[j];
                }
                if rr_sig.is_finite() && rr_sig > 0.0 {
                    let w = 1.0 / (rr_sig * rr_sig);
                    info_v += w * (los * los.transpose());
                }
                row_kind[row] = ROW_RR;
                innov[row] = rr_meas - rr_pred;
                r_diag[row] = rr_sig * rr_sig;
                row += 1;
            }
        }

        Some(Self {
            rows,
            row_kind,
            h,
            innov,
            r_diag,
            info_v,
        })
    }

    fn robustify_pseudorange(&mut self, p_rv: &Matrix6<f64>) {
        const ROBUST_K: f64 = 4.0;
        const ROBUST_CAP: f64 = 100.0;

        for i in 0..self.rows {
            if self.row_kind[i] != ROW_PR {
                continue;
            }
            let mut s_ii = self.r_diag[i];
            if !(s_ii.is_finite() && s_ii > 0.0) {
                continue;
            }
            for a in 0..6 {
                for b in 0..6 {
                    s_ii += self.h[(i, a)] * p_rv[(a, b)] * self.h[(i, b)];
                }
            }
            if !(s_ii.is_finite() && s_ii > 0.0) {
                continue;
            }
            let u = self.innov[i] / s_ii.sqrt();
            if !u.is_finite() {
                continue;
            }
            let abs = u.abs();
            if abs <= ROBUST_K {
                continue;
            }
            let scale = (abs / ROBUST_K).min(ROBUST_CAP);
            self.r_diag[i] *= scale * scale;
        }
    }

    fn r_matrix(&self) -> DMatrix<f64> {
        let mut r = DMatrix::<f64>::zeros(self.rows, self.rows);
        for i in 0..self.rows {
            r[(i, i)] = self.r_diag[i];
        }
        r
    }
}

struct GpsRawLin {
    p_dyn: DMatrix<f64>,
    r: DMatrix<f64>,
    d2: f64,
    k: DMatrix<f64>,
    dx_raw: DVector<f64>,
}

impl GpsRawLin {
    fn compute(p: &MatN, model: &GpsRawModel) -> Option<Self> {
        let p_dyn = DMatrix::<f64>::from_column_slice(STATE_DIM, STATE_DIM, p.as_slice());
        let r = model.r_matrix();
        let s = &model.h * &p_dyn * model.h.transpose() + &r;
        let s_chol = nalgebra::linalg::Cholesky::new(s.clone())?;
        let s_inv = s_chol.inverse();
        let d2 = (&model.innov.transpose() * &s_inv * &model.innov)[(0, 0)];
        let k = &p_dyn * model.h.transpose() * s_inv;
        let dx_raw = &k * &model.innov;
        Some(Self {
            p_dyn,
            r,
            d2,
            k,
            dx_raw,
        })
    }
}

impl Estimator {
    pub fn update_gps_raw(&mut self, meas: GpsRawMeas) {
        self.process_gps_raw(meas, true, true);
    }

    pub fn update_gps_raw_ref(&mut self, t: f64, sats: &[GpsSatMeas]) {
        self.process_gps_raw(
            GpsRawMeas {
                t,
                sats: sats.to_vec(),
            },
            true,
            true,
        );
    }

    pub(super) fn process_gps_raw(&mut self, meas: GpsRawMeas, record: bool, allow_rewind: bool) {
        let replay_events = if allow_rewind {
            self.rewind_if_needed(meas.t)
        } else {
            None
        };

        {
            let t = meas.t;
            let sats = meas.sats.as_slice();

            self.last_gnss_accepted = Some((t, false));
            if sats.len() < 4 {
                self.replay_tail(replay_events);
                return;
            }

            self.maybe_bias_cal_observe_gps_raw(t, sats);

            let p = self.cov();
            self.last_gnss_prior_t = Some(t);
            self.last_gnss_prior_rv_p = Some(p.fixed_view::<6, 6>(0, 0).into_owned());
            self.last_gnss_prior_r = Some(self.state.r);
            self.last_gnss_prior_v = Some(self.state.v);

            let Some(mut model) = GpsRawModel::build(self.state.r, self.state.v, sats) else {
                self.replay_tail(replay_events);
                return;
            };
            model.robustify_pseudorange(&p.fixed_view::<6, 6>(0, 0).into_owned());
            let Some(base_lin) = GpsRawLin::compute(&p, &model) else {
                self.replay_tail(replay_events);
                return;
            };
            let raw_scale =
                self.robust_covariance_scale(ROBUST_SENSOR_GPS_RAW, model.rows, base_lin.d2);
            if raw_scale > 1.0 {
                model.r_diag *= raw_scale;
            }
            let Some(mut lin) = GpsRawLin::compute(&p, &model) else {
                self.replay_tail(replay_events);
                return;
            };

            self.last_gnss_nis = Some((t, lin.d2, model.rows));
            let gate_sq =
                if self.noise.gps_raw_nis_gate.is_finite() && self.noise.gps_raw_nis_gate > 0.0 {
                    self.noise.gps_raw_nis_gate
                } else {
                    self.robust_gate_sq(model.rows, self.noise.nis_alpha)
                };
            if !lin.d2.is_finite() || lin.d2 > gate_sq {
                self.replay_tail(replay_events);
                return;
            }
            self.last_gnss_accepted = Some((t, true));
            self.update_iae_scale(ROBUST_SENSOR_GPS_RAW, model.rows, lin.d2);

            self.apply_gps_raw_active_set(t, model.rows, &lin.dx_raw, &mut lin.k);
            suppress_pseudorange_velocity(&model.row_kind, model.rows, &mut lin.k);
            maybe_project_velocity_weak_direction(&mut lin.k, model.info_v, model.rows);

            let dx_dyn = &lin.k * &model.innov;
            self.apply_state_error(dvec_to_vecn(&dx_dyn));

            let i_kh = DMatrix::<f64>::identity(STATE_DIM, STATE_DIM) - &lin.k * &model.h;
            let p_new_dyn =
                &i_kh * &lin.p_dyn * i_kh.transpose() + &lin.k * &lin.r * lin.k.transpose();
            let mut p_new = MatN::from_column_slice(p_new_dyn.as_slice());
            enforce_pd(&mut p_new);
            self.l_cov = chol_lower_with_enforce(p_new, "update_gps_raw");

            self.update_aero_learning(self.state.v, t);
            self.state.t = t;
            self.maybe_checkpoint(true);
        }

        if record {
            self.history.push_event(Event::GpsRaw(meas));
        }
        self.replay_tail(replay_events);
    }

    fn maybe_bias_cal_observe_gps_raw(&mut self, t: f64, sats: &[GpsSatMeas]) {
        if !self.bias_cal.enabled {
            return;
        }
        let Some(sol) = gps_raw_wls_rv(self.state.r, self.state.v, sats) else {
            return;
        };
        let gate_sq = chi_square_gate(sol.dof, 0.01);
        if !(sol.nis.is_finite() && sol.nis <= gate_sq) {
            return;
        }
        self.bias_cal.from_pvt = false;
        self.bias_cal_observe_rv(t, sol.pos, sol.vel, &sol.r_cov);
    }

    fn apply_gps_raw_active_set(
        &self,
        t: f64,
        rows: usize,
        dx_raw: &DVector<f64>,
        k: &mut DMatrix<f64>,
    ) {
        let (mut bias_updates_enabled, max_ba_update) = self.bias_update_profile(t);
        if self.bias_cal.enabled {
            bias_updates_enabled = false;
        }
        let bias_gain = self.noise.gps_ba_update_gain.clamp(0.0, 1.0);

        for i in 0..STATE_DIM {
            let active_rv = i < 6;
            let active_ba = bias_updates_enabled && (15..18).contains(&i);
            let scale = if active_rv {
                1.0
            } else if active_ba {
                let delta = bias_gain * dx_raw[i];
                let abs = delta.abs();
                if abs > max_ba_update && abs > 0.0 {
                    bias_gain * (max_ba_update / abs)
                } else {
                    bias_gain
                }
            } else {
                0.0
            };
            if scale == 0.0 {
                for j in 0..rows {
                    k[(i, j)] = 0.0;
                }
                continue;
            }
            if (scale - 1.0).abs() > 0.0 {
                for j in 0..rows {
                    k[(i, j)] *= scale;
                }
            }
        }
    }
}

fn count_rows(sats: &[GpsSatMeas]) -> Option<usize> {
    let mut rows = 0usize;
    for sat in sats {
        rows += 1;
        if sat.range_rate_mps.is_some() {
            rows += 1;
        }
    }
    if rows > GPS_RAW_MAX_ROWS {
        return None;
    }
    Some(rows)
}

fn suppress_pseudorange_velocity(
    row_kind: &[u8; GPS_RAW_MAX_ROWS],
    rows: usize,
    k: &mut DMatrix<f64>,
) {
    for j in 0..rows {
        if row_kind[j] != ROW_PR {
            continue;
        }
        for i in 3..6 {
            k[(i, j)] = 0.0;
        }
    }
}

fn maybe_project_velocity_weak_direction(k: &mut DMatrix<f64>, info_v: Matrix3<f64>, rows: usize) {
    let Some(ch) = Cholesky::new(info_v) else {
        return;
    };
    let cov_v = ch.inverse();
    let eig = SymmetricEigen::new(cov_v);
    let mut weak_idx = 0usize;
    let mut weak_var = eig.eigenvalues[0];
    for i in 1..3 {
        if eig.eigenvalues[i] > weak_var {
            weak_var = eig.eigenvalues[i];
            weak_idx = i;
        }
    }
    let weak_sig = weak_var.max(0.0).sqrt();
    if weak_sig <= 0.25 {
        return;
    }
    let weak = eig.eigenvectors.column(weak_idx).into_owned();
    let proj = Matrix3::<f64>::identity() - weak * weak.transpose();
    for j in 0..rows {
        let kv = Vector3::new(k[(3, j)], k[(4, j)], k[(5, j)]);
        let kv = proj * kv;
        k[(3, j)] = kv.x;
        k[(4, j)] = kv.y;
        k[(5, j)] = kv.z;
    }
}

fn dvec_to_vecn(dx: &DVector<f64>) -> super::VecN {
    let mut out = super::VecN::zeros();
    for i in 0..STATE_DIM {
        out[i] = dx[i];
    }
    out
}

fn gps_raw_wls_rv(r0: Vector3<f64>, v0: Vector3<f64>, sats: &[GpsSatMeas]) -> Option<WlsRv> {
    if sats.len() < 4 {
        return None;
    }

    let (r_est, cov_r) = wls_position(r0, sats)?;
    let (v_est, cov_v) = wls_velocity(v0, r_est, sats)?;
    let (d2, dof) = wls_nis_and_dof(r_est, v_est, sats)?;

    let mut r_cov = Matrix6::<f64>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            r_cov[(i, j)] = cov_r[(i, j)];
            r_cov[(i + 3, j + 3)] = cov_v[(i, j)];
        }
    }
    Some(WlsRv {
        pos: r_est,
        vel: v_est,
        r_cov,
        nis: d2,
        dof,
    })
}

fn wls_position(r0: Vector3<f64>, sats: &[GpsSatMeas]) -> Option<(Vector3<f64>, Matrix3<f64>)> {
    let mut r_est = r0;
    let mut info_r = Matrix3::<f64>::zeros();
    let mut chol_r: Option<Cholesky<f64, nalgebra::U3>> = None;

    for _ in 0..2 {
        info_r.fill(0.0);
        let mut b_r = Vector3::<f64>::zeros();
        let mut used = 0usize;

        for sat in sats {
            let sigma = sat.pr_sigma;
            if !(sigma.is_finite() && sigma > 0.0) {
                continue;
            }

            let dr = sat.pos_ecef - r_est;
            let rho = dr.norm();
            if !(rho.is_finite() && rho > 1e-6) {
                return None;
            }
            let los = dr / rho;
            let res = sat.pr_m - rho;
            if !res.is_finite() {
                return None;
            }

            let w = 1.0 / (sigma * sigma);
            info_r += w * (los * los.transpose());
            b_r += -w * los * res;
            used += 1;
        }

        if used < 4 {
            return None;
        }

        let chol = Cholesky::new(info_r)?;
        let delta = chol.solve(&b_r);
        if !(delta.x.is_finite() && delta.y.is_finite() && delta.z.is_finite()) {
            return None;
        }
        r_est += delta;
        chol_r = Some(chol);
        if delta.norm() < 1e-3 {
            break;
        }
    }

    let chol_r = chol_r?;
    let cov_r = chol_r.inverse();
    Some((r_est, cov_r))
}

fn wls_velocity(
    v0: Vector3<f64>,
    r_est: Vector3<f64>,
    sats: &[GpsSatMeas],
) -> Option<(Vector3<f64>, Matrix3<f64>)> {
    let mut v_est = v0;
    let mut info_v = Matrix3::<f64>::zeros();
    let mut b_v = Vector3::<f64>::zeros();
    let mut used_rr = 0usize;

    for sat in sats {
        let Some(rr_meas) = sat.range_rate_mps else {
            continue;
        };
        let Some(rr_sig) = sat.rr_sigma else {
            continue;
        };
        if !(rr_sig.is_finite() && rr_sig > 0.0) {
            continue;
        }

        let dr = sat.pos_ecef - r_est;
        let rho = dr.norm();
        if !(rho.is_finite() && rho > 1e-6) {
            return None;
        }
        let los = dr / rho;
        let rr_pred = (sat.vel_ecef - v_est).dot(&los);
        let res = rr_meas - rr_pred;
        if !res.is_finite() {
            return None;
        }

        let w = 1.0 / (rr_sig * rr_sig);
        info_v += w * (los * los.transpose());
        b_v += -w * los * res;
        used_rr += 1;
    }

    if used_rr < 4 {
        return None;
    }

    let chol_v = Cholesky::new(info_v)?;
    let delta_v = chol_v.solve(&b_v);
    if !(delta_v.x.is_finite() && delta_v.y.is_finite() && delta_v.z.is_finite()) {
        return None;
    }
    v_est += delta_v;
    let cov_v = chol_v.inverse();
    Some((v_est, cov_v))
}

fn wls_nis_and_dof(
    r_est: Vector3<f64>,
    v_est: Vector3<f64>,
    sats: &[GpsSatMeas],
) -> Option<(f64, usize)> {
    let mut d2 = 0.0f64;
    let mut n_pr = 0usize;
    let mut n_rr = 0usize;
    for sat in sats {
        let sigma = sat.pr_sigma;
        if !(sigma.is_finite() && sigma > 0.0) {
            continue;
        }

        let dr = sat.pos_ecef - r_est;
        let rho = dr.norm();
        if !(rho.is_finite() && rho > 1e-6) {
            continue;
        }
        let res = sat.pr_m - rho;
        if res.is_finite() {
            d2 += (res * res) / (sigma * sigma);
        }
        n_pr += 1;

        let Some(rr_meas) = sat.range_rate_mps else {
            continue;
        };
        let Some(rr_sig) = sat.rr_sigma else {
            continue;
        };
        if !(rr_sig.is_finite() && rr_sig > 0.0) {
            continue;
        }
        let los = dr / rho;
        let rr_pred = (sat.vel_ecef - v_est).dot(&los);
        let res_rr = rr_meas - rr_pred;
        if res_rr.is_finite() {
            d2 += (res_rr * res_rr) / (rr_sig * rr_sig);
        }
        n_rr += 1;
    }

    if n_pr < 4 || n_rr < 4 {
        return None;
    }
    let dof = (n_pr - 3) + (n_rr - 3);
    if dof == 0 {
        return None;
    }
    Some((d2, dof))
}
