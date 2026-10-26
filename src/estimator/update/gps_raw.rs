use super::{
    gating::chi_square_gate, history::Event, math::chol_lower_with_enforce,
    robust::ROBUST_SENSOR_GPS_RAW, Estimator, GpsRawMeas, GpsSatMeas, MatN, MeasMat,
};
use crate::types::STATE_DIM;
use crate::util::enforce_pd;
use nalgebra::{
    linalg::Cholesky, DMatrix, DVector, Matrix3, Matrix6, SVector, SymmetricEigen, Vector3,
};

const GPS_RAW_MAX_ROWS: usize = 64;
const ROW_PR: u8 = 0;
const ROW_RR: u8 = 1;
const RECOVERY_WINDOW_SEC: f64 = 2.0;
const FLOAT_TIME_TOL_SEC: f64 = 1.0e-9;

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
        let predicted_end_state = if self.seg_has_data {
            // RTS needs the propagation prior at this timestamp after the update
            // mutates `self.state`, so the pre-update state must be owned here.
            Some(self.state.clone())
        } else {
            None
        };

        {
            let t = meas.t;
            let sats = meas.sats.as_slice();

            self.last_gnss_accepted = Some((t, false));
            if sats.len() < 4 {
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
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
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
                self.replay_tail(replay_events);
                return;
            };
            let base_d2 = GpsRawLin::compute(&p, &model).map_or(f64::NAN, |lin| lin.d2);
            self.last_gnss_base_nis = Some((t, base_d2, model.rows));
            if let Some(accepted) = self.acquire_raw_gnss(&p, t, sats) {
                self.last_gnss_accepted = Some((t, accepted));
                self.record_gnss_posterior(t);
                self.maybe_checkpoint_with_prediction(true, predicted_end_state);
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
                self.replay_tail(replay_events);
                return;
            }

            // Gate the original innovation before any residual-dependent covariance
            // inflation. Otherwise an extreme Cauchy epoch can certify itself.
            if !base_d2.is_finite()
                || base_d2 > self.robust_gate_sq(model.rows, self.noise.nis_alpha)
            {
                self.gnss_reject_streak = self.gnss_reject_streak.saturating_add(1);
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
                self.replay_tail(replay_events);
                return;
            }
            model.robustify_pseudorange(&p.fixed_view::<6, 6>(0, 0).into_owned());
            let Some(base_lin) = GpsRawLin::compute(&p, &model) else {
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
                self.replay_tail(replay_events);
                return;
            };
            let raw_scale =
                self.robust_covariance_scale(ROBUST_SENSOR_GPS_RAW, model.rows, base_lin.d2);
            if raw_scale > 1.0 {
                model.r_diag *= raw_scale;
            }
            let Some(mut lin) = GpsRawLin::compute(&p, &model) else {
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
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
                if record {
                    self.history.push_event(Event::GpsRaw(meas));
                }
                self.replay_tail(replay_events);
                return;
            }
            self.last_gnss_accepted = Some((t, true));
            self.gnss_reject_streak = 0;
            self.update_iae_scale(ROBUST_SENSOR_GPS_RAW, model.rows, lin.d2);

            self.apply_gps_raw_active_set(t, model.rows, &lin.dx_raw, &mut lin.k);
            if self.bias_cal.enabled {
                suppress_pseudorange_velocity(&model.row_kind, model.rows, &mut lin.k);
            }
            if self.bias_cal.enabled && self.should_project_weak_velocity_direction(t) {
                maybe_project_velocity_weak_direction(&mut lin.k, model.info_v, model.rows);
            }

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
            self.record_gnss_posterior(t);
            self.maybe_checkpoint_with_prediction(true, predicted_end_state);
        }

        if record {
            self.history.push_event(Event::GpsRaw(meas));
        }
        self.replay_tail(replay_events);
    }

    fn acquire_raw_gnss(&mut self, p: &MatN, t: f64, sats: &[GpsSatMeas]) -> Option<bool> {
        if self.bias_cal.enabled
            || (self.gnss_reject_streak < 8
                && p.fixed_view::<3, 3>(0, 0).trace() < 3.0 * 100.0_f64.powi(2)
                && p.fixed_view::<3, 3>(3, 3).trace() < 3.0)
        {
            return None;
        }
        // The Doppler/range model is nonlinear for a kilometre-scale prior.
        // Solve redundant observations first, checking their post-fit residual,
        // then assimilate that solution once with its joint RV covariance.
        // Rank-deficient constellations retain the ordinary raw update path.
        let solution = gps_raw_wls_rv(self.state.r, self.state.v, sats)?;
        if std::env::var_os("TRACE_RAW").is_some() {
            eprintln!(
                "raw acquisition t={t:.3} fit_nis={:.3} fit_dof={} pos_step={:.3} vel_step={:.3}",
                solution.nis,
                solution.dof,
                (solution.pos - self.state.r).norm(),
                (solution.vel - self.state.v).norm()
            );
        }
        if !solution.nis.is_finite()
            || solution.nis > chi_square_gate(solution.dof, self.noise.nis_alpha)
        {
            self.gnss_reject_streak = 0;
            return Some(false);
        }
        let mut innovation = SVector::<f64, 6>::zeros();
        innovation
            .fixed_rows_mut::<3>(0)
            .copy_from(&(solution.pos - self.state.r));
        innovation
            .fixed_rows_mut::<3>(3)
            .copy_from(&(solution.vel - self.state.v));
        let mut h = MeasMat::<6>::zeros();
        for i in 0..6 {
            h[(i, i)] = 1.0;
        }
        let s = p.fixed_view::<6, 6>(0, 0).into_owned() + solution.r_cov;
        let d2 = s
            .cholesky()
            .map_or(f64::NAN, |c| innovation.dot(&c.solve(&innovation)));
        self.last_gnss_nis = Some((t, d2, 6));
        if std::env::var_os("TRACE_RAW").is_some() {
            eprintln!(
                "raw acquisition t={t:.3} prior_nis={d2:.3} prior_sigma_r={:.3}",
                p.fixed_view::<3, 3>(0, 0).trace().sqrt()
            );
        }
        let accepted =
            self.measurement_update(ROBUST_SENSOR_GPS_RAW, innovation, h, solution.r_cov, 6);
        if accepted {
            self.gnss_reject_streak = 0;
        } else {
            // A Gaussian initial prior has tails outside any finite gate. Do
            // not let that one rare draw reject every subsequent good epoch.
            // Internally inconsistent WLS fits above break this streak. Widen
            // only after repeated independently validated solutions; the next
            // epoch must still pass acquisition and its ordinary prior gate.
            self.gnss_reject_streak = self.gnss_reject_streak.saturating_add(1);
            if self.gnss_reject_streak >= 8 {
                self.inflate_rv(1000.0, 100.0);
                self.gnss_reject_streak = 0;
            }
        }
        Some(accepted)
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

    fn should_project_weak_velocity_direction(&self, t_meas: f64) -> bool {
        let Some((_, blackout_end)) = self.blackout_window else {
            return false;
        };
        t_meas >= blackout_end - FLOAT_TIME_TOL_SEC
            && t_meas <= blackout_end + RECOVERY_WINDOW_SEC + FLOAT_TIME_TOL_SEC
    }

    fn apply_gps_raw_active_set(
        &self,
        t: f64,
        rows: usize,
        dx_raw: &DVector<f64>,
        k: &mut DMatrix<f64>,
    ) {
        if !self.bias_cal.enabled && self.noise.gps_ba_update_gain == 1.0 {
            return;
        }
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

    let (r_est, _) = wls_position(r0, sats)?;
    let (v_est, _) = wls_velocity(v0, r_est, sats)?;
    let (d2, dof) = wls_nis_and_dof(r_est, v_est, sats)?;
    let model = GpsRawModel::build(r_est, v_est, sats)?;
    let mut information = Matrix6::zeros();
    for row in 0..model.rows {
        let h = SVector::<f64, 6>::from_fn(|i, _| model.h[(row, i)]);
        information += h * h.transpose() / model.r_diag[row];
    }
    let r_cov = information.cholesky()?.inverse();
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

    for _ in 0..5 {
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
