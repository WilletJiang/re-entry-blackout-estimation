use super::{history::Event, history::SnapshotInternal, Estimator, MatN, MeasVec};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, SMatrix, SVector, Vector3};

#[derive(Clone)]
pub(super) struct BiasCal {
    pub enabled: bool,
    pub ins_inited: bool,
    pub ins_r: Vector3<f64>,
    pub ins_v: Vector3<f64>,
    pub t0: f64,
    pub from_pvt: bool,
    pub sens_r_ba: Matrix3<f64>,
    pub sens_v_ba: Matrix3<f64>,
    pub normal_a: SMatrix<f64, 9, 9>,
    pub normal_b: SVector<f64, 9>,
    pub obs_count: u32,
}

impl BiasCal {
    pub(super) fn new(enabled: bool) -> Self {
        Self {
            enabled,
            ins_inited: false,
            ins_r: Vector3::zeros(),
            ins_v: Vector3::zeros(),
            t0: 0.0,
            from_pvt: false,
            sens_r_ba: Matrix3::zeros(),
            sens_v_ba: Matrix3::zeros(),
            normal_a: SMatrix::<f64, 9, 9>::zeros(),
            normal_b: SVector::<f64, 9>::zeros(),
            obs_count: 0,
        }
    }

    pub(super) fn reset_interval(&mut self) {
        self.sens_r_ba.fill(0.0);
        self.sens_v_ba.fill(0.0);
    }
}

impl Estimator {
    pub fn set_precalibration(&mut self, enabled: bool, end_time: f64) {
        self.allow_bias_precalibration = enabled;
        self.precalibration_end_time = end_time;
    }

    pub fn apply_bias_calibration(&mut self) -> bool {
        let trace = std::env::var("TRACE").ok().is_some();
        if !self.bias_cal.enabled {
            if trace {
                eprintln!("[BA-CAL] skip: disabled");
            }
            return false;
        }
        if self.bias_cal.obs_count < 10 {
            if trace {
                eprintln!("[BA-CAL] skip: obs_count={}", self.bias_cal.obs_count);
            }
            return false;
        }

        let mut a = self.bias_cal.normal_a;
        for i in 0..9 {
            for j in (i + 1)..9 {
                let avg = 0.5 * (a[(i, j)] + a[(j, i)]);
                a[(i, j)] = avg;
                a[(j, i)] = avg;
            }
        }
        for i in 0..9 {
            a[(i, i)] += 1.0e-9;
        }

        let sigma_r0 = 5_000.0_f64;
        let sigma_v0 = 50.0_f64;
        if sigma_r0.is_finite() && sigma_r0 > 0.0 {
            let lambda = 1.0 / (sigma_r0 * sigma_r0);
            for i in 0..3 {
                a[(i, i)] += lambda;
            }
        }
        if sigma_v0.is_finite() && sigma_v0 > 0.0 {
            let lambda = 1.0 / (sigma_v0 * sigma_v0);
            for i in 3..6 {
                a[(i, i)] += lambda;
            }
        }

        let is_pvt = self.bias_cal.from_pvt;
        let sigma_ba = if is_pvt {
            150e-6_f64 * 9.80665
        } else {
            300e-6_f64 * 9.80665
        };
        if sigma_ba.is_finite() && sigma_ba > 0.0 {
            let lambda_ba = 1.0 / (sigma_ba * sigma_ba);
            for i in 6..9 {
                a[(i, i)] += lambda_ba;
            }
        }

        let mut a_fact = a;
        let mut chol = Cholesky::new(a_fact.clone());
        if chol.is_none() {
            for i in 0..9 {
                a_fact[(i, i)] += 1.0e-6;
            }
            chol = Cholesky::new(a_fact.clone());
        }
        let Some(chol) = chol else {
            if trace {
                let mut non_finite = 0usize;
                for i in 0..9 {
                    for j in 0..9 {
                        if !a_fact[(i, j)].is_finite() {
                            non_finite += 1;
                        }
                    }
                }
                eprintln!(
                    "[BA-CAL] skip: normal_a not SPD after jitter; non_finite={}",
                    non_finite
                );
            }
            return false;
        };
        let theta = chol.solve(&self.bias_cal.normal_b);
        let ba_post_cov = chol.inverse().fixed_view::<3, 3>(6, 6).into_owned();
        let mut delta = theta.fixed_rows::<3>(6).into_owned();
        let max_norm = if is_pvt {
            400e-6_f64 * 9.80665
        } else {
            800e-6_f64 * 9.80665
        };
        let n = delta.norm();
        if n.is_finite() && n > max_norm && max_norm > 0.0 {
            delta *= max_norm / n;
        }

        self.state.ba += delta;
        self.apply_bias_calibration_covariance(ba_post_cov);
        self.bias_cal.enabled = false;
        self.bias_cal.reset_interval();
        true
    }

    pub fn apply_bias_calibration_record(&mut self) -> bool {
        let t = self.state.t;
        let predicted_end_state = if self.seg_has_data {
            // RTS needs the propagation prior at this timestamp after the
            // calibration mutates `self.state`.
            Some(self.state.clone())
        } else {
            None
        };
        let applied = self.apply_bias_calibration();
        self.maybe_checkpoint_with_prediction(true, predicted_end_state);
        self.history.push_event(Event::BiasCalApply { t });
        applied
    }

    fn apply_bias_calibration_covariance(&mut self, ba_post_cov: Matrix3<f64>) {
        let mut p = self.cov();
        for row in 0..3 {
            let state_row = 15 + row;
            for col in 0..24 {
                p[(state_row, col)] = 0.0;
                p[(col, state_row)] = 0.0;
            }
        }
        for row in 0..3 {
            for col in 0..3 {
                p[(15 + row, 15 + col)] = ba_post_cov[(row, col)];
            }
        }
        crate::util::enforce_pd(&mut p);
        self.l_cov = super::math::chol_lower_with_enforce(p, "apply_bias_calibration_covariance");
    }

    #[inline]
    pub(super) fn snapshot_internal(&self) -> SnapshotInternal {
        SnapshotInternal {
            bias_cal: self.bias_cal.clone(),
            gnss_reject_streak: self.gnss_reject_streak,
            robust_gate_alpha: self.robust_gate_alpha,
            robust_r_scale_cap: self.robust_r_scale_cap,
            robust_student_nu: self.robust_student_nu,
            iae_gamma: self.iae_gamma,
            iae_scales: self.iae_scales,
            aero_learning_enabled: self.aero_learning_enabled,
            aero_learning_rate: self.aero_learning_rate,
            last_v_meas: self.last_v_meas,
            learned_aero_correction: self.learned_aero_correction,
        }
    }

    #[inline]
    pub(super) fn restore_internal(&mut self, internal: SnapshotInternal) {
        self.bias_cal = internal.bias_cal;
        self.gnss_reject_streak = internal.gnss_reject_streak;
        self.robust_gate_alpha = internal.robust_gate_alpha;
        self.robust_r_scale_cap = internal.robust_r_scale_cap;
        self.robust_student_nu = internal.robust_student_nu;
        self.iae_gamma = internal.iae_gamma;
        self.iae_scales = internal.iae_scales;
        self.aero_learning_enabled = internal.aero_learning_enabled;
        self.aero_learning_rate = internal.aero_learning_rate;
        self.last_v_meas = internal.last_v_meas;
        self.learned_aero_correction = internal.learned_aero_correction;
    }

    #[inline]
    pub(super) fn bias_cal_observe_rv(
        &mut self,
        t: f64,
        pos: Vector3<f64>,
        vel: Vector3<f64>,
        r_cov: &Matrix6<f64>,
    ) {
        if !self.bias_cal.enabled {
            return;
        }

        if !self.bias_cal.ins_inited {
            self.bias_cal.ins_inited = true;
            self.bias_cal.ins_r = pos;
            self.bias_cal.ins_v = vel;
            self.bias_cal.t0 = t;
            self.bias_cal.sens_r_ba.fill(0.0);
            self.bias_cal.sens_v_ba.fill(0.0);
            self.bias_cal.normal_a.fill(0.0);
            self.bias_cal.normal_b.fill(0.0);
            self.bias_cal.obs_count = 0;
            return;
        }

        let mut res = MeasVec::<6>::zeros();
        res.fixed_rows_mut::<3>(0)
            .copy_from(&(pos - self.bias_cal.ins_r));
        res.fixed_rows_mut::<3>(3)
            .copy_from(&(vel - self.bias_cal.ins_v));

        let dt_rel = (t - self.bias_cal.t0).max(0.0);
        let mut s = SMatrix::<f64, 6, 9>::zeros();
        for i in 0..3 {
            s[(i, i)] = 1.0;
            s[(i, i + 3)] = dt_rel;
            s[(i + 3, i + 3)] = 1.0;
            for j in 0..3 {
                s[(i, j + 6)] = self.bias_cal.sens_r_ba[(i, j)];
                s[(i + 3, j + 6)] = self.bias_cal.sens_v_ba[(i, j)];
            }
        }

        const POS_OUTLIER_M: f64 = 100.0;
        const ROBUST_CAP: f64 = 1_000.0;
        for i in 0..3 {
            let abs = res[i].abs();
            if !(abs.is_finite() && abs > POS_OUTLIER_M) {
                continue;
            }
            let scale = (POS_OUTLIER_M / abs).max(1.0 / ROBUST_CAP);
            res[i] *= scale;
            for j in 0..9 {
                s[(i, j)] *= scale;
            }
        }

        let mut r_cov = r_cov.clone();
        for i in 0..6 {
            for j in (i + 1)..6 {
                let avg = 0.5 * (r_cov[(i, j)] + r_cov[(j, i)]);
                r_cov[(i, j)] = avg;
                r_cov[(j, i)] = avg;
            }
        }
        let r_chol = match Cholesky::new(r_cov) {
            Some(c) => c,
            None => return,
        };
        let w_res = r_chol.solve(&res);
        let w_s = r_chol.solve(&s);
        self.bias_cal.normal_b += s.transpose() * w_res;
        self.bias_cal.normal_a += s.transpose() * w_s;
        self.bias_cal.obs_count += 1;
    }

    #[inline]
    pub(super) fn bias_update_profile(&self, t_meas: f64) -> (bool, f64) {
        let base = self.noise.gps_raw_max_ba_update;
        if !base.is_finite() || base <= 0.0 {
            return (false, 0.0);
        }
        if self.noise.gps_ba_update_gain <= 0.0 {
            return (false, 0.0);
        }

        let in_precal =
            self.allow_bias_precalibration && t_meas <= self.precalibration_end_time + 1e-9;
        let in_postcal = if let Some((bs, be)) = self.blackout_window {
            let window = (be - bs).max(10.0);
            t_meas >= be && t_meas <= be + window
        } else {
            false
        };

        let max_step = if in_precal || in_postcal {
            base
        } else {
            0.5 * base
        };
        (true, max_step)
    }
}

pub(super) fn initial_cov() -> MatN {
    let mut p = MatN::zeros();
    for i in 0..3 {
        p[(i, i)] = (5_000.0_f64).powi(2);
    }
    for i in 3..6 {
        p[(i, i)] = (50.0_f64).powi(2);
    }
    let att_var = (5_f64.to_radians()).powi(2);
    for i in 6..9 {
        p[(i, i)] = att_var;
    }
    let omg_var = (2_f64.to_radians()).powi(2);
    for i in 9..12 {
        p[(i, i)] = omg_var;
    }
    let bg_var = (0.3_f64.to_radians() / 3600.0).powi(2);
    for i in 12..15 {
        p[(i, i)] = bg_var;
    }
    let ba_var = (300e-6_f64 * 9.80665).powi(2);
    for i in 15..18 {
        p[(i, i)] = ba_var;
    }
    let a_res_var = 0.05_f64.powi(2);
    for i in 18..21 {
        p[(i, i)] = a_res_var;
    }
    let t_res_var = 25.0;
    for i in 21..24 {
        p[(i, i)] = t_res_var;
    }
    p
}
