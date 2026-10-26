use super::{
    metrics, MonteCarloResult, SimConfig, SimRun, CHI2_3_HI, CHI2_3_LO, CHI2_6_HI, CHI2_6_LO,
};
use crate::sim::BlackoutTrace;
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, SVector, Vector3};

pub(super) struct SimFinish {
    pub result: MonteCarloResult,
    pub trace: Option<BlackoutTrace>,
}

impl SimRun {
    pub(super) fn accumulate_recovery_window_stats(&mut self, step: usize) {
        let t = self.est.state.t;
        if !in_recovery_window(t, &self.cfg) {
            return;
        }
        self.recovery_coverage.record(t, self.cfg.dt);
        let pos_err = self.est.state.r - self.truth.r;
        let vel_err = self.est.state.v - self.truth.v;
        self.pos_err2_acc += pos_err.norm_squared();
        self.vel_err2_acc += vel_err.norm_squared();
        let dq = self.est.state.q * self.truth.q.inverse();
        self.att_err2_acc += dq.angle().powi(2);
        self.omega_err2_acc += (self.est.state.omega - self.truth.omega).norm_squared();
        if self.trace && (step % 200 == 0) {
            eprintln!(
                "[t={:.3}] |pos_err|={:.3} m |vel_err|={:.3} m/s",
                t,
                pos_err.norm(),
                vel_err.norm()
            );
        }
        let bg_err = self.est.state.bg - self.truth.bg;
        let ba_err = self.est.state.ba - self.truth.ba;
        self.bg_err2_acc += bg_err.norm_squared();
        self.ba_err2_acc += ba_err.norm_squared();
        self.err_samples += 1.0;

        let p = self.est.cov();
        self.recovery
            .push(crate::sim::campaign::RecoverySample::new(
                &self.truth,
                &self.est.state,
                &p,
            ));
        self.recovery_numerically_valid &= state_is_finite(&self.est.state)
            && state_is_finite(&self.truth)
            && p.iter().all(|v| v.is_finite())
            && [
                self.pos_err2_acc,
                self.vel_err2_acc,
                self.att_err2_acc,
                self.omega_err2_acc,
                self.bg_err2_acc,
                self.ba_err2_acc,
            ]
            .iter()
            .all(|v| v.is_finite());
        self.pos_sigma3_max_m = self.pos_sigma3_max_m.max(sigma3_norm_from_cov(
            &p.fixed_view::<3, 3>(0, 0).into_owned(),
        ));
        self.vel_sigma3_max_mps = self.vel_sigma3_max_mps.max(sigma3_norm_from_cov(
            &p.fixed_view::<3, 3>(3, 3).into_owned(),
        ));
        self.att_sigma3_max_deg = self
            .att_sigma3_max_deg
            .max(sigma3_norm_from_cov(&p.fixed_view::<3, 3>(6, 6).into_owned()).to_degrees());
        self.omega_sigma3_max_dps = self
            .omega_sigma3_max_dps
            .max(sigma3_norm_from_cov(&p.fixed_view::<3, 3>(9, 9).into_owned()).to_degrees());
        self.bg_sigma3_max_dps = self
            .bg_sigma3_max_dps
            .max(sigma3_norm_from_cov(&p.fixed_view::<3, 3>(12, 12).into_owned()).to_degrees());
        self.ba_sigma3_max_ug = self.ba_sigma3_max_ug.max(
            sigma3_norm_from_cov(&p.fixed_view::<3, 3>(15, 15).into_owned()) / 9.80665 * 1.0e6,
        );
    }

    pub(super) fn capture_recovery_at_deadline(&mut self) {
        let deadline = self.cfg.blackout_end + 2.0;
        if (self.est.state.t - deadline).abs() > 1e-8 {
            return;
        }
        let mut samples = Vec::with_capacity(201);
        let success = self.est.visit_smoothed(2.0, |state, covariance| {
            if !in_recovery_window(state.t, &self.cfg)
                || (state.t * 100.0 - (state.t * 100.0).round()).abs() > 1e-6
            {
                return;
            }
            // All checkpoint times are on the 400 Hz truth grid. The same
            // tolerance used for delayed sensors prevents an off-by-one frame.
            if let Some((lo, hi, _)) = metrics::truth_bracket(&self.truth_times, state.t) {
                if lo == hi {
                    samples.push(crate::sim::campaign::RecoverySample::new(
                        &self.truth_hist[lo],
                        state,
                        covariance,
                    ));
                }
            }
        });
        if success.is_some() {
            samples.reverse();
            self.recovery_at_deadline = samples;
        }
    }

    pub(super) fn accumulate_smoothed_recovery_window_stats(&mut self) {
        if !self.cfg.enable_smoothing_diagnostics {
            return;
        }
        let current_t = self.est.state.t;
        if current_t + 1.0e-9 < self.cfg.history_sec {
            return;
        }
        let target_t = current_t - self.cfg.history_sec;
        if !in_recovery_window(target_t, &self.cfg) {
            return;
        }

        let Some((smoothed_state, _smoothed_cov)) = self.est.rts_smooth(self.cfg.history_sec)
        else {
            return;
        };
        if !in_recovery_window(smoothed_state.t, &self.cfg) {
            return;
        }
        let Some(truth_at) =
            truth_state_at_or_before(&self.truth_times, &self.truth_hist, smoothed_state.t)
        else {
            return;
        };
        let pos_err = smoothed_state.r - truth_at.r;
        let vel_err = smoothed_state.v - truth_at.v;
        self.sm_pos_err2_acc += pos_err.norm_squared();
        self.sm_vel_err2_acc += vel_err.norm_squared();
        self.sm_samples += 1.0;
    }

    pub(super) fn record_health_peak(&mut self) {
        let health = self.est.health();
        self.history_peak_snapshots = self.history_peak_snapshots.max(health.history_snapshots);
        self.history_peak_events = self.history_peak_events.max(health.history_events);
        self.history_peak_props = self.history_peak_props.max(health.history_props);
        self.history_peak_capacity_bytes = self
            .history_peak_capacity_bytes
            .max(health.history_capacity_bytes);
    }

    pub(super) fn finish(self) -> SimFinish {
        let (nees_stats, nees_phase_stats) = compute_nees_stats(self.nees_samples);

        let nees_mean = nees_stats.mean;
        let nees_covered_run = nees_stats.covered;
        let nis_mean = if self.nis_ct > 0 {
            self.nis_acc / (self.nis_ct as f64)
        } else {
            f64::NAN
        };
        let nis_covered_run = if self.nis_ct > 0 {
            (self.nis_hit as f64) / (self.nis_ct as f64)
        } else {
            f64::NAN
        };

        let rms_pos = rms_from_accumulator(self.pos_err2_acc, self.err_samples);
        let rms_vel = rms_from_accumulator(self.vel_err2_acc, self.err_samples);
        let rms_att_deg = if self.err_samples > 0.0 {
            (self.att_err2_acc / self.err_samples).sqrt().to_degrees()
        } else {
            f64::NAN
        };

        let bg_rms_si = rms_from_accumulator(self.bg_err2_acc, self.err_samples);
        let ba_rms_si = rms_from_accumulator(self.ba_err2_acc, self.err_samples);
        let bg_rms_dps = bg_rms_si.to_degrees();
        let ba_rms_ug = if ba_rms_si > 0.0 {
            ba_rms_si / 9.80665 * 1.0e6
        } else {
            ba_rms_si
        };
        let sm_rms_pos = rms_from_accumulator(self.sm_pos_err2_acc, self.sm_samples);
        let sm_rms_vel = rms_from_accumulator(self.sm_vel_err2_acc, self.sm_samples);

        let recovery_window_complete = self.recovery_coverage.complete(&self.cfg, self.est.state.t);
        let pass = recovery_window_complete
            && self.recovery_numerically_valid
            && passes_covariance_limits([
                self.pos_sigma3_max_m,
                self.vel_sigma3_max_mps,
                self.att_sigma3_max_deg,
                self.omega_sigma3_max_dps,
                self.bg_sigma3_max_dps,
                self.ba_sigma3_max_ug,
            ]);

        if std::env::var("TRACE").ok().is_some() {
            eprintln!(
                "[WIN-SUMMARY] samples={} rms_pos={:.3} rms_vel={:.3} pos_3sigma_max={:.3} vel_3sigma_max={:.3}",
                self.err_samples, rms_pos, rms_vel, self.pos_sigma3_max_m, self.vel_sigma3_max_mps
            );
        }

        let (nees_mean_phase, nees_cov_phase) = phase_means_and_cov(
            nees_phase_stats.sum,
            nees_phase_stats.ct,
            nees_phase_stats.hit,
        );
        let (nis_mean_phase, nis_cov_phase) =
            phase_means_and_cov(self.nis_phase_acc, self.nis_phase_ct, self.nis_phase_hit);

        let result = MonteCarloResult {
            recovery: self.recovery,
            recovery_at_deadline: self.recovery_at_deadline,
            aux_nis: self.aux_nis,
            nis_count: self.nis_ct,
            nis_hits: self.nis_hit,
            seed: self.seed,
            trajectory: self.trajectory,
            recovery_samples: self.recovery_coverage.samples,
            recovery_window_complete,
            recovery_numerically_valid: self.recovery_numerically_valid,
            rms_omega_dps: rms_from_accumulator(self.omega_err2_acc, self.err_samples).to_degrees(),
            sigma3_omega_max_dps: self.omega_sigma3_max_dps,
            rms_pos,
            rms_vel,
            rms_att_deg,
            sigma3_pos_max_m: self.pos_sigma3_max_m,
            sigma3_vel_max_mps: self.vel_sigma3_max_mps,
            sigma3_att_max_deg: self.att_sigma3_max_deg,
            sigma3_bg_max_dps: self.bg_sigma3_max_dps,
            sigma3_ba_max_ug: self.ba_sigma3_max_ug,
            nees_mean,
            nis_mean,
            nees_covered: nees_covered_run,
            nis_covered: nis_covered_run,
            nees_r_mean: nees_stats.r_mean,
            nees_v_mean: nees_stats.v_mean,
            nees_ba_mean: nees_stats.ba_mean,
            nees_r_covered: nees_stats.r_covered,
            nees_v_covered: nees_stats.v_covered,
            nees_ba_covered: nees_stats.ba_covered,
            nees_mean_pre: nees_mean_phase[0],
            nees_mean_blk: nees_mean_phase[1],
            nees_mean_post: nees_mean_phase[2],
            nees_mean_tail: nees_mean_phase[3],
            nees_covered_pre: nees_cov_phase[0],
            nees_covered_blk: nees_cov_phase[1],
            nees_covered_post: nees_cov_phase[2],
            nees_covered_tail: nees_cov_phase[3],
            nis_mean_pre: nis_mean_phase[0],
            nis_mean_blk: nis_mean_phase[1],
            nis_mean_post: nis_mean_phase[2],
            nis_mean_tail: nis_mean_phase[3],
            nis_covered_pre: nis_cov_phase[0],
            nis_covered_blk: nis_cov_phase[1],
            nis_covered_post: nis_cov_phase[2],
            nis_covered_tail: nis_cov_phase[3],
            post_gnss_updates: self.post_gnss_updates,
            post_gnss_accepted: self.post_gnss_accepted,
            post_sat_mean: if self.cfg.use_tight_coupling && self.post_sat_ct > 0 {
                (self.post_sat_sum as f64) / (self.post_sat_ct as f64)
            } else {
                f64::NAN
            },
            post_sat_min: if self.post_sat_min == usize::MAX {
                0
            } else {
                self.post_sat_min
            },
            smoothed_rms_pos: sm_rms_pos,
            smoothed_rms_vel: sm_rms_vel,
            bias_bg_rms_dps: bg_rms_dps,
            bias_ba_rms_ug: ba_rms_ug,
            history_snapshots: self.history_peak_snapshots,
            history_events: self.history_peak_events,
            history_props: self.history_peak_props,
            history_capacity_bytes: self.history_peak_capacity_bytes,
            passed: pass,
        };
        let trace = self.blackout_trace.map(|t| t.finish());
        SimFinish { result, trace }
    }
}

fn in_recovery_window(t: f64, cfg: &SimConfig) -> bool {
    t >= cfg.blackout_end - 1.0e-9 && t <= cfg.blackout_end + 2.0 + 1.0e-9
}

fn rms_from_accumulator(sum_sq: f64, samples: f64) -> f64 {
    if samples > 0.0 {
        (sum_sq / samples).sqrt()
    } else {
        f64::NAN
    }
}

fn sigma3_norm_from_cov(cov: &Matrix3<f64>) -> f64 {
    if !cov.iter().all(|v| v.is_finite()) || cov.diagonal().iter().any(|v| *v < 0.0) {
        return f64::INFINITY;
    }
    3.0 * cov.trace().sqrt()
}

fn truth_state_at_or_before<'a>(
    truth_times: &[f64],
    truth_hist: &'a [crate::types::State],
    t_target: f64,
) -> Option<&'a crate::types::State> {
    let idx = truth_times
        .iter()
        .enumerate()
        .rev()
        .find(|(_, &ti)| ti <= t_target)
        .map(|(i, _)| i)
        .or_else(|| (!truth_hist.is_empty()).then_some(0))?;
    truth_hist.get(idx)
}

fn phase_means_and_cov(sum: [f64; 4], ct: [usize; 4], hit: [usize; 4]) -> ([f64; 4], [f64; 4]) {
    let mut mean = [f64::NAN; 4];
    let mut cov = [f64::NAN; 4];
    for ph in 0..4 {
        if ct[ph] == 0 {
            continue;
        }
        mean[ph] = sum[ph] / (ct[ph] as f64);
        cov[ph] = (hit[ph] as f64) / (ct[ph] as f64);
    }
    (mean, cov)
}

struct NeesStats {
    mean: f64,
    covered: f64,
    r_mean: f64,
    v_mean: f64,
    ba_mean: f64,
    r_covered: f64,
    v_covered: f64,
    ba_covered: f64,
}

struct PhaseStats {
    sum: [f64; 4],
    ct: [usize; 4],
    hit: [usize; 4],
}

fn compute_nees_stats(samples: Vec<metrics::NeesSample>) -> (NeesStats, PhaseStats) {
    let mut nees_acc = 0.0f64;
    let mut nees_ct = 0usize;
    let mut nees_hit = 0usize;
    let mut nees_phase_sum = [0.0f64; 4];
    let mut nees_phase_ct = [0usize; 4];
    let mut nees_phase_hit = [0usize; 4];
    let mut nees_r_acc = 0.0f64;
    let mut nees_v_acc = 0.0f64;
    let mut nees_ba_acc = 0.0f64;
    let mut nees_r_ct = 0usize;
    let mut nees_v_ct = 0usize;
    let mut nees_ba_ct = 0usize;
    let mut nees_r_hit = 0usize;
    let mut nees_v_hit = 0usize;
    let mut nees_ba_hit = 0usize;

    for s in samples {
        let ph = s.phase as usize;

        if let Some(v) = nees_6(&s.err_rv, &s.p_rv) {
            nees_phase_sum[ph] += v;
            nees_phase_ct[ph] += 1;
            let in_gate6 = v >= CHI2_6_LO && v <= CHI2_6_HI;
            if in_gate6 {
                nees_phase_hit[ph] += 1;
            }
            if s.gnss_accepted {
                nees_acc += v;
                nees_ct += 1;
            }
            if s.gnss_accepted && in_gate6 {
                nees_hit += 1;
            }
        }

        if !s.gnss_accepted {
            continue;
        }

        let p_rr = s.p_rv.fixed_view::<3, 3>(0, 0).into_owned();
        let p_vv = s.p_rv.fixed_view::<3, 3>(3, 3).into_owned();
        let err_r = s.err_rv.fixed_rows::<3>(0).into_owned();
        let err_v = s.err_rv.fixed_rows::<3>(3).into_owned();

        if let Some(vr) = nees_3(&err_r, &p_rr) {
            nees_r_acc += vr;
            nees_r_ct += 1;
            if vr >= CHI2_3_LO && vr <= CHI2_3_HI {
                nees_r_hit += 1;
            }
        }
        if let Some(vv) = nees_3(&err_v, &p_vv) {
            nees_v_acc += vv;
            nees_v_ct += 1;
            if vv >= CHI2_3_LO && vv <= CHI2_3_HI {
                nees_v_hit += 1;
            }
        }
        if let Some(vb) = nees_3(&s.err_ba, &s.p_ba) {
            nees_ba_acc += vb;
            nees_ba_ct += 1;
            if vb >= CHI2_3_LO && vb <= CHI2_3_HI {
                nees_ba_hit += 1;
            }
        }
    }

    let nees_mean = if nees_ct > 0 {
        nees_acc / (nees_ct as f64)
    } else {
        f64::NAN
    };
    let nees_cov = if nees_ct > 0 {
        (nees_hit as f64) / (nees_ct as f64)
    } else {
        f64::NAN
    };
    let nees_r_mean = if nees_r_ct > 0 {
        nees_r_acc / (nees_r_ct as f64)
    } else {
        f64::NAN
    };
    let nees_v_mean = if nees_v_ct > 0 {
        nees_v_acc / (nees_v_ct as f64)
    } else {
        f64::NAN
    };
    let nees_ba_mean = if nees_ba_ct > 0 {
        nees_ba_acc / (nees_ba_ct as f64)
    } else {
        f64::NAN
    };
    let nees_r_cov = if nees_r_ct > 0 {
        (nees_r_hit as f64) / (nees_r_ct as f64)
    } else {
        f64::NAN
    };
    let nees_v_cov = if nees_v_ct > 0 {
        (nees_v_hit as f64) / (nees_v_ct as f64)
    } else {
        f64::NAN
    };
    let nees_ba_cov = if nees_ba_ct > 0 {
        (nees_ba_hit as f64) / (nees_ba_ct as f64)
    } else {
        f64::NAN
    };

    (
        NeesStats {
            mean: nees_mean,
            covered: nees_cov,
            r_mean: nees_r_mean,
            v_mean: nees_v_mean,
            ba_mean: nees_ba_mean,
            r_covered: nees_r_cov,
            v_covered: nees_v_cov,
            ba_covered: nees_ba_cov,
        },
        PhaseStats {
            sum: nees_phase_sum,
            ct: nees_phase_ct,
            hit: nees_phase_hit,
        },
    )
}

fn nees_6(err: &SVector<f64, 6>, p: &Matrix6<f64>) -> Option<f64> {
    let ch = Cholesky::new(*p)?;
    let x = ch.solve(err);
    let v = err.dot(&x);
    v.is_finite().then_some(v)
}

fn nees_3(err: &Vector3<f64>, p: &Matrix3<f64>) -> Option<f64> {
    let ch = Cholesky::new(*p)?;
    let x = ch.solve(err);
    let v = err.dot(&x);
    v.is_finite().then_some(v)
}

fn passes_covariance_limits(maxima: [f64; 6]) -> bool {
    maxima
        .into_iter()
        .zip([40.0, 0.5, 0.15, 0.05, 0.005, 120.0])
        .all(|(value, limit)| value.is_finite() && value >= 0.0 && value <= limit)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_six_limits_are_enforced_and_nonfinite_values_fail() {
        let limits = [40.0, 0.5, 0.15, 0.05, 0.005, 120.0];
        assert!(passes_covariance_limits(limits));
        for i in 0..6 {
            for bad in [limits[i] * 1.01, f64::NAN, f64::INFINITY, -1.0] {
                let mut values = limits;
                values[i] = bad;
                assert!(!passes_covariance_limits(values));
            }
        }
    }

    #[test]
    fn invalid_covariance_cannot_be_hidden_by_maximum_accumulation() {
        for bad in [f64::NAN, f64::INFINITY, -1.0] {
            let mut cov = Matrix3::identity();
            cov[(0, 0)] = bad;
            assert!(0.0_f64.max(sigma3_norm_from_cov(&cov)).is_infinite());
        }
        let mut cov = Matrix3::identity();
        cov[(0, 1)] = f64::NAN;
        assert!(sigma3_norm_from_cov(&cov).is_infinite());
    }
}

pub(super) fn state_is_finite(state: &crate::types::State) -> bool {
    state.t.is_finite()
        && state.q.quaternion().coords.iter().all(|v| v.is_finite())
        && [
            &state.r,
            &state.v,
            &state.omega,
            &state.bg,
            &state.ba,
            &state.aero_res_force,
            &state.aero_res_torque,
        ]
        .into_iter()
        .all(|vector| vector.iter().all(|v| v.is_finite()))
}
