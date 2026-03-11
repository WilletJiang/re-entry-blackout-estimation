use super::{
    metrics, MonteCarloResult, SimConfig, SimRun, CHI2_3_HI, CHI2_3_LO, CHI2_6_HI, CHI2_6_LO,
};
use crate::estimator::Estimator;
use crate::sim::BlackoutTrace;
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, SVector, Vector3};

pub(super) struct SimFinish {
    pub result: MonteCarloResult,
    pub trace: Option<BlackoutTrace>,
}

impl SimRun {
    pub(super) fn accumulate_recovery_errors(&mut self, t: f64, step: usize) {
        let pos_err = self.est.state.r - self.truth.r;
        let vel_err = self.est.state.v - self.truth.v;
        if t < self.cfg.blackout_end || t > self.cfg.blackout_end + 2.0 {
            return;
        }
        self.pos_err2_acc += pos_err.norm_squared();
        self.vel_err2_acc += vel_err.norm_squared();
        let dq = self.est.state.q * self.truth.q.inverse();
        self.att_err2_acc += dq.angle().powi(2);
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
    }

    pub(super) fn finish(self) -> SimFinish {
        let (nees_stats, nees_phase_stats) = compute_nees_stats(self.nees_samples);

        let (sm_rms_pos, sm_rms_vel) =
            smoother_metrics(&self.est, &self.cfg, &self.truth_times, &self.truth_hist);

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

        let rms_pos = if self.win_samples > 0.0 {
            (self.win_pos_err2_acc / self.win_samples).sqrt()
        } else {
            0.0
        };
        let rms_vel = if self.win_samples > 0.0 {
            (self.win_vel_err2_acc / self.win_samples).sqrt()
        } else {
            0.0
        };
        let rms_att_deg = if self.err_samples > 0.0 {
            (self.att_err2_acc / self.err_samples).sqrt().to_degrees()
        } else {
            0.0
        };

        let bg_rms_si = if self.win_samples > 0.0 {
            (self.win_bg_err2_acc / self.win_samples).sqrt()
        } else {
            0.0
        };
        let ba_rms_si = if self.win_samples > 0.0 {
            (self.win_ba_err2_acc / self.win_samples).sqrt()
        } else {
            0.0
        };
        let bg_rms_dps = bg_rms_si.to_degrees();
        let ba_rms_ug = if ba_rms_si > 0.0 {
            ba_rms_si / 9.80665 * 1.0e6
        } else {
            0.0
        };

        let pass = rms_pos <= 40.0 / 3.0
            && rms_vel <= 0.5 / 3.0
            && rms_att_deg <= 0.05
            && bg_rms_dps <= 0.005
            && ba_rms_ug <= 120.0;

        if std::env::var("TRACE").ok().is_some() {
            eprintln!(
                "[WIN-SUMMARY] samples={} rms_pos={:.3} rms_vel={:.3}",
                self.win_samples, rms_pos, rms_vel
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
            rms_pos,
            rms_vel,
            rms_att_deg,
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
            passed: pass,
        };
        let trace = self.blackout_trace.map(|t| t.finish());
        SimFinish { result, trace }
    }
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

fn smoother_metrics(
    est: &Estimator,
    cfg: &SimConfig,
    truth_times: &[f64],
    truth_hist: &[crate::types::State],
) -> (f64, f64) {
    let Some((xs, _ps)) = est.rts_smooth(2.0) else {
        return (0.0, 0.0);
    };
    let target_t = (cfg.t_final - 2.0).max(0.0);
    let idx = truth_times
        .iter()
        .enumerate()
        .rev()
        .find(|(_, &ti)| ti <= target_t)
        .map(|(i, _)| i)
        .unwrap_or(0);
    let truth_at = truth_hist
        .get(idx)
        .unwrap_or_else(|| truth_hist.first().unwrap());
    let pos_err_s = xs.r - truth_at.r;
    let vel_err_s = xs.v - truth_at.v;
    (pos_err_s.norm(), vel_err_s.norm())
}
