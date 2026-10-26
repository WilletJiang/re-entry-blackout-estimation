use crate::estimator::MatN;
use crate::sim::SimConfig;
use nalgebra::{Matrix3, Matrix6, SVector, UnitQuaternion, Vector3};

#[derive(Clone)]
pub(super) struct NeesSample {
    pub phase: u8,
    pub gnss_accepted: bool,
    pub err_rv: SVector<f64, 6>,
    pub p_rv: Matrix6<f64>,
    pub err_ba: Vector3<f64>,
    pub p_ba: Matrix3<f64>,
}

#[inline]
pub(super) fn phase_id(t: f64, cfg: &SimConfig) -> usize {
    if t < cfg.blackout_start {
        return 0;
    }
    if t < cfg.blackout_end - 1e-9 {
        return 1;
    }
    if t <= cfg.blackout_end + 2.0 + 1e-9 {
        return 2;
    }
    3
}

#[derive(Default)]
pub(super) struct RecoveryCoverage {
    pub samples: usize,
    first: f64,
    last: f64,
    broken: bool,
}

impl RecoveryCoverage {
    pub fn record(&mut self, t: f64, dt: f64) {
        if self.samples == 0 {
            self.first = t;
        } else if (t - self.last - dt).abs() > 1e-9 {
            self.broken = true;
        }
        self.broken |= !t.is_finite();
        self.last = t;
        self.samples += 1;
    }

    pub fn complete(&self, cfg: &SimConfig, final_t: f64) -> bool {
        if self.broken || self.samples == 0 || !cfg.dt.is_finite() || cfg.dt <= 0.0 || cfg.dt > 2.0
        {
            return false;
        }
        let start_step = ((cfg.blackout_end - 1e-9) / cfg.dt).ceil().max(1.0);
        let end_step = ((cfg.blackout_end + 2.0 + 1e-9) / cfg.dt).floor();
        self.samples as f64 == end_step - start_step + 1.0
            && (self.first - start_step * cfg.dt).abs() <= 1e-9
            && (self.last - end_step * cfg.dt).abs() <= 1e-9
            && final_t >= cfg.blackout_end + 2.0 - 1e-9
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_rejects_partial_missing_duplicate_and_nonfinite_samples() {
        let cfg = SimConfig {
            dt: 0.5,
            blackout_end: 2.0,
            ..SimConfig::default()
        };
        for times in [
            vec![],
            vec![2.0, 2.5],
            vec![2.0, 2.5, 3.5, 4.0],
            vec![2.0, 2.5, 2.5, 3.5, 4.0],
            vec![2.0, 2.5, f64::NAN, 3.5, 4.0],
        ] {
            let mut coverage = RecoveryCoverage::default();
            for t in times {
                coverage.record(t, cfg.dt);
            }
            assert!(!coverage.complete(&cfg, 4.0));
        }
        let mut coverage = RecoveryCoverage::default();
        for t in [2.0, 2.5, 3.0, 3.5, 4.0] {
            coverage.record(t, cfg.dt);
        }
        assert!(coverage.complete(&cfg, 4.0));
        assert!(!coverage.complete(&cfg, 3.5));
    }

    #[test]
    fn recovery_phase_includes_blackout_exit() {
        let cfg = SimConfig::default();
        assert_eq!(phase_id(160.0, &cfg), 2);
        assert_eq!(phase_id(162.0, &cfg), 2);
        assert_eq!(phase_id(159.99, &cfg), 1);
    }
}

#[inline]
pub(super) fn cov_rv_from_l(l: &MatN) -> Matrix6<f64> {
    let mut p = Matrix6::<f64>::zeros();
    for i in 0..6 {
        for j in 0..=i {
            let mut s = 0.0f64;
            for k in 0..=j {
                s += l[(i, k)] * l[(j, k)];
            }
            p[(i, j)] = s;
            p[(j, i)] = s;
        }
    }
    p
}

#[inline]
pub(super) fn cov_ba_from_l(l: &MatN) -> Matrix3<f64> {
    let mut p = Matrix3::<f64>::zeros();
    for ii in 0..3 {
        let i = 15 + ii;
        for jj in 0..=ii {
            let j = 15 + jj;
            let mut s = 0.0f64;
            for k in 0..=j {
                s += l[(i, k)] * l[(j, k)];
            }
            p[(ii, jj)] = s;
            p[(jj, ii)] = s;
        }
    }
    p
}

pub(super) fn find_truth_q_at(
    times: &[f64],
    hist: &[crate::types::State],
    t_target: f64,
) -> Option<UnitQuaternion<f64>> {
    let (lower, upper, weight) = truth_bracket(times, t_target)?;
    Some(hist.get(lower)?.q.slerp(&hist.get(upper)?.q, weight))
}

pub(super) fn truth_bracket(times: &[f64], t: f64) -> Option<(usize, usize, f64)> {
    if !t.is_finite()
        || times.is_empty()
        || t < times[0] - 1e-9
        || t > times[times.len() - 1] + 1e-9
    {
        return None;
    }
    let upper = times.partition_point(|&time| time < t).min(times.len() - 1);
    if (times[upper] - t).abs() < 1e-9 || upper == 0 {
        return Some((upper, upper, 0.0));
    }
    let lower = upper - 1;
    if (t - times[lower]).abs() < 1e-9 {
        return Some((lower, lower, 0.0));
    }
    Some((
        lower,
        upper,
        (t - times[lower]) / (times[upper] - times[lower]),
    ))
}
