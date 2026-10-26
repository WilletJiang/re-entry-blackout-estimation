//! Recovery evidence sampled independently of measurement acceptance.
use super::MonteCarloResult;
use crate::types::State;
use crate::MatN;
use nalgebra::SVector;
use std::io::Write;

pub const LIMITS: [f64; 6] = [40.0, 0.5, 0.15, 0.05, 0.005, 120.0];
const NAMES: [&str; 6] = [
    "position_m",
    "velocity_mps",
    "attitude_deg",
    "rate_dps",
    "gyro_bias_dps",
    "accel_bias_ug",
];
const CHI18: [f64; 2] = [8.230746194756668, 31.526378440386626];

#[derive(Clone)]
pub struct RecoverySample {
    pub t: f64,
    pub squared_error: [f64; 6],
    pub variance: [f64; 6],
    pub nees: f64,
    pub block_nees: [f64; 6],
}

impl RecoverySample {
    pub(crate) fn new(truth: &State, estimate: &State, p: &MatN) -> Self {
        let errors = [
            truth.r - estimate.r,
            truth.v - estimate.v,
            (truth.q * estimate.q.inverse()).scaled_axis(),
            truth.omega - estimate.omega,
            truth.bg - estimate.bg,
            truth.ba - estimate.ba,
        ];
        let units = [
            1.0,
            1.0,
            180.0 / std::f64::consts::PI,
            180.0 / std::f64::consts::PI,
            180.0 / std::f64::consts::PI,
            1e6 / 9.80665,
        ];
        let mut error = SVector::<f64, 18>::zeros();
        let mut squared_error = [0.0; 6];
        let mut variance = [0.0; 6];
        let mut block_nees = [0.0; 6];
        for i in 0..6 {
            error.fixed_rows_mut::<3>(3 * i).copy_from(&errors[i]);
            squared_error[i] = errors[i].norm_squared() * units[i].powi(2);
            let block = p.fixed_view::<3, 3>(3 * i, 3 * i).into_owned();
            variance[i] = block.trace() * units[i].powi(2);
            block_nees[i] = block
                .cholesky()
                .map_or(f64::INFINITY, |c| errors[i].dot(&c.solve(&errors[i])));
        }
        let nees = p
            .fixed_view::<18, 18>(0, 0)
            .into_owned()
            .cholesky()
            .map_or(f64::INFINITY, |c| error.dot(&c.solve(&error)));
        Self {
            t: estimate.t,
            squared_error,
            variance,
            nees,
            block_nees,
        }
    }
}

pub fn write_assessment(results: &[MonteCarloResult], out_dir: &str) -> anyhow::Result<bool> {
    write_selected_assessment(results, out_dir, false)
}

pub fn write_deadline_assessment(
    results: &[MonteCarloResult],
    out_dir: &str,
) -> anyhow::Result<bool> {
    write_selected_assessment(results, out_dir, true)
}

fn selected_samples(run: &MonteCarloResult, deadline: bool) -> &[RecoverySample] {
    if deadline {
        &run.recovery_at_deadline
    } else {
        &run.recovery
    }
}

fn write_selected_assessment(
    results: &[MonteCarloResult],
    out_dir: &str,
    deadline: bool,
) -> anyhow::Result<bool> {
    anyhow::ensure!(!results.is_empty(), "empty campaign");
    anyhow::ensure!(
        results
            .iter()
            .map(|r| r.seed)
            .collect::<std::collections::HashSet<_>>()
            .len()
            == results.len(),
        "duplicate Monte Carlo seeds"
    );
    let samples = selected_samples(&results[0], deadline).len();
    anyhow::ensure!(samples > 0, "empty recovery window");
    anyhow::ensure!(
        results
            .iter()
            .all(|r| selected_samples(r, deadline).len() == samples),
        "mismatched recovery grids"
    );
    let mut curve = csv::Writer::from_path(format!("{out_dir}/recovery_ensemble.csv"))?;
    curve.write_record([
        "t",
        "quantity",
        "empirical_3rms",
        "reported_3sigma_rms",
        "limit",
    ])?;
    let mut maxima = [0.0_f64; 6];
    for k in 0..samples {
        let t = selected_samples(&results[0], deadline)[k].t;
        anyhow::ensure!(
            results
                .iter()
                .all(|r| (selected_samples(r, deadline)[k].t - t).abs() < 1e-9),
            "mismatched timestamps"
        );
        for i in 0..6 {
            let mse = results
                .iter()
                .map(|r| selected_samples(r, deadline)[k].squared_error[i])
                .sum::<f64>()
                / results.len() as f64;
            let var = results
                .iter()
                .map(|r| selected_samples(r, deadline)[k].variance[i])
                .sum::<f64>()
                / results.len() as f64;
            let empirical = 3.0 * mse.sqrt();
            maxima[i] = if empirical.is_finite() {
                maxima[i].max(empirical)
            } else {
                f64::INFINITY
            };
            curve.serialize((t, NAMES[i], empirical, 3.0 * var.sqrt(), LIMITS[i]))?;
        }
    }
    curve.flush()?;
    let total = results.len() * samples;
    let covered = results
        .iter()
        .flat_map(|r| selected_samples(r, deadline))
        .filter(|s| s.nees >= CHI18[0] && s.nees <= CHI18[1])
        .count() as f64
        / total as f64;
    let mean_nees = results
        .iter()
        .flat_map(|r| selected_samples(r, deadline))
        .map(|s| s.nees)
        .sum::<f64>()
        / total as f64;
    let nis_count = results.iter().map(|r| r.nis_count).sum::<usize>();
    let nis_covered = results.iter().map(|r| r.nis_hits).sum::<usize>() as f64 / nis_count as f64;
    let mut auxiliary_coverage = [0.0; 2];
    for (i, coverage) in auxiliary_coverage.iter_mut().enumerate() {
        let count = results.iter().map(|r| r.aux_nis[i].count).sum::<usize>();
        let hits = results.iter().map(|r| r.aux_nis[i].hits).sum::<usize>();
        *coverage = hits as f64 / count as f64;
    }
    let run_passes = |r: &MonteCarloResult| {
        if deadline {
            let values = selected_samples(r, true);
            r.recovery_numerically_valid
                && r.recovery_window_complete
                && values.len() == 201
                && values
                    .iter()
                    .enumerate()
                    .all(|(i, s)| (s.t - r.recovery[0].t - i as f64 * 0.01).abs() < 1e-8)
                && values.iter().all(|s| {
                    s.variance
                        .iter()
                        .zip(LIMITS)
                        .all(|(&v, limit)| v.is_finite() && v >= 0.0 && 3.0 * v.sqrt() <= limit)
                })
        } else {
            r.passed
        }
    };
    let covariance_passed = results.iter().all(run_passes);
    let mut per_run = csv::Writer::from_path(format!("{out_dir}/per_run.csv"))?;
    per_run.write_record([
        "seed",
        "samples",
        "covariance_passed",
        "position_m",
        "velocity_mps",
        "attitude_deg",
        "rate_dps",
        "gyro_bias_dps",
        "accel_bias_ug",
        "nees18_sum",
        "nees18_hits",
        "gnss_nis_samples",
        "gnss_nis_hits",
        "star_nis_samples",
        "star_nis_hits",
        "radar_nis_samples",
        "radar_nis_hits",
    ])?;
    for run in results {
        let mut max_sigma = [0.0_f64; 6];
        for sample in selected_samples(run, deadline) {
            for i in 0..6 {
                let variance = sample.variance[i];
                max_sigma[i] = max_sigma[i].max(if variance.is_finite() && variance >= 0.0 {
                    3.0 * variance.sqrt()
                } else {
                    f64::INFINITY
                });
            }
        }
        let mut row = vec![
            run.seed.to_string(),
            selected_samples(run, deadline).len().to_string(),
            run_passes(run).to_string(),
        ];
        row.extend(max_sigma.map(|v| v.to_string()));
        row.push(
            selected_samples(run, deadline)
                .iter()
                .map(|s| s.nees)
                .sum::<f64>()
                .to_string(),
        );
        row.push(
            selected_samples(run, deadline)
                .iter()
                .filter(|s| s.nees >= CHI18[0] && s.nees <= CHI18[1])
                .count()
                .to_string(),
        );
        row.extend(
            [
                run.nis_count,
                run.nis_hits,
                run.aux_nis[0].count,
                run.aux_nis[0].hits,
                run.aux_nis[1].count,
                run.aux_nis[1].hits,
            ]
            .map(|n| n.to_string()),
        );
        per_run.write_record(row)?;
    }
    per_run.flush()?;
    let empirical_passed = maxima
        .iter()
        .zip(LIMITS)
        .all(|(&v, limit)| v.is_finite() && v <= limit);
    let passed = results.len() >= 100
        && covariance_passed
        && empirical_passed
        && covered >= 0.9
        && nis_covered >= 0.9
        && mean_nees <= 18.0 * 1.3
        && auxiliary_coverage.iter().all(|&v| v >= 0.9);
    let mut summary = std::fs::File::create(format!("{out_dir}/campaign.txt"))?;
    writeln!(
        summary,
        "assessment={}",
        if deadline {
            "retrospective checkpoints, data cutoff at exit+2s, 100Hz"
        } else {
            "causal, 400Hz"
        }
    )?;
    writeln!(summary, "Numerical reference campaign (not vehicle or ARM certification)\nruns={} recovery_samples_per_run={}\ncovariance_limits_passed={} empirical_limits_passed={}\nrecovery_nees18_coverage={:.6} recovery_nees18_mean={:.6} normalized_nees={:.6}\nbase_gnss_nis_coverage={:.6} nis_samples={}\nstatistical_screen_passed={}",results.len(),samples,covariance_passed,empirical_passed,covered,mean_nees,mean_nees/18.0,nis_covered,nis_count,passed)?;
    writeln!(
        summary,
        "base_star_nis_coverage={:.6} base_radar_nis_coverage={:.6}",
        auxiliary_coverage[0], auxiliary_coverage[1]
    )?;
    for i in 0..6 {
        let block_mean = results
            .iter()
            .flat_map(|r| selected_samples(r, deadline))
            .map(|s| s.block_nees[i])
            .sum::<f64>()
            / total as f64;
        let max_reported = results
            .iter()
            .flat_map(|r| selected_samples(r, deadline))
            .map(|s| {
                if s.variance[i].is_finite() && s.variance[i] >= 0.0 {
                    3.0 * s.variance[i].sqrt()
                } else {
                    f64::INFINITY
                }
            })
            .fold(0.0_f64, f64::max);
        writeln!(
            summary,
            "{}: max_empirical_3rms={:.9} max_reported_3sigma={:.9} limit={} mean_block_nees={:.6}",
            NAMES[i], maxima[i], max_reported, LIMITS[i], block_mean
        )?;
    }
    println!(
        "Assessment: {}",
        if deadline {
            "retrospective at exit+2s (100Hz checkpoints)"
        } else {
            "causal (400Hz)"
        }
    );
    println!("Empirical recovery maxima (3×ensemble RMS): {maxima:?}");
    println!("Recovery NEES(18) coverage={covered:.4}, mean={mean_nees:.3}; unconditioned base GNSS NIS coverage={nis_covered:.4}; statistical screen={passed}");
    println!(
        "Unconditioned star NIS coverage={:.4}, radar={:.4}",
        auxiliary_coverage[0], auxiliary_coverage[1]
    );
    Ok(passed)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_nees_includes_cross_covariance() {
        let mut truth = State::new();
        let estimate = State::new();
        truth.r.x = 1.0;
        truth.v.x = 1.0;
        let mut p = MatN::identity();
        p[(0, 3)] = 0.5;
        p[(3, 0)] = 0.5;
        let sample = RecoverySample::new(&truth, &estimate, &p);
        assert!((sample.nees - 4.0 / 3.0).abs() < 1e-12);
        assert_eq!(sample.squared_error[0], 1.0);
        assert_eq!(sample.block_nees[0] + sample.block_nees[1], 2.0);
    }
}
