use entry_blackout_estimation::{
    data::load_vehicle_props,
    data::{AeroTable, RcsTorqueSeries},
    plot::{
        plot_bias_histograms, plot_blackout_envelope, plot_blackout_nees_time, plot_histogram,
        plot_nees_nis, plot_pass_rate, write_blackout_envelope_csv, write_nees_blocks_csv,
        write_nees_phase_csv, write_post_window_csv, write_summary_csv, write_trajectory_csv,
    },
    sim::{
        run_monte_carlo, run_monte_carlo_with_trace, BiasJumpConfig, BlackoutTraceConfig, SimConfig,
    },
    types::NoiseParams,
};
use std::env;

const DEFAULT_RUNS: usize = 100;
const BYTES_PER_MIB: f64 = 1_048_576.0;

fn main() -> anyhow::Result<()> {
    let runs = parse_runs_arg()?;

    let out_dir = env::var("OUT_DIR").unwrap_or_else(|_| "output".to_string());
    std::fs::create_dir_all(&out_dir)?;

    let dyn_params = load_vehicle_props("data/Vehicle Properties Preview.csv")?;
    let mut aero = AeroTable::load("data/Aero Coeff Table.csv")?;
    let rcs = RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv")?;
    let mut cfg = SimConfig::default();
    let scenario = std::env::var("SCENARIO").unwrap_or_else(|_| "reference-entry".to_string());
    if scenario == "reference-entry" {
        let (initial, reference_aero) =
            entry_blackout_estimation::sim::scenario::reference_entry(&dyn_params);
        cfg.initial_truth = Some(initial);
        cfg.initial_bias_jitter = Some((0.3_f64.to_radians() / 3600.0, 150e-6 * 9.80665));
        aero = reference_aero;
    } else if scenario != "legacy" {
        anyhow::bail!("SCENARIO must be reference-entry or legacy");
    }
    eprintln!(
        "scenario={scenario}; reference-entry uses an explicitly synthetic aerodynamic model"
    );
    cfg.use_tight_coupling = std::env::var("TIGHT").ok().is_some_and(|v| v != "0");
    if std::env::var("BIAS_JUMP").ok().is_some_and(|v| v != "0") {
        cfg.bias_jump = Some(BiasJumpConfig {
            t: 0.5 * (cfg.blackout_start + cfg.blackout_end),
            bg_step_dps: 0.005,
            ba_step_ug: 80.0,
        });
    }
    if let Ok(offset) = env::var("SEED_OFFSET") {
        cfg.seed_offset = offset.parse()?;
    }
    let mode = if cfg.use_tight_coupling {
        "tight"
    } else {
        "pvt"
    };
    let noise = NoiseParams::tuned();
    let assessment = env::var("ASSESSMENT").unwrap_or_else(|_| "causal".to_string());
    anyhow::ensure!(
        matches!(assessment.as_str(), "causal" | "deadline"),
        "ASSESSMENT must be causal or deadline"
    );
    let revision = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    let dirty = std::process::Command::new("git")
        .args([
            "diff",
            "--quiet",
            "HEAD",
            "--",
            "src",
            "Cargo.toml",
            "Cargo.lock",
        ])
        .status()
        .map_or(true, |s| !s.success());
    std::fs::write(format!("{out_dir}/configuration.txt"), format!(
        "source_revision={revision}\nsource_code_dirty={dirty}\nscenario={scenario}\nmode={mode}\nseed_schema=2\nseed_offset={}\nruns={runs}\ndt_s={}\nt_final_s={}\nhistory_s={}\noutlier_rate={}\ndropout_rate={}\nassessment={assessment}\narch={}\nstar_convergence_count={}\nfilter_noise={noise:?}\n",
        cfg.seed_offset, cfg.dt, cfg.t_final, cfg.history_sec, cfg.outlier_rate, cfg.dropout_rate, std::env::consts::ARCH, cfg.star_convergence_count))?;

    println!(
        "Running {runs} Monte Carlo cases (mode={mode}, history={:.1}s)...",
        cfg.history_sec
    );
    let record_trace = std::env::var("BLACKOUT_TRACE")
        .ok()
        .is_some_and(|v| v != "0");
    let (results, traces, trace_cfg) = if record_trace {
        let trace_cfg = BlackoutTraceConfig {
            t_start: (cfg.blackout_start - 10.0).max(0.0),
            t_end: (cfg.blackout_end + 10.0).min(cfg.t_final),
            dt: 0.1,
            blackout_start: cfg.blackout_start,
            blackout_end: cfg.blackout_end,
        };
        let out = run_monte_carlo_with_trace(runs, dyn_params, noise, aero, rcs, cfg, trace_cfg);
        let (results, traces): (Vec<_>, Vec<_>) =
            out.into_iter().map(|o| (o.result, o.trace)).unzip();
        (results, Some(traces), Some(trace_cfg))
    } else {
        let results = run_monte_carlo(runs, dyn_params, noise, aero, rcs, cfg);
        (results, None, None)
    };
    println!("Done. Writing summary...");
    let causal_passed =
        entry_blackout_estimation::sim::campaign::write_assessment(&results, &out_dir)?;
    let deadline_dir = format!("{out_dir}/deadline");
    std::fs::create_dir_all(&deadline_dir)?;
    let deadline_passed = entry_blackout_estimation::sim::campaign::write_deadline_assessment(
        &results,
        &deadline_dir,
    )?;
    let statistical_passed = if assessment == "deadline" {
        deadline_passed
    } else {
        causal_passed
    };

    // CSV (extended summary with bias and pass flag)
    write_summary_csv(&results, &format!("{out_dir}/mc_summary.csv"))?;
    write_nees_blocks_csv(&results, &format!("{out_dir}/mc_nees_blocks.csv"))?;
    write_nees_phase_csv(&results, &format!("{out_dir}/mc_nees_phase.csv"))?;
    write_post_window_csv(&results, &format!("{out_dir}/mc_post_window.csv"))?;
    write_summary_csv(&results, &format!("{out_dir}/mc_summary_{mode}.csv"))?;
    write_nees_blocks_csv(&results, &format!("{out_dir}/mc_nees_blocks_{mode}.csv"))?;
    write_nees_phase_csv(&results, &format!("{out_dir}/mc_nees_phase_{mode}.csv"))?;
    write_post_window_csv(&results, &format!("{out_dir}/mc_post_window_{mode}.csv"))?;

    write_trajectory_csv(&results, &format!("{out_dir}/mc_trajectory_{mode}.csv"))?;

    // Plot NEES/NIS coverage
    plot_nees_nis(&results, &format!("{out_dir}/mc_nees_nis.svg"))?;
    plot_histogram(&results, &format!("{out_dir}/mc_rmspos_hist.svg"))?;
    // New plots
    plot_bias_histograms(
        &results,
        &format!("{out_dir}/mc_gyro_bias_hist.svg"),
        &format!("{out_dir}/mc_accel_bias_hist.svg"),
    )?;
    plot_pass_rate(&results, &format!("{out_dir}/mc_pass_rate.svg"))?;
    plot_nees_nis(&results, &format!("{out_dir}/mc_nees_nis_{mode}.svg"))?;
    plot_histogram(&results, &format!("{out_dir}/mc_rmspos_hist_{mode}.svg"))?;
    plot_bias_histograms(
        &results,
        &format!("{out_dir}/mc_gyro_bias_hist_{mode}.svg"),
        &format!("{out_dir}/mc_accel_bias_hist_{mode}.svg"),
    )?;
    plot_pass_rate(&results, &format!("{out_dir}/mc_pass_rate_{mode}.svg"))?;

    if let (Some(traces), Some(trace_cfg)) = (traces.as_ref(), trace_cfg) {
        write_blackout_envelope_csv(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_envelope.csv"),
        )?;
        write_blackout_envelope_csv(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_envelope_{mode}.csv"),
        )?;
        plot_blackout_envelope(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_envelope.svg"),
        )?;
        plot_blackout_nees_time(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_nees_time.svg"),
        )?;
        plot_blackout_envelope(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_envelope_{mode}.svg"),
        )?;
        plot_blackout_nees_time(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_nees_time_{mode}.svg"),
        )?;
    }

    // Pass rate summary
    let passed = results.iter().filter(|r| r.passed).count();
    let rate = if results.is_empty() {
        0.0
    } else {
        100.0 * usize_to_f64(passed, "passed_runs")? / usize_to_f64(results.len(), "total_runs")?
    };
    let history_mean_bytes = if results.is_empty() {
        0.0
    } else {
        let history_sum_bytes = results.iter().try_fold(0.0, |acc, run_result| {
            Ok::<f64, anyhow::Error>(
                acc + usize_to_f64(run_result.history_capacity_bytes, "history_capacity_bytes")?,
            )
        })?;
        history_sum_bytes / usize_to_f64(results.len(), "total_runs")?
    };
    let history_max_bytes = results
        .iter()
        .map(|r| r.history_capacity_bytes)
        .max()
        .unwrap_or(0);
    println!(
        "Causal recovery covariance screen: {:.1}% ({} / {})",
        rate,
        passed,
        results.len()
    );
    println!("This is not full task acceptance; consistency, trajectory validity, and target runtime require separate evidence.");
    let ascending = results
        .iter()
        .filter(|r| r.trajectory.initial_radial_velocity_mps >= 0.0)
        .count();
    let no_radar = results
        .iter()
        .filter(|r| r.trajectory.radar_eligible_samples == 0)
        .count();
    println!("Truth audit: {ascending} initially ascending runs, {no_radar} runs never below the radar altitude ceiling.");
    println!("Truth diagnostics: {out_dir}/mc_trajectory_{mode}.csv (seed schema 2)");
    println!(
        "History capacity: mean {:.3} MiB, max {:.3} MiB",
        history_mean_bytes / BYTES_PER_MIB,
        usize_to_f64(history_max_bytes, "history_capacity_bytes")? / BYTES_PER_MIB
    );

    println!("Outputs: {out_dir}/mc_summary.csv, {out_dir}/mc_nees_blocks.csv, {out_dir}/mc_nees_phase.csv, {out_dir}/mc_post_window.csv, {out_dir}/mc_nees_nis.svg, {out_dir}/mc_rmspos_hist.svg, {out_dir}/mc_gyro_bias_hist.svg, {out_dir}/mc_accel_bias_hist.svg, {out_dir}/mc_pass_rate.svg");
    println!("Outputs ({mode}): {out_dir}/mc_summary_{mode}.csv, {out_dir}/mc_nees_blocks_{mode}.csv, {out_dir}/mc_nees_phase_{mode}.csv, {out_dir}/mc_post_window_{mode}.csv, {out_dir}/mc_nees_nis_{mode}.svg, {out_dir}/mc_rmspos_hist_{mode}.svg, {out_dir}/mc_gyro_bias_hist_{mode}.svg, {out_dir}/mc_accel_bias_hist_{mode}.svg, {out_dir}/mc_pass_rate_{mode}.svg");
    if record_trace {
        println!("Outputs (trace): {out_dir}/mc_blackout_envelope.csv, {out_dir}/mc_blackout_nees_time.svg");
        println!("Outputs (trace/{mode}): {out_dir}/mc_blackout_envelope_{mode}.csv, {out_dir}/mc_blackout_nees_time_{mode}.svg, {out_dir}/mc_blackout_envelope_{mode}.svg");
    }
    if env::args().any(|arg| arg == "--check") {
        anyhow::ensure!(
            statistical_passed,
            "{assessment} statistical reference screen failed; see {out_dir}/campaign.txt and {out_dir}/deadline/campaign.txt"
        );
    }
    Ok(())
}

fn parse_runs_arg() -> anyhow::Result<usize> {
    let Some(raw_runs) = env::args().nth(1) else {
        return Ok(DEFAULT_RUNS);
    };
    raw_runs
        .parse::<usize>()
        .map_err(|err| anyhow::anyhow!("invalid runs argument '{raw_runs}': {err}"))
}

fn usize_to_f64(value: usize, field_name: &str) -> anyhow::Result<f64> {
    let value_u32 = u32::try_from(value).map_err(|_| {
        anyhow::anyhow!(
            "{field_name}={value} exceeds u32 range needed for loss-aware f64 conversion"
        )
    })?;
    Ok(f64::from(value_u32))
}
