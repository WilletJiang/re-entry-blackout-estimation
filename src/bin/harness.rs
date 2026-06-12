use entry_blackout_estimation::{
    data::load_vehicle_props,
    data::{AeroTable, RcsTorqueSeries},
    plot::{
        plot_bias_histograms, plot_blackout_envelope, plot_blackout_nees_time, plot_histogram,
        plot_nees_nis, plot_pass_rate, write_blackout_envelope_csv, write_nees_blocks_csv,
        write_nees_phase_csv, write_post_window_csv, write_summary_csv,
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
    let aero = AeroTable::load("data/Aero Coeff Table.csv")?;
    let rcs = RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv")?;
    let mut cfg = SimConfig::default();
    cfg.use_tight_coupling = std::env::var("TIGHT").ok().is_some_and(|v| v != "0");
    if std::env::var("BIAS_JUMP").ok().is_some_and(|v| v != "0") {
        cfg.bias_jump = Some(BiasJumpConfig {
            t: 0.5 * (cfg.blackout_start + cfg.blackout_end),
            bg_step_dps: 0.005,
            ba_step_ug: 80.0,
        });
    }
    let mode = if cfg.use_tight_coupling {
        "tight"
    } else {
        "pvt"
    };
    let noise = NoiseParams::tuned();

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

    // CSV (extended summary with bias and pass flag)
    write_summary_csv(&results, &format!("{out_dir}/mc_summary.csv"))?;
    write_nees_blocks_csv(&results, &format!("{out_dir}/mc_nees_blocks.csv"))?;
    write_nees_phase_csv(&results, &format!("{out_dir}/mc_nees_phase.csv"))?;
    write_post_window_csv(&results, &format!("{out_dir}/mc_post_window.csv"))?;
    write_summary_csv(&results, &format!("{out_dir}/mc_summary_{mode}.csv"))?;
    write_nees_blocks_csv(&results, &format!("{out_dir}/mc_nees_blocks_{mode}.csv"))?;
    write_nees_phase_csv(&results, &format!("{out_dir}/mc_nees_phase_{mode}.csv"))?;
    write_post_window_csv(&results, &format!("{out_dir}/mc_post_window_{mode}.csv"))?;

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
    println!("Pass rate: {:.1}% ({} / {})", rate, passed, results.len());
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
