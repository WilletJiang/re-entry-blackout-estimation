use entry_blackout_estimation::{
    data::{AeroTable, RcsTorqueSeries},
    plot::{plot_histogram, plot_nees_nis, plot_bias_histograms, plot_pass_rate, plot_blackout_envelope, plot_blackout_nees_time, write_blackout_envelope_csv, write_nees_blocks_csv, write_nees_phase_csv, write_post_window_csv, write_summary_csv},
    sim::{run_monte_carlo, run_monte_carlo_with_trace, BiasJumpConfig, BlackoutTraceConfig, SimConfig},
    types::NoiseParams,
    data::load_vehicle_props,
};
use std::env;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    let runs: usize = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(100);

    let out_dir = env::var("OUT_DIR").unwrap_or_else(|_| "output".to_string());
    std::fs::create_dir_all(&out_dir)?;

    let dyn_params = load_vehicle_props("data/Vehicle Properties Preview.csv")?;
    let aero = AeroTable::load("data/Aero Coeff Table.csv")?;
    let rcs = RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv")?;
    let mut cfg = SimConfig::default();
    cfg.use_tight_coupling = std::env::var("TIGHT")
        .ok()
        .is_some_and(|v| v != "0");
    if std::env::var("BIAS_JUMP").ok().is_some_and(|v| v != "0") {
        cfg.bias_jump = Some(BiasJumpConfig {
            t: 0.5 * (cfg.blackout_start + cfg.blackout_end),
            bg_step_dps: 0.005,
            ba_step_ug: 80.0,
        });
    }
    let mode = if cfg.use_tight_coupling { "tight" } else { "pvt" };
    let noise = NoiseParams::tuned();

    println!("Running {runs} Monte Carlo cases...");
    let record_trace = std::env::var("BLACKOUT_TRACE").ok().is_some_and(|v| v != "0");
    let (results, traces, trace_cfg) = if record_trace {
        let trace_cfg = BlackoutTraceConfig {
            t_start: (cfg.blackout_start - 10.0).max(0.0),
            t_end: (cfg.blackout_end + 10.0).min(cfg.t_final),
            dt: 0.1,
            blackout_start: cfg.blackout_start,
            blackout_end: cfg.blackout_end,
        };
        let out = run_monte_carlo_with_trace(runs, dyn_params, noise, aero, rcs, cfg, trace_cfg);
        let (results, traces): (Vec<_>, Vec<_>) = out.into_iter().map(|o| (o.result, o.trace)).unzip();
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
        write_blackout_envelope_csv(traces, &trace_cfg, &format!("{out_dir}/mc_blackout_envelope.csv"))?;
        write_blackout_envelope_csv(
            traces,
            &trace_cfg,
            &format!("{out_dir}/mc_blackout_envelope_{mode}.csv"),
        )?;
        plot_blackout_envelope(traces, &trace_cfg, &format!("{out_dir}/mc_blackout_envelope.svg"))?;
        plot_blackout_nees_time(traces, &trace_cfg, &format!("{out_dir}/mc_blackout_nees_time.svg"))?;
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
    let rate = (passed as f64) / (results.len() as f64) * 100.0;
    println!("Pass rate: {:.1}% ({} / {})", rate, passed, results.len());

    println!("Outputs: {out_dir}/mc_summary.csv, {out_dir}/mc_nees_blocks.csv, {out_dir}/mc_nees_phase.csv, {out_dir}/mc_post_window.csv, {out_dir}/mc_nees_nis.svg, {out_dir}/mc_rmspos_hist.svg, {out_dir}/mc_gyro_bias_hist.svg, {out_dir}/mc_accel_bias_hist.svg, {out_dir}/mc_pass_rate.svg");
    println!("Outputs ({mode}): {out_dir}/mc_summary_{mode}.csv, {out_dir}/mc_nees_blocks_{mode}.csv, {out_dir}/mc_nees_phase_{mode}.csv, {out_dir}/mc_post_window_{mode}.csv, {out_dir}/mc_nees_nis_{mode}.svg, {out_dir}/mc_rmspos_hist_{mode}.svg, {out_dir}/mc_gyro_bias_hist_{mode}.svg, {out_dir}/mc_accel_bias_hist_{mode}.svg, {out_dir}/mc_pass_rate_{mode}.svg");
    if record_trace {
        println!("Outputs (trace): {out_dir}/mc_blackout_envelope.csv, {out_dir}/mc_blackout_nees_time.svg");
        println!("Outputs (trace/{mode}): {out_dir}/mc_blackout_envelope_{mode}.csv, {out_dir}/mc_blackout_nees_time_{mode}.svg, {out_dir}/mc_blackout_envelope_{mode}.svg");
    }
    Ok(())
}
