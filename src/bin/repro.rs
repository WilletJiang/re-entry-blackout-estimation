use entry_blackout_estimation::{
    data::{load_vehicle_props, AeroTable, RcsTorqueSeries},
    sim::{simulate_single, SimConfig},
    types::NoiseParams,
};
use std::sync::Arc;

const BYTES_PER_MIB: f64 = 1_048_576.0;

fn main() -> anyhow::Result<()> {
    let seed = parse_seed_arg()?;

    let dyn_params = load_vehicle_props("data/Vehicle Properties Preview.csv")?;
    let mut aero = Arc::new(AeroTable::load("data/Aero Coeff Table.csv")?);
    let rcs = Arc::new(RcsTorqueSeries::load(
        "data/RCS Torque Series 100Hz Preview.csv",
    )?);

    let mut cfg = SimConfig::default();
    let scenario = std::env::var("SCENARIO").unwrap_or_else(|_| "reference-entry".to_string());
    if scenario == "reference-entry" {
        let (initial, reference_aero) =
            entry_blackout_estimation::sim::scenario::reference_entry(&dyn_params);
        cfg.initial_truth = Some(initial);
        cfg.initial_bias_jitter = Some((0.3_f64.to_radians() / 3600.0, 150e-6 * 9.80665));
        aero = Arc::new(reference_aero);
    } else if scenario != "legacy" {
        anyhow::bail!("SCENARIO must be reference-entry or legacy");
    }
    eprintln!(
        "scenario={scenario}; reference-entry uses an explicitly synthetic aerodynamic model"
    );
    cfg.use_tight_coupling = std::env::var("TIGHT").ok().is_some_and(|v| v != "0");
    let noise = NoiseParams::tuned();

    let mode = if cfg.use_tight_coupling {
        "tight"
    } else {
        "pvt"
    };
    eprintln!("mode={mode} seed={seed}");
    let r = simulate_single(seed, dyn_params, noise, aero, rcs, cfg);
    println!(
        "passed={} rms_pos_m={:.3} rms_vel_mps={:.3} rms_att_deg={:.4} sigma3_pos_max_m={:.3} sigma3_vel_max_mps={:.3} sigma3_att_max_deg={:.4} ba_rms_ug={:.1} bg_rms_dps={:.6} sigma3_ba_max_ug={:.1} sigma3_bg_max_dps={:.6} hist_mib={:.3} nees_mean={:.3} nis_mean={:.3}",
        if r.passed { 1 } else { 0 },
        r.rms_pos,
        r.rms_vel,
        r.rms_att_deg,
        r.sigma3_pos_max_m,
        r.sigma3_vel_max_mps,
        r.sigma3_att_max_deg,
        r.bias_ba_rms_ug,
        r.bias_bg_rms_dps,
        r.sigma3_ba_max_ug,
        r.sigma3_bg_max_dps,
        usize_to_f64(r.history_capacity_bytes, "history_capacity_bytes")? / BYTES_PER_MIB,
        r.nees_mean,
        r.nis_mean
    );
    let mut smoothed_max = [0.0_f64; 6];
    for sample in &r.recovery_at_deadline {
        for i in 0..6 {
            smoothed_max[i] = smoothed_max[i].max(3.0 * sample.variance[i].sqrt());
        }
    }
    println!(
        "retrospective at exit+2s: samples={} max_3sigma={smoothed_max:?}",
        r.recovery_at_deadline.len()
    );
    println!("angular_rate: rms_dps={:.6} sigma3_max_dps={:.6}; recovery: samples={} complete={} finite={}",
        r.rms_omega_dps, r.sigma3_omega_max_dps, r.recovery_samples,
        r.recovery_window_complete, r.recovery_numerically_valid);
    println!("truth (seed schema 2): initial_radial_mps={:.3} altitude_m={:.1}->{:.1} Mach={:.2}..{:.2} outside_aero_table={}/{} radar_eligible={}",
        r.trajectory.initial_radial_velocity_mps, r.trajectory.initial_altitude_m,
        r.trajectory.final_altitude_m, r.trajectory.min_mach, r.trajectory.max_mach,
        r.trajectory.outside_aero_table_samples, r.trajectory.samples, r.trajectory.radar_eligible_samples);
    Ok(())
}

fn parse_seed_arg() -> anyhow::Result<u64> {
    let Some(raw_seed) = std::env::args().nth(1) else {
        return Err(anyhow::anyhow!(
            "missing seed argument; usage: repro <seed>"
        ));
    };
    raw_seed
        .parse::<u64>()
        .map_err(|err| anyhow::anyhow!("invalid seed argument '{raw_seed}': {err}"))
}

fn usize_to_f64(value: usize, field_name: &str) -> anyhow::Result<f64> {
    let value_u32 = u32::try_from(value).map_err(|_| {
        anyhow::anyhow!(
            "{field_name}={value} exceeds u32 range needed for loss-aware f64 conversion"
        )
    })?;
    Ok(f64::from(value_u32))
}
