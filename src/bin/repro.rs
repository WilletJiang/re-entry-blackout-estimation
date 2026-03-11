use entry_blackout_estimation::{
    data::{load_vehicle_props, AeroTable, RcsTorqueSeries},
    sim::{simulate_single, SimConfig},
    types::NoiseParams,
};
use std::env;
use std::sync::Arc;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    let seed: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or_else(|| {
        eprintln!("usage: repro <seed>");
        std::process::exit(2);
    });

    let dyn_params = load_vehicle_props("data/Vehicle Properties Preview.csv")?;
    let aero = Arc::new(AeroTable::load("data/Aero Coeff Table.csv")?);
    let rcs = Arc::new(RcsTorqueSeries::load(
        "data/RCS Torque Series 100Hz Preview.csv",
    )?);

    let mut cfg = SimConfig::default();
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
        "passed={} rms_pos_m={:.3} rms_vel_mps={:.3} rms_att_deg={:.4} ba_rms_ug={:.1} bg_rms_dps={:.6} nees_mean={:.3} nis_mean={:.3}",
        r.passed as u8,
        r.rms_pos,
        r.rms_vel,
        r.rms_att_deg,
        r.bias_ba_rms_ug,
        r.bias_bg_rms_dps,
        r.nees_mean,
        r.nis_mean
    );
    Ok(())
}
