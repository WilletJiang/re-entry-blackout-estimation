use entry_blackout_estimation::{
    data::{AeroTable, RcsTorqueSeries},
    sim::{run_monte_carlo, SimConfig},
    types::{DynamicsParams, NoiseParams},
};

fn main() -> anyhow::Result<()> {
    #[cfg(debug_assertions)]
    eprintln!("WARNING: Running bench in debug mode; use `--release` for accurate performance.");
    let aero = AeroTable::load("data/Aero Coeff Table.csv")?;
    let rcs = RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv")?;
    let noise = NoiseParams::tuned();
    let cfg = SimConfig::default();
    let results = run_monte_carlo(
        100,
        DynamicsParams::default(),
        noise,
        aero,
        rcs,
        cfg,
    );
    let mut rms_pos = 0.0;
    let mut rms_vel = 0.0;
    let mut rms_att = 0.0;
    let mut nees = 0.0;
    let mut nis = 0.0;
    let mut nees_cov = 0.0;
    let mut nis_cov = 0.0;
    for r in &results {
        rms_pos += r.rms_pos;
        rms_vel += r.rms_vel;
        rms_att += r.rms_att_deg;
        nees += r.nees_mean;
        nis += r.nis_mean;
        nees_cov += r.nees_covered;
        nis_cov += r.nis_covered;
    }
    let n = results.len() as f64;
    println!("MC runs: {}", results.len());
    println!("RMS pos (m): {:.2}", rms_pos / n);
    println!("RMS vel (m/s): {:.3}", rms_vel / n);
    println!("RMS att (deg): {:.4}", rms_att / n);
    println!("NEES mean: {:.3}, coverage: {:.2}%", nees / n, 100.0 * nees_cov / n);
    println!("NIS  mean: {:.3}, coverage: {:.2}%", nis / n, 100.0 * nis_cov / n);
    Ok(())
}
