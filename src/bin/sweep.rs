use entry_blackout_estimation::{
    data::{AeroTable, RcsTorqueSeries},
    sim::{run_monte_carlo, RobustTuning, SimConfig},
    types::{DynamicsParams, NoiseParams},
};

fn main() -> anyhow::Result<()> {
    let aero = AeroTable::load("data/Aero Coeff Table.csv")?;
    let rcs = RcsTorqueSeries::load("data/RCS Torque Series 100Hz Preview.csv")?;
    let mut configs = Vec::new();
    let gates = [0.9990, 0.9995];
    let caps = [1.0e3, 1.0e4];
    let student_nus = [0.0, 6.0];
    let iags: [Option<f64>; 2] = [None, Some(0.35)];
    let reacq: [(f64, f64); 2] = [(30.0, 10.0), (50.0, 15.0)];
    for &gp in &gates {
        for &cap in &caps {
            for &student_nu in &student_nus {
                for &ig in &iags {
                    for &(ps, vs) in &reacq {
                        let mut cfg = SimConfig::default();
                        cfg.tuning = Some(RobustTuning {
                            gate_p: gp,
                            r_scale_cap: cap,
                            student_nu,
                            iae_gamma: ig,
                            reacq_pos_scale: ps,
                            reacq_vel_scale: vs,
                        });
                        configs.push(cfg);
                    }
                }
            }
        }
    }
    println!("Sweeping {} configs...", configs.len());
    for (i, cfg) in configs.into_iter().enumerate() {
        let runs = 16usize;
        let noise = NoiseParams::tuned();
        let results = run_monte_carlo(
            runs,
            DynamicsParams::default(),
            noise,
            aero.clone(),
            rcs.clone(),
            cfg.clone(),
        );
        let mut rms_pos = 0.0;
        let mut rms_vel = 0.0;
        let mut rms_att = 0.0;
        let mut nees = 0.0;
        let mut nis = 0.0;
        let mut nees_cov = 0.0;
        let mut nis_cov = 0.0;
        let mut pass = 0usize;
        for r in &results {
            rms_pos += r.rms_pos;
            rms_vel += r.rms_vel;
            rms_att += r.rms_att_deg;
            nees += r.nees_mean;
            nis += r.nis_mean;
            nees_cov += r.nees_covered;
            nis_cov += r.nis_covered;
            if r.passed {
                pass += 1;
            }
        }
        let n = results.len() as f64;
        let rp = rms_pos / n;
        let rv = rms_vel / n;
        let ra = rms_att / n;
        let nm = nis / n;
        let ne = nees / n;
        let nc = 100.0 * nis_cov / n;
        let ec = 100.0 * nees_cov / n;
        let pr = 100.0 * (pass as f64) / n;
        println!(
            "cfg#{:02}: gate_p={:.5} cap={:.1e} nu={:.1} iae={:?} reacq=({:.0},{:.0}) | pass={:.1}% rms(m)={:.2} rmsv={:.3} att(deg)={:.3} | NIS={:.2}({:.1}%) NEES={:.2}({:.1}%)",
            i, cfg.tuning.unwrap().gate_p, cfg.tuning.unwrap().r_scale_cap, cfg.tuning.unwrap().student_nu, cfg.tuning.unwrap().iae_gamma, cfg.tuning.unwrap().reacq_pos_scale, cfg.tuning.unwrap().reacq_vel_scale,
            pr, rp, rv, ra, nm, nc, ne, ec
        );
    }
    Ok(())
}
