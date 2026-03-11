mod finish;
mod gnss;
mod init;
mod propagation;
mod radar;
mod star;
mod trace;

use super::{
    gps, metrics, BiasJumpConfig, BlackoutTrace, BlackoutTraceConfig, MonteCarloResult, SimConfig,
};
use crate::data::{AeroTable, RcsTorqueSeries};
use crate::estimator::Estimator;
use crate::physics::nominal_body_dynamics_with_residual;
use crate::types::{DynamicsParams, NoiseParams, State};
use nalgebra::{Matrix3, Vector3};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};
use std::sync::Arc;

const CHI2_6_LO: f64 = 1.237_344_245_791_202_5_f64;
const CHI2_6_HI: f64 = 14.449_375_335_447_918_f64;
const CHI2_3_LO: f64 = 0.215_795_282_623_897_94_f64;
const CHI2_3_HI: f64 = 9.348_403_604_496_145_f64;
const Z_95: f64 = 1.959_963_984_540_054_f64;

fn bias_jump_vectors(
    cfg: Option<BiasJumpConfig>,
    standard_normal: &Normal<f64>,
    rng: &mut StdRng,
) -> (Vector3<f64>, Vector3<f64>) {
    let Some(cfg) = cfg else {
        return (Vector3::zeros(), Vector3::zeros());
    };

    let mut dir = Vector3::new(
        standard_normal.sample(rng),
        standard_normal.sample(rng),
        standard_normal.sample(rng),
    );
    let n = dir.norm();
    if n.is_finite() && n > 0.0 {
        dir /= n;
    } else {
        dir = Vector3::new(1.0, 0.0, 0.0);
    }

    let bg_step = cfg.bg_step_dps.to_radians();
    let ba_step = cfg.ba_step_ug * 1.0e-6 * 9.80665;
    (dir * bg_step, dir * ba_step)
}

pub fn run_monte_carlo(
    runs: usize,
    dyn_params: DynamicsParams,
    noise: NoiseParams,
    aero: AeroTable,
    rcs: RcsTorqueSeries,
    cfg: SimConfig,
) -> Vec<MonteCarloResult> {
    use rayon::prelude::*;
    let aero = Arc::new(aero);
    let rcs = Arc::new(rcs);
    (0..runs)
        .into_par_iter()
        .map(|seed| {
            simulate_single(
                seed as u64,
                dyn_params.clone(),
                noise.clone(),
                aero.clone(),
                rcs.clone(),
                cfg.clone(),
            )
        })
        .collect()
}

pub struct MonteCarloTraceResult {
    pub result: MonteCarloResult,
    pub trace: BlackoutTrace,
}

pub fn run_monte_carlo_with_trace(
    runs: usize,
    dyn_params: DynamicsParams,
    noise: NoiseParams,
    aero: AeroTable,
    rcs: RcsTorqueSeries,
    cfg: SimConfig,
    trace_cfg: BlackoutTraceConfig,
) -> Vec<MonteCarloTraceResult> {
    use rayon::prelude::*;
    let aero = Arc::new(aero);
    let rcs = Arc::new(rcs);
    (0..runs)
        .into_par_iter()
        .map(|seed| {
            simulate_single_with_trace(
                seed as u64,
                dyn_params.clone(),
                noise.clone(),
                aero.clone(),
                rcs.clone(),
                cfg.clone(),
                trace_cfg,
            )
        })
        .collect()
}

pub fn simulate_single(
    seed: u64,
    dyn_params: DynamicsParams,
    noise: NoiseParams,
    aero: Arc<AeroTable>,
    rcs: Arc<RcsTorqueSeries>,
    cfg: SimConfig,
) -> MonteCarloResult {
    let mut run = SimRun::new(seed, dyn_params, noise, aero, rcs, cfg, None);
    run.run();
    run.finish().result
}

pub fn simulate_single_with_trace(
    seed: u64,
    dyn_params: DynamicsParams,
    noise: NoiseParams,
    aero: Arc<AeroTable>,
    rcs: Arc<RcsTorqueSeries>,
    cfg: SimConfig,
    trace_cfg: BlackoutTraceConfig,
) -> MonteCarloTraceResult {
    let mut run = SimRun::new(seed, dyn_params, noise, aero, rcs, cfg, Some(trace_cfg));
    run.run();
    let out = run.finish();
    MonteCarloTraceResult {
        result: out.result,
        trace: out.trace.expect("trace enabled"),
    }
}

struct SimRun {
    rng: StdRng,
    standard_normal: Normal<f64>,
    trace: bool,
    no_star: bool,
    dyn_params: DynamicsParams,
    noise: NoiseParams,
    aero: Arc<AeroTable>,
    rcs: Arc<RcsTorqueSeries>,
    cfg: SimConfig,
    est: Estimator,
    truth: State,
    next_gnss: f64,
    next_star: f64,
    next_radar: f64,
    gnss_reacq_done: bool,
    bias_cal_applied: bool,
    star_update_count: usize,
    nis_acc: f64,
    nis_ct: usize,
    nis_hit: usize,
    nis_phase_acc: [f64; 4],
    nis_phase_ct: [usize; 4],
    nis_phase_hit: [usize; 4],
    pos_err2_acc: f64,
    vel_err2_acc: f64,
    att_err2_acc: f64,
    err_samples: f64,
    bg_err2_acc: f64,
    ba_err2_acc: f64,
    win_pos_err2_acc: f64,
    win_vel_err2_acc: f64,
    win_bg_err2_acc: f64,
    win_ba_err2_acc: f64,
    win_samples: f64,
    nees_samples: Vec<metrics::NeesSample>,
    truth_times: Vec<f64>,
    truth_hist: Vec<State>,
    post_gnss_updates: usize,
    post_gnss_accepted: usize,
    post_sat_sum: usize,
    post_sat_ct: usize,
    post_sat_min: usize,
    bias_jump: Option<BiasJumpConfig>,
    bias_jump_applied: bool,
    bias_jump_bg: Vector3<f64>,
    bias_jump_ba: Vector3<f64>,
    phi_f_const: f64,
    phi_t_const: f64,
    inertia_inv: Matrix3<f64>,
    total_steps: usize,
    blackout_trace: Option<trace::TraceRecorder>,
}

impl SimRun {
    fn new(
        seed: u64,
        dyn_params: DynamicsParams,
        noise: NoiseParams,
        aero: Arc<AeroTable>,
        rcs: Arc<RcsTorqueSeries>,
        cfg: SimConfig,
        trace_cfg: Option<BlackoutTraceConfig>,
    ) -> Self {
        let mut rng = StdRng::seed_from_u64(0xBADA55 + seed);
        let standard_normal = Normal::new(0.0, 1.0).unwrap();
        let trace = std::env::var("TRACE").ok().is_some();
        let no_star = std::env::var("NO_STAR").ok().is_some();

        let mut est = Estimator::new(noise.clone(), dyn_params.clone(), 6.0);
        est.set_aero_table_shared(aero.clone());
        if let Some(blk_f) = noise.aero_force_rw_blackout {
            est.set_blackout_residual(
                cfg.blackout_start,
                cfg.blackout_end,
                blk_f,
                noise.aero_torque_rw_blackout,
            );
        }
        if let Some(tune) = cfg.tuning {
            est.set_robust_tuning(tune.gate_p, tune.r_scale_cap, tune.student_nu);
            if let Some(g) = tune.iae_gamma {
                est.enable_iae(g);
            } else {
                est.disable_iae();
            }
        }
        if cfg.enable_aero_learning {
            est.enable_aero_learning(0.02);
        }
        est.set_precalibration(true, cfg.blackout_start);

        let truth = init::init_truth_state(&dyn_params, &standard_normal, &mut rng);
        init::init_est_state(&mut est, &truth, &standard_normal, &mut rng);

        let bias_jump = cfg.bias_jump;
        let (bias_jump_bg, bias_jump_ba) = bias_jump_vectors(bias_jump, &standard_normal, &mut rng);

        let total_steps = (cfg.t_final / cfg.dt) as usize;
        let nees_cap = (cfg.t_final * cfg.gnss_rate).ceil() as usize + 16;

        let phi_f_const = (-cfg.dt / dyn_params.tau_aero_force).exp();
        let phi_t_const = (-cfg.dt / dyn_params.tau_aero_torque).exp();
        let inertia_inv = dyn_params
            .inertia
            .try_inverse()
            .expect("Inertia matrix must be invertible");

        let blackout_trace = trace_cfg.map(|tc| trace::TraceRecorder::new(&cfg, tc));

        Self {
            rng,
            standard_normal,
            trace,
            no_star,
            dyn_params,
            noise,
            aero,
            rcs,
            cfg,
            est,
            truth,
            next_gnss: 0.0,
            next_star: 0.0,
            next_radar: 0.0,
            gnss_reacq_done: false,
            bias_cal_applied: false,
            star_update_count: 0,
            nis_acc: 0.0,
            nis_ct: 0,
            nis_hit: 0,
            nis_phase_acc: [0.0; 4],
            nis_phase_ct: [0; 4],
            nis_phase_hit: [0; 4],
            pos_err2_acc: 0.0,
            vel_err2_acc: 0.0,
            att_err2_acc: 0.0,
            err_samples: 0.0,
            bg_err2_acc: 0.0,
            ba_err2_acc: 0.0,
            win_pos_err2_acc: 0.0,
            win_vel_err2_acc: 0.0,
            win_bg_err2_acc: 0.0,
            win_ba_err2_acc: 0.0,
            win_samples: 0.0,
            nees_samples: Vec::with_capacity(nees_cap),
            truth_times: Vec::with_capacity(total_steps + 1),
            truth_hist: Vec::with_capacity(total_steps + 1),
            post_gnss_updates: 0,
            post_gnss_accepted: 0,
            post_sat_sum: 0,
            post_sat_ct: 0,
            post_sat_min: usize::MAX,
            bias_jump,
            bias_jump_applied: false,
            bias_jump_bg,
            bias_jump_ba,
            phi_f_const,
            phi_t_const,
            inertia_inv,
            total_steps,
            blackout_trace,
        }
    }

    fn run(&mut self) {
        for step in 0..self.total_steps {
            self.step(step);
        }
    }

    fn step(&mut self, step: usize) {
        let t = step as f64 * self.cfg.dt;
        self.maybe_apply_bias_jump(t);
        self.update_truth_residuals(t);
        let (tau_b, f_b) = self.rcs.sample(t);
        let nominal_truth = nominal_body_dynamics_with_residual(
            &self.aero,
            &self.dyn_params,
            self.truth.r,
            self.truth.v,
            &self.truth.q,
            self.truth.omega,
            f_b,
            tau_b,
            self.truth.aero_res_force,
            self.truth.aero_res_torque,
        );
        self.propagate_estimator(t, tau_b, f_b, &nominal_truth);
        self.propagate_truth(&nominal_truth);
        self.record_truth(t + self.cfg.dt);
        self.maybe_apply_bias_calibration(t + self.cfg.dt);
        self.maybe_gnss(t);
        self.maybe_star(t);
        self.maybe_radar(t);
        self.accumulate_recovery_errors(t, step);
        self.maybe_record_trace(step);
    }

    fn maybe_record_trace(&mut self, step: usize) {
        let Some(trace) = &mut self.blackout_trace else {
            return;
        };
        trace.maybe_record(step, &self.est, &self.truth);
    }

    fn maybe_apply_bias_jump(&mut self, t: f64) {
        if self.bias_jump_applied {
            return;
        }
        let Some(cfg) = self.bias_jump else {
            return;
        };
        if t + 1e-12 < cfg.t {
            return;
        }

        self.truth.bg += self.bias_jump_bg;
        self.truth.ba += self.bias_jump_ba;
        self.bias_jump_applied = true;

        if self.trace {
            let bg = self.bias_jump_bg.map(|v| v.to_degrees());
            let ba_ug = self.bias_jump_ba / 9.80665 * 1.0e6;
            eprintln!(
                "[FAULT] bias jump at t={:.3}: Δbg_dps=({:.3},{:.3},{:.3}) Δba_ug=({:.1},{:.1},{:.1})",
                t, bg.x, bg.y, bg.z, ba_ug.x, ba_ug.y, ba_ug.z
            );
        }
    }

    fn maybe_apply_bias_calibration(&mut self, t: f64) {
        if self.bias_cal_applied {
            return;
        }
        if t < self.cfg.blackout_start {
            return;
        }
        let applied = self.est.apply_bias_calibration_record();
        self.bias_cal_applied = true;
        if applied && std::env::var("TRACE").ok().is_some() {
            let est_ba = self.est.state.ba / 9.80665 * 1.0e6;
            eprintln!(
                "[BA-CAL] applied at t={:.3} est_ba_ug=({:.1},{:.1},{:.1})",
                t, est_ba.x, est_ba.y, est_ba.z
            );
        }
    }
}
