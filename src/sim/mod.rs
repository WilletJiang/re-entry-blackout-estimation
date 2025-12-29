mod gps;
mod metrics;
mod run;

pub use run::{run_monte_carlo, simulate_single};
pub use run::{run_monte_carlo_with_trace, simulate_single_with_trace, MonteCarloTraceResult};

#[derive(Clone)]
pub struct SimConfig {
    pub dt: f64,
    pub t_final: f64,
    pub blackout_start: f64,
    pub blackout_end: f64,
    pub gnss_rate: f64,
    pub star_rate: f64,
    pub radar_rate: f64,
    pub radar_max_alt: f64,
    pub outlier_rate: f64,
    pub dropout_rate: f64,
    pub tuning: Option<RobustTuning>,
    pub use_tight_coupling: bool,
    pub star_convergence_count: usize,
    pub enable_aero_learning: bool,
    pub bias_jump: Option<BiasJumpConfig>,
}

#[derive(Clone, Copy)]
pub struct BlackoutTraceConfig {
    pub t_start: f64,
    pub t_end: f64,
    pub dt: f64,
    pub blackout_start: f64,
    pub blackout_end: f64,
}

impl BlackoutTraceConfig {
    pub fn sample_count(&self) -> usize {
        let n = ((self.t_end - self.t_start) / self.dt).round() as isize;
        (n.max(0) as usize) + 1
    }

    pub fn time_at(&self, sample_idx: usize) -> f64 {
        self.t_start + (sample_idx as f64) * self.dt
    }
}

pub struct BlackoutTrace {
    pub pos_err_m: Vec<f32>,
    pub pos_sigma3_m: Vec<f32>,
    pub nees_rv: Vec<f32>,
}

impl Default for SimConfig {
    fn default() -> Self {
        SimConfig {
            dt: 0.0025,
            t_final: 180.0,
            blackout_start: 120.0,
            blackout_end: 160.0,
            gnss_rate: 10.0,
            star_rate: 2.0,
            radar_rate: 20.0,
            radar_max_alt: 25_000.0,
            outlier_rate: 0.05,
            dropout_rate: 0.20,
            tuning: None,
            use_tight_coupling: false,
            star_convergence_count: 4,
            enable_aero_learning: false,
            bias_jump: None,
        }
    }
}

#[derive(Clone, Copy)]
pub struct RobustTuning {
    pub gate_p: f64,
    pub r_scale_cap: f64,
    pub student_nu: f64,
    pub iae_gamma: Option<f64>,
    pub reacq_pos_scale: f64,
    pub reacq_vel_scale: f64,
}

#[derive(Clone, Copy)]
pub struct BiasJumpConfig {
    pub t: f64,
    pub bg_step_dps: f64,
    pub ba_step_ug: f64,
}

#[derive(Clone)]
pub struct MonteCarloResult {
    pub rms_pos: f64,
    pub rms_vel: f64,
    pub rms_att_deg: f64,
    pub nees_mean: f64,
    pub nis_mean: f64,
    pub nees_covered: f64,
    pub nis_covered: f64,
    pub nees_r_mean: f64,
    pub nees_v_mean: f64,
    pub nees_ba_mean: f64,
    pub nees_r_covered: f64,
    pub nees_v_covered: f64,
    pub nees_ba_covered: f64,
    pub nees_mean_pre: f64,
    pub nees_mean_blk: f64,
    pub nees_mean_post: f64,
    pub nees_mean_tail: f64,
    pub nees_covered_pre: f64,
    pub nees_covered_blk: f64,
    pub nees_covered_post: f64,
    pub nees_covered_tail: f64,
    pub nis_mean_pre: f64,
    pub nis_mean_blk: f64,
    pub nis_mean_post: f64,
    pub nis_mean_tail: f64,
    pub nis_covered_pre: f64,
    pub nis_covered_blk: f64,
    pub nis_covered_post: f64,
    pub nis_covered_tail: f64,
    pub post_gnss_updates: usize,
    pub post_gnss_accepted: usize,
    pub post_sat_mean: f64,
    pub post_sat_min: usize,
    pub smoothed_rms_pos: f64,
    pub smoothed_rms_vel: f64,
    pub bias_bg_rms_dps: f64,
    pub bias_ba_rms_ug: f64,
    pub passed: bool,
}
