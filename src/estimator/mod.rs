use crate::data::AeroTable;
use crate::types::{DynamicsParams, NoiseParams, State, STATE_DIM};
use nalgebra::{Matrix3, Matrix6, SMatrix, SVector, UnitQuaternion, Vector3};
use std::sync::Arc;

pub type MatN = SMatrix<f64, STATE_DIM, STATE_DIM>;
pub type VecN = SVector<f64, STATE_DIM>;
pub type MeasMat<const M: usize> = SMatrix<f64, M, STATE_DIM>;
pub type MeasVec<const M: usize> = SVector<f64, M>;

const STATE_DIM2: usize = STATE_DIM * 2;
type Mat2NByN = SMatrix<f64, STATE_DIM2, STATE_DIM>;

#[path = "calibration/bias_cal.rs"]
mod bias_cal;
#[path = "filter/core.rs"]
mod core;
#[path = "support/diag.rs"]
mod diag;
#[path = "update/gating.rs"]
mod gating;
#[path = "history/mod.rs"]
mod history;
#[path = "support/math.rs"]
mod math;
#[path = "support/noise.rs"]
mod noise;
#[path = "filter/propagate.rs"]
mod propagate;
#[path = "history/rewind.rs"]
mod rewind;
#[path = "filter/rts.rs"]
mod rts;
#[path = "filter/state_ops.rs"]
mod state_ops;
#[path = "update/common.rs"]
mod update_common;
#[path = "update/gnss.rs"]
mod update_gnss;
#[path = "update/gps_raw.rs"]
mod update_gps_raw;
#[path = "update/radar.rs"]
mod update_radar;
#[path = "update/star.rs"]
mod update_star;

#[cfg(test)]
#[path = "tests/mod.rs"]
mod tests;

#[derive(Clone)]
pub struct GnssMeas {
    pub t: f64,
    pub pos: Vector3<f64>,
    pub vel: Vector3<f64>,
    pub r_cov: Matrix6<f64>,
}

/// Single satellite measurement for tight-coupled GPS/INS
#[derive(Clone)]
pub struct GpsSatMeas {
    pub pos_ecef: Vector3<f64>,
    pub vel_ecef: Vector3<f64>,
    pub pr_m: f64,
    pub pr_sigma: f64,
    pub range_rate_mps: Option<f64>,
    pub rr_sigma: Option<f64>,
}

/// Raw GPS measurements for tight coupling
#[derive(Clone)]
pub struct GpsRawMeas {
    pub t: f64,
    pub sats: Vec<GpsSatMeas>,
}

#[derive(Clone)]
pub struct StarMeas {
    pub t: f64,
    pub q_bi: UnitQuaternion<f64>,
    pub r_cov: Matrix3<f64>,
}

#[derive(Clone)]
pub struct RadarMeas {
    pub t: f64,
    pub alt: f64,
    pub sigma: f64,
}

pub struct Estimator {
    pub state: State,
    pub l_cov: MatN, // lower Cholesky so P = L L^T
    noise: NoiseParams,
    pub dyn_params: DynamicsParams,
    history: history::HistoryBuffer,
    history_snapshot_dt: f64,
    last_snapshot_t: f64,
    seg_phi: MatN,
    seg_q: MatN,
    seg_has_data: bool,
    pub aero_table: Option<Arc<AeroTable>>,
    blackout_window: Option<(f64, f64)>,
    aero_force_rw_blackout: Option<f64>,
    aero_torque_rw_blackout: Option<f64>,
    pub last_gnss_nis: Option<(f64, f64, usize)>,
    pub last_gnss_prior_t: Option<f64>,
    pub last_gnss_prior_rv_p: Option<Matrix6<f64>>,
    pub last_gnss_prior_r: Option<Vector3<f64>>,
    pub last_gnss_prior_v: Option<Vector3<f64>>,
    pub last_gnss_post_t: Option<f64>,
    pub last_gnss_post_rv_p: Option<Matrix6<f64>>,
    pub last_gnss_post_r: Option<Vector3<f64>>,
    pub last_gnss_post_v: Option<Vector3<f64>>,
    pub last_gnss_accepted: Option<(f64, bool)>,
    gnss_reject_streak: u32,
    aero_learning_enabled: bool,
    last_v_meas: Option<(f64, Vector3<f64>)>,
    pub learned_aero_correction: Vector3<f64>,
    aero_learning_rate: f64,
    pub allow_bias_precalibration: bool,
    pub precalibration_end_time: f64,
    bias_cal: bias_cal::BiasCal,
}

#[derive(Debug, Clone)]
pub struct EstimatorHealth {
    pub t: f64,
    pub last_gnss_nis: Option<(f64, f64, usize)>,
    pub last_gnss_accepted: Option<(f64, bool)>,
    pub gnss_reject_streak: u32,
    pub history_snapshots: usize,
    pub history_events: usize,
    pub history_props: usize,
    pub history_capacity_bytes: usize,
}
