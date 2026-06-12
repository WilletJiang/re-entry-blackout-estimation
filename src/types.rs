use nalgebra::{Matrix3, UnitQuaternion, Vector3};

pub const STATE_DIM: usize = 24;

#[derive(Debug, Clone)]
pub struct DynamicsParams {
    pub mass: f64,
    pub inertia: Matrix3<f64>,
    pub s_ref: f64,
    pub b_ref: f64,
    pub c_ref: f64,
    pub radius_earth: f64,
    pub mu_earth: f64,
    pub tau_aero_force: f64,
    pub tau_aero_torque: f64,
}

impl Default for DynamicsParams {
    fn default() -> Self {
        DynamicsParams {
            mass: 8.0e4,
            inertia: Matrix3::from_diagonal(&Vector3::new(3.5e5, 4.2e5, 3.0e5)),
            s_ref: 25.0,
            b_ref: 5.5,
            c_ref: 3.0,
            radius_earth: 6_378_137.0,
            mu_earth: 3.986004418e14,
            tau_aero_force: 5.0,
            tau_aero_torque: 5.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NoiseParams {
    pub gyro_noise: f64,
    pub accel_noise: f64,
    pub gyro_bias_rw: f64,
    pub accel_bias_rw: f64,
    pub torque_noise: f64,
    pub accel_residual: f64,
    /// Dimensionless scaling applied to the translational process noise (RV block)
    /// before/within/after GNSS blackout.
    pub q_scale_rv_pre: f64,
    pub q_scale_rv_blk: f64,
    pub q_scale_rv_post: f64,
    /// Dimensionless scaling applied to accelerometer-bias RW before/within/after blackout.
    pub q_scale_ba_pre: f64,
    pub q_scale_ba_blk: f64,
    pub q_scale_ba_post: f64,
    pub aero_force_rw: f64,
    pub aero_torque_rw: f64,
    pub aero_force_rw_blackout: Option<f64>,
    pub aero_torque_rw_blackout: Option<f64>,
    /// Generic NIS significance level (α) for chi-square gates
    pub nis_alpha: f64,
    /// Star tracker NIS gate (χ² threshold) in attitude space
    pub star_nis_gate_sq: f64,
    /// GPS raw NIS gate for tight coupling (very permissive by default)
    pub gps_raw_nis_gate: f64,
    /// Maximum accelerometer-bias update per GPS raw update (m/s²)
    pub gps_raw_max_ba_update: f64,
    /// Extra gain on GNSS-driven accelerometer-bias updates (0..1)
    pub gps_ba_update_gain: f64,
    /// Maximum aerodynamic correction magnitude learned online (m/s²)
    pub aero_learning_max_correction: f64,
}

impl Default for NoiseParams {
    fn default() -> Self {
        NoiseParams {
            gyro_noise: 0.005_f64.to_radians(),
            accel_noise: 120e-6 * 9.80665,
            gyro_bias_rw: 4e-4_f64.to_radians(),
            // Reduced accel-bias RW to prevent excessive cross-covariance buildup
            // Truth ba is drawn from N(0, 150ug), so we use a smaller process noise
            accel_bias_rw: 50e-6 * 9.80665,
            torque_noise: 20.0,
            // Translational residual specific-force spectral density (m/s^2/√Hz)
            // Moderate value to balance position noise vs covariance growth
            accel_residual: 2.0,
            q_scale_rv_pre: 1.0,
            q_scale_rv_blk: 1.0,
            q_scale_rv_post: 1.0,
            q_scale_ba_pre: 1.0,
            q_scale_ba_blk: 1.0,
            q_scale_ba_post: 1.0,
            aero_force_rw: 500.0,
            aero_torque_rw: 160.0,
            aero_force_rw_blackout: None,
            aero_torque_rw_blackout: None,
            nis_alpha: 0.001,                      // 99.9% gate by default
            star_nis_gate_sq: 100.0,               // ~10σ on attitude angle-axis
            gps_raw_nis_gate: 0.0, // <=0 means "use chi-square gate(rows, nis_alpha)"
            gps_raw_max_ba_update: 5e-6 * 9.80665, // 5 µg in SI units
            gps_ba_update_gain: 1.0,
            aero_learning_max_correction: 0.5, // cap aero learning at 0.5 m/s²
        }
    }
}

impl NoiseParams {
    pub fn tuned() -> Self {
        let mut n = Self::default();
        n.gyro_bias_rw = 7e-5_f64.to_radians(); // 0.07 mdps/√Hz
                                                // Bias is nearly constant over the short mission window; keep RW small so
                                                // GNSS can actually converge it instead of it random-walking away.
        n.accel_bias_rw = 2e-6_f64 * 9.80665; // 2 µg/√Hz in SI units
                                              // Keep translational process noise in-family with the IMU spec so that
                                              // GNSS velocity can actually calibrate accelerometer bias (otherwise Q dominates).
        n.accel_residual = 0.003; // m/s²/√Hz (~300 µg)
                                  // Phase-dependent inflation: modest before/after blackout, stronger in blackout.
        n.q_scale_rv_pre = 1.0;
        n.q_scale_rv_blk = 10.0;
        n.q_scale_rv_post = 2.0;
        n.q_scale_ba_pre = 1.0;
        n.q_scale_ba_blk = 1.0;
        n.q_scale_ba_post = 1.0;
        // Tight-coupling NIS gating: use chi-square gate based on the row count.
        n.gps_raw_nis_gate = 0.0;
        // Allow measurable but bounded accelerometer-bias corrections during
        // GNSS availability windows (pre/post blackout).
        n.gps_raw_max_ba_update = 20e-6_f64 * 9.80665; // 20 µg per update
        n.gps_ba_update_gain = 0.0;
        // Increase aero residual within blackout to ~1000 N for better covariance growth
        n.aero_force_rw_blackout = Some(1000.0);
        // Also increase torque residual during blackout for attitude uncertainty
        n.aero_torque_rw_blackout = Some(300.0);
        n
    }
}
#[derive(Debug, Clone)]
pub struct ImuSample {
    pub t: f64,
    pub dt: f64,
    pub accel_meas: Vector3<f64>,
    pub gyro_meas: Vector3<f64>,
    pub torque_cmd: Vector3<f64>,
    pub force_cmd: Vector3<f64>,
}

#[derive(Debug, Clone)]
pub struct State {
    pub t: f64,
    pub r: Vector3<f64>,
    pub v: Vector3<f64>,
    pub q: UnitQuaternion<f64>,
    pub omega: Vector3<f64>,
    pub bg: Vector3<f64>,
    pub ba: Vector3<f64>,
    pub aero_res_force: Vector3<f64>, // body-frame residual force [N]
    pub aero_res_torque: Vector3<f64>, // body-frame residual torque [Nm]
}

impl State {
    pub fn new() -> Self {
        Self {
            t: 0.0,
            r: Vector3::zeros(),
            v: Vector3::zeros(),
            q: UnitQuaternion::identity(),
            omega: Vector3::zeros(),
            bg: Vector3::zeros(),
            ba: Vector3::zeros(),
            aero_res_force: Vector3::zeros(),
            aero_res_torque: Vector3::zeros(),
        }
    }
}
