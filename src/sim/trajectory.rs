use crate::{
    data::AeroTable,
    physics::isa_atmosphere,
    types::{DynamicsParams, State},
};

/// Truth-only evidence about whether a campaign exercises the intended scenario.
/// These diagnostics do not certify the aerodynamic model or the trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryDiagnostics {
    pub initial_altitude_m: f64,
    pub initial_radial_velocity_mps: f64,
    pub final_altitude_m: f64,
    pub min_altitude_m: f64,
    pub max_altitude_m: f64,
    pub min_mach: f64,
    pub max_mach: f64,
    pub samples: usize,
    pub outside_aero_table_samples: usize,
    pub radar_eligible_samples: usize,
    pub finite: bool,
}

impl TrajectoryDiagnostics {
    pub(super) fn new(state: &State, params: &DynamicsParams) -> Self {
        let altitude = state.r.norm() - params.radius_earth;
        Self {
            initial_altitude_m: altitude,
            initial_radial_velocity_mps: state.r.dot(&state.v) / state.r.norm(),
            final_altitude_m: altitude,
            min_altitude_m: altitude,
            max_altitude_m: altitude,
            min_mach: f64::INFINITY,
            max_mach: f64::NEG_INFINITY,
            samples: 0,
            outside_aero_table_samples: 0,
            radar_eligible_samples: 0,
            finite: true,
        }
    }

    pub(super) fn record(
        &mut self,
        state: &State,
        params: &DynamicsParams,
        table: &AeroTable,
        radar_max_alt: f64,
    ) {
        self.samples += 1;
        let altitude = state.r.norm() - params.radius_earth;
        let (_, sound_speed) = isa_atmosphere(altitude);
        let mach = state.v.norm() / sound_speed;
        let body_velocity = state.q.inverse_transform_vector(&state.v);
        let alpha = body_velocity.z.atan2(body_velocity.x).to_degrees();
        let beta = (body_velocity.y / body_velocity.norm())
            .clamp(-1.0, 1.0)
            .asin()
            .to_degrees();
        if ![altitude, mach, alpha, beta].iter().all(|v| v.is_finite()) {
            self.finite = false;
            return;
        }
        self.final_altitude_m = altitude;
        self.min_altitude_m = self.min_altitude_m.min(altitude);
        self.max_altitude_m = self.max_altitude_m.max(altitude);
        self.min_mach = self.min_mach.min(mach);
        self.max_mach = self.max_mach.max(mach);
        if !within_axis(mach, &table.mach)
            || !within_axis(alpha, &table.alpha)
            || !within_axis(beta, &table.beta)
        {
            self.outside_aero_table_samples += 1;
        }
        if altitude >= 0.0 && altitude < radar_max_alt {
            self.radar_eligible_samples += 1;
        }
    }
}

fn within_axis(value: f64, axis: &[f64]) -> bool {
    match (axis.first(), axis.last()) {
        (Some(min), Some(max)) => value >= *min && value <= *max,
        _ => false,
    }
}
