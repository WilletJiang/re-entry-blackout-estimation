use crate::estimator::Estimator;
use crate::types::{DynamicsParams, State};
use nalgebra::{UnitQuaternion, Vector3};
use rand::rngs::StdRng;
use rand_distr::{Distribution, Normal};

pub(super) fn init_truth_state(
    dyn_params: &DynamicsParams,
    standard_normal: &Normal<f64>,
    rng: &mut StdRng,
) -> State {
    let mut truth = State::new();
    truth.r = Vector3::new(dyn_params.radius_earth + 60_000.0, 0.0, 0.0);
    truth.v = Vector3::new(7_600.0 * 0.7660, 0.0, -7_600.0 * 0.6428);
    truth.bg = Vector3::new(
        0.3_f64.to_radians() / 3600.0 * standard_normal.sample(rng),
        0.3_f64.to_radians() / 3600.0 * standard_normal.sample(rng),
        0.3_f64.to_radians() / 3600.0 * standard_normal.sample(rng),
    );
    truth.ba = Vector3::new(
        150e-6 * 9.80665 * standard_normal.sample(rng),
        150e-6 * 9.80665 * standard_normal.sample(rng),
        150e-6 * 9.80665 * standard_normal.sample(rng),
    );
    truth
}

pub(super) fn init_est_state(
    est: &mut Estimator,
    truth: &State,
    standard_normal: &Normal<f64>,
    rng: &mut StdRng,
) {
    est.state.t = 0.0;
    est.state.r = truth.r
        + Vector3::new(
            Normal::new(0.0, 5_000.0).unwrap().sample(rng),
            Normal::new(0.0, 5_000.0).unwrap().sample(rng),
            Normal::new(0.0, 5_000.0).unwrap().sample(rng),
        );
    est.state.v = truth.v
        + Vector3::new(
            50.0 * standard_normal.sample(rng),
            50.0 * standard_normal.sample(rng),
            50.0 * standard_normal.sample(rng),
        );
    let misalign_axis = Vector3::new(
        1.0_f64.to_radians() * standard_normal.sample(rng),
        1.0_f64.to_radians() * standard_normal.sample(rng),
        1.0_f64.to_radians() * standard_normal.sample(rng),
    );
    est.state.q = UnitQuaternion::new_normalize(
        (UnitQuaternion::from_scaled_axis(misalign_axis) * truth.q).into_inner(),
    );
    est.state.omega = truth.omega;
}
