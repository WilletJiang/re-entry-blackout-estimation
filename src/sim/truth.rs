use crate::types::State;
use nalgebra::Vector3;
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, StandardNormal};

/// Physical simulation assumptions, independent of the estimator's NoiseParams.
/// White-noise densities use Var(integral over dt) = density^2 * dt.
#[derive(Debug, Clone)]
pub struct TruthNoise {
    /// rad / sqrt(s); interval-average rate noise has sigma / sqrt(dt).
    pub gyro_noise: f64,
    /// (m/s) / sqrt(s); interval-average acceleration noise has sigma / sqrt(dt).
    pub accel_noise: f64,
    /// (rad/s) / sqrt(s).
    pub gyro_bias_rw: f64,
    /// (m/s^2) / sqrt(s).
    pub accel_bias_rw: f64,
    /// Stationary OU standard deviations, in N and Nm (not diffusion densities).
    pub aero_force_sigma: f64,
    pub aero_torque_sigma: f64,
    pub blackout_force_sigma: f64,
    pub blackout_torque_sigma: f64,
}

impl Default for TruthNoise {
    fn default() -> Self {
        Self {
            gyro_noise: 0.005_f64.to_radians(),
            accel_noise: 120e-6 * 9.80665,
            gyro_bias_rw: 7e-5_f64.to_radians(),
            accel_bias_rw: 2e-6 * 9.80665,
            aero_force_sigma: 500.0,
            aero_torque_sigma: 160.0,
            blackout_force_sigma: 1000.0,
            blackout_torque_sigma: 300.0,
        }
    }
}

impl TruthNoise {
    pub(super) fn sample_imu_noise(
        &self,
        dt: f64,
        rng: &mut StdRng,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let inv_sqrt_dt = 1.0 / dt.sqrt();
        (
            normal_vector(rng) * (self.accel_noise * inv_sqrt_dt),
            normal_vector(rng) * (self.gyro_noise * inv_sqrt_dt),
        )
    }

    pub(super) fn advance_biases(&self, state: &mut State, dt: f64, rng: &mut StdRng) {
        state.bg += normal_vector(rng) * (self.gyro_bias_rw * dt.sqrt());
        state.ba += normal_vector(rng) * (self.accel_bias_rw * dt.sqrt());
    }
}

pub(super) fn normal_vector(rng: &mut StdRng) -> Vector3<f64> {
    Vector3::new(
        StandardNormal.sample(rng),
        StandardNormal.sample(rng),
        StandardNormal.sample(rng),
    )
}

// Domain-separated seeds keep physical truth and other sensors unchanged when
// a GNSS mode consumes a different number of random draws. This is seed schema 2.
pub(super) fn seeded_stream(seed: u64, domain: u64) -> StdRng {
    let mut bytes = [0u8; 32];
    bytes[..8].copy_from_slice(&seed.to_le_bytes());
    bytes[8..16].copy_from_slice(&domain.to_le_bytes());
    bytes[16..24].copy_from_slice(&0xBADA55_u64.to_le_bytes());
    bytes[24..].copy_from_slice(&2u64.to_le_bytes());
    StdRng::from_seed(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrated_white_noise_and_bias_variance_do_not_depend_on_sample_rate() {
        let noise = TruthNoise::default();
        let trials = 4000;
        for rate in [100, 400] {
            let dt = 1.0 / rate as f64;
            let mut rng = seeded_stream(19, rate);
            let mut squares = [0.0; 4];
            let mut sums = [0.0; 4];
            for _ in 0..trials {
                let mut state = State::new();
                let mut delta_v = 0.0;
                let mut delta_angle = 0.0;
                for _ in 0..rate {
                    let (accel, gyro) = noise.sample_imu_noise(dt, &mut rng);
                    delta_v += accel.x * dt;
                    delta_angle += gyro.x * dt;
                    noise.advance_biases(&mut state, dt, &mut rng);
                }
                let normalized = [
                    delta_v / noise.accel_noise,
                    delta_angle / noise.gyro_noise,
                    state.bg.x / noise.gyro_bias_rw,
                    state.ba.x / noise.accel_bias_rw,
                ];
                for i in 0..4 {
                    sums[i] += normalized[i];
                    squares[i] += normalized[i].powi(2);
                }
            }
            for i in 0..4 {
                let mean = sums[i] / trials as f64;
                let variance = squares[i] / trials as f64 - mean.powi(2);
                assert!(mean.abs() < 0.07, "rate={rate}, component={i}, mean={mean}");
                assert!(
                    (variance - 1.0).abs() < 0.10,
                    "rate={rate}, component={i}, variance={variance}"
                );
            }
        }
    }

    #[test]
    fn zero_noise_is_exactly_deterministic() {
        let noise = TruthNoise {
            gyro_noise: 0.0,
            accel_noise: 0.0,
            gyro_bias_rw: 0.0,
            accel_bias_rw: 0.0,
            ..TruthNoise::default()
        };
        let mut rng = seeded_stream(u64::MAX, 1);
        let mut state = State::new();
        state.bg.fill(0.1);
        state.ba.fill(0.2);
        let (accel, gyro) = noise.sample_imu_noise(0.0025, &mut rng);
        noise.advance_biases(&mut state, 0.0025, &mut rng);
        assert_eq!(accel, Vector3::zeros());
        assert_eq!(gyro, Vector3::zeros());
        assert_eq!(state.bg, Vector3::repeat(0.1));
        assert_eq!(state.ba, Vector3::repeat(0.2));
    }
}
