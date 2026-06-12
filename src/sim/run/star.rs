use super::{metrics, SimRun};
use nalgebra::{Matrix3, UnitQuaternion, Vector3};
use rand::Rng;
use rand_distr::Distribution;

const STAR_BLACKOUT_DELAY_START_SEC: f64 = 10.0;
const STAR_BLACKOUT_DELAY_END_SEC: f64 = 30.0;

impl SimRun {
    pub(super) fn maybe_star(&mut self, t: f64) {
        if self.no_star || t < self.next_star {
            return;
        }
        self.next_star += 1.0 / self.cfg.star_rate;
        let star_blackout_start = self.cfg.blackout_start + STAR_BLACKOUT_DELAY_START_SEC;
        let star_blackout_end = self.cfg.blackout_start + STAR_BLACKOUT_DELAY_END_SEC;
        if (t >= star_blackout_start && t <= star_blackout_end)
            || self.rng.gen_range(0.0..1.0) <= self.cfg.dropout_rate
        {
            return;
        }
        let star_delay = 0.12;
        let meas_time = t - star_delay;
        if meas_time < 0.0 {
            return;
        }
        let q_true_at_meas =
            metrics::find_truth_q_at(&self.truth_times, &self.truth_hist, meas_time)
                .unwrap_or(self.truth.q);
        let noise_axis = Vector3::new(
            20e-6 * self.standard_normal.sample(&mut self.rng),
            20e-6 * self.standard_normal.sample(&mut self.rng),
            20e-6 * self.standard_normal.sample(&mut self.rng),
        );
        let q_meas = UnitQuaternion::new_normalize(
            (UnitQuaternion::from_scaled_axis(noise_axis) * q_true_at_meas).into_inner(),
        );
        self.est.update_star(crate::estimator::StarMeas {
            t: meas_time,
            q_bi: q_meas,
            r_cov: Matrix3::identity() * (60e-6 * 60e-6),
        });
        self.star_update_count += 1;
    }
}
