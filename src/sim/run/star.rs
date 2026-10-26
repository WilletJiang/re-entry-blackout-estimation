use super::{metrics, SimRun};
use nalgebra::{Matrix3, UnitQuaternion, Vector3};
use rand::Rng;
use rand_distr::{Cauchy, Distribution};

const STAR_BLACKOUT_DELAY_START_SEC: f64 = 10.0;
const STAR_BLACKOUT_DELAY_END_SEC: f64 = 30.0;

impl SimRun {
    pub(super) fn maybe_star(&mut self, t: f64) {
        if self.no_star || t + 1e-9 < self.next_star {
            return;
        }
        self.next_star += 1.0 / self.cfg.star_rate;
        let star_blackout_start = self.cfg.blackout_start + STAR_BLACKOUT_DELAY_START_SEC;
        let star_blackout_end = self.cfg.blackout_start + STAR_BLACKOUT_DELAY_END_SEC;
        let star_delay = 0.12;
        let meas_time = t - star_delay;
        if (meas_time >= star_blackout_start && meas_time <= star_blackout_end)
            || self.star_rng.gen_range(0.0..1.0) <= self.cfg.dropout_rate
        {
            return;
        }
        if meas_time < 0.0 {
            return;
        }
        let q_true_at_meas =
            metrics::find_truth_q_at(&self.truth_times, &self.truth_hist, meas_time)
                .unwrap_or(self.truth.q);
        let noise_axis = if self.star_rng.gen_range(0.0..1.0) < self.cfg.outlier_rate {
            let outlier = Cauchy::new(0.0, 0.01).unwrap();
            Vector3::from_fn(|_, _| outlier.sample(&mut self.star_rng))
        } else {
            Vector3::new(
                20e-6 * self.standard_normal.sample(&mut self.star_rng),
                20e-6 * self.standard_normal.sample(&mut self.star_rng),
                20e-6 * self.standard_normal.sample(&mut self.star_rng),
            )
        };
        let q_meas = UnitQuaternion::new_normalize(
            (UnitQuaternion::from_scaled_axis(noise_axis) * q_true_at_meas).into_inner(),
        );
        self.est.update_star(crate::estimator::StarMeas {
            t: meas_time,
            q_bi: q_meas,
            r_cov: Matrix3::identity() * (20e-6 * 20e-6),
        });
        let d2 = self
            .est
            .last_star_base_nis
            .filter(|(t, _, _)| (*t - meas_time).abs() < 1e-8)
            .map_or(f64::NAN, |(_, d, _)| d);
        self.aux_nis[0].record(d2, 3);
        self.star_update_count += 1;
    }
}
