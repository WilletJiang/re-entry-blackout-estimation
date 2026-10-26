use super::SimRun;
use rand::Rng;
use rand_distr::{Cauchy, Distribution};

const RADAR_OUTLIER_SCALE_M: f64 = 10.0;

impl SimRun {
    pub(super) fn maybe_radar(&mut self, t: f64) {
        if t + 1e-9 < self.next_radar {
            return;
        }
        self.next_radar += 1.0 / self.cfg.radar_rate;
        let radar_delay = 0.07;
        let meas_time = t - radar_delay;
        if meas_time < 0.0 {
            return;
        }
        let alt_true = truth_alt_at(
            &self.truth_times,
            &self.truth_hist,
            meas_time,
            self.dyn_params.radius_earth,
            &self.truth,
        );
        if alt_true >= self.cfg.radar_max_alt
            || self.radar_rng.gen_range(0.0..1.0) <= self.cfg.dropout_rate
        {
            return;
        }
        let alt_meas = alt_true
            + if self.radar_rng.gen_range(0.0..1.0) < self.cfg.outlier_rate {
                cauchy_with_scale(RADAR_OUTLIER_SCALE_M, "radar altitude outlier scale")
                    .sample(&mut self.radar_rng)
            } else {
                1.0 * self.standard_normal.sample(&mut self.radar_rng)
            };
        self.est.update_radar(crate::estimator::RadarMeas {
            t: meas_time,
            alt: alt_meas,
            sigma: 1.0,
        });
        let d2 = self
            .est
            .last_radar_base_nis
            .filter(|(t, _, _)| (*t - meas_time).abs() < 1e-8)
            .map_or(f64::NAN, |(_, d, _)| d);
        self.aux_nis[1].record(d2, 1);
    }
}

fn cauchy_with_scale(scale: f64, label: &str) -> Cauchy<f64> {
    match Cauchy::new(0.0, scale) {
        Ok(dist) => dist,
        Err(err) => {
            panic!("{label} must be finite and > 0; got scale={scale}: {err}");
        }
    }
}

fn truth_alt_at(
    times: &[f64],
    hist: &[crate::types::State],
    t_target: f64,
    radius_earth: f64,
    truth: &crate::types::State,
) -> f64 {
    let Some((lower, upper, weight)) = crate::sim::metrics::truth_bracket(times, t_target) else {
        return truth.r.norm() - radius_earth;
    };
    let position = hist[lower].r * (1.0 - weight) + hist[upper].r * weight;
    position.norm() - radius_earth
}
