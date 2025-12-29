use super::{history::Event, Estimator, MeasMat, RadarMeas};
use nalgebra::{SMatrix, SVector};

impl Estimator {
    pub fn update_radar(&mut self, meas: RadarMeas) {
        self.process_radar(meas, true, true);
    }

    pub(super) fn process_radar(&mut self, meas: RadarMeas, record: bool, allow_rewind: bool) {
        let replay_events = if allow_rewind {
            self.rewind_if_needed(meas.t)
        } else {
            None
        };
        let r_norm = self.state.r.norm();
        let alt_pred = r_norm - self.dyn_params.radius_earth;
        let z = SVector::<f64, 1>::new(meas.alt - alt_pred);
        let mut h = MeasMat::<1>::zeros();
        if r_norm > 1e-6 {
            let radial = self.state.r / r_norm;
            for i in 0..3 {
                h[(0, i)] = radial[i];
            }
        }
        let r_cov = SMatrix::<f64, 1, 1>::new(meas.sigma * meas.sigma);
        self.measurement_update(z, h, r_cov, 1);
        self.maybe_checkpoint(true);
        if record {
            self.history.push_event(Event::Radar(meas));
        }
        if let Some(ev) = replay_events {
            for e in ev {
                self.replay_event(e);
            }
        }
    }
}
