use super::{history::Event, robust::ROBUST_SENSOR_RADAR, Estimator, MeasMat, RadarMeas};
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
        let predicted_end_state = if self.seg_has_data {
            // RTS needs the propagation prior at this timestamp after the update
            // mutates `self.state`, so the pre-update state must be owned here.
            Some(self.state.clone())
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
        let accepted = self.measurement_update(ROBUST_SENSOR_RADAR, z, h, r_cov, 1);
        if !accepted {
            self.replay_tail(replay_events);
            return;
        }
        self.maybe_checkpoint_with_prediction(true, predicted_end_state);
        if record {
            self.history.push_event(Event::Radar(meas));
        }
        self.replay_tail(replay_events);
    }
}
