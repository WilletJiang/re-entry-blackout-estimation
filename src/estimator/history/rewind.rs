use super::{history::Event, Estimator, MatN};

impl Estimator {
    pub(super) fn replay_tail(&mut self, replay_events: Option<Vec<Event>>) {
        let Some(events) = replay_events else {
            return;
        };
        for event in events {
            self.replay_event(event);
        }
    }

    pub(super) fn rewind_if_needed(&mut self, t_meas: f64) -> Option<Vec<Event>> {
        if t_meas >= self.state.t - 1e-9 {
            return None;
        }
        let idx = self.history.find_snapshot_index(t_meas)?;
        let snap = self.history.snapshots[idx].clone();
        let snap_t = snap.state.t;

        self.history.snapshots.truncate(idx + 1);
        self.history.props.truncate(idx);

        let mut replay_all: Vec<Event> = Vec::new();
        while let Some(ev) = self.history.events.back() {
            if ev.time() <= snap_t + 1e-9 {
                break;
            }
            replay_all.push(self.history.events.pop_back().unwrap());
        }
        replay_all.reverse();

        self.state = snap.state;
        self.l_cov = snap.l_cov;
        self.restore_internal(snap.internal);
        self.last_snapshot_t = snap_t;
        self.seg_phi = MatN::identity();
        self.seg_q = MatN::zeros();
        self.seg_has_data = false;

        let mut tail: Vec<Event> = Vec::new();
        for ev in replay_all {
            if ev.time() <= t_meas + 1e-9 {
                self.replay_event(ev);
                continue;
            }
            tail.push(ev);
        }
        debug_assert!(
            (self.state.t - t_meas).abs() < 1e-6,
            "rewind alignment drift: state.t={} t_meas={}",
            self.state.t,
            t_meas
        );
        Some(tail)
    }

    pub(super) fn replay_event(&mut self, ev: Event) {
        match ev {
            Event::Imu(log) => self.propagate_internal(log.sample, true),
            Event::Gnss(m) => self.process_gnss(m, true, false),
            Event::GpsRaw(m) => self.process_gps_raw(m, true, false),
            Event::Star(m) => self.process_star(m, true, false),
            Event::Radar(m) => self.process_radar(m, true, false),
            Event::InflateRv {
                t: _,
                pos_scale,
                vel_scale,
            } => {
                self.inflate_rv(pos_scale, vel_scale);
                self.maybe_checkpoint(true);
            }
            Event::BiasCalApply { t: _ } => {
                let predicted_end_state = if self.seg_has_data {
                    // RTS needs the propagation prior at this timestamp after
                    // the calibration mutates `self.state`.
                    Some(self.state.clone())
                } else {
                    None
                };
                self.apply_bias_calibration();
                self.maybe_checkpoint_with_prediction(true, predicted_end_state);
            }
        }
    }
}
