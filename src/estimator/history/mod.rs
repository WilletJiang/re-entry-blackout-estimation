use super::{bias_cal::BiasCal, GnssMeas, GpsRawMeas, MatN, RadarMeas, StarMeas};
use crate::types::ImuSample;
use std::collections::VecDeque;

const SNAPSHOT_RESERVE_PAD: usize = 32;
const EVENT_RESERVE_PAD: usize = 128;

#[derive(Clone)]
pub(super) struct Snapshot {
    pub state: crate::types::State,
    pub l_cov: MatN,
    pub internal: SnapshotInternal,
}

#[derive(Clone)]
pub(super) struct PropRecord {
    pub phi: MatN,
    pub q: MatN,
    pub predicted_end_state: crate::types::State,
}

#[derive(Clone)]
pub(super) struct LoggedImu {
    pub sample: ImuSample,
}

#[derive(Clone)]
pub(super) struct SnapshotInternal {
    pub bias_cal: BiasCal,
    pub gnss_reject_streak: u32,
    pub robust_gate_alpha: Option<f64>,
    pub robust_r_scale_cap: f64,
    pub robust_student_nu: f64,
    pub iae_gamma: Option<f64>,
    pub iae_scales: [f64; super::robust::ROBUST_SENSOR_COUNT],
    pub aero_learning_enabled: bool,
    pub aero_learning_rate: f64,
    pub last_v_meas: Option<(f64, nalgebra::Vector3<f64>)>,
    pub learned_aero_correction: nalgebra::Vector3<f64>,
}

#[derive(Clone)]
pub(super) enum Event {
    Imu(LoggedImu),
    Gnss(GnssMeas),
    GpsRaw(GpsRawMeas),
    Star(StarMeas),
    Radar(RadarMeas),
    InflateRv {
        t: f64,
        pos_scale: f64,
        vel_scale: f64,
    },
    BiasCalApply {
        t: f64,
    },
}

impl Event {
    pub(super) fn time(&self) -> f64 {
        match self {
            Event::Imu(l) => l.sample.t,
            Event::Gnss(m) => m.t,
            Event::GpsRaw(m) => m.t,
            Event::Star(m) => m.t,
            Event::Radar(m) => m.t,
            Event::InflateRv { t, .. } => *t,
            Event::BiasCalApply { t } => *t,
        }
    }
}

pub(super) struct HistoryBuffer {
    horizon: f64,
    pub snapshots: VecDeque<Snapshot>,
    pub props: VecDeque<Option<PropRecord>>,
    pub events: VecDeque<Event>,
}

impl HistoryBuffer {
    pub(super) fn new(horizon: f64, snapshot_dt: f64) -> Self {
        let imu_rate_hz_est = 400.0;
        let snapshot_dt = snapshot_dt.clamp(0.0025, 0.5);
        let snap_rate_hz = 1.0 / snapshot_dt;

        let snap_cap = ((horizon * snap_rate_hz).ceil() as usize + 8).clamp(64, 10_000);
        let ev_cap = ((horizon * imu_rate_hz_est).ceil() as usize + 128).clamp(256, 50_000);
        Self {
            horizon,
            snapshots: VecDeque::with_capacity(snap_cap),
            props: VecDeque::with_capacity(snap_cap),
            events: VecDeque::with_capacity(ev_cap),
        }
    }

    pub(super) fn push_snapshot(&mut self, snap: Snapshot) {
        reserve_small_pad_if_full(&mut self.snapshots, SNAPSHOT_RESERVE_PAD);
        self.snapshots.push_back(snap);
        self.pad_missing_prop_intervals();
        self.prune();
    }

    pub(super) fn push_prop(&mut self, prop: PropRecord) {
        reserve_small_pad_if_full(&mut self.props, SNAPSHOT_RESERVE_PAD);
        self.props.push_back(Some(prop));
    }

    pub(super) fn push_event(&mut self, ev: Event) {
        debug_assert!(
            self.events
                .back()
                .map(|last| ev.time() + 1e-9 >= last.time())
                .unwrap_or(true),
            "event time went backwards: new={} last={}",
            ev.time(),
            self.events.back().map(|e| e.time()).unwrap_or(f64::NAN)
        );
        reserve_small_pad_if_full(&mut self.events, EVENT_RESERVE_PAD);
        self.events.push_back(ev);
        self.prune();
    }

    pub(super) fn prune(&mut self) {
        let Some(newest) = self.snapshots.back() else {
            return;
        };
        let cutoff = newest.state.t - self.horizon;
        while self
            .snapshots
            .front()
            .map(|s| s.state.t < cutoff)
            .unwrap_or(false)
        {
            self.snapshots.pop_front();
            if !self.props.is_empty() {
                self.props.pop_front();
            }
        }
        while self
            .events
            .front()
            .map(|e| e.time() < cutoff)
            .unwrap_or(false)
        {
            self.events.pop_front();
        }
        while self.props.len() + 1 > self.snapshots.len() && !self.props.is_empty() {
            self.props.pop_front();
        }
    }

    pub(super) fn find_snapshot_index(&self, t: f64) -> Option<usize> {
        self.snapshots
            .iter()
            .enumerate()
            .rev()
            .find(|(_, s)| s.state.t <= t)
            .map(|(i, _)| i)
    }

    fn pad_missing_prop_intervals(&mut self) {
        while self.props.len() + 1 < self.snapshots.len() {
            reserve_small_pad_if_full(&mut self.props, SNAPSHOT_RESERVE_PAD);
            self.props.push_back(None);
        }
    }
}

fn reserve_small_pad_if_full<T>(buffer: &mut VecDeque<T>, reserve_pad: usize) {
    if buffer.len() < buffer.capacity() {
        return;
    }
    buffer.reserve_exact(reserve_pad);
}
