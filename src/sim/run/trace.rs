use crate::estimator::Estimator;
use crate::sim::{metrics, BlackoutTrace, BlackoutTraceConfig, SimConfig};
use crate::types::State;
use nalgebra::{linalg::Cholesky, Matrix6, SVector, Vector3};

pub(super) struct TraceRecorder {
    start_step: usize,
    end_step: usize,
    stride_steps: usize,
    pos_err_m: Vec<f32>,
    vel_err_mps: Vec<f32>,
    pos_sigma3_m: Vec<f32>,
    vel_sigma3_mps: Vec<f32>,
    nees_rv: Vec<f32>,
}

impl TraceRecorder {
    pub fn new(sim_cfg: &SimConfig, cfg: BlackoutTraceConfig) -> Self {
        let start_step = time_to_step(sim_cfg.dt, cfg.t_start);
        let end_step = time_to_step(sim_cfg.dt, cfg.t_end);
        let stride_steps = time_to_step(sim_cfg.dt, cfg.dt).max(1);
        let cap = (end_step.saturating_sub(start_step) / stride_steps) + 1;
        Self {
            start_step,
            end_step,
            stride_steps,
            pos_err_m: Vec::with_capacity(cap),
            vel_err_mps: Vec::with_capacity(cap),
            pos_sigma3_m: Vec::with_capacity(cap),
            vel_sigma3_mps: Vec::with_capacity(cap),
            nees_rv: Vec::with_capacity(cap),
        }
    }

    pub fn maybe_record(&mut self, step: usize, est: &Estimator, truth: &State) {
        if step < self.start_step || step > self.end_step {
            return;
        }
        if (step - self.start_step) % self.stride_steps != 0 {
            return;
        }

        let pos_err = est.state.r - truth.r;
        let vel_err = est.state.v - truth.v;
        let p_rv = metrics::cov_rv_from_l(&est.l_cov);

        let pos_sigma3 = pos_sigma3_from_p(&p_rv);
        let vel_sigma3 = vel_sigma3_from_p(&p_rv);
        let nees_rv = nees_rv(&pos_err, &vel_err, &p_rv).unwrap_or(f32::NAN);

        self.pos_err_m.push(f64_to_f32(pos_err.norm(), "pos_err_m"));
        self.vel_err_mps
            .push(f64_to_f32(vel_err.norm(), "vel_err_mps"));
        self.pos_sigma3_m.push(pos_sigma3);
        self.vel_sigma3_mps.push(vel_sigma3);
        self.nees_rv.push(nees_rv);
    }

    pub fn finish(self) -> BlackoutTrace {
        BlackoutTrace {
            pos_err_m: self.pos_err_m,
            vel_err_mps: self.vel_err_mps,
            pos_sigma3_m: self.pos_sigma3_m,
            vel_sigma3_mps: self.vel_sigma3_mps,
            nees_rv: self.nees_rv,
        }
    }
}

fn time_to_step(dt: f64, t: f64) -> usize {
    let ratio = t / dt;
    let rounded = ratio.round();
    let step = if (ratio - rounded).abs() <= 1.0e-9 {
        rounded
    } else {
        (ratio + 1.0e-9).floor()
    };
    (step.max(0.0)) as usize
}

fn pos_sigma3_from_p(p_rv: &Matrix6<f64>) -> f32 {
    let p_rr = p_rv.fixed_view::<3, 3>(0, 0).into_owned();
    let trace = p_rr.trace().max(0.0);
    f64_to_f32(3.0 * trace.sqrt(), "pos_sigma3_m")
}

fn vel_sigma3_from_p(p_rv: &Matrix6<f64>) -> f32 {
    let p_vv = p_rv.fixed_view::<3, 3>(3, 3).into_owned();
    let trace = p_vv.trace().max(0.0);
    f64_to_f32(3.0 * trace.sqrt(), "vel_sigma3_mps")
}

fn nees_rv(pos_err: &Vector3<f64>, vel_err: &Vector3<f64>, p_rv: &Matrix6<f64>) -> Option<f32> {
    let ch = Cholesky::new(*p_rv)?;
    let mut e = SVector::<f64, 6>::zeros();
    e.fixed_rows_mut::<3>(0).copy_from(pos_err);
    e.fixed_rows_mut::<3>(3).copy_from(vel_err);
    let x = ch.solve(&e);
    let v = e.dot(&x);
    v.is_finite().then_some(f64_to_f32(v, "nees_rv"))
}

fn f64_to_f32(value: f64, field_name: &str) -> f32 {
    if value.is_nan() {
        return f32::NAN;
    }
    if value == f64::INFINITY {
        return f32::INFINITY;
    }
    if value == f64::NEG_INFINITY {
        return f32::NEG_INFINITY;
    }

    let min_value = f64::from(f32::MIN);
    let max_value = f64::from(f32::MAX);
    if value < min_value || value > max_value {
        panic!(
            "{field_name} must fit in f32; got value={value} outside [{min_value}, {max_value}]"
        );
    }

    match value.to_string().parse::<f32>() {
        Ok(converted) => converted,
        Err(err) => {
            panic!("{field_name} conversion to f32 failed for value={value}: {err}");
        }
    }
}
