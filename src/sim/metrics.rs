use crate::estimator::MatN;
use crate::sim::SimConfig;
use nalgebra::{Matrix3, Matrix6, SVector, UnitQuaternion, Vector3};

#[derive(Clone)]
pub(super) struct NeesSample {
    pub phase: u8,
    pub gnss_accepted: bool,
    pub err_rv: SVector<f64, 6>,
    pub p_rv: Matrix6<f64>,
    pub err_ba: Vector3<f64>,
    pub p_ba: Matrix3<f64>,
}

#[inline]
pub(super) fn phase_id(t: f64, cfg: &SimConfig) -> usize {
    if t < cfg.blackout_start {
        return 0;
    }
    if t <= cfg.blackout_end {
        return 1;
    }
    if t <= cfg.blackout_end + 2.0 {
        return 2;
    }
    3
}

#[inline]
pub(super) fn cov_rv_from_l(l: &MatN) -> Matrix6<f64> {
    let mut p = Matrix6::<f64>::zeros();
    for i in 0..6 {
        for j in 0..=i {
            let mut s = 0.0f64;
            for k in 0..=j {
                s += l[(i, k)] * l[(j, k)];
            }
            p[(i, j)] = s;
            p[(j, i)] = s;
        }
    }
    p
}

#[inline]
pub(super) fn cov_ba_from_l(l: &MatN) -> Matrix3<f64> {
    let mut p = Matrix3::<f64>::zeros();
    for ii in 0..3 {
        let i = 15 + ii;
        for jj in 0..=ii {
            let j = 15 + jj;
            let mut s = 0.0f64;
            for k in 0..=j {
                s += l[(i, k)] * l[(j, k)];
            }
            p[(ii, jj)] = s;
            p[(jj, ii)] = s;
        }
    }
    p
}

pub(super) fn find_truth_q_at(
    times: &[f64],
    hist: &[crate::types::State],
    t_target: f64,
) -> Option<UnitQuaternion<f64>> {
    if hist.is_empty() {
        return None;
    }
    let idx = times
        .iter()
        .enumerate()
        .rev()
        .find(|(_, &ti)| ti <= t_target)
        .map(|(i, _)| i)?;
    Some(hist[idx].q)
}
