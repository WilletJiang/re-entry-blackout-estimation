#[cfg(feature = "diagnostics")]
use nalgebra::linalg::Cholesky;
#[cfg(feature = "diagnostics")]
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(feature = "diagnostics")]
use super::MatN;

#[cfg(feature = "diagnostics")]
static COV_CHECK_COUNTER: AtomicUsize = AtomicUsize::new(0);

#[cfg(feature = "diagnostics")]
#[inline]
pub(super) fn maybe_check_covariance(p: &MatN, context: &str) {
    const STRIDE: usize = 1024;
    let k = COV_CHECK_COUNTER.fetch_add(1, Ordering::Relaxed);
    if k % STRIDE != 0 {
        return;
    }

    let mut max_asym = 0.0f64;
    for i in 0..super::STATE_DIM {
        for j in (i + 1)..super::STATE_DIM {
            max_asym = max_asym.max((p[(i, j)] - p[(j, i)]).abs());
        }
    }
    if max_asym > 1.0e-8 {
        eprintln!(
            "[WARN] covariance symmetry drift in {}: max_asym={:.3e}",
            context, max_asym
        );
    }
    if Cholesky::new(*p).is_none() {
        eprintln!("[WARN] covariance not SPD in {} (Cholesky failed)", context);
    }
}
