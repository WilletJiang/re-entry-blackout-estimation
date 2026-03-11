use super::{Mat2NByN, MatN};
use crate::types::STATE_DIM;
use nalgebra::{linalg::Cholesky, Matrix3, SymmetricEigen};

#[inline]
pub(super) fn qr_r_upper_in_place(a: &mut Mat2NByN) -> MatN {
    let m = super::STATE_DIM2;
    let n = super::STATE_DIM;

    for k in 0..n {
        let mut norm_sq = 0.0f64;
        for i in k..m {
            let v = a[(i, k)];
            norm_sq += v * v;
        }
        let norm_x = norm_sq.sqrt();
        if norm_x == 0.0 {
            continue;
        }

        let x0 = a[(k, k)];
        let alpha = if x0 >= 0.0 { -norm_x } else { norm_x };

        let v0 = x0 - alpha;

        let mut v_norm_sq = v0 * v0;
        for i in (k + 1)..m {
            let vi = a[(i, k)];
            v_norm_sq += vi * vi;
        }
        if v_norm_sq == 0.0 {
            a[(k, k)] = alpha;
            continue;
        }

        let tau = 2.0 / v_norm_sq;

        for j in (k + 1)..n {
            let mut dot = v0 * a[(k, j)];
            for i in (k + 1)..m {
                dot += a[(i, k)] * a[(i, j)];
            }
            let w = tau * dot;
            a[(k, j)] -= w * v0;
            for i in (k + 1)..m {
                a[(i, j)] -= w * a[(i, k)];
            }
        }

        a[(k, k)] = alpha;
        for i in (k + 1)..m {
            a[(i, k)] = 0.0;
        }
    }

    let mut r = MatN::zeros();
    for i in 0..n {
        for j in i..n {
            r[(i, j)] = a[(i, j)];
        }
    }
    r
}

pub(super) fn insert_block(
    dst: &mut MatN,
    block: &Matrix3<f64>,
    row_offset: usize,
    col_offset: usize,
) {
    for r in 0..3 {
        for c in 0..3 {
            dst[(row_offset + r, col_offset + c)] = block[(r, c)];
        }
    }
}

pub(super) fn chol_lower_with_enforce(mut p: MatN, context: &str) -> MatN {
    if let Some(ch) = Cholesky::new(p.clone()) {
        return ch.l().into_owned();
    }

    let jitter = 1.0e-6;
    for i in 0..STATE_DIM {
        p[(i, i)] += jitter;
    }
    if let Some(ch) = Cholesky::new(p.clone()) {
        return ch.l().into_owned();
    }

    let mut eig = SymmetricEigen::new(p.clone());
    let mut min_eig = eig.eigenvalues[0];
    for i in 1..STATE_DIM {
        min_eig = min_eig.min(eig.eigenvalues[i]);
    }
    let floor = if min_eig.is_finite() {
        min_eig.abs().max(1.0e-8)
    } else {
        1.0e-6
    };
    for i in 0..STATE_DIM {
        if !eig.eigenvalues[i].is_finite() || eig.eigenvalues[i] < floor {
            eig.eigenvalues[i] = floor;
        }
    }
    let mut p_spd = eig.recompose();
    if let Some(ch) = Cholesky::new(p_spd.clone()) {
        return ch.l().into_owned();
    }

    let jitter2 = 1.0e-6;
    for i in 0..STATE_DIM {
        p_spd[(i, i)] += jitter2;
    }
    if let Some(ch) = Cholesky::new(p_spd.clone()) {
        return ch.l().into_owned();
    }

    let mut min_diag = p_spd[(0, 0)];
    for i in 1..STATE_DIM {
        min_diag = min_diag.min(p_spd[(i, i)]);
    }
    let eig2 = SymmetricEigen::new(p_spd.clone());
    let mut min_eig2 = eig2.eigenvalues[0];
    for i in 1..STATE_DIM {
        min_eig2 = min_eig2.min(eig2.eigenvalues[i]);
    }

    let diag_val = min_diag.abs().max(1.0e-6);
    eprintln!(
        "[WARN] chol_lower_with_enforce spectral repair failed in {}: min_diag={:.3e}, min_eig={:.3e}, using diag {:.3e}",
        context, min_diag, min_eig2, diag_val
    );
    let mut fallback = MatN::zeros();
    for i in 0..STATE_DIM {
        fallback[(i, i)] = diag_val;
    }
    if let Some(ch) = Cholesky::new(fallback.clone()) {
        return ch.l().into_owned();
    }
    panic!(
        "chol_lower_with_enforce: fallback also failed in {}; this should be unreachable",
        context
    );
}
