use crate::types::STATE_DIM;
use nalgebra::{linalg::Cholesky, Matrix3, SMatrix, SymmetricEigen, Vector3};

pub fn bracket(grid: &[f64], x: f64) -> (usize, usize) {
    debug_assert!(!grid.is_empty());
    let n = grid.len();
    if n == 0 {
        return (0, 0);
    }
    if x <= grid[0] {
        return (0, 0);
    }
    if x >= grid[n - 1] {
        let last = n - 1;
        return (last, last);
    }

    // Binary search for the smallest hi such that grid[hi] >= x,
    // then (hi-1, hi) brackets x.
    let mut lo = 0usize;
    let mut hi = n - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if grid[mid] < x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo, hi)
}

pub fn unique_sorted(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    v
}

pub fn skew(v: Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0)
}

pub fn enforce_pd(p: &mut SMatrix<f64, STATE_DIM, STATE_DIM>) {
    // Make symmetric first (cheap, avoids Cholesky failures due to asymmetry drift).
    for i in 0..STATE_DIM {
        for j in (i + 1)..STATE_DIM {
            let a = p[(i, j)];
            let b = p[(j, i)];
            let avg = 0.5 * (a + b);
            p[(i, j)] = avg;
            p[(j, i)] = avg;
        }
    }

    // Fast check: if already SPD, do nothing (do NOT reconstruct p = L Lᵀ).
    if Cholesky::new(p.clone()).is_some() {
        return;
    }

    // Small diagonal jitter, then re-check.
    let jitter = 1.0e-9;
    for i in 0..STATE_DIM {
        p[(i, i)] += jitter;
    }
    if Cholesky::new(p.clone()).is_some() {
        return;
    }

    enforce_pd_slow(p);
}

#[cold]
fn enforce_pd_slow(p: &mut SMatrix<f64, STATE_DIM, STATE_DIM>) {
    // Spectral repair with eigenvalue flooring.
    let mut eig = SymmetricEigen::new(p.clone());
    let mut min_eig = eig.eigenvalues[0];
    for i in 1..STATE_DIM {
        min_eig = min_eig.min(eig.eigenvalues[i]);
    }
    let floor = if min_eig.is_finite() {
        min_eig.abs().max(1.0e-12)
    } else {
        1.0e-9
    };
    for i in 0..STATE_DIM {
        if !eig.eigenvalues[i].is_finite() || eig.eigenvalues[i] < floor {
            eig.eigenvalues[i] = floor;
        }
    }
    *p = eig.recompose();
    // Preserve symmetry after recomposition.
    for i in 0..STATE_DIM {
        for j in (i + 1)..STATE_DIM {
            let a = p[(i, j)];
            let b = p[(j, i)];
            let avg = 0.5 * (a + b);
            p[(i, j)] = avg;
            p[(j, i)] = avg;
        }
    }
}
