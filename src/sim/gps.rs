use nalgebra::Vector3;

const PLANES: usize = 6;
const SATS_PER_PLANE: usize = 4;

// Idealized circular six-plane constellation, not a broadcast ephemeris.
// All vectors use the simulator's inertial frame. Public GpsSatMeas field
// names retain their legacy `_ecef` spelling; callers must use one frame.
fn satellite(t: f64, plane: usize, slot: usize) -> (Vector3<f64>, Vector3<f64>) {
    let radius: f64 = 6_378_137.0 + 20_200_000.0;
    let motion = (3.986004418e14_f64 / radius.powi(3)).sqrt();
    let node = plane as f64 * std::f64::consts::TAU / PLANES as f64;
    let phase = slot as f64 * std::f64::consts::TAU / SATS_PER_PLANE as f64
        + plane as f64 * std::f64::consts::TAU / (PLANES * SATS_PER_PLANE) as f64
        + motion * t;
    let (si, ci) = 55.0_f64.to_radians().sin_cos();
    let (sn, cn) = node.sin_cos();
    let (su, cu) = phase.sin_cos();
    let position = Vector3::new(cn * cu - sn * su * ci, sn * cu + cn * su * ci, su * si) * radius;
    let velocity =
        Vector3::new(-cn * su - sn * cu * ci, -sn * su + cn * cu * ci, cu * si) * (radius * motion);
    (position, velocity)
}

pub(super) fn generate_gps_sats(
    t: f64,
    user_pos: Vector3<f64>,
) -> Vec<(Vector3<f64>, Vector3<f64>)> {
    let mut sats = Vec::with_capacity(PLANES * SATS_PER_PLANE);
    let up = user_pos.normalize();
    for plane in 0..PLANES {
        for slot in 0..SATS_PER_PLANE {
            let (pos, vel) = satellite(t, plane, slot);
            let los = (pos - user_pos).normalize();
            if los.dot(&up) > 5.0_f64.to_radians().sin() {
                sats.push((pos, vel));
            }
        }
    }
    sats
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn distinct_orbital_planes_and_velocity_matches_position_derivative() {
        let mut normals: Vec<Vector3<f64>> = Vec::new();
        for plane in 0..PLANES {
            let (r, v) = satellite(123.0, plane, 0);
            let normal = r.cross(&v).normalize();
            assert!((normal.z - 55.0_f64.to_radians().cos()).abs() < 1e-12);
            for previous in &normals {
                assert!((normal - previous).norm() > 0.1);
            }
            normals.push(normal);
            let derivative = (satellite(123.01, plane, 0).0 - satellite(122.99, plane, 0).0) / 0.02;
            assert!((v - derivative).norm() < 1e-5);
        }
    }
    #[test]
    fn entry_geometry_observes_all_three_position_axes() {
        for t in [0.0, 120.0, 160.0, 180.0] {
            let user = Vector3::new(6_400_000.0, 0.0, 400_000.0);
            let sats = generate_gps_sats(t, user);
            assert!(sats.len() >= 4);
            let mut info = nalgebra::Matrix3::zeros();
            for (pos, _) in sats {
                let los = (pos - user).normalize();
                info += los * los.transpose();
            }
            let eig = info.symmetric_eigen().eigenvalues;
            assert!(eig.min() > 0.1 && eig.max() / eig.min() < 20.0);
        }
    }
}
