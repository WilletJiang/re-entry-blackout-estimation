use nalgebra::Vector3;

pub(super) fn generate_gps_sats(t: f64, user_pos: Vector3<f64>) -> Vec<(Vector3<f64>, Vector3<f64>)> {
    const PLANES: usize = 6;
    const SATS_PER_PLANE: usize = 4;

    let mut sats = Vec::with_capacity(PLANES * SATS_PER_PLANE);
    let gps_alt = 20_200_000.0;
    let omega_gps = std::f64::consts::TAU / (12.0 * 3600.0);

    let inclination = 55.0_f64.to_radians();
    let (sin_i, cos_i) = inclination.sin_cos();
    let r_orbit = 6_378_137.0 + gps_alt;
    let v_orbit = (3.986004418e14_f64 / r_orbit).sqrt();
    let user_up = user_pos.normalize();

    for plane in 0..PLANES {
        let plane_offset = (plane as f64) * std::f64::consts::TAU / (PLANES as f64);
        for sat_idx in 0..SATS_PER_PLANE {
            let phase =
                (sat_idx as f64) * std::f64::consts::TAU / (SATS_PER_PLANE as f64) + omega_gps * t;
            let ang = phase + plane_offset;
            let (sin_a, cos_a) = ang.sin_cos();

            let sat_pos = Vector3::new(
                r_orbit * cos_a * cos_i,
                r_orbit * sin_a,
                r_orbit * cos_a * sin_i,
            );
            let sat_vel = Vector3::new(
                -v_orbit * sin_a * cos_i,
                v_orbit * cos_a,
                -v_orbit * sin_a * sin_i,
            );

            let los = sat_pos - user_pos;
            let elevation = (los.dot(&user_up) / los.norm()).asin();
            if elevation > 5.0_f64.to_radians() {
                sats.push((sat_pos, sat_vel));
            }
        }
    }
    sats
}
