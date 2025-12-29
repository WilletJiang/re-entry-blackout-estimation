use super::{Estimator, VecN};
use nalgebra::UnitQuaternion;

impl Estimator {
    pub(super) fn apply_state_error(&mut self, dx: VecN) {
        self.state.r += dx.fixed_rows::<3>(0).into_owned();
        self.state.v += dx.fixed_rows::<3>(3).into_owned();
        let dq = UnitQuaternion::from_scaled_axis(dx.fixed_rows::<3>(6).into_owned());
        self.state.q = UnitQuaternion::new_normalize((dq * self.state.q).into_inner());
        self.state.omega += dx.fixed_rows::<3>(9).into_owned();
        self.state.bg += dx.fixed_rows::<3>(12).into_owned();
        self.state.ba += dx.fixed_rows::<3>(15).into_owned();
        self.state.aero_res_force += dx.fixed_rows::<3>(18).into_owned();
        self.state.aero_res_torque += dx.fixed_rows::<3>(21).into_owned();
    }
}
