use entry_blackout_estimation::estimator::Estimator;
use entry_blackout_estimation::types::{DynamicsParams, ImuSample, NoiseParams};
use nalgebra::Vector3;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

struct CountingAlloc;

thread_local! {
    static ALLOC_COUNT: Cell<usize> = Cell::new(0);
}

#[inline]
fn alloc_count() -> usize {
    ALLOC_COUNT.with(|c| c.get())
}

unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOC_COUNT.with(|c| c.set(c.get() + 1));
        System.alloc(layout)
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOC_COUNT.with(|c| c.set(c.get() + 1));
        System.alloc_zeroed(layout)
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        System.dealloc(ptr, layout)
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOC_COUNT.with(|c| c.set(c.get() + 1));
        System.realloc(ptr, layout, new_size)
    }
}

#[global_allocator]
static A: CountingAlloc = CountingAlloc;

#[test]
fn propagate_has_no_heap_allocations_after_warmup() {
    let noise = NoiseParams::default();
    let dyn_params = DynamicsParams::default();
    let mut est = Estimator::new(noise, dyn_params, 1.0);

    // Make the state physically sane so any optional aero/atmos paths remain well-behaved.
    est.state.r = Vector3::new(est.dyn_params.radius_earth + 1.0, 0.0, 0.0);

    let dt = 0.0025;
    let warm = 600usize;
    for k in 1..=warm {
        est.propagate(ImuSample {
            t: k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }

    let before = alloc_count();
    let iters = 600usize;
    let t0 = est.state.t;
    for k in 1..=iters {
        est.propagate(ImuSample {
            t: t0 + k as f64 * dt,
            dt,
            accel_meas: Vector3::zeros(),
            gyro_meas: Vector3::zeros(),
            torque_cmd: Vector3::zeros(),
            force_cmd: Vector3::zeros(),
        });
    }
    let after = alloc_count();
    assert_eq!(
        after, before,
        "allocations detected: before={} after={}",
        before, after
    );
    let before_smooth = alloc_count();
    let smoothed = est.rts_smooth(0.5);
    let after_smooth = alloc_count();
    assert!(smoothed.is_some());
    assert_eq!(
        before_smooth, after_smooth,
        "RTS must use bounded stack workspace"
    );
}
