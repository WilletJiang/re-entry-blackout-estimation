//! Single-thread estimator resource probe. Run on the target for target evidence.
use entry_blackout_estimation::{
    estimator::{Estimator, GnssMeas, GpsSatMeas, RadarMeas, StarMeas},
    sim::scenario::reference_entry,
    types::{DynamicsParams, ImuSample, NoiseParams},
};
use nalgebra::{Matrix3, Matrix6, UnitQuaternion, Vector3};
use std::{
    alloc::{GlobalAlloc, Layout, System},
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};

struct Meter;
static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
static ALLOCS: AtomicUsize = AtomicUsize::new(0);
fn added(bytes: usize) {
    let live = LIVE.fetch_add(bytes, Ordering::Relaxed) + bytes;
    PEAK.fetch_max(live, Ordering::Relaxed);
    ALLOCS.fetch_add(1, Ordering::Relaxed);
}
unsafe impl GlobalAlloc for Meter {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc(layout);
        if !p.is_null() {
            added(layout.size());
        }
        p
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let p = System.alloc_zeroed(layout);
        if !p.is_null() {
            added(layout.size());
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
        System.dealloc(p, layout);
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let next = System.realloc(p, layout, size);
        if !next.is_null() {
            LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
            added(size);
        }
        next
    }
}
#[global_allocator]
static ALLOCATOR: Meter = Meter;
const STACK: usize = 512 * 1024;
const CYCLES: usize = 2000;

// Linux clock ID from <bits/time.h>. CPU time excludes descheduling; wall
// time below remains the relevant end-to-end deadline observation.
#[cfg(target_os = "linux")]
fn thread_cpu_ms() -> f64 {
    #[repr(C)]
    struct Timespec {
        seconds: std::ffi::c_long,
        nanoseconds: std::ffi::c_long,
    }
    unsafe extern "C" {
        fn clock_gettime(clock: std::ffi::c_int, value: *mut Timespec) -> std::ffi::c_int;
    }
    let mut value = Timespec {
        seconds: 0,
        nanoseconds: 0,
    };
    // SAFETY: the Linux ABI takes a writable timespec of two C long values.
    let status = unsafe { clock_gettime(3, &mut value) };
    if status == 0 {
        value.seconds as f64 * 1000.0 + value.nanoseconds as f64 / 1e6
    } else {
        f64::NAN
    }
}
#[cfg(not(target_os = "linux"))]
fn thread_cpu_ms() -> f64 {
    f64::NAN
}

struct ProbeResult {
    heap: usize,
    total: usize,
    wall_p99: f64,
    wall_max: f64,
    cpu_p99: f64,
    cpu_max: f64,
    wall_over_budget: usize,
    cpu_over_budget: usize,
    rts10_ms: f64,
    rts2_ms: f64,
    rts_allocations: usize,
}

fn probe(tight: bool) -> ProbeResult {
    let baseline = LIVE.load(Ordering::Relaxed);
    PEAK.store(baseline, Ordering::Relaxed);
    let params = DynamicsParams::default();
    let (_, aero) = reference_entry(&params);
    let radius = params.radius_earth + 10_000.0;
    let truth_r = Vector3::new(radius, 0.0, 0.0);
    let accel = Vector3::new(params.mu_earth / radius.powi(2), 0.0, 0.0);
    let mut est = Estimator::new(NoiseParams::tuned(), params, 10.0);
    est.set_precalibration(false, 120.0);
    est.set_aero_table(aero);
    est.state.r = truth_r;
    let sats: Vec<_> = (0..8)
        .map(|i| {
            let a = i as f64 * std::f64::consts::TAU / 8.0;
            let pos = truth_r + Vector3::new(0.6, a.cos(), a.sin()).normalize() * 20_200_000.0;
            GpsSatMeas {
                pos_ecef: pos,
                vel_ecef: Vector3::zeros(),
                pr_m: (pos - truth_r).norm(),
                pr_sigma: 1.5,
                range_rate_mps: Some(0.0),
                rr_sigma: Some(0.02),
            }
        })
        .collect();
    let mut cycles = [0.0_f64; CYCLES];
    let mut cpu_cycles = [0.0_f64; CYCLES];
    for cycle in 0..CYCLES {
        let cpu_start = thread_cpu_ms();
        let start = Instant::now();
        for sub in 1..=4 {
            let t = (cycle * 4 + sub) as f64 * 0.0025;
            est.propagate(ImuSample {
                t,
                dt: 0.0025,
                accel_meas: accel,
                gyro_meas: Vector3::zeros(),
                force_cmd: Vector3::zeros(),
                torque_cmd: Vector3::zeros(),
            });
        }
        let t = (cycle + 1) as f64 * 0.01;
        if (cycle + 1) % 10 == 0 {
            if tight {
                est.update_gps_raw_ref(t, &sats);
            } else {
                est.update_gnss(GnssMeas {
                    t,
                    pos: truth_r,
                    vel: Vector3::zeros(),
                    r_cov: Matrix6::from_diagonal(&nalgebra::SVector::from_row_slice(&[
                        3.61, 3.61, 3.61, 0.0064, 0.0064, 0.0064,
                    ])),
                });
            }
        }
        if (cycle + 1) % 50 == 0 && t >= 0.12 {
            est.update_star(StarMeas {
                t: t - 0.12,
                q_bi: UnitQuaternion::identity(),
                r_cov: Matrix3::identity() * 4e-10,
            });
        }
        if (cycle + 1) % 5 == 0 && t >= 0.07 {
            est.update_radar(RadarMeas {
                t: t - 0.07,
                alt: 10_000.0,
                sigma: 1.0,
            });
        }
        cycles[cycle] = start.elapsed().as_secs_f64() * 1000.0;
        cpu_cycles[cycle] = thread_cpu_ms() - cpu_start;
    }
    let before_smooth = ALLOCS.load(Ordering::Relaxed);
    let start = Instant::now();
    assert!(est.rts_smooth(10.0).is_some());
    let smooth_ms = start.elapsed().as_secs_f64() * 1000.0;
    let smooth_allocs = ALLOCS.load(Ordering::Relaxed) - before_smooth;
    let start = Instant::now();
    let mut count = 0;
    assert!(est.visit_smoothed(2.0, |_, _| count += 1).is_some());
    assert_eq!(count, 201);
    let recovery_smooth_ms = start.elapsed().as_secs_f64() * 1000.0;
    assert_eq!(
        ALLOCS.load(Ordering::Relaxed),
        before_smooth,
        "streaming RTS allocated"
    );
    let heap = PEAK.load(Ordering::Relaxed) - baseline;
    let wall_over_budget = cycles.iter().filter(|&&ms| ms > 5.0).count();
    let cpu_over_budget = cpu_cycles.iter().filter(|&&ms| ms > 5.0).count();
    #[cfg(target_os = "linux")]
    assert!(cpu_cycles.iter().all(|ms| ms.is_finite() && *ms >= 0.0));
    cycles.sort_by(f64::total_cmp);
    cpu_cycles.sort_by(f64::total_cmp);
    assert!(est.cov().cholesky().is_some());
    ProbeResult {
        heap,
        total: heap + STACK,
        wall_p99: cycles[(CYCLES * 99) / 100],
        wall_max: cycles[CYCLES - 1],
        cpu_p99: cpu_cycles[(CYCLES * 99) / 100],
        cpu_max: cpu_cycles[CYCLES - 1],
        wall_over_budget,
        cpu_over_budget,
        rts10_ms: smooth_ms,
        rts2_ms: recovery_smooth_ms,
        rts_allocations: smooth_allocs,
    }
}
fn main() {
    let tight = std::env::args().any(|arg| arg == "--tight");
    let result = std::thread::Builder::new()
        .stack_size(STACK)
        .spawn(move || probe(tight))
        .unwrap()
        .join()
        .unwrap();
    println!(
        "arch={} mode={} cycles={} history_s=10 reserved_stack_bytes={STACK}",
        std::env::consts::ARCH,
        if tight { "tight" } else { "pvt" },
        CYCLES
    );
    println!("peak_heap_bytes={} heap_plus_reserved_stack_bytes={} cycle_p99_ms={:.6} cycle_max_ms={:.6} rts_ms={:.6} rts_allocations={}", result.heap, result.total, result.wall_p99, result.wall_max, result.rts10_ms, result.rts_allocations);
    println!("cycle_cpu_p99_ms={:.6} cycle_cpu_max_ms={:.6} wall_cycles_over_5ms={} cpu_cycles_over_5ms={}", result.cpu_p99, result.cpu_max, result.wall_over_budget, result.cpu_over_budget);
    println!("recovery_window_rts_ms={:.6}", result.rts2_ms);
    assert!(
        result.total <= 20_000_000,
        "estimator exceeds 20 MB including thread stack"
    );
    assert_eq!(result.rts_allocations, 0, "RTS allocated on heap");
    if std::env::args().any(|arg| arg == "--check-time") {
        assert!(
            result.wall_max <= 5.0,
            "wall-clock cycle budget exceeded; target and scheduling evidence required"
        );
    }
}
