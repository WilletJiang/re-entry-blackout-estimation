use super::{gps, metrics, SimRun};
use crate::estimator::{GnssMeas, GpsRawMeas, GpsSatMeas};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, SVector, Vector3};
use rand::Rng;
use rand_distr::{Cauchy, Distribution};

const POS_SIGMA_M: f64 = 1.9;
const PVT_VEL_SIGMA_MPS: f64 = 0.08;
const TIGHT_RANGE_RATE_SIGMA_MPS: f64 = 0.02;
const PVT_OUTLIER_SCALE_M: f64 = 200.0;
const TIGHT_PSEUDORANGE_OUTLIER_SCALE_M: f64 = 50.0;
const RECOVERY_WINDOW_SEC: f64 = 2.0;
const FLOAT_TIME_TOL_SEC: f64 = 1.0e-9;

impl SimRun {
    pub(super) fn maybe_gnss(&mut self, t: f64) {
        let mut gnss_updated = false;
        let mut gnss_update_time = None;
        let mut sat_count_meas = 0usize;
        let mut sat_sig_v_ls = None;

        let t = t + self.cfg.dt;
        if t + FLOAT_TIME_TOL_SEC < self.next_gnss {
            return;
        }
        self.next_gnss += 1.0 / self.cfg.gnss_rate;

        let t_meas = t;
        let in_blackout = t_meas >= self.cfg.blackout_start - FLOAT_TIME_TOL_SEC
            && t_meas + FLOAT_TIME_TOL_SEC < self.cfg.blackout_end;
        let allow_gnss = self.star_update_count >= self.cfg.star_convergence_count
            || t_meas >= self.cfg.blackout_start;

        let meas_considered =
            allow_gnss && !in_blackout && self.gnss_rng.gen_range(0.0..1.0) > self.cfg.dropout_rate;
        if meas_considered {
            {
                let (updated, sat_count, sat_sig) = self.gnss_update_normal(t_meas);
                gnss_updated = updated;
                gnss_update_time = updated.then_some(t_meas);
                sat_count_meas = sat_count;
                sat_sig_v_ls = sat_sig;
            }
        }

        self.finish_gnss_step(
            gnss_updated,
            gnss_update_time,
            t_meas,
            sat_count_meas,
            sat_sig_v_ls,
        );
    }

    fn finish_gnss_step(
        &mut self,
        gnss_updated: bool,
        gnss_update_time: Option<f64>,
        nees_time: f64,
        sat_count_meas: usize,
        sat_sig_v_ls: Option<Vector3<f64>>,
    ) {
        let gnss_accepted = gnss_update_time.is_some_and(|update_time| {
            self.est
                .last_gnss_accepted
                .is_some_and(|(t_acc, ok)| ok && (t_acc - update_time).abs() < 1e-6)
        });

        if let Some(update_time) = gnss_update_time {
            self.record_post_blackout_gnss_stats(
                gnss_updated,
                update_time,
                gnss_accepted,
                sat_count_meas,
                sat_sig_v_ls,
            );
            self.accumulate_nis(update_time);
            self.trace_bias_near_boundaries(update_time);
            self.record_nees_sample(update_time, gnss_accepted);
            return;
        }

        self.record_nees_sample(nees_time, false);
    }

    fn record_post_blackout_gnss_stats(
        &mut self,
        gnss_updated: bool,
        update_time: f64,
        gnss_accepted: bool,
        sat_count_meas: usize,
        sat_sig_v_ls: Option<Vector3<f64>>,
    ) {
        if gnss_updated
            && update_time >= self.cfg.blackout_end
            && update_time <= self.cfg.blackout_end + RECOVERY_WINDOW_SEC
        {
            self.post_gnss_updates += 1;
            if gnss_accepted {
                self.post_gnss_accepted += 1;
            }
            if self.cfg.use_tight_coupling && sat_count_meas > 0 {
                self.post_sat_sum += sat_count_meas;
                self.post_sat_ct += 1;
                self.post_sat_min = self.post_sat_min.min(sat_count_meas);
            }
        }

        if std::env::var("TRACE_GNSS").ok().is_some()
            && gnss_updated
            && update_time >= self.cfg.blackout_end
            && update_time <= self.cfg.blackout_end + RECOVERY_WINDOW_SEC
        {
            self.trace_post_gnss(update_time, gnss_accepted, sat_count_meas, sat_sig_v_ls);
        }
    }

    fn gnss_update_normal(&mut self, t_meas: f64) -> (bool, usize, Option<Vector3<f64>>) {
        if self.cfg.use_tight_coupling {
            return self.gnss_update_tight(t_meas);
        }
        self.gnss_update_pvt(t_meas);
        (true, 0, None)
    }

    fn gnss_update_tight(&mut self, t_meas: f64) -> (bool, usize, Option<Vector3<f64>>) {
        let sats = gps::generate_gps_sats(t_meas, self.truth.r);
        if sats.len() < 4 {
            return (false, 0, None);
        }
        let is_outlier = self.gnss_rng.gen_range(0.0..1.0) < self.cfg.outlier_rate;
        let mut sat_meas = Vec::with_capacity(sats.len());
        let mut info_v = Matrix3::<f64>::zeros();
        for (sat_pos, sat_vel) in sats {
            let true_range = (sat_pos - self.truth.r).norm();
            let true_rr = (sat_vel - self.truth.v).dot(&(sat_pos - self.truth.r).normalize());

            let rr_meas = true_rr
                + TIGHT_RANGE_RATE_SIGMA_MPS * self.standard_normal.sample(&mut self.gnss_rng);
            let pr_meas = if is_outlier {
                let cauchy = cauchy_with_scale(
                    TIGHT_PSEUDORANGE_OUTLIER_SCALE_M,
                    "tight pseudorange outlier scale",
                );
                true_range + cauchy.sample(&mut self.gnss_rng)
            } else {
                let pr_noise = 1.5 * self.standard_normal.sample(&mut self.gnss_rng);
                true_range + pr_noise
            };

            sat_meas.push(GpsSatMeas {
                pos_ecef: sat_pos,
                vel_ecef: sat_vel,
                pr_m: pr_meas,
                pr_sigma: 1.5,
                range_rate_mps: Some(rr_meas),
                rr_sigma: Some(TIGHT_RANGE_RATE_SIGMA_MPS),
            });

            let dr = sat_pos - self.est.state.r;
            let rho = dr.norm();
            if rho > 1e-6 {
                let los = dr / rho;
                let w = 1.0 / (TIGHT_RANGE_RATE_SIGMA_MPS * TIGHT_RANGE_RATE_SIGMA_MPS);
                info_v += w * (los * los.transpose());
            }
        }

        let sat_sig_v_ls = if std::env::var("TRACE_GNSS").ok().is_some()
            && t_meas >= self.cfg.blackout_end
            && t_meas <= self.cfg.blackout_end + 2.0
        {
            Cholesky::new(info_v).map(|ch| {
                let cov_v = ch.inverse();
                Vector3::new(
                    cov_v[(0, 0)].abs().sqrt(),
                    cov_v[(1, 1)].abs().sqrt(),
                    cov_v[(2, 2)].abs().sqrt(),
                )
            })
        } else {
            None
        };

        let sat_count = sat_meas.len();
        self.est.update_gps_raw(GpsRawMeas {
            t: t_meas,
            sats: sat_meas,
        });
        (true, sat_count, sat_sig_v_ls)
    }

    fn gnss_update_pvt(&mut self, t_meas: f64) {
        let meas = self.sample_pvt_measurement(t_meas, PVT_VEL_SIGMA_MPS, true);
        self.est.update_gnss(meas);
    }

    fn sample_pvt_measurement(
        &mut self,
        t_meas: f64,
        vel_sigma_mps: f64,
        allow_pos_outlier: bool,
    ) -> GnssMeas {
        let mut pos_meas = self.truth.r;
        let mut vel_meas = self.truth.v;
        let is_pos_outlier =
            allow_pos_outlier && self.gnss_rng.gen_range(0.0..1.0) < self.cfg.outlier_rate;
        if is_pos_outlier {
            let cauchy = cauchy_with_scale(PVT_OUTLIER_SCALE_M, "PVT position outlier scale");
            pos_meas += Vector3::new(
                cauchy.sample(&mut self.gnss_rng),
                cauchy.sample(&mut self.gnss_rng),
                cauchy.sample(&mut self.gnss_rng),
            );
        } else {
            pos_meas += Vector3::new(
                POS_SIGMA_M * self.standard_normal.sample(&mut self.gnss_rng),
                POS_SIGMA_M * self.standard_normal.sample(&mut self.gnss_rng),
                POS_SIGMA_M * self.standard_normal.sample(&mut self.gnss_rng),
            );
        }
        vel_meas += Vector3::new(
            vel_sigma_mps * self.standard_normal.sample(&mut self.gnss_rng),
            vel_sigma_mps * self.standard_normal.sample(&mut self.gnss_rng),
            vel_sigma_mps * self.standard_normal.sample(&mut self.gnss_rng),
        );
        let r_cov = pvt_covariance(vel_sigma_mps);
        GnssMeas {
            t: t_meas,
            pos: pos_meas,
            vel: vel_meas,
            r_cov,
        }
    }

    fn trace_post_gnss(
        &mut self,
        t_meas: f64,
        gnss_accepted: bool,
        sat_count: usize,
        sat_sig_v_ls: Option<Vector3<f64>>,
    ) {
        let pe = self.est.state.r - self.truth.r;
        let ve = self.est.state.v - self.truth.v;
        let p_rv = metrics::cov_rv_from_l(&self.est.l_cov);
        let sig_v = Vector3::new(
            p_rv[(3, 3)].abs().sqrt(),
            p_rv[(4, 4)].abs().sqrt(),
            p_rv[(5, 5)].abs().sqrt(),
        );
        let mut sig_att_deg = [0.0f64; 3];
        for ii in 0..3 {
            let i = 6 + ii;
            let mut s = 0.0f64;
            for k in 0..=i {
                let v = self.est.l_cov[(i, k)];
                s += v * v;
            }
            sig_att_deg[ii] = s.sqrt().to_degrees();
        }
        let p_full = self.est.cov();
        let sig_att_deg_full = [
            p_full[(6, 6)].abs().sqrt().to_degrees(),
            p_full[(7, 7)].abs().sqrt().to_degrees(),
            p_full[(8, 8)].abs().sqrt().to_degrees(),
        ];
        let sig_bg_dps = [
            p_full[(12, 12)].abs().sqrt().to_degrees(),
            p_full[(13, 13)].abs().sqrt().to_degrees(),
            p_full[(14, 14)].abs().sqrt().to_degrees(),
        ];
        let sig_ba_ug = [
            p_full[(15, 15)].abs().sqrt() / 9.80665 * 1.0e6,
            p_full[(16, 16)].abs().sqrt() / 9.80665 * 1.0e6,
            p_full[(17, 17)].abs().sqrt() / 9.80665 * 1.0e6,
        ];
        let nis = self.est.last_gnss_nis.and_then(|(t_nis, d2, dof)| {
            if (t_nis - t_meas).abs() < 1e-6 {
                Some((d2, dof))
            } else {
                None
            }
        });
        if let Some((d2, dof)) = nis {
            eprintln!(
                "[POST-GNSS] t={:.3} ok={} sats={} |pe|={:.3} |ve|={:.3} sig_v=({:.3},{:.3},{:.3}) sig_att_deg=({:.4},{:.4},{:.4}) sig_att_deg_full=({:.4},{:.4},{:.4}) sig_bg_dps=({:.6},{:.6},{:.6}) sig_ba_ug=({:.1},{:.1},{:.1}) gnss_sig_v=({:.3},{:.3},{:.3}) NIS={:.2} dof={}",
                t_meas,
                gnss_accepted as u8,
                sat_count,
                pe.norm(),
                ve.norm(),
                sig_v.x,
                sig_v.y,
                sig_v.z,
                sig_att_deg[0],
                sig_att_deg[1],
                sig_att_deg[2],
                sig_att_deg_full[0],
                sig_att_deg_full[1],
                sig_att_deg_full[2],
                sig_bg_dps[0],
                sig_bg_dps[1],
                sig_bg_dps[2],
                sig_ba_ug[0],
                sig_ba_ug[1],
                sig_ba_ug[2],
                sat_sig_v_ls.map(|v| v.x).unwrap_or(f64::NAN),
                sat_sig_v_ls.map(|v| v.y).unwrap_or(f64::NAN),
                sat_sig_v_ls.map(|v| v.z).unwrap_or(f64::NAN),
                d2,
                dof
            );
        } else {
            eprintln!(
                "[POST-GNSS] t={:.3} ok={} sats={} |pe|={:.3} |ve|={:.3} sig_v=({:.3},{:.3},{:.3}) sig_att_deg=({:.4},{:.4},{:.4}) sig_att_deg_full=({:.4},{:.4},{:.4}) gnss_sig_v=({:.3},{:.3},{:.3})",
                t_meas,
                gnss_accepted as u8,
                sat_count,
                pe.norm(),
                ve.norm(),
                sig_v.x,
                sig_v.y,
                sig_v.z,
                sig_att_deg[0],
                sig_att_deg[1],
                sig_att_deg[2],
                sig_att_deg_full[0],
                sig_att_deg_full[1],
                sig_att_deg_full[2],
                sat_sig_v_ls.map(|v| v.x).unwrap_or(f64::NAN),
                sat_sig_v_ls.map(|v| v.y).unwrap_or(f64::NAN),
                sat_sig_v_ls.map(|v| v.z).unwrap_or(f64::NAN)
            );
        }
    }

    fn accumulate_nis(&mut self, t_meas: f64) {
        let (_, d2, dof) = self
            .est
            .last_gnss_base_nis
            .filter(|(t, _, _)| (*t - t_meas).abs() < 1e-6)
            .unwrap_or((t_meas, f64::NAN, 0));

        let covered = crate::sim::innovation::in_interval(d2, dof);
        let ph = metrics::phase_id(t_meas, &self.cfg);
        self.nis_acc += d2;
        self.nis_ct += 1;
        if covered {
            self.nis_hit += 1;
        }
        self.nis_phase_acc[ph] += d2;
        self.nis_phase_ct[ph] += 1;
        if covered {
            self.nis_phase_hit[ph] += 1;
        }
    }

    fn trace_bias_near_boundaries(&mut self, t_meas: f64) {
        if std::env::var("TRACE").ok().is_none() {
            return;
        }
        let near_pre_end =
            t_meas >= self.cfg.blackout_start - 0.5 && t_meas <= self.cfg.blackout_start;
        let near_post_start =
            t_meas >= self.cfg.blackout_end && t_meas <= self.cfg.blackout_end + 0.05;
        if !(near_pre_end || near_post_start) {
            return;
        }
        let est_ba = self.est.state.ba / 9.80665 * 1.0e6;
        let tru_ba = self.truth.ba / 9.80665 * 1.0e6;
        eprintln!(
            "[BA] t={:.3} est_ba_ug=({:.1},{:.1},{:.1}) truth_ba_ug=({:.1},{:.1},{:.1})",
            t_meas, est_ba.x, est_ba.y, est_ba.z, tru_ba.x, tru_ba.y, tru_ba.z
        );
    }

    fn record_nees_sample(&mut self, t_meas: f64, gnss_accepted: bool) {
        let ph = metrics::phase_id(t_meas, &self.cfg) as u8;
        let pe = self.est.state.r - self.truth.r;
        let ve = self.est.state.v - self.truth.v;
        let err_rv = SVector::<f64, 6>::new(pe.x, pe.y, pe.z, ve.x, ve.y, ve.z);
        let p_rv = metrics::cov_rv_from_l(&self.est.l_cov);
        let err_ba = self.est.state.ba - self.truth.ba;
        let p_ba = metrics::cov_ba_from_l(&self.est.l_cov);
        self.nees_samples.push(metrics::NeesSample {
            phase: ph,
            gnss_accepted,
            err_rv,
            p_rv,
            err_ba,
            p_ba,
        });
    }
}

fn cauchy_with_scale(scale: f64, label: &str) -> Cauchy<f64> {
    match Cauchy::new(0.0, scale) {
        Ok(dist) => dist,
        Err(err) => {
            panic!("{label} must be finite and > 0; got scale={scale}: {err}");
        }
    }
}

fn pvt_covariance(vel_sigma_mps: f64) -> Matrix6<f64> {
    let mut r_cov = Matrix6::zeros();
    for i in 0..3 {
        r_cov[(i, i)] = POS_SIGMA_M * POS_SIGMA_M;
        r_cov[(i + 3, i + 3)] = vel_sigma_mps * vel_sigma_mps;
    }
    r_cov
}
