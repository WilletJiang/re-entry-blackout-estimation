use super::{gps, metrics, SimRun, CHI2_6_HI, CHI2_6_LO, Z_95};
use crate::estimator::{GnssMeas, GpsRawMeas, GpsSatMeas};
use nalgebra::{linalg::Cholesky, Matrix3, Matrix6, SVector, Vector3};
use rand::Rng;
use rand_distr::{Cauchy, Distribution, Normal};

impl SimRun {
    pub(super) fn maybe_gnss(&mut self, t: f64) {
        if t < self.next_gnss {
            return;
        }
        self.next_gnss += 1.0 / self.cfg.gnss_rate;

        let t_meas = t + self.cfg.dt;
        let in_blackout = (t_meas >= self.cfg.blackout_start) && (t_meas <= self.cfg.blackout_end);
        let allow_gnss = self.star_update_count >= self.cfg.star_convergence_count
            || t_meas >= self.cfg.blackout_start;

        let mut gnss_updated = false;
        let mut sat_count_meas = 0usize;
        let mut sat_sig_v_ls: Option<Vector3<f64>> = None;

        let meas_considered =
            allow_gnss && !in_blackout && self.rng.gen_range(0.0..1.0) > self.cfg.dropout_rate;
        if meas_considered {
            if !self.gnss_reacq_done && t_meas >= self.cfg.blackout_end {
                self.gnss_reacquire(t_meas);
                gnss_updated = true;
                self.gnss_reacq_done = true;
            } else {
                let (updated, sat_count, sat_sig) = self.gnss_update_normal(t_meas);
                gnss_updated = updated;
                sat_count_meas = sat_count;
                sat_sig_v_ls = sat_sig;
            }
        }

        let gnss_accepted = self
            .est
            .last_gnss_accepted
            .is_some_and(|(t_acc, ok)| ok && (t_acc - t_meas).abs() < 1e-6);

        if gnss_updated && t_meas >= self.cfg.blackout_end && t_meas <= self.cfg.blackout_end + 2.0
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
            && t_meas >= self.cfg.blackout_end
            && t_meas <= self.cfg.blackout_end + 2.0
        {
            self.trace_post_gnss(t_meas, gnss_accepted, sat_count_meas, sat_sig_v_ls);
        }

        if gnss_accepted {
            self.accumulate_nis(t_meas);
        }
        if meas_considered {
            self.accumulate_window_gnss_rms(t_meas);
            self.trace_bias_near_boundaries(t_meas);
        }

        self.record_nees_sample(t_meas, gnss_accepted);
    }

    fn gnss_reacquire(&mut self, t_meas: f64) {
        let mut pos_meas = self.truth.r;
        let mut vel_meas = self.truth.v;
        let n = Normal::new(0.0, 1.9).unwrap();
        pos_meas += Vector3::new(
            n.sample(&mut self.rng),
            n.sample(&mut self.rng),
            n.sample(&mut self.rng),
        );
        let vel_sigma = if self.cfg.use_tight_coupling {
            0.02
        } else {
            0.125
        };
        let nv = Normal::new(0.0, vel_sigma).unwrap();
        vel_meas += Vector3::new(
            nv.sample(&mut self.rng),
            nv.sample(&mut self.rng),
            nv.sample(&mut self.rng),
        );
        let mut r_cov = Matrix6::zeros();
        for i in 0..3 {
            r_cov[(i, i)] = 1.9_f64.powi(2);
        }
        for i in 0..3 {
            r_cov[(i + 3, i + 3)] = vel_sigma.powi(2);
        }

        if let Some(tune) = self.cfg.tuning {
            self.est
                .inflate_rv_record(tune.reacq_pos_scale, tune.reacq_vel_scale);
        } else {
            self.est.inflate_rv_record(1000.0, 100.0);
        }
        self.est.update_gnss(GnssMeas {
            t: t_meas,
            pos: pos_meas,
            vel: vel_meas,
            r_cov,
        });
        if std::env::var("TRACE").ok().is_some() {
            let est_ba = self.est.state.ba / 9.80665 * 1.0e6;
            let tru_ba = self.truth.ba / 9.80665 * 1.0e6;
            eprintln!(
                "[BA] t={:.3} est_ba_ug=({:.1},{:.1},{:.1}) truth_ba_ug=({:.1},{:.1},{:.1})",
                t_meas, est_ba.x, est_ba.y, est_ba.z, tru_ba.x, tru_ba.y, tru_ba.z
            );
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
        let is_outlier = self.rng.gen_range(0.0..1.0) < self.cfg.outlier_rate;
        let mut sat_meas = Vec::with_capacity(sats.len());
        let mut info_v = Matrix3::<f64>::zeros();
        for (sat_pos, sat_vel) in sats {
            let true_range = (sat_pos - self.truth.r).norm();
            let true_rr = (sat_vel - self.truth.v).dot(&(sat_pos - self.truth.r).normalize());

            let (pr_meas, rr_meas) = if is_outlier {
                let cauchy = Cauchy::new(0.0, 50.0).unwrap();
                (true_range + cauchy.sample(&mut self.rng), true_rr)
            } else {
                let pr_noise = 1.5 * self.standard_normal.sample(&mut self.rng);
                let rr_noise = 0.02 * self.standard_normal.sample(&mut self.rng);
                (true_range + pr_noise, true_rr + rr_noise)
            };

            sat_meas.push(GpsSatMeas {
                pos_ecef: sat_pos,
                vel_ecef: sat_vel,
                pr_m: pr_meas,
                pr_sigma: 1.5,
                range_rate_mps: Some(rr_meas),
                rr_sigma: Some(0.02),
            });

            let dr = sat_pos - self.est.state.r;
            let rho = dr.norm();
            if rho > 1e-6 {
                let los = dr / rho;
                let w = 1.0 / (0.02 * 0.02);
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
        let mut pos_meas = self.truth.r;
        let mut vel_meas = self.truth.v;
        if self.rng.gen_range(0.0..1.0) < self.cfg.outlier_rate {
            let cauchy = Cauchy::new(0.0, 200.0).unwrap();
            pos_meas += Vector3::new(
                cauchy.sample(&mut self.rng),
                cauchy.sample(&mut self.rng),
                cauchy.sample(&mut self.rng),
            );
        } else {
            let n = 1.9 * self.standard_normal.sample(&mut self.rng);
            pos_meas += Vector3::new(
                n,
                1.9 * self.standard_normal.sample(&mut self.rng),
                1.9 * self.standard_normal.sample(&mut self.rng),
            );
            vel_meas += Vector3::new(
                0.08 * self.standard_normal.sample(&mut self.rng),
                0.08 * self.standard_normal.sample(&mut self.rng),
                0.08 * self.standard_normal.sample(&mut self.rng),
            );
        }
        let mut r_cov = Matrix6::zeros();
        for i in 0..3 {
            r_cov[(i, i)] = 1.9_f64.powi(2);
        }
        for i in 0..3 {
            r_cov[(i + 3, i + 3)] = 0.08_f64.powi(2);
        }
        self.est.update_gnss(GnssMeas {
            t: t_meas,
            pos: pos_meas,
            vel: vel_meas,
            r_cov,
        });
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
        let nis = self.est.last_gnss_nis.and_then(|(t_nis, d2, dof)| {
            if (t_nis - t_meas).abs() < 1e-6 {
                Some((d2, dof))
            } else {
                None
            }
        });
        if let Some((d2, dof)) = nis {
            eprintln!(
                "[POST-GNSS] t={:.3} ok={} sats={} |pe|={:.3} |ve|={:.3} sig_v=({:.3},{:.3},{:.3}) sig_att_deg=({:.4},{:.4},{:.4}) sig_att_deg_full=({:.4},{:.4},{:.4}) gnss_sig_v=({:.3},{:.3},{:.3}) NIS={:.2} dof={}",
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
        let Some((t_nis, d2, dof)) = self.est.last_gnss_nis else {
            return;
        };
        if dof == 0 || (t_nis - t_meas).abs() >= 1e-6 || !d2.is_finite() {
            return;
        }

        let (chi2_lo, chi2_hi) = if dof == 6 {
            (CHI2_6_LO, CHI2_6_HI)
        } else {
            let nu = dof as f64;
            let a = 1.0 - 2.0 / (9.0 * nu);
            let b = (2.0 / (9.0 * nu)).sqrt();
            let lo = nu * (a - Z_95 * b).max(0.0).powi(3);
            let hi = nu * (a + Z_95 * b).powi(3);
            (lo, hi)
        };
        let ph = metrics::phase_id(t_meas, &self.cfg);
        self.nis_acc += d2;
        self.nis_ct += 1;
        if d2 >= chi2_lo && d2 <= chi2_hi {
            self.nis_hit += 1;
        }
        self.nis_phase_acc[ph] += d2;
        self.nis_phase_ct[ph] += 1;
        if d2 >= chi2_lo && d2 <= chi2_hi {
            self.nis_phase_hit[ph] += 1;
        }
    }

    fn accumulate_window_gnss_rms(&mut self, t_meas: f64) {
        if t_meas < self.cfg.blackout_end || t_meas > self.cfg.blackout_end + 2.0 {
            return;
        }
        let pe = self.est.state.r - self.truth.r;
        let ve = self.est.state.v - self.truth.v;
        self.win_pos_err2_acc += pe.norm_squared();
        self.win_vel_err2_acc += ve.norm_squared();
        let bge = self.est.state.bg - self.truth.bg;
        let bae = self.est.state.ba - self.truth.ba;
        self.win_bg_err2_acc += bge.norm_squared();
        self.win_ba_err2_acc += bae.norm_squared();
        self.win_samples += 1.0;
        if std::env::var("TRACE").ok().is_some() {
            eprintln!(
                "[WIN] t={:.3} |pos_err|={:.3} |vel_err|={:.3}",
                t_meas,
                pe.norm(),
                ve.norm()
            );
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
