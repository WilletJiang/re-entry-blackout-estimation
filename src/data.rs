use crate::types::DynamicsParams;
use crate::util::unique_sorted;
use anyhow::Result;
use nalgebra::Vector3;
use serde::Deserialize;
use std::path::Path;

#[derive(Debug, Deserialize)]
#[allow(non_snake_case)]
struct VehiclePropsRow {
    #[serde(rename = "S_ref_m2")]
    s_ref_m2: f64,
    #[serde(rename = "c_ref_m")]
    c_ref_m: f64,
    #[serde(rename = "b_ref_m")]
    b_ref_m: f64,
    #[serde(rename = "Jx_kgm2")]
    jx_kgm2: f64,
    #[serde(rename = "Jy_kgm2")]
    jy_kgm2: f64,
    #[serde(rename = "Jz_kgm2")]
    jz_kgm2: f64,
    #[serde(rename = "mass_kg")]
    mass_kg: f64,
}

pub fn load_vehicle_props(path: impl AsRef<Path>) -> Result<DynamicsParams> {
    let mut rdr = csv::Reader::from_path(path)?;
    let row: VehiclePropsRow = rdr.deserialize().next().ok_or_else(|| {
        anyhow::anyhow!("vehicle properties file empty or parse error")
    })??;
    Ok(DynamicsParams {
        mass: row.mass_kg,
        inertia: nalgebra::Matrix3::from_diagonal(&Vector3::new(
            row.jx_kgm2, row.jy_kgm2, row.jz_kgm2,
        )),
        s_ref: row.s_ref_m2,
        b_ref: row.b_ref_m,
        c_ref: row.c_ref_m,
        ..DynamicsParams::default()
    })
}

#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
struct AeroRow {
    Mach: f64,
    alpha_deg: f64,
    beta_deg: f64,
    Cx: f64,
    Cy: f64,
    Cz: f64,
    Cl: f64,
    Cm: f64,
    Cn: f64,
    Cl_p: f64,
    Cm_q: f64,
    Cn_r: f64,
}

#[derive(Debug, Clone)]
pub struct AeroCoeff {
    pub cx: f64,
    pub cy: f64,
    pub cz: f64,
    pub cl: f64,
    pub cm: f64,
    pub cn: f64,
    pub cl_p: f64,
    pub cm_q: f64,
    pub cn_r: f64,
}

#[derive(Debug, Clone)]
pub struct AeroTable {
    pub mach: Vec<f64>,
    pub alpha: Vec<f64>,
    pub beta: Vec<f64>,
    pub coeffs: Vec<AeroCoeff>,
    pub n_a: usize,
    pub n_b: usize,
}

impl AeroTable {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let mut rdr = csv::Reader::from_path(path)?;
        let mut rows = Vec::new();
        for rec in rdr.deserialize() {
            let row: AeroRow = rec?;
            rows.push(row);
        }
        rows.sort_by(|a, b| {
            a.Mach
                .partial_cmp(&b.Mach)
                .unwrap()
                .then(a.alpha_deg.partial_cmp(&b.alpha_deg).unwrap())
                .then(a.beta_deg.partial_cmp(&b.beta_deg).unwrap())
        });
        let mach = unique_sorted(rows.iter().map(|r| r.Mach).collect());
        let alpha = unique_sorted(rows.iter().map(|r| r.alpha_deg).collect());
        let beta = unique_sorted(rows.iter().map(|r| r.beta_deg).collect());
        let n_m = mach.len();
        let n_a = alpha.len();
        let n_b = beta.len();
        let mut coeffs = vec![
            AeroCoeff {
                cx: 0.0,
                cy: 0.0,
                cz: 0.0,
                cl: 0.0,
                cm: 0.0,
                cn: 0.0,
                cl_p: 0.0,
                cm_q: 0.0,
                cn_r: 0.0,
            };
            n_m * n_a * n_b
        ];
        for r in rows {
            let mi = mach.iter().position(|&m| (m - r.Mach).abs() < 1e-8).unwrap();
            let ai = alpha.iter().position(|&a| (a - r.alpha_deg).abs() < 1e-8).unwrap();
            let bi = beta.iter().position(|&b| (b - r.beta_deg).abs() < 1e-8).unwrap();
            let idx = (mi * n_a + ai) * n_b + bi;
            coeffs[idx] = AeroCoeff {
                cx: r.Cx,
                cy: r.Cy,
                cz: r.Cz,
                cl: r.Cl,
                cm: r.Cm,
                cn: r.Cn,
                cl_p: r.Cl_p,
                cm_q: r.Cm_q,
                cn_r: r.Cn_r,
            };
        }
        Ok(AeroTable {
            mach,
            alpha,
            beta,
            coeffs,
            n_a,
            n_b,
        })
    }

    pub fn interp(&self, mach: f64, alpha_deg: f64, beta_deg: f64) -> AeroCoeff {
        // Locate bracketing indices along each axis using a binary search.
        let (mi0, mi1) = crate::util::bracket(&self.mach, mach);
        let (ai0, ai1) = crate::util::bracket(&self.alpha, alpha_deg);
        let (bi0, bi1) = crate::util::bracket(&self.beta, beta_deg);

        let m0 = self.mach[mi0];
        let m1 = self.mach[mi1];
        let a0 = self.alpha[ai0];
        let a1 = self.alpha[ai1];
        let b0 = self.beta[bi0];
        let b1 = self.beta[bi1];

        // 1D weights along each axis – no heap allocation, no cloning.
        let (wm0, wm1) = if mi0 == mi1 {
            (1.0, 0.0)
        } else {
            let w = (mach - m0) / (m1 - m0);
            (1.0 - w, w)
        };
        let (wa0, wa1) = if ai0 == ai1 {
            (1.0, 0.0)
        } else {
            let w = (alpha_deg - a0) / (a1 - a0);
            (1.0 - w, w)
        };
        let (wb0, wb1) = if bi0 == bi1 {
            (1.0, 0.0)
        } else {
            let w = (beta_deg - b0) / (b1 - b0);
            (1.0 - w, w)
        };

        let idx_m = [mi0, mi1];
        let w_m = [wm0, wm1];
        let idx_a = [ai0, ai1];
        let w_a = [wa0, wa1];
        let idx_b = [bi0, bi1];
        let w_b = [wb0, wb1];

        let mut accum = AeroCoeff {
            cx: 0.0,
            cy: 0.0,
            cz: 0.0,
            cl: 0.0,
            cm: 0.0,
            cn: 0.0,
            cl_p: 0.0,
            cm_q: 0.0,
            cn_r: 0.0,
        };

        // Trilinear interpolation over the 2×2×2 cell around (mach, alpha, beta).
        for im in 0..2 {
            let w_mi = w_m[im];
            if w_mi == 0.0 {
                continue;
            }
            let m_idx = idx_m[im];
            for ia in 0..2 {
                let w_ai = w_a[ia];
                if w_ai == 0.0 {
                    continue;
                }
                let a_idx = idx_a[ia];
                for ib in 0..2 {
                    let w_bi = w_b[ib];
                    if w_bi == 0.0 {
                        continue;
                    }
                    let b_idx = idx_b[ib];
                    let w = w_mi * w_ai * w_bi;
                    let c = &self.coeffs[(m_idx * self.n_a + a_idx) * self.n_b + b_idx];
                    accum.cx += w * c.cx;
                    accum.cy += w * c.cy;
                    accum.cz += w * c.cz;
                    accum.cl += w * c.cl;
                    accum.cm += w * c.cm;
                    accum.cn += w * c.cn;
                    accum.cl_p += w * c.cl_p;
                    accum.cm_q += w * c.cm_q;
                    accum.cn_r += w * c.cn_r;
                }
            }
        }
        accum
    }
}

#[derive(Debug, Clone, Deserialize)]
#[allow(non_snake_case)]
pub struct RcsTorqueRow {
    pub time_s: f64,
    pub tau_x_Nm: f64,
    pub tau_y_Nm: f64,
    pub tau_z_Nm: f64,
    pub Fx_N: f64,
    pub Fy_N: f64,
    pub Fz_N: f64,
}

#[derive(Debug, Clone)]
pub struct RcsTorqueSeries {
    pub dt: f64,
    pub samples: Vec<RcsTorqueRow>,
}

impl RcsTorqueSeries {
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let mut rdr = csv::Reader::from_path(path)?;
        let mut samples = Vec::new();
        for rec in rdr.deserialize() {
            let row: RcsTorqueRow = rec?;
            samples.push(row);
        }
        samples.sort_by(|a, b| a.time_s.partial_cmp(&b.time_s).unwrap());
        let dt = if samples.len() >= 2 {
            samples[1].time_s - samples[0].time_s
        } else {
            0.01
        };
        Ok(RcsTorqueSeries { dt, samples })
    }

    pub fn sample(&self, t: f64) -> (Vector3<f64>, Vector3<f64>) {
        if self.samples.is_empty() {
            return (Vector3::zeros(), Vector3::zeros());
        }
        let idx = ((t - self.samples[0].time_s) / self.dt).round() as isize;
        if idx < 0 || idx as usize >= self.samples.len() {
            return (Vector3::zeros(), Vector3::zeros());
        }
        let s = &self.samples[idx as usize];
        (
            Vector3::new(s.tau_x_Nm, s.tau_y_Nm, s.tau_z_Nm),
            Vector3::new(s.Fx_N, s.Fy_N, s.Fz_N),
        )
    }
}
