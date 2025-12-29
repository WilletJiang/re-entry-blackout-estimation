use crate::sim::{BlackoutTrace, BlackoutTraceConfig, MonteCarloResult};
use plotters::element::PathElement;
use plotters::prelude::*;
use std::fs::File;
use std::io::{Write, BufWriter};

pub fn plot_nees_nis(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let root = SVGBackend::new(path, (900, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let max_run = results.len() as f64;
    let max_val = results
        .iter()
        .map(|r| r.nees_mean.max(r.nis_mean))
        .fold(0.0f64, f64::max)
        .max(15.0);
    let mut chart = ChartBuilder::on(&root)
        .caption("NEES/NIS per run", ("sans-serif", 24))
        .margin(15)
        .x_label_area_size(35)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..max_run, 0.0..max_val)?;
    chart.configure_mesh().x_desc("run").y_desc("value").draw()?;
    // 95% chi2 lines for dof 6 and 3
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 12.592), (max_run, 12.592)],
        &BLUE.mix(0.3),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 7.815), (max_run, 7.815)],
        &RED.mix(0.3),
    )))?;
    chart
        .draw_series(LineSeries::new(
            results.iter().enumerate().map(|(i, r)| (i as f64, r.nees_mean)),
            &BLUE,
        ))?
        .label("NEES mean")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    chart
        .draw_series(LineSeries::new(
            results.iter().enumerate().map(|(i, r)| (i as f64, r.nis_mean)),
            &RED,
        ))?
        .label("NIS mean")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    Ok(())
}

pub fn plot_histogram(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let root = SVGBackend::new(path, (900, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let vals: Vec<f64> = results.iter().map(|r| r.rms_pos).collect();
    let max_v = vals.iter().cloned().fold(0.0f64, f64::max).max(1.0);
    let bins = 30usize;
    let bin_w = max_v / bins as f64;
    let mut counts = vec![0u32; bins];
    for v in vals {
        let idx = ((v / max_v).min(0.9999) * bins as f64) as usize;
        counts[idx] += 1;
    }
    let y_max = counts.iter().cloned().max().unwrap_or(1) + 1;
    let mut chart = ChartBuilder::on(&root)
        .caption("RMS Position Error Histogram (3σ≤40 m target)", ("sans-serif", 24))
        .margin(15)
        .x_label_area_size(35)
        .y_label_area_size(50)
        .build_cartesian_2d(0f64..max_v, 0u32..y_max)?;
    chart.configure_mesh().x_desc("RMS pos (m)").y_desc("count").draw()?;
    // acceptance line at 40/3 ≈ 13.33 m (RMS corresponding to 3σ 40 m)
    let accept = 40.0 / 3.0;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(accept, 0u32), (accept, y_max as u32)],
        &BLUE.mix(0.4),
    )))?;
    chart.draw_series(
        (0..bins).map(|i| {
            let x0 = i as f64 * bin_w;
            let x1 = x0 + bin_w;
            let y = counts[i];
            Rectangle::new([(x0, 0u32), (x1, y)], RED.mix(0.5).filled())
        }),
    )?;
    Ok(())
}

// Bias RMS histograms (gyro in deg/s, accel in µg)
pub fn plot_bias_histograms(results: &[MonteCarloResult], path_bg: &str, path_ba: &str) -> anyhow::Result<()> {
    // Gyro bias (deg/s)
    {
        let vals: Vec<f64> = results.iter().map(|r| r.bias_bg_rms_dps).collect();
        let root = SVGBackend::new(path_bg, (900, 500)).into_drawing_area();
        root.fill(&WHITE)?;
        let max_v = vals.iter().cloned().fold(0.0f64, f64::max).max(0.02);
        let bins = 30usize;
        let bin_w = max_v / bins as f64;
        let mut counts = vec![0u32; bins];
        for v in vals {
            let idx = ((v / max_v).min(0.9999) * bins as f64) as usize;
            counts[idx] += 1;
        }
        let y_max = counts.iter().cloned().max().unwrap_or(1) + 1;
        let mut chart = ChartBuilder::on(&root)
            .caption("Gyro bias RMS (deg/s) histogram (≤0.005 target)", ("sans-serif", 24))
            .margin(15)
            .x_label_area_size(35)
            .y_label_area_size(50)
            .build_cartesian_2d(0f64..max_v, 0u32..y_max)?;
        chart.configure_mesh().x_desc("deg/s").y_desc("count").draw()?;
        let accept = 0.005;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(accept, 0u32), (accept, y_max as u32)],
            &BLUE.mix(0.4),
        )))?;
        chart.draw_series(
            (0..bins).map(|i| {
                let x0 = i as f64 * bin_w;
                let x1 = x0 + bin_w;
                let y = counts[i];
                Rectangle::new([(x0, 0u32), (x1, y)], RED.mix(0.5).filled())
            }),
        )?;
    }
    // Accel bias (µg)
    {
        let vals: Vec<f64> = results.iter().map(|r| r.bias_ba_rms_ug).collect();
        let root = SVGBackend::new(path_ba, (900, 500)).into_drawing_area();
        root.fill(&WHITE)?;
        let max_v = vals.iter().cloned().fold(0.0f64, f64::max).max(300.0);
        let bins = 30usize;
        let bin_w = max_v / bins as f64;
        let mut counts = vec![0u32; bins];
        for v in vals {
            let idx = ((v / max_v).min(0.9999) * bins as f64) as usize;
            counts[idx] += 1;
        }
        let y_max = counts.iter().cloned().max().unwrap_or(1) + 1;
        let mut chart = ChartBuilder::on(&root)
            .caption("Accel bias RMS (µg) histogram (≤120 target)", ("sans-serif", 24))
            .margin(15)
            .x_label_area_size(35)
            .y_label_area_size(50)
            .build_cartesian_2d(0f64..max_v, 0u32..y_max)?;
        chart.configure_mesh().x_desc("µg").y_desc("count").draw()?;
        let accept = 120.0;
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(accept, 0u32), (accept, y_max as u32)],
            &BLUE.mix(0.4),
        )))?;
        chart.draw_series(
            (0..bins).map(|i| {
                let x0 = i as f64 * bin_w;
                let x1 = x0 + bin_w;
                let y = counts[i];
                Rectangle::new([(x0, 0u32), (x1, y)], RED.mix(0.5).filled())
            }),
        )?;
    }
    Ok(())
}

// Pass-rate plot (single line showing pass percentage)
pub fn plot_pass_rate(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let root = SVGBackend::new(path, (800, 200)).into_drawing_area();
    root.fill(&WHITE)?;
    let total = results.len() as f64;
    let passed = results.iter().filter(|r| r.passed).count() as f64;
    let rate = if total > 0.0 { passed / total * 100.0 } else { 0.0 };
    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Pass rate: {:.1}%", rate), ("sans-serif", 24))
        .margin(15)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0.0..100.0, 0.0..1.0)?;
    chart.configure_mesh().x_desc("percent").disable_y_mesh().y_labels(0).draw()?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, 0.5), (rate, 0.5)],
        &GREEN,
    )))?;

    chart.draw_series(std::iter::once(PathElement::new(
        vec![(95.0, 0.0), (95.0, 1.0)],
        &BLUE.mix(0.5),
    )))?;
    Ok(())
}

pub fn write_blackout_envelope_csv(
    traces: &[BlackoutTrace],
    cfg: &BlackoutTraceConfig,
    path: &str,
) -> anyhow::Result<()> {
    if traces.is_empty() {
        return Ok(());
    }
    let n = cfg.sample_count();
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "t_s,pos_err_p50_m,pos_err_p95_m,pos_sigma3_mean_m,nees_rv_mean"
    )?;

    for i in 0..n {
        let t = cfg.time_at(i);
        let (p50, p95, sigma3, nees) = summarize_trace_sample(traces, i);
        writeln!(w, "{:.6},{:.6},{:.6},{:.6},{:.6}", t, p50, p95, sigma3, nees)?;
    }
    Ok(())
}

pub fn plot_blackout_envelope(
    traces: &[BlackoutTrace],
    cfg: &BlackoutTraceConfig,
    path: &str,
) -> anyhow::Result<()> {
    if traces.is_empty() {
        return Ok(());
    }

    let n = cfg.sample_count();
    let mut t_vec = Vec::with_capacity(n);
    let mut p50_vec = Vec::with_capacity(n);
    let mut p95_vec = Vec::with_capacity(n);
    let mut sigma3_vec = Vec::with_capacity(n);
    for i in 0..n {
        t_vec.push(cfg.time_at(i));
        let (p50, p95, sigma3, _) = summarize_trace_sample(traces, i);
        p50_vec.push(p50);
        p95_vec.push(p95);
        sigma3_vec.push(sigma3);
    }

    let mut y_max = 1.0f64;
    for i in 0..n {
        y_max = y_max.max(p95_vec[i]).max(sigma3_vec[i]);
    }
    y_max *= 1.2;

    let root = SVGBackend::new(path, (900, 450)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Blackout crossing: position error envelope vs covariance",
            ("sans-serif", 22),
        )
        .margin(15)
        .x_label_area_size(35)
        .y_label_area_size(60)
        .build_cartesian_2d(cfg.t_start..cfg.t_end, 0.0..y_max)?;
    chart.configure_mesh().x_desc("t (s)").y_desc("m").draw()?;

    // blackout markers
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.blackout_start, 0.0), (cfg.blackout_start, y_max)],
        &BLACK.mix(0.25),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.blackout_end, 0.0), (cfg.blackout_end, y_max)],
        &BLACK.mix(0.25),
    )))?;

    chart
        .draw_series(LineSeries::new(
            (0..n).map(|i| (t_vec[i], p95_vec[i])),
            &RED,
        ))?
        .label("pos err p95")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], &RED));
    chart
        .draw_series(LineSeries::new(
            (0..n).map(|i| (t_vec[i], p50_vec[i])),
            &BLUE,
        ))?
        .label("pos err p50")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], &BLUE));
    chart
        .draw_series(LineSeries::new(
            (0..n).map(|i| (t_vec[i], sigma3_vec[i])),
            &GREEN,
        ))?
        .label("mean 3σ (trace(P_rr)/3)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], &GREEN));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;
    Ok(())
}

pub fn plot_blackout_nees_time(
    traces: &[BlackoutTrace],
    cfg: &BlackoutTraceConfig,
    path: &str,
) -> anyhow::Result<()> {
    if traces.is_empty() {
        return Ok(());
    }

    let n = cfg.sample_count();
    let mut t_vec = Vec::with_capacity(n);
    let mut nees_mean = Vec::with_capacity(n);
    for i in 0..n {
        t_vec.push(cfg.time_at(i));
        let (_, _, _, nees) = summarize_trace_sample(traces, i);
        nees_mean.push(nees);
    }

    let mut y_max = 15.0f64;
    for v in &nees_mean {
        y_max = y_max.max(*v);
    }
    y_max *= 1.1;

    let root = SVGBackend::new(path, (900, 450)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Blackout crossing: NEES(RV) vs time", ("sans-serif", 22))
        .margin(15)
        .x_label_area_size(35)
        .y_label_area_size(60)
        .build_cartesian_2d(cfg.t_start..cfg.t_end, 0.0..y_max)?;
    chart
        .configure_mesh()
        .x_desc("t (s)")
        .y_desc("NEES (dof=6)")
        .draw()?;

    // blackout markers
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.blackout_start, 0.0), (cfg.blackout_start, y_max)],
        &BLACK.mix(0.25),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.blackout_end, 0.0), (cfg.blackout_end, y_max)],
        &BLACK.mix(0.25),
    )))?;

    // 95% chi2 interval for dof=6
    let chi2_lo = 1.237_344_245_791_202_5_f64;
    let chi2_hi = 14.449_375_335_447_918_f64;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.t_start, chi2_lo), (cfg.t_end, chi2_lo)],
        &RED.mix(0.35),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.t_start, chi2_hi), (cfg.t_end, chi2_hi)],
        &RED.mix(0.35),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cfg.t_start, 6.0), (cfg.t_end, 6.0)],
        &BLUE.mix(0.25),
    )))?;

    chart.draw_series(LineSeries::new(
        (0..n).map(|i| (t_vec[i], nees_mean[i])),
        &BLACK,
    ))?;
    Ok(())
}

fn summarize_trace_sample(traces: &[BlackoutTrace], sample_idx: usize) -> (f64, f64, f64, f64) {
    let mut errs: Vec<f32> = Vec::with_capacity(traces.len());
    let mut sigma3_sum = 0.0f64;
    let mut sigma3_ct = 0usize;
    let mut nees_sum = 0.0f64;
    let mut nees_ct = 0usize;

    for t in traces {
        if let Some(v) = t.pos_err_m.get(sample_idx) {
            if v.is_finite() {
                errs.push(*v);
            }
        }
        if let Some(v) = t.pos_sigma3_m.get(sample_idx) {
            if v.is_finite() {
                sigma3_sum += *v as f64;
                sigma3_ct += 1;
            }
        }
        if let Some(v) = t.nees_rv.get(sample_idx) {
            if v.is_finite() {
                nees_sum += *v as f64;
                nees_ct += 1;
            }
        }
    }

    errs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let p50 = quantile_sorted(&errs, 0.50) as f64;
    let p95 = quantile_sorted(&errs, 0.95) as f64;
    let sigma3_mean = if sigma3_ct > 0 {
        sigma3_sum / (sigma3_ct as f64)
    } else {
        f64::NAN
    };
    let nees_mean = if nees_ct > 0 {
        nees_sum / (nees_ct as f64)
    } else {
        f64::NAN
    };
    (p50, p95, sigma3_mean, nees_mean)
}

fn quantile_sorted(vals: &[f32], q: f64) -> f32 {
    if vals.is_empty() {
        return f32::NAN;
    }
    let q = q.clamp(0.0, 1.0);
    let idx = ((vals.len() - 1) as f64 * q).round() as usize;
    vals[idx.min(vals.len() - 1)]
}


pub fn write_summary_csv(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "run,rms_pos_m,rms_vel_mps,rms_att_deg,bias_bg_rms_dps,bias_ba_rms_ug,nees_mean,nis_mean,nees_covered,nis_covered,sm_rms_pos,sm_rms_vel,passed"
    )?;
    for (i, r) in results.iter().enumerate() {
        writeln!(
            w,
            "{},{:.6},{:.6},{:.6},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{}",
            i,
            r.rms_pos,
            r.rms_vel,
            r.rms_att_deg,
            r.bias_bg_rms_dps,
            r.bias_ba_rms_ug,
            r.nees_mean,
            r.nis_mean,
            r.nees_covered,
            r.nis_covered,
            r.smoothed_rms_pos,
            r.smoothed_rms_vel,
            if r.passed { 1 } else { 0 }
        )?;
    }
    Ok(())
}

pub fn write_post_window_csv(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "run,post_gnss_updates,post_gnss_accepted,post_accept_rate,post_sat_mean,post_sat_min"
    )?;
    for (i, r) in results.iter().enumerate() {
        let rate = if r.post_gnss_updates > 0 {
            (r.post_gnss_accepted as f64) / (r.post_gnss_updates as f64)
        } else {
            f64::NAN
        };
        writeln!(
            w,
            "{},{},{},{:.6},{:.6},{}",
            i,
            r.post_gnss_updates,
            r.post_gnss_accepted,
            rate,
            r.post_sat_mean,
            r.post_sat_min
        )?;
    }
    Ok(())
}

pub fn write_nees_blocks_csv(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "run,nees_r_mean,nees_v_mean,nees_ba_mean,nees_r_covered,nees_v_covered,nees_ba_covered"
    )?;
    for (i, r) in results.iter().enumerate() {
        writeln!(
            w,
            "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            i,
            r.nees_r_mean,
            r.nees_v_mean,
            r.nees_ba_mean,
            r.nees_r_covered,
            r.nees_v_covered,
            r.nees_ba_covered,
        )?;
    }
    Ok(())
}

pub fn write_nees_phase_csv(results: &[MonteCarloResult], path: &str) -> anyhow::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "run,nees_mean_pre,nees_mean_blk,nees_mean_post,nees_mean_tail,nees_covered_pre,nees_covered_blk,nees_covered_post,nees_covered_tail,nis_mean_pre,nis_mean_blk,nis_mean_post,nis_mean_tail,nis_covered_pre,nis_covered_blk,nis_covered_post,nis_covered_tail"
    )?;
    for (i, r) in results.iter().enumerate() {
        writeln!(
            w,
            "{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
            i,
            r.nees_mean_pre,
            r.nees_mean_blk,
            r.nees_mean_post,
            r.nees_mean_tail,
            r.nees_covered_pre,
            r.nees_covered_blk,
            r.nees_covered_post,
            r.nees_covered_tail,
            r.nis_mean_pre,
            r.nis_mean_blk,
            r.nis_mean_post,
            r.nis_mean_tail,
            r.nis_covered_pre,
            r.nis_covered_blk,
            r.nis_covered_post,
            r.nis_covered_tail,
        )?;
    }
    Ok(())
}
