# GN&C Task Specification: SE-2317 "Entry Blackout State Estimation"

## Mission Context

* Vehicle: Reusable orbital-class stage, **6-DOF** rigid-body dynamics, belly-first blunt-body re-entry, **RCS** attitude control + movable flaps.
* Conditions: **Plasma blackout** (GNSS & comm loss) around peak heating, strongly coupled nonlinear aerodynamics ($M\in[2,12]$, high $q_\mathrm{dyn}$), thermal-induced **IMU drift** and **scale-factor nonlinearity**.
* Objective: Vehicle must enter a controllable guidance corridor within **≤ 2 s** after blackout window exit, with state error and covariance consistency meeting specified constraints.

---

## Sensor Configuration (Asynchronous Multi-Sensor)

* **IMU** (400 Hz): 3-axis gyro & accelerometer; noise (1σ): $\sigma_\omega=0.01^\circ/\mathrm{s}/\sqrt{\mathrm{Hz}}$, $\sigma_a=150\,\mu g/\sqrt{\mathrm{Hz}}$; **thermal drift** as random walk.
* **GNSS** (10 Hz): PVT solution; **complete signal loss** in $[t_0+120\,\mathrm{s},\, t_0+160\,\mathrm{s}]$.
* **Star Tracker** (2 Hz): quaternion attitude; **fixed delay** $\delta t_\star=120\,\mathrm{ms}$; sun glare **saturation** (data loss) in $[t_0+130,\, t_0+150]$.
* **Radar Altimeter** (20 Hz): terrain-normal altitude, **valid only below $h<25\,\mathrm{km}$**; delay $\delta t_h=70\,\mathrm{ms}$.
* **Measurement anomalies**: 5% probability of **heavy-tail/outlier** (Gaussian mixture: 95% $\mathcal N$, 5% Cauchy); data may arrive **out-of-order**.

---

## State Vector Definition

$$
x = \big[r,\ v,\ q,\ \omega,\ b_g,\ b_a\big]^\top
$$

* $r\in\mathbb R^3$ position; $v\in\mathbb R^3$ velocity; $q\in\mathbb H$ unit quaternion; $\omega\in\mathbb R^3$ angular velocity;
* $b_g, b_a$ are gyro/accelerometer biases (random walk).

**Continuous-time model:**

$$
\begin{aligned}
\dot r &= v \\
\dot v &= R(q)\big(a_\text{meas}-b_a-n_a\big) + g(r) + a_\text{aero}(r,v) + a_\text{RCS}(t) \\
\dot q &= \tfrac{1}{2}\,\Omega(\omega)\,q \\
\dot \omega &= J^{-1}\!\left(\tau_\text{meas}-b_g-n_g - \omega\times J\omega + \tau_\text{aero}(r,v) + \tau_\text{RCS}(t)\right) \\
\dot b_g &= w_g,\quad \dot b_a = w_a
\end{aligned}
$$

* $R(q)$: direction cosine matrix (body-to-inertial); $\Omega(\omega)$: quaternion rate matrix;
* $a_\text{aero}, \tau_\text{aero}$: **unknown/time-varying aerodynamic residuals** (primary source of model mismatch, strongly nonlinear with $M,\,q_\mathrm{dyn}$);
* $w_g,w_a$: zero-mean Gaussian process noise (random walk).

**Measurement model:**

$$
\begin{aligned}
z_\mathrm{GNSS}(t_k) &= h_\mathrm{GNSS}(x(t_k)) + \nu_g \\
z_\star(t_k-\delta t_\star) &= h_\star(x(t_k-\delta t_\star)) + \nu_\star \\
z_h(t_k-\delta t_h) &= h_h(x(t_k-\delta t_h)) + \nu_h
\end{aligned}
$$

* Note: **delayed observations** and **out-of-order arrival** require time calibration and **stochastic cloning / delay compensation** mechanisms.

---

## Quantitative Targets

Within **2 s after blackout exit** ($t_\mathrm{exit}=t_0+160\,\mathrm{s}$), including fixed-lag smoothing:

1. **Position error (3σ) ≤ 40 m**; **velocity error (3σ) ≤ 0.5 m/s**.
2. **Attitude error (3σ) ≤ 0.15°** (equivalent cone angle); **angular rate (3σ) ≤ 0.05°/s**.
3. **Bias convergence**: $\|b_g\|_{3\sigma} \le 0.005^\circ/\mathrm{s}$, $\|b_a\|_{3\sigma} \le 120\,\mu g$.
4. **Consistency**: NEES/NIS 95% confidence interval coverage ≥ 90% over 100 Monte Carlo runs.
5. **Robustness**: No divergence under 5% outlier + 20% random dropout stress test; covariance not overly optimistic (underestimate by no more than 30%).
6. **Real-time**: **100 Hz main loop**, **single-cycle compute budget ≤ 5 ms** (ARM A78AE single-core equivalent), **memory ≤ 20 MB**; support **fixed-lag smoothing window ≤ 10 s**.
