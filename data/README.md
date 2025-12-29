# Data Files

Numerical ranges and trends reference NASA open-source literature on **blunt-body/lifting-entry vehicles** (Apollo/Orion class): hypersonic $C_D \sim 1.1\text{–}1.6$, $C_L$ increasing with angle of attack and peaking near $35^\circ\!-\!45^\circ$. The database uses $(M,\alpha,\beta)$ as independent variables to output $C_x,C_y,C_z, C_l,C_m,C_n$ and damping derivatives ($C_{l_p}, C_{m_q}, C_{n_r}$), along with **RCS layout and blackout-window pulse sequence** (100 Hz synthesized torque). These constructions and typical ranges align qualitatively/quantitatively with Apollo/Orion aerodynamic databases and wind-tunnel/flight data analyses. ([NASA Technical Reports Server][1])

### Files

* [vehicle\_properties.csv]
* [aero\_coeff\_table.csv]
  (Grid: Mach 2–12 (step ≈0.5), $\alpha$ −10°–50° (5°), $\beta$ −10°–10° (5°))
* [rcs\_layout.csv]
  (12 thrusters with position/direction/thrust/min pulse width/rise time)
* [rcs\_pulses.csv]
  (On/off pulse sequence during 120–160 s blackout window, with throttle ratio)
* [rcs\_torque\_timeseries\_100Hz.csv]
  (Aggregated 100 Hz synthesized torque and force)

---

### Data Structure and Physical Definitions

* **Aerodynamics (Coefficients)**

  * Input: `Mach, alpha_deg, beta_deg`
  * Output: `Cx, Cy, Cz, Cl, Cm, Cn, Cl_p, Cm_q, Cn_r`
  * Conventions: body-axes, forward $+x$, right $+y$, down $+z$. Lift up → $C_z<0$. Damping derivatives use standard non-dimensional rates $\hat p=\frac{pb}{2V},\,\hat q=\frac{qc}{2V},\,\hat r=\frac{rb}{2V}$.
  * Reference quantities: $S_\text{ref},\,b,\,c$ given in `vehicle_properties.csv`. Aero force/moment calculation:

    $$
    F = q_\infty S_\text{ref}\,[C_x,C_y,C_z],\quad
    M = q_\infty S_\text{ref}\,[b\,C_l,\,c\,C_m,\,b\,C_n].
    $$

* **RCS (Reaction Control System)**

  * `rcs_layout.csv`: body-frame position $(x,y,z)$, direction unit vector $\hat d$, max thrust, min pulse width, rise time for each thruster.
  * `rcs_pulses.csv`: `time_on_s, thruster_id, duration_s, throttle`.
  * `rcs_torque_timeseries_100Hz.csv`: aggregates thruster forces at 100 Hz, computing $\tau = r \times (F_\text{max}\cdot \text{throttle}\cdot \hat d)$ for net torque and force.

---

### Sources and References

* **Blunt-body hypersonic aerodynamics**: Cd magnitude, Cl growth with AoA, L/D range, and Orion/Apollo database construction methods can be verified in NASA open-source materials. ([NASA Technical Reports Server][1])
* **Wind-tunnel/flight comparison and uncertainty**: Apollo/Orion literature shows matching and uncertainty bands between wind-tunnel and flight data at high AoA and high $M$, which informed the Mach/AoA-dependent coefficients and damping derivative ranges used here. ([NASA Technical Reports Server][2])
* **General blunt-body trends**: For $M\!\gtrsim\!5$, Cd is relatively flat, Cl follows geometric amplitude ($\sin 2\alpha$) with a peak, and neutral sideslip derivatives are small—see NASA technical report surveys. ([NASA Technical Reports Server][3])

> **Disclaimer**: The above data is **synthetic data for research/simulation purposes only** and does not represent any manufacturer or mission configuration.

[1]: https://ntrs.nasa.gov/api/citations/20110013644/downloads/20110013644.pdf "Development of the Orion Crew Module Static ..."
[2]: https://ntrs.nasa.gov/api/citations/19670027745/downloads/19670027745.pdf "Entry flight aerodynamics from apollo mission AS-202"
[3]: https://ntrs.nasa.gov/api/citations/19650003700/downloads/19650003700.pdf "Aerodynamic Characteristics of Blunt Bodies"

---

# Extended Configuration Notes

- **Blackout aerodynamic mismatch scaling**: `config/sim.json` field `blackout.aero_mismatch_scale` (default 0.30) applies proportional scaling to aero coefficients during GNSS blackout to simulate model mismatch.

- **Out-of-order / dropout rate configuration matrix**: Each sensor config can optionally define a `schedule` field for time-phased policies, e.g.:

```json
"gnss": {
  "rate_hz": 10.0, "delay_s": 0.0, "jitter_std": 0.010,
  "noise": {"sigma": [2.5,2.5,2.5]},
  "outlier_rate": 0.05, "dropout_rate": 0.2, "ooo_rate": 0.05,
  "schedule": [
    {"t0": 0.0, "t1": 120.0, "dropout_rate": 0.1, "ooo_rate": 0.02},
    {"t0": 120.0, "t1": 160.0, "dropout_rate": 1.0, "ooo_rate": 0.10},
    {"t0": 160.0, "t1": 9999.0, "dropout_rate": 0.2, "ooo_rate": 0.05}
  ]
}
```

When `schedule` is not provided, the constant `dropout_rate/ooo_rate` values are used.
