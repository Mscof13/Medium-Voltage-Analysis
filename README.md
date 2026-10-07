# Medium-Voltage Line Short-Circuit Experiment

Modeling of a 35 kV medium-voltage line, first at no load and then during a short circuit, driven by three different sources. Each case is solved analytically, solved numerically in MATLAB, and simulated in Simulink with both Simulink and Simscape components.

Course project "Application of Computer Tools in Power Engineering", Faculty of Technical Sciences, University of Novi Sad (2025).

**Authors:** Nebojša Todorović, Nikola Pantelić, Đorđe Dimitrić

> The report is written in Serbian (`docs/report_sr.pdf`). This README summarizes it in English.

## Line model

A lumped-parameter model (series R and L, shunt C and G) of an aluminium overhead line.

| Quantity | Value |
|---|---|
| Nominal voltage | 35 kV (RMS) |
| Rated power | 12.1 MVA, cos φ = 0.9 |
| Conductor | Aluminium, 240 mm² |
| Length | 19.1 km |
| Per-km parameters | r = 0.1175 Ω/km, l = 0.884 mH/km, c = 0.02 µF/km, g = 1 µS/km |
| Line totals | R = 2.24 Ω, L = 16.9 mH, C = 0.382 µF, G = 19.1 µS |

## The three cases

### 1. Constant current source (I0 = 200 A)

At no load the line is described by a first-order differential equation. Because an ideal current source keeps pushing current into the open line, the voltage rises toward I0/G ≈ 10.5 MV, an idealized value. After the short circuit the receiving-end voltage is zero and the current equals I0. The analytical solution, the MATLAB numerical solution, and the Simulink and Simscape models all agree.

### 2. Sinusoidal voltage source (35 kV RMS, 50 Hz)

At no load the line is a second-order system with an underdamped transient. The steady-state receiving-end voltage is about 36.15 kV RMS, higher than the source voltage. This is the Ferranti effect, which appears on lightly loaded high-voltage lines.

For the short circuit, the current is the sum of a decaying transient and a steady-state part. The size of the short-circuit current depends on the instant of switching, so the report derives the instants that give the smallest and largest current. The transient dies out after about five time constants L/R.

### 3. Piecewise-linear voltage source

Source parameters: U0 = 10, t1 = 1 s, t2 = 2 s, t3 = 4 s. The analytical solution uses the Laplace transform and Heaviside step functions. It is compared with the MATLAB numerical solution (`ode23`) and with the Simulink and Simscape models, and two peaks are examined in detail.

## Files

- `part1_constant_current.m`
- `part2_sinusoidal_voltage.m`
- `part3_arbitrary_voltage.m`
- `docs/report_sr.pdf`: full report (in Serbian)

## Simulink / Simscape models

The `pod_a/b/c` in the file names stands for subtask a/b/c of the project (Serbian "podzadatak"), one per source type:

| Subtask | Source | Simulink | Simscape |
| ------- | ------ | -------- | -------- |
| a | Constant current (I0 = 200 A) | `Projekat_pod_a_simulink.slx` | `Projekat_pod_a_simscape.slx` |
| b | Sinusoidal voltage (35 kV RMS, 50 Hz) | `Projekat_pod_b_simulink.slx` | `Projekat_pod_b_simscape.slx` |
| c | Piecewise-linear voltage | `Projekat_pod_c_simulink.slx` | `Projekat_pod_c_simscape.slx` |

Saved in MATLAB R20XX.

## Tools

MATLAB (`ode45`, `ode23`, Symbolic Math Toolbox), Simulink, Simscape
