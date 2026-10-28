Radiation Dose Accounting: MATLAB Reproduction and Bin-Width Sensitivity

Overview

This repository contains an independent MATLAB reproduction and synthetic sensitivity analysis of time-binned radiation dose accounting.

The work was conducted in response to a research task proposed by Abhimanyu Singh through the NASA Open Science Data Repository (OSDR) Analysis Working Groups (AWG) community.

The objective is to investigate how bin width, pulse timing, and partial-window allocation affect cumulative radiation dose estimates.

Methods

1. Reference Reproduction

Four synthetic reference cases were reproduced using duration-weighted dose integration:

Case	Reference Dose (µGy)
Constant	5/9
Linear	5/9
Adaptive pulse	19/36
Narrow pulse	5/18

2. Pulse and Bin-Width Sensitivity

A rectangular pulse was defined with:

* Dose rate: 100 µGy/h
* Pulse duration: 10 seconds
* Pulse start times: 40, 45, and 50 seconds
* Analysis interval: [0, 100) seconds
* Uniform bin widths: 5, 10, 20, 50, and 100 seconds
* Window endpoints: 48, 52, and 100 seconds

Exact bin means were calculated using pulse-bin overlap durations.

For partial bins, proportional allocation was evaluated against the analytically known pulse overlap.

Results

All four reference totals were reproduced.

Across 45 pulse/bin/window combinations, the full-interval dose was conserved at 5/18 µGy.

Partial-window proportional allocation produced absolute errors of up to 0.133333 µGy.

The error was not necessarily monotonic with bin width and depended on the relative alignment of the pulse, bin boundaries, and window endpoint.

Repository Contents

* run_interval_study.m — MATLAB analysis script
* control_reproduction.csv — Reference reproduction results
* pulse_bin_error_table.csv — Sensitivity analysis results
* pulse_bin_error.png — Absolute error visualization

Reproducibility

Run run_interval_study.m in MATLAB.

The script reproduces the synthetic reference cases, evaluates the pulse/bin-width configurations, and generates the result tables and figure.

Scope and Limitations

This study uses synthetic signals with analytically known reference totals.

The observed partial-window errors arise from proportional allocation assumptions applied to time-averaged bins. They should not be interpreted as demonstrated errors in physical radiation instruments.

Real DosTel measurement intervals were not analyzed in this initial study.

Attribution

Original research problem, reference implementation, and synthetic test specifications: Abhimanyu Singh, Space Radiation Measurement Reliability project.

Independent MATLAB reproduction and sensitivity analysis: Yeongho Lee.

This repository is an individual analysis contribution and is not an official NASA software release.
