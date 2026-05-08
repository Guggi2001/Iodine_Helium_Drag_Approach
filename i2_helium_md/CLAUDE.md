# CLAUDE.md

This project is a Python port of a legacy MATLAB molecular-dynamics codebase
for iodine / iodine-ion dynamics in helium nanodroplets.

Python package:

```text
i2_helium_md
```

Legacy MATLAB reference:

```text
legacy_matlab_repository/
```

## Current State

The transfer of the legacy MATLAB simulation code is mainly completed. The
neutral and ion propagation drivers are implemented, the main single-pulse
input presets are migrated, and MATLAB/Python cross-reference validation has
covered the important propagation and bookkeeping paths.

The remaining work is mainly post-processing plot recreation and refinement:
matching legacy MATLAB figures, separating plots that use different simulation
settings, and keeping those plotting scripts connected to the correct run
directories and reference data.

Do not broaden into new physics, pump-probe support, effusive dynamics, Abel
inversion, or full experimental VMI interpretation unless the user explicitly
asks.

## Implemented Core Path

Main public run script:

```text
scripts/run_single_pulse.py
```

Available input presets:

```text
single_pulse_N2000              # 9 A HeDFT comparison case
single_pulse_N2000_18Angst      # 18 A HeDFT comparison case
single_pulse_droplet_distribution
```

The run script writes self-contained run directories through `RunDirectory`:

```text
cfg.json
neutral.npz
ion.npz
```

Core implemented modules include:

```text
i2_helium_md/config.py
i2_helium_md/presets.py
i2_helium_md/physics/
i2_helium_md/sampling/
i2_helium_md/simulation/
i2_helium_md/postprocess/
```

## Current Plotting Scripts

The post-processing plots are now intentionally split because the legacy
combined MATLAB figure mixed outputs from different run settings.

HeDFT trajectory comparison:

```text
scripts/post_processing/plot_hedft_comparison.py
```

This script uses the HeDFT-comparison run and produces:

- I-I distance trajectory comparison,
- velocity-vs-time trajectory comparison.

Experimental velocity-distribution comparison:

```text
scripts/post_processing/plot_experimental_comparison.py
```

This script uses the realistic experimental-condition run and produces:

- gas-phase experimental VMI reference,
- droplet experimental VMI reference,
- simulated mass-selected final-velocity histograms for `I+He` and `I+He2`.

Keep these workflows separate. The top velocity-vs-time HeDFT comparison should
not be plotted from the experimental-condition run, and the experimental VMI
distribution should not be plotted from the HeDFT-comparison run.

Legacy MATLAB live-debug + paper figure reproduction:

```text
scripts/post_processing/plot_neutral_energy_balance.py
scripts/post_processing/plot_ion_energy_balance.py
scripts/post_processing/plot_ion_temperature_diagnostic.py
scripts/post_processing/plot_paper_figure.py
```

These run post-hoc from a finished `RunDirectory` and reproduce the
legacy figures from `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (energy
balance), `vmi_sim_3d_ion_propa.m` (energy balance + temperature
diagnostic), and the simulation-side panels of
`post_process_single_pulse_paper_v3.m` (radial v + phi histogram +
mass spectrum, exported as `compare_simulation_and_measurement.pdf`).
The polar-VMI panels (cos^2 anisotropy fit, beta(v), 3-D surf of
polar VMI image) are intentionally deferred -- the polar VMI image
data is not in `data/reference/` and CLAUDE.md scopes "full
experimental VMI interpretation" out by default.

## Post-Processing Data Contracts

Expected normalized reference inputs:

```text
data/reference/9A_All_Data.csv
data/reference/18A_All_Data.csv
data/reference/vmi_iplus_he.csv
data/reference/vmi_iplus_gas.csv
```

`9A_All_Data.csv` and `18A_All_Data.csv` use the 8-column header:

```text
Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance
```

The earlier split legacy 9 A files (`data_vabs2.csv`, `R1-R2.csv`) are
provenance, not the current Python loader contract.

Do not use hardcoded absolute MATLAB paths. Before adding or repurposing data
files, verify whether existing `data/reference/*.csv` files are experimental
VMI data or HeDFT trajectory data.

## Main Working Mode

Use the current Python package and run-directory APIs first. Prefer this
workflow:

1. Load an existing run with `RunDirectory`.
2. Load HeDFT references with `load_hedft_trajectory`.
3. Compute numerical diagnostics with `compare_distance` and
   `compare_velocity_magnitude` before changing trajectory plots.
4. Use `velocity_distribution.py` for VMI reference loading and mass-selected
   final-velocity histograms.
5. Keep plot changes local to `scripts/post_processing/` unless a package API
   change is actually needed.
6. Add or update focused pytest coverage when changing behavior.

The goal now is reliable reproduction of the legacy post-processing figures,
not broad refactoring or new simulation features.

## Relevant Files

Post-processing Python:

```text
i2_helium_md/postprocess/hedft_loader.py
i2_helium_md/postprocess/compare_trajectories.py
i2_helium_md/postprocess/velocity_distribution.py
i2_helium_md/postprocess/energy_balance.py
i2_helium_md/postprocess/_smoothing.py
i2_helium_md/postprocess/__init__.py
scripts/post_processing/plot_hedft_comparison.py
scripts/post_processing/plot_experimental_comparison.py
scripts/post_processing/plot_neutral_energy_balance.py
scripts/post_processing/plot_ion_energy_balance.py
scripts/post_processing/plot_ion_temperature_diagnostic.py
scripts/post_processing/plot_paper_figure.py
```

Simulation I/O:

```text
i2_helium_md/simulation/run_directory.py
i2_helium_md/simulation/checkpoint.py
```

Preset/run configuration:

```text
i2_helium_md/presets.py
scripts/run_single_pulse.py
```

Relevant MATLAB provenance:

```text
legacy_matlab_repository/inputfiles_dft_comparison/single_pulse_N2000.m
legacy_matlab_repository/inputfiles_dft_comparison/single_pulse_N2000_18Angst.m
legacy_matlab_repository/inputfiles_dft_comparison/single_pulse_droplet_distribution.m
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image_only_trajectories.m
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image.m
legacy_matlab_repository/vmi_sim_3d_neutral_propa_HeDFT_mimic.m
legacy_matlab_repository/vmi_sim_3d_ion_propa.m
legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m
```

Do not read the full legacy plotting stack unless a specific numerical or
visual discrepancy requires it.

## Known Plotting Details

For the velocity-vs-time HeDFT panel, the MATLAB code sampled roughly 15
molecules and plotted both iodine atoms, so Python caps the overlay near 30
velocity traces.

For the experimental velocity-distribution panel, MATLAB used:

```text
edges_velocity = 0:0.04:26
vd_ion = movmean(h, 15)
xlim([0, 28])
```

Python should preserve that scale when matching the legacy figure: fine
`0.04 A/ps` bins, 15-bin moving mean on simulation histograms, and a displayed
velocity range out to `28 A/ps`.

## Editing Limits

Before editing, run:

```bash
git status
```

If the working tree is dirty, do not revert or overwrite unrelated changes.
Proceed only if the user explicitly scopes the edit or the dirty files are your
own current work.

Avoid changing these unless a focused test reveals a real bug:

```text
checkpoint schema
physical constants
neutral propagation physics
ion propagation physics
collision physics
single-pulse run script behavior
```

For current work, prefer changes in:

```text
scripts/post_processing/
i2_helium_md/postprocess/
focused tests
status/docs files
```

## Testing Expectations

At minimum for post-processing changes:

- loader tests for normalized HeDFT CSV format,
- comparison tests using tiny synthetic `IonCheckpoint` objects,
- validation for missing files and non-overlapping time axes,
- VMI reference loader and mass-filtered histogram tests,
- plotting smoke coverage with `matplotlib` in non-interactive mode.

For preset/run-script changes, run the focused preset and run-script tests.

Do not generate figures or production-sized checkpoints in tests.

## Reporting Format

For post-processing work, report:

1. MATLAB files inspected,
2. Python files changed,
3. data files or run directories used,
4. numerical or plotting behavior changed,
5. tests run and results,
6. remaining risks or deferred behavior.

After a focused post-processing change is complete, stop and report. Do not
automatically move to Abel inversion, full experimental VMI interpretation, or
new physics branches.
