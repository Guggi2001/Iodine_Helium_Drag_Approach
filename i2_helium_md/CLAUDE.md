# CLAUDE.md

This project is a Python port of a legacy MATLAB molecular-dynamics codebase for iodine / iodine-ion dynamics in helium nanodroplets.

Python package:

```text
i2_helium_md
```

Legacy MATLAB reference:

```text
legacy_matlab_repository/
```

## Current task

The neutral and ion propagation drivers are implemented, ion-stage
MATLAB/Python cross-reference validation is complete, and the public
single-pulse run script is implemented:

```text
scripts/run_single_pulse.py
```

Step 13 now has a first implemented post-processing path in `postprocess/`.
The current task is to validate, document, and carefully extend that path only
where needed.

Implemented Step 13 files:

```text
i2_helium_md/postprocess/hedft_loader.py
i2_helium_md/postprocess/compare_trajectories.py
i2_helium_md/postprocess/velocity_distribution.py
scripts/plot_hedft_comparison.py
scripts/post_processing_comparison/compare.py
```

`scripts/post_processing_comparison/compare.py` is imported VMI-reference
verification context. The package-level APIs in `i2_helium_md/postprocess/`
are the main code path for HeDFT loading, trajectory comparison, VMI reference
loading, and final-velocity histogram calculations.

Do not continue with Abel inversion, full experimental VMI interpretation,
pump-probe support, effusive dynamics, analytical-force cleanup, broad
refactors, or out-of-scope MATLAB branches unless the user explicitly asks.

## Main working mode

Work data-contract first, and keep post-processing changes narrow.

Prefer this workflow when continuing from the current state:

1. Load the existing production run through `RunDirectory`.
2. Load `data/reference/9A_All_Data.csv` with `load_hedft_trajectory`.
3. Compute explicit numerical diagnostics with `compare_distance` and
   `compare_velocity_magnitude`.
4. Use `velocity_distribution.py` only for the current VMI reference overlay
   and mass-selected final-velocity histograms.
5. Add or update focused pytest coverage before changing public behavior.
6. Only then adjust `scripts/plot_hedft_comparison.py` if the numerical API
   proves the plot needs a change.

The goal is reliable loading, numerical comparison, and a reproducible first
plotting path. Avoid turning this into final publication plotting or full VMI
analysis unless requested.

## Expected Step 13 outputs

Implemented files:

```text
i2_helium_md/postprocess/hedft_loader.py
i2_helium_md/postprocess/compare_trajectories.py
i2_helium_md/postprocess/velocity_distribution.py
i2_helium_md/postprocess/__init__.py
scripts/plot_hedft_comparison.py
tests/test_hedft_loader.py
tests/test_compare_trajectories.py
tests/test_velocity_distribution.py
tests/test_plot_hedft_comparison_smoke.py
docs/hedft_loader_module.md
docs/compare_trajectories_module.md
docs/velocity_distribution_module.md
docs/plot_hedft_comparison_script.md
```

The trajectory comparison computes:

- mean MD I-I distance vs. time from the ion checkpoint,
- interpolation of MD distance onto the normalized HeDFT time grid,
- overlap interval,
- RMSE in Å,
- mean MD/HeDFT distance ratio.

It also computes analogous I1/I2 velocity-magnitude comparisons. The VMI
helper currently loads `vmi_iplus_he.csv` and `vmi_iplus_gas.csv` and computes
plain mass-selected final-velocity histograms for simulation overlays.

This mirrors the numerical block in `simulation_image_only_trajectories.m`:

```text
dR_mean = mean(sqrt(dx^2 + dy^2 + dz^2), axis=particles)
tmax = min(max(tR), max(time_i))
dR_mean_on_tR = interp1(t_md, d_md, tR_use)
rmse = sqrt(mean((dR_mean_on_tR - R_use)^2))
ratio = mean(dR_mean_on_tR / R_use)
```

## Data policy

Do not use hardcoded absolute MATLAB paths.

Expected normalized reference inputs are:

```text
data/reference/9A_All_Data.csv
data/reference/18A_All_Data.csv
data/reference/vmi_iplus_he.csv
data/reference/vmi_iplus_gas.csv
```

`9A_All_Data.csv` and `18A_All_Data.csv` use the 8-column header
`Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`. The earlier split
legacy 9 Å files (`data_vabs2.csv`, `R1-R2.csv`) are provenance, not the
current Python loader contract.

Before adding data files, verify whether existing `data/reference/*.csv` files
are experimental VMI data or HeDFT trajectory data. Do not silently repurpose
files with ambiguous names.

## Relevant files

Start with these Python files:

```text
i2_helium_md/postprocess/hedft_loader.py
i2_helium_md/postprocess/compare_trajectories.py
i2_helium_md/postprocess/velocity_distribution.py
scripts/plot_hedft_comparison.py
scripts/post_processing_comparison/compare.py
```

Then inspect simulation I/O only as needed:

```text
i2_helium_md/simulation/run_directory.py
i2_helium_md/simulation/checkpoint.py
```

Relevant MATLAB files:

```text
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image_only_trajectories.m
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image.m
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/importfile_v2.m
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/importfile_R1_R2.m
```

Do not read the full legacy plotting stack unless a numerical comparison
depends on it.

## Editing limits

Before editing, run:

```bash
git status
```

If the working tree is dirty, do not revert or overwrite unrelated changes.
Proceed only if the user explicitly scopes the edit or the dirty files are your
own current work. Otherwise stop and ask.

During post-processing work, do not change unless a focused comparison test
reveals a real bug:

```text
checkpoint schema
physical constants
neutral propagation physics
ion propagation physics
collision physics
single-pulse run script behavior
```

Keep code changes local to `postprocess/`, `scripts/plot_hedft_comparison.py`,
focused tests, small reference CSV copies if needed, and status documentation.

## Testing expectations

At minimum:

- loader tests for the normalized HeDFT CSV format,
- comparison tests using a tiny synthetic `IonCheckpoint`,
- validation for missing files and non-overlapping time axes,
- VMI reference loader and mass-filtered histogram tests,
- plotting smoke coverage with `matplotlib` in non-interactive mode,
- relevant existing checkpoint/run-directory tests if loader code uses them.

Do not generate figures or production-sized checkpoints in tests.

## Reporting format

For Step 13 work, report:

1. data files inspected or copied,
2. MATLAB files inspected,
3. Python files written or changed,
4. numerical comparison API,
5. tests run and results,
6. remaining risks or deferred behavior.

After a focused post-processing change is complete, stop and report. Do not
automatically move to Abel inversion, full experimental VMI interpretation, or
new physics branches.
