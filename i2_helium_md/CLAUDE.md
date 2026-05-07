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

The current task is **Step 13: HeDFT loading and trajectory comparison in
`postprocess/`**.

Some post-processing comparison code has already been imported in:

```text
scripts/post_processing_comparison/compare.py
```

That script currently verifies the exported VMI reference CSVs in
`data/reference/` (`vmi_iplus_he.csv` and `vmi_iplus_gas.csv`) by recreating a
comparison plot. Treat it as existing context for post-processing work, but do
not confuse it with the pending HeDFT trajectory loader/comparison API for
Step 13.

Do not continue with plotting-heavy workflows, Abel inversion, experimental VMI
comparison, pump-probe support, 18 Å HeDFT support, analytical-force cleanup,
broad refactors, or out-of-scope MATLAB branches unless the user explicitly
asks.

## Main working mode

Work data-contract first.

Prefer this workflow:

1. Inventory the available reference data under `data/reference/` and the
   legacy 9 Å HeDFT files under
   `legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/`.
2. Inspect the MATLAB comparison script only for the numerical comparison
   contract:
   - `legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image_only_trajectories.m`
   - `legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/importfile_v2.m`
   - `legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/importfile_R1_R2.m`
3. Implement the smallest loader for the 9 Å trajectory reference files.
4. Implement a numerical comparison against an existing `RunDirectory` ion
   checkpoint.
5. Add focused pytest coverage with tiny synthetic CSV/checkpoint data.
6. Only after the numerical comparison is clear, decide whether plotting or a
   public command-line wrapper is useful.

Keep this step narrow. The goal is reliable loading and explicit numerical
comparison values, not final publication plots.

## Expected Step 13 outputs

Preferred new files:

```text
i2_helium_md/postprocess/hedft_loader.py
i2_helium_md/postprocess/compare_trajectories.py
tests/test_hedft_loader.py
tests/test_compare_trajectories.py
docs/hedft_loader_module.md
docs/compare_trajectories_module.md
```

Only create `i2_helium_md/postprocess/__init__.py` if the package directory
does not already exist.

The first comparison should compute, at minimum:

- mean MD I-I distance vs. time from the ion checkpoint,
- interpolation of MD distance onto the HeDFT `R1-R2.csv` time grid,
- overlap interval,
- RMSE in Å,
- mean MD/HeDFT distance ratio.

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

Expected 9 Å reference inputs are:

```text
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/data_vabs2.csv
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/R1-R2.csv
```

If corresponding normalized copies are missing under `data/reference/`, add
small documented copies there using stable ASCII filenames, for example:

```text
data/reference/hedft_9A_velocity.csv
data/reference/hedft_9A_distance.csv
```

Before adding data files, verify whether existing `data/reference/*.csv` files
are experimental VMI data or HeDFT trajectory data. Do not silently repurpose
files with ambiguous names.

## Relevant files

Start with these Python files:

```text
scripts/run_single_pulse.py
i2_helium_md/simulation/run_directory.py
i2_helium_md/simulation/checkpoint.py
i2_helium_md/config.py
```

Then inspect or create:

```text
i2_helium_md/postprocess/
```

Also inspect the existing imported post-processing comparison context if the
work touches VMI reference-data verification:

```text
scripts/post_processing_comparison/compare.py
```

Relevant MATLAB files:

```text
legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/simulation_image_only_trajectories.m
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

During Step 13, do not change unless a focused comparison test reveals a real
bug:

```text
checkpoint schema
physical constants
neutral propagation physics
ion propagation physics
collision physics
single-pulse run script behavior
```

Keep code changes local to `postprocess/`, focused tests, small reference CSV
copies if needed, and status documentation.

## Testing expectations

At minimum:

- loader tests for the 9 Å velocity and distance CSV formats,
- comparison tests using a tiny synthetic `IonCheckpoint`,
- validation for missing files and non-overlapping time axes,
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

After the first numerical comparison path is complete, stop and report. Do not
automatically move to plotting or experimental VMI comparison.
