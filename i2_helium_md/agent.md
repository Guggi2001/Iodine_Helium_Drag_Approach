# Agent Guide

This file bundles the project guidance that agents need before editing or
reviewing this repository. It summarizes the live rules from `CLAUDE.md`,
`migration_log.md`, `current_state.md`, `testing.md`,
`agent_protocol.md`, and `docs/post_processing_port_plan.md`.

## Project Snapshot

This repository is a Python port of a legacy MATLAB molecular-dynamics
codebase for iodine / iodine-ion dynamics in helium nanodroplets.

Main package:

```text
i2_helium_md
```

Legacy MATLAB reference:

```text
legacy_matlab_repository/
```

The core MATLAB-to-Python simulation transfer is mainly complete:

- neutral propagation is implemented,
- ion propagation is implemented,
- the main single-pulse presets are migrated,
- run directories write `cfg.json`, `neutral.npz`, and `ion.npz`,
- important MATLAB/Python propagation and bookkeeping paths have been
  cross-reference validated,
- in-scope post-processing helpers and scripts are implemented,
- the full test suite was verified at 422 passing tests.

The current phase is authentic post-processing porting and cleanup. Compare
the generated Python PDFs against the legacy MATLAB figures, then tighten
binning, smoothing, normalization, filtering, fitting, and plotting
conventions where the behavior still differs. Keep the goal faithful
reproduction of in-scope legacy diagnostics, not new analysis scope.

## Current Scope

In scope:

- single-pulse neutral and ion dynamics,
- 9 A and normalized 18 A HeDFT comparison inputs already present in
  `data/reference/`,
- VMI reference loading and final-velocity histogram helpers,
- consolidated post-processing diagnostics from finished run directories,
- authentic reproduction of legacy post-processing figures where reference
  data and run outputs are available,
- focused MATLAB/Python reference validation.

Out of scope unless the user explicitly asks:

- pump-probe support,
- effusive / gas-phase dynamics,
- Abel inversion,
- full experimental VMI image interpretation,
- broad experimental VMI analysis beyond current overlays,
- new physics branches,
- broad refactors,
- live-debug 3D animations and visualization-only MATLAB utilities.

## Project Quality Principles

Use these principles for every porting decision, code review, and cleanup:

1. No duplicate implementations of the same physics. Shared conversions,
   formulas, and constants belong in `constants.py` or the appropriate shared
   module.
2. No dead code. Remove unused imports, commented-out blocks, and speculative
   branches.
3. Encode units and conventions in names: `mass_kg`, `time_ps`,
   `T_particles_K`, `R0_GS_angstrom`, etc.
4. Validate early and fail loudly. Wrong shape, unsupported collision mode,
   invalid type, or non-overlapping time axes should raise clear errors.
5. Public functions need docstrings with units, shapes, edge cases, inputs,
   and outputs.
6. Organize modules by concern: `physics/` is science, `sampling/` is
   randomness, `simulation/` is orchestration, `postprocess/` is analysis of
   finished runs.
7. Tests document intended behavior, units, tolerances, and known
   MATLAB/Python deviations.
8. Audit after refactors for dead imports, duplicate physics, and drift
   between docs and code.
9. For any known reference output, literal transliteration comes before clean
   refactor. First reproduce the MATLAB behavior, then refactor once the
   numerical result is verified.
10. Do not preserve bad legacy behavior merely for byte identity. Python uses
    corrected modern constants and fixes known MATLAB bookkeeping bugs unless
    the user explicitly asks for legacy behavior.

## Architecture Rules

Use `SimConfig` instead of globals. MATLAB global settings were consolidated
into `i2_helium_md/config.py`; functions that need parameters should receive
`cfg: SimConfig` explicitly.

Use preset functions instead of scripts:

```python
from i2_helium_md import (
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
    single_pulse_droplet_distribution,
)
```

Use `RunDirectory` for simulation artifacts. A run directory is
self-describing and should contain:

```text
cfg.json
neutral.npz
ion.npz
figures/      optional post-processing outputs
```

Checkpoint I/O rules:

- checkpoints are explicit dataclasses,
- every checkpoint has `schema_version`,
- incompatible versions fail at load time,
- no constants are saved in checkpoints,
- config belongs in `cfg.json`,
- load with `allow_pickle=False`,
- shape validation should use `cfg` when available.

Avoid changing checkpoint schema, physical constants, neutral propagation,
ion propagation, collision physics, or `scripts/run_single_pulse.py` behavior
unless a focused test reveals a real bug or the user explicitly asks.

## Current Working Mode

Prefer current Python APIs over ad hoc scripts:

1. Load a finished run with `RunDirectory`.
2. Load HeDFT references with `load_hedft_trajectory`.
3. Compute numerical diagnostics with `compare_distance` and
   `compare_velocity_magnitude` before changing trajectory plots.
4. Use `velocity_distribution.py` for VMI reference loading and mass-selected
   final-velocity histograms.
5. Use focused post-processing helpers instead of rolling new histograms:
   `energy_balance.py`, `polar_velocity.py`, `velocity_2d.py`,
   `pair_correlation.py`, `time_resolved.py`, `boltzmann_overlay.py`.
6. Use `scripts/post_processing/plot_run_summary.py` for every in-scope
   diagnostic from a finished run.
7. Keep plot changes local to `scripts/post_processing/` unless a package API
   change is actually needed.
8. Add or update focused pytest coverage when behavior changes.

For investigation, audit, inspect, compare, or explain requests:

- do not edit files,
- read relevant docs and tests,
- report files inspected,
- report conclusions and uncertainties,
- suggest the smallest safe next edit.

For implementation or fix requests:

- run `git status` first,
- preserve unrelated user changes,
- make the smallest coherent change,
- do not touch unrelated files,
- add or update tests when behavior changes,
- run relevant tests,
- show changed files and remaining risks.

## Post-Processing Workflow

The legacy combined MATLAB figure mixed outputs from different run settings.
Keep these workflows separate:

- HeDFT trajectory comparison uses a HeDFT-comparison run.
- Experimental VMI distribution comparison uses the realistic
  experimental-condition run.
- Do not plot the HeDFT velocity-vs-time panel from an experimental-condition
  run.
- Do not plot experimental VMI distributions from a HeDFT-comparison run.

Current script entry points:

```text
scripts/post_processing/plot_hedft_comparison.py
scripts/post_processing/plot_experimental_comparison.py
scripts/post_processing/plot_neutral_energy_balance.py
scripts/post_processing/plot_ion_energy_balance.py
scripts/post_processing/plot_ion_temperature_diagnostic.py
scripts/post_processing/plot_paper_v2.py
scripts/post_processing/plot_paper_v3.py
scripts/post_processing/plot_paper_v4.py
scripts/post_processing/plot_run_summary.py
```

Build the 9 A HeDFT summary:

```bash
python scripts/post_processing/plot_run_summary.py data/runs/9A_hedft_comparison --hedft-ref data/reference/9A_All_Data.csv --no-show
```

Build the experimental-condition droplet summary:

```bash
python scripts/post_processing/plot_run_summary.py data/runs/single_pulse_droplet --vmi-ref-he data/reference/vmi_iplus_he.csv --vmi-ref-gas data/reference/vmi_iplus_gas.csv --no-show
```

The next post-processing work is authentic comparison against legacy MATLAB
figures. For any mismatch:

1. inspect the exact MATLAB recipe,
2. literal-port the normalization, binning, smoothing, filtering, or fit
   behavior,
3. verify numerically or visually,
4. only then refactor the Python helper for clarity.

Do not read the full legacy plotting stack unless a specific numerical or
visual discrepancy requires it.

## Data Contracts

Normalized reference data lives in `data/reference/`.

Expected files:

```text
9A_All_Data.csv
18A_All_Data.csv
vmi_iplus_he.csv
vmi_iplus_gas.csv
```

HeDFT trajectory CSVs use this header:

```text
Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance
```

The earlier split 9 A files (`data_vabs2.csv`, `R1-R2.csv`) are provenance,
not the current loader contract.

Do not use hardcoded absolute MATLAB paths. Before adding or repurposing data,
verify whether existing `data/reference/*.csv` files are experimental VMI data
or HeDFT trajectory data.

When a legacy MATLAB post-processing script loads or processes experimental
data, treat that MATLAB path as the legacy source of truth for the data
processing recipe. First run or adapt the MATLAB processing path to export the
processed experimental result into a small, inspectable reference format,
preferably CSV. Save the exported file under `data/reference/`, and add or keep
the MATLAB export script under `data/reference/scripts/`, following the
existing `data/reference/scripts/export_vmi_reference_data.m` precedent.

For processed 2-D VMI image references, use matrix data plus a JSON sidecar
when CSV would make the artifact large or awkward. MATLAB exporters should
prefer `.mat` files so they do not depend on MATLAB's Python bridge; Python may
also accept `.npz` with the same fields for manually converted references. The
matrix file should store calibrated axis arrays and intensity separately, for
example `vx_mps`, `vy_mps`, and `intensity` (m/s on disk; Python loaders convert
to A/ps internally). Normalize exported image grids for Matplotlib
`pcolormesh(X, Y, C)`: `vx_mps` should be the plot x-grid, `vy_mps`
the plot y-grid, and `intensity` the color array. If the MATLAB source plots
full coordinate matrices, export full 2-D coordinate grids rather than slicing
constant-looking row or column vectors. The sidecar should document the MATLAB
source, measurement or MAT-file source, center, velocity factor, axis
equations, units, and external toolbox requirement. Keep 1-D radial or angular
curves as CSV whenever practical.

For these experimental-data exports, document the provenance: original MATLAB
script or function, measurement IDs or input files, processing steps,
calibration and scaling factors, output columns and units, and any external
toolbox requirement. Python should load the exported reference data and
reproduce or compare against it rather than reimplementing opaque extraction
from raw lab/toolbox inputs. Full experimental VMI interpretation and
Abel/image-processing expansion remain out of scope unless explicitly
requested.

Reference data should be small, inspectable, and reproducible. Prefer small
CSV, JSON, NPZ, or text files. Avoid committing large checkpoints, large MAT
files, generated figures, temporary debugging outputs, or full simulation
output directories.

If MATLAB reference data is generated, document the MATLAB script/command,
input parameters, random seed, molecule count, timestep count, enabled or
disabled physics features, and whether the data contains a known MATLAB bug or
an intentional Python correction.

## Known Plotting Conventions

Velocity-vs-time HeDFT panel:

- MATLAB sampled roughly 15 molecules and plotted both iodine atoms.
- Python should cap the overlay near 30 velocity traces.

Experimental velocity distribution:

```text
edges_velocity = 0:0.04:26
vd_ion = movmean(h, 15)
xlim([0, 28])
```

Preserve the fine `0.04 A/ps` bins, 15-bin moving mean, and displayed range to
`28 A/ps` when matching the legacy figure.

Polar histogram and anisotropy:

- fit model: `f(phi) = a + b * cos(phi - phi0)^2`,
- beta recovery: `beta = 2*b / (2*a + b)`,
- beta range: `[-1, 2]`,
- `beta(v)` skips bins with fewer than 50 counts by default.

Angular pair covariance:

- `theta = arctan2(vx, vy) + pi`,
- zero the covariance diagonal by default to mirror
  `cov_angular - diag(...)`.

Boltzmann overlay:

- use `physics.droplet_potential`,
- use `cfg.potential_steepness_molecule`,
- convert `cfg.binding_energy_molecule_K` to eV,
- normalize by trapezoidal integration on the chosen radial grid.

## Scientific-Code Caution

Clean code is not automatically correct physics.

Before changing a formula, unit conversion, force sign, random sampler,
normalization convention, indexing convention, or draw order:

1. locate the corresponding MATLAB source or previous Python test,
2. explain the convention,
3. add a regression test or focused numerical check,
4. then edit.

Do not start validation with a full stochastic trajectory comparison. Too many
effects are entangled. Validate in this order:

1. direct formula comparison,
2. shape and unit checks,
3. one-step deterministic comparison,
4. multi-step deterministic comparison,
5. energy bookkeeping comparison,
6. stochastic statistical comparison,
7. full driver smoke test.

For deterministic tests, disable stochastic features when possible:

- collisions disabled,
- mass attachment disabled,
- fixed initial state,
- fixed molecule count,
- few timesteps,
- fixed droplet radius,
- fixed random seed.

For stochastic tests, prefer distribution and moment checks over exact
trajectory matching unless RNG identity is guaranteed.

## Known MATLAB Bugs Not To Reproduce

Do not force Python to match these known MATLAB bookkeeping bugs:

- neutral-stage `E_pot` at `t=0` omitted partner Morse contribution,
- ion-stage `E_kin` at `t=0` used an incorrect velocity expression,
- ion-stage `E_kin` at `t=0` omitted `vz`,
- ion-stage `E_pot` at `t=0` omitted the `z` coordinate,
- ion-stage `E_pot` at `t=0` omitted the partner Coulomb term.

Tests should state whether Python is expected to match MATLAB, match within
known constant/unit differences, statistically match, or intentionally differ
because a MATLAB bug was corrected.

## Testing

Default test command:

```bash
pytest
```

In this environment, Python may not be on PATH. The absolute interpreter that
has been used successfully is:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

After editing code:

1. run the narrowest relevant test first,
2. run broader tests if multiple modules are affected,
3. report exactly which tests were run,
4. report failures honestly,
5. do not claim correctness without tests or a focused numerical check.

For post-processing changes, minimum coverage should include whichever apply:

- loader tests for normalized HeDFT CSV format,
- comparison tests using tiny synthetic checkpoints,
- validation for missing files and non-overlapping time axes,
- VMI reference loader tests,
- mass-filtered final-velocity histogram tests,
- plotting smoke coverage with non-interactive matplotlib.

Do not generate figures or production-sized checkpoints in tests.

Numerical tolerances must be justified:

- analytical formula ports should use tight tolerances,
- finite-difference forces need tolerances consistent with the FD step,
- Monte Carlo samplers need statistical tolerances based on sample size,
- MATLAB/Python comparisons may need tolerances for known constant updates,
- loose tolerances need an explanatory test comment.

Do not update expected values blindly. First determine whether the code is
wrong, the test encoded a legacy bug, the tolerance is unreasonable, the model
changed, a constant difference is expected, the stochastic test is
under-sampled, or the reference data came from the wrong MATLAB path.

## Forbidden Without Explicit User Approval

- deleting reference data,
- changing physical constants,
- changing checkpoint schema,
- changing random-number draw order,
- changing default simulation scope,
- changing neutral or ion propagation physics without a focused bug,
- changing collision physics without a focused bug,
- broad refactors,
- optimizing performance by changing numerical behavior,
- implementing out-of-scope MATLAB paths.

## Reporting

For post-processing work, report:

1. MATLAB files inspected,
2. Python files changed,
3. data files or run directories used,
4. numerical or plotting behavior changed,
5. tests run and results,
6. remaining risks or deferred behavior.

After a focused post-processing change, stop and report. Do not automatically
continue into Abel inversion, full experimental VMI interpretation, new
physics branches, or broad cleanup.
