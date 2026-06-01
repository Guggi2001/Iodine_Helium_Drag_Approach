# CLAUDE.md

This file bundles the project guidance that Claude-style agents need before
editing or reviewing this repository. It intentionally mirrors `agent.md` so
both agent entry points carry the same live rules.

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

The MATLAB-to-Python transfer is complete:

- neutral propagation is implemented,
- ion propagation is implemented,
- the main single-pulse presets are migrated,
- run directories write `cfg.json`, `neutral.npz`, and `ion.npz`,
- important MATLAB/Python propagation and bookkeeping paths have been
  cross-reference validated,
- the full in-scope post-processing surface is implemented: focused
  legacy-debug and paper scripts (`plot_neutral_energy_balance.py`,
  `plot_ion_energy_balance.py`, `plot_ion_temperature_diagnostic.py`,
  `plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`,
  `plot_paper_cov.py`), the consolidated `plot_run_summary.py` driver,
  and the supporting helpers under `i2_helium_md/postprocess/`,
- the experimental reference exports under `data/reference/paper_v2/`,
  `paper_v3/`, `paper_v4/`, `paper_cov/`, and `vmi_summary/` are
  populated and frozen.

The current phase is **drag-model physics**: deprecate the hard-sphere
collision model in `physics/collisions.py` and replace it with a
TDDFT-calibrated drag-force model for I⁺ in a helium bubble. The
post-processing surface is the comparison layer for the new physics and
must stay stable. This is the first scope item that explicitly overrides
the "do not change collision physics" forbidden-list rule below; the
exception is scoped to the drag-model port only.

**Planning is complete.** All architectural choices are locked in
`DRAG_PORT_DESIGN_DECISIONS.md` (noise, mass attachment, drag form,
integrator, spatial gating, validation). The frozen MD baseline is
`PHYSICS_BASELINE.md`; the upstream drag-law extraction pipeline is
`Drag_extraction_code.md`. The port is being implemented in dependency-ordered
**slices**. **Slice 1 (the pure gated-drag physics module) is complete and
reviewed; the active task is Slice 2** — the BAOAB ion-stage stepper —
specified in `SLICE2_GOALS_baoab_ion_stepper.md`. See the "Drag-Model Port"
section below for the working rules that apply to this phase.

## Current Scope

In scope:

- single-pulse neutral and ion dynamics,
- 9 A and normalized 18 A HeDFT comparison inputs already present in
  `data/reference/`,
- VMI reference loading and final-velocity histogram helpers,
- consolidated post-processing diagnostics from finished run directories,
- authentic reproduction of legacy post-processing figures where reference
  data and run outputs are available,
- focused MATLAB/Python reference validation,
- drag-model physics calibrated against TDDFT for I⁺ in the helium
  bubble (replacing the hard-sphere collision model).

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

## Drag-Model Port

This section governs the active drag-model phase. It sits above the
post-processing rules below, which now describe the stable **comparison
layer**, not active development.

### Working method

- **Physics-Definition vs. Software-Implementation boundary is strict.**
  Default to mathematical formulation (LaTeX) and logic trees. Do **not**
  write implementation code until the user gives the explicit trigger
  `[PROCEED TO IMPLEMENTATION]`. Discussion and specification come first.
- **Every model choice stays interchangeable behind a `SimConfig` enum** so
  the empirical cross-check against the TDDFT / VMI references can select
  among hypotheses. No drag form, noise model, mass scenario, or gate is
  hard-wired; each lives behind its own enum surface documented in
  `DRAG_PORT_DESIGN_DECISIONS.md`.
- **Before proposing a governing equation, perform and display a strict
  dimensional analysis.** Define the units of every parameter (forces in
  `amu·Å/ps²`, velocity in `Å/ps`, etc.) and reject any formulation whose
  units do not balance.
- When brainstorming a model transition (discrete collisions → continuous
  drag), list the physical trade-offs explicitly: lost degrees of freedom,
  violated conservation laws, thermodynamic changes (e.g. the loss of the
  Langevin fluctuation channel).

### Friction convention (unified — do not reintroduce the second one)

- **`γ(v)` is a force coefficient**, units **amu/ps**, defined
  `γ(v) = |F_drag(v)| / v`. The friction force is `γ(v)·v` — **no leading
  `m`**.
- The friction **rate** is `γ/m` [1/ps] and appears **only** inside the
  BAOAB damping exponent `e^(−γ·dt/m)`.
- The FDT noise amplitude `√(2·γ·k_B·T_eff)` uses `γ` [amu/ps] directly.
- Consequence: the drag-physics module is **mass-agnostic**. It never takes
  `m`. Mass enters only at the integrator's O-step, as one explicit division
  by `m(t)`. Do not pass mass into the drag module.

### Slice plan (dependency-ordered)

The Tier-0 critical path is the only path runnable with coefficients in hand
(`linear_cubic` + `mass_scenario=fixed` at `m_eff`, noise off). It is built in
four slices:

1. **Slice 1 — pure gated-drag physics module.** *Complete.* No
   dependencies. Delivered: `physics/drag.py` (three pure mass-free functions +
   `DragCoefficients` bundle type). Coefficients frozen under
   `data/reference/drag/<case>/linear_and_cubic/`. See
   `SLICE1_GOALS_gated_drag_module.md` and `drag_module.md`.
2. **Slice 2 — BAOAB ion-stage stepper.** *Active task.* Consumes Slice 1's
   `γ(v)`; new `physics/baoab.py`; noise amplitude pinned to zero at Tier 0.
3. **Slice 3 — `SimConfig` enum surface + §6.5 mass↔coefficient consistency
   guard.**
4. **Slice 4 — ion-driver rewiring + O-step energy accounting.**

Mass dynamics (§2), the `IonCheckpoint` v6 rename
(`E_mass_attach_defect_eV` → `E_mass_transfer_eV`), and the noise machinery
(§1) stay **stubbed behind their enums and inert** until their validation
tier comes up. None are on the Tier-0 path.

### Slice 1 — complete

Delivered `physics/drag.py`: three pure, config-free, mass-free functions —
`drag_force(v, depth, …)` [amu·Å/ps²], `drag_gamma(v, depth, …)` [amu/ps,
closed form `g·(a+b·v²)`], `spatial_gate(depth, steepness)` [dimensionless
erf-complement] — plus the `DragCoefficients` bundle type (form tag,
coefficients, `extraction_mass_model` + `extraction_mass_amu` provenance).
`linear_cubic` realised; `linear_quadratic` / `threshold` / `power_law`
reserved and raise `NotImplementedError`. `physics/collisions.py` left intact
and importable (additive, parallel — not a deletion). Specs:
`SLICE1_GOALS_gated_drag_module.md`, `drag_module.md`.

Two upstream verifications were flagged and remain the user's to confirm at the
extraction source (a self-consistent refit cannot detect either):

- the stamped `extraction_mass_amu ≈ 202.954` is the mass the force balance
  *actually ran under*, not a relabelled value (an earlier literal was
  `179.912`);
- `drag_data.csv`'s `F_drag` column stores the **positive magnitude** (signed
  force balance is `F_drag = m_eff·a − F_C`, negative during the explosion).

Post-extraction empirical finding folded into the decisions doc: the
`power_law` exponent is `n ≈ +2` (not the anticipated `−2`), so the form is
regular at `v→0` and `drag_low_v_floor` is inert for the real coefficients
(retained only for a hypothetical `n<0` re-extraction).

### Slice 2 — active task

Specified in full in `SLICE2_GOALS_baoab_ion_stepper.md`. New module
`physics/baoab.py`: the BAOAB operator-split ion-stage stepper (decision §4.6),
replacing `velocity_verlet_step` for the ion stage only. Summary of intent:

- **Scheme B–A–O–A–B:** B/A are the baseline kick/drift (conservative force =
  Coulomb + droplet via the existing `_ion_accel_fn`); O is the new physics —
  drag as multiplicative velocity damping `v ↦ e^(−γ·dt/m)·v`, plus a dormant
  Langevin-noise site.
- **Tier-0 reduction:** mass fixed at `m_eff`, `T_eff = 0` (noise off),
  `linear_cubic` drag. The stepper is built and tested fully deterministically
  and is still the production Tier-0 integrator.
- **Asymmetric γ-freeze (intentional, documented in the module docstring):**
  γ's velocity argument frozen at the O-step input velocity; γ's depth/gate
  argument at the current O-step position (freshly known after the first
  half-drift). This keeps the never-adds-energy dissipativity exact.
- **Placement:** new `physics/baoab.py`, with `_kick`/`_drift` helpers
  **extracted from `leapfrog.py`** to avoid duplicate physics. This is the one
  knowing touch to the frozen baseline integrator — accepted because §4.6 added
  an integrator at all, and the extraction is self-verified by the anchor test.
- **Energy:** the O-step returns dissipated energy in **amu·Å²/ps²** (Slice 4
  converts to eV); returned now though consumed only by Slice 4's §2.9 invariant.
- **Per-step closure rebuild:** `make_ion_baoab_step` mirrors
  `make_ion_step` and is rebuilt per step (mass changes under future
  scenarios), matching `ion_propagation_step.py:184-193`.
- **Anchor (killer) test:** at `γ=0` with noise off, BAOAB ≡ baseline
  `velocity_verlet_step` to round-off — proving both integrator correctness and
  that the kick/drift extraction left the frozen baseline behaviourally
  unchanged.

Slice 2 scope fence — does **not** touch: `SimConfig` fields; driver wiring
(`ion_propagation_step.py:211-256`); checkpoint schema / energy rename; active
noise; mass dynamics; the eV conversion; the neutral stage (the only
`leapfrog.py` change is the behaviour-preserving kick/drift extraction).

### Validation hierarchy (sequential, not simultaneous)

The undetermined parameters are entangled; validate in tier order, fixing
each tier's winner before introducing the next unknown:

- **Tier 0** — drag form, deterministic, fixed mass, in-window TDDFT traces.
- **Tier 1** — mass scenario (A/B/biphasic), deterministic.
- **Tier 2** — terminal I⁺(He)ₙ size distribution vs. experimental detector
  data (the only observable that separates the mass scenarios).
- **Tier 3** — ensemble second moments (noise) vs. VMI references.

A `mass_scenario`↔`drag_coefficients` consistency guard is enforced at
config-load (§6.5): constant-mass coefficients are self-consistent **only**
with `mass_scenario=fixed`; the inconsistent pair is a hard refuse unless
`allow_inconsistent_mass_pairing=True`. Histogram comparisons default to the
Wasserstein metric.

## Post-Processing Comparison Layer

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

## General Request Handling

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

During the drag-model phase, "implementation" requests still require the
`[PROCEED TO IMPLEMENTATION]` trigger before any code is written (see
Drag-Model Port → Working method).

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
scripts/post_processing/plot_paper_cov.py
scripts/post_processing/plot_run_summary.py
```

`plot_run_summary.py` is configured through a `# USER SETTINGS` block at the
top of the script (the same pattern as the focused scripts). Edit the
constants and run the script — either from PyCharm's run button or from a
shell:

```bash
python scripts/post_processing/plot_run_summary.py
```

For the 9 A HeDFT summary, set:

```python
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_hedft_comparison"
HEDFT_REF_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
VMI_REF_HE_PATH = None
VMI_REF_GAS_PATH = None
```

For the experimental-condition droplet summary, set:

```python
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
HEDFT_REF_PATH = None
VMI_REF_HE_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he.csv"
VMI_REF_GAS_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"
```

The post-processing port is complete. Further changes are bug-fix-only:
if a panel's numerical or visual behavior disagrees with the legacy
MATLAB figure, follow the legacy-first rule:

1. inspect the exact MATLAB recipe,
2. literal-port the normalization, binning, smoothing, filtering, or fit
   behavior,
3. verify numerically or visually,
4. only then refactor the Python helper for clarity.

Do not read the full legacy plotting stack unless a specific numerical or
visual discrepancy requires it. Do not add new post-processing analysis
scope (Abel inversion, pump-probe, effusive comparison, full experimental
VMI image interpretation remain out of scope unless explicitly requested).

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
`pcolormesh(X, Y, C)`: `vx_mps` should be the plot x-grid, `vy_mps` the plot
y-grid, and `intensity` the color array. If the MATLAB source plots full
coordinate matrices, export full 2-D coordinate grids rather than slicing
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
edges_velocity = 0:0.04:26  (A/ps internally; 0:4:2600 m/s on display)
vd_ion = movmean(h, 15)
xlim([0, 2800])  m/s on display (= 28 A/ps in the legacy MATLAB figure)
```

Preserve the fine `0.04 A/ps` (4 m/s) bins, 15-bin moving mean, and displayed
range to `2800 m/s` when matching the legacy figure. The Python plotting
scripts (`plot_experimental_comparison.py`, `plot_paper_v2.py`) display m/s on
the axis; simulation histograms still bin in A/ps internally and expose a
`bin_centers_mps` field for plotting (multiply by 100).

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

For the drag-model port the deterministic configuration is the Tier-0 setup:
`mass_scenario=fixed` at `m_eff`, noise amplitude zero, `linear_cubic` drag.
The seven-step validation order above maps onto the drag-port tier hierarchy
(Drag-Model Port → Validation hierarchy): formula/shape/one-step/multi-step
deterministic checks are Tier 0–1; energy bookkeeping closes the §2.9
invariant; statistical comparison is Tier 2–3.

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

Active exception (scoped, time-limited): the user has explicitly
approved replacing the hard-sphere collision model with a
TDDFT-calibrated drag model for I⁺ in the helium bubble. This relaxes
the "changing collision physics" rule only for that work. The neutral
driver, checkpoint schema, RNG draw order, default simulation scope,
and physical-constants table remain off-limits.

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
