# Current migration state

## Completed

- `physics/constants.py`
- `config.py`
- `presets.py`
  - `single_pulse_N2000()`
  - `single_pulse_droplet_distribution()`
- `physics/potentials.py`
- `physics/interactions.py`
- `physics/leapfrog.py`
- `sampling/droplet_sizes.py`
- `sampling/radial_positions.py`
- `simulation/checkpoint.py`
- `simulation/run_directory.py`
- `sampling/orientations.py`
- `physics/collisions.py`
- `simulation/initial_state.py`
- `simulation/propagation_step.py`
- `simulation/neutral.py`
- Ion steps 11a, 11b, 11c, 11d, and 11e:
  - velocity-dependent cross-section helper
  - ion initial-state builder
  - pure ion propagation step
  - full ion propagation driver in `simulation/ion.py`
  - `E_mass_attach_defect_eV` diagnostic ported (IonCheckpoint schema
    v3 → v4); closes the per-side conservation invariant
    `E_kin + E_pot + E_dissip + E_mass_attach_defect ≈ const` when
    helium attachment is enabled.
- Ion-stage MATLAB/Python cross-reference validation:
  - target 1: ion `t=0` state from a tiny neutral checkpoint
  - targets 2-4: deterministic one-step, multi-step, and energy
    bookkeeping checks with collisions disabled
  - target 5: forced stochastic collision/mass-attachment behavior
- `scripts/run_single_pulse.py`:
  - public single-pulse neutral + ion pipeline entry point
  - supports selecting `single_pulse_N2000` or
    `single_pulse_droplet_distribution` as the input preset
  - writes `cfg.json`, `neutral.npz`, and `ion.npz` through `RunDirectory`
  - includes clobber protection unless `OVERWRITE_EXISTING_RUN = True`
- Step 13 post-processing path:
  - `i2_helium_md/postprocess/hedft_loader.py` loads normalized 8-column
    HeDFT trajectory CSVs (`9A_All_Data.csv`, `18A_All_Data.csv`)
  - `i2_helium_md/postprocess/compare_trajectories.py` computes numerical
    MD/HeDFT distance and velocity-magnitude comparisons on the overlap grid
  - `i2_helium_md/postprocess/velocity_distribution.py` loads VMI reference
    CSVs and computes mass-selected final-velocity histograms
  - `scripts/post_processing/plot_hedft_comparison.py` recreates the focused
    HeDFT comparison plotting workflow from an existing `RunDirectory`
- Legacy MATLAB live-debug and paper post-processing reproduction:
  - `i2_helium_md/postprocess/energy_balance.py` -- recipe helpers for
    sum-over-atoms neutral / per-molecule ion energy traces, simulated
    azimuthal phi histogram, and final ion mass spectrum.
  - `i2_helium_md/postprocess/_smoothing.py` -- shared MATLAB-style
    `movmean` and trace normaliser, reused by both
    `plot_experimental_comparison.py` and the new paper figure script.
  - `IonCheckpoint` schema bumped v4 -> v5: adds
    `temperature_diagnostic: (num_steps, 3)` capturing the per-step
    legacy MATLAB `[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>]`
    accumulator from `vmi_sim_3d_ion_propa.m:683`. Older v4 files must
    be regenerated.
  - `physics/collisions.py` exposes `CollisionDiagnostics` and
    `temperature_diagnostic_from_collision`; `apply_collision`
    accepts a `return_diagnostics=False` keyword for opt-in capture
    without changing existing call sites.
  - Focused post-processing scripts under `scripts/post_processing/`:
    `plot_neutral_energy_balance.py`, `plot_ion_energy_balance.py`,
    `plot_ion_temperature_diagnostic.py`, `plot_paper_v2.py`,
    `plot_paper_v3.py`, `plot_paper_v4.py`.
    - `plot_paper_v2.py` fully ports the active non-effusive
      `post_process_single_pulse_paper_IplusHe_comparison.m` branch for
      selected-run Python checkpoints, including processed experimental 2-D
      VMI reference loading, simulated physical `v_x/v_y` VMI map, radial
      comparison, and separate phi comparison. MATLAB exports are under
      `data/reference/scripts/export_paper_v2_reference_data.m`.
    - `plot_paper_v3.py` fully ports the active non-effusive
      `post_process_single_pulse_paper_v3.m` branch for selected-run Python
      checkpoints: radial and phi comparison plus separate ion mass histogram.
    - `plot_paper_v4.py` fully ports the active non-effusive
      `post_process_single_pulse_paper_v4.m` branch for selected-run Python
      checkpoints: radial comparison, simulated angular pair covariance, and
      separate ion mass histogram.
    - `plot_paper_cov.py` ports the active non-effusive
      `post_process_single_pulse_paper_IplusHe_comparison_cov.m` branch as
      six split figures: VMI comparison (exp + sim), 1-D velocity
      distribution comparison (exp I+ gas, exp I+He, exp v-cov trace from
      cov_radial diagonal, sim v-cov trace), 1-D phi(angle) distribution
      overlay (exp from precomputed MATLAB CSV
      `data/reference/paper_cov/iplus_he_phi.csv` -- literal MATLAB
      `mean(res_Iplus_He.image_polar(:, b_r), 2) / max(...)` for the same
      three measurement IDs as the covariance reference. The phi pipeline
      uses the HARDCODED `plot_processed_VMI` center `[524.5297, 380.8430]`
      from `_cov.m` line 100 (the covariance pipeline still uses the
      auto-detected center per `_cov.m` line 356); these two centers
      differ and reusing the auto-detected center for the phi pipeline
      rotates the resulting angular distribution out of phase with the
      live MATLAB figure. Sim phi from
      `atan2(vy, vx) + pi`), angular pair covariance (exp + sim
      side-by-side), radial pair-speed covariance (exp + sim side-by-side),
      and a 1-D pair-cov axis-sum trace figure (angular + radial traces,
      sim vs exp overlaid; `movmean(., 3)` -> `- min` -> `/ max`
      normalisation per the legacy recipe). The experimental covariance
      matrices and phi CSV live as frozen MATLAB-exported reference data
      under `data/reference/paper_cov/` (`iplus_he_covariance.mat` + JSON
      sidecar + `iplus_he_phi.csv`); the exporter is
      `data/reference/scripts/export_paper_cov_reference_data.m`. Python
      helpers in `i2_helium_md/postprocess/paper_cov.py`:
      `radial_pair_speed_covariance` (sim radial-speed pair covariance with
      MATLAB-style diag-zero and 2x2 movmean), `radial_covariance_trace`
      (1-D trace from cov_radial over the [4, 22] A/ps band),
      `load_paper_cov_experimental_reference`,
      `simulated_phi_distribution`, and
      `covariance_axis_sum_normalised`. The phi CSV reuses the generic
      `load_paper_v2_phi_reference` loader.
  - `scripts/post_processing/plot_run_summary.py` consolidates the in-scope
    diagnostics into one multi-page PDF plus per-panel PNGs.
  - Additional post-processing helpers cover polar velocity histograms,
    cos^2 anisotropy fits, beta(v), 2D velocity density, pair correlation,
    time-resolved radial distributions, Boltzmann overlays, bimodal Gaussian
    fits, and neutral-side HeDFT comparison.
  - `tests/test_energy_balance.py` covers the recipe helpers, the
    schema-v5 round trip, and the v4-reject path.
  - `tests/test_plot_legacy_debug_smoke.py` runs each new script in
    non-interactive mode against an existing run directory.
  - `tests/test_plot_run_summary_smoke.py` covers the consolidated summary
    driver with non-interactive matplotlib.
  - `tests/test_paper_cov.py` covers the `paper_cov` postprocess helpers
    (radial pair-speed covariance binning, mass + outside filters, 2x2
    movmean smoothing, the legacy `[4, 22] A/ps` v-cov trace recipe, and
    the `.mat` / `.npz` reference round-trip).
  - `tests/test_plot_paper_cov_smoke.py` covers the `plot_paper_cov.py`
    driver with and without the experimental covariance reference
    available (4 figures vs 2 figures fallback path).

## Current phase

The MATLAB → Python port is complete:

- neutral and ion propagation drivers,
- ion-stage MATLAB/Python cross-reference validation,
- public single-pulse run script,
- focused legacy-debug and paper post-processing scripts
  (`plot_neutral_energy_balance.py`, `plot_ion_energy_balance.py`,
  `plot_ion_temperature_diagnostic.py`, `plot_paper_v2.py`,
  `plot_paper_v3.py`, `plot_paper_v4.py`, `plot_paper_cov.py`),
- the consolidated `plot_run_summary.py` driver,
- experimental reference exports under `data/reference/paper_v2/`,
  `paper_v3/`, `paper_v4/`, `paper_cov/`, and `vmi_summary/`,
- visual comparison of the generated Python figures against the legacy
  MATLAB figures, with any recipe mismatches literal-ported.

The current phase is **drag-model physics**: deprecate the hard-sphere
collision model in `physics/collisions.py` and replace it with a
TDDFT-calibrated drag-force model for I⁺ in a helium bubble. Treat the
post-processing surface as the comparison layer for the new physics --
keep it stable. This is the first scope item that explicitly overrides
the `CLAUDE.md` "do not change collision physics" rule; the exception is
scoped to the drag-model port only.

## Currently pending

1. Record numerical MD/HeDFT comparison values for the 9 A HeDFT run
   (`data/runs/9A_hedft_comparison`) and decide which outputs should be kept
   as documented reference diagnostics.
2. Keep post-processing tests focused on loader contracts, overlap
   interpolation, VMI reference loading, final-velocity histogram filters, and
   plotting smoke coverage.
3. Keep Abel inversion, pump-probe, effusive dynamics, MATLAB multi-start
   matrix behavior, and full experimental VMI image interpretation out of
   scope unless explicitly requested. (Note: paper-cov experimental
   pair-covariance is in scope via the frozen MATLAB-exported reference
   under `data/reference/paper_cov/`.)

Items previously listed here that are now complete:

- `data/reference/scripts/export_paper_cov_reference_data.m` was re-run
  with the corrected `plot_processed_VMI` center `[524.5297, 380.8430]`
  for the phi pipeline; `data/reference/paper_cov/iplus_he_covariance.mat`,
  `iplus_he_covariance.json`, and `iplus_he_phi.csv` are regenerated and
  frozen (dated 2026-05-19). `plot_paper_cov.py` now produces all four
  comparison figures.
- The two `run_summary.pdf` outputs (`data/runs/9A_hedft_comparison/figures/`
  and `data/runs/single_pulse_droplet/figures/`) were visually reviewed
  against the legacy MATLAB figures.
- Recipe mismatches found during the review were literal-ported into the
  active Python helpers before any cosmetic cleanup.

## Recommended next task

Begin the drag-model port. As a first step, *survey* before editing:

1. enumerate every active call site of `apply_collision` in
   `physics/collisions.py`, `simulation/propagation_step.py`, and
   `simulation/ion.py`, including the `IonStepState` /
   `IonCheckpoint.temperature_diagnostic` data path that consumes the
   collision diagnostic;
2. identify which TDDFT reference dataset will calibrate the I⁺ /
   helium-bubble drag coefficient (velocity-dependent friction γ(v) or
   equivalent), and whether the existing 9 A / 18 A HeDFT trajectories
   under `data/reference/` are reusable or whether a new export is
   required;
3. design the new module as `physics/drag.py` parallel to
   `physics/collisions.py` rather than mutating the existing collision
   code; keep the collision module importable until the drag model is
   validated, so the two can be A/B compared against the same
   post-processing surface;
4. plan the deterministic / stochastic validation order matching
   `CLAUDE.md` "Scientific-Code Caution": formula check first, then
   one-step deterministic, then multi-step, then energy bookkeeping,
   then stochastic statistics, then full driver smoke test.

Do not broaden into Abel inversion or full experimental VMI interpretation
unless explicitly requested.

## Files to inspect for the current phase

Drag-model survey scope. Start with the directly relevant Python files:

- `i2_helium_md/physics/collisions.py` (existing hard-sphere model;
  source of the per-step `CollisionDiagnostics` and
  `temperature_diagnostic_from_collision` helpers — the contract the
  drag model must replace)
- `i2_helium_md/physics/interactions.py` (where ion / atom interaction
  potentials live; check whether the drag force should plug in here or
  in a new module)
- `i2_helium_md/simulation/propagation_step.py`
- `i2_helium_md/simulation/ion.py` (consumes the collision return value
  and writes `IonCheckpoint.temperature_diagnostic`)
- `i2_helium_md/config.py` (collision-related flags such as
  `collision_mode`, `mass_attach_enabled`)
- `i2_helium_md/physics/constants.py` (helium / I⁺ masses, accommodation
  / cross-section constants that may be replaced or retained)

Relevant tests:

- `tests/test_collisions.py`
- `tests/test_ion_propagation_step.py`
- `tests/test_ion.py`

Relevant MATLAB provenance:

- `legacy_matlab_repository/.../vmi_sim_3d_ion_propa.m` (current collision
  call sites, for reference only)
- whatever TDDFT-output files (forces, friction vs. velocity, or
  equivalent) end up being designated as the drag-calibration source —
  to be identified during the survey step.

The post-processing surface stays untouched in this phase; its files are
inspected only when validating the drag model output against existing
HeDFT / VMI references.

## Known MATLAB bugs not to reproduce

Future MATLAB/Python checks must not force Python to match known MATLAB
bookkeeping bugs.

Known intentional Python corrections include:

- neutral-stage `E_pot` at `t=0` includes the partner Morse contribution,
- ion-stage `E_kin` at `t=0` fixes the MATLAB velocity-expression bug,
- ion-stage `E_pot` at `t=0` fixes the MATLAB radial-coordinate bug,
- ion-stage `E_pot` at `t=0` includes the partner Coulomb contribution,
- Python uses more accurate physical constants than rounded MATLAB constants.

Tests should explicitly state whether Python is expected to match MATLAB or
intentionally differ.

## Do not do yet

- Do not refactor the neutral or ion drivers unless a test reveals a real bug.
- Do not implement out-of-scope MATLAB features.
- Do not expand the VMI plotting helpers into Abel inversion or full
  experimental analysis without an explicit request.

Drag-model exception: the user has explicitly approved replacing the
hard-sphere collision model with a TDDFT-calibrated drag model. This
is the one scoped exception to the `CLAUDE.md` forbidden-list item on
collision physics; it applies only to the drag-model work. The neutral
driver, the propagation cadence, the checkpoint schema, the RNG draw
order, and the physical-constants table remain off-limits.
