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

## Current phase

The neutral and ion propagation drivers are implemented and the ion-stage
MATLAB/Python cross-reference validation is complete. The public single-pulse
run script is implemented. The in-scope post-processing port now includes both
focused plotting scripts and the consolidated `plot_run_summary.py` driver.

The current phase is authentic post-processing porting and cleanup: compare
the generated Python PDFs against the legacy MATLAB figures, tighten visual
and numerical conventions where behavior still differs, and keep changes
narrowly focused on faithful reproduction rather than new analysis scope.

## Currently pending

1. Review `data/runs/9A_hedft_comparison/figures/run_summary.pdf` and
   `data/runs/single_pulse_droplet/figures/run_summary.pdf` against the
   corresponding legacy MATLAB post-processing figures.
2. When mismatches are found, port the relevant MATLAB post-processing recipe
   more authentically before refactoring the Python version for clarity.
3. Record numerical MD/HeDFT comparison values for the 9 A HeDFT run
   (`data/runs/9A_hedft_comparison`) and decide which outputs should be kept
   as documented reference diagnostics.
4. Keep post-processing tests focused on loader contracts, overlap
   interpolation, VMI reference loading, final-velocity histogram filters, and
   plotting smoke coverage.
5. Keep Abel inversion, pump-probe, effusive dynamics, and full experimental
   VMI image interpretation out of scope unless explicitly requested.

## Recommended next task

Perform an authentic post-processing pass on the generated summary PDFs:

- inspect the legacy MATLAB recipe for any panel whose Python output differs,
- port the MATLAB normalization, binning, smoothing, filtering, or fit recipe
  literally enough to reproduce the intended behavior,
- then document the resulting numerical diagnostics, especially
  `compare_distance(ion, hedft_9A)` RMSE/ratio and
  `compare_velocity_magnitude(..., atom="I1"|"I2")` RMSE/ratio.

Do not broaden into Abel inversion or full experimental VMI interpretation
unless explicitly requested.

## Files to inspect for the current phase

Start with the directly relevant Python files:

- `i2_helium_md/postprocess/hedft_loader.py`
- `i2_helium_md/postprocess/compare_trajectories.py`
- `i2_helium_md/postprocess/velocity_distribution.py`
- `i2_helium_md/postprocess/energy_balance.py`
- `i2_helium_md/postprocess/polar_velocity.py`
- `i2_helium_md/postprocess/velocity_2d.py`
- `i2_helium_md/postprocess/pair_correlation.py`
- `i2_helium_md/postprocess/time_resolved.py`
- `i2_helium_md/postprocess/boltzmann_overlay.py`
- `scripts/post_processing/plot_run_summary.py`
- `scripts/post_processing/plot_hedft_comparison.py`
- `scripts/post_processing/plot_experimental_comparison.py`

Relevant tests:

- `tests/test_hedft_loader.py`
- `tests/test_compare_trajectories.py`
- `tests/test_velocity_distribution.py`
- `tests/test_plot_hedft_comparison_smoke.py`
- `tests/test_plot_run_summary_smoke.py`
- `tests/test_polar_velocity.py`
- `tests/test_velocity_2d.py`
- `tests/test_pair_correlation.py`
- `tests/test_time_resolved.py`
- `tests/test_boltzmann_overlay.py`

Relevant MATLAB provenance files:

- `simulation_image_only_trajectories.m`
- `simulation_image.m`
- files under `single_pulse_simulation/HeDFT_comparison/9Angström/`

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
