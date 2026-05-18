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

The neutral and ion propagation drivers are implemented and the ion-stage
MATLAB/Python cross-reference validation is complete. The public single-pulse
run script is implemented. The in-scope post-processing port now includes both
focused plotting scripts and the consolidated `plot_run_summary.py` driver.

The current phase is authentic post-processing cleanup and review: the focused
v2/v3/v4 paper scripts are ported for the active droplet branches, but visual
comparison against regenerated MATLAB references should continue to catch
remaining convention mismatches. Keep changes narrowly focused on faithful
reproduction rather than new analysis scope.

## Currently pending

1. Re-run `data/reference/scripts/export_paper_cov_reference_data.m`
   once in MATLAB (with the legacy VMI toolbox on the path) to
   regenerate `data/reference/paper_cov/iplus_he_covariance.mat` +
   JSON sidecar AND `data/reference/paper_cov/iplus_he_phi.csv`. The
   exporter was recently corrected to use the hardcoded
   `plot_processed_VMI` center `[524.5297, 380.8430]` for the phi
   pipeline (matching `_cov.m` line 100); the previous version reused
   the auto-detected covariance center and produced a phi curve out of
   phase with the live MATLAB figure. The CSV must be regenerated so
   `plot_paper_cov.py` produces the two side-by-side covariance
   figures, the 1-D pair-cov trace figure, and the experimental
   overlay on the phi distribution. Without the covariance reference,
   the three cov-derived figures are skipped (with a warning). Without
   the phi CSV the phi figure draws the simulated curve alone (with a
   warning).
2. Review `data/runs/9A_hedft_comparison/figures/run_summary.pdf` and
   `data/runs/single_pulse_droplet/figures/run_summary.pdf` against the
   corresponding legacy MATLAB post-processing figures.
3. When mismatches are found in already-ported scripts, correct the relevant
   MATLAB post-processing recipe detail before refactoring the Python version
   for clarity.
4. Record numerical MD/HeDFT comparison values for the 9 A HeDFT run
   (`data/runs/9A_hedft_comparison`) and decide which outputs should be kept
   as documented reference diagnostics.
5. Keep post-processing tests focused on loader contracts, overlap
   interpolation, VMI reference loading, final-velocity histogram filters, and
   plotting smoke coverage.
6. Keep Abel inversion, pump-probe, effusive dynamics, MATLAB multi-start
   matrix behavior, and full experimental VMI image interpretation out of
   scope unless explicitly requested. (Note: paper-cov experimental
   pair-covariance is now in scope via the frozen MATLAB-exported reference
   under `data/reference/paper_cov/`.)

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
