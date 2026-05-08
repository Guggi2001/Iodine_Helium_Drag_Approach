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
  - `scripts/plot_hedft_comparison.py` recreates the first HeDFT comparison
    plotting workflow from an existing `RunDirectory`
  - `scripts/post_processing_comparison/compare.py` remains as imported
    VMI-reference verification context
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
  - Four new post-processing scripts under `scripts/post_processing/`:
    `plot_neutral_energy_balance.py`, `plot_ion_energy_balance.py`,
    `plot_ion_temperature_diagnostic.py`, `plot_paper_figure.py`.
  - `tests/test_energy_balance.py` covers the recipe helpers, the
    schema-v5 round trip, and the v4-reject path.
  - `tests/test_plot_legacy_debug_smoke.py` runs each new script in
    non-interactive mode against an existing run directory.

## Current phase

The neutral and ion propagation drivers are implemented and the ion-stage
MATLAB/Python cross-reference validation is complete. The public single-pulse
run script is implemented. A first post-processing comparison path is also
implemented.

The current phase is post-processing validation and cleanup: use the new
comparison API on real run directories, document the numerical outputs, and
keep any further plotting or analysis changes narrowly scoped.

## Currently pending

1. Record numerical MD/HeDFT comparison values for the current production run
   (`data/runs/single_pulse_N_2000`) and decide which outputs should be kept
   as documented reference diagnostics.
2. Keep post-processing tests focused on loader contracts, overlap
   interpolation, VMI reference loading, final-velocity histogram filters, and
   plotting smoke coverage.
3. Regenerate existing ion checkpoints at schema v5 so the new
   temperature-diagnostic figures and the data-gated tests
   (`test_compare_trajectories::TestEndToEndReal`,
   `test_plot_*_smoke`) can run without manual deselection.
4. Polar-VMI panels of `post_process_single_pulse_paper_v3.m`
   (cos^2 angular anisotropy fit, beta(v) function, 3-D surf of
   polar VMI image) remain deferred. They require a 2-D polar VMI
   image not currently in `data/reference/`.

## Recommended next task

Run the implemented post-processing path on the existing production run and
write down explicit numerical diagnostics:

- `compare_distance(ion, hedft_9A)` RMSE and mean ratio,
- `compare_velocity_magnitude(..., atom="I1")` and `"I2"` RMSE/ratio,
- any mass-selected histogram caveats, especially if no atoms pass a target
  mass/outside filter.

Do not broaden into Abel inversion or full experimental VMI interpretation
unless explicitly requested.

## Files to inspect for the current phase

Start with the directly relevant Python files:

- `i2_helium_md/postprocess/hedft_loader.py`
- `i2_helium_md/postprocess/compare_trajectories.py`
- `i2_helium_md/postprocess/velocity_distribution.py`
- `scripts/plot_hedft_comparison.py`
- `scripts/post_processing_comparison/compare.py`

Relevant tests:

- `tests/test_hedft_loader.py`
- `tests/test_compare_trajectories.py`
- `tests/test_velocity_distribution.py`
- `tests/test_plot_hedft_comparison_smoke.py`

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
