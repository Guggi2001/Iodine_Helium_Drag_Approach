# Current migration state

## Completed

- `physics/constants.py`
- `config.py`
- `presets.py`
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
  - writes `cfg.json`, `neutral.npz`, and `ion.npz` through `RunDirectory`
  - includes clobber protection unless `--force` is used

## Current phase

The neutral and ion propagation drivers are implemented and the ion-stage
MATLAB/Python cross-reference validation is complete. The public single-pulse
run script is implemented.

The next migration phase is Step 13: HeDFT loading and trajectory comparison in
`postprocess/`.

## Currently pending

1. Step 13: HeDFT loading and trajectory comparison in `postprocess/`

## Recommended next task

Design the first HeDFT loading and trajectory-comparison path in
`postprocess/`, starting with an inventory of the available reference data.
The legacy 9 Å files are present under
`legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/9Angström/`;
verify whether normalized copies already exist under `data/reference/` before
adding or using them.

Keep the first comparison narrow: load the reference data, load an existing
single-pulse run directory, and compute explicit numerical comparison values.
Avoid plotting-heavy workflows until the data-loading and numerical comparison
contracts are clear.

## Files to inspect for the current phase

Start with only the directly relevant Python files:

- `scripts/run_single_pulse.py`
- `i2_helium_md/presets.py`
- `i2_helium_md/simulation/run_directory.py`
- `i2_helium_md/simulation/checkpoint.py`
- `i2_helium_md/config.py`

Then inspect the relevant MATLAB orchestration files in
`legacy_matlab_repository/`, especially:

- `simulation_image_only_trajectories.m`
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

- Do not implement plotting before the HeDFT data-loading and numerical
  comparison path is clear.
- Do not refactor the neutral or ion drivers unless a test reveals a real bug.
- Do not implement out-of-scope MATLAB features.
