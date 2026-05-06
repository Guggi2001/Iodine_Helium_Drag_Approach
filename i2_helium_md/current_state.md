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

## Current phase

The ion driver is implemented. The project is now in the MATLAB/Python
cross-reference phase for the ion stage.

The immediate goal is not to implement a new simulation feature. The immediate
goal is to validate the completed Python ion driver against the corresponding
legacy MATLAB ion-propagation behavior using the smallest useful deterministic
reference case first.

## Currently pending

1. Ion-stage MATLAB/Python cross-reference validation. Targets 1–4
   (deterministic) are done in
   `scripts/cross_reference/ion_t0_state/` and
   `scripts/cross_reference/ion_multistep_no_collision/`. Target 5
   (stochastic forced-event driver-level comparison) is planned but
   not yet implemented.
2. Step 12: `scripts/run_single_pulse.py`
3. Step 13: HeDFT loading and trajectory comparison in `postprocess/`

## Recommended next task

Design and implement the first focused MATLAB/Python cross-reference test for
the ion stage.

The first validation should be small and deterministic. Prefer, in order:

- ion `t=0` state copied from a tiny neutral checkpoint,
- deterministic one-step ion propagation with collisions disabled,
- deterministic multi-step ion propagation with collisions disabled,
- energy-bookkeeping sanity in a controlled deterministic case.

Do not start with a full stochastic simulation comparison. Too many effects are
entangled at once.

## Files to inspect for the current phase

Start with only the directly relevant Python files:

- `i2_helium_md/simulation/ion.py`
- `i2_helium_md/simulation/ion_initial_state.py`
- `i2_helium_md/simulation/ion_propagation_step.py`
- `i2_helium_md/simulation/checkpoint.py`
- `tests/test_ion.py`
- `tests/test_ion_initial_state.py`
- `tests/test_ion_propagation_step.py`

Then inspect the relevant MATLAB ion-stage files in
`legacy_matlab_repository/`, especially:

- `vmi_sim_3d_ion_propa.m`
- `frog_step_ion.m`
- `ion_interaction_potential.m`
- `add_partner_interaction_ion.m`
- `droplet_potential.m`

## Known MATLAB bugs not to reproduce

The cross-reference phase must not force Python to match known MATLAB
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

- Do not implement `scripts/run_single_pulse.py` before the first ion-stage
  MATLAB/Python cross-reference is planned or completed.
- Do not implement HeDFT comparison before the ion stage has been validated.
- Do not implement plotting before the full simulation path is executable and
  validated.
- Do not refactor the neutral or ion drivers unless a test reveals a real bug.
- Do not implement out-of-scope MATLAB features.
