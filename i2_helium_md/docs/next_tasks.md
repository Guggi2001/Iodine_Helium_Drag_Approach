# Next tasks

This file defines the intended task order for the remaining MATLAB-to-Python
migration. Claude Code should follow this order unless the user explicitly
changes priorities.

## Task 1 — Implement full ion propagation driver

### Goal

Implement the full ion-stage driver in `simulation/ion.py`.

The driver should mirror the structure of `simulation/neutral.py`, but use
the ion-specific state builder and propagation step.

### Likely files affected

- `i2_helium_md/simulation/ion.py`
- `i2_helium_md/simulation/__init__.py`
- `tests/test_ion.py` or `tests/test_ion_driver.py`
- optional smoke test

### Required building blocks

- `build_initial_ion_state`
- `ion_propagation_step`
- `IonCheckpoint`
- `RunDirectory`

### Constraints

- Preserve internal timestep even when storage is strided.
- Store the final state even if the final internal step is not aligned with the stride.
- Correctly populate `mass_history_kg`.
- Correctly populate `E_kin_eV`, `E_pot_eV`, and `E_dissip_eV`.
- Preserve reproducibility with fixed RNG seed.
- Do not implement out-of-scope ion features.
- Do not change checkpoint schema unless explicitly requested.
- Do not change physical constants.
- Do not refactor neutral propagation unless a test failure proves it is necessary.

### Acceptance criteria

- A small ion run completes.
- Output shapes match the `IonCheckpoint` schema.
- `time_ps` is consistent with stored columns.
- Final state is stored.
- `mass_history_kg` is populated consistently.
- Fixed-seed runs are reproducible.
- RunDirectory save/load roundtrip works.
- Relevant ion tests pass.

---

## Task 2 — Design MATLAB cross-reference tests for the ion stage

### Goal

After the Python ion driver exists, design focused tests that compare the
Python ion-stage behavior against MATLAB reference output or MATLAB-derived
reference values.

This should be analogous to the validation strategy used for the neutral
stage: compare the new Python implementation against the legacy MATLAB model
where appropriate, while explicitly documenting intentional deviations caused
by fixed MATLAB bugs or more accurate constants.

### Important principle

Do not blindly force Python to match MATLAB if MATLAB is known to be wrong.

Known fixed MATLAB bookkeeping bugs include:

- ion-stage `E_kin` at `t=0` used an incorrect velocity expression,
- ion-stage `E_pot` at `t=0` omitted the `z` coordinate and omitted the partner Coulomb term,
- neutral-stage `E_pot` at `t=0` omitted the partner Morse contribution.

The cross-reference tests should distinguish between:

1. behavior that should match MATLAB,
2. behavior that should intentionally differ because Python fixes a known MATLAB bug,
3. behavior that may differ slightly because Python uses more accurate physical constants.

### Likely files affected

- `tests/test_ion_matlab_reference.py`
- `tests/reference/`
- possibly `scripts/export_matlab_reference_*` or documentation describing how reference data was generated

### Possible reference data

Use one or more small reference cases exported from MATLAB:

- tiny `N`, very few timesteps,
- fixed seed if possible,
- collision-disabled deterministic case,
- mass-attachment-disabled case,
- single-pulse case,
- ion propagation from a stored neutral checkpoint.

The reference data should be small enough to commit to the repository if it is text or compact NumPy data. Large generated checkpoints should not be committed.

### Suggested test categories

#### 1. Initial ion state comparison

Check quantities derived from the neutral checkpoint:

- initial positions,
- initial velocities,
- initial masses,
- initial kinetic energy,
- initial ion-droplet potential,
- initial partner Coulomb contribution.

For known MATLAB t=0 bookkeeping bugs, write tests that document the intentional Python correction instead of matching the wrong MATLAB value.

#### 2. Deterministic no-collision propagation

Use a configuration where collisions and mass attachment are disabled.

Compare:

- positions after a small number of steps,
- velocities after a small number of steps,
- kinetic energy,
- potential energy,
- total energy drift.

This is the cleanest test of the leapfrog + force model.

#### 3. Collision probability sanity check

For collision-active runs, exact trajectory matching may be difficult if RNG streams differ between MATLAB and Python.

Instead compare controlled statistical quantities:

- collision fraction,
- velocity-dependent cross-section behavior,
- mean energy loss per collision,
- mass attachment rate if enabled.

#### 4. Checkpoint continuity

Verify that the last neutral state used as the ion initial state is copied
correctly into the ion checkpoint at `t=0`.

#### 5. Energy bookkeeping

Check that:

```text
E_total = E_kin + E_pot + E_dissip