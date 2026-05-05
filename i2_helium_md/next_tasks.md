# Next tasks

This file defines the active task order for the current MATLAB-to-Python migration phase.

Claude Code should focus on the current task only unless the user explicitly changes priorities.

## Current priority — Step 11d: Implement full ion propagation driver

### Goal

Implement the full ion-stage driver in:

```text
i2_helium_md/simulation/ion.py
```

The driver should mirror the structure of:

```text
i2_helium_md/simulation/neutral.py
```

but use the ion-specific state builder and propagation step.

Do not move on to later tasks automatically. After Step 11d is implemented and tested, stop and report the result.

## Relevant MATLAB source

Before editing Python code, inspect the corresponding MATLAB ion-stage implementation in:

```text
legacy_matlab_repository/
```

Likely relevant MATLAB files include:

```text
vmi_sim_3d_ion_propa.m
frog_step_ion.m
ion_interaction_potential.m
add_partner_interaction_ion.m
droplet_potential.m
```

If the exact file paths differ, search inside `legacy_matlab_repository/` using function names, variable names, or distinctive formulas.

Useful search terms:

```text
vmi_sim_3d_ion_propa
E_kin_ion
E_pot_ion
E_dissip
mass_i
mass_attachment
p_attach
sigma_dependent_on_v
ion_propagation_checkpoint
```

## Relevant Python files to inspect first

Before implementing, inspect:

```text
i2_helium_md/simulation/neutral.py
i2_helium_md/simulation/initial_state.py
i2_helium_md/simulation/propagation_step.py
i2_helium_md/simulation/ion_initial_state.py
i2_helium_md/simulation/ion_propagation_step.py
i2_helium_md/simulation/checkpoint.py
i2_helium_md/simulation/run_directory.py
```

Also inspect the relevant tests for:

```text
neutral propagation
ion initial state
ion propagation step
checkpoint I/O
run directory I/O
```

## Likely Python files affected

The implementation should be narrow.

Likely affected files:

```text
i2_helium_md/simulation/ion.py
tests/test_ion.py
```

Possibly affected files:

```text
i2_helium_md/simulation/__init__.py
tests/smoke_test_ion.py
```

Do not touch unrelated modules unless a real blocker is found and explained.

## Required existing building blocks

Reuse the existing implementation pieces:

```text
build_initial_ion_state
ion_propagation_step
IonCheckpoint
RunDirectory
```

The ion driver should not duplicate ion-step physics that already lives in `ion_propagation_step`.

## Required behavior

The ion driver should:

- build an initial `IonCheckpoint` from a completed `NeutralCheckpoint`,
- run internal ion timesteps using `ion_propagation_step`,
- preserve the internal timestep even when storage is strided,
- write stored states into the checkpoint,
- store the final state even if the final internal step is not aligned with the storage stride,
- correctly populate `mass_history_kg`,
- correctly populate `E_kin_eV`,
- correctly populate `E_pot_eV`,
- correctly populate `E_dissip_eV`,
- preserve reproducibility with a fixed RNG seed,
- optionally save via `RunDirectory`,
- avoid `os.chdir()`,
- avoid hardcoded absolute paths.

## Constraints

Do not:

- change checkpoint schema,
- change physical constants,
- change random-number draw order unless clearly justified and tested,
- refactor neutral propagation unless a real blocker is found,
- implement unsupported ion features,
- implement MATLAB cross-reference tests yet,
- implement plotting,
- implement `scripts/run_single_pulse.py`,
- implement HeDFT postprocessing,
- add large binary reference files,
- perform broad architecture cleanup.

## Known MATLAB bugs not to reproduce

Do not reintroduce known MATLAB bookkeeping bugs just to match legacy output.

Known intentional Python corrections include:

- neutral-stage `E_pot` at `t=0` includes the partner Morse contribution,
- ion-stage `E_kin` at `t=0` fixes the MATLAB velocity-expression bug,
- ion-stage `E_pot` at `t=0` fixes the MATLAB radial-coordinate bug,
- ion-stage `E_pot` at `t=0` includes the partner Coulomb contribution,
- Python uses more accurate physical constants than rounded MATLAB constants.

The ion driver should preserve these existing Python corrections.

## Test expectations for Step 11d

Add or update focused ion-driver tests.

At minimum, tests should cover:

- a small ion run completes,
- output shapes match `IonCheckpoint`,
- `time_ps` is consistent with stored columns,
- final state is stored,
- `mass_history_kg` is populated consistently,
- fixed-seed runs are reproducible,
- RunDirectory save/load roundtrip works,
- unsupported scope flags raise clear errors if relevant,
- energy bookkeeping is sane for a small deterministic or controlled case.

Suggested test file:

```text
tests/test_ion.py
```

If a smoke test is useful, add:

```text
tests/smoke_test_ion.py
```

## Suggested implementation plan

Use this as a starting plan, but inspect the code before finalizing it.

1. Inspect `run_neutral_propagation` in `i2_helium_md/simulation/neutral.py`.
2. Inspect the ion step-state API in `i2_helium_md/simulation/ion_propagation_step.py`.
3. Inspect the ion checkpoint schema in `i2_helium_md/simulation/checkpoint.py`.
4. Inspect `RunDirectory` save/load conventions.
5. Inspect the MATLAB ion propagation loop.
6. Define a public ion driver API that is consistent with the neutral driver.
7. Implement the ion driver with minimal new logic.
8. Add focused tests.
9. Run the narrowest relevant tests first.
10. Run broader tests if the narrow tests pass.
11. Report files inspected, MATLAB files inspected, files changed, tests run, and remaining risks.

## Acceptance criteria

Step 11d is complete when:

- the ion driver exists,
- a small neutral-to-ion pipeline can run through the ion stage,
- output checkpoint fields have correct shapes,
- stored time columns are consistent,
- final state storage is guaranteed,
- mass history is populated,
- fixed-seed behavior is reproducible,
- RunDirectory roundtrip works,
- focused ion-driver tests pass,
- no out-of-scope features were implemented,
- no checkpoint schema or physical constants were changed.

## Explicit stop point

After Step 11d is complete, stop.

Do not proceed to:

- MATLAB cross-reference tests,
- `scripts/run_single_pulse.py`,
- HeDFT loader,
- trajectory comparison,
- plotting,
- analytical force cleanup.

The user will explicitly approve the next task.
