# CLAUDE.md

This project is a Python port of a legacy MATLAB molecular-dynamics codebase for iodine / iodine-ion dynamics in helium nanodroplets.

The Python package is `i2_helium_md`.

The original MATLAB code is available inside this repository under:

```text
legacy_matlab_repository/
```

The MATLAB repository is the primary reference implementation for migration, debugging, and cross-checking.

## Immediate priority

The current task is Step 11d only:

```text
Implement the full ion propagation driver in i2_helium_md/simulation/ion.py.
```

Do not proceed to MATLAB cross-reference tests, run scripts, HeDFT postprocessing, plotting, or broader cleanup unless the user explicitly asks.

After Step 11d is implemented and tested, stop and report the result. The user will decide the next task.

## Mandatory rule

Before changing migrated Python code, inspect the relevant MATLAB source file whenever a MATLAB source exists.

Do not rely only on documentation if the original MATLAB implementation is available. Documentation summarizes decisions; MATLAB shows the legacy behavior.

When Python intentionally differs from MATLAB, document and test the difference.

## Read before editing

For any coding task, read:

1. `README.md`
2. `current_state.md`
3. `next_tasks.md`
4. `testing.md`
5. `agent_protocol.md`
6. The relevant module documentation in `docs/`
7. The relevant Python tests in `tests/`
8. The corresponding MATLAB source in `legacy_matlab_repository/`

Use `migration_log.md` when historical decisions, known deviations, or unresolved issues matter.

For the immediate Step 11d task, especially inspect:

```text
i2_helium_md/simulation/neutral.py
i2_helium_md/simulation/initial_state.py
i2_helium_md/simulation/propagation_step.py
i2_helium_md/simulation/ion_initial_state.py
i2_helium_md/simulation/ion_propagation_step.py
i2_helium_md/simulation/checkpoint.py
i2_helium_md/simulation/run_directory.py
```

Also inspect the relevant tests and MATLAB ion propagation source.

## Current state

- Neutral propagation is complete.
- Ion propagation has 11a, 11b, and 11c complete.
- Step 11d, the full ion propagation driver, is pending.
- `scripts/run_single_pulse.py` is still pending.
- HeDFT postprocessing/comparison is still pending.

## Current scope

In scope now:

- implementing the full ion propagation driver,
- adding focused ion-driver tests,
- checking consistency with the existing neutral driver design,
- checking consistency with the relevant MATLAB ion propagation code,
- updating documentation only if needed to reflect Step 11d.

Out of scope for the current task unless explicitly requested:

- MATLAB cross-reference test implementation,
- full single-pulse run script,
- HeDFT loader,
- HeDFT trajectory comparison,
- plotting,
- pump-probe,
- effusive dynamics,
- experimental VMI comparison,
- Abel inversion,
- 18 Å HeDFT comparison,
- image-processing utilities,
- broad architecture rewrites,
- hardcoded absolute Windows paths.

## Unit conventions

Unless stated otherwise:

- Length: Å
- Time: ps
- Velocity: Å/ps
- Energy: eV
- Mass: kg in config/checkpoints unless explicitly `_u`
- Force: eV/Å
- Acceleration: Å/ps²
- Temperature: K

Unit conversions must be explicit and centralized in `i2_helium_md/physics/constants.py`.

## MATLAB cross-reference protocol

When editing a migrated module, identify and inspect the corresponding MATLAB source.

For Step 11d, inspect the MATLAB ion propagation implementation first.

Likely relevant MATLAB files include:

```text
vmi_sim_3d_ion_propa.m
frog_step_ion.m
ion_interaction_potential.m
add_partner_interaction_ion.m
droplet_potential.m
```

If the exact file path differs, search inside `legacy_matlab_repository/` by function name, variable name, or formula.

Useful search terms for Step 11d:

```text
vmi_sim_3d_ion_propa
E_kin_ion
E_pot_ion
mass_i
sigma_dependent_on_v
p_attach
E_dissip
ion_propagation_checkpoint
```

Common mappings:

```text
physical_constants.m -> i2_helium_md/physics/constants.py
droplet_potential.m, get_morse_potential_*.m -> i2_helium_md/physics/potentials.py
atom_interaction_potential.m, ion_interaction_potential.m -> i2_helium_md/physics/interactions.py
frog_step_neutral.m, frog_step_ion.m -> i2_helium_md/physics/leapfrog.py
generate_droplet_sizes.m, get_dropletsize.m -> i2_helium_md/sampling/droplet_sizes.py
generate_radial_samples_3d.m -> i2_helium_md/sampling/radial_positions.py
vmi_sim_3d_neutral_propa_HeDFT_mimic.m -> i2_helium_md/simulation/neutral.py and related neutral modules
vmi_sim_3d_ion_propa.m -> i2_helium_md/simulation/ion.py and related ion modules
run_simulation.m -> scripts/run_single_pulse.py
simulation_image*.m -> i2_helium_md/postprocess/
```

Always report which MATLAB files were inspected.

## Known intentional Python deviations

Do not reintroduce known MATLAB bugs for byte-identical output.

Known intentional fixes include:

- neutral-stage `E_pot` at `t=0` includes partner Morse contribution,
- ion-stage `E_kin` at `t=0` fixes the MATLAB velocity-expression bug,
- ion-stage `E_pot` at `t=0` fixes the MATLAB radial-coordinate bug and includes partner Coulomb contribution,
- Python uses more accurate physical constants than rounded MATLAB constants.

For Step 11d, preserve these existing Python corrections. Do not change them unless the user explicitly asks.

## Step 11d implementation constraints

The ion driver should mirror the design of `run_neutral_propagation` where appropriate.

Use the existing building blocks:

```text
build_initial_ion_state
ion_propagation_step
IonCheckpoint
RunDirectory
```

Implementation constraints:

- Do not change checkpoint schema.
- Do not change physical constants.
- Do not refactor neutral propagation unless a real blocker is found.
- Do not implement unsupported ion features.
- Do not implement MATLAB cross-reference tests yet.
- Do not implement plotting.
- Do not implement `scripts/run_single_pulse.py` yet.
- Do not implement HeDFT postprocessing.
- Preserve internal timestep even if checkpoint storage is strided.
- Store the final state even if the final internal step is not aligned with the storage stride.
- Correctly populate `mass_history_kg`.
- Correctly populate `E_kin_eV`, `E_pot_eV`, and `E_dissip_eV`.
- Preserve reproducibility with a fixed seed.
- Use `RunDirectory` for optional run-directory saving.
- Do not use `os.chdir()`.
- Do not hardcode absolute Windows paths.

## Step 11d test expectations

Add or update focused tests for the ion driver.

At minimum, tests should cover:

- small ion run completes,
- output shapes match `IonCheckpoint`,
- `time_ps` is consistent with stored columns,
- final state is stored,
- `mass_history_kg` is populated consistently,
- fixed-seed runs are reproducible,
- RunDirectory save/load roundtrip works,
- unsupported scope flags raise clear errors if relevant,
- energy bookkeeping is sane for a small deterministic or controlled case.

Run the narrowest relevant tests first, then broader tests if needed.

## Testing

Before editing:

```bash
git status
```

If the working tree is not clean, stop and ask before editing.

After editing, run the narrowest relevant tests first:

```bash
pytest tests/test_<relevant_module>.py
```

For Step 11d, likely relevant tests include:

```bash
pytest tests/test_ion_initial_state.py
pytest tests/test_ion_propagation_step.py
pytest tests/test_checkpoint.py
pytest tests/test_run_directory.py
pytest tests/test_ion.py
```

Not every listed test file must already exist. Add the focused ion-driver test file if needed.

Then run broader tests if needed:

```bash
pytest
```

If pytest is unavailable, inspect and run the relevant smoke test if one exists.

Do not claim correctness unless tests passed, a focused numerical check was run, or the task was inspection-only.

Numerical tolerances must be justified by formula agreement, MATLAB/Python deviation, finite-difference error, statistical uncertainty, physical approximation, or updated constants.

## Working modes

Investigation mode:

- do not edit files,
- inspect docs, tests, and MATLAB source,
- report findings and uncertainties,
- suggest the smallest safe next step.

Planning mode:

- do not edit files,
- inspect all relevant Python and MATLAB sources,
- identify ambiguities,
- ask clarification questions if needed,
- propose the smallest safe implementation plan.

Edit mode:

- make the smallest coherent change,
- avoid unrelated files,
- inspect MATLAB source first,
- add or update tests,
- run relevant tests,
- report files changed and remaining risks.

## Forbidden without explicit approval

Do not:

- delete reference data,
- delete the MATLAB repository,
- change physical constants,
- change checkpoint schema,
- change random-number draw order,
- perform broad refactors,
- optimize by changing numerical behavior,
- implement out-of-scope features,
- add large binary reference files,
- implement MATLAB cross-reference tests before the ion driver is complete and explicitly approved,
- move on to Step 12 or Step 13 automatically.

## Reporting format

For coding tasks, report:

1. files inspected,
2. MATLAB files inspected,
3. implementation plan or root cause,
4. files changed,
5. tests run,
6. numerical checks, if any,
7. intentional MATLAB/Python deviations,
8. remaining risks.

For Step 11d, explicitly state whether the ion driver is complete and what still remains before moving to the next task.
