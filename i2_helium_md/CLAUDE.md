# CLAUDE.md

This project is a Python port of a legacy MATLAB molecular-dynamics codebase for iodine / iodine-ion dynamics in helium nanodroplets.

The Python package is `i2\_helium\_md`.

The original MATLAB code is available inside this repository under:

```text
legacy\_matlab\_repository/
```

The MATLAB repository is the primary reference implementation for migration, debugging, and cross-checking.

## Mandatory rule

Before changing migrated Python code, inspect the relevant MATLAB source file whenever a MATLAB source exists.

Do not rely only on documentation if the original MATLAB implementation is available. Documentation summarizes decisions; MATLAB shows the legacy behavior.

When Python intentionally differs from MATLAB, document and test the difference.

## Read before editing

For any coding task, read:

1. `README.md`
2. `current\_state.md`
3. `docs/next\_tasks.md`
4. `docs/testing.md`
5. `docs/agent\_protocol.md`
6. The relevant module documentation in `docs/`
7. The relevant Python tests in `tests/`
8. The corresponding MATLAB source in `legacy\_matlab\_repository/`

Use `migration\_log.md` when historical decisions, known deviations, or unresolved issues matter.

## Current state

* Neutral propagation is complete.
* Ion propagation has 11a, 11b, and 11c complete.
* Step 11d, the full ion propagation driver, is pending.
* After the ion driver, the next priority is MATLAB cross-reference validation for the ion stage.
* `scripts/run\_single\_pulse.py` and HeDFT postprocessing are still pending.

## Scope

In scope:

* single-pulse neutral and ion dynamics,
* 9 Å HeDFT comparison,
* checkpoint-based run directories,
* MATLAB/Python cross-reference validation,
* small reference tests.

Out of scope unless explicitly requested:

* pump-probe,
* effusive dynamics,
* experimental VMI comparison,
* Abel inversion,
* 18 Å HeDFT comparison,
* image-processing utilities,
* broad architecture rewrites,
* hardcoded absolute Windows paths.

## Unit conventions

Unless stated otherwise:

* Length: Å
* Time: ps
* Velocity: Å/ps
* Energy: eV
* Mass: kg in config/checkpoints unless explicitly `\_u`
* Force: eV/Å
* Acceleration: Å/ps²
* Temperature: K

Unit conversions must be explicit and centralized in `i2\_helium\_md/physics/constants.py`.

## MATLAB cross-reference protocol

When editing a migrated module, identify and inspect the corresponding MATLAB source.

Common mappings:

```text
physical\_constants.m -> i2\_helium\_md/physics/constants.py
droplet\_potential.m, get\_morse\_potential\_\*.m -> i2\_helium\_md/physics/potentials.py
run\_simulation.m -> scripts/run\_single\_pulse.py
simulation\_image\*.m -> i2\_helium\_md/postprocess/
```

If a file name differs, search the MATLAB repository by function name, variable name, or formula.

Always report which MATLAB files were inspected.

## Known intentional Python deviations

Do not reintroduce known MATLAB bugs for byte-identical output.

Known intentional fixes include:

* neutral-stage `E\_pot` at `t=0` includes partner Morse contribution,
* ion-stage `E\_kin` at `t=0` fixes the MATLAB velocity-expression bug,
* ion-stage `E\_pot` at `t=0` fixes the MATLAB radial-coordinate bug and includes partner Coulomb contribution,
* Python uses more accurate physical constants than rounded MATLAB constants.

## Testing

Before editing:

```bash
git status
```

If the working tree is not clean, stop and ask before editing.

After editing, run the narrowest relevant tests first:

```bash
pytest tests/test\_<relevant\_module>.py
```

Then run broader tests if needed:

```bash
pytest
```

If pytest is unavailable, inspect and run the relevant smoke test if one exists.

Do not claim correctness unless tests passed, a focused numerical check was run, or the task was inspection-only.

Numerical tolerances must be justified by formula agreement, MATLAB/Python deviation, finite-difference error, statistical uncertainty, physical approximation, or updated constants.

## Working modes

Investigation mode:

* do not edit files,
* inspect docs, tests, and MATLAB source,
* report findings and uncertainties,
* suggest the smallest safe next step.

Edit mode:

* make the smallest coherent change,
* avoid unrelated files,
* inspect MATLAB source first,
* add or update tests,
* run relevant tests,
* report files changed and remaining risks.

## Forbidden without explicit approval

Do not:

* delete reference data,
* delete the MATLAB repository,
* change physical constants,
* change checkpoint schema,
* change random-number draw order,
* perform broad refactors,
* optimize by changing numerical behavior,
* implement out-of-scope features,
* add large binary reference files.

## Reporting format

For coding tasks, report:

1. files inspected,
2. MATLAB files inspected,
3. files changed,
4. tests run,
5. numerical checks, if any,
6. intentional MATLAB/Python deviations,
7. remaining risks.

