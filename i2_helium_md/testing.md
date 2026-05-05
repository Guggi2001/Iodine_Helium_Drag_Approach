# Testing guide

This project is a scientific MATLAB-to-Python migration.

Tests are not only software checks. They document physical conventions, unit choices, numerical tolerances, MATLAB/Python deviations, and known legacy bugs.

The original MATLAB source is available in:

```text
`legacy\_matlab\_repository/`
```

Use it when designing reference tests or changing migrated behavior.

## Default command

Run the full test suite with:

```bash
pytest
```

## Module-specific tests

Prefer narrow tests while developing.

Examples:

```bash
pytest tests/test\_constants.py
pytest tests/test\_potentials.py
pytest tests/test\_interactions.py
pytest tests/test\_leapfrog.py
pytest tests/test\_ion\_propagation\_step.py
pytest tests/test\_ion.py
```

Not every listed file must exist at all times.

Use the closest matching test file for the module being edited.



## Before editing

Before editing any code, check repository state:

```bash
git status
```

If the working tree is not clean, stop and ask the user before editing.

## After editing

After editing code:

1. Run the narrowest relevant test first.
2. Run broader tests if the change affects multiple modules.
3. Report exactly which tests were run.
4. Report failures honestly.
5. Do not claim correctness without tests or a focused numerical check.

## Numerical tolerance rules

Do not choose tolerances arbitrarily.

A test tolerance must be justified by the numerical or physical situation:

* direct analytical formula ports should use tight tolerances,
* finite-difference forces need tolerances consistent with the finite-difference step,
* Monte-Carlo samplers need statistical tolerances based on sample size,
* MATLAB-vs-Python comparisons may need tolerances for known constant updates,
* known legacy MATLAB bugs should not be preserved unless explicitly required,
* stochastic tests should check distributions or moments, not exact samples unless RNG identity is guaranteed.

If a tolerance is loose, explain why in a test comment.

## MATLAB cross-reference tests

When comparing against MATLAB reference results, tests should state clearly whether they check:

1. exact numerical agreement,
2. qualitative physical behavior,
3. agreement within expected constant/unit differences,
4. agreement after intentionally fixing a known MATLAB bug,
5. statistical agreement for stochastic routines.

The test name and comments should make the comparison target clear.

## MATLAB reference-data rules

Reference data should be small, inspectable, and reproducible.

Preferred formats:

* small `.csv`,
* small `.json`,
* small `.npz`,
* short text files documenting MATLAB output.

Avoid committing:

* large full-run checkpoints,
* large `.mat` files,
* generated figures,
* temporary debugging outputs,
* full simulation-output directories.

If reference data is generated from MATLAB, document:

* MATLAB script or command used,
* relevant input parameters,
* random seed if applicable,
* number of molecules,
* number of timesteps,
* disabled/enabled physics features,
* whether the data includes a known MATLAB bug or an intentional correction.

## Recommended MATLAB/Python cross-reference sequence

For a newly migrated simulation component, validate in this order:

1. Direct formula comparison.
2. Shape and unit checks.
3. One-step deterministic comparison.
4. Multi-step deterministic comparison.
5. Energy bookkeeping comparison.
6. Stochastic statistical comparison.
7. Full driver smoke test.

Do not start with a full stochastic simulation comparison. Too many effects are entangled.

## Deterministic tests

When possible, disable stochastic features first.

Useful deterministic settings include:

* collisions disabled,
* mass attachment disabled,
* fixed initial state,
* fixed molecule count,
* very few timesteps,
* fixed droplet radius,
* fixed random seed.

Deterministic tests are preferred for:

* force calculation,
* leapfrog stepping,
* energy bookkeeping,
* checkpoint continuity,
* initial-state transfer from neutral to ion stage.

## Stochastic tests

For stochastic code, exact trajectory matching is usually fragile unless MATLAB and Python share identical random-number streams.

Prefer tests of:

* empirical collision rate,
* mean energy loss,
* mass-attachment rate,
* sampled distribution moments,
* isotropy measures,
* reproducibility within Python for a fixed seed.

Use enough samples for statistical tests, but keep the test suite fast.

## After changing simulation flow

For changes to drivers such as neutral or ion propagation, test at least:

* output shapes,
* time-axis consistency,
* reproducibility with fixed seed,
* checkpoint save/load roundtrip,
* final state storage,
* energy bookkeeping,
* behavior with small `num\_molecules`,
* behavior with strided storage if applicable.

## Known MATLAB bugs to guard against

Do not reintroduce these bugs:

* neutral-stage `E\_pot` at `t=0` omitted partner Morse contribution,
* ion-stage `E\_kin` at `t=0` used an incorrect velocity expression,
* ion-stage `E\_pot` at `t=0` omitted the `z` coordinate,
* ion-stage `E\_pot` at `t=0` omitted the partner Coulomb term.

Tests should explicitly guard against these if the relevant code is touched.

## What not to do

Do not update expected values blindly.

If a test fails after a code change, first determine whether the failure means:

* the code is wrong,
* the test encoded a legacy MATLAB bug,
* the tolerance was physically unreasonable,
* the intended model changed and documentation must be updated,
* a MATLAB/Python constant difference is expected,
* a stochastic test is under-sampled,
* or the reference data was generated from the wrong MATLAB code path.

