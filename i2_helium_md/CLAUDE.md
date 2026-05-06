# CLAUDE.md

This project is a Python port of a legacy MATLAB molecular-dynamics codebase for iodine / iodine-ion dynamics in helium nanodroplets.

Python package:

```text
i2_helium_md
```

Legacy MATLAB reference:

```text
legacy_matlab_repository/
```

## Current task

The neutral and ion propagation drivers are implemented.

The current task is **MATLAB/Python cross-reference validation of the ion driver**.

Do not continue with broad migration work, plotting, HeDFT loading, full single-pulse scripts, analytical-force cleanup, or refactors unless the user explicitly asks.

## Main working mode

Work script-first.

The goal is to produce small, deterministic comparison results between MATLAB and Python.

Prefer this workflow:

1. Inspect the relevant MATLAB ion code.
2. Inspect the relevant Python ion code.
3. Identify the smallest deterministic comparison case.
4. Write a minimal MATLAB script that exports reference values.
5. Write a minimal Python script that produces the corresponding values.
6. Write a comparison script that loads both outputs and reports differences.
7. Only after the comparison is clear, decide whether a formal pytest regression test is useful.

Do not start by adding large test infrastructure.

Do not start with a full stochastic trajectory comparison.

## First validation targets

Validate in this order:

1. Ion `t=0` state copied from a tiny neutral checkpoint.
2. One deterministic ion step with collisions disabled.
3. Several deterministic ion steps with collisions disabled.
4. Energy bookkeeping in deterministic mode.
5. Collision/statistical behavior only after deterministic comparisons are stable.

## Relevant files

Start with these Python files:

```text
i2_helium_md/simulation/ion.py
i2_helium_md/simulation/ion_initial_state.py
i2_helium_md/simulation/ion_propagation_step.py
i2_helium_md/simulation/checkpoint.py
```

Start with these MATLAB files or search for them under `legacy_matlab_repository/`:

```text
vmi_sim_3d_ion_propa.m
frog_step_ion.m
ion_interaction_potential.m
add_partner_interaction_ion.m
droplet_potential.m
```

Useful search terms:

```text
E_kin_ion
E_pot_ion
E_dissip
mass_i
p_attach
sigma_dependent_on_v
ion_propagation_checkpoint
```

Read other files only if they are needed for the current comparison.

Do not read the full `migration_log.md` unless a documented deviation or historical decision is needed.

## MATLAB/Python deviation policy

MATLAB is the legacy reference, but known MATLAB bugs should not be reproduced blindly.

Known intentional Python corrections:

- neutral-stage `E_pot` at `t=0` includes the partner Morse contribution,
- ion-stage `E_kin` at `t=0` fixes the MATLAB velocity-expression bug,
- ion-stage `E_pot` at `t=0` fixes the MATLAB radial-coordinate bug,
- ion-stage `E_pot` at `t=0` includes the partner Coulomb contribution,
- Python uses more accurate physical constants than rounded MATLAB constants.

Every comparison result must state whether Python is expected to:

1. match MATLAB,
2. intentionally differ because Python fixes a MATLAB bug,
3. differ slightly because constants differ,
4. differ statistically because stochastic streams differ.

## Script rules

Comparison scripts should be small, readable, and disposable if needed.

Preferred locations:

```text
scripts/cross_reference/
tests/reference/
```

Preferred outputs:

```text
small .csv
small .json
small .npz
```

Avoid:

```text
large .mat files
large .npz checkpoints
generated figures
full simulation-output directories
temporary debug dumps
```

Each reference script should document:

- command used to run it,
- number of molecules,
- number of timesteps,
- random seed if applicable,
- enabled/disabled physics features,
- whether known MATLAB bugs are present in the exported data.

## Editing limits

Before editing, run:

```bash
git status
```

If the working tree is not clean, stop and ask.

During this phase, do not change unless a comparison reveals a real bug:

```text
checkpoint schema
physical constants
neutral propagation
plotting or postprocessing
single-pulse public script
HeDFT loader
```

Keep changes local to the comparison scripts or the smallest required ion-driver fix.

## Reporting format

For each cross-reference attempt, report:

1. MATLAB files inspected,
2. Python files inspected,
3. scripts written or changed,
4. comparison case,
5. expected match or intentional deviation,
6. numerical comparison result,
7. files changed,
8. remaining risks.

After one comparison task is complete, stop and report. Do not automatically move to the next migration step.
