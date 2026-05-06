# i2_helium_md

Modern Python port of Michael Stadlhofer's MATLAB molecular-dynamics code
for simulating iodine in helium nanodroplets (single-pulse,
HeDFT-comparison scope).

## Migration progress

| # | Python module | MATLAB source | Status |
|---|---------------|---------------|--------|
| 1 | `physics/constants.py` | `physical_constants.m` | done |
| 2 | `config.py` (SimConfig) | ~36 MATLAB globals | done |
| 3 | `presets.py` | `inputfiles_dft_comparison/single_pulse_N2000.m` | done |
| 4 | `physics/potentials.py` | `droplet_potential.m`, `get_morse_potential_X.m`, `get_morse_potential_I2plus.m`, `morse_potential_I2plus_state_select.m` | done |
| 5 | `physics/interactions.py` | `atom_interaction_potential.m`, `ion_interaction_potential.m`, `add_partner_interaction.m`, `add_partner_interaction_ion.m` | done |
| 6 | `physics/leapfrog.py` | `frog_step_neutral.m`, `frog_step_ion.m` | done |
| 7 | `sampling/droplet_sizes.py` | `generate_droplet_sizes.m`, `get_dropletsize.m` | done |
| 8 | `sampling/radial_positions.py` | `generate_radial_samples_3d.m` | done |
| 9 | `simulation/checkpoint.py`, `simulation/run_directory.py` | `save('neutral_propagation_checkpoint', ...)` | done |
| 10 | `simulation/neutral.py`, `sampling/orientations.py`, `physics/collisions.py`, `simulation/initial_state.py`, `simulation/propagation_step.py` | `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` | done |
| 11 | `simulation/ion.py` | `vmi_sim_3d_ion_propa.m` | done; MATLAB/Python cross-reference complete |
| 12 | `scripts/run_single_pulse.py` | `run_simulation.m` | done |
| 13 | `postprocess/hedft_loader.py` + `compare_trajectories.py` | `simulation_image_only_trajectories.m` | pending |

## Current phase

The neutral and ion propagation drivers are implemented. Ion-stage
MATLAB/Python cross-reference validation is complete.

The public single-pulse run script is implemented. The current focus is Step
13: HeDFT loading and trajectory comparison in `postprocess/`.

Completed ion cross-reference artifacts:

- `scripts/cross_reference/ion_t0_state/`
- `scripts/cross_reference/ion_multistep_no_collision/`
- `scripts/cross_reference/ion_stochastic_forced/`

## Documentation

- `docs/physics_background.md` — physical model, potentials, and design rationale
- `docs/constants_module.md` — walkthrough of the `constants.py` module
- `docs/interactions_module.md` — walkthrough of the `interactions.py` module
- `docs/leapfrog_module.md` — walkthrough of the `leapfrog.py` integrator
- `docs/droplet_sizes_module.md` — walkthrough of the `droplet_sizes.py` sampler
- `docs/droplet_sizes_diagnostics_module.md` — debugging plots for the pickup simulation
- `docs/radial_positions_module.md` — walkthrough of the `radial_positions.py` sampler
- `docs/checkpoint_module.md` — walkthrough of the `checkpoint.py` I/O module
- `docs/run_directory_module.md` — walkthrough of the `RunDirectory` convention layer
- `docs/orientations_module.md` — walkthrough of the `orientations.py` angular sampler
- `docs/initial_state_module.md` — walkthrough of `build_initial_state` (Step 10c-i)
- `docs/propagation_step_module.md` — walkthrough of `neutral_propagation_step` (Step 10c-ii)
- `docs/neutral_module.md` — walkthrough of `run_neutral_propagation` driver (Step 10c-iii)
- `docs/ion_initial_state_module.md` — walkthrough of `build_initial_ion_state` (Step 11b)
- `docs/ion_propagation_step_module.md` — walkthrough of `ion_propagation_step` (Step 11c)
- `docs/collisions_module.md` — walkthrough of the `collisions.py` hard-sphere physics
- `migration_log.md` — chronological record of decisions, deviations, and open questions
- `current_state.md` — completed modules and current validation phase
- `next_tasks.md` — task list and acceptance criteria for upcoming work
- `testing.md` — testing conventions, tolerances, and MATLAB cross-reference rules
- `agent_protocol.md` — investigation vs. edit-mode rules for collaborators

## Project layout

```text
i2_helium_md_py/
├── data/reference/              data files copied from legacy repo
├── docs/                        physics background and module docs
├── i2_helium_md/                Python package
│   ├── config.py                SimConfig dataclass
│   ├── presets.py               preset builders
│   ├── physics/                 constants, potentials, interactions, integrators
│   ├── sampling/                random samplers
│   ├── simulation/              neutral and ion propagation
│   └── postprocess/             HeDFT comparison utilities
├── legacy_matlab_repository/    original MATLAB reference implementation
├── scripts/                     entry points
├── tests/                       smoke tests and pytest suite
└── pyproject.toml
```

## Data files needed in `data/reference/`

Copy these three files from the legacy MATLAB repo:

| Legacy path | → | New path |
|-------------|---|----------|
| `HeDFT_MD_comparison_neutral/custom_start_interpolating_functions.mat` | → | `data/reference/hedft_custom_start.mat` |
| `single_pulse_simulation/HeDFT_comparison/9Angström/data_vabs2.csv` | → | `data/reference/hedft_9A_velocity.csv` |
| `single_pulse_simulation/HeDFT_comparison/9Angström/R1-R2.csv` | → | `data/reference/hedft_9A_distance.csv` |

Once copied, they are referenced via `SimConfig.data_dir`; no hardcoded paths
should be used.

## Quickstart

Run a small single-pulse pipeline from the repository root:

```bash
python scripts/run_single_pulse.py --run-dir results/single_pulse_test --num-molecules 10 --seed 123
```

This writes `cfg.json`, `neutral.npz`, and `ion.npz` into the run directory.
Use `--force` to intentionally overwrite an existing run directory.

Direct function calls are also supported:

```python
from i2_helium_md import single_pulse_N2000
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.ion import run_ion_propagation

cfg = single_pulse_N2000(num_molecules=500, seed=123)
neutral_result = run_neutral_propagation(cfg)
ion_result = run_ion_propagation(cfg, neutral_result)
```

If the exact public ion-driver function name differs, update this quickstart
to match `i2_helium_md/simulation/ion.py`.

## Scope decisions agreed with user

In scope:

- single-pulse neutral and ion dynamics,
- 9 Å HeDFT comparison,
- MATLAB/Python reference validation.

Out of scope:

- pump-probe,
- effusive dynamics,
- experimental VMI comparison,
- 18 Å HeDFT comparison,
- Abel inversion,
- image-processing utilities.
