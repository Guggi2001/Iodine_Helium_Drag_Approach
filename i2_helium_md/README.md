# i2_helium_md

Modern Python port of Michael Stadlhofer's MATLAB molecular-dynamics code
for simulating iodine in helium nanodroplets (single-pulse,
HeDFT-comparison scope).

## Migration progress

| # | Python module | MATLAB source | Status |
|---|---------------|---------------|--------|
| 1 | `physics/constants.py` | `physical_constants.m` | done |
| 2 | `config.py` (SimConfig) | ~36 MATLAB globals | done |
| 3 | `presets.py` | `inputfiles_dft_comparison/single_pulse_N2000.m`, `inputfiles_dft_comparison/single_pulse_droplet_distribution.m` | done |
| 4 | `physics/potentials.py` | `droplet_potential.m`, `get_morse_potential_X.m`, `get_morse_potential_I2plus.m`, `morse_potential_I2plus_state_select.m` | done |
| 5 | `physics/interactions.py` | `atom_interaction_potential.m`, `ion_interaction_potential.m`, `add_partner_interaction.m`, `add_partner_interaction_ion.m` | done |
| 6 | `physics/leapfrog.py` | `frog_step_neutral.m`, `frog_step_ion.m` | done |
| 7 | `sampling/droplet_sizes.py` | `generate_droplet_sizes.m`, `get_dropletsize.m` | done |
| 8 | `sampling/radial_positions.py` | `generate_radial_samples_3d.m` | done |
| 9 | `simulation/checkpoint.py`, `simulation/run_directory.py` | `save('neutral_propagation_checkpoint', ...)` | done |
| 10 | `simulation/neutral.py`, `sampling/orientations.py`, `physics/collisions.py`, `simulation/initial_state.py`, `simulation/propagation_step.py` | `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` | done |
| 11 | `simulation/ion.py` | `vmi_sim_3d_ion_propa.m` | done; MATLAB/Python cross-reference complete |
| 12 | `scripts/run_single_pulse.py` | `run_simulation.m` | done |
| 13 | `postprocess/hedft_loader.py`, `postprocess/compare_trajectories.py`, `postprocess/velocity_distribution.py`, `scripts/plot_hedft_comparison.py` | `simulation_image_only_trajectories.m`, parts of `simulation_image.m` | done; first plotting path and VMI helpers present |

## Current phase

The neutral and ion propagation drivers are implemented. Ion-stage
MATLAB/Python cross-reference validation is complete.

The public single-pulse run script is implemented. Step 13 now has a first
post-processing path: normalized HeDFT trajectory loading, numerical
distance/velocity comparison against an ion checkpoint, final-velocity
histogram helpers, and a plotting script for the HeDFT comparison figures.

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
- `docs/run_single_pulse_script.md` — usage guide for the public single-pulse script
- `docs/collisions_module.md` — walkthrough of the `collisions.py` hard-sphere physics
- `docs/hedft_loader_module.md` — walkthrough of normalized HeDFT reference loading
- `docs/compare_trajectories_module.md` — walkthrough of numerical MD/HeDFT trajectory comparison
- `docs/velocity_distribution_module.md` — walkthrough of VMI reference and final-velocity histogram helpers
- `docs/plot_hedft_comparison_script.md` — usage guide for the HeDFT comparison plotting script
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

## Reference data in `data/reference/`

The current post-processing code expects normalized reference CSVs:

| File | Purpose |
|------|---------|
| `9A_All_Data.csv` | normalized 8-column 9 A HeDFT trajectory reference |
| `18A_All_Data.csv` | normalized 8-column 18 A HeDFT trajectory reference |
| `vmi_iplus_gas.csv` | exported experimental gas-phase VMI reference |
| `vmi_iplus_he.csv` | exported experimental droplet VMI reference |

The HeDFT loader reads the common columns
`Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`. The earlier split
legacy 9 A files (`data_vabs2.csv` and `R1-R2.csv`) remain useful provenance,
but the Python API now uses the normalized all-data files.

## Quickstart

Run the single-pulse pipeline from the repository root:

```bash
python scripts/run_single_pulse.py
```

Edit the `USER SETTINGS` section at the top of
`scripts/run_single_pulse.py` to choose `INPUT_PRESET`, `RUN_DIR`,
`RUN_SIZE`, `NUM_MOLECULES`, `SEED`, and `ION_TIME_PS`. Supported input
presets currently mirror `single_pulse_N2000.m` and
`single_pulse_droplet_distribution.m`. The script writes `cfg.json`,
`neutral.npz`, and `ion.npz` into the run directory.

See `docs/run_single_pulse_script.md` for the full usage workflow.

Direct function calls are also supported:

```python
from i2_helium_md import single_pulse_N2000, single_pulse_droplet_distribution
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.ion import run_ion_propagation

cfg = single_pulse_N2000(num_molecules=500, seed=123)
# or:
cfg = single_pulse_droplet_distribution(num_molecules=500, seed=123)
neutral_result = run_neutral_propagation(cfg)
ion_result = run_ion_propagation(cfg, neutral_result)
```

Plot an existing production run against the default 9 A HeDFT reference:

```bash
python scripts/plot_hedft_comparison.py
```

Or use the numerical comparison API directly:

```python
from i2_helium_md.postprocess import compare_distance, load_hedft_trajectory
from i2_helium_md.simulation.run_directory import RunDirectory

run = RunDirectory("data/runs/single_pulse_N_2000")
ion = run.load_ion()
hedft = load_hedft_trajectory("data/reference/9A_All_Data.csv")
distance_result = compare_distance(ion, hedft)
print(distance_result.rmse, distance_result.mean_ratio)
```

## Scope decisions agreed with user

In scope:

- single-pulse neutral and ion dynamics,
- 9 Å HeDFT comparison,
- normalized 18 Å HeDFT reference loading and smoke coverage,
- first VMI reference-data loading and final-velocity histogram helpers,
- MATLAB/Python reference validation.

Out of scope:

- pump-probe,
- effusive dynamics,
- full experimental VMI analysis beyond the current reference overlays,
- 18 Å production interpretation beyond normalized reference loading,
- Abel inversion,
- image-processing utilities.
