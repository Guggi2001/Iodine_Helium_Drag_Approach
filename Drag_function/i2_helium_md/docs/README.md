# i2_helium_md

Modern Python port of Michael Stadlhofer's MATLAB molecular-dynamics code
for simulating iodine in helium nanodroplets (single-pulse, HeDFT-comparison
scope).

## Migration progress

| # | Python module | MATLAB source | Status |
|---|---------------|---------------|--------|
| 1 | `physics/constants.py` | `physical_constants.m` | ✅ done |
| 2 | `config.py` (SimConfig) | ~36 MATLAB globals | ✅ done |
| 3 | `presets.py` | `inputfiles_dft_comparison/single_pulse_N2000.m` | ✅ done |
| 4 | `physics/potentials.py` | `droplet_potential.m`, `get_morse_potential_X.m`, `get_morse_potential_I2plus.m`, `morse_potential_I2plus_state_select.m` | ✅ done |
| 5 | `physics/interactions.py` | `atom_interaction_potential.m`, `ion_interaction_potential.m`, `add_partner_interaction.m`, `add_partner_interaction_ion.m` | ✅ done |
| 6 | `physics/leapfrog.py` | `frog_step_neutral.m`, `frog_step_ion.m` | ✅ done |
| 7 | `sampling/droplet_sizes.py` | `generate_droplet_sizes.m`, `get_dropletsize.m` | ✅ done |
| 8 | `sampling/radial_positions.py` | `generate_radial_samples_3d.m` | ✅ done |
| 9 | `simulation/checkpoint.py` | `save('neutral_propagation_checkpoint', ...)` | ⏳ |
| 10 | `simulation/neutral.py` | `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` | ⏳ |
| 11 | `simulation/ion.py` | `vmi_sim_3d_ion_propa.m` | ⏳ |
| 12 | `scripts/run_single_pulse.py` | `run_simulation.m` | ⏳ |
| 13 | `postprocess/hedft_loader.py` + `compare_trajectories.py` | `simulation_image_only_trajectories.m` | ⏳ |

## Documentation

- `docs/physics_background.md` — physical model, potentials, and design rationale
- `docs/constants_module.md` — walkthrough of the `constants.py` module
- `docs/interactions_module.md` — walkthrough of the `interactions.py` module
- `docs/leapfrog_module.md` — walkthrough of the `leapfrog.py` integrator
- `docs/droplet_sizes_module.md` — walkthrough of the `droplet_sizes.py` sampler
- `docs/radial_positions_module.md` — walkthrough of the `radial_positions.py` sampler
- `docs/migration_log.md` — chronological record of decisions, deviations, and open questions

## Project layout

```
i2_helium_md_py/
├── data/reference/              data files copied from legacy repo (see below)
├── docs/                        physics background + migration log
├── i2_helium_md/                Python package
│   ├── config.py                SimConfig dataclass
│   ├── presets.py               preset builders
│   ├── physics/                 constants, potentials, interactions, integrators
│   ├── sampling/                random samplers
│   ├── simulation/              neutral + ion propagation
│   └── postprocess/             HeDFT comparison plots
├── scripts/                     entry points
├── tests/                       smoke tests + pytest suite
└── pyproject.toml
```

## Data files needed in `data/reference/`

Copy these three files from the legacy MATLAB repo:

| Legacy path | → | New path |
|-------------|---|----------|
| `HeDFT_MD_comparison_neutral/custom_start_interpolating_functions.mat` | → | `data/reference/hedft_custom_start.mat` |
| `single_pulse_simulation/HeDFT_comparison/9Angstr*ö*m/data_vabs2.csv` | → | `data/reference/hedft_9A_velocity.csv` |
| `single_pulse_simulation/HeDFT_comparison/9Angstr*ö*m/R1-R2.csv` | → | `data/reference/hedft_9A_distance.csv` |

Once copied they're referenced via `SimConfig.data_dir`, no hardcoded paths
anywhere in the code.

## Quickstart (after all steps done)

```python
from i2_helium_md import single_pulse_N2000
from i2_helium_md.simulation import run_neutral, run_ion

cfg = single_pulse_N2000(num_molecules=500, seed=123)
neutral_result = run_neutral(cfg)
ion_result = run_ion(cfg, neutral_result)
```

## Scope decisions (agreed with user)

* **In scope:** single-pulse neutral + ion dynamics; 9 Å HeDFT comparison
* **Out of scope:** pump-probe, effusive, VMI experimental comparison, 18 Å HeDFT
  (no data in repo), Abel inversion, all image-processing utilities

## Open items

* `g()` function referenced in `get_morse_potential_X.m` is not defined in
  the legacy repo. We implement it as a Gaussian
  `g(sigma, A, r) = A * exp(-(r - sigma)^2 / (2 * width^2))`
  pending user confirmation of the original definition.
