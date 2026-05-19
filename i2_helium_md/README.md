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
| 13 | `postprocess/hedft_loader.py`, `postprocess/compare_trajectories.py`, `postprocess/velocity_distribution.py`, `postprocess/energy_balance.py`, `postprocess/polar_velocity.py`, `postprocess/velocity_2d.py`, `postprocess/pair_correlation.py`, `postprocess/time_resolved.py`, `postprocess/boltzmann_overlay.py`, `scripts/post_processing/plot_run_summary.py` | `simulation_image_only_trajectories.m`, parts of `simulation_image.m`, legacy live-debug and in-scope paper post-processing scripts | done; consolidated post-processing summary and helper APIs present |
| 14 | `postprocess/paper_v2.py`, `postprocess/paper_v2_plotting.py`, `postprocess/paper_v3.py`, `postprocess/paper_v4.py`, `postprocess/paper_cov.py`, `postprocess/paper_cov_plotting.py`, `postprocess/_smoothing.py`, `scripts/post_processing/plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`, `plot_paper_cov.py` + reference exports under `data/reference/paper_v2/`, `paper_v3/`, `paper_v4/`, `paper_cov/` | `post_process_single_pulse_paper_IplusHe_comparison.m`, `post_process_single_pulse_paper_v3.m`, `post_process_single_pulse_paper_v4.m`, `post_process_single_pulse_paper_IplusHe_comparison_cov.m` | done; active droplet branches fully ported, experimental references frozen |

## Current phase

The MATLAB → Python port is complete:

- neutral and ion propagation drivers,
- ion-stage MATLAB/Python cross-reference validation,
- public single-pulse run script writing `cfg.json`, `neutral.npz`,
  and `ion.npz` via `RunDirectory`,
- full in-scope post-processing surface: focused legacy-debug and
  paper scripts (`plot_neutral_energy_balance.py`,
  `plot_ion_energy_balance.py`, `plot_ion_temperature_diagnostic.py`,
  `plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`,
  `plot_paper_cov.py`), the consolidated `plot_run_summary.py`
  driver, and the helpers under `i2_helium_md/postprocess/`,
- experimental reference exports under `data/reference/paper_v2/`,
  `paper_v3/`, `paper_v4/`, `paper_cov/`, and `vmi_summary/`,
- visual comparison of the generated Python figures against the legacy
  MATLAB figures, with any recipe mismatches literal-ported.

The next phase is the **drag-model physics port**: deprecate the
hard-sphere collision model in `physics/collisions.py` and replace it
with a TDDFT-calibrated drag-force model for I⁺ in the helium bubble.
The post-processing surface is the comparison layer for the new
physics and stays stable.

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
- `docs/post_process/hedft_loader_module.md` — walkthrough of normalized HeDFT reference loading
- `docs/post_process/compare_trajectories_module.md` — walkthrough of numerical MD/HeDFT trajectory comparison
- `docs/post_process/velocity_distribution_module.md` — walkthrough of VMI reference and final-velocity histogram helpers
- `docs/post_process/energy_balance_module.md` — walkthrough of neutral / ion energy-balance helpers
- `docs/post_process/plot_hedft_comparison_script.md` — usage guide for the focused HeDFT comparison plotting script
- `docs/post_process/plot_legacy_debug_scripts.md` — usage guide for the focused legacy-debug plotting scripts
- `docs/post_process/post_processing_strategy.md` — overall strategy notes for the post-processing port
- `docs/post_process/post_processing_port_plan.md` — final inventory for the post-processing port
- `docs/post_process/scripts/plot_paper_v2.md` — walkthrough of `plot_paper_v2.py`
- `docs/post_process/scripts/plot_paper_v3.md` — walkthrough of `plot_paper_v3.py`
- `docs/post_process/scripts/plot_paper_v4.md` — walkthrough of `plot_paper_v4.py`
- `docs/post_process/scripts/plot_paper_cov.md` — walkthrough of `plot_paper_cov.py`
- `migration_log.md` — chronological record of decisions, deviations, and open questions
- `current_state.md` — completed modules and current validation phase
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
│   └── postprocess/             HeDFT, VMI, energy, and summary diagnostics
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
| `vmi_summary/vmi_iplus_gas.csv` | exported experimental gas-phase VMI reference |
| `vmi_summary/vmi_iplus_he.csv` | exported experimental droplet VMI reference |
| `vmi_summary/vmi_iplus_he_high_snr.csv` | exported high-SNR experimental droplet VMI reference |
| `paper_v2/*` | radial / phi / 2-D image references for `plot_paper_v2.py` (high-SNR + 160/300/600 mW droplet curves, processed VMI image MAT + JSON sidecars) |
| `paper_v3/*` | radial / phi references for `plot_paper_v3.py` (high-SNR + timescan curves) |
| `paper_v4/*` | radial references for `plot_paper_v4.py` (droplet + gas variants at multiple powers) |
| `paper_cov/iplus_he_covariance.mat` (+ JSON sidecar) | experimental pair-covariance matrices for `plot_paper_cov.py` |
| `paper_cov/iplus_he_phi.csv` | experimental phi distribution for `plot_paper_cov.py` |

The HeDFT loader reads the common columns
`Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`. The earlier split
legacy 9 A files (`data_vabs2.csv` and `R1-R2.csv`) remain useful provenance,
but the Python API now uses the normalized all-data files. The per-paper
exporters under `data/reference/scripts/export_paper_*.m` are the MATLAB
sources of truth for the corresponding subdirectory contents.

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

Build the consolidated post-processing PDF for the 9 A HeDFT comparison run:

```bash
python scripts/post_processing/plot_run_summary.py data/runs/9A_hedft_comparison --hedft-ref data/reference/9A_All_Data.csv --no-show
```

Build the consolidated PDF for the experimental-condition droplet run:

```bash
python scripts/post_processing/plot_run_summary.py data/runs/single_pulse_droplet --vmi-ref-he data/reference/vmi_iplus_he.csv --vmi-ref-gas data/reference/vmi_iplus_gas.csv --no-show
```

The focused paper-figure scripts are configured through `# USER SETTINGS`
blocks at the top of each file and are launched with no CLI arguments —
edit and run:

```bash
python scripts/post_processing/plot_paper_v2.py
python scripts/post_processing/plot_paper_v3.py
python scripts/post_processing/plot_paper_v4.py
python scripts/post_processing/plot_paper_cov.py
```

Or use the numerical comparison API directly:

```python
from i2_helium_md.postprocess import compare_distance, load_hedft_trajectory
from i2_helium_md.simulation.run_directory import RunDirectory

run = RunDirectory("data/runs/9A_hedft_comparison")
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
- VMI reference-data loading and final-velocity histogram helpers,
- consolidated in-scope post-processing diagnostics,
- MATLAB/Python reference validation,
- authentic post-processing reproduction of legacy figures where reference
  data and run outputs are available,
- drag-model physics calibrated against TDDFT for I⁺ in the helium
  bubble (replacing the hard-sphere collision model).

Out of scope:

- pump-probe,
- effusive dynamics,
- full experimental VMI analysis beyond the current reference overlays,
- 18 Å production interpretation beyond normalized reference loading,
- Abel inversion,
- image-processing utilities.
