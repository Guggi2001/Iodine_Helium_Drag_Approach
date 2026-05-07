# How to Run a Single-Pulse Simulation

This project is organized as a pipeline:

1. choose simulation settings,
2. run neutral propagation,
3. run ion propagation,
4. save a self-contained run directory,
5. later load that run for comparison or post-processing.

The user-facing file for steps 1-4 is:

```text
scripts/run_single_pulse.py
```

This script is intentionally not a command-line tool. Open it, edit the user
settings at the top, and run it.

## First Run

From the repository root:

```bash
python scripts/run_single_pulse.py
```

For a first run, use a small smoke configuration:

```text
INPUT_PRESET = "single_pulse_N2000"
RUN_SIZE = "smoke"
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_smoke"
NUM_MOLECULES = 10
SEED = 123
ION_TIME_PS = 0.1
```

This is the right first run. It checks that the whole project works on your
machine without starting a large production simulation.

The script writes:

```text
i2_helium_md/data/runs/single_pulse_smoke/
├── cfg.json
├── neutral.npz
└── ion.npz
```

## What to Edit

Open `scripts/run_single_pulse.py` and edit only the `USER SETTINGS` section
for normal use.

### `INPUT_PRESET`

Choose which migrated MATLAB input file should provide the base configuration:

```python
INPUT_PRESET = "single_pulse_N2000"
INPUT_PRESET = "single_pulse_droplet_distribution"
```

`single_pulse_N2000` mirrors:

```text
legacy_matlab_repository/inputfiles_dft_comparison/single_pulse_N2000.m
```

It uses the fixed 2000-atom droplet setup with `R0_GS = 9 A`.

`single_pulse_droplet_distribution` mirrors:

```text
legacy_matlab_repository/inputfiles_dft_comparison/single_pulse_droplet_distribution.m
```

It uses the droplet-size sampler (`use_single_droplet_size = False`), starts
from `R0_GS = 2.666 A`, and uses `num_molecules = 8000` in production mode.

### `RUN_SIZE`

Choose one of:

```python
RUN_SIZE = "smoke"
RUN_SIZE = "custom"
RUN_SIZE = "production"
```

Use `smoke` for checking that the pipeline works.

Use `custom` when you want to control molecule count, seed, and ion simulation
time yourself.

Use `production` for the selected `INPUT_PRESET` exactly, except for any values
you set in the production override section.

### `RUN_DIR`

This is where the run is saved:

```python
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "my_run_name"
```

Use a new folder name for each run you want to keep. The run directory is the
unit of output in this project. Anchoring the path to `PROJECT_ROOT` makes the
script write to the project-level `data/runs` folder even when your editor
launches the script with `scripts/` as the working directory.

### `OVERWRITE_EXISTING_RUN`

This protects you from accidentally replacing data:

```python
OVERWRITE_EXISTING_RUN = False
```

If `RUN_DIR` already contains `cfg.json`, `neutral.npz`, or `ion.npz`, the
script stops. To intentionally rerun into the same folder:

```python
OVERWRITE_EXISTING_RUN = True
```

### `NUM_MOLECULES`, `SEED`, `ION_TIME_PS`

These control smoke/custom runs:

```python
NUM_MOLECULES = 10
SEED = 123
ION_TIME_PS = 0.1
```

For a quick check, use a small molecule count and short ion time.

For a more meaningful small run, try:

```python
RUN_SIZE = "custom"
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "N100_seed123"
NUM_MOLECULES = 100
SEED = 123
ION_TIME_PS = 20.0
```

### Production Runs

For the full fixed-N HeDFT preset:

```python
INPUT_PRESET = "single_pulse_N2000"
RUN_SIZE = "production"
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_N2000"
```

Production mode uses:

```text
num_molecules = 2000
ion_simulation_time = 20 ps
seed = None
```

For the full droplet-distribution preset:

```python
INPUT_PRESET = "single_pulse_droplet_distribution"
RUN_SIZE = "production"
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet_distribution"
```

Production mode then uses:

```text
num_molecules = 8000
R0_GS_angstrom = 2.666
use_single_droplet_size = False
ion_simulation_time = 20 ps
seed = None
```

If you want a production-like run with one override, use the production
override section:

```python
PRODUCTION_SEED = 123
```

Leave production overrides as `None` to use the preset exactly.

## How to Load Results

After a run finishes:

```python
from i2_helium_md.simulation.run_directory import RunDirectory

run = RunDirectory("data/runs/single_pulse_smoke")
cfg = run.load_cfg()
neutral = run.load_neutral()
ion = run.load_ion()
```

Important checkpoint fields:

- `neutral.positions_x`, `neutral.positions_y`, `neutral.positions_z`
- `neutral.velocities_x`, `neutral.velocities_y`, `neutral.velocities_z`
- `ion.positions_x`, `ion.positions_y`, `ion.positions_z`
- `ion.velocities_x`, `ion.velocities_y`, `ion.velocities_z`
- `ion.positions_final_x`, `ion.positions_final_y`, `ion.positions_final_z`
- `ion.mass_final_kg`
- `ion.b_ion_outside`

The first dimension of trajectory arrays is `2 * num_molecules`: first all atom
1 entries, then all atom 2 entries.

## Effective Project Workflow

Recommended workflow:

1. Run the default smoke script.
2. Confirm `cfg.json`, `neutral.npz`, and `ion.npz` were written.
3. Load the run with `RunDirectory`.
4. Increase `NUM_MOLECULES` or `ION_TIME_PS` gradually.
5. Use a new `RUN_DIR` for every run you want to keep.
6. Only use `OVERWRITE_EXISTING_RUN = True` for disposable test folders.

Avoid starting with production settings. A full run is much larger and slower,
and it is harder to debug if basic environment or path setup is wrong.

## What the Script Does Not Do

The script does not:

- plot trajectories,
- compare to HeDFT,
- compare to experimental VMI images,
- run pump-probe simulations,
- perform Abel inversion.

Those are separate stages. The purpose of this script is to create a clean,
loadable simulation run directory.

## MATLAB Note

The legacy MATLAB `run_simulation.m` runs neutral propagation first, then ion
propagation. This Python script follows the same ordering.

The MATLAB script toggles `sigma_dependent_on_v` off for neutral propagation
and back on for ion propagation. In this Python port, neutral propagation uses
the neutral cross-section setting directly, while ion propagation reads
`cfg.sigma_dependent_on_v`. Therefore no stage-specific config mutation is
needed in `scripts/run_single_pulse.py`.
