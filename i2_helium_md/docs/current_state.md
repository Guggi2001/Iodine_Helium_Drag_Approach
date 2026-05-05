# Current migration state

Last updated: YYYY-MM-DD

## Completed

- `physics/constants.py`
- `config.py`
- `presets.py`
- `physics/potentials.py`
- `physics/interactions.py`
- `physics/leapfrog.py`
- `sampling/droplet_sizes.py`
- `sampling/radial_positions.py`
- `simulation/checkpoint.py`
- `simulation/run_directory.py`
- `sampling/orientations.py`
- `physics/collisions.py`
- `simulation/initial_state.py`
- `simulation/propagation_step.py`
- `simulation/neutral.py`
- Ion steps 11a, 11b, 11c:
  - velocity-dependent cross section helper
  - ion initial-state builder
  - pure ion propagation step

## Currently pending

1. Step 11d: full ion propagation driver in `simulation/ion.py`
2. Step 12: `scripts/run_single_pulse.py`
3. Step 13: HeDFT loading and trajectory comparison in `postprocess/`

## Recommended next task

Implement Step 11d: full ion propagation driver.

The driver should mirror the structure of `simulation/neutral.py`:

- decide internal ion step count
- decide storage stride if needed
- build initial ion state from the neutral checkpoint
- loop with `ion_propagation_step`
- write stored states into `IonCheckpoint`
- optionally save via `RunDirectory`
- add focused tests for shape, final-state storage, mass history, energy bookkeeping, reproducibility, and run-directory roundtrip

## Do not do yet

- Do not implement HeDFT comparison before the full ion driver exists.
- Do not implement plotting before the full simulation path is executable.
- Do not refactor the neutral driver unless a test fails.
- Do not implement out-of-scope MATLAB features.