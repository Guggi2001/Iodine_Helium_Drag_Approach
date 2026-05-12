# The `simulation/neutral.py` module

## What problem does this file solve?

This is the top-level driver for neutral propagation. Given a
`SimConfig`, it runs a complete simulation and produces a
`NeutralCheckpoint`. It orchestrates the three pieces we already built:

1. **`build_initial_state`** — sets up t=0 and pre-allocates the
   trajectory arrays.
2. **`neutral_propagation_step`** (pure) — advances one timestep.
3. **`RunDirectory`** — optional checkpoint serialization.

It also handles two run-time concerns the lower-level pieces don't:

- **Single-pulse mode** runs only 2 internal steps (matching MATLAB
  `t_max = dt*2`).
- **Auto-stride** for memory-bounded checkpoints: if the full-resolution
  trajectory would exceed ~300 MB, only every K-th internal step is
  stored. Internal stepping still happens at `cfg.dt_neutral` so the
  physics (collision rate, leapfrog stability) is unchanged.

## Public API

```python
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.run_directory import RunDirectory

# Simple call: returns a NeutralCheckpoint
ckpt = run_neutral_propagation(cfg, rng=np.random.default_rng(42))

# With saving:
rd = RunDirectory("runs/my-run-1")
ckpt = run_neutral_propagation(cfg, run_dir=rd, verbose=True)

# Custom memory budget (default 300 MB):
ckpt = run_neutral_propagation(cfg, max_bytes=100_000_000)
```

## Single-pulse vs long runs

The `cfg.single_pulse` flag controls how long the run is:

- **`single_pulse=True`** (default for HeDFT comparison snapshots):
  exactly **2 internal steps** are taken. The simulation effectively
  freezes time at the moment of photoexcitation. We retain a 2-column
  checkpoint so downstream plotters and analysis code can treat
  single-pulse and time-resolved runs uniformly.
- **`single_pulse=False`**: the full
  `num_timesteps_neutral = ceil(t_max_neutral / dt_neutral)`
  internal steps are taken.

This matches the MATLAB:

```matlab
if single_pulse
    t_max = dt*2;
end
num_timesteps = ceil(t_max/dt);
```

## Auto-stride for memory budget

If the full-resolution checkpoint would exceed `max_bytes`, the driver
chooses an integer storage stride `K` such that storing every K-th
internal step keeps the checkpoint within the budget:

```
stride = ceil(num_internal_steps / max_storable_steps)
num_stored_steps = ceil(num_internal_steps / stride)
```

The internal loop **always** runs at full resolution (`dt_neutral`).
The collision sampler relies on the previous internal step's
displacement, which would be dramatically wrong if we tried to skip
internal steps. Only **storage** is downsampled.

The state itself is carried as a `NeutralStepState` dataclass between
internal calls; only every K-th state is copied into the checkpoint
arrays.

### Memory estimates

| Run | N | t_max | Internal steps | Full size | Stride | Stored | Saved size |
|---|---|---|---|---|---|---|---|
| Single-pulse | 2000 | 0.02 ps | 2 | 1 MB | 1 | 2 | 1 MB |
| HeDFT 9 Å | 500 | 20 ps | 2000 | 160 MB | 1 | 2000 | 160 MB |
| HeDFT 9 Å | 2000 | 20 ps | 2000 | 640 MB | 3 | 667 | 214 MB |
| Long | 2000 | 200 ps | 20000 | 6400 MB | 22 | 910 | 291 MB |

## DFT pre-fill stub

If `cfg.custom_DFT_start = True`, the driver raises
`NotImplementedError`. The TD-HeDFT pre-fill block in MATLAB
(replacing the first few timesteps with TD-DFT trajectory data) is
not yet ported. When needed, add a function `apply_dft_prefill(state, cfg)`
and call it after `build_initial_state`.

## RunDirectory integration

If `run_dir` is provided:

1. The cfg is saved to `cfg.json` if not already present.
2. The checkpoint is saved to `neutral.npz`.

The driver does NOT re-validate the cfg before saving — that's the
job of `load_neutral_checkpoint(path, cfg=...)` on read.

## What's inside

```
run_neutral_propagation(cfg, *, rng, run_dir, max_bytes, verbose)
├─ _internal_step_count(cfg)                 -- single-pulse override
├─ _decide_stride(N, internal_steps, cap)    -- auto-stride math
├─ build_initial_state(cfg, num_steps, rng)  -- column 0 populated
├─ raise NotImplementedError if cfg.custom_DFT_start
├─ inner loop:
│    state = state_from_checkpoint_column(ckpt, 0)
│    for internal_id in 1..num_internal_steps-1:
│        new = neutral_propagation_step(state, ..., prev_distance, rng)
│        prev_distance = |new - state|
│        if internal_id % stride == 0:
│            write_state_to_checkpoint_column(new, ckpt, next_storage_idx++)
│        state = new
│    if next_storage_idx < num_stored_steps:    -- ensure final state stored
│        write_state_to_checkpoint_column(state, ckpt, next_storage_idx)
├─ run_dir.save_neutral(ckpt)  if run_dir is not None
└─ return ckpt
```

## Tests

- **15 pytest tests** in `test_neutral.py` covering: API contract,
  single-pulse step count, long-run consistency, reproducibility,
  DFT-prefill stub, stride math (small/large cases),
  end-to-end strided run, RunDirectory round trip, energy bookkeeping,
  internal-helper bounds.
- **23 sandbox checks** in `smoke_test_neutral.py` exercising the
  same paths plus a real save/load through `RunDirectory`.

The most informative regression guard is **energy bookkeeping**:
`E_kin + E_pot + E_dissip` should be approximately conserved over a
50-step collision-active run (~3% drift in our test, within 10%
tolerance).

## Departures from MATLAB

1. **Auto-stride is new.** MATLAB stores every step unconditionally;
   that doesn't scale to N=2000 + 200 ps runs. Auto-stride keeps the
   physics identical but caps memory.
2. **`single_pulse` enforced inside the driver, not in cfg.** MATLAB
   sets `t_max = dt*2` if `single_pulse=true`; we keep `t_max_neutral`
   in cfg as the user's *intent* and use `_internal_step_count(cfg)`
   to decide what actually runs.
3. **`custom_DFT_start` stub.** Raises a clear error rather than
   silently doing nothing.
4. **Final state always stored.** When `(internal_steps - 1)` is not
   divisible by `stride`, MATLAB would lose the final state. We
   explicitly write it into the last reachable column so trajectories
   end at the actual end-time.
