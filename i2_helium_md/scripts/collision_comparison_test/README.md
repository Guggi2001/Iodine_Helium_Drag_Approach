# Hard-sphere collision: Python vs MATLAB statistical comparison

This directory contains a statistical cross-check of the hard-sphere
collision implementation (Mode 3) between our Python port and the legacy
MATLAB code.

## Why a statistical test?

Unlike the leapfrog test in `compare_with_matlab/` (which is fully
deterministic), the collision code uses random numbers — and Python's
`numpy.random.default_rng` produces a fundamentally different sequence
from MATLAB's `rand`. Trajectories will diverge atom-by-atom even if
the physics is identical.

What we **can** test is whether the **distributions** match: collision
rate, energy-loss-per-event statistics, ensemble averages over many
atoms, etc.

## Test setup

- **1000 non-interacting atoms** placed on a deterministic grid inside a
  single He droplet of radius 27.97 Å (well clear of the surface).
- **No Morse, no droplet potential**: atoms free-stream between
  collisions. This isolates the collision sampler and the elastic-
  scattering kinematics from any other physics.
- Each atom starts with **|v| = 5.0 Å/ps** (~0.4 eV) in a deterministic
  direction sampled from a Halton sequence (so directions are
  reproducible across implementations, no RNG involved in init).
- **200 timesteps × 0.01 ps** = 2 ps of simulation.
- Both implementations use the same physical parameters: σ = 30 Å²,
  ρ_droplet = 0.8 × 0.0219 = 0.01752 atoms/Å³, m_He = 4 amu, no Gaussian
  scattering-angle smearing.

This setup yields ~5 collisions per atom on average → ~5000 collision
events per run, plenty of statistics for first/second moments.

## Files

- `generate_init_state.py` — writes `init_state.csv` (1000 atoms,
  positions + velocities, no RNG involved)
- `run_python_collisions.py` — Python driver; writes
  `python_summary.csv` and `python_collision_events.csv`
- `run_matlab_collisions.m` — MATLAB driver; writes
  `matlab_summary.csv` and `matlab_collision_events.csv`
- `compare_collisions.py` — performs the statistical comparison

## How to run

### 1. Generate the initial state

```bash
python scripts/compare_collisions_with_matlab/generate_init_state.py
```

This writes `init_state.csv` next to the script. It's deterministic —
running it again produces the same file.

### 2. Run the Python side

```bash
python scripts/compare_collisions_with_matlab/run_python_collisions.py
```

Writes `python_summary.csv` (per-step ensemble statistics) and
`python_collision_events.csv` (one row per collision).

### 3. Run the MATLAB side

Edit `MATLAB_LEGACY_DIR` at the top of `run_matlab_collisions.m` to
point to your copy of the legacy `Iodine_Helium_Simulation/`
directory. Then:

```bash
matlab -batch "run_matlab_collisions"
```

Writes `matlab_summary.csv` and `matlab_collision_events.csv`.

### 4. Compare

```bash
python scripts/compare_collisions_with_matlab/compare_collisions.py
```

The script reports:
1. **Total collision count** (Poisson noise σ ~ 70 for ~5000 events)
2. **Time-resolved `<E_kin>(t)` and `<E_dissip>(t)`** at sample times
3. **Per-collision ΔE/E₀ distribution**: mean, min, max, χ² test
4. **Theoretical predictions** for elastic scattering with ρ = m_I/m_He = 31.75:
   - `⟨ΔE/E₀⟩ = 2ρ/(1+ρ)² ≈ 5.92%`
   - `max ΔE/E₀ = 1 - ((1-ρ)/(1+ρ))² ≈ 11.84%`

Each test produces a "deviation in σ" number; passing means < 3σ
(consistent with statistical noise).

## What "agreement" means here

For ~5000 events, statistical noise is ~1% on the means. So we expect:
- Collision count agreement to a few percent (Poisson)
- Mean energy agreement to ~1% (sample size)
- Distribution shape agreement (χ² per bin ~ 1)

If any test exceeds 3σ deviation, that **strongly** suggests a real
difference between the implementations — not just noise.

## A note on the MATLAB script

The legacy collision code is embedded inside the ~900-line
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m`, deeply intertwined with the
rest of the simulation. To make this a focused unit test, I
**extracted just the Mode-3 collision logic** into the standalone
`run_matlab_collisions.m`. I did this by copying the relevant lines
(~600-820 of the original) verbatim where possible.

This means **the MATLAB script in this folder is my own assembly**
from legacy parts, not a turnkey legacy script. If the comparison
shows large statistical deviations, check first whether I extracted
the collision kernel faithfully — particularly:

- Order of `rand()` calls (the legacy draws random numbers for *every*
  atom every step, not just the colliders, and we mirror that)
- The `b_collision = trial < p_scatter & depth < 0` condition
- The `cos(Θ)` ↔ lab-frame transformation
- The Landau cutoff

The line-for-line annotations in `run_matlab_collisions.m` reference
the corresponding MATLAB legacy line numbers.

## What this DOESN'T test

- The full `propagation_step` integration with collisions enabled
  (this isolates the collision sampler; integration is tested
  separately in `compare_with_matlab/`)
- Mode 1 and Mode 2 collision modes (Python explicitly raises for
  these; they were never used in production)
- Velocity-dependent cross sections (`sigma_dependent_on_v` in MATLAB —
  not used in production)
- The Gaussian angle smearing (we set it to zero for this test)
