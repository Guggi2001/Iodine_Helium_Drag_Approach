# The `leapfrog.py` module — a walkthrough

## What problem does this file solve?

Given the forces on every atom (from `interactions.py`) and the droplet
potential (from `potentials.py`), we need a rule for **moving atoms forward
in time**. This module implements that rule.

## Position in the dependency chain

```
config.py
   ↓
physics/constants.py      (universal constants)
   ↓
physics/potentials.py     (U(r): Morse, Coulomb, droplet)
   ↓
physics/interactions.py   (pair forces: partner_interaction_*)
   ↓
physics/leapfrog.py       ← THIS MODULE
   ↓
simulation/neutral.py     (full simulation driver — next)
```

The file answers the question: **"Given positions, velocities, and a
force model, where are the atoms one timestep later?"**

## The algorithm: velocity-Verlet (not classical leapfrog)

Despite the legacy file being named `frog_step_neutral.m`, the algorithm is
actually **velocity-Verlet** in its "kick-drift-kick" form:

```
# 1. Compute acceleration at current position
a0 = accel(x0)

# 2. Drift: update position using current velocity and acceleration
x1 = x0 + v0*dt + 0.5*a0*dt²

# 3. Compute acceleration at new position
a1 = accel(x1)

# 4. Kick: update velocity with the average acceleration
v1 = v0 + 0.5*(a0 + a1)*dt
```

### Why this algorithm?

Velocity-Verlet is the **gold standard** for molecular dynamics because:

1. **Symplectic** — it preserves the symplectic structure of Hamiltonian
   mechanics, which means energy drift is bounded over time rather than
   growing unboundedly (as would happen with naive Euler integration).

2. **Second-order accurate** — error in positions scales as O(dt²), so
   halving dt quarters the error.

3. **Time-reversible** — running the simulation backward recovers the
   original trajectory (up to floating-point precision). This is physics:
   classical mechanics is time-reversible, and a good integrator should
   be too.

4. **Only one force evaluation per step** in the "one-pass" variant (we
   reuse `a1` of this step as `a0` of the next — though our current
   implementation recomputes it for simplicity).

### What it's **not**

Classical leapfrog stores velocity at half-timesteps (`v(t+dt/2)`). Our
algorithm keeps position and velocity co-located at integer timesteps,
which is what "velocity-Verlet" means. The two algorithms are
mathematically equivalent for conservative forces, but velocity-Verlet is
easier to start, stop, and checkpoint.

## Architectural design: the "pluggable acceleration" pattern

The biggest refactor from MATLAB is the **separation of algorithm from
physics**. The MATLAB `frog_step_neutral.m` was a monolithic
~150-line function that interleaved:

- The velocity-Verlet math
- Droplet force computation
- Partner interaction computation
- Optional velocity scattering
- Optional ionization branching

In Python we split these into layers:

```
velocity_verlet_step(pos, vel, acc_fn, dt)
    │
    └─ calls acc_fn(pos)  which is one of:
        ├─ neutral:  _neutral_accel_fn  →  droplet + partner_interaction_neutral
        └─ ion:      _ion_accel_fn      →  droplet + partner_interaction_ion
```

The integrator itself is **physics-agnostic**. It takes any function of the
form `positions → (accelerations, potential_energy)` and steps forward.
This means:

- Testing the integrator becomes trivial (feed it constant or zero forces
  and check the algebra — see `TestVelocityVerlet`).
- Adding a new force model means writing a new `_accel_fn` — no changes to
  the integrator.
- If we ever want to swap velocity-Verlet for a higher-order integrator
  (RK4, Gear predictor-corrector), the force assemblers stay the same.

## Public API

```python
from i2_helium_md.physics.leapfrog import (
    velocity_verlet_step,  # pure algorithm, takes an acc_fn
    make_neutral_step,     # returns a step function for neutral propagation
    make_ion_step,         # returns a step function for ion propagation
)
```

The **factory functions** (`make_*_step`) are what simulation drivers will
use. Each one binds the context (mass, droplet radii, charge, etc.) and
returns a closure that the main loop calls with just
`(positions, velocities, dt)`.

### Usage pattern

```python
from i2_helium_md import single_pulse_N2000
from i2_helium_md.physics.leapfrog import make_neutral_step

cfg = single_pulse_N2000()

# Setup (done once, before the main loop)
step = make_neutral_step(cfg, mass, droplet_radii)

# Main loop
pos = (x0, y0, z0)
vel = (vx0, vy0, vz0)
for t in range(num_timesteps):
    pos, vel, E_pot = step(pos, vel, dt=cfg.dt_neutral)
    # ... record trajectory, check for escape conditions, etc.
```

For ion propagation, identical pattern with `make_ion_step`:

```python
step = make_ion_step(cfg, mass, droplet_radii, charge, state_ids=None)
pos, vel, E_pot = step(pos, vel, dt=cfg.dt_ion)
```

## Internal walkthrough

### `velocity_verlet_step(pos, vel, acc_fn, dt)`

The pure algorithm. Takes tuples of positions and velocities (each a 3-tuple
of 2N arrays), a callable that maps positions to accelerations, and a
timestep. Returns updated positions, velocities, and the potential energy
evaluated at the new positions.

**Why return E_pot?** For diagnostics. After each step we want to track
energy conservation (KE + PE) to verify the integrator is behaving well.

**Why take `acc_fn` as an argument?** So unit tests can pass in analytical
force fields (zero, constant, harmonic) and verify the algorithm
independently of the physics.

### `_droplet_acceleration(x, y, z, mass, droplet_radii, cfg, use_ion_binding)`

Computes the force from the droplet solvation potential and converts it to
acceleration. This is the same computation for both neutral and ion runs
— only the binding energy differs.

Force chain (conceptual):
```
   r      = |(x, y, z)|                                             (Angstrom)
   depth  = r - droplet_radii                                       (Angstrom)
   dU/dr  = droplet_force(depth, steepness, binding_energy)         (eV/A)
   F      = -dU/dr  * r_hat                                         (eV/A, radial)
   F[N]   = F * EV / 1e-10                                          (Newtons)
   a[m/s²]= F[N] / mass                                             (m/s²)
   a[A/ps²] = a[m/s²] * 1e10 * (1e-12)^2                            (A/ps²)
```

In practice the three unit conversions collapse into a single multiplication
by the shared constant `EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 = EV * 1e-4 ≈ 1.602e-23`:

```python
a_mag = -dU_dr / mass * EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2   # A/ps^2
```

The result is projected onto the radial unit vector and returned as 3
Cartesian components.

### `_neutral_accel_fn(pos, ctx, cfg)`

Sums the neutral acceleration contributions:

1. Droplet force (always on, uses I atom binding energy).
2. I–I partner Morse force (only if `cfg.partner_interaction=True`).

Returns `((ax, ay, az), E_pot_per_pair)` in the format the integrator
expects.

### `_ion_accel_fn(pos, ctx, cfg)`

Sums the ion acceleration contributions:

1. Droplet force (always on, uses **ion** binding energy).
2. I⁺–I⁺ partner force — either pure Coulomb or Coulomb+Morse depending
   on config.

Converts the per-atom `E_pot` from `partner_interaction_ion` back to
per-pair form (length N) so neutral and ion interfaces match.

### `_StepContext` dataclass

A small bundle of per-simulation inputs that the acceleration functions
need but the integrator doesn't care about: mass, droplet_radii, charge,
state_ids. Packaged together so the factory functions can bind them once
and hide them from the integrator.

## Unit flow

This is the single most bug-prone area of MD simulations, so call it out
explicitly:

| Quantity   | MATLAB/Python | Units       |
|------------|---------------|-------------|
| Position   | x, y, z       | Angstrom    |
| Velocity   | vx, vy, vz    | Angstrom/ps |
| Acceleration | ax, ay, az  | Angstrom/ps² |
| Mass       | mass          | kg          |
| Time       | dt            | picosecond  |
| Force      | dU/dr         | eV/Angstrom |
| Potential  | U             | eV          |

**Conversion factor used everywhere:**
```
a[A/ps²] = F[eV/A] / mass[kg] * EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2
                                = 1.602e-23
```

Defined once in `physics/constants.py` and imported by both `leapfrog.py`
(droplet acceleration) and `interactions.py` (pair force acceleration).
Previously these two files had two slightly different numerical
implementations of the same conversion — see `migration_log.md` for the
audit.

## Features intentionally skipped

The MATLAB ion integrator had two extra features that are out of scope:

1. **`he_direction_scattering`**: random perturbation of velocity
   direction each step. Was set to 0 in the input files we target, so
   it never activated. Skipped.

2. **`additional_droplet_charges`**: models a charged droplet background
   via `add_helium_interaction_coulomb` or `add_helium_interaction_charged_droplet`.
   Used for the pump-probe scenarios which are out of scope. Skipped.

If these become relevant later, they're just extra terms added to
`_ion_accel_fn`, and the integrator won't need to change.

## Testing strategy

The test file `tests/test_leapfrog.py` is organized into three layers,
each testing one abstraction level:

### 1. `TestVelocityVerlet` — the pure algorithm

- Zero force → position advances as x + v·dt (trivial kinematics).
- Constant force → matches analytical `x = v·dt + 0.5·a·dt²` exactly.

These tests never involve any physics; they only verify the math of
velocity-Verlet is correctly implemented.

### 2. `TestNeutralStep` / `TestIonStep` — force assembly

- Atoms at equilibrium barely move.
- Compressed molecules dissociate (atom 1 gains +vx, atom 2 gains -vx).
- Ions with same charge Coulomb-explode.
- Missing charge raises `ValueError` early.

### 3. `TestDropletForce` — the droplet contribution

- A stationary atom at the droplet surface gets pulled inward.

### The energy conservation test

The crown jewel of the suite: simulate a vibrating I₂ molecule for 200 fs
and check that total energy (KE + PE) drifts less than 1%.

This is the single test that would catch:

- Wrong sign on a force.
- Wrong unit conversion.
- Integrator bugs (using `v0` instead of `v1` anywhere).
- Subtle asymmetries in how atom 1 and atom 2 are updated.

Current result: **~0.18% drift over 200 steps** — well within the 1%
tolerance. Most of this drift is the finite-difference truncation in the
Morse force (see `interactions.py` docs); switching to analytical
derivatives would drop it further.

## Common pitfalls

1. **Don't pass `mass` in amu.** It's in kg. Construct via
   `np.full(2*N, MASS_I_AMU * U)` where `U` is the amu-to-kg factor from
   `constants.py`.

2. **The `acc_fn` must return a pair** `(accelerations_tuple, E_pot_array)`.
   Returning just the accelerations will crash the integrator.

3. **`E_pot` returned by `_ion_accel_fn` is per-pair (length N)**, converted
   from `partner_interaction_ion`'s per-atom output. Don't accidentally sum
   it over both atoms.

4. **Timestep too large** → energy conservation fails. The Morse vibration
   period is ~160 fs, so dt should be ≤ 2–5 fs for clean oscillations.
   `cfg.dt_neutral = 0.01 ps = 10 fs` is right at the edge; our energy-
   conservation test uses `dt = 0.001 ps = 1 fs` for a cleaner check.

## Future improvements

1. **Analytical partner forces.** Currently Morse forces use finite
   differences (preserving the MATLAB convention). Swapping to
   `dU/dr = 2·a·D_e·(1 − exp(−a(r−R_e)))·exp(−a(r−R_e))` (with the Xdip
   derivative added when active) would halve runtime and eliminate FD
   truncation error. Plan to do this once the full pipeline runs end-to-end.

2. **Adaptive timestep.** Short-range ionic Coulomb forces can require
   dt ~ 10⁻³ ps when atoms are close; out in free space dt ~ 10⁻¹ ps is
   fine. An adaptive dt strategy could speed long runs substantially.

3. **Vectorized multi-timestep storage.** The main simulation loop will
   want to call `step` ~20000 times and store the results. Pre-allocating
   a (2N, num_timesteps) array and writing into slices is faster than
   appending to Python lists.
