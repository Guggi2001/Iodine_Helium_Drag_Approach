# The `interactions.py` module — a walkthrough

## What problem does this file solve?

The simulation's integrator (next step) will ask, at every timestep, a simple
question:

> "Given where every atom is right now, what is the acceleration on each atom?"

The potentials module (`potentials.py`) gives us the *energy* U(r) for a pair
of atoms at distance r. But an integrator needs *acceleration*, not energy.
Going from one to the other requires several mechanical steps:

1. Compute the pair distance r from 3D coordinates.
2. Differentiate U(r) to get a scalar force F.
3. Apply the force along the pair axis (unit vector).
4. Split F into Cartesian (x, y, z) components.
5. Divide by mass and apply correct unit conversions.
6. Apply Newton's third law: atom 1 and atom 2 get opposite accelerations.

All of this is mechanical bookkeeping — not physics. Keeping it in the
potentials file would clutter clean math. Putting it in the integrator would
tangle the integrator with geometry. So it lives in its own module:
`interactions.py`.

## Position in the dependency chain

```
config.py
   ↓
physics/constants.py      (universal constants: eV, u, k_B, ...)
   ↓
physics/potentials.py     (U(r): Morse, Coulomb, droplet)
   ↓
physics/interactions.py   ← THIS MODULE
   ↓
physics/leapfrog.py       (time stepping)
   ↓
simulation/neutral.py     (full simulation)
```

Each layer answers exactly one question. This module answers:
**"Given 3D atom positions, what are the 3D accelerations?"**

## The "2N layout" convention

One convention you must internalize to read this file:

> **Coordinate arrays in this simulation have length 2N for N molecules.**
>
> - Indices `0 .. N-1` hold atom 1 (the "first" atom) of each molecule.
> - Indices `N .. 2N-1` hold atom 2 (the "twin") of each molecule.

This mirrors the original MATLAB code, which used the same scheme to
vectorize pair operations. If N = 5 molecules, the coordinate array `x`
looks like:

```
index:  [  0    1    2    3    4     5    6    7    8    9  ]
         ┌──────────────────────┐ ┌──────────────────────────┐
         │  atom-1 of molecule  │ │  atom-2 of molecule      │
         │  0   1   2   3   4   │ │  0   1   2   3   4       │
         └──────────────────────┘ └──────────────────────────┘
```

Molecule 3's two atoms are at indices 3 and 8. Their positions are
`(x[3], y[3], z[3])` and `(x[8], y[8], z[8])`.

**Why this layout?** It lets you compute all pair distances with one
vectorized subtraction:

```python
dr_vec = np.stack([x[:N] - x[N:], y[:N] - y[N:], z[:N] - z[N:]], axis=-1)
```

No Python loop, no per-molecule overhead.

**Where does this convention live?**
`interactions.py` is the *only* module that exposes this layout. Everything
above (configs, potentials) works with scalar pair quantities. Everything
below (leapfrog, simulation) receives clean per-atom arrays. This is
encapsulation: the messy convention is isolated behind a clean interface.

## Public API (what other modules will use)

```python
from i2_helium_md.physics.interactions import (
    atom_interaction_potential,    # scalar U for neutral I-I pairs
    ion_interaction_potential,     # scalar U for ion pairs (Coulomb ± Morse)
    partner_interaction_neutral,   # a_x, a_y, a_z, E_pot for neutral pairs
    partner_interaction_ion,       # a_x, a_y, a_z, E_pot for ion pairs
)
```

**The two force functions (`partner_interaction_neutral` / `_ion`) are what
the integrator calls every timestep.** The two potential functions are used
internally and also available as stand-alone numerical tools.

## Function-by-function walkthrough

### `atom_interaction_potential(dr, cfg)`

A thin wrapper around `morse_X`. The MATLAB code had a `atom_interaction_potential.m`
file separate from `get_morse_potential_X.m`, and both were used — so we
preserve the name for readability and to mirror the legacy structure.

**Input:** array of pair distances, `SimConfig`
**Output:** array of neutral I-I potential energies in eV

### `ion_interaction_potential(dr, q1, q2, cfg, *, state_ids=None)`

Implements the full ion interaction logic from `ion_interaction_potential.m`:

- **Pure Coulomb term**: `E_coulomb_scale · q1 · q2 · 14.4 / r` — always included.
- **Morse I2+ term**: added only when `cfg.single_charge_ionization_allowed`
  is True AND the pair is asymmetrically ionized (q1 + q2 == 1).

The `E_coulomb_scale` knob lets the author empirically tune the Coulomb
strength without changing charges. This reflects partial screening by the
helium environment and is a known free parameter in the model.

See `docs/physics_background.md` §2 for *why* the usual double-ionization
case uses pure Coulomb and not Morse — the short answer is length-scale
and energy-scale arguments.

### `partner_interaction_neutral(x, y, z, mass, cfg)`

This is what the **neutral** integrator calls every timestep.

**Inputs (all length 2N):**
- `x, y, z` — atom positions in Angstrom
- `mass` — atom masses in kg
- `cfg` — simulation config

**Outputs:**
- `ax, ay, az` — acceleration in Å/ps² (length 2N)
- `E_pot` — potential energy per pair in eV (length N)

**Internal flow:**
```
(x, y, z)                       length 2N
    ↓   _split_pair_coordinates
(r1, r2)                        two arrays length N
    ↓   _pair_geometry
(dr, dr_unit)                   length N each
    ↓   atom_interaction_potential  (uses morse_X)
E_pot(dr)                       length N
    ↓   _force_from_potential_fd
F(dr)                           length N, eV/Å
    ↓   _acceleration_from_force  (Newton's 3rd law + unit conversion)
(ax, ay, az)                    length 2N, Å/ps²
```

### `partner_interaction_ion(x, y, z, mass, charge, cfg, *, state_ids=None)`

Same pipeline as `partner_interaction_neutral`, except:

1. Extracts per-molecule charges from the 2N-long `charge` array.
2. Calls `ion_interaction_potential` instead of `atom_interaction_potential`.
3. Follows MATLAB's convention of returning `E_pot` **per atom** (length 2N)
   rather than per pair (length N). Each atom receives half the pair
   energy, so `sum(E_pot)` equals the total pair energy. This matters for
   energy-conservation diagnostics.

## Private helpers (read if you want to understand the internals)

These are implementation details. You don't call them directly, but
reading them clarifies what `partner_interaction_*` actually does.

### `_split_pair_coordinates(x, y, z)`

Splits a 2N array into two N arrays. Raises `ValueError` if the input has
odd length (a defensive check that would silently produce garbage in MATLAB).

### `_pair_geometry(x, y, z)`

Given 2N coordinate arrays, returns:
- `dr` (length N) — scalar distance between the two atoms of each molecule
- `dr_unit` (shape (N, 3)) — unit vector pointing from atom 2 toward atom 1

This is used once at the start of every force computation. The unit vector
encodes the direction along which the pair force acts.

### `_force_from_potential_fd(potential_fn, dr, h=1e-4)`

The finite-difference force:

```
F(r) = (U(r) - U(r+h)) / h   ≈   -dU/dr
```

**Why this sign convention?** It gives a positive F when the potential
decreases with increasing r — i.e. when the force wants to push r larger.
For a repulsive short-range potential (like the Morse at r < R_e or
Coulomb), this means F > 0 and the atoms fly apart. For attractive regions
(Morse at r just > R_e), F < 0 and the atoms are pulled together.

We kept the MATLAB formula verbatim for byte-for-byte compatibility.

**Limitation:** Finite-difference forces are first-order accurate and have
truncation error proportional to h × U''(r). At a minimum of U (where
U' = 0), the FD method returns a tiny-but-nonzero force — `~0.04 Å/ps²`
for I2 X state at R_e. This isn't a bug, it's a known artifact of FD. See
`docs/migration_log.md` for the discussion.

A future improvement is to replace FD with analytical derivatives (like
we already did for `droplet_force`). This is worth doing if energy
conservation becomes a concern.

### `_acceleration_from_force(F, dr_unit, mass)`

Converts a scalar pair force into per-atom 3D accelerations.

The physics:
- Atom 1 gets force `+F · dr_unit` (pushed in the direction of dr_unit).
- Atom 2 gets force `−F · dr_unit` (Newton's 3rd law).
- Acceleration = force / mass, with unit conversion `eV/(Å·u) → Å/ps²`.

The unit conversion factor `9648.533` was derived in a comment at the top
of the module. It comes from:

```
1 eV/Å / 1 u = 1.602e-19 J/m / 1.66054e-27 kg = 9.648e17 m/s²
             = 9.648e17 × 1e10 Å/m × (1e-12 s/ps)²
             = 9648.5 Å/ps²
```

## Example usage

```python
import numpy as np
from i2_helium_md import single_pulse_N2000
from i2_helium_md.physics.constants import MASS_I_AMU, U
from i2_helium_md.physics.interactions import partner_interaction_neutral

# Three molecules, all aligned along x at equilibrium distance
N = 3
R_e = 2.666
x = np.concatenate([np.full(N, +R_e / 2), np.full(N, -R_e / 2)])
y = np.zeros(2 * N)
z = np.zeros(2 * N)
mass = np.full(2 * N, MASS_I_AMU * U)

cfg = single_pulse_N2000()
ax, ay, az, E_pot = partner_interaction_neutral(x, y, z, mass, cfg)

# Shapes
assert ax.shape == (6,)       # 2N per-atom accelerations
assert E_pot.shape == (3,)    # N per-pair energies

# At equilibrium, force ~ 0 (up to FD residual)
assert abs(ax[0]) < 1.0       # Å/ps²

# Newton's 3rd law
assert np.allclose(ax[:N], -ax[N:])
```

## Common pitfalls for future contributors

1. **Don't assume length N — coordinate arrays are 2N.** The docstring of
   every function states the layout.

2. **The `mass` argument is in kg**, not atomic mass units. Multiply by
   `U` (from `constants.py`) when constructing.

3. **The `charge` argument is a 2N array of per-atom charges**, not a list
   of (q1, q2) tuples. The function splits it internally.

4. **`state_ids` is a keyword-only argument** (`*,`). You can't pass it
   positionally — this is enforced to prevent accidental misuse when
   `single_charge_ionization_allowed` is False.

5. **FD force at a minimum is not exactly zero.** See the note under
   `_force_from_potential_fd`.

## Testing

See `tests/test_interactions.py`. Key tests:

- Shape checks (output dimensions match input).
- Newton's 3rd law (atom 1 and 2 accelerations sum to zero for equal masses).
- Force direction (repulsive below R_e, attractive above).
- `E_coulomb_scale` linearly scales the ion force.
- Invalid input shapes raise `ValueError` rather than silently broadcasting.
- Single-charge mode without `state_ids` raises rather than silently
  returning wrong numbers.
