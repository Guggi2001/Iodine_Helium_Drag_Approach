# The `constants.py` module — a walkthrough

## What problem does this file solve?

Every physical constant the simulation needs lives here. The MATLAB code
defined these in a script (`physical_constants.m`) that polluted the global
workspace. The Python version is a clean module of named constants — no
global state, no script execution.

## Position in the dependency chain

```
physics/constants.py        ← THIS MODULE
   ↓
config.py                   (uses K_B, EV, MASS_I_AMU for derived properties)
   ↓
everything else             (physics, sampling, simulation modules)
```

Constants is the lowest level of the codebase — it imports nothing from
within the project, and almost everything else depends on it. Keep it
small, well-documented, and free of any logic.

## What's inside

### 1. Universal physics (port of `physical_constants.m`)

```python
E_CHARGE: float = -1.602e-19          # C
EPSILON_0: float = 8.85418781e-12     # F/m
U: float = 1.66053907e-27             # kg
EV: float = 1.602e-19                 # J
K_B: float = 1.380649e-23             # J/K
HC: float = 1240.0                    # eV*nm
EV_PER_WAVENUMBER: float = 1.0 / 8065.54429
```

Each maps one-to-one to a line in the MATLAB source. The `ALL_CAPS` naming
convention is Python's standard signal for "fixed value, don't reassign."

### 2. Iodine-specific constant

```python
MASS_I_AMU: float = 127.0
```

The MATLAB code wrote `127*u` inline everywhere. Pulling it into a named
constant makes the species explicit at every call site.

### 3. Helium droplet density

```python
BULK_DENSITY_HELIUM: float = 0.0219            # atoms / A^3
DENSITY_DROPLET: float = 0.8 * BULK_DENSITY_HELIUM
```

The 0.8 factor accounts for surface effects and is a standard assumption
in helium-droplet literature (Phys. Rev. B 58, 3341).

`DENSITY_DROPLET` is automatically recomputed if you ever change
`BULK_DENSITY_HELIUM`. A small thing MATLAB couldn't easily express.

### 4. Coulomb helpers

```python
def coulomb_energy(r_angstrom: np.ndarray | float) -> np.ndarray | float
def coulomb_velocity(r_angstrom, mass_kg) -> np.ndarray | float
```

These replace MATLAB's anonymous-function-handle style:
```matlab
coulomb_energy = @(x) e_charge.^2 ./ (4*pi*epsilon_0*x*1E-10);
coulomb_velocity = @(x, m) sqrt(coulomb_energy(x)./m);
```

Used by ion-related code paths. Distances are in Angstrom (typical
simulation units), but the output is in SI Joules — convert with
`/EV` if you need eV.

### 5. The unit-conversion constant (added during a post-Step-7 refactor)

```python
EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2: float = EV * 1e-4    # = 1.602e-23
```

This is the **single source of truth** for converting

> "force in eV/Å acting on a mass in kg → acceleration in Å/ps²"

It exists because both `interactions.py` (pair-force acceleration) and
`leapfrog.py` (droplet-force acceleration) need the same conversion. They
used to compute it independently with two slightly different forms:

- `interactions.py`: `F * (u/mass_kg) * 9648.5` (Å/ps² per eV/(Å·u))
- `leapfrog.py`:     `F * 1.602e-9 / mass_kg * 1e-14` (two-stage conversion)

Both produced numerically equal results, but the rounding in `1.602e-9` vs
`9648.5` made them disagree in the 4th significant digit. Maintaining two
parallel implementations of the same physics is a common source of drift,
so we consolidated them here.

### Derivation
```
Start:  F[eV/A] acting on m[kg]

Step 1: Convert F to Newtons.
        F[N] = F[eV/A] * EV[J/eV] / 1e-10[m/A]

Step 2: Get acceleration in SI.
        a[m/s²] = F[N] / m[kg]

Step 3: Convert acceleration to A/ps².
        a[A/ps²] = a[m/s²] * 1e10[A/m] * (1e-12[s/ps])²
                 = a[m/s²] * 1e-14

Combined:
        a[A/ps²] = F[eV/A] / m[kg] * EV * 1e10 * 1e-14
                 = F[eV/A] / m[kg] * EV * 1e-4
                 = F[eV/A] / m[kg] * 1.602e-23
```

## Departures from MATLAB

### 1. No global namespace pollution
MATLAB's `physical_constants.m` was a *script*, not a function — running
it dumped names like `eV`, `u`, `k_B` into the caller's workspace. That
made it impossible to track where a constant came from when reading code
20 lines down. Python's `from .constants import EV` is unambiguous.

### 2. Anonymous functions become real functions
`coulomb_energy = @(x) ...` becomes `def coulomb_energy(r_angstrom):`. The
benefits: docstrings, type hints, recursion, debugger support, importable.

### 3. Unit suffixes in helper signatures
`coulomb_velocity(r_angstrom, mass_kg)` — the units are in the parameter
names. MATLAB had no equivalent; you had to read the comment block above
the function.

### 4. New constant for unit conversion
`EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2` doesn't exist in the MATLAB code at
all. It emerged when we noticed two parallel implementations in Python and
chose to consolidate.

## Testing

`tests/test_foundations.py` covers:

- The numerical values of the elementary constants (cross-checked vs the
  MATLAB source).
- That `MASS_I_AMU = 127.0` exactly.
- That `binding_energy_I_atom_eV` (a derived `SimConfig` property) matches
  the MATLAB formula `318.43 * k_B / eV`.
- That `EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 ≈ 1.602e-23` and that it
  cross-checks against the alternative `u * 9648.5` form to within
  rounding.

## Style guide for additions

If you find yourself wanting to add a new constant:

1. **Is it a universal physical constant?** → put it here, ALL_CAPS.
2. **Is it specific to one molecule species (e.g. iodine)?** → put it
   here with a species suffix (`MASS_I_AMU`, not `MASS_AMU`).
3. **Is it a tunable simulation parameter?** → it belongs in
   `config.py` as a `SimConfig` field, not here.
4. **Is it a unit-conversion factor used in multiple places?** → put it
   here. Single source of truth beats parallel implementations.
5. **Is it only used inside one module?** → keep it private (`_LOCAL_NAME`)
   in that module's file. Only promote to `constants.py` if a second
   module starts needing it.

## References

- `physical_constants.m` (legacy MATLAB) — primary source for steps 1–4.
- Phys. Rev. B 58, 3341 (1998) — bulk liquid helium density.
- CODATA 2018 — modern values of `K_B`, `e`, `u`, `hc` (current values
  are slightly more precise than the MATLAB ones, but we keep MATLAB's
  values to match legacy results bit-for-bit).
