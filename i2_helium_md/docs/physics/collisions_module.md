# The `collisions.py` module

## What problem does this file solve?

Inside a helium droplet, an iodine atom (or molecule) traveling through
the liquid undergoes hard-sphere collisions with helium atoms. Each
collision changes the iodine's velocity direction and removes some
kinetic energy. This module provides the physics:

1. **Stochastically decide** which particles collide on a given timestep.
2. **Apply** an elastic 2-body scattering event: sample an impact
   parameter, compute the new energy and direction.

The implementation is a clean port of **Mode 3** of the legacy
MATLAB code's hard-sphere collision block
(`vmi_sim_3d_neutral_propa_HeDFT_mimic.m`, lines ~595-820), which is
itself based on **Andreas Braun's PhD thesis, section D.2**.

## Position in the dependency chain

```
physics/constants.py
   ↓
physics/collisions.py    ← THIS MODULE   (no SimConfig dependency)
   ↓
simulation/neutral.py     (passes cfg fields in as kwargs)
simulation/ion.py
```

Deliberately a leaf in the dependency graph: no I/O, no plotting, no
config-coupled state.

## Public API

```python
from i2_helium_md.physics.collisions import (
    sample_collision_events,
    apply_collision,
    velocity_dependent_cross_section,
    CollisionDiagnostics,
    temperature_diagnostic_from_collision,
)
```

### `sample_collision_events(...)` — Mode 3

Returns a boolean array marking which particles collide this timestep:

```
p_scatter = distance_travelled * sigma * rho_droplet
collide   = (uniform() < p_scatter) AND (depth < 0) AND (E0 >= E_min)
```

### `apply_collision(...)` — scatter velocities

For colliders, samples `b/R = sqrt(uniform)`, computes COM angle from
hard-sphere geometry, transforms to lab frame, optionally smears
Gaussian, and assembles a new 3D velocity vector. Non-colliders are
returned with their incoming velocity exactly (and `delta_E = 0`).

By default returns the 4-tuple `(vx_new, vy_new, vz_new, delta_E_eV)`.
Pass `return_diagnostics=True` to additionally get a
`CollisionDiagnostics` namedtuple
`(b_collision, COSTHETA, COStheta_lab, rho, E0_eV, E1_eV)` -- the
per-atom internals the ion-propagation driver needs to record the
legacy MATLAB temperature-diagnostic accumulator
(`vmi_sim_3d_ion_propa.m:683`). `COStheta_lab` is the **lab-frame**
post-smearing cosine, which is what the MATLAB diagnostic uses;
`COSTHETA` is kept alongside for completeness. The default keeps the
4-tuple shape so existing neutral-side callers are unaffected.

### `temperature_diagnostic_from_collision(diag)` — recipe helper

Reduces a `CollisionDiagnostics` to the three-element legacy MATLAB
``diagnostic_array`` row::

    [ mean(E1[b]/E0[b]),
      mean((1 + rho[b]**2) / (1 + rho[b])**2),
      mean(arccos(clip(COStheta_lab[b], -1, 1))) ]

where ``b = diag.b_collision``. The third entry uses the **lab-frame**
scattering angle (matching MATLAB ``vmi_sim_3d_ion_propa.m:561``,
where ``theta = acos(COStheta(b_collision)) + smearing``). For heavy
projectile on light target the lab-frame cone is very narrow
(max ``asin(1/rho) ~ 1.81 deg`` for I+ on He), so the angle trace
sits near 0-2 deg rather than spanning the full 0-180 deg COM range.
Returns all-NaN when no atom collided in the step. Used by
``ion_propagation_step`` to write each step's diagnostic row into the
v5 ion checkpoint's ``temperature_diagnostic`` field.

### `velocity_dependent_cross_section(...)` — for ions

Helper for the ion-stage v-dependent model:

```
sigma_per_particle = sigma_0 * v ** exponent
```

Production setting is `exponent = -2` (so `sigma ~ 1/v²`). Pass the
result as the `sigma_angstrom_sq` argument to `sample_collision_events`.

**Edge case at v=0**: with negative exponent this returns `+inf`, not
NaN. That propagates cleanly through `sample_collision_events` because
`trial < +inf` is always True — i.e. an infinitely-slow ion always
collides in any finite step. This is mathematically defensible (an
infinitely-slow particle has infinite mean collision time) and avoids
the need for an arbitrary clamp. The Landau cutoff (`E0 < E_min`) in
`sample_collision_events` provides the physical low-velocity floor.

## Physics details (the parts worth understanding)

### Why `b/R = sqrt(uniform)`?

For a hard sphere, the probability of an impact parameter `b` is
proportional to `2π b db` (annulus area). The CDF is
`F(b) = (b/R)²` and the inverse is `b/R = sqrt(u)`. This means
**grazing collisions are more likely than head-on**, which is
geometrically correct.

### Why `cos(θ_COM) = 2(b/R)² − 1`?

For hard-sphere elastic scattering,
`sin(θ_COM/2) = sqrt(1 − (b/R)²)`. Using the double-angle identity:

```
cos(θ_COM) = 1 − 2 sin²(θ_COM/2) = 1 − 2(1 − (b/R)²) = 2(b/R)² − 1
```

So `b/R = 0` (head-on) gives `cos θ_COM = -1` (backscatter), and
`b/R = 1` (grazing) gives `cos θ_COM = +1` (no deflection).

### Energy transfer formula

For elastic 2-body scattering with incoming speed `v0`, target initially
at rest, mass ratio `ρ = m_proj / m_target`:

```
E1/E0 = (1 + 2ρ cos θ_COM + ρ²) / (1 + ρ)²
```

Limits worth remembering:
- `ρ → ∞` (heavy projectile): `E1/E0 → 1`, no energy loss.
- `ρ = 1` (equal masses): `E1/E0 = (1 + cos θ_COM)/2`, head-on
  (`cos θ_COM = -1`) gives full transfer `E1 = 0`.
- For iodine (127 amu) on helium (4 amu): `ρ = 31.75`,
  max fractional loss is `4ρ/(1+ρ)² ≈ 0.119` (12% per head-on hit).
  Mean loss with uniform `cos θ_COM` is half of that ≈ 6%.

### COM → lab transformation

When the target is a free particle (not the lab frame),

```
cos θ_lab = (cos θ_COM + ρ) / sqrt(1 + 2ρ cos θ_COM + ρ²)
```

Same heavy-projectile limit: `θ_lab → 0`.

## The "azimuthal smearing" convention (worth flagging)

The MATLAB code samples the scattering azimuth with

```matlab
COSBETA = (rand() - 0.5) * 2;       % uniform in [-1, 1]
SINBETA = sqrt(1 - COSBETA^2);
```

This is **NOT** a uniform azimuth. A true uniform azimuth would be

```python
phi = uniform(0, 2π)
COSBETA = cos(phi)
SINBETA = sin(phi)
```

The MATLAB convention concentrates mass at `|β| ≈ π/2` and gives zero
density at `β = 0, π`.

**Why we mirror it exactly:**

- The legacy MATLAB is the project's reference. Diverging silently
  would create a discrepancy between the two implementations that
  would be invisible until someone tried to compare results.
- The `velocity_normal_1` and `velocity_normal_2` axes are themselves
  defined by a **fresh random reference direction per particle per
  step**. The non-uniformity in `(cos β, sin β)` rotates uniformly
  with the random axes, which we hypothesize washes out.

**Verified:** `test_perpendicular_plane_isotropic` (in
`test_collisions.py`) checks that the azimuth `φ = atan2(vz, vy)` of
the scattered velocity for incoming `+x` is uniform in `[-π, π]` to
within 1% on `<cos φ>`, `<sin φ>`, and `<cos 2φ>`. The 2-fold
moment is the canary — if the COSBETA bias did NOT wash out, it
would show up as a non-zero `<cos 2φ>`. Empirically `<cos 2φ> ≈
−0.002` with 200,000 samples, well within statistical noise.

**If we ever want to switch to uniform-φ:** replace
the lines

```python
COSBETA = (rng.uniform(0.0, 1.0, size=n) - 0.5) * 2.0
SINBETA = np.sqrt(np.clip(1.0 - COSBETA**2, 0.0, 1.0))
```

with

```python
phi_az = rng.uniform(0.0, 2.0 * np.pi, size=n)
COSBETA = np.cos(phi_az)
SINBETA = np.sin(phi_az)
```

The unit tests should still pass; comparison against the legacy
MATLAB output would diverge but the physics would be more clearly
correct.

## Why no SimConfig dependency

The module takes plain function parameters: `sigma_angstrom_sq`,
`scatter_mass_amu`, `E_min_eV`, etc. The simulation driver pulls these
from `cfg` and passes them in. This:

- Lets us test the physics in isolation without constructing a full
  `SimConfig`.
- Allows reuse for both neutral and ion stages, which use different
  cross sections and scatter masses.
- Makes the function signature a self-documenting contract for what
  the physics actually depends on.

## What's NOT here

- **Modes 1 and 2** of the legacy collision logic. Mode 1 uses a
  fixed per-step probability; Mode 2 tracks the mean free path. Both
  are considered legacy code paths; we'll add them if needed.
- **Velocity-dependent cross section** (`sigma_lookup(v)`). Not
  needed for the neutral stage; will be added for ions in Step 11.
- **Mass attachment** (helium atoms accreting onto the projectile).
  Not used in neutral propagation per the legacy code's
  `attach_he = 0` setting.

## Statistical regression-test signatures

Locked in by `test_collisions.py`:

| Quantity | Expected | Tolerance |
|---|---|---|
| Empirical collision rate | `d × σ × ρ` | ±0.005 (n=200k) |
| Mean fractional energy loss (I/He) | ~5–7% | range [0, 10%] |
| Max fractional energy loss (equal masses) | → 1.0 | > 0.98 (n=50k) |
| `<cos θ_lab>` | analytic integral | ±0.005 (n=100k) |
| `<cos φ>`, `<sin φ>`, `<cos 2φ>` (perpendicular plane) | 0 | ±0.01 (n=200k) |
| `\|v_new\|² == 2E1/m` | exact | rtol 1e-9 |
| Energy strictly non-increasing for colliders | yes | exact |
| Non-colliders unchanged | byte-identical | exact |
