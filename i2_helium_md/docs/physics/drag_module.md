# The `drag.py` module

## What problem does this file solve?

The hard-sphere collision model (`collisions.py`) treats the helium as a
gas of discrete scatterers that randomly kick the projectile. For the
**ion stage** we are replacing that with a continuous, TDDFT-calibrated
**drag force** for I⁺ moving through the helium bubble: instead of
sampling individual collisions, the ion feels a smooth velocity-dependent
retarding force whose law was *extracted* from TD-DFT trajectories by
force balance.

**Signed force-balance convention (matches upstream).** From the equation
of motion `m_eff · a = F_C + F_drag`, the extracted signed drag is

```
F_drag(t) = m_eff · a(t) − F_C(t)
```

as written in `Drag_extraction_code.md` and §6.5 of the decisions doc.
During the Coulomb explosion the ion accelerates outward (`F_C` outward,
`m_eff·a < F_C`), so the *signed* `F_drag` is negative — it opposes the
motion, as a drag must. The module functions below return the **positive
magnitude** `g·(a·v + b·v³)` (with extracted `a,b > 0`); the integrator
restores the sign and applies `−F_drag` along `v̂`. The earlier
`F_drag = F_C − m_eff·a` phrasing was the magnitude with the sign folded
in and is *not* the signed force-balance equation — do not propagate it.

> **Verify before locking:** confirm the `F_drag_amuAps2` column in
> `drag_data.csv` is stored as the positive magnitude. If it carries the
> signed (negative) value, the lstsq refit would recover sign-flipped
> coefficients and the killer test would silently encode the wrong sign.

This module is **Slice 1** of that port: the pure-physics core that
occupies the collision-model swap point. It provides three functions and
one data type, and nothing else — no integrator, no config, no mass
dynamics, no noise. `collisions.py` is **left intact and importable**;
this module is *additive and parallel*, not a replacement (removing the
collision model needs explicit approval).

Authoritative specs:
- `DRAG_PORT_DESIGN_DECISIONS.md` — all architectural choices.
- `SLICE1_GOALS_gated_drag_module.md` — this slice's contract.
- upstream extraction: `../Drag_function/drag_function/drag_calculation.py`.

## Position in the dependency chain

```
physics/constants.py
   ↓
physics/drag.py        ← THIS MODULE   (no SimConfig dependency, no mass)
   ↓
(Slice 2) BAOAB ion-stage O-step       — consumes γ(v); applies −F_drag/m(t)
(Slice 3) SimConfig enum surface + §6.5 mass↔coefficient guard
(Slice 4) ion-driver rewiring + O-step energy accounting
```

Deliberately a leaf: no I/O, no plotting, no config-coupled state — the
same clean seam as `collisions.py`.

## The one convention that drives everything: γ is a *force coefficient*

There is a single friction convention across the whole drag port. Get
this right and the rest follows:

- **`γ(v)` has units `amu/ps`** and is defined `γ(v) = |F_drag(v)| / v`.
  The friction *force* is `γ(v)·v` — **no leading mass `m`**.
- The friction **rate** is `γ/m` `[1/ps]`; it appears **only** inside the
  BAOAB damping exponent `e^(−γ·dt/m)` (Slice 2), never as a multiplier on
  the force.
- The FDT noise amplitude `√(2·γ·k_B·T_eff)` (Slice ≥3) uses `γ`
  `[amu/ps]` directly.

**Consequence — this module is mass-agnostic.** It never takes `m` as an
argument. Mass enters the simulation in exactly one place: the
integrator's O-step, as one explicit division by `m(t)`. This removes the
factor-of-`m` hazard from the module boundary entirely.

## Governing equations (primary form `linear_cubic`)

With `depth = r_atom − r_droplet` (negative **inside** the droplet,
positive outside — the same convention as
`potentials.droplet_potential`):

```
g(depth)        = ½ · (1 − erf(depth / steepness))     # gate, dimensionless ∈ [0,1]
F_drag(v,depth) = g(depth) · (a·v + b·v³)              # amu·Å/ps²
γ(v,depth)      = g(depth) · (a + b·v²)                # amu/ps  (CLOSED FORM)
```

The functions return the **positive magnitude** form; the integrator
restores the sign and applies `−F_drag` along `v̂` (drag opposes motion).

### Dimensional check (units balance to a force)

| quantity | units | check |
|---|---|---|
| `v` | Å/ps | — |
| `a` | amu/ps | `a·v` → amu·Å/ps² ✓ |
| `b` | amu·ps/Å² | `b·v³` → amu·Å/ps² ✓ |
| `g`, `depth/steepness` | dimensionless | erf argument dimensionless ✓ |
| `F_drag` | amu·Å/ps² | a force; matches Coulomb/droplet accel after ÷m ✓ |
| `γ` | amu/ps | `b·v²` = amu/ps, same as `a` ✓ |

## Public API

```python
from i2_helium_md.physics.drag import (
    DragCoefficients,     # the coefficient-bundle TYPE (frozen dataclass)
    spatial_gate,         # g(depth)
    drag_force,           # F_drag(v, depth)
    drag_gamma,           # γ(v, depth)  -- closed form
    LINEAR_CUBIC, LINEAR_QUADRATIC, THRESHOLD, POWER_LAW,  # form tags
)
```

### `spatial_gate(depth, steepness)` → dimensionless ∈ [0, 1]

`g(depth) = ½(1 − erf(depth/steepness))`: **1** deep inside, **½** at the
nominal surface (`depth = 0`), **0** outside. Smooth and `C¹` — the
continuity the discarded sharp boolean gate lacked. It reuses the same
erf/`steepness` machinery as the confining potential
(`potentials.droplet_potential`), only **complemented** so drag turns
*off* outside the droplet where there is no helium to drag against.
Raises `ValueError` for non-positive `steepness`.

### `drag_force(v, depth, coeffs, steepness)` → amu·Å/ps²

`g(depth)·(a·v + b·v³)` for `linear_cubic`. The Coulomb/droplet
accelerations and this force share units after dividing by `m(t)`.

### `drag_gamma(v, depth, coeffs, steepness)` → amu/ps

`g(depth)·(a + b·v²)`. Later consumed by **both** the O-step rate `γ/m`
and the FDT noise amplitude, and it carries the **same** `g(depth)` that
scales `drag_force` (hard FDT coupling: the noise must be gated wherever
the friction is, or the ion would get thermal kicks in vacuum where it
feels no friction).

### `DragCoefficients` — the coefficient-bundle type

Frozen dataclass; Slice 1 owns the **type** (the §6.5 consistency guard
that *uses* the metadata is Slice 3). Fields:

| field | meaning |
|---|---|
| `form` | drag-form tag (`linear_cubic`, …) |
| `coefficients` | form-tagged, variable-arity, e.g. `{"a","b"}` for `linear_cubic` |
| `extraction_mass_model` | `"constant"` or `"time_resolved"` |
| `extraction_mass_amu` | the effective mass the law was extracted under — **provenance only** |

It validates on construction (unknown form, missing coefficients, bad
mass model, non-positive mass all raise `ValueError`). The module
*consumes* a bundle and never builds one from config.

> **Deliberately deferred fields (Slice 3+), flagged so the type is not
> frozen too narrowly.** Two known additions are *not* in the Slice 1
> bundle and are fine to omit at constant-mass deterministic Tier 0, but
> the type should leave room for them:
> - **`a_err`, `b_err`** — the 1σ fit uncertainties are present in
>   `fit_parameters.json` but not carried on the bundle. They feed the
>   Tier 3 ensemble uncertainty band (§6 of the decisions doc), so Slice 3
>   will want them on the type.
> - **A time-resolved `m(t)` handle** — `extraction_mass_model` already
>   admits `"time_resolved"`, but only the scalar `extraction_mass_amu` is
>   stored. A `"time_resolved"` bundle needs a reference to the `m(t)`
>   profile it was extracted under for the §6.5 guard to match it against
>   an evolving-mass scenario. Out of Slice 1 scope (constant-mass only),
>   but the field shape is a Slice 3 decision.

## The closed-form γ (a physics point, not a coding nicety)

`γ` is exposed via its **closed form** `g·(a + b·v²)`, **never** via
`|F_drag|/v`. The two are analytically identical, but the division
manufactures a `0/0` singularity at `v → 0` that the `linear_cubic` form
does **not** have — there `γ → g·a`, finite. The O-step and the FDT noise
both need `γ` near rest, where the physics is well-behaved, so computing
it by division would invent a singularity the law does not contain.

This is exactly the `v → 0` boundary where the **`power_law`** form (with
`n < 0`) *genuinely* diverges and would need a low-velocity floor — which
is why the closed-form-vs-division choice is made **explicit per form**
rather than globally.

## Form dispatch — one realised, three reserved

```
coeffs.form?
├─ linear_cubic     → implemented fully
├─ linear_quadratic → raise NotImplementedError  (no fit pass yet)
├─ threshold        → raise NotImplementedError  (no fit pass yet)
└─ power_law        → raise NotImplementedError  (needs the low-v floor;
                                                   out of Tier-0 scope)
```

Only `linear_cubic` is realised — it is the form with coefficients in
hand for both bubble sizes. The other three are reserved behind the same
dispatch so adding them later is a *branch*, not a signature change.
Refusals are **explicit** (`NotImplementedError` with a reason), never a
silent wrong answer.

> **Note on `power_law`.** The design docs anticipated `n ≈ −2` (a
> hard-sphere `σ ∝ v⁻²` artifact, singular at rest). The actual exported
> extraction gives `n ≈ +2` (9 Å: 18 Å similar), i.e. *no* low-`v`
> divergence. This does not affect Slice 1 (the form is deferred either
> way), and the low-`v` "contrast against divergence" is asserted against
> the *hypothetical* `n < 0` form, not the real export.

## Calibration data and provenance

The `linear_cubic` coefficients come from the upstream TD-DFT extraction
and are frozen under `data/reference/drag/<case>/linear_and_cubic/`:

```
data/reference/drag/9A/linear_and_cubic/{fit_parameters.json, drag_data.csv}
data/reference/drag/18A/linear_and_cubic/{fit_parameters.json, drag_data.csv}
```

- `fit_parameters.json`: `{a, b, a_err, b_err, meff_amu}`.
- `drag_data.csv`: `t_ps, v_spline_Aps, F_drag_amuAps2` — the **full**
  extraction window (the trusted interior is recovered by dropping the
  first/last 500 points, matching the extraction's truncation).

| case | a (amu/ps) | b (amu·ps/Å²) | m_eff (amu) |
|---|---|---|---|
| 9 Å | 13.860 | 2.581 | 202.954 (I + 19 He) |
| 18 Å | 14.556 | 2.053 | 202.954 (I + 19 He) |

`m_eff ≈ 203 amu` is the **mass the drag law was extracted under**, the
window-representative ~19-He shell — *not* "the true mass of the ion."
The module never uses it; it is carried as provenance so the Slice 3
`mass_scenario ↔ coefficients` guard can later check consistency.

> **Verify before locking — the provenance check the killer test cannot
> catch.** The stamped `meff_amu = 202.954` must be the mass the force
> balance *actually used* when `{a,b}` and `drag_data.csv` were produced.
> An earlier extraction literal was `179.912`. Because `drag_data.csv` and
> `fit_parameters.json` come from the *same* run, the lstsq killer test is
> self-consistent regardless of which mass was used — it cannot detect a
> relabelled mass. If the coefficients were fit at 179.912 and only the
> JSON field was bumped to 203, the stamped provenance is false, and the
> Slice 3 guard would validate `fixed`-mass runs against a mass the law
> was never extracted under (a ~13% error feeding
> `F_drag = m_eff·a − F_C`). **Confirm the extraction was re-run at
> 202.954, not relabelled.** The value itself is consistent:
> I-127 (126.90) + 19·He (19 × 4.0026 = 76.05) ≈ 202.95 amu.

## Why no SimConfig dependency

`coeffs` and `steepness` arrive as plain function arguments (the same
kwarg pattern as `collisions.py`). The future driver pulls them from
`cfg` and passes them in. This keeps the physics testable in isolation,
keeps the function signature a self-documenting contract for what the law
depends on, and keeps every model choice behind its own swappable surface
rather than hard-wired here.

## What's NOT here (Slice 1 scope fence)

- **No integrator.** The BAOAB O-step, the damping exponent `e^(−γ·dt/m)`,
  and the OU update are Slice 2.
- **No `SimConfig` fields.** Only the bundle *type* lives here; the enum
  surfaces (`drag_form`, `drag_spatial_gate`, …) and the §6.5 guard are
  Slice 3.
- **No driver wiring.** Replacing the collision call sites in
  `ion_propagation_step` is Slice 4.
- **No checkpoint changes / energy rename.** The `IonCheckpoint` v6 bump
  and `E_mass_attach_defect_eV → E_mass_transfer_eV` come with mass
  dynamics.
- **No noise.** The FDT amplitude is defined-but-inactive; Slice 1 only
  ensures `γ` is exposed in the gated, closed form the noise will consume.
- **No mass dynamics.** The module is mass-free by construction; `m(t)`
  lives entirely in the integrator.

## Regression-test signatures

Locked in by `tests/test_drag.py` (parametrized over the 9 Å and 18 Å
cases):

| Quantity | Expected | Tolerance |
|---|---|---|
| **Killer test** — lstsq refit of `a·v + b·v³` to the trusted-interior scatter recovers stored `{a,b}` | the module *is* the extracted law | rel 1e-4 |
| `drag_force` vs. raw `(v, F_drag)` scatter | within the genuine fit residual (9 Å R² ≈ 0.974) | rel-RMS < 7% |
| Dissipativity `F_drag·v ≥ 0` over 0–30 Å/ps | drag never adds energy | exact |
| Guard not vacuous: in-hand `a>0`, `b≥0` | no turnover speed `v† = √(−a/b)` | exact |
| Gate limits `g(−∞), g(0), g(+∞)` | `1, ½, 0` | approx |
| Gate monotone-decreasing in depth; `∈ [0,1]` | yes | exact |
| Gate `C¹`: central FD vs. analytic slope | match | atol 1e-7 |
| FDT carrier: `drag_force` and `drag_gamma` share the same `g(depth)` (so `F = γ·v`) | identical gate factor | rtol 1e-12 |
| Low-`v` regularity: `γ(0⁺) → g·a`, finite | no floor needed for `linear_cubic` | approx |
| Mass-agnosticism: no function takes a mass argument | enforced | exact |
| Reserved forms (`linear_quadratic`, `threshold`, `power_law`) | raise `NotImplementedError` | exact |

> **Max-trajectory-speed open item is dormant, not resolved.** The
> `linear_cubic` turnover guard `v† = √(−a/b)` only bites when `b < 0`.
> Both extracted cases have `b > 0`, so there is no real turnover and the
> §6.10 "source a max trajectory speed" requirement is moot *for these
> coefficients*. A future re-extraction (different case, charge state, or
> the `time_resolved`-m(t) refit) that yields `b < 0` reactivates it —
> the guard then needs a real max-speed ceiling to test `v†` against.
