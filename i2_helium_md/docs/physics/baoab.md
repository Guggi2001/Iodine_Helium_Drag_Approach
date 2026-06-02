# The `baoab.py` module

## What problem does this file solve?

The baseline ion stage integrates with `velocity_verlet_step`
(`leapfrog.py`), which is kick-drift-kick and assumes a **position-only**
force `F = F(x)`. The drag-model port breaks that assumption twice over:
the drag force is **velocity-dependent** `F_drag(v)`, and (later) the
Langevin noise turns the equation of motion into a stochastic
differential equation. An explicit Verlet treatment of either is only
*conditionally* stable and has no statistically-correct site for the
noise (decisions doc §4.2).

This module is **Slice 2** of the drag port: the **BAOAB operator-split
ion-stage stepper** that replaces `velocity_verlet_step` *for the ion
stage only*. It splits one timestep into a conservative Verlet
kick/drift (the **B** and **A** sub-steps, handling Coulomb + droplet via
the existing acceleration callable) and an Ornstein–Uhlenbeck sub-step
(the **O** sub-step, handling drag — and, when activated, noise),
composed symmetrically as **B–A–O–A–B**. It consumes Slice 1's `γ(v)`
(`drag.py`) as its only new physics input.

The neutral stage is untouched and keeps `velocity_verlet_step` (the
decisions-doc §4.6 accepted asymmetry: the neutral stage has no drag and
no noise, so BAOAB buys nothing there).

Authoritative specs:
- `DRAG_PORT_DESIGN_DECISIONS.md` §4 — why operator splitting (Option III).
- `SLICE2_GOALS_baoab_ion_stepper.md` — this slice's contract.
- `docs/physics/drag_module.md` — the `γ(v)` this module consumes.
- `docs/physics/leapfrog_module.md` — the baseline integrator and the
  `_kick`/`_drift` primitives shared with this module.

## Position in the dependency chain

```
physics/constants.py
   ↓
physics/drag.py          ← Slice 1: γ(v), mass-free
   ↓
physics/baoab.py         ← THIS MODULE (Slice 2): B–A–O–A–B; mass enters here
   ↑  (reuses)
physics/leapfrog.py        _kick / _drift primitives, _ion_accel_fn
   ↓
(Slice 3) SimConfig enum surface + §6.5 mass↔coefficient guard
(Slice 4) ion-driver rewiring + O-step energy accounting (eV conversion)
```

## The BAOAB scheme

One symmetric (Strang) step over `dt`, ion stage only. Let
`acc_fn(x) → (a_cons, E_pot)` be the *existing* conservative ion
acceleration (Coulomb + droplet, `leapfrog._ion_accel_fn`), unchanged.

```
state (x0, v0)
  B   v ← v + (dt/2)·a_cons(x)         half-kick, conservative force only
  A   x ← x + (dt/2)·v                 half-drift
  O   v ← e^(−γ·dt/m)·v   (+ noise)    drag damping (+ dormant noise site)
  A   x ← x + (dt/2)·v                 half-drift
  B   v ← v + (dt/2)·a_cons(x)         half-kick, conservative force only
return (x1, v1, E_pot, ΔE_dissip)
```

- **B (kick)** and **A (drift)** are *exactly* velocity-Verlet's kick and
  drift — this is the fact the anchor test exploits.
- **O** is the only genuinely new physics. Conservative forces never enter
  O; drag never enters B.
- `acc_fn` is evaluated **exactly twice per step** — once at the entry
  position (first B) and once at the post-step position (last B) — matching
  the two evaluations `velocity_verlet_step` makes. This count is **asserted**
  in the anchor test (`calls == 2·n`, equal to Verlet's), so the eval-for-eval
  parity is a tested contract, not a comment. The two evaluations are
  re-computed each step (no cross-step force caching); the choice is a
  performance/cleanliness one the anchor test guards either way (§11
  open item, resolved to *re-evaluate*). A future final-B→next-first-B cache
  would drop the count to ≈`n+1` and update that assertion deliberately.

## The one convention that drives everything: mass enters *here*, in amu

`drag.py` is mass-agnostic — `γ(v)` is a force coefficient in `amu/ps`,
defined `γ(v) = |F_drag(v)|/v`, with **no leading mass**. This module is
**the single place mass enters the drag model**, and it does so in **amu**:

- **`m` arrives in amu** (not kg). The damping exponent `γ·dt/m` is then
  dimensionless directly: `(amu/ps)·(ps)/amu = 1`.
- Mass appears in exactly two spots: the O-step damping exponent
  `e^(−γ·dt/m)` and the dissipated-energy bookkeeping `½·m·(…)`. The
  conservative inertia (`F = m·a`) is already handled inside `acc_fn`.
- **Locked unit contract:** the module is pure-mechanical. `ΔE_dissip` is
  returned in **amu·Å²/ps²**, *not* eV — Slice 4 owns both the kg↔amu
  plumbing at the driver boundary and the eV conversion of the
  dissipated-energy accumulator. No eV factor and no kg appear in this file.

## Governing O-step

### Tier-0 (noise off, the production deterministic integrator)

```
v_out = e^(−γ(v_in, depth)·dt/m) · v_in
ΔE_dissip = ½·m·(‖v_in‖² − ‖v_out‖²)  =  ½·m·‖v_in‖²·(1 − e^(−2γ·dt/m))   ≥ 0
```

per atom. Because Slice 1 guarantees `γ ≥ 0` for `linear_cubic`
(`a, b > 0`), the damping factor `e^(−γ·dt/m) ∈ (0, 1]`, so
`‖v_out‖ ≤ ‖v_in‖` — **the O-step can never add kinetic energy, exactly**,
and `ΔE_dissip ≥ 0` by construction.

### Full form (for reference — noise dormant at Tier 0)

```
v_out = e^(−γ·dt/m)·v_in + √( (k_B·T_eff/m)·(1 − e^(−2γ·dt/m)) ) · ξ ,   ξ ~ N(0, 1)
```

`ξ` is a dimensionless standard normal. At `T_eff = 0` the second term
vanishes and this reduces to the Tier-0 form. The gate `g(depth)` lives
*inside* `γ` (Slice 1: `γ = g·(a + b·v²)`), so it scales **both** the
damping and the noise amplitude — the FDT coupling, automatic and free.

### Dimensional check

| expression | units | balance |
|---|---|---|
| damping exponent `γ·dt/m` | `(amu/ps)(ps)/amu` | dimensionless ✓ |
| damped velocity `e^(−γ·dt/m)·v` | (dimensionless)·Å/ps | Å/ps ✓ |
| half-kick `(dt/2)·a_cons` | (ps)(Å/ps²) | Å/ps ✓ |
| half-drift `(dt/2)·v` | (ps)(Å/ps) | Å ✓ |
| `ΔE_dissip = ½·m·(v_in² − v_out²)` | (amu)(Å/ps)² | amu·Å²/ps², an energy; ≥ 0 ✓ |

## The asymmetric γ-freeze (intentional — do not "fix" it)

The nonlinear O-step is frozen to an `O(dt)` approximation by evaluating
`γ` **once**, at a deliberately *mixed* point:

- **Velocity → frozen at the O-step input velocity** `v_in^O` (the
  post-first-B velocity; the first A does not change velocity). This is
  the only choice available without an implicit solve, and it makes the
  returned `ΔE_dissip` *exact for the frozen γ*: the energy removed is
  computed with the same `γ` that damped the velocity.
- **Depth (gate) → evaluated at the current O-step position** (after the
  first half-drift). In B–A–O–A–B the position is *already* updated when O
  runs, so `depth = r_atom − r_droplet` is freshly and legitimately known.

**Consequence (intentional):** `γ` is evaluated at a mixed
`(v_old, depth_new)` point. This is correct — they are independent
arguments to a pure function. Freezing depth at step-entry instead would
let a fast surface-crossing gate carry a stale wrong-side `g` and forfeit
the exact never-adds-energy property for no gain. This is documented in
the module docstring so it is not "corrected" later.

## Public API

```python
from i2_helium_md.physics.baoab import make_ion_baoab_step
```

### `make_ion_baoab_step(m, droplet_radii, acc_fn, gamma_fn, *, T_eff=0.0, rng=None)` → `step`

Builds one BAOAB ion-stage step closure. Mirrors
`leapfrog.make_ion_step`, **including its call convention**: `dt` is *not*
bound in the factory — it is passed to `step` on each call, sourced from
`SimConfig.dt_ion` by the Slice-4 driver (as `ion_propagation_step.py`
does `dt = cfg.dt_ion`).

| parameter | meaning |
|---|---|
| `m` | per-atom physical mass **in amu**, shape `(2N,)`. Enters only the O-step exponent and the energy bookkeeping. |
| `droplet_radii` | per-atom droplet radius (Å), shape `(2N,)`; forms `depth = r_atom − droplet_radii` for the gate. |
| `acc_fn` | conservative ion acceleration `(x,y,z) → ((ax,ay,az), E_pot)` (the existing `_ion_accel_fn`); evaluated twice per step. |
| `gamma_fn` | friction coefficient `(speed, depth) → γ` in amu/ps, shape `(2N,)`. Tier-0: Slice 1 `drag_gamma` closed over `coeffs` + `steepness`. |
| `T_eff` | effective noise temperature. Tier-0 default `0.0` (noise dormant). `> 0` raises `NotImplementedError` (Slice ≥3). |
| `rng` | reserved for the future noise draw; **never consumed at `T_eff = 0`**. |

Returns the closure:

```python
step(pos, vel, dt) -> (new_pos, new_vel, E_pot_per_pair, dE_dissip)
```

- `dt`: timestep in picoseconds, passed on the call (`SimConfig.dt_ion`,
  default 0.01 ps) — mirroring `make_ion_step`'s `step(pos, vel, dt)`.
- `E_pot_per_pair`: shape `(N,)`, eV — from the final B's `acc_fn`
  evaluation at the new position, matching `velocity_verlet_step`.
- `dE_dissip`: shape **`(2N,)`**, **amu·Å²/ps²**, per atom, `≥ 0`.

**Per-step closure rebuild.** Slice 4 will rebuild this closure every step
because **mass** changes under future mass scenarios, exactly as
`ion_propagation_step.py` rebuilds `make_ion_step`; `dt` itself rides on
the per-call argument, not the rebuild. The factory is built to support
that even though Tier-0 mass is fixed, so mass dynamics drop in later with
no restructuring.

## The shared `_kick`/`_drift` extraction (the one baseline touch)

BAOAB's **B** and **A** *are* velocity-Verlet's kick and drift. To avoid a
duplicate integrator implementation (CLAUDE.md rule 1), the kick and drift
were factored out of the monolithic `velocity_verlet_step` into
`leapfrog._kick` / `leapfrog._drift`, and both integrators now call them:

```
_kick(vel, acc, c)  →  v + c·a      (component-wise)
_drift(pos, vel, c) →  x + c·v      (component-wise)
```

`velocity_verlet_step` is rewritten as half-kick / full-drift / half-kick
through these primitives — analytically the identical KDK update.

This is the **one knowing touch to the frozen baseline integrator**
(`PHYSICS_BASELINE.md` §5 noted the drag port "should not need to touch
`leapfrog.py`"). It is accepted because decisions-doc §4.6 — the newer,
more specific decision — added an ion-stage integrator at all, and because
the extraction is a pure refactor made **self-verifying by the anchor
test**: if extracting kick/drift changed `velocity_verlet_step`'s
behaviour, the anchor test (below) fails loudly. Recorded here and in the
module docstring.

## Noise — dormant at Tier 0

- Noise is **off** at `T_eff = 0` (the Tier-0 default): the second O-step
  term is identically zero and **no random number is drawn**.
- The O-step is *structured* so the Langevin term slots in at the single
  marked site (`_o_step`) without touching B, A, or the energy-return
  shape — and its injected energy will be tracked in a **separate**
  accumulator from `ΔE_dissip` (the noise *injects* energy to the thermal
  floor; drag *removes* it).
- Activating it (`T_eff > 0`) is Slice ≥3 work: the FDT amplitude, the RNG
  draw-order pinning, and the second noise-energy accumulator. Until then
  `T_eff > 0` raises `NotImplementedError` rather than silently running an
  unvalidated path (`T_eff < 0` raises `ValueError`).
- The second-energy-channel return slot is **deferred** (§11 open item):
  the Slice 2 return stays a 4-tuple; adding the noise-energy term is a
  one-time signature bump at Slice ≥3.

## What's NOT here (Slice 2 scope fence)

- **No `SimConfig` fields.** Parameters arrive as function arguments; the
  enum surfaces and the §6.5 guard are Slice 3.
- **No driver wiring.** Replacing the collision call sites in
  `ion_propagation_step.py` and rebuilding the closure per step in the
  *driver* is Slice 4. Slice 2 provides the factory; Slice 4 calls it.
- **No checkpoint changes / energy rename.** `IonCheckpoint` v6 and
  `E_mass_attach_defect_eV → E_mass_transfer_eV` come with mass dynamics.
- **No active noise.** Dormant at `T_eff = 0`; structured for later.
- **No mass dynamics.** Mass fixed at Tier 0; the per-step rebuild pattern
  is in place so dynamics drop in later.
- **No eV conversion.** The stepper stays mechanical (amu·Å²/ps²); Slice 4
  converts.
- **Neutral stage untouched.** The only `leapfrog.py` change is the
  behaviour-preserving kick/drift extraction.

## Regression-test signatures

Locked in by `tests/test_baoab.py`:

| Quantity | Expected | Tolerance |
|---|---|---|
| **Anchor (killer) test** — `γ = 0`, noise off: BAOAB step ≡ baseline `velocity_verlet_step` on identical state + `acc_fn`, over many steps (positions, velocities, `E_pot`) | exact agreement (no-drag limit recovers the frozen baseline; also proves the kick/drift extraction is behaviour-preserving) | rtol 1e-12 |
| **Eval-for-eval parity** (asserted inside the anchor test) — `acc_fn` is called **exactly twice per step** (entry + post-step position), the same count `velocity_verlet_step` makes | `calls == 2·n == verlet_calls` (present no-cache contract; the trip-wire for the future final-B caching, which will drop it to ≈`n+1`) | exact |
| Analytic decay — constant linear `γ`, no conservative force: `v(t) = v₀·e^(−γt/m)` over many steps | matches closed form (O-step exact for constant `γ`) | rtol 1e-12 |
| Mass enters only via `γ/m` — doubling `m` only halves the damping rate; `acc_fn` output never scaled by `m` | `v = e^(−γ·dt/m)·v₀` | rtol 1e-12 |
| Dissipation identity + units — returned `ΔE_dissip = ½·m·(‖v_in‖² − ‖v_out‖²)`, shape `(2N,)`, in amu·Å²/ps² | exact (raw mechanical value; an eV conversion would be ~1e-4× and break the match) | rtol 1e-12 |
| Dissipativity — O-step never raises `‖v‖` for `γ ≥ 0`; `ΔE_dissip ≥ 0` | `‖v_out‖ ≤ ‖v_in‖` | exact |
| Gate-off vacuum — real `drag_gamma`, atoms far outside droplet (`g → 0 ⇒ γ → 0`): O-step is identity, `ΔE_dissip = 0` | no drag, no dissipation outside the droplet | rtol 1e-12 |
| Noise dormant — `T_eff = 0` with a real seeded `rng` passed: `rng.bit_generator.state` is **snapshotted before/after** and unchanged (the dormancy proof — *not* the output equality, which holds at `T_eff = 0` regardless) | no RNG draw consumed | exact |
| Active noise refused — `T_eff > 0` | raises `NotImplementedError` | exact |
