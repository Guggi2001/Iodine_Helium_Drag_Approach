# Slice 4 — Ion-driver rewiring + O-step energy accounting

## What problem does this slice solve?

Slices 1–3 of the drag-model port were **additive**: they delivered a pure
physics module (`physics/drag.py`), a BAOAB stepper (`physics/baoab.py`), and
the `SimConfig` drag surface + config-load guard + coefficient loader — but
none of it was ever *called*. A drag-enabled preset still ran the hard-sphere
collision path.

Slice 4 is the **wiring slice**. It makes a drag-enabled preset run end-to-end,
deterministically, at **Tier 0** (`mass_scenario=fixed`, noise off,
`linear_cubic`), producing the first runnable drag trajectory. It is the first
slice that *exercises* the scoped CLAUDE.md "do not change collision physics"
exception rather than merely declaring it.

**Deliverable: mechanical correctness, not a TDDFT match.** Energy closes, the
checkpoint round-trips, and a deterministic smoke harness passes. The TDDFT
comparison (`compare_distance` / `compare_velocity_magnitude` vs
`9A_All_Data.csv`) is the *separate next task* that sets the deferred §6.10
thresholds; it is deliberately **not** part of this slice.

Authoritative specs:
- `SLICE4_GOALS_ion_driver_rewiring.md` — this slice's contract.
- `DRAG_PORT_DESIGN_DECISIONS.md` — §2.9 (energy invariant), §4 (BAOAB),
  §5.5 (gate collapse), §6.4/§6.5 (Tier 0, mass↔coefficient consistency).
- `docs/config_and_preset.md` — Slice 3 surface this slice consumes.

## Position in the dependency chain

```
physics/leapfrog.py  ── make_ion_accel_fn ──┐
physics/baoab.py     ── make_ion_baoab_step ┤
physics/drag.py      ── drag_gamma ─────────┤
physics/constants.py ── U, EV ──────────────┤
                                            ↓
simulation/ion_propagation_step.py   baoab_propagation_step + _check_drag_scope + lifts
                                            ↓
simulation/ion.py    dispatch (drag_coefficients is not None) + per-step closure build
```

No import cycles: `baoab.py` depends only on `leapfrog.py`;
`ion_propagation_step.py` imports the `BaoabStep` type from `baoab.py`; `ion.py`
imports the factories from `physics/` and the new step + guard from
`ion_propagation_step.py`.

---

## The architecture: a parallel path, not an in-place switch

The central risk-containment choice (spec §2): the existing collision per-step
function is **not rewritten in place**. Slice 4 adds a *third sibling* per-step
function next to `neutral_propagation_step` and `ion_propagation_step`, and the
driver picks between them **once per run**:

```
ion.py  (run_ion_propagation)
  ├─ cfg.drag_coefficients is not None  →  build BAOAB closure, loop baoab_propagation_step
  └─ else                               →  existing collision path, loop ion_propagation_step  (UNCHANGED)
```

- **Dispatch predicate: `cfg.drag_coefficients is not None`.** The one condition
  that is both *necessary* (BAOAB cannot run without coefficients) and *already
  validated* (a non-`None` bundle passed the Slice 3 `check_drag_config`:
  consistency + dissipativity + form agreement). No separate `use_drag` flag
  that could disagree with the coefficients' presence.
- Evaluated **once**, before the step loop — not re-checked per step.
- Consequence: "drag preset" ≡ "has coefficients." A/B scenario comparison uses
  *different presets*, not the same preset toggled.

Because the hard-sphere path stays bit-identical and uncalled-for-drag, the
regression surface the drag model will eventually validate against (the §6
hard-sphere-variance calibration) survives untouched.

---

## What was implemented, file by file

### 1. `physics/leapfrog.py` — expose the bare conservative acceleration

`make_ion_baoab_step` needs the *position-only* conservative ion acceleration
(`acc_fn`) for its B/A kicks. Previously that lived only inside `make_ion_step`,
bundled into a velocity-Verlet closure and never exposed.

Added **`make_ion_accel_fn(cfg, mass, droplet_radii, charge, state_ids=None)
-> AccelFn`** and refactored `make_ion_step` to consume it, so there is a single
source of the ion force assembly (CLAUDE.md rule 1). This is a *behavior-
preserving* touch to the frozen integrator file — the same discipline as the
Slice 2 `_kick`/`_drift` extraction — guarded by the leapfrog + ion-step
regression tests staying green.

> **Scope note.** `SLICE4_GOALS §11` literally fences off `physics/` changes.
> This factory is a `physics/` edit; it was taken deliberately (user-approved)
> as the clean rule-1 option, by precedent of the Slice 2 leapfrog touch, and
> because `make_ion_step`'s output is byte-identical after the refactor.

### 2. `simulation/ion_propagation_step.py` — new sibling, lifts, scope guard

**Behavior-preserving lifts (spec §3.1).** Three physics-free helpers extracted
and the existing collision step routed through them (gated on
`test_ion_propagation_step.py` staying green):

| helper | meaning |
|---|---|
| `_depth(x, y, z, droplet_radii)` | per-atom `r_atom − r_droplet` (negative inside) |
| `_E_kin_eV(mass_kg, vx, vy, vz)` | `½ m v²` in eV (`v` in m/s via `×100`) |
| `_E_pot_per_atom(depth, E_pot_per_pair, cfg)` | ion-droplet binding + half partner Coulomb |

**`baoab_propagation_step(state, *, step, cfg, droplet_radii) -> IonStepState`.**
The Tier-0 drag-branch analog of `ion_propagation_step` — symmetric, but with
the *middle* replaced: no collision sampling, no `apply_collision`, no mass
attachment. Per call it:

1. runs the pre-built BAOAB closure `step(pos, vel, dt)` →
   `(pos', vel', E_pot_per_pair, ΔE_dissip)`;
2. computes `depth` and the eV `E_kin`/`E_pot` via the shared lifts;
3. converts `ΔE_dissip` (amu·Å²/ps²) to eV and accumulates `E_dissip_eV`
   (see "The eV conversion" below);
4. applies the Tier-0 checkpoint fills.

It is **thin**: closure construction lives in the driver (§3 below); this
function does per-step accounting only.

**`_check_drag_scope(cfg)`.** The drag-branch analog of the collision path's
`_check_scope`. It asserts the Tier-0 *runnability envelope* (distinct from
Slice 3's `check_drag_config`, which validates config *internal consistency*).
It raises `NotImplementedError` when any of these hold:

```
noise_form != "none"            (active Langevin noise is Tier 3)
mass_scenario != "fixed"        (mass dynamics is Tier 1)
drag_form != "linear_cubic"     (other forms NotImplemented in drag.py)
effusive_dynamics / single_charge_ionization_allowed / additional_droplet_charges > 0
```

> The spec phrases the first reject as "`T_eff > 0`". There is no literal
> `T_eff` field; it maps to `noise_form != "none"`. The driver passes
> `T_eff=0.0` into `make_ion_baoab_step` at Tier 0.

### 3. `simulation/ion.py` — dispatch + per-step closure construction

After the existing driver scope check, `run_ion_propagation`:

- sets `use_drag = cfg.drag_coefficients is not None`;
- if drag: calls `_check_drag_scope(cfg)`, resolves the gate steepness via the
  new `_drag_gate_steepness(cfg)`, and builds
  `gamma_fn = partial(drag_gamma, coeffs=…, steepness=…)`;
- in the step loop, on the drag branch, **rebuilds the BAOAB closure every
  step** (matching the `make_ion_step` rebuild pattern; Tier-1-ready though
  Tier-0 mass is fixed):

```python
acc_fn = make_ion_accel_fn(cfg, state.mass_kg, droplet_radii, charge)
step   = make_ion_baoab_step(state.mass_kg / U, droplet_radii, acc_fn, gamma_fn, T_eff=0.0)
new    = baoab_propagation_step(state, step=step, cfg=cfg, droplet_radii=droplet_radii)
```

The collision branch (`ion_propagation_step` + `prev_distance` tracking) is
unchanged and draws RNG exactly as before; the drag branch needs no
`prev_distance` and draws no RNG.

**`_drag_gate_steepness(cfg)` — the §5.5 G4→G2 collapse.** The erf gate itself
lives inside `drag_gamma` (via its `steepness` argument); this helper only
chooses *which* steepness to pass:

| `drag_spatial_gate` | steepness used | rationale |
|---|---|---|
| `density_proportional` (default) | `cfg.potential_steepness` | no density profile yet → collapses to the erf complement (G2) |
| `erf_tied` | `cfg.potential_steepness` | the explicit erf-complement gate — *identical today* |
| `erf_independent` | `cfg.drag_gate_steepness` | separate steepness (defaults to `potential_steepness`) |
| `sharp` | — (raises) | discarded G1: a discontinuous force breaks the BAOAB O-step |

The trap this avoids: handling only `erf_tied` would let the *default*
`density_proportional` preset fall through with no gate built.

### The eV conversion (spec §4)

`make_ion_baoab_step`'s O-step returns `ΔE_dissip` in **amu·Å²/ps²** (Slice 2's
deliberate pure-mechanical handoff). Slice 4 converts to eV using the **same**
idiom as the baseline mass-attach defect (`ion_propagation_step.py`):

```
dE_dissip_eV = dE_dissip * U * 100**2 / EV
#              amu→kg (U) · (Å/ps→m/s)² (100²) · J→eV (EV)
```

This keeps `E_dissip_eV` in a consistent eV with `E_kin_eV`/`E_pot_eV`. A
factor off by `U` or `100²` would make energy *appear* to close within the drag
channel while drifting against `E_kin` — guarded directly by a test asserting
the conversion against a known `½ m v²` in both unit systems (rtol 1e-12).

### Checkpoint — v5 retained (spec §5)

No schema bump (the v6 rename `E_mass_attach_defect_eV → E_mass_transfer_eV` is
tied to Tier-1 mass dynamics). Under the drag branch the v5 fields take:

| field | Tier-0 drag value |
|---|---|
| `E_dissip_eV` | the live, cumulative drag dissipation |
| `E_mass_attach_defect_eV` | `0` (no mass attachment — true, not a placeholder) |
| `temperature_diagnostic` | all-`NaN` `(3,)` (no collision events to diagnose) |
| `mass_history_kg` | constant (fixed mass) |
| `number_of_collisions` | `0` |

The checkpoint round-trips through `save_ion_checkpoint` / `load_ion_checkpoint`
with these fills intact (the NaN diagnostic survives the `.npz` round-trip).

---

## The energy-closure invariant (the real acceptance criterion)

At Tier 0, `E_mass_transfer = 0`, so the §2.9 invariant reduces to

```
E_kin(t) + E_pot(t) + E_dissip(t) ≈ const.
```

Closure is **tight, not loose**. The O-step removes kinetic energy and books it
to `E_dissip` *analytically* (`½m(‖v_in‖²−‖v_out‖²)`), so the dissipative part
is exact; the only drift is the baseline conservative Verlet drift on
`E_kin + E_pot`. The smoke test therefore asserts closure to the same tolerance
as baseline Verlet (sub-1%), making it a sharp instrument rather than a loose
sanity check.

---

## Tests

**`tests/test_baoab_propagation_step.py`** — focused unit coverage of the
Slice-4-specific pieces:
- eV conversion vs a known `½ m v²` in both unit systems (rtol 1e-12); cumulative
  accumulation;
- Tier-0 fills (`E_mass_attach_defect_eV`/`number_of_collisions` stay 0, mass
  fixed, `temperature_diagnostic` all-NaN `(3,)`);
- gate collapse (`density_proportional` and `erf_tied` resolve to the same
  steepness; `erf_independent` uses its own; `sharp` raises);
- per-step closure rebuild ≡ build-once at fixed mass (round-off identical).

**`tests/test_ion_drag_smoke.py`** — the spec §8 "thorough debug" artifact: a
tiny fixed-seed drag run on a **synthetic `NeutralCheckpoint`** (decoupled from
the neutral stage), `N=2`, 20 internal steps. Asserts dispatch to the drag
branch, finite trajectories, tight energy closure, monotone positive
dissipation, v5 checkpoint round-trip, and scope-guard rejection of
`noise_form`/`mass_scenario`. It asserts **nothing** about TDDFT.

### Results at implementation

- New Slice 4 tests: **14 passed**.
- Regression (`test_ion_propagation_step`, `test_leapfrog`, `test_baoab`,
  `test_ion`, `test_ion_initial_state`): **85 passed** — the factory refactor
  and §3.1 lifts are behavior-preserving.
- Full suite: **545 passed, 4 skipped, 6 failed**. The 6 failures are
  **pre-existing and unrelated** to Slice 4 — they come from a committed,
  stale `data/runs/single_pulse_droplet/ion.npz` written at schema v4 while the
  code expects v5, loaded by post-processing smoke tests. Verified to reproduce
  on the base branch with the Slice 4 changes stashed.

---

## What Slice 4 does NOT touch (scope fence)

- No in-place rewrite of `ion_propagation_step` (only the §3.1 lifts).
- No checkpoint schema bump (v5 retained; v6 is Tier 1).
- No mass dynamics — the mass-attach block stays in the uncalled-for-drag
  collision path (not deleted, not made dormant).
- No active noise (`T_eff=0`; the scope guard rejects otherwise).
- No TDDFT validation / no threshold setting.
- No new `SimConfig` fields (Slice 4 *reads* the Slice 3 fields).
- The neutral stage is untouched.

---

## Remaining risks and deferred items

1. **Initial ion mass in the production drag path — FIXED**
   (`SLICE4_FIX_initial_mass_consistency.md`). The original defect:
   `build_initial_ion_state` inherited `mass_kg` from the neutral checkpoint
   (bare iodine, ~127 amu) and `baoab_propagation_step` integrated at that mass,
   while `DRAG_PORT_DESIGN_DECISIONS.md §2.2/§6.5` require the `fixed` scenario to
   run at `m_eff ≈ 203 amu` (the law's extraction mass). Resolved by three
   coordinated changes:
   - **Change A** — `build_initial_ion_state` now fills the ion mass array
     uniformly with `mass_initial_amu × U` when
     `drag_coefficients is not None AND mass_scenario == "fixed"` (else inherits
     the neutral mass bit-identically). The override site carries a §6.5 `why`
     comment.
   - **Change B** — `_check_drag_scope(cfg, initial_mass_kg)` now reads the
     **realized** initial ion mass (`ckpt.mass_kg`, downstream of the override)
     and raises if it is more than `_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU`
     (~8 amu) from `m_eff_amu` — a trip-wire against a future refactor or bypass
     leaving the mass at ~127 amu. Reads the realized array, not the config
     field, so it cannot merely echo Change A.
   - **Change C** — direct tests: override-fires (127 → `m_eff`), guard-rejects-127,
     and override-does-not-fire on non-drag / non-`fixed` configs.

   The production `single_pulse_N2000_drag → neutral → ion` path now integrates
   at `m_eff`, the precondition for a meaningful Tier-0 TDDFT comparison. The
   fix's "Risk 3" — the smoke harness's `m_eff` sidestep being a silent blind
   spot — is **retired** by Change C: the override and the realized-mass guard
   now have direct coverage, so a regression that breaks the override is caught.

2. **Tier-0 TDDFT comparison is the next task, not this slice.** It sets the
   deferred §6.10 acceptance thresholds by inspecting the first drag run against
   `9A_All_Data.csv` inside the extraction window. The mechanical smoke harness
   here is explicitly not that check.

3. **`temperature_diagnostic` all-NaN under the drag branch.** Confirm the
   post-processing plot `plot_ion_temperature_diagnostic.py` tolerates an
   all-NaN array rather than erroring. Untested here; flagged for the task that
   first reads a drag run's diagnostics.

4. **`_E_pot_per_atom` lift.** Confirmed byte-identical between the two paths
   and shared. If a future divergence appears, the rule is: duplicate
   *structure* is acceptable, duplicate *physics* is not.

5. **Per-step closure rebuild cost.** Negligible at fixed mass but retained for
   Tier-1 readiness; if profiling later flags it, the recorded trap (from
   Slice 2) is to cache the conservative *force* `F_cons`, never the
   *acceleration* `a = F/m` (which goes stale under mass dynamics).
