# Slice 4 — Goals: Ion-Driver Rewiring + O-Step Energy Accounting

**Status:** Specification. Driver wiring + energy accounting + scope guard +
debug harness contract only. No implementation code until
`[PROCEED TO IMPLEMENTATION]`.

**Scope:** Slice 4 *only* — wire the Slice 2 BAOAB stepper into the ion stage so
a drag-enabled preset runs end-to-end, deterministically, at Tier 0
(`mass_scenario=fixed`, `T_eff=0`, `linear_cubic`). It adds a parallel
`baoab_propagation_step`, a config dispatch in `ion.py`, the amu·Å²/ps²→eV
conversion of the drag dissipation, a drag-branch scope guard, and a
deterministic smoke-run harness. It consumes Slices 1–3 and produces the first
runnable drag trajectory.

**Deliverable (wiring-only — D1).** The drag preset runs end-to-end, energy
closes, the checkpoint round-trips, and the smoke harness passes. The deliverable
is *mechanical correctness*, **not** a TDDFT match. The first
`compare_distance` / `compare_velocity_magnitude` run against `9A_All_Data.csv`
is a **separate Tier-0 analysis task** that *sets* the deferred validation
thresholds (§6.10) by inspecting the result — it is not gated by a threshold
that does not exist yet. Folding it in would make "done" subjective.

---

## 1. The framing fact: this is the only slice that runs, and the riskiest

Slices 1–3 were additive (new files, new fields) and produced no runnable drag
trajectory. Slice 4 is where everything executes — and it is the slice that
**exercises the scoped CLAUDE.md "do not change collision physics" exception**,
not just declares it. The risk is contained by the central architectural choice
(§2): the existing collision path is **not rewritten in place**; it stays
bit-identical and uncalled-for-drag, so the hard-sphere regression surface — the
very thing the drag model validates against (§6 hard-sphere-variance noise
calibration depends on the old path staying runnable) — survives untouched.

---

## 2. Architecture — parallel per-step function, dispatch in `ion.py` (D3 + B)

The codebase already has **two sibling per-step functions**,
`neutral_propagation_step` (`propagation_step.py:99`) and `ion_propagation_step`
(`ion_propagation_step.py:141`), with the stage drivers picking one. Slice 4
adds a **third sibling**, `baoab_propagation_step`, in `ion_propagation_step.py`,
and `ion.py` chooses between it and `ion_propagation_step`.

```
ion.py  (stage driver)
  ├─ drag_coefficients is not None  →  build BAOAB closure, loop baoab_propagation_step
  └─ else                          →  existing collision path, loop ion_propagation_step  (UNCHANGED)
```

- **Dispatch predicate: `cfg.drag_coefficients is not None`** (D3, confirmed).
  The one condition that is both *necessary* (BAOAB cannot run without
  coefficients) and *already validated* (a non-`None` bundle has passed the
  Slice 3 guard: consistency + dissipativity + form agreement). No separate
  `use_drag` flag — that would be a second source of truth that could disagree
  with the coefficients' presence.
- **Dispatch evaluated once per run**, in `ion.py`, before the step loop — not
  re-checked per step.
- **Consequence (accepted):** "drag preset" ≡ "has coefficients." A/B scenario
  comparison uses *different presets* (drag vs. non-drag), not the same preset
  toggled — which is the cleaner design anyway.

### 2.1 `ion_propagation_step` is not touched (beyond behavior-preserving lifts)

The existing collision per-step function gets **no mode switch grafted in**. It
stays single-purpose and bit-identical, so the hard-sphere regression surface
(`test_ion_propagation_step.py`, the HeDFT comparisons) stays pristine. The only
edits to it are the small behavior-preserving *extractions* in §3.1, each
guarded by that suite staying green — the same discipline as the Slice 2
`_kick`/`_drift` lift.

### 2.2 `ion.py` owns closure construction (B)

`ion.py` assembles the `gamma_fn` and the spatial gate and builds the BAOAB
closure via `make_ion_baoab_step` — consistent with where the baseline already
rebuilds `make_ion_step` (`ion_propagation_step.py:184-193`).
`baoab_propagation_step` *receives* the ready `step` plus the conserved-energy
pieces and does per-step orchestration only (call step, convert dissipation,
book energy, populate checkpoint fields). This keeps the new per-step function
thin and keeps the gate/`gamma_fn` assembly in the orchestration layer.

**Per-step closure rebuild (confirmed — rebuild every step).** `ion.py` rebuilds
the BAOAB closure *every step*, matching the baseline `make_ion_step` rebuild
pattern, even though Tier-0 mass is fixed and a build-once would be correct
today. Rationale: consistency with the established pattern and Tier-1-readiness
(mass dynamics drop in with no restructuring) outweigh the negligible
fixed-mass cost — and it is exactly what `make_ion_baoab_step` was shaped to
support.

### 2.3 The gate assembly implements the §5.5 G4→G2 collapse

Because `ion.py` builds `gamma_fn`, it reads `drag_spatial_gate` and
`drag_gate_steepness` and assembles the gate. At Tier 0
`drag_spatial_gate=density_proportional` (the default) **collapses to the erf
complement** (§5.5: no `helium_density_profile` exists yet), *identical* to
`erf_tied`. The gate assembly must therefore treat **both `density_proportional`
and `erf_tied` as the erf-complement gate today**, with a recorded "collapses
until a density profile exists" note. Trap avoided: if the assembly handled only
`erf_tied`, the default `density_proportional` preset would fall through with no
gate built.

---

## 3. `baoab_propagation_step` — the new per-step function

Per-step orchestration, symmetric with `ion_propagation_step` but with the
*middle* replaced: no collision sampling, no `apply_collision`, no mass
attachment — just the BAOAB `step` and the dissipation booking.

| # | Action | vs. `ion_propagation_step` |
|---|---|---|
| 1 | receive rebuilt BAOAB `step` from `ion.py` | replaces the per-step `make_ion_step` build |
| 2 | `step(pos, vel, dt) → (pos', vel', E_pot_per_pair, ΔE_dissip)` | replaces leapfrog + sampling + apply + mass-attach (`:211-256`) |
| 3 | per-atom depth `r' − droplet_radii` | shared scaffolding (§3.1) |
| 4 | `E_kin` at new velocity (fixed mass), eV | shared scaffolding (§3.1) |
| 5 | `E_pot` from `step`'s returned `E_pot_per_pair`, eV | the conservative ion potential, as the baseline |
| 6 | convert `ΔE_dissip` amu·Å²/ps² → eV (§4); accumulate `E_dissip_eV` | replaces the collision `E_dissip` accumulation |
| 7 | `E_mass_attach_defect_eV` = 0; `temperature_diagnostic` = NaN | Tier-0 fills (§5) |

### 3.1 Shared scaffolding (Sub-decision A) — behavior-preserving lifts only

Extract *only* the genuinely-identical, physics-free pieces, each as a
behavior-preserving lift out of `ion_propagation_step`, **gated on
`test_ion_propagation_step.py` staying green** (Slice-2 `_kick`/`_drift`
discipline). Candidate extractions:

- **depth computation** `r − droplet_radii` — pure geometry, identical.
- **eV kinetic-energy assembly** `½ m v² → eV` — identical at fixed mass.
- **eV potential assembly** from `E_pot_per_pair` — *if* it lifts cleanly
  (both paths use the conservative ion potential; confirm the assembly is
  byte-identical before sharing).

**Rule:** if a piece cannot be proven clean by the regression suite,
`baoab_propagation_step` writes its own rather than forcing a share that
couples the two functions or risks perturbing the protected one. Duplicate
*structure* is acceptable where the *meaning* differs; duplicate *physics* is
not (rule 1) — these candidates are physics-free, so lifting them is the
rule-1-correct move when it is clean.

---

## 4. The eV conversion (4b) — reuse the baseline idiom, do not reinvent

`ΔE_dissip` returns from the stepper in **amu·Å²/ps²** (Slice 2's deliberate
pure-mechanical handoff). Slice 4 converts to eV for `E_dissip_eV`, reusing the
**exact** baseline conversion path — the mass-attach defect already does this:

```
dE_defect_eV = -0.5 * mass_diff_kg * (v_post_sq * 100**2) / EV
```

i.e. amu→kg (via `U`), Å/ps→m/s (×100), ½mv² in J, ÷`EV`. The drag conversion
must use the **same** amu→kg / ×100² / ÷`EV` path so `E_dissip_eV` is in a
consistent eV with `E_kin_eV` and `E_pot_eV`.

**Trap:** a hand-rolled factor off by the `U` (amu→kg) or the `100²` (Å/ps→m/s)
would make energy *appear* to close *within* the drag channel while drifting
against `E_kin`. Test asserts the conversion against a known `½ m v²` value in
both unit systems.

---

## 5. Checkpoint — v5 retained at Slice 4 (D4)

No schema bump. The `IonCheckpoint` v6 (rename
`E_mass_attach_defect_eV → E_mass_transfer_eV`) is tied to **mass dynamics
(Tier 1)**, not Slice 4. At Tier 0 the v5 fields take these values under the
drag branch:

| field | Tier-0 drag value | reason |
|---|---|---|
| `E_dissip_eV` | real drag dissipation (cumulative, from §4) | the live channel |
| `E_mass_attach_defect_eV` | `0` | no mass attachment at fixed mass (true, not a placeholder) |
| `temperature_diagnostic (T,3)` | `NaN` | no collision events → nothing to diagnose; documented sentinel |
| `mass_history_kg` | constant (= initial) | fixed mass |
| `number_of_collisions` | `0` | no collisions in the drag branch |

The checkpoint **round-trips** (writes v5, loads back, NaN-filled
`temperature_diagnostic` survives the `.npz` round-trip, `E_mass_attach_defect_eV`
is 0). Loader is unchanged (still enforces v5 `expected_fields ⊆ npz.files`).

---

## 6. `_check_drag_scope` — the drag-branch scope guard (D5)

The drag-branch analog of `_check_scope` (`ion_propagation_step.py:308`, which
demands `hard_sphere_collision_mode == 3` — irrelevant under drag). Asserts the
**Tier-0 envelope** so an out-of-scope drag config fails at the driver, not deep
in a half-implemented path. Mirrors how `make_ion_baoab_step` already raises
`NotImplementedError`:

```
reject  T_eff > 0                              (active noise — Slice ≥4 / Tier 3)
reject  mass_scenario != "fixed"               (mass dynamics — Tier 1)
reject  drag_form != "linear_cubic"            (other forms NotImplemented in drag.py)
reject  single_charge_ionization_allowed, additional_droplet_charges > 0,
        effusive_dynamics                      (inherited from _check_scope)
```

Runs in `ion.py` (or at the head of `baoab_propagation_step`) before the loop.
The Slice 3 `check_drag_config` validates the *config's internal consistency*;
`_check_drag_scope` validates the *Tier-0 runnability envelope* — distinct
checks (config-valid ≠ Tier-0-runnable).

---

## 7. Energy-closure invariant (the real acceptance criterion)

At Tier 0, `E_mass_transfer = 0` (no mass dynamics), so the §2.9 invariant is:

$$E_\text{kin}(t) + E_\text{pot}(t) + E_\text{dissip}(t) \approx \text{const}$$

**Closure is tight, not loose.** The O-step removes KE and books it to
`E_dissip` *analytically* (Slice 2's exact `½m(‖v_in‖²−‖v_out‖²)`), so the
dissipative part is exact and the only drift is the baseline's *conservative*
Verlet drift on `E_kin + E_pot`. The test therefore asserts closure to the
**same tolerance as the baseline conservative Verlet drift** — drag dissipation
is exact, not approximate. This makes the closure test a sharp instrument, not a
loose sanity check.

---

## 8. Smoke-run harness (the "thorough debug" artifact, D1 meta)

A **deterministic mechanical-correctness instrument**, separate from any TDDFT
comparison and *not* gated on the deferred thresholds. Config:
`single_pulse_N2000_drag`, tiny (N = 2–4 molecules, ~10–50 steps), fixed seed,
`T_eff=0`. Asserts:

- **Finite trajectories** — no NaN/Inf in positions/velocities (catches a blown
  exponent, a sign error in drag direction, a bad gate).
- **Energy closes** to the §7 tight tolerance (`E_kin+E_pot+E_dissip` conserved
  to Verlet-drift level).
- **`E_dissip_eV` monotone non-decreasing and `> 0`** (drag actually removed
  energy — catches a gate or `gamma_fn` wired to zero).
- **Checkpoint round-trips** — writes v5, loads back; `temperature_diagnostic`
  NaN-filled, `E_mass_attach_defect_eV` = 0, `number_of_collisions` = 0.
- **Scope guard rejects** an out-of-Tier-0 drag config (`T_eff>0`,
  `mass_scenario≠fixed`).

**Explicitly does NOT assert anything about matching TDDFT** — that is the
separate Tier-0 analysis task that *sets* thresholds. The smoke run proves the
machine runs *correctly*, which is the gap between "wiring-only done" and
"validated."

---

## 9. Acceptance criteria

| Quantity | Expected | Tolerance |
|---|---|---|
| Dispatch — `drag_coefficients is not None` routes to `baoab_propagation_step`; `None` routes to `ion_propagation_step` | enforced | exact |
| Hard-sphere path bit-identical — a non-drag preset run is byte-identical to pre-Slice-4 | no regression | exact / round-off |
| `test_ion_propagation_step.py` green after the §3.1 lifts | behavior-preserving | round-off |
| Energy closure (§7) — `E_kin+E_pot+E_dissip` conserved over a run | baseline Verlet-drift tolerance | tight |
| eV conversion (§4) — drag `ΔE_dissip` → eV matches a known ½mv² in both unit systems | exact | rtol ~1e-12 |
| Gate collapse (§2.3) — `density_proportional` and `erf_tied` build the identical erf-complement gate | identical | exact |
| Per-step rebuild — BAOAB closure rebuilt each step; result identical to a build-once reference at fixed mass | identical | round-off |
| Checkpoint v5 round-trip (§5) — NaN `temperature_diagnostic`, 0 defect, 0 collisions survive write/load | enforced | exact |
| Scope guard (§6) — `T_eff>0` / `mass_scenario≠fixed` / `drag_form≠linear_cubic` rejected | raises | exact |
| Smoke harness (§8) — all five assertions | pass | per §8 |

---

## 10. Open items surfaced (recorded, not blockers)

- **Tier-0 TDDFT comparison is the *next* task, not this slice.** It sets the
  deferred §6.10 thresholds by inspecting the first drag run against
  `9A_All_Data.csv` inside the extraction window. Recorded as the immediate
  successor to Slice 4.
- **§3.1 potential-assembly lift** — confirm the eV `E_pot` assembly is
  byte-identical between the two paths before sharing; if not, `baoab` writes
  its own. Decide at implementation against the regression suite.
- **`temperature_diagnostic` NaN vs. sentinel** — NaN chosen (D4); confirm the
  post-processing plots (`plot_ion_temperature_diagnostic.py`) tolerate an
  all-NaN array under the drag branch rather than erroring. Flag for the
  Tier-0-comparison task that first reads these.

---

## 11. Scope fence — what Slice 4 does NOT touch

- **No in-place rewrite of `ion_propagation_step`** beyond the behavior-
  preserving §3.1 lifts — it stays the single-purpose collision step.
- **No checkpoint schema bump.** v5 retained; v6 comes with mass dynamics
  (Tier 1).
- **No mass dynamics.** Mass fixed; the mass-attach block stays in the
  uncalled-for-drag collision path (D2 — not deleted, not dormant-in-new-path).
- **No active noise.** `T_eff=0`; the scope guard rejects otherwise.
- **No TDDFT validation / no threshold setting.** Wiring-only; the comparison is
  the next task.
- **No new `SimConfig` fields.** Slice 4 *reads* the Slice 3 fields; it adds
  none.
- **No `physics/` change.** `drag.py` / `baoab.py` are consumed unchanged; the
  only edits are in `ion_propagation_step.py` (new function + §3.1 lifts) and
  `ion.py` (dispatch + closure construction).
- **Neutral stage untouched.**

---

## 12. Definition of done

- `baoab_propagation_step` added to `ion_propagation_step.py`; the existing
  `ion_propagation_step` bit-identical beyond the §3.1 behavior-preserving lifts.
- `ion.py` dispatches on `drag_coefficients is not None`, builds the
  `gamma_fn`/gate (with the §2.3 collapse) and the BAOAB closure, rebuilds per
  step, loops `baoab_propagation_step`.
- `ΔE_dissip` converted to eV via the baseline idiom (§4) and accumulated into
  `E_dissip_eV`.
- Checkpoint v5 retained with the §5 Tier-0 fills; round-trips.
- `_check_drag_scope` rejects the out-of-Tier-0 envelope (§6).
- The smoke harness (§8) passes; energy closes to baseline-Verlet tolerance.
- Hard-sphere path and `test_ion_propagation_step.py` regression-green.
- The §10 open items recorded in-repo; the Tier-0 TDDFT comparison flagged as
  the immediate next task.
