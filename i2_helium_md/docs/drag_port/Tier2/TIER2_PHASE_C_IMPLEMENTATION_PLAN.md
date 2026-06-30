# Tier 2 — Phase C Implementation Plan (Schema + Generative Integrator)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no pseudo-code until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations from
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (LaTeX + dimensional checks);
> module descriptions are *interface contracts*, not implementations.
>
> **Entry docs:** `TIER2_IMPLEMENTATION_PLAN.md` §3/§4 (Phase C in the dependency
> graph); `TIER2_PHASE_A_IMPLEMENTATION_PLAN.md` (L/K/U energetics) and
> `TIER2_PHASE_B_IMPLEMENTATION_PLAN.md` (ρ/P/Q channels) — the layers Phase C
> *composes*; `MASS_DYNAMICS_LOCKED_*.md` §6 (the budget S1/S2/K1/K2 + the **5-term
> invariant**, eq. line 935), §7 (schema implications), A13 (jump operators,
> one-event-per-step, jump-then-O), §8 R4 (the double-count guard). **Reuse spine:**
> `simulation/checkpoint.py` (v6 + the v5→v6 shim), `postprocess/energy_balance.py`
> (`ion_energy_totals`, `ion_ledger_closure` — the 4-term gate),
> `simulation/ion_propagation_step.py` (`IonStepState`, the `shed_step` jump-then-O
> seam, `baoab_propagation_step`), `physics/baoab.py` (the SQ1 O-step).

---

## 0. Status and intent

Tier 0 / Tier 1a are delivered; **Phase A** (energetics primitives L/K/U) and **Phase
B** (stochastic channels ρ/P/Q) are delivered *as plans*. **Phase C composes A + B into
the production loop**: it persists the new `E_int` reservoir (Slice **X** — checkpoint
**v6→v7** + the **5-term invariant** gate) and assembles the `biphasic` generative
integrator (Slice **G** — per-step cooling + one stochastic mass event + jump-then-O +
`E_int` evolution). This is the first slice where the channels stop being isolated
primitives and become a *trajectory*.

Phase C asks **two** questions: *how is the internal-energy reservoir persisted and how
is closure extended from four terms to five (X); and how do the per-step conservative
kicks, Newton cooling, the at-most-one pickup/evaporation event, the jump-then-O seam,
and the `E_int` update compose into one deterministic-given-RNG step (G)?* X maps onto
the **energy-bookkeeping** validation step (the §2.9→§6 invariant closes); G maps onto
the **full-driver smoke** step — but always behind the invariant as the cross-cutting
correctness gate.

**Central engineering fact.** Phase C is **composition, not new physics**: the only new
physics already landed in A/B (the rungs, cooling, channels, resets). G is a per-step
*orchestration* that wires accepted modules together; X is a *schema + closure* extension
that mirrors the v5→v6 / 4-term precedent exactly. The bar here is **the 5-term invariant
closing to Verlet drift** on analytic limits *before* any real run, plus the `fixed`
regression staying untouched.

**Locked decisions (2026-06-29, this session).**
- **Slice G placement — extend `simulation/ion_propagation_step.py`** with a
  `biphasic_step` seam alongside the existing `shed_step` / `baoab_propagation_step`,
  reusing the `IonStepState` dataclass (extended with `E_int_eV`) and the SQ1 O-step in
  place. Mirrors the Tier-1a `shed_step` precedent.
- **Production scenario name — keep `biphasic`.** The `MassScenario` literal stays
  `{fixed, biphasic, anchored_discrete}`; `biphasic` *is* the production value (full name
  `biphasic_energy_gated`, MASS §11, recorded in docs/log). No config-enum rename, no
  guard-set churn.
- **Schema bump — v6→v7** (adds `E_int_eV (2N,T)`), with a **v6→v7 back-compat shim**
  (E_int absent → synthesize zeros + flag) mirroring the delivered v5→v6 shim.
- **5-term closure — extend in place.** Add `E_int_eV` to `EnergyTotals` /
  `ion_energy_totals` / `ion_ledger_closure` as an optional term (exactly how
  `E_mass_transfer_eV` was added at v6), not a parallel function.

---

## 1. Scope lock — what Phase C does and does not do

**Does:**
- **Slice X:** the `IonCheckpoint` **v6→v7** bump (new `E_int_eV (2N,T)`, `mass_scenario`
  tag `biphasic`), the **v6→v7 migration shim**, the extension of `IonStepState` with the
  `E_int_eV` per-atom field, and the **5-term invariant closure**
  `E_kin + E_pot + E_dissip + E_mass_transfer + E_int ≈ const` in `energy_balance.py`;
  the **RNG draw-order lock** (forbidden-list item) is *specified and frozen* here.
- **Slice G:** the `biphasic_step` generative integrator in
  `simulation/ion_propagation_step.py` — per step: conservative kicks → **K2 Newton
  cooling** of `E_solv.struct` → **at most one mass event** (shed via Q else pickup via
  P; **shed-then-pickup** if both Bernoulli draws fire) applied **jump-then-O** (reuse the
  `shed_step`/capture seam + the SQ1 O-step rebuilt at post-jump `m⁺`) → **E_int update**
  (S1 / S2-onset / K1 / K2 via U); the config wiring (`mass_scenario=biphasic` enters the
  §6.5 `time_resolved` arm; the `s≥1` guard; noise inert).

**Does NOT (fails review if it leaks in):**
- any **new channel physics** — pickup/evaporation/cooling/ladder live in A/B and are
  consumed, not re-derived here;
- the **bridge / comparison / calibration** (Phases D/E/F): no generative-vs-anchored run,
  no size-distribution loader, no Wasserstein, no campaign;
- the **noise model** (Tier 3 — stays inert behind its enum; G reads no noise amplitude);
- any change to **neutral propagation**, the **Tier-0-locked drag law**, or the
  `IonCheckpoint` field *semantics* beyond the additive v6→v7 bump.

**Double-count guard (R4, load-bearing).** Drag work → `E_dissip` *only*; binding release
on pickup splits **`+f_ret·D_0` → `E_int`** and **`(1−f_ret)·D_0` → `E_dissip`** (bath);
cold-shed binding return books to `E_pot`/`E_dissip` per the K1 neutrality (§6 K2). G must
never let the same joule enter two reservoirs — this is the chief energy-balance risk and
the thing the 5-term closure is built to catch.

---

## 2. Locked physics — the invariant and schema (the test oracles)

### 2.1 The 5-term invariant (MASS §6, eq. line 935)

$$
E_\text{kin}+E_\text{pot}+E_\text{dissip}+E_\text{mass\_transfer}+E_\text{int}\approx\text{const}
\quad(\text{closes to Verlet drift}).
$$

It extends the delivered Tier-1a **4-term** invariant
(`E_kin+E_pot+E_dissip+E_mass_transfer`, `energy_balance.ion_ledger_closure`) by the new
`E_int` reservoir. **Per-channel closure (each sums to zero, MASS §6 audit):**

| Channel | Bookings | Closes because |
|---|---|---|
| drag (continuous) | `E_kin ↓`, `E_dissip ↑` | exact SQ1 dissipation defect (Tier-0) |
| pickup S1 | `E_int += f_ret·D_0`, `E_dissip += (1−f_ret)·D_0`, `E_mass_transfer +=` capture defect | `−D_0 + f_ret·D_0 + (1−f_ret)·D_0 = 0` + reduced-mass KE defect |
| cold-shed K1 | `E_int −= D_0`, `E_pot += D_0`, A8 marginal → bath | binding returned as the rung leaves |
| cooling K2 | `dE_int = −dE_dissip` at fixed N | Newton drain books to `E_dissip` |

**Reduced-mass capture defect (the +He counterpart booked into `E_mass_transfer`, A13):**
$$
\Delta E_\text{cap}=\tfrac12\frac{m\,m_\text{He}}{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2
\xrightarrow{u_\text{He}=0}\tfrac12\frac{m\,m_\text{He}}{m+m_\text{He}}\|v^-\|^2\ (>0),
$$
*not* the heavy-ion `½·m_He·v²` (closes **exactly**, not to ~3%, under the §4 momentum
reset). The cold-shed defect (negative) is the delivered `mass_jump.cold_shed` form.
*Dim:* all terms eV (after the `amu·Å²/ps² → eV` conversion `mass_jump`/`baoab` already
own). ✓

### 2.2 Schema delta (MASS §7; `checkpoint.py` v6 → v7)

- **ADD** `E_int_eV (2N, T)` — the internal-energy trajectory (per-atom, eV).
- **TAG** `mass_scenario = "biphasic"` (production); the field already exists at v6.
- **KEEP** every v6 field (`E_mass_transfer_eV`, `n_shell`, `mass_history_kg` with its
  dropped-monotonicity note, …). `n_shell` becomes a **genuine state** under `biphasic`
  (pickup/evap change it directly), not a quantity derived from `mass_history_kg` — the
  writer sets it from the channel, and the v6→v7 shim leaves the v6 derivation intact.
- **Migration v6→v7:** `E_int_eV` absent → synthesize a zeros array (+ a load-time flag
  that the file predates the reservoir), mirroring the v5→v6 `_migrate_ion_checkpoint`
  pattern (`checkpoint.py:312`). v5 still migrates v5→v6→v7 transitively or is rejected —
  decide at build, but the **v6→v7 shim is mandatory** so existing Tier-1a `ion.npz`
  load.
- `_ION_SCHEMA_VERSION` 6 → 7; add `E_int_eV` to the `trajectory_2N_T_fields` shape-check
  list in `_validate_against_cfg` (`checkpoint.py:499`).

### 2.3 Step composition (MASS §4 + A13; the G contract)

Per BAOAB sub-step `dt = dt_ion = 0.01 ps`, in **fixed order**:
1. conservative **B/A** kicks (positions/velocities under the conservative force);
2. **K2 cooling** of `E_solv.struct` toward `E_∞(N)` (closed-form `e^{−dt/τ}`, Slice K);
3. **at most one mass event** (A13): draw shed (Q) and pickup (P) Bernoullis; if **both
   fire, apply shed then pickup** (fixed order, joint prob. `~λν dt² ~ 4e−4`);
4. apply the event **jump-then-O**: the jump (capture/`cold_shed` reset, `m⁺`) at the seam,
   then the **SQ1 O-step rebuilt at `m⁺`** (friction + the inert noise amplitude both read
   `m⁺`/`γ(v⁺)`, SQ3);
5. **E_int update**: S2 onset deposit at `t=0`; S1 heat on pickup; K1 drain on shed; K2
   drain folded continuously (the §6 jump-consistency keeps S1/K1 unchanged under the
   variable swap).

*Dim/limit checks:* force-free + no events ⇒ the step reduces to the Tier-1a path and the
5-term sum is constant; one isolated pickup/shed ⇒ the matching `E_int`/`E_mass_transfer`
increments make the residual flat to machine precision. ✓

---

## 3. New-work slices

Each slice carries the template: **Purpose / Module + interface / Encoded form / Knobs /
Oracle / Independence / Test spec / Acceptance**.

### Slice X — Checkpoint v6→v7 + 5-term invariant gate *(schema + closure; composes nothing stochastic)*

**Modules.** `simulation/checkpoint.py` (schema), `postprocess/energy_balance.py`
(closure), `simulation/ion_propagation_step.py` (`IonStepState` field add).

**Purpose.** Persist the `E_int` reservoir and extend the correctness gate from four terms
to five — the detector that catches a G miswire from S onward.

**Interface (proposed).**
- *checkpoint:* add `E_int_eV: np.ndarray` to `IonCheckpoint`; `_ION_SCHEMA_VERSION = 7`;
  extend `_migrate_ion_checkpoint` with the **v6→v7** arm (synthesize+flag); add `E_int_eV`
  to the `trajectory_2N_T_fields` validation tuple.
- *energy_balance:* add `E_int_eV` to `EnergyTotals` (optional, `None` for neutral / a
  v6-origin run); `ion_energy_totals` sums it into `E_system_eV`; `ion_ledger_closure`
  then reports the **5-term** residual unchanged in shape (reuses the same machinery —
  this is the "extend in place" decision).
- *step state:* add `E_int_eV: np.ndarray (2N,)` to `IonStepState` (the per-step carrier G
  reads/writes; defaults to zeros so the `fixed`/Tier-1a paths are untouched).

**Encoded form.** §2.1 (invariant) + §2.2 (schema delta).

**Knobs (config §4):** none new here — X is schema/closure; the `mass_scenario=biphasic`
tag is read from the existing field. The **RNG draw-order lock** is *recorded* here as a
forbidden-list item (pickup-vs-shed sequence fixed; documented in `DRAG_PORT_DESIGN_
DECISIONS.md` and the log).

**Oracle values.**
- v7 round-trips (`save`→`load`) bit-for-bit on a tiny synthetic checkpoint.
- v6 and v5 still load (shim) — `E_int_eV` synthesized zeros + flag for v6-origin files.
- 5-term closure: on a synthetic run with a known `E_int` trajectory, `residual_eV` =
  `E_system[t] − E_system[0]` is zero to machine precision when the channels balance, and
  **diverges** under a deliberately-miswired `E_int` increment (relabel-fault still caught,
  now in five terms).

**Independence.** Pure schema + arithmetic on stored arrays; no RNG, no integrator. Built
before G so G has its gate ready.

**Test spec (`tests/test_checkpoint.py`, `tests/test_energy_balance.py`).**
- v7 round-trip; v6→v7 and v5→…→v7 migration (synthesize+flag); missing-`E_int_eV` on a
  genuine v7 raises (not silently zero-filled).
- `_validate_against_cfg` rejects a wrong-shape `E_int_eV`.
- `ion_ledger_closure` 5-term: flat residual on a balanced synthetic stream; divergent on a
  miswire; neutral/`E_int=None` path unaffected.

**Acceptance.** v7 round-trips; v6/v5 still load; 5-term closure reuses `ion_energy_totals`;
relabel-fault still caught; `fixed` runs produce an all-zero `E_int_eV` and the residual is
unchanged from the 4-term baseline.

---

### Slice G — Generative driver `biphasic_step` *(composes A + B; the production loop)*

**Module.** `simulation/ion_propagation_step.py` (new `biphasic_step` seam).

**Purpose.** The production `biphasic` integrator: assemble the accepted A/B modules into
one per-step loop with the **one-event-per-step**, **jump-then-O**, **5-term-closing**
contract.

**Interface (proposed).**
- `biphasic_step(state, *, rng, cfg, droplet_radii, conservative_force_fn, ...) ->
  IonStepState` — runs §2.3 steps 1–5 and returns the advanced `IonStepState` (with updated
  `vx/vy/vz`, `mass_kg`, `n_shell`-carrying state, `E_int_eV`, `E_mass_transfer_eV`,
  `E_dissip_eV`, `time_ps`). Composes:
  - `physics/internal_energy_cooling.newton_cool_step` (K2),
  - `physics/evaporation.evaporation_step` (Q) and `physics/pickup.pickup_step` (P) under
    the **shed-then-pickup** order, each calling the `mass_jump` reset,
  - the existing `shed_step`/capture **seam** + `make_ion_baoab_step` **SQ1 O-step rebuilt
    at `m⁺`** (SQ3),
  - `physics/internal_energy_budget` (U) for the S1/S2/K1 `E_int` increments.
- Reuse the existing driver loop (the function that calls `baoab_propagation_step` /
  `shed_step` per stored step) — dispatch to `biphasic_step` when
  `cfg.mass_scenario == "biphasic"`.

**Encoded form.** §2.3 (step composition); the channel forms are MASS §4 (Phase B).

**Knobs (config §4):** consumes the full Phase-A/B surface (κ, picture, |S|, τ, f_int,
f_ret, λ_0, p, ν, s, the cap/form/reset enums). New **integrator** flags:
`one_mass_event_per_step=true`, `jump_o_step_ordering=jump_then_O`,
`mass_jump_velocity_reset=momentum_conserving` (forbid `label_only` in production),
`he_capture_velocity=at_rest`.

**Config guards (§6.5 / A11 / Tier 3).**
- `mass_scenario=biphasic` is already in `_EVOLVING_MASS_SCENARIOS` → requires
  `time_resolved` coefficients → the production pairing **trips the §6.5 guard
  structurally** → runs under **`allow_inconsistent_mass_pairing=True`** (the §6.6 mid-window
  `m≈19 He = m_eff` defence; warns, does not raise).
- `evap_rrk_dof` **`s≥1`** load guard (Phase B) is enforced before G runs.
- **Noise stays inert** (Tier 3): G reads no FDT noise amplitude; the O-step noise term is
  zero behind its enum.

**Oracle values / analytic limits (checked before any real run).**
- **`fixed` regression unaffected:** with `mass_scenario=fixed` the dispatch never enters
  `biphasic_step`; the Tier-0/1a trajectories are bit-identical.
- **Force-free, no-event run:** the 5-term invariant is constant to machine precision.
- **One isolated pickup / one isolated shed:** the `E_int`/`E_mass_transfer`/`E_dissip`
  increments make the residual flat (the per-channel closure of §2.1).
- **One-event-per-step honored:** even with both Bernoullis forced to fire, exactly one
  shed **and** one pickup are applied in the fixed shed-then-pickup order, `n` net unchanged,
  bookkeeping closed.
- **Draw order honored + reproducible:** a fixed seed reproduces the trajectory; swapping the
  documented draw order changes results (proving the lock matters).
- **E_int finite + gate correct:** `E_int ≥ 0`, no evaporation while `E_int > Σ(n)`, and the
  gate opens at `t_×` as cooling drains `E_int` below the integrated ladder.

**Independence.** Composes only **accepted** modules (A delivered, B delivered, Tier-0/1a
seam delivered). For its own tests it can mock the conservative force and inject a seeded
RNG; the channel internals are already unit-tested in A/B, so G's suite asserts the
**composition** (order, one-event, closure), not the channel physics again.

**Test spec (`tests/test_ion_propagation_step.py` / a focused `test_biphasic_step.py`).**
- `fixed` dispatch unchanged (regression); `biphasic` dispatch enters the new seam.
- Force-free 5-term closure to machine precision over many steps (mock force = 0).
- Isolated-event closure (one pickup, one shed) and the both-fire one-event-per-step order.
- Seeded reproducibility; draw-order sensitivity.
- `E_int` non-negativity + gate behaviour on a short synthetic cooling ramp.
- The §6.5 guard trips and runs under `allow_inconsistent_mass_pairing=True`; the `s≥1`
  guard blocks a bad override before the loop.

**Acceptance.** `fixed` regression bit-identical; force-free run conserves the 5-term
invariant; one-event-per-step + locked draw order honored; `E_int` finite with correct gate
behavior; the pairing guard + `s≥1` guard fire as specified; **no noise amplitude read**.

---

## 4. Config contract (fields landing / activated in Phase C — MASS §11)

Most knobs are declared in Phase A/B (rule-2 exception); Phase C **activates** them by
wiring G, and adds the integrator flags. Validators reuse `check_drag_config`.

- **Integrator (G):** `one_mass_event_per_step` (default `true`), `jump_o_step_ordering`
  (fixed `jump_then_O`), `mass_jump_velocity_reset ∈ {momentum_conserving(default),
  label_only}` (forbid `label_only` in production), `he_capture_velocity ∈ {at_rest(default),
  thermal}` (thermal deferred — Tier 3 noise coupling).
- **Scenario:** `mass_scenario` keeps the literal **`biphasic`** as the production value
  (full name `biphasic_energy_gated`, MASS §11 — documented, not renamed). It already sits
  in `_EVOLVING_MASS_SCENARIOS`, so the §6.5 `time_resolved` requirement + the
  `allow_inconsistent_mass_pairing` override apply unchanged.
- **Schema (X):** no config field; `_ION_SCHEMA_VERSION = 7`.
- **RNG draw order:** *frozen* at X (forbidden-list: changing pickup-vs-shed draw order
  needs explicit approval), recorded in `DRAG_PORT_DESIGN_DECISIONS.md` + the log.

**Naming note (build-time).** The `mass_rate_*`→`pickup_*` rename flagged in Phase B lands
with P; Phase C only consumes the renamed fields. Record the final names in
`drag_migration_log_tier2.md`.

---

## 5. Dependency / build order (within Phase C)

```
X (checkpoint v6->v7 + IonStepState.E_int + 5-term closure)   ── built & tested first
        │   (gate ready before the driver exists)
        ▼
G (biphasic_step: composes A {K,U} + B {P,Q,ρ} + Tier-1a {seam, SQ1 O-step})
```

X is independent of G (pure schema + closure) and **must precede** G so the 5-term invariant
exists as G's acceptance gate. G composes only **accepted** modules — every Phase-A/B slice
and the Tier-0/1a seam — and adds no new physics. Nothing in Phase C depends on Phases D/E/F.

---

## 6. Testing methodology

- **Schema (X):** round-trip + migration on **tiny synthetic** checkpoints; shape-validation
  rejects; **no production-sized checkpoints** in pytest.
- **Closure (X):** 5-term residual on synthetic streams — flat when balanced (machine
  precision), divergent on a deliberate miswire.
- **Driver (G):** analytic limits **before** any real run — force-free conservation,
  isolated-event closure, one-event order, seeded reproducibility, draw-order sensitivity,
  gate behaviour on a short cooling ramp — with a **mocked conservative force** and seeded
  RNG. The full-driver smoke (a few-step real run) is the last check, not the first.
- **No figures, no production-sized runs, no noise** in pytest. Justify tolerances (machine
  precision for closure identities; sample-size bands only where a channel draw is involved).

Suites: `tests/test_checkpoint.py`, `tests/test_energy_balance.py`,
`tests/test_ion_propagation_step.py` (+ a focused `test_biphasic_step.py`). Run the narrowest
first, then the full suite. Interpreter (Python may not be on PATH):

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

---

## 7. Acceptance criteria

**Slice X:** v7 round-trips; v6/v5 still load via the shim; the 5-term closure reuses
`ion_energy_totals` and catches a relabel-fault in five terms; `fixed` runs carry an all-zero
`E_int_eV` and leave the 4-term residual unchanged.

**Slice G:** the `fixed` regression is bit-identical; a force-free run conserves the 5-term
invariant to Verlet drift; one-event-per-step + the locked shed-then-pickup draw order are
honored; `E_int` stays finite with correct gate behavior; the §6.5 pairing guard and the
`s≥1` guard fire as specified.

**Phase C overall:** full suite stays green (current baseline all-green on
`drag_implementation`); the invariant is the cross-cutting gate from X onward.

**Out-of-scope guard (fails review if it leaks in):** any Phase-C code path that introduces
new channel physics, reads a **noise amplitude**, changes the **Tier-0 drag law**, changes
the **RNG draw order** without explicit approval, runs a bridge/comparison/calibration step
(D/E/F), or touches neutral propagation. Those are later phases / Tier 3.

---

## 8. Risks / notes carried

- **R4 double-count is the chief energy risk** (§6 / §8). Drag→`E_dissip` only; pickup
  binding release splits `E_int` (`f_ret`) / `E_dissip` (remainder); cold-shed neutrality
  keeps K1∩K2 disjoint. The 5-term closure is built to catch a violation — treat a divergent
  residual as a miswire, never as "loosen the tolerance."
- **`E_int` is a *constructed* reservoir (R2, HIGH).** It has no native MD source; G builds it
  from S2/S1/K1/K2. A wrong construction mistimes evaporation — but the cascade is bounded by
  the `D_0` ladder regardless of partition (energetics cap the shed count), so errors are
  bounded, not divergent; the partition is pinned at Tier 2 (Phase F), not here.
- **MASS §7 says "v6"; the live schema is already v6** (Tier-1a used the slot). Tier 2 bumps
  **v6→v7** to add `E_int_eV` — the §7 text predates the Tier-1a v6 use; the v6→v7 shim is the
  reconcile. Record this in the log so the doc/code version story stays straight.
- **`biphasic` vs `biphasic_energy_gated` — resolved (keep `biphasic`).** Config literal
  unchanged; the full MASS §11 name is documentation only. No guard-set or test churn.
- **RNG draw order: locked here.** Phase B *specified* the shed-then-pickup order; Slice X
  *freezes* it as a forbidden-list item so results are reproducible and reviewable.
- **Mechanism locked, values open** — G wires forms and exposes knobs; no fitting in Phase C
  (that is the Phase F campaign).
