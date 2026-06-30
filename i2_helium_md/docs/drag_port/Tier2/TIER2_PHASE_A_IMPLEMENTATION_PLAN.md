# Tier 2 — Phase A Implementation Plan (Energetics Primitives)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no pseudo-code until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations from
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (LaTeX + dimensional checks);
> module descriptions are *interface contracts*, not implementations.
>
> **Entry docs:** `TIER2_IMPLEMENTATION_PLAN.md` §3/§4 (Phase A in the dependency graph);
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §4 (mechanism), §6 (budget
> S1/S2/K1/K2 + invariant), §8 R3/A5 (Form U ladder), §11 (config); `CALIBRATION_MAP.md`
> (parameter classes). Predecessor / structural template: `../Tier1/TIER1A_IMPLEMENTATION_PLAN.md`
> and the delivered `physics/shell_schedule.py`.

---

## 0. Status and intent

Tier 0 (drag form locked: `shared_pure_cubic`, `γ=g·b·v²`) and Tier 1a (anchored
kinematic mass-dynamics validation) are **delivered**. Tier 2 builds the *generative*
`biphasic_energy_gated` mass mechanism. **Phase A is its foundation layer**: the three
**pure, stateless energetics primitives** that every later slice consumes.

Phase A asks **one** question per slice and answers it with a closed-form function
checked against an oracle: *what is the per-rung dissociation cost `D_0(n)` and the
self-bound gate threshold (L); how does the solvation-structure energy cool and where does
it cool to (K); and what are the per-channel `E_int` heating/draining increments (U)?*
These map onto **validation steps 1–2** (direct formula + shape/unit) of the
scientific-code-caution order — the cheapest, most-entangled-free checks, done first.

**Central engineering fact.** Phase A is *all new* but *structurally trivial*: no state,
no RNG, no integrator. It is the Tier-2 analog of Slice S (`shell_schedule.py`) — the
"first, fully-independent build unit." Everything stateful (channels P/Q, Phase B),
composed (integrator G, Phase C), and observable (H/W, Phase E) **consumes** Phase A and
lives later. The correctness bar here is *exact agreement with the MASS-doc oracles*, not
trajectory fidelity.

**Locked decisions (2026-06-29).**
- **Module layout:** three separate one-concern modules (mirrors the Tier-1a
  `shell_schedule.py` / `mass_jump.py` precedent and the L/K/U independence in the Tier-2
  dependency graph):
  - `physics/dissociation_ladder.py` — Slice L,
  - `physics/internal_energy_cooling.py` — Slice K,
  - `physics/internal_energy_budget.py` — Slice U.
- New sourced constants land in `physics/constants.py`; new knobs in `config.py`.

---

## 1. Scope lock — what Phase A does and does not do

**Does:**
- the Form U dissociation ladder `D_0(n)` and its cumulative self-bound **gate threshold**
  `Σ(n) = Σ_{i≤n} D_0(i)` (L);
- the Newton cooling **driver** `dE_solv.struct/dt` and its occupancy-resolved asymptote
  `E_∞(N) = −|S(N)|` with the **pair + electrostriction binding split** (K);
- the per-channel **`E_int` budget rules** — S2 onset deposit, S1 pickup heating, K1 shed
  drain — plus the **post-t× reconstruction** identity (U).

**Does NOT (fails review if it leaks in):**
- any **RNG / Bernoulli / Poisson** draw — that is the pickup/evaporation *channels*
  (P, Q; Phase B);
- any **integrator / O-step** coupling, jump-then-O seam, or `m(t)` plumbing — the
  generative driver (G; Phase C);
- any **`E_int` state persistence**, checkpoint schema bump, or the **5-term invariant
  closure** — that is the schema slice (X; Phase C);
- the **noise model** (Tier 3 — stays inert behind its enum);
- reading the **Tier-0 drag law** or any `γ`/force-coefficient quantity.

**Mass-agnostic guard.** No Phase-A function takes the ion mass `m(t)` or a velocity. The
only place a velocity-defect *form* is named is U's documentation of the S1/K1 sign
convention (the reduced-mass capture defect is **booked by P at Phase B**, not computed
here). This preserves the unified friction convention (CLAUDE.md "Friction convention"):
the drag-physics module stays mass-agnostic and Phase A never touches it.

---

## 2. Locked physics — constants and equations (the test oracles)

All numbers below are sourced from the MASS doc and verified by direct conversion
(`EV_PER_WAVENUMBER = 1/8065.543937 = 1.23984e-4 eV/cm⁻¹`, already in `constants.py`).
They are the oracles the Phase-A pytest suites assert against.

### 2.1 Constants (land in `physics/constants.py`)

| Symbol | Value | Units | Source / provenance |
|---|---|---|---|
| `D0_1_X2` | 106.9  ( = **0.013254 eV** = 13.3 meV) | cm⁻¹ | IHe05 EPAPS exact J=0 ZPE; MASS R3 (pinned ±3 cm⁻¹) |
| `D0_1_MIX` | 74.4  ( = **0.009224 eV** = 9.23 meV) | cm⁻¹ | IHe05 statistical SO mixture (X₂+I₁+I₀)/3; MASS A10 |
| `D0_1_COOLING_RELAXED` | ∈ (74.4, 106.9) | cm⁻¹ | relaxation-weighted blend; MASS rev 2026-06-21 #3 |
| `D_FLOOR` | 4.97  ( = **6.16e-4 eV** = 7.15 K) | cm⁻¹ | bulk-He chemical potential \|μ_He^bulk\|; MASS R3 |
| `N_STAR` | 21 | count | I2-notes first-shell cation; MASS OQ8 (`21⁺⁰₋₁`) |
| `NU_EVAP_PER_PS` | 2.42 | ps⁻¹ | IHe05 X₂ well curvature ω_e=80.6 cm⁻¹; MASS A11 (read by Slice Q, Phase B) |
| `S_ABS_eV` (\|S\|) | 0.308  ( = 2484 cm⁻¹; \|S\|/n* = **118 cm⁻¹** = 170 K @ n*=21) | eV | DFT first-shell solvation S_I⁺; MASS K2 |
| `E_BIND_DRAG_eV` | 0.1168  ( = 942 cm⁻¹) | eV | §6.5.1 drag effective binding (integrated cross-check, *not* a rung-sum target) |

Scenario-keyed (already present as `coulomb_available_eV`, `config.py:210`):

| `E_avail^ion` | Value | Regime |
|---|---|---|
| validation (d = 9 Å, ½·14.40/9) | **0.80 eV** | the budget the drag + 21→19→14 shell refs were generated under |
| production (Rₑ = 2.666 Å, ½·14.40/2.666) | **2.70 eV** | vertical double-ionization onset (3.375× hotter) |

Cooling time band (read by Slice K):
`τ_dissip ∈ [2.6, 16.5] ps` — MASS R8 sweep band (Na⁺ experiment 2.6±0.4; TDDFT 7.3–16.5),
a ±factor-3 prior, **not** a fit target in Phase A.

### 2.2 Encoded forms (LaTeX + dimensional check)

**(L) Form U ladder + cumulative gate.**
$$
D_0^{\,\mathrm{I^+}}(n)=D_\text{floor}+\big(D_0(1)-D_\text{floor}\big)\frac{1-\sigma(n)}{1-\sigma(1)},
\qquad
\sigma(n)=\big[1+e^{-\kappa(n-n^*-\tfrac12)}\big]^{-1},
\qquad
\Sigma(n)=\sum_{i=1}^{n}D_0(i).
$$
*Dim:* `σ` dimensionless (κ per-unit-`n`, `(n−n*−½)` dimensionless); bracket dimensionless;
`D_0(n)`, `Σ(n)` in eV. ✓ Cliff centered at `n*+½ = 21.5` (between last in-shell atom `n*`
and first shell-2 atom).

**(K) Newton cooling + occupancy-resolved asymptote + binding split.**
$$
\left.\frac{dE_\text{solv.struct}}{dt}\right|_{K2}=-\frac{E_\text{solv.struct}-E_\infty(N)}{\tau_\text{dissip}},
\qquad
E_\infty(N)=-|S(N)|,\qquad |S(N)|=|S|\,\frac{\Sigma(N)}{\Sigma(n^*)},
$$
$$
E_\text{solv.struct}(N)=\underbrace{-\Sigma(N)}_{E_\text{bind}^\text{pair}(N)}
+\underbrace{-\big(|S(N)|-\Sigma(N)\big)}_{E_\text{elec}(N)\le 0}
+\,E_\text{int}(N).
$$
*Dim:* every term eV; `(E_solv.struct − E_∞)/τ` in eV/ps = power. ✓ `E_∞(N) → 0` as
`N → 0` (occupancy-resolved → total strip reachable, OQ6).

**(U) `E_int` budget rules.**
$$
\text{S2 onset: } E_\text{int}(0)=f_\text{int}\,E_\text{avail}^\text{ion};\qquad
\text{S1 pickup: } \Delta E_\text{int}=+f_\text{ret}\,D_0(n{+}1)\ \ (\text{bath gets }(1-f_\text{ret})D_0);
$$
$$
\text{K1 shed: } \Delta E_\text{int}=-D_0(n);\qquad
\text{reconstruction (post-}t_\times\text{ only): } E_\text{int}=E_\text{solv.struct}-E_\text{bind}^\text{pair}(N)-E_\text{elec}(N).
$$
*Dim:* `f_int`, `f_ret` dimensionless; all energies eV. ✓ All Phase-A quantities are
energies (eV) or a rate (eV/ps); **no force-coefficient `γ` appears anywhere in Phase A.**

---

## 3. New-work slices

Each slice carries the Tier-1a template: **Purpose / Module + interface / Encoded form /
Knobs / Oracle / Independence / Test spec / Acceptance**. Interfaces are *proposed
contracts* (final names settle at build), vectorized scalar-in/array-in following
`shell_schedule.py`.

### Slice L — Dissociation ladder + integrated gate *(pure; fully independent)*

**Module.** `physics/dissociation_ladder.py`

**Purpose.** The Form U sigmoid ladder `D_0(n)` and its cumulative sum `Σ(n)` (the
self-bound gate threshold). No physics state. The Tier-2 analog of Slice S — built and
tested first.

**Interface (proposed).**
- `d0_of_n(n, *, picture, kappa) -> eV` — single-rung dissociation energy.
- `sigma(n, *, kappa)` — the centered sigmoid (helper).
- `ladder_cumsum(n, *, picture, kappa) -> eV` — `Σ(n)`, the gate threshold (consumed by K
  and U).
- `gate_threshold(n, ...)` — readable alias for `ladder_cumsum`.
- Optional `tabulated_ladder(...)` **declared fallback** (a non-monotone / second-cliff
  histogram would trigger it; not the default).
- `picture` defaults to `statistical_mixture`; `D_floor`, `n*` pulled from `constants.py`.

**Encoded form.** §2.2 (L).

**Knobs (config §5):** `ladder_electronic_picture ∈ {statistical_mixture(default),
x2_only, cooling_relaxed}` (sets `D_0(1)` = 74.4 / 106.9 / between); `ladder_steepness` κ
(Free knob); `D_floor`, `n*` sourced; `dissociation_ladder` selects Form U vs tabulated.

**Oracle values.**
- First rung: `x2_only` → 106.9 cm⁻¹ = 0.013254 eV; `statistical_mixture` → 74.4 cm⁻¹ =
  0.009224 eV; `cooling_relaxed` strictly between.
- Floor: 4.97 cm⁻¹ = 6.16e-4 eV; cliff at `n*+½ = 21.5`.
- Integrated first-shell sum `Σ(21)`: **X₂ 0.25–0.28 eV / mixture 0.17–0.19 eV**, varying
  only ~11% over κ∈[0.3,5] — the **near-κ-independence / decoupling** check (the gate
  threshold and the S2 floor are picture-set and pinnable *ahead* of the κ fit).
- Drag-binding cross-check: a moderate-κ ladder meets `E_BIND_DRAG_eV = 0.1168 eV`
  naturally (representative sums 835–1571 cm⁻¹); `|S|` is an **upper bound**, never a
  rung-sum target (rungs must **not** be calibrated to `|S|/n*`).

**Independence.** Depends only on constants + κ + picture. Built and tested first; K and U
mock it.

**Test spec (`tests/test_dissociation_ladder.py`).**
- Rungs/sums match oracle to 4 figures across κ∈[0.3,5] and all 3 pictures.
- `Σ(21)` lands in the X₂ / mixture bands and stays within ~11% across the κ range
  (decoupling assertion).
- Cliff geometry: `D_0(n*) ≈ D_0(1)` (full-depth in-shell), `D_0(n*+1) → D_floor` for
  large κ.
- Monotone non-increasing in `n`; `D_0(n) ∈ [D_floor, D_0(1)]`.
- Vectorized scalar-in/array-in shape parity.
- Fail-loud: `D_0(1) > D_floor` enforced; unknown `picture` rejected.

**Acceptance.** All pass; gate threshold ≈κ-independent; tabulated fallback round-trips a
hand-built ladder.

---

### Slice K — Newton cooling + occupancy-resolved asymptote *(pure; mocks L)*

**Module.** `physics/internal_energy_cooling.py`

**Purpose.** The cooling **driver** for `E_solv.struct` and the binding split it cools
toward. A pure analytic step — **no integrator state**; the per-step closed-form
relaxation only.

**Interface (proposed).**
- `s_collective_eV(N, *, picture, kappa)` — `|S(N)| = |S|·Σ(N)/Σ(n*)`.
- `e_infinity_eV(N, ...)` — `E_∞(N) = −|S(N)|`.
- `e_bind_pair_eV(N, ...)` — `−Σ(N)` (thin wrapper over L's `ladder_cumsum`).
- `e_electrostriction_eV(N, ...)` — `E_elec(N) = −(|S(N)| − Σ(N)) ≤ 0`.
- `newton_cool_step(E_solv_struct_eV, N, *, tau_ps, dt_ps, ...)` — one closed-form
  exponential relaxation toward `E_∞(N)` (exact `e^{−dt/τ}` form so it is `dt`-robust and
  needs no integrator). Consumes L's `ladder_cumsum`.

**Encoded form.** §2.2 (K).

**Knobs (config §5):** `internal_energy_cooling_tau_ps` τ (Bounded sweep [2.6,16.5]);
`solv_struct_asymptote` |S| = 0.308 eV; `electrostriction_binding` (derived, surfaced for
inspection).

**Oracle values.**
- `E_∞(N) → 0` as `N → 0` (OQ6: total strip stays dynamically reachable; a *fixed*
  full-shell `E_∞` would mechanically halt shedding — this slice must reproduce the
  occupancy-resolved form, not the fixed one).
- `|S|/n* = 118 cm⁻¹` (170 K) at `n* = 21` (per-atom collective marginal).
- `E_elec(N) ≤ 0` for all N; electrostriction is the *dominant* binding term
  (≈4.7× the pair-at-radius ~25 cm⁻¹/atom).
- Cold-shed neutrality: `ΔE_solv.struct = 0` across a shed (the K1∩K2 non-overlap
  identity, modulo the A8 marginal-electrostriction bath booking) — so K1 (discrete) and
  K2 (continuous) never double-count and `τ_GAH25` transplants directly.
- `E_int^eq = 0` exactly under the split (R12 systematic eliminated).
- Rate units eV/ps.

**Independence.** Mocks L's `ladder_cumsum` with a known stub; no integrator, no state.

**Test spec (`tests/test_internal_energy_cooling.py`).**
- `e_infinity_eV` monotone in N and `→ 0` at `N = 0`.
- `e_electrostriction_eV ≤ 0` everywhere; equals `−(|S(N)| − Σ(N))`.
- `s_collective_eV(n*) = |S| = 0.308 eV`; `|S(n*)|/n* = 118 cm⁻¹` to 4 figures.
- `newton_cool_step` relaxes toward `E_∞` with the exact factor `e^{−dt/τ}`; fixed point
  at `E_solv.struct = E_∞`; `dt`-robust (same endpoint via one big step vs many small).
- Cold-shed neutrality identity holds to machine precision on a synthetic
  `(N → N−1, E_int −= D_0)` micro-case.

**Acceptance.** Newton step + split match oracle; neutrality identity holds; `E_∞`
monotone and reaches 0.

---

### Slice U — `E_int` budget bookkeeping rules *(pure; mocks L, K)*

**Module.** `physics/internal_energy_budget.py`

**Purpose.** The per-channel ΔE_int helpers + the post-t× reconstruction — the bricks the
Phase-C driver and the 5-term ledger consume. Pure; no reservoir state lives here (the
*state* is Phase C).

**Interface (proposed).**
- `e_int_onset_eV(*, f_int, e_avail_eV)` — S2: `f_int·E_avail`.
- `dE_int_pickup_eV(n, *, f_ret, picture, kappa)` — S1: `+f_ret·D_0(n+1)`; companion
  `pickup_bath_release_eV(...)` returns `(1−f_ret)·D_0(n+1)`.
- `dE_int_shed_eV(n, *, picture, kappa)` — K1: `−D_0(n)`.
- `reconstruct_e_int_eV(E_solv_struct_eV, N, *, picture, kappa, post_crossing: bool)` —
  the A9 reconstruction; **raises** loudly if `post_crossing` is False (pre-t× use is
  invalid). Consumes L (rungs) and K (binding split).
- `f_int_floor(*, e_avail_eV, picture, kappa)` — the derived self-unbound floor
  `Σ(n*)/E_avail` (diagnostic helper for the config bound).

**Encoded form.** §2.2 (U).

**Knobs (config §5):** `internal_energy_partition_fraction` f_int (Bounded; scenario-keyed
floor), `internal_energy_retained_fraction` f_ret (Bounded), `coulomb_available_eV`
E_avail (exists; 0.80 validation / 2.70 production).

**Oracle values.**
- f_int floor `= Σ(n*)/E_avail`:
  - @ **0.80 eV**: X₂ **0.31–0.35**, mixture **0.21–0.24**;
  - @ **2.70 eV**: X₂ **0.09–0.10**, mixture **~0.065**.
  (Verified numerically; the self-unbound onset is robust across both budgets but the
  headroom shrinks from ~14× at 2.70 eV to ~4.5× at 0.80 eV.)
- S1/K1 sign + magnitude: pickup heats by `+f_ret·D_0(n+1)`, shed drains by `−D_0(n)`,
  bath gets the S1 remainder.
- Jump-consistency: the `E_solv.struct ↔ E_int` variable swap leaves S1/K1 unchanged
  (MASS §6 "jump-consistency check") — at fixed N, `dE_solv.struct = dE_int`.

**Independence.** Mocks L (rungs) and K (split) with stubs; pure arithmetic on known
inputs.

**Test spec (`tests/test_internal_energy_budget.py`).**
- Each rule reproduces the MASS §6 identity exactly (sign + magnitude).
- S1 split closes: `f_ret·D_0 + (1−f_ret)·D_0 = D_0` to machine precision.
- f_int floor matches the scenario-keyed oracle at both budgets and both pictures.
- `reconstruct_e_int_eV` **raises** when `post_crossing=False`; returns the split-exact
  value when True (`E_int^eq = 0` recovered on an equilibrium micro-case).
- Jump-consistency: a synthetic swap leaves the S1/K1 increments invariant.

**Acceptance.** Each rule reproduces the §6 identity; reconstruction guard fires pre-t×;
floor matches the oracle.

---

## 4. Golden-oracle fixtures (single source of locked values)

The consolidated table the three suites assert against (all derived from §2 + verified
conversions):

| Quantity | Value | Used by |
|---|---|---|
| `D_0(1)` X₂ / mix | 106.9 cm⁻¹ = 0.013254 eV / 74.4 cm⁻¹ = 0.009224 eV | L |
| `D_floor` | 4.97 cm⁻¹ = 6.16e-4 eV (7.15 K) | L, K |
| cliff center | `n*+½ = 21.5` | L |
| `Σ(21)` band (κ∈[0.3,5]) | X₂ 0.25–0.28 eV / mix 0.17–0.19 eV (~11% spread) | L, U |
| `\|S\|` and `\|S\|/n*` | 0.308 eV = 2484 cm⁻¹; 118 cm⁻¹ @ n*=21 | K |
| `E_bind` drag cross-check | 0.1168 eV = 942 cm⁻¹ (upper bound, not a target) | L |
| f_int floor @ 0.80 eV | X₂ 0.31–0.35 / mix 0.21–0.24 | U |
| f_int floor @ 2.70 eV | X₂ 0.09–0.10 / mix ~0.065 | U |
| cooling neutrality | `ΔE_solv.struct = 0` across a shed | K |
| `E_int^eq` | 0 (exact, under the split) | K, U |

---

## 5. Config contract (new `SimConfig` fields landing in Phase A)

Per the rule-2 declared-but-unread exception, each field lands with its owning slice and
is removed from the exception table when its slice activates it.

- **L:** `dissociation_ladder` (Form U selector / tabulated fallback), `ladder_steepness`
  κ, `ladder_electronic_picture ∈ {statistical_mixture, x2_only, cooling_relaxed}`,
  `solv_struct_asymptote` |S|, `electrostriction_binding` (derived).
- **K:** `internal_energy_cooling_tau_ps` τ.
- **U:** `internal_energy_partition_fraction` f_int,
  `internal_energy_retained_fraction` f_ret. (`coulomb_available_eV` already exists,
  `config.py:210` — scenario-keyed 0.80 / 2.70.)

**Validators.** Reuse the existing config-load guard scaffolding (`check_drag_config`,
`config.py:~314+`) for any Phase-A checks: the `ladder_electronic_picture` enum reject arm
(mirroring the `mass_scenario` reject arm), and the `D_0(1) > D_floor` assertion. These
are load-time fail-loud guards, not silent clamps.

**Out-of-Phase-A flag.** The `MassScenario` literal / `biphasic` vs `biphasic_energy_gated`
naming reconcile (Tier-2 plan §5) is a **Phase C** concern — note it, do not touch it here.

---

## 6. Dependency / build order (within Phase A)

```
L (ladder D_0(n), gate Σ(n))   ── independent; built & tested first
        │
        ├──►  K (cooling; consumes Σ(n) via ladder_cumsum)
        │
        └──►  U (budget; consumes rungs from L and the split from K)
```

L is fully independent. K and U each **mock** L's `ladder_cumsum` (and, for U, K's split)
with known stubs for their own unit tests, then compose against the real L/K. No Phase-A
slice depends on Phase B (channels) or Phase C (driver/schema). Phase A is parallelizable
to the extent its mocks allow, but the natural order is L → K → U.

---

## 7. Testing methodology

Mirror Tier-1a's pure-function discipline (validation steps 1–2):
- closed-form **oracle** asserts with **tight analytical tolerances** (4 figures) for
  rungs, floors, conversions;
- **band-membership** asserts for the κ-range integrated sums and f_int floors;
- **identity** asserts (`= 0` to machine precision) for cooling neutrality, the S1 split,
  and `E_int^eq`.
- **No figures, no checkpoints, no RNG** in pytest.

Suites: `tests/test_dissociation_ladder.py`, `tests/test_internal_energy_cooling.py`,
`tests/test_internal_energy_budget.py`. Run the narrowest first, then the full suite.
Interpreter (Python may not be on PATH):

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

---

## 8. Acceptance criteria

**Per slice:** the slice's own suite green with neighbours mocked.

**Phase A overall:** full suite stays green (current baseline is all-green on
`drag_implementation`); L/K/U reproduce every §2/§4 oracle to the stated tolerance; the
decoupling check (gate threshold ≈κ-independent) holds; all fail-loud guards fire on bad
input.

**Out-of-scope guard (fails review if it leaks in):** any Phase-A code path that draws
RNG, references the integrator / O-step / jump seam, persists `E_int` *state*, computes a
5-term closure, reads a noise amplitude, or touches the Tier-0 drag law. Those are Phase
B / C / Tier 3.

---

## 9. Risks / notes carried

- **OQ1 (electronic-picture provenance) — carry-both by design (not a Phase-A action
  item).** L exposing all three pictures (`statistical_mixture` default, `x2_only`,
  `cooling_relaxed`) as the Tier-2 co-fit knob *is* the settled design response to the
  open provenance question — deliberate, not a loose end. The default stays
  `statistical_mixture` (do not hard-collapse), and `E_bind` is not cited as independent
  support for the mixture. Provenance remains pending author contact; revisit only if a
  resolution arrives.
- **Mechanism locked, values open** — Phase A encodes *forms* and *sourced* anchors only;
  κ, f_int, f_ret, τ stay knobs. **No fitting in Phase A** (that is the Phase F campaign).
- **`|S|` is an upper bound, not a rung-sum target** (R3/R12) — rungs are never calibrated
  to `|S|/n*`; the reachable integrated cross-check is the drag binding `E_BIND_DRAG_eV`.
- **CLAUDE.md status (reconciled 2026-06-29):** `CLAUDE.md` already records **Tier 1a —
  DELIVERED** and **Tier 2 — ACTIVE (current goal)** (Project Snapshot + validation
  hierarchy). The earlier "snapshot frames Tier 1a as ACTIVE" drift note is **stale and
  dropped** — no reconcile outstanding.
- **Reconstruction is post-t× only (A9)** — the U guard must enforce this loudly so no
  downstream slice reconstructs `E_int` before the gate crossing.
