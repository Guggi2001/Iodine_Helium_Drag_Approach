# Tier 2 — Phase B Implementation Plan (Stochastic Mass Channels)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no pseudo-code until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations from
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (LaTeX + dimensional checks);
> module descriptions are *interface contracts*, not implementations.
>
> **Entry docs:** `TIER2_IMPLEMENTATION_PLAN.md` §3/§4 (Phase B in the dependency
> graph); `TIER2_PHASE_A_IMPLEMENTATION_PLAN.md` (the energetics layer Phase B
> consumes — L rungs, K split, U ΔE_int rules); `MASS_DYNAMICS_LOCKED_*.md` §1/§4
> (two-channel mechanism + dimensional table), §5 (density-only confirmed), §6 S1/K1
> (the budget channels Phase B emits into), §6.11 (Π order parameter, t×, regime),
> §11 (config), A11 (ν/s sourcing, n=1 direct dissociation, the s≥1 guard), A12
> (Langmuir occupancy cap), A13 (jump operators, one-event-per-step, jump-then-O).
> **Reset spine:** the delivered Tier-1a `physics/mass_jump.py`
> (`cold_shed`, `continuous_velocity_shed`, `kick_factor`,
> `_reduced_mass_defect_coeff`, and their `*_components` variants).

---

## 0. Status and intent

Tier 0 (drag form locked: `shared_pure_cubic`, `γ=g·b·v²`) and Tier 1a (anchored
kinematic mass-dynamics validation) are **delivered**; **Phase A is delivered as a
plan** (`TIER2_PHASE_A_IMPLEMENTATION_PLAN.md`) — the three pure, stateless energetics
primitives L/K/U. **Phase B builds the two stateful, RNG-driven channels** that turn
the locked mechanism into a *generative* integer-`n` process: the Poisson **pickup**
gain channel and the energy-gated, RRK-rate-limited **evaporation** loss channel, plus
the helium **density profile** that gates pickup.

Phase B asks **one** question per slice and answers it with a stochastic primitive
checked against a closed-form moment/identity oracle under a *seeded/mock RNG*:
*what is the local He density gate `ρ_He(depth)/ρ_bulk` (ρ); what is the per-step
He-capture probability, its momentum-conserving velocity reset, and the S1 heat it
deposits (P); and what is the per-step gated RRK shed probability, its cold-shed reset,
and the K1 drain it produces (Q)?* These map onto **validation step 6** (statistical)
of the scientific-code-caution order, but at the **unit** level — single-channel
distribution/moment checks with neighbours mocked, *not* a full stochastic trajectory.

**Central engineering fact.** Phase B is the first place RNG, per-event mass change, and
velocity resets appear. But every channel stays a **pure-ish primitive**: it takes the
current `(n, E_int, v, m)` plus an injected RNG and returns the post-event
`(n', E_int', v', m', ΔE_mass_transfer)` — it performs **no integration, no `E_int`
state persistence, no schema I/O, and no 5-term closure** (those are Phase C: X schema +
G driver). The channels consume Phase A (L rungs `D_0(n)`, K cooling/split, U ΔE_int
rules) and reuse the Tier-1a reset spine; only the **+He capture counterpart** of
`mass_jump` and the two channel modules are new physics here.

**Locked decisions (2026-06-29, this session).**
- **Module layout — three one-concern modules + capture in `mass_jump.py`** (mirrors
  the Phase-A 3-module split and the ρ/P/Q independence in the Tier-2 dependency graph):
  - `physics/helium_density.py` — Slice ρ,
  - `physics/pickup.py` — Slice P,
  - `physics/evaporation.py` — Slice Q,
  - `physics/mass_jump.py` — **extended** with the symmetric `capture` reset (the +He
    counterpart of `cold_shed`/`continuous_velocity_shed`), so the reduced-mass defect
    formula stays in one place (CLAUDE.md rule 1).
- **ρ sourcing — reuse the erf-complement gate.** Slice ρ uses the same
  `0.5·(1−erf(depth/steepness))` form as `drag.py::spatial_gate` (steepness 14.2 Å) for
  the `ρ_He/ρ_bulk` ratio; a **sourced baseline/TDDFT density profile is the declared
  fallback** (rule-2 convention), *not* built in Phase B. The drag spatial gate stays
  **G2** (no G2→G4 promotion this phase; CALIBRATION_MAP row 5).
- New knobs land in `config.py`; the one new sourced constant (`λ_0` prior) is a config
  default, not a `constants.py` entry (it is Bounded/calibrated, not pinned).

---

## 1. Scope lock — what Phase B does and does not do

**Does:**
- the helium **density gate** `ρ_He(depth)/ρ_bulk ∈ [0,1]` feeding pickup (ρ);
- the **pickup channel**: the per-step Bernoulli draw `P_attach = 1−e^{−λ·dt}` with
  `λ_attach = λ_0·(ρ_He/ρ_bulk)·(1−n/n*)_+^p`, the **momentum-conserving capture reset**
  `v⁺ = m/(m+m_He)·v⁻` (He at rest), and the **S1** heat `+f_ret·D_0(n+1)` deposited via
  U (P);
- the **+He capture reset** in `mass_jump.py` (symmetric to `cold_shed`), exposing the
  reduced-mass capture defect that becomes `E_mass_transfer` at Phase C;
- the **evaporation channel**: the **self-bound gate** (suppress while `E_int > Σ(n)`),
  the saturating **RRK rate** `k = ν·(1−D_0(n)/E_int)^{s−1}` with `P_shed = 1−e^{−k·dt}`
  (`s = 3n−3`, `n≥2`; **`n=1` direct dissociation `k=ν`**), reuse of
  `mass_jump.cold_shed` for the reset, and the **K1** drain `−D_0(n)` via U (Q).

**Does NOT (fails review if it leaks in):**
- any **integrator / O-step** coupling, the jump-then-O seam, or `m(t)` plumbing — the
  generative driver (G; Phase C);
- any **`E_int` state persistence**, checkpoint schema bump (v6→v7), or the **5-term
  invariant closure** — that is the schema slice (X; Phase C);
- the **noise model** (Tier 3 — stays inert behind its enum);
- reading the **Tier-0 drag law** or any `γ`/force-coefficient quantity;
- *fitting* λ_0/f_ret/ν/s/κ — Phase B encodes forms and exposes knobs; calibration is
  the Phase F campaign.

**Mass-agnostic guard.** The channels take the ion mass `m` only as the **shed/capture
quantity** in the reset (exactly as `mass_jump` already does), never as a drag-law input.
The unified friction convention (CLAUDE.md "Friction convention") is preserved: no
Phase-B function evaluates `γ`, and the post-jump O-step at `m⁺` is Phase C, not here.

**RNG contract.** Every stochastic primitive takes an **injected** `numpy.random.Generator`
(no module-level global RNG, no implicit seeding). The **draw order** (pickup vs shed) is
*specified* here for unambiguous bookkeeping but is **locked** as a forbidden-list item at
Slice X (Phase C) — Phase B documents it, Phase C freezes it in the schema/driver.

---

## 2. Locked physics — constants and equations (the test oracles)

All numbers are sourced from the MASS doc and `CALIBRATION_MAP.md`; energies convert via
`EV_PER_WAVENUMBER = 1/8065.543937 = 1.23984e-4 eV/cm⁻¹` (already in `constants.py`).
The Phase-A rung/floor/`n*`/ν constants (`D0_1_X2`, `D0_1_MIX`, `D_FLOOR`, `N_STAR`,
`NU_EVAP_PER_PS`) land with Phase A; Phase B **consumes** them.

### 2.1 Constants and priors

| Symbol | Value | Units | Class / source | Read by |
|---|---|---|---|---|
| `MASS_HE_AMU` | 4.0026 | amu | exists, `constants.py:109` | P, Q (capture/shed) |
| `N_STAR` | 21 | count | Phase-A const; [I2-notes] cation first shell | P (cap), Q (gate via Σ) |
| `NU_EVAP_PER_PS` (ν) | 2.42 | ps⁻¹ | Phase-A const; **pinned** ([IHe05] curvature, A11) | Q |
| `s` (RRK DOF) | `3n−3` (n≥2); `n=1` direct | — | **Derived**, mode-counted (A11); guarded `s≥1` | Q |
| `λ_0` (pickup coeff.) | central ~0.7–1.1, ±factor-2 | ps⁻¹ | **Sourced+Bounded** (GAH25 Rb⁺/Cs⁺; *not* 2.0 Na⁺; CALIBRATION_MAP row 7) | P |
| `p` (occupancy exponent) | default tied to κ | — | **Derived** (A12 `p↔κ`); split only if size dist. demands | P |
| `f_ret` | small prior, `[0,1]` | — | **Bounded** (CALIBRATION_MAP row 13) | P (S1 heat) |
| `ρ_He/ρ_bulk` | erf-complement, steepness 14.2 Å | — | **Sourced** (row 8); Phase-B reuse of `spatial_gate` | ρ → P |
| `dt_ion` | 0.01 | ps | integrator step (MASS §4) | P, Q (per-step prob.) |

### 2.2 Encoded forms (LaTeX + dimensional check)

**(ρ) Helium density gate.** Reuse of the drag erf-complement, read as a density ratio:
$$
\frac{\rho_\text{He}(\text{depth})}{\rho_\text{bulk}}
= \tfrac12\big(1-\operatorname{erf}(\text{depth}/\text{steepness})\big)\in[0,1],
\qquad \text{depth}=r_\text{atom}-r_\text{droplet}.
$$
*Dim:* dimensionless ratio; `depth`, `steepness` in Å. ✓ →1 deep inside (`depth≪0`),
0.5 at the nominal surface, →0 outside (`depth≫0`, ion exit → λ→0 → termination).

**(P) Pickup — Poisson→Bernoulli per step, capture reset, S1 heat.**
$$
P_\text{attach}(dt)=1-e^{-\lambda_\text{attach}\,dt},\qquad
\lambda_\text{attach}=\lambda_0\,\frac{\rho_\text{He}}{\rho_\text{bulk}}\,
\Big(1-\tfrac{n}{n^*}\Big)_+^{\,p},
$$
$$
\text{on fire: } n\to n+1,\ \ m\to m+m_\text{He},\ \
v^+=\frac{m\,v^-+m_\text{He}u_\text{He}}{m+m_\text{He}}\xrightarrow{u_\text{He}=0}\frac{m}{m+m_\text{He}}v^-,\ \
\Delta E_\text{int}=+f_\text{ret}D_0(n{+}1).
$$
Capture KE defect (the +He counterpart of the shed defect, books to
`E_mass_transfer` at Phase C):
$$
\Delta E_\text{mass\_transfer}^\text{capture}
=\tfrac12\frac{m\,m_\text{He}}{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2\ \ (>0,\ u_\text{He}=0).
$$
*Dim:* `λ` ps⁻¹ (`[λ dt]=1` ✓); occupancy factor dimensionless (`n,n*` counts, `p` pure);
defect `amu·Å²/ps²` = energy ✓. The S1 split closes: `f_ret·D_0 + (1−f_ret)·D_0 = D_0`,
remainder to bath (U).

**(Q) Evaporation — self-bound gate + saturating RRK + cold-shed + K1 drain.**
$$
\text{gate: suppress all sheds while } E_\text{int}>\Sigma(n)=\sum_{i\le n}D_0(i);
$$
$$
\text{else } k(E_\text{int},n)=\nu\Big(1-\frac{D_0(n)}{E_\text{int}}\Big)^{s-1}\ (n\ge2),\quad
k=\nu\ (n{=}1,\ \text{direct, gated }E_\text{int}>D_0(1)),\quad
P_\text{shed}=1-e^{-k\,dt};
$$
$$
\text{on fire: } n\to n-1,\ \ E_\text{int}\mathrel{-}=D_0(n),\ \
v^+=\frac{m}{m-m_\text{He}}v^-\ (\text{cold shed, }u_\text{He}\approx0),\ \
\Delta E_\text{mass\_transfer}=-\tfrac12\frac{m\,m_\text{He}}{m-m_\text{He}}\|v^-\|^2.
$$
*Dim:* `s=3n−3` dimensionless (full complex, A11; `s−1≥1` for `n≥2`); bracket
`1−D_0/E_int` dimensionless (eV/eV); `k` ps⁻¹, **bounded `k∈[0,ν)`** for `s>1` and `k=ν`
at `n=1`; per-step shed prob. `≤ 1−e^{−ν dt} ≈ ν·dt` (**no gate-open avalanche**).
Below threshold (`E_int ≤ D_0(n)`) the rate is exactly zero.

**(Π) Regime order parameter (derived diagnostic, not a channel output).**
$$
\Pi(n)\equiv\lambda(n)\,f_\text{ret}\,\tau\quad(\text{dimensionless: }\text{ps}^{-1}\!\cdot1\cdot\text{ps}\ \checkmark);
\qquad \Pi>1\ \text{shedding persists},\ \Pi<1\ \text{freeze}.
$$
Surfaced as a small pure helper (post-hoc reconstructable; the actual regime reconstruction
is Slice D2, Phase E). Included here only so the pickup/cooling priors can be sanity-checked
against the §6.11 stability picture.

No term fails dimensional balance.

---

## 3. New-work slices

Each slice carries the Tier-1a/Phase-A template: **Purpose / Module + interface / Encoded
form / Knobs / Oracle / Independence / Test spec / Acceptance**. Interfaces are *proposed
contracts* (final names settle at build), vectorized scalar-in/array-in following
`mass_jump.py` and `shell_schedule.py`.

### Slice ρ — Helium density gate *(pure; fully independent)*

**Module.** `physics/helium_density.py`

**Purpose.** The `ρ_He(depth)/ρ_bulk ∈ [0,1]` field that gates pickup. Kept a distinct
sourced quantity feeding P, even though it reuses the drag erf form, because its *role*
(occupancy/capture gate) is independent of drag and it carries the declared sourced-profile
fallback.

**Interface (proposed).**
- `rho_he_ratio(depth, *, steepness=POTENTIAL_STEEPNESS_ANGSTROM) -> ndarray` — the
  erf-complement ratio in `[0,1]`. Thin, deliberately delegates to (or mirrors) the same
  erf machinery as `drag.spatial_gate` so the two never drift (CLAUDE.md rule 1: a shared
  `_erf_complement(depth, steepness)` helper is the single source).
- Optional `tabulated_density_profile(...)` — **declared fallback** (a sourced baseline/
  TDDFT `ρ_He(r)`); not the default, returns the same `[0,1]` ratio contract.

**Encoded form.** §2.2 (ρ).

**Knobs (config §6):** `helium_density_profile` (existing placeholder, `config.py:220`) —
selector `erf_complement(default)` / `tabulated`; reuses `potential_steepness` (14.2 Å).

**Oracle values.** ratio = 1 at `depth≪0`, **0.5 at `depth=0`**, →0 at `depth≫0`;
monotone non-increasing in `depth`; identical to `drag.spatial_gate(depth, 14.2)` to
machine precision (shared-form check).

**Independence.** Depends only on `depth` + steepness. Built and tested first; P mocks it.

**Test spec (`tests/test_helium_density.py`).**
- ratio ∈ [0,1], =0.5 at depth 0, monotone non-increasing; →0 / →1 at the asymptotes.
- **Shared-form parity:** equals `drag.spatial_gate` on a `depth` grid to machine precision
  (guards against the two erf forms diverging).
- Vectorized scalar-in/array-in shape parity; fail-loud on `steepness ≤ 0`.
- Tabulated fallback round-trips a hand-built `[0,1]` profile.

**Acceptance.** All pass; ratio matches `spatial_gate`; fallback round-trips.

---

### Slice P — Pickup channel + `mass_jump` capture reset *(stateful primitive; mock RNG)*

**Modules.** `physics/pickup.py` (the channel) and `physics/mass_jump.py` (**extended**
with the capture reset).

**Purpose.** The gain channel: one independent Bernoulli draw per ion per step; on fire,
apply the momentum-conserving **capture** reset and report the S1 heat. The reset is the
symmetric +He counterpart of `cold_shed` and lives in `mass_jump.py` so the reduced-mass
defect coefficient stays in one place.

**Interface (proposed).**
- *In `mass_jump.py`:* `capture(v_minus, m_minus_amu, *, m_he_amu=MASS_HE_AMU, u_he=0.0)`
  → `ShedResult`-shaped tuple with `v_plus = m/(m+m_He)·v_minus`,
  `m_plus_amu = m + m_He`, `dE_mass_transfer = +½·(m·m_He)/(m+m_He)·|v⁻|²` (reuses
  `_reduced_mass_defect_coeff` with the **+** sign convention); plus the `*_components`
  vectorized form mirroring `cold_shed_velocity_components`.
- *In `pickup.py`:*
  - `attach_probability(lambda_attach, dt_ps) -> P_attach` (`1−e^{−λ dt}`).
  - `lambda_attach(rho_ratio, n, *, lambda0, n_star, p, cap)` — the rate with the Langmuir
    factor (`cap="langmuir"` applies `(1−n/n*)_+^p`; `"none"` recovers density-only).
  - `pickup_step(state, *, rng, rho_ratio, lambda0, f_ret, ...)` — draw one Bernoulli; on
    fire compose `mass_jump.capture` (reset/mass/defect) + `internal_energy_budget.
    dE_int_pickup_eV` (S1 heat); return post-event `(n', m', v', ΔE_int, ΔE_mass_transfer,
    fired)`. Consumes ρ (gate), U (S1).

**Encoded form.** §2.2 (P).

**Knobs (config §6):** `pickup_rate_coefficient` λ_0 (**repurpose/rename**
`mass_rate_coefficient`, `config.py:218`); `pickup_rate_form ∈ {density_only(default),
sweeping, dwell_time}` (existing `mass_rate_form`, `config.py:217`);
`pickup_occupancy_cap ∈ {langmuir(default), none}`; `pickup_occupancy_exponent` p
(default tied to κ); `he_capture_velocity ∈ {at_rest(default), thermal}` (thermal
deferred — Tier 3 noise coupling).

**Oracle values.**
- **Momentum/defect exact:** capture reset reproduces `v⁺ = m/(m+m_He)·v⁻` and the
  reduced-mass defect to machine precision; round-trip `capture` then `cold_shed` returns
  the original mass and a net-zero defect identity at fixed `|v|` only up to the order of
  operations (documented, not asserted as exact across the asymmetric resets).
- **Poisson moments:** over many seeded steps at constant `λ`, the count of fires per
  window matches Poisson `P_n(t)` mean `λt` and variance `λt` within a sample-size band
  (A1).
- **Cap behaviour:** `(1−n/n*)_+^p → 0` at `n=n*`; `cap="none" ≡ "langmuir"` for the
  ejection case where `n ≪ n*` (cap inert).
- **S1 sign/magnitude:** heat `= +f_ret·D_0(n+1) ≥ 0`; bath remainder `(1−f_ret)·D_0`.

**Independence.** Mocks ρ (known ratio stub), U's `dE_int_pickup_eV` (known stub), and the
RNG (seeded `Generator` / mock). No integrator, no `E_int` state.

**Test spec (`tests/test_pickup.py`, + capture asserts in `tests/test_mass_jump.py`).**
- `capture` reset: `v⁺`, `m⁺`, defect exact; sign opposite to `cold_shed`; fail-loud on
  non-finite / non-positive masses (reuse `_check_masses`).
- `attach_probability` matches `1−e^{−λ dt}`; `→ λ dt` for small `λ dt`.
- Poisson moment check (seeded, sample-size tolerance); independence across ions.
- Langmuir cap →0 at `n*`; density-only/langmuir equivalence for `n ≪ n*`.
- S1 heat sign + magnitude via the mocked U stub.

**Acceptance.** Momentum/defect exact; Poisson moments within the sample-size band; cap →0
at `n*`; S1 sign correct.

---

### Slice Q — Evaporation channel *(stateful primitive; mock RNG; reuse `cold_shed`)*

**Module.** `physics/evaporation.py`

**Purpose.** The loss channel: the parameter-free self-bound gate, then the saturating RRK
shed of the top rung; on fire reuse `mass_jump.cold_shed` and report the K1 drain.

**Interface (proposed).**
- `rrk_rate(E_int_eV, n, *, nu, picture, kappa) -> k_per_ps` — the RRK rate with the
  `s=3n−3` bracket for `n≥2`, the **`n=1` direct `k=ν`** branch, and `k=0` below threshold;
  consumes L's `d0_of_n`.
- `shed_probability(k_per_ps, dt_ps) -> P_shed` (`1−e^{−k dt}`).
- `is_self_bound(E_int_eV, n, *, picture, kappa) -> bool` — the gate `E_int > Σ(n)`
  (suppress while True); consumes L's `ladder_cumsum`.
- `evaporation_step(state, *, rng, nu, ...)` — if self-unbound, return no-shed; else draw
  one Bernoulli, and on fire compose `mass_jump.cold_shed` (reset/mass/defect) +
  `internal_energy_budget.dE_int_shed_eV` (K1 drain `−D_0(n)`); return post-event tuple.
- `effective_dof(n)` — `3n−3` (n≥2), with the **config-load `s≥1` guard** applied to any
  effective-scalar override (rejects `s<1`, the divergent-rate regime).

**Encoded form.** §2.2 (Q).

**Knobs (config §6):** `evap_rate_prefactor_per_ps` ν (Sourced, pinned 2.42);
`evap_rrk_dof` s (Derived `3n−3`; optional effective-scalar override, **guarded `s≥1`** at
config-load). `evap_gate_onset_eV` is **retired → optional override** (the gate is the
integrated ladder `Σ(n)`, parameter-free — no free onset knob).

**Oracle values.**
- `k ∈ [0, ν)` for `s>1` (`n≥2`); `k=ν` exactly at `n=1`; `k=0` for `E_int ≤ D_0(n)`.
- **Gate suppresses all sheds** while `E_int > Σ(n)` (net self-unbound) — zero fires
  regardless of RNG.
- **No avalanche:** per-step shed probability `≤ 1−e^{−ν dt}` even as `E_int→∞`
  (saturation), so at most ~`ν·dt` per step; **one rung per step**.
- **`s≥1` guard fires** on an override `s<1` at config-load (clear error, no silent clamp).
- K1 drain `= −D_0(n)` exactly (via the mocked U stub).

**Independence.** Mocks L (`d0_of_n`, `ladder_cumsum`), U (`dE_int_shed_eV`), and the RNG
with known stubs. Reuses the *real* `mass_jump.cold_shed` (already delivered + tested) for
the reset. No integrator, no `E_int` state.

**Test spec (`tests/test_evaporation.py`).**
- `rrk_rate` matches the oracle: bounded `[0,ν)`, `=ν` at `n=1`, `=0` below threshold,
  monotone increasing in `E_int`.
- Gate: `is_self_bound` flips at `E_int = Σ(n)`; suppresses all sheds while self-unbound.
- Avalanche guard: per-step prob `≤ 1−e^{−ν dt}` across a wide `E_int` sweep; exactly one
  rung drawn per fire.
- `effective_dof` = `3n−3`; the config-load `s≥1` guard rejects `s<1` (load-time fail-loud).
- K1 drain sign + magnitude via the mocked U stub; shed reset delegates to `cold_shed`
  (no second copy of the reduced-mass formula).

**Acceptance.** Rate matches oracle; gate + avalanche guarantees hold; `s≥1` guard fires;
K1 drain exact; reset reuses `cold_shed`.

---

## 4. Golden-oracle fixtures (single source of locked values)

| Quantity | Value | Used by |
|---|---|---|
| density gate at depth 0 | 0.5; →1 inside, →0 outside | ρ |
| density ↔ `spatial_gate` parity | machine precision @ steepness 14.2 Å | ρ |
| `P_attach`, `P_shed` | `1−e^{−x dt}`; `→ x dt` small-`x dt` | P, Q |
| capture reset | `v⁺=m/(m+m_He)·v⁻`; defect `+½(m·m_He)/(m+m_He)·|v⁻|²` | P |
| cold-shed reset (reuse) | `v⁺=m/(m−m_He)·v⁻`; defect `−½(m·m_He)/(m−m_He)·|v⁻|²` | Q |
| Poisson moments | mean=var=`λt` (sample band) | P |
| Langmuir cap | `(1−n/n*)_+^p`; →0 at `n=n*=21` | P |
| RRK rate | `k∈[0,ν)` (`s>1`); `k=ν` (`n=1`); `k=0` (`E_int≤D_0`) | Q |
| `ν` | 2.42 ps⁻¹ (pinned) | Q |
| `s` | `3n−3` (`n≥2`); guarded `s≥1` | Q |
| self-bound gate | suppress while `E_int>Σ(n)` | Q |
| K1 / S1 increments | `−D_0(n)` / `+f_ret·D_0(n+1)` | Q / P (via U) |
| `Π(n)` | `λ(n)·f_ret·τ`; >1 shed / <1 freeze | (diagnostic) |

---

## 5. Config contract (new/renamed `SimConfig` fields landing in Phase B — MASS §11)

Per the rule-2 declared-but-unread exception, each field lands with its owning slice and
is removed from the exception table when its slice activates it. Validators reuse the
existing `check_drag_config` scaffolding (`config.py:~314+`) — load-time fail-loud guards,
not silent clamps.

- **ρ:** `helium_density_profile` (existing placeholder `config.py:220`) → selector
  `erf_complement(default)` / `tabulated`.
- **P:** `pickup_rate_coefficient` λ_0 (**rename/repurpose** `mass_rate_coefficient`),
  `pickup_rate_form` (existing `mass_rate_form`), `pickup_occupancy_cap ∈ {langmuir, none}`,
  `pickup_occupancy_exponent` p, `he_capture_velocity ∈ {at_rest, thermal}`.
- **Q:** `evap_rate_prefactor_per_ps` ν, `evap_rrk_dof` s (**guarded `s≥1`**),
  `evap_gate_onset_eV` (retired → optional override).
- **Validators:** `pickup_rate_form` and `pickup_occupancy_cap` enum-reject arms (mirror
  the `mass_scenario` reject arm); the `evap_rrk_dof` `s≥1` load guard.

**Rename note (build-time item, not done in this doc).** `config.py` currently exposes
`mass_rate_coefficient` / `mass_rate_form` (Tier-1 generic names). Phase B repurposes them
to `pickup_*`; keep a back-compat alias or do the rename atomically in the P build, and
record it in `drag_migration_log_tier2.md`. The `biphasic` vs `biphasic_energy_gated`
`MassScenario` naming reconcile stays a **Phase C** concern (Tier-2 plan §5) — note, do
not touch here.

---

## 6. Dependency / build order (within Phase B)

```
ρ (helium_density)            ── independent; built & tested first
        │
        ├──►  P (pickup; consumes ρ gate + mass_jump.capture + U-S1; mocks RNG)
        │
        └──►  Q (evaporation; consumes L gate/rungs + mass_jump.cold_shed + U-K1; mocks RNG)
```

ρ is fully independent. P and Q each compose Phase-A primitives (mocked for their own unit
tests, real for the composed checks) and the `mass_jump` resets, and each mocks the RNG. **P
and Q are mutually independent** (parallelizable) — they meet only in the Phase-C driver G
(one-event-per-step, shed-then-pickup order). No Phase-B slice depends on Phase C (schema /
driver) or Tier 3 (noise). Natural order: ρ → {P, Q}.

---

## 7. Testing methodology

Mirror the Tier-1a stochastic-unit discipline (validation step 6, unit level):
- **seeded/mock RNG** injected into every channel; **distribution/moment** asserts with
  **sample-size-justified tolerances** (Poisson mean/variance; Bernoulli rates);
- **identity / machine-precision** asserts for the capture & cold-shed reset defects and
  the S1/K1 increments;
- **bound/guarantee** asserts for the RRK rate range, the self-bound gate, the
  no-avalanche per-step cap, and the `s≥1` load guard;
- **shared-form parity** for ρ vs `drag.spatial_gate`;
- **No figures, no checkpoints, no real integrator, no production-sized runs** in pytest.

Suites: `tests/test_helium_density.py`, `tests/test_pickup.py`, `tests/test_evaporation.py`,
plus capture-reset asserts added to `tests/test_mass_jump.py`. Run the narrowest first,
then the full suite. Interpreter (Python may not be on PATH):

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

---

## 8. Acceptance criteria

**Per slice:** the slice's own suite green with neighbours (ρ/L/U) and the RNG mocked.

**Phase B overall:** full suite stays green (current baseline all-green on
`drag_implementation`); ρ/P/Q reproduce every §2 oracle to the stated tolerance; the
capture and cold-shed resets are exact and share the one reduced-mass formula; Poisson and
RRK moment/bound checks hold; all fail-loud guards (`s≥1`, enum rejects, mass checks) fire.

**Out-of-scope guard (fails review if it leaks in):** any Phase-B code path that runs the
integrator / O-step / jump-then-O seam, persists `E_int` *state*, bumps the checkpoint
schema, computes a 5-term closure, reads a noise amplitude, touches the Tier-0 drag law, or
*fits* λ_0/f_ret/ν/s/κ. Those are Phase C / Phase F / Tier 3.

---

## 9. Risks / notes carried

- **λ_0 is a Na⁺→I⁺ ±factor-2 prior, not a Phase-B fit target.** Central ~0.7–1.1/ps
  (GAH25 Rb⁺/Cs⁺; I⁺ is Rb⁺-like), *not* the 2.0/ps Na⁺ datum. Phase B encodes the form and
  exposes the knob; the size distribution arbitrates it at Tier 2 (CALIBRATION_MAP R1, R7).
- **Density-only ≡ Langmuir for the production ejection problem (cap inert).** The ion exits
  (`ρ_He→0`) before the shell saturates, so `Π` crosses 1 via density first (§6.11). The cap
  earns its keep only for the resting/slow-ion (A12) and the Tier-1 tail — kept for
  correctness and reviewer-defensibility, not because production needs it.
- **`n=1` direct dissociation is the boundary of validity of the statistical picture**
  (A11): a diatomic has a single vibrational mode, so RRK is degenerate (`s=3n−3=0`); the
  `k=ν` direct branch keeps the rate bounded through the last atom, and it is exactly where
  cold-shed neutrality (A8) is weakest — flagged, not papered over.
- **RNG draw order: specified here, locked at X.** Phase B documents the shed-then-pickup
  order for one-event-per-step bookkeeping (A13); the actual *lock* (forbidden-list:
  draw-order changes need explicit approval) lands with the schema/driver at Phase C.
- **ρ reuses the G2 erf gate; no G2→G4 promotion this phase.** The sourced baseline/TDDFT
  density profile is the declared fallback (CALIBRATION_MAP rows 5/8); promoting the drag
  spatial gate to G4 is explicitly out of Phase-B scope.
- **Mechanism locked, values open** — Phase B encodes *forms* and *sourced* anchors only;
  λ_0, f_ret, ν/s, κ, p stay knobs. **No fitting in Phase B** (that is Phase F).
