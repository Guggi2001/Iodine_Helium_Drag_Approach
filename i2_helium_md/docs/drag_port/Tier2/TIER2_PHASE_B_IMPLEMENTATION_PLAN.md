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
| `NU_EVAP_PER_PS` (ν) | 2.42 | ps⁻¹ | **added by Q** (constants anchor); **pinned** ([IHe05] curvature, A11) | Q |
| `s` = `effective_dof(n)` | **`n=2→4`** (linear), `3n−3` (`n≥3`); `n=1` direct | — | **Derived**, mode-counted (A11); guarded `s≥1` | Q |
| `λ_0` (pickup coeff.) | central ~0.7–1.1, ±factor-2 | ps⁻¹ | **Sourced+Bounded** (GAH25 Rb⁺/Cs⁺; *not* 2.0 Na⁺; CALIBRATION_MAP row 7) | P |
| `p` (occupancy exponent) | **fixed `1.0`** (not `p=κ`; tie is inverse) | — | held fixed; free at Phase F if size-dist demands (A12; ⚠ §3.2) | P |
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
*Dim:* `s = effective_dof(n)` dimensionless — **`n=2 → 4`** (linear 3-atom `3N−5`;
CALIBRATION row 10 "n=2 linear +1"), **`3n−3` for `n≥3`** (nonlinear full complex, A11);
`s−1 ≥ 1` throughout (n=2 → `s−1 = 3`); bracket `1−D_0/E_int` dimensionless (eV/eV); `k`
ps⁻¹, **bounded `k∈[0,ν)`** for `s>1` and `k=ν` at `n=1`; per-step shed prob.
`≤ 1−e^{−ν dt} ≈ ν·dt` (**no gate-open avalanche**). Below threshold (`E_int ≤ D_0(n)`) the
rate is exactly zero. Signed gate margin `G = E_int − Σ(n)` is the diagnostic readout.

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

**Interface (decided — refinement 2026-06-30, see §3.1).**
- `_erf_complement(depth, steepness) -> ndarray` — the **single source** for
  `½(1−erf(depth/steepness))`, **extracted** to a neutral location
  (**locked `physics/_gates.py`** — user 2026-07-01; a gate-specific neutral home so ρ
  never imports `drag.py`, preserving ρ's drag-independent role). Fail-loud on `steepness ≤ 0`.
  `drag.spatial_gate` is **rewired to call it** — a Tier-0 *no-behavior-change* refactor
  (identical output), guarded by a parity regression test (CLAUDE.md rule 1: one formula,
  one place; not a physics change to the locked drag law).
- `rho_he_ratio(depth, *, steepness) -> ndarray` — the erf-complement ratio in `[0,1]`,
  routed through `_erf_complement`. **No default steepness** (mirrors `drag.spatial_gate`,
  which has none); the driver passes **`_drag_gate_steepness(cfg)`** (the `simulation/ion.py`
  resolver — user 2026-07-01) so density and drag see the **one consistent surface** *across
  all gate modes*. NB: the resolver returns `cfg.potential_steepness` under the production
  `density_proportional`/`erf_tied` gate and `cfg.drag_gate_steepness` only under G3
  `erf_independent`; threading the **resolver** (not the raw `drag_gate_steepness` field)
  guarantees ρ tracks whatever surface drag actually uses — and matches CALIBRATION_MAP row 5
  (source = confining-potential steepness = `potential_steepness`). Module stays
  config-agnostic (Phase-A precedent).
- `tabulated_density_profile(...)` — interpolation machinery returning the same `[0,1]`
  ratio contract. **Built + round-trip tested in Phase B** on a hand-built profile (mirrors
  the Slice-L `tabulated_ladder` precedent); the **sourced baseline/TDDFT `ρ_He(r)` data
  array is a deferred rule-2 carry** (CALIBRATION row 8; pinned when the profile exists).
  **Interface locked (user 2026-07-01):** takes a hand-built `(depth_grid, ratio_grid)`,
  **linear interpolation on `depth`** (not `r`), round-trip exact at the grid nodes, and
  **clamps out-of-range `depth` to the tail values** (`1.0` inside for `depth < grid.min`,
  `0.0` outside for `depth > grid.max`) rather than raising — matching the erf-complement
  asymptotes so the two `HeliumDensityProfile` arms share the same `[0,1]` boundary contract.

**Encoded form.** §2.2 (ρ).

**Knobs (config §6):** **new** `HeliumDensityProfile = Literal["erf_complement","tabulated"]`
and field `helium_density_profile: HeliumDensityProfile = "erf_complement"` — **repurpose**
the existing `Optional[object]` placeholder (`config.py:261`; its `# future G4 density
profile` comment is **rewritten** — ρ stays **G2**, this is the surface-density gate, not a
G2→G4 promotion). Reject-arm guard mirrors `mass_scenario` / `check_ladder_config`. Steepness
is **not** a new ρ field — the driver passes `_drag_gate_steepness(cfg)` (14.2 Å; the
resolver, = `potential_steepness` in production — see the interface note above).

**Oracle values.** ratio = 1 at `depth≪0`, **0.5 at `depth=0`**, →0 at `depth≫0`;
monotone non-increasing in `depth`; identical to `drag.spatial_gate(depth, 14.2)` (both
route through `_erf_complement`, so equality is by construction — the parity test is a
regression lock that neither re-inlined the formula, not a two-forms drift check).

**Independence.** Depends only on `depth` + steepness. Built and tested first; P mocks it.

**Test spec (`tests/test_helium_density.py`, + a parity assert added to
`tests/test_drag.py`).**
- ratio ∈ [0,1], =0.5 at depth 0, monotone non-increasing; →1 / →0 at the asymptotes
  (erf saturates *exactly* to ±1 in float at large `|depth|`, so assert exact `1.0`/`0.0`
  at the far tail — no epsilon needed, per the Slice-K tail-saturation precedent).
- **Parity / single-source regression:** `rho_he_ratio` equals `drag.spatial_gate` on a
  `depth` grid to machine precision (locks that both still route through `_erf_complement`).
- Vectorized scalar-in/array-in shape parity; `_erf_complement` (and through it both
  callers) fail-loud on `steepness ≤ 0`.
- Tabulated machinery round-trips a hand-built `[0,1]` profile.

**Acceptance.** All pass; `rho_he_ratio` == `spatial_gate`; both route through the extracted
`_erf_complement`; `spatial_gate` output unchanged (Tier-0 parity); tabulated round-trips.

#### §3.1 Slice ρ refinement decisions (user, 2026-06-30)

Cross-checked against MASS §4/§5 + dim-table and CALIBRATION_MAP rows 5/8 — no
contradictions (row 5 = erf-tied **G2** gate, source 14.2 Å; row 8 = Sourced ρ_He profile =
deferred fallback). Docs-only; `[PROCEED TO IMPLEMENTATION]` boundary holds.

1. **Single-source erf via an extracted `_erf_complement` helper** (not delegate-to-drag,
   not an independent mirror). Both `drag.spatial_gate` and `rho_he_ratio` route through it.
   Truest rule-1; keeps ρ from importing the drag law (preserves the ρ-independent role);
   the `spatial_gate` rewire is a **Tier-0 no-behavior-change refactor** guarded by a parity
   regression test (build-time edit to a locked module, output identical).
2. **No default steepness** on `rho_he_ratio` — caller passes it, exactly as
   `drag.spatial_gate` does. Driver supplies **`_drag_gate_steepness(cfg)`** (14.2 Å) — the
   `simulation/ion.py` resolver, **not** the raw `drag_gate_steepness` field — so density and
   drag share the **identical surface across all gate modes** (the resolver returns
   `potential_steepness` under the production `density_proportional`/`erf_tied` gate,
   `drag_gate_steepness` only under G3 `erf_independent`). This matches CALIBRATION_MAP row 5
   (source = confining-potential steepness = `potential_steepness`). **Corrected 2026-07-01**
   (was "`cfg.drag_gate_steepness`", which only matches the drag surface in the G3 branch or
   while the two fields happen to be equal). No new `POTENTIAL_STEEPNESS_ANGSTROM` constant
   (avoids a 4th home for 14.2; `potential_steepness` / `potential_steepness_molecule` /
   `drag_gate_steepness` already exist).
3. **Explicit `HeliumDensityProfile` Literal selector** (`erf_complement` default /
   `tabulated`) with an enum reject-arm guard, repurposing the `Optional[object]`
   placeholder (`config.py:261`) and rewriting its stale `# future G4` comment. The
   tabulated **machinery is built + round-trip tested** in Phase B; only the **sourced
   TDDFT data array** is deferred (rule-2). Tabulated profile data is passed as a **module
   arg**, not a config field, in Phase B.

**Drift fixed by this refinement:** stale line ref (`config.py:220`→`:261`); misleading
`# future G4 density profile` field comment (ρ stays G2); non-existent
`POTENTIAL_STEEPNESS_ANGSTROM` default; the §2.1/interface "tabulated not built vs.
round-trip tested" tension (machinery built, data deferred); and the self-contradictory
"shared helper" vs "two forms diverging" parity-test framing (now a single-source
regression lock).

---

### Slice P — Pickup channel + `mass_jump` capture reset *(stateful primitive; mock RNG)*

**Modules.** `physics/pickup.py` (the channel) and `physics/mass_jump.py` (**extended**
with the capture reset).

**Purpose.** The gain channel: one independent Bernoulli draw per ion per step; on fire,
apply the momentum-conserving **capture** reset and report the S1 heat. The reset is the
symmetric +He counterpart of `cold_shed` and lives in `mass_jump.py` so the reduced-mass
defect coefficient stays in one place.

**Interface (decided — refinement 2026-07-01, see §3.2).**
- *In `mass_jump.py`:* `capture(v_minus, m_minus_amu, *, m_he_amu=MASS_HE_AMU, u_he=0.0)`
  → a **new `CaptureResult`** dataclass `(v_plus, m_plus_amu, dE_mass_transfer)` with
  `v_plus = m/(m+m_He)·v_minus`, `m_plus_amu = m + m_He` (a **gain**),
  `dE_mass_transfer = +½·(m·m_He)/(m+m_He)·|v⁻|²` (**> 0**, into `E_mass_transfer`); plus the
  `*_components` vectorized form mirroring `cold_shed_velocity_components`. A distinct
  `CaptureResult` (not `ShedResult`) is used so the type name matches the gain semantics.
  **Shared-coeff refactor:** `_reduced_mass_defect_coeff` is changed to return the **unsigned
  magnitude** `0.5·(m·m_He)/m_plus`; `cold_shed` applies the `−`, `capture` the `+` — one
  formula, one place (rule 1), signs at the call sites. `cold_shed` output stays byte-identical
  (regression test locks it).
- *In `pickup.py`:*
  - `attach_probability(lambda_attach, dt_ps) -> P_attach` (`1−e^{−λ dt}`).
  - `lambda_attach(rho_ratio, n, *, lambda0, n_star, p, cap)` — the rate with the Langmuir
    factor (`cap="langmuir"` applies `(1−n/n*)_+^p`; `"none"` recovers density-only).
  - `pickup_step(state, *, rng, rho_ratio, lambda0, f_ret, picture, kappa, ...)` — draw one
    Bernoulli; on fire compose `mass_jump.capture` (reset/mass/defect) +
    `internal_energy_budget.dE_int_pickup_eV(n, f_ret=…, picture=…, kappa=…)` (S1 heat, `n` =
    **pre-pickup**); return post-event `(n', m', v', ΔE_int, ΔE_mass_transfer, fired)`. Consumes
    ρ (gate), U (S1). **`picture`+`kappa` are threaded** (both required by U's `d0_of_n`).
  - **`pickup_rate_form`: only `density_only` is built**; `sweeping`/`dwell_time` are
    **declared-but-unread rule-2 arms** (enum present; raise a clear `NotImplementedError`
    until Tier-1/2 demands them — MASS §5 locks density-only, CALIBRATION R7). Same for
    `he_capture_velocity`: `at_rest` built, `thermal` a rule-2 arm.

**Encoded form.** §2.2 (P).

**Knobs (config §6):** `pickup_rate_coefficient` λ_0 (**atomic rename**, no alias, of
`mass_rate_coefficient` `config.py:259` — zero readers, verified); `pickup_rate_form ∈
{density_only(default), sweeping, dwell_time}` (rename of `mass_rate_form` `config.py:258`);
`pickup_occupancy_cap ∈ {langmuir(default), none}`; **`pickup_occupancy_exponent` p — fixed
default `1.0`, NOT `p=κ`** (the A12 p↔κ tie is **physically inverse**: rigid shell = large κ =
sharper cutoff = *smaller* p, so a literal `p=κ` inverts it; free p at Phase F only if the
size-dist first-shell edge can't be met with p=1 + κ, and if ever tied, match cutoff *slopes*,
not values); `he_capture_velocity ∈ {at_rest(default), thermal}` (thermal deferred — Tier 3;
`at_rest` = `u_he=0`) — **declared here (owns the `capture` primitive) but declared-but-unread
in Phase B**: the `capture` reset takes `u_he` as a param, and the config field is read by the
**Phase-C driver** (the f_int/f_ret precedent). **Retire** the dead `mass_relaxation_tau_ps`
(`config.py:260`, zero readers, superseded by the Slice-K `internal_energy_cooling_tau_ps`).

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
- `capture` reset returns a `CaptureResult`: `v⁺`, `m⁺ = m+m_He`, defect exact and **> 0**;
  defect **sign opposite** to `cold_shed`; fail-loud on non-finite / non-positive masses
  (reuse `_check_masses`).
- **Coeff-refactor regression:** `_reduced_mass_defect_coeff` returns the unsigned magnitude;
  `cold_shed` output is **byte-identical** to the delivered Tier-1a values (lock — no physics
  drift from moving the sign to the call site).
- `attach_probability` matches `1−e^{−λ dt}`; `→ λ dt` for small `λ dt`.
- Poisson moment check (seeded, sample-size tolerance); independence across ions.
- Langmuir cap →0 at `n*`; density-only/langmuir equivalence for `n ≪ n*`; **p=1 default**
  gives `(1−n/n*)` linear.
- S1 heat sign + magnitude via the mocked U stub (`picture`/`kappa` threaded).
- `pickup_rate_form ∈ {sweeping, dwell_time}` and `he_capture_velocity="thermal"` raise
  `NotImplementedError` (rule-2 arms).

**Acceptance.** Momentum/defect exact; Poisson moments within the sample-size band; cap →0
at `n*`; S1 sign correct.

---

### Slice Q — Evaporation channel *(stateful primitive; mock RNG; reuse `cold_shed`)*

**Module.** `physics/evaporation.py`

**Purpose.** The loss channel: the parameter-free self-bound gate, then the saturating RRK
shed of the top rung; on fire reuse `mass_jump.cold_shed` and report the K1 drain.

**Interface (decided — refinement 2026-07-01, see §3.2).**
- `rrk_rate(E_int_eV, n, *, nu, picture, kappa) -> k_per_ps` — the RRK rate with the
  `effective_dof(n)` bracket for `n≥2`, the **`n=1` direct `k=ν`** branch, and `k=0` below
  threshold; consumes L's `d0_of_n`.
- `shed_probability(k_per_ps, dt_ps) -> P_shed` (`1−e^{−k dt}`).
- `is_self_bound(E_int_eV, n, *, picture, kappa) -> bool` — the gate `E_int > Σ(n)`
  (suppress while True); consumes L's `ladder_cumsum`.
- `gate_margin_eV(E_int_eV, n, *, picture, kappa) -> ndarray` — the **signed diagnostic**
  `G = E_int − Σ(n)` (`is_self_bound` is `G > 0`). Pure/cheap; the mechanism's natural
  readout. Q exposes the helper here; the *trajectory record* of `Σ(n,t)` / `G(t)` (alongside
  `t×` and `Π`) is written by the Phase-C driver / Phase-E diagnostics, **not** in Phase B.
- `evaporation_step(state, *, rng, nu, picture, kappa, ...)` — if self-unbound, return
  no-shed; else draw one Bernoulli, and on fire compose `mass_jump.cold_shed` (reset/mass/
  defect) + `internal_energy_budget.dE_int_shed_eV(n, picture=…, kappa=…)` (K1 drain `−D_0(n)`,
  `n` = **pre-shed**); return post-event tuple.
- `effective_dof(n)` — **`n==2 → 4`** (linear 3-atom `3N−5`; I⁺He₂ = ion + 2 He = 3 atoms),
  **`3n−3` for `n≥3`** (nonlinear full complex), `n=1` is the direct branch (no bracket).
  The **config-load `s≥1` guard** applies to any effective-scalar override (rejects `s<1`, the
  divergent-rate regime).

**Encoded form.** §2.2 (Q).

**Knobs (config §6):** `evap_rate_prefactor_per_ps` ν (Sourced, pinned 2.42) defaulting to a
**new `NU_EVAP_PER_PS = 2.42` `constants.py` anchor** (the Slice-K `|S|` single-source
precedent — the plan's "ν lands with Phase A" was inaccurate; Phase A added only `N_STAR` /
`D0_*` / `D_FLOOR` / `S_ABS_EV`, so **Q adds `NU_EVAP_PER_PS`**). `evap_rrk_dof` s
(`Optional[float] = None`: None → per-`n` `effective_dof`; a set float → effective-scalar
override, **guarded `s≥1`** *when set*, the f_int=None precedent). **Gate onset (3 facets):**
(1) **no calibration knob** — the gate is `Σ(n)`, parameter-free (MASS R9), stays **Derived**;
(2) the signed **gate margin `G`** is a *diagnostic output* (helper above; recorded at C/E);
(3) an optional **`gate_onset_override_eV: Optional[float] = None`** for sensitivity/
falsification — None → computed `Σ(n)`; a float forces a fixed threshold **only behind a loud
provenance guard** (mirrors `allow_unvalidated_binding_pairing`: config-load **refuses unless**
an explicit `allow_gate_onset_override=True`), so it can never silently enter a production /
Tier-2-lock run (R10 diagnostic-lever philosophy, not a knob).

**Oracle values.**
- `k ∈ [0, ν)` for `s>1`; `k=ν` exactly at `n=1`; `k=0` for `E_int ≤ D_0(n)`.
- **`effective_dof`: `n=2 → 4`** (linear), **`n=3 → 6`, `n=21 → 60`** (nonlinear `3n−3`);
  `n=1` uses the direct branch. So the `n=2` bracket exponent is `s−1 = 3` (was 2 under the
  uniform count).
- **Gate suppresses all sheds** while `E_int > Σ(n)` (net self-unbound) — zero fires
  regardless of RNG. `gate_margin_eV = E_int − Σ(n)` is `< 0` exactly when self-bound.
- **No avalanche:** per-step shed probability `≤ 1−e^{−ν dt}` even as `E_int→∞`
  (saturation), so at most ~`ν·dt` per step; **one rung per step**.
- **`s≥1` guard fires** on an override `s<1` at config-load (clear error, no silent clamp).
- **Gate-onset override guard fires:** a non-None `gate_onset_override_eV` without
  `allow_gate_onset_override=True` is refused at config-load (mirrors the binding-pairing guard).
- K1 drain `= −D_0(n)` exactly (via the mocked U stub).

**Independence.** Mocks L (`d0_of_n`, `ladder_cumsum`), U (`dE_int_shed_eV`), and the RNG
with known stubs. Reuses the *real* `mass_jump.cold_shed` (already delivered + tested) for
the reset. No integrator, no `E_int` state.

**Test spec (`tests/test_evaporation.py`, + config guards in `tests/test_*_config.py`).**
- `rrk_rate` matches the oracle: bounded `[0,ν)`, `=ν` at `n=1`, `=0` below threshold,
  monotone increasing in `E_int`.
- Gate: `is_self_bound` flips at `E_int = Σ(n)`; suppresses all sheds while self-unbound;
  `gate_margin_eV = E_int − Σ(n)` sign + zero-crossing.
- Avalanche guard: per-step prob `≤ 1−e^{−ν dt}` across a wide `E_int` sweep; exactly one
  rung drawn per fire.
- `effective_dof`: **`n=2 → 4`**, `n=3 → 6`, `n=21 → 60`; the config-load `s≥1` guard rejects
  a set override `s<1` (load-time fail-loud), no-op when None.
- **`gate_onset_override_eV` guard:** non-None without `allow_gate_onset_override` refused at
  config-load; None → parameter-free `Σ(n)`.
- K1 drain sign + magnitude via the mocked U stub; shed reset delegates to `cold_shed`
  (no second copy of the reduced-mass formula).
- `NU_EVAP_PER_PS = 2.42` anchor exists; `evap_rate_prefactor_per_ps` defaults to it.

**Acceptance.** Rate matches oracle; gate + avalanche guarantees hold; `s≥1` guard fires;
gate-onset override guard fires; K1 drain exact; reset reuses `cold_shed`.

#### §3.2 Slice P/Q refinement decisions (user, 2026-07-01)

Cross-checked against MASS §4/§5 + dim-table, §11 config, A11/A12, and CALIBRATION rows
7/8/10/13. One wording discrepancy surfaced (see ⚠ below). Docs-only; the
`[PROCEED TO IMPLEMENTATION]` boundary holds.

**Slice P**
1. **`CaptureResult` dataclass** (new), not a reused `ShedResult` — the type name matches the
   +He *gain* (`m⁺ = m+m_He`, `dE_mass_transfer > 0`). `_reduced_mass_defect_coeff` is
   refactored to the **unsigned magnitude** `0.5·(m·m_He)/m_plus`; `cold_shed` applies `−`,
   `capture` `+` (rule 1, one formula; `cold_shed` byte-identical, regression-locked).
2. **`pickup_occupancy_exponent` p = fixed `1.0`, not `p=κ`.** ⚠ MASS §11 says "p default
   *tied to* κ," which reads as `p=κ`. The tie is **physically inverse** (rigid shell = large
   κ = sharper cutoff = *smaller* p), so a literal `p=κ` is wrong. Decision: **hold p=1** for
   calibration; free p at Phase F only if the size-dist first-shell edge can't be met with
   p=1 + κ (and if ever tied, match cutoff *slopes*). **MASS §11 + CALIBRATION row 7 wording
   flagged for a clarifying annotation** (not silently rewritten).
3. **Rename is atomic, no alias; `mass_relaxation_tau_ps` retired.** `mass_rate_form` /
   `mass_rate_coefficient` / `mass_relaxation_tau_ps` have **zero readers** in the package
   (grep-verified) → `→ pickup_*` atomically; the τ field is dead (Tier-2 cooling is the
   Slice-K `internal_energy_cooling_tau_ps`) → removed.
4. **`density_only` + `at_rest` built; `sweeping`/`dwell_time`/`thermal` are rule-2 arms**
   (raise `NotImplementedError`). MASS §5 locks density-only; CALIBRATION R7 defers the rest.

**Slice Q**
5. **`effective_dof`: `n=2 → 4`** (linear 3-atom `3N−5`, CALIBRATION row-10 "n=2 linear +1"
   applied), `3n−3` for `n≥3`, `n=1` direct `k=ν`. The uniform `3n−3` is no longer used at n=2.
6. **Gate onset — 3 facets:** (a) **no calibration knob** (Derived, gate = `Σ(n)`,
   parameter-free, MASS R9); (b) **diagnostic** signed margin `gate_margin_eV = E_int − Σ(n)`
   exposed as a pure Q helper (trajectory record of `Σ(n,t)`/`G(t)` alongside `t×`/`Π` is
   Phase C/E); (c) **`gate_onset_override_eV: None`** with a **loud provenance guard**
   (`allow_gate_onset_override`, mirroring `allow_unvalidated_binding_pairing`) so a forced
   threshold can never silently enter production (R10 diagnostic-lever).
7. **`NU_EVAP_PER_PS = 2.42` is added by Q** (constants.py anchor + `evap_rate_prefactor_per_ps`
   defaulting to it) — the plan's "ν lands with Phase A" was inaccurate; Phase A shipped no
   evap prefactor. `evap_rrk_dof` lands `Optional[float] = None` (per-`n` default; guarded
   `s≥1` only when set).

**Drift fixed:** stale config line refs (`:217/:218` → `:258/:259`); the "ν lands with Phase
A" gap; the reused-`ShedResult`-for-a-gain semantic mismatch; the `p=κ` ambiguity.

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
| Langmuir cap | `(1−n/n*)_+^p`; **p=1 default**; →0 at `n=n*=21` | P |
| RRK rate | `k∈[0,ν)` (`s>1`); `k=ν` (`n=1`); `k=0` (`E_int≤D_0`) | Q |
| `ν` | 2.42 ps⁻¹ (pinned; new `NU_EVAP_PER_PS` anchor) | Q |
| `s` = `effective_dof(n)` | **`n=2→4`** (linear), `3n−3` (`n≥3`); `n=1` direct; guarded `s≥1` | Q |
| self-bound gate | suppress while `E_int>Σ(n)`; margin `G=E_int−Σ(n)` | Q |
| K1 / S1 increments | `−D_0(n)` / `+f_ret·D_0(n+1)` | Q / P (via U) |
| `Π(n)` | `λ(n)·f_ret·τ`; >1 shed / <1 freeze | (diagnostic) |

---

## 5. Config contract (new/renamed `SimConfig` fields landing in Phase B — MASS §11)

Per the rule-2 declared-but-unread exception, each field lands with its owning slice and
is removed from the exception table when its slice activates it. Validators reuse the
existing `check_drag_config` scaffolding (`config.py:~314+`) — load-time fail-loud guards,
not silent clamps.

- **ρ:** new `HeliumDensityProfile = Literal["erf_complement","tabulated"]`; field
  `helium_density_profile` **repurposes** the `Optional[object]` placeholder (`config.py:261`,
  stale `# future G4` comment rewritten — ρ stays **G2**) → default `erf_complement`. Steepness
  is **not** a new ρ field; the driver passes `_drag_gate_steepness(cfg)` (14.2 Å; the resolver,
  = `potential_steepness` in production — CALIBRATION row 5). See §3.1.
- **P:** `pickup_rate_coefficient` λ_0 (**atomic rename** of `mass_rate_coefficient`
  `:259`), `pickup_rate_form` (rename of `mass_rate_form` `:258`), `pickup_occupancy_cap ∈
  {langmuir, none}`, `pickup_occupancy_exponent` p (**fixed default `1.0`**, not `p=κ`),
  `he_capture_velocity ∈ {at_rest, thermal}`. **Retire** `mass_relaxation_tau_ps` (`:260`,
  dead).
- **Q:** `evap_rate_prefactor_per_ps` ν (defaults to the **new** `NU_EVAP_PER_PS = 2.42`
  `constants.py` anchor added by Q), `evap_rrk_dof` s (`Optional=None`; **guarded `s≥1`**
  when set), `gate_onset_override_eV` (`Optional=None`, **parameter-free by default**) +
  `allow_gate_onset_override` (loud provenance guard).
- **Validators:** `helium_density_profile`, `pickup_rate_form`, `pickup_occupancy_cap`, and
  `he_capture_velocity` enum-reject arms (mirror `mass_scenario` / `check_ladder_config`); the
  `evap_rrk_dof` `s≥1` load guard (when set); the `gate_onset_override_eV` ↔
  `allow_gate_onset_override` provenance refuse (mirrors `allow_unvalidated_binding_pairing`).

**A13 integrator-policy fields are Phase C, not Phase B (MASS §11 reconcile).** The A13 knobs
`mass_jump_velocity_reset ∈ {momentum_conserving, label_only}`, `one_mass_event_per_step`, and
`jump_o_step_ordering` are **driver/integrator policy** and land + activate with **Slice G
(Phase C)** — see `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §4 "Integrator (G)". They are **not**
Phase-B config fields: Phase-B `capture`/`cold_shed` are momentum-conserving **by
construction**, so `label_only` (the non-closing §6-invariant diagnostic) is a **Phase-C driver
branch** Phase B never implements; `one_mass_event_per_step` / `jump_o_step_ordering` compose
the two channels + the O-step, which do not exist until G. Only `he_capture_velocity` is
declared here (Slice P, owning `capture`) and *activated* at C — the one field shared across the
B/C boundary, resolved by the declare-at-P / read-at-C split above.

**Rename note (build-time).** `mass_rate_coefficient` / `mass_rate_form` /
`mass_relaxation_tau_ps` have **zero readers** (grep-verified) → **atomic rename to
`pickup_*`, no back-compat alias**, and **`mass_relaxation_tau_ps` removed** (superseded by
Slice-K `internal_energy_cooling_tau_ps`); record in `drag_migration_log_tier2.md`. The
`biphasic` vs `biphasic_energy_gated` `MassScenario` naming reconcile stays a **Phase C**
concern (Tier-2 plan §5) — do not touch here.

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
