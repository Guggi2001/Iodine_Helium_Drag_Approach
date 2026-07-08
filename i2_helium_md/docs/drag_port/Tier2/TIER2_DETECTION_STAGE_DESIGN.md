# Tier 2 — Detection-Time Continuation Stage (event-driven post-ejection cascade)

> **Boundary.** This is a *design document*, not code. The strict
> Physics-Definition / Software-Implementation boundary holds: no Python until
> the explicit `[PROCEED TO IMPLEMENTATION]` trigger. Module descriptions are
> *interface contracts*; the validation plan is a *specification*.
>
> **Status:** DESIGN APPROVED in discussion (user, 2026-07-07) — the four
> design forks (§1) were adjudicated by the user. The Sourced input is
> **supplied (user, 2026-07-07): t_detect = 8.53 µs = 8.53·10⁶ ps**, from the
> experimental-setup publication (CALIBRATION_MAP row 24). **Amended
> 2026-07-07 (post-approval, user):** decision 5 (§1) — E2 retained but
> re-scoped (bridge + diagnostic), gated caps shortened by config, and the
> stage made **skippable** (direct `ion.npz` seeding under the §2.1 guard);
> one connected code review covers DS + the E2 re-frame. All §7 pre-trigger
> items are resolved (slice placement: standalone pre-F5 **Slice DS**).
>
> **DELIVERED 2026-07-07** under the `[PROCEED TO IMPLEMENTATION]` trigger
> (delivery record: `drag_migration_log_tier2.md`, "Slice DS DELIVERED").
> **As-built deviations (all recorded in place below):**
> (1) **"v unchanged" was wrong** — the delivered `mass_jump.cold_shed`
> applies the momentum-conserving kick `v⁺ = m/(m−m_He)·v` per fire (He
> leaves cold; the KE rise is cancelled by the booked `E_mass_transfer`
> defect). The (E_int, n) jump chain never reads `v`, so the exactness
> claim is untouched; §2.2/§2.4/§2.5 corrected in place.
> (2) The P3 bound's density factor is **gate-conditional**: the local
> erfc ρ/ρ_bulk under `density_scaled`, exactly **1 under `none`** (ungated
> cooling acts everywhere — a live ungated ion can never hand over).
> (3) `check_detection_config` additionally requires **ν > 0** (with ν = 0
> in-band ions sit at k = 0 without being frozen/suppressed — the §2.3
> taxonomy would be unsound) and bounds `detection_time_ps` beyond the
> nominal seed-window end; the stage re-checks against the realized t_h.
> (4) `DetectionResult`/`detection.npz` carry, beyond the §3.2 list, the
> terminal velocities and the per-ion terminal E_kin/E_pot/E_dissip/
> E_mass_transfer, so the 5-term closure is checkable across the stage
> boundary from the artifact alone (fork-4 spirit).
> Constants as built: `EPS_DRAIN = 1e-6`, `DETECTION_STREAM_KEY = 0xD5_2026`.
>
> **CONNECTED CODE REVIEW EXECUTED + FIXES APPLIED (2026-07-07/08)** — the
> §1-item-5 review ran (multi-agent, high effort; 10 confirmed findings,
> none refuted). The three correctness fixes are folded into this document
> in place: the **review-hardened P1–P3 guard** (§2.1 NB — cooling bound
> extended to suppressed ions; a helium-exposure bound for *all* ions incl.
> frozen, defending P1/P2 on the skip path), the **RRK-underflow frozen
> lane** (§2.3 — previously a crash), and the **report's skip-path support**
> (`relaxation.npz` optional, plus a stale-`detection.npz` coherence guard).
> Cleanups: point-of-use `tabulated` refusals; bit-exact `fold = −dE_int`;
> shared seed-guard + stage-RNG helpers in `simulation/checkpoint.py`
> (E2 rewired to them, zero behavior change); E1 doc-drift fixes. Delivery
> record: `drag_migration_log_tier2.md` "Slice DS code review". Suite after
> fixes: 2176 passed.
>
> **Entry docs:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4b/§6 (I11: the gated
> terminal read is flight-time dependent — the problem this stage resolves);
> `drag_migration_log_tier2.md` ("Wave-5 gated landing re-location EXECUTED").
> **Parent:** `TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (the F2/F5 observables
> this stage makes well-defined).

---

## 0. Problem statement

The Tier-2 arbitration observable is the **terminal** I⁺Heₙ size distribution
(`data/reference/integrated_i_he_abundance.csv`). Wave 5 established (I11)
that under the intended production arm (`cooling_spatial_gate =
density_scaled`) the cascade does **not** converge within any tractable
fixed-dt window at the staircase-landing s_eff: at s_eff = 30 the 1000 ps
relaxation cap ends with `frac_frozen = 0`, n̄ still descending (13.80 →
6.27), the energetic floor n ≈ 2 unreached. The relaxed read is then a
*snapshot of a live cascade*, not a terminal — cap-artifactual.

Fixed-dt integration cannot reach the experimental timescale structurally:
the one-event-per-step bias guard `ν·dt ≤ 0.1` with ν = 2.42 ps⁻¹ (Sourced,
`NU_EVAP_PER_PS`) caps `dt_relax ≤ 0.041 ps`, so a ~8.5 µs flight needs
≥ 2·10⁸ steps per run. `check_relaxation_config` already records that no
sourced experimental flight time was available (MASS §R5: "hundreds of ps").

This stage resolves the observable by **changing the integration scheme, not
the physics**: in the post-ejection regime the locked mass mechanism reduces
exactly to a Markov jump chain, which an event-driven solver integrates to
*any* detection time in ≤ ~21 events per ion.

## 1. Design decisions (user, 2026-07-07)

Four forks were presented with alternatives; the adjudicated choices:

1. **Detection-time semantics: Sourced TOF flight time.** The distribution is
   read at detector arrival; `t_detect` becomes a **Sourced** constant
   (apparatus-derived: flight geometry/fields → CALIBRATION_MAP row).
   **Resolved 2026-07-07:** t_detect = 8.53 µs, from the experimental-setup
   publication (user-confirmed accurate). Rejected alternatives: an earlier
   mass-freeze point (extraction region); an unsourced Bounded pin with a
   sensitivity band (the band survives as a free report-side read, §3.3).
2. **Realization: stochastic event-driven (Gillespie).** Exponential waiting
   times per shed; each ion yields an integer n — drop-in for E1's per-ion
   `terminal_n` and the Wasserstein comparison; per-event ledger booking
   stays possible. The exact hypoexponential closed form becomes the *pytest
   oracle*, not a production arm (rejected: analytic probability-vector
   observable — changes the E1 input type, no per-ion realization; rejected:
   both-behind-an-enum — duplicate physics for a check a test oracle covers).
3. **Architecture: additive stage after E2.** The delivered
   `relaxation_stage` runs unmodified (to `relaxation_time_ps` or
   all-frozen); the detection stage seeds from its final state, with an
   exactness guard (§2.1) at handover. Rejected: shortening E2's contract by
   design; integrating event-driven logic into the reviewed E2 loop.
   Campaign economy (a shorter E2 cap once the continuation owns the tail)
   stays a later *config* decision, not a design fork.
4. **Bookkeeping: full per-event ledger.** Each shed books its K1 drain,
   e_bind fold, and cold-shed mass-transfer defect, so the 5-term invariant
   closes across the new stage; per-ion event times are stored (enables the
   detection-time sensitivity re-read without re-running). Rejected: a
   minimal `(n, E_int)`-only terminal read (breaks the invariant chain at
   the stage boundary; loses event-time diagnostics).
5. **E2 relationship (post-approval amendment, user, 2026-07-07 — amends
   fork 3's "additive after E2" to "additive after E2 *or* the ion stage").**
   Adjudicated after an explicit "is E2 dead code?" investigation of
   `relaxation_stage.py`:
   - **E2 is retained, not deleted.** It is the only correct integrator
     while K2 cooling is live (DS's exactness requires E_int constant
     between sheds — cooling breaks exactly that): the entire ungated
     `none` arm (a working-method-mandated sensitivity leg; ions arrive at
     ion-end neither frozen nor P3-compliant) and any in-/near-bubble ion
     at handover. It also produced the 102 existing probe artifacts Wave 7
     re-reads.
   - **E2's contract is re-framed** from "reach the terminal by matched-time
     integration" (the R5 mission, unachievable on the gated arm per Wave 5)
     to: ungated — terminal solver to freeze-out (DS no-ops); gated — a
     short **handover bridge** to a P3-clean state. The relaxed read is
     demoted to a convergence diagnostic (§3.5). Docstring mission statement
     updated in the slice.
   - **Gated caps are shortened by config** (future-facing USER SETTINGS,
     e.g. ~10 ps instead of 1000 ps; the §2.1 guard fails loud if shortened
     past ejection/decoupling). Existing probe dirs are untouched.
   - **E2 becomes skippable.** With `relaxation_stage_enabled=False`, DS
     seeds **directly from `ion.npz`**; the §2.1 P1–P3 guard is the sole
     defense (its failure remedy on this path: enable E2 / extend the cap).
     Both artifacts are v7 `IonCheckpoint`s, so this is one seeding
     contract, not a second code path. On the skip path there is no
     `frac_frozen` / relaxed diagnostic — the DS `state_reason` fractions
     take over.
   - **One connected code review** covers the DS build + the E2 re-frame
     together (review question: is the stage boundary and the demoted E2
     contract coherent).

## 2. Physics definition

### 2.1 Exactness preconditions (the handover guard)

The stage is an **exact solver** of the delivered mechanism valid iff, at the
handover state (the seed checkpoint's final column — E2's, or the ion
stage's on the skip path (§1 item 5) — time t_h), for every ion:

> **Review-hardened (code review, 2026-07-07/08).** The approved draft
> exempted every `k = 0` (permanent-state) ion from the bound and checked
> only cooling. The connected review confirmed two holes: (a) **suppression
> is cooling-reversible** (K2 drains E_int below Σ(n) and re-opens the gate
> — `biphasic_step` runs cooling *before* the gate for exactly this reason),
> so suppressed ions need the cooling bound like live ions; (b) pickup and
> drag are only *omitted* by the stage, not proven absent — a frozen
> in-bubble ion has live pickup (each capture re-heats E_int by `f_ret·D₀`,
> which can **unfreeze** the cascade) and live drag. The as-built guard
> below reflects the fix.

- **(P1/P2)** pickup and drag dead **by position**, for **every** ion
  (frozen included): `ρ_He(r_i)/ρ_bulk · max(λ₀, 1/τ) · (t_detect − t_h) ≤
  ε_drain`. The stage never calls those channels; this bound is what makes
  the omission exact — it is the skip path's actual defense.
- **(P3)** cooling negligible over the whole remaining flight, for every
  ion **not at the energetic floor** (live in-band *and* suppressed ions;
  only `n = 0` / `E_int < D₀(n)` ions are exempt — cooling cannot unfreeze):

  `ρ_cool(r_i) · (t_detect − t_h)/τ ≤ ε_drain`

  with `ρ_cool` the local erfc density under `density_scaled` and exactly
  **1 under `none`** (ungated cooling acts everywhere — a non-frozen ungated
  ion can never hand over; it must freeze in E2). `ε_drain` is a documented
  module constant (order 10⁻⁶ — a *defensive bound*, not a knob: gated
  ejected ions have erfc-suppressed ρ/ρ_bulk that underflows to 0, so the
  guard passes with enormous margin or fails loudly on a genuinely wrong
  input). Dimensionless: [1]·[ps]/[ps] = 1. The P3 bound is the fractional
  E_int error committed by treating cooling as zero (|ΔE_int|/E_int ≤
  ε_drain), t_detect-aware by construction.

Violations raise with the per-ion list and the remedy (extend
`relaxation_time_ps`; on the skip path, enable the relaxation stage; or the
run is on the ungated arm and must freeze in E2 — ungated runs arrive
`frac_frozen = 1` and the stage no-ops per ion).

Under P1–P3, E_int is **constant between sheds**: the per-ion state (E_int, n)
evolves only by discrete shed events at the delivered rate (`rrk_rate`
verbatim, all gating encoded):

```
k(E_int, n) = ν·(1 − D₀(n)/E_int)^(s−1)   for D₀(n) < E_int < Σ(n), n ≥ 2
            = ν                            for n = 1, E_int > D₀(1)  (direct)
            = 0                            otherwise (frozen / suppressed / n ≤ 0)
```

with s = `effective_dof(n)` or the `cfg.evap_rrk_dof` override — whatever the
run config says; the stage adds no rate physics.

### 2.2 Event loop (per ion, exact)

```
t ← t_h
loop:
    k ← rrk_rate(E_int, n)                         [ps⁻¹]
    if k == 0:            permanent state — record reason, stop
    Δt ← −ln(u)/k,  u ~ U(0,1)                     [ps]
    if t + Δt > t_detect: time exhausted — stop (n unchanged)
    t ← t + Δt;  fire:
        n     → n − 1
        E_int → E_int − D₀(n_pre)                  (K1, pre-shed rung)
        m     → m − m_He;  v → m/(m−m_He)·v        (cold shed, as-built:
                                                    the delivered kick)
        book per-event ledger terms (§2.5)
```

Termination is unconditional: n strictly decreases on every fire, so the loop
runs ≤ n_h ≈ 21 iterations regardless of t_detect (µs, ms — free).

**No re-entry pathology:** between fires E_int is constant, so the gate status
cannot change spontaneously; across a fire the self-bound margin
G = E_int − Σ(n) is shed-invariant (E_int and Σ drop by the same D₀). An
in-band ion therefore stays in-band until the energetic floor
E_int < D₀(n); a suppressed or frozen ion is suppressed/frozen **forever**.

### 2.3 Permanent states (k = 0 is final under P1–P3)

| Reason | Condition | Physical reading |
|---|---|---|
| `frozen` | n = 0, or E_int < D₀(n), or the RRK-bracket **underflow band** (in-band with E_int within ~10⁻⁶ relative of D₀(n): `(1−D₀/E_int)^(s−1)` underflows to k = 0; expected waiting time exceeds any flight — review fix, previously a crash) | energetic floor — converged terminal |
| `suppressed` | E_int > Σ(n), n ≥ 2 | ejected self-unbound complex; rides to the detector at its handover n (the gated τ ≥ 16.5 class, n = 21). Permanent **only because the guard verified zero cooling** — suppression is cooling-reversible (§2.1 NB) |
| `time_exhausted` | t_detect reached with k > 0 | cascade live at the detector — the *physical* in-flight snapshot (this is the read Wave 5 could not compute) |

The `suppressed` semantics are the **delivered gate semantics faithfully
extended** — whether an ejected net-self-unbound complex should instead
fragment (a candidate bare-I⁺ channel, OQ-B) is explicitly *not* resolved by
this stage; cross-linked, not built.

### 2.4 Translation is out of scope (and why that is exact)

Cold shed kicks v by the delivered `m/(m−m_He)` reset at each fire
(**as-built correction** — the approved draft said "v unchanged"), so
post-ejection v is **piecewise constant**: changed only at fires, never
between them, and never read by the (E_int, n) jump chain. Positions never
feed back (P1–P3 remove every position-dependent channel, verified at
handover). The stage therefore propagates the **mass subsystem plus the
v/ledger passengers** — no positions, no forces, no BAOAB. E_kin at
detection is recomputed from the terminal (m, v) by the existing path; this
is the same "decoupling fact" already documented and audited in
`relaxation_stage.py`.

### 2.5 Dimensional analysis (strict, per the working method)

| Quantity | Units | Balance check |
|---|---|---|
| ν, k | ps⁻¹ | RRK bracket (1 − D₀/E_int)^(s−1) dimensionless: eV/eV |
| u | 1 | uniform draw |
| Δt = −ln(u)/k | ps | [1]/[ps⁻¹] = ps |
| t_h, t_detect, event times | ps | one absolute axis (continues the ion-stage/E2 axis) |
| E_int, D₀(n), Σ(n) | eV | K1 drain: eV − eV |
| m | amu (stored kg via existing conversion) | fire: m − m_He |
| v | Å/ps | cold-shed kick m/(m−m_He) per fire (as-built); constant between fires |
| ΔE_mass_transfer | amu·Å²/ps² | `mass_jump` mechanical convention; converted at booking by the existing path |
| ε_drain criterion | 1 | [1]·[ps]/[ps] |

There is **no timestep**, hence no `ν·dt ≤ 0.1` guard: exponential waiting
times have zero discretization bias. The formal equivalence claim: the E2
per-step Bernoulli chain `P = 1 − exp(−k·dt)` at constant k converges in
distribution to this exact process as dt → 0 — i.e. the delivered fixed-dt
path is the biased approximation of *this* stage, not the reverse. Test 6
(§5) verifies the claim statistically.

### 2.6 Composition (rule-1 reuse, no duplicate physics)

The fire path composes the **same delivered primitives** as
`evaporation_step`: `rrk_rate` (all gating), `mass_jump.cold_shed`
(velocity/mass/KE-defect reset), `internal_energy_budget.dE_int_shed_eV` (K1,
pre-shed n), `solvation_cooling.e_bind_pair_eV` (the E_pot fold). Only the
*when* changes: the per-dt Bernoulli draw is replaced by the exponential
waiting-time draw. No formula is reimplemented.

## 3. Interface contracts (no code until trigger)

### 3.1 Module and entry point

`i2_helium_md/simulation/detection_stage.py` ::
`run_detection_stage(seed_ckpt, cfg, *, rng=None, save_path=None)
-> DetectionResult`

- Seeds from the seed checkpoint's **final** stored column with E2's
  seed-coherence guards (biphasic-scenario check; `mass_history_kg[:, -1] ==
  mass_final_kg` stride guard). `seed_ckpt` is `relaxation.npz` when E2 ran,
  or `ion.npz` on the skip path (§1 item 5) — both are v7 `IonCheckpoint`s,
  one loading contract; the pipeline decides which file to hand over
  (`relaxation_stage_enabled` selects it).
- Enforces §2.1 P3 per ion; fail-loud with the violator list.
- **RNG:** fresh stream via `SeedSequence((cfg.seed, DETECTION_STREAM_KEY))`
  (a fixed module-level key, the `RELAXATION_STREAM_KEY` pattern). No
  existing stream is touched; draws are per-ion sequential in ion-major
  order. The draw count is data-dependent — acceptable because the stream is
  stage-private; there is no cross-stage draw-order contract to preserve.

### 3.2 `DetectionResult` and artifact

Per-ion arrays: `n_detected`, `E_int_detected_eV`, `mass_detected_kg`,
`state_reason` (`frozen` / `suppressed` / `time_exhausted`), plus **ragged
per-ion event records**: absolute shed times [ps], pre-shed rungs, per-event
ledger terms (K1 drain, e_bind fold, mass-transfer defect).

Artifact: **`detection.npz`** in the run directory — a new small file
(ragged lists as flat arrays + per-ion offsets; `allow_pickle=False`; shape
validation on load). It is **not** an `IonCheckpoint`; **schema v7 is
untouched** (no forbidden-list contact).

### 3.3 Free sensitivity re-read

Because every shed time is stored, n(t) for **any** t ∈ [t_h, t_detect] is
reconstructable from the artifact alone — but only *up to* the generated
t_detect (events beyond it were never sampled). The report therefore emits
the detection-time sensitivity band (e.g. n̄ at t_detect/4, t_detect/2,
t_detect) as a pure post-processing read, no re-run; if an upper band edge
beyond the Sourced value is wanted, the run must be *generated* at that
largest time and the Sourced value read as an intermediate point. This
subsumes the rejected "Bounded pin + band" option.

### 3.4 Config surface (additive; each field fail-loud at load)

- `detection_stage_enabled: bool = False` — opt-in, the
  `relaxation_stage_enabled` pattern. Not a physics enum: the stage is an
  exact solver of the existing mechanism, not a model choice, so the
  working-method "every model choice behind an enum" rule is satisfied by
  the enabled flag alone (documented rationale, not an omission).
- `detection_time_ps: Optional[float] = None` — **Sourced; no default.**
  `check_detection_config`: enabled ⇒ set, > 0, > the seed window end on the
  absolute axis; biphasic scenario. The relaxation stage is **not** required
  (§1 item 5 skip path): with `relaxation_stage_enabled=False` the stage
  seeds from `ion.npz` and the §2.1 P1–P3 guard is the sole defense.
  **Sourced value (resolved 2026-07-07): 8.53·10⁶ ps (8.53 µs)** — the TOF
  flight time from the experimental-setup publication (user-confirmed;
  CALIBRATION_MAP row 24 carries the provenance). The config field still
  requires the explicit value (no baked-in default — the constant is
  calibration data, not code).
- `ε_drain`: documented module constant (§2.1), **not** a config field.

### 3.5 E1 / report integration

`compute_terminal_shell_distribution` gains the `source_tag="detected"` read
through the same duck-typing hook as `"relaxed"` (a `terminal_n`-bearing
input). Probe/campaign reports gain `n_detect_*` columns and the state-reason
fractions. **The detected read becomes the Tier-2 arbitration observable**;
the relaxed read is demoted to a convergence diagnostic alongside
`frac_frozen`. (Report-side changes ride the F3 idiom; the campaign-scorer
wiring belongs to the B.5/F2 re-scope, not this stage.)

## 4. Scope boundaries — what this stage is NOT

- **Not new physics.** Rates, gates, ladder, K1, cold shed: delivered code,
  reused (§2.6). The stage changes the integration scheme in a regime where
  that change is exact.
- **Model-scope caveat (new OQ-E, documented not built):** over the µs
  flight the only active channel is RRK evaporation — no radiative cooling,
  no electronic relaxation, no residual-gas collisions. A domain-expert
  question, recorded for the production discussion.
- **OQ-B untouched:** `suppressed` ions ride to the detector intact at their
  handover n (= 21 for the gated τ ≥ 16.5 class). Whether such ejected
  self-unbound complexes should instead fragment to bare I⁺ is the standing
  bare-peak-channel discussion, cross-linked, not resolved here.
- **Nothing here discharges the F5 gate** — the stage makes the gated
  terminal *computable*; adjudication stays with the 0.80 eV N=500 campaign
  (reporting-gate stance unchanged).

## 5. Validation plan (the seven-step order)

1. **Formula:** seed-fixed Δt = −ln(u)/k against hand values; `rrk_rate`
   reuse means no new rate tests (they exist).
2. **Shape/units:** `detection.npz` round-trip; ragged offsets; kg/amu
   conversions at the boundary.
3. **One-step deterministic:** k = 0 lanes (`frozen`, `suppressed`, n ≤ 0)
   → zero events, zero draws consumed, n preserved to t_detect, reason
   recorded. Skip-path seeding: an `ion.npz`-seeded call produces the same
   result as a `relaxation.npz`-seeded call whose E2 window did nothing
   (all-frozen / already-P3-clean seed); guard failure on a
   non-P3-compliant ion seed raises with the per-ion list.
4. **Multi-event:** constructed ladder with a known rung sequence → occupancy
   at t vs the **hypoexponential closed form** (the analytic oracle);
   sample-size-justified tolerance.
5. **Energy bookkeeping:** 5-term closure across the E2 → detection boundary;
   residual in the existing ~2·10⁻⁵ eV class.
6. **Statistical:** same handover state, short horizon — event-driven vs a
   fixed-dt E2 extension; distribution/moment agreement (the §2.5
   equivalence claim).
7. **Driver smoke:** tiny-N end-to-end run dir producing `detection.npz`;
   report reads it. No figures, no production-sized artifacts in pytest.

**Production wiring oracle (manual, non-pytest):** continue the Wave-5 gated
s_eff = 30 run to t_detect — the first physically meaningful gated terminal
read; immediately answers whether that point reaches the n ≈ 2 floor by the
detector or is caught mid-cascade (`time_exhausted`).

**Review scope (§1 item 5):** one connected code review covers the DS build
*and* the E2 contract re-frame (docstring mission update, relaxed-read
demotion, skip-path config guard) — the review question is the coherence of
the stage boundary and the demoted E2 contract.

## 6. Recorded consequences (not built now)

- The **B.5 / F2 campaign re-scope is unblocked**: the gated terminal
  observable is well-defined at a Sourced detection time; the campaign
  scorer reads `detected`, reports `relaxed` + `frac_frozen` as convergence
  diagnostics.
- **Campaign economy — now adjudicated (§1 item 5):** with the continuation
  owning the tail, gated E2 caps are shortened by config (future-facing
  USER SETTINGS; the §2.1 guard fails loud if shortened past
  ejection/decoupling), and gated runs may skip E2 entirely
  (`relaxation_stage_enabled=False`, `ion.npz` seeding). The ungated arm
  keeps its full E2 (it is the terminal solver there); on the skip path the
  `frac_frozen` convergence diagnostic does not exist — the DS
  `state_reason` fractions are the replacement read.
- The staircase-probe findings' boundary items on cap-truncated gated reads
  (§6 of the findings doc) become historical once the detected read exists.

## 7. Open items before `[PROCEED TO IMPLEMENTATION]`

1. ~~`detection_time_ps` value + provenance~~ **Resolved 2026-07-07:**
   8.53 µs, setup publication (user-confirmed) → CALIBRATION_MAP row 24.
2. ~~Slice placement in the Phase-F program~~ **Resolved 2026-07-07:**
   standalone pre-F5 **Slice DS** (the staircase-probe insertion pattern),
   scope = DS build + probe-report `n_detect_*` extension + the §1-item-5
   E2 re-frame/skip path, one connected code review. The B.5/F2 plan update
   stays parked with the campaign (post-Wave-6 decision 3) — explicitly
   *not* bundled into this slice.

None remain — the design awaits `[PROCEED TO IMPLEMENTATION]`.

**Cross-links:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` (I11/I14, §6
boundaries); `TIER2_STAIRCASE_PROBE_PLAN.md` Addendum B;
`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (F2/F5); `CALIBRATION_MAP.md`
(pending t_detect row); `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`
(§R5, A11); `drag_migration_log_tier2.md` (decision/delivery record).
