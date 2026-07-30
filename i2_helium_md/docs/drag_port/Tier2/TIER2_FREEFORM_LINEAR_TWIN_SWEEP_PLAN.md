# Tier 2 — Free-Form Linear Drag Twin Sweep (atlas counterfactual arm)

**Status: DESIGNED 2026-07-30 (discussion-adjudicated in-session; awaiting
user review). Zero code exists; all implementation sits behind
`[PROCEED TO IMPLEMENTATION]`.**

Companion documents: `TIER2_SENSITIVITY_ATLAS_PLAN.md` (§3.5c is the
gate/grid precedent reused here), `TIER2_SENSITIVITY_ATLAS_FINDINGS.md`
(§6.6/§6.7 — the lq counterfactual precedent and the form-blindness
record), `TIER2_PARAMETER_INFLUENCE.md` §1/§2 (the form-entry and
capped-tail ledgers this arm feeds), `drag_migration_log_tier2.md`
(2026-07-30 entries: velocity anatomy, tail-force invariance, the
mechanism correction).

---

## 0. Stance, motivation, and what this arm is *not*

**Motivating premise (user, 2026-07-30, recorded verbatim in spirit):**
at the corrected geometry the production cap sits at v_c = 5.5, at the
top of the TDDFT band (4.95 Å/ps) — so for the fast class (n = 1/n = 2
exit velocities ≈ 9.7–12.3 Å/ps) the operative drag curve is already
**fully experiment-fitted**. Sweeping a free-form linear coefficient is
therefore "not worse than fitting v_c": it trades the in-band Tier-0
anchor for one fewer structural assumption (no cap, no kink), and lets
the experimental observables arbitrate a maximally rigid one-parameter
family.

**Stance.** This is an **atlas counterfactual arm** in the §6.6 lq
tradition, not a production-form change and not a Tier-0 revision:

- Tier-0 remains the recorded in-band form authority
  (`shared_pure_cubic`, held-out transfer test; the lq FAIL 0.699 and
  the linear-family in-band evidence — a → 0 whenever a higher power
  competes, n̂ = 2.927 — stand unedited).
- The histogram landing can never be cited as evidence *for* any form
  (form-blindness is measured, W₁-only — findings §6.6/§6.7). What
  *can* discriminate forms at Tier-2 level are the form-sensitive
  observables (KE₁/KE₂, fate split, midHot/n̄ at matched well). If the
  linear system were to beat cubic there **at MD level**, that opens an
  adoption discussion with real standing; nothing short of that does.
- Any landing found here is an *effective-theory* landing. The arm's
  primary product is understanding: it is the sharpest available test
  of the plateau-pinning claim ("no landed system carries sub-plateau
  tail force" — currently n = 1 form families at the corrected
  geometry), because a constant-γ form has **no second lever**: one
  parameter sets the force at every velocity, so the fast class cannot
  be softened without committing the slow class to the consequences.

**What this arm is not:** not a p_tail re-open (§3.5h stands), not a
Tier-0 re-fit (the pure-linear Method-B rows are a separate, smaller
proposal — see §9), not a source-side/(C) revisit.

## 1. Forms and dimensional analysis

Friction convention unchanged (CLAUDE.md): γ(v) is a force coefficient
[amu/ps], F_drag = γ·v [amu·Å/ps²], drag module mass-agnostic, density
gate ρ̂(r) (shared erf steepness) multiplies γ exactly as in production.

**Arm 1 — pure linear (`lin`):**

    γ(v) = ρ̂(r) · a            a [amu/ps]
    F(v) = ρ̂ · a · v           [amu/ps]·[Å/ps] = [amu·Å/ps²]  ✓

One drag parameter. Constant friction rate a/m: every ion decelerates
on the same timescale regardless of speed — the maximal-rigidity limit.

**Arm 2 — linear + quadratic, linear part frozen (`linq`):**

    γ(v) = ρ̂(r) · (a* + c·v)   a* [amu/ps] frozen (§5 R5), c [amu/Å]
    F(v) = ρ̂ · (a*·v + c·v²)   [amu/Å]·[Å²/ps²] = [amu·Å/ps²]  ✓

c ≥ 0 adds ballistic ram growth at the top; ram reading
σ_eff = c/(ρ_He·m_He) with ρ_He = 0.0218 Å⁻³, m_He = 4.0026 amu
(c = 6.13 ⇒ σ_eff ≈ 70 Å², r_eff ≈ 4.7 Å — the user's 18 Å per-case lq
fit, D4 form table). If the landing *wants* c > 0, that is the data
asking for the ram term — a measured statement, recorded as such.

**Tail-force coordinate (the kill-3 metric).** For every cell:

    Φ ≡ F(9.7 Å/ps) / 418.5

418.5 amu·Å/ps² = b·v_c³ at the standing corrected-geometry point
(b = 2.5153508541760052, v_c = 5.5) — the h405 plateau, the absolute
anchor. Arm 1: Φ = 9.7·a/418.5 (a ≈ 43.1 ⇒ Φ = 1; a ≈ 30.2 ⇒ Φ = 0.7).

## 2. Step 0 — twin KE₁ ranking-authority measurement (zero MD, runs FIRST)

The observable this arm most cares about (KE₁) has **no measured twin
licensure** (G4 Block 0 licensed W₁ ρ +0.82 and midHot ρ +1.00,
rejected deepKE ρ +0.33; KE₁ was never in the box). It is measurable at
zero MD:

- **Instrument:** the twin's existing `ke_curve` (per-bin mean
  KE = ½·m(n_det)·v_inf², mass at detected n, Coulomb in v_inf,
  co-moving-consistent vacuum cascade). KE₁ = bin-1 mean; KE₂ = bin-2;
  above-1.15 share computable per-fragment. Surfacing these as scan
  columns is the only code.
- **Reference set:** the committed `atlas_ke_lown_scan.csv` (52 rows,
  MD KE1_mean/KE2_mean with full (v_c, well, τ, E₀) pins). Replay set =
  the corrected-geometry cells with distinct pins; standing-geometry
  reference rows excluded; pooled rows preferred over members where
  both exist. **The exact row list is frozen in the scorer before any
  ρ is computed.**
- **Measurement:** twin run at each row's pins (capped_cubic, p_tail
  −1 — the family the MD rows belong to); Spearman ρ(twin KE₁, MD KE₁)
  across the replay set; ρ(KE₂) reported as secondary.
- **Pre-registered outcome bands:**
  - ρ ≥ 0.8 → **KE₁-ranking LICENSED** (the W₁-precedent level): the
    sweep may *rank* cells by twin KE₁ and report a best-KE₁ basin.
  - 0.5 ≤ ρ < 0.8 → **directional only**: twin KE₁ reported per cell,
    never used to select or rank; KE₁ verdict is MD-only.
  - ρ < 0.5 → **UNLICENSED**: the sweep gates histogram-side only
    (deepKE precedent); twin KE₁ columns still recorded for the
    post-hoc authority ledger, with the unlicensed stamp.
- **Transfer caveat (pre-registered):** licensure is measured on
  capped-cubic cells; its transfer to the linear family is an
  assumption. Backstop: the conditional MD ring (§6) delivers the KE
  verdict regardless of what the twin claimed.
- KE₁ reference convention throughout: the §3.5g amendment — scored
  against the experimental n = 1 **peak ≈ 1.00 eV**; older anchors
  (median 1.128, mean 1.302) reported, not chased. KE₂ vs mean 0.706.

## 3. Sweep design (zero MD)

Structure and machinery are §3.5c verbatim (nested Route A/B): trap and
the KE curve depend only on the chord integration — here the
**(a, E_bind)** surface (arm 2: **(c, E_bind)** at frozen a*); τ (exact
K rescale) and E₀ (fate map) are free post-processing. Twin caveat
carried forward: in MD, τ feeds back on trap through the mass
mechanism — the clean decomposition is itself a twin approximation
(authority box, D0 §14.4).

**Fixed:** the committed corrected master draw (seed 20260727,
m = 20000), rq4graded ladder, p = 1, detection = the non-trapped chord
read, E_bind grid = the Tier-0 extracted spread (the §6.5/§6.7
documented joint-pairing exception, stamped per table).

**Arm-1 grids:**

- Chord surface: **a ∈ [15, 60] amu/ps, step 2.5** (19 values; the
  plateau-match hypothesis a ≈ 43 and the sub-plateau hypothesis
  a ≈ 30 both sit mid-grid) × **E_bind ∈ {0.048, 0.1168, 0.154} eV**
  → 57 chord integrations (≈ 1.5–2× the g3scan chord cost; npz-cached).
- Free surface: **τ ∈ {2.4, 3.2, 4.8, 6.4, 9.6, 12.8} ps** (any landing
  requiring τ ≫ the sourced 6.55 is flagged against its calibration
  class, never adopted silently); **E₀ ∈ [0.17, 0.52] step 0.01**
  (§3.5c values unchanged) → 57 × 216 ≈ 1.2·10⁴ scored cells.

**Arm-2 grids (conditional, §5 R5):** a* frozen;
**c ∈ {0, 1, 2, 3, 4.5, 6.13, 9, 12.79} amu/Å** (0 = arm-1 anchor
column; 6.13 = the user's 18 Å per-case fit; 12.79 = the shared-lq
landmark) × the same E_bind/τ/E₀ grids. 24 chord integrations.

**Per-cell outputs:** the committed observable vector + det_yield + the
detected-subset R quantiles (§3.5c convention) + **KE₁_mean, KE₂_mean,
above-1.15, Φ** (this arm's additions).

## 4. Gates — §3.5c frozen bands, reused verbatim

**Hard gate (unchanged):** n₁_solv ∈ [0.19, 0.30] AND n̄ ∈ [4.4, 7.1]
(target 4.07 + the residence-conditional twin bias bracket [+0.3, +3]).

**Reported, not gating (unchanged):** W₁ (twin-bias-loaded); supp soft
(RQ3 mixture); trap reported as a **floor** with composition beside it.

**This arm's reported additions (never gating):** KE₁/KE₂/above-1.15
(status per Step-0 licensure), Φ per cell, and **deepKE with the
standing blind-spot disclosure** (ρ +0.33 — the twin cannot see the
deep bins; a twin-gated linear basin can still crater deep-bin KE in
MD; pre-registered as the known MD failure mode so a later MD kill is
not goalpost-moving).

No band is re-frozen, tightened, or loosened for this arm.

## 5. Pre-registered readings and outcome semantics

- **R1 — existence.** Does any (a, E_bind, τ, E₀) cell pass the hard
  gate? **NO cell** ⇒ the constant-γ family is dead **by experiment at
  the corrected geometry** — the direct, TDDFT-independent answer to
  "why is linear ruled out", and a strengthening of the plateau-pinning
  record from outside the capped family. (Zero MD spent.)
- **R2 — where the basin sits (the kill-3 read).** For gated cells,
  read Φ:
  - **Φ ∈ [0.85, 1.15] (plateau-convergent):** a third independent
    form family lands on the same fast-class force — plateau pinning
    upgraded to form-universal at the corrected geometry; KE₁ stays
    cold by construction; the KE₁ lever stays source-side per the (C)
    record.
  - **Φ ≤ 0.7 (sub-plateau):** the histogram can be held with
    materially less fast-class drag — the pinning claim breaks at twin
    level, the drag-side KE₁ door reopens, and the MD ring (§6)
    decides.
  - **Intermediate (0.7 < Φ < 0.85):** reported as a measured boundary;
    MD ring at user discretion.
- **R3 — E_bind preference.** Pre-registered prediction: the basin (if
  any) sits at the **shallow end (0.048)** — the lq precedent (a
  low-power form buys back its low-v over-drag with a shallow well,
  §6.7 item-2). Confirmation/refutation is a mechanism statement about
  the (form, well) co-compensation, recorded either way.
- **R4 — KE₁ ranked read** (only if Step 0 licensed): does any gated
  cell's twin KE₁ beat the h405-family twin KE₁ at matched gates? If
  licensure came back directional-only or unlicensed, R4 is skipped and
  its question moves wholesale to the MD ring.
- **R5 — arm-2 trigger.** Arm 2 runs only if arm 1 produces a gated
  basin (a* = the arm-1 basin optimum by the §3.5c ranking convention)
  OR the user explicitly asks for the c-continuum anyway (a* then = the
  best-W₁ non-gated a). The c-read: does the landing want c > 0 (the
  ram term) and does c move Φ/KE₁ at fixed histogram?

## 6. Conditional MD confirmation ring (behind its own trigger)

Runs only on a §5-qualifying basin (R2 sub-plateau, or R4 positive, or
user call on an intermediate). Design frozen at that point, not now;
its committed conventions are fixed here:

- ≤ 10 cells × N = 1000, fresh seed, **CRN-paired vs the h405 pins**
  (the §3.5k/§6.7 pairing convention).
- KE/fate acceptance bands **frozen pre-launch** in the ring design
  (the G3 Step 3 precedent), covering: KE₁_mean/KE₂_mean/above-1.15,
  W₁/n₁_solv/n̄, supp/trap, **deepKE** (the twin-blind observable gets
  its verdict here).
- Scorer extends `tier2atlas_g3ring_table.py` conventions; **the full
  observable vector (all KE columns included) is committed to CSV**
  (the §6.7 regeneration-cost lesson, memory rule).
- The ring is the backstop for every twin-authority assumption in this
  plan (Step-0 transfer caveat, deepKE blind spot, n̄ bias bracket).

## 7. Oracles and non-regression (run before any new number)

1. **Landmark oracle:** `h2b_g3_corrected_row.csv` re-derived bit-exact
   at the standing cell (g3scan Block-0 convention) + the S6 machinery
   oracle — proves the twin is unchanged where it was measured.
2. **Non-regression:** with the linear branch present but unselected,
   the g3scan cached chord surface reproduces bit-exact (the capped
   path must be byte-identical; `lin`/`linq` live behind their own form
   switch, default off).
3. **Step-0 internal oracle:** the twin replay at the standing h405
   battery pins reproduces the twin numbers recorded in the Step-0
   table on re-run (drift guard for the licensure measurement itself).
4. **Arm-2 seam oracle:** `linq` at c = 0 reproduces the arm-1 `lin`
   row at the same (a*, E_bind, τ, E₀) bit-exact.

## 8. Records and bookkeeping

- Results → a new **"Free-form linear twin sweep"** section in
  `TIER2_SENSITIVITY_ATLAS_FINDINGS.md`; influence deltas → D0
  (`TIER2_PARAMETER_INFLUENCE.md`) §1 form-entry bullet **first** (the
  D0-first memory rule), §2 tail ledger if Φ moves anything;
  chronology → `drag_migration_log_tier2.md`.
- The Step-0 licensure result also back-fills the D0 §14.4 authority
  box (KE₁ joins W₁/midHot/deepKE with a measured ρ) — it is reusable
  program capital independent of this arm's outcome.
- Nothing in this arm moves `finc1v725`/h405, the Tier-0 record, or
  any production preset. Adoption semantics: §0.
- CSVs land under `data/runs/h2b_forward_model/` and are committed
  (`atlas_ke1_authority.csv`, `atlas_linsweep.csv`, `atlas_linqsweep.csv`).

## 9. Code scope (all behind `[PROCEED TO IMPLEMENTATION]`)

Twin-side + scorer-side only; **zero production-path changes**:

1. `scripts/tier2_h2b_forward_model.py`: `lin`/`linq` branch in the
   chord accel (one `gamma` expression each, behind a form switch
   defaulting to the existing capped path), KE₁/KE₂/above-1.15/Φ
   columns, three new stages (`ke1auth`, `linscan`, `linqscan`).
2. A Step-0 scorer (`tier2atlas_ke1_authority.py` or a stage in 1):
   the frozen replay-row list, the ρ computation, the outcome stamp.
3. Focused tests: the §7 oracles as pytest where cheap (form-switch
   default-off byte-identity; linq c = 0 seam).

**Separate smaller proposal, not part of this arm (own trigger):** the
pure-linear Method-B rows in the D4 Step 1 form table (18 Å-only /
9 Å-only / shared / held-out) — closes "linear was never fit to the
traces" as a measured Tier-0-side row set; pre-registered expectation:
in-band objective fails on the 9 Å curvature, held-out worse than lq's
0.699.

---

*Design discussion record (2026-07-30, this session): the four-kill
audit (in-band rejection / slow-ion arithmetic / plateau pinning /
Tier-0 authority), the concession that the fast class is already fully
experiment-fitted at v_c 5.5 ≈ band top, the band-matched-linear ≈
plateau arithmetic (a ≈ 43 ⇒ F(9.7) ≈ 418), the §6.6 counterfactual
precedent, and the Step-0 licensure idea (KE₁'s authority is measurable
against the committed 52-row KE table before the sweep leans on it).*
