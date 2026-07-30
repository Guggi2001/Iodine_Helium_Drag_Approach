# Tier 2 — Free-Form Linear Drag Twin Sweep (atlas counterfactual arm)

**Status: Step 0 + BOTH ARMS EXECUTED 2026-07-30 (findings
in-chapter): §2.1 twin KE₁ LICENSED ρ 0.998; §3.1/§5.1 — arm-1 basin
35 cells incl. 6 sub-plateau (kill-3 universality breaks at twin
level), R4 fires (twin KE₁ 0.892 vs anchor 0.624); arm 2 — the
landing rejects the ram term, c = 0 optimal, pure linear selected
in-family; §7.1 oracles all green. NEXT GATE: the §6 CRN MD ring vs
h405 (user trigger) — candidate set = the arm-1 sub-plateau core
a 27.5–35.**

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

### 2.1 Findings — Step 0 EXECUTED 2026-07-30 (zero MD): twin KE₁ ranking **LICENSED**, ρ = 0.998

Stage `ke1auth` in `scripts/tier2_h2b_forward_model.py` (committed);
outputs `atlas_ke1_authority.csv` + `_summary.csv` (committed); helper
tests `tests/test_tier2atlas_ke1auth.py` (11 passed).

- **Oracle:** the committed `h2b_g3_corrected_row.csv` re-derived
  string-identically at the standing cell before any new number (the
  twin is unchanged where its authority was measured).
- **Frozen replay set:** 26 distinct-pin corrected-ensemble cells kept
  out of the 52-row table; 26 exclusions, every one recorded with its
  reason in the CSV (6 incumbent = standing geometry; 11 geogrid =
  controlled single-R cells; 5 h405 battery members → their pooled row;
  3 duplicate pins deduped per pooled N = 5000 > g4finals N = 1000 >
  g3ring N = 500; 1 f725 = MD n₁ = 0, KE₁ unmeasured). Pin
  spread: v_c 5.0–6.5, all three wells, τ 4.4–6.4, E₀ 0.29–0.415; MD
  KE₁ dynamic range 0.279–0.956 eV.
- **Result: Spearman ρ(KE₁) = +0.9979 over 26 cells (0 twin-NaN);
  ρ(KE₂) = +0.9966. Verdict vs the pre-registered bands (0.8/0.5):
  LICENSED** — twin KE₁ may rank cells and report a best-KE₁ basin
  (§5 R4 active). The n = 1 bin is the twin's *best* KE observable
  (deepKE stays +0.33 unlicensed): short chords, mechanism-light —
  the prior's favorable branch confirmed.
- **Level bias (levels stay unlicensed, standing rule):** twin KE₁
  uniformly cold — ratio twin/MD 0.897–0.984 (median 0.975), additive
  −0.013…−0.037 eV (median −0.017); 3/25 adjacent rank inversions,
  all among near-tie neighbours. Any twin KE₁ number quoted from the
  sweep carries the ~−2 % level stamp.
- **Carried caveat (unchanged):** licensure measured on the
  capped-cubic family; transfer to `lin`/`linq` is an assumption the
  conditional MD ring (§6) backstops.
- D0 §14.4 authority box back-filled (KE₁ row beside W₁/midHot/deepKE).

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

### 3.1 Findings — arm 1 EXECUTED 2026-07-30 (zero MD): 35/12 312 cells gate; the basin spans Φ 0.64–1.39

Stage `linscan` (57 chord integrations ≈ 2.3 h, npz-cached; committed
outputs `atlas_linsweep.csv` / `_chords.csv` / `_gated_ke.csv` /
`_summary.csv`). Basin anatomy:

- **35 gated cells**: τ = 4.8 carries 28/35 (τ = 6.4 the rest; **zero
  τ-flags** — no landing leans on an off-calibration cooling clock);
  all three wells gate (eb0482 19 / eb1168 9 / eb154 7 — shallow
  *favored*, not required: the R3 shallow-end-only prediction was too
  strong); a spans 27.5–60 with the coherent core 27.5–42.5 at
  E₀ 0.35–0.39.
- **Φ anatomy (R2):** 6 sub-plateau (Φ ≤ 0.7: the a = 27.5 and a = 30
  columns at τ 4.8, all three wells), 7 intermediate, 17
  plateau-convergent, 5 above-window; **0 sub-bare** — the entire
  gated set stays above the bare-ram floor, so no cell needs the
  unphysical-flag discount.
- **KE₁ along the basin (twin-ranked, licensed):** monotone in Φ —
  a = 27.5/Φ 0.64 → twin KE₁ 0.892; a = 30/0.69 → 0.817; a = 35/0.81 →
  0.703; a = 40/0.93 → 0.575 (eb0482, τ 4.8 column). Softer fast-class
  force buys n = 1 heat *inside* the histogram gate.
- **Twin-level flags (pre-registered blind spots, MD-undecidable
  here):** deepKE at the sub-plateau cells runs **2.0–3.4× HOT**
  (unlicensed observable — recorded, carries no authority, and is the
  first place the MD ring can kill); above-1.15 = 0.0 everywhere (the
  twin has no source-KER width by construction — the width axis stays
  (B)-side and MD-only); twin trap ≤ 0.011 but trap is a floor and
  the lin form's low-v over-drag (×8–18 at v = 1) makes MD trap the
  single most exposed gate.
- W₁ (reported, bias-loaded): best gated 0.586–0.673 — same range the
  capped twin basin reported at this geometry.

**Arm 2 (`linqscan`, a\* = 35.0, EXECUTED 2026-07-30; 24 chord
families ≈ 1 h, 5 184 cells; committed `atlas_linqsweep*.csv`; §7.4
seam oracle PASSED — linq c = 0 bit-identical to the lin family):
the landing does NOT want the ram term.** 13 cells gate and the c = 0
seam column carries every best cell (W₁ 0.586/0.643/0.677 = the arm-1
a = 35 cells reproduced); everything degrades monotonically with c —
W₁ 0.586 (c = 0) → 0.608 (c = 1) → 0.675 (c = 2) → 0.959+ (c ≥ 3);
twin KE₁ 0.703 → 0.531 → 0.345 → ≤ 0.252; trap explodes (0.18/0.29 at
c = 2/3 bundle-well up to 0.93–0.97 at c = 12.79, where nothing gates
at any well — the shared-lq landmark is dead in this family); the
only far-c gated cells (c = 6.13/9) need τ-flagged off-calibration
cooling (τ 9.6/12.8) and sit at KE₁ 0.11–0.15. c > 0 buys nothing on
any observable at fixed a\*: within the free-form family the
experimental arbitration prefers **pure linear**, and softer-is-hotter
(arm 1) extends to softer-is-*better-landing* against the quadratic
admixture.

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

### 5.1 Findings — readings R1–R5 (arm 1, 2026-07-30)

- **R1 — basin EXISTS** (35 cells). The constant-γ family is *not*
  dead by experiment at twin level; "linear is ruled out" remains
  true only in-band (Tier-0), not as a system statement.
- **R2 — sub-plateau class POPULATED: kill-3's universality claim
  BREAKS at twin level.** Six gated cells carry Φ ≤ 0.7 (down to
  0.64), none sub-bare. The histogram gate does not pin the
  fast-class force inside the linear family — the plateau-pinning
  record (958.5/990.6 convergence, §3.5h) does not extend to this
  family at the corrected geometry. Twin-level statement only; the
  §6 ring is the arbiter.
- **R3 — shallow-well preference CONFIRMED as a preference, REFUTED
  as a requirement** (19/9/7 across the wells; the lq-precedent
  co-compensation shows up as ordering, not as a gate).
- **R4 — FIRES: best gated twin KE₁ 0.892 (a = 27.5, Φ 0.64) vs the
  h405 twin anchor 0.624**, monotone softer-tail → hotter-n = 1 along
  the basin. Under the ~2 % cold level stamp this forecasts an MD
  KE₁ ≈ 0.91 *if* the capped→lin transfer holds — against h405's MD
  0.637, the (A)-ceiling 0.708, and the 1.00 reference peak. This is
  the strongest drag-side KE₁ signal the program has produced;
  everything rides on the transfer caveat and the MD ring.
- **R5 — arm-2 trigger MET; a\* = 35.0 frozen** (the gated-basin W₁
  optimum, 0.586 at eb0482/τ 4.8/E₀ 0.37, per the §3.5c ranking
  convention). Noted for the record: the KE₁ optimum sits at
  a = 27.5 — the a\* freeze follows the pre-registered rule, not the
  KE₁-favorable choice.
- **R5 c-read (arm 2 EXECUTED): the landing wants c = 0.** At frozen
  a\* = 35.0 every observable degrades monotonically with the ram
  term (W₁, KE₁, trap; §3.1 arm-2 block); c = 12.79 (the shared-lq
  landmark) gates nowhere; the far-c stragglers need τ-flags. Within
  the free-form family the data selects **pure linear** — the c-axis
  closes at twin level, and the MD-ring candidate set stays the
  arm-1 basin (sub-plateau core a 27.5–35).
- **§6 proceed criterion: MET** (contiguous gated sub-plateau basin,
  not sub-bare-only, R4 positive). The CRN MD ring vs h405 is now the
  decision point — design to be frozen at ring time, **behind its own
  trigger** (user gate).

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

### 7.1 Findings — oracle outcomes (2026-07-30)

1. Landmark oracle: PASSED before Step 0, arm 1, and arm 2 (row + KE
   files string-identical each time).
2. Non-regression: **measured, not assumed** — after the form-switch
   edit, `ke1auth` re-ran end-to-end: landmark row bit-exact AND both
   committed Step-0 CSVs reproduced string-identically (the capped
   path is byte-identical with the `lin`/`linq` branch present).
   Plus 22/22 focused tests (`test_tier2atlas_ke1auth.py` +
   `test_tier2_linsweep.py` — validation, seam bit-identity on a toy
   ensemble, dissipation direction, Φ/sub-bare/R2-class arithmetic).
3. Step-0 drift oracle: exercised on the re-run (green, see 2).
4. Seam oracle: run at arm-2 launch (result recorded in §3.1's arm-2
   block).

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
