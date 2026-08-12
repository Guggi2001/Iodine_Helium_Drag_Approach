# Tier 2 — Free-Form Linear Drag Twin Sweep (atlas counterfactual arm)

**Status: Step 0 + BOTH ARMS EXECUTED 2026-07-30 (findings
in-chapter): §2.1 twin KE₁ LICENSED ρ 0.998; §3.1/§5.1 — arm-1 basin
35 cells incl. 6 sub-plateau (kill-3 universality breaks at twin
level), R4 fires (twin KE₁ 0.892 vs anchor 0.624); arm 2 — the
landing rejects the ram term, c = 0 optimal, pure linear selected
in-family; §7.1 oracles all green. Ring design FROZEN 2026-07-31
(§6.1: 8 × N = 500 CRN, kills per §4.1 — trap kill removed,
policy-robustness added); outcome-space edges ADJUDICATED same day
(§6.2). **MD ring EXECUTED 2026-07-31 (§6.3): the ring LANDS — 6/7
cells gate both-policy-clean, K-KE does not fire, SUCCESS at five
cells (best KE₁ 0.903 at a 27.5 vs h405p 0.637, ΔKE₁ +0.266
CRN-paired); capped→lin twin transfer HELD (~1–2 % cold); trap dead
(≤ 0.011, lin frees h405p's trapped ions); costs = mid-band
overheating (midHot 1.4–2.0) + χ²; nothing adopted. NEXT USER GATE:
the adoption discussion + escalation battery (§6.3).**

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
  unphysical-flag discount. **The plateau-convergent class is the
  family's SECOND personality — read out in §6.4: it clones the h405
  landing with no cap at all.**
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

### 4.1 Ring-level kill adjudication (2026-07-31, discussion record): the trap kill is REMOVED

The draft ring design carried "trap ≥ 0.15 at every core cell → family
dead." **User challenge sustained; the number was incumbent-pattern-
matched (f725's 0.437), not experiment-derived, and it does not
survive scrutiny:**

- f725's falsification ran through the histogram gates (n₁ = 0); its
  trap number was the *mechanism* of the failure, not the verdict.
- Trap is not constrained by the scored surface at all: a trapped ion
  exits the detected window entirely (experimentally: an unscored
  heavy ion-doped-droplet channel). The n = 20–21/suppressed lesson
  applies verbatim — fate bookkeeping outside the detected window has
  no experimental anchor. A high-trap system that lands the detected
  gates is making an unscored *prediction*, not failing.

**What replaces it (frozen for the ring):**

1. **Kills are experiment-anchored only:** the §3.5c hard gate
   (n₁_solv/n̄, MD-side) and K-KE (CRN-paired ΔKE₁ vs h405, §6.1).
2. **Policy-robustness condition:** every cell is scored under both
   ends of the retained-policy bracket (the open G4
   `exclude_all_coupled` question); a cell whose gate verdict flips
   with policy is **policy-blocked** (neither passed nor killed —
   escalates to the user). A wholly policy-blocked basin escalates,
   never auto-resolves.
3. **Trap is reported, decomposed** (bound/marginal, g3ring column
   convention), plus the **CRN fate-flow read**: same-seed per-ion
   comparison vs the h405 partner — who traps under lin that escaped
   under capped, from which (R, depth, v₀) stratum. This is the
   landing-by-amputation diagnostic (the §4ee precedent in reverse):
   reported, never gating; it determines what a landing *means*
   (cascade-carried vs selection-carried), not whether it counts.
4. **deepKE stays interpretation-banded, not auto-killing** (same
   logic applied consistently): [0.6, 1.6] reads neutral-to-improved
   (the NB-RQ11-12 mid-band-softening direction); ≥ 2.0 at every core
   cell reads as the twin's hot flag confirmed → **form-questioned,
   user adjudicates**.

Consequence, stated openly: the ring's genuine kill surface is narrow
(histogram gates + ΔKE₁) — which is honest, because those are the
only observables the experiment actually pins. Everything else the
ring produces is calibration (the lin-family authority box) and
mechanism understanding.

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

### 6.1 Ring design FROZEN (2026-07-31, discussion-adjudicated; awaiting `[PROCEED TO IMPLEMENTATION]`)

**N choice (user-adjudicated): N = 500, not 1000** — the G3-Step-3
precedent, motivated by the recalibration strategy (MD may shift the
basin; the ring then measures the lin-family twin↔MD bias vector and
the twin re-scans with it, rather than burning MD hunting). Measured
resolution at N = 500: trap SE ≈ 0.01, KE₁ seed-SE ≈ 0.005 vs an
expected signal +0.1…+0.27, n₁ SE ≈ 0.02 vs gate half-width 0.055,
deepKE direction-only (sufficient for its band read), W₁ ≈ ±0.13
single-seed (disclosed; W₁ is nowhere gating).

**Cells: 8 × N = 500, one fresh shared seed 20260731, CRN across all
cells including the incumbent partner.** Axes deliberately spread so
the ring doubles as the lin-family authority box (∂(twin−MD)/∂knob
measurable even if every cell individually misses):

| # | cell | role |
|---|---|---|
| 1 | a 27.5 / eb0482 / τ 4.8 / E₀ 0.35 | KE₁ optimum, Φ 0.64 |
| 2 | a 27.5 / eb1168 / τ 4.8 / E₀ 0.35 | well axis at the optimum |
| 3 | a 30 / eb0482 / τ 4.8 / E₀ 0.36 | sub-plateau W₁-best |
| 4 | a 30 / eb1168 / τ 4.8 / E₀ 0.36 | well axis, second point |
| 5 | a 32.5 / eb0482 / τ 4.8 / E₀ 0.36 | mid-core |
| 6 | a 35 / eb0482 / τ 4.8 / E₀ 0.37 | W₁ optimum (a\*), Φ 0.81 |
| 7 | a 27.5 / eb0482 / τ 4.8 / E₀ 0.37 | E₀ axis at fixed chord |
| 8 | h405 (capped, v_c 5.5 / τ 4.4 / E₀ 0.405) | CRN incumbent partner |

**Frozen bands and semantics (kills per §4.1 — trap kill removed):**

- Hard gate: §3.5c bands MD-side (n₁_solv [0.19, 0.30], n̄ [4.4, 7.1]).
- **K-KE kill:** CRN-paired ΔKE₁ < +0.05 eV vs h405 at every gated
  cell → co-selection wins in real MD; the drag-side KE₁ question
  closes with the full instrument chain.
- **Success:** any gated cell with KE₁ ≥ 0.75 (beats the measured
  (A)-ceiling 0.708).
- Policy-robustness, trap reporting + CRN fate-flow, deepKE
  interpretation band: §4.1 verbatim.
- Reported: W₁ (±0.13 disclosure), above-1.15 + KE₁_SD (the width
  read — the twin predicts nothing here), supp, full observable
  vector to committed CSV (the regeneration-cost rule).
- **Escalation pre-registered:** survivors → one N = 1000 × 5-seed
  battery at the single best cell (the G4 pattern); no wider ring
  without a new adjudication. (Survivor defined in §6.2 item 2;
  all-miss semantics in §6.2 item 1.)

**Code scope (production path — this is the first non-twin code of
the program arm):** `lin` drag form in `physics/drag.py` behind its
own `SimConfig` enum value (architecture rule: every form behind an
enum; mass-agnostic γ(v) = ρ̂·a, same gate convention), config-load
guard posture as biphasic (`allow_inconsistent_mass_pairing=True`,
§6.6 mid-window defense), generator `gen_tier2atlas_linring.py` (g3ring
pattern, CRN seed plumbing), scorer extending
`tier2atlas_g3ring_table.py` (KE columns + trap decomposition + both
retained-policy scores + the CRN fate-flow read). Cost ≈ 5.5 h MD +
scoring. All behind `[PROCEED TO IMPLEMENTATION]`.

### 6.2 Outcome-space edges ADJUDICATED (2026-07-31, user; pre-registered before launch)

Six semantic gaps the frozen §6.1 design left open, settled by
discussion — none touches the frozen bands, the cell list, or the
§4.1 kill surface:

1. **All-miss is recalibration, not a kill.** If all 7 lin cells miss
   the §3.5c hard gate, the family is **not** declared dead — the
   N = 500 choice was made *for* this branch. One bias-corrected twin
   re-scan is licensed: the ring's measured lin-family twin↔MD bias
   vector (the authority-box product) is applied and the §3 twin scan
   re-run with it. The family dies at this stage only if (a) the
   corrected re-scan finds no gateable basin inside the a ∈ [15, 60]
   grid, or (b) the bias vector is incoherent (sign-flipping across
   cells, so no coherent correction exists). Any second recalibration
   iteration (more MD) is a **new user gate**, never automatic.
2. **Survivor ≠ success.** Escalation eligibility (the N = 1000 ×
   5-seed battery) = hard-gated + not policy-blocked + CRN-paired
   ΔKE₁ ≥ +0.05 eV vs h405. The KE₁ ≥ 0.75 band is the *interpretive*
   success stamp (beats the measured (A)-ceiling 0.708), **not** the
   escalation threshold — a gated cell at e.g. KE₁ 0.70 with
   ΔKE₁ ≥ +0.05 escalates.
3. **Retained-policy bracket ends NAMED** (the §4.1 condition made
   concrete): **Arm A** = `exclude_all_coupled` (headline scoring,
   marginals excluded) vs **Arm B** = marginal-injection at handover
   n with conservative asymptotic KE — the §3.5b bracket construction
   (`tier2atlas_retained_bracket.py`: `marginal_dossier` /
   `arm_b_detection`, the sanctioned-reuse precedent; that file's
   consumer note gets extended when the ring scorer imports it). A
   cell whose hard-gate verdict differs between Arm A and Arm B is
   policy-blocked per §4.1.
4. **deepKE mixed outcome defaults to reported-no-flag.** The
   form-questioned escalation fires only on the frozen condition
   (≥ 2.0 at *every* core cell). Mixed cells and the 1.6–2.0 strip
   are reported with the band read, carry no flag, and trigger
   nothing.
5. **h405 partner is Δ-only.** Lin-cell gate verdicts are absolute
   against the frozen bands. The fresh N = 500 h405 partner serves
   the paired Δ reads (K-KE, fate-flow) exclusively; if seed noise
   puts the partner itself marginally outside a gate band, that is
   disclosed as seed noise and never contaminates any lin verdict.
6. **Twin forecast columns join the ring CSV.** The scorer commits
   each cell's twin-predicted observables beside the MD values, so
   the bias vector ∂(twin−MD)/∂knob is a pure CSV read (full-vector
   memory rule extended to the authority-box product).

**Sequencing (user):** the §9 pure-linear Method-B row set stays
deferred until after the ring — it never blocks the ring and is
revisited only if the ring outcome makes the Tier-0-side answer
relevant.

### 6.3 Findings — MD ring EXECUTED 2026-07-31: the ring LANDS — 6/7 cells gate both-policy-clean, K-KE does not fire, the success band is met at five cells (best KE₁ 0.903 vs h405p 0.637); the capped→lin twin transfer HELD

Instrument (behind `[PROCEED TO IMPLEMENTATION]`, first production-path
code of the arm): `pure_linear` form in `physics/drag.py` +
`SimConfig`/guard surface (a > 0; new `free_form` extraction-method
vocabulary — a free-form parameter may not carry a §6.5.1 binding
stamp, enforced), generator `gen_tier2atlas_linring.py` (8 × N = 500,
CRN seed 20260731, LR-P1 committed-artifact oracle, cfg-diff guard),
scorer `tier2atlas_linring_table.py` (Arm A/Arm B via the §3.5b
bracket construction, KE + width columns, trap decomposition, CRN
fate flow, twin join). Full suite 2992 passed with the form present
(default path inert); focused tests 38. MD ≈ 2 h wall at
concurrency 3, all 8 cells clean. Committed CSV:
`atlas_linring_table.csv`.

**Gate-band clarification (documented at scoring, before any verdict
was read):** §6.1's "hard gate: §3.5c bands MD-side" quoted the
*twin*-side n̄ band [4.4, 7.1] (target 4.07 + twin-hot bias bracket);
the MD-side realization is the committed g3ring/ke_lown precedent
**n̄ ∈ [3.77, 4.37]** (4.07 ± 0.3) — applying the twin band MD-side
would fail the h405 partner itself (pooled MD n̄ 3.955). The scorer
uses [3.77, 4.37]; n₁ band unchanged [0.19, 0.30].

**Results (Arm A; Arm B verdict-identical at every lin cell —
policy-blocked = 0, marginal class = 0 in-family):**

| cell | a/well/E₀ | n₁ | n̄ | gate | KE₁ | ΔKE₁ vs h405p | deepKE | midHot | trap |
|---|---|---|---|---|---|---|---|---|---|
| lr1 | 27.5/eb0482/0.35 | 0.203 | 3.98 | ✓ | **0.903** | **+0.266** | 2.37 | 1.98 | 0.000 |
| lr2 | 27.5/eb1168/0.35 | 0.202 | 3.98 | ✓ | 0.867 | +0.230 | 1.92 | 1.84 | 0.004 |
| lr3 | 30/eb0482/0.36 | 0.209 | 3.89 | ✓ | 0.830 | +0.192 | 1.67 | 1.74 | 0.003 |
| lr4 | 30/eb1168/0.36 | 0.202 | 3.87 | ✓ | 0.795 | +0.158 | 1.46 | 1.61 | 0.009 |
| lr5 | 32.5/eb0482/0.36 | 0.198 | 4.11 | ✓ | 0.786 | +0.149 | 1.29 | 1.59 | 0.007 |
| lr6 | 35/eb0482/0.37 | 0.202 | 4.03 | ✓ | 0.715 | +0.077 | 0.95 | 1.39 | 0.011 |
| lr7 | 27.5/eb0482/0.37 | 0.235 | 3.38 | ✗ (twin-predicted) | 0.847 | — | 2.18 | 1.84 | 0.000 |
| h405p | capped 5.5/4.4/0.405 | 0.205 | 3.75 | (Δ-only) | 0.637 | 0 | 0.47 | 0.95 | 0.063 |

- **Hard gate:** 6/7 lin cells pass under BOTH policy ends; the gate
  pattern transferred from the twin **7/7** (incl. lr7's miss, twin
  n̄ 3.838 → MD 3.376). h405p itself sits just under the n̄ floor this
  seed (3.748; Arm B flips it in → the *partner* is the only
  policy-blocked row) — disclosed Δ-only per §6.2 item 5, lin
  verdicts absolute and unaffected.
- **K-KE does NOT fire** (every gated cell ΔKE₁ ≥ +0.077, CRN-paired,
  seed-SE ≈ 0.005). **Survivors = all six gated cells; SUCCESS band
  (KE₁ ≥ 0.75, beats the measured (A)-ceiling 0.708) met at five**
  (lr1–lr5). Best: **KE₁ 0.903 at lr1** vs the 1.00 reference peak —
  the twin's R4 forecast (~0.91 if transfer holds) confirmed at MD
  level.
- **Authority box (the §6.2 item-6 product): the capped→lin transfer
  HELD.** Twin KE₁ uniformly ~1–2 % cold in-family (twin−MD
  −0.010…−0.015 eV — the Step-0 level stamp transfers); n₁ transfer
  ≤ 0.006; n̄ twin-hot bias uniform +0.46…+0.59 (small end of the
  bracket); the recalibration branch (§6.2 item 1) was never needed.
- **Both pre-registered kill axes benign:** trap ≤ 0.011 (all bound,
  ZERO marginal — the lin low-v over-drag fear is dead); the CRN fate
  flow runs in reverse: new-trapped = 0 at every lin cell while lin
  *frees* the 52–63 ions h405p traps — the landing is
  cascade-carried, not selection-carried (the §4ee amputation
  diagnostic finds nothing amputated). deepKE mixed
  (hot/strip/neutral along the a-axis, monotone in Φ) → per §6.2
  item 4 reported-no-flag; the twin's 2.0–3.4× hot flag was itself
  ~30 % hot (MD 2.37 at lr1).
- **Honest costs (reported, none gating):** midHot 1.39–1.98 vs soft
  band [0.85, 1.15] and χ²_med 795–2160 vs h405p 267 — the n = 1/2
  heat is bought with **mid-band overheating** (the NB-RQ11-12
  mid-band coordinate measured from the lin side; constant γ
  under-drags the 5–9 Å/ps class relative to cubic). W₁ 0.78–0.97 vs
  h405p 0.767 (±0.13 single-seed disclosure). above-1.15 = 0
  everywhere incl. h405p — the KER width axis stays (B)-side, exactly
  as the (C) record predicts.

**Standing per §0:** an MD-level win on the form-sensitive observables
(KE₁, fate, trap) at matched histogram gates is measured — this
**opens the adoption discussion with real standing**; nothing is
adopted, Tier-0's in-band record and h405/finc1v725 stand untouched.
**NEXT USER GATE:** the adoption discussion + the pre-registered
escalation (one N = 1000 × 5-seed battery at the *single best cell*;
the KE₁-vs-W₁ best-cell convention — lr1 vs lr6 — is part of that
call), and whether the §9 Method-B Tier-0-side row set now becomes
relevant (the ring outcome makes it so).

### 6.4 Adoption-discussion reads (2026-08-10, zero MD, re-read of committed arm-1 rows): the family has TWO measured personalities — the KE₁-buyer AND an uncapped h405 clone

Discussion record from the post-ring adoption thread (user question:
*"can we reach the relatively good result of h405 without capping —
this is so unphysical — maybe by giving up KE₁?"*). Answered by
re-reading the committed `atlas_linsweep.csv`; **no new computation,
nothing adopted.** Recorded because the second personality was
previously durable only as the bare count "17 plateau-convergent" —
the reading itself would have been lost.

**1. The family is a DIAL, not two corners.** a is a one-knob traverse
of a continuous frontier (all rows gated, eb0482 unless noted):

| cell | Φ | twin KE₁ | midHot | deepKE | W₁ | reading |
|---|---|---|---|---|---|---|
| a 27.5, τ 4.8, E₀ 0.35 | 0.64 | **0.892** | 2.16 | 3.38 | 0.673 | KE₁-buyer |
| a 35, τ 4.8, E₀ 0.37 | 0.81 | 0.703 | 1.53 | 1.86 | 0.586 | mid |
| a 42.5, τ 6.4, E₀ 0.31 | **0.985** | 0.615 | **1.003** | 0.644 | 0.686 | **h405 clone** |
| a 45, τ 4.8, E₀ 0.40 | 1.04 | 0.487 | 0.914 | 0.879 | 0.880 | past the clone |
| h405 (capped) twin | 1.00 | 0.624 | 1.069 | 0.844 | 0.703 | incumbent |

**2. The user's proposal is CONFIRMED at twin level.** `a = 42.5 /
eb0482 / τ 6.4 / E₀ 0.31` reproduces the h405 twin vector near
observable-for-observable — KE₁ 0.615 vs 0.624, midHot 1.003 vs
1.069, W₁ 0.686 vs 0.703, n₁ 0.197 vs 0.208, n̄ 4.58 vs 4.48 (deepKE
0.644 vs 0.844, the lin cell slightly colder) — from a **pure linear
form with no cap and no kink**, at the cost of the KE₁ gain, exactly
as proposed. Standing: the §6.3 ring measured the lin-family twin↔MD
box (KE₁ ~1–2 % cold, gate pattern 7/7, n̄ bias +0.46…+0.59), so this
forecast is licensed-transfer, not speculation; an MD check would be
a 1–2 cell CRN spot-check, not a new ring. **Arm 2 already closed the
quadratic question** (c = 0 measured optimal, §3.1/§5.1): the landing
wants maximal flatness, so a ram term only tilts the curve away.

**3. Why this is NOT a reason to retire the cap (the deviation-
accounting argument).** The plateau-matched linear does not remove
the effective element — it *relocates* it into the region where the
program has ab-initio authority. Crossover √(a/b) = 4.11 Å/ps
(b = 2.5153, Tier-0): the a = 42.5 line runs **×4.2 the Tier-0 cubic
at v = 2, ×1.9 at v = 3, and 31 % SOFT at the band top 4.95** — i.e.
wrong by factors precisely inside the TDDFT-calibrated window, which
is the Tier-0 linear rejection (n̂ = 2.927; a → 0 whenever a higher
power competes) restated in system terms. `capped_cubic` by contrast
is *exact* where the traces speak and parks its whole free element
above 4.95, where TDDFT is permanently infeasible (D0 §2 domain
input) and only the experiment has authority. **Both forms carry a
fiction; the cap hides it where nothing can see it, the smooth line
hides it where the traces can.** By the validation hierarchy that is
a worse trade, bought with aesthetics.

**4. Cap-physics assessment (recorded so "the cap is unphysical" is
not re-litigated from scratch).** v_c 5.5 Å/ps ≈ **Mach 2.3** against
the He sound speed (~2.4 Å/ps); constant force above it = constant
energy loss per Å, the signature of wave-drag / vortex-shedding
regimes — saturating drag on a supersonic bubble is respectable
physics, not a fudge. What *is* idealized is the sharp corner: any
smooth saturation family (e.g. `F = b·v³/(1 + (v/v_c)^m)`, cubic at
low v, saturating at high v) departs from the Tier-0 cubic near the
band top unless m is large — and **m → ∞ is exactly the cap**. The
kink's sharpness is therefore *forced* by the two authorities
(in-band trace lock + the experiment's demand for a flat fast-class
force), not chosen.

**5. POSTED, NOT TAKEN — the smooth-saturation twin read.** If the
kink itself is the objection, the physically-motivated instrument is
not the polynomial family but a Hill-type smooth saturation (b
Tier-0-locked; free (v_c, m)), twin-scanned at zero MD to measure
*how sharp the corner has to be* — i.e. the largest smoothing the
in-band lock and the histogram jointly tolerate. Turns an aesthetic
into a measurement. Not adjudicated; behind its own trigger. (Note:
this is NOT a §3.5h p_tail re-open — p_tail softens the *tail
exponent* at fixed corner; this softens the *corner* at fixed tail.)

**Standing:** nothing adopted; both personalities are counterfactual
readings of a committed twin artifact. Tier-0, h405, and finc1v725
stand.

### 6.5 h405-clone MD battery — DESIGN FROZEN (2026-08-12, discussion-adjudicated; `[PROCEED TO IMPLEMENTATION]` given)

The §6.4 item-2 reading is a **twin** forecast. This battery measures it
in MD: does the uncapped, kink-free pure-linear cell actually reproduce
the h405 landing? The §6.3 ring never went there — it covered
a ∈ [27.5, 35] at τ 4.8, and the clone sits at **a 42.5 / τ 6.4**,
outside the measured twin↔MD bias box in *both* swept axes. That
extrapolation is the reason MD is required and, per CL-P6 below, is a
program product regardless of the verdict.

**Cell (one) × seeds (three), N = 500.** `pure_linear a = 42.5 amu/ps /
E_bind 0.0482 (eb0482) / τ 6.4 ps / E₀ 0.31 eV`, corrected geometry,
standing pins otherwise — the committed `atlas_linsweep.csv` row
verbatim (Φ 0.985; twin KE₁ 0.6153, midHot 1.0144, deepKE 0.6444,
W₁ 0.6857, n₁ 0.1966, n̄ 4.576, trap-floor 0.0476, supp 0.0456).
Seeds **20260731** (the §6.3 ring seed — the clone thereby slots into
the committed ring table as a same-N, same-seed 8th cell directly
comparable to lr1–lr7, *and* CRN-pairs for free against the committed
`h405p` partner) **+ 20260812 + 20260813** (fresh, audited unused —
committed generators end at 20260734). ≈ 45 min at concurrency 3.

**No new h405 MD is needed, and that is measured, not assumed.** The
committed g4 Step-2 battery (h405, N = 1000, seeds 20260730–34) has a
seed spread small enough that the comparison target is a constant:

| | KE₁ | midHot | n₁ | n̄ | W₁ | trap | deepKE |
|---|---|---|---|---|---|---|---|
| h405 5-seed mean | 0.6406 | 0.9464 | 0.2094 | 3.877 | 0.765 | 0.078 | 0.503 |
| per-seed SD | **0.0024** | 0.0074 | 0.0102 | 0.047 | 0.030 | 0.0076 | 0.028 |

Cross-N licence: the committed N = 500 `h405p` agrees with the N = 1000
pooled battery on every intensive observable within ~2σ (KE₁ 0.6371 vs
0.6406, midHot 0.9466 vs 0.9464, W₁ 0.7675 vs 0.7653, n₁ 0.2054 vs
0.2094; the two ~2σ entries are trap 0.063 vs 0.078 and n̄ 3.748 vs
3.877). **Exception — χ²_med is N-extensive** (267 at N = 500 vs 461 at
N = 1000): the clone's χ² is compared to `h405p` and per-seed only,
never to the pooled N = 1000 number.

**Measured resolution at 3 × N = 500 (pooled 1500 ions).** Per-seed SDs
are taken from the battery table above, scaled ×√2 to N = 500, and
SE = SD/√3:

| read | pooled SE | band | quality |
|---|---|---|---|
| CL-P2 ΔKE₁ | ≈ 0.003 | ≤ 0.05 | decisive |
| CL-P3 midHot | ≈ 0.006 | ≤ 1.15 (forecast ≈ 0.92) | decisive |
| CL-P1 n̄ | ≈ 0.038 | [3.77, 4.37] (forecast 3.99–4.12) | decisive |
| CL-P1 n₁ | ≈ 0.008 | ≥ 0.19 (forecast **0.19–0.20**) | **marginal — the failure mechanism** |
| CL-P4 trap | ≈ 0.006 | reported | fine |
| CL-P5 W₁ | ≈ 0.024 (ring disclosure ±0.13/seed) | reported-only | fine |

The N cut costs nothing on the two headline reads; it costs only on n₁,
which is where the run was already most likely to fail.

**Frozen reads (pre-registered before launch):**

- **CL-P1 — hard gate**, pooled ensemble, both retained-policy arms:
  n₁_solv ∈ [0.19, 0.30] ∧ n̄_det ∈ [3.77, 4.37] (the §6.3 MD-side
  realization). **Three-way outcome, forced by the N cut:** a pooled
  value within 1 SE of a band edge returns **gate-marginal** — neither
  pass nor fail — where SE is the *measured* per-seed SD/√3 of this
  battery, not a forecast. Per-seed verdicts are reported with their
  spread; no single seed decides.
- **CL-P2 — equivalence:** |ΔKE₁| ≤ 0.05 eV vs h405. The ring's K-KE
  kill is **inverted here** — this run is an equivalence test, the
  clone *wants* ΔKE₁ ≈ 0, and a large positive ΔKE₁ would mean the cell
  is not a clone, not that it won.
- **CL-P3 — the cost read:** midHot ≤ 1.15 (the soft band top). The
  §6.3 ring's honest cost was midHot 1.39–1.98; the twin says the clone
  returns to 1.014. This is the read that decides whether mid-band
  overheating is an a-dial artifact or intrinsic to constant γ.
- **CL-P4 — trap:** reported decomposed (bound/marginal); no kill
  (§4.1). Flagged if the marginal class is non-zero (the ring measured
  zero in-family) or Δtrap vs h405 > +0.05. Twin floor here is 0.0476,
  4× the ring cells' — the a-axis raises the low-v over-drag.
- **CL-P5 — W₁: reported-only** (user-adjudicated 2026-08-12). W₁ never
  gated anywhere in this arm and the histogram-level landing was
  measured form-blind in §6.6; a pooled loss > +0.10 vs h405 is recorded
  as an honest cost line beside midHot/χ², and cannot overturn
  CL-P1..P3. Standing expectation: the lin family's twin W₁ ran ~0.2
  cold in the ring, so a clone W₁ ≈ 0.85 vs h405's 0.766 is the
  *forecast*, not a surprise.
- **CL-P6 — authority-box extension:** the measured twin−MD deltas at
  a 42.5 / τ 6.4 are committed beside the MD values and compared to the
  ring-measured in-family bias ranges (KE₁ twin−MD −0.010…−0.015; n₁
  ≤ 0.006; n̄ +0.46…+0.59). Outside those ranges ⇒ the bias vector is
  declared a-/τ-dependent. **This product lands regardless of the
  verdict** — it is the first bias data outside the ring's box.

**Failure branch, pre-registered (no improvisation).** Any of CL-P1
miss / CL-P2 / CL-P3 miss routes to **§6.2 item 1's recalibration
machinery, not a new ring**: CL-P6's measured bias vector is applied to
the committed `atlas_linsweep.csv` and the twin re-scanned
bias-corrected — **zero MD** — which *names* the corrected clone cells.
A second MD iteration is a **new user gate**, never automatic.

**Instrument.** `scripts/gen_tier2atlas_linclone.py` — 3 runs built
through `gen_tier2atlas_linring.build_cell` with only the seed replaced
(rule 1: no second copy of the cell construction), reusing that module's
`verify_corrected_geometry` / `verify_against_reference` guards and free-form
provenance posture (`extraction_method="free_form"`, no §6.5.1 binding
stamp, `allow_unvalidated_binding_pairing`). Scorer
`scripts/post_processing/tier2atlas_linclone_table.py` extends the
§6.3 ring scorer (Arm A/Arm B via the §3.5b bracket, KE + width columns,
trap decomposition, CRN fate flow, twin join) and pools the three seeds
with `pool_confirmation_reads`. Artifact:
`atlas_linclone_table.csv` — per-seed + pooled, full observable vector
incl. every KE column, twin-forecast columns beside MD (memory rule).

**Oracles, all before any new number:**

1. **LC-P1** — the frozen clone twin row reproduces **string-exact**
   from the committed `atlas_linsweep.csv`, plus the §6.3 ring's own
   LR-P1 (7 lin rows + h405 authority anchor) re-run unchanged.
2. **CRN guard (the pairing *is* the guard).** The seed-20260731 cell's
   cfg diff vs the committed `h405p` cfg must contain neither `seed` nor
   `num_molecules` — their identity is what makes the pair CRN. The two
   fresh-seed cells must diff against the 20260731 cell in exactly
   `{"seed"}` (the g4 Step-2 precedent).
3. **Partner anchor** — the committed `atlas_linring_table.csv` `h405p`
   row is re-derived by this scorer before it is used as the Δ target
   (gate-on-committed-artifacts rule).
4. Scorer-drift + KE-path oracles inherited from the ring scorer;
   default path inert (`pure_linear` already suite-covered).

**Outcome semantics (§0 unchanged — nothing is adopted by this run).**
CL-P1 ∧ CL-P2 ∧ CL-P3 all passing = *an uncapped, kink-free pure-linear
cell reproduces the h405 landing in real MD*. That **opens** the
form-choice on the §6.4 item-3 deviation-accounting grounds; it does
not settle it, because Tier-0's in-band rejection (n̂ = 2.927) is
untouched by any Tier-2 landing. A gate miss on the n₁ floor means the
twin transfer breaks outside the ring's measured box. midHot > 1.15
means the mid-band overheating is **form-intrinsic to constant γ**, not
an a-dial artifact — which is the cleanest available answer to the
question that opened §6.4.

**Findings — EXECUTED 2026-08-12 (3 × N = 500; full record in
`TIER2_SENSITIVITY_ATLAS_FINDINGS.md` "§6.5 h405-clone MD battery" and
D0 §1; artifact `atlas_linclone_table.csv`).** All four oracles green
including the CRN guard. **The clone is a KE-equivalent, not a
landing-equivalent:**

- **CL-P2 PASS** — pooled KE₁ 0.6348 vs the h405 battery 0.6406,
  ΔKE₁ **−0.0059 eV** (CRN pair −0.0016) against a ±0.05 band.
- **CL-P3 PASS, and it inverts the §6.3 cost** — midHot **0.891**,
  *cooler* than h405's 0.946 and far below the ring's 1.39–1.98;
  per-seed χ² 173/234/241 vs h405p's 267. **Mid-band overheating is an
  a-coordinate, not a property of constant γ.**
- **CL-P1 `gate-marginal`, both policy ends** — n̄ 3.993 passes; n₁
  0.1921 sits *on* the 0.19 floor (gap 0.0021, per-seed SD 0.0181 ⇒
  ~72 seeds of N = 500 to separate). Boundary is the measurement, not a
  resolution shortfall.
- **CL-P4 benign** — trap 0.062 *below* h405's 0.078, marginal 2 ions,
  CRN flow new-trapped 1 / freed 9; supp 0.051 vs 0.216.
- **CL-P5 the cost that names the result** — W₁ **1.079 vs 0.765**
  (+0.314) while both histogram *moments* sit at/inside band: the
  moments match, the shape does not.
- **CL-P6 the box does NOT extend here** — n₁/n̄ transfer inside the
  ring ranges, KE₁ (−0.0195) outside, and **twin W₁ collapses (0.686 →
  MD 1.079, −0.39 bias vs the ring's ≈ −0.2)**. The §6.4 clone claim
  rested on a 0.017 twin-W₁ gap; the transfer error here is ~20× it.

**Standing:** nothing adopted. §6.2 item 1 is pre-registered for a
*miss* and this is not one, so routing (accept-marginal / re-site off
the n₁ floor / recalibrate) is a **USER GATE**. Note the target has
changed shape: a re-sited search is no longer for an h405 *clone* — the
twin coordinate that promised one is now measured untrustworthy at this
corner — but for a linear cell that lands the histogram *shape*, which
no twin column currently predicts here.

### 6.6 The (a, τ) joint ring — DESIGN FROZEN (2026-08-12, discussion-adjudicated; `[PROCEED TO IMPLEMENTATION]` given)

User question after §6.5: *"what sweep can improve histogram and kinetic
energy simultaneously (that KE₁ doesn't land fully is accepted)"*. This
ring answers it, and it exists because the §6.5 anatomy identified a
mechanism the program had never isolated.

**What the §6.5 anatomy measured (all zero-MD reads of committed
artifacts; the reason this ring is well-posed):**

1. **The clone's W₁ damage is entirely a tail deficit.** Signed CDF-gap
   decomposition: lead mass (sim short in the tail) 0.880 for the clone
   vs 0.439 h405p / 0.482 lr6 / 0.725 lr1; the clone's *low-n* side is
   the best of the four (0.200). At n = 12/13/14 it holds
   0.0116/0.0064/0.0041 vs reference 0.0222/0.0197/0.0187.
2. **The tail was evaporated, not trapped.** lr6, h405p, clone-s1 all
   share seed 20260731, so the per-ion CRN flow is exact: of the 237
   ions at n ≥ 10 under lr6, only **43** are trapped in the clone; 194
   survive detected and the heavily-dressed group falls from mean
   n **18.3 → 9.6**. Low-n ions *gain* (+0.49/+0.25/+0.12 for lr6-n of
   0–2/2–5/5–10). The clone compresses the distribution from both ends
   at a correct mean — right n̄, wrong width.
3. **E₀ is NOT the width knob.** lr1 → lr7 (a 27.5, E₀ 0.35 → 0.37)
   doubles suppression 0.093 → 0.191 and moves n̄ 3.98 → 3.38, while W₁
   moves 0.9477 → 0.9459 — i.e. nothing.
4. **`a` improves W₁ monotonically at τ 4.8** (0.948 → 0.877 → 0.856 →
   0.785 over a 27.5 → 35), taking midHot 1.976 → 1.393 and χ²
   2160 → 795 with it, and costing only KE₁ (0.903 → 0.715).
5. **By elimination, τ owns the clone's +0.35 W₁.** Extrapolating (4) to
   a 42.5 predicts W₁ ≈ 0.72; the clone measured 1.079, and τ is the
   only remaining moved knob.
6. **τ ↓ is not available** (checked before pinning any cell, and it
   killed the first version of this design): at a 42.5, τ 3.2 gives
   twin n̄ **6.7–16.0** across the whole E₀ grid (τ 2.4: 10.9–17.7).
   Shorter τ = faster cooling = *less* evaporation = far more retained
   helium. With the validated n̄ transfer (−0.55 MD-side) those cells
   land n̄ ≈ 6+ against a [3.77, 4.37] band — not marginal, a different
   distribution.

**Cells: 6 × N = 500, seed 20260731 — deliberately the §6.3 ring seed**,
so every new cell is CRN-paired ion-for-ion with lr1–lr7, `h405p` **and**
clone-s1. Cell-to-cell differences are therefore common-random-number
paired (far tighter than absolute noise), while absolute values carry the
per-seed SD measured by the §6.5 battery in this exact family at N = 500:
**W₁ 0.031, n₁ 0.0181, n̄ 0.099, KE₁ 0.0022, midHot 0.0024, deepKE 0.040,
trap 0.0066.** All cells `eb0482` (better than eb1168 on both W₁ and KE₁
throughout the ring).

| # | cell | role | projected MD (measured biases KE₁ +0.02, midHot ×0.88, n̄ −0.55) |
|---|---|---|---|
| `s37` | a 37.5 / τ 4.8 / E₀ 0.38 | a-curve | KE₁ 0.66, midHot 1.18, n̄ 4.01, n₁ 0.203 |
| `s40` | a 40 / τ 4.8 / E₀ 0.39 | a-curve | KE₁ 0.60, midHot 1.03, n̄ 3.96, n₁ 0.203 |
| `s425` | a 42.5 / τ 4.8 / E₀ 0.39 | **decomposition** — the clone's a at the ring's clock | KE₁ 0.56, midHot 0.94, n̄ 4.17, n₁ 0.193 |
| `s45` | a 45 / τ 4.8 / E₀ 0.40 | does the a-trend saturate? | KE₁ 0.51, midHot 0.82, n̄ 4.09, n₁ 0.194 |
| `d2` | a 35 / τ 6.4 / E₀ 0.30 | **τ mirror** — lr6's a at the clone's clock | KE₁ 0.79, midHot 1.26, n̄ 3.87, n₁ 0.192 |
| `j1` | a 40 / τ 6.4 / E₀ 0.31 | **the joint candidate** | KE₁ 0.67, midHot 0.98, n̄ 3.84, n₁ 0.211 |

With lr6 (a 35, τ 4.8) and the clone (a 42.5, τ 6.4) already committed,
this gives a **5-point a-curve at τ 4.8** and a **3-point τ contrast at
a 35 / 40 / 42.5**. (`j1` misses the *twin*-side n̄ band at 4.392 — the
twin band's floor is 4.4 — but lands inside the MD-side band under the
validated −0.55 transfer; recorded here so its twin `gate = 0` is not
later misread as a prediction of failure.)

**Frozen reads (pre-registered before launch):**

- **JR-P1 — the decomposition (this ring's primary question).** τ owns
  the W₁ damage ⟺ `s425` W₁ ≤ 0.85 **and** `d2` W₁ ≥ 0.95. Killed if
  `s425` ≥ 1.0 (then a = 42.5 itself broke the a-trend and the dial
  saturates) or `d2` ≤ 0.85 (then τ 6.4 is harmless at a 35 and the
  clone's damage is an a×τ interaction). Any other combination is
  reported as mixed and claims nothing.
- **JR-P2 — the a-curve at τ 4.8**: W₁, KE₁, midHot, χ² and the tail
  diagnostic vs a ∈ {35, 37.5, 40, 42.5, 45}; the saturation point (if
  any) is reported.
- **JR-P3 — the joint target (the user's actual question).** SUCCESS =
  any cell beating h405 on **both** W₁ (< 0.767) **and** KE₁ (> 0.637)
  with midHot ≤ 1.15. `j1` is the pre-registered favourite and it wins
  only if JR-P1 comes back "a×τ interaction" rather than "τ".
- **JR-P4 — the fallback is a result, not a failure.** If no cell meets
  JR-P3, the Pareto front of (W₁, KE₁) over the whole measured family is
  reported and named as the boundary: KE₁ and histogram are then
  **not simultaneously improvable in this family**, with lr6 (equal W₁,
  KE₁ +0.078) as its ceiling.
- **JR-P5 — the tail diagnostic** (the mechanism read that made this
  ring possible): scored n ≥ 10 fraction, max n, and the per-ion CRN
  shell flow vs lr6, per cell.
- **Gate: reported, never selecting** (user call, 2026-08-12) — pass /
  gate-marginal / fail under both retained-policy arms, with the §6.5
  three-way band verdict. Rationale: every cell in this family sits at
  n₁ 0.19–0.21 against a 0.19 floor because of the shelved **(C)**
  source-side deficit that h405 shares, so gating on it discriminates
  nothing *within* the family and only manufactures marginal verdicts.
  Selection is on W₁ + KE₁, with midHot as the cost line.

**Honest prior, on record before launch:** the τ 4.8 cells that are
projected to fix W₁ land KE₁ ≈ 0.56–0.60, i.e. **below** h405's 0.637 —
on that line the trade is real and not broken. `j1` is the only cell
projected to break it. A JR-P4 outcome is therefore a live and
respectable possibility, not a disappointment.

**Instrument.** `scripts/gen_tier2atlas_linjoint.py` — same construction
path as §6.5 (cells built through `gen_tier2atlas_linring.build_cell`,
its geometry/free-form guards reused; rule 1), LJ-P1 oracle freezing all
six committed twin rows string-exact, and the CRN guard extended: with
the seed *shared* with `h405p`, `seed` and `num_molecules` must both be
absent from every cell's diff. Scorer
`scripts/post_processing/tier2atlas_linjoint_table.py` joins the
committed lr6 / clone-pool / h405p rows for the curve and the contrast,
adds the JR-P5 tail diagnostic, and commits the full vector to
`atlas_linjoint_table.csv`.

**Standing:** nothing adopted; instrument runs only. Tier-0's in-band
rejection (n̂ = 2.927), h405 and `finc1v725` stand.

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
