# Tier 2 — Parameter Influence Reference (compact)

> **What this is.** The D0 deliverable of `TIER2_SENSITIVITY_ATLAS_PLAN.md`:
> a **pure distillation** of `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee,
> I1–I100) into per-knob influence statements. **No new analysis** — every
> number here is a pointer-backed quote of a measured result. The findings
> doc stays the archive/provenance; this doc is the navigation layer.
>
> **How to read an entry.** *Role* = where the knob enters the mechanism.
> *Class* = CALIBRATION_MAP calibration class (Sourced / Bounded / Free /
> Derived). *Influence* = measured effect on the observables, with
> magnitudes and insight pointers (I-numbers → the findings-doc register;
> § → findings-doc sections). *Couplings* = structural entanglements with
> other knobs. *Status* = locked / standing / open / **GAP** (no measured
> influence yet — atlas target).
>
> **Two chapters, two questions.** §1–§16 answer *what each knob does*
> (influence). **§17 — the physical-sensibility ledger** — answers *which
> parts of the model are physics and which are scaffolding*, with what
> would retire each. The program's purpose is model understanding and
> physical soundness in general; RQ11 is one symptom, not the target
> (plan §0, user 2026-07-26).
>
> **Living doc.** Each sensitivity-atlas stage merges its results in and
> closes its GAP marker. Last built: 2026-07-23 (distillation of
> I1–I100); merged 2026-07-24: D4 Step 1 + §6.6 counterfactual
> (twin + MD spot-check — the form-blindness result, §1); §6.7 item 1
> (lq sanity battery — over-suppression, §1) + item 2 (E_bind scan —
> over-suppression is the FORM; **E_bind GAP closed, §9**); merged
> 2026-07-26: the **D2b as-built provenance audit** (§15 — which
> sampler actually ran in which production, code-verified and
> distribution-verified) + the **parent-document geometry anchor**
> (§15.4/§15.5 — solvation-depth reproduction and the source-condition
> droplet range; the influence GAPs in §14/§15 stay open).
>
> **Standing point (context for every "standing" value):** finc1v725 —
> `capped_cubic` v_c 7.25 / p_tail −1, τ 3.2 ps, E₀ 0.27 eV, `rq4graded`
> ladder, Landau v_L 0.58 Å/ps, locked Tier-0 b, E_bind 0.1168 eV,
> margin 3 Å, N = 5000 pooled battery reference (§4cc).

---

## 0. Summary table

| knob | class | headline influence | status |
|---|---|---|---|
| drag form + b (in-band) | Derived (Method-B/TDDFT) | owns KE scale + fragment composition; **the landed observables do NOT identify the form — cubic and quadratic systems both land in full MD (§6.6)**; deep-KE lever re-pointed at the 5–9 Å/ps mid-band | locked in-band by Tier-0 traces ONLY; landing form-blind (MD-measured) |
| v_c, p_tail (capped tail) | Free | v_c owns the whole mid-bin KE curve; joint-landing basin v_c ∈ [7.25, 7.5]; p_tail is not a second lever | standing 7.25 / −1 |
| τ (cooling clock) | Bounded | descent clock; joint closure only as a (v_c, τ) pair; race coordinate with f_int under the gate | standing 3.2 ps |
| E₀ = E_int(0) | Bounded [0.2, 0.5] eV | fate-cliff position (bare↔shell split); re-landed at solvation scale 0.22–0.27 under full geometry | standing 0.27 eV |
| κ (ladder steepness) | Free | near-dead staircase lever (inverted + normalization-capped) | standing (rq4graded supersedes) |
| picture (electronic) | Free (selection) | magnitude-degenerate ≤ 4 %; discrimination belongs to the size distribution | standing `statistical_mixture` family |
| f_int | Derived coordinate | superseded by absolute E₀; gated: only the race margin Δ× is physical | retired as free knob |
| f_ret | Bounded (prior small) | pinned 0.1 throughout; never swept | **GAP** |
| λ₀ pickup (+ p cap) | Sourced + Bounded (p conditional-Free) | re-filling ≈ 7 ions/100; p = 1 de-suppresses to the twin and lifts n₁ KE toward experiment | standing p = 1 |
| s_eff (RRK dof) | Bounded | cascade *rate*, arm-conditional (8 ungated / 30 gated); detector-compressed except the n = 1 bin | standing 8-family |
| ladder D₀(n), Σ | rungs Sourced; shape RQ4-graded | owns n₁/small-n shape; rq4graded > floor1 ≫ 2σ; deep tail Σ-coupled ×20 too stiff to move deep KE | standing rq4graded |
| Landau v_L | Sourced 0.58 Å/ps | bit-flat on scored surface across 0.30–0.58; acts only on the retained class | standing 0.58, quiet |
| per-shed ε | Bounded small (NB-RQ23-1) | refuted as RQ11 owner; ε ≈ 1–2 meV/shed sub-dominant compatible | ε = 0 standing |
| shed momentum convention | convention (two-valued) | cold-shed injects ×1.611 KE over a full strip; histogram convention-blind; bare-bin KE reads the fragmentation convention | co-moving basis for twin parity; RQ3-coupled |
| birth margin | pinned convention (3 Å) | **the sensitive robustness lever**: 6 Å moves n₁_solv −0.171 (3.4σ) | standing 3 Å; in-tier open item (I88) |
| E_bind (ion–droplet well) | Derived (joint Method-B) | **swept §6.7 item 2**: trap +0.058/0.1168-step (clean well lever), n̄ −0.50, midHot −0.076; over-suppression is the FORM, not the well | measured (§9) |
| droplet geometry (R × r) | controlled (Axis A G1 11-cell grid); **size externally anchored** (§15.5) | birth depth is the physics knob, R a selection knob; the landing needs the size *distribution* (pinning R̄ alone: W₁ 0.571 → 0.813); ≈ 95 % of detected-size variance geometry-inherited; deepKE crosses 1 at birth depth ≈ 11–15 Å; trap → 0.40–0.57 at the anchored radii | **GAP closed** (§14); **G2 ADOPTED 2026-07-27** — the corrected geometry is the target, re-arbitration pending (G3) |
| sampling laws (size + position) | theory-laden legacy ports, **partly bypassed in the drag branch** | provenance audited (§15): production uses the analytic ⟨N⟩ = 2000 prior + uniform_volume, *not* the legacy pickup MC + Boltzmann; E_solv 14 vs 30 meV discrepancy is inert here; **⟨N⟩-pin influence measured at ensemble level by grid re-weighting (§15.7, zero MD): the corrected ensemble breaks the landing** (trap 0.31–0.42, W₁ ≈ 9.0–9.8, deepKE 1.80–1.90) | ⟨N⟩ pin **measured** (§15.7); distribution A/B remainder open (D2b) |

---

## 1. Drag law — in-band form + coefficient b

| | |
|---|---|
| role | continuous friction γ(v) = ρ̂·b·v² (force ρ̂·b·v³) inside the bubble; owns all in-window dissipation |
| class | Derived — Method-B extraction from 9/18 Å TDDFT traces (bands: 18 Å 2.54–3.02, 9 Å 2.83–4.95 Å/ps); locked Tier 0 |

**Influence (measured):**

- Owns the KE *scale* and the fragment *composition*: n₁/bare occupants
  are near-edge births losing 0.8–1.2 eV to drag vs 2.6–2.8 eV for deep
  bins (I60). Geometry owns the KE *shape*; the envelope [current-law,
  ballistic] brackets experiment everywhere (I43).
- Over-dissipates at production speeds: detected model speeds undershoot
  the VMI I⁺He peak ~2–2.5× (I40); measured K₂.₇₀ = 0.746 with 41.7 % of
  exposure beyond the calibrated band (I37).
- The RQ11 deep-bin cold tail (sim/ref 0.81 → 0.38 over n = 10–17, a
  slope not a scale, I92) accrues in the **in-band pure-cubic phase at
  v ≈ 2–7 Å/ps, t ≈ 2–14 ps, in-bubble** (I90); with candidates
  (i)/(ii)/(iii) refuted, the in-band law at low v is the surviving
  owner-candidate (I96) — ~1.5 % of the deep survivors' dissipation.
  **Revised by the §6.6 MD spot-check (NB-RQ11-12): the velocity window
  that actually moves deep-bin KE is the upper-mid band v ≈ 5–9 Å/ps,
  not the theorized sub-2.5 low band.** Measured: the lq cell at
  standing τ (γ lower than capped-cubic above ~5 Å/ps, cap at 9)
  lifts deep bins ~+0.1–0.2 (qcb 1.02…0.45 vs pooled cubic
  0.81…0.38) *despite* carrying ×2.5 the drag at v = 2 — the sub-2.5
  region barely registers on deep KE, the 5–9 window dominates. The
  deficit's slope survives every form tested (form-robust).
  Direction-only at N = 500 deep-bin counts; N = 1000 confirmation is
  the escalation path.
- Low-v form space mapped in the twin (§4ee): hard/erf floors are
  perfectly speed-selective and act by un-trapping the marginal retained
  class (trap 0.042 → 0); the sub form at v_f ≈ 1.5–2.0 is the twin's
  landing window (deep → 1.06…0.83, midHot ≤ 1.14); the parameter-free
  shifted-cubic at locked b over-corrects globally (all bins ×1.2–1.4,
  midHot 1.40 — though W₁ improves 0.678 → 0.502) (I98/I99).
- Literature: above ~1 Å/ps in a droplet the drag *shape* is
  unconstrained — TDDFT is the only authority; the (v−v_L)³ law is
  near-threshold/negative-ion/mean-of-sawtooth and does not transfer
  (I97, I100, NB-RQ11-1..9).
- Atlas D4 Step 1 (2026-07-23, artifact reuse, zero MD): the three-mode
  form table stands as compiled — lq rejected (shared Δobj +0.0339,
  held-out 9 Å FAIL) with its low-v limb ×2.54 the production drag at
  v = 2 (wrong direction for RQ11); power law cubic-equivalent
  (n̂ = 2.927). Tail re-inspection: **no drag-dominated trace samples
  below v ≈ 2.5 Å/ps** — the sub-2.5 shape is TDDFT-blind pending the
  collaborator ask. The subtractive gated form is trace-constrainable
  in-band (×0.65–0.38 at 2.54 for v_f 1.5–2.0 vs the whole 18 Å window
  at 2.54–3.02); its re-fit is a pending decision. Full tables:
  `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` (NB-RQ11-10).
- **The landed *histogram* does not identify the in-band form — the
  histogram landing is form-blind, measured at battery scale (§6.6 +
  the §6.7 N = 1000 paired battery, 2026-07-24).** W₁, midHot, n̄, and
  n₁_solv are all form-comparable (pooled lq W₁ 0.548 vs cubic 0.571,
  midHot 1.023 vs 1.014). "Form-blind" is a **W₁-only** statement, and
  the §6.7 item-2 E_bind scan sharpens why: item-1's midHot/n̄ *match* is
  a **(form, well) co-compensation**, not form-insensitivity. At a
  *matched* well (0.1168) the forms diverge strongly — lq is colder
  (midHot Δ−0.070, ~35σ) and smaller (n̄ Δ−0.57, ~20σ) than cubic (§6.7
  item-2 firm-up); item-1 saw them equal only because lq's shallower
  0.048 well cancelled its intrinsic coldness. W₁ alone stays form-blind
  even at matched well (Δ+0.009, ns), but the fate split (supp/trap) and
  the KE/size (midHot/n̄) *do* separate the forms — so this is not a
  claim that every observable is form-insensitive.
  A quadratic-form system with its own Tier-0 artifacts (a ≈ 0,
  c = 12.79, E_bind 0.0482) and re-arbitrated downstream knobs lands
  the full Tier-2 observable vector in *full MD*: three cells
  (v_c 8.8–9.0 × τ 3.2–3.5), all landing, including the cell the twin
  predicted to fail by ~5 seed-SD (the twin's needle-thin lq basin was
  a frozen-chord artifact; the MD basin is wide — mechanism feedback
  smooths the (τ, E₀) knife-edge). Consequences, program-level: the
  experimental landing can never be cited as evidence for any drag
  form (cubic included); form authority lives exclusively in the
  Tier-0 trace instruments (held-out above all, where lq fails 0.699);
  and a large share of the landing quality is owned by the downstream
  mechanism's compensation capacity (through K, I47), not by the drag
  law. This is the measured closure of the §6.1 anti-circularity
  question.
- **The §6.6 read that a quadratic law lands *better* on the KE surface
  does NOT replicate — REFUTED by the §6.7 N = 1000 battery (2026-07-24).**
  The §6.6 same-N same-seed comparison (three lq cells halving χ²_med,
  64.7 / 61.3 / 82.6 vs cubic finc1v725 125.7) was a **single-seed
  artifact** (seed 20260721, N = 500). Paired across the 5 fresh battery
  seeds at N = 1000, lq's median-anchored KE χ²_med is *worse* than cubic
  on **every** seed (paired Δ+13.3±10.0; pooled 277 vs 242). The KE
  advantage was not even seed-robust — it evaporated *before* the E_bind
  confound is addressed.
- **lq over-suppresses — a resolved drag-form influence on the fate
  split (§6.7 battery + item-2 scan).** supp 0.213 vs cubic 0.187 at the
  battery (paired Δ+0.0255 ± 0.0041, ~6σ), with under-trapping. **The
  §6.7 item-2 E_bind scan DISCHARGES the confound: the over-suppression
  is the FORM, not the well** — at a *matched* well (0.1168) lq still
  over-suppresses by Δ+0.034 ± 0.004 (~8σ, N = 1000 × 3 seeds); the well
  adds only +0.010 more. lq also runs colder/smaller at matched well
  (midHot Δ−0.070, n̄ Δ−0.57) — the intrinsic quadratic-form signature
  that item-1's shallow co-calibrated well masked. (χ²_med stays
  seed-noise-dominated and does not resolve form vs well even at N = 1000.
  Deep well 0.154 over-retains: 2/3 cells trip the handover guard as
  ~2/2000 ions never decouple in 8000 ps.)
  **Atlas stance: no adoption** — the Tier-0 traces reject this form
  (held-out 0.699 FAIL); these entries document what the observables
  *reward/penalize*, not a form change.

**Couplings:** jointly extracted with E_bind (§9 below); KE↔histogram
anti-correlate through the exposure integral K (I47) — any drag change
re-opens the (v_c, τ, E₀) arbitration (I99).

**Status:** locked in-band **by the Tier-0 trace instruments only** —
after §6.6 (twin + MD, 2026-07-24) the Tier-2 landing carries zero
form-discrimination power (see the form-blindness influence bullet
above), so the traces are not merely the *best* authority but the
*sole* one; lq stays rejected there (held-out 0.699 FAIL, n̂ = 2.927).
RQ11's drag-shape axis is re-pointed at the 5–9 Å/ps mid-band
(NB-RQ11-12). Open items: subtractive-form re-fit decision
(Step-1 gate fired), the variable-mass Method-B re-run
(§6.6 pre-registered follow-up), the collaborator ask
(sub-2.5 + production-kinematics traces). Chronology + full tables:
`drag_migration_log_tier2.md`, `TIER2_SENSITIVITY_ATLAS_FINDINGS.md`.

## 2. Capped tail — v_c and p_tail

| | |
|---|---|
| role | production arbitration's above-band extension: cap at v_c with tail exponent p_tail |
| class | Free (the program's genuinely fitted knobs, arbitrated by experiment) |

**Influence (measured):**

- **v_c owns the entire mid-bin KE curve**: v6.5 → v7.5 at p−1 moves
  midHot 1.74 → 0.95 with n2–n8 within ~10 % (I78). The KE χ² argmin was
  right-censored at the bracket top until the taper was understood (I71).
- **Joint landing is a basin, not a knife-edge**: v_c ∈ [7.25, 7.5] at
  (τ3.2, E₀0.27, rq4graded); W₁_solv bowls at 7.25 (0.424), n₁_solv
  crosses 0.31 at ≈ 7.25, midHot crosses 1.0 at ≈ 7.2; closes above
  ~7.75 (whole curve ×0.68–0.88 at v8.0) (I80/I81).
- **p_tail is the wrong second lever**: a more-negative exponent re-heats
  the whole KE curve *and* re-strips — KE and histogram move together,
  not orthogonally (I78). p = 1 tails excluded outright (I51).
- Three observables prefer three v_c (KE χ² ≥ 8.5, n₁ anchor ≈ 6.2–6.3,
  W₁ ≈ 7.25) — a drag-form-*shape* statement, not a calibration error
  (I72).
- Hard cutoff (tail removed) excluded: 99.7 % suppression → bare at
  v_c ≤ 10 (I48).

**Couplings:** with τ (joint closure exists only as a pair, I49); with
the ladder — orthogonal (ladder is KE-neutral: midHot flat across the
ladder family at every v_c, I79).

**Status:** standing v_c 7.25 / p_tail −1; every sensitivity-ring
perturbation at the blessed point degrades W₁ (+0.18..+0.33) — the
standing cell sits at the basin optimum (I88). The candidate smooth
cap-replacement (Padé saturating cubic, one knob v_s) is **excluded by
zero-cost arithmetic** (atlas D4 Step 1): in-band cubicity forces
v_s ≳ ~15 Å/ps, which forfeits the p_tail = −1 saturation at the
production peak — the principled cap-remover remains extending TDDFT
authority to production kinematics (the collaborator ask), not a
different closed form. **Mid-band γ magnitude (band top → cap,
v ≈ 5–9 Å/ps) is the newly-identified deep-KE lever** (§6.6 MD
spot-check, NB-RQ11-12; escalated by the §6.7 item-1 lq battery to
~300 counts/bin, 2026-07-24): the lq system's softer γ in exactly this
window lifts deep bins ~+0.1–0.2 at n = 10–12 (direction confirmed at
battery scale), a region where cubic-vs-quadratic differ most and TDDFT
authority is absent (the collaborator-ask window). **CORRECTED by the
§6.7 battery:** the §6.6 read that this mid-band softening *"halves the
same-N χ²_med / lands the KE better"* **does NOT replicate paired** —
across 5 fresh seeds at N = 1000, lq's median-anchored KE χ²_med is
*worse* than cubic on every seed (paired Δ+13.3±10.0; pooled 277 vs
242). The §6.6 halving was a favourable-seed artifact (seed 20260721,
N = 500). So the mid-band warming is real (direction) but does **not**
yield a better KE landing; and the lq system additionally
**over-suppresses** (supp 0.213 vs 0.187, paired Δ+0.026 ~6σ) — a
distinct resolved fate-split influence. Histogram form-blindness (W₁,
midHot) still holds at battery scale. The cap arbitration (I72's
three-observables-three-v_c tension) should be re-read in this light at
synthesis. **Confound pending (§6.5 E_bind scan):** the lq cells carry
their co-extracted E_bind 0.048 eV (~0.07 eV shallower exit well), so
the residual supp/deep-direction attribution is provisional until the
E_bind-only scan separates the jointly-extracted pair.
**Corrected geometry (G3 Step 2 + Step 3 ring, §14.5):** the gate
basin sits at **v_c 5.5–6.0, MD-CONFIRMED** (a037/b031/d030/e154 land;
c50 shows the basin nearly reaches 5.0; c65 fails on n₁) — the
standing 7.25 is MD-measured broken at the corrected geometry (f725:
trap 0.437, n₁ 0). The shift is exactly the mid-band-softening
direction NB-RQ11-12 identified.

## 3. τ — Newton cooling clock

| | |
|---|---|
| role | E_int Newton-cooling time constant; sets the evaporative-descent clock and the in-bubble leak |
| class | Bounded (GAH25 band; re-classed from the 6.55 pin, §I.8) |

**Influence (measured):**

- Alone, cannot fix the KE↔histogram anti-correlation; **(v_c, τ)
  jointly close** what neither does alone (18-cell basin at τ ≈ 4.1;
  refined to τ 3.2 at the c1 winner, τ 3.8 at the c4 basin) (I49, I82).
- Under the gate, τ and f_int are **not independent** — only the race
  margin Δ× = t_eject − t× is physical; the gated arm is binary in the
  t×↔ejection race (cliff, re-opening at τ 16.5 × f_int 0.25) (I7, I16,
  I18).
- s_eff↔τ separate through first-shed timing (I4).

**Couplings:** v_c (joint basin); f_int/E₀ (race coordinate); the
cooling exposure K is the single integral through which KE and stripping
trade off (I47).

**Status:** standing 3.2 ps. **Corrected geometry (G3 Step 2 + Step 3
ring, §14.5):** the basin wants τ 4.8–6.4 — *toward* the sourced GAH25
6.55 — and the ring MD-confirms **both** τ values land at v5.5
(a037 @ 4.8, d030 @ 6.4); the τ ≫ 6.55 arms are not needed.

## 4. E₀ = E_int(0) — internal-energy budget

| | |
|---|---|
| role | absolute internal energy at onset [eV]; feeds the RRK gate via the Σ(21) crossing |
| class | Bounded — RQ1 band [0.2, 0.5] eV, floor ≈ 0.22 (row 16; replaced f_int as the physical variable, I25) |

**Influence (measured):**

- Sets the **fate cliff**: closed-form bare ⇔ K < K* = ln(E₀/Σ(21));
  the (E₀, K) two-parameter fate map reproduces measured columns to
  0.84 He (I34). At a congruent single-point ensemble the suppressed
  fraction is a delta step — coexistence requires ensemble heterogeneity
  (I26).
- The experimental histogram inverts into a smooth unimodal p(E₀)
  coinciding with the RQ1 band untuned (W₁ = 0.086 at the probe basis)
  (I29–I31).
- **Scale is geometry-dependent**: pinned-droplet reads wanted
  E₀ ≈ 0.38–0.41 (I39); full position geometry re-lands the optimum at
  the solvation scale 0.22–0.27 eV (I45). Standing 0.27.
- Budget-corner regime change: c1's suppressed channel closes entirely
  at E₀ = 0.23 (supp ≈ 0), inverting orderings (I82).
- The scenario *budget* (0.80 vs 2.70 eV) is bookkeeping-only in the
  delivered model — the gated map is budget-invariant in absolute E₀
  (I28).

**Couplings:** τ (race), droplet-K axis (the cliff position moves with
exposure: E*(K) = Σ(21)·e^K, I39), sampling width δ (bare-fraction
trade-off, I35/I36).

**Status:** standing 0.27 eV. **Corrected geometry (G3 Step 2 + Step 3
ring, §14.5):** MD lands at E₀ 0.30–0.37 (gated cells 0.30/0.31/
0.36/0.37 depending on (v_c, τ)) — still inside the RQ1 band; the
geometry-dependence chain gains its third point (pinned 0.38–0.41 →
standing mixture 0.22–0.27 → corrected ~0.30–0.37).

## 5. κ — ladder steepness (Form U)

| | |
|---|---|
| role | sigmoid steepness of the D₀(n > 1) ladder |
| class | Free (one of the two original Free knobs) |

**Influence:** near-dead: inverted *and* normalization-capped for the
21→14 staircase (Form-U floors D₀(21) at ≈ 0.53·D₀(1)) (I2); Σ(21) is
nearly κ-independent (~11 % over the pure-κ range, CALIBRATION_MAP row
21). Superseded in practice by the rq4graded ladder family (§8).

**Status:** standing within rq4graded; not an active lever.

## 6. picture — electronic-state selection

| | |
|---|---|
| role | selects the D₀ curve family (statistical_mixture / x2_only / cooling_relaxed) |
| class | Free (selection) |

**Influence:** magnitude-degenerate — ≤ 4 % spread from x-invariance;
its only staircase signal is gate-open timing via Σ(21), degenerate with
f_int (I3); landing band picture-robust (Wave 3). Discrimination belongs
to the terminal size distribution, i.e. Tier-2 arbitration, where the
mixture family stands.

**Status:** standing `statistical_mixture` family (c1 bundle).

## 7. f_int, f_ret, λ₀, p — budget fraction, retention, pickup

- **f_int** (Derived coordinate = E₀/E_avail): parametrization broken in
  both directions (no mechanical partition; literal coupling ~1 % would
  be inert) → replaced by absolute E₀ (I25). Under the gate it was a
  joint timing + budget knob (I15); only the race margin Δ× is physical
  — independent τ × f_int grids sample the race incidentally (I18).
  Retired as a free knob.
- **f_ret** (Bounded, prior small): pinned 0.1 through the whole probe +
  campaign chain; **never swept — GAP** (atlas candidate if ever
  motivated; no measured influence on record).
- **λ₀ pickup** (Sourced + Bounded, ~0.7–1.1/ps): live Langmuir
  re-filling is real and measurable — MD suppression sits ≈ 7 ions/100
  below the frozen twin at leg B, ≈ 1–3/100 at p = 1 (I63, I65).
- **p occupancy-cap exponent** (conditional-Free, standing p = 1):
  p = 0 → 1 de-suppresses ≈ 3× to the twin's prediction (best twin↔MD
  W₁ 0.31–0.47 of the chain) and lifts solvated n₁ KE to ≈ 1.1 eV —
  closest to the experimental n₁ on the solvated branch; touches only
  the onset, never chord dynamics (trapped p-invariant) (I65/I66).

## 8. Ladder — D₀(n) rungs, tail, Σ

| | |
|---|---|
| role | per-shell binding ladder: RRK gate, suppression criterion (E_ej > Σ(n₀)), descent bookkeeping, e_bind E_pot fold |
| class | D₀(1) Sourced (±3 cm⁻¹, IHe05); shape RQ4-graded (2.2 : 1.5 : 1.3 diagnostic); floor Sourced (bulk µ_He) |

**Influence (measured):**

- **Owns the small-n histogram shape**: n = 1 is a one-rung-wide E_ej
  window — the taper controls n₁; equal-deep-rung families revert n₁/n₂
  to the flat chord-K ratio (I42). The solvated histogram independently
  demands the RQ4 target ratios through a purely geometric forward model
  (I42).
- **rq4graded > floor1, decisively**: N = 500 verdict on every histogram
  read (W₁ 0.524/0.657 vs 0.843, ≫ 2σ; c4's KE curve also fails
  structurally) (I86). Cheap-bottom ladders strip *past* n = 1 into bare
  and *lower* n₁_solv — rq4graded's 20 meV D₀(1) barrier is what piles
  weight at n = 1 (I79).
- The ladder is **KE-neutral** (midHot flat across the family at every
  v_c) — orthogonal to the drag axis (I79).
- **The deep tail cannot move deep-bin KE**: deep-KE and histogram are
  rigidly Σ-coupled, ~×20 too stiff — ±5 % deep rungs ↔ ∓3.5 pts
  suppressed weight; ×1.05 already bleeds supp 0.202 → 0.167, ×1.2
  wrecks midHot (1.46); Σ-preserving transfers blow W₁ to 1.1–2.0.
  Corollary: the seed-robust histogram landing is itself evidence the
  deep rungs are approximately right (I93, §4dd).

**Couplings:** Σ(21) sets the gate threshold (κ-insensitive); the
suppression criterion and KE re-mapping read the *same* tail (the I93
stiffness); picture selects the curve family (§6).

**Status:** standing rq4graded; RQ4 external arbitration remains the
blocking authority for the taper's *physics* (I42).

## 9. E_bind — ion–droplet surface well (0.1168 eV)

| | |
|---|---|
| role | depth of the mean-field droplet exit barrier for the ion (`binding_energy_I_ion_eV`); decides eject vs trapped |
| class | Derived — jointly extracted with the drag coefficients (Tier-0 Method-B); §6.5.1 guard enforces the exact pairing |

**Influence (measured — DIRECT OAT sweep; §6.7 item-2 E_bind scan,
2026-07-24; GAP CLOSED).** Swept on the lq system at N = 1000 × 3 paired
seeds (well overridden, drag stamp held → honest under
`allow_unvalidated_binding_pairing`). Δ per the 0.048 → 0.1168 eV step,
mean ± SD:

- **trap +0.058 ± 0.004 — the clean well lever** (deeper well retains
  more; single-seed OAT trap 0.038 / 0.075 / 0.100 across
  0.048 / 0.1168 / 0.154 eV).
- **n̄ −0.50 ± 0.05, midHot −0.076 ± 0.003** — deeper well → smaller,
  colder detected clusters (the large/slow ions are retained away).
- supp +0.010 ± 0.002 (small); W₁ +0.023 ± 0.012 (marginal);
  χ²_med not resolved (±35).
- **Deep well 0.154 over-retains** — 2/3 N = 1000 cells trip the P1–P3
  detection handover guard (~2/2000 ions never decouple from He in
  8000 ps; ρ̂ ≈ 0.93, live drag exposure).
- **Not the driver of lq's over-suppression** — that is the drag FORM:
  the matched-well form gap Δsupp +0.034 (~8σ) survives, the well adds
  only +0.010 (see §1). This scan's purpose was to separate that pair.

Prior indirect reads (retained):

- Defines the **trapped droplet-retained class**: 5–11 % of fragments
  (inward partners of off-center births, chord K up to ≈ 18) dissipated
  below the 0.117 eV barrier, never ejecting; includes centrifugal
  resonances; fragile at meV scale (1–3 ions/config flip at a
  convention change) (I44, I56, I61).
- The §4ee low-v twin sweep moved this boundary *from the drag side*:
  un-trapping the marginal class (trap 0.042 → 0) was the floors'
  dominant effect — an in-model echo of RQ11 candidate (iv).

**Couplings:** drag law (jointly calibrated pair — "correct drag traps
the ions", TIER0_FINDINGS); *not* connected to the ladder D₀(n)/Σ
(different bookkeeping surfaces: droplet exit well vs per-atom shell
energetics). **Corrected geometry (G3 Step 2, §14.5):** the trap lever
steepens dramatically with v_c at long chords (eb0482 holds trap
0.01–0.04 everywhere; eb154 reaches 0.55–0.61 at v_c ≥ 8 with K_q50
17–18) — the (v_c, well) coupling is the §3.5b over-dissipation
distinguishing evidence; the gate itself is only weakly well-sensitive
(gated cells exist at all three wells, v_c window shifting ~0.5).
**Ring MD (arm E, §14.5):** trap ordering CONFIRMED in MD at the basin
chord — 0.026 / 0.087 / 0.115 for 0.0482 / 0.1168 / 0.154 eV (≈ 0.8/eV,
matching the item-2 lq-system 0.85/eV — the lever transfers), and the
deep well *gates* (e154), so the basin tolerates the full Tier-0 well
spread.

### 9.1 Is E_bind R-dependent? — bounded analytically (2026-07-26, zero cost)

The well was co-extracted in an ⟨N⟩ = 2000 droplet, so Axis A's large-R
cells would run it outside its calibration geometry. The bound splits the
well by range:

- **Local term (snowball / electrostriction)** — the compressed first
  shells within a few Å of the ion. Nonlinear, dominates the absolute
  depth, and **R-independent**: it does not know where the surface is.
- **Far-field term** — the dielectric continuum out to the surface. The
  only part that carries R, and weak because helium is nearly
  non-polarizable.

Born for an ion of radius a at the center of a dielectric sphere of
radius R, `W(R) = −c·(1/a − 1/R)` with
`c = (q²/8πε₀)(1 − 1/ε) = 7.20 eV·Å × 0.0541 = 0.390 eV·Å`
(liquid He ε = 1.0572). The R-term is `+c/R` — a smaller droplet is a
*shallower* well:

| R [Å] | c/R [eV] | ΔE_bind vs R1 | implied Δtrap (at 0.85/eV) |
|---|---|---|---|
| 26.6 (R1, standing) | 0.0146 | — | — |
| 34.0 (R2) | 0.0115 | +0.0032 | +0.003 |
| 49.4 (R3) | 0.0079 | +0.0068 | +0.006 |
| 68.3 (largest cell) | 0.0057 | +0.0089 | +0.008 |

The Δtrap column uses this section's own measured lever (trap
+0.058 ± 0.004 per +0.0686 eV ⇒ ≈ 0.85 per eV).

- **Deepening over the entire grid span ≤ 0.009 eV = 7.6 % of 0.1168** —
  about ⅛ of the smallest step measured to matter, moving trapped by
  ≈ 0.006 against the 0.0003 → 0.158 R-selection effect the grid
  measures. One-sided and small.
- **Free sanity check on the extraction:** at a = 3.5 Å the same formula
  gives W(26.6 Å) = 0.096 eV vs the jointly-extracted 0.1168 eV — within
  ~20 %, the first independent number of any kind on E_bind.
- **Caveat:** Born underestimates the *absolute* depth (no
  electrostriction), so this validates nothing about 0.1168 eV. It bounds
  the **difference** only — legitimately, because the local term cancels
  there.

**G0 decision (2026-07-26):** defer — all Axis A cells run at 0.1168 eV;
the R-dependence is carried as this bound. An MD bracket cell at R3 is
**armed as a rider**, firing only if G1 shows the trapped channel near a
boundary (trapped > ~0.4, or a detection-handover-guard trip as at
0.154 eV). Modelling E_bind(R) as a derived quantity is explicitly not
done in this program.

## 10. Landau v_L — dissipation threshold

| | |
|---|---|
| role | γ = 0 gate below v_L (0.58 Å/ps); also a latent shared reader (sets the neutral-stage E_min threshold, I84) |
| class | Sourced — bulk roton Landau value, droplet-persistent down to ~1000 atoms (NB-RQ11-1..3) |

**Influence (measured):** **bit-flat on every scored observable** across
v_L ∈ {0.30, 0.40, 0.58} (Δn₁_solv = ΔW₁ = 0.0000) while demonstrably
acting on the retained class (retained KE 0.277 → 0.005–0.023 eV,
monotone) — the scored surface is structurally insulated, not the arm
inert (I84, re-confirmed at the blessed point I88). The Landau floor
(≈ 0.003 eV) sits an order of magnitude below the scored deep bins —
part of why RQ11 candidate (i) fails (I89).

**Status:** standing 0.58, quiet on the scored surface; literature
deep-read confirms the value but shows the *gate shape* above v_L is
unconstrained (I100).

## 11. Per-shed ε — evaporative recoil

| | |
|---|---|
| role | translational energy carried per shed He (currently ε = 0; OQ-F origin, I23) |
| class | Bounded small (NB-RQ23-1: ~0.5–0.6 meV prescription) |

**Influence (measured):** required ε_req(n) to close RQ11 is
sign-changing (−3…−4 meV/shed at n = 2–3, +3 → +9 over n = 10–17); flat
ε ≈ 6–8 meV closes the deep slope (χ²_med 242 → 66) only by breaking
midHot (1.25–1.33 vs [0.84, 1.14]); physical drag-retention is
*inverted* (deep-bin sheds happen in-bubble, 43 % retention vs 95 % at
n₁) — retention-weighted forward model floors at χ²_med ≈ 198.
**Refuted as RQ11 owner; ε ≈ 1–2 meV sub-dominant stays compatible**
(I94, §4dd).

## 12. Shed momentum convention (cold-shed vs co-moving)

| | |
|---|---|
| role | what the evaporated He carries away: at-rest-in-lab (cold-shed, momentum-conserving) vs co-moving (ε ≈ 0) |
| class | convention — RQ2/RQ3-coupled (OQ-J) |

**Influence (measured):** cold-shed injects +0.91 eV over a full
post-exit strip (×1.611 = m₂₁/m₁ in KE) — it *was* the twin−MD n₁ KE
divergence, closed entirely by the co-moving counterfactual (I59, I61).
The histogram is convention-blind (W₁(A″, A′) = 0.14–0.20); the
suppressed-class *fragmentation* read is convention-decided: pooled
conventions give intact 1.96 / co-moving 1.36 / cold-shed 2.82 eV vs
ref bare mean 3.706 — direction momentum-conserving, all 6–10× too
narrow in σ, and the bare mean is kinematically un-sourceable from the
production channel (cap ≈ 3.37 eV pre-drag) → n = 0 is a two-channel
mixture; the suppressed class belongs to the low-KE shoulder (I60, I64,
I95).

## 13. Birth margin — 3 Å pinned convention

| | |
|---|---|
| role | minimum birth depth below the droplet surface for the sampled I₂ center |
| class | pinned convention (I88) |

**Influence (measured):** **the sensitive robustness lever** at the
blessed point: margin 3 → 4.67 / 6.0 Å moves n₁_solv −0.094 / −0.171
(≈ 1.9σ / 3.4σ), collapses the ratio, halves suppression — the
histogram landing depends on this pin (I88). In-tier open item.

**Couplings:** position sampling (§15) — the margin truncates the
Boltzmann radial law near the surface; the position axis is a two-sided
race lever (I55).

## 14. Droplet geometry — size R × birth position r **(GAP CLOSED — all 11 grid cells scored; Axis A G1 2026-07-26/27 + the retained-class arm 2026-07-27)**

> **Caveat that travels with the R ≥ 49 Å rows only.** Those cells exclude a
> `droplet_retained_marginal` class (ions energetically able to escape but
> still helium-coupled at handover) whose size is conditional on a cubic drag
> law extrapolated ~3× beyond its 9/18 Å calibration band. §14.2 measures
> what that exclusion does to each observable: it is **resolved, not
> negligible**, on `deepKE`, `n̄` and `W₁`. The R ≤ 34 Å rows have an
> **empty** marginal class and are unaffected.

Never controlled independently; everything known is inference from
sampled ensembles:

- The droplet-radius axis has a **quantified demand**: ≥ 9.6 % of the
  ensemble below the 0.22 eV solvation floor requires a larger-K
  (larger-droplet) tail (I30); width attribution: the ~21-bin span must
  be carried by the E₀ and/or droplet-K axes — per-point reads are
  ≤ 1 bin wide (I33).
- The **position axis alone opens the two-sided race** (first MD weight
  at n = 0/1 with no dressing; suppressed/bare ordering transfers
  exactly from the twin) (I55); center-pinning was half of the T3 park
  (I52).
- The **droplet prior and the KE axis are decoupled**: flipping the
  kornilov prior moves the histogram (n̄ ≈ −1 He) while solvated n₁ KE
  moves ≤ 0.12 eV (I68); prior axis asymmetric at the blessed point
  (δ 0.40 / pickup ≈ 1.5σ; δ 0.80 noise-level) (I88).
- **Fate is birth-dressing + droplet-size ordered** (deep survivors:
  larger droplets 30.1 Å, fuller dressing 18.4; retained: 34.0 Å)
  (I90); frozen-geometry twin bias (channel (d)) scales with in-droplet
  residence, 0 → −0.8 He across the v_c bracket (I73), (τ, E₀)-blind
  (I74).

**Measured influence — entangled battery read (2026-07-26, zero MD;
full tables in `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` "Axis A pre-read").**
The pooled N = 5000 battery binned by its own sampled geometry. Sampled,
not controlled — R and depth co-vary and the T5 `density_tied` dressing
rides on depth (⟨n₀⟩ = 14.0 / 16.4 / 19.3 across depth terciles) — so
this orders and bounds, it does not attribute.

- **Birth depth is the physics knob** (terciles 4.4 / 7.7 / 14.6 Å):
  supp 0.381 → 0.182 → 0.003, n̄ 2.36 → 3.52 → 6.29, n₁_solv
  0.358 → 0.109, midHot 0.703 → 1.209, **deep-KE 0.316 → 0.792**.
- **Droplet radius is a selection knob**: at *fixed* depth, R is nearly
  inert on the scored ensemble (supp and deep-KE flat to ±0.03 across R
  terciles within each depth tercile); what R does move is trapped,
  0.0003 → 0.042 → 0.158. Mechanism: the escaping fragment's path is its
  birth depth (R-independent at fixed depth); R sets the *inward*
  partner's path ≈ 2R − depth, and those fragments leave the ensemble.
- **The landing is a mixture property**: pooled W₁ 0.571 beats every
  geometry bin (best tercile 0.750, best quintile 0.639). No single
  geometry lands the histogram.
- **RQ11 relevance**: the deep-bin cold tail closes by ×2.5 with birth
  depth — the first candidate owner this program has found. Confounded
  with the dressing. **The control that would decompose it does not
  exist** (2026-07-26): the three density-keyed legs — T5 birth dressing,
  ρ̂-scaled cooling, drag gate — all read the single shared surface
  `rho_he_ratio(depth, …)`, so they are *one* physical fact (local He
  density at birth) expressed three ways, and every proposed control
  (`initial_shell_model="full"`, `internal_energy_partition_law=
  "constant"`, `cooling_spatial_gate="none"`) requires the less physical
  arm of its pair. The correction achieves what they were meant to
  measure: at the anchored geometry ρ̂ saturates and T5/T6 go inert on
  their own.

### 14.1 Measured influence — the CONTROLLED grid (Axis A G1, all 11 cells)

**All 11 cells scored.** Six landed on the first pass (2026-07-26); the five
anchored-radius cells were completed on 2026-07-27 by the retained-class arm
(`detection_droplet_retained_policy="exclude_all_coupled"`) — a
**detection-only re-run from the stored `relaxation.npz`, zero new MD**.
N = 500, **one shared seed 20260727** (common random numbers), fixed droplet
size per cell, everything else at the standing point; committed scorer,
pooled-battery oracle reproduced (all 7 columns within 0.002). Re-running
detection on the six original cells under the new arm reproduced them
**bit-for-bit** (the arm oracle — the marginal class is empty at R ≤ 34 Å).
Per the §3.1 mixture rule these cells are read against **each other**, never
against the acceptance.

`trap` is decomposed: **`t_b`** = `droplet_retained` (**physics** — total
energy below the effective-potential barrier including angular momentum, so
the ion cannot escape at any relaxation length) and **`t_m`** =
`droplet_retained_marginal` (**modelling exclusion** — can escape, still
helium-coupled at handover). They are never summed in a reading; `trap` is
their total.

| cell | R [Å] | law | ⟨depth⟩ | trap | t_b | t_m | supp | n̄_det | n₁_solv | W₁ | midHot (bins) | deepKE (bins) | χ²_med (npts) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r1l1 | 26.6 | center-pin | 26.6 | 0 | 0 | 0 | 0 | 8.39 | 0 | 4.98 | 1.04 (3) | 1.83 (1) | 406 (4) |
| r1l2 | 26.6 | parent Boltzmann | 24.9 | 0 | 0 | 0 | 0 | 8.35 | 0 | 4.76 | 1.25 (4) | 1.58 (2) | 117 (6) |
| r1l3 | 26.6 | `uniform_volume` m3 | 9.0 | 0.059 | 0.059 | 0 | 0.167 | 4.46 | 0.235 | 0.81 | 1.02 (7) | 0.51 (8) | 266 (17) |
| r2l1 | 34.0 | center-pin | 34.0 | 0 | 0 | 0 | 0 | 12.58 | 0 | 7.96 | NaN (0) | 1.49 (6) | 339 (6) |
| r2l2 | 34.0 | parent Boltzmann | 29.3 | 0 | 0 | 0 | 0 | 12.30 | 0 | 7.53 | 1.42 (2) | 1.38 (8) | 46 (10) |
| r2l3 | 34.0 | `uniform_volume` m3 | 10.9 | 0.196 | 0.196 | 0 | 0.140 | 5.33 | 0.184 | 1.40 | 1.11 (7) | 0.75 (8) | 127 (17) |
| r3l1 | 49.4 | center-pin | 49.4 | 0.538 | 0.324 | **0.214** | 0 | 18.35 | 0 | 13.46 | NaN (0) | 1.94 (2) | 39 (2) |
| r3l2 | 49.4 | parent Boltzmann | 34.8 | 0.476 | 0.470 | 0.006 | 0 | 15.12 | 0 | 10.24 | NaN (0) | 2.19 (8) | 311 (8) |
| r3l3 | 49.4 | `uniform_volume` m3 | 14.8 | 0.404 | 0.402 | 0.002 | 0.123 | 6.30 | 0.151 | 2.31 | 1.31 (7) | 1.59 (8) | 109 (17) |
| r4l2 | 68.3 | parent Boltzmann | 40.5 | 0.800 | 0.200 | **0.600** | 0 | 15.28 | 0 | 10.39 | NaN (0) | 1.93 (8) | 173 (8) |
| r4l3 | 68.3 | `uniform_volume` m3 | 19.5 | 0.570 | 0.396 | **0.174** | 0.126 | 6.71 | 0.128 | 2.83 | 1.33 (7) | 1.97 (8) | 363 (17) |
| *pooled battery* | *sampled* | *m3* | *9.0* | *0.067* | *0.067* | *0* | *0.187* | *4.07* | *0.243* | *0.571* | *1.011* | *0.631* | *242 (17)* |

**Read χ²_med only with `ke_npts`.** The deep-birth cells vacate the low-n
bins, so their χ² sums over 4–10 points against the pooled 17; r2l2's 46 and
r1l2's 117 are fewer terms, **not** better KE agreement. Same for midHot /
deepKE at L1/L2 — r1l1's deepKE rests on one bin, r2l1's midHot band is
empty (NaN, correctly reported as 0 bins).

1. **The birth law is the dominant axis; R at fixed law is second-order.**
   At R1, going from the standing shallow law to either anchored deep law:
   supp 0.167 → **0**, n₁_solv 0.235 → **0**, n̄ 4.46 → 8.39, W₁ 0.81 → 4.98.
   Under L3, R 26.6 → 34 Å moves n̄ 4.46 → 5.33, trap 0.059 → 0.196,
   deepKE 0.51 → 0.75 (with mean depth also rising 9.0 → 10.9 Å, since
   `uniform_volume` depth scales with R — the two are not separable on this
   law). Confirms the pre-read's "depth is the physics knob, R the selection
   knob" with controlled cells.
2. **Width decomposition — the headline (§3.3 Q5 answered).** The detected
   size distribution's width tracks the *birth-depth* spread almost one-for-one,
   on top of a small irreducible mechanism term:

   | cell | SD(birth depth) [Å] | SD(n_det) [He] |
   |---|---|---|
   | r1l1 / r2l1 (center-pin) | **0** | 0.63 / 0.98 |
   | r1l2 / r2l2 (Boltzmann) | 0.86 / 1.78 | 0.93 / 1.74 |
   | r1l3 / r2l3 (m3) | 4.48 / 5.88 | 4.34 / 4.88 |

   So SD(n) ≈ **1 He per Å of birth-depth spread**, with the zero-spread
   cells exposing the mechanism-only floor (Poisson pickup + RRK evaporation,
   0.6–1.0 He). At the standing law that makes **≈ 95 % of the size variance
   geometry-inherited** (0.63²/4.34² ≈ 2 % mechanism). The L1 cells earn
   their place precisely here — the near-degeneracy with L2 in the *means*
   (n̄ 8.39 vs 8.35) is what makes them a clean zero-spread reference.
   No ceiling pile-up anywhere (share at n = 21 is 0.000 in all six cells),
   so the n̄ climb is not saturation.
3. **The landing needs the size *distribution*, not the right mean size.**
   r1l3 is the standing geometry with only the droplet-size spread removed
   (R fixed at its own realized mean): W₁ degrades **0.571 → 0.813** while
   supp / n₁_solv / trap / n̄ stay close. A controlled sharpening of the
   pre-read's mixture finding.
4. **RQ11: the deep-bin deficit is geometry-reachable and it inverts.**
   deepKE 0.51 (depth 9 Å) → 0.75 (10.9) → 1.38 / 1.49 (29–34 Å), the
   deep-birth values now on 6–8 occupied bins rather than 1–2. So the cold
   tail does not merely close with depth — it crosses 1.0 and overshoots to
   ≈ +40 % too hot, implying some intermediate depth matches the deep-bin KE
   exactly. First controlled evidence that RQ11 sits on an axis geometry can
   traverse. **Not** a claim that these cells are better models — their
   histograms are destroyed (W₁ 7.5–8.0).
5. **Trapping requires a shallow, off-centre birth, then scales with R.**
   trap is exactly 0 in all four L1/L2 cells and 0.059 → 0.196 across L3.
   The long path belongs to the *inward* Coulomb partner (≈ 2R − depth);
   at L1/L2 depth ≈ R so both fragments travel only R.
6. **The frozen anchored-cell prediction is confirmed on every clause**
   (supp ≈ 0 → exactly 0; n̄ > 6.3 → 8.4; n₁_solv ≲ 0.1 → exactly 0;
   W₁ ≳ 1.6 → 4.8–5.0). The anchored geometry does not shift the histogram,
   it **evacuates its low-n half**.
7. **Along the production birth law, every observable degrades monotonically
   with R** — the completed rows make this the axis's cleanest single trend
   (L3 column, R 26.6 → 34.0 → 49.4 → 68.3 Å, mean depth 9.0 → 10.9 → 14.8
   → 19.5 Å):

   | R [Å] | trap | supp | n̄_det | n₁_solv | W₁ | midHot | deepKE |
   |---|---|---|---|---|---|---|---|
   | 26.6 | 0.059 | 0.167 | 4.46 | 0.235 | 0.81 | 1.02 | 0.51 |
   | 34.0 | 0.196 | 0.140 | 5.33 | 0.184 | 1.40 | 1.11 | 0.75 |
   | 49.4 | 0.404 | 0.123 | 6.30 | 0.151 | 2.31 | 1.31 | 1.59 |
   | 68.3 | 0.570 | 0.126 | 6.71 | 0.128 | 2.83 | 1.33 | 1.97 |

   Trapping and W₁ rise, n₁_solv falls, and the KE observables heat
   monotonically. Note `uniform_volume` depth scales with R, so this column
   moves both axes together by construction — it is the *production law's* R
   response, not an R-at-fixed-depth read (which the pre-read showed is
   nearly inert).
8. **The deepKE = 1 crossing is located, on the production law.** Along L3,
   deepKE crosses 1.0 between R 34.0 and 49.4 Å, i.e. **mean birth depth
   ≈ 11–15 Å**, on 8 occupied bins throughout — no longer the 1–2-bin read
   finding 4 had to caveat. This sharpens finding 4: RQ11's deficit and its
   overshoot bracket a depth in that window. **The R ≥ 49 Å deepKE values
   are, however, bracket-sensitive** (§14.2): at r4l3 the marginal exclusion
   alone accounts for +0.42 of the 1.97, and at r3l1 for +1.61 of 1.94. The
   *crossing* survives because r3l3/r4l3 carry marginal fractions of 0.002 /
   0.174 and the r3l3 bracket is tight, but any quantitative deep-KE claim at
   the anchored radii must be quoted with its bracket.

### 14.2 The retained class at the anchored radii — decomposed, and the exclusion bracketed

All five R ≥ 49.4 Å cells originally failed the detection-stage P1–P3 handover
guard. Not a bug and not a tolerance artifact — the measured mechanism:

| cell | R | still helium-coupled at 8 ns | **bound** (`droplet_retained`) | **marginal** (`…_marginal`) |
|---|---|---|---|---|
| r3l1 | 49.4 | 538/1000 | 324 | 214 |
| r3l2 | 49.4 | 476 | 470 | 6 |
| r3l3 | 49.4 | 404 | 402 | 2 |
| r4l2 | 68.3 | 800 | 200 | 600 |
| r4l3 | 68.3 | 570 | 396 | 174 |

**The composition flips with R, which is why the two classes are never
summed.** At R = 49.4 Å the coupled population is essentially all *bound* —
98–99 % of it at r3l2/r3l3, a physics statement. At R = 68.3 Å on the parent
law it is 75 % *marginal* (600 of 800), i.e. three quarters of that cell's
exclusion is a modelling convention. A single `trap` column would have read
0.476 and 0.800 as the same kind of number.

- **A third to a half of all ions never leave the droplet** at the parent's
  radii — plan §3.6's "does the trapped channel explode at large R?" answered
  emphatically yes.
- Those still inside sit 20–30 Å **below** the surface with |v_rad| ≈ 0.01–0.04
  Å/ps. Cubic drag over a 25–50 Å path takes them below the Landau threshold
  (v_L = 0.58 Å/ps), where `landau_gated_drag` switches dissipation off; they
  then **oscillate** conservatively in the well. Median net radial progress
  over the window's last 4 ns is **−0.5 Å** (r3l1) / **−0.9 Å** (r4l2);
  only ~35 % still move outward, ~18 % have flatly asymptoted.
- The escaping subset needs **~80–690 ns** more at its observed late drift
  rate — **10–90×** the 8 ns relaxation window, not the ~5× a naive read
  suggests.
- **The timescale separation the staging assumes does not exist at the
  corrected geometry.** 30 ps drag-active MD → 8 ns conservative E2 →
  8.53 µs free flight is calibrated to R ≈ 27 Å, where ejected ions leave at
  several Å/ps and cover ~48 000 Å (≈ 5 µm) inside the E2 window — decoupled
  by an enormous margin. At R ≥ 49 Å escape time and flight time become
  comparable, with helium still present. Recorded as a ledger row (§17): a
  structural consequence the correction **creates**, the first found.

**Resolution (2026-07-27): count the coupled class instead of integrating it.**
`detection_droplet_retained_policy="exclude_all_coupled"` excludes every
still-coupled ion and counts it, keeping bound/marginal decomposed. Adopted
after the ~µs staging extension was rejected on the grounds that it would
spend ~100× the relaxation cost integrating trajectories that are an artifact
of a law extrapolated past its calibration band. Recovery cost **zero MD** —
the cells had failed at detection, after the MD.

**The exclusion is NOT observable-neutral — measured, per the pre-registered
§3.5b item-8 bracket.** Each cell was scored twice: **Arm A** with marginals
excluded (the headline row above) and **Arm B** with them injected at their
handover `n` and their exact conservative asymptotic KE. Arm B costs nothing
new: E2 is zero-gamma, so an ion clearing its barrier arrives at infinity with
exactly `E_tot − U(∞)`, a quantity the classifier already computes. Seed-SDs
measured from the five N = 1000 battery members and scaled by
√(N_member/N_cell):

| cell | t_m | Δn̄ | Δn₁_solv | ΔmidHot | ΔdeepKE | Δχ²_med | ΔW₁ | verdict |
|---|---|---|---|---|---|---|---|---|
| r3l1 | 0.214 | −0.44 | 0 | — | **−1.61** | +128 | −0.44 | **WIDE** |
| r3l2 | 0.006 | +0.01 | 0 | — | −0.02 | −23 | +0.01 | tight |
| r3l3 | 0.002 | +0.04 | −0.001 | 0.000 | −0.010 | −2.6 | +0.04 | tight |
| r4l2 | 0.600 | **+2.50** | 0 | — | **−0.41** | −38 | **+2.50** | **WIDE** |
| r4l3 | 0.174 | **+3.38** | **−0.040** | 0.000 | **−0.42** | −123 | **+3.35** | **WIDE** |

Direction exactly as pre-registered — deepKE **down**, n̄ **up**, midHot
**flat**, supp unchanged — and the magnitude is 10–38 seed-SD on deepKE where
the class is large. **So the R ≥ 49 Å deepKE and n̄ values are conditional on
the exclusion, and the R ≤ 34 Å rows are not** (empty marginal class,
bit-for-bit reproduced).

**The dossier is the physical finding, and it settles what the class is.**
Marginal ions' asymptotic energetics, read off the stored handover states:

| cell | n_marg | median handover n | median depth [Å] | median KE_asym [eV] | KE_ceiling [eV] | ref ⟨KE⟩ at that n [eV] |
|---|---|---|---|---|---|---|
| r3l1 | 214 | 17 | −12.9 | **0.000** | 0.094 | 0.066 |
| r3l2 | 6 | 16 | −12.2 | **0.000** | 0.093 | 0.069 |
| r3l3 | 2 | 17 | −11.8 | **0.000** | 0.091 | 0.067 |
| r4l2 | 600 | 19 | −15.9 | **0.000** | 0.031 | 0.066 |
| r4l3 | 174 | 19 | −15.5 | **0.000** | 0.031 | 0.066 |

Two readings, both important:

1. **They would arrive essentially at rest, carrying a large shell.** Median
   handover n = 16–19 with asymptotic KE clipping to zero. That is precisely
   the deep-n / low-KE corner, which is why injecting them craters deepKE
   (1.94 → 0.33 at r3l1). The intuition that the retained class is "slow,
   high-n ions" is **confirmed quantitatively**; the accompanying intuition
   that it therefore "doesn't shift the result" is **refuted** — slow and
   high-n is not a null direction, it is the RQ11 direction.
2. **Most of the class is probably bound, not marginal.** `KE_asym` uses the
   *physical* half-credit pair-Coulomb split; `_conservatively_bound` credits
   the **full** pair energy to both fragments (a deliberate, documented
   over-estimate that keeps a `bound` verdict certain). The median marginal
   ion has `E_tot − ½E_coul − E_bind ≤ 0`, i.e. it clears its barrier *only*
   on the over-crediting convention. "Marginal" therefore means **"not
   provably bound"**, and the marginal fraction is an **upper bound** on
   genuine slow escapers. The full-credit ceiling (0.03–0.09 eV) sits at or
   below the experimental mean KE at the same n (0.066 eV) — these fragments
   would be a cold population the reference does not show. That is not a
   detector-acceptance test (the reference is a first-moment table, not an
   acceptance function), but it is the direction that matters.

**What the bracket does not close.** Arm B holds `n` at its handover value and
lets the ion coast out unchanged. Over the ~80–690 ns these ions actually
need, pickup is live (λ₀ρ) and would both raise `n` and mass-load them —
plausibly into the bound class, which would mean the conservative split
*understates* trapping. That is the kinetic-forward-model follow-up
(conservative orbit + Poisson pickup + RRK), still open.

**As-built geometry the standing point actually sits at** (audit
2026-07-26, §15): `kornilov_lognormal` prior pinned at ⟨N⟩ = 2000,
δ = 0.625, truncated [250, 16000] He; births `uniform_volume` with a
3 Å margin. Realized in `bigc1v725s1` (N = 1000): **R q10/med/mean/q90
= 20.1 / 26.1 / 26.6 / 33.5 Å**; closed-form prior quantiles
N = 742 / 1647 / 3665 He → R = 20.1 / 26.2 / 34.2 Å (truncation
discards 0.14 %). These are the numbers the §3 grid should be built
from — no sampler draw needed.

**Structural note (why the two axes may not separate).** R and r enter
the physics *only* through `depth = r − R`: the erf solvation
potential/force, the erf-complement drag gate `g(depth)`
(`physics/drag.py`), the `density_tied` birth dressing, and the
escape/handover test. He density is a constant. To first order both
axes therefore act through the single scalar chord `R − r₀`, so the
3 × 3 factorial is expected to be **near-degenerate along the exposure
diagonal**; the discriminating read is whether the n = 1 velocity
distribution breaks that degeneracy (plan §3.3 Q2). Worth
pre-registering before the grid runs.

**Atlas target — EXECUTED, all 11 cells (see §14.1):** the §3 factorial grid,
adjudicated **option (B)**
2026-07-26 — R ∈ {26.6, 34.0, 49.4} Å × birth law ∈ {center-pin,
parent-consistent Boltzmann, `uniform_volume` m 3}: a
standing-vs-anchored-geometry test spanning birth depths 8.9 → 49.4 Å
(the isotropic-chord coordinate is *withdrawn* — the pre-read showed R is
inert at fixed depth). **Funded at 11 cells** (2026-07-26): the 3 × 3
plus two R = 68.3 Å cells on L2/L3, the parent's largest droplet, which
also close the D2b re-weighting extrapolation. Plus the D2b re-weighting
(the ⟨N⟩ = 12794 variant sits inside the grid support).

**Design caution carried into the grid** (from §9's own scan, finding 4):
at N = 500 with a single seed, W₁ swings ≈ 0.1 on seed alone. Per-cell W₁
is therefore *reported, not read as a discriminator* below |Δ| ≈ 0.15;
the grid's discriminating observables are supp, trap, n̄, n₁_solv, midHot
and deep-bin KE, which resolved at 8–35σ in that scan. All cells share
one seed (common random numbers).

**Near-degeneracy noted, kept deliberately:** L1 (center-pin) and L2
(parent Boltzmann) have almost the same mean depth and chord at every R
(26.6/26.6, 34.0/33.8, 49.4/47.9 Å) and both saturate the dressing, so
they differ essentially only in birth-depth *spread*. That contrast is
the axis's handle on how much ensemble spread is geometry-inherited
(§3.3 Q5); all three L1 cells were kept for a balanced factorial.

### 14.3 Why birth depth is the lever — two-channel decomposition (synthesis, 2026-07-27 discussion; measured components §4p + G1)

The birth-depth influence is two separable geometry-derived channels, both
already measured, and the distinction decides which one downstream knobs
can compensate:

1. **Partial initial dressing** (T5 `density_tied`; probe findings **§4p**,
   the one-lever A/B at the standing geometry): flipping n₀ = 21-for-all to
   depth-dressed (mean n₀ ≈ 17.2 at 9 Å births) **triples the
   suppressed/bare class** (0.096 → 0.322 at c1; ordering transfers from
   the twin), drops n̄_det by 1.0–1.9 He, truncates the deep tail at the
   dressed band, and pulls n₁ KE from ≈ 1.0 to ≈ 0.63 eV — W₁(lever) =
   1.0–1.9 bins, one of the largest single levers of the whole build, and
   the best twin↔MD agreement of the oracle chain (W₁ 0.41–0.51). *(A
   session recollection that this activation was low-impact was checked
   against §4p and corrected.)* **Geometry-locked:** the dressing
   saturates by depth ≈ 25 Å (§G) and shares the one `rho_he_ratio`
   surface with the cooling and drag gates (§3.1c rejection) — no
   downstream knob can un-saturate it at realistic depths.
2. **Short transit** (exposure): shallow births exit fast, so pickup
   re-filling is small (measured ≈ 7 ions/100 at shallow chords, §4p) and
   the cascade freezes early; deep births transit long at ρ̂ ≈ 1 and
   re-fill toward n ≈ 20 — the origin of the ≈ 1 He/Å depth slope
   (§14.1). **Parameter-accessible:** transit time, cascade duration,
   budget and exit boundary are exactly what (v_c, τ, E₀, E_bind) touch.

Consequence for the G3 question ("can the realistic geometry land with
the right parameters"): channel 1 is dead at realistic depths regardless
of parameters; channel 2 is open. Two candidate routes, both twin-scannable:
**Route A (cascade)** — shed ~10 more He/ion (≈ 0.05–0.15 eV extra E_int,
same order as E₀) and cool deep arrivals (the deepKE sign flip and its
crossing at 11–15 Å mean the arbitration passes through the right value);
the small in-tier lever slopes were all measured at *short* transits and do
not bound long-transit behavior. **Route B (selection)** — make only
shallow-born ions detectable: the trap/retention boundary (E_bind lever
0.85/eV measured, §9; v_c) plus the marginal class sitting at/below the
experimental KE support (§14.2) can carve the detected ensemble toward the
shallow tail — the emergent version of what the ⟨N⟩ = 2000 prior was
hand-doing. Acceptance note: supp is a *softer* target than n₁_solv/n̄/W₁
(the experimental n = 0 bin is an RQ3 two-channel mixture; the solvated
histogram renormalizes to n ≥ 1), so losing the suppression channel at
saturated dressing is not by itself disqualifying.

### 14.4 Twin authority at the corrected geometry (G3 Step 1, 2026-07-27 — the scan instrument's error model; **ring-validated at Step 3**: n̄ bias 0.24–2.66 He small-end residence-scaled at 14/14 cells, n₁ transfer ≤ 0.036, trap floor +0.012…+0.13)

The S6 twin was re-issued at the corrected geometry and measured against
all 11 G1 MD cells (stage `g3landmarks`; full record findings "G3
Step 1"). The per-observable authority box, superseding the
standing-geometry channel-(d) numbers for G3 use: **supp / n₁_solv
near-quantitative** (within 0.02 / 0.016 on L3, exact 0 on L1/L2);
**n̄ hot with a residence-scaled bias** (+0.2…+3.2 He across R1→R4, W₁
inheriting it) — twin scans aim at target + bias; **trap a twin floor**
(−0.03…−0.14; the r3l1 center-pin trap channel, MD 0.538, is
structurally twin-invisible — mechanism-made by pickup mass-loading +
Landau freeze); **KE direction-only** (+15–30 % hot); all orderings
preserved incl. the deepKE = 1 crossing bracket. The corrected-ensemble
twin row at finc1v725 parameters confirms the §15.7 forecast
cross-instrument (trap 0.308, supp 0, n̄ 17.0, W₁ 12.1, deepKE 2.12 —
each discrepancy = the measured bias). Landmark table re-issued and
law-tagged (the recorded center-pin K 0.74460 is pure-cubic; under the
capped tail the anchored radii are traversable at K_capped 1.7–20.2).

### 14.5 The corrected-geometry landing surface — (v_c, τ, E₀, E_bind) re-arbitrated by the G3 Step 2 twin scan (2026-07-27, zero MD)

The §3.5c nested Route A/B factorial: 30 chord families (v_c × E_bind,
one integration of the committed corrected master each) × 216 free cells
(τ × E₀) = 6480 scored cells against the pre-registered hard gate
(n₁_solv ∈ [0.19, 0.30] ∧ n̄ ∈ [4.4, 7.1]; W₁/KE/supp non-gating; trap a
floor). Oracles bit-exact (S6 rows + the Step-1 corrected landmark). Full
record: findings "G3 Step 2"; stage `g3scan`.

**Result — the basin exists and is Tier-0-legitimate: 24/6480 cells
gate, every one inside the calibrated v_c range** (the pre-registered
failure criterion did NOT fire; the sub-band diagnostic arms land
nothing). Four cells are clean of *both* caveat stamps (standing well
eb1168, τ ≤ the sourced 6.55): **(v_c 5.5, τ 4.8, E₀ 0.36–0.37) and
(v_c 6.0, τ 6.4, E₀ 0.32–0.33)** — twin W₁ 0.50–0.74, trap floor
0.053–0.129, supp 0.08–0.13, deepKE 0.46–0.80 (direction-only). Versus
finc1v725 the corrected geometry re-arbitrates **all three knobs in
physically comfortable directions**: v_c 7.25 → 5.5–6.0 (softer mid-band
γ — the NB-RQ11-12 deep-KE lever direction), τ 3.2 → 4.8–6.4 (*toward*
the sourced GAH25 6.55 — the calibration-class tension relaxes), E₀
0.27 → 0.32–0.37 (inside the RQ1 band; extends the §4 "scale is
geometry-dependent" chain: pinned 0.38–0.41 → standing-mixture
0.22–0.27 → corrected ~0.35).

**Mechanism reads:**

- **Cascade-carried, not selection-carried.** At the gated families
  det_yield is 0.87–0.99 and the detected-subset R quantiles sit within
  1–2 Å of the source — the Route-B droplet-size selection axis is
  barely exercised. The landing comes from the chord law itself:
  K655_q50 1.334 → 0.60–0.73 (≈ ×2 exposure reduction) *plus* the trap
  composition change (0.308 → 0.05–0.13).
- **Pure exposure rescaling cannot land it** (pre-scan): uniform
  X → f·X of the *standing* chord reaches the n₁ band at 141/216 free
  cells (Route-A kill NOT fired) but the **full joint gate at 0/216·20
  scalings** — the standing chord's K-shape + frozen trap 0.308 block
  n̄. The gate opens only when the chord *reshapes* (trap + K tail),
  i.e. the landing genuinely needs the chord surface, not a scale
  factor.
- **The standing drag point does not survive the corrected geometry:**
  the (v_c 7.25, eb1168) family gates at 0/216 — no (τ, E₀) rescue.
- **E_bind × v_c coupling (over-dissipation regime):** the trap ladder
  steepens with v_c — eb0482 keeps trap 0.01–0.04 everywhere, while at
  v_c ≥ 8 the deep well 0.154 sends K_q50 to 17–18 and trap to
  0.55–0.61. The helium-coupled class is strongly (v_c, well)-sensitive
  — evidence for the §3.5b "over-dissipation artifact" branch of the
  retained-policy sub-decision (G2 record).
- **Why the sub-band diagnostic arms fail:** at v_c 3.5/4.25 the n₁
  band is reachable only at n̄ 3.0–3.6 (< the 4.4 gate floor) —
  over-stripped. The corrected-geometry landing does *not* point below
  the TDDFT band top; the collaborator-ask branch did not fire.

**Authority caveats (§14.4 applied):** twin-level statements only — n̄
gate already carries the residence bias bracket; trap is a floor (MD
adds pickup mass-loading + the Landau freeze — the twin-invisible
center-pin channel); KE direction-only; W₁ bias-loaded. Nothing is
adopted: finc1v725 stands until G4; the designed next step is the MD
confirmation ring (~10–20 × N = 500) behind its own trigger.

**MD-CONFIRMED (G3 Step 3 ring, 2026-07-27/28 — 14 × N = 500, seed
20260728, plan §3.5d; findings "G3 Step 3").** The basin is real in MD:
**four cells land the bias-free MD acceptance** (n₁ ∈ [0.19, 0.30] ∧
n̄ ∈ [3.77, 4.37]) — **a037 (v5.5, τ4.8, E₀0.37; n̄ 4.097), b031
(v6.0, τ6.4, E₀0.31), d030 (v5.5, τ6.4, E₀0.30), e154 (v5.5 + 0.154
well, τ4.8, E₀0.36)** — so the MD basin spans both τ values at v5.5
and tolerates the deep well. GR-P2..P6 all confirmed (P5 c50 split:
the control fails as predicted but on the n₁ clause, not n̄ — its n̄
4.35 is *inside* the band, so the basin nearly reaches v_c 5.0). The
instrument itself is now MD-calibrated at the corrected geometry: n̄
twin-hot at 14/14 with Δ 0.24–2.66 He (small-end, residence-scaled),
n₁ transfer ≤ 0.036 everywhere. **f725 (finc1v725 at the corrected
geometry) measured broken: trap 0.437 (bound 0.278 + marginal 0.159),
supp 0, n₁ 0, n̄ 14.3.** Policy evidence: the marginal (modelling-
exclusion) class is ≈ 0 at the basin (≤ 0.018) and 0.159 at the
standing chord — the coupled-class question largely *evaporates* at
the re-arbitrated point (§3.5b over-dissipation reading, MD-grade).
KE trade-off at the gated cells (reported): midHot 0.94/0.50/0.81/0.96,
deepKE 0.43/0.72/0.27/0.43 (a037/b031/d030/e154; standing pooled
reference 1.011/0.631) — no gated cell holds both axes at once; the
mid-vs-deep KE tension persists at the corrected geometry (RQ11
successor question, now at v_c 5.5–6.0). G4 (successor-point
adjudication + re-baseline) decides among these on user authority.

## 15. Sampling laws — size distribution + radial positions **(provenance AUDITED 2026-07-26; ⟨N⟩-pin influence MEASURED by grid re-weighting 2026-07-27 — §15.7; distribution-level A/B remainder open)**

### 15.1 The laws as implemented

- **Size law** (`sampling/droplet_sizes.py`): nozzle correlation
  ⟨N⟩ = k₁·p^0.97·T^−3.88·d² (Lackner); log-normal δ = 0.625 (Kornilov);
  pickup + evaporation chain (mean free path 4 Å, ε = 0.04/collision,
  **E_solv = 14 meV vs thesis-text 30 meV — documented discrepancy**);
  analytic D4 prior family (`droplet_size_prior`, Slice T8).
- **Position law** (`sampling/radial_positions.py`): thermal-equilibrium
  ansatz p(r) ∝ r²·exp(−U_drop(r−R)/k_B T) at T_particles_K
  (`binding_energy_molecule_meV`, steepness), or the twin's L1 law
  `uniform_volume` + hard margin (Slice T7).

### 15.2 Which sampler actually ran, per production (code + distribution verified)

| production | size | position | realized geometry |
|---|---|---|---|
| legacy MATLAB (thesis) | pickup MC, source-condition ⟨N⟩ | Boltzmann | N̄ ≈ 16.4k He, R̄ ≈ 54 Å |
| `main`-branch Python port | pickup MC (`mode="post_pickup"`) | Boltzmann | N̄ 16.4k–16.7k, R̄ 53.9–54.2 Å |
| **Tier-2 drag (standing point)** | **analytic `kornilov_lognormal`, ⟨N⟩ pinned 2000** | **`uniform_volume`, margin 3 Å** | **N̄ 2.0k, R̄ 26.6 Å** |
| Tier-0 / HeDFT presets | fixed N = 2000 | center-pinned | R = 28.0 Å, r₀ ≡ 0 |

- **Legacy production did use both samplers.** `run_simulation.m:66` runs
  `inputfiles_dft_comparison/single_pulse_droplet_distribution.m`
  (`use_single_droplet_size=false`, `single_initial_position=false`,
  8000 molecules, 40 mbar / 14 K); `vmi_sim_3d_neutral_propa_HeDFT_mimic.m`
  calls `generate_droplet_sizes(...)` at `:141` and
  `generate_radial_samples_3d(...)` at `:196`. No later assignment
  overwrites `droplet_radii` (the only constant-N assignment is `:145`,
  inside the `use_single_droplet_size` branch). Trap, harmless:
  `generate_droplet_sizes.m:16-17` hardcodes `p0=40; T0=14`, clobbering
  the globals with the same values.
- **Distribution-level identification** (flags alone cannot separate raw
  from post-pickup): correlation ⟨N⟩ = 12794; sampler `raw` → mean
  12589 / med 10384 / R̄ 49.4 Å; sampler `post_pickup` → mean 16267 /
  med 13400 / R̄ 53.8 Å (the +30 % mean shift is the N^(2/3) pickup
  weighting, exp(⅔δ²) = 1.30, net of evaporation). Stored runs:
  `single_pulse_droplet` 16392 / 13571, `single_pulse_droplet_long`
  16658 / 13697 → **post-pickup**; `single_pulse_droplet_log_droplet`
  12843 / 10547 → raw (the named log-normal comparison run).
- **The drag branch bypasses both legacy laws.** The leg-D generators
  (e.g. `gen_tier2atlas_lqbattery.py:94-103`) pin
  `SIZE_PRIOR="kornilov_lognormal"`, `BIRTH_LAW="uniform_volume"`,
  `BIRTH_MARGIN_ANGSTROM=3.0`, `SINGLE_INITIAL_POSITION=False`.
  `single_droplet_size=2000` from the 9 Å preset is inert; the number
  that replaced it is *also* 2000 because no generator ever overrides
  the config default `droplet_prior_mean_N=2000` — the twin's Wave-10/11
  D4 pin, inherited from the 9 Å TDDFT droplet, **never re-derived from
  the source conditions**. Consequence: the legacy pickup MC and the
  E_solv 14-vs-30 meV discrepancy are *inert* at the standing point.

### 15.3 N → R convention (one convention on the propagation path)

Propagation converts with **bulk** density:
`droplet_radius_bulk_angstrom` = (3N/4πn_He)^(1/3) at n_He = 0.0219 Å⁻³
= 2.2173·N^(1/3) (`simulation/initial_state.py:86`), matching legacy
MATLAB's rounded `2.22·N^(1/3)` (the 0.8·ρ_bulk line is commented out at
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m:151-153`). The 0.8·ρ_bulk
convention (2.39·N^(1/3), **+7.7 % in R**) lives only *inside* the
size sampler — its pickup cross-section, `mean_droplet_size`,
`droplet_radius_from_N` — and is never seen by propagation. Any grid
or re-weighting that converts N → R must use the bulk form.

### 15.4 The Boltzmann position law is strongly R-dependent — the parent model's own behavior

`U` depends on `r − R` only, with a fixed erf width 14.3324 Å and well
573.3 K = 49.4 meV at T = 0.4 K, so the allowed shell is `r ≲ R − 35 Å`
at any R: an interior shell at legacy R̄ ≈ 53 Å (r/R q10/med/q90
0.131 / 0.286 / 0.450) but effectively center-pinned at R ≈ 28 Å
(median r₀ = 1.36 Å, V0-3 `birthlaw` oracle). That is why T7 introduced
`uniform_volume` — the legacy law was structurally inert at the
production droplet size, not rejected on physics.

**Parent-document anchor.** The parent computes the molecule position
from *"the thermal velocity of the iodine molecule and its droplet
solvation potential"* at 0.4 K (its Fig. 6.9), quoting mean solvation
depth 29 Å at R = 34 Å → 40 Å at R = 68.3 Å:

| well depth used | depth @ R = 34 Å | depth @ R = 68.3 Å |
|---|---|---|
| parent document | 29.0 | 40.0 |
| ours at β₂ = 26.99 meV (the DFT fit) | **29.34** | **40.43** |
| ours at 573.3 K = 49.4 meV (as-built) | 30.35 | 41.80 |

**Epistemic status: port-fidelity reproduction, not validation** — the
parent computes the same ansatz we implement, HeDFT entering only as the
potential *shape*; depths stay Derived, 0.4 K a stated assumption. The
parent used the fit's own β₂ = 26.99 meV, so the code's 573.3 K carries
a ~1.3 Å systematic (open item; inert today).

**Measured departure at the standing point** (`bigc1v725s1`, N = 1000):

| birth law on the production droplets (R̄ 26.6 Å) | ⟨r₀⟩ | mean depth | mean isotropic chord |
|---|---|---|---|
| production `uniform_volume`, m = 3 Å | 17.6 Å | **9.0 Å** | **21.7 Å** |
| Boltzmann (the parent model) on the same droplets | 1.6 Å | 25.0 Å | 26.5 Å |

`uniform_volume` **cannot** reach the parent depth at this R
(`⟨depth⟩ = 0.25·R + 0.75·m` would need m = 24.4 Å > R), so at R ≈ 27 Å
only a centered law is parent-consistent — but the chord differs by only
22 %, i.e. birth depth is the secondary term and droplet size the
primary one (§15.5).

### 15.5 Consequence for D2b — the variant that matters most

Legacy production vs standing point: **R̄ 54 → 26.6 Å (×2.0), N̄ 16.4k
→ 2.0k (×8.2)**. Since exposure ≈ the chord `R − r₀` (§14), that gap is
larger than the whole proposed §3 grid span (20 → 34 Å). The
first-class D2b variant is therefore **the ⟨N⟩ pin itself** (2000 vs
the source-condition 12794 and its pickup-weighted 16.3k), ahead of
δ perturbations and the E_solv 14/30 meV pair, which are inert at the
standing point. Whether the standing landing survives at the
experiment's own source-condition droplet size is an open model risk,
not a settled choice.

**External confirmation of the source-condition ensemble (2026-07-26).**
The parent document quotes its droplet range as R = 34 Å (smallest) to
68.3 Å (largest) ⇒ N = 3605 / 29227, ratio 8.11. Against the
source-condition ln-normal (⟨N⟩ = 12794 at 40 mbar / 14 K, δ = 0.625)
those sit at quantiles **0.043 and 0.949** — i.e. the ~5–95 % range of
exactly that distribution. Its **raw** quantiles (q05 34.9 / q95
68.7 Å) match the quoted pair; post-pickup would give 37.7 / 74.1 Å. So
the parent ensemble is the source-condition distribution, **raw, not
pickup-weighted** — which both fixes the external anchor at
⟨N⟩ ≈ 12.8k and pre-answers one §4.2 A/B item.

**Exposure ledger (mean isotropic chord, the quantity drag integrates):**

| geometry | R | mean depth | mean chord |
|---|---|---|---|
| standing point (⟨N⟩ 2000, `uniform_volume` m 3) | 26.6 | 9.0 | **21.7 Å** |
| same droplets, parent birth law | 26.6 | 25.0 | 26.5 Å |
| parent small droplet | 34.0 | 30.3 | 33.8 Å |
| parent large droplet | 68.3 | 41.7 | 64.5 Å |
| legacy production (R̄, Boltzmann) | 54.0 | 37.8 | 52.2 Å |

**≈ 2.4–3× less helium traversed than the parent geometry.** The drag
law itself is unaffected (b is calibrated per unit density behind the
ρ̂ gate; interior density is R-independent, so no Tier-0 re-extraction),
but (v_c, τ, E₀) were arbitrated at this exposure and trade against it
through K (I47/I72) — an exposure change of that size is outside every
knob in this atlas.

**Atlas target (remaining):** the grid re-weighting is **EXECUTED**
(§15.7 — the ⟨N⟩-pin variant measured at ensemble level, zero MD); left
open are the §4.2 distribution-level A/B plots, the ≤ 2 MD confirmations
(the below-support R ≈ 20 Å × L3 cell is now the designed first
candidate, §15.7 oracle read), and the salvageability verdict, which is
G2-coupled: the standing ⟨N⟩ = 2000 pin is *known wrong* against the
source conditions, so "landed under the standing pin" is no longer a
defensible-variant statement.

### 15.6 Reachability of the corrected ensemble — the G0-1 decision (2026-07-26)

Code fact established this session: **`mode="raw"` is unreachable from
config.** `simulation/initial_state.py:83` hardcodes
`sample_droplet_sizes(cfg, mode="post_pickup")` for the `legacy` arm; the
selector exists only as a function argument
(`sampling/droplet_sizes.py:141`). The three candidate corrected-ensemble
arms therefore differ in cost, not just in physics:

| arm | reachable today | realized ensemble | cost |
|---|---|---|---|
| legacy `post_pickup` | yes | ⟨N⟩ 16267, R̄ 53.8 Å — the *pickup-weighted* ensemble, which §15.5 showed is **not** the parent's | zero |
| **legacy `raw`** ← adopted | **no** | ⟨N⟩ 12589, R̄ 49.4 Å; q05/q95 34.9 / 68.7 Å vs the parent's quoted 34 / 68.3 Å | one config field + plumb |
| analytic prior re-pinned at ⟨N⟩ 12794 | needs `DROPLET_PRIOR_N_HI` raised from 16000 | truncated family; ≈ 25 % of the mass sits above the present window at δ = 0.625 | constant change + re-pin of the twin's D4 family |

**Decisive point:** under `legacy` the mean is not a knob — it comes from
the nozzle correlation evaluated on the preset's **own** p = 40 mbar /
T = 14 K, which *is* ⟨N⟩ = 12794. The corrected ensemble therefore needs
**no new number**: no re-pinned mean (the `legacy` guard at
`config.py:1152` forbids an off-default `droplet_prior_mean_N` anyway), no
widened truncation window. It also retires the analytic prior family
along with the ⟨N⟩ pin — one **S** row more than §17 credited.

**Adopted:** a `droplet_size_sampler_mode ∈ {raw, post_pickup}` selector
defaulting to `post_pickup`, so every existing run is bit-for-bit
unchanged (no default-scope change). Not required for G1 (fixed-R cells);
required for G3/G4 and the D2b confirmations. Built with the stage-2b
generator under its trigger.

### 15.7 Grid re-weighting — the ⟨N⟩-pin influence at ensemble level (D2b §4.3 EXECUTED 2026-07-27, zero MD)

`scripts/post_processing/tier2atlas_geometry_reweight.py` re-weights the
§14 grid with candidate size densities. Method conventions (frozen
before the first run): **column-matched 1-D re-weighting** — the
candidate position laws coincide with grid columns exactly (production
`uniform_volume` m3 ≡ L3; corrected Boltzmann 313.2 K ≡ L2), so the
conditional birth depth at each R is the cell's own and only the R
marginal is interpolated (nearest-midpoint vs piecewise-linear, clamped
outside support with the clamp reported); the mixture is scored by the
**committed scorer** through an exact weighted sufficient-statistics
read (test-locked against literal pooling); `χ²_med` is intentionally
absent (its sim-SE convention has no weighted analogue).

**Oracle verdicts.** The pre-registered §4.3 oracle — reconstruct the
pooled N = 5000 row from the standing density — is **INADMISSIBLE, 3/7
columns pass**: **53.8 % of the standing density lies below the grid's
R = 26.6 Å support edge** (the support starts at the standing *mean*),
and the clamp overshoots trap (+0.027, 9 SD), n̄ (+0.58) and W₁ (+0.34).
A **post-hoc in-support oracle** (declared post-hoc; the battery's own
R ≥ 26.6 Å sub-ensemble, 46.2 % of its ions, support-covered by
construction) passes **7/7 on both conventions** (`nearest` adopted;
worst error W₁ −0.18 ≈ 1.9 SD, all others ≤ ~1 SD). Reading: the
*method* is valid where the support covers the density; the full-density
failure is the support hole. The corrected density is covered to ~95 %
(clamp 0.2 % below / 5.1 % above), so its forecast inherits the
in-support error scale, not the full-density failure. The designed
remedy for the hole — only needed if a standing-mixture reconstruction
is ever load-bearing — is a below-support R ≈ 20 Å × L3 cell (the §4.3
"≤ 2 confirmations" budget).

**Forecast (nearest convention; Arm A = marginals excluded / Arm B =
marginals injected at exact conservative asymptotic KE — the §14.2
NOT-TIGHT bracket propagated to ensemble level; quote deepKE and n̄ as
A–B ranges):**

| mixture | trap (marg) | det_yield | supp | n̄ | n₁_solv | W₁_solv | midHot | deepKE |
|---|---|---|---|---|---|---|---|---|
| pooled recorded | 0.067 | — | 0.187 | 4.07 | 0.243 | 0.571 | 1.011 | 0.631 |
| std × L3 (oracle recon) | 0.094 | 0.91 | 0.161 | 4.65 | 0.224 | 0.912 | 1.04 | 0.599 |
| std × L2 | 0.005 | 1.00 | 0 | 9.29 | 0 | 4.67 | 1.252 (4 bins) | 1.47 |
| corr × L3 | 0.35–0.38 (0.032) | 0.62–0.65 | 0.12–0.13 | 6.0–6.7 | 0.15–0.16 | 2.07–2.69 | 1.254 | 1.35–1.45 |
| **corr × L2 (corrected geometry)** | **0.31–0.42 (0.110)** | **0.58–0.69** | **0** | **13.9–14.7** | **0** | **9.0–9.8** | **1.256 (4 bins)** | **1.80–1.90** |

**Reads.**

1. **The corrected geometry breaks the landing at ensemble level,
   through the §3.1-predicted channel:** suppression → 0 and the low-n
   half evacuates (n₁_solv 0, W₁ ≈ 9–10), with trap 0.31–0.42. What the
   anchored *cells* showed (§14.1), the anchored *ensemble* confirms —
   this is the quantitative form of the G2 clause "the landing is
   expected to break".
2. **Decomposition:** the **birth law owns the histogram breakage**
   (std × L2 already gives W₁ 4.67, n̄ 9.3, supp 0 at unchanged sizes,
   with trap ≈ 0.005) and the **size distribution owns trapping**
   (corr × L3: trap 0.35–0.38 at supp 0.13); each alone pushes deepKE
   past 1 (1.35–1.47), together 1.80–1.90 — the two axes are separately
   fatal, not one compound effect.
3. **RQ11 at ensemble level:** deepKE crosses 1.0 between the standing
   and corrected geometry along *either* axis — the ensemble echo of
   §14.1's depth-crossing at ≈ 11–15 Å. The corrected geometry
   *overshoots* the deep-bin KE (too hot) by ×1.8–1.9.
4. **Detected-subset bias:** det_yield 0.58–0.69 at the corrected
   geometry — the detected ensemble is a small-R / shallow-birth biased
   subset of the source ensemble, so every detection-conditional
   observable there must be read jointly with the fate split.
5. **Caveats:** grid cells are N = 500 single-seed; midHot at the L2
   mixtures rests on 4 of 7 band bins; the marginal fraction (0.11 of
   source ions at corr × L2 Arm A) is conditional on the cubic law
   extrapolated ~3× past its calibration band (§14.2); the in-support
   interpolation error travels with every row.

## 16. Quiet / structural surfaces (for completeness)

- **s_eff (RRK dof, Bounded):** a *rate*, not an amount (I13);
  arm-conditional landing (ungated ≈ 8, gated ≈ 30 — the ungated prior
  does not carry, I10); the detector compresses it (arrival state, not
  arrival n) *except* the n = 1 bin, which selects the kinetics band
  (I20, I32). Standing: the 8-family via the production bundle.
- **E2 relaxation stage:** KE-frozen at handover on every solvated bin
  (E_dissip gain 0.0000; Landau-gated E2 drag touches only the retained
  class); detection-stage sheds KE-negligible — the post-window chain is
  causally disconnected from the KE curve (I89); E2 checkpoint stride is
  scored-read-neutral (I67).
- **Detection time:** Sourced 8.53 µs; ungated reads are
  detector-converged (n_detect ≡ n_relaxed); the gated long-time
  structure is the Klots evaporative-ensemble regime (I19, I23).
- **Scenario budget (0.80/2.70 eV):** bookkeeping-only in the delivered
  model; the gated map is budget-invariant in absolute E₀ (I28).
- **Guard flags:** biphasic runs under `allow_inconsistent_mass_pairing`
  (§6.6 mid-window defense); E_bind sweeps will use
  `allow_unvalidated_binding_pairing` (§6.5.1). Deliberate, logged.

---

## 17. Physical-sensibility ledger (added 2026-07-26)

Influence answers *what a knob does*. This chapter answers a different
question, the one the atlas exists for: **which parts of the model are
physics, and which are scaffolding?** Role codes:

- **P** — physics-constrained: external evidence pins it (TDDFT traces,
  literature, source conditions, experiment).
- **C** — convention: a choice that must be made; documented, defensible,
  not free.
- **E** — effective: no external constraint; fitted, and absorbing
  whatever physics is missing there.
- **S** — scaffolding: exists to compensate another element's wrong
  value; would retire if that element were fixed.
- **M** — missing: known-absent physics, tracked so it is not mistaken
  for zero.

| element | provenance | external constraint? | role | what would retire / confirm it |
|---|---|---|---|---|
| in-band drag form + b | Tier-0 Method-B on 9/18 Å TDDFT | **yes**, 2.54–4.95 Å/ps only | **P** (band-limited) | extended TDDFT (collaborator ask). Note §1: the Tier-2 landing is form-blind — it is *not* evidence |
| capped tail v_c, p_tail | fit to experiment above the TDDFT band | no | **E** | TDDFT at production kinematics; the "S-shape" framing (§6.1) says a cap is a stand-in, not a law |
| τ cooling clock | arbitration basin §4w | no | **E** | a physical cooling model for the bubble |
| E₀ = E_int(0) | solvation-scale bracket + arbitration | weak (bracket only) | **E** | Axis B curve + an independent onset argument |
| κ, picture | Form-U parametrisation | no | **E** (both near-dead levers) | superseded by the tabulated `rq4graded` ladder |
| ladder D₀(n) rungs | RQ4 literature | **yes** (rungs) | **P** rungs / **E** graded shape | shape: an independent size-resolved binding source |
| λ₀ pickup, detection time | source conditions / flight geometry | **yes** | **P** | — |
| Landau v_L 0.58 | literature bulk value | **yes** | **P** (measured quiet) | — |
| per-shed ε | bounded small (NB-RQ23-1) | weak | **C** (set 0) | a recoil model |
| shed convention (cold/co-moving) | two-valued bookkeeping choice | no | **C** | RQ3 resolution |
| E_bind 0.1168 eV | joint Method-B with the drag pair | **yes** (jointly) | **P** (paired) | measured to matter (§9); pairing is the constraint, not the value. R-dependence **bounded** ≤ 0.009 eV over the Axis A span (§9.1, Born far-field; the local snowball term cancels) |
| **⟨N⟩ = 2000 droplet pin** | Tier-0 TDDFT droplet → S2-D4 freeze → twin D4 family | **contradicted**: parent ensemble is ⟨N⟩ ≈ 12794 (§15.5) | **S** | adopting the source-condition prior (Axis A / D2b) — which under `legacy` needs no pin at all, the nozzle correlation supplies it from p/T (§15.6) |
| **analytic `kornilov_lognormal` prior family** | T8 twin-parity device carrying the ⟨N⟩ pin; window [250, 16000] He | no — and the window cannot hold the source-condition distribution | **S** | retires with the pin: the `legacy` + `raw` sampler is the parent's own law (§15.6) |
| **`uniform_volume` birth law** | T7 twin parity, because Boltzmann is center-pinning at R ≈ 28 Å | **contradicted**: the parent model is the thermal law (§15.4) | **S** | adopting the parent birth law — which is *also* our own ported sampler |
| **margin 3 Å** | pinned convention (I88) | no | **E**, and load-bearing (n₁_solv 3.4σ at 6 Å) | a birth law that needs no margin (i.e. the parent law) |
| **T5 `density_tied`** | dressing for under-dressed shallow births | no | **S** — inert at the parent geometry (§G) | adopting the parent geometry. *No control can decompose it* — it shares one `rho_he_ratio` surface with the cooling gate and the drag gate (§14) |
| **T6 `sigma_proportional`** | onset scaled by Σ(n₀)/Σ(n*) | no | **S** — same; and it silently imposes `E_int/Σ = 1.309` at every birth | same: the geometry correction, not a control cell (C1/C2 rejected 2026-07-26 — each needs the less physical arm of its pair) |
| `cooling_spatial_gate="density_scaled"` | ρ̂-scaled Newton drain | no | **E**, but the alternative is unphysical | nothing available: `"none"` means cooling with no bath to cool into. Third leg of the same shared ρ̂ surface |
| Boltzmann well 573.3 K | code value (legacy MATLAB); the DFT fit it cites is 26.99 meV = 313.2 K | **mismatch** (~1.3 Å in depth) | **provenance defect** | **G0 decision 2026-07-26:** Axis A L2 cells override to 313.2 K **per run**; the config default stays 573.3 K (legacy port fidelity — a default change buys nothing while the arm is unused in production). 573.3 K is a robustness check, not a cell |
| mass / binding pairing hatches | §6.5–§6.6 documented exceptions | n/a | **C** (logged per run) | — |
| ensemble second moments | not modelled | experiment says under-dispersed | **M** | Tier 3 (discrete-emission drag is the standing candidate) |
| **three-stage timescale separation** (30 ps drag-active MD → 8 ns conservative E2 → 8.53 µs free flight) | staging calibrated at R ≈ 27 Å, where ejection is effectively instantaneous | **contradicted at the anchored radii**: at R ≥ 49 Å escape takes ~0.1–1 µs with He still present, and ⅓–½ of ions never leave (§14.2) | **S** — scaffolding for the wrong geometry; **no longer a blocker** | **partly retired 2026-07-27** by the second option: the µs-orbiting class *is* defined as retained (`exclude_all_coupled`), decomposed into bound (physics) and marginal (convention), so the grid is readable — §14.2. What remains open is whether the µs residence is real at all, i.e. whether cubic drag over-dissipates on 25–50 Å paths (the collaborator ask; a *geometric* retained fraction is insensitive to (v_c, b), an over-dissipation artifact is not). **The first ledger row the correction created rather than retired** |

**What the ledger says.** The **S** rows share one root: the droplet is
~2× too small and births ~3× too shallow, and several arms exist to
describe or compensate for that. Adopting the anchored geometry retires
**six rows at once** — the ⟨N⟩ pin, the analytic prior family that
carries it (§15.6), the `uniform_volume` birth law, the margin knob, and
both dressing arms — i.e. it *simplifies* the model rather than
complicating it. That is the strongest
physical-sensibility argument the program has produced, and it is the
reason Axis A was re-cut as a geometry test (plan §3.1b).

The **E** rows are where the model is fitted rather than known: the
capped tail, τ, E₀, the graded ladder shape, the margin, the cooling
gate. Each is a candidate for an external constraint we do not yet have;
the collaborator ask (extended TDDFT) covers the first, Axis B measures
the next two, and the geometry decision removes the fourth and fifth.

**Ledger discipline.** Role codes are claims about *evidence*, not
quality — an **E** entry is not wrong, it is unconstrained. Nothing in
this chapter adopts anything; adoption stays a separate pre-registered
decision (plan §0). The ledger's job is to be the dossier that decision
would read.

---

## Cross-references

- Archive/provenance: `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee,
  I1–I100, boundaries §6, open questions §7, run inventory §8).
- Program that fills the GAPs: `TIER2_SENSITIVITY_ATLAS_PLAN.md`.
- Calibration classes + identifiability: `CALIBRATION_MAP.md`.
- Open physics questions: `RESEARCH_QUESTIONS.md` (RQ1–RQ11 + NB
  registers).
- Decision history: `drag_migration_log_tier2.md`.
