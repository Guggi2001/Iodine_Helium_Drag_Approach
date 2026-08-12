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
| E_bind (ion–droplet well) | Derived (joint Method-B) | **swept §6.7 item 2**: trap +0.058/0.1168-step (clean well lever), n̄ −0.50, midHot −0.076; over-suppression is the FORM, not the well. **KE transfer function measured (§9.2, 2026-08-10):** affine, **dKE₁/dE_bind = −0.542 (h405) / −0.549 (lin)** eV/eV — a *fixed ~45 % refund*, form- and depth-invariant; trap lever form-SPLIT (0.620 vs 0.033 /eV); whole axis worth ≤ +0.064 eV on KE₁. **MD-CONFIRMED (§9.3): −0.5308 measured, 2 % from the twin**. **MECHANISM SOLVED (§9.6, 2026-08-11): the "45 % refund" is a MASS-FRAME PARTITION, not drag** — the ion pays **98.1 %** of the well; the shortfall is 38.0 % `m(1)/m(21)` + 6.7 % cascade + **1.9 % genuine drag refund** (form-split 3.3× vs the lin arm's 6.2 %). Density-width scan **NULL** (c 0.5422 → 0.5597 over a 4.5× sharpening ⇒ +0.002 eV on KE₁) | measured (§9); KE half of the GAP **closed** (§9.2) + MD-confirmed (§9.3) + mechanism + width GAP closed (§9.6) |
| droplet geometry (R × r) | controlled (Axis A G1 11-cell grid); **size externally anchored** (§15.5) | birth depth is the physics knob, R a selection knob; the landing needs the size *distribution* (pinning R̄ alone: W₁ 0.571 → 0.813); ≈ 95 % of detected-size variance geometry-inherited; deepKE crosses 1 at birth depth ≈ 11–15 Å; trap → 0.40–0.57 at the anchored radii | **GAP closed** (§14); **G2 ADOPTED 2026-07-27** — the corrected geometry is the target, re-arbitration pending (G3) |
| sampling laws (size + position) | theory-laden legacy ports, **partly bypassed in the drag branch** | provenance audited (§15): production uses the analytic ⟨N⟩ = 2000 prior + uniform_volume, *not* the legacy pickup MC + Boltzmann; E_solv 14 vs 30 meV discrepancy is inert here; **⟨N⟩-pin influence measured at ensemble level by grid re-weighting (§15.7, zero MD): the corrected ensemble breaks the landing** (trap 0.31–0.42, W₁ ≈ 9.0–9.8, deepKE 1.80–1.90) | ⟨N⟩ pin **measured** (§15.7); distribution A/B remainder open (D2b) |
| drag state coupling s(n) (R_core, ρ_shell) | Bounded geometric closure (design doc) | **MEASURED DEAD in-window (§18): gate-clipped** — no ion reaches n ≤ 8 while inside the droplet (min-n-inside ≥ 9 for 100 %, mean exit n = 19.0), so the low-n regime of s(n) is structurally unreachable; in-window s ∈ [0.85, 1.06] ≈ ρ-independent; net effect KE₁ **−0.04 (sign-inverted)**, needle SD unmoved (0.038), all SC predictions refuted | probe EXECUTED 2026-07-29; axis stopped (SC-P2 signature failure); code stays behind `off` default |

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
- **lq runs COLDER at n = 1/n = 2 — the direct KE₁ read (§6.7 item 3,
  2026-07-30, committed `tier2atlas_lq_ke1_table.py`,
  `atlas_lq_ke1_table.csv`).** Measured on the regenerated §6.7 battery
  (the lq run dirs had been cleaned locally; regenerated seed-exact from
  the committed generator and gated by a printed-precision reproduction
  oracle on the committed §6.7 rows — all 10 members + the pooled row
  reproduce). Paired (lq − cubic), 5 × N = 1000: **ΔKE₁_mean
  −0.047 ± 0.005 eV (lower on every seed, ~10σ paired)**, ΔKE₂_mean
  −0.060 ± 0.004, Δabove-1.15 −0.093 ± 0.018 (pooled 0.095 vs 0.188;
  reference 0.487), ΔKE₁_SD +0.005 ± 0.002 (no width gain), Δn₁_solv
  +0.001 ± 0.006 (no population gain). Pooled KE₁: lq 0.987 ± 0.132 vs
  cubic 1.034 ± 0.128 (reference n = 1 mean 1.302). **This closes the
  "could a lower-power form raise the low-n KE" question with the
  direct observable: below the TDDFT-calibrated band a lower power
  decays slower ⇒ MORE drag on slow ions — quadratic moves KE₁/KE₂ the
  wrong way, buys no width and no n₁, and keeps its over-suppression.**
  Consistent with (and sharper than) the item-1 χ² read; the drag-form
  KE lever at low n is measured DOWNHILL-ONLY toward lq. Pooled lq
  figure set rendered (`..._lq_N5000_tier2atlas_conf270_qccbigpooled`,
  pair-preserving container via the committed
  `build_pooled_detection_container.py`, builder-oracled against
  `bigc1v725pooled`). **Provenance caveat (user-raised, 2026-07-30):
  this battery sits at the PRE-CORRECTION geometry (leg-D
  uniform-volume birth, old droplet ensemble) — the corrected-geometry
  lq-vs-cubic KE ordering is UNMEASURED, and the Δ is a whole-path
  integral (lq is softer in the 5–9 mid-band yet colder at n = 1/2
  net), so the ordering does not transfer by argument alone. Candidate
  resolving test posted (log 2026-07-30): lq twin scan at the
  corrected geometry for lq's own basin, then a CRN-paired MD ring vs
  h405 with frozen KE/fate bands — NOT adjudicated.**

- **The pure-linear counterfactual holds the corrected-geometry twin
  gate INCLUDING at sub-plateau tail force — kill-3's universality
  BREAKS at twin level (free-form linear sweep arm 1 `linscan`,
  2026-07-30, zero MD; `TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md`
  §3.1/§5.1, `atlas_linsweep*.csv`).** γ = ρ̂·a swept a ∈ [15, 60] × 3
  wells × the §3.5c free surface (12 312 cells, frozen gate verbatim):
  **35 cells gate** (R1), τ = 4.8 dominant (28/35, zero τ-flags), all
  three wells (shallow-favored 19/9/7), a core 27.5–42.5. **Φ =
  F(9.7)/418.5 spans 0.64–1.39 with 6 sub-plateau cells (Φ ≤ 0.7),
  none below the bare-ram floor** — the histogram gate does NOT pin
  the fast-class force inside the linear family. **R4 (twin KE₁
  ranked under the Step-0 licensure, levels ~2 % cold): softer tail →
  hotter n = 1 monotonically — best gated twin KE₁ 0.892 at a = 27.5
  (Φ 0.64) vs the h405 twin anchor 0.624 (+0.27).** Twin-level
  caveats, pre-registered: deepKE (unlicensed, ρ 0.33) runs 2.0–3.4×
  HOT at the sub-plateau cells — the declared MD-undecidable axis;
  above-1.15 = 0 (the twin carries no source-KER width by
  construction); trap is a twin floor and the lin form's low-v
  over-drag makes MD trap the risk axis. **Nothing adopted; the MD
  ring (plan §6) is the decision instrument and sits behind its own
  trigger.** **Arm 2 (`linqscan`, a* = 35.0 W₁-optimum, R5 freeze;
  seam oracle bit-exact): the landing REJECTS the ram term** — at
  fixed a* every observable degrades monotonically with c (W₁ 0.586 →
  1.02, twin KE₁ 0.703 → ≤ 0.25, trap → 0.93–0.97 by c = 12.79 where
  nothing gates — the shared-lq landmark is dead in-family; far-c
  stragglers τ-flagged). Within the free-form family the data selects
  **pure linear**; the MD-ring candidate set is the arm-1 sub-plateau
  core a 27.5–35.

- **The MD ring CONFIRMS the twin landing — the drag-side KE₁ door is
  MD-OPEN (free-form linear ring, 2026-07-31, 8 × N = 500 CRN seed
  20260731 vs the h405 partner; plan §6.3, `atlas_linring_table.csv`).**
  6/7 lin cells pass the MD-side hard gate under BOTH retained-policy
  ends (zero policy-blocked; lr7's miss was twin-predicted — the gate
  pattern transfers 7/7). **K-KE does NOT fire: ΔKE₁ = +0.077…+0.266 eV
  at every gated cell; five cells beat the (A)-ceiling success band
  (KE₁ ≥ 0.75), best KE₁ 0.903 at a = 27.5/eb0482/E₀ 0.35 vs h405p
  0.637 — the first MD-measured drag-side KE₁ result above the exit-
  stripping ceiling.** The capped→lin twin transfer HELD: twin KE₁
  uniformly ~1–2 % cold in-family (d(twin−MD) −0.010…−0.015 eV), n₁
  transfer ≤ 0.006, n̄ twin-hot bias uniform ≈ +0.46…+0.59 — the
  lin-family authority box is a near-copy of the capped one. **The two
  pre-registered kill axes both came back benign:** trap ≤ 0.011 (all
  bound, zero marginal; the low-v over-drag fear is DEAD — the CRN
  fate flow runs the other way, lin *frees* the ~60 ions h405p traps,
  new-trapped = 0 at every cell); deepKE mixed/neutral (2.37 at the
  KE₁-best cell down to 0.95 at a = 35 — the twin's 2.0–3.4× hot flag
  was itself ~30 % hot; no form-question). **Honest costs, reported:**
  midHot 1.39–1.98 (the soft band is 0.85–1.15) and χ²_med 795–2160 vs
  h405p's 267 — the n = 1/2 heat is bought with mid-band overheating
  (the NB-RQ11-12 coordinate, now measured from the lin side); W₁
  0.78–0.97 vs h405p 0.767 (±0.13 single-seed, not gating);
  above-1.15 = 0 everywhere (source-KER width stays (B)-side, as
  designed). h405p itself sat just under the n̄ band this seed (3.748
  vs 3.77; Arm-B flips it in) — disclosed Δ-only per §6.2 item 5.
  **Nothing adopted; per plan §0 an MD-level win on the form-sensitive
  observables opens the adoption discussion — that discussion + the
  pre-registered escalation (one N = 1000 × 5-seed battery at the
  single best cell; KE₁-best lr1 vs W₁-convention lr6 is part of the
  call) are the NEXT USER GATE.**

- **The linear family has TWO personalities, and the second one
  CLONES h405 with no cap — the choice between forms is *where the
  deviation is hidden*, not whether there is one (adoption-discussion
  re-read of the committed arm-1 rows, 2026-08-10, zero MD; plan
  §6.4).** a is a continuous dial: Φ 0.64 buys KE₁ 0.892 at midHot
  2.16/deepKE 3.38 (the ring's MD-confirmed corner), while **Φ 0.985
  (a = 42.5 / eb0482 / τ 6.4 / E₀ 0.31) reproduces the h405 twin
  vector near observable-for-observable — KE₁ 0.615 vs 0.624, midHot
  1.003 vs 1.069, W₁ 0.686 vs 0.703, n₁ 0.197 vs 0.208, n̄ 4.58 vs
  4.48 — from a kink-free pure-linear law**, at the price of the KE₁
  gain (the user's proposed trade, confirmed). **But the cap is not
  thereby retired:** crossover √(a/b) = 4.11 Å/ps means that same
  line runs ×4.2 the Tier-0 cubic at v = 2, ×1.9 at v = 3 and 31 %
  soft at the band top — the fiction moves *into* the TDDFT-validated
  window (the Tier-0 linear rejection, n̂ = 2.927, restated at system
  level), whereas `capped_cubic` is exact in-band and parks its free
  element above 4.95 where TDDFT is permanently infeasible. Cap
  physics recorded: v_c 5.5 ≈ **Mach 2.3** (He sound ~2.4 Å/ps),
  constant force = constant loss per Å = wave-drag/vortex-shedding
  signature; the *sharpness* of the corner is forced (any smooth
  saturation `b·v³/(1+(v/v_c)^m)` breaks the in-band lock unless
  m → ∞, which is the cap). Posted-not-taken: a Hill-type
  smooth-saturation twin scan to measure how sharp the corner must be
  (plan §6.4 item 5; NOT a §3.5h p_tail re-open — corner vs tail
  exponent). Arm 2 already closed the ram term (c = 0 optimal), so
  the flatness is what the landing wants.

- **The clone is MD-measured a KE-equivalent, NOT a landing-equivalent
  — and the ring's mid-band overheating is a position on the a dial,
  not a property of constant γ (h405-clone battery, 2026-08-12,
  3 × N = 500 seeds 20260731/0812/0813, CRN vs the committed h405p;
  plan §6.5, `atlas_linclone_table.csv`).** At `pure_linear a 42.5 /
  eb0482 / τ 6.4 / E₀ 0.31` the KE side of the §6.4 forecast transfers
  and the histogram side does not:
  - **KE₁ equivalence confirmed** — pooled 0.6348 vs the committed
    h405 battery 0.6406 (Δ **−0.0059 eV**; the same-seed same-N CRN
    pair gives **−0.0016**), against a ±0.05 band. χ²_med per seed
    173/234/241 vs h405p's 267 — the KE *curve* fit is as good or
    better (pooled χ² is N-extensive and not comparable).
  - **midHot 0.891 — cooler than h405 itself (0.946), against the
    ring's 1.39–1.98 at a 27.5–35.** The §6.3 "honest cost" is
    therefore **not form-intrinsic**: constant γ can be mid-band
    neutral; where it sits on the NB-RQ11-12 coordinate is set by a,
    not by the form. (deepKE 0.376, i.e. cold, vs h405p 0.472.)
  - **Trap stays dead at 1.5× the ring's chord:** 0.062 (bound 0.0613,
    marginal 0.0007 = 2 ions) *below* h405's 0.078; CRN fate flow
    new-trapped 1 / freed 9. The low-v over-drag fear does not revive
    at a = 42.5. Suppression collapses to 0.051 vs h405p's 0.216
    (τ 6.4 + E₀ 0.31).
  - **But W₁ 1.079 vs h405 0.765 (+0.31)** while n₁ 0.192 and n̄ 3.99
    sit at/inside their gate bands — the two histogram *moments* match
    and the *shape* does not. Reported-only per the §6.5 user
    adjudication; it cannot overturn the KE reads, and it is what makes
    "clone" the wrong word for anything but the KE observables.
  - **n₁ sits ON the gate floor**: pooled 0.1921 vs the 0.19 edge, gap
    0.0021, per-seed 0.1734/0.1937/0.2095 (SD 0.0181) ⇒ **~72 seeds of
    N = 500** would be needed to separate it from the edge. Verdict
    `gate-marginal` under **both** retained-policy ends (not
    policy-blocked). That is a measurement, not a resolution
    shortfall — this cell is boundary, and no affordable MD makes it
    otherwise.
  - **CL-P6 — the in-family twin↔MD box does NOT extend to this
    corner.** n₁ (twin−MD +0.0045) and n̄ (+0.583) transfer inside the
    §6.3 ring ranges, but KE₁ (−0.0195) falls outside [−0.015, −0.010],
    and **W₁ transfer collapses: twin 0.686 → MD 1.079, a −0.39 bias
    against the ring's ≈ −0.2.** The §6.4 clone claim rested on the
    twin W₁ being within 0.02 of h405's; that specific number is now
    measured wrong by an order of magnitude more than the gap it was
    asserting. Twin W₁ is not transferable at a 42.5 / τ 6.4.

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
different closed form. **Domain input (user, 2026-07-30): the
collaborator ask is INFEASIBLE — TDDFT itself breaks down at
production kinematics. The principled cap-remover is therefore
unavailable: the cap is a PERMANENT effective element (§17 class
unchanged), and no velocity band outside the existing traces can ever
gain ab-initio authority.**

**Velocity anatomy + tail-force invariance (zero-MD arithmetic,
2026-07-30, from the two pooled containers + committed coefficient
bundles through the committed drag module; triggered by the user's
form-logic challenge):** detected |v| by terminal n (old-geometry
pools, co-moving shed preserves v through the vacuum cascade, minus
residual post-exit Coulomb gain): **n = 1 mean 12.3 Å/ps (p10 11.3),
n = 2 10.4** — the n = 1/n = 2 ions live entirely ABOVE both caps
(user-argued premise CONFIRMED; at the corrected geometry KE₁ 0.637 ⇒
≈ 9.7 Å/ps vs v_c 5.5, even deeper in the tail); n = 5 at 6.5
(band top), n = 10 at 3.5 (mid/low band). Production force curves:
cubic tail plateau **958.5** amu·Å/ps² (b 2.515, v_c 7.25) vs lq tail
plateau **990.6** (c 12.79, v_c 8.8) — **the two independently
arbitrated systems, one power apart, land within 3 % of the SAME tail
force, with the LOWER power ending HIGHER** (its landing pushed v_c up
to recover dissipation). The form redistributes drag only BELOW the
cap (lq ×5.1 harder at v = 1, ×0.70 softer at 7.25) — where the
n = 1/n = 2 ions never are. **Reading: at n = 1/n = 2 velocities the
drag is plateau-set, not power-set, and the plateau is
histogram-pinned (the §3.5h p_tail kill is its direct lowering test) —
a lower-power form does NOT deliver lower drag to the fast ions in any
landed system.** **Mechanism correction (user-caught, 2026-07-30):
drag work never enters E_int (it books to E_dissip; E_int sources are
the E₀ seed + the f_ret·D₀ pickup heat, sinks are the density-gated τ
cooling + the per-shed D₀ — `biphasic_step` K1/K2/S1). The
drag↔shedding coupling is a TIME coupling, not an energy pipe:
E_int^exit ≈ E₀·exp(−t_res/τ_eff), and the plateau sets t_res for the
fast class. Reaching n = 1 requires exiting hot, and the ions that
exit hot are exactly the ones that exit fast — the KE↔histogram
anti-correlation (I47) is a CO-SELECTION on residence time. This is
why (v_c, τ) close only jointly (the race coordinate), why §3.5h
softening un-damps the cascade (more ions exit fast-and-hot,
over-shedding small-n while mid bins starve), and why the source-KER
lever is the unique KE₁ mover (it shortens residence: raises exit KE
AND preserves E_int simultaneously).** **Mid-band γ magnitude (band top → cap,
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
**Deep-KE ceiling MEASURED (G4 ladder v525 probe, N = 1000):** dropping
to **v_c 5.25** reaches **deepKE 0.656 in a GATED cell** — above the
standing pooled 0.633 — but at midHot 1.394 / χ²_med 912. So the
mid-band lever has genuine deep-KE headroom and it is paid for in
midHot: ≈ 0.66 at midHot ≈ 1.4, or ≈ 0.60 at midHot ≈ 0.96 (v_c 5.5).
v_c 5.25 also traps least (0.048 vs 0.087 at 5.5). **`p_tail` was not
needed:** the pre-registered Block-2 trigger did not fire — the deep-KE
axis proved parameter-accessible on (v_c, τ, E₀) at fixed form
(§14.5 G4F-P4), so the "p_tail is not a second lever" verdict stands
un-retested at v_c 5.5–6.0 rather than being overturned.

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
**Ordering MEASURED at N = 1000 (G4 ladder, three CRN-paired τ arms at
v_c 5.5, each swept in E₀):** on W₁ at matched gate status,
**τ 4.4 (0.669–0.709) < τ 4.8 (0.733–0.792) < τ 5.2 (0.872–0.892)** —
shorter τ is better at the corrected geometry, i.e. back *toward* the
standing 3.2, not toward 6.55. deepKE follows the same ordering
(0.55–0.61 / 0.53–0.59 / 0.43–0.50). **The twin ranked τ 5.2 best**
(twin W₁ 0.42–0.49) — the second twin-W₁ inversion (§14.4). Successor
candidate τ = 4.4; τ 4.0 is untested at v_c 5.5 (the twin gates nothing
there, but its n₁ gate is the mis-calibrated one — §14.4).

## 4. E₀ = E_int(0) — internal-energy budget

| | |
|---|---|
| role | absolute internal energy at onset [eV]; feeds the RRK gate via the Σ(21) crossing |
| class | Bounded — RQ1 band [0.2, 0.5] eV, floor ≈ 0.22 (row 16; replaced f_int as the physical variable, I25) |

**Influence (measured):**

- **Bin-level signature (G4 Step 2 Block D, CRN-paired τ 4.4 arm, per
  +0.005 eV):** the E₀ lever closes the solvated CDF gaps **only at
  n = 1–4** (Δgap +0.002…+0.012 there; |Δgap| < 0.003 and sign-unstable
  for n ≥ 5). **Bins 7–20 are E₀-inaccessible** — they carry ≈ 45 % of
  the pooled-h405 W₁ residual (five of its nine top-70 % bins), which is
  the measured identity of the W₁ floor's carriers (findings "G4 Step 2
  — Block D").
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

**Slope MEASURED — the cleanest one-knob response in the atlas (G4
Block 3 ladder, 2026-07-28; 15 × N = 1000 CRN-paired, §14.5).** Inside
the corrected-geometry ridge, per **+0.005 eV** of E₀, monotone in all
three τ arms (4.4 / 4.8 / 5.2):

| observable | Δ per +0.005 eV | direction |
|---|---|---|
| W₁_solv | **−0.02 … −0.04** | improves |
| n₁_solv | **+0.004 … +0.007** | improves (toward exp. 0.243) |
| n̄_det | **−0.11 … −0.13 He** | degrades (away from exp. 4.07) |
| midHot | **−0.02 … −0.03** | degrades |
| supp | +0.02 … +0.03 | — |

Each τ arm therefore gates over a window of only **~0.01 eV in E₀** —
E₀ is a sharp knob at the corrected geometry, not a soft one, and it is
the knob that trades histogram quality against n̄ and midHot.
**Consequence (the G4 headline, §14.5 item 3):** because n₁ and n̄ move
in opposite directions under E₀, they cannot be matched simultaneously —
n₁ = 0.243 needs E₀ ≈ 0.435 where n̄ ≈ 3.3, 0.77 He below experiment.
This is what sets the W₁ floor ≈ 0.67 and localizes the corrected-
geometry residual to the *mechanism*, not the drag surface. Corrected-
geometry successor candidate: **E₀ 0.405** at (v_c 5.5, τ 4.4) — the
E₀-band chain's fourth point, back up near the pinned-droplet scale.

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

### 9.2 The KE transfer function — MEASURED (atlas §6.5 Step 2 twin scan, 2026-08-10, zero MD)

§9 carried trap / n̄ / midHot but never the **KE** response, which is the
one the user's first-order argument is about. Measured on the corrected
geometry at two arms × 7 wells (`atlas_ebind_twin.csv`; stage
`ebindscan`): **arm H** = h405 (`capped_cubic` v_c 5.5 / τ 4.4 / E₀ 0.405),
**arm L** = the MD-measured lin chord (`pure_linear` a 27.5 / τ 4.8 /
E₀ 0.35). Both oracles green (arm L reproduces three committed
`atlas_linsweep.csv` wells string-exact; arm H reproduces the committed
h405 twin anchor string-exact).

**The first-order argument, and why it is wrong by a factor of 2.** The
droplet well and the He density gate are the *same* erf at the *same*
width — identically `U(r) = E_bind·(1 − ρ̂(r))`, `cfg.potential_steepness`
= 14.2 Å. Births sit ~35 Å inside a ~48 Å droplet, so `U(birth) = 0` and
for a **fixed trajectory** E_bind is a purely additive, velocity-independent
per-fragment exit toll: `dKE_∞/dE_bind = −1` exactly, uniformly for every
peak. Measured, it is **≈ −0.545**.

> **THE MISSING 45 % IS SOLVED — see §9.6 (2026-08-11). It is not a drag
> refund and there is no defect.** The ion pays **98.1 %** of the well
> (measured two independent ways: a dressed-frame derivative, and energy
> closure at 0.1167 eV against the nominal 0.1168 eV). The scored gap is
> a **mass-frame partition**: the toll is paid by the 211 amu *dressed*
> complex, while KE₁ scores the 131 amu *bare* ion at that complex's final
> velocity, so only `m(1)/m(21) = 0.6205` of the toll reaches the
> observable. Decomposition of the measured 0.5422: **38.0 % mass
> partition + 6.7 % cascade + 1.9 % genuine drag refund**. Every
> statement below about a "refund" being 45 % should be read as the 1.9 %
> figure; the *measured numbers* in this section are all unaffected —
> only the mechanism attributed to them moved. The consequences for the
> flight mass are a separate, larger thread: `TIER2_MASS_SCENARIOS.md`.

| quantity | arm H (h405, capped) | arm L (lin a 27.5) |
|---|---|---|
| dKE₁/dE_bind [eV/eV] | **−0.5421** | **−0.5486** |
| max fit residual [eV] | 0.00008 | 0.00007 |
| dKE₂/dE_bind | −0.5513 | −0.5619 |
| d⟨KE⟩(n2–8, geo)/dE_bind | −0.5617 | −0.6077 |
| d⟨KE⟩(n10–17)/dE_bind | −0.4171 | −0.6813 |
| dtrap/dE_bind [1/eV] | **+0.6199** | **+0.0326** |
| KE₁ at E_bind = 0 | 0.6876 | 0.9189 |
| KE₁ at the bundle well | 0.6236 | 0.8546 |

Fits are over the five provenance wells only (0.048 / 0.071 / 0.113 /
0.1168 / 0.154); 0.0 and 0.2168 eV are labelled out-of-provenance
diagnostics (`in_provenance = 0`) and enter no fit.

**Four things this establishes.**

1. **The response is affine over the whole span, not just locally.**
   Residuals against the provenance line stay ≤ 0.0007 eV out to
   *both* diagnostic wells — 0.0 and 0.2168 eV, a 4.5× range in depth.
   There is no saturation and no curvature: the refund is a **fixed
   fraction**, not a depth-dependent one.
2. **It is drag-form invariant.** 0.5421 vs 0.5486 — 1.2 % apart across
   two genuinely different laws (Φ 0.64 lin vs the h405 plateau). A
   refund driven by drag *magnitude* would not do that.

   > **Mechanism REVISED THREE TIMES. The final form is §9.6, not this
   > block.** Form-invariance is not evidence about the drag refund at
   > all: it is the signature of a term that does not involve the drag
   > law — the `m(1)/m(21)` mass partition, identical on both arms by
   > construction. Once that is removed, the genuine drag refund *is*
   > form-split as this block's reasoning expected, 1.9 % (capped) vs
   > 6.2 % (linear), a 3.3× ratio in the predicted direction — it was
   > simply 30× smaller than the effect it was being read from. The
   > product-of-two-factors physics below is **correct as physics** and is
   > what governs the 1.9 %/6.2 % residual and the R = 9 Å calibration
   > end; it was mis-scaled, not mis-derived. Retained for that reason.
   >
   > From the energy balance
   > `c = 1 + d(∫F_drag ds)/dδ`, the refund is a **product of two
   > necessary factors**:
   >
   > 1. **overlap** — the extra deceleration must occur where drag is
   >    live (if the well is paid outside the helium, the ion is only
   >    slowed downstream and `∫F ds` is untouched);
   > 2. **velocity sensitivity** — the drag force must respond to being
   >    slowed (a saturated, v-independent force does the same work over
   >    the same path however slow the ion is).
   >
   > Either factor at zero kills the refund. This is why a centre-born
   > ion under the h405 cap measures c = 1.000 flat across R = 9–47.8 Å
   > (factor 2 = 0, so varying factor 1 does nothing), while under the
   > **uncapped** Tier-0 cubic c ≈ 0.48–0.73 and **sharpening the density
   > halves the refund** (factor 1 is the manipulable one). The two
   > measured facts (fixed fraction, form-invariance) stand throughout;
   > only the explanation moved.
3. **The trap lever is form-SPLIT — the genuinely new number.** 0.620/eV
   on capped_cubic (inside the D0 §9 band 0.85 ± , measured 0.6199) but
   **0.033/eV on the linear arm**: over the entire 0 → 0.2168 eV span the
   lin family moves trap 0.0003 → 0.0085. The lin family is
   ~19× less well-sensitive on the trapped channel — the ring's "trap
   axis dead" single-point read, now a slope over seven wells.
4. **The axis is closed as a KE lever.** Even at `E_bind = 0` — the
   physically impossible limit — h405 reaches KE₁ 0.6876, short of the
   0.75 (A)-ceiling band and far short of the 1.00 eV reference peak. The
   whole well axis is worth **≤ +0.064 eV** on KE₁.

**Cross-instrument check.** The committed lin MD ring measures the same
coefficient on CRN well pairs: −0.5268 (lr1→lr2), −0.5062 (lr3→lr4),
mean **−0.5165** vs the twin's −0.5486 on the same chord — the twin runs
≈ 6 % steep there. **MD-confirmed on the capped arm (§9.3): −0.5308 vs
the twin's −0.5421, a 2 % transfer.**

**W₁ side effect (twin, non-gating).** On the capped arm W₁ is
*non-monotone* in the well: 1.388 (0) → 0.703 (bundle) → **0.600**
(0.154) → 0.676 (0.2168) — the 0.154 cell sits below the G4 W₁ floor
≈ 0.67, but it **fails the twin n̄ gate** (4.247 < 4.4), i.e. it buys W₁
by pushing n̄ out of band; the known trade, not a new landing. The lin
arm runs the other way (0.665 → 0.734, monotone worsening) and stays
gated at every well. Gate verdicts here are the §3.5c **twin-convention**
band (n₁ [0.19, 0.30], n̄ [4.4, 7.1]); the ring recorded that this differs
from the MD-side realization band.

**Pre-registered predictions (frozen before execution).** EB-P1 **PASS**
on its threshold (|slope_H| 0.542 < 0.56) but the *mechanism* claim
behind it — that capped_cubic would refund visibly more — is only
marginally supported (1.2 %); read the result as form-invariance, not as
a form effect. EB-P2 **PASS** (residual 8e-5 ≪ 0.005). EB-P3 **FAIL,
both arms** — the size grading is real and correctly signed (heavier
clusters pay closer to the nominal toll: +0.020 H, +0.059 L) but 2.5–7×
smaller than the predicted 0.10–0.20 band. EB-P4 **PASS** capped /
**FAIL** lin — the failure *is* finding 3. EB-P5 **PASS**. EB-P6 **PASS
on sign only** (+0.0007 eV excess at E_bind = 0): the predicted
saturation is absent, which is finding 1.

**Provenance caveat, on every row.** E_bind is *Derived* — jointly
extracted with the drag pair (§6.5.1). Every cell here overrides it
alone and so deliberately breaks that pairing (the §6.7 item-2
precedent). These are sensitivity reads, never candidate points; §9.1
bounds the physically defensible variation at ≤ 0.009 eV, so the
0.2168 eV diagnostic is ~14× outside what any geometry argument licenses.

**Retires:** the KE half of this section's GAP. **Leaves open:** the
steepness-decoupling test of the geometric-overlap explanation. (The
capped-arm MD confirmation named here as open was executed same-day —
§9.3.)

### 9.3 MD confirmation of the transfer function — EXECUTED (atlas §6.5 Step 3, 2026-08-10, 1 × N = 500)

One MD cell (`ebmdh405s`: h405 pins, well 0.0482 eV), CRN-paired against
the **committed** `linrh405p` ring run — cfg diff is exactly
`{binding_energy_I_ion_eV, allow_unvalidated_binding_pairing}`, so the
pair differs in the well and nothing else. Partner oracle green (nine
committed `atlas_linring_table.csv` observables reproduced).
Artifact `atlas_ebind_md.csv`.

| observable | bundle 0.11676 | shallow 0.0482 | MD slope [eV/eV] | twin |
|---|---|---|---|---|
| ⟨KE⟩ n=1 | 0.6371 | 0.6735 | **−0.5308** | −0.5421 |
| ⟨KE⟩ n=2 | 0.5511 | 0.5848 | −0.4914 | −0.5513 |
| ⟨KE⟩ n2–8 (geo) | 0.3205 | 0.3581 | −0.5493 | −0.5617 |
| ⟨KE⟩ n10–17 | 0.0789 | 0.1098 | −0.4521 (3 bins — fragile) | −0.4171 |
| trap | 0.0630 | 0.0120 | **+0.7439 /eV** | +0.6199 |
| n̄ | 3.748 | 4.245 | −7.25 /eV | −7.36 /eV |

**The transfer function is confirmed.** MD −0.5308 vs twin −0.5421 —
**2 %**, better than the lin arm's 6 %. The twin's error is a near-constant
**level** offset on KE₁ (−0.0135 eV at the bundle well, −0.0127 at the
shallow one), not a slope error: the *slope* is the licensed quantity.
The n̄ lever agrees to 1.5 % across instruments and matches the §6.7
item-2 lq value (−0.50 per 0.0686 eV step) — three independent systems.
The trap lever lands at 0.744/eV, inside D0 §9's band and between the
twin's 0.620 and the §6.7 0.85.

**Pre-registered MD-P1..P6:** P1 (refund exists) **PASS**; P2 (twin
transfer, band [−0.55, −0.47]) **PASS**; P3 (rigid translation, KE2
within 0.10 of KE1) **PASS** at 0.039; P4 (trap ≤ 0.025) **PASS** at
0.0120; P6 (ceiling holds) **PASS** at 0.6735 — extrapolating the MD
slope to `E_bind = 0` gives 0.699, still below 0.75. **P5 FAILED**, and
it corrects a §9.2 claim — see below.

**CORRECTION to §9.2 (the twin-grading caveat is withdrawn as stated).**
§9.2 recorded, from the lin arm alone, that "the twin under-reports
cross-bin KE grading by ~2.5×". MD-P5 predicted that would reproduce on
the capped arm (band [0.05, 0.25]). It does not: MD grading
(mid-band − n₁) is **+0.0184** against the twin's **+0.0277** — the twin
slightly *over*-reports here, and both are far below the lin arm's
values (MD +0.13, twin +0.059). The defensible statement is narrower:
**cross-bin KE grading is small and drag-form-dependent** (≈ 0.02 on the
capped arm, ≈ 0.06–0.13 on the linear arm); the twin tracks it to ±0.01
on the capped arm and misses it by ≈ 0.07 on the linear arm. It is *not*
a systematic twin under-report, and §9.2's EB-P3 band [0.10, 0.20] was
mis-calibrated by generalizing from one lin measurement. Slopes transfer
across instruments; differences of slopes are form-specific and should
be measured per arm, not carried over.

### 9.4 RQ12 probe — the refund is a *velocity-sensitivity* effect, and the Tier-0 calibration is NOT biased by it (2026-08-10, zero MD)

Built to test a specific hypothesis (user): if the model cashes only ~55 %
of the well, then Tier-0's co-fit of `E_bind` — anchored to real TDDFT
traces at R = 9/18 Å — must have returned an *inflated* `E_bind` ≈
`E_bind^true`/c, and correcting the density width would deflate it and
raise production KE₁. The transfer factor is `T = c(R_prod)/c(R_cal)`;
`T = 1` would make the whole refund an unobservable re-parametrisation.

Instrument: twin stage `refundscan` (+ a byte-inert `rho_steepness`
override on `integrate_pairs`, the twin-side mirror of production's
`erf_independent` gate). Controlled single trajectories — monodisperse R,
fully dressed fixed mass, radial launch, two wells, per-fragment finite
difference; no fate map, no ensemble. Artifacts
`atlas_refund_geometry{,_summary}.csv`.

**Correction applied during execution (not a refinement).** The raw
finite difference conflated the refund with a birth-geometry term: at
small R the ion is born partway *up* the well, so it can never pay the
full depth. That term alone reproduced the raw c to four decimals at
R = 9 and 18 Å (0.7774 predicted vs 0.7771 measured). All values below
are divided by the available fraction `ρ̂(depth_birth)` at the potential
width; the first-pass numbers were discarded.

Birth-corrected cashed fraction c (centre-born arm):

| s_ρ [Å] | R 9 | R 18 | R 26.6 | R 34.4 | R 47.8 | R 68.3 |
|---|---|---|---|---|---|---|
| **14.2** (standing) | **1.000** | **1.000** | 0.535 | 0.556 | 0.671 | −0.088 |
| 7.1 | 1.000 | 1.000 | 0.652 | 0.661 | 0.739 | −0.114 |
| 3.14 (Harms DFT) | 1.000 | 1.000 | 0.703 | 0.700 | 0.762 | −0.123 |

(The `frac027` birth arm runs 1.061 / 0.880 / 0.798 / 0.612 / 0.200 /
0.249 at s_ρ = 14.2 — see the birth-sensitivity point below.)

**Second correction, same session — the drag LAW.** The first working
build measured both sides of the transfer question under a *capped* law,
and hardcoded `G3_STANDING[1]` = v_c **7.25** (the superseded chord)
rather than h405's 5.5. Both are wrong for this question: Method B
extracted {a, b, E_bind} under the **uncapped** `shared_pure_cubic`,
while production runs the h405 cap. Above a cap the drag force is
v-independent, so a refund cannot exist *by construction* — which is
where a spurious `c = 1.000` at the calibration radii came from. The
probe now carries both laws explicitly (`REFUND_LAWS`), and the transfer
factor is **cross-law by construction**:

$$T = \frac{c_\text{prod}(\text{h405 cap},\, R \approx 48)}
          {c_\text{cal}(\text{uncapped cubic},\, R = 9)}$$

Birth-corrected c under the **Tier-0 (uncapped) law**, centre-born:

| s_ρ [Å] | R 9 | R 18 | R 26.6 | R 34.4 |
|---|---|---|---|---|
| **14.2** (standing) | **0.482** | 0.494 | 0.553 | 0.612 |
| 7.1 | 0.663 | 0.625 | 0.659 | 0.700 |
| 3.14 (Harms DFT) | **0.754** | 0.681 | 0.700 | 0.732 |

**1. The mechanism — a PRODUCT of two necessary factors.** From
`c = 1 + d(∫F_drag ds)/dδ`, a refund requires **both** (i) that the extra
deceleration happen where drag is live — *overlap*, the factor the
density width controls — and (ii) that the drag force respond to being
slowed — *velocity sensitivity*. Either at zero kills it. Under the h405
cap a centre-born ion stays above v_c from 9 to 47.8 Å and measures
c = 1.000 flat: factor (ii) is zero, so varying (i) does nothing. Under
the uncapped cubic, (ii) is live everywhere and c ≈ 0.48–0.73, and
**sharpening the density halves the refund** (0.518 → 0.246 at R = 9)
— factor (i) is the manipulable one. In the strongly over-dissipated
corner (R = 68.3, or R = 47.8 uncapped + sharpened) c goes **negative**:
a deeper well slows the ion enough that it loses *less* to cubic drag
than it gained in toll — a > 100 % refund. Physical, not a failure.

*(Two earlier framings in §9.2 — "geometric overlap" alone, then
"velocity sensitivity, not overlap" — were each wrong as stated; the
second was a false dichotomy. Corrected there too.)*

**2. Sharpening makes E_bind's effect STRONGER.** c rises toward 1 at
every radius as s_ρ sharpens (0.482 → 0.754 at R = 9; 0.553 → 0.700 at
26.6). Less refund, more of the well actually paid — the predicted
direction (RF-P2, which holds everywhere off the over-dissipated cell).

**3. E_bind IS inflated by the refund — and the transfer error runs
toward more KE, not less.** With c_cal = 0.482, a fit anchored to the
TDDFT traces returns `E_bind^fit ≈ E_bind^true / 0.482 ≈ 2.1 ×
E_bind^true`. Measured transfer factors:

| arm | T (standing) | T (sharpened) |
|---|---|---|
| centre | 2.073 | 1.327 |
| production-like birth | 1.488 | 1.032 |

**T > 1: production cashes more than the calibration did, i.e. it
over-pays the exit toll.** Correcting the density width lowers the
effective toll and **raises** KE₁ — and RF-P3 **PASSES** on both arms
(sharpening drives T toward 1, removing the transfer error).

**Magnitude — T is NOW a number (§9.6, 2026-08-11): T = 2.04,
ensemble-confirmed.** The apparent disagreement between the probe's
production-side c (0.80–1.00) and the MD ensemble 0.5308 was **not** an
ensemble-vs-trajectory discrepancy — it was the mass frame. Scored in the
frame the toll is actually paid in (dressed, m(21)), the production
ensemble measures **c_prod = 0.9815** on 6427 ions, i.e. it agrees with
the single-trajectory probe's 1.000 to 2 %. Both terms of

$$T = \frac{c_\text{prod}}{c_\text{cal}} = \frac{0.9815}{0.482} = 2.04$$

are now fixed-mass trajectory-level quantities, so the category error that
sank the first estimate cannot recur. This lands on the probe's own
centre-arm value 2.073. **`E_bind^fit ≈ 2× E_bind^true` stands** — the
robust part was always `T > 1`, and it is now sized.

The "+0.006–0.02 eV" payoff estimate is **reinstated at ≈ +0.02 eV**, but
for different and now level-consistent reasons: it is entirely
**calibration-side** arithmetic. Sharpening lifts c_cal 0.482 → 0.754 at
R = 9 Å, so a re-fit returns `E_bind` ≈ 0.1168 × (0.482/0.754) ≈ 0.075 eV
(≈ 0.085 via R = 18), and propagating that through the **MD-confirmed**
system slope −0.531 gives **ΔKE₁ ≈ +0.017…+0.022 eV**. Caveat: this
assumes the re-fit rescales `E_bind` alone, whereas Method B fits
{a, b, E_bind} jointly — an order-of-magnitude forecast, not a
prediction. **The production-side leg is dead** (§9.6): at R ≈ 48 Å the
drag refund is 1.9 %, so no width correction can act there.

**4. What IS falsified — RF-P5, and only that.** RF-P5 predicted
c(18)/c(9) ≈ 2.17, the ratio of the Tier-0 per-case fitted E_binds
(0.154 at 9 Å / 0.071 at 18 Å), if that split were the refund. Measured
**1.023** (centre) / 0.994 (frac027) under the correct Tier-0 law: c is
essentially *equal* at the two calibration radii, so the refund cannot
produce a 2.2× difference in fitted E_bind. **The 0.154/0.071 split has
another, still-unknown cause** — a live loose end.

**5. Birth position matters as much as radius.** Under the h405 cap at
R = 47.8, c = 1.000 centre-born vs 0.800 at the production-median birth
fraction; under the Tier-0 law, 0.706 vs 0.164. The refund is not a
geometric constant — it is a strong function of birth depth through the
speed history, coupling it to the §14.3 birth-depth lever.

**Verdicts.** RF-P1 **PASS**; RF-P2 **PASS** (production-like arm; the
centre arm's FAIL is the over-dissipated R = 68.3 cell only); RF-P3
**PASS** both arms; RF-P5 **FAIL** — the headline. RF-P4 was withdrawn as
mis-specified: a single trajectory cannot reproduce an ensemble mean, and
the two birth arms *bracket* it instead.

**Consequence for RQ12.** Sharpening the density width is right physics
**and** correctly signed for KE₁ — it removes a real
calibration-to-production transfer error. Its size is not yet established
(see the magnitude caveat above). It is not a solution to the KE₁
deficit. Doing it properly requires a Tier-0 re-extraction under the
corrected width, and the production-side c must first be re-measured on
the ensemble rather than on single trajectories.

### 9.5 Is the measured response dynamical or cascade? — DECOMPOSED (2026-08-11, zero MD, zero new integrations)

The §9.4 probe measures c = 1.000 for a supercritical single trajectory
while the ensemble measures 0.5422 (twin) / 0.5308 (MD). Since the probe
runs at **fixed mass with no fate map**, the hypothesis was that much of
the ensemble "refund" is not a drag refund at all but **cascade
repopulation** of the n = 1 bin — which the density width could not
reach, making the whole RQ12 leverage argument mis-attributed.

Method: KE₁ depends on the well through two bundles — **T** (trajectory:
`v_inf`, `trapped`) and **C** (cascade: `n_det` from the fate map). Score
all four cross-combinations KE₁(T_x, C_y) and difference. The split is
exact by construction; the ordering-dependence is reported as the
interaction term. Stage `ebinddecomp`, artifact `atlas_ebind_decomp.csv`;
both un-crossed corners **anchored** against the committed
`atlas_ebind_twin.csv` bin means before any crossed value was read.
Because KE₁ conditions on n = 1, the complex mass is pinned at m(1) — the
cascade enters **only** through bin membership, not through mass.

| cell | c_total | c_traj (dynamical) | c_casc (cascade) | cascade share | interaction |
|---|---|---|---|---|---|
| **h405** | 0.5422 | **+0.6090** | −0.0669 | **12.3 %** | −0.0000 |
| **lr1** (lin) | 0.5488 | **+0.5817** | −0.0330 | 6.0 % | +0.0001 |

**The hypothesis is REFUTED — DC-P1 FAIL** (registered: cascade ≥ 0.15;
measured 0.067 / 0.033). In refund terms the total refund 0.458 splits
**0.391 dynamical (85 %) / 0.067 cascade (15 %)**. So the response *is*
overwhelmingly the drag channel, and **the density width can reach it** —
RQ12's leverage argument is not mis-attributed on this axis.

**It also resolves the 1.000-vs-0.53 tension, and not as guessed.** The
ensemble's *dynamical* value is 0.6090 (DC-P2 PASS: a real refund even
with membership frozen), so the gap from the single trajectory's 1.000
down to 0.542 is ≈ 0.39 **ensemble geometry** (birth-position spread and
non-radial paths putting real ions below v_c during the crossing) and
only ≈ 0.07 cascade. The probe's defect was representativeness, not the
missing mass mechanism.

> **THIS ATTRIBUTION IS WITHDRAWN (§9.6, 2026-08-11).** The 0.39 is not
> ensemble geometry. It is the **mass frame**: `refundscan` scores at the
> same fixed dressed mass it integrates with (so m_ratio = 1 → c = 1.000),
> while the ensemble scores KE₁ at m(n_det) = m(1). Re-scored in the
> dressed frame the ensemble gives **c_traj = 0.9815**, i.e.
> 0.6090 = 0.9815 × 0.6205 exactly. Birth spread and non-radial paths
> contribute the 1.9 % residual, not 0.39. The **DC-P1 verdict is
> unaffected** — the cascade really is small — but the "85 % dynamical /
> 15 % cascade" split is better stated as *38 % mass partition,
> 6.7 % cascade, 1.9 % drag*.
>
> **The open T1 caveat below is also retired, and not by measurement:**
> `trapped` cannot contaminate c_traj through the n = 1 bin, because a
> trapped ion parks inside the droplet for the full 150 ps, accumulates a
> large cooling exposure K, and therefore lands at n_det = 21 — never
> n = 1. The `& ~trapped` mask is near-inert for KE₁. Splitting trapping
> into its own bundle is **not worth running**; the split that was
> actually missing was the mass frame.

**Cascade sign.** c_casc is small and **negative**: deepening the well
shifts n = 1 membership toward slightly *faster* ions, partially
offsetting the dynamical loss.

**Caveat on attribution (open).** `trapped` was bundled with T, being a
trajectory output, so **c_traj = velocity response + trapping
selection**. Trap moves 0.012 → 0.063 across this well pair, so part of
the 0.609 is selection rather than pure dynamics. Splitting trapping into
its own bundle is one further crossing and should be done **before the
85 % figure is leaned on**. DC-P3 (closure) and DC-P4 (interaction
≤ 0.10) both PASS — the split is well-conditioned.

### 9.6 THE 45 % IS SOLVED — mass partition, not a refund; and the density width is measured NULL at production (2026-08-11, zero MD)

This section supersedes the mechanism statements in §9.2, §9.4 and §9.5.
None of the *measured numbers* in those sections change; the physics
attributed to them does. **There is no defect and no missing energy: the
model is working correctly.**

#### 9.6.1 The question

Every earlier reading said the ion cashes only ~55 % of the solvation
well and that ~45 % is "refunded as un-incurred drag". That was always
uncomfortable — an ion cannot escape a droplet while paying only half
its binding energy — and three measured invariances made it untenable:
the "refund" was invariant to the drag **form** (0.5421 capped vs 0.5486
linear, §9.2), invariant to birth **depth** (§9.2), and — measured here —
invariant to the density **width**. A quantity indifferent to everything
about the drag is not a drag effect.

#### 9.6.2 The answer: a mass-frame partition

$$\boxed{\;c \;=\; \frac{m(1)}{m(21)} \;=\; \frac{130.903}{210.955} \;=\; 0.6205\;}$$

The twin integrates every trajectory at the **dressed** mass
`complex_mass_amu(ne_mol)`, and at the corrected geometry *every* fragment
carries `n₀ = 21` (measured: 20000/20000; the shallowest birth in the
committed master sits at depth −19.95 Å, so the pickup gate is saturated).
`_bin1_mean_ke` then scores `kinetic_energy_eV(complex_mass_amu(n_det),
v_inf)` — the **bare** mass. Hence

$$\Delta\!\left(\tfrac12 m_{21} v^2\right) = -E_\text{bind}
\;\Longrightarrow\;
\Delta KE_1 = \tfrac12 m_1 \Delta(v^2) = -E_\text{bind}\cdot\frac{m_1}{m_{21}}$$

**The toll is paid in full. The observable only ever sees the iodine's
share of it.** The other 38 % was paid by helium that evaporates
afterwards and is no longer part of the scored particle — under the
Tier-1a velocity-preserving shed each departing atom leaves at the
complex velocity and carries its own ½m_He v² away. Nothing is lost and
nothing is double-counted.

Confirmed by re-scoring the committed chord families in both frames
(both un-crossed corners anchored string-exact against
`atlas_ebind_twin.csv` before any crossed value was read):

| arm | frame | c_total | c_traj | c_casc | toll actually paid |
|---|---|---|---|---|---|
| **h405** (capped) | bare (production) | +0.5422 | +0.6090 | −0.0669 | |
| **h405** | **dressed** | +0.8737 | **+0.9815** | −0.1078 | **98.1 %** |
| **lr1** (lin) | bare (production) | +0.5488 | +0.5817 | −0.0330 | |
| **lr1** | **dressed** | +0.8843 | **+0.9375** | −0.0532 | **93.8 %** |

#### 9.6.3 Independent confirmation by energy closure

A completely different route — the ledger of the 6375 scored n = 1 ions
at h405 — gives the same answer:

| term | eV | |
|---|---|---|
| source budget (14.3996 / 2.666 / 2) | 2.7006 | per fragment |
| − solvation toll actually paid | **0.1167** | **99.9 % of E_bind = 0.1168** |
| − drag dissipation | 1.5790 | 58 % of the budget |
| = exit KE, dressed frame | 1.0050 | |
| × m(1)/m(21) | 0.6205 | |
| = scored KE₁ | **0.6236** | matches the committed value |

The well line is the direct answer to "where is that energy coming from":
**nowhere — it is paid, 99.9 % of it.**

#### 9.6.4 The corrected decomposition

| component | h405 (capped) | lr1 (lin) |
|---|---|---|
| mass partition `m(1)/m(21)` | 38.0 % | 38.0 % |
| cascade (bin membership) | 6.7 % | 3.3 % |
| **genuine drag refund** | **1.9 %** | **6.2 %** |
| = observed shortfall | 45.8 % | 45.1 % |

The genuine drag refund is **form-split 3.3×** in the direction §9.2's
mechanism block predicted — it was simply buried under a constant 30×
larger. It also explains why n-graded slopes rise (|slope| KE₁ < KE₂ <
mid on both arms, since m(n)/m(21) rises with n) and **retro-explains the
EB-P3 FAIL**: the mass ratio predicts a grading increment ≈ 0.057, not the
registered 0.10–0.20 band. That band was mis-specified, not the model.

#### 9.6.5 The density-width scan — atlas §6.8 T2, EXECUTED (NULL)

`c` re-measured at five widths across the Harms bracket, both arms,
two-well finite difference (licensed by EB-P2 affinity). Anchors: the
re-integration at s_ρ = 14.2 is **bit-identical** to the committed chord
npz, and all four (arm × well) `n1_ke_eV` reproduce
`atlas_ebind_twin.csv` **string-exact**.

| s_ρ [Å] | 14.2 (standing) | 7.10 | 4.40 | 3.50 | 3.14 (Harms DFT) |
|---|---|---|---|---|---|
| h405 (capped) | 0.5422 | 0.5548 | 0.5583 | 0.5590 | **0.5597** |
| lr1 (lin) | 0.5488 | 0.5635 | 0.5686 | 0.5693 | **0.5695** |

**A 4.5× sharpening moves `c` by 3 %**, converged by s_ρ = 4.4. In KE₁
terms **+0.002 eV**, against a 0.36 eV deficit. This is exactly the size
§9.6.4 predicts (the whole recoverable drag refund is 1.9 %), so the null
is *structural*, not bad luck — two independent routes to the same
number.

Why: at R ≈ 48 Å, ∫ρ̂ dr = R to 0.0 % (§17), so sharpening only
redistributes exposure across the surface shell; and a deep-born ion
(median birth depth 33.8 Å) is already at terminal speed when it crosses
that shell, so there is no velocity gradient for the redistribution to
exploit. The width bites only where **s ∼ R** — i.e. at the 9 Å Tier-0
calibration droplet (c_cal 0.482 → 0.754), and nowhere in production.

Side observables, capped arm only (twin, non-gating, hatch reads): n̄
4.484 → 4.322, W₁ 0.7025 → **0.6318**, trap 0.0535 → 0.0637, midHot
1.0685 → 1.0557. The lin arm is flat (n̄ 4.507 → 4.507). The W₁ movement
is below the G4 floor ≈ 0.67 and is worth noting, but n̄ 4.32 falls out of
the twin gate band [4.4, 7.1], i.e. it buys W₁ the known way. **Mechanism
not verified** — the plausible route is the cooling exposure
`K = ∫ρ dt` losing its ~30 Å of post-surface accumulation, which raises
`E_ej = E₀e^(−K)` and evaporates more. One integration would settle it;
not run.

#### 9.6.6 Consequences

- **RQ12's production-side leg is retired.** Sharpening the density width
  is right physics and remains a model-correctness fix, but it cannot
  move production observables: there is nothing there to recover.
- **RQ12's calibration-side leg is intact and now sized**: T = 2.04,
  ΔKE₁ ≈ +0.02 eV via a Tier-0 re-extraction (§9.4, T5).
- **T1 is not worth running** (see the §9.5 note); **T3's four-consumer
  confound is moot on the twin side** — pickup is measured inert and the
  twin has no detection stage, so "change the sharpness" is one argument
  to `integrate_pairs`.
- **A structural ceiling on every source-side lever**, new:

  $$\frac{\partial KE_1}{\partial(\text{any pre-evaporation energy})}
    \;\le\; \frac{m(1)}{m(21)} = 0.6205$$

  §3.5k's measured `S_k = 0.389 eV/eV` is **63 % of that ceiling**, not a
  free number that happened to land low; and the top **32 %** of its
  registered band [0.35, 0.75] was structurally unreachable. The BP-KILL
  verdict is unaffected (0.389 sits inside, above 0.35).
- **The large consequence is not in this section.** The same mass frame
  says KE₁ = E_exit · m(1)/m_flight with E_exit measured *mass-invariant*,
  which turns the entire KE₁ deficit into a question about what mass
  actually flies. That thread — §3.5g cross-check, the flight-mass
  measurement, the KE-vs-n shape, and the speed-gated co-moving
  hypothesis — is **`TIER2_MASS_SCENARIOS.md`**.

**Instruments:** scratchpad probes (`rq12_width_slope.py`,
`rq12_mass_frame.py`, `rq12_g35g_crosscheck.py`, `rq12_flight_mass.py`),
each anchored on committed artifacts before any new number was read. No
committed stage was modified and no committed artifact was rewritten.
Promotion of any of these to a committed stage is a separate decision.

**Retires:** the RQ12 magnitude GAP and the §9.2/§9.4/§9.5 mechanism
ambiguity. **Stance unchanged:** nothing adopted; Tier-0, h405 and
`finc1v725` all stand.

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

**Low-n KE is the same transit channel, measured (G4 Step 2, user-found
2026-07-28):** the corrected geometry **halves the n = 1 mean KE**
(pooled h405 battery 0.641 ± 0.003 eV vs the incumbent's 1.034 against
ref mean/median/mode 1.302/1.128/0.891) and takes ~20 % off n = 2
(0.549 ± 0.002 vs ref 0.706), while n ≥ 3 lands (h405 n = 3: 0.471 vs
0.496). The fastest, most-stripped fragments pay the full ~35 Å
transit-drag toll — the incumbent's n = 1 KE match was partly a
shallow-birth artifact, the KE-side twin of its W₁ 0.571. The scored
surface was blind (n = 1 KE only inside reported-not-scored χ²_med,
which did register 242 → 350+; midHot's n = 2–8 geo-mean dilutes n = 2);
now **pre-registered as a hard axis** for every future ranking (plan
§3.5f: KE₁ vs ref median 1.128, KE₂ vs ref mean 0.706,
|ln ratio|/ln(1.15) score terms). p_tail's "not a second lever" verdict
was earned on deep-KE and does **not** carry over to this axis.
*(Anchor amendment 2026-07-29, plan §3.5g: KE₁ is scored against the
experimental n = 1 **peak 1.00 eV** — the wide upper tail is read as a
second process; the median anchor stays reported.)*

**The low-n KE in-surface freedom, MEASURED (§3.5g retro-scan,
2026-07-29, zero MD, 52 committed rows — the per-knob KE₁ influence
GAP closes):**

- **∂KE₁/∂E₀ = −2.8…−3.9 eV/eV** (arm-resolved, CRN-paired), gate-bounded
  at **+0.038 eV** (the n₁ floor binds after ΔE₀ ≈ −0.014).
- **∂KE₁/∂τ = +0.083 eV/ps** along the gated ridge chord (the first τ
  read on this axis), gate-bounded at **+0.094 eV** (n₁ floor at
  τ ≈ 5.5).
- **∂KE₁/∂v_c = −0.43 eV per Å/ps**, paying **∂midHot/∂v_c = −1.48**:
  c50 (v_c 5.0) reaches KE₁ 0.956 / KE₂ 0.881 at midHot **1.923**, v525
  KE₁ 0.774 at midHot 1.394 — the KE₁↔midHot trade is measured across
  the whole v_c range, not a linearization artifact. In-gate bound
  **+0.054 eV**.
- **R5 verdict: the additive in-gate bound is +0.19 eV against the
  required +0.31 eV (target 0.95 = 1.00 − 1σ) — the (v_c, τ, E₀)
  surface CANNOT reach the KE₁ anchor.** The p_tail axis is motivated
  with this record as its evidence (the E₀-arm/W₁-floor structure).
- **The n = 1 channel is bifurcated (the R4 geometry read, standing
  chord):** the deep-born cells (Boltzmann/center, depth ≥ 25 Å) produce
  **zero** n = 1 fragments at every R (n₁ = 0.000, 6/6 cells; f725 at
  the corrected ensemble likewise n₁ = 0, n̄ 14.3, trap 0.437 —
  under-stripped, not over-stripped); **every** standing-chord n = 1 is
  shallow-born (depth 9–19.5 Å) with **KE₁ 1.042–1.091 eV** — the
  reference peak scale — and shallow n₁_solv falling with R
  (0.235/0.184/0.151/0.128 at R 26.6/34.0/49.4/68.3). The corrected
  basin instead makes its n = 1 deep-born and slow (0.64). The
  experimental peak ≈ 1.0 numerically coincides with the model's
  shallow-birth channel — the mixture/birth-law lever (a shallow-birth
  or small-droplet weight the pure center-weighted Boltzmann law lacks)
  is a *second* live candidate beside p_tail, and it is the only one of
  the two that also feeds the Block-D n = 1 *weight* deficit (−0.101).
  *(REJECTED on physics by the user 2026-07-29 — solvated I₂ sits
  centrally; the R4 bifurcation re-reads as evidence the p_tail = −1
  tail over-drags real deep-born fragments.)*

**p_tail single-knob, MEASURED AND KILLED (§3.5h ring, 2026-07-29,
3 × N = 1000 CRN-paired at the h405 pins):**

- **The KE₁ placement works:** pt15 (−1.5) 0.917 ∈ [0.79, 0.99], pt20
  (−2.0) 1.151 ∈ [0.96, 1.16]; pt30 (−3.0) 1.428 over-band (the 1-D
  model breaks at far softening). ∂KE₁/∂p_tail ≈ −0.26 eV per unit
  near −1.5.
- **But the tail is NOT an exit-toll-only knob — it re-times the whole
  cascade.** At pt15/pt20/pt30: midHot 1.92 / 3.02 / 4.41, deepKE
  0.88 / 4.32 / 12.5, n̄_det 3.08 / 2.39 / 1.84, supp 0.26 / 0.34 /
  0.41, n₁ 0.259 / 0.295 / 0.346, trap 0.011 / 0 / 0, χ²_med up to
  4.6e4. The 5.5–9.5 Å/ps band is where mid/deep ions do their late
  deceleration (and where the E2 Landau-gated relaxation thermalizes
  them): relief there keeps the *entire* ensemble fast and extends
  stripping. **PT-P4 sign prediction refuted** (n̄ falls, n₁/supp rise
  — the residence/relaxation channel dominates the heating channel),
  **PT-P5 orthogonality refuted** (deep bins are the most affected,
  not the least).
- **PT-P3 kill criterion FIRED** (every cell with KE₁ ≥ 0.95 has
  midHot ≫ 1.15) → per the frozen registration, **no further
  single-knob tail cells**; the axis moves to a joint re-tune or the
  honest-residual branch (user).
- **The selectivity numbers any successor shape must beat:** required
  relief n = 1 +0.36 eV vs n = 2 ≤ +0.16 vs mid ≲ +0.10 — a ≳ 2.3:1
  differential between exit speeds 9.9 and 9.1 Å/ps. Uniform tail
  softening measured a ≈ 1.3:1 differential (pt15: KE₁ +0.28, KE₂
  +0.30 — KE₂ overshoots ref before KE₁ reaches target). A
  KE₁-selective lever must confine relief to **v ≳ 9.5 Å/ps** (the
  band only the n = 1 exiters occupy at late times, once-through for
  everyone else) — a band-limited tail (second cap v_c2 with softening
  only above it) is the shape this measurement points at; **not
  registered, user adjudication pending**.

### 14.4 Twin authority at the corrected geometry (G3 Step 1, 2026-07-27 — the scan instrument's error model; **ring-validated at Step 3**: n̄ bias 0.24–2.66 He small-end residence-scaled at 14/14 cells, n₁ transfer ≤ 0.036, trap floor +0.012…+0.13; **ranking authority MEASURED at G4 Block 0, 2026-07-28: W₁ ρ +0.82 and midHot ρ +1.00 LICENSED, deepKE ρ +0.33 NOT — the deep-KE axis is not twin-scannable**)

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

**Ranking authority — MEASURED (G4 Block 0, 2026-07-28, zero MD).** The
14 paired (twin, MD) ring cells were read as a transfer measurement for
the first time (report `tier2atlas_g4_transfer.py`, CSV
`h2b_g4_transfer.csv`; both §1.4 oracles passed first). Spearman ρ,
twin vs MD, with 95 % bootstrap CI, against the pre-registered
permission gate ρ ≥ 0.7:

| observable | ρ (n = 14) | 95 % CI | verdict |
|---|---|---|---|
| W₁_solv | **+0.824** | [+0.44, +0.98] | **LICENSED** — the twin may *rank* on W₁ |
| midHot | **+0.996** | [+0.93, +1.00] | **LICENSED** — near-perfect rank transfer |
| deepKE | +0.327 | [−0.28, +0.82] | **NOT LICENSED** — gate-only; MD must rank it |

Consequences, both first-order for the atlas:

- The §14.4 header caveats "W₁ bias-loaded" and "KE direction-only" are
  now **split**: W₁ and midHot are *rank*-faithful (the bias is a level
  shift, not a reordering), while **deepKE is not** — the twin's deep-bin
  ratio does not order MD cells. **The deep-KE axis is not twin-scannable
  at all**; any RQ11 statement must come from MD. This is measured, not
  assumed, and it is the reason the G4 MD finalists must span the ridge
  rather than sit at the twin's top-scoring cell.
- **n̄ bias model** (replacing the crude [+0.3, +3] bracket):
  `Δn̄ = a + b·n̄_twin`, Δn̄ = MD − twin, `a = +0.58, b = −0.238`
  (13 cells, f725 excluded as leverage) or `a = +0.23, b = −0.173`
  (all 14). **Residual SD 0.28 He either way** — the fit's R² is carried
  by f725 but the *band* is robust, so predicted-MD n̄ is good to ±0.3 He
  and the twin gate can be run at the MD band itself.

**Level transfer is a REGRESSION, not an offset (same session).** Ranking
licence says nothing about level, and the score is a *distance* from the
experiment, so the level model matters for cell selection. Over the 13
non-f725 paired cells:

- `MD_W₁ = 0.671 + 0.412 · twin_W₁` (R² 0.60, resid SD 0.18),
- `MD_midHot = 0.001 + 0.897 · twin_midHot` (R² **0.994**, resid SD 0.035).

A constant shift is refuted: the mean W₁ shift over all 13 cells is +0.03 but
over the four MD-*gated* cells it is **+0.32 ± 0.10**, and it is negative at
the badly-landing cells (d036 −0.62) — the twin over-disperses W₁ and the
transfer regresses toward the middle. Two consequences:

1. **The W₁ intercept 0.671 is a predicted floor.** Even a twin cell at
   W₁ = 0 maps to MD W₁ ≈ 0.67, *above* the standing pooled 0.571. Extrapolated
   from cells that mostly land badly (R² 0.60), so it is a lead, not a verdict
   — but it is the quantitative form of "the corrected geometry may not reach
   the old histogram quality".
2. **midHot must be ranked on the corrected value.** |ln x| is V-shaped about
   1, so a multiplicative bias does *not* preserve the ordering of the raw
   values even though the rank correlation is +0.996: a twin midHot of 1.00
   lands at 0.90 in MD, while a twin 1.07 lands at 0.96. Ranking cells on raw
   twin midHot systematically picks the wrong ones.

**Seed-SD normalizers (same session, five N = 1000 battery members).**
W₁ 0.5791 ± 0.0954, midHot 1.0108 ± 0.0222, deepKE 0.6327 ± 0.0209,
χ²_med 131.5 ± 15.4. Two **provenance corrections to plan §1.2**
(recorded, not silently overwritten):

1. §1.2's "W₁_solv SD ≈ 0.04" is the **SEM of the pooled mean**
   (0.0954/√5 = 0.043), not the per-seed SD. The per-seed N = 1000
   scatter is 0.095, so an N = 500 single-seed W₁ carries ≈ 0.13 — the
   G3-ring W₁ column is even softer than it was read as.
2. §1.2's "deep-bin KE read 0.0603 ± 0.0033" does **not** reconcile with
   the deepKE ratio (0.633 ± 0.021). It is not the ratio's mean/SD; the
   most likely reading is the deep-bin *population weight*. Flagged, not
   resolved; the G4 score normalizes by the measured 0.0209.

**KE₁ joins the authority box — ranking LICENSED (free-form linear
sweep Step 0 `ke1auth`, 2026-07-30, zero MD;
`TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md` §2/§2.1,
`atlas_ke1_authority.csv`).** Twin KE₁ (the `g3_score` n = 1 bin mean)
vs MD KE1_mean over the 26 distinct-pin corrected-ensemble cells of the
committed `atlas_ke_lown_scan.csv` (dedupe: pooled N = 5000 > finals
N = 1000 > ring N = 500; f725 excluded — MD n₁ = 0): **Spearman ρ
+0.9979** (KE₂ +0.9966), against the pre-registered bands 0.8/0.5 —
far above the licensure line, the twin's best KE bin (deepKE stays
+0.33 NOT licensed; the n = 1 class is short-chord and
mechanism-light, the opposite of the deep bins). Levels are *not*
licensed (the standing rule): twin KE₁ runs uniformly cold, ratio
0.897–0.984 (median 0.975), additive −0.013…−0.037 eV (median
−0.017); 3/25 adjacent rank inversions, all among near-ties.
Transfer caveat pre-registered: measured on capped-cubic cells;
transfer to any new form family is an assumption backstopped by that
family's MD ring.

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

**REFINED — the fine ridge scan (G4 Block 1, 2026-07-28, zero MD; stage
`g4scan`, 5544 cells, CSVs `h2b_g4scan_{predictions,gated_ke,ridge}.csv`).**
v_c step 0.25, τ step 0.4, E₀ step 0.005, E_bind unchanged; the Step-2
grids are exact sub-lattices and **G4-P1 passed — 408 shared cells
reproduce the committed Step-2 values string-exactly**.

- **574/5544 cells gate (172 clean of both caveat stamps), against 24/6480
  at Step 2.** The dominant cause is *not* the finer grid: it is the gate.
  Step 2 had to carry the crude bias bracket (twin n̄ ∈ [4.4, 7.1]); the
  Block-0 bias model gates on predicted MD n̄ ∈ [3.77, 4.37] ± 0.28, i.e.
  twin n̄ ∈ [3.8, 5.4] — *narrower*, but positioned where the n₁ band
  actually lives. n₁ and n̄ are anti-correlated under stripping, so the old
  n̄ floor of 4.4 was fighting the n₁ band. **The binding constraint at
  Step 2 was the instrument's error model, not the parameter surface.**
- **G4-P2 CONFIRMED — the ridge is connected.** The clean gated set forms
  one 4-neighbour-connected diagonal band from (v_c 5.5, τ 4.8) to
  (v_c 6.0, τ 6.4): at τ ≤ 4.4 only v_c ≤ 5.5 gates, at τ ≥ 5.2 the band
  reaches v_c 6.25. The two Step-2 "sub-basins" were opposite corners of a
  single ridge, exactly as hypothesized.
- **Clean optimum: v_c 5.5, τ 4.8–5.2, E₀ 0.34–0.37** (E₀ now resolved to
  0.005). a037's corner is confirmed near-optimal but slightly off — the
  best clean cells are (5.5, 4.8, **0.365**) and (5.5, 5.2, **0.34**).
- **No clean cell is predicted to beat the incumbent.** On the licensed
  axes with the transfer regression applied, best clean S_pred **1.59**
  (per-cell uncertainty ≈ 0.57) vs the pooled battery's measured
  **1.091**; best of all, including the off-bundle-well stamp, is
  (5.5, eb0482, 5.6, 0.33) at 1.47. In-sample check at the four ring
  cells: predicted S 1.82/6.89/2.66/1.79 vs measured 1.90/6.83/3.52/2.05.
- **G4-P3 is NOT EVALUABLE at twin level.** 75 gated cells hold
  midHot ∈ [0.85, 1.15] ∧ twin deepKE ≥ 0.6, but deepKE is precisely the
  axis Block 0 measured as rank-untransferable (ρ +0.33) — a037's twin
  deepKE 0.755 landed at MD 0.43. The pre-registered Block-2 (p_tail)
  trigger therefore **cannot be decided by this scan** and moves to the MD
  finalists. This is a pre-registration defect discovered by G4's own
  Block 0, recorded rather than reinterpreted.
- **Open confound on the W₁ comparison:** the ring cells are N = 500
  (≈ 900 scored ions) while the incumbent reference is pooled N = 5000.
  Sampling noise inflates W₁, so part of the 0.83–1.17 vs 0.571 gap may be
  statistical. The five battery members (N = 1000, W₁ 0.48–0.73, mean
  0.579) bracket the pooled 0.571, suggesting the N = 1000 → 5000 step is
  small; the N = 500 → 1000 step is untested. The mandatory a037/b031
  N = 1000 replicates in Block 3 are the control that settles it.
  **SETTLED — it is not sample size (below).**

**MD-CONFIRMED at N = 1000 (G4 Block 3, 2026-07-28; 6 cells, fresh shared
seed 20260729, table `atlas_g4finals_table.csv`).** Oracles passed
(frozen twin rows string-exact, pooled-battery drift).

| cell | (v_c, τ, E₀) | n̄ | n₁ | W₁ | midHot | deepKE | χ²_med | gate |
|---|---|---|---|---|---|---|---|---|
| f1 | 5.5, 4.8, 0.365 | 4.169 | 0.1878 | 0.821 | 0.990 | 0.523 | 411 | 0 |
| f2 | 5.5, 5.2, 0.34 | 4.194 | 0.1876 | 0.916 | 0.980 | 0.438 | 440 | 0 |
| **f3** | 5.5, 4.4, 0.395 | 4.205 | 0.1867 | **0.769** | **1.023** | **0.614** | 426 | 0 |
| **a037** | 5.5, 4.8, 0.37 | 4.024 | 0.1984 | 0.792 | 0.960 | 0.526 | 360 | **1** |
| b031 | 6.0, 6.4, 0.31 | 4.263 | 0.1771 | 1.113 | 0.520 | 0.732 | 68 | 0 |
| x345 | 5.5, 4.8, 0.345 | 4.782 | 0.1491 | 0.966 | 1.155 | 0.499 | 603 | 0 |

1. **The n̄ bias model is quantitative.** Predicted − measured MD n̄ =
   −0.005 / +0.000 / +0.004 at f1/f2/f3 (fitted residual SD 0.28 He).
   Bias-corrected, the twin is an n̄ instrument, not just a ranker.
2. **n₁ carries transfer scatter of ±0.02 with inconsistent sign, and it
   decides the gate.** f1/f2/f3 land n₁ 0.1867–0.1878 — missing the 0.19
   floor by 0.001–0.003 — because MD came in 0.011–0.015 *below* twin,
   whereas at the G3 ring MD ran *above* twin (+0.001…+0.026). Block 1
   bias-corrected n̄ but treated n₁ as unbiased on the strength of its
   ≤ 0.036 tolerance; that tolerance is scatter, not zero bias. **Any
   future twin gate on n₁ must carry a ±0.02 band, and a ridge optimum
   selected on twin n₁ sits ≈ half an E₀ step too low.**
3. **The W₁ deficit is REAL, not sampling noise (G4F-P3).** a037
   N = 500 → 1000: W₁ 0.826 → 0.792 (Δ −0.034, against a per-seed SD of
   0.095); b031 moved the other way (1.061 → 1.113). The corrected
   geometry sits at **W₁ ≈ 0.77–0.92 against the incumbent's 0.579**, and
   the measured floor agrees with the regression intercept 0.671.
4. **The mid-vs-deep KE tension is BROKEN (G4F-P4).** Twin deepKE
   0.691/0.776/0.887 → MD 0.438/0.523/0.614 at f2/f1/f3, ordering
   preserved, span 0.176 across cells that are near-identical in n̄/n₁.
   **f3 holds midHot 1.023 *and* deepKE 0.614 simultaneously** — the
   incumbent's own KE profile (1.011 / 0.633). So the deep-KE axis *is*
   parameter-accessible on the (v_c, τ, E₀) surface at fixed drag form;
   the G3-ring reading "no cell holds both axes" was a resolution artifact
   of its coarse E₀ ladder. **Block 2 (p_tail) is not required by this
   evidence.** Note the twin ranked f3 *worst* of the three; MD ranks it
   best — the unlicensed-deepKE finding in operation.
5. **a037 replicates its landing** on a fresh seed at double N (n̄ 4.024,
   n₁ 0.1984) — the corrected-geometry basin is seed- and N-robust.
6. **b031's corner is out:** trap 0.203, midHot 0.520, S 6.63. The
   v_c 6.0 / τ 6.4 sub-basin does not survive N = 1000.
7. **Fate/trap at the ridge:** trap 0.088–0.093 with the marginal
   (modelling-exclusion) class at **0.002–0.003** — the retained-policy
   question is empirically empty at the successor point (G4 input).

**Net:** at the corrected geometry the residual deficit is **one axis —
W₁ — not a KE trade-off**. Best S on the licensed axes: f1 1.511 ≈ f3
1.512 < a037 1.683, all above the incumbent's 1.091, and only a037 gates.
The measured n₁ bias plus f3's KE profile point the successor at
**E₀ ≈ 0.38–0.41 on the (v_c 5.5, τ 4.4–4.8) chord** — untested.

**E₀ LADDER — the successor located and the W₁ floor measured (G4 Block 3
ladder arm, 2026-07-28; 9 further cells × N = 1000, same seed 20260729,
CRN-paired with the six above).** Three τ arms at v_c 5.5 plus one v_c
probe. The E₀ direction is the atlas's cleanest one-knob response:

| τ | E₀ | n̄ | n₁ | W₁ | midHot | deepKE | gate |
|---|---|---|---|---|---|---|---|
| 4.4 | 0.395 (f3) | 4.205 | 0.1867 | 0.769 | 1.023 | 0.614 | 0 |
| 4.4 | **0.405 (h405)** | 3.955 | 0.2078 | **0.709** | 0.957 | 0.603 | **1** |
| 4.4 | 0.410 (h410) | 3.862 | 0.2119 | 0.674 | 0.931 | 0.570 | 1 |
| 4.4 | 0.415 (h415) | 3.753 | 0.2147 | 0.669 | 0.906 | 0.548 | 0 |
| 4.8 | 0.370 (a037) | 4.024 | 0.1984 | 0.792 | 0.960 | 0.526 | 1 |
| 4.8 | 0.375 (h375) | 3.889 | 0.2141 | 0.747 | 0.926 | 0.585 | 1 |
| 4.8 | 0.380 (h380) | 3.769 | 0.2202 | 0.733 | 0.895 | 0.577 | 0 |
| 5.2 | 0.345–0.355 | 4.04→3.73 | 0.194→0.222 | 0.892→0.872 | 0.945→0.859 | 0.43→0.50 | 1/1/0 |

1. **E₀ influence (the knob's own entry, §4):** inside the ridge, raising
   E₀ by 0.005 eV moves **W₁ −0.02…−0.04, n₁ +0.004…+0.007, n̄
   −0.11…−0.13 He, midHot −0.02…−0.03**, monotonically in every arm. It
   is a single-parameter trade of histogram quality against n̄ and midHot,
   and each τ arm gates over a window of only ~2 E₀ steps (0.01 eV).
2. **The W₁ floor is real and ≈ 0.67.** The τ 4.4 arm flattens:
   0.769 → 0.709 → 0.674 → **0.669**, against the incumbent's 0.579.
   Block 0's transfer regression predicted the MD floor as its intercept
   **0.671** — an independent estimate agreeing to 0.002. The corrected
   geometry cannot reach the standing point's histogram quality by any
   (v_c, τ, E₀) setting.
3. **Why — n₁ and n̄ cannot be matched simultaneously.** W₁ improves
   because n₁ climbs toward the experimental 0.243, but n̄ falls as it
   does. Extrapolating the measured τ 4.4 slopes, n₁ = 0.243 needs
   E₀ ≈ 0.435, where n̄ ≈ 3.3 — **0.77 He below the experimental 4.07**.
   The standing point achieved n₁ 0.243 *and* n̄ 4.068 together only at
   the wrong geometry. **This is the precise, quantified form of the
   corrected-geometry deficit: the detected size distribution has the
   wrong shape, not the wrong scale** — which localizes the residual to
   the mechanism (ladder shape, pickup, per-shed ε), exactly the branch
   §3.5c's failure criterion pointed at, and not to the drag surface.
4. **τ ordering is the twin's reverse.** On W₁ at matched gate status:
   τ 4.4 (0.669–0.709) < τ 4.8 (0.733–0.792) < τ 5.2 (0.872–0.892). The
   twin ranked τ 5.2 *best* (twin W₁ 0.42–0.49). **Twin W₁ is not usable
   for fine discrimination** — a second inversion after f1/f3, despite
   ρ = 0.824. Read: the licensing was measured over cells spanning W₁
   0.5–12 and does not survive down to separations of ~0.05.
5. **Deep-KE ceiling (v525 probe, v_c 5.25/τ4.4/E₀0.39):** deepKE
   **0.656 > the incumbent's 0.633**, and the cell *gates* (n₁ 0.208,
   n̄ 3.994, W₁ 0.727) — but midHot 1.394 and χ²_med 912. So the deep-KE
   axis does have headroom above the standing point, reachable by
   lowering v_c, and it is paid for in midHot. RQ11 has a ceiling answer:
   ~0.66 with midHot ≈ 1.4, or ~0.60 with midHot ≈ 0.96.
6. **n̄ model validity range:** |pred − MD| ≤ 0.08 He everywhere inside
   the ridge (systematically −0.04…−0.08 at the high-E₀ end of each arm),
   but −0.30 at b031 and +0.11 at v525/x345 — the model is an
   in-ridge instrument, not a global one.
7. **Successor candidate: h405** (v_c 5.5, τ 4.4, E₀ 0.405) — gated,
   best gated S 1.559, and better than a037 on **every** axis (W₁ 0.709
   vs 0.792, deepKE 0.603 vs 0.526, midHot 0.957 vs 0.960, n₁ 0.208 vs
   0.198). Its supp 0.183 also lands on the standing 0.187 (RQ3), its
   trap is 0.087 and its marginal class 0.002.

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
| **`potential_steepness` 14.2 Å reused as the He *density* width** | DFT fit of the **solvation potential** (legacy: `% from ernesto dft result beta = [14.3324 26.9916 34.4431]`, "from fit of solvation potential DFT result") — then reused verbatim by `rho_he_ratio` via `drag_gate_steepness` | **mismatch**: 14.2 Å ⇒ a 10–90 % interface of **25.7 Å**, vs **5.7 Å** (DFT) / **6–8 Å** (experiment) for ⁴He droplets at N = 10³–10⁴ [Harms/Toennies/Dalfovo PRB 58, 3341 (1998)] ⇒ $s_\rho \approx 3.1$–4.4 Å, a factor ≈ 4. The potential is $\rho\otimes V_\text{I–He}$ + cavity, so ≈ 14 Å is *right for U* and too wide for ρ — the defect is the **reuse**, not the value | **convention (C)**, quantitatively bounded below | **RQ12** — the DFT *density* profile from the same calculation. Instrument already exists and is bit-inert: `drag_spatial_gate="erf_independent"` + `cfg.drag_gate_steepness` (DESIGN §5.7 G3, held for exactly this) |
| ↳ *its measured leverage* (2026-08-10/11, §9.4–§9.5) | (a) analytic: for a radial exit $\int_0^\infty\hat\rho(r-R)dr = R$, so a **saturated** drag is width-blind on total dissipation (0.0 % at R = 47.8 Å, +11 % at the 9 Å Tier-0 droplet); (b) the **E_bind transfer factor** $T = c_\text{prod}/c_\text{cal}$, measured cross-law; (c) the dynamical/cascade split of the ensemble response | (a) leaves the drag channel nearly inert on *total* dissipation; (b) $c_\text{cal} = 0.482$ under the **uncapped** Tier-0 cubic ⇒ the co-fitted $E_\text{bind}$ is inflated ≈ 2.1×, and **T > 1** (1.66–2.07 taken consistently at trajectory level): production **over**-pays. Sharpening drives T → 1 and **halves the refund** (0.518 → 0.246 at R = 9); (c) the ensemble response is **85 % dynamical / 15 % cascade** (§9.5) — so the density width *can* reach it | — | **consequence:** correcting $s_\rho$ is right physics and correctly signed for KE₁, but its **size is not yet established** (the "+0.006–0.02 eV" estimate mixed trajectory- and system-level derivatives — withdrawn). Needs: the trapping bundle split out of $c_\text{traj}$, an **ensemble** $c(s_\rho)$, then a **Tier-0 re-extraction**. Third channel, untouched: **birth dressing** ($\hat\rho(10\,\text{Å})$ 0.84 → ~1.00, $n_0 \approx 18 \to 21$), the §14.3 lever |
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

## 18. Drag state coupling s(n) — MEASURED DEAD in-window: gate-clipped (probe 2026-07-29)

| | |
|---|---|
| knob | `drag_state_coupling` ∈ {off, shell_area} + (R_core [Å], ρ_shell [Å⁻³]) |
| class | Bounded geometric closure (design `TIER2_DRAG_STATE_COUPLING_DESIGN.md`) |
| instrument | 3-cell probe sa22/sa30/sa44 (ρ_shell bulk/prior/2×bulk, R_core 3.2) × N = 1000 at the h405 pins, seed 20260729 CRN-paired; oracles O1/O2 + cfg-diff + unit oracle all pass |
| status | **probe EXECUTED; every registered prediction REFUTED; axis STOPPED** (SC-P2 signature failure kills the central claim per the registration). Code delivered and regression-clean; production default stays `off`. |

**Measured influence (vs the CRN h405 baseline KE₁ 0.637, needle SD 0.038):**

- KE₁ 0.596 / 0.599 / 0.601 at sa22/sa30/sa44 — **sign-inverted** (−0.04,
  the design predicted +0.2..+0.5) and **ordering-inverted** (strongest
  coupling = lowest KE₁).
- **Needle SD unmoved**: 0.038 in every cell (SC-P2 floor was 0.08) — the
  n–KE lock survives the state coupling entirely.
- Cells nearly indistinguishable across the full ρ Bounded range
  (factor 2): the swept parameter has almost no in-window authority.
- Exploratory (SC-P5): n̄ +0.26..+0.30, n₁ −0.012..−0.017, supp −0.014,
  trap +0.005, W₁ +0.046..+0.062 — all cells strictly worse than baseline
  on the landing surface; deepKE −0.08..−0.11.

**Why (the mechanism finding — measured on sa30's ion checkpoint, zero MD):**
the coupling is **gate-clipped**. Stripping is slow relative to transit:
**no ion reaches n ≤ 8 while still inside the droplet** (fraction 0.0000;
n ≤ 14 inside: 7.4 %), the mean shell at first surface crossing is
**19.0 He (median 20)** — ions exit essentially fully dressed — and 93 %
have exited by sim-end. All stripping to the design's large-|effect|
regime (n ≤ 2, s ≈ 0.32–0.43) happens outside, where the spatial gate
g(d) = 0 has already switched the drag off. In-window the shell stays in
[14, 21] → s ∈ [0.85, 1.06], and there s is nearly ρ-independent by
construction (ρ differentiates s only where R_core³ competes with n/ρ,
i.e. at low n). The small negative net KE₁ comes from the early-transit
s(21) ≈ 1.06 over-drag at the highest speeds plus back-reaction (less
mid-window drag → n̄ rises, the §14.5 residence channel sign). The
needle cannot break because per-ion in-window s variance is minute —
the strip-timing decorrelation the design §2 posited does not exist
inside the gate.

**Consequence:** being at low n and being inside the drag gate are
mutually exclusive states of the delivered mechanism — *no* γ(v,n)
coupling of this family (state factor × gated γ) can move KE₁, because
the exit toll is paid dressed. Any successor low-n KE lever must either
act outside the gate (a different force surface entirely) or change the
strip-timing itself (the OQ-F cooling-contact discussion, or upstream
E_int dynamics). Options on the table per the registration: the
honest-residual branch extended to KE₁, or OQ-F — user adjudication
either way. Nothing adopted; `finc1v725` and the h405 successor
candidacy stand as before the probe.

Records: findings "§3.5i s(n) probe" (scorer table + verdicts), scorer
CSV `atlas_sn_probe.csv`, design doc §8 (registration), log entries
2026-07-29.

---

## 19. Coulomb budget (E_coulomb_scale / per-fragment KER) — MEASURED ALIVE: S_k = 0.389 eV/eV; suppression is E_int-driven (probe 2026-07-29)

| | |
|---|---|
| knob | `E_coulomb_scale` (pair Coulomb ×scale at fixed R₀ 2.666) + `coulomb_available_eV` stamp; `internal_energy_partition_fraction` used as the E₀ pin |
| class | **Sourced** (user paper input 2026-07-29: gas phase reproduced at 0.8·E_C with Q = 2/Q = 3 → channels 2.16 / 4.32 eV per I⁺; RQ7/RQ8 NBs). The standing production 2.70 (scale 1.0) is ~25 % high vs this calibration. |
| instrument | §3.5k probe bud226k/bud411k/bud411 × N = 1000 at the h405 pins, seed 20260729 CRN-paired; oracles O1/O2 + cfg-diff + kinematic unit oracle all pass |
| status | **probe EXECUTED; BP-KILL NOT fired (S_k = 0.389 ≥ 0.2) — the source-side lever is ALIVE; the (C) design discussion proceeds per the registration** |

**Measured influence (vs the CRN h405 baseline, budget 2.7006):**

- **Kinematic transfer slope S_k = dKE₁/dbudget = 0.389 eV/eV**
  (E₀-pinned line 2.26 → 2.7006 → 4.11: KE₁ 0.486 → 0.637 → 1.200;
  BP-P2 band [0.35, 0.75] CONFIRMED). The 1-D placement slope 0.587
  over-predicts by the back-reaction factor ≈ 0.66 — this also puts
  bud226k 0.006 eV above its BP-P1 band edge (0.486 vs [0.29, 0.48]);
  the registered diagnose-obligation is discharged by the S_k reading
  itself (the placement's constant-toll picture is ~⅓ too steep).
  bud411k lands in-band (1.200 ∈ [1.13, 1.70]). KE₂ tracks (0.404 →
  0.551 → 1.081).
- **The suppression mechanism, measured (BP-P3 sign-INVERTED — the
  probe's mechanism yield):** supp(bud411) 0.909 > supp(bud411k) 0.546
  at identical kinematics — raising E₀ 0.405 → 0.6165 adds +36 pp
  suppressed (∂supp/∂E₀ ≈ +1.7 /eV). Suppression is **E_int-driven
  (self-unbound at detection)**, not freeze-driven: more deposited or
  onset energy → hotter E_int → more of the ensemble relabeled bare.
  Under the full (B) proportional wiring the fast ions go to
  suppressed-bare and **KE₁ stays at baseline** (bud411: KE₁ 0.623,
  n = 1 count 44) — the Q3 → n = 1-upper-half leg therefore requires
  (A) exit stripping (rebin at existing KE) or partial E₀ decoupling.
  The (A)/(B) complementarity is now measured from both sides.
- **Needle persists** (BP-P4 ✓): KE₁ SD 0.034 / 0.047 at the k-cells —
  a single-valued budget cannot widen the n = 1 KE; width requires the
  channel mixture.
- **Uniform-budget landing damage (exploratory — NOT the (C)
  channel-weighted forecast):** bud226k n̄ 5.33 / n₁ 0.126 / trap 0.179
  (2×) / W₁ 1.365 — a naive re-anchor to the calibrated 2.16–2.26
  **breaks the standing basin** (re-tune required, RQ7 NB); bud411k
  n̄ 1.63 / supp 0.546 / midHot 2.301 / deepKE 1.093 / W₁ 1.294;
  bud411 n̄ 0.34 / supp 0.909 / trap → 0.001.

**Consequence:** the KER → KE₁ transfer is physically real in-model at
≈ 0.39 eV/eV, the honest-residual branch is NOT forced, and the (C)
design (paper-anchored channel mixture {single, Q2, Q3} + exit
stripping) proceeds with this section as its slope authority. Per-channel
E₀ coupling (f_int channel-dependence) is a live design axis with
measured leverage. Nothing adopted; `finc1v725` and the h405 candidacy
stand.

**Addendum (2026-07-29, post-§3.5l):** the 0.8 calibration's
provenance is **CLOSED** — Hatherly et al, J. Phys. B 27 (1994)
2993–3003 (doi:10.1088/0953-4075/27/14/032; consulted, not
repo-kept): gas-phase finite-pulse CE, fraction
channel-independent (0.75/0.65·E_C at 90/200 fs), intrinsic channel
widths (FWHM 2.6/5.8 eV at 200 fs); NOT droplet screening; the Abel
gas export confirms 2.26/4.11 in the repo's own frame (a 2-D-rim
"frame split" rider was posted and withdrawn same-session). Under the
adjudicated (C) design (`TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`
OQ-B) **this knob's scalar form retires**: the per-molecule sampled
channel budget E_m (means 2.16/4.32 × shared Bounded f, widths
sampled) replaces the single `coulomb_available_eV` stamp; the
uniform-re-anchor damage row above stays the measured warning against
any single-budget re-anchor.

Records: findings "§3.5k budget probe" (table + verdicts), scorer CSV
`atlas_budget_probe.csv`, plan §3.5k (registration), RQ7/RQ8 NBs
2026-07-29, log entries 2026-07-29; provenance + (C)-design records:
findings "§3.5l reads (i)+(iii)", `TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`.

---

## 20. CE channel mixture (B) + depth-graded exit strip (A) — MEASURED at the (C) probe (2026-07-29): registration FAILED, all three kills fired; the strip arm over-tolls (ε-dominated) and FEEDS the suppressed gate; the mixture alone PLACES KE₁

**Status: probe-measured (4 × N = 1000 at the h405 pins, seed 20260729,
CRN; `atlas_ce_probe.csv`; findings "(C) probe EXECUTED"). Knobs built
2026-07-29 (config surface `ce_*` / `exit_strip_*`, checkpoint v8, both
enums default off = bit-identical). Frozen inputs and CP bands: design
doc §8. NOTHING adopted; the §11 PC-3 failure semantics are triggered
and await user adjudication.**

Measured influence (deltas vs the committed h405 baseline row
trap 0.087 / supp 0.183 / n̄ 3.955 / n₁ 0.208 / W₁ 0.709 / midHot 0.957 /
KE₁ 0.637 / KE₁ SD 0.038):

| knob / arm | measured influence |
|---|---|
| `ce_channel_mode = sampled` alone (bonly) | **KE₁ 0.920** (± the 1.00 anchor; above-1.15 share 0.42 ∈ the CP-3 band; KE₁ SD 0.286 — needle broken): the source-side placement works as forecast. **But** trap 0.087 → 0.465 (the E_single ≈ 0.53 channel, s_m ≈ 0.2, cannot climb out of realistic droplets — a class P3's registered single-channel assumption never priced), n̄ 6.69, n₁ 0.080, W₁ 2.97 (CP-7 range). supp 0.138 (mild w_Q3-diluted parking, not the full §3.5k complementarity). |
| `exit_strip_mode = depth_graded` alone (aonly) | **Catastrophic over-toll.** 90 % of ions strip, knock count mean 13.7 / p90 20: the box-mid ε_carry 0.025 × the measured knock count = 0.34–0.50 eV, plus Σ-outer 0.14–0.20 eV → total toll ≈ 0.5–0.7 eV per exiter, 2.5–3.5× the design's "≈ 0.2 eV" anchor. KE₁ 0.637 → **0.119**, midHot 0.96 → 0.30, W₁ 0.709 → 0.947, trap 0.087 → 0.167 (toll converts marginal escapers to retained). **supp 0.183 → 0.551**: the strip lowers Σ(n) faster than E_int drains — post-strip survivors at n 2–12 with E_int ≈ 0.15 eV sit above the tiny post-strip Σ(n) and park suppressed forever. **(A) feeds the suppressed scaffolding instead of retiring it** — the §17 ledger entry inverts. |
| C-full (both on) | The two failures compose: trap 0.564, supp 0.368, KE₁ 0.180, W₁ 1.887, midHot 0.370, slow-bare 0.446. CP-1..4 FAIL, CP-5/6/7 kills FIRED. |
| `f_int,Q3` 0.0985 → 0.15 (coupling arm, in-mixture) | Nearly flat: midHot 0.370 → 0.350, supp 0.368 → 0.361, KE₁ 0.180 → 0.213 — the coupling is second-order behind the strip toll; **CP-6 fired at both bracket ends**. |

Structural reads (why the P1/P3 forecasts missed):

1. **ε-box internal inconsistency, now measured.** ε_carry ∈ [0, 0.05]
   was anchored on "total ≈ 0.2 eV" assuming few knocks; the measured
   knock count (median ~14, p90 20 at the pinned box: P₀ ≈ 0.25–1 for
   ordinary exiters × G ≈ 1 on all rungs above j₀ ≈ 2) makes the box
   mid alone worth 0.35–0.50 eV. Consistency with the 0.2 eV anchor at
   the measured knock counts needs ε ≲ 0.005 eV/He.
2. **P1's evidence base stripped the suppressed class only** (the CF-3
   counterfactual re-labels), with no energy toll; the live §3.3 strip
   acts on EVERY outbound crossing of EVERY ion — the general solvated
   population that carried the good histogram is stripped too (W₁/midHot
   kills), and multi-crossing (OQ-J) ratchets the count to p90 = 20.
3. **The E_int/Σ(n) gate inversion.** CP-1's "conversion by
   construction" was occupancy arithmetic; live, a stripped survivor
   keeps its E_int while Σ(n) collapses → gate-closed → suppressed
   (aonly: 918 suppressed at n_det 2–12, E_int ≈ 0.15 eV). Emptying the
   class needs the strip to touch E_int (or the E_int at crossing to be
   below the post-strip Σ) — a design change, not a knob setting.

Class/identifiability: the (C) knobs stay Bounded-with-anchor as
registered; no value inside the frozen boxes can undo mechanisms 1–3
(the kills are family/structure-level, per the design §11 PC-3
pre-commitment). The scalar `coulomb_available_eV` budget did NOT
retire — it stands with `finc1v725`/h405 (nothing adopted).

**ADJUDICATED CLOSED (2026-07-30, user).** The PC-3 consequences are
ratified ((A) strip design v1 stopped; the current f_int wiring
closed) and **all three candidate follow-ups are DECLINED** (ε → ~0 +
suppressed-class-gated strip; E_int co-strip; single-channel
trap-class pricing) — the (C) line is SHELVED as a measured boundary,
re-openable only on new external evidence, not a re-parameterization.
Structural reading recorded at closure (arithmetic estimate from the
committed probe-table scored counts, NOT a per-ion read — the v8
`ce_channel`/`ce_E_m_eV` fields support an exact zero-MD
decomposition if ever needed): under the covariance-anchored weights
the n = 1 feeder pool is essentially the partner-masked Q3 channel
(~200 of ~1800 scored-eligible ions per 1000 molecules ≈ 0.11 of the
scored ensemble, a hard ceiling); even converting every bonly
suppressed ion to n = 1 gives n₁ ≈ (78 + 134)/970 ≈ 0.22 < 0.31 —
**the anchored mixture cannot populate n = 1 at the reference level**,
and re-splitting the provisional w_single/w_Q2 pair drains the trap
class without feeding n = 1 (Q2 at 2.16 eV lands dressed). The n = 1
residual therefore has exactly three live attributions: (1)
exit-conversion physics (the strip family — v1 measured dead here);
(2) the soft anchors (w_single/w_Q2 — measurably insufficient for
n₁); (3) the reference n = 1 bin partly not owned by the modeled
droplet ensemble (small-droplet / surface / gas-side feed — the
unmeasured branch). Record: findings "(C) adjudication CLOSED", log
2026-07-30.

Records: findings "(C) probe EXECUTED" (full table + verdicts +
decomposition), scorer CSV `atlas_ce_probe.csv`, generator
`gen_tier2atlas_ce_probe.py`, scorer `tier2atlas_ce_probe_table.py`,
design doc §8/§9/§11, log entries 2026-07-29 + 2026-07-30.

---

## Cross-references

- Archive/provenance: `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee,
  I1–I100, boundaries §6, open questions §7, run inventory §8).
- Program that fills the GAPs: `TIER2_SENSITIVITY_ATLAS_PLAN.md`.
- Calibration classes + identifiability: `CALIBRATION_MAP.md`.
- Open physics questions: `RESEARCH_QUESTIONS.md` (RQ1–RQ11 + NB
  registers).
- Decision history: `drag_migration_log_tier2.md`.
