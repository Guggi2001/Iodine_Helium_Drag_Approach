L# Tier 2 — Sensitivity Atlas: Findings

> Created 2026-07-23 at first execution (stage 2a). Protocol, yardstick,
> and stance: `TIER2_SENSITIVITY_ATLAS_PLAN.md` (§1, §8) — reported, not
> adjudicated; `finc1v725` stays the standing point. Findings are recorded
> parameter→influence-style; chronology lives in
> `drag_migration_log_tier2.md`. Axes A (droplet geometry), A/D2b
> (sampling laws), and B (E₀/τ) have not started — their sections are
> created when their stages fire.

## D4 Step 1 — Method-B form table + trace-tail re-inspection (2026-07-23, zero MD)

Executed per plan §6.2: **existing Tier-0 fit artifacts reused verbatim,
no new fits run.** Machinery: `scripts/extraction/atlas_d4_step1_report.py`
(table + γ(v) + tails) and `scripts/extraction/atlas_d4_step1_oracle.py`
(machinery verification; both read-only w.r.t. `data/reference/`). Dual
purpose: the RQ11 lever-hierarchy step-1 reading is recorded as
**NB-RQ11-10** in `RESEARCH_QUESTIONS.md`.

### Oracle (machinery verification)

- **Part A — point evaluation: PASS, bit-exact.** The joint objective at
  the stored `shared_pure_cubic` point reproduces the committed record
  exactly: objective 0.1292649398514104 Å/ps, RMSE₁₈ 0.13440214672133738,
  RMSE₉ 0.12412773298148343 (drag ion path RNG-free, neutral stage
  seeded — bitwise determinism confirmed on today's environment).
- **Part B — full refit: PASS, bit-exact.** The complete §9 Stage-2
  joint Nelder-Mead re-run (Method-A anchors a₀ = 14.5556 /
  b₀ = 2.0534, maxfev 400/start, 161 joint evaluations — the same count
  as the committed run, 420 s) reproduces the stored minimum
  identically: b = 2.5153508541760052, E_bind = 0.11675778353879479 eV,
  objective 0.1292649398514104 Å/ps. The Method-B machinery is verified
  end-to-end (optimizer path included); every artifact-derived number
  above rides on a reproducible instrument.

### Per-form × per-mode table (Tier-0 artifacts, reused)

Incumbent: `shared_pure_cubic`, objective 0.12926 Å/ps; equivalence band
±0.005; held-out band ≤ 0.45 Å/ps. Tier-0 verdicts stand as
**pre-registered priors** (plan §6.1).

**Shared joint fit** (both traces, one coefficient set + E_bind):

| form | coefficients | E_bind [eV] | obj [Å/ps] | Δobj | standing verdict |
|---|---|---|---|---|---|
| pure cubic | b = 2.5154 | 0.1168 | 0.12926 | — | **locked incumbent** |
| linear+cubic (a free) | a → 0 (2.5·10⁻⁴), b = 2.5159 | 0.1168 | 0.12927 | +0.0000 | a unidentifiable → pure cubic (T_a0) |
| linear+quadratic | a → 0 (10⁻⁴), c = 12.79 | 0.0482 | 0.16315 | **+0.0339** | **REJECTED** (outside ±0.005) |
| power law | C = 2.835, n̂ = 2.927 | 0.1126 | 0.12917 | −0.0001 | equivalent; n ≈ 3 unforced |

**Held-out** (fit 18 Å alone → predict untouched 9 Å; band ≤ 0.45 Å/ps):

| form | 9 Å prediction RMSE | pass | note |
|---|---|---|---|
| linear+cubic (§8 bundle) | 0.2685 | ✓ | the designed anti-circularity test |
| linear+quadratic | 0.6987 | **✗ FAIL** | the trap-catcher result |
| power law | 0.4110 | ✓ (0.41 of 0.45 — marginal) | 18 Å-only fit rails n = 4 (bound) |

**Per-case single-curve** (diagnostic only, never preset-wired;
identifiability/transferability probe):

| case | form | coefficients | E_bind [eV] | obj [Å/ps] |
|---|---|---|---|---|
| 9 Å | linear+cubic | a = 7.50, b = 1.751 | 0.1543 | (§8 bundle) |
| 18 Å | linear+cubic | a = 1.456, b = 3.316 | 0.0709 | (§8 bundle) |
| 9 Å | linear+quadratic | a = 0, c = 10.32 | 0.1303 | 0.0446 |
| 18 Å | linear+quadratic | a = 16.63, c = 6.13 | 0.0532 | 0.0994 |
| 9 Å | power law | n = 2.6505 ± 0.0264 | 0.1566 | 0.0406 |
| 18 Å | power law | n = 4.0 (railed) | 0.0825 | 0.0929 |

Reads carried from Tier 0, now in one place: identifiability is
case-asymmetric (9 Å pins n, 18 Å rails); the per-case-vs-shared
coefficient gaps (cubic: b 1.75 vs 3.32; lq: c 10.3 vs 6.1) are the
shared-form-ansatz stress signal; the 9 Å per-case coefficients absorb
the Tier-0 Finding-2 transverse contamination. The E_bind spread across
these fits {0.048, 0.071, 0.113, 0.117, 0.154} eV is the §6.5 E_bind
scan grid.

### γ(v) diagnostic table (gate = 1, i.e. ρ̂ = 1; [amu/ps])

Band 2.54–4.95 Å/ps = TDDFT-calibrated; 7.25 = production cap; 10.5 =
production peak. Parked forms at locked b = 2.5154 with **implied**
parameters (v_s ≙ arbitrated cap; v_f from the §4ee twin window) — the
plan-§6.1 zero-cost entry gate, not fits.

| form | v=2 | 2.54 | 4.95 | 7.25 | 10.5 | F(10.5)/F(7.25) |
|---|---|---|---|---|---|---|
| capped_cubic v_c 7.25 (production) | 10.06 | 16.23 | 61.63 | 132.21 | 91.29 | **1.00** |
| pure cubic (uncapped) | 10.06 | 16.23 | 61.63 | 132.21 | 277.32 | 3.04 |
| lin+quad shared | 25.58 | 32.49 | 63.32 | 92.74 | 134.32 | 2.10 |
| power law shared | 10.78 | 17.09 | 61.81 | 128.94 | 263.22 | 2.96 |
| lin+cubic 9 Å per-case | 14.50 | 18.80 | 50.41 | 99.56 | 200.59 | 2.92 |
| lin+cubic 18 Å per-case | 14.72 | 22.85 | 82.71 | 175.76 | 367.05 | 3.02 |
| PARKED Padé v_s 7.25 | 9.85 | 15.56 | 46.75 | 66.11 | 68.68 | 1.50 |
| PARKED subtractive v_f 1.5 | 4.40 | 10.57 | 55.97 | 126.55 | 271.66 | 3.11 |
| PARKED subtractive v_f 2.0 | 0.00 | 6.17 | 51.57 | 122.15 | 267.26 | 3.17 |

Low-v limb vs the RQ11 deficit (drag at v ≈ 2 must go *down*, not up):
lq carries **×2.54** the production drag at v = 2 — wrong direction,
confirming the §6.1 framing; the subtractive forms bend the right way
(×0.44 / ×0.0); Padé and power law are cubic-like (×0.98–1.07).
Cap need: every polynomial/power form still grows ×2.1–3.2 in force
between cap and production peak — none removes the cap; only the Padé
saturates by construction (1.50 at implied v_s, → 1 asymptotically).

### Parked-form entry-gate verdicts (plan §6.1)

- **Saturating (Padé) cubic — gate does NOT fire for its intended
  purpose; stays parked.** At the implied v_s = 7.25 the form is *not*
  cubic-equivalent in-band: γ is −4 % at 2.54, **−24 % at the 4.95 band
  top, −50 % at the cap** — a fit would be pushed to larger v_s to
  restore in-band cubicity. But the deviation scales as (v/v_s)³, so
  holding it to few-% at 4.95 needs v_s ≳ ~15 Å/ps — at which the
  high-v saturation the form exists to provide is forfeited
  (F(10.5)/F(7.25) ≈ 2.6, vs the arbitrated 1.0). The order-3 Padé
  cannot simultaneously be cubic in-band and constant-force at the
  production peak: **the one-knob (v_s) replacement of (v_c, p_tail) is
  arithmetically over-constrained.** A Step-1 fit would only rail v_s
  upward (traces give a lower bound, tail stays Free) — no fit run.
  Tail authority still requires extended/production-kinematics traces
  (the standing collaborator ask).
- **Subtractive (gated) cubic — gate FIRES; promotion to a fit is a
  separate logged decision.** The suppression persists into the
  calibrated band (v_f = 1.5: ×0.65 at 2.54, ×0.91 at 4.95; v_f = 2.0:
  ×0.38 / ×0.84), and the **entire 18 Å window lives at
  v = 2.54–3.02 Å/ps** — exactly where the two candidates differ by
  ×0.65–0.38 from cubic — so the traces genuinely constrain v_f. The
  incumbent's tight in-window residual (RMSE₁₈ 0.134) is prima facie
  evidence *against* suppression that strong, but the joint fit's
  compensation freedom (b, E_bind) means only a re-fit measures it.
  Implementation note for that decision: in-window speeds never drop
  below ~2.5 Å/ps > v_f ≤ 2.0, so the zero-clamp never activates
  in-window — the existing `linear_cubic` engine with the a ≥ 0 bound
  relaxed (a = −b·v_f² < 0) *is* the subtractive form for fitting
  purposes; the clamp matters only for a production enum (Step 3).
- **Shifted cubic — stays parked** (dropped by user decision 2026-07-23;
  re-add condition was "genuine threshold structure in the trace tails" —
  not found, next section).

### Low-v trace-tail re-inspection (smoothed 9/18 Å references)

- **18 Å** (window 4.54–14.77 ps): monotonic deceleration to the window
  floor v = 2.542 Å/ps at window end; no coast (mean dv/dt −0.041 Å/ps²
  over the last 1 ps). The force-balance read γ_emp ≈ 3.2–3.4 amu/ps at
  v ≈ 2.6–2.75 sits far below the ρ̂ = 1 cubic value (~17) —
  the late window is surface/exit-well mixed (gate < 1 + barrier
  climb), **not a clean low-v drag probe**; the local-exponent read
  n_loc ≈ 1.0 over v < 3 is attribution-confounded for the same reason.
- **9 Å** (window 2.67–14.08 ps): reaches min v = 2.833 Å/ps at
  12.41 ps, then *slightly re-accelerates* (+0.033 Å/ps² over the last
  1 ps, v(end) = 2.872). A flattening exists — but its attribution is
  ambiguous three ways (escape/ρ̂ → 0; genuine threshold behavior; the
  known Tier-0 Finding-2 transverse contamination), and the local
  exponent is unusable there (n_loc 38–288 → the tail is not
  drag-dominated).
- **Conclusion:** the existing traces contain **no drag-dominated
  samples below v ≈ 2.5 Å/ps** and no resolvable coast/threshold
  structure. The sub-2.5 Å/ps drag shape — where the candidate v_f
  lives and the RQ11 deficit's lower half sits — is **unconstrained by
  the existing traces**. The collaborator ask (TDDFT tails / extended
  traces / production-kinematics traces) is the only instrument that
  reaches it; the shifted-cubic re-add condition is not met.

## D4 §6.6 — Quadratic counterfactual arbitration (2026-07-24, zero MD, twin)

**Question (pre-registered, plan §6.6):** could the lq form have landed
the Tier-2 observables if the downstream knobs had been re-arbitrated
around it? **Answer: yes — outcome (b), measured: the Tier-2
observables cannot discriminate the drag forms.** (The initial
coarse-grid read was (a); the user-caught fine (τ, E₀) rescan below
overturned it — recorded in full because the correction is itself an
atlas result about basin width.)

**Design as executed.** lq with its **own** Tier-0 artifacts (shared
a = 9.805·10⁻⁵ ≈ 0, c = 12.792 amu/Å, co-extracted E_bind = 0.0482 eV)
through the reconstructed S6 twin (leg-D draw discipline, m = 20000;
oracle reproduced the committed `h2b_s6_final` v7.25 row bit-for-bit
before any lq cell). Downstream freedom: v_c ∈ {∞, 10.5, 9.0, 8.9, 8.8,
8.65, 8.0, 7.25, 6.5, 6.0, 5.5} with the constant-force tail
(p_tail = −1 analog; 8.65 = the cell whose tail force equals the
production cubic's b·v_c³) × τ ∈ {1.6, 2.4, 3.2, 4.8, 6.4} ps × E₀ ∈
{0.17, 0.22, 0.27, 0.32, 0.37} eV — **275 cells**. Pre-registered
comparability vs the base cubic twin row: |ΔW₁| ≤ 0.08 ∧ midHot ∈
[0.84, 1.14] ∧ |Δsupp| ≤ 0.04 ∧ |Δn₁s| ≤ 0.03; deep bins
direction-only. Machinery: `quad_counterfactual_twin.py` +
chord workers + `quad_score_all.py`, scratchpad-tier; CSV
`quad_counterfactual_results.csv` (session scratchpad).

**Coarse-grid pass (275 cells): 0 comparable.** At the standing
(τ, E₀), the v_c continuation seemed to show the two shape observables
demanding disjoint v_c windows:

| lq v_c (τ 3.2, E₀ 0.27) | supp | n₁s | W₁ | midHot | d10 | d17 |
|---|---|---|---|---|---|---|
| 8.65 (tail force ≡ cubic) | 0.201 | 0.239 | **0.696 ✓** | 1.225 ✗ | 1.22 | 0.76 |
| 8.8 | 0.199 | 0.236 | 0.777 ✗ | 1.187 ✗ | 1.20 | 0.79 |
| 8.9 | 0.198 | 0.234 | 0.833 ✗ | 1.164 ✗ | 1.19 | 0.80 |
| 9.0 | 0.197 | 0.231 | 0.890 ✗ | **1.144 (✓)** | 1.18 | 0.80 |
| cubic base (7.25) | 0.202 | 0.243 | 0.678 ✓ | 1.067 ✓ | 0.99 | 0.71 |

That reading was a **grid-coarseness artifact** (user-caught, same
day): the coarse τ grid stepped 3.2 → 4.8, and at fixed chord the
τ-leg moves midHot ~5× faster than supp — leaving room for an
intermediate τ to pull midHot into band before supp leaves it. Since τ
and E₀ are post-chord knobs, the refinement was free:

**Fine (τ, E₀) rescan — 11 chords × τ 2.6–5.2 (step 0.1) × E₀
0.20–0.36 (step 0.01) = 4,959 cells: SIX four-way passes**, at
v_c ∈ {8.8, 8.9, 9.0} × τ ∈ {3.3, 3.4, 3.5} × E₀ = 0.27:

| cell (τ, E₀ 0.27) | supp | n₁s | W₁ | midHot | n₁KE | d10…d17 |
|---|---|---|---|---|---|---|
| vc9.0 τ3.5 (best) | 0.238 | 0.246 | **0.641** | **0.990** | 0.872 | 1.00…0.67 |
| vc9.0 τ3.4 | 0.224 | 0.243 | 0.708 | 1.038 | 0.898 | 1.06…0.70 |
| vc8.9 τ3.4 | 0.226 | 0.245 | 0.663 | 1.056 | 0.914 | 1.07…0.70 |
| vc8.9 τ3.3 | 0.212 | 0.240 | 0.735 | 1.107 | 0.941 | 1.12…0.75 |
| vc8.8 τ3.4 | 0.227 | 0.247 | 0.615 | 1.076 | 0.929 | 1.07…0.69 |
| vc8.8 τ3.3 | 0.213 | 0.242 | 0.690 | 1.130 | 0.956 | 1.14…0.72 |
| cubic base (7.25, 3.2) | 0.202 | 0.243 | 0.678 | 1.067 | 0.972 | 0.99…0.71 |

The best lq cell beats the cubic base on W₁ *and* midHot centering,
with deep bins comparable (not inverted — the standing-(τ, E₀) hot
deep bins were a slice artifact too). Residual visible difference:
n₁KE runs ~10 % colder (0.87 vs 0.97). Structural observation that
survives the correction: the lq landing is a **needle, not a basin** —
passes exist only at E₀ = 0.27 exactly (grid step 0.01) and a τ window
~0.3 ps wide, vs the cubic §4w basin (v_c 7.25–7.5 with whole-grid
tolerance); seed-robustness of the needle is untested (fixed twin
draws).

**Pre-registered interpretation applied: outcome (b).** The Tier-2
observable vector **cannot discriminate the drag forms** — an
lq-based system re-arbitrated to (v_c ≈ 8.8–9.0, τ ≈ 3.3–3.5,
E₀ 0.27) lands the histogram and KE observables comparably to the
production cubic. Per the pre-registration, the form choice therefore
rests **solely and explicitly on the Tier-0 trace instruments** —
where lq remains rejected (shared Δobj +0.0339; held-out 9 Å 0.699
FAIL; free exponent n̂ = 2.927). The §6.1 anti-circularity note is
**upheld and sharpened**: the Tier-2 landing must never be cited as
evidence for cubic — that is now a *measured* statement, not a
caution. **Consequences (pre-registered):** (1) a twin (b) needs an
**MD spot-check** before any strong claim (lq enum build behind
`[PROCEED TO IMPLEMENTATION]`, N = 500 at an lq-passing cell);
(2) the **conditional follow-up fires**: Method-B form discrimination
re-run under the Tier-1a anchored variable-mass m(t), turning the
~0.2-in-n exponent-bias estimate into a measurement — both are
user decisions, neither moves the standing point (atlas stance).
**Authority box:** frozen-chord twin (I69/I73, channel (d)
under-expressed); deep-bin columns direction-only; m = 20000 twin
ensemble with fixed draws, not MD; the needle-vs-basin contrast is a
twin-geometry statement pending MD/seed confirmation.

### §6.6 MD spot-check (2026-07-24, 3 × N = 500, triggered) — all three cells LAND; (b) is MD-confirmed and the twin needle was a frozen-chord artifact

**Build:** `capped_linear_quadratic` enum (`physics/drag.py` +
`config.py` guard, tail conventions mirrored from `capped_cubic`,
in-band byte-identical to `linear_quadratic`);
`scripts/gen_tier2atlas_spotcheck.py` (atlas namespace
`tier2atlas_conf270_qc{a,b,c}`, substring locks; cfg verified
field-by-field against the standing finc1v725 `cfg.json` — diff exactly
{drag_form, drag_coefficients, binding, τ}); cells built from the
`shared_lq` Tier-0 bundle so drag ↔ E_bind (0.0482 eV) is a consistent
§6.5.1 pair. Same seed as finc1v725 (20260721) — the same-N same-seed
cubic row is the natural comparator. Scoring: committed scorer, oracle
rows (finc launch record + pooled §4cc) reproduced before any new read.

| cell (all E₀ 0.27) | supp | n₁s | W₁ | midHot | χ²_med | vs twin prediction |
|---|---|---|---|---|---|---|
| **qca** lq v_c 9.0 / τ 3.5 | 0.202 | 0.2715 | **0.5231** | 0.919 | **64.7** | twin PASS → **MD lands** |
| **qcb** lq v_c 9.0 / τ 3.2 (control) | 0.158 | 0.2621 | **0.5597** | 1.069 | **61.3** | twin FAIL (W₁ 0.89) → **MD lands anyway** |
| **qcc** lq v_c 8.8 / τ 3.4 | 0.191 | 0.2725 | **0.4956** | 1.001 | **82.6** | twin PASS → **MD lands** |
| finc1v725 (cubic, same N + seed) | 0.161 | 0.2718 | 0.5238 | 0.990 | 125.7 | — |
| pooled battery (cubic, N = 5000) | 0.187 | 0.2433 | 0.5713 | 1.014 | 242.0 | — |

**Reads (pre-registered grid):**

1. **A lands + B lands → the twin needle was a frozen-chord artifact.**
   The off-needle control (τ 3.2, twin-predicted ~5-seed-SD W₁ failure)
   lands *better* than the pooled cubic reference. Mechanism feedback
   (generative pickup/evaporation/E_int dynamics) blurs the twin's
   knife-edge into a wide MD basin — the needle-vs-basin contrast of the
   twin table is **not** a physical discriminator.
2. **Outcome (b) is MD-confirmed:** the Tier-2 observable surface cannot
   discriminate the drag forms *in full MD* — an lq-based system lands
   the histogram (W₁ 0.50–0.56 vs cubic 0.52/0.57) and the KE bands at
   multiple (v_c, τ) points. The form choice rests solely on the Tier-0
   trace instruments (lq held-out FAIL stands); the experimental landing
   supports **no** form claim, cubic included — now an MD-measured
   statement, closing the §6.1 anti-circularity question at every level.
3. **The lq cells halve the median-anchored KE χ²** (61–83 vs the
   same-N-same-seed cubic 125.7) with midHot staying in-band — under the
   atlas stance this is *reported, not adjudicated*: no adoption path
   exists here (plan §8), and Tier-0 trace authority still rejects lq.
   What it demonstrates is methodological: downstream-mechanism
   flexibility, not drag-form correctness, owns much of the landing
   quality.
4. **RQ11 signal (direction-only, thin bins):** deep-bin ratios per cell
   — qca 0.85…0.55, qcb 1.02…0.45, qcc 0.83…—, vs pooled cubic
   0.81…0.38. The deficit's *slope* persists under lq (form-robust),
   but qcb (standing τ) runs ~+0.1–0.2 warmer across n = 10–15 —
   consistent with lq's *lower* γ in the 5–9 Å/ps mid-band (its cap sits
   at 9.0 and its γ crosses below capped-cubic above ~5), pointing the
   deep-KE lever at the **upper-mid band (v ≈ 5–9)** rather than the
   sub-2.5 region. Recorded to `RESEARCH_QUESTIONS.md` as NB-RQ11-12.
   Counts n = 10–17 are ~2–37 per cell — direction-only; an N = 1000
   confirmation is the escalation if this read is to become
   headline-bearing.

**Authority box:** single fixed seed per cell at N = 500 (seed-SD ×√2
vs battery members); χ²_med compared same-N only (it scales with scored
count); deep bins thin; suppressed-class spread (0.158–0.202) within
the known cell-to-cell scatter. Nothing here moves the standing point;
`finc1v725` remains production (atlas stance).

### Status of the D4 form axis after Step 1

Registered priors unchanged: lq stays rejected (shared +0.0339, held-out
FAIL), power law stays cubic-equivalent, linear+cubic stays a → 0. New
knowledge from the zero-cost gate: Padé excluded as a cap-replacement
*by arithmetic* (no fit needed); subtractive v_f is trace-constrainable
and awaits a promotion decision; sub-2.5 Å/ps stays TDDFT-blind. Step 2
(twin sweeps) can proceed on the standing form list; Step 3 (MD
spot-checks / enum build) has no candidate yet that survives to it.

## §6.7 item 1 — lq sanity battery (5 × N = 1000, paired-by-seed) — PRE-REGISTRATION (frozen 2026-07-24, before launch)

**Frozen before the first run** per plan §6.7 (BN-style pre-registration
discipline). This block is the pre-commitment; the results subsection is
appended below it *after* the runs land. Nothing here moves `finc1v725`
(atlas stance).

**Purpose.** Escalate the §6.6 MD spot-check (3 × N = 500, single seed) to
the N = 1000 battery scale with a **paired-by-seed** design, so the
form-blindness statement (outcome (b), MD-confirmed) carries a paired-SD
error bar rather than a √2-inflated single-seed read; and escalate the
NB-RQ11-12 deep-bin mid-band read (~300 fragments at n = 10 pooled) from
direction-only to headline-bearing.

**Design.** Base cell **qcc** (chosen on MD merit at N = 500:
W₁ 0.496, midHot 1.001) — `capped_linear_quadratic`, `shared_lq` Tier-0
bundle (a = 9.805e-05, c = 12.792, jointly-extracted E_bind = 0.0482 eV;
a *consistent* §6.5.1 pair, **no** binding escape hatch), v_c 8.8 / p_tail
−1, τ 3.4, E₀ 0.27, standing finc1v725 pins otherwise (leg-D uniform-volume
birth, Landau-gated E2 arm). **5 members, N = 1000, seeds 20260722–20260726
= the cubic battery's exact seeds** → each member pairs by seed with
`bigc1v725s{1..5}` (same neutral draws), so every observable becomes a
per-seed paired Δ and seed noise cancels. Namespace
`tier2atlas_conf270_qccbigs{1..5}`.

**cfg discipline.** Each member is verified field-by-field against **its
paired cubic member's** `cfg.json` (matched N + seed) — the diff must be
exactly {`drag_form`, `drag_coefficients`, `binding_energy_I_ion_eV`,
`internal_energy_cooling_tau_ps`} (the form swap + the lq pair's E_bind +
τ 3.2 → 3.4). Any other diff aborts the build (battery precedent: the
run-dir cfg IS the spec).

**FROZEN predicted pass bands** (met → landing reproduces at battery scale
and form-blindness holds paired; deviation → a real finding, reported not
adjudicated):

| # | observable | frozen prediction |
|---|---|---|
| 1 | pooled W₁_solv | ∈ ~[0.46, 0.57] |
| 2 | midHot (pooled) | in the seed-robust band [0.84, 1.14] |
| 3 | deep-bin KE slope | deficit persists (form-robust); the NB-RQ11-12 mid-band warm read escalates, not reverses |
| 4 | per-member χ²_med | below its paired cubic member (same-N paired) |
| 5 | suppressed weight | within cell scatter of the standing 0.187 |

**Confound caveat (carried in every results table).** qcc's E_bind 0.0482
is ~0.07 eV shallower than the cubic pair's 0.1168 — large on the deep-bin
KE scale. This battery tests **reproducibility of the landing and the
deep-KE slope at scale**; it does **not** separate the drag-form effect
from the E_bind effect. That separation is the deferred §6.7 item 2
E_bind pair-separation scan; until it runs, any "lq lands better" / KE-χ²
advantage stays attribution-provisional.

**Authority box (pre-committed).** Paired-SD error bar from 5 seeds;
χ²_med compared same-N paired only (scales with scored count); deep bins
still thin per member but pooled ~300 at n = 10. Twin authority does not
enter (this is MD). Adoption remains outside this program (plan §8);
Tier-0 trace authority still rejects lq.

**Status: FROZEN 2026-07-24, pre-launch. MD pending.**

### §6.7 item 1 — lq sanity battery RESULTS (2026-07-24, 5 × N = 1000 executed) — histogram form-blindness CONFIRMED; the §6.6 KE advantage does NOT replicate (paired); lq over-suppresses

**Runs.** `9A_drag_shared_lq_N1000_tier2atlas_conf270_qccbigs{1..5}` at seeds
20260722–26, paired 1:1 with `bigc1v725s{1..5}` (same neutral draws).
Scored with the §6.6 committed instrument reproduced verbatim
(`atlas_spotcheck_score.py` → the battery driver `lqbattery_score.py`;
midHot = mean of sim/ref over n = 2–8, mean-anchored; χ²_med = median-anchored
KE chi-square). **Oracles PASS bit-exact before any new read:** finc1v725
(943 / supp 0.161 / W₁ 0.5238 / midHot 0.9905 / χ²_med 125.7) and the pooled
§4cc cubic row (9330 / W₁ 0.5713 / midHot 1.0139 / χ²_med 242.0).

**Per-member, paired by seed** (cubic / lq; W₁ · midHot · supp · χ²_med):

| seed | cubic W₁·midHot·supp·χ² | lq W₁·midHot·supp·χ² | χ²_med |
|---|---|---|---|
| s1 20260722 | 0.728 · 1.036 · 0.192 · 155 | 0.692 · 1.042 · 0.217 · 176 | lq **>** cub |
| s2 20260723 | 0.606 · 1.023 · 0.185 · 114 | 0.615 · 1.029 · 0.213 · 120 | lq **>** cub |
| s3 20260724 | 0.519 · 1.021 · 0.188 · 129 | 0.504 · 1.026 · 0.206 · 143 | lq **>** cub |
| s4 20260725 | 0.560 · 1.011 · 0.192 · 136 | 0.481 · 1.008 · 0.219 · 161 | lq **>** cub |
| s5 20260726 | 0.482 · 0.981 · 0.180 · 123 | 0.487 · 1.007 · 0.209 · 124 | lq **>** cub |

**Paired Δ (lq − cubic), mean ± sample SD (n = 5):**

| observable | Δ (lq − cubic) | read |
|---|---|---|
| supp | **+0.0255 ± 0.0041** | ~6σ paired — lq over-suppresses, consistent |
| trap | **−0.0305 ± 0.0019** | ~16σ paired — lq under-traps (weight → suppressed/detected) |
| χ²_med | **+13.3 ± 10.0** | lq KE chi² *worse* on all 5 seeds (pooled 277 vs 242) |
| n̄_det | −0.059 ± 0.028 | marginally colder size |
| W₁ | −0.023 ± 0.036 | comparable (lq mildly better, not significant) |
| midHot | +0.008 ± 0.011 | comparable |
| n₁_solv | +0.001 ± 0.006 | comparable |

**Pooled (concat 5 × N = 1000):**

| pool | scored | trap | supp | n̄ | n₁s | W₁ | midHot | χ²_med |
|---|---|---|---|---|---|---|---|---|
| lq | 9635 | 0.036 | 0.213 | 4.008 | 0.2447 | **0.5484** | **1.0232** | 277.3 |
| cubic (=§4cc) | 9330 | 0.067 | 0.187 | 4.068 | 0.2433 | 0.5713 | 1.0139 | 242.0 |

deep strip (n = 10–17, sim/ref, mean-anchored, counts): lq
`0.84(312) 0.78 0.76 0.64 0.53 0.42 0.44 0.40(41)`; cubic
`0.81(301) 0.75 0.74 0.67 0.59 0.58 0.55 0.38(17)`.

**Frozen-band scorecard: 3 MET / 2 MISSED.**

- **1 W₁ MET** (0.548 ∈ [0.46,0.57]); **2 midHot MET** (1.023 ∈ [0.84,1.14]);
  **3 deep slope MET** (both forms decline n10→n17, deficit persists,
  form-robust). → **Histogram-level form-blindness is confirmed at battery
  scale** (outcome (b) holds for W₁/midHot/n̄/n₁_solv; lq lands the histogram
  as well as — mildly better than — cubic on W₁).
- **4 per-member χ²_med < paired cubic — MISSED, REVERSED.** Predicted lq
  lower; measured lq **higher on all 5 seeds** (Δ+13.3±10.0; pooled 277 vs
  242). **The §6.6 "lq halves χ²_med / lands better on KE" does NOT
  replicate.** That read was a single-seed artifact: §6.6 compared qcc to
  finc1v725 at seed 20260721/N = 500 (χ²_med 82.6 vs 125.7 — lq favoured);
  across the 5 *different* battery seeds the sign flips. The form-χ²
  difference is within seed-noise range and does **not** robustly favour lq.
  (lq scores ~3 % more ions, which inflates χ² modestly — but the ~10–15 %
  gap and 5/5 consistency exceed that; the direction is robust.)
- **5 supp within cell scatter of 0.187 — MISSED.** lq pooled supp 0.213 vs
  0.187; paired Δ+0.0255±0.0041 (member SD 0.0053) — a real, consistent
  drag-form effect on the fate split: **lq shifts weight trapped → suppressed**
  (Δtrap −0.031). This is a new resolved-influence entry, not noise.

**Synthesis.** The paired battery cleanly separates two claims the §6.6 N = 500
read had fused: (i) lq **lands the histogram** — confirmed and robust (W₁,
midHot, n̄, n₁_solv all comparable to cubic; form-blindness (b) upheld at
N = 1000); (ii) lq **lands the KE distribution better** — **refuted**: paired,
lq's median-anchored KE chi² is worse on every seed, and lq over-suppresses.
The §6.6 "lands better" was a favourable-seed artifact, exposed by the
paired-by-seed design *before* the E_bind confound is even addressed.

**Deep mid-band (NB-RQ11-12 escalation, ~300 counts/bin).** lq runs warmer
than cubic at n = 10–12 (0.84/0.78/0.76 vs 0.81/0.75/0.74) and comparable/
colder at n = 14–17 — the mid-band (v ≈ 5–9 Å/ps) warming signature, now at
headline count. But the overall KE chi² is worse, so the mid-band warming
does *not* translate into a better KE landing.

**Confound status.** Both forms here carry lq's E_bind 0.0482; this battery
does **not** separate form from E_bind. It does show the KE advantage was not
even seed-robust, independent of E_bind — which the §6.7 item-2 E_bind
pair-separation scan (deferred) will still probe for the supp/deep effects.

**Atlas stance.** Nothing here moves `finc1v725`; Tier-0 trace authority still
rejects lq (held-out FAIL). Reported, not adjudicated. **Status: EXECUTED
2026-07-24; pre-registration outcome 3 MET / 2 MISSED (bands 4, 5).**

## §6.7 item 2 — E_bind pair-separation scan RESULTS (2026-07-24) — the over-suppression is a FORM effect; item-1's histogram match is a (form, well) CO-COMPENSATION

**Question.** Item 1's lq/cubic differences (over-suppression; colder KE) are
entangled with lq's co-extracted **shallower exit well** (E_bind 0.0482 vs the
cubic pair's 0.1168). This scan holds the lq drag law fixed and sweeps **only
the well** so the form-vs-well attribution becomes a measurement. The drag
bundle's stamped `effective_binding_energy_I_ion_eV` stays 0.0482 (provenance);
only the run's climbed well is overridden, honestly under
`allow_unvalidated_binding_pairing` (§6.5.1). Machinery:
`gen_tier2atlas_ebindscan.py` (initial, N = 500) + `gen_tier2atlas_ebindseeds.py`
(firm-up, N = 1000); scored with the committed row (`ebindscan_score.py` /
`ebindseeds_score.py`), oracle A (finc) bit-exact first.

### Initial single-seed scan (N = 500, seed 20260721) — directional

lq E_bind OAT (well →): supp 0.191 / 0.194 / 0.196 (flat); trap 0.038 / 0.075 /
0.100 (rises); W₁ 0.496 / 0.638 / 0.793; midHot 1.001 / 0.932 / 0.887 at
E_bind 0.0482 / 0.1168 / 0.154. Read as directional only — the W₁ numbers proved
seed-driven (see firm-up).

### Firm-up (N = 1000, seeds 20260722/23/24, paired to the item-1 battery)

Each cell pairs 1:1 with `bigc1v725s{k}` (cubic @ 0.1168) and `qccbigs{k}`
(lq @ 0.0482) at matched seed + N, so **form** and **well** are each isolated
across 3 paired seeds. (`eb1168` all landed; `eb154` at 0.154 tripped the
detection handover guard on 2/3 seeds — see deep-well note.)

**FORM at matched well 0.1168 — Δ(lq@0.1168 − cubic@0.1168), mean ± SD (n = 3):**

| observable | Δ | read |
|---|---|---|
| midHot | **−0.070 ± 0.002** | ~35σ — lq fragments much colder |
| n̄_det | **−0.568 ± 0.029** | ~20σ — lq clusters much smaller |
| supp | **+0.034 ± 0.004** | ~8σ — **lq over-suppresses (FORM, at matched well)** |
| trap | **+0.027 ± 0.003** | ~9σ — lq traps more |
| n₁_solv | +0.020 ± 0.003 | more n = 1 |
| W₁ | +0.009 ± 0.015 | **not significant — comparable** |
| χ²_med | +3.9 ± 38.6 | **inconclusive** (seed-noise dominates even at N = 1000) |

**WELL on lq — Δ(lq@0.1168 − lq@0.0482), mean ± SD (n = 3):** trap
**+0.058 ± 0.004** (the clean well lever — deeper well retains more); n̄
−0.503 ± 0.052; midHot −0.076 ± 0.003; supp +0.010 ± 0.002 (small); W₁
+0.023 ± 0.012 (marginal); χ²_med −9.7 ± 34.7 (ns).

### Resolved conclusions

1. **The over-suppression is a FORM effect, not the well** — at the *matched*
   well 0.1168 lq still over-suppresses vs cubic by +0.034 (~8σ); the well
   adds only +0.010 more (item-1's Δ+0.026 was ~⅓ well, ⅔ form). The item-1
   confound caveat is **discharged: it is the quadratic form.**
2. **Item-1's histogram form-blindness is a (form, well) CO-COMPENSATION.** At
   matched well, lq produces **distinctly colder (midHot −0.070) and smaller
   (n̄ −0.57) fragments** than cubic. The reason item-1 saw matching midHot
   (~1.02) / n̄ / W₁ is that lq's intrinsic coldness was cancelled by its
   shallower 0.0482 well. Match the wells → the form's coldness is exposed.
   The forms are physically **distinguishable**, not interchangeable.
3. **W₁ specifically is form-blind even at matched well** (Δ+0.009, ns): the
   solvated-histogram Wasserstein does not see the form, though midHot / n̄ /
   supp / trap do (8–35σ). So W₁ is a weak form discriminator; the KE / size /
   fate observables are strong ones.
4. **The single-seed W₁ reads were seed-noise.** The N = 500 seed-20260721
   scan suggested "cubic lands W₁ much better at matched well" (Δ+0.114) and
   "lq's W₁ degrades strongly with the well" — both **refuted** by the firm-up
   (matched-well ΔW₁ +0.009; well effect +0.023). Recorded as a caution: at
   N = 500 single seed W₁ swings ~0.1 on seed alone.
5. **χ²_med stays inconclusive** for form-vs-well even at N = 1000 × 3 seeds
   (SD ±38) — not chased further.
6. **Deep-well over-retention (physical, not a bug).** At E_bind 0.154, 2/3
   N = 1000 cells hit the P1–P3 detection handover guard: ~2/2000 ions never
   decoupled from He in the 8000 ps window (ρ̂ ≈ 0.93, live drag exposure, not
   energetically-bound so `retained_policy=exclude` does not catch them). The
   one seed that landed (s3) shows trap 0.125 / n̄ 3.28 — the deep well
   over-retains, to the point of numerical non-decoupling. eb154 is therefore
   a 1-seed point; the deep arm is not pursued further.

**Atlas stance.** Nothing here moves `finc1v725`; adoption stays outside the
program; Tier-0 still rejects lq. The scan's role was diagnostic — it
discharges the item-1 E_bind confound (→ form) and reframes the histogram
form-blindness as co-compensation. **Status: EXECUTED 2026-07-24 (initial
N = 500 single-seed + firm-up N = 1000 × 3 seeds).**

## Axis A pre-read — the pooled N = 5000 battery binned by its own sampled geometry (2026-07-26, zero MD)

Plan §3.4's "optional follow-up", executed **before** the Axis A grid because
it is free: the pooled battery already samples R ∈ [14, 54] Å and birth depth
∈ [3, 38] Å, so binning it gives the geometry response of the standing system
at zero cost. Scratchpad read; committed §4cc scorer reused verbatim
(`postprocess/tier2_confirmation`), no repo code added.

**Oracle (pooled row vs the recorded §4cc values):** n₁_solv **0.2433**
(recorded 0.243), W₁ **0.5713** (0.571), supp **0.1872** (0.187), n̄ 4.068,
trapped 0.067, deep-KE 0.631. midHot **1.0110** vs recorded 1.014 — a 0.3 %
definitional residue, not a convention choice (min_count 1 and 2 give the
identical 1.0110); flagged, everything else reproduces at recorded precision.

**Entanglement caveat, stated once and applying to every number below.** These
bins are *sampled*, not controlled: R and birth depth co-vary, and the T5
`density_tied` dressing ties the birth shell directly to depth
(n₀ = round(n*·ρ̂(depth)) — measured **n₀ = 14.0 / 16.4 / 19.3** across the
depth terciles). This read orders and bounds; it does not separate exposure
from dressing. That is the controlled grid's job.

### A. Binned by droplet radius (terciles)

| | R̄ 21.0 | R̄ 26.2 | R̄ 32.9 |
|---|---|---|---|
| trapped | 0.0003 | 0.042 | **0.158** |
| supp | 0.232 | 0.182 | 0.141 |
| n̄ | 2.90 | 4.23 | 5.26 |
| n₁_solv | 0.302 | 0.229 | 0.196 |
| W₁ | 1.127 | 0.750 | 1.360 |
| midHot | 0.941 | 0.990 | 1.106 |
| deep-KE | 0.526 | 0.529 | 0.757 |

### B. Binned by birth depth (terciles) — the stronger ordering

| | depth 4.4 Å | depth 7.7 Å | depth 14.6 Å |
|---|---|---|---|
| supp | **0.381** | 0.182 | **0.003** |
| n̄ | 2.36 | 3.52 | 6.29 |
| n₁_solv | 0.358 | 0.322 | 0.109 |
| W₁ | 1.083 | 0.588 | 1.593 |
| midHot | 0.703 | 0.892 | 1.209 |
| deep-KE | **0.316** | 0.417 | **0.792** |
| ⟨n₀⟩ (dressing) | 14.0 | 16.4 | 19.3 |

### C. 2-D (R × depth) — **the registered "R dominates" expectation is REFUTED**

At **fixed birth depth, droplet radius is nearly inert** on the scored
ensemble (supp 0.392/0.377/0.366 across R at depth 1; 0.201/0.175/0.161 at
depth 2; 0.005/0.001/0.004 at depth 3 — deep-KE likewise flat at
0.34/0.33/0.31, 0.51/0.40/0.41, 0.84/0.69/0.87). Only n̄ retains an R slope at
the deepest births (4.11 → 6.07 → 7.76).

**Mechanism.** The escaping fragment's He path is its *birth depth*, not the
isotropic chord: at fixed depth the radial path to the surface is depth,
independent of R. R acts on the *inward*-going Coulomb partner, whose path is
≈ 2R − depth — and those fragments leave the scored ensemble instead of
changing it (trapped 0.0003 → 0.158 across the R terciles). So **R is a
selection knob, depth is the physics knob.** The §3.1 pre-registration
("R moves the chord ×2.3, the law only ×1.24, so expect an R-dominated
response") used the isotropic chord as the exposure coordinate and is
withdrawn; birth depth is the right coordinate.

### D. The histogram landing is a MIXTURE property

Pooled W₁ = 0.571 is **better than every single geometry bin** (best tercile
0.750; best quintile 0.639; edges 1.13–1.64). No single (R, depth) lands the
solvated histogram — the ensemble average does. **Consequence for Axis A
scoring:** fixed-R/fixed-law cells must be compared to each other and to a
re-weighted mixture reconstruction, *never* judged directly against the
committed acceptance. Scoring nine fixed cells against the acceptance would
report nine spurious failures.

### E. RQ11 signal — the deep-bin deficit largely closes with birth depth

deep-KE ratio 0.316 → 0.792 across the depth terciles (0.63 pooled). Geometry
— or the dressing that rides on it — is the first candidate this program has
found that moves the RQ11 cold tail by a factor of 2.5. Confounded per the
caveat above; disentangling requires the `initial_shell_model="full"` control
(see the plan §3.1 amendment proposal).

### F. Pre-registered prediction for the anchored cells

At the parent birth law (mean depth ≈ 25 Å at R1 = 26.6 Å) the dressing
saturates: ρ̂ = 0.994 → **n₀ = 20.9 of 21**. Extrapolating trend B: suppression
→ ≈ 0, n̄ well above 6.3, n₁_solv ≈ 0.1 or below, W₁ ≳ 1.6. **The prediction is
that the landing breaks at the anchored geometry, and that it breaks through
the dressing/suppression channel rather than through drag exposure.** Frozen
here before the grid runs.

### G. Initial shell *energy* — the dressing is self-similar, and it saturates (2026-07-26, zero MD)

Follow-on analysis of §B/§C, prompted by the question of what the birth
dressing does to the *energy* budget rather than the shell count. Numbers at
the standing ladder (`rq4graded`, Σ(21) = 0.20629 eV, E₀ = 0.27 eV):

| birth depth | ρ̂ | n₀ | Σ(n₀) [eV] | Σ(n₀)/Σ(21) | E_int(0) [eV] | G = E_int − Σ |
|---|---|---|---|---|---|---|
| 4.4 | 0.669 | 14 | 0.1476 | 0.715 | 0.1932 | +0.0456 |
| 7.7 | 0.778 | 16 | 0.1660 | 0.805 | 0.2173 | +0.0513 |
| **8.9** (standing) | 0.812 | **17** | 0.1751 | 0.849 | **0.2292** | +0.0541 |
| 14.6 | 0.927 | 19 | 0.1927 | 0.934 | 0.2522 | +0.0595 |
| **25.0** (parent law) | 0.994 | **21** | 0.2063 | **1.000** | **0.2700** | +0.0637 |

1. **The as-built dressing is energetically self-similar.** Under T6
   `sigma_proportional` (p = 1) the onset is tied to the same n₀ that T5 sets,
   so **E_int(0)/Σ(n₀) = E₀/Σ(n*) = 1.309 at every depth**. The self-unbound
   margin G = 0.3088·Σ(n₀) is positive everywhere and **never changes sign
   with dressing** — only its absolute scale moves (0.046 → 0.064 eV). The
   depth ordering measured in §B is therefore *not* static energetics: it is
   the dynamic race (RRK reads absolute E_int; and Newton cooling is itself
   ρ̂-gated under `cooling_spatial_gate="density_scaled"`).
2. **Birth depth is a 4-leg bundle**, three of them keyed to the same
   ρ̂(depth): shell count (T5), onset energy (T6, via n₀), cooling rate (the
   cooling gate), plus the geometric drag exposure. No single-cell contrast
   can attribute the §B ordering; hence the C1/C2 controls (plan §3.1c).
3. **The dressing saturates by depth ≈ 25 Å** (ρ̂ 0.994 → n₀ = 21, ratio 1).
   At the grid's L1 (center-pin, depth 26.6 Å) and L2 (parent law, 24.8 Å)
   the T5 and T6 arms are **structurally inert**. Controls placed there would
   have measured nothing — they belong on L3, the shallow standing law, the
   only row where the dressing is live. (The original control proposal was
   corrected on this finding before launch.)
4. **Corollary — T5 and T6 are artifacts of the shallow birth law.** The birth
   heterogeneity those two delivered arms encode exists only because
   production births molecules ~9 Å below the surface; at the parent model's
   ~25 Å they do nothing. Recorded as a structural observation, not a defect
   claim: the arms were built for twin parity, and this is what they turn out
   to describe.

## Axis A / G1 — the controlled geometry grid (2026-07-26/27, MD: 11 cells launched, 6 scored)

Stage 2b executed under its trigger. `gen_tier2atlas_geometry.py`, N = 500 per
cell, **one shared seed 20260727** (common random numbers, so cross-cell
differences are partly paired), fixed droplet size per cell (size *sampling*
off — the radius is the controlled variable), every other knob at the standing
production point `finc1v725`, E_bind held at 0.1168 eV (G0-3). Scored with
`tier2atlas_geometry_table.py`; the pooled N = 5000 battery oracle reproduced
all seven recorded columns within 0.002 before any cell was read.

**Outcome: the R1/R2 rows (6 cells) are measured; the five anchored-radius
cells (R 49.4 / 68.3 Å) are blocked by a staging limit — see §G1.3.** The
influence tables and the width decomposition live in
`TIER2_PARAMETER_INFLUENCE.md` §14.1 (D0 is the influence reference; not
duplicated here). This section records what the axis *decided*.

### G1.1 The pre-registered prediction was confirmed on every clause

Frozen before launch (plan §3.1, from the pre-read §F): at the anchored birth
law the dressing saturates ⇒ supp ≈ 0, n̄ > 6.3, n₁_solv ≲ 0.1, W₁ ≳ 1.6, and
the landing breaks through the **dressing/suppression** channel rather than
through drag exposure. Measured at R1: supp **exactly 0**, n̄ **8.39**,
n₁_solv **exactly 0**, W₁ **4.98**. The anchored geometry does not shift the
histogram — it evacuates its low-n half. A confirmation, not a surprise; and
per the pre-registered G2 decision a broken landing at the anchored geometry
is not a reason to keep the wrong geometry.

### G1.2 Three results the grid establishes that the pre-read could not

1. **Size-distribution width is ~95 % geometry-inherited.** SD(n_det) tracks
   SD(birth depth) at ≈ 1 He per Å, on top of a mechanism-only floor of
   0.6–1.0 He exposed by the zero-spread center-pin cells. This is §3.3 Q5
   answered with controlled cells, and it is what the L1 column was for: its
   near-degeneracy with L2 in the *means* is exactly what makes it a clean
   zero-spread reference. (The design note that L1 and L2 "differ essentially
   only in birth-depth spread" is now a measurement, not a prediction.)
2. **The landing needs the size *distribution*, not the right mean size.**
   `r1l3` removes only the droplet-size spread from the standing configuration
   (R pinned at its own realized mean 26.6 Å) and W₁ degrades 0.571 → 0.813
   while supp / n₁_solv / trap / n̄ stay close — a controlled sharpening of the
   pre-read's mixture finding.
3. **RQ11 is on a geometry-traversable axis, and geometry overshoots it.**
   deepKE 0.51 → 0.75 → 1.38/1.49 across depth 9 → 11 → 29–34 Å, the
   deep-birth values resting on 6–8 occupied bins (the R1 row's 1–2 bins were
   not readable). The deficit does not merely close with depth: it crosses 1.0
   and overshoots to ≈ +40 % too hot, so some intermediate birth depth
   reproduces the deep-bin KE exactly. Reported under atlas stance — these
   cells' histograms are destroyed (W₁ 7.5–8.0), so this is a statement about
   *which axis the deficit lives on*, not a candidate point.

### G1.3 The anchored radii are blocked by a staging limit, not a bug

All five R ≥ 49.4 Å cells failed the detection-stage P1–P3 handover guard.
Diagnosed from the stored relaxation trajectories (read-only):

- **⅓–½ of all ions never leave the droplet** (324–470 of 1000 reclassified
  `droplet_retained` at R = 49.4 Å). Plan §3.6's "does the trapped channel
  explode at large R?" is answered **yes, emphatically**.
- The blocking ions sit **20–30 Å inside** the surface with |v_rad| ≈
  0.01–0.04 Å/ps. Cubic drag over a 25–50 Å path takes them below the Landau
  threshold (0.58 Å/ps), where `landau_gated_drag` switches dissipation off;
  they then oscillate conservatively in the droplet well. Median net radial
  progress over the window's last 4 ns is **−0.5 Å** (r3l1) / **−0.9 Å**
  (r4l2); ~35 % still move outward, ~18 % have flatly asymptoted.
- The escaping subset needs **~80–690 ns** more at its observed late drift
  rate: **10–90× the 8 ns E2 window**, and it needs that time *with drag and
  pickup live*, which E2 switches off below v_L. A modestly longer conservative
  window is therefore the wrong instrument.
- **The staging's timescale separation does not exist at the corrected
  geometry.** 30 ps MD → 8 ns E2 → 8.53 µs free flight is calibrated to
  R ≈ 27 Å, where ejected ions leave at several Å/ps and cover ~48 000 Å
  (≈ 5 µm) inside the E2 window. At R ≥ 49 Å escape time and flight time
  become comparable. Recorded as D0 §17's first ledger row that the geometry
  correction **creates** rather than retires.
- **The E_bind rider was NOT fired.** Its trigger clause included "a
  detection-handover-guard trip", written as a proxy for the trapped channel
  pressing against a boundary (the eb154 precedent, where a *deeper well*
  over-retained). Here the trip comes from undecided-fate ions in an
  undersized drag-active window, so an E_bind bracket cell would not address
  it. Treated as not met **in substance**, with the reason recorded, rather
  than executed mechanically.

**Consequence for the G-plan.** G1 is complete for R ≤ 34 Å and *blocked*
above it. G2 (the adoption decision) cannot be taken on the R3/R4 evidence it
was designed to read, because the model cannot currently produce that
evidence. The open question this hands forward is whether the ~µs in-droplet
residence at large R is a real prediction of the drag model or a sign that the
law over-dissipates at these path lengths — the Tier-0 calibration never saw a
25–50 Å path. Nothing adopted; `finc1v725` stands; F5 undischarged.

> **SUPERSEDED for the blocking claim only (2026-07-27, §G1.4):** the block is
> lifted by the retained-class arm and all 11 cells are scored. The physics
> diagnosis above stands unchanged; "the model cannot produce that evidence"
> is now false.

### G1.4 The grid completed — retained-class arm + the marginal bracket (2026-07-27, zero MD)

`detection_droplet_retained_policy="exclude_all_coupled"` (plan §3.5b) built
and run. The five blocked cells were completed **from their stored
`relaxation.npz` by re-running only the detection stage — zero new MD**, since
they had failed *after* the MD. Influence tables live in
`TIER2_PARAMETER_INFLUENCE.md` §14.1–§14.2 (D0 is the influence reference, not
duplicated here); this section records what the step *decided*.

**Arm oracle passed (§3.5b item 10).** Re-running detection on the six
already-scored R ≤ 34 Å cells under the new policy reproduced every field of
every `detection.npz` **bit-for-bit** — as it must, since every violator there
was provably bound, so the marginal class is empty at those radii. The
generator refuses to overwrite a differing result unless explicitly told to, so
this is enforced rather than asserted.

**The decomposition was load-bearing, not bookkeeping.** The coupled class's
*composition* inverts across the grid: at R = 49.4 Å it is 98–99 % **bound**
(physics — provably cannot escape), while at R = 68.3 Å on the parent law it is
75 % **marginal** (600 of 800 — a modelling convention). A single `trap` column
would have presented 0.476 and 0.800 as the same kind of quantity. Both
fractions are now printed separately, read through one shared constant.

**The pre-registered bracket came out NOT TIGHT — with the direction predicted
correctly.** Arm A (marginals excluded) vs Arm B (marginals injected at
handover `n` and their exact conservative asymptotic KE): tight at r3l2/r3l3
(marginal fractions 0.006 / 0.002), **wide at r3l1, r4l2, r4l3** — deepKE moves
−0.41 to −1.61 (10–38 seed-SD), n̄ up to +3.38, W₁ up to +3.35, while midHot is
flat to 0.000 and supp is unchanged, exactly as frozen. So the R ≥ 49 Å deepKE
and n̄ values are **conditional on the exclusion**; the R ≤ 34 Å rows — and
therefore all of §G1.2 — are not.

**Adjudication of the interpretation the arm was built under.** The recorded
read was that the retained class is slow, high-n, and does not move the result.
Measured: the first two clauses are **confirmed** (median handover n = 16–19,
median asymptotic KE clipping to **0.000 eV**, full-credit ceiling
0.031–0.094 eV against an experimental mean of 0.066 eV at the same n); the
third is **refuted**, because slow-and-high-n is not a null direction — it is
the RQ11 direction, so it lands squarely on the deep-bin observables and
nowhere else. Both halves come from the same table.

**A convention-level result worth its own line.** `KE_asym` uses the physical
half-credit pair-Coulomb split, whereas `_conservatively_bound` credits the
**full** pair energy to both fragments — a deliberate, documented over-estimate
that keeps a `bound` verdict certain. The median marginal ion clears its
barrier *only* under that over-crediting. So "marginal" means **"not provably
bound"**, and the marginal fraction is an **upper bound** on genuine slow
escapers rather than an estimate of them. This is why excluding the class stays
defensible even with a wide bracket: the population is physically dubious, and
the honest statement is that the anchored-radius KE observables carry a stated
range instead of a value.

**Still open, deliberately.** Arm B is the conservative corner — it holds `n`
fixed and lets the ion coast. Over the ~80–690 ns these ions need, live pickup
would raise `n` and mass-load them, plausibly into the bound class, which would
mean the conservative split *understates* trapping. That is the kinetic forward
model (conservative orbit + Poisson pickup + RRK), unfunded. The arm also does
not touch whether cubic drag over-dissipates on 25–50 Å paths (the collaborator
ask; the cheap distinguishing test remains that a *geometric* retained fraction
is insensitive to (v_c, b) while an over-dissipation artifact is not).

**Atlas stance intact:** nothing adopted, `finc1v725` stands, F5 undischarged,
G2 not taken.

## D2b §4.3 — grid re-weighting: oracle verdicts + the corrected-ensemble forecast (2026-07-27, zero MD)

`scripts/post_processing/tier2atlas_geometry_reweight.py` (tests:
`tests/test_tier2atlas_reweight.py`). The G1 grid used as the transfer
function; candidate size densities re-weighted onto it. Scorer-drift oracle
reproduced before anything was read.

### Method conventions (fixed by this session; recorded per plan §4.3)

- **Column-matched 1-D re-weighting.** Both candidate position laws coincide
  with grid columns exactly (production `uniform_volume` m3 ≡ L3; corrected
  Boltzmann 313.2 K ≡ L2), so conditional depth|R is the cell's own by
  construction and only the R marginal is interpolated. 2-D (R, depth)
  interpolation across columns was deliberately rejected: the L1/L2 columns
  are near-degenerate in mean depth but differ in spread, which makes a
  depth coordinate ill-posed across laws.
- **Interpolation convention: `nearest` (midpoint binning in R), adopted by
  the oracle** over `linear` (hat weights) — pre-registered tie-break: more
  gate columns passed, then smaller Σ|err|/tol. Mass outside the support is
  clamped to the end cells and the clamped fraction is always reported.
- **Mixture scoring**: exact weighted sufficient statistics (fraction
  vector, per-bin KE means, fate fractions at source-ion level with the
  scored share folded in) fed through the committed scorer via a duck-typed
  read; test-locked to equal literal pooling when weights ∝ ion counts.
  `χ²_med` intentionally absent (its sim-SE widening is an ion-count
  convention with no exact weighted analogue).
- **Densities**: standing = the five battery members' realized per-ion
  droplet radii; corrected = `legacy`+`raw` draw (200k samples, seed
  20260727) from the r3l2 cfg's own source conditions — drawn ⟨N⟩ 12769 vs
  nozzle-correlation 12794 (D0 §15 anchor reproduced).

### Oracle 1 (pre-registered, plan §4.3): INADMISSIBLE — 3/7

Gate (frozen before the run): per observable |recon − recorded| ≤
max(2 seed-SD, 10 % of recorded), seed-SDs measured from the 5 members;
adopted convention needs ≥ 6/7. Result: 3/7 (both conventions). Cause is
visible in the weights: **53.8 % of the standing density lies below the
R = 26.6 Å support edge** (the grid support starts at the standing *mean*),
clamps into r1l3, and overshoots trap (+0.027 ≈ 9 SD), n̄ (+0.58) and W₁
(+0.34 ≈ 3.6 SD) while undershooting supp (−0.026). midHot and deepKE pass.

### Oracle 2 (in-support; **declared post-hoc** after oracle 1 failed): 7/7

Target = the battery's own R ≥ 26.6 Å sub-ensemble (46.2 % of its ions;
support-covered by construction; scored with the committed conventions).
Both conventions pass 7/7; `nearest` errors: trap +0.002 (0.6 SD), supp
+0.002, n̄ −0.151 (1.3 SD), n₁_solv +0.008, W₁ −0.184 (1.9 SD), midHot
−0.018, deepKE −0.021 (1.0 SD). The target and the grid cells are
independent runs, so these errors contain interpolation bias *and* one
N = 500 seed's scatter — they are the stated interpolation error carried by
the forecast.

**Verdict structure:** the *method* (column re-weighting) is validated where
the support covers the density; the pre-registered failure is the support
hole below R = 26.6 Å. The corrected density is covered to ~95 % (clamp
0.2 % below / 5.1 % above R = 68.3 Å), so its forecast inherits the
in-support error scale. The below-support R ≈ 20 Å × L3 cell is the designed
§4.3 confirmation candidate, needed only if a standing-mixture
reconstruction ever becomes load-bearing.

### The forecast (nearest; Arm A/B = §G1.4 bracket at ensemble level)

| mixture | trap (marg) | det_yield | supp | n̄ | n₁_solv | W₁_solv | midHot | deepKE |
|---|---|---|---|---|---|---|---|---|
| pooled recorded | 0.067 | — | 0.187 | 4.07 | 0.243 | 0.571 | 1.011 | 0.631 |
| std × L3 (recon) | 0.094 | 0.91 | 0.161 | 4.65 | 0.224 | 0.912 | 1.04 | 0.599 |
| std × L2 | 0.005 | 1.00 | 0 | 9.29 | 0 | 4.67 | 1.252 (4 bins) | 1.47 |
| corr × L3 | 0.35–0.38 (0.032) | 0.62–0.65 | 0.12–0.13 | 6.0–6.7 | 0.15–0.16 | 2.07–2.69 | 1.254 | 1.35–1.45 |
| **corr × L2** | **0.31–0.42 (0.110)** | **0.58–0.69** | **0** | **13.9–14.7** | **0** | **9.0–9.8** | **1.256 (4 bins)** | **1.80–1.90** |

Reads (full statement in D0 §15.7): the corrected geometry **breaks the
landing at ensemble level through the pre-registered dressing/suppression
channel** (supp → 0, low-n evacuated, W₁ ≈ 9–10, trap 0.31–0.42); the
**birth law owns the histogram breakage** (std × L2 alone: W₁ 4.67 at
trap ≈ 0.005) while the **size distribution owns trapping** (corr × L3:
0.35–0.38); each axis alone pushes deepKE past 1 (ensemble echo of the
§G1.2/§G1.4 depth-crossing at ≈ 11–15 Å), together 1.80–1.90 — the
corrected geometry *overshoots* RQ11's deficit direction. det_yield
0.58–0.69 means the corrected detected ensemble is a small-R/shallow-birth
biased subset — a G2/G3-relevant structural fact, and possibly a physical
one (the experiment may itself only see the shallow subset).

**Caveats travelling with every forecast row:** N = 500 single-seed cells;
midHot at the L2 mixtures rests on 4 of 7 band bins; the marginal fraction
is conditional on the ~3× extrapolated cubic law (§G1.4); in-support
interpolation error as measured above.

**Atlas stance intact:** nothing adopted, `finc1v725` stands, F5
undischarged. G2 can now be taken on this forecast + the direct G1 rows.

## G3 Step 1 — twin landmark re-issue at the corrected geometry (2026-07-27, zero MD)

`scripts/tier2_h2b_forward_model.py g3landmarks` (new committed stage;
helper tests in `tests/test_tier2_h2b_forward_model.py`), executed under
its own `[PROCEED TO IMPLEMENTATION]` after G2. Three parts, oracle first;
outputs `h2b_g3_grid_twin{,_ke}.csv`, `h2b_g3_corrected_landmarks.csv`,
`h2b_g3_corrected_{row,ke}.csv` next to the committed twin CSVs.

### Oracles (all PASSED before any new number was read)

- **G3-P1 (S6 wiring):** all three committed `h2b_s6_final` rows
  (predictions **and** KE tables) re-derived **string-identically** from
  the leg-D draw discipline + the S6 scorer block — the reconstruction is
  certified against the frozen record.
- **G3-P2 (continuity anchor):** the recorded production center-pin
  landmark K = 0.74460 is reproduced — under the **pure-cubic** law, which
  is what it was recorded under (an earlier draft asserted it against the
  capped tail and failed at K = 0.48875; the landmark family is
  law-tagged from now on).
- **Test-level:** the twin's L2 arm (repo Boltzmann sampler at 313.2 K)
  reproduces the plan §3.1 pre-registered exposure-table mean depth
  29.3 Å at R = 34 Å.

### The re-issued center-pin landmarks (production kinematics, n₀ = 21)

| pin | N [He] | R [Å] | K_cubic (τ6.55) | K_capped (τ6.55) | t_exit [ps] | v_inf [Å/ps] | n_det |
|---|---|---|---|---|---|---|---|
| R2000 anchor | 2000 | 27.94 | **0.74460** | 0.48875 | 2.95 | 5.22 | 10 |
| r1 | 1727 | 26.60 | 0.68520 | 0.44618 | 2.69 | 5.56 | 9 |
| r2 | 3605 | 34.00 | 1.05191 | 0.72305 | 4.41 | 4.08 | 14 |
| r3 | 11059 | 49.40 | 14.93149 | 1.70360 | 10.04 | 2.33 | 20 |
| corr ⟨N⟩ | 12794 | 51.86 | 17.01518 | 2.04271 | 11.20 | 2.10 | 21 |
| r4 | 29227 | 68.30 | 20.43799 | 20.23274 | 21.85 | 2.35 | 21 |

Reads: (i) **the capped tail is what makes the anchored radii traversable
at all** — under pure cubic the R ≥ 49 Å center pin saturates the 150 ps
window (K 14.9–20.4 against the ≈ 22.9 window ceiling), under the cap the
fragment exits in 10–22 ps; (ii) even under the cap, R = 68.3 Å is in the
crawl regime (K_capped 20.2 — the fragment escapes but arrives at
n_det = 21, fully cold); (iii) the old center-pin landmarks
(K_prod 0.74460, K_9Å 0.89767) are ⟨N⟩ = 2000 numbers and are hereby
**superseded as ensemble anchors** — at the corrected mean size the
capped-law center-pin K is 2.04.

### Twin ↔ MD transfer at the G1 grid — the re-issued authority box

All 11 G1 cells run through the twin at the standing point (m = 20000 per
cell, CRN seed 20260727, density-tied dressing, MD rung tables); compared
against the MD rows of D0 §14.1 (committed conventions; twin trapped vs
the MD **t_b + t_m total**). Full table `h2b_g3_grid_twin.csv`; the
authority statements, which **supersede the standing-geometry
channel-(d) numbers for all G3 use**:

1. **Suppression/selection transfers near-quantitatively.** supp within
   0.02 of MD on every L3 cell (0.186/0.167, 0.147/0.140, 0.115/0.123,
   0.103/0.126) and exactly 0 on every L1/L2 cell, matching MD; n₁_solv
   within 0.016 on L3 and exactly 0 on L1/L2. Route B's selection
   observables are the twin's best-transferred ones.
2. **n̄ carries a residence-scaled hot bias.** Twin − MD: +0.2…+0.6 He at
   R1 → +1.4…+1.6 at R2 → +1.9/+2.8 at R3 → +2.5/+3.2 at R4 — the
   channel-(d) frozen-chord bias re-measured at long chords, same sign as
   the standing-geometry −0.32 (MD − twin) and growing ≈ linearly with
   in-droplet residence, as I73 predicted. W₁ inherits it
   (+0.26 → +3.2). **Any G3 twin scan must aim n̄ at target + bias, and
   twin-landing cells are conservative candidates (MD lands lower).**
3. **trap is a lower bound, and the center-pin trap channel is invisible.**
   Twin under-traps every off-center cell by 0.03–0.14 (e.g. r3l3
   0.305/0.404, r4l2 0.704/0.800) — consistent with no pickup
   mass-loading and no Landau freeze. At **r3l1 the twin traps 0.000
   against MD 0.538**: at an exact center pin the 1D chord has no
   inward-partner asymmetry, so the MD trap channel there is
   **mechanism-made** (pickup mass-loading during the long transit +
   relaxation-stage Landau off-switch), structurally outside the twin.
   Route B trap reads are therefore twin **floors**, not estimates.
4. **KE observables are direction-only, twin ≈ 15–30 % hot** (midhot_geo
   +0.06…+0.28 where both defined; deepKE +0.09…+0.44 except the
   band-edge r3l3/r4l3 at ±0.1). Orderings are preserved everywhere,
   including the deepKE = 1 crossing bracketed between r2l3 and r3l3 by
   both instruments.
5. The window-saturated K tail (K_q95 ≈ 44–46 at R ≥ 34 Å off-center
   cells) is the twin image of the MD retained class; the twin cannot
   decompose bound vs marginal.

### The corrected-ensemble twin row (finc1v725 parameters, corrected geometry)

Master draw m = 20000, seed 20260727, draw order N → births → cosines:
`legacy`+`raw` sizes at the preset's own 40 mbar / 14 K (drawn ⟨N⟩ 12750
vs nozzle correlation 12794; R q05/q50/q95 = 34.6/48.5/68.4 Å — the
parent's own quantiles 34/68.3 reproduced) × Boltzmann 313.2 K births
(depth q05/q50/q95 = 27.6/33.8/45.9 Å). Dressing saturates (n₀ ≡ 21).
K(τ3.2) q05/q50/q95 = 1.26/2.73/44.0; t_exit q50/q95 = 7.7/18.7 ps.

| | standing (twin) | corrected (twin) | D2b §4.3 corr × L2 forecast |
|---|---|---|---|
| trap | 0.0416 | 0.3080 | 0.31–0.42 |
| supp | 0.2019 | **0** | 0 |
| n̄_det | 4.387 | 16.99 | 13.9–14.7 |
| n₁_solv | 0.2429 | **0** | 0 |
| W₁_solv | 0.678 | 12.11 | 9.0–9.8 |
| midHot (geo) | 1.067 | 1.588 (1 bin) | 1.256 (4 bins) |
| deepKE | 0.898 | 2.119 (8 bins) | 1.80–1.90 |

**Cross-instrument agreement — the headline.** Two independent zero-MD
instruments (grid re-weighting over MD cells vs the twin forward model)
now give the same corrected-geometry picture, and **every twin−forecast
discrepancy carries exactly the sign and magnitude of the measured twin
bias above**: trap at the forecast's low edge (twin under-traps), n̄ high
by ≈ +2.5 (the residence bias), W₁ high by the same channel, KE hot by
the KE bias. The G2 dossier's ensemble forecast is thereby independently
confirmed: at the corrected geometry the standing parameters break the
landing through the dressing/suppression channel (supp → 0, low-n
evacuated), trap ≈ ⅓, and the KE observables overshoot hot.

**Caveats:** twin numbers, frozen-chord authority as boxed above; the
corrected midHot rests on 1 occupied band bin; no bound/marginal
decomposition; nothing here is an MD result.

### Consequences for the G3 scan (design inputs, not adjudication)

- The twin is **certified as the G3 scan instrument** with the
  per-observable authority box above; the scan reads Route B's selection
  observables at near-quantitative authority and Route A's cascade
  observables with the stated hot bias applied.
- The n̄ gap the scan must close is enormous at the standing parameters
  (twin 17.0 vs target ≈ 4.1) and n̄ moves toward the target as trap
  rises (detected-subset selection) — consistent with Route B carrying
  more of the load than Route A energetics alone.
- Landmark continuity: any future twin session at the corrected geometry
  oracles against `h2b_g3_corrected_row.csv` (bit-exact re-derivation,
  same seed/draw order) the way S6 sessions oracle against
  `h2b_s6_final`.

**Atlas stance intact:** nothing adopted, `finc1v725` stands, F5
undischarged; the G3 scan itself (Route A/B (v_c, τ, E₀)[+ E_bind] grid)
stays behind its own trigger.

## G3 Step 2 — the Route A/B twin scan at the corrected geometry (2026-07-27, zero MD, triggered)

Stage `g3scan` (plan §3.5c executed as frozen): 30 chord families
(v_c ∈ {5.0…10.0} ∪ {3.5, 4.25 diagnostic} × E_bind
∈ {0.0482, bundle 0.1168, 0.154} — exact extracted values, §6.7 item-2
tags; one integration of the committed corrected master each,
npz-cached) × 216 free cells (τ ∈ {2.4…12.8} × E₀ ∈ [0.17, 0.52]) =
**6480 scored cells** against the pre-registered hard gate
(n₁_solv ∈ [0.19, 0.30] ∧ n̄ ∈ [4.4, 7.1]). `integrate_pairs` gained a
byte-inert `e_bind_ev` override (array-equal-tested). Outputs:
`h2b_g3scan_{prescan,chords,predictions,gated_predictions,gated_ke}.csv`.

### Oracles (all PASSED before any new number was read)

- G3-P1: the three committed `h2b_s6_final` rows re-derived bit-exact
  (predictions + KE).
- Landmark continuity: `h2b_g3_corrected_row.csv` + `_ke.csv` re-derived
  string-identically at the standing cell (through the value-identical
  `_g3_corrected_ensemble` factoring — the refactor is oracle-covered).

### Block 1 — zero-integration pre-scan (Route-A kill criterion)

Uniform exposure scalings X → f·X of the *stored standing chord*
(f ∈ 0.05…1.0): the n₁_solv gate band is reachable at some f ≥ 0.25 in
**141/216** (τ, E₀) cells ⇒ **the analytic Route-A kill criterion did
NOT fire.** But the **full joint gate fires at 0** of all (τ, E₀, f)
combinations — under pure rescaling the standing chord's K-shape and its
frozen trap 0.308 block n̄. Registered reading: a scale factor cannot
land the corrected geometry; the landing must come from chord
*reshaping* (trap composition + K tail) — measured next.

### Block 2 — the verdict: **24/6480 cells gate; the basin exists and is Tier-0-legitimate**

The pre-registered failure criterion did **not** fire. Class breakdown:
all 24 in the tier0 v_c range (diagnostic arms: 0); 17 at τ ≤ the
sourced 6.55, 7 τ-flagged; by well: eb0482 16 / eb1168 5 / eb154 3.
The full gated set (twin numbers; trap is a floor, KE
direction-only, W₁ bias-loaded — §14.4 authority box applied):

| v_c | well | τ | E₀ | n₁_solv | n̄ | W₁ | trap | supp | midHotG | deepKE | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5.5 | eb0482 | 4.8 | 0.37 | 0.202 | 5.05 | 0.935 | 0.011 | 0.133 | 1.18 | 1.10 | E |
| 5.5 | eb0482 | 4.8 | 0.38 | 0.215 | 4.77 | 0.899 | 0.011 | 0.176 | 1.11 | 1.06 | E |
| 5.5 | eb0482 | 4.8 | 0.39 | 0.224 | 4.51 | 0.893 | 0.011 | 0.221 | 1.03 | 1.02 | E |
| 5.5 | eb0482 | 6.4 | 0.30 | 0.199 | 4.82 | 0.383 | 0.011 | 0.029 | 1.09 | 0.79 | E |
| 5.5 | eb0482 | 6.4 | 0.31 | 0.241 | 4.43 | 0.323 | 0.011 | 0.074 | 0.97 | 0.75 | E |
| **5.5** | **eb1168** | **4.8** | **0.36** | 0.193 | 4.88 | 0.592 | 0.053 | 0.092 | 1.15 | 0.80 | — |
| **5.5** | **eb1168** | **4.8** | **0.37** | 0.210 | 4.56 | 0.503 | 0.053 | 0.134 | 1.07 | 0.76 | — |
| 5.5 | eb154 | 4.8 | 0.36 | 0.196 | 4.65 | 0.641 | 0.076 | 0.092 | 1.09 | 0.66 | E |
| 6.0 | eb0482 | 6.4 | 0.33–0.37 (5 cells) | 0.20–0.21 | 4.49–5.59 | 1.48–1.71 | 0.025 | 0.12–0.32 | 0.39–0.52 | 0.61–0.70 | E |
| 6.0 | eb0482 | 9.6 | 0.27 | 0.253 | 4.65 | **0.259** | 0.025 | 0.034 | 0.42 | 0.42 | T E |
| **6.0** | **eb1168** | **6.4** | **0.32** | 0.208 | 4.84 | 0.744 | 0.129 | 0.081 | 0.49 | 0.47 | — |
| **6.0** | **eb1168** | **6.4** | **0.33** | 0.226 | 4.49 | 0.679 | 0.129 | 0.132 | 0.45 | 0.46 | — |
| 6.0 | eb154 | 4.8 | 0.40–0.41 (2 cells) | 0.19–0.19 | 4.41–4.66 | 1.13–1.15 | 0.186 | 0.19–0.23 | 0.47–0.49 | 0.57 | E |
| 6.5 | eb0482 | 9.6 | 0.30–0.32 (3 cells) | 0.19–0.20 | 4.51–5.17 | 1.16–1.18 | 0.035 | 0.14–0.25 | 0.27–0.31 | 0.42–0.45 | T E |
| 6.5 | eb0482 | 12.8 | 0.26 | 0.226 | 4.71 | 0.515 | 0.035 | 0.046 | 0.28 | 0.33 | T E |
| 6.5 | eb1168 | 9.6 | 0.28 | 0.194 | 4.46 | 0.991 | 0.213 | 0.047 | 0.30 | 0.53 | T |
| 7.25 | eb0482 | 12.8 | 0.31 | 0.200 | 4.49 | 0.877 | 0.040 | 0.193 | 0.22 | 0.33 | T E |

(flags: T = τ > 6.55 calibration-class flag; E = off-bundle-well
joint-pairing stamp. Bold = the four cells clean of both.)

**The fully-unflagged basin: (v_c 5.5, τ 4.8, E₀ 0.36–0.37) and
(v_c 6.0, τ 6.4, E₀ 0.32–0.33) at the standing well.** Versus finc1v725
all three knobs re-arbitrate in physically comfortable directions:
v_c 7.25 → 5.5–6.0 (softer mid-band γ — the NB-RQ11-12 lever
direction), τ 3.2 → 4.8–6.4 (*toward* the sourced GAH25 6.55), E₀
0.27 → 0.32–0.37 (inside the RQ1 band). The vc5.5/τ4.8 pair holds the
better midHot (1.07–1.15); the vc6.0/τ6.4 pair the better W₁ — the two
clean sub-basins trade the KE axis against the histogram axis.

### Mechanism reads

- **Cascade-carried, not selection-carried.** Gated families run
  det_yield 0.87–0.99 and the detected-subset R quantiles sit within
  1–2 Å of the source (e.g. vc5.5/eb1168: R_det_q50 47.8 vs src 48.5) —
  the Route-B droplet-size selection axis is barely exercised. The work
  is done by the chord law: K655_q50 1.334 → 0.60–0.73 plus the trap
  drop 0.308 → 0.05–0.13. **The Step-1 "Route B carries more of the
  load" expectation is hereby corrected** — n̄ closes through exposure
  + fate energetics at near-full yield, not through detected-subset
  R-selection.
- **The standing chord (v_c 7.25, bundle well) gates at 0/216** — no
  (τ, E₀) rescues finc1v725's drag point at the corrected geometry; it
  passes from under-strip (n̄ 17) to over-strip without entering the
  joint gate.
- **The sub-band diagnostic arms land nothing** (as-designed check):
  at v_c 3.5/4.25 the n₁ band is reachable only at n̄ 3.0–3.6, below
  the 4.4 floor — over-stripped. The corrected geometry does **not**
  point below the TDDFT band top; the collaborator-ask escape branch
  stays un-fired.
- **E_bind × v_c over-dissipation regime:** the trap ladder steepens
  with v_c — eb0482 keeps trap ≤ 0.043 at every v_c, while eb154 at
  v_c ≥ 8 reaches trap 0.55–0.61 with K655_q50 17–18 (the
  helium-coupled class explodes). The coupled fraction being strongly
  (v_c, well)-sensitive is the §3.5b distinguishing evidence pointing
  at the over-dissipation-artifact branch — direct input to the
  retained-policy sub-decision deferred at G2.
- KE (direction-only): the clean basin runs deepKE 0.46–0.80 and the
  eb0482/τ4.8 corner reaches ≈ 1.0; with the twin's +15–30 % hot bias
  the MD values will sit lower — the RQ11 deep-bin question transfers
  to the corrected geometry and must be read at the MD ring, not here.

### Consequences (reported, not adjudicated)

- Designed next step per §3.5c: the **MD confirmation ring**
  (~10–20 × N = 500 around the clean basin) behind its own
  `[PROCEED TO IMPLEMENTATION]`; the natural ring covers the two clean
  sub-basins ± one step in each knob at the standing well, plus one
  well-bracket cell.
- The **retained-policy sub-decision** (G2) and the **Axis-B deferral
  question** now have their scan evidence (over-dissipation coupling
  above; E₀/τ influence curves partially covered by the free-surface
  maps in `h2b_g3scan_predictions.csv`).
- **Atlas stance intact:** nothing adopted, `finc1v725` stands, F5
  undischarged; form authority remains Tier-0's (the scan arbitrates
  (v_c, τ, E₀, E_bind) *at fixed form* — capped cubic p_tail −1).

## G3 Step 3 — the MD confirmation ring (2026-07-27/28, 14 × N = 500, plan §3.5d)

Generator `gen_tier2atlas_g3ring.py`, scorer
`tier2atlas_g3ring_table.py` (table CSV committed:
`data/runs/h2b_forward_model/atlas_g3ring_table.csv`); corrected
geometry (legacy+raw / Boltzmann 313.2 K; realized R_q50 49.5 Å, mean
birth depth 34.6 Å), fresh shared seed 20260728, run dirs
`…tier2atlas_conf270_g3r*`. Launch provenance: two harness-managed
background launches were externally killed mid-E2 (no Python error;
cells died pre-`relaxation.npz` and were rebuilt); the completed ring
ran as a detached OS process. Both §1.4 oracles passed before any MD
number was read (GR-P1 twin rows string-exact; pooled-battery scorer
drift ≤ 0.002).

### Verdicts (all six pre-registered predictions land)

- **GR-P2 CONFIRMED — the twin basin is real in MD.** Cells landing
  the bias-free MD acceptance (n₁ ∈ [0.19, 0.30] ∧ n̄ ∈ [3.77, 4.37]):
  **a037** (v5.5, τ4.8, E₀0.37 — n̄ 4.097, within 0.03 of the 4.07
  target), **b031** (v6.0, τ6.4, E₀0.31), and beyond the strict A/B
  clause also **d030** (v5.5, τ6.4, E₀0.30) and **e154** (v5.5, deep
  well 0.154, τ4.8, E₀0.36). The MD basin spans both τ values at
  v5.5 and tolerates the full Tier-0 well spread.
- **GR-P3 CONFIRMED:** n̄ twin-hot at 14/14 cells, Δ 0.24–2.66 He —
  the small end of the +0.2…+3.2 bracket, residence-scaled (arm A
  ≈ 0.5, arm B ≈ 0.9, f725 2.66).
- **GR-P4 CONFIRMED:** n₁ transfer ≤ 0.036 at every cell (band 0.05) —
  the twin's n₁ authority is MD-grade at the corrected geometry.
- **GR-P5 SPLIT-CONFIRMED:** c65 fails on n₁ (0.165 < 0.19) as
  predicted; f725 reproduces the broken landing (n₁ 0, n̄ 14.3, trap
  0.437 ≥ 0.17); **c50 fails as predicted but on the *other* clause**
  — its n̄ 4.35 is inside the band and it misses on n₁ by 0.005
  (0.1849), so the basin nearly reaches v_c 5.0 (a §6.6-style
  narrow-gate lesson, this time caught by the pre-registered control).
- **GR-P6 CONFIRMED:** trap ordering 0.026 (eb0482) < 0.087 (eb1168)
  < 0.115 (eb154) at the basin chord, ≈ 0.8/eV — consistent with the
  §6.7 item-2 lever (0.85/eV).

### The table (committed CSV holds full columns)

| cell | (v_c, well, τ, E₀) | trap (bnd+marg) | supp | n̄ | n₁ | W₁ | midHot | deepKE | χ²_med | gate |
|---|---|---|---|---|---|---|---|---|---|---|
| a031 | 5.5, eb1168, 4.8, 0.31 | 0.086 (.084+.002) | 0.002 | 6.072 | 0.031 | 1.72 | 1.49 | 0.56 | 997 | 0 |
| a033 | 5.5, eb1168, 4.8, 0.33 | 0.087 (.085+.002) | 0.009 | 5.361 | 0.102 | 1.23 | 1.28 | 0.50 | 818 | 0 |
| a035 | 5.5, eb1168, 4.8, 0.35 | 0.087 (.085+.002) | 0.058 | 4.692 | 0.181 | 0.95 | 1.10 | 0.45 | 659 | 0 |
| **a037** | 5.5, eb1168, 4.8, 0.37 | 0.085 (.084+.001) | 0.137 | **4.097** | **0.211** | 0.83 | 0.94 | 0.43 | 317 | **1** |
| b029 | 6.0, eb1168, 6.4, 0.29 | 0.207 (.189+.018) | 0.005 | 5.061 | 0.082 | 1.39 | 0.62 | 0.61 | 78 | 0 |
| **b031** | 6.0, eb1168, 6.4, 0.31 | 0.204 (.188+.016) | 0.039 | **4.314** | **0.203** | 1.06 | 0.50 | 0.72 | 73 | **1** |
| b033 | 6.0, eb1168, 6.4, 0.33 | 0.197 (.185+.012) | 0.138 | 3.676 | 0.244 | 0.97 | 0.41 | 0.63 | 72 | 0 |
| c50 | 5.0, eb1168, 4.8, 0.34 | 0.023 (.023+0) | 0.048 | 4.346 | 0.185 | 1.00 | 1.92 | 1.06 | 1322 | 0 |
| c65 | 6.5, eb1168, 6.4, 0.33 | 0.322 (.256+.066) | 0.063 | 4.695 | 0.165 | 1.32 | 0.42 | 0.81 | 98 | 0 |
| **d030** | 5.5, eb1168, 6.4, 0.30 | 0.089 (.089+0) | 0.030 | **3.943** | **0.212** | 1.17 | 0.81 | 0.27 | 400 | **1** |
| d036 | 6.0, eb1168, 4.8, 0.36 | 0.194 (.179+.015) | 0.052 | 5.326 | 0.161 | 1.28 | 0.64 | 0.68 | 44 | 0 |
| e0482 | 5.5, eb0482, 4.8, 0.36 | 0.026 (.026+0) | 0.086 | 4.920 | 0.197 | 0.67 | 1.13 | 0.68 | 308 | 0 |
| **e154** | 5.5, eb154, 4.8, 0.36 | 0.115 (.115+0) | 0.092 | **4.140** | **0.208** | 0.99 | 0.96 | 0.43 | 279 | **1** |
| f725 | 7.25, eb1168, 3.2, 0.27 | 0.437 (.278+.159) | 0.000 | 14.34 | 0.000 | 9.45 | 1.40 | 1.76 | 140 | 0 |

### Readings

1. **finc1v725 is MD-measured broken at realistic droplets** (f725:
   44 % of ions helium-coupled at handover, zero suppressed, zero n₁,
   n̄ 14.3) — the twin/D2b forecast chain confirmed end-to-end in real
   MD (trap 0.437 vs the D2b Arm A–B 0.31–0.42 + the twin-floor bias).
2. **Retained-policy evidence (G4 input):** the marginal
   (modelling-exclusion) class is ≈ 0 at the basin (≤ 0.018, mostly
   exactly 0) and 0.159 at the standing chord — the §3.5b
   over-dissipation reading now MD-grade: the coupled-class question
   largely evaporates at the re-arbitrated point.
3. **KE trade-off persists (RQ11 successor):** no gated cell holds
   both KE axes — a037/e154 hold midHot ≈ 0.94–0.96 with deepKE
   ≈ 0.43; b031 holds deepKE 0.72 (better than the standing pooled
   0.631) at midHot 0.50; d030 is coldest on both. The mid-vs-deep
   tension moves with (v_c, τ) inside the basin, i.e. it remains a
   *drag-shape* question (NB-RQ11-12) at the corrected geometry, now
   posed at v_c 5.5–6.0.
4. **W₁ at gated cells 0.83–1.17** vs the standing pooled 0.571 —
   single-seed N = 500 numbers (±0.1–0.15 scatter), reported not
   gated; a pooled battery at the G4 winner is the honest W₁ read.
5. Suppressed fraction at gated cells 0.03–0.14 (vs standing pooled
   0.187) — the RQ3 mixture read moves but stays plausible; f725's
   supp 0 shows the corrected geometry kills the suppressed channel at
   the standing chord entirely.

### Consequences (reported, not adjudicated)

**G4 is now live**: the successor-point choice among the gated cells
(a037 / b031 / d030 / e154 — or a small interpolation battery), the
retained-policy final call (evidence above), and the ledger re-issue
(retire the five scaffolding rows) are the user's adjudications. A
pooled battery (§4cc pattern) at the chosen point is the natural
verification step before re-baselining. Atlas stance intact:
`finc1v725` stands until G4 is taken.

## G4 Step 1 — Block 0 (twin ranking authority) + Block 1 (the fine ridge scan) (2026-07-28, zero MD, plan §3.5e)

Report `tier2atlas_g4_transfer.py` (Block 0, pure scorer) and stage
`g4scan` in `tier2_h2b_forward_model.py` (Block 1). Committed artifacts:
`h2b_g4_transfer.csv`, `h2b_g4scan_{predictions,gated_ke,ridge}.csv`.
All oracles passed before any number was read: the frozen twin
pre-registration, the pooled-battery scorer-drift oracle, the S6
machinery oracle, the corrected-row landmark oracle, and **G4-P1 — 408
shared sub-lattice cells reproduce the committed Step-2 values
string-exactly**.

### Block 0 — the twin's ranking authority, measured for the first time

The 14 paired (twin, MD) ring cells read as a transfer measurement.
Pre-registered permission gate ρ ≥ 0.7:

| observable | ρ (n = 14) | 95 % CI | verdict |
|---|---|---|---|
| W₁_solv | +0.824 | [+0.44, +0.98] | LICENSED |
| midHot | +0.996 | [+0.93, +1.00] | LICENSED |
| deepKE | +0.327 | [−0.28, +0.82] | **NOT LICENSED** |

**The deep-KE axis is not twin-scannable.** Every RQ11 statement must be
MD-measured; no twin scan, at any resolution, can rank it.

Level transfer is a **regression**, not an offset (13 cells, f725
dropped): `MD_W₁ = 0.671 + 0.412·twin_W₁` (R² 0.60, resid SD 0.18) and
`MD_midHot = 0.001 + 0.897·twin_midHot` (R² 0.994, resid SD 0.035). The
mean W₁ shift is +0.03 over all cells but **+0.32 ± 0.10 over the four
MD-gated ones** and negative at the bad ones. Two consequences carried
into the ranking: the W₁ intercept 0.671 is a predicted MD floor *above*
the standing pooled 0.571; and because |ln x| is V-shaped about 1, a
multiplicative midHot bias does not preserve the ordering of raw twin
values even at ρ = +0.996 — cells must be ranked on the corrected value.

n̄ bias model `Δn̄ = a + b·n̄_twin`, **residual SD 0.28 He**, robust to
dropping f725. In-sample: a037 predicted 4.06 vs MD 4.097, d030 3.91 vs
3.943, e154 4.13 vs 4.140, b031 4.56 vs 4.314.

Seed-SD normalizers (five N = 1000 battery members): W₁ 0.5791 ± 0.0954,
midHot 1.0108 ± 0.0222, deepKE 0.6327 ± 0.0209, χ²_med 131.5 ± 15.4.
Two **plan §1.2 provenance corrections**: its "W₁ SD ≈ 0.04" is the SEM
of the pooled mean (0.0954/√5 = 0.043), not the per-seed SD — so an
N = 500 single-seed W₁ carries ≈ 0.13 and the G3 ring never distinguished
its cells on W₁ at all; and its deep-bin read "0.0603 ± 0.0033" does not
reconcile with the deepKE ratio (0.633 ± 0.021), most likely being the
deep-bin population weight. Flagged, not overwritten.

Score pre-registration reproduces exactly: pooled battery
S_provisional = **4.366** vs the design's 4.37.

### Block 1 — the fine ridge (5544 cells: v_c step 0.25, τ step 0.4, E₀ step 0.005)

- **574/5544 gate (172 clean of both stamps)** against 24/6480 at Step 2.
  The cause is **the gate, not the grid**: Step 2 carried the crude
  [+0.3, +3] bias bracket (twin n̄ ∈ [4.4, 7.1]); the Block-0 model gates
  on predicted MD n̄, i.e. twin n̄ ∈ [3.8, 5.4] — *narrower*, but placed
  where the n₁ band actually lives. n₁ and n̄ are anti-correlated under
  stripping, so the old n̄ floor of 4.4 was fighting the n₁ band. **The
  binding constraint at Step 2 was the instrument's error model.**
- **G4-P2 CONFIRMED — connected ridge.** The clean gated set is one
  4-neighbour-connected diagonal band from (5.5, 4.8) to (6.0, 6.4):

  | v_c \ τ | 4.0 | 4.4 | 4.8 | 5.2 | 5.6 | 6.0 | 6.4 |
  |---|---|---|---|---|---|---|---|
  | 5.0 | 5 | 6 | 4 | 3 | 2 | 2 | 1 |
  | 5.25 | 2 | 7 | 6 | 5 | 4 | 3 | 2 |
  | 5.5 | 0 | 7 | 8 | 6 | 5 | 4 | 3 |
  | 5.75 | 0 | 2 | 10 | 9 | 8 | 6 | 5 |
  | 6.0 | 0 | 0 | 0 | 10 | 10 | 9 | 8 |
  | 6.25 | 0 | 0 | 0 | 0 | 0 | 3 | 7 |

  (clean gated cells per lattice point). The two Step-2 "sub-basins" were
  opposite corners of one ridge — the registered hypothesis holds.
- **Clean optimum v_c 5.5, τ 4.8–5.2, E₀ 0.34–0.37.** a037's corner is
  near-optimal but slightly off: the best clean cells are
  (5.5, 4.8, **0.365**) S_pred 1.59 and (5.5, 5.2, **0.34**) 1.60.
- **G4-P4 NOT MET on the licensed axes.** Best clean S_pred 1.59
  (per-cell uncertainty ≈ 0.57) vs the incumbent's measured **1.091**;
  best including the off-bundle-well stamp is (5.5, eb0482, 5.6, 0.33) at
  1.47. In-sample predicted/measured S at the four ring cells:
  1.82/1.90, 6.89/6.83, 2.66/3.52, 1.79/2.05.
  *(The stage's own printed "G4-P4 MET, best S 1.506" is the raw-twin,
  three-axis provisional score — superseded here by the bias-corrected
  two-axis reading, which is the defensible one.)*
- **G4-P3 NOT EVALUABLE.** 75 gated cells hold midHot ∈ [0.85, 1.15] ∧
  twin deepKE ≥ 0.6, but deepKE is the axis Block 0 measured as
  rank-untransferable (a037: twin 0.755 → MD 0.43). The pre-registered
  Block-2 (p_tail) trigger cannot be decided by this scan and **moves to
  the MD finalists**. A pre-registration defect found by G4's own Block 0;
  recorded, not reinterpreted.

### Open confound (Block 3 settles it)

The ring cells are N = 500 (≈ 900 scored ions), the incumbent reference is
pooled N = 5000. Sampling noise inflates W₁, so part of the 0.83–1.17 vs
0.571 gap may be statistical rather than physical. The five N = 1000
battery members (W₁ 0.48–0.73, mean 0.579) bracket the pooled 0.571, so
the N = 1000 → 5000 step is small; the N = 500 → 1000 step is untested.
The mandatory a037/b031 N = 1000 replicates are the control.

### Consequences (reported, not adjudicated)

Block 2 (p_tail) neither fires nor is dismissed — it is deferred to the MD
evidence. Block 3 (6 × N = 1000) is the designed next step and carries its
own trigger. Atlas stance intact: nothing adopted, `finc1v725` stands.

### Block 3 — the MD finalists (2026-07-28, 6 × N = 1000, seed 20260729)

Generator `gen_tier2atlas_g4finals.py`, scorer
`tier2atlas_g4finals_table.py`, table `atlas_g4finals_table.csv`. Both
oracles passed before any number was read (frozen twin rows string-exact
from the committed fine scan; pooled-battery scorer drift). Wall clock
~90 min at concurrency 3 (two batches of three).

| cell | (v_c, τ, E₀) | n̄ | n₁ | W₁ | midHot | deepKE | χ²_med | trap (marg) | gate | S |
|---|---|---|---|---|---|---|---|---|---|---|
| f1 | 5.5, 4.8, 0.365 | 4.169 | 0.1878 | 0.821 | 0.990 | 0.523 | 411 | 0.091 (.003) | 0 | 1.511 |
| f2 | 5.5, 5.2, 0.34 | 4.194 | 0.1876 | 0.916 | 0.980 | 0.438 | 440 | 0.093 (.003) | 0 | 1.746 |
| **f3** | 5.5, 4.4, 0.395 | 4.205 | 0.1867 | **0.769** | **1.023** | **0.614** | 426 | 0.089 (.002) | 0 | 1.512 |
| **a037** | 5.5, 4.8, 0.37 | 4.024 | 0.1984 | 0.792 | 0.960 | 0.526 | 360 | 0.091 (.003) | **1** | 1.683 |
| b031 | 6.0, 6.4, 0.31 | 4.263 | 0.1771 | 1.113 | 0.520 | 0.732 | 68 | 0.203 (.026) | 0 | 6.631 |
| x345 | 5.5, 4.8, 0.345 | 4.782 | 0.1491 | 0.966 | 1.155 | 0.499 | 603 | 0.093 (.003) | 0 | 2.723 |

(S = the §3.5e joint score on the Block-0-licensed axes W₁ + midHot,
provisional norms; the incumbent `finc1v725` scores **1.091**.)

**Verdicts.**

- **G4F-P2 REFUTED — but by 0.001–0.003.** None of f1/f2/f3 gates: all
  three land n₁ 0.1867–0.1878 against the 0.19 floor. The cause is a
  design gap, not a physics failure: Block 1 bias-corrected n̄ but treated
  n₁ as unbiased on its ≤ 0.036 tolerance. That tolerance is **scatter
  with inconsistent sign** — MD ran 0.011–0.015 *below* twin here and
  +0.001…+0.026 *above* twin at the G3 ring. A twin n₁ gate needs a
  ±0.02 band, and a ridge optimum picked on twin n₁ sits ≈ half an E₀
  step too low.
- **G4F-P3 — the W₁ deficit is REAL.** a037 N = 500 → 1000: W₁
  0.826 → 0.792 (Δ −0.034 against a per-seed SD of 0.095); b031 moved the
  other way (+0.052). The corrected geometry sits at W₁ ≈ 0.77–0.92
  against the incumbent's 0.579, consistent with the transfer
  regression's 0.671 intercept. **The "lands worse" reading is confirmed
  and localized to the histogram axis alone.**
- **G4F-P4 — the mid-vs-deep KE tension is BROKEN.** Twin deepKE
  0.691/0.776/0.887 → MD 0.438/0.523/0.614 (f2/f1/f3), ordering
  preserved, span 0.176 at near-identical n̄/n₁. **f3 holds midHot 1.023
  and deepKE 0.614 at once** — the incumbent's own KE profile
  (1.011/0.633). The deep-KE axis is parameter-accessible on the
  (v_c, τ, E₀) surface at fixed drag form, so the G3-ring reading "no
  gated cell holds both axes" was a resolution artifact of its coarse E₀
  ladder, and **Block 2 (p_tail) is not required by this evidence**. The
  twin ranked f3 *worst* of the three and MD ranks it best — the
  unlicensed-deepKE result in operation.
- **G4F-P5 CONFIRMED:** x345 missed on n₁ (0.149 < 0.19) as
  pre-registered, with n̄ 4.78 above the ceiling.
- **a037 replicates** its ring landing on a fresh seed at double N
  (n̄ 4.024, n₁ 0.1984) — the basin is seed- and N-robust.
- **b031's corner is out:** trap 0.203, midHot 0.520, S 6.63 at N = 1000.
- **Retained policy (G4 input):** at the ridge the marginal
  (modelling-exclusion) class is **0.002–0.003**, i.e. empirically empty;
  only b031 shows 0.026. The G3-ring reading holds at N = 1000.

**Net.** At the corrected geometry the residual deficit is **one axis —
W₁ — not a KE trade-off**. No gated cell beats the incumbent's 1.091
(best: f1 1.511 ≈ f3 1.512, a037 1.683). The measured n₁ bias plus f3's
KE profile point the successor at **E₀ ≈ 0.38–0.41 on the (v_c 5.5,
τ 4.4–4.8) chord**, which is untested — 2–3 cells, ~45 min.

**Atlas stance:** nothing adopted; `finc1v725` stands; the G4
successor-point choice, the retained-policy final call and the ledger
re-issue remain user adjudications.

### Block 3 ladder arm — the E₀ direction (2026-07-28, 9 further cells × N = 1000, seed 20260729, CRN-paired)

User-requested extension of Block 3 along the E₀ direction the six
finalists pointed at. Placed on the **measured** transfer (predicted MD n̄
inside [3.80, 4.06]; twin n₁ ≥ 0.208 to clear 0.19 after the measured
−0.012 offset), not on twin values. Registered post-hoc — the cells were
chosen after the first six were scored, and the generator says so.

| cell | v_c | τ | E₀ | trap (marg) | supp | n̄ | n₁ | W₁ | midHot | deepKE | χ²_med | gate | S |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **h405** | 5.5 | 4.4 | 0.405 | 0.087 (.002) | 0.183 | 3.955 | 0.2078 | **0.709** | 0.957 | 0.603 | 350 | **1** | **1.559** |
| h410 | 5.5 | 4.4 | 0.41 | 0.086 (.002) | 0.200 | 3.862 | 0.2119 | 0.674 | 0.931 | 0.570 | 248 | 1 | 1.696 |
| h415 | 5.5 | 4.4 | 0.415 | 0.086 (.002) | 0.221 | 3.753 | 0.2147 | 0.669 | 0.906 | 0.548 | 261 | 0 | 1.878 |
| h375 | 5.5 | 4.8 | 0.375 | 0.089 (.002) | 0.146 | 3.889 | 0.2141 | 0.747 | 0.926 | 0.585 | 331 | 1 | 1.861 |
| h380 | 5.5 | 4.8 | 0.38 | 0.089 (.002) | 0.170 | 3.769 | 0.2202 | 0.733 | 0.895 | 0.577 | 283 | 0 | 2.075 |
| h345 | 5.5 | 5.2 | 0.345 | 0.092 (.003) | 0.093 | 4.043 | 0.1942 | 0.892 | 0.945 | 0.432 | 373 | 1 | 1.965 |
| h350 | 5.5 | 5.2 | 0.35 | 0.091 (.003) | 0.110 | 3.877 | 0.2095 | 0.884 | 0.897 | 0.441 | 320 | 1 | 2.328 |
| h355 | 5.5 | 5.2 | 0.355 | 0.090 (.001) | 0.133 | 3.732 | 0.2224 | 0.872 | 0.859 | 0.500 | 299 | 0 | 2.619 |
| v525 | 5.25 | 4.4 | 0.39 | 0.048 (.000) | 0.148 | 3.994 | 0.2084 | 0.727 | **1.394** | **0.656** | 912 | 1 | 3.648 |

**1. The E₀ lift works exactly as designed.** Five of the nine cells gate,
converting the three 0.001–0.003 near-misses into real landings. Per
+0.005 eV of E₀, inside the ridge: **W₁ −0.02…−0.04, n₁ +0.004…+0.007,
n̄ −0.11…−0.13 He, midHot −0.02…−0.03**, monotone in all three arms. Each
τ arm gates over a window of only ~0.01 eV in E₀.

**2. The W₁ floor is ≈ 0.67 and is now measured twice, independently.**
The τ 4.4 arm flattens: 0.769 → 0.709 → 0.674 → **0.669**. Block 0's
twin→MD regression predicted the floor as its intercept **0.671**, from
completely different data. The incumbent's 0.579 is unreachable on the
(v_c, τ, E₀) surface at the corrected geometry.

**3. Why the floor exists — n₁ and n̄ are not simultaneously matchable.**
W₁ falls because n₁ climbs toward the experimental 0.243, and n̄ falls
with it. On the measured τ 4.4 slopes, n₁ = 0.243 needs E₀ ≈ 0.435, where
n̄ ≈ 3.3 — **0.77 He below the experimental 4.07**. The standing point hit
n₁ 0.243 *and* n̄ 4.068 together only at the wrong geometry. So the
corrected-geometry deficit is a **shape** mismatch in the detected size
distribution, not a scale error, and it localizes to the mechanism
(ladder shape, pickup, per-shed ε) rather than the drag surface — the
branch §3.5c's failure criterion named in advance.

**4. τ ordering is the reverse of the twin's.** On W₁: τ 4.4 (0.669–0.709)
< τ 4.8 (0.733–0.792) < τ 5.2 (0.872–0.892), while twin W₁ ranked τ 5.2
best (0.42–0.49). Second inversion after f1/f3. **Twin W₁ cannot
discriminate cells separated by ≲ 0.05**; its ρ = 0.824 licence was
earned on cells spanning W₁ 0.5–12. Recorded as a limit on the
twin-first cost ladder, not a scoring fix.

**5. Deep-KE ceiling (v525).** deepKE **0.656**, above the incumbent's
0.633, in a cell that *gates* (n₁ 0.208, n̄ 3.994, W₁ 0.727) — but midHot
1.394, χ²_med 912. Lowering v_c to 5.25 buys deep-bin KE and pays in
midHot. RQ11 now has a measured ceiling: ≈ 0.66 at midHot ≈ 1.4, or
≈ 0.60 at midHot ≈ 0.96.

**6. Successor candidate: h405** (v_c 5.5, τ 4.4, E₀ 0.405). Gated, best
gated S (1.559), and better than a037 on **every** axis: W₁ 0.709 vs
0.792, deepKE 0.603 vs 0.526, midHot 0.957 vs 0.960, n₁ 0.208 vs 0.198.
Its supp 0.183 lands on the standing 0.187 (RQ3 mixture read intact),
trap 0.087, marginal class 0.002.

**7. n̄ model validity:** |pred − MD| ≤ 0.08 He inside the ridge
(systematically −0.04…−0.08 at each arm's high-E₀ end), but −0.30 at b031
and +0.11 at v525/x345 — an in-ridge instrument, not a global one.

**Atlas stance:** nothing adopted; `finc1v725` stands. The G4
adjudications (successor point, retained policy, ledger re-issue) are the
user's; the natural verification at h405 is a pooled 5 × N = 1000 battery.

## G4 Step 2 — Block F: mechanism fingerprints (FROZEN before any Block-D read) (2026-07-28)

Plan §3.5f Block F. Derived from the built mechanism (`physics/evaporation.py`,
`physics/pickup.py`, `internal_energy_budget`) and the committed influence
records (D0 §7/§8/§11, I42/I63/I65/I79/I93/I94, NB-RQ23-1/2) — **no Block-D
number had been computed when this section was committed** (the commit
timestamp is the freeze evidence; the anatomy scorer had not yet run).

Common derivation target (from the committed G4 Step 1 record, plan §3.5f):
at matched n₁ (0.243, reached near E₀ ≈ 0.435 on the τ 4.4 arm) the model
runs n̄ ≈ 3.3 vs the experimental 4.07 — and at h405 (n₁ 0.208, n̄ 3.955)
**both** ends are light while the E₀ lever moves them in opposite
directions. The missing weight must therefore enter as *dispersion* (heavier
n = 1 **and** heavier n ≥ 5–6 shoulder relative to the n = 2–4 core), not as
a slide of the cascade-depth coordinate. Each fingerprint below answers
whether its knob can do that.

### F1 — Ladder shape (c1 grading; RQ4)

1. **Mechanism route.** Terminal n is where the RRK band
   `D₀(n) < E_int < Σ(n)` closes against Newton cooling: the rung table
   `D₀(n)` (rq4graded = Form-U with rungs 1–3 × (2.2, 1.5, 1.3)) is the
   *per-bin* gate. n = 1 is a one-rung-wide E_ej window — the taper owns
   n₁ (I42); raising a mid rung (n ≈ 4–8) closes the gate earlier there
   and stalls part of the flux at that n. The deep tail (n ≳ 10) is
   Σ-locked (I93: ±5 % deep rungs ↔ ∓3.5 pts suppressed weight) and is
   NOT part of this knob's playground.
2. **Sign on n₁_solv:** − for a bare mid-rung raise (it starves the flux
   reaching n ≤ 3); ~0 for a *compensated* move (mid rungs up, taper
   retuned) — the two sub-knobs are separable by I42's measured n₁
   ownership.
3. **Sign on n̄_det:** + (stalled flux accumulates at the raised rungs).
4. **Bin-pattern:** lowers F_sim over n ≈ 2–5 (mass moved upward into the
   raised-rung bins), leaves F_sim at n = 1 controlled independently by
   the taper, leaves the n ≥ 10 tail pinned (Σ-lock). Distinctive
   signature: a *localized* CDF-gap change whose edges sit at the edited
   rungs — the only knob with per-bin-resolved leverage.
5. **§3.5f constraint (add shoulder mass without paying n₁)?**
   **CONDITIONAL YES** — uniquely among the three. A compensated two-part
   move (raise mid rungs, retune the taper) adds n ≥ 5 weight while
   holding the n = 1 window, i.e. it *creates dispersion*. Measured
   support: the ladder is KE-neutral (midHot flat across the family,
   I79), so it would not disturb the h405 KE landing — the axes that are
   already right stay right. Caveat: the compensation is exactly the
   RQ4-external-arbitration territory (I42) — physics authority for the
   taper remains open.

### F2 — Pickup (λ₀)

1. **Mechanism route.** `λ_attach = λ₀ · (ρ_He/ρ_bulk) · (1 − n/n*)₊^p`
   (p = 1 standing): the occupancy cap makes pickup *anti-select* high n
   — it refills the stripped low-n end while the ion is still in dense
   helium, and each fire deposits only `+f_ret·D₀(n+1)` (f_ret 0.1) of
   S1 heat. Measured: live refilling shifts suppression by ≈ 1–3
   ions/100 at p = 1 (I63/I65).
2. **Sign on n₁_solv:** − (the refill flux feeds precisely on the n = 1–2
   population; the occupancy factor is largest there).
3. **Sign on n̄_det:** + (refilled ions land at n ≈ 2–5); supp − (the
   measured de-suppression side effect).
4. **Bin-pattern:** lowers F_sim at n = 1–3 by *draining those very bins*
   into n ≈ 3–6; no n ≥ 8 leverage (late-flight ρ_He ≈ 0 gates pickup
   off; the occupancy cap suppresses it at high n even in-droplet).
5. **§3.5f constraint?** **NO.** Pickup adds shoulder mass by consuming
   the n = 1–2 bins directly — it pays n₁ one-for-one and is therefore
   another slide along the measured n₁↔n̄ anti-correlation, with a supp
   side effect that would break the landed supp 0.183 ≈ 0.187 (the RQ3
   mixture read h405 preserves).

### F3 — Per-shed ε (RQ2)

1. **Mechanism route.** The K1 drain per shed becomes −(D₀(n) + ε)
   (currently ε = 0, adjudicated 2026-07-10; the Klots–Hansen
   prescription is ε ≈ c·D/G ≈ 0.5–0.6 meV, finite-heat-bath corrected
   *smaller*, NB-RQ23-1/2). ε > 0 shortens every cascade by the same
   relative amount — a uniform re-scaling of the cascade-depth
   coordinate, mechanically degenerate with (−E₀, +τ) on the histogram.
2. **Sign on n₁_solv:** − (shorter cascades reach n = 1 less often).
3. **Sign on n̄_det:** + (every terminal shifts up together).
4. **Bin-pattern:** a near-uniform rightward translation of the whole
   solvated distribution — the same signed-gap pattern as the measured
   E₀ lever run backwards (per −0.005 eV: n̄ +0.11…+0.13, n₁
   −0.004…−0.007), with **no dispersion change**. No bin-localized
   structure: ε carries no n-dependence the ladder doesn't already own.
5. **§3.5f constraint?** **NO** — twice over. (a) Shape: ε is
   E₀/τ-degenerate on the histogram, and the E₀ direction is exactly the
   lever the ladder arm measured to a floor — a degenerate knob cannot
   break a floor. (b) Magnitude: the physically allowed ε (≈ 0.5–2 meV,
   D0 §11: 6–8 meV already breaks midHot) moves n̄ by ≪ the 0.77 He the
   deficit requires.

### The frozen match rule and decision table (verbatim from plan §3.5f)

**Match rule:** a knob is a live candidate iff its signature moves, with
the correct sign, the bins carrying ≥ 70 % of the residual W₁, without a
wrong-signed prediction on n₁ or n̄.

**Decision table:**

- exactly one knob matches → it becomes the next designed axis (own plan
  section, twin-authority caveat carried: twin W₁ cannot discriminate
  ≲ 0.05, so that axis budgets MD from the start);
- multiple match → cheapest-first discrimination (analytic/twin
  fingerprint sharpening) before any MD;
- none match → the W₁ floor is recorded as the corrected geometry's
  **honest residual** (D0 §17 + findings), the G4 adjudications
  complete, and the tier proceeds to the open items (RQ3/RQ5 reads,
  margin-3 Å pin I88, D2b A/B remainder, then Tier-3).

**Pre-registered reading aid (falsifiable):** if the Block-D residual at
the n₁-matched end is a *dispersion* pattern (F_sim too high over the
n = 2–4 core, too low at n = 1 and over the n ≥ 5–6 shoulder, tail
pinned), F1 (ladder shape) matches and F2/F3 do not — their patterns
cannot lower the core without paying an end. If instead the residual is
a *uniform-shift* pattern, **no knob matches** (the E₀ lever already
exhausts that direction at its measured floor) and the honest-residual
branch fires.

Frozen before any Block-D number was read.

## G4 Step 2 — Block V: the h405 pooled battery (2026-07-28, 5 × N = 1000, seeds 20260730–34)

Generator `gen_tier2atlas_g4step2_battery.py` (each member =
`build_cell(h405)` with only the seed swapped; cfg-vs-committed-`g4fh405`
oracle: exactly `{"seed"}`), scorer
`tier2atlas_g4step2_battery_table.py`, artifact
`atlas_g4step2_battery.csv`. All oracles passed before any number was
read (h405 twin row string-exact; scorer-drift; committed Block-3 row to
4 decimals), on every launch.

**Execution event:** the first launch lost 4 of 5 members to disk-full
(T: at 100 %, `OSError 28`); recovered by the **user-approved trajectory
strip** — `ion.npz`/`relaxation.npz` deleted from the oldest 194 scored
tier-2 run dirs (detection-stage runs only; Tier-0/1a trajectory runs,
the standing battery members and the new members untouched; manifest
retained in the session scratchpad; 20 GB freed). Every committed scorer
re-verified after the strip (finals table + twin oracles bit-exact).
Relaunch on the same seeds — deterministic, no physics impact.

| member | seed | trap (marg) | supp | n̄ | n₁ | W₁ | midHot | deepKE | χ²_med | gate | S |
|---|---|---|---|---|---|---|---|---|---|---|---|
| s1 | 20260730 | 0.087 (.001) | 0.203 | 3.826 | 0.214 | 0.801 | 0.934 | 0.576 | 237 | 1 | 1.889 |
| s2 | 20260731 | 0.068 (.000) | 0.192 | 3.819 | 0.221 | 0.729 | 0.943 | 0.539 | 305 | 1 | 1.695 |
| s3 | 20260732 | 0.073 (.003) | 0.191 | 3.892 | 0.199 | 0.798 | 0.948 | 0.524 | 271 | 1 | 1.782 |
| s4 | 20260733 | 0.076 (.001) | 0.200 | 3.906 | 0.195 | 0.763 | 0.957 | 0.506 | 328 | 1 | 1.651 |
| s5 | 20260734 | 0.087 (.001) | 0.185 | 3.940 | 0.218 | 0.736 | 0.949 | 0.498 | 438 | 1 | 1.663 |
| **pooled** | — | 0.078 (.001) | 0.194 | **3.877** | **0.209** | **0.765** | **0.946** | 0.504 | 461 | **1** | **1.734** |

**Verdicts.**

- **GV-P1 CONFIRMED** — pooled gates (n₁ 0.2094 ∈ [0.19, 0.30], n̄ 3.877
  ∈ [3.77, 4.37]); every member gates individually. h405 is
  seed-robust on the gate.
- **GV-P2 CONFIRMED** — pooled W₁ 0.7653 ∈ [0.64, 0.78]: **the floor,
  measured a third time independently** (arm asymptote 0.669; Block-0
  intercept 0.671; now the pooled 0.765 at the gated point). Per-seed SD
  0.0337 — 2.8× tighter than the standing battery's 0.0954: the
  corrected geometry is markedly more seed-stable, so N = 1000 W₁
  differences ≥ ~0.07 are now 2σ-significant.
- **GV-P3 REFUTED** — pooled S 1.734 > a037's 1.683. The Block-3
  single-seed 1.559 was a favorable draw (seed 20260729: W₁ 0.709 vs
  the battery's 0.729–0.801). Caveat recorded: 1.683 is itself a
  single-seed number, so this comparison carries ± the same scatter —
  but the pre-registered clause was pooled < 1.683 and it did not hold.
  **No adjudication fires automatically; the successor question returns
  to the user with this evidence.**
- **Low-n KE axis (pre-registered mid-stage, plan §3.5f):** pooled KE₁
  0.641 ± 0.003 eV (ratio 0.568 vs ref median 1.128), KE₂ 0.549 ± 0.002
  (ratio 0.778 vs ref mean 0.706); mode 0.656/0.569. The deficit is
  physics, not noise (per-seed SD ≈ 0.003 eV). Score terms would add
  **5.84** to S — the largest single defect on the surface. **Provenance
  note on the frozen gate band:** the registered 2×seed-SD half-width
  (ratio ± 0.005/0.007) is far narrower than the reference's own
  systematic bands (calib 4 %, condition 6 %); recorded as the
  *instrument resolution*, while any operational gate should carry the
  reference systematics (≈ ±10 %) — flagged, not overwritten (the plan
  §1.2 correction precedent).

## G4 Step 2 — Block D: the W₁ residual anatomy (2026-07-28, zero MD)

Scorer `tier2atlas_g4step2_w1_anatomy.py`, artifact
`atlas_g4step2_w1_anatomy.csv`. Oracles: scorer-drift; committed arm W₁
to 4 decimals; **identity Σ|gap| == w1_solvated to 1e-12 on all five
distributions**. Block F was committed (8580a8f) before any number here
was read.

**Provenance correction (flagged, not overwritten).** The Step-1
narrative's "experimental n₁ 0.243 / n̄ 4.07" are the **incumbent's sim
values** (the scorer-drift oracle numbers), not the reference: the
committed abundance reference gives **solvated n₁ 0.3103 and solvated
n̄ 4.889** (raw n = 1 fraction 0.1753 with bare 0.4352 — the 0.243-like
number arises only on the full support incl. the RQ3 bare mixture bin).
The §3.5d/e gate bands are therefore *incumbent-anchored conventions*,
not experimental readings. W₁ and the floor are unaffected (always
computed against the true reference); the "n₁ and n̄ not simultaneously
matchable" conclusion **strengthens** (the true targets are further out
on both ends).

**The residual shape (pooled h405, gap = F_sim − F_ref).** Deficit at
**both ends**, excess in the core: PMF n = 1 −0.101; core n = 2–12
**+0.156** (peak CDF gap +0.056 at n = 12); tail n ≥ 14 −0.055 (ref
0.0757 vs sim 0.0203). Top-70 % bins {1–4, 10–14}. The detected solvated
distribution is **under-dispersed** — mass must leave the n = 2–12 core
for both n = 1 and n ≥ 14 at once. The wrong-geometry incumbent shows
the same signature smaller (n = 1 −0.067; its shallow geometry was
partly *supplying* the dispersion).

**The E₀-lever signature (CRN-paired arm, per +0.005 eV):** closes gaps
only at n = 1–4 (+0.002…+0.012); n ≥ 5 sign-unstable < 0.003; **bins
7–20 are E₀-inaccessible** and carry ≈ 45 % of the residual — the
measured identity of the floor's carriers (D0 §4 updated).

**The frozen match rule, applied.**

- **F2 pickup: NO MATCH** (as pre-registered) — wrong sign at n = 1
  (drains the bin the residual needs filled).
- **F3 per-shed ε: NO MATCH** (as pre-registered) — uniform shift,
  wrong-signed at one end whichever direction it runs; E₀-degenerate.
- **F1 ladder shape: PARTIAL, not a match under the rule.** The taper /
  mid-rung sub-knobs can close the n = 1–4 side (correct-signed,
  ≈ 40 % of the residual, KE-neutral per I79) — but the n ≥ 14 tail is
  **Σ-locked** (I93: deep-rung edits pay in supp/midHot), so F1 cannot
  source the tail, and its achievable move is wrong-signed on n̄ unless
  paired with an independent tail source. Under the frozen ≥ 70 % rule:
  **no knob matches.**
- **Decision-table branch: NONE → the W₁ floor is recorded as the
  corrected geometry's honest residual**, with the partial-F1 caveat on
  record (the low side is mechanism-addressable; the tail side is not,
  by any of the three fingerprinted knobs).

**Data observation (reported, not adjudicated):** the missing tail mass
(0.055–0.076 over n ≥ 14) numerically shadows the **excluded trapped
class** (pooled trap_bound 0.0768), which §3.5b measured as
slow-and-high-n. If part of the physically-retained class is in fact
detected experimentally (or the retained boundary sits differently at
realistic droplets), the tail deficit is partly **policy, not
mechanism** — this puts a quantitative stake on the open retained-policy
adjudication (user), and is the natural first probe of the honest
residual's decomposition.

## §3.5g — the low-n KE retro-scan (2026-07-29, zero MD, 52 committed rows)

Scorer `tier2atlas_ke_lown_scan.py`, artifact `atlas_ke_lown_scan.csv`.
Pre-registration committed before the scorer ran (plan §3.5g, incl. the
**KE₁ anchor amendment to the experimental peak 1.00 eV** — user
adjudication; the median 1.128 stays reported). All three oracles passed
before any new number was read: **O1** scorer-drift (standing pooled
battery), **O2** committed-KE (all five h405 members + pooled reproduce
the committed `atlas_g4step2_battery.csv` KE columns to 1e-9), **O3**
incumbent-KE (pooled 1.034 / 0.754 string-exact at 3 decimals).

### R1 (τ) — the first τ read on the KE axis

Along the gated ridge chord h405→h375→h345 (E₀ co-moving to hold the
gate): **∂KE₁/∂τ = +0.083 eV/ps**, ∂KE₂/∂τ +0.054. In-gate headroom
Δτ ≈ +1.14 ps (the n₁ floor binds at τ ≈ 5.5) → **ΔKE₁ ≤ +0.094 eV**.
Real but ≈ 3× too small on its own. Within each fixed-τ arm KE₁ falls
monotonically with E₀ (τ 4.4: 0.668→0.612 over E₀ 0.395→0.415).

### R2 (v_c) — the strongest in-surface lever, and its price

| v_c | cell | KE₁ | KE₂ | midHot | gate |
|---|---|---|---|---|---|
| 5.00 | c50 | 0.956 | 0.881 | **1.923** | 0 |
| 5.25 | v525 | 0.774 | 0.696 | **1.394** | 1 |
| 5.50 | h405 | 0.637 | 0.551 | 0.957 | 1 |
| 6.00 | b031 | 0.523 | 0.364 | 0.520 | 0 |
| 6.50 | c65 | 0.279 | 0.222 | 0.415 | 0 |
| 7.25 | f725 | — (n₁ = 0) | — | 1.398 | 0 |

∂KE₁/∂v_c = −0.43 eV per Å/ps (v525↔f3 pair) at **∂midHot/∂v_c =
−1.48**: lowering the cap buys low-n KE and pays mid-n heat across the
*whole* measured range — c50 nearly lands both KE₁ and KE₂ (0.956/0.881
vs 1.00/0.706) at a disqualifying midHot 1.92. The KE₁↔midHot trade is
structural, not a linearization artifact. In-gate bound **+0.054 eV**
(midHot ceiling binds at v_c ≈ 5.37). Note v525 *gates* and lands KE₂
0.696 ≈ ref 0.706.

### R3 (E₀)

Arm-resolved, CRN-paired: ∂KE₁/∂E₀ −2.79/−3.59/−3.86 eV/eV (τ
4.4/4.8/5.2), ∂KE₂/∂E₀ −2.59…−3.50; gate slopes ∂n₁/∂E₀ +1.4…+2.4,
∂n̄/∂E₀ −22.5…−31.0. In-gate headroom ΔE₀ ≈ −0.014 (n₁ floor) →
**ΔKE₁ ≤ +0.038 eV**.

### R4 (geometry at the fixed standing chord) — the n = 1 channel is bifurcated

The 11-cell grid + incumbent battery + f725, all at (7.25, 3.2, 0.27):

- **Every deep-born cell produces zero n = 1 fragments.** Boltzmann and
  center cells (depth 24.9–49.4 Å): n₁_solv = 0.000, KE₁ undefined, at
  every R (6/6). f725 (standing chord at the corrected *ensemble*):
  likewise n₁ = 0, with n̄ 14.3 and trap 0.437 — at deep birth the
  standing chord **under**-strips and traps; it does not produce slow
  n = 1, it produces none.
- **Every standing-chord n = 1 is shallow-born and reference-fast.** The
  shallow cells (depth 9.0/10.9/14.8/19.5 Å at R 26.6/34.0/49.4/68.3):
  KE₁ **1.042 / 1.062 / 1.064 / 1.091 eV** (n = 184/127/79/48), KE₂
  0.782–0.827, with shallow n₁_solv falling with R
  (0.235/0.184/0.151/0.128). The incumbent pooled (sampled ⟨R⟩ 26.6,
  depth ≈ 9): KE₁ 1.034 ± 0.010 (members 1.020–1.045), KE₂ 0.754.
- Reading: the model holds **two disjoint n = 1 channels** — a
  shallow-birth channel whose KE sits exactly at the experimental peak
  (≈ 1.0–1.09 eV, mildly rising with R), and the corrected-basin
  deep-cascade channel at 0.64. The 2026-07-29 discussion's provenance
  question is settled *measured*: at the standing chord the user's
  "low shells from further out" is exactly right; at the corrected
  basin n = 1 is deep-born because the ensemble has no shallow births
  to offer. The experimental peak coinciding with the shallow channel
  makes the **birth-law/mixture lever** (shallow-birth or small-droplet
  weight absent from the pure center-weighted Boltzmann law) a second
  live candidate beside p_tail — and the only one of the two that also
  feeds the Block-D n = 1 *weight* deficit (−0.101) and the wide
  experimental n = 1 upper tail (a fast second channel is literally
  what the model's shallow branch is).

### R5 — the reachability verdict

From pooled h405 (KE₁ 0.641, required ΔKE₁ ≥ +0.309 eV to reach
0.95 = 1.00 − 1σ): E₀ ≤ +0.038, τ ≤ +0.094, v_c ≤ +0.054; **additive
single-knob bound +0.186 eV — the (v_c, τ, E₀) surface CANNOT reach the
KE₁ anchor. The in-surface freedom is exhausted on the low-n KE axis.**
The far-field cells confirm the bound is conservative in form but right
in structure: every distant cell that lifts KE₁ (c50, v525) exits
through the midHot ceiling first. Per the pre-registration, the p_tail
axis is now motivated **with this record as its evidence** (the same
structure by which the E₀ arm established the W₁ floor); the R4
mixture lever stands beside it as the second candidate, with the two
being distinguishable by design — p_tail acts on the *deep-born*
channel's toll (KE-first, weight-second), the birth-law lever *adds*
the fast shallow channel (weight-first, KE-free).

**Atlas stance:** nothing adopted; `finc1v725` stands; the G4
adjudications (successor, retained policy, ledger re-issue) remain
open; the next-axis choice (p_tail study vs birth-law/mixture probe vs
both) is a user adjudication.

## §3.5h — the p_tail ring (2026-07-29, 3 × N = 1000, seed 20260729, CRN-paired)

Generator `gen_tier2atlas_ptail_ring.py` (dry-run oracle-verified before
launch; one exterior kill mid-first-launch, same-seed relaunch —
deterministic, no physics impact), scorer `tier2atlas_ptail_table.py`,
artifact `atlas_ptail_ring.csv`. O1 (scorer drift) and O2 (committed
finals h405 row to 4 decimals) passed before any new number was read.

| cell | p_tail | n̄ | n₁ | supp | trap | W₁ | midHot | deepKE | KE₁ | KE₂ | χ²_med | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h405 | −1.0 | 3.955 | 0.208 | 0.183 | 0.087 | 0.709 | 0.957 | 0.603 | 0.637 | 0.551 | 350 | 1 |
| pt15 | −1.5 | 3.080 | 0.259 | 0.261 | 0.011 | 0.862 | 1.917 | 0.881 | 0.917 | 0.849 | 2093 | 0 |
| pt20 | −2.0 | 2.386 | 0.295 | 0.335 | 0 | 1.331 | 3.021 | 4.317 | 1.151 | 1.122 | 8676 | 0 |
| pt30 | −3.0 | 1.841 | 0.346 | 0.410 | 0 | 1.769 | 4.414 | 12.47 | 1.428 | 1.447 | 4.6e4 | 0 |

**Verdicts (pre-registered, plan §3.5h):**

- **PT-P1 partially refuted:** KE₁ monotone and the placement bands hold
  at pt15 (0.917 ∈ [0.79, 0.99]) and pt20 (1.151 ∈ [0.96, 1.16]); pt30
  1.428 over-band — the anchored 1-D model under-predicts at far
  softening (the cascade/mass history it freezes changes).
- **PT-P2 confirmed:** KE₂ rises monotonically — but past the reference
  (0.849 at pt15 vs ref 0.706): **KE₂ overshoots before KE₁ reaches its
  target.** The measured n = 1 : n = 2 relief differential of uniform
  tail softening is ≈ 1.3 : 1; the deficit profile needs ≳ 2.3 : 1.
- **PT-P3 — THE KILL CRITERION FIRED.** Both cells with KE₁ ≥ 0.95
  (pt20, pt30) carry midHot 3.0 / 4.4 ≫ 1.15. Per the frozen decision
  rule: **no further single-knob tail cells**; the axis moves to a
  joint re-tune or the honest-residual branch (user).
- **PT-P4 refuted with a sign inversion worth the price of the ring:**
  n̄ *falls* (3.96 → 1.84), n₁/supp *rise* (0.21 → 0.35 / 0.18 → 0.41),
  trap → 0. The pre-registered picture (less drag work → less heating →
  shallower cascades) is wrong at this knob: the dominant channel is
  residence/relaxation-timing — ions that stay fast through the 5.5–9.5
  band and into the Landau-gated E2 stage keep stripping instead of
  thermalizing. The tail does not just set the exit toll; it sets *how
  long every ion stays in the dissipative regime*.
- **PT-P5 refuted, maximally:** deep bins are the MOST affected axis
  (deepKE 0.88 / 4.3 / 12.5 vs the h405 battery's 0.50 ± 0.02) — the
  exit-speed orthogonality picture is invalid; the §3.5g slope
  arithmetic cannot be reused under a live p_tail.

**Reading.** The p_tail band as implemented (all v > 5.5) overlaps the
late-deceleration band of every ion class; softening it converts the
model into a globally under-damped system long before KE₁ lands. The
axis is not dead on selectivity grounds *per se* — the measurement
localizes what a KE₁-selective drag change must look like: relief
confined to **v ≳ 9.5 Å/ps** (above the n = 2 exit speed 9.1, below the
n = 1 exit 9.9; the band the n = 1 exiters hover in at late times and
every other class only crosses once, early and briefly). A band-limited
tail — a second cap v_c2 ≈ 9.5 with softening only above it, saturated
(p = −1) between v_c and v_c2 — is the shape this ring points at. That
is a new form-surface member (new enum value, new physics branch) and
therefore a fresh design + adjudication, not a parameter cell.

**Atlas stance:** nothing adopted; `finc1v725` stands; the §3.5h
registered options after the kill — joint (p_tail × v_c/τ/E₀) re-tune,
honest-residual branch, or the (unregistered) band-limited tail form —
are the user's call. The G4 adjudications stay open.

## §3.5i — the s(n) drag-state-coupling probe (2026-07-29, 3 × N = 1000, seed 20260729, CRN-paired)

**Instrument.** `gen_tier2atlas_sn_probe.py` (build S1–S4 delivered behind
`[PROCEED TO IMPLEMENTATION]`; design + registration in
`TIER2_DRAG_STATE_COUPLING_DESIGN.md` §8/§10): sa22/sa30/sa44 = ρ_shell
at bulk 0.0218 / prior 0.030 / 2×bulk 0.0436 Å⁻³, R_core 3.2 Å, at the
committed h405 pins, seed 20260729 (CRN-paired to the finals row). All
oracles pass: cfg-diff exactly the three coupling fields (pinned-copy
post-reference precedent), unit oracle s(19) = 1 exact + registered
s-table to 3 dp, O1 scorer-drift, O2 committed-row to 4 decimals.
Operational note: the first launch was externally killed mid-E2 (not
the user); E2 + detection were resumed from the completed `ion.npz`
per cell (stage-stream RNG is seed-derived and ion-stage-independent,
so the continuation is bit-equivalent — the G1 recovery precedent);
stored cfg.json asserted equal to the rebuilt config before resuming.

**Scored table** (`atlas_sn_probe.csv`; baseline = CRN h405 row first):

| cell | ρ_shell | KE₁ | KE₁_sd | KE₂ | KE₃ | n̄ | n₁ | supp | trap | midHot | deepKE | W₁ | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h405 | – | 0.637 | 0.038 | 0.551 | 0.471 | 3.96 | 0.208 | 0.183 | 0.087 | 0.957 | 0.603 | 0.709 | 1 |
| sa22 | 0.0218 | 0.596 | 0.038 | 0.512 | 0.437 | 4.25 | 0.191 | 0.168 | 0.092 | 0.897 | 0.522 | 0.771 | 1 |
| sa30 | 0.030 | 0.599 | 0.038 | 0.514 | 0.438 | 4.23 | 0.191 | 0.168 | 0.092 | 0.899 | 0.497 | 0.765 | 1 |
| sa44 | 0.0436 | 0.601 | 0.038 | 0.516 | 0.441 | 4.21 | 0.196 | 0.169 | 0.092 | 0.905 | 0.507 | 0.755 | 1 |

**Registered verdicts (design §8):**

- **SC-P1 REFUTED** — KE₁ moves −0.04 (predicted +0.2..+0.5), the
  ordering is *inverted* (sa22 < sa30 < sa44, i.e. strongest coupling =
  lowest KE₁), and no cell is above baseline at all.
- **SC-P2 REFUTED (the signature)** — per-bin n = 1 KE SD is 0.038 in
  every cell vs the ≥ 0.08 floor: **the needle does not break.** Per the
  registration this kills the design's central claim regardless of means.
- **SC-P3 REFUTED** — no ΔKE grading (all Δ negative, wrong order);
  deepKE −0.08..−0.11 (sa30 outside the ±0.10 band).
- **SC-P4 not evaluable** — no cell reaches KE₁ ≥ 0.95.
- **SC-P5 (exploratory):** n̄ +0.26..+0.30, n₁ −0.012..−0.017, supp
  −0.014, trap +0.005, W₁ +0.046..+0.062 — every cell strictly worse on
  the landing surface. Success-shape cells: NONE.

**§3.5i.1 The mechanism finding — the coupling is GATE-CLIPPED (measured,
zero MD, sa30 ion checkpoint).** Min-n-while-inside: **no ion reaches
n ≤ 8 inside the droplet** (fraction 0.0000; n ≤ 14 inside 7.4 %); mean
shell at first surface crossing **19.0 He, median 20** (ions exit fully
dressed); 93 % exited by sim-end; ion-stage-end n̄ 9.69 with all n ≤ 2
population produced *outside*. Stripping is slow relative to transit,
so the design §2 premise (early strippers transit as small objects) is
factually absent: being at low n and being inside the drag gate are
**mutually exclusive states** of the delivered mechanism. In-window the
shell stays in [14, 21] → s ∈ [0.85, 1.06], nearly ρ-independent by
construction (ρ differentiates s only where R_core³ competes with n/ρ —
the low-n regime the gate never sees), which explains the near-identical
cells across the factor-2 Bounded range. The small negative KE₁ is the
early-transit s(21) ≈ 1.06 over-drag at the highest speeds plus the
residence back-reaction (n̄ up — the §14.5 channel sign). The needle
survives because per-ion in-window s variance is minute: the
strip-timing decorrelation the design posited does not exist inside the
gate. **Consequence: no state-factor × gated-γ coupling of this family
can move KE₁ — the exit toll is paid dressed.** A successor low-n KE
lever must act outside the gate (a different force surface) or move the
strip timing itself (OQ-F cooling contact / upstream E_int dynamics).

**Disposition.** The axis STOPS on the SC-P2 signature failure (per
registration; no parameter chase). Options revert to the registered
pair — the honest-residual branch extended to KE₁, or the OQ-F
cooling-contact discussion — **user adjudication either way**. The
delivered code stays (default `off` = bit-identical; full suite 2848
passed); nothing adopts; `finc1v725` stands; the h405 successor
candidacy and all G4 adjudications are untouched. D0 §18 carries the
compact influence record.

## §3.5i.2 — Post-probe kinematic decomposition: the fast n = 1 ions already exist and are parked in a spurious fully-dressed bin (2026-07-29, zero MD, h405 baseline)

Discussion-triggering reads on the committed `g4fh405` checkpoints after
the s(n) refutation (user adjudication: cooling/OQ-F rejected as the
route; "find another explanation how n = 1 becomes 1 eV"). All numbers
measured; nothing adopted.

**(1) The 0.36 eV KE1 deficit is dominantly a kinematic identity, not
missing energy.** The n = 1 enders cross the surface at v_exit = 10.16
A/ps = **1.086 eV dressed KE** -- the snowball already exits above 1 eV.
Detection reads 0.637 eV because the physical co_moving cascade
(production convention, verified in the committed cfg) conserves
velocity, so the fragment KE is (1/2)*m(IHe1)*v^2 -- the mass ratio
131/203 plus a ~5 % post-exit slowdown (v_det/v_exit = 0.953). Reaching
1.0 eV in-bin needs v = 12.14 A/ps (+25 %) -- beyond every in-gate
lever (the measured closure of the whole gamma-surface program).

**(2) The measured v_exit -> terminal-n map is non-monotone and exposes
a spurious population:**

| v_exit [A/ps] | terminal n-bar (detected) | frac n <= 2 | count |
|---|---|---|---|
| 8-9 | 3.60 | 0.000 | 284 |
| 9-10 | 1.85 | 0.888 | 365 |
| 10-11 | 9.16 | 0.592 | 365 |
| 11-13 | **21.00 exactly** | 0.000 | 185 |

The FASTEST ions keep the FULL shell: E0 = 0.405 sits below the
full-shell gate, terminal n is pickup-driven (E_int accumulates with
residence), so the fastest transits never unfreeze the cascade. The
needle is the narrow [9, 10) selection slice. Consequence:
**352/2000 = 17.6 % of detected fragments sit at n = 20-21** -- a
delta-spike the smoothly-decaying experimental abundance cannot
contain (an independent falsification hook already inside the W1/chi2
residual). Their as-if-n = 1 KE: p10-p90 **0.69-0.87 eV** (full
n >= 20 group, v-bar 10.15); the fast [11, 13) subset reaches
**0.82-1.15 eV** -- exactly the experimental n = 1 band, width for free.

**(3) Candidate explanations left standing (next-session discussion
agenda; NO design, NO code -- physics definition first):**

- **(A) Exit stripping / surface scrubbing (front-runner):** the
  21-He shell should not survive a 10-13 A/ps (~1 km/s) interface
  crossing intact; an impulsive shell-loss channel at the surface
  crossing rebins the fast dressed exiters into low n AT THEIR
  EXISTING KE (co-moving kinematics: the later cascade cannot change
  per-fragment KE). One physical statement fixes both anomalies: the
  n = 1 position/width AND the spurious n = 20-21 spike. Needs a
  defensible closure (ram energy per shell He vs binding; shell-depth
  argument) + an enum on the exit event. Prediction: n = 1 KE inherits
  the current fast population (0.69-1.15 eV), roughly budget-blind.
- **(B) Per-ion source KER spread:** the experimental Coulomb-explosion
  KER is a channel distribution, not the single 2.70 eV; E0 =
  f_int x budget is ALREADY wired proportionally, so a high-KER ion
  exits faster AND carries E0 above the gate -> strips -> fast n = 1.
  The p_tail = -1 constant-force tail makes the KER -> exit-KE
  transfer slope ~ 1 (real authority). Checkables first: the
  `vmi_iplus_gas.csv` channel structure, the 2.70 eV provenance;
  then a 2-3-cell budget probe (existing knob, zero new physics code).
  Prediction: n = 1 KE tracks the source channel energy.
- **(C) A + B combine; discriminator = budget dependence of the n = 1
  KE** (A: budget-blind at 0.69-1.15; B: slope ~ 1).

**(4) Dead ends closed by measurement (do not reopen):** the entire
gamma(v)/s(n) drag surface (Tier-0 lock + §3.5g/h/i); cooling contact
(user adjudication + the gate-clipping timing bar); the shed convention
(co_moving is the physical path; cold-shed is the only KE-raising
branch and is the unphysical diagnostic bound); per-shed recoil
(meV-scale, NB-RQ23-1).

**Open verification items before any design:** (i) confirm the
experimental abundance reference has ~zero weight at n >= 15-21 (the
spike falsification); (ii) read the gas-phase VMI I+ reference for KER
channel structure; (iii) audit the 2.70 eV budget provenance
(legacy MATLAB path). Reads only -- all zero MD.

## §3.5i.3 — The three pre-design verification reads EXECUTED (2026-07-29, zero MD, read-only; cold-shed reintroduction DECLINED en route)

Adjudication recorded first: the user asked whether reintroducing
**cold shed** would solve KE₁ before any §3.5j work; declined on three
measured/physical grounds — (1) evaporation is isotropic in the
snowball frame (co_moving + meV recoil, NB-RQ23-1); cold shed requires
perfectly collimated retro-exhaust with ~30 meV *directed* energy per
He, (2) the boost factor is fixed at (m_exit/m_frag)², giving KE₁
0.637 → ≈ 1.53 eV and KE₂ 0.551 → ≈ 1.25 eV — untunable overshoots
past both anchors with midHot rising alongside, (3) the n = 20–21
spike population never sheds at all, so shed kinematics cannot rebin
it. The §3.5i.2 item-4 dead-end list stands unchanged.

**(i) Abundance weight at n ≥ 15 — the spike falsification is
CONFIRMED (nuanced).** `integrated_i_he_abundance.csv` (n = 0–20 mass
windows): smooth monotone decay; n ≥ 15 carries **3.21 %** of all
detected (5.69 % of solvated) — "~zero above 15" was too strong. At
the spike bins the verdict is unambiguous: **n = 19–20 together carry
0.59 %** (0.33 % → 0.26 %, no upturn; the table ends at n = 20,
extrapolated n = 21 ≲ 0.25 %) vs the model's 17.6 % at n = 20–21 —
a **~30× excess**. The dressed-freeze delta is independently
falsified, W₁-independent.

**(ii) Gas-phase VMI I⁺ KER structure — a multi-channel distribution,
VERIFIED.** `vmi_summary/vmi_iplus_gas.csv` (Abel-inverted 3-D speed
distribution) resolves ≥ 5 channels (per-fragment KE at I mass; pair
KER = 2× for equal-mass partners): 0.21 / 0.27 (rel. h ≈ 0.5), 0.49
(0.84), **0.70 (main)**, 2.26 (0.80), 4.11 eV (0.36); main-peak FWHM
alone spans KE 0.31–1.10 eV. The model's single 2.70 eV/fragment
(v = 20.3 Å/ps) sits in the **valley between the two fast channels**,
~20 % above the measured CE channel (2.26 eV ⇒ pair 4.52 eV ⇒
point-Coulomb distance ≈ 3.2 Å); the dominant experimental channel
(0.70 eV) is ~4× below the model budget. Cross-check for free: the
**droplet** reference's main peak sits at 12.30 Å/ps = **0.996 eV at
bare-I mass** — within 1.3 % of the v = 12.14 Å/ps the §3.5i.2
decomposition says the model's n = 1 fragment needs in-bin; the 1.00 eV
KE₁ anchor *is* this peak. Caveats for candidate (B): gas-phase channel
*weights* need not transfer into the droplet, and the droplet curve is
post-drag; a dominant 0.70 eV source could never reach the 1.0 eV
target after any positive toll — (B) requires the droplet signal to be
fed by the fast CE channels.

**(iii) The 2.70 eV budget provenance — a geometric idealization, not
a measurement.** The budget is the point-charge Coulomb pair energy at
the *neutral* I₂ ground-state separation: 14.4 eV·Å / 2.666 Å =
5.401 eV pair, symmetric split ⇒ 2.70 eV/fragment (documented in
`RESEARCH_QUESTIONS.md` RQ7; stamped `tier2_h2b_forward_model.py`;
realized dynamically — production births both ions at
`R0_GS_angstrom = 2.666`, `presets.py`, and the propagated I⁺–I⁺
Coulomb repulsion delivers the KE; `coulomb_available_eV` is a
provenance stamp, no hard refuse). The legacy MATLAB `D_e = 2.7 eV`
(I₂⁺ X-state Morse, doi 10.1063/1.475194 Table V,
`I2+_potential_energy_curve_check/`) is a **numerical coincidence**,
not the source. The model's source term is therefore a delta at an
idealized upper-channel value — exactly the single-valuedness candidate
(B) targets. The RQ register's bare-peak note (3.706 eV > the 2.71 eV
ballistic ceiling) reads as the same gap from the other side: the
measured 4.11 eV channel is the missing supply.

**Consequences for the §3.5j discussion:** (A)'s falsification hook is
armed (a real 30× spike excess, not just KE residual); (B)'s premise is
verified but sharpened (channel structure real; in-droplet weights
unknown and only the fast channels can feed KE₁); the KE₁ anchor is
identified with the droplet VMI main peak. Nothing adopted;
`finc1v725` stands; h405 candidacy + G4 adjudications open.

## §3.5j — the (A)/(B) discussion round (2026-07-29): CE-variant DECLINED (bulk-refill); OQ-J..N adjudicated; the two placement reads EXECUTED — **(A) measured NECESSARY-NOT-SUFFICIENT (bracket ceiling KE₁ ≈ 0.71); the (C) combination is REQUIRED**

**CE-instantaneous stripping declined (user query, recorded per
request).** The variant "the shell is stripped at the Coulomb explosion
itself" (≙ born-bare, replacing the T5 birth-dressing convention) fails
for KE₁ on the **bulk-refill argument**: *fast collisions knock shell He
off wherever ram > binding; the bath refills the shell wherever there is
bath* — so net stripping exists only where refill ends, the outbound
surface crossing. At t = 0 the ion sits at maximum bath depth
(impulse-stripping is maximally self-healing; the model's own pickup
channel is the refill, cf. the T5 comment "under-dressed ions may
re-fill via the live pickup channel"); a bare-transiting ion would in
any case detect at n = 0 (RQ3 territory), not n = 1; and born-bare
rewires m(t) globally (the whole landing re-arbitrated). The same
statement makes dressed-in-bulk, the TDDFT dressed-object calibration,
the §3.5i gate-clipping result and exit stripping mutually consistent.
The legitimate residue — is `full`/`density_tied` birth dressing the
right convention — already lives in the T5/T6 enum surface + D0 §17.

**Adjudications (user 2026-07-29: "I follow your recommendations"):**
OQ-J closure class = stochastic per-He knockout (amended by the
counterfactual below: a **depth-graded survival element is required**,
not a pure binding-threshold or pure Poisson-count rule); OQ-K
bookkeeping = re-labeled surface-crossing drag work with an explicit
strip term (no new energy source/sink; ~0.2 eV scale ≈ the measured
post-exit slowdown); OQ-L bare outcomes **allowed** (no artificial
n = 1 floor; feeds RQ3); OQ-M zero-MD placement reads first (executed
below), the 2–3-cell budget probe stays the registered discriminator;
OQ-N "(B) parked" — **SUPERSEDED same-session by the reads** (see
verdict).

**Read 1 — the experimental n = 1 KED upper tail (zero MD,
`ihe_ked/IHe_KED_curves_n1.csv`, 3-D P(E)):** weight above 1.15 eV =
**48.7 %** (above 1.0: 59.4 %; above 1.5: 32.7 %); reference row: mean
1.302, mode 0.891, median 1.128, σ 0.697 eV. The (A)-unreachable
region (above the ≈ 1.15 eV cap set by the existing exit-velocity
distribution) is **half the distribution**, not a tail. Caveat: n = 1
is the one fragment whose ⟨E⟩ is background-choice sensitive
(`bgOffShift`, COLUMNS.md) — the peak region is firm, the far tail
less so.

**Read 2 — the (A) exit-stripping counterfactual rescore (zero MD,
committed `g4fh405` checkpoints, detection-only; scratchpad
instrument).** Assumptions: CF-1 strip at first outbound crossing
(terminal n_cf = min(n_detected, n_strip); post-exit cascade not
re-run), CF-2 velocity history unchanged (KE_cf = ½·m(n_cf)·v_det²),
CF-3 knockout = equal-mass max transfer vs the production `rq4graded`
rungs. **Convention gate passed** before any counterfactual number:
local conventions reproduce the committed row exactly (n̄ 3.955, supp
0.183, n₁ 0.208, KE₁ 0.637 ± 0.038, W₁ 0.709).

- **Structural identification (a finding on its own):** the §3.5i.2
  "n = 20–21 spike" IS the **`suppressed` fate class** — all 334
  suppressed fragments carry mechanical shell n = 21 (relaxation-end
  n ≥ 19: 381 = 334 supp + 47 trapped), cross the surface at v
  10.66–11.69 Å/ps (p10–p90), and are scored as **bare-equivalents**
  (supp = 334/1826 = 0.183; n̄ = Σn/1826 = 3.955 with suppressed as 0).
  The model already relabels its fast dressed exiters at scoring level;
  (A) would replace that scaffolding with physics. D0 §17 relevance:
  the suppression rule is the scaffolding this axis would retire.
- **Variants (detected ensemble 1826; W₁_loc = this read's solvated
  n = 1–20 convention):**

| variant | n̄ | supp | bare | n₁_solv | KE₁ | KE₁ p10–p90 | W₁_loc |
|---|---|---|---|---|---|---|---|
| committed (gate) | 3.955 | 0.183 | 0.000 | 0.208 | 0.637 ± 0.038 | 0.59–0.69 | 0.709 |
| sharp knockout | 2.901 | 0.000 | 0.329 | 0.445 | 0.484 ± 0.070 | 0.39–0.58 | 0.787 |
| stoch η = 0.3 | 6.198 | 0.042 | 0.000 | 0.177 | 0.637 | 0.59–0.69 | 1.613 |
| stoch η = 1.0 | 4.980 | 0.000 | 0.019 | 0.195 | 0.618 ± 0.075 | 0.53–0.69 | 0.796 |
| **bracket: supp → n = 1** | 4.138 | 0.000 | 0.000 | 0.353 | **0.708 ± 0.087** | 0.60–0.83 | 0.751 |

- **Mechanism readings:** (a) the fast exiters sit ABOVE the 9.9 Å/ps
  full-strip threshold, so a pure binding-threshold rule sends them to
  **bare, not n = 1** — and also strips the *current* n = 1 enders
  (v_exit 10.16) to bare, dropping KE₁ to 0.484; (b) Poisson-limited
  stripping (N_swept ≈ 14 collisions vs a 21-He shell) cannot take
  21 → 1 and parks the spike at n ≈ 7–17, exploding W₁ (1.5–1.6);
  (c) hence the OQ-J amendment: the closure needs **depth-graded
  survival** (outer shell eager, inner 1–3 shadowed/protected) — a
  real, physically-motivated but shape-carrying design element;
  (d) the **bracket** (maximal conversion, every suppressed → n = 1 at
  existing KE) measures the (A) ceiling: **KE₁ 0.708 mixed / 0.774
  rebinned-only (p10–p90 0.70–0.87), 0.3 % above 1.0 eV** — n₁ 0.353
  lands at/above the reference solvated share 0.310, the needle breaks
  (SD 0.038 → 0.087), supp → 0.

**Verdict (measured, superseding the discussion's front-runner
framing):** (A) fixes the spike, the n₁ *weight*, and the needle width,
and is the only candidate that retires the suppression scaffolding —
but its KE₁ ceiling is ≈ 0.71–0.77, ~0.3 eV below the 1.00 eV peak
anchor, and the 48.7 % KED weight above 1.15 eV is (A)-unreachable by
construction. **(B) is therefore necessary for the KE₁ position and
upper half, not merely the tail** — the natural division of labor is
the (C) combination: (B) supplies the source spread (fast CE channels →
faster exiters), (A) converts fast + dressed into small-at-existing-KE
at the boundary; the n = 1 KED becomes the drag-tolled image of the
fast KER channels (its σ 0.70 eV from channel structure, its
multi-eV tail from the 4.11 eV channel — unreachable from any
single 2.70 eV source). **Next (user adjudication): the (C) physics
definition (both mechanisms, one design doc) + registration of the
budget-slope MD probe (existing `coulomb_available_eV`/R₀ knob, 2–3
cells) as (B)'s authority measurement and the registered
discriminator.** Nothing adopted; `finc1v725` stands; h405 candidacy +
G4 adjudications open.

## §3.5k — the budget-slope probe EXECUTED (2026-07-29, 3 × N = 1000, seed 20260729, CRN-paired): **BP-KILL NOT FIRED — S_k = 0.389 eV/eV, the source-side lever is ALIVE**; BP-P3 sign-inverted (suppression is E_int-driven); needle persists; the (C) design proceeds per registration

Instruments `gen_tier2atlas_budget_probe.py` /
`tier2atlas_budget_table.py` (registration plan §3.5k; dry-run + scorer
oracles verified green before launch). O1 (scorer drift) and O2
(committed finals h405 row to 4 decimals) passed before any new number
was read. Artifact `atlas_budget_probe.csv`. Compact influence record:
**D0 §19** (updated first, per the standing rule). Mid-probe domain
input (user, from the paper): the gas spectrum is reproduced at
**0.8·E_C with Q = 2/Q = 3 from 2.666 Å** — channels re-anchored to
2.16 / 4.32 eV per I⁺, the fast peak assigned I⁺–I²⁺; full record RQ7/
RQ8 NBs 2026-07-29 (the probe's 2.26/4.11 brackets remain valid as
slope points).

| cell | budget | E₀ | n̄ | n₁ | supp | trap | W₁ | midHot | deepKE | KE₁ | KE₁_sd | KE₂ | gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h405 | 2.7006 | 0.405 | 3.955 | 0.208 | 0.183 | 0.087 | 0.709 | 0.957 | 0.603 | 0.637 | 0.038 | 0.551 | 1 |
| bud226k | 2.26 | 0.405 | 5.325 | 0.126 | 0.073 | 0.179 | 1.365 | 0.702 | 0.554 | 0.486 | 0.034 | 0.404 | 0 |
| bud411k | 4.11 | 0.405 | 1.632 | 0.329 | 0.546 | 0.001 | 1.294 | 2.301 | 1.093 | 1.200 | 0.047 | 1.081 | 0 |
| bud411 | 4.11 | 0.6165 | 0.344 | 0.243 | 0.909 | 0.001 | 1.222 | 0.861 | 0.409 | 0.623 | 0.036 | 0.525 | 0 |

**Registered verdicts:**

- **BP-P1 near-miss down / hold up:** bud226k 0.486 sits 0.006 eV above
  its [0.29, 0.48] band; bud411k 1.200 ∈ [1.13, 1.70]; ordering rises
  monotonically. The registered diagnose-obligation is discharged by
  BP-P2 itself: the 1-D constant-toll placement (slope 0.587) is ~⅓
  too steep — the cascade back-reaction eats the difference.
- **BP-P2 CONFIRMED — the probe's deliverable: S_k = 0.389 eV/eV**
  (E₀-pinned line 2.26/2.7006/4.11 → KE₁ 0.486/0.637/1.200; registered
  band [0.35, 0.75]).
- **BP-KILL NOT FIRED** (0.389 ≥ 0.2): the source-side lever is alive;
  per the pre-registration the honest-residual branch is *not* forced
  and **the (C) design discussion proceeds with S_k as its measured
  authority**.
- **BP-P3 REFUTED with a sign inversion (the mechanism yield):**
  supp(bud411) 0.909 > supp(bud411k) 0.546 — raising E₀ at fixed
  kinematics *adds* +36 pp suppressed (∂supp/∂E₀ ≈ +1.7 /eV).
  Suppression is **E_int-driven (self-unbound at detection)**, not
  freeze-driven. Corollary measured in the same rows: under the full
  proportional (B) wiring the fast ions all exit through the
  suppressed-bare door and **KE₁ stays at baseline** (bud411 KE₁ 0.623,
  n = 1 count 44) — Q3 can feed the n = 1 upper half only through (A)
  exit stripping (rebin at existing KE) or partial E₀ decoupling. The
  (A)/(B) complementarity is now measured from both directions.
- **BP-P4 CONFIRMED:** KE₁ SD 0.034/0.047 at the k-cells — the needle
  survives any single-valued budget; n = 1 width requires the mixture.
- **Exploratory (uniform-budget — explicitly NOT the channel-weighted
  (C) forecast):** a uniform re-anchor to the calibrated 2.16–2.26
  breaks the standing basin (bud226k: n₁ 0.126, trap 0.179, W₁ 1.365);
  a uniform 4.11 point overheats globally (midHot 2.301, deepKE 1.093).
  Both are the expected price of moving *every* ion; the (C) mixture
  applies the fast channel to a small weight.

**Atlas stance:** nothing adopted; `finc1v725` stands; h405 candidacy +
G4 adjudications open. **Next (user): the (C) physics-definition design
doc** — paper-anchored channel set {single-ionization, Q2 2.16, Q3
4.32 at 0.8·E_C}, per-channel E₀ coupling as a design axis (measured
leverage), (A) exit stripping with the depth-graded survival closure,
weights anchored on bare peak / I²⁺ / power series — plus the named
zero-MD verification reads (bare-KED bimodality, I²⁺ presence, power-
matched weights).

## §3.5l reads (i)+(iii) EXECUTED (2026-07-29, zero MD, read-only) — bare KED is FAST-FED and (weakly) two-lump; gas Q3 fast-share rises 10→17→23 % with power; **NEW: a ~24 % calibration-frame split between the repo's canonical exports and the paper's 0.8·E_C channel values (user adjudication required)**

Scratchpad instrument (`read_i_iii_instrument.py`, session-local;
trapezoid moments of `signal_3d_PE` on the native E grid, negatives
clipped). **Convention gate:** the §3.5j Read-1 scratchpad weights are
not bit-reproducible (session-local code), so the gate was moved to the
committed artifact: the `IHe_KED_reference.csv` n = 1 row is reproduced
exactly (mean 1.302 / σ 0.697 at the diagnosed signal-mask edge
E ≤ 3.5 eV; mode 0.891 at smoothing window 7–9; median 1.151 vs 1.128
— the one residual, the row's median uses the pipeline's own binning).
This instrument's full-range weights above 1.0/1.15/1.5 eV =
60.3/51.1/33.8 % (masked: 59.5/50.1/32.4) bracket the §3.5j quotes
59.4/48.7/32.7 within ≤ 1.4 pp — the "half the distribution above the
(A) cap" conclusion is convention-robust.

**Read (i) — bare-KED structure (`IHe_KED_curves_n0.csv`, 3-D P(E),
valid E > 0.4 eV Abel-center cut).** Committed row: mean **3.706** /
mode **4.758** / median 3.776 / σ 1.356 eV; `bgOffShift = 0.000`
(the bare row does **not** hinge on the background choice, unlike
n = 1's 0.378 eV); N_eff 193 725. Instrument full-range agrees
(3.711/1.351). Band weights: E < 1.0 eV **1.9 %**, 1.0–1.6 4.1 %,
1.6–3.05 **27.8 %**, > 3.05 **66.2 %** (45.9 % in 3.05–5.0 +
20.3 % in 5.0–6.2); above 2.16 / 3.7 / 4.32 eV = 84.6 / 52.7 /
38.7 %.

- **Bimodality verdict: WEAKLY PRESENT.** Two broad lumps — lower
  ≈ 2.0–2.7 eV (maxima 2.3–2.5 at windows 9–15) and dominant upper
  ≈ 3.5–5.1 eV (mode 4.76) — separated by a shallow dip near
  3.1–3.3 eV. A two-component decomposition is admissible, not
  compelled by shape alone; the *weights* are the firm content.
- **The toll arithmetic is mean-level, not mode-level:** mean 3.706 =
  4.32 − 0.61 (≈ the measured post-exit slowdown scale), but the upper
  lump's mode 4.758 sits *above* the paper-frame Q3 4.32 (see the
  calibration-frame item below, which resolves this naturally).
- **The bare row is fast-fed — a NEW (C) kill-direction.** The current
  model's bare-equivalents (the `suppressed` class) cross the surface
  at 10.66–11.69 Å/ps → KE 0.75–0.90 eV, a region holding **≤ 2 %**
  of the experimental bare KED (band 0.4–1.0 eV: 1.9 %). The
  experimental bare peak is dominated by fast under-dressed
  Q3/Q2-transit ions, NOT by tolled single-budget ions. Any (C)
  closure that routes its slow suppressed/stripped ions into n = 0
  at ~0.8 eV overpopulates a nearly-empty region — slow-bare weight
  is a scoreable kill axis for the (C) probe. (Caveat: the n = 0
  3-D curve is blind below 0.4 eV — genuinely *slow* bare ions
  < 0.4 eV are invisible here; RQ3 territory.)

**Read (iii) — power-matched gas weights (radial exports
`iplus_gas_{160,300,600}mw`, 2-D angle-integrated projections, NOT
Abel-inverted; 300 mW from paper_v2 center [524.5, 380.8], 160/600 mW
paper_v4 center [509.4, 387.6] — README known issue; center fill
inflates the slow band, so rankings are the robust content).** Bands
per I⁺ (m 126.90 u): slow < 1.0 / Q2 1.0–3.05 / Q3 > 3.05 eV (3.05 =
geometric mean of 2.16 and 4.32):

| power | slow | Q2-band | Q3-band | Q3/(Q2+Q3) |
|---|---|---|---|---|
| 160 mW | 74.6 % | 22.8 % | 2.6 % | **10.2 %** |
| 300 mW | 70.0 % | 24.9 % | 5.0 % | **16.8 %** |
| 600 mW | 52.7 % | 36.3 % | 10.9 % | **23.1 %** |

The Q3 fast-share rises monotonically with power and at 600 mW
(**23.1 %**) matches the I²⁺ covariance anchor (~20 % at
2.94×10¹⁴ W/cm² under the provisional 600 mW mapping) — the two
independent instruments agree at the reference condition. At 300 mW
the band read (16.8 %) sits above the covariance ~7 %; expected
direction (2-D Q2-tail leakage into the Q3 band inflates the band
read; the covariance counts only true coincident I²⁺). The slow band
(single-ionization + center fill) *decreases* 75 → 53 % as multiple
ionization grows. **Weight-anchoring consequence (design axis 4):
the power series confirms the covariance ordering and the 600 mW
Q3 ≈ 20 % anchor; the bare-row KED anchors the *joint*
(source × retention) prediction, not the source weights directly —
it is retention-filtered (fast Q3 escapes under-dressed → bare,
slow Q2/single dress → I⁺Heₙ), which is exactly what a (C) forward
model predicts and read (i) measures.**

**Rim-position rider — POSTED AND WITHDRAWN SAME-SESSION (record
kept per the withdrawn-readings discipline).** The 600 mW 2-D gas rims
sit at ≈ 2.67/4.88 eV — ~24 % above the paper's 0.8·E_C channel values
— which was briefly posted as a possible calibration-frame split
("repo frame at full E_C, no re-anchor owed"). **Refuted by the
committed Abel instrument the same session** (user prompt: the 0.8
factor is paper physics, not a frame): the trusted 3-D gas export
`vmi_summary/vmi_iplus_gas.csv` (43632, Abel-inverted I(v), canonical
vf) peaks at **2.26 / 4.11 eV** — the paper's 0.8·E_C channels
2.16/4.32 confirmed to ~5 % *in the repo's own calibration* (inside
the declared 4 % `calibSyst_frac`, and identical to the RQ8-NB session
readings). The 2-D radial exports are the biased instrument for
*absolute* positions (center ambiguity per the paper_v2/v4 READMEs,
projection, toolbox radial-profile convention; an r-Jacobian removal
does not reconcile them) — they remain valid for the *power ordering*
(the read-(iii) table above), not for channel energies. Standing
consequences (RQ7 NB unchanged): the production budget (scale 1.0 →
2.70) stays **known-high ~25 % vs the experiment's own calibration**;
any re-anchor still owes the bud226k-measured basin price; the
paper-frame 2.16/4.32 ARE the repo-frame channel energies.

**Provenance question ANSWERED same-session (user supplied the source
paper for this session's read; consulted 2026-07-29, then removed —
not repo-kept, key numbers extracted below).** The 0.8 factor's source
is **Hatherly, Stankiewicz, Codling, Frasinski & Cross, J. Phys. B:
At. Mol. Opt. Phys. 27 (1994) 2993–3003
(doi:10.1088/0953-4075/27/14/032)** — the multielectron dissociative
ionization of I₂ in
intense laser fields, a **pure gas-phase experiment (no helium
anywhere)**. The user's "reduced in helium nanodroplets" recall is
thereby corrected: the reduction is **finite-pulse CE physics** —
higher charge states are created sequentially while the ions already
separate (equivalently: vertical transitions from laser-induced
stretched states at ~3.5–4.1 Å), so peak channel KERs are a
**channel-independent fraction of the point-charge E_C from
2.666 Å**: measured ≈ 0.75·E_C at 90 fs and ≈ 0.65·E_C at 200 fs
pulse length (their Table 1). Consequences for the record:
(1) the 0.8 is **universal source physics** — legacy applying
`E_coulomb_scale = 0.8` in both the gas *and* droplet presets is
correct implementation of this calibration, and the in-droplet source
channels are **Q2 2.16 / Q3 4.32 eV per I⁺** with no additional
droplet-screening factor claimed anywhere in the record;
(2) the exact fraction is pulse-length-tied (0.65–0.75 in Hatherly's
conditions; the group's own 0.8 fit on the *own-experiment* gas
spectrum is the operative calibration, and the Abel export confirms
it at 2.26/4.11); (3) the standing scale-1.0 production budget
(2.70) remains known-high vs the source physics — RQ7 unchanged;
(4) **two (C) design supports fall out of the paper:** the fraction
is measured *channel-independent* at fixed pulse → ONE shared
Bounded fraction f ≈ 0.8 across {single, Q2, Q3}, not per-channel
factors; and the channels carry **intrinsic KER width** (total-KER
FWHM 2.6 eV for I⁺+I⁺, 5.8 eV for I⁺+I²⁺ at 200 fs) — the (B)
source spread is paper-measured, per-channel KER must be *sampled
with width*, not fixed at the channel mean (feeds the n = 1 KED
σ 0.70 eV requirement directly).

**Atlas stance:** nothing adopted; `finc1v725` stands; h405 candidacy +
G4 adjudications open. The (C) design discussion now has all three
§3.5l pre-design reads on record ((ii) was answered 2026-07-29).

## (C) pre-steps P1 + P2 EXECUTED (2026-07-29, zero MD) — P1: the strip prior box PINNED (a in [1.5, 2.5], j0 in [1.5, 2], w_j in [0.5, 1]; a = 1 disfavored on W1-parking + slow-bare); P2: measured channel widths are NARROWER than the Hatherly priors (sigma_Q2 0.31, sigma_Q3 0.55 eV upper bounds; E_single ~ 0.53); OQ-K CLOSED by user confirmation

**Context.** Design `TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` §9; user
trigger "continue with P1+P2" (2026-07-29). Same session, the OQ-K
condition mapping was **user-confirmed: 300 mW ↔ 1.47×10¹⁴ / 600 mW ↔
2.94×10¹⁴ W/cm²** — the 600 mW `ihe_ked` reference pairs with the
~20 % Q3 covariance share; w_Q3 ≈ 0.20 is no longer provisional. The
m/q-127-row covariance ask (the clean w_single/w_Q2 closer) is
**postponed by the user until the (C) method proves successful**; the
split rests on the gas slow band + Q2 partner structure at
registration.

### P1 — strip-form prior calibration (committed instrument `tier2atlas_ce_strip_prior.py`, artifact `atlas_ce_strip_p1.csv`)

Detection-only rescore of the committed `g4fh405` checkpoints with the
actual design form P_knock(j) = min(1, (v_x/9.9)^a)·G(j),
G = 1/(1+exp(−(j−j₀)/w_j)); §3.5j CF-1/CF-2 conventions unchanged
(first outbound crossing; velocity history unchanged; toll NOT
applied); 48 Bernoulli replicas, seed 20260729; grid 5×5×4 = 100
cells (a ∈ {1…3}, j₀ ∈ {1…3}, w_j ∈ {0.5…2}).

**Oracles (all passed before any new number):** O1 scorer drift; O2
committed finals h405 row to 4 decimals; O3 kinematics (stored KE =
½·m(n)·v² bit-tight; suppressed carry the full n = 21 mechanical
mass 210.955 amu); **O4 crossing extraction reproduces the committed
§3.5j suppressed band p10–p90 10.66–11.69 Å/ps and the n = 1-ender
mean v_x 10.16 exactly**; **O5 the CF-3 sharp variant AND the
supp→n=1 bracket reproduce the committed §3.5j variants table**
(n̄/n₁/KE₁/bare/W₁ within 0.005). 1951/2000 fragments cross
(every scored fragment crosses; the 49 non-crossers are retained).

**Measured map (w_j = 1 slice; full 100-cell CSV committed):**

| (a, j₀) | P(0)/P(1)@10.2 | n₁ | KE₁ | KE₁ SD | bin0 | n̄ | W₁ | midHot |
|---|---|---|---|---|---|---|---|---|
| 1, 1 | 0.63 | 0.304 | 0.650 | 0.115 | 0.127 | 3.47 | 0.932 | 0.829 |
| 1, 2 | 0.23 | 0.261 | 0.654 | 0.096 | 0.033 | 3.82 | 1.041 | 0.895 |
| 2, 1.5 | 0.38 | 0.263 | 0.673 | 0.091 | 0.064 | 4.09 | 0.643 | 0.981 |
| 2, 2 | 0.23 | 0.242 | 0.667 | 0.083 | 0.030 | 4.19 | 0.733 | 1.002 |
| 2, 3 | 0.09 | 0.193 | 0.649 | 0.061 | 0.003 | 4.40 | 0.850 | 1.059 |
| 3, 2 | 0.23 | 0.234 | 0.671 | 0.080 | 0.029 | 4.23 | 0.722 | 1.008 |

(baseline gate row: n₁ 0.208, KE₁ 0.637, SD 0.038, W₁ 0.709.)

**Readings (the prior-box pin — values moved inside the OQ-D box,
family untouched):**

1. **The retention-twin anchor selects (j₀, w_j), independent of a**
   (at the twin speed 10.2 > v_strip 9.9, P₀ = 1 for every a):
   P(0)/P(1) ≈ ⅓ lands at **(j₀ 1.5, w_j 0.5–1)** [0.32–0.38] or
   **(j₀ 2, w_j 1)** [0.23]; j₀ = 3 kills the bare twin (0.09),
   j₀ = 1 over-produces it (0.63–0.86).
2. **a = 1 is disfavored** — the soft velocity gate strips slow
   exiters too: W₁ degrades above baseline (0.93–1.08, the CP-7
   mid-tail-parking direction) and bin0 doubles (slow-bare, the CP-5
   kill direction). a ∈ [1.5, 2.5] behaves ram-like: W₁ *improves*
   (0.58–0.73) and bin0 stays ≤ 0.07. a beyond 2 is nearly
   degenerate with 2 (P₀ saturates over the fast exiters).
3. **Pinned prior box for registration: a = 2 (range [1.5, 2.5]),
   j₀ ∈ [1.5, 2], w_j ∈ [0.5, 1].** Inside it: n₁ 0.24–0.31 (CP-4
   direction, reference 0.31), KE₁ 0.66–0.69 with SD 0.08–0.10
   (needle broken, CP-2 threshold met at the box edge), bin0
   0.03–0.07, midHot 0.96–1.01 (CP-6 flat), W₁ 0.58–0.73 (CP-7
   safe).
4. **KE₁ stays at the (A) ceiling** (0.65–0.69 ≈ the §3.5j bracket
   0.708; converted-suppressed n = 1 members mean 0.774 = the
   committed rebinned-only value — instrument-consistency check).
   Position remains (B)'s job, exactly per the design's division of
   labor.
5. Stated limitations: single-budget counterfactual (no mixture — the
   fast Q3 population that CP-5 needs for the bare row does not exist
   here, so bin0 here is pure slow-bare risk); CF-1 survivor counts at
   high n are ceiling reads (post-strip cascade not re-run); toll not
   applied (§3.5j limitation carried; MD probe measures it).

### P2 — Abel gas channel-width read (committed instrument `tier2atlas_ce_gas_widths.py`)

Jacobian-correct KED P(E) ∝ I(v)/v from the committed
`vmi_summary/vmi_iplus_gas.csv` (43632, pyabel rIbeta 3-D I(v));
**peak oracle passed** (I(v) maxima reproduce the committed
2.26/4.11 eV within one grid step); bounded three-Gaussian fit
(single/Q2/Q3), RMS residual 4.4 % of peak, window-stable
(σ moves ≤ 0.02 eV across four fit windows):

| component | µ [eV] | σ [eV] | area share | design prior |
|---|---|---|---|---|
| single | 0.534 ± 0.011 | 0.416 ± 0.013 | 0.655 | E_single band [0.3, 0.8] |
| Q2 | 2.258 ± 0.022 | **0.311 ± 0.023** | 0.221 | Hatherly σ 0.55 |
| Q3 | 3.913 ± 0.103 | **0.553 ± 0.113** | 0.123 | Hatherly σ 1.23 |

**Readings:**

1. **The experiment's own channels are markedly narrower than the
   Hatherly 200 fs priors** (σ_Q2 0.31 vs 0.55; σ_Q3 0.55 vs 1.23) —
   and the measured values are *upper bounds* on the intrinsic widths
   (instrument response + Abel movmean smoothing + 4 % calibSyst
   folded in). Design consequence: the §3.5 table priors are replaced
   by **σ_Q2 = 0.31, σ_Q3 = 0.55 eV** at registration. The n = 1 KED
   σ 0.697 eV must then come mostly from the **channel separation +
   strip selection**, not from intrinsic per-channel width — a
   sharper, more falsifiable mixture forecast for P3.
2. **E_single sharpens to ≈ 0.53 eV (σ ≈ 0.42)** inside the design
   band — one effective Gaussian; the 0.27/0.49/0.70 sub-structure is
   below this 72-point export's resolving power (recorded, not fit).
3. µ readings are KED-frame (Jacobian-shifted below the I(v) peaks:
   Q3 3.91 vs 4.11) — **not a re-anchor**; the channel means stay
   paper-anchored at f·E_C (2.16/4.32), per OQ-B.
4. Area shares (0.66/0.22/0.12) are detection-filtered and are NOT
   channel weights — OQ-C stays anchored on the covariance + power
   series (Q3/(Q2+Q3) here = 0.36, above the 23 % band read;
   projection-free but acceptance-weighted).

**Atlas stance:** nothing adopted; `finc1v725` stands; h405 candidacy +
G4 adjudications open. Next per the design: **P3 composition forecast
+ the probe registration freeze** (weights now un-gated by OQ-K).

## (C) pre-step P3 EXECUTED (2026-07-29, zero MD) — composition forecast on record: C-full lands KE1 mean 0.92–1.04 / SD 0.51–0.55 / above-1.15 43–50 % (reference 48.7 %) with Q3 at 53–61 % of the n = 1 bin; TWO REGISTERED TENSIONS — the blend is two-lump (mode parks at the Q2 toll lump 0.46 vs reference single-mode 0.891) and n1 forecasts 0.18–0.23, short of the 0.31 reference

**Instrument** `tier2atlas_ce_p3_forecast.py` (committed; artifact
`atlas_ce_p3_forecast.csv`). Pure arithmetic over committed inputs:
the §3.5k budget rows as channel proxies (Q2 ← bud226k, Q3 ← bud411k;
KE₁ S_k-corrected to the channel means 2.16/4.32; occupancies
uncorrected — stated proxy bias 4–5 %), the P1 analytic
Poisson-binomial head at the pinned box corners as the
suppressed-conversion kernel (a-independent above v_strip), the P2
measured widths (σ_KE₁|c = √((S_k·σ_c)² + needle²)), converted-KE
ratio 0.774/0.637 (h405-measured, budget-independence assumed),
single channel as a registered assumption (n₁ ∈ [0, 0.05], no row
exists — its §3.5k extrapolated KE₁ < 0). **Oracles passed:** O1
(§3.5k committed rows reproduced to 3 decimals), O2 (the P1 CSV
kernel head equals the analytic recomputation to 1e-9). Weight
vector entering registration (OQ-C, OQ-K closed): **(w_single, w_Q2,
w_Q3) = (0.30, 0.50, 0.20)**.

**Forecast table (band ends over kernel corners × single-channel
assumption; full CSV committed):**

| variant | n₁ | KE₁ mean | KE₁ SD | mode | >1.15 eV | supp | bare | slow-bare | Q3@n=1 |
|---|---|---|---|---|---|---|---|---|---|
| B-only (strip off) | 0.13–0.14 | 0.79–0.87 | 0.45–0.49 | 0.447 | 0.33–0.37 | 0.146 | 0 | 0 | 0.46–0.51 |
| C-full (box) | 0.18–0.23 | 0.92–1.04 | 0.51–0.55 | 0.46 | 0.43–0.50 | 0 | 0.007–0.027 | 0.002–0.007 | 0.53–0.61 |

(reference n = 1 KED: mean 1.302 / mode 0.891 / σ 0.697 / above-1.15
48.7 %; reference n₁_solv 0.31; experimental bare < 1 eV: 1.9 %.)

**Readings:**

1. **The (C) division of labor is confirmed in-arithmetic.** B-only
   reproduces the §3.5k complementarity (n₁ 0.13, residual supp
   0.146); C-full breaks the needle (SD 0.51–0.55 vs 0.034–0.047),
   moves the blend mean into the 1.0 eV anchor band, and lands the
   upper-half weight 0.43–0.50 against the reference 0.487 — the
   selection amplification is now computed, not asserted: the 20 %
   Q3 source share carries 53–61 % of the n = 1 bin.
2. **CP-1/CP-5/CP-6-directions are safe at forecast level:** supp → 0
   by construction of the conversion; ensemble bare 0.7–2.7 % with
   slow-bare ≤ 0.7 % (the Q2-conversion lump at ≈ 0.53 eV is the
   only slow-bare feeder and the pinned kernel keeps it small).
3. **REGISTERED TENSION 1 — the mode.** The forecast n = 1 KED is
   **two-lump** (Q2 toll lump 0.45–0.55 carrying the density mode
   0.46; Q3 lump ≈ 1.28–1.56) while the experimental KED is
   single-mode at 0.891 between the lumps. CP-3 is therefore frozen
   at **blend-mean level** ([0.85, 1.20]); the mode/bimodality is a
   named probe discriminator: a resolved two-lump MD KED against the
   single-mode reference is a (C) shape failure — unless the strip
   toll + in-MD spread (both absent from this forecast) merge the
   lumps, which the probe measures.
4. **REGISTERED TENSION 2 — n₁ magnitude.** The graded kernel
   converts only P(1) ≈ 0.35 of the suppressed class to n = 1, so
   n₁ forecasts 0.18–0.23 — barely above the h405 baseline 0.208 and
   short of the reference 0.31 (the §3.5j all-to-n=1 bracket 0.353
   was the ceiling, not the expectation). CP-4 is frozen at the
   forecast band, not at the reference.
5. Caveats carried: toll not applied (KE forecast toll-high by up to
   ~0.2 eV at full strip); proxy-budget occupancy bias; KE_conv
   budget-independence; the single-channel row is an assumption.

**Registration-band proposal — FROZEN 2026-07-29 (user: "I approve
the bands"); the authoritative frozen table is design doc §8:** CP-1 supp ≤ 0.05; CP-2 KE₁ SD ∈ [0.30, 0.70];
CP-3 KE₁ mean ∈ [0.85, 1.20] AND above-1.15 ∈ [0.35, 0.55];
CP-4 n₁_solv ∈ [0.18, 0.26] (forecast band + one seed-SD);
CP-5 (kill) ensemble slow-bare ≤ 0.01 (bare-row share below 1 eV
≤ 5 %); CP-6 (kill) midHot ∈ [0.85, 1.15]; CP-7 (kill) W₁ ≤ 0.80
(baseline 0.709 + one seed-SD); CP-8 bare KED two-lump qualitative.
MD cells per the design §9: C-full / A-only / B-only / coupling arm
(f_int,Q3 0.0985 vs 0.15), N = 1000, CRN, `drag_state_coupling off`.

**Atlas stance:** nothing adopted; `finc1v725` stands; h405 candidacy +
G4 adjudications open. P1–P3 are complete and the **registration is
FROZEN** (design §8) — the probe build (channel sampler, exit strip,
checkpoint v8) is next, behind `[PROCEED TO IMPLEMENTATION]`.

## (C) probe EXECUTED (2026-07-29, 4 × N = 1000 at the h405 pins, seed 20260729) — **REGISTRATION FAILED: CP-1..4 FAIL, ALL THREE KILLS (CP-5/6/7) FIRED on C-full; CP-6 fired at BOTH coupling-arm ends.** Decomposition: the mixture alone PLACES KE₁ (0.920, above-1.15 0.42 in-band, needle broken) but traps the slow single channel; the strip AS SPECIFIED over-tolls (ε-dominated, 0.5–0.7 eV at the measured knock counts) and FEEDS the suppressed gate (supp 0.18 → 0.55 in A-only). PC-3 failure semantics TRIGGERED — user adjudication next

**Provenance.** Build + probe both this session (2026-07-29), behind the
user's `[PROCEED TO IMPLEMENTATION]`: sampler `sampling/ce_channels.py`
(dedicated stream 0xCE1_2026, frozen three-draw order), strip forms
`physics/exit_strip.py` + step operator `exit_strip_step` (co-moving
carry → E_mass_transfer; count-consistent top-rung D₀ toll → the E_pot
e_bind fold; ε_carry + full-strip E_int residual → E_dissip, OQ-H;
5-term closure bit-tight, unit-tested), config surface `ce_*` /
`exit_strip_*` (biphasic-only guards; partner mask non-optional),
checkpoint **v8** (exactly the OQ-E three fields: `ce_channel`,
`ce_E_m_eV`, `ce_strip_count`; silent exact v7→v8 shim), pair-scale
seam threaded through ion driver / relaxation / detection escape
energetics / t0 (off = `None` = byte-identical), partner-mask scoring
seam (`ce_partner_include_mask` + `include_mask` on the confirmation
read). Off-mode: zero new draws, defaults bit-identical (56 new tests +
full suite green, 2881 passed). Generator `gen_tier2atlas_ce_probe.py`
(unit + cfg-diff oracles passed pre-launch; A-only CRN-paired to the
committed h405 per PC-5), scorer `tier2atlas_ce_probe_table.py`
(O1 pooled-battery + O2 committed-h405-row oracles passed before any
new number), CSV `atlas_ce_probe.csv`.

**The table (partner-masked where channels on; h405 = committed
baseline rescored):**

| cell | scored | trap | supp | n̄ | n₁ | W₁ | midHot | KE₁ | KE₁ SD | >1.15 | slow-bare | strip frac / mean knocks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h405 | 1826 | 0.087 | 0.183 | 3.955 | 0.208 | 0.709 | 0.957 | 0.637 | 0.038 | 0.000 | 0.000 | — |
| cfull | 791 | 0.564 | 0.368 | 3.197 | 0.148 | 1.887 | 0.370 | 0.180 | 0.169 | 0.000 | 0.446 | 0.578 / 11.1 |
| aonly | 1666 | 0.167 | 0.551 | 1.428 | 0.210 | 0.947 | 0.302 | 0.119 | 0.049 | 0.000 | 0.740 | 0.903 / 13.7 |
| bonly | 970 | 0.465 | 0.138 | 6.690 | 0.080 | 2.968 | 1.167 | **0.920** | 0.286 | **0.418** | 0.000 | 0 / 0 |
| cq3hi | 792 | 0.563 | 0.361 | 3.221 | 0.165 | 1.937 | 0.350 | 0.213 | 0.207 | 0.000 | 0.442 | 0.573 / 11.3 |

**Frozen verdicts (design §8, scored on C-full):** CP-1 supp 0.368 >
0.05 FAIL; CP-2 KE₁ SD 0.169 below [0.30, 0.70] FAIL; CP-3 KE₁ mean
0.180 far below [0.85, 1.20] and above-1.15 = 0 FAIL; CP-4 n₁ 0.148
below [0.18, 0.26] FAIL; **CP-5 KILL** (slow-bare 0.446 vs ≤ 0.01);
**CP-6 KILL** (midHot 0.370, and 0.350 at the other coupling end —
fired at BOTH ends); **CP-7 KILL** (W₁ 1.887 vs ≤ 0.80). CP-8 moot at
this level (the bare row is toll-fed, mean 0.45 eV — not the
experimental fast bare lump).

**Decomposition (the attribution cells carry the physics):**

1. **(B) works where the design put it.** B-only KE₁ = 0.920 near the
   1.00 anchor with above-1.15 share 0.418 (CP-3 band [0.35, 0.55])
   and KE₁ SD 0.286 — the needle is broken by source spread alone, and
   the placement confirms the §3.5k S_k lever live in-mixture. The
   fast branch parks only mildly (supp 0.138 — w_Q3-diluted, partial
   §3.5k complementarity).
2. **The mixture's slow tail traps.** trap 0.465 (B-only): the
   E_single ≈ 0.53 eV channel (s_m ≈ 0.2) cannot climb out of
   realistic droplets — ~30 % of molecules land droplet-retained, n₁
   collapses to 0.080 and W₁ blows to 2.97. P3's registered
   single-channel assumption (n₁ ∈ [0, 0.05]) never priced this class.
3. **The strip as specified over-tolls, ε-dominated.** rq4graded outer
   rungs are tiny (Σ(21)−Σ(7) = 0.136 eV), so the D₀ part of the toll
   is small — but the measured knock counts (mean 13.7, p90 = 20 in
   A-only: P₀ ≈ 0.25–1 for ordinary exiters, G ≈ 1 above j₀, OQ-J
   multi-crossing ratchet) make the box-mid ε_carry = 0.025 eV/He
   worth 0.34–0.50 eV alone → total toll 0.5–0.7 eV ≈ 2.5–3.5× the
   design's "≈ 0.2 eV" anchor. KE₁ 0.637 → 0.119 (A-only); the toll
   also converts marginal escapers to retained (trap 0.087 → 0.167).
   **The ε box is internally inconsistent with its own 0.2 eV anchor
   at the measured knock counts** (consistency needs ε ≲ 0.005).
4. **The strip FEEDS the suppressed gate instead of retiring it.**
   A-only supp 0.183 → 0.551: 918 ions land suppressed at n_det 2–12
   with E_int ≈ 0.15 eV — the strip collapses Σ(n) under a survivor
   whose E_int it never touches, closing the self-bound gate at small
   n. CP-1's "conversion by construction" was occupancy arithmetic;
   the live E_int/Σ(n) gate inverts it. Structural: no value in the
   frozen boxes undoes it.
5. **Why P1/P3 missed:** P1's counterfactual stripped the suppressed
   class only, with no energy toll; P3 composed P1's kernel with
   toll-free channel proxies. The live §3.3 strip acts on every
   outbound crossing of every ion — the first time the general
   population, the toll, and the gate ran together.

**PC-3 (design §11) consequences — TRIGGERED, adjudication with the
user:** CP-5 fired at a = 2 and CP-7 fired ⇒ family-level kills — per
the pre-commitment the **(A) strip design v1 stops** and the question
returns to the §3.5j alternatives; CP-6 fired at both bracket ends ⇒
the **current (C) f_int wiring closes** (the mixture itself is NOT
killed — read 1 is its measured success). Candidate follow-ups for
adjudication (NOT taken): ε → ~0 + suppressed-class-gated strip
variant; an E_int co-strip term (the gate-inversion fix); pricing the
single channel's trap class into the mixture (weights are
covariance-anchored, so this is a fate question, not a weight refit).
PC-4: the m/q-127 ask does NOT fire (no CP pass). Tension-2 note
moot at this outcome level.

**Atlas stance:** nothing adopted, nothing moved — `finc1v725` stands,
h405 candidacy + G4 adjudications open, the scalar budget did NOT
retire (it stands with the standing points). The (C) run dirs are
instrument runs (`tier2atlas_conf270_ce*`).
