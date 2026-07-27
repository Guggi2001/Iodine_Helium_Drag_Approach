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
