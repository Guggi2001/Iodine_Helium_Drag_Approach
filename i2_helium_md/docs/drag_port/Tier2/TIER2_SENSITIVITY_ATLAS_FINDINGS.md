# Tier 2 — Sensitivity Atlas: Findings

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

### Status of the D4 form axis after Step 1

Registered priors unchanged: lq stays rejected (shared +0.0339, held-out
FAIL), power law stays cubic-equivalent, linear+cubic stays a → 0. New
knowledge from the zero-cost gate: Padé excluded as a cap-replacement
*by arithmetic* (no fit needed); subtractive v_f is trace-constrainable
and awaits a promotion decision; sub-2.5 Å/ps stays TDDFT-blind. Step 2
(twin sweeps) can proceed on the standing form list; Step 3 (MD
spot-checks / enum build) has no candidate yet that survives to it.
