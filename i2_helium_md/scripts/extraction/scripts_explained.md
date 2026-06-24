# The Method-B extraction scripts, explained

The scripts in this directory are the **drivers** of the Method-B
trajectory-matching extraction (`METHOD_B_trajectory_matching_extraction.md`);
all heavy lifting lives in `i2_helium_md/extraction/trajectory_matching.py`
(documented in `docs/extraction/trajectory_matching_module.md`). They were
run to completion in June 2026 and their verdicts are **recorded** — this
file explains what each script does, what it wrote, and what to watch out
for if one is ever run again.

```
method_b_extraction.py                    per-case joint {a, b, E_bind} fit   (plan B4)
method_b_cross_case_check.py              provisional cross-case held-out     (plan B5/B7)
method_b_shared_refit.py                  §9 shared-form joint refit, 2 stages (§9)
form_phase_common.py                      §10 shared helpers + LOCKED constants
method_b_form_refit_linear_quadratic.py   §10 linear_quadratic family driver  (§10)
method_b_form_refit_power_law.py          §10 power_law family driver          (§10)
```

Intended run order (as it happened historically): extraction per case →
cross-case check → (after the §8 failure and the §9 decisions) the shared
refit.

## Status box — read this first

- **All production runs are complete** (§8/§9 on 2026-06-11, the §10
  form-discrimination phase on 2026-06-14) and the verdicts are recorded in
  the artifacts and in `drag_migration_log_tier0.md`. Per the **first-runs rule**
  the pre-registered thresholds were fixed before the runs and the recorded
  verdicts stand — do not re-run a fit to "check" it.
- **Outcome in one line each:** the per-case fits were in-window excellent
  but **failed** the provisional cross-case bands (weak-`a`
  identifiability; bundles flagged not-yet-usable); the shared refit
  **passed** every pre-registered §9.4 band — the law is effectively
  pure-cubic, Tier 1 ungated; the §10 form phase then **confirmed the
  incumbent** — the free-`n` `power_law` fit recovers `n̂ ≈ 2.93 ≈ 3` and
  forced-`v²` `linear_quadratic` is measurably worse, so no alternative form
  beats `shared_pure_cubic` (no escalation, presets unchanged).
- **The §9.4 bands supersede the cross-case check's provisional bands.**
  `method_b_cross_case_check.py`'s check (i) (the 18 Å A↔B `a`-ratio band)
  is **dead** under the weak-`a` finding — it would auto-fail any `a → 0`
  fit and is not re-applied to the shared fit (2026-06-11 decision).
- **Post-re-wiring caveat (2026-06-12).** At delivery time the drag
  presets carried the per-case **Method-A** bundles, and two pieces of the
  scripts leaned on that: `build_case_setup` reads its normalization
  anchors `a0/b0` from the preset's coefficients, and
  `method_b_shared_refit.py` takes its anchors from `setups["18A"].a0/b0`
  ("the 18A drag preset wires exactly that bundle" — true then). Since the
  preset re-wiring (METHOD_B §10.1) the presets load the shared
  `shared_pure_cubic` bundle, so a fresh setup yields `a0 = 0.0`: the
  shared refit's positive-anchor check would raise, and the per-case fit's
  `[0.1, 10]·a0` bounds would collapse. **Resolution (user decision
  2026-06-12, METHOD_B §10.4):** future fit drivers carry their anchors as
  **pre-registered constants** — values read off the Method-A bundles once
  (e.g. `data/reference/drag/18A/linear_and_cubic/fit_parameters.json`),
  locked with provenance in the USER SETTINGS block at the pre-run session
  — never the preset-derived `setup.a0/b0` and never a runtime bundle read.

## Common patterns

All three scripts follow the repo's script conventions:

- **`# USER SETTINGS` block** at the top — edit constants, then run the
  file directly (no CLI arguments).
- **`PROJECT IMPORT SETUP`** — inserts the repo root into `sys.path` so
  the scripts run from any working directory.
- **`DRY_RUN = True`** (the two fit drivers) evaluates the objective twice
  at a fixed point, **asserts the two results are bitwise identical**, and
  writes nothing. This is the determinism guard: the drag ion path must be
  RNG-free, and any stochasticity leaking in fails loudly here before an
  expensive fit is started. The recorded production runs were preceded by
  a passing DRY_RUN (~1.6 s per joint evaluation = 2 ion runs).
- **Provenance stamps** — every written artifact carries date + git branch
  (via the module's `_git_branch`, which degrades to `"unknown"` rather
  than crash) and the thresholds it was judged against, so a later
  re-derivation can re-judge the same numbers.
- Shared numerical settings: `N = 50` ensemble (the committed Tier-0
  evidence point — N=50 vs N=2000 identical to 0.001 Å/ps), neutral seed
  `20260604`, ion stage 20 ps @ 0.01 ps.

---

## 1. `method_b_extraction.py` — the per-case fit driver (plan B4)

Runs the joint `{a, b, E_bind}` fit for **one case** (`CASE = "9A"` or
`"18A"`; fit 18 Å first — the trustworthy clean-radial case) against the
same-smoothed long reference over the full post-dynamic-start window.

Flow:

1. `build_case_setup` — one neutral run, cached for every evaluation;
   prints the reference window and anchors.
2. (`DRY_RUN`: two evaluations at the Method-A-era start point +
   determinism assert, then exit without writing.)
3. `fit_trajectory_matching` — the multi-start Nelder-Mead fit (`E_bind`
   pre-scan + two ridge-probing extra starts), printing each start's
   minimum and the best point.
4. **Abort on trapping:** if the best fit still carries an escape penalty,
   the script exits 1 and writes nothing — a binding↔drag pair that traps
   ions must never become a bundle.
5. **Seed sweep** — warm-started refits (reduced `maxfev = 120`) on four
   fresh neutral seeds; the parameter spread across seeds is one half of
   the uncertainty band. (Empirical §8 finding: the HeDFT-comparison
   presets are seed-insensitive — std ≈ 0 — so the band is carried by the
   sensitivity scan in practice.)
6. **Sensitivity scan** — `sensitivity_halfwidths` (parameter change for a
   10% objective rise).
7. **Band = `max(seed std, sensitivity half-width)` per parameter**,
   labeled `uncertainty_model = "seed_sweep_std_plus_rmse_sensitivity"` —
   a deliberately conservative *sensitivity band*, not a statistical CI.
8. **Ridge warning:** if the multi-start minima disagree beyond 3× the
   band, a loud warning flags the live drag↔binding degeneracy ridge (the
   bundle is still written — the spread is recorded in
   `optimizer.start_minima` — but it must be reviewed before any preset
   re-wiring).
9. `write_fit_parameters` →
   `data/reference/drag/<CASE>/trajectory_matching/fit_parameters.json`
   (sibling of `linear_and_cubic/`; Method-A artifacts untouched), with
   the 9 Å transverse-contamination flag set by `CASE`.

The closing note is part of the method: a fit produced here is
**calibrated, not validated** — the held-out validation is the separate
cross-case script (and ultimately the post-Tier-1 VMI tier).

**Historical outcome (2026-06-11, recorded in METHOD_B §8):** 18 Å
`a = 1.456` (pinned at the then-`0.1·a0` lower bound), `b = 3.316`,
`E_bind = 0.0709 eV`, RMSE 0.0968 Å/ps; 9 Å `a = 7.50` (multi-start ridge
spanning 3.68–7.50 at Δobjective ≈ 2e-4), `b = 1.751`, `E_bind = 0.154 eV`,
RMSE 0.0409 Å/ps — both fully escaping, both later **failing** the
provisional cross-case bands.

## 2. `method_b_cross_case_check.py` — the provisional held-out check (plan B5/B7)

Pure post-processing of the two written bundles (no MD): the load-bearing
held-out validation that the full-window calibration left in scope before
the VMI tier. Three checks against **named provisional thresholds** (the
§6.10 first-runs rule — set before the first run, recorded with the
verdict):

| check | threshold | rationale |
|---|---|---|
| (i) 18 Å A↔B coefficient agreement | `a_B/a_A`, `b_B/b_A` ∈ [0.5, 2.0] | the clean-radial fit should land near its Method-A region; large drift = overfit signature. The 9 Å A↔B check is **dropped** — its Method-A artifact is provenance-broken (2026-06-11 decision). |
| (ii) residual-quality asymmetry | RMSE 9 Å / 18 Å ≤ 3.0 | 9 Å is expected worse (model-dimensionality cost); an extreme ratio means transverse contamination dominated. |
| (iii) cross-case `E_bind` consistency | \|E_9A − E_18A\| ≤ 0.05 eV | the binding stands in for the same I⁺/He ejection physics in both cases. |

Writes `held_out_validation.json` next to each case's bundle (keeping
`fit_parameters.json` immutable post-fit), with one-line verdicts per case
and the `"vmi": "pending"` stamp. Degrades gracefully if the 9 Å bundle is
absent (checks (ii)/(iii) deferred).

**Historical outcome (2026-06-11):** (i) **FAIL** (`a_B/a_A = 0.10`; the
`b` ratio 1.62 passed), (ii) PASS (ratio 0.42 — note the 9 Å fit was
*better* than clean-radial 18 Å, itself the §5 warning signature),
(iii) **FAIL** (0.083 eV). Honest reading: dominated by weak-`a`
identifiability, not proven overfit — which led to the §9 shared refit.
**Superseded:** the §9.4 bands replace these provisional bands for judging
the shared fit; `method_b_shared_refit.py` appends its `shared_form_refit`
verdict block into these same `held_out_validation.json` files (the
provisional records are left intact above it).

## 3. `method_b_shared_refit.py` — the §9 shared-form joint refit driver

The script that resolved the §8 failure. Two stages, run **in order** so
the held-out content is recorded *before* the joint fit consumes the
cross-case axis. All thresholds are the **pre-registered §9.4 bands**,
fixed before any run as named constants in the USER SETTINGS block
(`S1_PRED_RMSE_MAX_APS = 0.45`, `S2_RMSE_18A_MAX_APS = 0.19`,
`S2_RMSE_9A_MAX_APS = 0.45`, escape = 1.0, `T_A0_APS = 0.005`).

Flow:

1. Build both case setups (one neutral run each); take the anchors from
   the 18 Å setup (see the **post-re-wiring caveat** in the status box —
   at delivery time this was the 18 Å Method-A bundle, `a0 = 14.5556`,
   `b0 = 2.0534`).
2. (`DRY_RUN`: one joint evaluation ×2 + bitwise-determinism assert.)
3. **Stage 1 — the untouched 9 Å prediction** from the delivered §8 18 Å
   bundle (read straight from its `fit_parameters.json`; **no refit**, per
   the pre-run clarification). One `evaluate_objective` call on the 9 Å
   setup + band scoring → `stage1_prediction.json`. **Recorded-only,
   non-gating** — only the Stage-2 bands carry verdict power.
4. **Stage 2 — two joint fits:** the primary `shared_3param` (`a` free,
   lower bound 0) and the `shared_pure_cubic` variant (`a ≡ 0`), each via
   `fit_shared_trajectory_matching` with the same multi-start scheme.
5. **`T_a0` classification:** `obj(a≡0) − obj(a free) ≤ 0.005 Å/ps` ⇒ the
   variants are equivalent and pure-cubic is recorded as the empirical
   conclusion on `a`'s identifiability (a classification, not a gate; the
   sign convention is the corrected §9.4 one — nested models make the
   difference ≥ 0 up to optimizer noise).
6. **Pre-registered Stage-2 band checks** on both variants
   (`_band_checks`); the primary's `all_pass` is the verdict bit.
7. **Stage-1↔Stage-2 parameter shift** — qualitative, recorded-only (a
   large shift when 9 Å enters the objective would flag 9 Å dragging the
   law toward its transverse contamination).
8. **`diagnostic_4param` only on primary failure** (shared `{a, b}`,
   per-case `E_bind` — localizes drag-law vs. binding failure; written to
   `diagnostic_4param.json`, never as a loadable bundle).
9. **Sensitivity + bundle writing for both shared variants** —
   `joint_sensitivity_halfwidths` fills the band
   (`uncertainty_model = "rmse_sensitivity_only_seed_sweep_omitted_per_s8"`;
   the §8 sweep showed seed-insensitivity), then
   `write_shared_fit_parameters` to the per-variant subdirectories. A
   trapped variant is skipped, not written.
10. **`verdict.json`** — the complete machine-readable record: thresholds
    + gating provenance, the full Stage-1 record, both Stage-2 summaries
    with their band checks, the `T_a0` conclusion, the shift record, the
    diagnostic (if run), `primary_pass`, bundle paths, the §9.4 verdict
    mapping, and the `"vmi": "pending"` stamp.
11. **Per-case `held_out_validation.json` update** — appends the
    `shared_form_refit` block (scored against §9.4; the provisional §8
    records above it are left intact, superseded).

Presets are **not** re-wired by this script (that was a separate user
decision, executed 2026-06-12 — see METHOD_B §10.1).

**Historical outcome (2026-06-11, METHOD_B §9.7): PASS on every band.**
Stage 1: 0.2685 Å/ps ≤ 0.45, full escape. Stage 2 `shared_3param`:
`a → 0` bound (all three starts agree), `b = 2.516`, `E_bind = 0.1168 eV`,
18 Å 0.1345 ≤ 0.19, 9 Å 0.1240 ≤ 0.45, escape 1.0/1.0;
`shared_pure_cubic` indistinguishable (`T_a0 = −2.3e-7`). Tier 1 ungated;
the diagnostic was not needed.

## 4. `form_phase_common.py` + the two `method_b_form_refit_*.py` drivers (§10)

The §10 **alternative-form discrimination** phase: realize the `power_law`
and `linear_quadratic` families and run each through the **same** §9
shared-form joint-refit machinery against the pure-cubic incumbent, to test
whether any alternative drag exponent fits the two cases as well or better.
This is **model selection on seen data** under a pre-registered protocol (the
§9 cross-case axis was already spent by the Stage-2 joint fit) — a legitimate
ranking, **not** fresh held-out validation; the winner's external test stays
VMI after Tier 1.

`form_phase_common.py` is the shared module (project rule 1 — no copy-paste
between the two siblings). It carries:

- the **LOCKED §10.4.1 pre-registered constants** as named constants with
  provenance: the reused §9.4 Stage-2 bands (18 Å ≤ 0.19, 9 Å ≤ 0.45, escape
  1.0), the two `T_form` thresholds (`T_FORM_EQUIV_APS = 0.005`,
  `T_FORM_BETTER_APS = 0.013`), the exact incumbent objective
  `INCUMBENT_OBJECTIVE_APS = 0.1292649398514104`, the `power_law` `n` bounds
  `(1, 4)` and pivot speed `V_REF_APS = 3.0`, and the 18 Å-only anchor
  constants `A0`/`C0_LQ`/`PL_C0`/`PL_N0` (read off the Method-A bundles once,
  locked here — **never** a runtime bundle read, since the re-wired presets
  carry `a0 = 0`);
- `band_checks` (the reused §9.4 bands, made form-generic — scores a
  `JointFormObjectiveResult`), `fit_summary`, `fit_is_trapped`;
- `classify_t_form` — the **two-threshold, asymmetric** classifier:
  `Δ = best-variant objective − incumbent` → `genuinely_better` (escalate,
  *and* must pass the bands) / `marginally_better` (recorded, no escalation) /
  `equivalent` (exponent degeneracy recorded, pure-cubic stays) / `worse`
  (rejected). The asymmetry keeps a Δ inside the objective's own flatness scale
  from being stamped "genuinely better";
- `update_form_comparison_verdict` — read-modify-write of the combined
  `form_comparison_verdict.json` (keyed by family) next to the §9
  `verdict.json`.

Each sibling driver (`method_b_form_refit_linear_quadratic.py`,
`method_b_form_refit_power_law.py`) mirrors the §9 driver flow: build both
setups → DRY_RUN (joint eval ×2 + bitwise-determinism assert) → **Stage-1
analog** (fit the family's full variant on 18 Å only, score the untouched 9 Å
prediction vs 0.45 Å/ps for comparability with pure-cubic's 0.2685 — recorded,
**non-gating**, record-only JSON, no bundle) → **Stage-2 shared joint fits**
(lq runs `lq_shared_3param` + `lq_shared_pure_quadratic` then the `T_a0`-analog
equivalence; pl runs `pl_shared_3param` only, in the locked **pivot**
parameterization `(γ_ref, n)`) → §9.4 band checks → `T_form` classification on
the family best variant → sensitivity (`form_sensitivity_halfwidths`; pl
converts the pivot half-width to `C_err = γ_ref_err / v_ref^(n̂−1)` and stamps
`n_err` directly — the sharpest single number this phase produces) → bundle
writes (skipping trapped fits) → per-family + combined verdicts. The anchors are
passed as the LOCKED constants, **never** `setup.a0/b0` (the post-re-wiring
caveat in the status box — a fresh setup reads `a0 = 0`).

**Production outcome (2026-06-14, METHOD_B §10.7): incumbent CONFIRMED.**

- `power_law` (free `n`): all 4 starts → **`n̂ = 2.927`** (C ≈ 2.835,
  E_bind ≈ 0.113 eV), objective 0.129171 → Δ = −0.0001 → **EQUIVALENT** to
  pure-cubic; bands PASS; `n_err` half-width 0.279. The free-exponent fit
  **independently recovers `n ≈ 3`**, not Method-A's `n ≈ 2.06`; Stage-1
  analog 9 Å predict 0.411 ≤ 0.45 PASS.
- `linear_quadratic` (forced `n = 2`): both variants pass the §9.4 bands but
  collapse to the pure-quadratic corner (`a → 0`, `c ≈ 12.8`, E_bind ≈ 0.048)
  with objective 0.16315 → Δ = +0.0339 → **WORSE, rejected-by-objective**;
  Stage-1 analog 9 Å predict 0.699 > 0.45 FAIL (non-gating).
- **No family beats `shared_pure_cubic`; no escalation; presets NOT re-wired.**
  The §10.2 exponent tension is resolved in favour of `n = 3` — the trajectory
  objective discriminates the exponent (unlike the linear term, §9).

---

## Artifact map

```
data/reference/drag/
├── 18A/
│   ├── linear_and_cubic/fit_parameters.json      Method-A (frozen; anchor source)
│   └── trajectory_matching/
│       ├── fit_parameters.json                   §8 per-case fit (flagged not-yet-usable)
│       └── held_out_validation.json              provisional verdicts + shared_form_refit block
├── 9A/
│   ├── linear_and_cubic/fit_parameters.json      Method-A (restored 2026-06-12; frozen)
│   └── trajectory_matching/                      (same two files; 9 Å flags set)
└── shared/trajectory_matching/
    ├── stage1_prediction.json                    §9 Stage-1 held-out record (non-gating)
    ├── verdict.json                              the complete §9 run record
    ├── shared_3param/fit_parameters.json         §9 variant bundle (a at the 0 bound)
    ├── shared_pure_cubic/fit_parameters.json     §9 variant — PRODUCTION CANDIDATE
    │                                             (wired into the presets 2026-06-12)
    ├── stage1_analog_linear_quadratic.json       §10 lq Stage-1 analog (non-gating)
    ├── stage1_analog_power_law.json              §10 pl Stage-1 analog (non-gating)
    ├── verdict_linear_quadratic.json             §10 lq family verdict (rejected: worse)
    ├── verdict_power_law.json                    §10 pl family verdict (equivalent → n≈3)
    ├── form_comparison_verdict.json              §10 combined verdict (keyed by family)
    ├── lq_shared_3param/fit_parameters.json      §10 rejected bundle (a→0, pure-quad)
    ├── lq_shared_pure_quadratic/fit_parameters.json  §10 rejected bundle (n=2 forced)
    └── pl_shared_3param/fit_parameters.json      §10 bundle: n̂=2.927 ≈ pure-cubic
                                                  (EQUIVALENT — incumbent confirmed)
```
