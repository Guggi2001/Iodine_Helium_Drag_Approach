# Per-case (individual 9A / 18A) Method-B fits for `linear_quadratic` and `power_law`

**Status:** planned, not yet implemented. Drag-phase implementation work — it
requires the explicit `[PROCEED TO IMPLEMENTATION]` trigger before any code is
written. Decisions below confirmed with the user 2026-06-14.

## Context

Method-B trajectory matching has so far produced, per drag form:

- **`linear_cubic`** — per-case bundles for 9A and 18A (§8) **and** the shared
  joint refit (§9, `shared_pure_cubic` = production candidate).
- **`linear_quadratic`** and **`power_law`** — **only** the shared joint refits
  (§10.7: `lq_shared_3param`, `lq_shared_pure_quadratic`, `pl_shared_3param`).

Goal: fill the **missing diagonal** — individual per-case (9A-only, 18A-only)
optimal fits for `linear_quadratic` and `power_law`, parallel to what §8 gave
`linear_cubic`.

**Why this is worth doing (diagnostic value).** Per §10.4.1 the exponent
leverage lives mainly in the *cross-case speed-scale difference*, not the
in-window shape. So per-case fits answer questions the shared fit cannot:
- per-case **`n̂`** for `power_law` — does 9A alone vs 18A alone prefer `n≈2`
  or `n≈3`, and how poorly is `n` identified *within a single curve* (the
  predicted large `n_err`)?
- does each case independently collapse to **pure-quadratic** (`a→0`) the way
  the shared `lq` fit did?

**Methodological status (must be stamped honestly).** These are
**calibrated-not-validated, single-case** fits. For these forms the §9 cross-case
held-out axis is already *spent*; a per-case fit trivially matches its own curve
(the §5/§8 warning). 9A additionally carries the transverse-contamination flag.
→ Diagnostic artifacts only; **never wired to presets**; the incumbent stays
`shared_pure_cubic`.

## Decisions (confirmed with user)

- **Artifact = loadable bundles** (parallel to §8), one `fit_parameters.json` per
  `(case, variant)`, loadable through `load_drag_coefficients`.
- **`linear_quadratic` = both variants** per case (`lq_shared_3param` +
  `lq_shared_pure_quadratic`). `power_law` = `pl_shared_3param` (free-`n`, only
  variant). → **6 bundles**: {9A, 18A} × {lq_3param, lq_pure_quadratic, pl_3param}.
- **Driver = thin sibling script** (`method_b_per_case_form.py`) reusing
  `form_phase_common`, **not** a generalization of the §8 `method_b_extraction.py`.

### Why a sibling, not a flag on `method_b_extraction.py`

The §8 driver's *structure* is the template (fit one case → write one loadable
per-case bundle), and the form-generic engine it would call already exists
(`fit_form_trajectory_matching`, `form_sensitivity_halfwidths`) — those are
reused either way. But `method_b_extraction.py` is wired to `linear_cubic` in
three load-bearing ways, so "set a variant" is not enough:

1. **Fit fn + param vector.** It calls `fit_trajectory_matching` →
   `TrajectoryMatchingFit` with `.best.a/.best.b` and an explicit `(a,b,E)`
   vector in `dry_run` and the seed-sweep `warm_start_refit`. Forms return
   `SharedFormTrajectoryMatchingFit` with `.best.coefficients` (`{C,n}`/`{a,c}`).
2. **Anchor source (the real blocker).** §8 conditions on `setup.a0/b0`
   (preset-derived). The form phase **forbids** that — presets were re-wired to
   `shared_pure_cubic`, so `setup.a0 = 0` — and **mandates** the LOCKED
   pre-registered anchors from `form_phase_common` (§10.4.1).
3. **Writer + uncertainty convention.** `write_fit_parameters` is the
   `linear_cubic` schema; §8 builds its band from a **seed sweep**, while the
   form phase deliberately uses **sensitivity-only** (`..._seed_sweep_omitted_per_s8`).

Folding forms into the §8 file therefore means form branches through the whole
fit/anchor/sensitivity/writer path **plus** a new writer anyway — while editing
the frozen §8 evidence script and mixing two anchor sources + two uncertainty
conventions behind conditionals. The sibling keeps §8 (and its committed
`linear_cubic` per-case bundles) untouched and reuses `form_phase_common` cleanly.
A form-capable single-case **writer is unavoidable in every option** (the form
writer `write_shared_form_fit_parameters` refuses single-case per the 2026-06-12
decision; the `linear_cubic` writer is the wrong schema).

## What already exists and is reused (no reinvention)

- `fit_form_trajectory_matching(setup, variant=, anchors=, ...)`
  (`i2_helium_md/extraction/trajectory_matching.py:1446`) — **the single-case
  form fit**. Currently used only as the 18A-only Stage-1 analog; call it per case.
- `evaluate_form_objective`, `form_sensitivity_halfwidths` — single-case-safe
  (the shared-`E_bind` check passes trivially for one case). For `power_law`,
  convert the `gamma_ref` half-width to `C_err = gamma_ref_err / v_ref**(n-1)`.
- `form_phase_common.py` — LOCKED §10.4.1 anchors (`A0`, `C0_LQ`, `PL_C0`,
  `PL_N0`, `V_REF_APS`, `N_BOUNDS`), `N`, `SEED`, `ION_TIME_PS`, `DT_ION_PS`,
  `UNCERTAINTY_MODEL`, `ANCHOR_PROVENANCE`, `stamp`, `dump_json`,
  `fit_is_trapped`. Anchors are conditioning-only → reused unchanged.

## Why new code is needed

`write_shared_form_fit_parameters` (`trajectory_matching.py:1733`) **deliberately
refuses single-case fits** (line 1759, the 2026-06-12 "Stage-1 analog is
record-only" decision). That refusal must stay intact. The per-case loadable
bundle therefore needs a **dedicated writer**, parallel to the `linear_cubic`
per-case writer `write_fit_parameters` (`:576`) but form-generic.

## Implementation

### 1. New per-case form writer — `i2_helium_md/extraction/trajectory_matching.py`

Add `write_per_case_form_fit_parameters(fit, out_dir, *, transverse_contaminated,
anchor_provenance, v_ref_Aps=None)`:

- Accept a single-case `SharedFormTrajectoryMatchingFit` (`len(fit.cases) == 1`).
- Refuse a trapped fit (`penalty != 0`) and unfilled `coeff_errs` /
  `uncertainty_model` — mirror both existing writers.
- Payload = the §6 loadable schema, form-generic: `extraction_method`,
  `form`, the form's raw coeff keys (`_REQUIRED_COEFF_KEYS[form]`) + per-key
  `<k>_err`, `effective_binding_energy_I_ion_eV` (+ err), single-case window,
  single `reference_file`, `variant`, `anchors`+provenance, optimizer/start
  minima, per-case RMSE/escape. `power_law` additionally stamps the pivot block
  (`v_ref_Aps`, `gamma_ref`) for re-derivability (mirror the shared writer).
- **`flags`**: `radial_projection_convention: True`,
  `transverse_contaminated_non_radial_reference: transverse_contaminated`,
  `full_window_heldout_window_axis_forfeited: True`, **plus new honesty flags**
  `per_case_calibrated_not_validated: True`,
  `cross_case_axis_not_applied: True` (single-case ⇒ no cross-case content;
  not a preset candidate).
- Leave `write_shared_form_fit_parameters` untouched.

### 2. New driver — `scripts/extraction/method_b_per_case_form.py`

Sibling of `method_b_form_refit_*.py`, importing `form_phase_common`. Pattern
mirrors `method_b_extraction.py` (the §8 per-case driver), not the joint refit:

- USER SETTINGS: `CASES = ("18A", "9A")`, `VARIANTS = (LQ_SHARED_3PARAM,
  LQ_SHARED_PURE_QUADRATIC, PL_SHARED_3PARAM)`, `DRY_RUN`, `MAXFEV_PER_START`.
- Build each case setup **once** (cache the neutral run), reuse across variants.
- Per `(case, variant)`: select anchors per variant from `form_phase_common`
  constants (`{a0,c0}` for lq, `{C0,n0,v_ref_Aps}` for pl); run
  `fit_form_trajectory_matching`; abort-and-skip if trapped; compute
  `form_sensitivity_halfwidths`; fill `fit.coeff_errs` / `fit.e_bind_err_eV` /
  `fit.uncertainty_model` (= `UNCERTAINTY_MODEL`, sensitivity-only, seed sweep
  omitted per §8); write the bundle to
  `data/reference/drag/<case>/trajectory_matching/<variant>/`.
- DRY_RUN: one form objective eval ×2 + bitwise-determinism assert, no writes
  (mirror the existing drivers).
- Print per-fit summary incl. `n̂` + `n_err` (pl) and the `a` value (lq collapse).

**Output path convention:** `data/reference/drag/<case>/trajectory_matching/<variant>/fit_parameters.json`
— extends the §8 per-case location (`<case>/trajectory_matching/…`) with a
variant subdir; no collision with the existing `linear_cubic`
`<case>/trajectory_matching/fit_parameters.json`.

### 3. Tests — `tests/test_extraction_trajectory_matching.py`

- New `TestPerCaseFormBundleWriter`: for each form, a single-case fit (stub
  `run_fn` / synthetic `SharedFormTrajectoryMatchingFit`, mirroring existing
  stub-law tests) → writer round-trips through `load_drag_coefficients`;
  bundle carries the per-case + transverse flags; `power_law` records the pivot
  block; writer refuses trapped fit and unfilled `coeff_errs`.
- Keep tests fast (no N=50 ion runs, no figures) per the testing rules.

### 4. Docs

- `METHOD_B_trajectory_matching_extraction.md`: new subsection (e.g. **§10.8
  Per-case alternative-form fits — diagnostic**) — motivation (per-case `n̂`,
  per-case pure-quadratic collapse), method, artifact paths, the
  calibrated-not-validated / spent-cross-case status, and the run outcomes.
- `drag_migration_log.md`: delivery + decision entry (loadable bundles, both lq
  variants, new writer, not-preset-wired).

### 5. Run the extraction (implementation-time, heavy compute)

After the trigger: DRY_RUN determinism check, then the full run producing the 6
bundles. Per-case is one ion run per objective eval (cheaper than the joint
fit's two), but still N=50 × multi-start × `MAXFEV_PER_START` per fit — a long
compute step. Record the resulting coeffs / RMSE / escape / `n̂`+`n_err` into the
doc + migration log.

## Critical files

- `i2_helium_md/extraction/trajectory_matching.py` — new writer (alongside `:576`
  / `:1733`); reuse `fit_form_trajectory_matching:1446`,
  `form_sensitivity_halfwidths:1642`.
- `scripts/extraction/method_b_per_case_form.py` — new driver (reuses
  `form_phase_common.py`).
- `tests/test_extraction_trajectory_matching.py` — new writer tests.
- `METHOD_B_trajectory_matching_extraction.md`, `drag_migration_log.md` — records.
- New artifacts: `data/reference/drag/{9A,18A}/trajectory_matching/{lq_shared_3param,lq_shared_pure_quadratic,pl_shared_3param}/fit_parameters.json`.

## Verification

1. `DRY_RUN=True` run → bitwise-determinism assert passes (ion drag path RNG-free).
2. Narrow tests: `python -m pytest tests/test_extraction_trajectory_matching.py -q`
   (new `TestPerCaseFormBundleWriter` + existing pass).
3. Each written bundle loads cleanly via
   `presets.load_drag_coefficients(<bundle_dir>, expected_m_eff_amu=202.953908)`.
4. Full suite green (705+/0) — no preset/loader regressions; confirm presets
   still resolve `shared_pure_cubic` (unchanged).
5. Sanity-read the 6 bundles: per-case `n̂`±`n_err` (pl), `a→0` collapse (lq),
   escape = 1.0, flags present.

## Out of scope / unchanged

- Presets stay on `shared_pure_cubic` (no re-wiring); incumbent unchanged.
- `write_shared_form_fit_parameters` single-case refusal stays intact.
- No Tier-1 start, no VMI, no new physics, no constants/RNG/schema changes.
- LOCKED §10.4.1 anchors reused as-is (not re-derived).
