# The `trajectory_matching.py` module

`i2_helium_md/extraction/trajectory_matching.py` — the Method-B
trajectory-matching drag extraction (per-case joint `{a, b, E_bind}` fit
**and** the §9 shared-form joint refit machinery).

## What problem does this file solve?

Method A (the original pipeline, `Drag_extraction_code.md`) regresses the
force-balance residual `|F_drag|` against speed — it fits an *intermediate*
quantity and never sees the forward-integrated trajectory. Method B
(`METHOD_B_trajectory_matching_extraction.md`) inverts that: the fit
objective **is** the trajectory match. This module implements Method B as a
named, provenance-stamped extraction method:

- fit the `linear_cubic` drag coefficients `{a, b}` **jointly with the
  effective droplet binding** `E_bind` (the binding↔drag coupled pair,
  METHOD_B §5.5 / `TIER0_FINDINGS.md` "Correct drag traps the ions"),
- by minimizing the **forward-integrated, from-onset, in-window RMSE** of
  the ensemble-mean clean-atom speed `|v2|` against the same-smoothed long
  reference (`data/reference/drag/<case>/velocity_smoothed/cleaned_data_long.csv`),
- plus an **escape penalty** that steers the optimizer away from binding
  depths that trap the ions.

The module also carries the **§9 shared-form joint refit**: one
size-independent drag law fit jointly across both HeDFT cases (9 Å + 18 Å),
which is the over-constrained (3 parameters / 2 trajectories), falsifiable
version of the fit. Its 2026-06-11 run **passed** every pre-registered §9.4
band and settled the law as effectively pure-cubic (`γ = g·b·v²`).

Authoritative specs and records:

- `METHOD_B_trajectory_matching_extraction.md` — the method (§1–§6), the
  per-case outcome (§8), the shared refit spec + outcome (§9, §9.7), and
  the follow-up phase (§10).
- `drag_migration_log_tier0.md` — delivery records and user-locked decisions.
- `docs/physics/drag_module.md`, `docs/physics/baoab.md` — the physics this
  module drives through (`γ(v)` and the BAOAB O-step).

## Position in the dependency chain

```
presets (single_pulse_N2000_drag / _18Angst_drag)   ← per-case SimConfig base
   ↓
simulation/neutral.run_neutral_propagation          ← run ONCE per (case, seed)
   ↓                                                  (cached in CaseSetup)
simulation/ion.run_ion_propagation                  ← once per objective
   ↓                                                  evaluation (BAOAB drag path)
postprocess/compare_trajectories.compare_speed_to_reference
postprocess/hedft_loader.load_smoothed_speed_reference
   ↓
THIS MODULE: objective → Nelder-Mead fit → provenance-stamped bundle
   ↑
driven by: scripts/extraction/method_b_extraction.py        (per-case, B4)
           scripts/extraction/method_b_shared_refit.py      (§9 shared)
           (see scripts/extraction/scripts_explained.md)
```

The written bundles load back through
`i2_helium_md.presets.load_drag_coefficients` — Method-B output is a
coefficient swap behind the same interchangeable surface; `drag.py` /
`baoab.py` consume it unchanged.

## Units and conventions

| symbol | unit | meaning |
|---|---|---|
| `a` | amu/ps | linear friction coefficient (`γ = g·(a + b·v²)`) |
| `b` | amu·ps/Å² | cubic-force coefficient |
| `E_bind` (`e_bind_eV`) | eV | effective droplet binding (calibrated, ≤ static 0.308) |
| objective, RMSE, penalty | Å/ps | all scores live on the speed axis |
| times, windows | ps | reference clock (explosion onset at t = 0) |

Module constants: `E_BIND_STATIC_EV = 0.308` (the physical **upper bound**
for the effective binding — the dynamical escape barrier is below the
static solvation energy, never above) and `_NONFINITE_SENTINEL_APS = 1e3`
(a large but finite objective for a non-finite trajectory, so Nelder-Mead
recovers instead of crashing).

## The objective

One evaluation at `θ = (a, b, E_bind)`:

```
objective(θ) = RMSE_window( mean_ensemble |v2|^MD(θ) , |v2|^smoothed-ref )
             + penalty_weight · (1 − escape_fraction)        [Å/ps]
```

- **From-onset:** each evaluation is a full ion-stage forward integration
  from the explosion onset, scored in-window only (via
  `compare_speed_to_reference`'s `window=`). The pre-t\* transient is part
  of what the coefficients must absorb — a locked design decision
  (production usage integrates from onset too).
- **Escape penalty:** `penalty_weight` defaults to 2.0 Å/ps — an order of
  magnitude above a good in-window RMSE (~0.1–0.4), so the trapping-basin
  boundary is unambiguous even where the windowed RMSE is locally flat. An
  accepted fit must have `penalty == 0` (full escape); the writers refuse
  anything else.
- **`_escape_fraction(ion, gate_steepness)`:** fraction of all 2N ions
  whose final depth `r − R_droplet` exceeds `+2·gate_steepness` (well past
  the erf gate, drag off) with non-negative radial velocity. A trapped
  ensemble (R(t) reversing inside the well) scores ~0.
- **Exact §6.5.1 pairing per evaluation:** each trial stamps the candidate
  `E_bind` into the coefficient bundle **and** sets
  `cfg.binding_energy_I_ion_eV` to the same value, so the config-load
  guard passes strictly — no escape hatch anywhere in the fit loop.
- **Deterministic:** the drag ion path is RNG-free (BAOAB at `T_eff = 0`,
  no collisions, no mass attachment), so the neutral checkpoint is computed
  **once** per (case, seed) and reused; the objective is bit-deterministic
  in `θ`. The driver scripts' DRY_RUN asserts exactly this.
- **VMI is never in the objective** — it stays the held-out downstream
  arbiter (METHOD_B §4).

## Data structures

| class | role |
|---|---|
| `ObjectiveResult` (frozen) | one evaluation: `(a, b, e_bind_eV)`, `objective`, `rmse_Aps`, `escape_fraction`, `penalty_Aps`, `num_overlap_points` |
| `CaseSetup` | per-case fit context: base `SimConfig`, the **cached neutral checkpoint**, reference arrays + window, normalization anchors `a0/b0`, reference path |
| `TrajectoryMatchingFit` | a completed per-case fit: `best`, per-start minima, eval count, convergence, provenance fields; uncertainty fields filled by the driver |
| `JointObjectiveResult` (frozen) | one §9 joint evaluation: shared `(a, b)`, per-case `e_bind_eV` dict, equal-weight-mean `objective`, per-case `ObjectiveResult`s |
| `SharedTrajectoryMatchingFit` | a completed §9 Stage-2 joint fit for one variant: `best`, starts, per-case windows/references, `anchors`, uncertainty fields |

## Per-case API (METHOD_B §6 / §8)

### `build_case_setup(case, *, num_molecules=50, seed=20260604, ion_time_ps=20.0, dt_ion_ps=0.01, reference_filename="cleaned_data_long.csv")`

Builds the per-case context and **runs the neutral stage once** (drag
coefficients and binding only enter the ion stage, so the checkpoint is
reusable across all evaluations). Defaults are the committed evidence
points: N=50 (Tier-0: N=50 vs N=2000 in-window RMSE identical to
0.001 Å/ps), seed 20260604 (the Tier-0 run-generation seed), 20 ps at
0.01 ps (covers the `[t*, ~14 ps]` window with escape margin — raises if
the reference window outruns the ion duration).

> **Post-re-wiring caveat (2026-06-12).** `CaseSetup.a0/b0` are read from
> the case preset's `drag_coefficients` — at delivery time that was the
> per-case **Method-A** bundle, which is what made them usable as
> optimizer normalization anchors. Since the preset re-wiring (METHOD_B
> §10.1) the presets load the shared `shared_pure_cubic` bundle, so a
> fresh `build_case_setup` now yields `a0 = 0.0` — unusable as an anchor
> (`fit_trajectory_matching`'s bounds collapse; the shared fit's anchor
> check raises). The recorded extraction runs are complete and are not to
> be re-run (first-runs rule). **Resolution (user decision 2026-06-12,
> METHOD_B §10.4):** future fits source their anchors as **pre-registered
> constants** — values read off the Method-A bundles once, locked with
> provenance in the driver's USER SETTINGS at the pre-run session — never
> via `setup.a0/b0` or a runtime bundle read.

### `evaluate_objective(theta, setup, *, penalty_weight_Aps=2.0, run_fn=run_ion_propagation)`

Evaluates the objective at `θ = (a, b, e_bind)`. Builds a config differing
from `setup.cfg_base` **only** in the ion-stage fields (candidate
coefficients + binding, stamped as an exact §6.5.1 pair), forward-integrates
from the cached neutral checkpoint, scores. Returns an `ObjectiveResult`;
a non-finite trajectory/score yields the finite sentinel objective (never
raises into the optimizer for that reason). `run_fn` is dependency-injected
so unit tests stub the MD entirely.

### `fit_trajectory_matching(setup, *, ...)`

The per-case joint 3-parameter fit. Mechanics:

- **Normalized coordinates** `x = (a/a0, b/b0, E/0.308)` — the
  Method-A-era anchors condition the optimizer.
- **Bounds:** `a ∈ [0.1, 10]·a0` (keeps the pre-§9.5 dissipativity `a > 0`
  satisfied; this is the "old `0.1·a0` lower bound" the §8 18 Å fit pinned
  against), `b ∈ [1e-3, 10]·b0`, `E ∈ [e_bind_min, e_bind_max]` with the
  static 0.308 eV as the default upper bound.
- **Initialization:** a coarse 1-D `E_bind` pre-scan at fixed `(a0, b0)`
  (default 8 points in [0.05, 0.305] eV) picks the best start; with
  `extra_starts` two more starts probe the documented drag↔binding
  degeneracy ridge (`(a0, b0, 0.9·E_max)` and `(1.3·a0, 0.7·b0, E_prescan)`).
- **Optimizer:** bounded Nelder-Mead per start (`fatol` 1e-3 Å/ps,
  `xatol` 1e-3 normalized, `maxfev` 400 per start).
- Returns a `TrajectoryMatchingFit`; `best` is the lowest-objective
  minimum across starts and `starts` records each start's minimum **so the
  degeneracy spread stays visible**. Uncertainty fields are left `None` —
  the driver script fills them.

### `sensitivity_halfwidths(setup, best, *, rel_rise=0.10, factors=(1.01…2.0), ...)`

RMSE-sensitivity half-widths for `(a, b, e_bind)` around the optimum,
via the shared `_objective_rise_halfwidths` core: for each parameter,
scan multiplicative factors in both directions (others held at the
optimum) and record the change at which the objective first rises 10%
above the best value; the half-width is the mean of the up/down crossing
distances (one side if only one crosses; the largest scanned change if
neither does — conservative and visible in review). A parameter exactly 0
(a fixed pure-cubic `a`) gets half-width 0.0 by convention — it was not a
free parameter. **This is a sensitivity band on a deterministic objective,
NOT a statistical confidence interval** (so labeled in the bundle's
`uncertainty_model`).

### `write_fit_parameters(fit, out_dir, *, transverse_contaminated)`

Writes the Method-B `fit_parameters.json` (METHOD_B §6 schema) — a
**superset of the Method-A file**, loadable by `load_drag_coefficients`
(which requires the Method-B provenance keys when
`extraction_method = "trajectory_matching"`). Carries: coefficients +
errors + `uncertainty_model`, `meff_amu`, the stamped
`effective_binding_energy_I_ion_eV` (+ error), the window
(`t_start`/`t_end`), `reference_file` (repo-relative), the objective tag
and score, seed-sweep records, the full optimizer record (per-start
minima), date/branch, and the standing flags
(`radial_projection_convention`,
`transverse_contaminated_non_radial_reference` — True for the genuinely
non-radial 9 Å case, False for clean-radial 18 Å —
`full_window_heldout_window_axis_forfeited`).

**Refusal paths (fail loudly, rule 4):** a trapped fit (`penalty != 0`)
is never written; unfilled uncertainty fields are refused. Method-A
artifacts are never touched (the bundle goes to the sibling
`trajectory_matching/` directory).

## Shared-form joint refit API (METHOD_B §9)

Variant tags: `SHARED_3PARAM` (fully shared `{a, b, E_bind}`, `a` free with
lower bound **0** — the §9.5 guard relaxation makes `a = 0` legal while
`b > 0`), `SHARED_PURE_CUBIC` (`a ≡ 0` fixed, fit `{b, E_bind}`), and
`DIAGNOSTIC_4PARAM` (shared `{a, b}`, per-case `E_bind` — run only on a
primary band failure to localize drag-law vs. binding generalization
failure; **never a production bundle**).

### `evaluate_joint_objective(a, b, e_bind_by_case_eV, setups, *, ...)`

Per-case `(RMSE + escape penalty)` **first**, then the **equal-weight mean**
across cases — the windows differ in length, so averaging the per-case
scalars prevents long-window dominance (§9.2). Case keys of the binding
dict and the setups must match exactly (raises otherwise). Equal binding
values realize the shared-`E_bind` variants; distinct values the
4-parameter diagnostic.

### `fit_shared_trajectory_matching(setups, *, variant, anchors, ...)`

The Stage-2 joint Nelder-Mead fit. Validation up front: at least 2 cases
(it is a cross-case fit by definition); all cases must share `seed`,
`num_molecules`, and `m_eff_amu` (one joint provenance stamp); anchors
must be positive. **`anchors` is an explicit keyword on purpose:** per
§9.3 the normalization anchors must come from the **18 Å Method-A bundle
only** — the 9 Å Method-A artifact was provenance-broken at the time the
rule was locked (hand-edited `a`; restored to the authentic values
2026-06-12, migration log) and was never allowed to condition the fit —
which is why the per-case `setup.a0/b0` are deliberately not used here. Parameter vectors per variant (normalized by
`(a0, b0, 0.308 eV)`):

| variant | vector | bounds |
|---|---|---|
| `shared_3param` | `(a, b, E)` | `a ∈ [0, 10]·a0` (lower bound 0), `b ∈ [1e-3, 10]·b0` |
| `shared_pure_cubic` | `(b, E)` | `a = 0` fixed |
| `diagnostic_4param` | `(a, b, E_case1…E_caseK)` | per-case bindings, setups insertion order |

Same pre-scan + multi-start scheme as the per-case fit (shared-`E_bind`
pre-scan at the anchor coefficients — `(0, b0)` under pure-cubic — plus
the high-binding and coefficient-tilted ridge-probing starts). Returns a
`SharedTrajectoryMatchingFit`; `n_evaluations` counts joint evaluations
(each = one ion run **per case**).

### `joint_sensitivity_halfwidths(setups, best, *, ...)`

Same band semantics as the per-case scan, on the joint objective, for
`(a, b, E_shared)`. Defined **only** for shared-`E_bind` variants — the
4-parameter diagnostic has no single shared `E` to scan and never produces
a production bundle, so it is refused. Under pure-cubic (`best.a == 0`)
the `a` half-width is 0.0 by convention.

### `write_shared_fit_parameters(fit, out_dir, *, stage="stage2_joint")`

Writes the shared-form bundle (METHOD_B §9.6) — a superset of the per-case
Method-B schema that **loads through `load_drag_coefficients` unchanged**.
§9.6 additions: `calibration_cases`, `stage`, `variant`, the per-case
`windows` / scores / `reference_file` list, and the anchor provenance
(including the "18A Method-A bundle only" note). The loader-required
`t_start`/`t_end` span the **union** of the per-case windows. The flags
block adds `cross_case_axis_consumed_by_joint_fit` and sets the
transverse-contamination flag True (9 Å is *inside* the Stage-2 joint
objective — the §3.5 accepted risk, carried by the bundle).

**Refusal paths:** the `diagnostic_4param` variant (per-case bindings are
a failure-localization device, never a loadable production pairing); a
trapped fit in **any** case; unfilled uncertainty fields.

## Form-discrimination API (METHOD_B §10)

The §10 layer runs the `linear_quadratic` and `power_law` families through
the **same** joint objective against the pure-cubic incumbent. It mirrors the
§9 API form-for-form, with the coefficients carried as the form's raw
closed-form **dict** (`{a, c}` / `{C, n}`) instead of the `linear_cubic`
scalars. Variant tags: `LQ_SHARED_3PARAM`, `LQ_SHARED_PURE_QUADRATIC`,
`PL_SHARED_3PARAM` (collected in `FORM_VARIANTS`). Data structures:
`FormObjectiveResult` / `JointFormObjectiveResult` (form-generic analogs of
the §9 results) and `SharedFormTrajectoryMatchingFit`.

- **`evaluate_form_objective(form, coefficients, e_bind_eV, setup, *, ...)`**
  and **`evaluate_joint_form_objective(form, coefficients, e_bind_by_case_eV,
  setups, *, ...)`** — the form-generic front-ends of `_run_and_score` (the
  core extracted from `evaluate_objective`, which now delegates to it
  bitwise-equally). Same objective semantics; `evaluate_form_objective`
  accepts `linear_cubic` too, for the cross-check against `evaluate_objective`.
- **`fit_shared_form_trajectory_matching(setups, *, variant, anchors,
  n_bounds=(1,4), ...)`** — the Stage-2 shared joint fit (≥2 cases), same
  normalized bounded Nelder-Mead + `E_bind` pre-scan + multi-start discipline
  as §9. The normalized vector is family-specific: lq `(a/a0, c/c0, E/0.308)`
  with `a` lower bound 0 (pure-quadratic reachable); pl
  `(γ_ref/γ_ref0, n, E/0.308)` in the **pivot** parameterization
  `γ_ref = C·v_ref^(n−1)` (§10.4.1 — axis-aligns the matching ridge so the
  `n̂` half-width is meaningful), with `n` direct-bounded by `n_bounds`. The
  `anchors` dict is the family's pre-registered constants (`{a0, c0}` /
  `{C0, n0, v_ref_Aps}`), validated positive — **never** a runtime bundle
  read.
- **`fit_form_trajectory_matching(setup, *, variant, anchors, **kwargs)`** —
  the single-case Stage-1 analog (fit on 18 Å only, predict 9 Å). Same engine;
  record-only output (the writer refuses single-case bundles).
- **`form_sensitivity_halfwidths(setups, best, *, v_ref_Aps=None, ...)`** —
  RMSE-sensitivity half-widths; lq scans `(a, c, E)` (fixed `a=0` → 0.0 by
  convention), pl scans the **pivot** `(γ_ref, n, E)` (requires `v_ref_Aps`;
  the driver converts the `γ_ref` half-width to `C_err = γ_ref_err/v_ref^(n̂−1)`
  at fixed `n` and stamps `n_err` directly).
- **`write_shared_form_fit_parameters(fit, out_dir, *, anchor_provenance,
  stage="stage2_joint")`** — the §10 bundle writer; same loader contract
  (raw coefficients + `<k>_err` + the §9.6 provenance superset + the variant
  tag + a new `model_selection_on_seen_data` flag). power_law bundles add a
  `pivot` block (`v_ref_Aps`, `gamma_ref`) for re-derivability. Refuses
  single-case (Stage-1) fits, trapped fits, and unfilled uncertainty.

Drivers: `scripts/extraction/form_phase_common.py` +
`method_b_form_refit_{linear_quadratic,power_law}.py`.

## Provenance helpers

`_reference_file_str` (repo-relative POSIX path for the stamp; absolute if
outside the repo, e.g. a test tmp dir) and `_git_branch` (current branch
or `"unknown"` — the provenance stamp must never crash a fit run).

## Status and test coverage

All production runs are complete and recorded (per-case §8: delivered but
failed the provisional cross-case bands, bundles flagged not-yet-usable;
shared §9.7: **PASS**, pure-cubic equivalent, Tier 1 ungated; §10
form-discrimination §10.7: **incumbent CONFIRMED** — the free-`n` `power_law`
fit recovers `n̂ = 2.927 ≈ 3` (EQUIVALENT, Δ = −0.0001) and forced-`v²`
`linear_quadratic` is measurably worse (Δ = +0.0339, rejected), so no
alternative form beats `shared_pure_cubic`). Per the first-runs rule the
recorded verdicts stand — do not re-run the fits to "check" them.

Tests: `tests/test_extraction_trajectory_matching.py` (60 tests; stub-MD
via the `run_fn` injection seam — fit recovery on stub laws incl. the §10
form variants and the pl pivot transform, joint objective semantics, the
bitwise `linear_cubic` cross-check of `evaluate_form_objective` vs
`evaluate_objective`, variant behavior, all refusal paths, writer/loader
round-trips, committed shared-bundle artifact checks; plus one real
N=2/0.2 ps wiring evaluation). Full suite **705 passed / 0 failed**.
