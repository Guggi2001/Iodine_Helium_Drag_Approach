# Tier-0 Drag-Comparison Scripts — Guide

This document describes the **three remaining Tier-0 scripts** and *what each
aims to do* at the current state of the drag-model port. Earlier helpers (the
t\*-seeded and same-smoothed comparison scripts) have been removed; their roles
are folded into the surviving three plus the package comparison APIs.

Read alongside:
- `METHOD_B_trajectory_matching_extraction.md` — the active extraction method
  (trajectory-matching) and the same-smoothed objective these scripts score.
- `TIER0_FINDINGS.md` — the settled verdict (18 Å clean pass, 9 Å non-radial)
  and the held-out-generalization role Tier 0 now plays.
- `DRAG_PORT_DESIGN_DECISIONS.md` §6 — the validation hierarchy.

## Current Tier-0 role

The production drag law is `shared_pure_cubic` (`γ = g·b·v²`), extracted by
**Method B** (trajectory matching) jointly with the effective droplet binding.
Under Method B the in-window trajectory match *is* the fit objective, so Tier 0's
old "consistency check" framing is retired: these scripts now exist to (a)
**re-run any catalog drag law / hand-tune** through the real pipeline and (b)
**score** the result against the very quantity the coefficients were fit to —
the CEEMDAN+SG same-smoothed `|v2|` reference (`cleaned_data_long.csv`).

The three scripts share one selection surface: a `CASE` (`9A`/`18A`) plus a
`VARIANT` (a key of `tier0_common.CATALOG`). The generator and the comparison
agree on the run-directory name through `tier0_common`, so the comparison always
finds exactly what the generator wrote.

```
tier0_common.py          -> shared catalog + cfg builder + naming (imported by the other two)
gen_tier0_runs.py        -> build + run a drag run for one CASE/VARIANT (data/runs, gitignored)
tier0_drag_comparison.py -> score that run vs the same-smoothed |v2| reference; plot; export means
```

---

## 1. `scripts/tier0_common.py`

**Aim.** Be the single source of truth for *which* Method-B drag bundle a Tier-0
run uses and *how* it is wired, so the generator and the scorer can never
disagree. It is a module (imported by the other two scripts), not a runnable
entry point.

**What it provides.**
- **`CATALOG`** — named map from a variant key to its on-disk bundle directory
  under `data/reference/drag`. Two families: **shared** joint-refit bundles
  (`shared_pure_cubic`, `shared_3param`, `shared_pl`, `shared_lq`,
  `shared_lq_pq`) that run on either geometry, and **per-case** single-curve
  bundles (`percase_linear_cubic`, `percase_pl`, `percase_lq`, `percase_lq_pq`)
  whose `{case}` placeholder resolves to the run's own geometry — a **geometry
  guard** that makes cross-geometry application impossible by construction.
  Method-A bundles are intentionally excluded (no jointly-validated binding).
- **`resolve_bundle_dir(case, variant)`** — catalog key → absolute bundle dir,
  with case/variant validation.
- **`build_drag_cfg(...)`** — load a bundle, optionally apply an inline
  hand-tuning override (`coeff_overrides`, `e_bind_override`), and wire the
  resulting `{form, coefficients, E_bind}` onto the case's **base** preset under
  the Tier-0 envelope (`mass_scenario="fixed"` at `m_eff`, no escape hatch). With
  no overrides the cfg is identical to the committed-bundle run. It wires onto
  the base presets (not the frozen production drag presets) precisely so the
  catalog stays interchangeable. The §6.5.1 binding-pairing identity holds by
  construction, so `validate()` passes with no hatch.
- **`run_dir_name(case, variant, n, run_tag)`** — the agreed `data/runs`
  basename; `run_tag` marks a hand-tuned run so it never overwrites the
  committed-variant run.
- **`window_source_dir(case)`** — the scoring window is a **case** property, not
  a variant property (a shared bundle stamps the 9 Å onset `t_start=2.67` for
  both cases), so the window is always read from the per-case
  `trajectory_matching` bundle (9 Å `[2.67, 14.081]`, 18 Å `[4.54, 14.768]`).

Covered by `tests/test_tier0_common.py`.

---

## 2. `scripts/gen_tier0_runs.py`

**Aim.** Produce one full-duration drag run for a chosen `CASE` + `VARIANT` (with
optional inline hand-tuning), to be scored by script 3.

**What it does.**
- Builds the cfg via `tier0_common.build_drag_cfg`, then runs the real pipeline
  (`run_neutral_propagation` → `run_ion_propagation`) in the Tier-0 envelope
  (full duration, `dt_ion=0.01`, fixed seed, noise off, fixed `m_eff`).
- **Catalog-driven, not preset-frozen.** The drag law is chosen by `VARIANT`
  independent of the production presets; with the default `shared_pure_cubic` and
  no overrides the run matches `single_pulse_N2000{,_18Angst}_drag`.
- **Optional hand-tuning.** `COEFF_OVERRIDES` (keys must be coefficients of the
  chosen form) and `E_BIND_OVERRIDE` allow a Method-B-style manual probe; setting
  either flags the run as tuned. `RUN_TAG` keeps a tuned run from overwriting the
  committed-variant run dir.
- **Reduced N to iterate, full N to confirm transfer.** `compare_*` reduce to
  means over molecules, so a small N sets the number and a full-N run confirms it
  transfers (18 Å N=50 vs N=2000 agree to 3 decimals).
- Output lands in `data/runs/<run_dir_name>/`, which is **gitignored** —
  checkpoints are large and the runs are deliberately ephemeral. The committed
  artifact is the tiny mean-series CSV exported by script 3.

**Run.** Edit the `# USER SETTINGS` block (`CASE`, `VARIANT`, `N`,
`ION_TIME_PS`, overrides, `RUN_TAG`), then `python scripts/gen_tier0_runs.py`.

---

## 3. `scripts/post_processing/tier0_drag_comparison.py`

**Aim.** Score a finished drag run against the reference, report the numbers,
plot, and export the (free-to-churn) ensemble-mean trajectory.

**What it does.**
- Loads the run (`RunDirectory`), the HeDFT reference (`load_hedft_trajectory`),
  the same-smoothed reference (`load_smoothed_speed_reference` on
  `cleaned_data_long.csv`), and the window + `m_eff` from the **per-case** bundle
  via `window_source_dir` — never from the (possibly shared) variant being run.
- **Time-origin safeguard.** Asserts MD ion-stage `t=0` (charge-switch)
  coincides with the HeDFT `t=0` and that the window lies inside both overlaps,
  so a silent offset can never score the wrong segment.
- **GATE = same-smoothed `|v2|` RMSE.** The scored gate is the in-window `|v2|`
  RMSE against `cleaned_data_long.csv` (column `cleaned_SG`) — *the exact
  quantity Method B minimised* (METHOD_B §6), so the comparison scores what the
  coefficients were fit to. The loader accepts both the 18 Å 2-column and the
  9 Å 3-column (`+IMF_cleaned`) layouts.
- **Diagnostics (reported, not gated):** raw-HeDFT distance RMSE, per-atom
  (I1/I2) `|v|` RMSE and their mean, the smoothed `mean_ratio`, the raw
  `mean_ratio`s, and the I1–I2 split (symmetry sanity). There is no smoothed
  distance or smoothed I1 reference, so the smoothed gate is intrinsically
  I2-velocity-only.
- **Plots full trajectory, scores in-window.** Distance + velocity panels show
  the *full* trajectory with the scored window shaded (the uncalibrated pre-`t*`
  transient and post-window divergence are visible but excluded). Optional
  per-atom kinetic-energy figure (`ENERGY_FIGURE`, with the binding-depth line)
  and radial-projected force-balance figure (`FORCE_FIGURE`), both reconstructed
  with the same physics functions the ion driver uses (no duplicated physics).
- **Exports the mean-series CSV** (`md_mean_trajectory.csv`): the four
  ensemble-mean series that losslessly reproduce every scored + diagnostic
  number. This export is a free-to-churn diagnostic, distinct from the committed
  `md_mean_trajectory_N50.csv` regression reference.

**Run.** Edit the `# USER SETTINGS` block (`CASE`, `VARIANT`, `N`, `RUN_TAG`,
figure toggles) to match the generated run, then
`python scripts/post_processing/tier0_drag_comparison.py`.
