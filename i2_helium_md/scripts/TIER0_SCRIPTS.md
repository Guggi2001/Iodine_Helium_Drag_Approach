# Tier-0 Drag-Comparison Scripts — Guide

This document describes the scripts and committed artifacts built for the
**Tier-0 TDDFT drag comparison** (the first validation tier of the drag-model
port). It explains *what each script does* and *with what intent*, so the
workflow can be re-run and the design choices are auditable.

Read alongside:
- `TIER0_COMPARISON_spec.md` — the task specification.
- `TIER0_FINDINGS.md` — the first-run verdict (18 Å passes, 9 Å is the
  documented different-regime case) and the diagnosis.
- `DRAG_PORT_DESIGN_DECISIONS.md` §6 — the validation hierarchy these scripts
  serve.

## What Tier-0 asks

Does the deterministic, fixed-`m_eff`, `linear_cubic` drag model reproduce the
TDDFT reference trajectory **inside the extraction window** `[t*, t_end]`? The
deliverable is a judgment plus a *manufactured threshold* committed as a
regression gate, not an objective pass/fail.

Every scored number is an **in-window RMSE** (distance and per-atom speed),
computed with the new `window=` parameter on
`i2_helium_md/postprocess/compare_trajectories.py`. `window=None` is the
pre-existing whole-overlap behaviour, bit-identical; `window=(t*, t_end)`
restricts scoring to the calibrated window. The window bounds and `m_eff` are
read from each case's `data/reference/drag/<case>/linear_and_cubic/
fit_parameters.json` — the *same provenance source* as the drag coefficients —
so the scored region can never drift from the extraction's own `[t*, t_end]`.

## The end-to-end workflow

```
1. gen_tier0_runs.py            -> produce a drag run (ephemeral, data/runs gitignored)
2. tier0_drag_comparison.py     -> score it from-onset; export the mean-series reference
   tier0_tstar_seeded_comparison.py -> the clean (transient-free) form test
3. tests/test_tier0_drag_comparison.py -> the committed regression gate
```

---

## 1. `scripts/gen_tier0_runs.py`

**Purpose.** Generate the full-duration drag runs that feed the comparison.

**Intent / design.**
- Runs the real neutral + ion pipeline (`run_neutral_propagation` →
  `run_ion_propagation`) for a chosen case (`9A`/`18A`) and molecule count `N`,
  in the strict **Tier-0 envelope**: full 20 ps, `dt_ion=0.01`, fixed seed,
  noise off, `mass_scenario=fixed` at `m_eff` (inherited from the
  `single_pulse_N2000{,_18Angst}_drag` presets).
- **Reduced N for the threshold-setting iteration, full N to confirm transfer.**
  `compare_*` are means over molecules of a single-trajectory-shaped quantity,
  so a small N suffices to set the threshold; a full-N run then confirms the
  number transfers (it does — 18 Å N=50 vs N=2000 are identical to 3 decimals).
- Outputs land in `data/runs/<case>_drag_tier0_N<N>/`, which is **gitignored**
  (checkpoints are large); the runs are deliberately ephemeral working
  artifacts. The committed regression reference is the tiny mean-series CSV
  (below), not the checkpoint.

**Run.**
```
python scripts/gen_tier0_runs.py 9A 50
python scripts/gen_tier0_runs.py 18A 50
python scripts/gen_tier0_runs.py 18A 2000     # transfer confirmation
```

---

## 2. `scripts/post_processing/tier0_drag_comparison.py`

**Purpose.** The **from-onset** comparison harness: score a finished drag run
against its TDDFT reference, report the Tier-0 numbers, plot, and export the
committed mean-series reference.

**Intent / design.**
- Loads the run (`RunDirectory.load_ion`), the HeDFT reference
  (`load_hedft_trajectory`), and the window from `fit_parameters.json`.
- **Time-origin safeguard (§3.2).** Asserts the MD ion-stage `t=0`
  (charge-switch) coincides with the HeDFT `t=0`, and that the window lies inside
  both overlaps. A silent half-ps offset would score the wrong segment; this is a
  hard precondition, cheap because the origins *do* align (both start at 0).
- **Scores in-window only:** distance RMSE (Å) and the mean of I1/I2
  velocity-magnitude RMSE (Å/ps) — the two Tier-0 gates. Reports `mean_ratio`
  and the I1–I2 split as *diagnostics, not gates* (symmetry sanity, §4).
- **"Plotted but not scored" (§3, D2).** The figure shows the *full* trajectory
  with the scored window shaded, so the uncalibrated pre-`t*` transient and the
  post-window divergence are visible but excluded from the scored scalar.
- **Exports the mean-series reference.** The four ensemble-mean series
  (`time_ps, mean_distance_A, mean_speed_I1_Aps, mean_speed_I2_Aps`) fully
  reproduce every scored + diagnostic number (the comparisons reduce to means
  over molecules), so this small CSV is the committed regression reference — no
  large checkpoint committed.

**Caveat.** This mode runs from the explosion onset, so it conflates the
*uncalibrated pre-`t*` transient* with the drag-form question. Use the t*-seeded
script (below) to isolate the form.

**Run.** Edit the `# USER SETTINGS` block (`RUN_DIR`, `HEDFT_PATH`,
`DRAG_COEFF_DIR`, `EXPORT_MEAN_SERIES_PATH`, `SHOW_FIGURE`), then
`python scripts/post_processing/tier0_drag_comparison.py`.

---

## 3. `scripts/post_processing/tier0_tstar_seeded_comparison.py`

**Purpose.** The **t\*-seeded clean-form test**: isolate the drag *form* from the
transient by seeding the ion at `t*` with the reference state and integrating the
drag law forward over `[t*, t_end]`.

**Intent / design.**
- Removes the from-onset transient contribution: if the integrated law tracks the
  reference *from a point on it*, the form is adequate and the from-onset
  divergence was the transient + off-trajectory compounding. This is the
  strictly-cleaner form isolation hinted at in `DRAG_PORT_DESIGN_DECISIONS` §6.4.
- Builds a 1-molecule synthetic `NeutralCheckpoint` at the reference state
  (same seeding pattern as the Slice-4 smoke harness) and reuses
  `run_ion_propagation`, so **no drag physics is re-implemented**. The ion clock
  is then shifted by `+t*` to align with the reference.
- **Seed convention — `|v|` placed radially (internal-consistency).** The
  reference CSVs now carry the **full 3D per-atom velocity *and* per-atom
  positions** (`V{1,2}_{x,y,z}`, `X/Y/Z{1,2}`), so the true I–I axis and the real
  radial/transverse split are recoverable. The seed nonetheless places each
  atom's *full speed* `|v_i|` along the constructed I–I axis (transverse zero),
  because the drag law was extracted with the scalar speed `v = |v2|` as the
  radial velocity — forward-integrating γ from `|v|` is the honest Tier-0
  internal-consistency check (does the BAOAB driver reproduce the speed γ was fit
  to?). This **retires the earlier `dR/dt`-reconstruction** (central difference on
  raw `R(t)` for the radial rate, leftover speed dumped into a fictitious
  transverse component); the finite difference is gone.
- **Real-split diagnostic + provenance cross-check.** `_report_real_radial_split`
  builds the true bond axis `R̂ = (r1−r2)/|r1−r2|` from the position columns at
  `t*` and **prints** each atom's `v·R̂` (radial) vs the transverse remainder —
  this *reports*, it does not feed the seed. It also asserts
  `|r1−r2| == R_distance` at `t*`. Result: **9 Å atom 2 is ~99.6% transverse**
  (radial −0.29, |v2|=4.90 — it genuinely co-translates); **18 Å atom 2 is
  ~99.9% radial** (radial −3.12).
- Reads the droplet radius from a from-onset run of the same case (so the spatial
  gate matches the real pipeline).

**Run.** Edit `CASE` / `ONSET_RUN_DIR` in `# USER SETTINGS`, then
`python scripts/post_processing/tier0_tstar_seeded_comparison.py`.

**Why two comparison modes.** From-onset = the production-representative path
(what the real simulation does). t\*-seeded = the diagnostic that separates
"transient error" from "form error". Comparing the two was decisive: at 9 Å the
t\*-seed cut the distance RMSE from 9.85 → 2.45 Å (so the transient is a large
contributor) while a velocity residual persisted. With the real 3D + position
data, that residual is now diagnosed by the real-split print above: **9 Å atom 2
is ~99.6% transverse at `t*` — it genuinely co-translates**, which a central-force
MD (Coulomb + radial droplet + radial drag) structurally cannot reproduce. The
residual is therefore a *different-regime* signal, **not** a harness artifact and
**not** a drag-law error (see `TIER0_FINDINGS.md` → "t\*-seed from real 3D
velocities").

---

## 4. `tests/test_tier0_drag_comparison.py` — the committed regression gate

**Purpose.** Pin the Tier-0 thresholds so later tiers cannot silently degrade the
in-window form match.

**Intent / design.**
- **Anchored on 18 Å**, which passes Tier-0 cleanly; **9 Å** is asserted
  *finite-only* (the documented different-regime case, `TIER0_FINDINGS.md`).
- Loads the committed mean-series CSV, rebuilds a 1-molecule `IonCheckpoint`
  whose molecule-means reproduce the stored series, reads the window from
  `fit_parameters.json`, and asserts the **windowed** RMSEs against named
  thresholds with a provenance comment:
  - distance RMSE ≤ **3.0 Å** (run value 2.512, N=50≈N=2000),
  - mean(I1,I2) |v| RMSE ≤ **0.25 Å/ps** (run value 0.163),
  - I1–I2 split ≤ 0.1 Å/ps (loose symmetry sanity, not a physics gate).
- **What it does / does not catch.** It locks the committed reference + the
  windowed-comparison machinery + the threshold numbers (the auditable home of
  the gate). Re-validating the drag *physics* after a change is done by
  regenerating the run (scripts 1–2), not in-test (no production checkpoint is
  generated in tests).

---

## Committed data artifacts

- `data/reference/drag/{9A,18A}/tier0/md_mean_trajectory_N50.csv` — the tiny
  (~200 KB) ensemble-mean references the regression gate consumes. Each carries a
  provenance header (preset, N, seed, duration, `m_eff`, window, date/branch).

## Package change that underpins all of the above

- `i2_helium_md/postprocess/compare_trajectories.py` gained the additive
  `window=(t_start, t_end)` keyword on `compare_distance` /
  `compare_velocity_magnitude` / `_compare_series`. `window=None` is
  bit-identical to the prior whole-overlap behaviour; windowing lives in one
  place (`_compare_series`) so the finite-mask / ratio-guard logic is not
  re-derived. Covered by `tests/test_compare_trajectories.py::TestWindow`.
