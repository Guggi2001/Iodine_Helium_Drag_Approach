# Tier-0 TDDFT Drag Comparison — Findings

**Status:** First-run findings from the §6.4 Tier-0 comparison (the
"read the first run" judgment, `TIER0_COMPARISON_spec.md`). Records what the
deterministic, fixed-`m_eff`, `linear_cubic` drag model reproduces against the
TDDFT references in-window, and *why the two cases diverge*. **Acceptance
thresholds remain deferred** pending the decision below (§6.10 permits this).

Date: 2026-06-04. Branch: `drag_implementation`.

## Setup

- Model: `drag_form=linear_cubic`, `mass_scenario=fixed` at `m_eff=202.954`,
  `noise=none` — the Tier-0 envelope. Coefficients/window from
  `data/reference/drag/<case>/linear_and_cubic/fit_parameters.json`
  (9 Å: a=13.86, b=2.58, window [2.67, 8.5]; 18 Å window [4.54, 8.0]).
- Runs: `single_pulse_N2000{,_18Angst}_drag`, reduced **N=50**, full 20 ps,
  `dt_ion=0.01`, scored in-window only via the new `window=` parameter on
  `compare_distance` / `compare_velocity_magnitude`.
- Two comparison modes:
  - **from-onset** — integrate from the ion-stage explosion onset (t=0), score
    `[t*, t_end]`. Conflates the uncalibrated pre-`t*` transient with the form.
    (`scripts/post_processing/tier0_drag_comparison.py`)
  - **t\*-seeded** — seed the ion at `t*` with the reference state and integrate
    forward; isolates the **form** from the transient (§6.4 strictly-cleaner).
    (`scripts/post_processing/tier0_tstar_seeded_comparison.py`)

## Results (in-window RMSE)

| case | mode | distance RMSE (Å) | mean(I1,I2) \|v\| RMSE (Å/ps) | \|v\| ratio I1/I2 |
|---|---|---|---|---|
| 9 Å  | from-onset | 9.85 | 1.55 | 0.59 / 0.72 |
| 9 Å  | t\*-seeded  | 2.45 | 1.09 | 0.70 / 0.88 |
| 18 Å | from-onset | 2.51 | **0.16** | 1.02 / 1.03 |
| 18 Å | t\*-seeded  | **0.30** | 0.16 | 0.99 / 1.00 |

## Verdict

- **18 Å: `linear_cubic` reproduces the TDDFT trace cleanly.** t\*-seeded
  distance tracks to **0.3 Å** and per-atom speed to ~0.16 Å/ps across the whole
  window (ratios ~1.0); from-onset is nearly as good (|v| RMSE 0.16). This is a
  genuine Tier-0 pass for 18 Å.
- **9 Å: poor lab-frame match, but *not* a form failure.** The residual is
  diagnosable, and `linear_quadratic` would not help (the issue is an over-large
  effective friction, not the high-`v` wing shape).

## Why the cases diverge (diagnosis)

The drag law is self-consistent with the extraction at `t*` by construction
(`drag_data.csv` `v_spline=4.90`, `F_drag=321` matches 9A_All_Data's I2 speed —
same trajectory). The divergence is in *forward integration*, and is a
**center-of-mass-drift / reference-frame effect**:

- At **18 Å** `t*` the explosion is clean-radial: per-atom `|v|` ≈ the radial
  separation velocity, COM drift ≈ 0. The extracted `γ(v)`, applied to lab-frame
  speed, *is* effectively applied to the radial speed → it reproduces the trace.
- At **9 Å** `t*` the reference per-atom velocity is dominated by a **~4 Å/ps COM
  drift** (atom 1 moves mostly in x, atom 2 mostly in y — nearly orthogonal; the
  I₂⁺ is drifting through the droplet, not just exploding). The per-atom drag
  `γ(|v_lab|)·v_lab` therefore damps the COM drift, collapsing MD speed
  (4.7→~2.0 across the window) while the reference coasts at ~4.2. The
  separation roughly tracks (t\*-seeded distance RMSE 2.45 Å) but the speed is
  over-damped.

**Physical reading:** drag should act on the ion-relative-to-helium velocity.
When the ion (or ion+bubble) co-drifts with the He, lab speed is high but
relative speed is low → little real drag. The extraction used **lab speed**
(`v_spline` ≈ lab `|v|`). For a clean radial explosion (18 Å) lab ≈ relative and
the law works; for a trajectory with large COM drift (9 Å) lab ≫ relative and
`γ(|v_lab|)` over-damps. This is the §3.6 "the cases are in different regimes"
finding, traced to an **extraction-frame** subtlety (§6.10 extraction-side item),
not the drag *form*.

The from-onset 9 Å number (9.85 Å) is additionally inflated by the
uncalibrated pre-`t*` transient (§6.7): t\*-seeding removes most of the distance
error (9.85→2.45 Å), confirming the transient is a large extra contributor at
9 Å but **not** the whole story — the velocity over-damping persists on-trajectory.

## Data limitation found (recorded)

The HeDFT CSVs store per-atom **speed magnitudes** + only the **x,z** velocity
components (no y). The full 3D per-atom velocity is therefore **not recoverable**
from the reference file, so the t\*-seed reconstructs only the rotation-invariants
the comparison scores (R, dR/dt, |v1|, |v2|). This is sufficient for the scored
quantities but means a *fully* faithful 3D seed (and a clean COM-vs-relative
decomposition of the 9 Å reference) needs richer reference data than the current
CSVs provide.

## Decisions taken (2026-06-04)

1. **Tier-0 anchored on 18 Å.** The regression gate
   (`tests/test_tier0_drag_comparison.py`) is set from the 18 Å in-window
   residual; 9 Å is recorded finite-only as the documented different-regime case.
   Thresholds (named, with provenance in the test):
   distance RMSE ≤ **3.0 Å** (run value 2.512), mean(I1,I2) |v| RMSE ≤
   **0.25 Å/ps** (run value 0.163).
2. **Transfer to full N confirmed.** 18 Å N=50 vs N=2000 in-window RMSE is
   essentially identical — distance 2.5124 vs 2.5120 Å, mean|v| 0.1629 vs
   0.1627 A/ps — so the committed reduced-N (N=50) reference is a faithful
   stand-in and the threshold transfers.
3. **9 Å COM-drift / frame question: recorded extraction-side item, not pursued
   here.** Confirm at the extraction source whether `v_spline` is lab speed vs
   ion–He relative speed, and whether the 9 Å reference's large COM drift should
   be removed (relative-velocity re-extraction, or gating drag on relative
   velocity). Likely needs richer reference data (full 3D velocities).

## Infrastructure delivered (independent of the physics outcome)

- `compare_trajectories.py`: additive `window=(t_start, t_end)` parameter
  (`window=None` bit-identical), + windowing tests. Spec §3, §9.
- `scripts/post_processing/tier0_drag_comparison.py`: from-onset harness
  (windowed scoring, origin-alignment safeguard, mean-series export, plot).
- `scripts/post_processing/tier0_tstar_seeded_comparison.py`: t\*-seeded
  clean-form test.
- `scripts/gen_tier0_runs.py`: reproducible run generation (runs are ephemeral;
  `data/runs/` is gitignored).
- `data/reference/drag/{9A,18A}/tier0/md_mean_trajectory_N50.csv`: committed
  tiny ensemble-mean references.
- `tests/test_tier0_drag_comparison.py`: the 18 Å-anchored regression gate +
  named thresholds (9 Å recorded finite-only).
