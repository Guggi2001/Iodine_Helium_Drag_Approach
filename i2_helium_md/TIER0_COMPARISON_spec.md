# Tier-0 TDDFT Comparison — Task Specification

**Status:** Task specification. Defines the comparison, the windowing, the
threshold-setting procedure, and the regression artifact. No implementation code
until `[PROCEED TO IMPLEMENTATION]`.

**What this task is.** The first scientific test of the drag port: does the
deterministic, fixed-`m_eff`, `linear_cubic` drag model reproduce the TDDFT
reference trajectory *inside the extraction window*? It is the §6.4 **Tier 0** of
the validation hierarchy, and it is a *different kind of task* than Slices 1–4:
its deliverable is a **judgment plus a manufactured threshold**, not an objective
pass/fail. §6.10 deferred the acceptance numbers precisely because they are set
*by reading the first run* — so this task (1) produces the windowed comparison,
(2) decides whether `linear_cubic` is adequate, and (3) **writes down the
threshold** that all later tiers inherit as a regression gate.

**Precondition (now satisfied).** The `SLICE4_FIX` mass-consistency fix makes the
`fixed` drag run integrate at `m_eff`, the mass the law was extracted under. A
comparison run at the wrong inertia would be confidently wrong; that hole is
closed.

---

## 1. What is compared (D-confirmed scope)

§6.4 Tier 0: `drag_form=linear_cubic`, `mass_scenario=fixed` at `m_eff`, noise
off, against the TDDFT distance and velocity-magnitude traces **inside the
extraction window only**. Instruments: the existing
`compare_distance` / `compare_velocity_magnitude` (`compare_trajectories.py`),
extended with windowing (§3). References: `9A_All_Data.csv`, `18A_All_Data.csv`.

This is **trajectory** comparison, not distribution comparison — so the metric is
the existing per-trajectory RMSE, **not** Wasserstein (D5; Wasserstein is the
histogram tiers 2/3).

---

## 2. The constant-mass caveat — read the residual correctly (D1)

This is the **first-pass** comparison (§6.6 "first pass adequate"), not the
strictly-clean form isolation. We run fixed-`m_eff` coefficients against a
reference whose true shell mass *declines* ~21→14 He across the window. §6.6:
the residual error tracks `|m(t) − m_eff|/m_eff`, near-zero mid-window and ~⅓ at
the ends.

**Therefore, record before reading any plot:** a **mid-window-good /
ends-poorer** residual is the *expected signature of the constant-mass
approximation*, **not** necessarily a `linear_cubic` form failure. Do not promote
`linear_quadratic` (§3.4) on the basis of end-of-window divergence alone — that
pattern is consistent with the known constant-mass residual. A *form* problem
would show as in-window-shape mismatch in the **mid-window** region where the
mass approximation is good. The strictly-clean version (time-resolved `m(t)`
re-extraction, §6.6) is the escalation if the mid-window residual is itself bad —
not a prerequisite for this first pass.

---

## 3. Windowing — extend `compare_trajectories.py` (C, approach i)

The existing `_compare_series` overlaps on the **full** series
(`[max(t0s), min(t_ends)]` — e.g. `[0.0, 14.08] ps` for the 9 Å smoke check) and
its RMSE is whole-overlap. That number is dominated by the pre-`t*` transient and
the post-window long-time divergence (MD ~112 Å vs HeDFT ~84 Å) — exactly the
regions D1/D2 say are uncalibrated and **must not be scored**.

**Change:** add an **additive optional** `window=(t_start, t_end)` parameter to
`compare_distance` / `compare_velocity_magnitude` / `_compare_series`. It
restricts the overlap (step 2 of the recipe) to
`[max(t_min, t_start), min(t_max, t_end)]`. Default `window=None` reproduces the
current whole-overlap behaviour **bit-identically** — the existing 26.19 Å
smoke-check number must remain exactly reproducible, and existing callers
(`plot_hedft_comparison.py`) are unaffected.

- Single-sourced scoring: windowing lives *inside* `_compare_series`, so the
  finite-mask / ratio-guard logic is not re-derived anywhere (rule 1).
- The dataclass already exposes `t_overlap_ps` / `md_on_hedft_grid` /
  `hedft_on_overlap`; with windowing these carry only the in-window samples, and
  `overlap_t_min_ps` / `overlap_t_max_ps` report the windowed bounds.

**"Plotted but not scored" (D2):** the plot shows the **full** trajectory
`[0, t_end_full]` (call `compare_*` with `window=None` for the plot arrays, or
plot the raw series); the **threshold reads only** the windowed RMSE
(`window=(t_start, t_end)`). The pre-`t*` transient and post-window divergence are
visible but excluded from the scored scalar. §6.7: transient tolerances are loose
for all scenarios.

### 3.1 Window-bound provenance (A — confirmed in JSON)

`t_start` / `t_end` are now exported in `fit_parameters.json` alongside
`{a, b, meff_amu}`. The comparison reads the window from that file — the **same
source** as the coefficients and `m_eff` — so the scored region can never drift
from the extraction's actual `[t*, t_end]`. No hard-coded `2.67` / `8.5`. This is
the same provenance discipline as the `extraction_mass_amu` stamp: the number
that defines the calibration travels with the coefficients.

### 3.2 Time-origin alignment (B — verify, do not assume)

Before windowing means anything, assert that the MD ion-stage `t=0` (charges
switch / explosion onset) coincides with the TDDFT `t=0` (ionization instant),
and that both axes are in **ps**. A half-ps offset would silently score the wrong
segment. This is a hard precondition check in the comparison harness — a
documented assertion, not an assumption. If the origins differ, the offset must
be resolved (or applied) before any threshold is read.

---

## 4. The threshold — two scalars on the windowed RMSE (C)

Tier 0 sets **two independent thresholds**, both on the **in-window** RMSE
(they fail independently — over-strong drag shows in velocity first, in distance
cumulatively):

- **distance RMSE** (Å) — `compare_distance(..., window=(t_start, t_end)).rmse`.
- **velocity-magnitude RMSE** (Å/ps) — the **mean of I1 and I2** windowed RMSE
  (new question, confirmed): gate on
  `½(rmse_I1 + rmse_I2)`.

**`mean_ratio` and the I1–I2 split are reported diagnostics, not gates:**

- `mean_ratio` near 1 with large RMSE = right-magnitude / wrong-shape; worth
  seeing, but RMSE is the gate.
- The **I1–I2 velocity split** (`|rmse_I1 − rmse_I2|`, and the per-atom ratios) is
  reported as a **symmetry sanity diagnostic**. In a symmetric Coulomb explosion
  the two iodines should be near mirror images; a large MD-side I1–I2 asymmetry
  flags a problem *upstream of the form question* (initial conditions, a
  per-atom bug), not a drag-law verdict. The smoke check already shows a split
  (I1 ratio 0.96, I2 1.17) — interpret, don't gate.

**Threshold-setting is the judgment (§6.10).** The numbers are not chosen a
priori; they are set by inspecting the first 9 Å run's windowed residual and
deciding "this agreement is good enough to call `linear_cubic` validated at
Tier 0." Record the chosen numbers *with the run that justified them* so the
threshold is auditable, not arbitrary.

---

## 5. The run that feeds the comparison (E)

**Threshold-setting iteration: reduced-N, full-duration, full-`dt`.** The smoke
harness (N=2, ~20 steps) is far too short to span `[t*, t_end]`. Tier 0 needs a
*physically meaningful* run: the `single_pulse_N2000_drag` 9 Å preset at
**reduced N** (small molecule count) but **full ion duration** (20 ps) and full
`dt_ion = 0.01 ps`, so the trajectory traverses the window. Reduced N is
adequate because `compare_distance` is a **mean over molecules** of a
single-trajectory-shaped quantity, and the TDDFT reference is itself one
trajectory.

**Confirm against a full run before committing the threshold.** `compare_distance`
means over molecules; that mean smooths differently at N=2000 than at small N.
The confirmation step checks **the threshold transfers** from reduced-N to
full-N (not merely that the full run executes) — if the windowed RMSE shifts
materially between reduced and full N, the threshold is set on the full-N number.

---

## 6. 9 Å first, then 18 Å as a recorded-finding consistency check (D3 + D)

**9 Å first** — the canonical case; the threshold-setting and debug run.

**18 Å second — cross-case consistency (§3.6 shared-form headline).** The signal
is "`linear_cubic` fits *both* cases with only the coefficients moving" — weak
evidence the law captures real transport physics rather than per-case
curve-fitting. The 18 Å criterion is therefore **not** "matches as well as 9 Å in
absolute terms" but "`linear_cubic` at the 18 Å coefficients achieves comparable
in-window agreement, with no `linear_quadratic` promotion and no dissipativity-
guard trip."

**18 Å is a recorded finding, not a hard gate (D, §3.6).** If 18 Å diverges, that
is the §3.6 "the cases are in different regimes" result — *information to record
with its reason*, not a Tier-0 failure. Per-case form divergence is an allowed,
documented outcome (the architecture supports per-preset `drag_form`); it does
not block Tier 0. Record:

- 18 Å windowed distance + mean-velocity RMSE, alongside 9 Å's.
- whether the agreement is in the *same class* as 9 Å (consistency holds) or
  materially worse (different-regimes finding).
- if worse: this is the recorded trigger to *consider* (not auto-promote)
  `linear_quadratic` or a per-case form at a later pass — flagged, not actioned
  here.

---

## 7. The regression artifact (D4)

Tier 0's manufactured threshold becomes a **committed regression assertion** so
later tiers cannot silently degrade the form match:
`tests/test_tier0_drag_comparison.py`.

- Runs the reduced-N full-duration 9 Å drag trajectory (or loads a committed
  reference `ion.npz` for it), calls the **windowed** `compare_*`, and asserts:
  - in-window distance RMSE ≤ the Tier-0 distance threshold,
  - in-window mean(I1,I2) velocity-magnitude RMSE ≤ the Tier-0 velocity
    threshold.
- The threshold constants are named and carry a provenance comment (the run /
  date / commit that set them, per §4).
- The I1–I2 split is asserted only as a *loose* symmetry sanity bound (or
  reported, not asserted) — it is a diagnostic, not the gate.
- 18 Å: assert it *runs* and produces finite windowed RMSEs; its agreement is
  **recorded** (logged / in the test's docstring), not hard-asserted against a
  pass threshold (D / §3.6).

This converts Tier 0's one-time judgment into a durable gate — the point of the
sequential hierarchy is that each tier's winner is *locked* before the next
unknown is introduced.

---

## 8. Open / deferred items

- **Strictly-clean time-resolved re-extraction (§6.6)** — the escalation if the
  *mid-window* residual is bad (not the end-of-window residual, which is the
  expected constant-mass signature, §2). Extraction-side; not done here.
- **`temperature_diagnostic` all-NaN under the drag branch** (slice_4 Risk 3
  successor) — this task is the **first to read a drag run's diagnostics**, so
  confirm `plot_ion_temperature_diagnostic.py` (and any diagnostic the comparison
  plotting touches) tolerates an all-NaN array rather than erroring. The first
  place the post-processing surface meets the drag branch.
- **`linear_quadratic` promotion** — explicitly *not* triggered by end-of-window
  divergence (§2). Triggered only by a mid-window form mismatch or a recorded
  18 Å different-regimes finding, and even then *considered*, not auto-promoted
  (§6.8 empirical lean is recorded, the hierarchy adjudicates).
- **Acceptance threshold numbers** — set by this task's first run (§4, §6.10);
  not pre-specified.

---

## 9. Scope fence — what this task does NOT do

- **No new physics, no drag-form change.** It *evaluates* `linear_cubic`; it does
  not modify `drag.py` / `baoab.py` or promote another form.
- **No mass-dynamics, no noise.** Strictly the Tier-0 envelope
  (`fixed` / `none` / `linear_cubic`); `_check_drag_scope` already enforces it.
- **No new comparison metric.** Trajectory RMSE only; Wasserstein is tiers 2/3.
- **The only code change is the additive `window=` parameter** to
  `compare_trajectories.py` (default `None` = current behaviour, bit-identical),
  plus the new comparison harness/script and the regression test. No change to
  the existing whole-overlap callers.
- **No Tier 1/2/3 work.** Mass scenarios, the I⁺(He)ₙ size distribution, and
  ensemble second moments are later tiers.

---

## 10. Definition of done

- `compare_trajectories.py` gains an additive `window=(t_start, t_end)`
  parameter; `window=None` is bit-identical to today (smoke-check 26.19 Å
  reproducible).
- The comparison harness reads `t_start`/`t_end` from `fit_parameters.json`
  (§3.1) and asserts MD↔TDDFT time-origin alignment (§3.2).
- A reduced-N full-duration 9 Å drag run produced; windowed distance RMSE and
  mean(I1,I2) velocity RMSE computed; the I1–I2 split + `mean_ratio` reported as
  diagnostics.
- Two Tier-0 thresholds set from that run, named with provenance, and confirmed
  to transfer to a full-N run (§5).
- 18 Å run + windowed RMSEs recorded as a consistency **finding** (not a gate);
  different-regimes outcome documented if it occurs (§6).
- `tests/test_tier0_drag_comparison.py` commits the thresholds as a regression
  gate (§7).
- The §8 items recorded; the `temperature_diagnostic` NaN-tolerance of the
  plotting path confirmed.
- A short written verdict: is `linear_cubic` validated at Tier 0, and does the
  9 Å↔18 Å consistency hold or reveal different regimes?
