# Tier-0 Comparison — Status: Implemented / Learned / Left To Do

**Purpose.** A single status page for the Tier-0 drag comparison after the first
run and its diagnosis. Consolidates what is built, what the investigation
concluded (and corrected along the way), and the small remaining task — so the
next session has an unambiguous picture without re-reading the whole thread.

**Read alongside:** `TIER0_FINDINGS.md` (the verdict + diagnosis),
`TIER0_COMPARISON_spec.md` (the task spec), `TIER0_SCRIPTS.md` (script guide),
`EXTRACTION_FRAME_FIX_milestone.md` (now demoted to a contingency, see below).

**One-line status:** Tier 0 is an **internal-consistency check** of the drag
implementation; 18 Å passes cleanly, 9 Å's residual is diagnosed (windowing +
bubble-mode, not a law/frame error), and the **only remaining task is the
same-smoothed consistency comparison** — data in hand, no external dependency.

---

## 1. What Tier 0 is (reframed)

An **internal-consistency check** that the Slice 1–4 drag implementation + BAOAB
driver **correctly forward-integrate the extracted γ**, reproducing the
(denoised) trajectory γ was fit to. It is **not** an independent physical
validation — γ was *fit* to the reference, so a forward-integration match is
consistency, not confirmation. The right comparison instrument is therefore the
reference **denoised by the extraction's own CEEMDAN+SG**, scored alongside the
raw reference.

---

## 2. Implemented (delivered, committed)

- **`compare_trajectories.py` `window=(t_start, t_end)` parameter** — additive on
  `compare_distance` / `compare_velocity_magnitude` / `_compare_series`;
  `window=None` bit-identical to the prior whole-overlap behaviour. Windowing
  lives in one place (no re-derived finite-mask/ratio logic). Tests:
  `tests/test_compare_trajectories.py::TestWindow`.
- **`scripts/post_processing/tier0_drag_comparison.py`** — from-onset harness:
  windowed scoring, MD↔TDDFT time-origin alignment safeguard, mean-series export,
  plot (full trajectory with the scored window shaded — "plotted but not
  scored").
- **`scripts/post_processing/tier0_tstar_seeded_comparison.py`** — t\*-seeded
  clean-form diagnostic: seeds the ion at `t*` with the reference state, reuses
  `run_ion_propagation` (no drag physics re-implemented), shifts the ion clock by
  `+t*`. The instrument that separated transient error from on-trajectory error.
- **`scripts/gen_tier0_runs.py`** — reproducible run generation; runs ephemeral
  (`data/runs/` gitignored).
- **`data/reference/drag/{9A,18A}/tier0/md_mean_trajectory_N50.csv`** — committed
  tiny ensemble-mean references (the regression gate consumes these, not large
  checkpoints), each with a provenance header.
- **`tests/test_tier0_drag_comparison.py`** — the committed regression gate,
  anchored on 18 Å: distance RMSE ≤ **3.0 Å** (run 2.512), mean(I1,I2) |v| RMSE ≤
  **0.25 Å/ps** (run 0.163), I1–I2 split ≤ 0.1 Å/ps (loose symmetry sanity).
  9 Å asserted finite-only (a NaN/crash tripwire, not a physics gate).
- **Window bounds + `m_eff` read from `fit_parameters.json`** (same provenance
  source as the coefficients); `t_start`/`t_end` exported there.

---

## 3. Learned (the investigation, including what was corrected)

### 3.1 Results in hand (in-window RMSE)

| case | mode | distance RMSE (Å) | mean(I1,I2) \|v\| RMSE (Å/ps) | \|v\| ratio I1/I2 |
|---|---|---|---|---|
| 9 Å  | from-onset            | 9.85 | 1.55 | 0.59 / 0.72 |
| 9 Å  | t\*-seeded            | 2.45 | 1.09 | 0.70 / 0.88 |
| 9 Å  | t\*-seeded, clean window (→6.2), existing γ | — | **0.86** | — |
| 9 Å  | t\*-seeded, clean window (→6.2), **re-extracted** γ | — | **0.88** | — |
| 18 Å | from-onset            | 2.51 | 0.16 | 1.02 / 1.03 |
| 18 Å | t\*-seeded            | 0.30 | 0.16 | 0.99 / 1.00 |

### 3.2 What we concluded

- **18 Å: clean internal-consistency pass.** The extraction window stays in the
  clean radial regime throughout; the MD reproduces it (0.16 Å/ps |v|).
- **9 Å: consistency holds in the clean regime; the inflated full-window numbers
  are explained by three non-law effects:**
  1. **Artifact atom in the comparison.** Atom 1's sideways motion is a TDDFT
     artifact, discarded at extraction (γ was fit to the single clean atom 2).
     But `9A_All_Data.csv` *contains both atoms*, so the mean-based comparison is
     scored partly against the artifact.
  2. **Window runs past the clean regime.** Atom 2 develops a
     consistent-direction **drift after ~6 ps**; a central-force MD (Coulomb +
     radial droplet + radial drag) has **no mechanism for a directional sideways
     drift**, so it cannot reproduce atom 2's late motion regardless of γ. The
     window's *upper* edge is unguarded (the transient-exclusion tool guards only
     `t*`). 18 Å, by contrast, stays clean across its whole window — which is the
     real reason it passed.
  3. **Bubble-mode oscillation.** γ was fit to a CEEMDAN+SG-**denoised** atom-2
     velocity (1.2 ps bubble mode removed). The MD integrates that smooth law →
     smooth trajectory. Scoring against the *raw* reference (large 9 Å bubble
     oscillation) penalises the MD for not reproducing the oscillation the
     extraction *defined as noise*. The eye averages through it; RMSE does not.
- **The fit is fine (quantitative):** clean-window re-extraction moved the |v|
  RMSE 0.86 → 0.88 (marginally *worse*), proving late-drift data did **not**
  contaminate the fit. γ is the same law on `[2.67,8.5]` or `[2.67,6.2]`.

### 3.3 What was corrected mid-investigation (audit trail)

- **Withdrawn: "9 Å is a different regime" (§3.6 first read).** Too generous; the
  residual is diagnosable, not a regime statement.
- **Withdrawn: "lab-vs-relative velocity *frame* systematic" (intermediate
  read).** This spawned `EXTRACTION_FRAME_FIX_milestone.md`. It was wrong because
  the molecular COM it rested on, `½(v₁+v₂)`, is built on the **artifact atom**
  (atom 1) — a COM was never in the extraction, and one input is distrusted. The
  real mechanism is windowing + bubble-mode (3.2), both data-in-hand. The
  milestone is **demoted to a contingency**.

---

## 4. Left to do

### 4.1 The one remaining Tier-0 task — same-smoothed consistency comparison

Score the MD against the reference processed through the **same** CEEMDAN+SG
denoising the extraction used (`Drag_extraction_code.md`: identical
`target_period_ps=1.2` / `noise_width` / SG window — **not** a hand-tuned
filter), and **report it alongside the raw-reference RMSE, never instead of it.**

- **Reuse the extraction pipeline**, do not hand-pick a smoother (a free
  smoothing knob would let RMSE be driven to anything and would launder a real
  residual).
- **Confirm the smoother preserves atom 2's late directional drift** (different
  timescale from the 1.2 ps oscillation) — so the late regime-change is not
  silently smoothed away.
- **Outcome branches:**
  - same-smoothed RMSE → ~18 Å-class (~0.2): the 0.86 was bubble-mode variance;
    Tier 0 is a clean internal-consistency pass; frame milestone fully retired;
    Tier 1 unblocks.
  - residual survives the denoised comparison: a real (small, well-quantified)
    residual; revive the *demoted* frame hypothesis (now "mild real drift of the
    single clean atom," not the dead molecular-COM story) or a magnitude
    recalibration — investigate then, not now.
- **Honest caveat to record with the result:** same-smoothed scoring is *almost*
  circular (it checks the MD reproduces the trajectory γ was fit to — internal
  consistency, not physical truth). That is what Tier 0 is for; the raw number is
  the harsher cross-check, which is why both are always reported.

### 4.2 Optional hardening (data in hand, do if cheap)

- **Principled window upper edge.** Extend `transient_exclusion.py` symmetrically:
  walk forward and end the 9 Å window where atom 2's residual *leaves* the
  stationary band (a `t_end_clean`), instead of the hard-coded ~6.2/8.5. Makes
  9 Å's window determined the same way 18 Å's is. Optional — the fit is already
  shown insensitive to the upper edge (0.86→0.88).
- **Score I2-only (clean atom) for 9 Å**, reported beside the mean, so the
  artifact-atom contribution is explicitly separated in the committed record.

### 4.3 Gated on the §4.1 outcome (not before)

- **Tier 1** — conditionally unblocked: proceed once the same-smoothed comparison
  confirms Tier 0 is a clean internal-consistency pass. The original frame-based
  block is lifted.
- **`EXTRACTION_FRAME_FIX_milestone.md`** — contingency only: pursue *iff* a
  residual survives the same-smoothed comparison *and* points to a genuine
  single-clean-atom drift effect. Its external 3D-velocity dependency is **no
  longer the project critical path**.

---

## 5. What is NOT in question (settled)

- `linear_cubic` is consistent with the extraction; **no form switch warranted**
  (`linear_quadratic`/`threshold` not promoted).
- The in-hand `{a,b}` are **not** frame-provisional; the fit is sound.
- The Slice 1–4 implementation forward-integrates the law correctly (the 18 Å
  pass and the clean-regime 9 Å match demonstrate the machine works).
- The Tier-0 infrastructure (windowing, both harnesses, the regression gate) is
  reusable as-is; nothing about it is provisional.
