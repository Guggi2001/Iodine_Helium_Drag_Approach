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

## Verdict (corrected 2026-06-04 — supersedes the same-day "frame" reading)

> **Retraction of the intermediate "frame systematic" reading.** During this
> investigation the 9 Å mismatch was first attributed to a lab-vs-relative
> *velocity-frame* error (a COM-drift effect), which prompted the
> `EXTRACTION_FRAME_FIX_milestone.md`. **That reading was wrong** and is
> withdrawn — see "Why the diagnosis changed" below. It is recorded here only
> for the audit trail. The corrected diagnosis is **windowing + bubble-mode
> oscillation**, both resolvable with data in hand, and Tier 0 is reframed as an
> **internal-consistency check**.

**What Tier 0 actually is.** An internal-consistency check that the new drag
implementation + BAOAB driver **correctly forward-integrate the extracted γ**,
reproducing the (denoised) trajectory γ was fit to. It is **not** an independent
physical validation of the drag law — γ was *fit* to the reference, so a
forward-integration match is consistency, not confirmation. This is a genuinely
useful check (it exercises Slices 1–4 end-to-end on real data) and it is what
Tier 0 can honestly deliver.

- **18 Å: clean consistency pass.** t\*-seeded distance 0.30 Å, mean |v| RMSE
  0.16 Å/ps across the window (ratios ~1.0); from-onset nearly as good. The
  extraction window stays in the clean radial regime throughout, so the MD
  reproduces it well. This is the committed regression floor.
- **9 Å: consistency holds in the clean regime; the residual is explained, not a
  law failure.** The t\*-seeded MD tracks the clean extraction atom (atom 2)
  **well from t\* until ~6.5 ps**, then diverges. The full-window numbers
  (from-onset distance 9.85 Å, t\*-seeded 2.45 Å, mean |v| 1.09–1.55) are
  inflated by **three identified, non-law effects**, none of which is a
  drag-form or frame problem:
  1. the **artifact atom** (atom 1) is in the mean-based comparison (see below);
  2. the extraction window runs **past ~6 ps** where atom 2 develops a
     directional drift the central-force MD *cannot* represent;
  3. the **1.2 ps bubble oscillation** in the *raw* reference, which the
     extraction removed before fitting γ.

**Consequence:** the in-hand `{a,b}` are **not** frame-provisional. `linear_cubic`
is consistent with the extraction; no form switch is warranted; the
`EXTRACTION_FRAME_FIX_milestone.md` is **demoted to a contingency**. The
remaining Tier-0 work is the **same-smoothed consistency comparison** (below),
data in hand. Tier 1's block is lifted to *conditional* — pending that
same-smoothed confirmation (see `tier0_comparison_tasks_left.md`).

## Why the diagnosis changed (frame → windowing + bubble-mode)

Three facts overturned the frame reading:

1. **Atom 1's sideways motion is a TDDFT artifact**, judged unphysical and
   **discarded at extraction**; γ was fit to the single clean atom (atom 2), not
   an average. So the "molecular COM" `½(v₁+v₂)` the frame story rested on was
   **never in the extraction**, and one of its two inputs (v₁) is distrusted. The
   COM-drift mechanism was built on the artifact — it is unsound.
2. **The break is local and directional, not a global bias.** The MD matches
   atom 2 cleanly from t\* to ~6.5 ps (at most minimal excess drag early), then
   breaks where **atom 2 itself begins a consistent-direction drift after ~6 ps**.
   A frame/magnitude error would bias the fit *throughout*; a sharp mid-window
   onset is a **regime change in the reference**. A central-force MD (Coulomb +
   radial droplet + radial drag) has **no mechanism to produce a directional
   sideways drift**, so it cannot reproduce atom 2's late motion regardless of γ
   — the late mismatch is **not a drag-quality signal**.
3. **18 Å stays clean across its whole window** `[4.54, 8.0]`. This — not "COM
   drift ≈ 0" — is why 18 Å passed: its window never leaves the clean radial
   regime, while 9 Å's `[2.67, 8.5]` runs past atom 2's ~6 ps drift onset. Same
   form, same physics; the difference is purely whether the window stays inside
   the clean regime.

**The fit is fine — confirmed quantitatively.** Re-extracting γ on the shortened
clean window `[2.67, ~6.2]` **barely changed** the early-band velocity RMSE
(0.86 → 0.88, marginally *worse*). If late-drift points had contaminated the
*fit*, re-fitting on the clean window would have *lowered* it. It didn't — so
γ is essentially the same law either way, and the late data only inflated the
*score*, not the *fit*. The "fit-contamination" hypothesis is dead.

**What the ~0.86 residual is.** With windowing and the artifact atom accounted
for, a residual ~0.86 Å/ps remains even in the clean band — but the MD curve
matches atom 2 *by eye*. The most likely explanation is the **1.2 ps bubble
oscillation**: γ was fit to a **CEEMDAN+SG-denoised** atom-2 velocity (bubble
mode removed by construction), and the MD integrates that smooth law to produce
a **smooth** trajectory. Scoring it against the **raw** `9A_All_Data.csv` —
which still carries the large 9 Å bubble oscillation — charges the MD for not
reproducing the very oscillation the extraction *defined as noise and removed*.
RMSE penalises that oscillation-variance; the eye averages through it. **This is
the next task to confirm** (below).

## The remaining Tier-0 task — same-smoothed consistency comparison

Score the MD against the reference processed through the **same** CEEMDAN+SG
denoising the extraction used (`Drag_extraction_code.md`, identical
`target_period_ps=1.2` / `noise_width` / SG window — **not** a hand-tuned
filter), and **report it alongside the raw-reference RMSE, never instead of it**.

- If the same-smoothed RMSE drops to ~18 Å-class (~0.2) → the 0.86 was
  bubble-mode variance the law never claimed to reproduce; Tier 0 is a clean
  **internal-consistency pass**, the frame milestone is fully retired, and Tier 1
  unblocks.
- If it stays elevated against the denoised reference → a genuine residual
  survives, and *then* the demoted frame hypothesis (now: mild real drift of the
  single clean atom, not the dead molecular-COM story) or a magnitude
  recalibration revives — as a small, well-quantified residual, not the original
  "9 Å fails badly" picture.

**Caveat (kept honest):** same-smoothed scoring is *almost* circular — it checks
the MD reproduces the denoised trajectory γ was fit to, i.e. internal
consistency, not physical correctness against the true trajectory. That is
exactly what Tier 0 is for; the raw-reference number remains the harsher honest
cross-check, which is why both are always reported. Confirm the smoother
**preserves atom 2's late directional drift** (different timescale from the
1.2 ps oscillation) so it is not silently smoothed away.

## Same-smoothed comparison — RESULT (2026-06-05)

The same-smoothed comparison ran (`compare_speed_to_reference` +
`load_smoothed_speed_reference` + `scripts/post_processing/
tier0_same_smoothed_comparison.py`), scoring MD |v2| against the CEEMDAN+SG
cleaned |v2| and the raw |v2|, both modes, over the cleaned window (9 Å
`[2.67, 6.0]`, 18 Å `[4.54, 8.0]`):

| case | mode | raw \|v2\| | smoothed \|v2\| | raw−smoothed |
|---|---|---|---|---|
| 9 Å  | from-onset | 1.164 | 0.815 | +0.349 |
| 9 Å  | t\*-seeded  | 0.878 | **0.402** | +0.476 |
| 18 Å | from-onset | 0.165 | 0.104 | +0.061 |
| 18 Å | t\*-seeded  | 0.156 | **0.091** | +0.065 |

**Outcome = the middle branch.** Same-smoothed scoring roughly **halves** the
9 Å clean-form residual (0.88 → 0.40), confirming the 1.2 ps bubble-mode
oscillation as the **dominant raw contributor** — but 0.40 does **not** collapse
to the 18 Å-class ~0.09. A real residual survives in the clean regime.

**Residual characterised (per-third decomposition, t\*-seeded):**

- **9 Å: a near-uniform ~10% magnitude deficit.** MD |v2| is low across the whole
  window (ratio 0.93/0.89/0.88 early/mid/late). Of the 0.40 total RMSE, **0.394 is
  constant bias, only 0.08 is shape**. Flat, not end-loaded, not a mid-window
  shape mismatch.
- **18 Å control: zero net bias** (+0.002); residual is **pure shape** with the
  textbook constant-`m_eff` §2 signature (mid-window 0.028/ratio 1.000; ends ~0.11).
  The cleaned instrument and the machine are sound.

**The surviving 9 Å residual is therefore NOT a form failure, NOT the constant-mass
end signature, and NOT the withdrawn frame story — it is a near-constant ~10%
over-damping** (the "magnitude recalibration" contingency). Prime suspect: the 9 Å
`a` nearly doubled in the clean-window re-extraction (`a`: 13.86 → **24.876**,
`b` ≈ unchanged 2.085), ~70% stronger than 18 Å's `a`=14.556 — strong enough that
forward-integration undershoots its own fit target.

**Decision (2026-06-05): investigate the residual before unblocking Tier 1.**
The next step is an **extraction-side audit of the 9 Å clean-window `a`** (why it
nearly doubled). Coefficients are not changed without the user. Tier 1 **stays
blocked** pending that audit (not on a frame re-extraction). The committed 9 Å
`md_mean_trajectory_N50.csv` was refreshed to the current γ (the prior one
predated the re-extraction). No same-smoothed regression assertion committed yet.

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
3. **9 Å mismatch → windowing + bubble-mode, NOT a frame systematic.** The
   intermediate "frame" reading (and the `EXTRACTION_FRAME_FIX_milestone.md` it
   spawned) is **withdrawn** — see "Why the diagnosis changed." The 9 Å residual
   is explained by (i) the artifact atom in the mean comparison, (ii) the window
   running past atom 2's ~6 ps directional drift onset (a central-force MD can't
   represent it), and (iii) the 1.2 ps bubble oscillation in the raw reference
   that the extraction removed before fitting γ. The clean-window re-extraction
   barely moved the RMSE (0.86→0.88), proving the *fit* is fine. **Remaining
   task: the same-smoothed consistency comparison** (score MD against the
   reference denoised by the extraction's own CEEMDAN+SG, reported alongside the
   raw number) — data in hand, no external dependency. `EXTRACTION_FRAME_FIX` is
   **demoted to a contingency**, triggered only if a residual survives the
   same-smoothed comparison.
4. **Tier 1 block lifted to conditional.** With the frame story retired, the
   §6.3 attribution objection (COM-drift contaminating Tier 1) no longer applies.
   Tier 1 is **gated only on the same-smoothed consistency confirmation**; once
   that shows Tier 0 is a clean internal-consistency pass, Tier 1 may proceed.
5. **18 Å regression gate stands** (item 1) as the committed floor; its
   "provisional / frame-null only" caveat from the intermediate reading is
   **removed** — 18 Å is a legitimate clean-regime consistency pass, not a
   luck-of-frame artifact.

(Items 1–2 below are unchanged and remain valid.)

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
