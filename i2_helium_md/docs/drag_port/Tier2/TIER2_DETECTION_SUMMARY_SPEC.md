# Tier-2 Detection-Stage Summary — Design Spec

**Status: DELIVERED 2026-07-22** (trigger given same day; §3 as amended
by the implementation-survey entries below; delivery record in
`drag_migration_log_tier2.md`).
Provenance: the post-production run_summary review (log entry 2026-07-22,
"run_summary stage-mismatch review") found the legacy detector-facing
sections rendering the 30 ps handover state against detector references
on confirmation runs. User adjudication: not a relabel-only fix — the
full detector-facing figure set is rebuilt from the detection stage in a
dedicated document.

## 1. Physics rationale

The legacy `plot_run_summary.py` design encodes a hard-sphere-era
assumption: **ejection is freeze-out** — after the ion leaves the
droplet, mass is fixed and velocity is ballistic, so the end-of-window
`ion.npz` state *is* the asymptotic observable and may be compared to
detector references directly.

Tier-2 physics breaks that assumption. After the in-window MD
(~30 ps handover) the ion still evolves through the E2 relaxation stage
(Coulomb, to ~8 ns) and the µs energy-gated RRK evaporation cascade;
the asymptotic state is reached only at `t_detect` (≈ 8.53 µs,
`detection.npz`). Measured on the blessed production run
(`finc1v725`): handover n̄ = 7.77 with zero bare ions, detected
n̄ = 4.00 with 16.1 % suppressed; per-n mean KE differs by up to 26×.
Detector comparisons therefore belong to the detection ensemble; the
handover state remains a legitimate *diagnostic*, not a comparison.

## 2. Architecture

Three pieces, all under `scripts/post_processing/` (workflow rule 7 —
plot changes stay local to the scripts layer; the one package change is
§3's read extension):

1. **`summary_sections.py` (new, shared).** The detector-facing section
   builders factored out of `plot_run_summary.py` — both existing
   `detected_*` builders and the legacy-recipe builders (size
   distribution, KED mean-energy, KED curves, mass-resolved, VMI,
   polar, cov) parameterized by ensemble. Rule 1: one implementation
   of every plotting recipe; both entry scripts import from here.
   Numeric recipes (binning, movmean, peak-normalization, polar fit,
   covariance conventions) are moved verbatim — no recipe changes in
   this work.
2. **`plot_detection_summary.py` (new entry point).** `# USER
   SETTINGS` block (run dir, reference paths) like every other
   post-processing script. Writes `detection_summary.pdf` + per-section
   PNGs into the run's `figures/`. **Fails loudly** (clear error, no
   output) on a run dir without `detection.npz` — this document has no
   legacy mode.
3. **`plot_run_summary.py` (modified, detection runs only).** Stays the
   MD-window diagnostic document. On a dir *with* `detection.npz`:
   - the four detector-facing legacy sections (`mass_spectrum`,
     `ihe_ked_mean_energy`, `ihe_ked_curves_3d`, `ihe_ked_curves_2d`)
     **drop the experimental overlay** and retitle:
     *"ion-stage handover (t = t_handover) — pre-detection diagnostic;
     detector comparison: detection_summary"*;
   - the legacy VMI/polar/cov and `mass_resolved_velocities` sections
     gain the same stage banner (content unchanged);
   - the three D4 `detected_*` sections **stay** (the D4 delivery is
     not reversed), now as thin calls into `summary_sections.py`.
   On a dir *without* `detection.npz`: **byte-identical output to
   today** — the existing gating contract extends over the refactor.

## 3. Data contracts

- **Detected ensemble view (amended 2026-07-22, implementation
  survey).** The recipe sections consume a checkpoint-shaped
  `DetectedEnsembleView` (`i2_helium_md/postprocess/detected_view.py`)
  duck-typing exactly the `IonCheckpoint` surface the shared
  mass-gated diagnostics read (`mass_final_kg`,
  `velocities_final_{x,y,z}`, `b_ion_outside`, `num_molecules`),
  populated from `detection.npz` (`vx/vy/vz_detected` [Å/ps],
  `mass_detected_kg`, `n_detected`, `state_reason`). Every legacy
  recipe then runs verbatim — rule 1 with zero recipe changes. The
  originally-specced `ConfirmationDetectionRead` velocity extension is
  superseded (the read keeps its scoring role). No schema change.
  Shape validation in the view factory (rule 4).
- **n-selection (amended).** Selection runs through the existing
  single mass-gate convention (`select_final_mass_gate`) applied to
  the view; no *new* gate machinery. This is provably equivalent to
  `n_detected` selection: for detected in-band ions
  `mass_detected_kg = m(n_detected)` exactly (discrete shed), pinned
  by test.
- **Retained ions** (`state_reason = droplet_retained`): their
  `detection.npz` rows hold the verbatim in-droplet handover state,
  never a detector arrival — the view forces their mass to NaN so
  they match **no** gate (per-fragment exclusion in every panel,
  including pair panels via the pair-AND).
- **Suppressed channel = the n = 0 ensemble.** Suppressed ions ride at
  their handover n in `detection.npz`; the view forces their mass to
  the bare gate `m(0)`, realizing the frozen suppressed → n = 0
  convention with their stored detected velocities. They populate the
  n = 0 panels (size distribution bin 0 as today; *newly*, the n = 0
  speed-distribution panel — the RQ3 bare-peak observable, which the
  legacy summary could never populate). Panels label the channel
  explicitly.
- **Two-detected-fragments convention (cov panels; amended —
  supersedes the partner-ballistic convention).** `detection.npz`
  arrays are shape (2N,): *both* iodine fragments of every molecule
  are independent I⁺Heₙ detections (the E1 convention carried through
  the detection stage). Cov panels pair the molecule's two **detected**
  rows; a pair enters only when both fragments are non-retained
  (automatic through the per-fragment NaN-mass exclusion).
  `plot_detection_summary.py` does **not** load `ion.npz`.
- **Tier-3 caveat label.** VMI/polar (and cov) comparisons render with
  an explicit *"ensemble second moments under-dispersed — Tier-3 noise
  stubbed"* annotation. They are previews, not Tier-3 validation.

## 4. Section roster (`detection_summary.pdf`, in order)

| # | Section | Source | Convention notes |
|---|---------|--------|------------------|
| 1 | `detection_metadata` | cfg.json + all three artifacts | stage identity, t_detect, scored/retained/suppressed counts, drag/ladder/prior knobs |
| 2 | `detected_size_distribution` | detection read | existing D4 builder (suppressed bin + solvated vs abundance, W1/n1/ratio) |
| 3 | `detected_ihe_ked_mean_energy` | detection read | existing D4 builder (committed error model, I77 median anchor, both χ²) |
| 4 | `detected_ke_anatomy` | detection read | existing D4 builder (signed per-bin z) |
| 5 | `detected_ihe_ked_curves_3d` | detection read | legacy recipe, n = 0 = suppressed channel; ion counts per the legacy recipe (sim curve labels carry N, legend on the n = 0 panel; an empty gate annotates "sim: no atoms in gate" on its own panel — amended 2026-07-22 post-review: the original "per-panel ion counts" wording overstated the verbatim-moved recipe, which renders one legend) |
| 6 | `detected_ihe_ked_curves_2d` | detection read | as 5, in-plane projection |
| 7 | `detected_mass_resolved_velocities` | detection read | legacy recipe, sim-only |
| 8 | `detected_paper_v2_vmi` | detection read | Tier-3 caveat label |
| 9 | `detected_paper_v2_polar` | detection read | Tier-3 caveat label |
| 10–14 | `detected_paper_cov_*` (5 panels) | detected view (both fragments) | two-detected-fragments pairing + Tier-3 caveat |

PNG naming: `detected_*` prefix throughout. The three D4 sections keep
their existing filenames; `run_summary.pdf` and `detection_summary.pdf`
both emit them via the same shared builders, so the files are
regenerated with identical content whichever script runs last.

## 5. Testing

- **Legacy invariance:** section-list + builder-identity regression for
  non-detection dirs across the refactor (extends the existing D4
  gating tests). No byte-comparison of PNGs; the lock is at the
  section-list/builder level as today.
- **New tests:** velocity-vector read extension (shapes, units,
  retained exclusion); suppressed-channel n = 0 panel convention; the
  partner-ballistic cov assembly; retitle/overlay-drop behavior of the
  four legacy sections on a detection dir; `plot_detection_summary.py`
  fail-loud on a dir without `detection.npz`; smoke render of every
  roster section on synthetic fixtures.
- Fixtures: the `test_tier2_confirmation` synthetic builders, Agg
  backend, figures closed, no production artifacts (testing rules).

## 6. Governance

- **No schema change, no physics change.** Post-processing extension
  under the drag-phase scoped exception; the "post-processing is
  bug-fix-only" freeze is explicitly relaxed for this spec by user
  adjudication (this document).
- **D4 follow-up adjudication (recorded here):** the D4 option-(a)
  wording "beside the existing ion-stage mass spectrum" is superseded
  on detection runs by §2.3's overlay-drop + retitle. Non-detection
  dirs are untouched, so the D4 lock tests keep their meaning.
- Rule-2: no declared-but-unread config fields are introduced.
- Implementation, as always, only on `[PROCEED TO IMPLEMENTATION]`.

## 7. Out of scope

- Any new scoring metric on detected speed distributions (display
  only; scoring stays `tier2_confirmation.py`).
- Tier-3 noise (stays stubbed behind its enum).
- Reordering or restyling `run_summary.pdf` beyond §2.3.
- Abel inversion / full VMI image interpretation (standing scope rule).
