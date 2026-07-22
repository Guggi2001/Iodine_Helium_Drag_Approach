# Tier-2 Detection-Stage Summary — Design Spec

**Status: DESIGN FROZEN 2026-07-22 — awaits its own `[PROCEED TO IMPLEMENTATION]`.**
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

- **Detection ensemble read.** `ConfirmationDetectionRead`
  (`i2_helium_md/postprocess/tier2_confirmation.py`) is extended to
  surface the per-ion terminal velocity vectors and detected masses
  already present in `detection.npz` (`vx/vy/vz_detected` [Å/ps],
  `mass_detected_kg`, `n_detected`, `state_reason`). No schema change —
  read-side only. Shape/unit validation on load (rule 4).
- **n-selection.** Sections select by integer `n_detected` directly.
  The legacy m(n) ± 0.5 amu mass-gate machinery does **not** come
  along (it exists only because `ion.npz` has no fragment label).
- **Retained ions** (`state_reason = droplet_retained`) are excluded
  from every panel by the standing detection contract.
- **Suppressed channel = the n = 0 ensemble.** Suppressed ions are
  bare I⁺ at detection; they populate the n = 0 panels (size
  distribution bin 0 as today; *newly*, the n = 0 speed-distribution
  panel — the RQ3 bare-peak observable, which the legacy summary could
  never populate). Panels label the channel explicitly.
- **Partner-ballistic convention (cov panels).** The neutral partner
  feels no forces after the MD window, so its handover velocity
  (`ion.npz`) *is* its `t_detect` velocity. Cov panels pair detected
  ion velocities with handover partner velocities and carry a one-line
  label recording this. `plot_detection_summary.py` therefore also
  loads `ion.npz` — for partner velocities only.
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
| 5 | `detected_ihe_ked_curves_3d` | detection read | legacy recipe, n = 0 = suppressed channel, per-panel ion counts (no silent caps) |
| 6 | `detected_ihe_ked_curves_2d` | detection read | as 5, in-plane projection |
| 7 | `detected_mass_resolved_velocities` | detection read | legacy recipe, sim-only |
| 8 | `detected_paper_v2_vmi` | detection read | Tier-3 caveat label |
| 9 | `detected_paper_v2_polar` | detection read | Tier-3 caveat label |
| 10–14 | `detected_paper_cov_*` (5 panels) | detection read + handover partners | partner-ballistic label + Tier-3 caveat |

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
