# Detection-Stage Summary — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans
> (inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`)
> syntax for tracking. Spec: `TIER2_DETECTION_SUMMARY_SPEC.md` (as amended by
> Task 1 of this plan).

**Goal:** Build `detection_summary.pdf` (full detector-facing figure set from
the detection stage) + de-mislead `run_summary.pdf` on detection runs, per the
frozen spec.

**Architecture:** A checkpoint-shaped `DetectedEnsembleView` (package) lets
every existing mass-gated plotting recipe run verbatim on the detection
ensemble; the section builders move to a shared `summary_sections.py`
consumed by both entry scripts. No postprocess recipe changes, no schema
change.

**Tech stack:** numpy / matplotlib (Agg in tests) / pytest; interpreter
`C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe`.

## Global constraints

- Non-detection run dirs render **byte-identically** through `plot_run_summary.py` (existing D4 gating contract extends over the refactor).
- Rule 1: zero duplicated plotting recipes or physics conventions.
- No checkpoint-schema, RNG, or physics changes. `detection.npz` is read-only input.
- TDD per task; tolerances justified; no figures or production checkpoints written by tests.
- Every task ends green on the narrow tests; full suite at close-out.
- Commit per task (project slice convention).

---

### Task 0: Checkpoint commit of pending work

The post-review fixes + spec + log entries are uncommitted. Commit them
as-is before new code.

- [ ] `git add` the six pending files (2 test files, `tier2_confirmation.py`, `plot_run_summary.py`, `tier2_confirmation_score.py`, tier-2 log) **plus** `TIER2_DETECTION_SUMMARY_SPEC.md` and this plan.
- [ ] Commit: `Post-production review fixes (median-anchor bands-on lock, nan/ratio guards) + detection-summary spec frozen + implementation plan`.

### Task 1: Spec amendments (facts learned from the code survey)

**Files:** Modify `docs/drag_port/Tier2/TIER2_DETECTION_SUMMARY_SPEC.md`.

Three amendments, each replacing the corresponding spec text:

- [ ] **Partner-ballistic convention → two-detected-fragments convention.** `DetectionResult` arrays are shape (2N,): *both* iodine fragments of every molecule are independent I⁺Heₙ detections (the E1 convention, `detection_stage.py:157-201`). Cov panels pair the molecule's two **detected** rows; a pair enters only when both fragments are non-retained (automatic through the per-fragment exclusion below). `plot_detection_summary.py` does **not** need `ion.npz`.
- [ ] **n-selection realized as a checkpoint-shaped view.** Selection runs through the existing single mass-gate convention (`select_final_mass_gate`) applied to a `DetectedEnsembleView` whose `mass_final_kg` is derived from `n_detected`/`state_reason`; exactness `mass == m(n)` is pinned by test. The "no mass-gate machinery" clause is amended to "no *new* gate machinery; the existing shared gate is the selection mechanism, provably equivalent to `n_detected` selection".
- [ ] **Suppressed/retained realization.** Suppressed ions ride at their handover n in `detection.npz` (`state_reason` docstring); the view forces their mass to the bare gate (`m(0)`), realizing the frozen suppressed = n = 0 convention with their stored detected velocities. Retained ions get mass = NaN (match no gate — excluded from every panel, per-fragment).

### Task 2: `DetectedEnsembleView` (package)

**Files:**
- Create: `i2_helium_md/postprocess/detected_view.py`
- Modify: `i2_helium_md/postprocess/__init__.py` (export)
- Test: `tests/test_detected_view.py`

**Produces (later tasks rely on):**
```python
@dataclass(frozen=True)
class DetectedEnsembleView:
    """Checkpoint-shaped view of the detected ensemble (duck-types the
    IonCheckpoint surface the mass-gated diagnostics read)."""
    num_molecules: int          # N (arrays are 2N, checkpoint layout)
    mass_final_kg: np.ndarray   # (2N,) see forcing rules
    velocities_final_x: np.ndarray  # (2N,) [A/ps]
    velocities_final_y: np.ndarray
    velocities_final_z: np.ndarray
    b_ion_outside: np.ndarray   # (N,) all True (exclusion is by mass)

def detected_ensemble_view(det: DetectionResult) -> DetectedEnsembleView: ...
```
Forcing rules in the factory: start from `mass_detected_kg`;
`state_reason == "suppressed"` → `MASS_I * U_KG` (bare gate);
`"droplet_retained"` → `np.nan`; `"frozen"`/`"time_exhausted"` → unchanged.
Velocities are `v{x,y,z}_detected` pass-through. Validate shapes (2N,)
against `num_molecules` (rule 4).

- [ ] **Failing tests first** (`tests/test_detected_view.py`, reuse `t2c.make_detection` builders):
  - `test_shapes_and_passthrough` — view arrays (2N,), velocities identical to `v*_detected`.
  - `test_frozen_mass_equals_m_of_n_exactly` — for frozen/time_exhausted rows, `mass_final_kg / U_KG == 126.90447 + 4.002602 * n_detected` exactly (the n_detected ↔ gate equivalence pin; use the constants imported from the mass model, not literals).
  - `test_suppressed_forced_to_bare_gate` — suppressed row passes `select_final_mass_gate(view, mass_amu=m(0))`, fails its handover-n gate.
  - `test_retained_matches_no_gate` — retained row is in no gate for any n 0..21.
  - `test_end_to_end_with_existing_recipes` — `compute_final_velocity_histogram(view, mass_amu=m(1), ...)` and `fragment_mean_kinetic_energy(view, 1)` run and count exactly the non-retained n = 1 rows.
- [ ] Run: expect FAIL (module missing). Implement. Run: PASS.
- [ ] Commit: `feat(postprocess): DetectedEnsembleView — checkpoint-shaped detected ensemble (spec Task 2)`.

### Task 3: Shared `summary_sections.py` + byte-identical refactor

**Files:**
- Create: `scripts/post_processing/summary_sections.py`
- Modify: `scripts/post_processing/plot_run_summary.py`
- Test: existing `tests/test_plot_run_summary_detection.py` and the `_nan_aware_moving_mean` pin test must pass **unchanged**.

Move verbatim (no recipe edits) from `plot_run_summary.py` into
`summary_sections.py`: `_SectionSkipped`, `_nan_aware_moving_mean`, the
recipe constants they read (`IHE_KED_*`, `HIST_*`, `MASS_SPECTRUM_MAX_AMU`,
`VELOCITY_PLOT_V_MAX_APS`, `PAPER_V2_MASS_AMU`, `EXPERIMENTAL_NOISE_FLOOR`,
`DETECTED_KE_*`), `_load_detection_read`, `_detected_ke_scores`,
`_required_paper_cov_reference`, and the builders: `_section_mass_spectrum`,
`_section_ihe_ked_mean_energy`, `_section_ihe_ked_curves`,
`_section_paper_v2_vmi`, `_section_paper_v2_polar`,
`_section_paper_cov_*` (5), `_section_mass_resolved`,
`_section_detected_size_distribution`, `_section_detected_ked_mean_energy`,
`_section_detected_ke_anatomy`. Window-only sections (metadata, energies,
temperature, radial, pair distance, boltzmann, hedft) stay in
`plot_run_summary.py`.

- [ ] Create module; `plot_run_summary.py` replaces the moved code with `from summary_sections import ...` aliases under the **same names** (the tests do `mod._section_detected_size_distribution` etc. on the loaded script module — aliases keep them and the section list working).
- [ ] Run: `pytest tests/test_plot_run_summary_detection.py tests/test_tier2_confirmation.py -q` and the movmean pin test → all PASS unchanged.
- [ ] Commit: `refactor(post-processing): shared summary_sections module (rule 1; run_summary behavior unchanged)`.

### Task 4: Stage notes + handover retitles in run_summary

**Files:**
- Modify: `scripts/post_processing/summary_sections.py`, `scripts/post_processing/plot_run_summary.py`
- Test: `tests/test_plot_run_summary_detection.py` (new class)

**Interfaces produced:**
```python
# summary_sections.py — added keyword-only params, all defaulting to
# legacy behavior (None/True), so Task 3 call sites are untouched:
_section_mass_spectrum(ion, abundance, *, stage_note=None)
_section_ihe_ked_mean_energy(ion, ked_ref, *, stage_note=None, include_reference=True)
_section_ihe_ked_curves(ion, ked_dir, ked_ref, representation, *, stage_note=None, include_reference=True)
_section_paper_v2_vmi(..., *, stage_note=None)          # likewise polar, cov x5,
_section_mass_resolved(ion, *, stage_note=None)          # mass_resolved
_apply_stage_note(fig, stage_note)  # fig.text banner, no-op on None
```
`include_reference=False`: mean-energy section plots sim points only (no
reference errorbars/bands, no `_SectionSkipped` on `ked_ref=None`); curves
section plots sim histograms + sim v(⟨E⟩) markers only. Titles gain the
handover retitle when `include_reference=False`.

In `plot_run_summary.py` (section list, detection runs only —
`detection_read is not None`):
- `mass_spectrum` → `_section_mass_spectrum(ion, None, stage_note=_HANDOVER_NOTE)` (sim-only branch);
- `ihe_ked_mean_energy`, both `ihe_ked_curves` → `include_reference=False, stage_note=_HANDOVER_NOTE`;
- VMI/polar/cov×5/`mass_resolved` → `stage_note=_HANDOVER_NOTE` (content unchanged);
- `_HANDOVER_NOTE = "ion-stage handover state (t = t_handover) — pre-detection diagnostic; detector comparison: detection_summary"`.

- [ ] **Failing tests first** (`TestHandoverRetitles` in `tests/test_plot_run_summary_detection.py`): on a synthetic detection run, mass_spectrum figure contains the stage note text and no experiment bars; `include_reference=False` mean-energy renders with `ked_ref=None`; curves figure has no "experiment" legend entries; a no-detection dir keeps legacy behavior (existing tests re-assert this).
- [ ] RED → implement → GREEN.
- [ ] Commit: `feat(run-summary): handover retitles + stage notes on detection runs (spec §2.3; legacy dirs unchanged)`.

### Task 5: `plot_detection_summary.py`

**Files:**
- Create: `scripts/post_processing/plot_detection_summary.py`
- Test: `tests/test_plot_detection_summary.py`

USER SETTINGS mirror `plot_run_summary.py` (run dir + the same reference
paths). Flow: `RunDirectory` → **fail loudly** (`FileNotFoundError` with the
spec sentence) if `detection.npz` missing → `load_detection_result` →
`read_confirmation_detection` (D4 sections) + `detected_ensemble_view`
(recipe sections) → roster (spec §4 order):

1. `detection_metadata` — cover: run dir name, t_handover_ps, detection_time_ps, state_reason counts/fractions, scored/retained totals, knob lines from cfg.json (drag_form, v_c, p_tail, tau_ps, ladder, f_int, prior — reuse `_knob_columns` idiom from `tier2_confirmation_score.py` by import? No: scripts must not import scripts — inline the 7 cfg.get lines, they are config-reads not physics).
2–4. the three D4 detected builders (from `summary_sections`).
5–6. `_section_ihe_ked_curves(view, ..., "3d"/"2d", stage_note=_DETECTED_NOTE)` — n = 0 panel is the suppressed/bare channel; `_DETECTED_NOTE` says so.
7. `_section_mass_resolved(view, stage_note=_DETECTED_NOTE)`.
8–9. VMI + polar with `stage_note=_TIER3_NOTE` ("detected ensemble; second moments under-dispersed — Tier-3 noise stubbed").
10–14. cov ×5 with `stage_note=_TIER3_NOTE` (pairs = the molecule's two detected fragments).

Output: `figures/detection_summary.pdf` + `detected_*.png` / `detection_metadata.png` per section (same save loop idiom as run_summary).

- [ ] **Failing tests first** (`tests/test_plot_detection_summary.py`, importlib-load the script like the run_summary tests): `test_fails_loudly_without_detection_npz` (tmp dir, `pytest.raises(FileNotFoundError, match="detection.npz")`); `test_metadata_section_renders` (synthetic detection + minimal cfg.json); `test_curves_sections_render_from_view` (synthetic detection, synthetic ked refs — suppressed ions appear in the n = 0 panel: assert the n = 0 axis has a sim line/N annotation); `test_reference_less_sections_skip` (VMI/cov refs None → `_SectionSkipped`); `test_section_list_order_matches_spec`.
- [ ] RED → implement → GREEN.
- [ ] Commit: `feat(post-processing): plot_detection_summary.py — detection-stage figure set (spec §4 roster)`.

### Task 6: Close-out

- [ ] Full suite: `python -m pytest -q` → green (report the count).
- [ ] Generate the real document once for the blessed run (`finc1v725`) — visual sanity only, artifacts stay uncommitted (data-contract rule).
- [ ] Spec status line → DELIVERED (with date); tier-2 log entry: DELIVERED record (what shipped, conventions realized, test counts, the Task-1 amendments called out).
- [ ] Commit: `Detection summary DELIVERED: shared sections + DetectedEnsembleView + detection_summary.pdf; run_summary handover retitles (log + spec updated)`.

## Self-review notes

- Spec §3 read-extension of `ConfirmationDetectionRead` is **superseded** by
  `DetectedEnsembleView` (Task 1 amendment covers this — the read keeps its
  scoring role; the view carries velocities).
- Spec coverage: §1 n/a; §2.1→Task 3, §2.2→Task 5, §2.3→Task 4; §3→Tasks 1+2;
  §4→Task 5; §5→each task's tests + Task 6; §6→Tasks 0/1/6; §7 respected
  (no new scoring, no noise, no reorder beyond §2.3).
- Type consistency: `view` duck-types `IonCheckpoint` for exactly the five
  attributes the shared gate reads (verified against
  `select_final_mass_gate`, `_paper_v2_atom_selection`, `pair_correlation`,
  `paper_v4` call sites).
