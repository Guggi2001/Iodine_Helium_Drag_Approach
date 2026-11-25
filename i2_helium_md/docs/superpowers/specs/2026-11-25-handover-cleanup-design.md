# Handover cleanup — design

Date 2026-11-25. Status: **general design approved by the user 2026-11-25**; all decisions of §8 closed; allowlist complete, awaiting sign-off. Planning only; code
needs `[PROCEED TO IMPLEMENTATION]`, and commits, pushes and tags each need an
explicit request. Once approved, this replaces the phase table and delete list of
`docs/POST_THESIS_CLEANUP.md` C5. The C5.1–C5.3 detail (presets, battery,
defaults), as amended by C5.7, stays the reference except where this file differs.

Markers: **[U]** = user decision, **[P]** = proposal for the user to check.

## 1. Intent

[U] A successor finds a clean codebase: a modular, well-documented Python
package. Unnecessary and outdated run scripts and documentation are deleted. What
remains is documentation that explains the physics or crucial scripts, and the
scripts that reproduce thesis findings.

**Goals**

- [U] On a new PC, 2 or 3 simple scripts regenerate all data from a fresh clone.
  No data archive is required. Thesis laws run by name (B2 / B3 / h405; B2 is "LQ"
  in the thesis).
- [P] The scripts: (1) drag batteries B2 / B3 / h405 with pool and score; (2)
  legacy hard-sphere and Tier-0 runs; (3) all thesis figures. The legacy runs get
  fixed seeds in their presets or driver: `single_pulse_droplet_long`,
  `9A_hedft_comparison` and `18A_hedft_comparison` were unseeded, so their new runs
  agree statistically with the submitted figures, not bit for bit.
- Everything in the repo is there for a stated reason (the allowlist, §4).
- Package: only model options the thesis uses or quotes ([U] option B).

**Non-goals.** New physics; Tier 3; performance work; re-running the production
batteries; changing constants, RNG draw order or checkpoint schema; open physics
questions (C3 stays an open question in the docs).

**Anti-goals.** A thesis number moves; a finished run stops loading; something is
deleted before its replacement is proven; dead branches survive because a test
still covers them.

## 2. State found (2026-11-25)

- `drag_implementation` is 413 commits ahead of `main`; `main` has 1 commit of its
  own (`80fc786`). `git ls-remote --tags origin` lists no tags, so
  `pre-drag-cleanup` and `pre-tddft-move` are local only.
- Outside git: `data/runs/` run directories (9.6 GB on disk in total; 102 files are
  tracked, 101 of them in `h2b_forward_model/`) and `legacy_matlab_repository/`
  (282 files).
- Python: package 67 files / 27,302 lines; `scripts/` 120 / 43,842;
  `tests/` 132 / 38,202 (3,342 tests collected, Python 3.13.0).
- Markdown under `docs/`: 61,061 lines. No README in `i2_helium_md/`; the
  repo-root README is a 306-byte stub; `PROJECT_OVERVIEW.md` describes
  `finc1v725` as production (2026-08-13).
- `data/figures/README.md` already maps all 43 thesis figure files to script and
  inputs. It is the seed of the allowlist.
- `tier2_b2_thesis_figure.py` imports one constant, `NBAR_LANDING`, from
  `tier2_g414_table.py`, and through it 10 `gen_*` modules. Four more figure
  scripts inherit this: `tier2_law_comparison_figure`, `tier0_law_figure`,
  `legacy_ke_size_figure`, `legacy_ked_shapes_figure` (plus `tier0_thesis_figure`
  via `tier0_law_figure`). `tier2_g6_chord_map.py` reaches 15 `gen_*` modules.
- Legacy figures read hard-sphere runs (`single_pulse_droplet_long`,
  `{9A,18A}_hedft_comparison`), so those are thesis inputs too.
  `18A_hedft_comparison` was missing and the user regenerated it on 2026-11-25.
  `fig_continuum_drag.*` renders again but does **not** reproduce the committed
  CSV, for 9 Å as well as 18 Å: both runs are unseeded (`seed: null`) and newer
  than the committed figure (2026-10-23; 9 Å run 2026-10-26). See D7.
  `legacy_ke_and_size_points.csv` and `legacy_ked_shapes_points.csv` reproduce
  exactly from the May `single_pulse_droplet_long` run (unseeded). That run was
  overwritten on 2026-11-25 and rerun by the user, so these two figures now agree
  statistically only; the committed files stay as submitted.
- `SimConfig` has 122 fields; B2 / B3 / h405 each store 122. Unknown fields make a
  `cfg.json` fail to load (C6 item 8).
- The repo root also tracks `Drag_function/` and `molecular_dynamics_lammps/`.

## 3. Approaches

1. **Prune in place, outermost first.** Small diffs; but thesis figures import the
   code being deleted, and leftovers are likely.
2. **New repository, copy in what is needed.** Cleanest; loses `git blame`, forces
   all tests to be re-sorted at once, no step-by-step comparison.
3. **Same repo: build the new path, prove it, delete everything off an allowlist.**
   [U] **Chosen.** Every deletion follows a proven replacement, and "unnecessary"
   is a checkable list.

Work happens on a new branch `handover` from tag `thesis-final`.
`drag_implementation` stays untouched, so thesis work can continue there.
[P] A thesis-driven fix on `drag_implementation` after the tag is cherry-picked to
`handover`, the tag moves only on user request, and P3/P5 gates are re-run.

## 4. Allowlist rule

A file or option stays only if it is needed by at least one of:

1. a file in `data/figures/thesis_*` (43 files);
2. a number quoted in the thesis;
3. running the model: presets, driver, package, hard-sphere baseline;
4. tests of 1–3;
5. documentation of §7.

Output of P1: `docs/handover/ALLOWLIST.csv`, columns
`item, kind (figure|number|script|module|option|data|doc|test), needed_by, inputs,
decision (keep|rewire|move|delete), note`. Every tracked file and every `SimConfig`
field has exactly one row or is covered by exactly one folder row. The user signs
it off before any deletion. It is removed in P6, once used.

History: [P] full history lives in git tag `thesis-final`, with a one-page pointer
in the README. No per-script removal manifest; this replaces C5's
`SCRIPT_REMOVAL_MANIFEST.md`.

## 5. Package (option B)

Initial reading, settled per option by P1:

| Area | Keep | Remove |
|---|---|---|
| `drag_form` | `linear_quadratic`, `capped_cubic` | other five, unless quoted |
| `mass_scenario` | `biphasic`, `fixed` | `anchored_discrete`, `anchor_*`, `shell_schedule` |
| `noise_*` | none | all stubs; Tier 3 becomes design notes |
| E2 relaxation | none | module, fields, guards (C1) |
| `ce_channel_mode`, `exit_strip_mode`, `drag_state_coupling` | off behaviour | on-branches and fields |
| `droplet_size_sampler_mode`, `initial_shell_model`, `self_bound_gate_release`, ladder / pickup / cooling options | thesis value plus quoted variants | rest |
| Hard-sphere path, neutral stage | keep, including the collision reader of `v_limit` (C1 warning) | none |
| Extraction | `lq_fixed_c` and what `method_b_g412_ridge.py` needs | force balance, old fit variants (C6.3, C6.4) |
| Post-processing | used by kept scripts | modules used only by their own tests |
| Twin, `tier2_h2b_forward_model.py` (5,437 lines) | none; the method is described in the docs | **[U]** all of it; cheap MD runs replace it |

Rules:

- **Single-valued enums disappear.** If one value remains, the field and its
  dispatch go and the code states the model directly. This deliberately ends
  CLAUDE.md's "every choice behind an enum" rule, which served the selection phase.
- **`cfg.json` migration.** One explicit table, `RETIRED_FIELDS: name -> inert
  value`. A stored retired field loads only at its inert value, otherwise the load
  fails naming field and value. Missing fields behave as today. Stored `cfg.json`
  files are never rewritten. `HISTORICAL_DRAG_BUNDLE_IDS` stays.
- **Retired but not inert.** The thesis members do store non-default values for
  fields that take no effect (for example `relaxation_dissipation =
  "landau_gated_drag"` with `relaxation_stage_enabled = false`). The table takes a
  predicate for these: "ignored when `<switch>` is off".
- **Drag off by default** (C5.3).
- **`config.py`** (2,529 lines) is split by concern only after pruning; public
  imports from `i2_helium_md.config` stay valid. No other refactor of module
  structure unless P1 finds duplicated physics.
- The `force_balance` value of `extraction_method` stays loadable (C6.4).

## 6. Phases and gates

One phase = one plan = one or more commits, each green. Baselines from P0.

| Phase | Work | Gate |
|---|---|---|
| **P0 Secure** | Tag `thesis-final` (thesis is submitted) and push it with the two existing tags, on explicit request. Branch `handover`. Record baseline: suite result with skips by name; hashes of all 43 thesis figure files and a fresh render of each. `data/runs/` inputs stay on disk untouched until P7 | baseline files written; each figure classified identical / differs / not rendered, with the reason |
| **P1 Allowlist** | §4 table | user sign-off; every tracked file and config field covered |
| **P2 New path** | Presets, battery module, driver: C5.1 + C5.2 as amended by C5.7. Shared basics move first (C5.7). Additive only | all 18 members' `cfg.json` rebuilt, every field; scores equal the decoupled rows (`B2`, `B3`, `h405`) of `g13_decoupled_table.csv` to 4 dp; pooled B2 rebuilt array for array; one B2 s1 rerun gives an array-identical `detection.npz`; N = 20 smoke run |
| **P3 Rewire** | Kept scripts read from package or driver, no `gen_*` import. `NBAR_LANDING` and similar constants move to their owner. Includes `g6_*` and `g13_*` renders, which C6 did not re-render | figures pixel-identical, CSVs byte-identical to P0; import check: no kept script reaches a non-allowlisted module |
| **P4 Delete scripts and data** | Scripts, tests, tables, figures and root files off the allowlist. Tests of package physics that import `tier2_common` are rewired first (C5) | suite green; P3 gate again |
| **P5 Prune package** | §5, one option family per commit; then `config.py` split | after every commit: suite green, all 18 members load, every rebuilt preset equals its stored cfg on surviving fields. At the end: P2 and P3 gates, including a second B2 s1 rerun |
| **P6 Docs** | §7: write first, then delete | every command in the docs run from a fresh clone; no dead links; no reference to a deleted symbol |
| **P7 Handover** | Integrate to `main` (decide on `80fc786`), tag `handover-v1`, cold-start test | someone who did not write it installs, runs the smoke preset and reproduces one figure using only the README |

Estimate by weight: P3 and P5 heavy, P2 and P6 medium, rest light. A gate that
cannot be met stops the phase; tolerances are not loosened.

MD cost: two B2 s1 reruns in total. Pixel identity holds for one machine and
matplotlib version; [P] P0 records the versions and the package gets a tested
lower bound plus a lock file for exact reproduction.

## 7. Documentation after P6

| File | Content | Main source |
|---|---|---|
| `i2_helium_md/README.md` | what it is, install, quick start, folder map, where history is | new |
| `docs/physics.md` | governing equations with units, stages, mass model, detection, conventions | `PHYSICS_BASELINE`, `MASS_DYNAMICS_LOCKED…`, `docs/physics/`, `PROJECT_OVERVIEW` |
| `docs/calibration.md` | Tier-0 ridge, parameter provenance, B2 / B3 / h405 | `CALIBRATION_MAP`, Tier-0 docs, G4.12 |
| `docs/reproduce.md` | presets, driver, data archive, figure table | `data/figures/README.md` |
| `docs/limits_and_open_questions.md` | known misfits, RQ / M questions, Tier-3 notes, C3 | `THESIS_DISCUSSION`, `RESEARCH_QUESTIONS`, `TIER2_MASS_SCENARIOS` |
| `docs/reference_data.md` | provenance of `data/reference` | existing sidecars, CLAUDE.md |
| `CLAUDE.md` | durable rules only, short | current file |

Module pages under `docs/` are kept only where they explain something the
docstring cannot; otherwise the content moves into the docstring. Logs, phase
plans, slice docs, atlas documents and the register are deleted (tag
`thesis-final`). Docs quote only committed numbers and name their source.

## 8. Open decisions (do not block approval; needed by the phase named)

| # | Question | Default [P] | Needed by |
|---|---|---|---|
| D1 | Source of quoted numbers | **[U] decided 2026-11-25:** the last numbers in `THESIS_DISCUSSION.md` are the thesis numbers; where sections disagree the latest wins | done |
| D2 | Data archive | **[U] 2026-11-25:** not needed; left open. Data must be regenerable (§1) | closed |
| D3 | `Drag_function/` | **[U] 2026-11-25:** outside `i2_helium_md`, left alone | closed |
| D4 | Other folders | **[U] 2026-11-25:** nothing outside `i2_helium_md` is touched (this includes `molecular_dynamics_lammps/` and the repo-root files; own reading of D3's answer). `old_data/`, `scripts/cross_reference/`: delete | closed |
| D5 | Legacy MATLAB source | **[U] 2026-11-25:** no action; the successor has access to the MATLAB git. Nothing at run time reads it (it appears only in comments and docstrings) | closed |
| D6 | Thesis submitted? | **[U] 2026-11-25:** yes; `thesis-final` marks the submitted state. The exact commit is confirmed when tagging | closed |
| D7 | `fig_continuum_drag` | **[U] 2026-11-25:** no dramatic change; re-render and compare. Done: the committed figure and the render from the runs on disk differ in 0.66 % of pixels; per-bin force chi² 5.8 / 10 (9 Å) and 9.2 / 8 (18 Å); same conclusion. The committed file stays as submitted | closed |

## 9. Left to the builder

File placement inside the package, test organisation, commit slicing within a
phase, helper names, format of the hash list. Not free: gates, migration rule,
allowlist sign-off, draw order, constants, schema.
