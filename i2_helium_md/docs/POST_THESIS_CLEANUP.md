# POST_THESIS_CLEANUP — deferred cleanup register

**Purpose.** Jobs that are deliberately **not** done during thesis writing,
because they would change code paths that produced thesis numbers, but that
must not be forgotten. Each entry: *what*, *why deferred*, *depends on*, *done
when*. This is a register, not a plan — an entry is executed only on an
explicit user decision (and, for code, `[PROCEED TO IMPLEMENTATION]`).

Created 2026-09-30 from the detector-stage-fix discussion
(`docs/drag_port/Tier2/TIER2_DetectorStageFix.md` §3a–§3d).

Status legend: **OPEN** · **BLOCKED** (dependency not met) · **DONE** (keep
the entry, add the commit / date).

---

## C1. Retire the E2 relaxation stage — BLOCKED

**What.** Remove the Stage-II / E2 relaxation stage once production runs on
the skip path (Stage I to ~0.5 ns → detection seeded from `ion.npz`):

- `i2_helium_md/simulation/relaxation_stage.py` (whole module) and its
  stream key `RELAXATION_STREAM_KEY` / checkpoint hooks in
  `simulation/checkpoint.py`;
- config fields `relaxation_stage_enabled`, `relaxation_time_ps`,
  `relaxation_dt_ps`, `relaxation_forces`, `relaxation_dissipation`, their
  enums (`RelaxationForces`, `RelaxationDissipation`) and the config-load
  guards (`check_relaxation_config` and the `landau_gated_drag` cross-guard);
- the `relaxation.npz` branch of the detection stage's seed handling
  (`detection_stage.py`: "relaxation.npz when E2 ran, or ion.npz");
- tests: `test_relaxation_stage.py`, `test_relaxation_dissipation_config.py`,
  `test_relaxation_landau_gamma.py`, and the E2 parts of
  `test_detection_stage.py`, `test_state_coupling.py`, `test_tier2_common.py`,
  `test_gen_tier2_*`, `test_gen_tier2atlas_*`,
  `test_tier2_size_distribution_table.py`, `test_tier2_staircase_probe.py`;
- `run_relaxation_stage` calls in the historical generators
  (`scripts/gen_tier2_*.py`, `scripts/gen_tier2atlas_*.py`, 20 files) —
  decide per file: delete, or freeze as a documented historical record.

**Why deferred.** Every committed Tier-2 atlas / finals number (incl. the
standing point and h405) was produced with E2. Until the thesis is final, the
E2 path must stay runnable to reproduce them. Until then it is kept and
**marked legacy** (docstring banner + config-field comments — part of the
pipeline implementation).

**Must NOT be removed with it.** `v_limit_m_per_s` / `v_limit_angstrom_per_ps`
also define `cfg.E_min_eV`, the Landau kinetic-energy cutoff of the hard-sphere
collision step, read by the **neutral stage** (`simulation/propagation_step.py`)
and the ion-stage collision path (`simulation/ion_propagation_step.py`). Only
the E2 reader (`landau_gated_drag`) goes.

**Depends on.** Thesis final; ~~production pipeline switched to the skip
path~~ (done 2026-09-30, `TIER2_DetectorStageFix.md` §5.7); decision on whether any E2-era
number must remain reproducible from code (vs. from committed artifacts).

**Done when.** No `relaxation` symbol left outside history docs; full
`pytest` green; `drag_migration_log_tier2.md` retirement entry; CLAUDE.md
rule-2 carries updated.

## C2. Stale "E2 is zero-gamma" docstrings — DONE (2026-09-30, with the pipeline implementation, `TIER2_DetectorStageFix.md` §5.4 D7)

**What.** `detection_stage.py` `escape_energetics`, `_conservatively_bound`
and the `EscapeEnergetics.asymptotic_ke_eV` field doc state that E2 is
zero-gamma (no drag). False since `landau_gated_drag` (2026-07-20) and
meaningless on the skip path. The *criterion* remains valid and conservative
(drag after handover could only make escape harder); only the justification
text is wrong.

**Why deferred.** Docstring-only; bundle it with the pipeline implementation
(the closure touches the same module) or with C1.

**Depends on.** Nothing (can be done with the pipeline implementation).

**Done when.** Docstrings describe the handover-time criterion without
reference to E2's dissipation arm.

## C3. Pair-aware escape criterion (replace the full-credit rule) — OPEN, optional

**What.** `escape_energetics` credits the full pair Coulomb energy to *both*
fragments (a deliberate conservative over-estimate). For molecules with
**both** fragments trapped (a doubly charged droplet) this double counting
labels them `droplet_retained_marginal` ("unbound") although they oscillate
in place indefinitely (`TIER2_DetectorStageFix.md` §2f: across 88 finished
runs, 1860 / 1867 marginal ions have a marginal partner, the other 7 a
bound partner). Replace with a molecule-level criterion (pair total energy vs
the joint barrier), so that bound / marginal reads as physics.

**Why deferred.** Score-inert (both classes are excluded). The thesis reports
"one fragment retained" vs "both fragments retained" instead of the split.

**Depends on.** A need to report the split as physics.

**Done when.** Criterion implemented behind an enum, validated on the
committed runs (marginal class should collapse into bound for trapped pairs).

## C4. `v_limit_m_per_s` does double duty — OPEN, note

**What.** One field sets two different physics: the neutral/ion collision
cutoff `E_min_eV` (legacy MATLAB, `E_min = m_I·v_L²/2`) and, until C1, the E2
drag gate. The finals generators override it to 58 m/s (config default 40),
which therefore also changes the neutral-stage collision cutoff. Decide
whether these should be separate fields.

**Why deferred.** Changing it alters neutral-stage behaviour for every run;
not a thesis-phase change.

**Depends on.** C1 (after which only the collision reader remains, and the
question may resolve itself).

**Done when.** Either documented as one intended parameter, or split into
two named fields with a migration note.

## C5. Drag presets, one battery driver, script cleanup — OPEN (design agreed 2026-10-14)

Design spec, agreed with the user section by section (2026-10-14). Executed
after thesis writing; code requires `[PROCEED TO IMPLEMENTATION]`; turn this
entry into a phase-by-phase implementation plan first.

**Goal.** The repo is handed over clean: a successor reproduces the thesis
laws B2 / B3 / h405 from named presets with one driver script, without the
34-module generator chain; no config default implies a drag law; study
scripts that only produced committed results are removed, with a manifest.

**Why deferred.** Touches the code paths that produced thesis numbers;
thesis first. Every phase below is gated by oracles against the committed
run artifacts, so no thesis number can move.

**Decisions (user, 2026-10-14).** Audience = handover, after submission ·
old study scripts deleted with a manifest (git history = archive) · driver
does run + pool + score (figures stay in the existing plot scripts) · drag
off by default · design A (explicit preset functions + a package battery
module + a thin driver; rejected: committed cfg files as presets — opaque,
drift with new fields; wrapping the generators — cleans nothing) · run-dir
naming lives in the driver, presets are pure physics · the twin moves into
the package (own phase) · Tier-0 extraction: keep only what the thesis cites.

**State found (2026-10-14, read-only).** B2 is built by
`gen_tier2_detfix_postpickup.py --b2` through 34 script modules (g414 →
g413 → atlas_g4finals → …), but differs from `single_pulse_N2000` in only
35 of 121 `SimConfig` fields; its six members differ only in `seed`.
B3 differs from B2 in 3 fields (E_bind, coefficients, f_int); h405 in 9
(sampler, law, gate, bundle, E_bind, τ, f_int). The h405 run dirs are gone;
the G8 `…g8bareh405s1_th500` cfg equals the h405 generator's s1 cfg except
`initial_shell_model` (the committed h405 oracle). The per-member run exists
in three near-identical copies (`run_one` ×2, `run_lq`). 47 tests import
scripts; 4 of them test package physics through `tier2_common`
(`test_detection_stage`, `test_evaporation_shed_convention`,
`test_initial_shell_model`, `test_internal_energy_partition_law`).

### C5.1 Presets (`presets.py`, section "Drag-model production presets")

```python
DRAG_BATTERY_SEEDS = {"s1": 20260730, "s2": 20260731, "s3": 20260732,
                      "s4": 20260733, "s5": 20260734, "r6": 20260729}
def drag_production_B2(member="s1", **overrides) -> SimConfig
def drag_effective_B3(member="s1", **overrides) -> SimConfig
def drag_contrast_h405(member="s1", **overrides) -> SimConfig
DRAG_PRESETS = {"B2": ..., "B3": ..., "h405": ...}
```

- Private `_drag_battery_base(member)` = `single_pulse_N2000()` + the 26
  shared fields, grouped and commented: source / geometry (R₀ 2.666,
  droplet distribution, N 1000, t_ion 500 ps), mass model (biphasic,
  n₀ 21, pickup 0.9/ps, RRK, tabulated ladder, co-moving shed), cooling
  (density-scaled gate, retained fraction 0.1), detection (stage on,
  8.53 µs, `exclude_all_coupled`, partner-aware closure),
  `allow_inconsistent_mass_pairing=True` (comment → §6.6 / G11.4 defence).
  Also the E2-era fields the runs carry (`relaxation_*`, `v_limit` 58) —
  needed for the oracle; C1 / C4 remove them later.
- Per preset only its law:

  | | B2 | B3 | h405 |
  |---|---|---|---|
  | sampler | `post_pickup` | `post_pickup` | `raw` |
  | law | lq a 18.007, c 1.6 | lq a 19.8, c 1.6 | capped cubic b 2.515, v_c 5.5 |
  | gate | erf 3.5 Å | erf 3.5 Å | density-proportional 14.2 Å |
  | E_bind [eV] | 0.11727 | 0.0482 | 0.11676 |
  | τ [ps] | 6.4 | 6.4 | 4.4 |
  | E₀ [eV] | 0.3297 | 0.3389 | 0.405 |

- E₀ written in eV; `internal_energy_partition_fraction = E₀ /
  coulomb_available_eV` derived and commented. Coefficients from the drag
  bundle catalog (B2 = `srho3p5A18A_linear_quadratic`; B3 overrides a and
  E_bind and states the departure in its docstring); ladder rungs from the
  ladder source, not literals. Unknown member → `KeyError`; `**overrides`
  as in the legacy presets.
- **Oracle tests (exact):** B2 every member = committed
  `g5lqppc160e3297{m}` cfg.json, every field; B3 likewise
  (`g5lqppc160a198w0482e3389{m}`); h405 s1 = the G8 cfg with
  `initial_shell_model` reset to `density_tied`.

### C5.2 Battery module + driver

`i2_helium_md/simulation/battery.py`:

```python
T_HANDOVER_PS = 500.0          # Stage I end; E2 skipped (TIER2_DetectorStageFix §5.7)
MAX_BYTES_ION, STORAGE_STRIDE  # moved from gen_tier2_detfix_battery
def run_member(cfg, run_dir, *, keep_ion=False) -> Path
def run_battery(preset, members, out_root, tag, *, concurrency=2,
                overrides=None) -> list[Path]
def pool_battery(run_dirs, out_dir) -> Path
```

- `run_member` = the one consolidated copy: handover-alignment check → save
  cfg → neutral → ion to 500 ps → detection (partner-aware closure) → drop
  `ion.npz` with the existing marker unless `keep_ion`. Resumable (complete
  dir skipped); an existing cfg.json that differs → error with the field
  diff (never overwritten).
- `run_battery`: ≤ 2 concurrent by default (process pool); first failure
  stops the batch.
- `pool_battery`: the pair-preserving block layout + pairing spot-check
  from `build_pooled_detection_container.py`; `README_POOLED.txt` with all
  member seeds.
- Naming: `{case}_drag_{bundle}_N{n}_{tag}{member}_th500`, pool
  `…N6000_{tag}pooled_th500`. Default tag per preset = today's
  (B2 `detfix_conf270_g5lqppc160e3297`, B3 `…g5lqppc160a198w0482e3389`,
  h405 `detfix_conf270_h405`) → existing dirs are recognised, nothing
  re-runs.
- Scoring: `score_pool(run_dirs) -> dict` in
  `postprocess/tier2_confirmation.py` (today `_row` in a post-processing
  script): n̄, n₁, W₁, midHot, deepKE, KE₁, KE₂, n ≥ 15, S, χ²_med,
  trapped, suppressed, size columns; one-row CSV.
- `scripts/run_single_pulse_drag.py` (style of `run_single_pulse.py`, no
  CLI, USER SETTINGS):

  ```python
  PRESET = "B2"                  # "B2" | "B3" | "h405"
  MEMBERS = "all"                # or ["s1"]
  RUN_SIZE = "production"        # "smoke": N = 20, short t_ion, no pool / score
  CONCURRENCY = 2
  OUTPUT_ROOT = PROJECT_ROOT / "data" / "runs"
  TAG = None                     # None -> the preset's committed tag
  OVERRIDES = {}                 # variants, e.g. {"initial_shell_model": "bare"} (G8)
  ```

  Flow: run → pool (only when all members complete) → score → summary table
  + CSV path. Errors: unknown preset / member → `KeyError` with choices;
  mismatched cfg → field diff; pool with missing members → names them.
  `OVERRIDES` + `TAG` replace the one-off probe generators.

### C5.3 Config defaults + guard

- Principle: no default implies a drag law; `SimConfig()` and the legacy
  presets stay the hard-sphere MATLAB port.
- To `None`: `drag_form` (was `"linear_cubic"`), `drag_bundle`
  (`"shared_pure_cubic"`), `drag_spatial_gate` (`"density_proportional"`),
  `drag_gate_steepness` (14.2). Unchanged: `mass_scenario="fixed"`,
  `allow_*`, `m_eff_amu` / `mass_initial_amu`, Tier-1a fields.
- Guard (config load), all-or-nothing: coefficients set → all four set +
  existing cross-checks (form, bundle stamp, gate pairing); coefficients
  `None` → all four `None`, else error naming the fields.
- `cfg.json` migration: coefficients null **and** the four fields at exactly
  the old defaults → set to `None`, noted in the cfg provenance; any other
  combination → error. Drag cfgs load unchanged. Checkpoint schema, RNG
  draw order, constants untouched (user approval for the default-scope
  change: 2026-10-14).
- Not here: bundle-catalog pruning (reference data stays), the Tier-1a
  `anchored_discrete` code path, `run_single_pulse.py`.

### C5.4 Cleanup + manifest

Rule: a script survives only if it is a general tool or regenerates a
thesis number / figure from data still on disk; committed results under
`data/` always stay.

- **Keep:** `run_single_pulse.py`; `plot_run_summary.py`,
  `plot_detection_summary.py`, `summary_sections.py`,
  `plot_experimental_comparison.py`, `plot_hedft_comparison.py`,
  `plot_histogram_comparison.py`, `plot_speed_distribution_comparison.py`,
  `plot_mean_kinetic_comparison.py`, `plot_ion_energy_balance.py`,
  `plot_neutral_energy_balance.py`, `plot_ion_temperature_diagnostic.py`,
  `plot_ihe_ked_reference_mean_energy.py`, `plot_paper_cov.py`,
  `plot_paper_v2–v4.py`; Tier 0: `gen_tier0_runs.py`, `tier0_common.py`,
  `tier0_thesis_figure.py`, `tier0_drag_comparison.py`, extraction
  `method_b_g412_ridge.py` + `form_phase_common.py` +
  `method_b_srho3p5_refit.py`, `method_b_18A_only_refit.py` (G4.9),
  `method_a_g411_cross_case.py` (G4.11 caveat).
- **Rewire:** `tier2_g6_chord_map.py` → presets + battery naming; the 4
  physics tests → configs from `drag_production_B2()` / a test fixture.
- **Delete** (≈ 85 files + ≈ 40 script-level tests): all `gen_tier2*` /
  `gen_tier2atlas_*` (28), `gen_tier1a*`, `gen_detstagefix_validation.py`,
  `tier1a_common.py`, `tier2_common.py`, `build_pooled_detection_container.py`
  (absorbed), every `tier2*` / `tier2atlas_*` / `tier1a_*` table / report
  script except G6, the other 11 extraction scripts, `plot_shell_schedule.py`.
- **Manifest** `docs/SCRIPT_REMOVAL_MANIFEST.md`: per deleted script →
  study → doc section → committed results → last commit containing it;
  header: `git checkout <sha> -- <path>`. Deleted tests listed with the
  reason "tested the deleted script". Historical docs keep their
  references (the manifest resolves them); `CLAUDE.md`, `TIER0_SCRIPTS.md`
  and the THESIS_DISCUSSION pointers are updated.

### C5.5 The twin into the package (own phase)

`scripts/tier2_h2b_forward_model.py` (5593 lines) = core forward model
(1-D chord ODE + fate map) plus study stages (`stage_scan`,
`stage_levers`, …) and re-implemented physics (own drag-tail γ, ladder sums,
KE). Move only the core into the package (placement decided in the plan:
`simulation/` or `postprocess/`), make it call package physics (no
duplicate physics, rule 1), delete the stages with the manifest; its 5
tests move with it. Oracle: the committed G4.12 twin outputs reproduce.

### C5.6 Phases, gates, documentation

| phase | content | gate |
|---|---|---|
| P1 | presets + seeds (C5.1) | cfg oracles B2 / B3 all members, h405 via G8 cfg |
| P2 | battery module + `score_pool` + driver (C5.2) | rescoring the B2 / B3 / A5 member dirs reproduces the committed `g5_winners_table` rows to 4 dp; pooled B2 container rebuilt array-for-array; smoke battery end to end |
| P3 | defaults → `None`, guard, migration (C5.3) | legacy-config snapshot unchanged; half-set configs rejected; every committed cfg.json under `data/runs` loads |
| P4 | cleanup + manifest + rewired tests + doc pointers (C5.4) | full `pytest` green; G6 regenerates its CSVs identically |
| P5 | twin core into the package (C5.5) | committed G4.12 twin outputs reproduce |
| — | optional, once | re-run B2 s1 → `detection.npz` identical (~1 MD run) |

One commit per phase, each green before the next. Tests: no production-size
runs (smoke N = 20, short t_ion); oracles against committed run dirs skip
cleanly when the dirs are absent; tolerances exact or 4 dp with the reason
stated. Docs: CLAUDE.md (presets, driver, battery module, drag-off rule,
manifest), a how-to header in the driver, one log entry per phase.

**Interplay.** Do C5 P1–P4 **before C1**: deleting the generators removes
most of C1's per-file decisions (its list of 20 historical generators and
their tests). The presets carry the E2-era fields (`relaxation_*`,
`v_limit` 58) for the oracle; C1 and C4 then remove or split them in the
presets as well.

**Depends on.** Thesis writing finished; `[PROCEED TO IMPLEMENTATION]`;
an implementation plan written from this entry.

**Done when.** P1–P5 gates passed; `scripts/` holds only the keep list +
the driver; `SimConfig()` has no drag default; manifest committed;
CLAUDE.md updated.

### C5.7 Amendment after G13 (2026-10-27) — the design above predates the release

The spec above was written on the `never` model. Since G13
(`THESIS_FINAL_GOALS.md` G13.5–G13.7, 2026-10-26) the thesis numbers come
from different run directories and generators. Fold the points below into
the implementation plan before P1.

**How the thesis runs are produced today** (three-script chain; each
imports from the one above it):

| script | origin | role today |
|---|---|---|
| `gen_tier2_detfix_battery.py` | 2026-09-30 | h405 `raw` / 14.2 Å battery (history). Still imported for the shared basics: `MEMBER_SEEDS`, `T_HANDOVER_PS`, `MAX_BYTES_ION`, `verify_handover_alignment`. |
| `gen_tier2_detfix_postpickup.py` | 2026-10-02 | Default `--cell h405`: sampler A/B (history). `--cell lq --b1/--b1tau/--b1c/--b1w`: G5 exploration (history). **`--cell lq --b2` / `--b3`: the thesis runs** (`build_thesis_lq` + `run_thesis_member`, `THESIS_GATE_RELEASE = "decoupled"`, `ion.npz` kept) → `data/runs/B2_N1000_s1…s6`, `B3_…`. |
| `gen_tier2_g13_decoupled.py` | 2026-10-26 | One-off G13 driver. `--oracle` (B2 s1 under `never` = committed `detection.npz`), `--run` (calls the builders above; same dirs), **`--pool` → `{B2,B3,h405}_N6000_pooled` + `g13_decoupled_table.csv`**. |

h405 thesis runs: `gen_tier2_h405_pp35.py --scan-pool` → `h405_N1000_s…`.
On 2026-10-26 the B2 / B3 dirs were written through `gen_tier2_g13_decoupled.py
--run`; `gen_tier2_detfix_postpickup.py --cell lq --b2` builds the identical
cfg (an existing dir with a different `cfg.json` is refused).

**Changes to the C5 design:**

- **Presets (C5.1):** set `self_bound_gate_release="decoupled"` in the
  three drag presets (the `SimConfig` default stays `never`). The h405
  preset is the **re-landed** contrast law (`post_pickup`, erf 3.5 Å,
  capped cubic b 2.515 / v_c 5.40, E_bind 0.1168, τ 4.4, E₀ 0.424;
  `THESIS_DISCUSSION.md` §15), not the `raw` / 14.2 Å row in the C5.1
  table. Oracles: every member = the `cfg.json` of `{B2,B3,h405}_N1000_s1…s6`
  (s6 = seed of r6); the `never` batteries stay as history only.
- **Battery / driver (C5.2):** naming = the plain thesis dirs
  (`{law}_N1000_s{k}`, pools `{law}_N6000_pooled`), not the long
  `…_th500` tags; `keep_ion=True` is the thesis default (the shed-timing
  and exit reads need `ion.npz`). Score oracle = the decoupled rows of
  `g13_decoupled_table.csv` (4 dp), not `g5_winners_table`.
- **Cleanup (C5.4):** add `gen_tier2_g13_decoupled.py` and
  `gen_tier2_h405_pp35.py` to the delete list (absorbed by the driver;
  the G13 `--oracle` check becomes a driver test or a manifest note).
  Plot scripts `tier2_b2_thesis_figure.py` and `tier2_g6_chord_map.py`
  (both read the decoupled dirs since 2026-10-26/27) go on the keep list,
  rewired to the driver's naming.
- **Shared basics first:** move `MEMBER_SEEDS`, `T_HANDOVER_PS`,
  `MAX_BYTES_ION` / `STORAGE_STRIDE` and the handover-alignment check out of
  `gen_tier2_detfix_battery.py` into the battery module **before** deleting
  any generator; the thesis path imports them today.
- Scale for the manifest: 35 `scripts/gen_*.py` files as of 2026-10-27.
