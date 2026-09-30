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

**Depends on.** Thesis final; production pipeline switched to the skip path
(`TIER2_DetectorStageFix.md` §3c item 3); decision on whether any E2-era
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
