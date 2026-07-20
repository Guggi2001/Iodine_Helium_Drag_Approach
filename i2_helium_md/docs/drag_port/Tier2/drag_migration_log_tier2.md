# Drag-Model Port — Migration Log: Tier 2

**Purpose.** This log is the **decision + delivery history for Tier 2** (generative
mass dynamics + size-distribution comparison). Tier-0 history lives in
`docs/drag_port/Tier0/drag_migration_log_tier0.md`; Tier-1a history in
`docs/drag_port/Tier1/drag_migration_log_tier1a.md`.

**Companion docs:** `TIER2_IMPLEMENTATION_PLAN.md` (the active build plan),
`DRAG_PORT_DESIGN_DECISIONS.md` (frozen design),
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (the live mass-model detail),
`CALIBRATION_MAP.md` (parameter classes + Tier anchors).

---

## Tier-2 planning — decision record (2026-06-29)

Study of `CLAUDE.md`, `DRAG_PORT_DESIGN_DECISIONS.md`,
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`, `CALIBRATION_MAP.md`, and the
Tier-1a plan/log, plus a code-surface audit. No code written (the
`[PROCEED TO IMPLEMENTATION]` boundary holds); these are documentation + scoping
decisions for the eventual build.

- **Tier-2 broken into 11 slices across 6 phases** (`TIER2_IMPLEMENTATION_PLAN.md`):
  Phase A energetics primitives (L ladder, K cooling, U E_int budget); Phase B
  stochastic channels (ρ density, P pickup, Q evaporation); Phase C schema +
  generative integrator (X v6→v7 + 5-term invariant, G generative driver); Phase D
  bridge (Z generative-vs-anchored); Phase E observable (H extraction, W Wasserstein,
  D2 diagnostics); Phase F campaign (R2). The 5-term invariant is the cross-cutting
  correctness gate from X onward.

- **Three scope decisions (user, 2026-06-29):** (1) **full program** —
  mechanism + comparison + calibration; (2) **validation-first** — first generative
  runs at the 0.80 eV (d=9 Å) budget cross-checked against the TDDFT 21→19→14 decline
  and the Tier-1a anchored results, then the 2.70 eV production onset; (3) **include
  the generative-vs-anchored bridge slice** (Z) to de-risk the mechanism before the
  experimental histogram.

- **Tier-2 reference data confirmed present** (no export slice needed):
  `data/reference/integrated_i_he_abundance.csv` (experimental I⁺Heₙ size
  distribution, the arbiter), `data/reference/vmi_summary/vmi_iplus_*.csv`
  (per-fragment velocities), `data/reference/9A_All_Data.csv` (TDDFT cross-check).

- **Reuse confirmed from Tier 1a:** `physics/mass_jump.py` resets, `physics/baoab.py`
  SQ1 O-step, the `shed_step` jump-then-O seam, the v6 checkpoint + v5 shim, the
  4-term `ion_ledger_closure`, `velocity_distribution.py`, and `tier1a_common`/
  `gen_tier1a_runs` orchestration. Tier 2 makes sheds generative (adds pickup + RRK
  evaporation + E_int reservoir + cooling), extends to the 5-term invariant, and adds
  the size-distribution comparison layer.

- **Naming item flagged for build:** MASS §11 names the production value
  `biphasic_energy_gated`; `config.py` currently has `biphasic`. Recommend keeping
  `biphasic` and aliasing in docs; resolve at the Phase-C build.

- **Open: OQ1 (drag electronic-state provenance)** stays open → carry both electronic
  pictures as the co-fit free knob; do not cite E_bind as independent support for the
  mixture. User is in author contact; update on resolution.

---

## Phase A detailed plan — delivered (2026-06-29)

The rough Phase A sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 was expanded into a
dedicated, Tier-1a-depth plan: **`TIER2_PHASE_A_IMPLEMENTATION_PLAN.md`**. No code
(the `[PROCEED TO IMPLEMENTATION]` boundary holds); planning doc only.

- **Module layout decided (user, 2026-06-29):** three separate one-concern pure
  modules — `physics/dissociation_ladder.py` (L), `physics/internal_energy_cooling.py`
  (K), `physics/internal_energy_budget.py` (U) — mirroring the Tier-1a
  `shell_schedule.py` / `mass_jump.py` precedent and the L/K/U independence in the
  dependency graph. New sourced constants land in `constants.py`; new knobs in
  `config.py`.
- **Content:** per-slice interface contracts, LaTeX-with-dimensional-checks encoded
  forms (Form U ladder + gate; Newton cooling + occupancy-resolved `E_∞(N)` +
  pair/electrostriction split; S1/S2/K1 ΔE_int rules + post-t× reconstruction), a
  consolidated golden-oracle fixture table, the config field map (MASS §11 names),
  intra-phase build order (L → K → U with mocks), test methodology, and the
  out-of-scope guard (no RNG / integrator / E_int state / 5-term closure / noise —
  those are Phase B/C/Tier 3).
- **Golden values verified by direct conversion:** `D_0(1)` 106.9 cm⁻¹ = 0.013254 eV
  (X₂) / 74.4 cm⁻¹ = 0.009224 eV (mix); `D_floor` 4.97 cm⁻¹; `Σ(21)` X₂ 0.25–0.28 eV /
  mix 0.17–0.19 eV; `|S|/n*` = 118 cm⁻¹ @ n*=21; f_int floor @ 0.80 eV X₂ 0.31–0.35 /
  mix 0.21–0.24, @ 2.70 eV X₂ 0.09–0.10 / mix ~0.065.

---

## Phase B detailed plan — delivered (2026-06-29)

The rough Phase B sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 (the two stochastic
channels) was expanded into a dedicated, Tier-1a/Phase-A-depth plan:
**`TIER2_PHASE_B_IMPLEMENTATION_PLAN.md`**. No code (the `[PROCEED TO IMPLEMENTATION]`
boundary holds); planning doc only.

- **Two design decisions locked (user, 2026-06-29):**
  1. **Module layout — three one-concern modules + capture in `mass_jump.py`:**
     `physics/helium_density.py` (ρ), `physics/pickup.py` (P), `physics/evaporation.py`
     (Q); the +He **capture reset** extends `physics/mass_jump.py` symmetric to
     `cold_shed`/`continuous_velocity_shed`, keeping the reduced-mass defect formula in
     one place. Mirrors the Phase-A 3-module split and the ρ/P/Q dependency-graph
     independence.
  2. **ρ sourcing — reuse the erf-complement gate:** Slice ρ uses the same
     `0.5·(1−erf(depth/steepness))` form as `drag.py::spatial_gate` (steepness 14.2 Å) as
     the `ρ_He/ρ_bulk` ratio; a **sourced baseline/TDDFT density profile is the declared
     fallback** (rule-2 convention), not built in Phase B. The drag spatial gate stays
     **G2** (no G2→G4 promotion this phase; CALIBRATION_MAP rows 5/8).
- **Content:** per-slice interface contracts (ρ/P/Q + the `mass_jump.capture`
  counterpart), LaTeX-with-dimensional-checks for the density gate, pickup
  (`P_attach`, `λ_attach` with Langmuir cap, momentum-conserving capture reset + defect,
  S1 heat) and evaporation (self-bound gate, saturating RRK rate with `n=1` direct
  dissociation, cold-shed reuse, K1 drain), the consolidated golden-oracle fixture table,
  the config field map (MASS §11 names + the `mass_rate_*`→`pickup_*` rename flagged as a
  build-time item), the intra-phase build order (ρ → {P, Q}, P/Q parallelizable), the
  seeded/mock-RNG test methodology, and the out-of-scope guard (no integrator / O-step /
  `E_int` state / schema bump / 5-term closure / noise / drag-law / fitting — those are
  Phase C / Phase F / Tier 3).
- **Reuse confirmed:** `mass_jump.py` reset spine (`cold_shed`,
  `_reduced_mass_defect_coeff`, `kick_factor`, component variants),
  `drag.py::spatial_gate` (the erf-complement gate, shared via a single helper),
  Phase-A L/U primitives (mocked at unit level), and `config.py` `check_drag_config`
  validator scaffolding. New physics is only the `mass_jump.capture` counterpart and the
  two channel modules.

---

## Phase C detailed plan — delivered (2026-06-29)

The rough Phase C sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 (Slice X schema + 5-term
invariant; Slice G generative driver) was expanded into a dedicated, Phase-A/B-depth
plan: **`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md`**. No code (the
`[PROCEED TO IMPLEMENTATION]` boundary holds); planning doc only.

- **Two design decisions locked (user, 2026-06-29):**
  1. **Slice G placement — extend `simulation/ion_propagation_step.py`** with a
     `biphasic_step` seam alongside `shed_step` / `baoab_propagation_step`, reusing the
     `IonStepState` dataclass (extended with `E_int_eV`) and the SQ1 O-step in place.
     Mirrors the Tier-1a `shed_step` precedent.
  2. **Production scenario name — keep `biphasic`.** The `MassScenario` literal stays
     `{fixed, biphasic, anchored_discrete}`; `biphasic` *is* the production value (full
     name `biphasic_energy_gated`, MASS §11, documentation only). No config-enum rename,
     no §6.5 guard-set churn.
- **Content:** the two slices with full interface contracts —
  **X** (checkpoint **v6→v7** adding `E_int_eV (2N,T)` + a v6→v7 back-compat shim
  mirroring the delivered v5→v6 shim; `IonStepState.E_int_eV` field; the **5-term
  invariant** `E_kin+E_pot+E_dissip+E_mass_transfer+E_int≈const` extended *in place* in
  `energy_balance.ion_energy_totals`/`ion_ledger_closure`; the **RNG draw-order lock**),
  and **G** (the `biphasic_step` per-step loop: conservative kicks → K2 cooling → at most
  one mass event, shed-then-pickup if both → jump-then-O with the SQ1 O-step rebuilt at
  `m⁺` → S1/S2/K1/K2 `E_int` update). Plus the §6 per-channel closure table, the
  reduced-mass capture defect, the R4 double-count guard, the config contract (integrator
  flags + the §6.5 `time_resolved`/`allow_inconsistent_mass_pairing` path), build order
  (X before G), analytic-limit test methodology, and the out-of-scope guard (no new
  channel physics / noise / drag-law change / draw-order change / D-E-F steps).
- **Reuse confirmed (code audited):** `checkpoint.py` v6 + the v5→v6
  `_migrate_ion_checkpoint` shim + `_validate_against_cfg` `trajectory_2N_T_fields`;
  `energy_balance.py` `ion_energy_totals` / `ion_ledger_closure` (the 4-term gate, with
  `E_mass_transfer` already added as an optional term — the pattern X follows for
  `E_int`); `ion_propagation_step.py` `IonStepState`, the `shed_step` SQ2/SQ3 seam, and
  `baoab_propagation_step`; `physics/baoab.py` SQ1 O-step; the §6.5 guard
  (`_EVOLVING_MASS_SCENARIOS` → `time_resolved` requirement → `allow_inconsistent_mass_
  pairing` override). G adds **no new physics** — it composes accepted A/B modules + the
  Tier-1a seam.
- **Schema version reconcile noted:** MASS §7 text says "v6" (predates the Tier-1a v6 use
  of the slot); the live schema is already v6, so Tier 2 bumps **v6→v7** for `E_int_eV` —
  the v6→v7 shim is the reconcile, recorded so the doc/code version story stays straight.

---

## Phase D detailed plan — delivered (2026-06-29)

The rough Phase D sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 (Slice Z, the
generative-vs-anchored bridge) was discussed and expanded into a dedicated plan:
**`TIER2_PHASE_D_IMPLEMENTATION_PLAN.md`**. No code (the `[PROCEED TO IMPLEMENTATION]`
boundary holds); planning doc only. Brainstormed live with the user before writing.

- **Three design decisions locked (user, 2026-06-29):**
  1. **Minimalistic existence demonstration** — Slice Z runs a **single representative
     priored parameter point** (λ₀, f_int, f_ret, τ chosen in-band, *not* fitted) to show a
     priored point *exists* that reproduces the 9 Å kinematics. **Not** a sweep, **not** a
     mean-field ODE layer, **not** a calibration (those stay in Phase F).
  2. **Minimal reconstruction self-contained in D** — the bridge builds only the thin
     pieces it needs (mean `n(t)` = average of the v7 `n_shell` array over the run's ions;
     tiny `t×` = first `E_int<Σ(n)` and `Π(t)=λ·f_ret·τ` helpers; reuse
     `compare_distance`/`compare_velocity_magnitude` for `R(t)`/`|v(t)|`). **Phase-E
     Slices H/W/D2 are NOT pulled forward.**
  3. **Targets loaded, not regenerated** — read the delivered Tier-1a anchored `n(t)`
     artifact + the frozen 9 Å TDDFT curve (`data/reference/9A_All_Data.csv`) as fixed
     comparison targets (no anchored-artifact churn). **Ensemble = the ions in one
     generative run** (no multi-run seeding).
- **Nature:** Phase D builds **no new physics and no new modules** — it is a **run +
  compare + report** slice composing the delivered Phase-C `biphasic` driver, the Tier-1a
  `tier1a_common`/`gen_tier1a_runs` scaffolding, the `compare_*` helpers, and the 5-term
  `ion_ledger_closure`, plus two tiny `t×`/`Π` reconstruction helpers (placed where
  Phase-E D2 later generalizes them).
- **Content:** the run (9 Å / 0.80 eV `biphasic`, scenario-stamped → §6.5 pairing guard →
  `allow_inconsistent_mass_pairing=True`), the loaded kinematic targets (anchored
  21→19→14 staircase, 9 Å TDDFT `R(t)`/`|v(t)|`, GAH25 `t×` 5–6.5 ps ±factor-2 prior, Π
  regime), the two reconstruction helpers + their unit tests, the few-step driver smoke
  (5-term invariant), and the **written diagnostic report** with the actionable-lever
  framing for any miss. **Reported, not auto-adjudicated** (mirrors the Tier-1a reporting
  stance; the bridge *localizes* a bug, it does not arbitrate fidelity — experiment does
  that at Phase E/F).

---

## Phase E detailed plan — delivered (2026-06-29)

The rough Phase E sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 (Slices H/W/D2 — the
Tier-2 observable / comparison layer) was expanded into a dedicated, Phase-A/B/C/D-depth
plan: **`PHASE_E_IMPLEMENTATION_PLAN.md`**. No code (the `[PROCEED TO IMPLEMENTATION]`
boundary holds); planning doc only.

- **Three scope decisions locked (user, 2026-06-29):**
  1. **R5 truncation → build a post-ejection relaxation stage.** Terminal `n` at 20 ps is
     only an upper bound (MASS §R5); Phase E adds a **reduced relaxation driver** (Slice E2,
     `simulation/relaxation_stage.py`) that propagates the energy-gated cascade forward to
     the experimental timescale with **pickup and drag off** (the locked Q/K/U channels in
     their free-flight, evaporation-only limit — no new physics), so the size distribution
     is read at matched time.
  2. **Size distribution only.** The `vmi_summary/*.csv` references are **aggregate**
     I⁺He / I⁺-gas distributions, not per-fragment per-`n`; Phase E **drops velocity
     comparison entirely** and commits to the integer-n size distribution
     (`integrated_i_he_abundance.csv`) as the sole observable. The original Slice-H
     per-fragment velocity histograms are **cut**.
  3. **Pure functions only.** No CLI / figure / report scripts; all orchestration and
     figures deferred to the Phase F campaign (R2).
- **Re-slice H/W/D2 → E1–E5:** **E1** terminal integer-n size-distribution extractor
  (`postprocess/size_distribution.py`); **E2** relaxation stage (above); **E3** abundance
  reference loader mirroring `hedft_loader` (`postprocess/abundance_loader.py`); **E4**
  pure-numpy integer-support Wasserstein `W₁=Σ|F_sim−F_ref|`
  (`postprocess/distribution_compare.py`, activates `validation_histogram_metric`);
  **E5** derived diagnostics `t×`/`Π(t)`/regime/total-strip reachability
  (`postprocess/derived_diagnostics.py`).
- **Hard prerequisites:** Slice X (v7 `E_int_eV`) gates E2/E5; L+P+K gate E5; G+Q+K+U gate
  E2. E1/E3/E4 are buildable now on synthetic v7 checkpoints. Real runs are wired in at
  Phase F (R2) — Phase E is tested entirely on synthetic checkpoints + the real CSV.
- **Reuse confirmed (code audited):** `shell_schedule.complex_mass_amu` (E1 rung labels),
  `RunDirectory`/`IonCheckpoint.n_shell`/`mass_final_kg` (E1), the Slice-G `biphasic_step`
  loop + `shed_step` seam + `mass_jump.cold_shed` (E2, pickup/drag disabled),
  `hedft_loader.load_hedft_trajectory` contract (E3 mirror),
  `config.validation_histogram_metric` (E4), `energy_balance.ion_ledger_closure` 5-term
  (E2/E5 invariant context). New physics is **none** — E2 only composes accepted channels.
- **Cross-cutting gates:** the 5-term invariant is E2's correctness detector (pickup off,
  drag frozen → `E_dissip` constant); integer-n support + pure-numpy Wasserstein
  throughout; the R5 matched-time/upper-bound caveat is carried in every E4 record; the
  locked RNG draw order is *extended* (shed-only) for the relaxation stage, default
  `relaxation_stage_enabled=False` so default scope is unchanged.

---

## Phase F detailed plan — delivered (2026-06-29)

The rough Phase F sketch in `TIER2_IMPLEMENTATION_PLAN.md` §4 (Slice R2 — the
calibration campaign + production switch) was expanded into a dedicated plan:
**`PHASE_F_IMPLEMENTATION_PLAN.md`**. No code (the `[PROCEED TO IMPLEMENTATION]`
boundary holds); planning doc only. Phase F is the stage that **absorbs everything
the earlier phases postponed** — the κ+picture co-fit, the f_int/f_ret/τ bounded
scalars, the 0.80→2.70 eV production switch, identifiability + regime reporting, the
Calvo24 total-strip runs, and all CLI/scoreboard/figure scripts Phase E deferred.

- **Three scope decisions locked (user, 2026-06-29):**
  1. **Staged (clean-fallback) calibration** (MASS §6.11): co-fit only the two Free
     knobs (κ + electronic picture); pin `f_int` at a fixed `τ`; then a
     `τ`-insensitivity sweep; `f_ret` via the 9/18 Å density contrast. The full
     κ×picture×f_int×f_ret×τ grid is a **documented escalation** only.
  2. **Larger-N, single-seed runs** (N≈500 → ~1000 fragments/run) for a stable
     size-distribution histogram / Wasserstein score.
  3. **Build the total-strip secondary-run path now**; **R6** 9-Å re-extraction and
     the **p↔κ** cap split stay **documented conditional triggers** (built on demand).
- **Re-slice R2 → F1–F6:** **F1** campaign harness (`scripts/tier2_common.py`,
  `build_biphasic_cfg` + run-dir convention mirroring `tier1a_common`); **F2** staged
  run-matrix generator (`scripts/gen_tier2_runs.py`, 0.80 eV first, 9 Å + 18 Å,
  ion→relaxation pipeline); **F3** scoreboard
  (`scripts/post_processing/tier2_size_distribution_table.py`, composes E1/E2/E4/E5 →
  per-run W₁ matched + upper-bound + diagnostics, mirrors `tier1a_rmse_table.py`);
  **F4** identifiability + regime report (`tier2_identifiability_report.py`); **F5**
  2.70 eV production switch (gated on 0.80 eV landing) + total-strip secondary runs;
  **F6** overlay figures (`SHOW_FIGURE`-gated).
- **Nature:** Phase F builds **no new physics and no new `SimConfig` physics fields**
  — it composes only **accepted A–E modules** + the Phase-D bridge gate, and only
  calibration *values* move via config knobs. New surfaces are orchestration scripts +
  report assemblers + figures (plus a `total_strip` run-tag variant).
- **Reuse confirmed (code audited):** `tier1a_common.build_anchored_cfg` /
  `tier1a_run_tag` / `tier1a_run_dir_name` (F1 mirror), `tier0_common.build_drag_cfg` /
  `run_dir_name`, `complex_mass_amu(21)`; `gen_tier1a_runs._run_one` neutral→ion
  pipeline (F2); `tier1a_rmse_table` `collect_*`/`score_*`/`format_table`/
  `write_rows_csv` + figure-builder pattern (F3/F6); E1–E5 + `load_he_abundance_
  reference` + 5-term `ion_ledger_closure` (F3/F4).
- **Cross-cutting gates:** validation-first (0.80 eV, Phase-D-bridge-cross-checked,
  before the 2.70 eV switch); **reported, not auto-adjudicated** (no code asserts a
  fidelity verdict); scenario-stamped `coulomb_available_eV` guards the budget↔drag
  pairing with `allow_inconsistent_mass_pairing=True` the normal production path (R6,
  §6.6); 5-term invariant residual carried per run; out of scope = noise/Tier 3
  (inert), Tier-0 drag law / neutral propagation / RNG draw order / mechanism, and the
  R6/p↔κ conditional triggers (document-only). Identifiability of 8+ quantities on one
  observable is the central reported dependency (CALIBRATION_MAP "Anchor coverage").

---

## Phase A — three open questions resolved (2026-06-30)

A review pass over `TIER2_PHASE_A_IMPLEMENTATION_PLAN.md`, cross-referenced against
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` and `CALIBRATION_MAP.md`, surfaced
five open questions; three needed a user call and are resolved (decision owner: user).
No code — the `[PROCEED TO IMPLEMENTATION]` boundary still holds; docs-only amendment.

- **1. `cooling_relaxed` ladder picture → rule-2 declared-but-unread stub.** The source
  leaves the third electronic-picture first rung **unpinned** (MASS rev-2026-06-21 #3:
  "e.g. a relaxation-weighted blend", *strictly between* mix 74.4 and X₂ 106.9 cm⁻¹;
  CALIBRATION row 20: "one selection, between mix & X₂"). Phase A ships the enum arm
  declared-but-unread, asserts **ordering only** (`statistical_mixture < cooling_relaxed
  < x2_only`), no 4-figure oracle; the concrete blend is pinned at Phase F (removed from
  the rule-2 exception table then). Defer-as-stub is faithful to source, not a shortcut.
- **2. Slice K module renamed `internal_energy_cooling.py` → `solvation_cooling.py`.**
  K2 cools the GAH25 variable **`E_solv.struct = E_bind + E_int`** (MASS lines
  334/837/848); `E_int` is only one additive term *inside* it, owned by Slice U
  (`internal_energy_budget.py`). The module is now named for the variable it relaxes
  (Quality Principle 6). The config knob name `internal_energy_cooling_tau_ps` is the
  sourced MASS §11 name and is **unchanged**. Sibling-doc filename references updated for
  drift (PHASE_C `newton_cool_step` path, PHASE_F Phase-A module list); the historical
  Phase-A-delivered entry above (2026-06-29) predates the rename.
- **3. Full K/U surface built in Phase A.** `newton_cool_step` (K) and
  `reconstruct_e_int_eV` (U) are built and oracle-tested in Phase A even though their
  first *composing* caller is the Phase-C driver (G) — keeps Phase A self-contained.

Two further points flagged for a clarifying comment only (no fork):
- **`E_int^eq = 0` is a split-consistency (tautology) check, not independent
  corroboration** — with `E_int` defined as the residual after subtracting the pair +
  electrostriction split, it is identically 0 at equilibrium by construction (MASS line
  881). Annotated in the Slice K test spec.
- **No default `f_int` lands in Phase A** — only the `f_int_floor` helper is exercised;
  the committed default is a Phase-F calibration output (the floor is scenario-split,
  0.21–0.35 @ 0.80 eV vs 0.065–0.10 @ 2.70 eV, so no single default serves both budgets).

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. Knob classes
match (κ + picture **Free** rows 19/20; τ/f_int/f_ret **Bounded** rows 11/13/14;
D₀(1)/D_floor/n\*/|S| **Sourced** rows 18/23/22/12; Σ(21) **Derived** row 21); config field
names match MASS §11; |S|=0.308 / drag-binding 0.1168 confirmed as "separate cross-checks,
not band ends" (row 21).

---

## Phase A — Slice L DELIVERED (2026-06-30)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN), oracle asserts against the plan §2/§4 golden values; pure,
stateless, mass-agnostic (no `γ`, no `m`, no RNG, no integrator).

**Build (4 files):**
- `physics/constants.py` — new **sourced** anchors (cm⁻¹, converted via
  `EV_PER_WAVENUMBER`): `N_STAR=21`, `D0_1_X2_WAVENUMBER=106.9`,
  `D0_1_MIX_WAVENUMBER=74.4`, `D0_1_COOLING_RELAXED_WAVENUMBER=90.65`,
  `D_FLOOR_WAVENUMBER=4.97`. Additions, **not** edits to the MD constants table.
- `physics/dissociation_ladder.py` — `sigma`, `d0_of_n`, `ladder_cumsum`,
  `gate_threshold` (alias), `first_rung_d0_eV` resolver, `tabulated_ladder` /
  `TabulatedLadder` fallback, module-level `D_FLOOR_EV` / `CLIFF_CENTER` /
  `_FIRST_RUNG_EV` (the single source of truth for picture keying).
- `config.py` — `LadderForm` / `LadderElectronicPicture` literals; live fields
  `dissociation_ladder` (`form_u`), `ladder_steepness` κ (`1.0`, Free),
  `ladder_electronic_picture` (`statistical_mixture`); `check_ladder_config` guard
  (selector reject + picture reject + `D_0(1) > D_floor`) wired into `validate()`.
- `tests/test_dissociation_ladder.py` (57) + `tests/test_ladder_config.py` (9).

**ExitPlanMode resolutions applied at build (user, 2026-06-30):**
- **(b) `cooling_relaxed` ships a provisional value = arithmetic mean** of mix and
  X₂ (90.65 cm⁻¹ = 0.011240 eV). This **refines** planning-decision #1 above: the
  arm is *not* unread — it returns the provisional blend, which the **ordering**
  oracle (`mix < cooling_relaxed < x2_only`) reads; there is no 4-figure value
  oracle for it, and the concrete relaxation-weighted blend is still pinned at
  Phase F (strict betweenness preserved; faithful to MASS rev-2026-06-21 #3).
- **τ knob name kept** `internal_energy_cooling_tau_ps` (decision #2 / Slice K).
- **`tabulated_ladder` fallback built + round-trip tested** in Phase A (not deferred).
- f_int default → `Optional[float] = None`, config-silent until Phase F — a **Slice
  U** concern, recorded here, not built at L.

**Oracle cross-validation.** First rungs match to 4 figures and are κ-independent
(Form-U normalisation pins `D_0(1)` = source exactly). `Σ(21)` lands in the plan
2-decimal bands (mix 0.17–0.19, X₂ 0.25–0.28 eV) across κ∈[0.3,5] with ~10–11%
spread (the decoupling property); the κ=5 mixture endpoint computes to **0.19306 eV
→ 0.19** and the κ=0.3 X₂ endpoint to **0.24973 → 0.25**, so band membership is
asserted at the plan's stated two-decimal precision (`round(value, 2)`), not an
ad-hoc epsilon — the rounded endpoints *are* the computed values.

**Rule-2 table.** L's three config fields are **born live** (read by the module /
guard on arrival) — nothing to remove from the Tier-0 exception table. The only
Phase-A rule-2 carry is the `cooling_relaxed` *concrete blend* (provisional-average
placeholder until Phase F), tracked in the plan §9 notes.

**Tests:** `test_dissociation_ladder.py` + `test_ladder_config.py` = 66 green; full
suite **1007 passed, 0 failed** (17 warnings = pre-existing `anchored_discrete`
mass-pairing, unrelated).

### Slice L — review + test-hardening pass (2026-06-30)

A high-recall review (8 angles, run inline — small self-authored diff) surfaced
three edge/robustness items; all fixed test-first. No silently-wrong physics in the
sourced κ∈[0.3,5] range was found (the oracles stand).

- **κ positivity guard (was missing).** `d0_of_n` divides by `(1 − σ(1))`; κ ≤ 0
  inverts/flattens the cliff → near-zero denominator → inf rungs (κ = 0 → degenerate
  flat ladder). `check_ladder_config` now refuses `ladder_steepness ≤ 0` (CLAUDE.md
  principle 4). This is also the field's **genuine config-level read** — before the
  guard, `ladder_steepness` was declared-but-unread (module took κ as a param), so
  the "LIVE at Slice L" claim is now actually true.
- **σ overflow → tanh form.** `sigma` reimplemented as the identical
  `0.5·(1 + tanh(x/2))`; saturates instead of forming an `inf` intermediate at sharp
  κ (no `RuntimeWarning`, no oracle change).
- **`ladder_cumsum` fractional-n guard.** Rejects a genuinely fractional occupancy
  loudly (was a cryptic `cum[float]` `IndexError`); integer-valued floats accepted.
- **Test-side fix:** the σ-monotonicity assertion was over-strict (σ saturates
  exactly to 0/1 in the tails → consecutive-equal); corrected to non-decreasing
  globally + strictly-increasing across the unsaturated cliff window.

**Tests after hardening:** Slice-L suites **104 green** (66 → +38); full suite
**1045 passed, 0 failed** (same 17 unrelated warnings). Next: Slice K
(`solvation_cooling.py`).

---

## Phase A — Slice K DELIVERED (2026-06-30)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (plan approved in
plan mode). TDD throughout (test→RED→GREEN); pure, stateless, mass-agnostic (no
`γ`, no `m`, no RNG, no integrator). Four open points were resolved by the user
*before* coding (see below).

**Build (4 files):**
- `physics/constants.py` — new **sourced** anchor `S_ABS_EV = 0.308` (eV-primary;
  = 2484 cm⁻¹ total, 118 cm⁻¹/atom @ n*=21; DFT first-shell solvation, MASS K2).
  Addition, **not** an edit to the MD constants table.
- `physics/solvation_cooling.py` — `s_collective_eV` (`|S(N)|=|S|·Σ(N)/Σ(n*)`),
  `e_infinity_eV` (`−|S(N)|`), `e_bind_pair_eV` (`−Σ(N)`, wraps L), `e_electrostriction_eV`
  (`−(|S(N)|−Σ(N))≤0`), `newton_cool_step` (exact `e^{−dt/τ}` relaxation, `dt`-robust,
  fail-loud on `τ≤0`). Consumes L's `ladder_cumsum`.
- `config.py` — fields `solv_struct_asymptote_eV` (`=S_ABS_EV`, single-source) and
  `internal_energy_cooling_tau_ps` (`6.55`, geometric mid of [2.6,16.5]); derived
  `electrostriction_binding_eV` **@property** (`E_elec(n*)` for inspection);
  `check_solvation_cooling_config` guard (`τ>0`, `|S|>0`) wired into `validate()`.
- `tests/test_solvation_cooling.py` + `tests/test_solvation_cooling_config.py`.

**Four open questions resolved (user, 2026-06-30, before coding):**
1. **Config surface = all three fields** (plan §5 literal): `solv_struct_asymptote_eV`,
   `electrostriction_binding_eV` (derived property), `internal_energy_cooling_tau_ps`.
2. **τ default = geometric mid 6.55 ps** (`√(2.6·16.5)`), in the R8 band.
3. **τ guard = `τ>0` only**; the [2.6,16.5] band stays a soft prior, not enforced.
4. **|S| = eV-primary `0.308`**, surfaced as a config field backed by one
   `constants.py` anchor (honors Sourced class, CALIBRATION row 12).

**Plan refinements recorded (as-built vs plan §2/§3/§5):**
- Plan §2.1 listed `S_ABS_eV` as a constant; it ships as a **config field defaulting
  to a single `constants.py` anchor** (the value lives once; field is the override/
  inspection surface). |S| stays **Sourced** — never tuned (the field is not a Free knob).
- Plan §5 filed `solv_struct_asymptote`/`electrostriction_binding` under **L's** bullet,
  but L did not land them — they are **K's** (where consumed). L delivered only its three
  ladder fields. `electrostriction_binding` ships as a **derived @property**, not a stored
  field (it is N-/picture-/κ-dependent).
- **Cold-shed neutrality clarified** to the **pair + `E_int` sub-sum** form
  (`Δ(−Σ + E_int) = +D₀ − D₀ = 0`, machine precision); the collective electrostriction
  marginal (`|S| > Σ(n*)`) is the A8 bath booking, intentionally excluded from the identity.
- **`E_int^eq = 0`** asserted but annotated as a split-consistency **tautology** (MASS
  line 881), not an independent anchor.

**Oracle cross-validation.** `s_collective_eV(n*) = |S| = 0.308 eV` exactly for any
picture/κ (ratio 1); `|S|/n* = 118.3 cm⁻¹` (3 figs); `E_∞(0) = 0` (OQ6); `E_elec ≤ 0`
everywhere; split closes `E_∞ = E_bind_pair + E_elec`; `newton_cool_step` exact decay
factor + fixed point + `dt`-robust; neutrality identity = 0 to 1e-15. Independence test
stubs `ladder_cumsum` (Σ(n)=n) to prove K composes L only through the gate.

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. |S| Sourced
(row 12); τ Bounded (row 11), default in-band; the config field name
`internal_energy_cooling_tau_ps` keeps the sourced MASS §11 name even though the module
is `solvation_cooling.py` (the 2026-06-30 rename kept the knob name).

**Rule-2 table.** K's two fields are **born live** (read by the module / guard on
arrival); nothing to remove from the exception table. The only Phase-A rule-2 carry
remains the `cooling_relaxed` *concrete blend* (Slice L; pinned at Phase F).

**Tests:** Slice-K suites **141 green**; full suite **1186 passed, 0 failed** (same 17
unrelated `anchored_discrete` warnings).

### Slice K — review + test-hardening pass (2026-06-30)

An inline review (independent numerical re-verification of every oracle on a fresh
interpreter, plus an 8-angle correctness/edge sweep) confirmed the physics and
surfaced one hardening gap + several coverage gaps; all addressed.

- **Stronger-than-planned result (locked as a regression).** `E_elec ≤ 0` is a
  **structural guarantee for all κ>0**, not just the sourced [0.3,5] band: `Σ(n*) ≤
  n*·D_0(1) = 0.278 eV (X₂) / 0.194 eV (mix)`, both always `< |S| = 0.308`. Verified to
  κ=100 and asserted via `test_sigma_at_nstar_below_S_abs_for_all_kappa` +
  `test_electrostriction_nonpositive_for_all_kappa`.
- **`newton_cool_step` `dt < 0` guard (was missing).** A negative step gives
  `exp(+|dt|/τ) > 1` — silent anti-cooling away from `E_∞`. Now fail-loud (principle 4),
  added test-first (RED→GREEN); `dt = 0` remains a valid no-op. This complements the
  existing `τ > 0` guard.
- **Coverage extended (+36 tests, regression locks on already-correct behavior):**
  picture-correct asymptote targeting, long-time convergence to `E_∞`, array-N/scalar-E
  and dt=0 edges, scalar→float / array→ndarray return discipline on the cooling step,
  `s_abs_eV` override through `e_infinity_eV`, and **fail-loud propagation** of L's
  negative/fractional-N guard through all four consumers + the cooling step.
- **Deliberate non-finding.** The plan §3 "electrostriction ≈4.7× pair-at-radius ~25
  cm⁻¹/atom" parenthetical stays **unencoded** — there is no sourced 25 cm⁻¹ anchor in the
  codebase, so asserting it would test a doc aside, not physics. The encoded electrostriction
  oracles are sign (`≤0`), the split-definition equality, and `|S|/n* = 118 cm⁻¹`.
- **Deliberate non-guard.** Module functions do **not** guard `s_abs_eV > 0` (only the
  config guard does) — the module stays config-agnostic by design (plan §3); `τ`/`dt` are
  guarded because they crash (div-by-zero) or silently invert the cooling, `|S|<0` only
  mis-signs and is caught at config-load.

**Tests after hardening:** Slice-K suites **177 green** (141 → +36); full suite **1222
passed, 0 failed** (same 17 unrelated warnings). Next: Slice U (`internal_energy_budget.py`).

---

## Phase A — Slice U DELIVERED (2026-06-30) — Phase A COMPLETE

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN); pure, stateless, mass-agnostic (no `γ`, no `m`, no RNG, no
integrator, no `E_int` *state* — the reservoir is Phase C). Oracle asserts against the
MASS §6 budget rules + the plan §2/§4 golden values. Four open questions (OQ-U1–U4)
were resolved by the user *before* coding, then a final MASS/CALIBRATION cross-check
ran clean.

**Build (4 files):**
- `physics/internal_energy_budget.py` — `e_int_onset_eV` (S2 `f_int·E_avail`),
  `dE_int_pickup_eV` (S1 `+f_ret·D_0(n+1)`) + `pickup_bath_release_eV`
  (`(1−f_ret)·D_0(n+1)`), `dE_int_shed_eV` (K1 `−D_0(n)`), `reconstruct_e_int_eV`
  (A9 `E_solv.struct − E_bind^pair − E_elec`, **raises** when `post_crossing=False`),
  `f_int_floor` (`Σ(n*)/E_avail`, `e_avail>0` guard). Consumes L (`d0_of_n` /
  `ladder_cumsum`) and K (`e_bind_pair_eV` / `e_electrostriction_eV`). No new constants.
- `config.py` — two fields `internal_energy_partition_fraction` (f_int) and
  `internal_energy_retained_fraction` (f_ret), both `Optional[float] = None`;
  `check_internal_energy_budget_config` guard (**bounded-when-set** `0 ≤ f ≤ 1`; None →
  no-op) wired into `validate()`.
- `tests/test_internal_energy_budget.py` (294) + `tests/test_internal_energy_budget_config.py` (21).

**OQ-U1–U4 resolved (user, 2026-06-30, before coding):**
1. **U1 — config surface.** Both `f_int` and `f_ret` land as `Optional[float] = None`
   declared-but-unread (the *module* takes them as kwargs; the Phase-C driver reads the
   config), MASS §11 names. Guard = **bounded-when-set** `[0,1]` (hard cap, CALIBRATION
   rows 14/13); the soft ~0.2 ceiling + the self-unbound floor stay advisory, NOT
   enforced.
2. **U2 — `reconstruct_e_int_eV`.** `post_crossing` is a **required** keyword (no
   default); `False`/omitted raises (A9 — static reconstruction errs pre-t×). Composes
   the **real** K split; tests stub it.
3. **U3 — `f_int_floor`.** Bare scalar `ladder_cumsum(N_STAR,…)/e_avail_eV`, takes
   `e_avail_eV` explicitly (config-agnostic), `e_avail>0` guard, no headroom extras.
4. **U4 — onset does not enforce the floor.** `e_int_onset_eV` is pure arithmetic; a
   sub-floor `f_int` is accepted (MASS S2 / CALIBRATION row 14: advisory, not a
   constraint). `f_int_floor` is diagnostic only.

**As-built refinements vs plan §3:**
- The plan's "fail-loud: L's negative/fractional-n guards reach U" maps onto the
  **`reconstruct`→K→`ladder_cumsum`** path (negative/fractional `N` raises). The
  d0-based `dE_int_pickup_eV`/`dE_int_shed_eV` inherit L's *smooth, deliberately
  unguarded* `d0_of_n` (only the cumulative/indexed path guards) — consistent with L's
  design, so no negative-n guard was added to the thin wrappers.
- f_int floor band asserts use **`round(floor, 2)`** band-membership (the Slice-L
  precedent): the κ-endpoints round to the MASS-stated edges (e.g. κ=5 mix
  Σ(21)=0.19306 → /0.80 = 0.2413 → 0.24). The 2.70 eV mix band is MASS's "~0.065" at
  two decimals → [0.06, 0.07].
- `e_int_onset_eV` is **config-agnostic** (no `f_int`/`e_avail` guard in the module),
  mirroring K's deliberate non-guard stance — bounds are caught at config-load.
- `E_int^eq = 0` asserted but annotated as a split-consistency **tautology** (MASS line
  881), not an independent anchor (carries the Slice-K annotation forward to U).

**Oracle cross-validation.** S2 onset exact (`0.3·0.80 = 0.24`); S1 split closes to
machine precision (`f_ret·D_0 + (1−f_ret)·D_0 = D_0`); K1 = `−D_0(n)`; shed pair+E_int
sub-sum neutral to 1e-15 (jump-consistency); reconstruction matches the split exactly
and recovers `E_int^eq = 0` at equilibrium; f_int floor lands in every scenario×picture
band (X₂ 0.31–0.35 / mix 0.21–0.24 @ 0.80 eV; X₂ 0.09–0.10 / mix ~0.065 @ 2.70 eV).
Independence test stubs `d0_of_n` to prove pickup/shed compose L only through the rung.

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. f_int/f_ret
**Bounded** (rows 14/13); config field names match MASS §11 (lines 1987/1970); the soft
~0.2 ceiling + scenario-keyed floor confirmed advisory (row 14); reconstruction post-t×
only (A9).

**Rule-2 table.** Slice U **adds two declared-but-unread carries** —
`internal_energy_partition_fraction` (f_int) and `internal_energy_retained_fraction`
(f_ret) — read by the Phase-C generative driver, not by Phase A. They join the
`cooling_relaxed` *concrete blend* (Slice L; pinned at Phase F) as the standing Phase-A
rule-2 entries. (L's and K's own fields were born-live; U's are the first Phase-A fields
that are genuinely deferred.)

**Tests:** Slice-U suites **315 green** (`test_internal_energy_budget.py` 294 +
`test_internal_energy_budget_config.py` 21); full suite **1537 passed, 0 failed** (same
17 unrelated `anchored_discrete` warnings).

### Slice U — review + test-hardening pass (2026-06-30)

An inline review (independent numerical re-verification of every §6 oracle on a fresh
interpreter — onset, S1 pickup/bath/split, K1 shed, the f_int floor across κ/picture/
budget, and `E_int^eq=0`, all matched to machine precision — plus a multi-angle
correctness/edge/vectorisation/independence sweep) **confirmed the physics with no bug**
and surfaced only coverage gaps; all closed test-first.

- **No physics bug, no code change.** The module was correct and complete; the
  hardening is regression/coverage locks (mirrors the Slice-K review stance).
- **De-speculated the `e_int_onset_eV` array branch (rule 2).** It was unread → now
  exercised by vectorised tests (array `f_int`; scalar-`f_int`/array-budget Phase-F
  sweep). Kept (not removed) — consistent with the module's vectorisation contract.
- **Coverage added (+16 tests):** vectorised reconstruction over `N` + scalar
  float-return discipline; a **stub-independence** lock proving `reconstruct` composes K
  *only* through `e_bind_pair_eV` + `e_electrostriction_eV`; the deliberate
  **module-level non-guard** lock (out-of-`[0,1]` `f_int` is *computed*, not rejected —
  the bound is a config-load concern, mirroring K); `f_int_floor` **negative**-budget
  guard (only zero was covered); return-`ndarray` dtype locks on the vectorised helpers;
  the `cooling_relaxed` picture-ordering carry on the S1 increment
  (mix < cooling_relaxed < x2); and config-guard both-set paths (both-valid pass,
  offender-named on a mixed set, `[0,1]` inclusive endpoints).
- **Deliberate non-findings (consistent with L/K, not gaps):** the d0-based
  `pickup`/`shed` wrappers stay **unguarded** for negative/fractional/zero `n` —
  they inherit L's *smooth* `d0_of_n` by design (only the cumulative/indexed path
  guards); a guard would diverge from L. `e_int_onset_eV` stays **config-agnostic**
  (no `f_int`/`e_avail` range check in the module) — bounds caught at config-load.

**Tests after hardening:** Slice-U suites **331 green** (315 → +16:
`test_internal_energy_budget.py` 307 + `test_internal_energy_budget_config.py` 24); full
suite **1553 passed, 0 failed** (same 17 unrelated `anchored_discrete` warnings).
**Phase A (L + K + U) is complete.** Next: Phase B (Slice ρ `helium_density.py`, then
P/Q) — `TIER2_PHASE_B_IMPLEMENTATION_PLAN.md`.

---

## Phase B — Slice ρ refinement decisions (2026-06-30)

A pre-build refinement pass over the Slice ρ spec in
`TIER2_PHASE_B_IMPLEMENTATION_PLAN.md`, cross-referenced against
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §4/§5 + dim-table and
`CALIBRATION_MAP.md` rows 5/8, plus a code-surface audit (`physics/drag.py::spatial_gate`,
`config.py`, `constants.py`). No code — the `[PROCEED TO IMPLEMENTATION]` boundary holds;
docs-only amendment (Slice ρ §3.1 record + §5 config bullet/validators). The spec carried
four drift items and three genuinely-open design calls; all resolved (decision owner: user).

- **1. Single-source erf → extract `_erf_complement`.** Both `drag.spatial_gate` and the
  new `rho_he_ratio` route through one extracted helper (proposed `physics/_gates.py`).
  Chosen over delegate-to-drag (would make ρ import the drag module, muddying its
  independent role + tautologizing the parity test) and over an independent mirror
  (duplicate formula). The `spatial_gate` rewire is a **Tier-0 no-behavior-change refactor**
  guarded by a parity regression test — a build-time edit to a locked module with identical
  output, *not* a drag-law physics change.
- **2. No default steepness on `rho_he_ratio`.** Mirrors `spatial_gate` (which has none);
  the Phase-C driver passes `cfg.drag_gate_steepness` (14.2 Å) so density and drag share the
  identical surface. No new `POTENTIAL_STEEPNESS_ANGSTROM` constant (the planned default
  referenced a constant that does not exist; 14.2 already lives as `potential_steepness` /
  `potential_steepness_molecule` / `drag_gate_steepness`). Module stays config-agnostic.
- **3. Explicit `HeliumDensityProfile` Literal selector.** New
  `Literal["erf_complement","tabulated"]`, default `erf_complement`, enum reject-arm guard
  (mirrors `mass_scenario` / `check_ladder_config`), repurposing the `Optional[object]`
  placeholder at `config.py:261` (stale `# future G4 density profile` comment rewritten — ρ
  stays **G2**, no G2→G4 promotion). Tabulated **machinery built + round-trip tested in
  Phase B** (Slice-L `tabulated_ladder` precedent); the **sourced TDDFT data array deferred
  rule-2** (CALIBRATION row 8); tabulated data passed as a module arg, not a config field.

**Drift fixed:** stale `config.py:220`→`:261`; misleading `# future G4` comment;
non-existent `POTENTIAL_STEEPNESS_ANGSTROM`; the "tabulated not built vs round-trip tested"
tension (machinery built, data deferred); the self-contradictory "shared helper" + "two
forms diverging" parity framing (now a single-source regression lock).

**Cross-reference verdict:** no contradictions. Gate is **Derived/G2** (row 5, erf-tied
until ρ_He profile exists), source = 14.2 Å steepness; ρ_He profile is **Sourced** (row 8),
the deferred fallback. ρ is purely geometric (depth + steepness only) — no `n`, no Langmuir
cap (that is P), no `γ`/drag-law read, no mass. **Next:** Slice P/Q refinement, then the
Phase-B build under the `[PROCEED TO IMPLEMENTATION]` trigger (ρ first).

---

## Phase B — Slice P/Q refinement decisions (2026-07-01)

A pre-build refinement pass over the Slice P (pickup + `mass_jump.capture`) and Slice Q
(evaporation) specs in `TIER2_PHASE_B_IMPLEMENTATION_PLAN.md`, cross-referenced against
MASS §4/§5 + dim-table + §11 config + A11/A12/R9/R10, CALIBRATION rows 7/8/10/13, and a
code-surface audit (`physics/mass_jump.py`, `internal_energy_budget.py`,
`dissociation_ladder.py`, `config.py`, `constants.py`). No code — the
`[PROCEED TO IMPLEMENTATION]` boundary holds; docs-only amendment (Slice P/Q Interface/Knobs/
Oracle/Test-spec, §2.1/§2.2, §3.2 record, §4 table, §5 config). Seven design calls resolved
(decision owner: user).

**Slice P**
- **1. New `CaptureResult` dataclass** (not a reused `ShedResult`) — the type name matches the
  +He *gain* (`m⁺ = m+m_He`, `dE_mass_transfer > 0` into `E_mass_transfer`).
  `_reduced_mass_defect_coeff` is refactored to the **unsigned magnitude** `0.5·(m·m_He)/m_plus`;
  `cold_shed` applies the `−`, `capture` the `+` (rule 1, one formula, signs at the call site);
  `cold_shed` output stays byte-identical (regression lock). The S1 heat (`+f_ret·D_0(n+1)` into
  `E_int`) and the capture KE defect (into `E_mass_transfer`) are **two distinct injections** —
  no double-count; the 5-term closure combining them is Phase C.
- **2. `pickup_occupancy_exponent` p = fixed `1.0`, NOT `p=κ`.** The A12 p↔κ tie is **physically
  inverse** (rigid shell = large κ = sharper cutoff = *smaller* p), so a literal `p=κ` inverts
  it. Hold p=1 for calibration; free p at Phase F only if the size-dist first-shell edge can't
  be met with p=1 + κ (if ever tied, match cutoff *slopes*). ⚠ **MASS §11 "p default tied to κ"
  + CALIBRATION row 7 flagged for a clarifying annotation** — decision recorded, locked docs not
  silently rewritten (pending user OK).
- **3. Atomic rename, no alias; `mass_relaxation_tau_ps` retired.** `mass_rate_form` /
  `mass_rate_coefficient` / `mass_relaxation_tau_ps` have **zero readers** in the package
  (grep-verified) → `→ pickup_*` atomically; the τ field is dead (Tier-2 cooling is the
  Slice-K `internal_energy_cooling_tau_ps`) → removed at the P build.
- **4. `density_only` + `at_rest` built; `sweeping`/`dwell_time`/`thermal` are rule-2 arms**
  (raise `NotImplementedError`) — MASS §5 locks density-only, CALIBRATION R7 defers the rest.
- **U-signature confirmed:** `dE_int_pickup_eV(n, *, f_ret, picture, kappa)` takes `n` =
  pre-pickup and returns `+f_ret·D_0(n+1)`; `picture`+`kappa` are threaded through `pickup_step`.

**Slice Q**
- **5. `effective_dof`: `n=2 → 4`** (linear 3-atom `3N−5`; CALIBRATION row-10 "n=2 linear +1"
  applied), `3n−3` for `n≥3`, `n=1` direct `k=ν`. The uniform `3n−3` is no longer used at n=2
  (n=2 bracket exponent `s−1 = 3`, was 2).
- **6. Gate onset — 3 facets:** (a) **no calibration knob** (Derived, gate = `Σ(n)`,
  parameter-free, MASS R9); (b) **diagnostic** signed margin `gate_margin_eV = E_int − Σ(n)`
  exposed as a pure Q helper — the record of `Σ(n,t)` / `G(t)` alongside `t×` / `Π` is written by
  the Phase-C driver / Phase-E diagnostics, not in Phase B; (c) **`gate_onset_override_eV:
  None`** with a **loud provenance guard** (`allow_gate_onset_override`, mirroring
  `allow_unvalidated_binding_pairing`) so a forced threshold can never silently enter a
  production / Tier-2-lock run (R10 diagnostic-lever, not a knob).
- **7. `NU_EVAP_PER_PS = 2.42` is added by Q** (constants anchor + `evap_rate_prefactor_per_ps`
  defaulting to it) — the plan's "ν lands with Phase A" was **inaccurate** (Phase A shipped
  `N_STAR` / `D0_*` / `D_FLOOR` / `S_ABS_EV` only). `evap_rrk_dof` lands `Optional[float] = None`
  (per-`n` default; `s≥1` guard fires only when set).
- **U-signature confirmed:** `dE_int_shed_eV(n, *, picture, kappa)` takes `n` = pre-shed and
  returns `−D_0(n)`; Q reuses the **real** `cold_shed` (no second copy of the reduced-mass form).

**Drift fixed:** stale config line refs (`:217/:218` → `:258/:259`); the "ν lands with Phase A"
gap; the reused-`ShedResult`-for-a-gain semantic mismatch; the `p=κ` ambiguity.

**Cross-reference verdict:** consistent with MASS/CALIBRATION **except** the flagged `p` wording
(§3.2 item 2) — a clarifying annotation to MASS §11 + CALIBRATION row 7 is the one carried
follow-up (awaiting user OK). ν pinned (row —/Sourced), λ_0 Sourced+Bounded (row 7), f_ret
Bounded (row 13), `s` Derived (row 10, n=2-linear now applied). **Next:** the Phase-B build
under the `[PROCEED TO IMPLEMENTATION]` trigger — ρ → {P, Q}.

### Phase B/C boundary — `mass_jump_velocity_reset` reconcile (2026-07-01)

A follow-up check on the MASS §11 config block found `mass_jump_velocity_reset ∈
{momentum_conserving, label_only}` (line 1940) absent from the Phase-B config map. **Verdict:
correctly so — it is Phase C, not a Phase-B omission.** It is A13 *driver/integrator policy*
and already owned by `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §4 ("Integrator (G)"), alongside
`one_mass_event_per_step` and `jump_o_step_ordering`. Phase-B `capture`/`cold_shed` are
momentum-conserving **by construction**, so `label_only` (the non-closing §6-invariant
diagnostic) is a Phase-C driver branch Phase B never implements.

- **Real drift found + fixed:** `he_capture_velocity` was **double-listed** (Phase-B Slice P
  *and* Phase-C §4 "integrator flags"). Resolved via the f_int/f_ret pattern — **declared at
  Slice P** (it owns the `capture` primitive; declared-but-unread in Phase B) and **activated
  by the Phase-C G driver** (reads the field, passes `u_he`). Both plan docs updated so the
  field has one declaring owner (P) and one activator (C), no double-add.
- **Docs touched:** Phase-B plan Slice P knobs (he_capture_velocity note) + §5 (new "A13
  integrator-policy fields are Phase C" boundary note); Phase-C plan §4 (split "added here"
  vs "activated here, declared at P"). No code; boundary holds.

### Phase B — Slice ρ pre-build open questions resolved (2026-07-01)

A final discussion pass over the Slice ρ spec before the `[PROCEED TO IMPLEMENTATION]`
trigger, grounded in a code-surface re-audit (`physics/drag.py::spatial_gate`,
`simulation/ion.py::_drag_gate_steepness`, `config.py:261`) and cross-checked against
CALIBRATION_MAP rows 5/8 + MASS §4/§5. No code — boundary holds; docs-only amendment to the
Phase-B plan (interface bullets, §3.1 item 2, §5). Three calls resolved (decision owner: user).

- **1. Steepness source corrected → `_drag_gate_steepness(cfg)`, not `drag_gate_steepness`.**
  The prior wording had the Phase-C driver thread the raw `cfg.drag_gate_steepness` field to
  `rho_he_ratio` "so density and drag share one surface." **But the live drag path does not
  read that field in production:** `simulation/ion.py::_drag_gate_steepness` returns
  `cfg.potential_steepness` under the default/production `density_proportional` (and
  `erf_tied`) gate, and `cfg.drag_gate_steepness` **only** under G3 `erf_independent`. Both
  default to 14.2 so it is numerically identical today, but threading the raw field would let
  ρ and drag silently diverge if the two were ever set unequal — the exact failure the "one
  surface" decision exists to prevent. **Fix:** the driver threads the **resolver**
  `_drag_gate_steepness(cfg)`, which tracks whatever surface drag actually uses across all
  gate modes. Confirmed against **CALIBRATION row 5** (gate source = "confining-potential
  steepness (14.2 Å)" = `potential_steepness`), which the resolver returns in production; the
  raw `drag_gate_steepness` field is the G3 override, not the sourced anchor. ρ module stays
  config-agnostic (takes `steepness` as an arg — unchanged).
- **2. `_erf_complement` home locked → `physics/_gates.py`.** The plan's "final placement
  settles at build" is resolved now: a gate-specific neutral module so ρ never imports
  `drag.py` (preserving ρ's drag-independent role — the whole reason extraction beat
  delegate-to-drag). `drag.spatial_gate` is rewired to call it (Tier-0 no-behavior-change
  refactor, parity-regression-locked).
- **3. `tabulated_density_profile` interface locked.** Takes a hand-built
  `(depth_grid, ratio_grid)`, **linear interpolation on `depth`** (not `r`), round-trip exact
  at grid nodes, and **clamps out-of-range `depth` to the tail values** (`1.0` inside, `0.0`
  outside) rather than raising — matching the erf-complement asymptotes so both
  `HeliumDensityProfile` arms share one `[0,1]` boundary contract. Machinery built + round-trip
  tested in Phase B (Slice-L `tabulated_ladder` precedent); the sourced TDDFT `ρ_He(r)` array
  stays the deferred rule-2 carry (CALIBRATION row 8).

**Cross-reference verdict:** no contradictions. Row 5 (Derived/G2, erf-tied, source
`potential_steepness` 14.2 Å) actively supports the resolver fix; row 8 (Sourced ρ_He, deferred
fallback) matches the tabulated-data deferral; MASS §4/§5 say nothing about the ρ gate steepness
(its `steepness` hits are all the ladder κ). **Slice ρ is build-ready** — next step is the
Phase-B build under the `[PROCEED TO IMPLEMENTATION]` trigger (ρ first, then {P, Q}).

---

## Phase B — Slice ρ DELIVERED (2026-07-01)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN); pure, stateless, config-agnostic, mass-agnostic (depth + steepness
only — no `n`, no Langmuir cap, no `γ`, no mass, no RNG, no integrator, no `E_int`
state). Oracle asserts against the plan §2.2 (ρ) / §3.1 / §4 golden values.

**Build (5 files):**
- `physics/_gates.py` — **new leaf module.** `_erf_complement(depth, steepness)`, the
  single source of `0.5·(1−erf(depth/steepness))` with the `steepness>0` fail-loud
  guard. A neutral home so `helium_density` never imports the drag law and `drag`
  never imports a Phase-B module — both import this leaf (rule 1). Depends only on
  numpy + `scipy.special.erf`.
- `physics/drag.py` — **rewired** `spatial_gate` to `return _erf_complement(depth,
  steepness)` (dropped the local `erf` import + inlined arithmetic + guard). Tier-0
  **no-behaviour-change** refactor: output byte-identical, locked by parity regression
  tests from *both* sides (`test_drag.py::TestSpatialGateSingleSource` +
  `test_helium_density.py::TestSpatialGateParity`).
- `physics/helium_density.py` — **new module.** `rho_he_ratio(depth, *, steepness)`
  (the erf-complement gate, no default steepness, scalar→float/array→ndarray);
  `TabulatedDensityProfile` + `tabulated_density_profile(depth_grid, ratio_grid)` (the
  declared sourced-profile fallback: `np.interp` linear-on-`depth`, tail-clamped to
  endpoint ratios, fail-loud on length mismatch / non-increasing depth / ratio∉[0,1] /
  empty). Data array deferred (rule-2; CALIBRATION row 8).
- `config.py` — `HeliumDensityProfile = Literal["erf_complement","tabulated"]`; field
  `helium_density_profile` **repurposed** from the `Optional[object]=None` G4
  placeholder → `"erf_complement"` default (stale `# future G4 density profile` comment
  rewritten — ρ stays **G2**); `check_helium_density_config` enum reject-arm guard
  (mirrors `check_ladder_config`) wired into `validate()`.
- `tests/test_helium_density.py` (19) + `tests/test_helium_density_config.py` (6) +
  2 parity asserts appended to `tests/test_drag.py`.

**As-built notes:**
- **Steepness resolver deferred to the caller.** `rho_he_ratio` has **no default
  steepness** (mirrors `spatial_gate`); the plan's `_drag_gate_steepness(cfg)` threading
  is a **Phase-C driver** concern — Slice ρ is config-agnostic and takes `steepness` as
  a kwarg, so nothing here reads `cfg`. The one config touchpoint is the profile-selector
  guard.
- **`ion.py::_drag_gate_steepness` docstring de-drifted** (docstring only, no behaviour):
  the stale "at Tier 0 no `helium_density_profile` exists" line now states the field is
  the *pickup* G2 occupancy gate, separate from the drag gate (no G4 promotion), sharing
  the resolver's steepness.
- **Scalar/array discipline** follows the `shell_schedule`/`dissociation_ladder` idiom
  (`float(out) if np.ndim(depth)==0 else out`); `_erf_complement` returns the raw np
  result (0-d for scalar) so `spatial_gate` stays byte-identical.

**Oracle cross-validation.** ratio = 0.5 at depth 0; **exact** 1.0 inside / 0.0 outside
at the far tails (erf saturates in float — no epsilon); monotone non-increasing; ρ ==
`spatial_gate` to machine precision on a 201-pt grid (single-source lock); both callers
fail-loud on `steepness ≤ 0`; tabulated round-trips grid nodes exactly, interpolates
linearly (0.5 at the midpoint), clamps tails, and rejects malformed tables.

**Rule-2 table.** `helium_density_profile` is **born live** (read by
`check_helium_density_config` on arrival) — nothing added to the exception table. The
sourced **TDDFT `ρ_He(r)` data array** is the one Slice-ρ rule-2 carry (machinery built +
round-trip tested; concrete array pinned when the profile exists; CALIBRATION row 8).

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. Gate stays
Derived/**G2** (row 5, source = confining-potential steepness); ρ_He profile **Sourced**
(row 8, deferred fallback). ρ is purely geometric (depth + steepness) — no `n`, no cap,
no `γ`, no mass.

**Tests:** Slice-ρ suites **27 green** (`test_helium_density.py` 19 +
`test_helium_density_config.py` 6 + 2 `test_drag.py` parity asserts); full suite **1580
passed, 0 failed** (was 1553; +27; same 17 unrelated `anchored_discrete` warnings).

### Slice ρ — review + test-hardening pass (2026-07-01)

An inline adversarial review (10-angle correctness / edge / vectorisation / independence /
import / doc-drift sweep over the small self-authored diff) confirmed the physics and
surfaced **two genuine robustness/placement gaps** plus coverage gaps; all closed test-first
(the guard gap test→RED→GREEN, the rest regression/characterisation locks).

- **Genuine gap 1 — `TabulatedDensityProfile` direct construction bypassed validation.**
  Only the `tabulated_density_profile` *builder* validated; constructing the frozen dataclass
  directly with an unsorted `depth_grid` (or a ratio ∉[0,1], length mismatch, empty table)
  passed silently, and `np.interp` on unsorted `xp` returns **silently-wrong** values (worse
  than a crash; principle 4). **Fix:** validation **centralised in `__post_init__`** (the
  single source for both the builder and any direct construction); the builder is now a thin
  coercion wrapper. Test-first (4 direct-construction guard tests → RED → GREEN).
- **Genuine gap 2 — config field placement drift.** `helium_density_profile` was born-live
  (read by `check_helium_density_config`) but sat under the `# -- Deferred (no Tier-0
  reader) --` block. **Fix:** relocated to a dedicated `# -- Tier-2 Phase-B helium density
  gate (Slice rho) --` section beside the Phase-A L/K/U live blocks. No-behaviour-change field
  reorder among defaulted fields; `test_run_directory` round-trips confirm serialisation
  unaffected.
- **Coverage added (+8 regression/characterisation locks):** `ratio()` array-in→ndarray;
  **deliberate non-guards locked** — a **non-monotone ratio_grid is accepted** (a sourced
  TDDFT profile may show a near-surface density lip; only the *depth* axis must be sorted for
  `np.interp`) and a **single-node table is a valid constant profile** (matches
  `tabulated_ladder`'s ≥1 contract); **NaN steepness rejected** by the single guard (both
  callers, since `not (nan > 0)` is True); scalar-return discipline across
  python-int / 0-d-array / list / ndarray inputs to `rho_he_ratio`.
- **`ion.py::_drag_gate_steepness` docstring de-drift** (from the delivery pass) retained —
  docstring only, no behaviour.
- **Deliberate non-findings:** NaN/inf `depth` is **not** guarded (an upstream bug if it
  occurs; the drag `spatial_gate` never guarded it either — consistency); `ratio()`
  rebuilding arrays from the stored tuples per call is an accepted micro-cost for a lookup
  helper; no `__init__.py` export added (Slice P imports the module directly).

**Import audit.** Removing `from scipy.special import erf` from `drag.py` (now routed through
`_gates`) strands no importer — every `physics.drag` consumer imports only the public symbols
(`spatial_gate` / `drag_force` / `drag_gamma` / `DragCoefficients` / form tags /
`REALIZED_FORMS` / `_REQUIRED_COEFF_KEYS`), grep-verified.

**Tests after hardening:** Slice-ρ suites **39 green** (`test_helium_density.py` 31 +
`test_helium_density_config.py` 6 + 2 `test_drag.py` parity asserts); full suite **1592
passed, 0 failed** (1580 → +12; same 17 unrelated `anchored_discrete` warnings). Next:
Slice P (`physics/pickup.py` + `mass_jump.capture`) / Slice Q (`physics/evaporation.py`).

---

## Phase B — Slice P pre-build interface decisions (2026-07-01)

A pre-build discussion pass over the Slice P spec before the `[PROCEED TO IMPLEMENTATION]`
trigger, grounded in a code-surface audit (`physics/mass_jump.py`,
`physics/internal_energy_budget.py` U signatures, `config.py:275–277` `mass_rate_*` fields).
No code — the boundary holds; docs-only amendment (Slice P interface bullets §3, §5 rename
note, `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §4 for the deferred `m↔n` guard). Five interface
calls resolved (decision owner: user) plus two code-audit corrections to the plan.

- **Q1 — `PickupResult` dataclass out, loose kwargs in.** `pickup_step` returns a
  `PickupResult(n_plus, m_plus_amu, v_plus, dE_int_eV, dE_mass_transfer, fired)` dataclass
  (symmetric to `ShedResult`/the new `CaptureResult`, not a bare 6-tuple). Inputs are loose
  scalar kwargs — the persistent `(n, E_int, v, m)` carrier is G's `IonStepState` (Phase C),
  not a Phase-B state container.
- **Q2 — vectorized-over-ensemble form built in Phase B.** A `pickup_step_components(...)`
  variant over the `(2N,)`/N-ion arrays (mirroring `cold_shed_velocity_components`), one
  vectorized Bernoulli draw per step. Chosen over scalar-only + driver-loop because the
  Poisson-moment + cross-ion-independence oracles exercise the ensemble directly; the scalar
  `pickup_step` is the single-ion oracle the components form is checked against.
- **Q3 — `m↔n` consistency deferred to G, documented there.** The Phase-B primitives advance
  mass (`capture`: `m⁺=m+m_He`) and integer `n` (`pickup_step`: `n→n+1`) independently;
  `m == m_I⁺ + n·m_He` is documented-but-not-enforced in the primitives (no Phase-B
  persistence/closure). Recorded in `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §4 (Slice G) as a
  per-step structural invariant G asserts alongside the 5-term energy closure — G owns the
  guard (first place a persistent `(n_shell, mass)` state co-evolves). Primitives stay
  mass-based so `capture`/`cold_shed` remain the single reduced-mass reset source (rule 1).
- **Q4 — unbuilt-arm `NotImplementedError` at point-of-use, not config-load.**
  `sweeping`/`dwell_time`/`thermal` are *valid* enum literals → `validate()` accepts them (a
  carrying config round-trips, the rule-2 declared-but-unread contract); the error fires lazily
  inside `lambda_attach`/`pickup_step`/`capture`. The enum reject-arm rejects only invalid
  typos, never a valid-but-unbuilt member. (Best-practice call, user-delegated.)
- **Q5 — RNG draw convention (documented here, locked at X).** One vectorized uniform draw per
  step, `rng.random(size=n_ions)`, fire on `draw < P_attach`; shed-then-pickup channel order
  documented for one-event-per-step (A13) but frozen only at the Phase-C schema/driver
  (forbidden-list draw-order lock lands at Slice X).

**Two code-audit corrections to the plan (were inaccurate, now fixed in the Phase-B plan):**
- **"zero readers" was package-only.** The atomic `mass_rate_*→pickup_*` rename (no alias)
  cascades to **two non-package readers**: `tests/test_drag_config.py:92`
  (`assert cfg.mass_rate_form == "density_only"` → migrate to `pickup_rate_form`) and
  `docs/config_and_preset.md:152`. Both must be touched at the P build.
- **The coeff refactor has two call sites.** `_reduced_mass_defect_coeff` is consumed by
  `cold_shed` (`mass_jump.py:150`) **and** `cold_shed_velocity_components` (`:234`); moving the
  sign to the call site means negating at both, plus rewriting the coeff docstring. Both stay
  byte-identical (regression-locked from both sides).

**Cross-reference verdict:** these are software-interface decisions (result-shape,
vectorisation, error placement, RNG draw), not physics-parameter changes — no MASS/CALIBRATION
value is touched. The `m↔n` (Q3) and draw-order (Q5) items are consistent with A13
(one-event-per-step, jump-then-O) and the 5-term invariant ownership already in the Phase-C
plan. **Slice P is build-ready** pending the `[PROCEED TO IMPLEMENTATION]` trigger (P, then Q).

### Phase B — Slice P pre-build open question resolved (2026-07-01)

A final discussion pass before the trigger, grounded in a code-surface re-audit
(`physics/mass_jump.py::_check_masses`, `_reduced_mass_defect_coeff` call sites,
`config.py:57/275–277`, `tests/test_drag_config.py:36/92/658`). One genuine physics-guard
fork surfaced; resolved (decision owner: user). No code at decision time — boundary held.

- **`capture` mass guard → dedicated positivity-only guard (not shared `_check_masses`).**
  The plan test-spec said "reuse `_check_masses`", but that guard imposes the *shed*
  precondition `m_minus − m_He > 0` (post-shed mass positive), which a **capture does not
  have** — a capture *adds* mass (`m⁺ = m + m_He`, always positive). Reuse never wrongly
  fires for a physical I⁺ (127 ≫ 4 amu He) but would fail loud for the *wrong reason* on a
  light ion. **Decision:** a new `_check_masses_gain` requiring only finite + positive
  `m_minus`/`m_He` (CLAUDE.md principle 4 — fail loud for the right reason). `cold_shed`
  keeps `_check_masses` unchanged.
- **Rename cascade correction.** The plan's §5 "two non-package readers" undercounted: the
  atomic `mass_rate_*→pickup_*` rename also renames the **type alias itself**
  (`config.py:57` `MassRateForm → PickupRateForm`), cascading to `test_drag_config.py:36`
  (import) and `:658` (members assert) in addition to `:92` (field assert) and
  `docs/config_and_preset.md`. All migrated at build.

---

## Phase B — Slice P DELIVERED (2026-07-01)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. Pure-ish stochastic
primitive: injected RNG, per-event mass/velocity resets, **no** integrator, no `E_int`
state persistence, no schema I/O, no 5-term closure (all Phase C). Mass-agnostic in the
drag sense — no Phase-P function evaluates `γ`; `m` enters only as the captured quantity in
the reset. Oracle asserts against the plan §2.2 (P) / §4 golden values.

**Build (files):**
- `physics/mass_jump.py` — **extended.** New `CaptureResult` dataclass (the +He gain
  counterpart of `ShedResult`); `capture(v_minus, m_minus_amu, *, m_he_amu, u_he=0.0)`
  (`v⁺ = (m·v⁻ + m_He·u_He)/(m+m_He)`, `m⁺ = m+m_He`, defect `+½(m·m_He)/(m+m_He)·|v⁻−u_He|²
  > 0`) + `capture_velocity_components` (per-ion array form, **per-ion mass** since fires
  diverge — unlike the scalar-mass `cold_shed_velocity_components`); `_check_masses_gain`
  (positivity-only guard, per the resolved open question). **Coeff refactor:**
  `_reduced_mass_defect_coeff` now returns the **unsigned magnitude** `0.5·(m·m_He)/m_plus`;
  the `−` is applied at both shed call sites (`cold_shed`, `cold_shed_velocity_components`),
  the `+` at both capture call sites — one reduced-mass formula, signs at the call site
  (rule 1). `cold_shed` / `cold_shed_velocity_components` output **byte-identical**
  (locked by the pre-existing Tier-1a reset/energy tests, all green).
- `physics/pickup.py` — **new module.** `attach_probability` (`1−e^{−λdt}`), `lambda_attach`
  (`λ₀·ρ·(1−n/n*)_+^p`, Langmuir/`none` cap), `pickup_step` → `PickupResult(n_plus,
  m_plus_amu, v_plus, dE_int_eV, dE_mass_transfer, fired)` (composes `mass_jump.capture` +
  `internal_energy_budget.dE_int_pickup_eV`, `n` = pre-pickup, picture/kappa threaded), and
  `pickup_step_components` (one vectorized `rng.random(size=M)` draw, per-ion independence,
  `np.where` fire-mask). Rule-2 arms (`sweeping`/`dwell_time` rate form, `thermal` capture
  velocity) raise `NotImplementedError` **lazily at point-of-use** (Q4), never at config-load.
- `config.py` — atomic rename (no alias) `MassRateForm→PickupRateForm`,
  `mass_rate_form→pickup_rate_form`, `mass_rate_coefficient→pickup_rate_coefficient`;
  **removed** dead `mass_relaxation_tau_ps` (superseded by `internal_energy_cooling_tau_ps`).
  New Literals `PickupOccupancyCap`/`HeCaptureVelocity`; new fields `pickup_occupancy_cap`
  (langmuir), `pickup_occupancy_exponent` (p=1.0), `he_capture_velocity` (at_rest). New
  `check_pickup_config` enum-reject guard (accepts valid-but-unbuilt rule-2 members, rejects
  only typos) wired into `validate()`.
- `tests/` — `test_pickup.py` (33) + `test_pickup_config.py` (9) new; capture asserts added
  to `test_mass_jump.py` (+10); rename cascade migrated in `test_drag_config.py` (import +
  field assert + members assert); `docs/config_and_preset.md` type table + deferred list
  updated (also de-drifted the stale `helium_density_profile` "None; future G4" doc line to
  its live G2 status).

**Oracle cross-validation.** Capture reset exact to machine precision (`v⁺`, `m⁺`, defect);
defect **> 0** and sign-opposite to `cold_shed`; reduced-mass form cross-checked by hand;
`u_He≠0` conserves centre-of-mass momentum. `attach_probability` matches `1−e^{−λdt}` and
`→ λdt` small; Langmuir cap `→0` at `n=n*=21` and clamps `≥ n*`; density-only/`none`
equivalence at `n=0`; `p=1` linear, `p>1` sharpens. Poisson moment: fire count over 40 000
seeded steps within a 5σ band of `λt`; ensemble fraction fired within 5σ of `P` over 200 000
ions; components form matches the scalar oracle ion-by-ion (all-fire), leaves non-fired ions
untouched (mixed pattern), and is identity under no-fire. S1 heat `= +f_ret·D_0(n+1)`
composed from the real U, picture/kappa threaded.

**Rule-2 table.** `pickup_rate_form` / `pickup_occupancy_cap` / `he_capture_velocity` are
**born live** (read by `check_pickup_config`). `pickup_rate_coefficient` (λ₀) and
`pickup_occupancy_exponent` (p) are the two Slice-P **declared-but-unread** carries (read by
the Phase-C driver — the f_int/f_ret precedent). `he_capture_velocity` is declared here (owns
`capture`) and *activated* at Phase-C G (passes `u_he`). The `sweeping`/`dwell_time`/`thermal`
enum arms and the sourced TDDFT `ρ_He` array (Slice ρ) remain the standing Phase-B rule-2
carries; the `cooling_relaxed` concrete blend (Slice L) is still pinned at Phase F.

**Carried follow-up — RESOLVED (user OK 2026-07-01).** The `p=κ` wording clarification (§3.2
item 2) is applied: every bare "default tied to κ" in MASS was annotated with the **inverse-tie
NB** (rigid shell = large κ = *smaller* p; **not** literal `p=κ`; Phase-B holds `p=1` fixed) —
MASS lines ~158 (A12 Form-B rationale), ~1737 (A12 coupling discussion), ~1806 (§11 knob-index
row); the §11 config block (~1934) and CALIBRATION row 7p / line 83 already carried it. Locked
values untouched (annotation only, not a rewrite); `p=1.0` remains the Phase-B default and is
freed only at Phase F if the size-dist first-shell edge demands it.

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. Capture defect is
the exact reduced-mass counterpart of the cold-shed defect (rule 1, one formula); λ₀
Sourced+Bounded (row 7), f_ret Bounded (row 13), p held fixed (A12). No `γ`/drag-law read, no
integrator, no `E_int` state — the out-of-scope guard holds.

**Tests:** Slice-P suites **52 green** (`test_pickup.py` 33 + `test_pickup_config.py` 9 +
`test_mass_jump.py` capture asserts 10); full suite **1644 passed, 0 failed** (1592 → +52;
same 17 unrelated `anchored_discrete` warnings).

### Slice P — review + test-hardening pass (2026-07-01)

An inline adversarial review (correctness / edge / vectorisation / independence / guard /
import angles over the small self-authored diff) **confirmed the physics with no bug** and
surfaced one genuine guard gap plus coverage gaps; all closed test-first (the guard gap
RED→GREEN, the rest characterisation/regression locks).

- **Genuine gap — `lambda_attach` silent inf/nan on pathological knobs.** The Langmuir
  factor `np.maximum(0, 1−n/n*)**p` gives `0**neg = inf` for `p<0` (at `n≥n*`) and `0/0 =
  nan` for `n_star<=0` (at `n=0`) — a silently-wrong *fire-every-step* rate. No config guard
  protects `p` (declared-but-unread until Phase C), so the module is the only defence. **Fix:**
  fail-loud on `p<0` / `n_star<=0` **inside the langmuir branch only** (both unused under
  `cap="none"`, which stays truly p/n\*-independent). Same div-by-zero/inf-producer philosophy
  as the Slice-K/L `τ>0` / `κ>0` guards (bounds that *crash or silently invert* are guarded;
  advisory bounds are caught at config-load).
- **Doc nuance (not a bug):** `attach_probability` can return exactly `1.0` (float underflow
  at large `λ·dt`), so the "[0,1)" docstring was corrected to "[0,1], saturating to 1 in the
  fire-every-step limit; non-decreasing in λ".
- **Coverage added (+15 characterisation/regression locks):** capture **energy-conservation
  identity** (`dE_mass_transfer == KE_before − KE_after` for the mechanical/COM channel, both
  `u_He=0` and a moving non-axis-aligned He); zero-relative-velocity → zero defect; components
  **scalar-mass broadcast** + per-ion momentum conservation with a moving He; `attach_probability`
  monotone/bounded + `dt=0` no-fire; `lambda_attach` monotone-non-increasing in `n` + zero-gate
  (`ρ=0`) → zero rate; and the physically-important **step-level gate properties** — a closed
  gate (`ρ=0`, the ion-exit termination) and a full shell under Langmuir **never fire even with
  a fire-forcing RNG**, while `cap="none"` lets a full shell still capture.
- **Deliberate non-guards (consistent with the sibling reset modules):** per-ion array-length
  mismatch in `pickup_step_components` is unguarded (parity with `cold_shed_velocity_components`,
  which trusts the caller); negative `dt_ps` / negative `n` are treated as upstream bugs (the
  drag/shed modules never guarded them either). `pickup_step` is documented single-ion (an array
  `rho_ratio`/`n` would trip `bool(array)` — a misuse, not a supported path).

**Import audit.** `pickup.py` imports only public symbols (`capture` /
`capture_velocity_components` from `mass_jump`, `dE_int_pickup_eV` from
`internal_energy_budget`, `MASS_HE_AMU`/`N_STAR` from `constants`); no cycle (config imports
physics lazily), no `__init__` export added (the Phase-C driver will import the module
directly — the Slice-ρ precedent).

**Tests after hardening:** Slice-P suites **67 green** (52 → +15); full suite **1659 passed,
0 failed** (1644 → +15; same 17 unrelated `anchored_discrete` warnings). Next: Slice Q
(`physics/evaporation.py`).

---

## Phase B — Slice Q pre-build interface decisions (2026-07-01)

A pre-build discussion pass over the Slice Q (evaporation) spec before the
`[PROCEED TO IMPLEMENTATION]` trigger, grounded in a code-surface audit (`physics/pickup.py`
as the delivered P pattern, `internal_energy_budget.dE_int_shed_eV`,
`dissociation_ladder.d0_of_n`/`ladder_cumsum`, `mass_jump.cold_shed`/`ShedResult`, and the
`config.py` `check_pickup_config` / `allow_unvalidated_binding_pairing` guards) and
cross-checked against MASS §R9 + the plan §2.2/§4. The Slice-Q interface text predated the
P delivery decisions and carried drift; **eight calls resolved** (decision owner: user). No
code — the boundary holds; docs-only amendment to the Phase-B plan (Slice Q interface bullets,
Knobs, §5 validators, new §3.3 record).

- **1. Gate boolean corrected — the one physics-critical fix.** The spec was
  self-contradictory: line 432 defined `is_self_bound(...) -> E_int > Σ(n)` ("suppress while
  True"), but §2.2/§4 + line 466 say `gate_margin_eV = E_int − Σ(n) < 0` "exactly when
  self-bound". **MASS §R9 settles it:** `G ≡ E_int − Σ(n)` is the **self-*un*bound margin**;
  evaporation is *suppressed* while net self-unbound (`E_int > Σ(n)`, `G > 0`) and *enabled*
  once self-bound (`E_int < Σ(n)`, `G < 0`). The **physics direction was always right** (it
  matches the encoded form); only the `is_self_bound` helper had its boolean inverted vs its
  name. **Resolved as recommended:** `is_self_bound` returns **`E_int < Σ(n)`** (True =
  self-bound = shed enabled), `gate_margin_eV = E_int − Σ(n)`, and `evaporation_step` sheds
  only where `is_self_bound` is True.
- **2. `EvaporationResult` dataclass**, not the plan's "post-event tuple" — symmetric to the
  delivered `PickupResult`/`ShedResult`/`CaptureResult`; keeps `n_plus`/`m_plus_amu`
  "post-event" names though evaporation is a loss (`n_plus = n−1`, `m_plus_amu = m − m_He`).
- **3. Loose scalar kwargs, not a `state` container** (the P Q1 precedent): `evaporation_step(*,
  rng, E_int_eV, n, v, m_amu, nu, picture, kappa, dt_ps, evap_rrk_dof=None, …)`; the persistent
  `(n, E_int, v, m)` carrier is G's `IonStepState` (Phase C).
- **4. Build `evaporation_step_components` in Phase B** (the P Q2 precedent): the vectorized
  `(M,)`-ensemble form (one `rng.random(size=M)` draw, per-ion gate + `draw < P_shed` fire
  mask, `cold_shed_velocity_components` with per-ion mass since fires diverge + `dE_int_shed_eV`).
- **5. Unconditional Bernoulli draw for scalar↔components RNG parity.** The scalar form draws
  one uniform every step (`k=0`/`P_shed=0` under suppression → never fires but still consumes
  the draw), so scalar and components consume the RNG identically and the scalar stays the
  exact ion-by-ion oracle. A minor deviation from the plan's literal "if self-unbound return
  no-shed before drawing"; chosen for the Slice-X draw-order lock (outcome unchanged).
- **6. `NU_EVAP_PER_PS = 2.42`** added to `constants.py` as a **rate-primary float [ps⁻¹]** (not
  wavenumber-converted like the ladder anchors); `evap_rate_prefactor_per_ps` defaults to it
  (Q7, locked).
- **7. `evap_rrk_dof` override scope = `n≥2` bracket only** (`n=1` stays the direct `k=ν`
  branch); the `s≥1` load guard fires only when the override is set (None → per-`n`
  `effective_dof`, naturally ≥1, unguarded).
- **8. New `check_evaporation_config` guard** (mirrors `check_pickup_config`): the `evap_rrk_dof`
  `s≥1`-when-set check + the `gate_onset_override_eV` ↔ `allow_gate_onset_override` provenance
  refuse (mirrors the `allow_unvalidated_binding_pairing` refuse→warn arm, `config.py:846`).
  **Knob-name reconcile flagged:** MASS §R9 names it `evap_gate_onset_eV`; the plan uses
  `gate_onset_override_eV` (conveys the diagnostic-override intent). Plan name at build; MASS
  annotated not silently rewritten (the `p=κ` annotation precedent).

**Cross-reference verdict:** these are software-interface + one boolean-direction fix; no
MASS/CALIBRATION *value* is touched. The gate-direction correction *aligns* the helper with
MASS §R9 (self-unbound margin `G`, suppress while `G>0`) — the encoded physics was already
correct. ν pinned (Sourced), `s` Derived (row 10, n=2-linear), gate parameter-free (R9,
Derived). **Slice Q is build-ready** pending the `[PROCEED TO IMPLEMENTATION]` trigger.

---

## Phase B — Slice Q DELIVERED (2026-07-01) — Phase B COMPLETE

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN): the suites were written first and confirmed RED (missing
`NU_EVAP_PER_PS` / `check_evaporation_config` / `physics.evaporation`), then the module
made them GREEN. Pure-ish stochastic primitive: injected RNG, per-event mass/velocity
reset, **no** integrator, no `E_int` state persistence, no schema I/O, no 5-term closure
(all Phase C). Mass-agnostic in the drag sense — no Phase-Q function evaluates `γ`; `m`
enters only as the shed quantity in the reset. Oracle asserts against the plan §2.2 (Q) /
§4 golden values, composed against the **real** L (`d0_of_n`/`ladder_cumsum`), U
(`dE_int_shed_eV`), and `mass_jump.cold_shed` neighbours.

**Physics resolved at build (the §3.3 gate-direction fix carried through):** MASS §4
locks the shed band `D_0(n) < E_int < Σ(n)` — the upper bound is the self-bound gate
(`is_self_bound == E_int < Σ(n)`, True = self-bound = shed enabled), the lower is
RRK-bracket positivity. **`n=1` is the degenerate diatomic** (the band collapses,
`s=3n−3=0`), so it is the **direct-dissociation** special case with the *inverted* gate
`E_int > D_0(1)`, `k=ν` (MASS A11 boundary-of-validity). The fire path therefore keys on
`rrk_rate` (which encodes *both* the n≥2 self-bound gate and the n=1 direct gate);
`is_self_bound` / `gate_margin_eV` are the n≥2 diagnostics only (documented as **not** the
fire predictor at n=1).

**Build (files):**
- `physics/constants.py` — new **rate-primary** anchor `NU_EVAP_PER_PS = 2.42` [ps⁻¹]
  (2.42 directly, *not* wavenumber-converted like the ladder anchors; Sourced/pinned,
  MASS A11). Addition, not an edit to the MD constants table.
- `physics/evaporation.py` — **new module.** `EvaporationResult` dataclass (post-event
  names, a loss); `effective_dof(n)` (`n=2→4` linear, `3n−3` for `n≥3`, **fail-loud
  `ValueError` for `n<2`** — the direct branch owns n=1); `shed_probability` (`1−e^{−k dt}`,
  the shed counterpart of `pickup.attach_probability`); `gate_margin_eV`
  (`G=E_int−Σ(n)`); `is_self_bound` (`E_int < Σ(n)`, Python-bool scalar / bool-ndarray);
  `rrk_rate` (the full gated rate — n≥2 saturating RRK in-band + self-bound gate, n=1
  direct, `k=0` for n≤0; bounded `[0,ν)`; guarded so **no `0**0` / `0**neg` /
  divide-by-zero** is ever formed; module-level `s≥1` defense on the `evap_rrk_dof`
  override mirroring pickup's `p<0` guard); `evaporation_step` (**loose kwargs**,
  **unconditional draw**, composes `cold_shed` + `dE_int_shed_eV`, `n`=pre-shed →
  `EvaporationResult`); `evaporation_step_components` (vectorized `(M,)` ensemble, one
  `rng.random(size=M)` draw, per-ion gate + `draw < P_shed` mask). Optional diagnostic
  `gate_onset_eV` threading (n≥2 gate only; n=1 keeps its `D_0(1)` onset).
- `physics/mass_jump.py` — **`cold_shed_velocity_components` generalized to per-ion mass**
  (the symmetric counterpart of the Slice-P `capture_velocity_components`): its scalar
  `_check_masses` guard swapped for a new array-aware `_check_masses_shed`; the arithmetic
  (`m+`, `kick`, `_reduced_mass_defect_coeff`) already broadcasts, so **the scalar path is
  byte-identical** (Tier-1a contract; all delivered mass_jump/reset tests stay green). This
  is the honest resolution of the §3.3 decision-#4 "per-ion mass, fires diverge" note (the
  function previously took a *uniform* scalar mass; one cold-shed-components function now
  serves both the Tier-1a uniform schedule and the Tier-2 per-ion ensemble — rule 1).
- `config.py` — new fields `evap_rate_prefactor_per_ps` (ν, defaults to the
  `NU_EVAP_PER_PS` anchor — single source), `evap_rrk_dof` (`Optional=None`),
  `gate_onset_override_eV` (`Optional=None`), `allow_gate_onset_override` (`False`); new
  `check_evaporation_config` guard (the `evap_rrk_dof` `s≥1`-when-set check + the
  `gate_onset_override_eV` ↔ `allow_gate_onset_override` **refuse→warn** provenance guard
  mirroring `allow_unvalidated_binding_pairing`) wired into `validate()`.
- `tests/` — `test_evaporation.py` (42 + 8 hardening = 50) + `test_evaporation_config.py`
  (16); `docs/config_and_preset.md` deferred-field list extended with the three Slice-Q
  driver-read fields.

**Oracle cross-validation.** `effective_dof` = 4/6/60 at n=2/3/21, fail-loud <2;
`shed_probability` = `1−e^{−k dt}` (→`k dt` small, saturates to `1.0`); `rrk_rate` `=0`
below threshold and while self-unbound, `∈(0,ν)` strictly across the whole in-band sweep
for n∈{2,3,7,15,21} (no avalanche), monotone↑ in `E_int`, `=ν` for the n=1 hot direct
branch and `=0` when cold; the `evap_rrk_dof` override applies to n≥2 only (n=1 stays
direct) and `s<1` fails loud; `gate_onset_eV` shifts the n≥2 gate but never n=1;
`gate_margin_eV`/`is_self_bound` sign + zero-crossing; the scalar step's fire path matches
`cold_shed` (v/m/defect) + `dE_int_shed_eV` (`−D_0(n)`) exactly, no-fire is a defensive
copy, suppressed/below-threshold never fire even under a fire-forcing RNG, one draw
consumed **unconditionally**; the components form matches the scalar oracle ion-by-ion
(all-fire), leaves non-fired ions untouched, draws once (size M) with **RNG consumption
parity** to the scalar loop, and reproduces the seeded Bernoulli fire fraction within 5σ
over 200 000 ions. The per-ion `cold_shed_velocity_components` matches the scalar
`cold_shed` ion-by-ion and fails loud on an underweight ion.

**Rule-2 table.** `evap_rate_prefactor_per_ps` (ν), `evap_rrk_dof` (s), and
`gate_onset_override_eV` are the three Slice-Q **declared-but-unread** carries (read by the
Phase-C driver / evaporation module; `allow_gate_onset_override` and the `s≥1` guard are
born-live via `check_evaporation_config`). They join the standing Phase-B carries
(`sweeping`/`dwell_time`/`thermal` pickup arms, the sourced TDDFT `ρ_He` array) and the
Slice-L `cooling_relaxed` concrete blend (pinned at Phase F).

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. ν pinned
(Sourced, A11), `s` Derived (row 10, n=2-linear applied), the self-bound gate
parameter-free (R9, Derived), the override a guarded R10 diagnostic-lever. No `γ`/drag-law
read, no integrator, no `E_int` state — the out-of-scope guard holds.

**Tests:** Slice-Q suites **66 green** (`test_evaporation.py` 50 + `test_evaporation_config.py`
16); full suite **1725 passed, 0 failed** (1659 → +66; same 17 unrelated `anchored_discrete`
warnings). **Phase B (ρ + P + Q) is complete.** Next: Phase C (Slice X checkpoint v6→v7 +
5-term invariant, then Slice G generative driver) — `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md`.

### Slice Q — extended review + test pass (2026-07-01)

An extended adversarial review (correctness / edge / vectorisation / independence / guard /
RNG-parity / doc-drift angles over the self-authored diff, plus a re-audit of the touched
`mass_jump` function's callers) **confirmed the physics with no correctness bug** and
surfaced **one genuine dead-code / wrong-comment item**; fixed, plus +10 deeper locks.

- **Dead exponent clamp removed (rule 2).** `rrk_rate` computed the RRK bracket exponent as
  `np.maximum(s - 1.0, 0.0)`, but `s` is **always ≥ 1** — per-`n` `effective_dof(max(n,2)) ≥ 4`,
  and the override is guarded `s ≥ 1` in the same function — so the clamp **never activated**,
  and its comment misdescribed *why* (it claimed the n<2 lanes get a negative `s`, but they get
  `effective_dof(2)=4` via the `max(n,2)` clamp). Simplified to `base_safe ** (s - 1.0)`; the
  `base > 0` mask remains the sole (sufficient) protection against `0**0`/`0**neg`, and `E_safe`
  keeps the division safe. Comment rewritten to state the true invariant. No behaviour change
  (the clamp was provably inert) — full suite unchanged.
- **Caller re-audit (mass_jump generalization).** The `cold_shed_velocity_components`
  per-ion-mass generalization has **no production callers** (the driver is Phase C) — only
  `tests/test_mass_jump.py` and `tests/test_ion_variable_mass.py`. The latter's guard-reject
  test (`match="m_he|pre-shed"`) still passes: `_check_masses_shed`'s message carries "pre-shed",
  and the scalar-mass reject path is unchanged. Scalar byte-identity is now **explicitly**
  locked (uniform-array == scalar-mass result) in addition to the pre-existing Tier-1a tests.
- **+10 extended locks (all green, several bug-catching):** RNG-stream parity **with a
  suppressed ion** (the load-bearing unconditional-draw property — a suppressed ion still
  consumes its components draw, matching the scalar loop); non-finite `E_int` characterisation
  (NaN → `k=0` no crash; `inf` → shell suppressed / `n=1` direct fires); `n=0` has no rung
  (`k=0`, scalar + vectorized); negative `n` fails loud via `ladder_cumsum`; the `s=1` override
  edge (flat `k=ν` in-band, still `0` below-threshold / suppressed); the override **respects the
  gate** (never un-suppresses); vectorized gate diagnostics with the `gate_onset_eV` override;
  the fire booking **two distinct channels** (K1 eV drain vs the mechanical cold-shed defect,
  no double-count); and the scalar↔uniform-array `cold_shed_velocity_components` equality.
- **Deliberate non-guards (consistent with pickup/drag):** non-finite `E_int`/`v` is treated as
  an upstream bug (characterised → `k=0`/no-fire, not guarded); `evaporation_step` is documented
  single-ion (an array `E_int`/`n` would trip `bool(...)` — a misuse, not a supported path); the
  components form trusts caller array-length parity (parity with `pickup_step_components`).

**Tests after the review pass:** Slice-Q suites **76 green** (`test_evaporation.py` 60 +
`test_evaporation_config.py` 16); full suite **1735 passed, 0 failed** (1725 → +10; same 17
unrelated `anchored_discrete` warnings). **Phase B (ρ + P + Q) is complete.** Next: Phase C
(Slice X checkpoint v6→v7 + 5-term invariant, then Slice G generative driver) —
`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md`.

---

## Phase C — Slice X pre-build decisions (2026-07-01)

A pre-build discussion pass over the Slice X spec (checkpoint **v6→v7** + `IonStepState`
extension + the **5-term** invariant closure) in `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md`,
grounded in a code-surface audit of `simulation/checkpoint.py`
(`_ION_SCHEMA_VERSION`, `_migrate_ion_checkpoint`, `_validate_against_cfg`
`trajectory_2N_T_fields`), `postprocess/energy_balance.py` (`EnergyTotals` /
`ion_energy_totals` / `ion_ledger_closure`), and `simulation/ion_propagation_step.py`
(`IonStepState`, the `ion_state_from_checkpoint_column` / `write_ion_state_to_checkpoint_column`
seams). No code — the `[PROCEED TO IMPLEMENTATION]` boundary holds; docs-only amendment. Three
design calls resolved (decision owner: user) plus five code-audit findings.

**Five code-audit findings.**
- **F1 — Slice X necessarily touches the fixed/anchored driver seams.**
  `ion_state_from_checkpoint_column` (`ion_propagation_step.py:636`) and
  `write_ion_state_to_checkpoint_column` (`:653`) are the read/write seams the *current*
  (non-biphasic) driver uses, and `IonCheckpoint` is allocated in `ion.py`. Once `E_int_eV`
  is a v7 field a `fixed` / `anchored_discrete` run must still emit a **valid all-zero v7**
  (plan §7 acceptance). So X is not purely additive on the dataclasses — it threads a zeros
  array through both seams + the driver allocation. The `E_int` *evolution* stays in G; only
  the zeros plumbing is X.
- **F2 — the shim is single-hop today.** `_migrate_ion_checkpoint` handles only `version==5→6`
  and is checked once; v7 needs either a cascade or a v7-only arm. Plan §2.2 left this "decide
  at build".
- **F3 — two distinct "no real E_int" representations.** Neutral → `None` (no field); a
  migrated v6-origin ion file → **zeros array** (not `None`). So `ion_energy_totals` sums
  `E_int_eV` unconditionally for any ion checkpoint (always present post-shim); the `None`
  branch is neutral-only. The plan's "None for a v6-origin run" wording is corrected here to
  **zeros, not None**.
- **F4 — `n_shell` "genuine state" is a G concern, not X.** The field already exists at v6; X
  changes nothing about it (the v6 derivation / writer rule stays intact).
- **F5 — draw-order lock is docs-only in X.** No driver exists yet to enforce it; CLAUDE.md
  already forbids "changing random-number draw order" generically, so X only *records* the
  specific shed-then-pickup order in `DRAG_PORT_DESIGN_DECISIONS.md` + this log.

**Three design calls resolved (user, 2026-07-01).**
1. **v5 handling → chain v5→v6→v7.** `_migrate_ion_checkpoint` is restructured to apply
   migrations **stepwise** (the existing v5→v6 arm returns 6, a new v6→v7 arm then runs) so all
   surviving Tier-0/1a v5 **and** v6 `ion.npz` stay loadable. Chosen over v7-only/reject-v5
   (which would strand pre-Tier-1a artifacts) — cost is trivial.
2. **`E_int_eV` discipline → required field, wire zeros now.** `E_int_eV` lands as a
   **required** field on both `IonStepState` (`(2N,)`) and `IonCheckpoint` (`(2N,T)`), mirroring
   `E_mass_transfer_eV` (no `Optional`/`None` debt). X updates the two column seams + the
   `IonCheckpoint` allocation in `ion.py` so `fixed` / `anchored_discrete` runs emit a valid
   all-zero v7 — satisfying the plan §7 acceptance ("fixed runs produce an all-zero `E_int_eV`")
   literally. Chosen over `Optional=None` (defers plumbing to G but scatters None-checks through
   `energy_balance` + the seams).
3. **Migration provenance → load-time `warnings.warn`.** The v6→v7 arm synthesizes
   `np.zeros_like(E_mass_transfer_eV)` and emits a **warning** that the file predates the
   `E_int` reservoir — faithful to the plan's "flag" without an unplanned persisted schema
   field. (The v5→v6 arm stays silent, as delivered.) Chosen over silent-zeros (drops the
   provenance signal) and a persisted marker field (unplanned v7 schema addition).

**Resulting Slice-X build spec.**
- `checkpoint.py`: `_ION_SCHEMA_VERSION 6→7`; add `E_int_eV (2N,T)` to `IonCheckpoint`; cascade
  shim (v5→v6 arm returns 6, new v6→v7 arm synthesizes zeros + warns); add `E_int_eV` to the
  `trajectory_2N_T_fields` shape-check tuple.
- `energy_balance.py`: `E_int_eV` added to `EnergyTotals` (`None` = neutral only);
  `ion_energy_totals` sums it into `E_system_eV` (per-molecule `/n`, additive `+` per MASS §6
  eq. line 935); `ion_ledger_closure` reports the 5-term residual with unchanged machinery
  (extend-in-place, the `E_mass_transfer` precedent).
- `ion_propagation_step.py`: `IonStepState.E_int_eV` required; both column seams read/write it.
- `ion.py`: allocate the zeros `E_int_eV` array in the checkpoint assembly.
- Docs: freeze the shed-then-pickup RNG draw order (forbidden-list item) in
  `DRAG_PORT_DESIGN_DECISIONS.md` + this log.
- Tests: v7 round-trip; v5→v6→v7 and v6→v7 migration (zeros + warning asserted); missing
  `E_int_eV` on a genuine v7 raises; wrong-shape `E_int_eV` rejected; 5-term closure flat on a
  balanced synthetic stream / divergent on a miswire / neutral-`None` unaffected; the `fixed`
  regression (all-zero `E_int_eV`, 4-term residual unchanged).

**Boundary holds — no code.** Next: the MASS / CALIBRATION_MAP cross-check of the Phase-C plan
(below), then the Slice-X build under the `[PROCEED TO IMPLEMENTATION]` trigger (X before G, §5).

### Phase C — MASS / CALIBRATION_MAP cross-check of the Phase-C plan (2026-07-01)

A cross-reference pass over `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` (X + G) against
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §6/§7/§11 + A13 + §8 R2/R4 and
`CALIBRATION_MAP.md` (parameter classes + row 169 A13). No code — docs-only. **Verdict: the
plan is faithful on all physics / schema / config substance; Phase C adds no new calibration
parameter** (the four integrator flags are non-parameters per CALIBRATION row 169). Verified:
the 5-term invariant (MASS §6 eq. line 938), the reduced-mass capture defect (§6 943–951, *not*
`½m_He·v²`), the per-channel closure table (drag / pickup S1 / cold-shed K1 / K2), the four
integrator flags (§11 1946–1958), jump-then-O / one-event shed-then-pickup (A13 1786–1789), the
`§6.5` pairing guard + `s≥1` guard, and the R2 (constructed-reservoir, HIGH, bounded-by-ladder)
/ R4 (double-count) framings — all match line-by-line.

**Four discrepancies surfaced; three actioned as NB annotations (annotation-not-rewrite, the
`p=κ` precedent — locked values/mechanism untouched), one no-action.**
- **D1 (version slot) — annotated.** MASS §7 line 1132 ("`IonCheckpoint → v6`") and §6.11 line
  983 ("already in the v6 schema") predate the Tier-1a v6 use of the slot; the live build bumps
  **v6→v7** for `E_int_eV`. NB added at both sites pointing to the Phase-C plan §2.2/§8 + this
  log (version numbers only; mechanism unchanged). Already flagged as a reconcile in the plan §8.
- **D2 (scenario literal) — annotated.** MASS §11 line 1920–1921 still names
  `biphasic_energy_gated` as the production literal with `scenario_A_accretion` /
  `scenario_B_stripping` / `biphasic` as legacy; the delivered `MassScenario` set is
  **`{fixed, biphasic, anchored_discrete}`** (production value = `biphasic`; the A/B scenarios
  were superseded by `anchored_discrete`, never implemented). NB added at §11 (main) and §7
  (scenario-metadata bullet, cross-pointer). Config literal names only; §11 knobs unchanged.
- **D3 (Tier-2 observable) — annotated.** CALIBRATION "Validation tiers" Tier-2 line still reads
  "size distribution (+ per-fragment velocity histograms)"; Phase E **cut** the velocity
  comparison (aggregate `vmi_summary/*.csv`, not per-n). NB added pointing to
  `PHASE_E_IMPLEMENTATION_PLAN.md`.
- **D4 (R4 wording) — no action.** The plan calls R4 the "chief energy risk" while §8 tags it
  MEDIUM, but §6 line 933 itself calls it "the chief energy-balance risk" — the plan is faithful
  to MASS's own framing.

**Files touched (annotations only):** `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §7
(D1 + D2 pointer), §6.11 (D1), §11 (D2); `CALIBRATION_MAP.md` Tier-2 validation line (D3). No
code; boundary holds. **Next:** the Slice-X build under the `[PROCEED TO IMPLEMENTATION]`
trigger (X before G, §5).

---

## Phase C — Slice X DELIVERED (2026-07-01)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN). Pure schema + closure extension — **no RNG, no integrator, no new
channel physics**; composes nothing stochastic (the pre-build spec, §2.2/§3). The three
Slice-X design calls (v5→v6→v7 cascade; `E_int_eV` **required** field, wire zeros now;
load-time `warnings.warn` provenance) and the five code-audit findings (F1–F5) carried
straight through to the build.

**Build (production, 4 files):**
- `simulation/checkpoint.py` — `_ION_SCHEMA_VERSION 6→7`; **add required
  `E_int_eV (2N,T)`** to `IonCheckpoint` (placed after `E_mass_transfer_eV`, before the
  defaulted `mass_scenario`/`schema_version`); `_migrate_ion_checkpoint` **restructured to a
  stepwise cascade** — the v5 arm now *falls through* (no early return) into a new v6→v7 arm
  that synthesizes `E_int_eV = zeros_like(E_mass_transfer_eV)` and emits a `UserWarning`
  (the v5→v6 arm stays silent, as delivered); `E_int_eV` added to the
  `trajectory_2N_T_fields` shape-check tuple. `import warnings` added. Docstring updated
  (field bullet + the v6/v7 schema paragraph).
- `postprocess/energy_balance.py` — **5-term closure extended in place.** `E_int_eV`
  added to `EnergyTotals` (`None` = neutral only); `ion_energy_totals` sums
  `e_int = Σ E_int_eV / n` into `E_system_eV` (per-molecule `/n`, additive `+`, MASS §6
  eq. line 935); `neutral_energy_totals` passes `E_int_eV=None`; `ion_ledger_closure`
  machinery **unchanged** (reads the now-5-term `E_system`). `LedgerClosure` /
  `ion_ledger_closure` docstrings updated four-term→five-term.
- `simulation/ion_propagation_step.py` — **required `E_int_eV (2N,)`** on `IonStepState`
  (after `E_mass_transfer_eV`); both column seams read/write it
  (`ion_state_from_checkpoint_column`, `write_ion_state_to_checkpoint_column`); the two
  fresh `IonStepState(...)` constructions (`baoab_propagation_step`, the Tier-0 drag step)
  carry `E_int_eV=state.E_int_eV` through untouched. `shed_step` uses `replace(...)` so it
  carries the field for free. Docstring updated.
- `simulation/ion_initial_state.py` — allocate the all-zero `E_int_eV (2N,T)` in the
  `build_initial_ion_state` checkpoint assembly (the plan's "ion.py allocation"; the actual
  production builder is `ion_initial_state.py`).

**Design calls as-built (unchanged from the 2026-07-01 pre-build record):**
1. **v5→v6→v7 cascade.** The shim reads the stored `schema_version` (via `_load_checkpoint`,
   called only when `version != expected`), so a v5 file walks v5→v6→v7, a v6 file runs the
   v6→v7 arm only, and a genuine v7 skips migration → a v7 file missing `E_int_eV` hits the
   `missing fields` raise (never zero-filled). The v5→v6 output stays byte-identical
   (regression-locked by the pre-existing v5-shim + Tier-1a reset tests).
2. **`E_int_eV` required (no `Optional`/`None` debt).** Mirrors `E_mass_transfer_eV`; cascaded
   to **every** `IonCheckpoint` / `IonStepState` construction site (1 production builder +
   ~17 test builders, incl. one dict-style `bad_dict` the kwarg grep missed).
3. **Load-time `warnings.warn` (UserWarning).** The v6→v7 arm flags synthesized zeros; the
   v5→v6 arm silent. A v6-origin file is 4-term-equivalent (all-zero reservoir), so the
   `EnergyTotals.E_int_eV` for a migrated file is **zeros, not `None`** (F3) — the `None`
   branch is neutral-only.

**RNG draw-order lock (F5, docs-only in X).** The shed-then-pickup order is **frozen** as a
forbidden-list item in `DRAG_PORT_DESIGN_DECISIONS.md` §2.9 (Slice-X DELIVERED block).
Enforcement lands with the Slice-G driver (none exists yet); X only records it.

**Oracle cross-validation.** v7 round-trips `E_int_eV` bit-for-bit; v6→v7 and v5→…→v7
synthesize all-zero `E_int_eV` + warn (`pytest.warns(UserWarning, match="E_int")`); a genuine
v7 missing the field raises `missing fields`; a wrong-shape `(N,T)` `E_int_eV` is rejected by
`_validate_against_cfg`. 5-term closure: flat residual (≤1e-12) on a balanced stream, reduces
to the 4-term baseline when `E_int≡0`, an `E_int` deposit offset by an equal `E_dissip` drain
closes, and an **un-offset** `E_int` deposit diverges (0.6 eV/molecule, caught in five terms).
`fixed` Tier-0 run round-trips a present, all-zero `E_int_eV` (`test_ion_drag_smoke`).

**Doc drift fixed.** Phase-C plan §3 F3 line (v6-origin → **zeros, not None**);
`test_ion_drag_smoke::test_checkpoint_v5_round_trip` assert `6→7` + an all-zero `E_int_eV`
check (it round-trips a freshly-produced current-schema checkpoint, not a real v5 file).

**Rule-2 table.** Slice X adds **no config field** (schema/closure only) — nothing to add to
or remove from the exception table. `_ION_SCHEMA_VERSION = 7`.

**New warnings (expected, not failures).** 8 `UserWarning`s now fire in the suite where real
**v6-origin** reference `ion.npz` fixtures are loaded (`test_compare_trajectories` real-9A,
the legacy-debug / paper_cov plot smokes) — the v6→v7 provenance shim working as designed. The
17 pre-existing `anchored_discrete` mass-pairing `RuntimeWarning`s are unchanged.

**Tests:** new coverage = `test_checkpoint.py::TestIonSchemaV7` (4) + the extended v5-shim
assert + `test_energy_balance.py::TestLedgerClosure` 5-term (3) + the totals-sum + smoke
assert. Full suite **1742 passed, 0 failed** (was 1735; +7 net), 25 warnings (17 pre-existing
`anchored_discrete` + 8 new v6→v7 shim). **Next:** Slice G — the `biphasic_step` generative
driver (`simulation/ion_propagation_step.py`), composing A {K,U} + B {ρ,P,Q} + the Tier-1a
seam behind the now-ready 5-term invariant (`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §3 Slice G).

### Slice X — review + test-hardening pass (2026-07-01)

An extended adversarial review (8 angles — line-by-line diff scan / removed-behavior /
cross-file tracer / reuse / simplification / efficiency / altitude / CLAUDE.md-conventions,
run **inline** over the small self-authored diff per the delivered-slice precedent)
**confirmed the physics/schema with no correctness bug** and surfaced only coverage gaps;
all closed test-first (characterisation/regression locks).

- **No production bug, no code change.** The v6→v7 bump is a faithful in-place mirror of the
  Tier-1a `E_mass_transfer_eV` precedent. Verified: the migration cascade preserves the
  v5→v6 sub-transform byte-for-byte (falls through instead of returning), a genuine v7
  missing `E_int_eV` raises rather than zero-fills (migrate never runs at `version==expected`),
  the two step-state carry-throughs (`E_int_eV=state.E_int_eV`, no copy) match the baoab
  `E_mass_transfer` precedent and are read-only within the step chain (the column read seam
  `.copy()`s), and the drop of the early `if version != 5: return` re-establishes the
  "unknown version → returned for the strict check to raise" invariant via the cascade tail.
- **Cross-file surface clean (verified).** `EnergyTotals` is only ever **keyword**-constructed
  (both sites in `energy_balance.py`) — the new defaulted `E_int_eV` field breaks no positional
  caller; **no** `IonCheckpoint(` / `IonStepState(` / `EnergyTotals(` construction exists under
  `scripts/` (production builds go through `build_initial_ion_state`), so the required-field
  cascade is fully covered by the 1 production + ~18 test builders already updated.
- **Coverage gaps closed (+3 locks):** (1) the per-molecule `/n` division of `E_int` in
  `ion_energy_totals` was only exercised **all-zero** (the closure tests zero `E_int`) → a
  nonzero value-lock on `totals.E_int_eV == ΣE_int/n` + its inclusion in `E_system`; (2) the
  two v7 **column seams** (`ion_state_from_checkpoint_column` /
  `write_ion_state_to_checkpoint_column`) had **no direct test** (only the all-zero driver
  path) → a nonzero round-trip lock on an **isolated** checkpoint copy (`replace(ion,
  E_int_eV=…copy())`, so the module-scoped fixture is never polluted) that also asserts the
  read-side `.copy()` independence; (3) `neutral_energy_totals` `E_int_eV is None` lock.
- **Deliberate non-finding.** `warnings.warn(..., stacklevel=2)` points at the
  `_load_checkpoint` migrate call rather than the user's `load_ion_checkpoint(...)` (4 wrapper
  frames deep); left as-is — a fixed larger stacklevel is brittle (breaks if `_load_checkpoint`
  is called directly), the delivered v5→v6 arm set no stacklevel at all, and the category +
  message are the load-bearing signal (asserted by `pytest.warns(UserWarning, match="E_int")`).
- **Deliberate non-guards (consistent with the delivered schema code):** the v6→v7 arm skips
  synthesis when `E_mass_transfer_eV` is absent (a corrupt v6) and lets the `missing fields`
  check raise, rather than KeyError-ing inside `np.zeros_like`; non-finite `E_int` is an
  upstream-bug concern, not guarded at the schema layer.

**Tests after hardening:** narrow suites (`test_energy_balance` + `test_ion_propagation_step` +
`test_checkpoint`) **62 green**; full suite **1745 passed, 0 failed** (1742 → +3), same 25
warnings (17 pre-existing `anchored_discrete` + 8 v6→v7 shim). **Next:** Slice G — the
`biphasic_step` generative driver (`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §3 Slice G).

---

## Phase C — Slice G pre-build decisions (2026-07-01)

A pre-build discussion pass over the Slice G spec (`biphasic_step` generative driver) in
`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md`, grounded in a full code-surface audit of every
composed surface: `simulation/ion_propagation_step.py` (`IonStepState`, `shed_step`,
`baoab_propagation_step`, the two column seams), `simulation/ion.py` (the `use_drag`
dispatch loop + the per-step BAOAB-closure rebuild + `_drag_gate_steepness`),
`simulation/ion_initial_state.py` (`build_initial_ion_state` mass/`n_shell`/E_int column-0
fills), `physics/baoab.py` (`make_ion_baoab_step` per-atom-mass O-step),
`physics/solvation_cooling.py` (`newton_cool_step`), `physics/pickup.py` /
`physics/evaporation.py` (the `*_components` ensemble forms + RNG contracts),
`physics/internal_energy_budget.py` (S1/S2/K1), and `config.py` (§6.5 guard,
`_EVOLVING_MASS_SCENARIOS`, the None-/0.0-defaulted knobs). No code — the
`[PROCEED TO IMPLEMENTATION]` boundary holds; docs-only amendment. Four design calls
resolved (decision owner: user) plus the audit-settled composition facts recorded below.

**Audit-settled (no fork — consistent with plan + delivered code):**
- **Dispatch.** Add a `biphasic` arm in `ion.py`'s `use_drag` block beside the
  `anchored_discrete`/`shed_step` arm; `biphasic` already sits in
  `_EVOLVING_MASS_SCENARIOS` → trips §6.5 → runs under `allow_inconsistent_mass_pairing=
  True`. The arm threads `rng` (the channel draws) — the current drag path does not.
- **Scope guard.** `_check_drag_scope` currently *rejects* `biphasic`; G flips it to accept
  (the `noise_form='none'` check stays; the `m_eff` band trip-wire already `return`s for
  non-`fixed`, so evolving biphasic mass is exempt like `anchored_discrete`).
- **Per-atom mass already supported.** `make_ion_baoab_step` takes a `(2N,)` mass vector and
  the O-step divides by `m(t)` per atom, so **non-uniform** per-ion biphasic mass integrates
  unchanged (unlike the uniform-mass `anchored_discrete` schedule).
- **Ensemble granularity.** `pickup_step_components` / `evaporation_step_components` run over
  the `(2N,)` ion ensemble, one unconditional `rng.random(size=2N)` draw each.
- **Units seam (G owns).** The channels' `dE_mass_transfer` is **amu·Å²/ps²**; the
  `E_mass_transfer_eV` field is eV → G converts with the delivered `×U×100²/EV` idiom (the
  `shed_step`/`baoab_propagation_step` path). Channel `dE_int` is already eV.
- **Cooling ≡ E_int decay.** `newton_cool_step` returns `E_solv.struct = E_∞(N) + E_int`, so
  at fixed N it is exactly `E_int·e^{−dt/τ}`; the drained `E_int·(1−decay)` books to
  `E_dissip` (§6 closure row K2). G **reuses `newton_cool_step`** (rule 1: single source of
  the cooling form), not an inline decay.
- **Ordering that the tests lock.** Cooling runs **before** the evaporation draw — the
  self-bound gate opens (`t×`) precisely as cooling drains `E_int` below `Σ(n)`. Forward
  `E_int` evolution is **additive** (`E_int += Σ dE_channel` after the K2 decay);
  `reconstruct_e_int_eV` (A9) stays a post-t× diagnostic, never called in the loop.
- **Draw order (frozen at X).** Evaporation draw, then pickup draw; if both fire on an ion →
  apply **shed then pickup** (pickup's rate reads the post-shed `n`; net `n` unchanged; both
  E_int/E_mass_transfer bookings applied).

**Four design calls resolved (user, 2026-07-01):**
1. **Integrator placement → pre-step seam (reuse), NOT in-step B/A→jump→O.** `biphasic_step`
   mirrors the Tier-1a `shed_step` pattern: K2 cooling + ≤1 mass event (shed-then-pickup) +
   the `E_int` update happen at the **step seam**, then the BAOAB closure is rebuilt at `m⁺`
   so the conservative kicks *and* the SQ1 O-step both read the post-jump mass. Reuses
   delivered code, leaves `baoab.py`/`make_ion_baoab_step` untouched, O(dt)-benign (the
   `shed_step` docstring's jump-measure-zero argument). **Resolves the plan §2.3 numbered-list
   tension:** that list is *logical* (the composed operations), not a literal in-BAOAB
   ordering — the seam does cooling+event *before* the B/A kicks, not between B/A and O.
2. **`n_shell` → genuine state on `IonStepState`.** Add `n_shell: np.ndarray (2N,)` to
   `IonStepState`; the P/Q `*_components` channels advance it (`n→n±1`) **independently** of
   the mass reset; the writer stores `state.n_shell` directly (stops deriving it via `rint`
   for the biphasic path); and `biphasic_step` asserts `m == m_I⁺ + n·m_He` (within float
   tol) per step — the **deferred Q3 guard** (Phase-B Slice P) now has real teeth (independent
   `n` to catch a channel/reset drift, not a `rint`-tautology). Touches the two column seams +
   `build_initial_ion_state` + `write_ion_state_to_checkpoint_column`. The `fixed`/
   `anchored_discrete` paths keep `n_shell` derived-from-mass at write (byte-identical).
   *(NB: `n_shell` is already a **checkpoint** field at v6 — this adds it only to the
   in-memory `IonStepState` carrier; no checkpoint-schema bump, no forbidden-list trip.)*
3. **Biphasic initial state → n₀=21, onset seeded in the builder.** `build_initial_ion_state`
   gets a `biphasic` branch: start at the full first shell **n₀ = `ANCHOR_N_START` = 21**
   (`mass = complex_mass_amu(21)`, matching the 21→19→14 validation target and the
   `anchored_discrete` start), and deposit the **S2 onset** `E_int[:,0] = f_int ·
   coulomb_available_eV` in the builder (it already owns column-0 physics: E_kin/E_pot/mass/
   n_shell). `E_avail ← cfg.coulomb_available_eV` (0.80 eV validation / 2.70 eV production).
   The onset is deposited **once** at t=0; only S1/K1/K2 modify `E_int` thereafter.
4. **Missing-knob guard → config-load `check_biphasic_config`.** A new `validate()` guard:
   when `mass_scenario=='biphasic'`, require `internal_energy_partition_fraction` (f_int) and
   `internal_energy_retained_fraction` (f_ret) **non-None** (both are None-defaulted and the
   onset/S1 cannot be computed without them) and flag `pickup_rate_coefficient` (λ₀) `> 0`
   (0.0-default = pickup structurally inert — warn or require, decide at build). Mirrors
   `check_pickup_config` / `check_evaporation_config`; fails at config-load (earliest), not
   deep in the driver. κ/picture/τ/ν/p/cap all carry live defaults, so they are not part of
   the required set.

**Config field deltas (Slice G):** the three A13 integrator flags land + activate here
(`one_mass_event_per_step`, `jump_o_step_ordering`, `mass_jump_velocity_reset`; plan §4);
`he_capture_velocity` (declared at Slice P) is *activated* (G passes `u_he`). The Slice-P/U
declared-but-unread carries (`pickup_rate_coefficient`, `pickup_occupancy_exponent`, f_int,
f_ret, `evap_rate_prefactor_per_ps`, `evap_rrk_dof`, `gate_onset_override_eV`) are **read by
the driver here** → removed from the rule-2 exception table at the G build.

**Cross-reference verdict:** these are software-composition + integrator-policy calls; no
MASS/CALIBRATION *value* is touched. Decision #1 (seam) is faithful to A13 (jump-then-O,
one-event-per-step) and the plan §0 locked "reuse the `shed_step` seam" decision; #2 realizes
the Phase-B Q3 deferred guard; #3 matches the validation-first 0.80 eV / n=21 target; #4
mirrors the existing channel config guards. **Slice G is build-ready** pending the
`[PROCEED TO IMPLEMENTATION]` trigger (the last slice of Phase C).

### Phase C — Slice G two remaining build calls resolved (2026-07-01, user)

A second pre-build pass (full code-surface audit of the composed modules) closed the two
items design call #4 / plan §4 left open. No code — the `[PROCEED TO IMPLEMENTATION]`
boundary holds; docs-only amendment.

- **5. λ₀ guard → warn, do not require.** `check_biphasic_config` **hard-requires**
  `internal_energy_partition_fraction` (f_int) and `internal_energy_retained_fraction`
  (f_ret) non-None when `mass_scenario=='biphasic'` (onset/S1 undefined without them), and
  emits a **load-time warning** — not a raise — when `pickup_rate_coefficient` (λ₀) `== 0.0`
  (pickup structurally inert). Rationale: an evaporation-only biphasic run seeded at n₀=21
  is a legitimate limit (the Phase-E relaxation stage is its close cousin), so λ₀=0 is not a
  config error; the validation-first bridge (Phase D) still uses λ₀>0. Fail-loud is reserved
  for the genuinely-undefined knobs (f_int/f_ret), advisory for the merely-inert one.
- **6. The three A13 integrator flags → fixed driver policy, NOT config toggles.** Plan §4
  listed `one_mass_event_per_step` (=true), `jump_o_step_ordering` (=jump_then_O), and
  `mass_jump_velocity_reset ∈ {momentum_conserving, label_only}` as fields landing at G. Each
  has exactly **one built/honored value** (`label_only` has no primitive — the Phase-B resets
  are momentum-conserving by construction), so live toggles would be **dead one-valued config
  surface** (rule 2). Resolution: encode all three as **structural policy** — the
  `biphasic_step` seam does one-event-per-step (evaporation-then-pickup, shed-then-pickup),
  jump-then-O (rebuild BAOAB at m⁺), and momentum-conserving resets **by construction**;
  documented in the driver docstring + `DRAG_PORT_DESIGN_DECISIONS.md` and **asserted by the
  Slice-G tests**, with **no new SimConfig fields**. This refines plan §4 (which pre-supposed
  the fields); the falsification-lever knobs that *do* land are the existing physics ones
  (κ, picture, |S|, τ, f_int, f_ret, λ₀, p, ν, s, cap/form, `gate_onset_override_eV`).
  Consequently the rule-2 carries removed at the G build are only the Slice-P/U physics
  fields (`pickup_rate_coefficient`, `pickup_occupancy_exponent`, f_int, f_ret,
  `evap_rate_prefactor_per_ps`, `evap_rrk_dof`, `gate_onset_override_eV`) — no integrator
  field is added or removed.

**Implementation note (no decision needed).** `write_ion_state_to_checkpoint_column` will
branch on `mass_scenario`: store `state.n_shell` directly on the `biphasic` path (genuine
state) and keep the `rint`-from-mass derivation for `fixed`/`anchored_discrete` (byte-identical
regression). With both calls resolved, Slice G is fully specified and build-ready.

---

## Phase C — Slice G DELIVERED (2026-07-02) — Phase C COMPLETE

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD throughout
(test→RED→GREEN), composing only accepted Phase-A/B modules + the Tier-0/1a seam; **no new
channel physics**. The four pre-build design calls + the two build calls (λ₀ warn; integrator
flags as fixed policy) carried straight through.

**Build (production, 4 files):**
- `simulation/ion_propagation_step.py` — new **`biphasic_step`** pre-step seam (K2 cooling →
  evaporation-then-pickup draws → S1/K1 `E_int` + `E_mass_transfer` (eV via `×U×100²/EV`) +
  bath/cooling `E_dissip` bookings → the `m == m_I⁺ + n·m_He` assert); `IonStepState` gains a
  defaulted `n_shell: np.ndarray | None = None` (genuine state on biphasic, `None` elsewhere)
  carried through the two return constructors + the read seam; `write_ion_state_to_checkpoint_
  column` gains a `mass_scenario` kwarg branching biphasic (store `state.n_shell`) vs the
  `rint`-from-mass derivation (byte-identical); `_check_drag_scope` flipped to **accept
  `biphasic`** (noise-`none` still required; m_eff trip-wire already exempts non-`fixed`).
  Module const `_M_N_CONSISTENCY_TOL_AMU = 1e-6`.
- `simulation/ion.py` — `biphasic` dispatch arm in the `use_drag` loop: `biphasic_step`
  (threading `rng`, passing the resolved `gate_steepness`) → rebuild the BAOAB closure at `m⁺`
  → `baoab_propagation_step` → **fold `e_bind_pair(n)` into `E_pot`**; both writer calls pass
  `mass_scenario=cfg.mass_scenario`.
- `simulation/ion_initial_state.py` — `biphasic` branch: start at `complex_mass_amu(21)`,
  seed the **S2 onset** `E_int[:,0] = f_int·coulomb_available_eV`, and seed the **binding fold**
  `E_pot[:,0] += e_bind_pair(21)` (consistent offset from t0).
- `config.py` — `check_biphasic_config` (require f_int/f_ret non-None under biphasic; **warn,
  not raise**, on λ₀==0) wired into `validate()`.

**Key as-built accounting decision — the E_pot binding fold (R4 closure).** The plan §2.1
closure table lists `E_pot += D_0` (shed) / `−D_0` (pickup) but the Slice-G *interface* omitted
`E_pot`, and MD `E_pot` is recomputed each step with no binding ladder. Reading MASS §6 (the
jump-consistency `E_bind ± D_0`) resolved it: the 5-term invariant closes iff the binding
release is booked to a tracked term. As-built, the **pair-binding potential
`e_bind_pair(n) = −Σ_{i≤n}D_0(i)` is folded into the stored `E_pot`** on the biphasic path (a
pure function of the genuine `n`, so an event shifts it by exactly `±D_0`, matching the table);
the electrostriction marginal stays the untracked "A8 → bath" collective term. Verified by the
isolated forced-shed / forced-pickup tests, which close the **full per-ion 5-term residual to
machine precision (atol 1e-12)** — the strongest evidence the R4 partition is right. Worked
oracle in the build scratchpad.

**Composition facts as-built.** Frozen draw order = evaporation `rng.random(2N)` then pickup
`rng.random(2N)`; both-fire ⇒ shed-then-pickup (pickup rate reads post-shed `n`); K2 reuses
`newton_cool_step` on `E_solv.struct = E_∞(n)+E_int` (asymptote cancels at fixed `n` → the
drain `E_int·(1−e^{−dt/τ})` books to `E_dissip`); cooling precedes the evaporation draw (gate
`t×` opens as `E_int` drains below `Σ(n)`); `E_int` evolves additively (`reconstruct` never
called); pickup density gate uses `rho_he_ratio` at the pre-step depth with the resolved
drag-gate steepness (shared surface).

**Integrator flags → fixed policy (no config fields).** One-event-per-step, jump-then-O, and
momentum-conserving resets are structural in `biphasic_step`/the driver arm and asserted by the
tests; no `one_mass_event_per_step` / `jump_o_step_ordering` / `mass_jump_velocity_reset`
fields were added (rule 2 — no dead one-valued surface).

**Rule-2 table.** The Slice-P/U declared-but-unread carries are now **read by the Slice-G
driver** and are removed from the exception table: `pickup_rate_coefficient`,
`pickup_occupancy_exponent`, `internal_energy_partition_fraction` (f_int),
`internal_energy_retained_fraction` (f_ret), `evap_rate_prefactor_per_ps`, `evap_rrk_dof`,
`gate_onset_override_eV`. `he_capture_velocity` is activated (G passes `u_he=at_rest`). No new
config field lands (integrator policy is structural).

**Tests (2 new files + 2 updated).** `tests/test_biphasic_config.py` (9: presence-require +
λ₀ warn), `tests/test_biphasic_step.py` (19: n_shell seams + scope flip; biphasic init;
cooling-drain exact closure; `m↔n` invariant over 50 event steps; **isolated forced-shed /
forced-pickup full 5-term closure to 1e-12**; seeded reproducibility + event firing; full
`run_ion_propagation` biphasic smoke — v7, finite, onset, `E_int≥0`, `m↔n` across columns,
5-term `ion_ledger_closure` residual < 1e-2 Verlet). Two tests updated for the legitimately
changed behavior: `test_ion_drag_smoke` (biphasic now **admitted**; active noise still refused)
and `test_ion_initial_state` (biphasic now starts at the n=21 complex mass).

**Regression.** `fixed` / `anchored_discrete` never enter the biphasic arm (dispatch on
`mass_scenario`); the writer's `rint` branch keeps `n_shell` byte-identical; the whole Tier-0/1a
drag + collision suite is unchanged. Full suite **1773 passed, 0 failed** (was 1745; +28), same
25 warnings (17 pre-existing `anchored_discrete` + 8 v6→v7 shim). **Phase C (Slice X + Slice G)
is complete.** Next: Phase D (Slice Z — the generative-vs-anchored bridge),
`TIER2_PHASE_D_IMPLEMENTATION_PLAN.md`.

---

## Phase C — Slice G post-delivery code review + fixes (2026-07-02)

An 8-angle post-delivery review of the Slice-G diff (3 correctness + reuse/simplification/
efficiency/altitude/conventions finders, 1-vote verify) surfaced 10 findings; all fixed under a
fresh `[PROCEED TO IMPLEMENTATION]`. The review **confirmed the physics composition**: every
energy booking (S1 split on the pre-pickup rung, K1 drain, capture/shed reduced-mass defects,
the eV conversion, the frozen evaporation-then-pickup draw order, cooling-before-draw) checks
out against the Phase-A/B primitives — the findings were dispatch/edge holes and
maintainability debt, not channel physics.

**Correctness fixes (4 production files + tests):**
1. **Biphasic-without-drag-bundle hole (the top finding).** `mass_scenario='biphasic'` with
   `drag_coefficients=None` passed every guard and silently dispatched onto the hard-sphere
   collision path — with the n=21 mass / S2 onset / E_pot fold seeded at column 0, the fold
   lost from column 1 (spurious +Σ(21) ledger jump), `E_int` frozen, and a frozen verbatim
   `n_shell` stored while attachment grows the mass. Now refused twice: `check_biphasic_config`
   requires the bundle at config-load; `ion._check_scope_ion_driver` re-checks for
   non-validated configs (NotImplementedError before any stepping).
2. **Writer contract hardened; as-built record corrected.** The delivery note above claimed
   "both writer calls pass `mass_scenario`" — **wrong as-built**: the post-loop tail write
   (`ion.py`) omitted it. Fixed, and `write_ion_state_to_checkpoint_column`'s `mass_scenario`
   kwarg is now **required** (an omitted scenario on any future caller is a TypeError, not a
   silent rint fallback), with a fail-loud ValueError on `biphasic` + `n_shell=None`. Analysis
   note: the tail branch itself is *unreachable defensive code*
   (`1 + floor((n−1)/s) ≡ ceil(n/s)`), so no stored column was ever actually mis-written — the
   fix closes the latent contract break, not a live corruption.
3. **Negative λ₀ sign hole.** `check_biphasic_config` warned only on `== 0.0`; a sign typo
   (λ₀ < 0) validated silently and structurally disabled pickup (`P_attach < 0`, draws in
   `[0,1)` never fire). Now a config-load ValueError under biphasic **and** a module-level
   guard in `pickup.lambda_attach` (mirroring the `p < 0` defense). The Slice-P "no load-time
   bound on λ₀" decision covers the value prior, not the sign. λ₀ == 0.0 stays advisory
   (decision #5 unchanged).
4. **`helium_density_profile='tabulated'` lazily refused.** `biphasic_step` hardcoded the erf
   gate and silently ignored the config-accepted `tabulated` member; it now raises
   NotImplementedError at point-of-use — the same rule-2 contract as
   `pickup_rate_form='sweeping'` / `he_capture_velocity='thermal'`.

**Cleanup fixes:**
5. **Stray script flags reverted.** The `tier0_drag_comparison.py` ENERGY_FIGURE/FORCE_FIGURE
   False→True flips were unrelated to Slice G; reverted to the committed state (re-flip
   locally if wanted — they are `# USER SETTINGS`).
6. **eV-conversion single-sourced (rule 1).** New `_amu_ang2_ps2_to_eV` helper in
   `ion_propagation_step.py` replaces the three inline `U*(100²)/EV` idioms
   (`baoab_propagation_step`, `shed_step`, `biphasic_step`). Its op order is the exact
   pre-refactor inline sequence, so the **Tier-0/1a outputs stay byte-identical**; the
   biphasic site previously pre-multiplied (different association) and changes in the last
   ulp — free pre-commit, tests updated in lockstep.
7. **Biphasic t0 physics grouped.** `build_initial_ion_state`: the byte-identical
   anchored/biphasic mass-init branches merged into one arm; the E_pot binding fold moved
   down to join the S2 onset in a single biphasic column-0 block (float-identical: the fold
   now applies to `E_pot_eV[:,0]` after the fill instead of to `E_pot_t0` before it).
8. **Ladder hot-path cache (bit-identical).** `dissociation_ladder.ladder_cumsum` now serves
   Σ from an `lru_cache`d, padded, read-only prefix table (`_sigma_prefix_table`; prefix sums
   are independent of later rungs, so one `(picture, κ)` entry serves the whole run). Kills
   the ~5 redundant per-step ladder rebuilds the biphasic driver paid (e_bind fold, e_inf ×2,
   gate threshold). No numerical-behavior change — regression-tested for exact float equality
   across call orders and for cached-table immutability.
9. **Test-helper dedup (rule 1).** `test_biphasic_step.py` now imports the production
   `_E_kin_eV` and `_amu_ang2_ps2_to_eV` instead of shadow copies, shares one
   `_biphasic_cfg` base with `test_biphasic_config.py` (single source of the "minimal valid
   biphasic cfg"; local wrapper only defaults pickup off), and the per-ion 5-term residual
   became one module-level `_five_term_residual`.
10. Driver-level unknown-`mass_scenario` reject coverage (lost with the old
    `test_scope_guard_rejects_mass_scenario`) restored.

**Review-extension tests (same session, before the fixes; all plan-§3 oracles previously
uncovered):** both-fire ⇒ shed-then-pickup with exact `−D₀(n₀)+f_ret·D₀(n₀)` E_int booking
(the swapped order would book the n₀+1 rung — numerically distinct); a bitwise manual-replay
lock of the frozen evaporation-then-pickup RNG stream; gate opening at t× on a cooling ramp
(no shed above Σ, certain shed on crossing, τ-independent construction); `E_int ≥ 0` + finite
across 50 event steps; `n_shell=None` ValueError; the m↔n guard's teeth (0.5 amu drift trips).

**Files.** Production: `config.py` (guard extended), `simulation/ion.py` (driver guard + tail
kwarg), `simulation/ion_propagation_step.py` (writer contract, tabulated refusal, conversion
helper), `simulation/ion_initial_state.py` (branch merge + t0 grouping),
`physics/pickup.py` (λ₀ sign guard), `physics/dissociation_ladder.py` (cached Σ table).
Docs: `docs/simulation/ion_module.md` pseudocode de-drifted (biphasic arm + writer kwarg).
Tests: `test_biphasic_config.py`, `test_biphasic_step.py`, `test_dissociation_ladder.py`,
`test_ion_propagation_step.py` (writer kwarg). Full suite **1790 passed, 0 failed**
(1773 → +7 review-extension, +10 fix tests), 26 warnings (25 pre-existing + 1 expected §6.5
pairing RuntimeWarning now that the shared biphasic test cfg carries a real constant-mass
bundle under `allow_inconsistent_mass_pairing=True`). Not addressed (recorded, low priority):
the E_pot-fold ownership stays in the driver per the as-built decision above; the per-step
`m↔n` allclose assert stays (safety net beats its ~µs cost).

---

## Phase A — post-delivery external review pass (2026-07-02)

An extensive post-delivery code review of the full Phase-A surface (the L/K/U
modules, the `constants.py` anchors, the config fields + three guards + derived
property, and the six suites), followed by fixes under the
`[PROCEED TO IMPLEMENTATION]` trigger. Independent numerical re-verification of
**all 32 plan §2/§4 golden oracles on a fresh interpreter** (first rungs +
ordering + floor + cliff geometry; Σ(21) bands + ~10% κ-spread decoupling;
|S(n*)| pinning + |S|/n*; E_∞ monotonicity + OQ6 zero; E_elec ≤ 0 incl. κ=100;
split closure; Newton exact factor + dt-robustness + fixed point; cold-shed
pair+E_int neutrality; S2/S1/K1 identities; reconstruction + pre-t× refusal;
f_int floors at both budgets × pictures) passed **32/32 — no physics bug**.
Scope guards verified: no RNG / integrator / mass / velocity / `γ` in Phase A;
import graph is exactly the L→K→U DAG. Four minor findings, behavioral ones
fixed test-first (watched RED → GREEN):

1. **`dissociation_ladder="tabulated"` was accepted but silently ignored.**
   `check_ladder_config` validates the enum, but `biphasic_step` composes the
   Form-U module functions directly and never read the selector — a `tabulated`
   config validated cleanly and silently ran Form U (principle-4 violation; the
   gap opened at Slice G, which gave the twin `helium_density_profile=
   'tabulated'` arm a lazy refusal but not this one). **Fixed:** `biphasic_step`
   now refuses `dissociation_ladder != 'form_u'` with `NotImplementedError` at
   point-of-use (same rule-2 lazy-refusal contract; the `TabulatedLadder`
   fallback has no config data path yet), + a `TestDriverGuards` regression.
2. **`TabulatedLadder` missed the fractional-n hardening its Form-U twin got**
   at the Slice-L review: `d0_of_n(2.5)` / `ladder_cumsum(2.5)` died with
   numpy's cryptic integer-index `IndexError`. **Fixed:** extracted the shared
   `_as_integer_occupancy` helper (rule 1 single source; the Form-U
   `ladder_cumsum` inline block now routes through it — message and behavior
   unchanged, bit-identity cache regression still green) and applied it to both
   tabulated lookups: fractional `n` raises loudly, integer-valued floats are
   accepted and cast (the `ladder_cumsum` contract). +2 tests (reject/accept).
3. **Cache-coherence hazard documented** (note only, no code change):
   `_sigma_prefix_table` (the Slice-G hot-path lru_cache) is keyed by
   `(picture, kappa, n_top)` but builds through the module-global `d0_of_n` /
   `_FIRST_RUNG_EV`, so monkeypatching those *in the ladder's own namespace*
   after a table is cached would serve stale values. Current suites are
   hygienic (stubs patch consumer namespaces only; the `"broken"`-picture test
   never calls `ladder_cumsum`); the hazard + the `cache_clear()` escape hatch
   are now stated in the cache docstring.
4. **`cooling_relaxed` docstring drift fixed:** the module docstring and the
   `constants.py` comment still called the arm "a rule-2 declared-but-unread
   stub"; per the Slice-L ExitPlanMode resolution (b) the *arm* is live
   (provisional arithmetic mean, read by the ordering oracle) and the rule-2
   carry is the concrete blend *value* (pinned at Phase F). Wording aligned in
   both places; no behavior change.

**Deliberate non-fixes (recorded, below the bar):** a NaN `dt_ps` slips past
`newton_cool_step`'s `dt < 0` guard and propagates NaN (a NaN `tau_ps` is
caught); each distinct `n > _SIGMA_TABLE_N_TOP` builds its own prefix-table
cache entry (`maxsize=64` bounds it; no correctness impact); boolean `n` is
accepted as occupancy 0/1 by the integer coercion.

**Files.** Production: `physics/dissociation_ladder.py` (shared
`_as_integer_occupancy` + tabulated guards + cooling_relaxed / cache-hazard
docstrings), `simulation/ion_propagation_step.py` (ladder-form point-of-use
refusal + Raises doc), `physics/constants.py` (cooling_relaxed comment).
Tests: `test_dissociation_ladder.py` (+2), `test_biphasic_step.py` (+1).

**Tests:** new tests watched RED for the diagnosed reasons (cryptic
`IndexError`; `DID NOT RAISE`) then GREEN; Phase-A suites + biphasic/
evaporation consumers 707 green; full suite **1793 passed, 0 failed**
(1790 → +3), 26 warnings (all pre-existing: §6.5 mass-pairing RuntimeWarnings
+ the v6→v7 `E_int` migration UserWarning).

## Phase B — post-delivery external review pass (2026-07-02)

An extensive post-delivery code review of the full Phase-B surface
(`b3187d5..ef1b2d0`: Slices ρ/P/Q + the config/constants/cross-cutting layer) by
four parallel review subagents, followed by fixes under the
`[PROCEED TO IMPLEMENTATION]` trigger. **Verdict: zero Critical findings; every
high-risk item independently verified correct** — the self-bound gate direction
(boundary-consistent across `is_self_bound`/`gate_margin_eV`/`rrk_rate`, the item
the plan itself once had inverted), the `_reduced_mass_defect_coeff` sign refactor
(**bit-exact** `cold_shed` regression over 2000 random cases vs the base commit),
the `spatial_gate` rewire (verbatim body move, scalar return type unchanged), the
S1/K1 pre-event-`n` off-by-one hotspots, the pickup/evaporation RNG contract
(empirically: M scalar calls ≡ one components call incl. post-call generator
state), the §8 out-of-scope leak guard, and the atomic `pickup_*` rename. Full
suite at ef1b2d0: 1726 passed + 9 env-dependent skips.

**Timing note / already-fixed findings.** The review ran against the Phase-B
commit range *after* Phase C (Slices X + G) had landed on the branch, so two
findings were already fixed by Slice G before this pass: the `lambda_0 < 0`
config-load + point-of-use guards ("Slice-G review fix"), and the
`helium_density_profile='tabulated'` point-of-use refusal in the driver.

### Fixes applied (test-first; watched RED `DID NOT RAISE` → GREEN)

1. **ν sign guard (the λ₀ class, both layers).** A negative
   `evap_rate_prefactor_per_ps` gave `k < 0 → P_shed < 0`: the evaporation
   channel *silently never fires* (draws live in `[0, 1)`). Added the
   config-load refuse to `check_biphasic_config` (now "non-negative channel
   rates": λ₀ **and** ν) + the mirroring module defense in
   `evaporation.rrk_rate`. The Slice-Q "ν carries no load-time bound" decision
   covered the *value* prior (Sourced 2.42), not the sign;
   `check_evaporation_config`'s docstring now says where the sign guard lives.
2. **`rho_ratio < 0` point-of-use guard in `pickup.lambda_attach`** — same
   silent-shut-off class; defense-in-depth (helium_density's two arms are
   `[0, 1]` by contract and cannot produce one).
3. **Pickup scalar↔components RNG-consumption parity lock with a *real*
   generator** (`TestRNGConsumptionParity`): M scalar `pickup_step` calls must
   equal one `pickup_step_components` call — outputs, fires (incl. a suppressed
   full-shell ion that still consumes its draw), **and the post-call PCG64
   stream state** (next-uniform bitwise equal). The stub-based
   `TestComponentsMatchScalar` checks outputs but sidesteps consumption; a
   refactor skipping the draw at `P == 0` (the original plan wording!) would
   have broken the Slice-X-frozen contract with all tests green. Evaporation
   already had this lock; pickup did not.
4. **`TabulatedDensityProfile` finiteness guard.** NaN passed the `[0, 1]` check
   vacuously (fails both comparisons) → `ratio()` returned NaN silently; ±inf
   depth nodes degenerate `np.interp` (probed: `ratio(-5) = 0.5` where ≈1 is
   right); a single-node `(nan,)` grid passed everything. `__post_init__` now
   refuses non-finite nodes (+3 reject tests). Realistic failure mode: one NaN
   row in the future sourced TDDFT CSV.
5. **Gate-override n-validation bypass closed.** With `gate_onset_eV` set,
   `_gate_threshold_eV` never called `ladder_cumsum`, so its `n >= 0` /
   integer-occupancy rejection was skipped (`rrk_rate(E, -1, gate_onset_eV=…)`
   silently returned 0.0 while the no-override path raises). Both threshold
   sources now validate `n` identically (covers `rrk_rate`, `is_self_bound`,
   `gate_margin_eV`). Reuses the Phase-A-review `_as_integer_occupancy` single
   source.
6. **`effective_dof` fractional-n rejection.** `effective_dof(2.5)` silently
   truncated on the scalar path (`int(4.5) = 4`) while the array path kept 4.5;
   now both raise via `_as_integer_occupancy` (integer-valued floats still cast).
7. **Direct golden lock on `_reduced_mass_defect_coeff`** (named in the plan's
   test spec but only indirectly covered): bit-exact against the literal
   `0.5*(m*m_He)/m_plus` in the implementation's operation order, both
   directions (shed `m−m_He`, capture `m+m_He`).
8. **`p = 0` Langmuir characterization locked + documented**: `0**0 == 1`, so
   the cap is structurally inert at `n ≥ n*` (≡ `cap="none"` under a langmuir
   label). Deliberately *allowed* (only `p < 0` is rejected); docstring +
   characterization tests make any future tightening a conscious choice.
9. **Golden value pin on `_erf_complement`**: the parity suites are routing
   locks (equality by construction); one new test pins the helper bitwise to
   the literal inline `0.5*(1−erf(depth/14.2))` so a numerically-different
   rewrite (e.g. `0.5*erfc`) cannot slip through green.
10. **Test hygiene:** `test_pickup` disjoint-halves tolerance tightened from
    0.02 (~17σ, nearly vacuous) to a sample-size-justified 5σ band on the
    difference of two Bernoulli means; dead `_ALWAYS` fixture removed from
    `test_evaporation`; stale `_ensemble` comment ("third is
    self-unbound-eligible" — all three are in-band) fixed.
11. **Docstring corrections:** the `k ∈ [0, ν)` overclaim qualified in
    `evaporation.py` (module + `rrk_rate`: strictly `< ν` only for `s > 1`;
    `= ν` at `n = 1` and under an `s = 1` override — asserted by the existing
    s=1 test); `_gates._erf_complement` / `rho_he_ratio` scalar-return
    annotations honest (`float | np.ndarray`; the erf ufunc collapses 0-d to
    `np.float64` — behaviour unchanged); `capture`/`capture_velocity_components`
    now flag that a scalar `u_he` means the He *vector* `(u, u, u)` and that the
    Tier-3 `thermal` arm needs per-ion velocity *components*, not a scalar speed.
12. **CLAUDE.md rule-2 pointer fixed:** it referenced a "rule-2 exception table
    in `drag_migration_log_tier0.md`" that does not exist; the carries are prose
    entries in the **active phase's** log (this file). The pointer now says so
    and states the nomenclature convention the review flagged as three-way
    inconsistent: a field whose only reader is a config-load *guard* is still a
    carry — **"guard-live" ≠ "physics-live"** (applies to `he_capture_velocity`,
    `gate_onset_override_eV`, `evap_rrk_dof`).

### Record corrections (review findings on this log itself)

- **Suite-count convention:** the Slice-Q record's "full suite 1735 passed, 0
  failed" conflates environment-dependent skips — a fresh checkout gives
  **1726 passed + 9 skipped (= 1735 collected)**; the 9 skips need
  experimental-run / VMI reference data. Counts below follow "N passed (+M
  env-dependent skips)".
- **Old run directories now fail loudly at load (intentional, previously
  unrecorded):** the atomic `mass_rate_* → pickup_*` rename means every
  pre-Phase-B `cfg.json` (which serialized the old field names via `asdict`)
  is refused by `run_directory.load_cfg`'s unknown-field check with the
  "different version" message. Correct fail-loud behaviour, not a bug —
  regenerate or hand-migrate archived run dirs.
- **Plan-deviation record:** the P/Q test suites use the *real* L/U modules
  where the plan's Independence sections said "mocked stubs". Stronger as
  composition checks (and picture-threading is asserted by cross-picture
  inequality), but a ladder bug would co-vary with the test oracle — recorded
  as the accepted trade, not silent drift.

### Deliberate non-fixes (recorded, below the bar)

- NaN `rho_ratio` / NaN `E_int` still pass silently (the documented deliberate
  non-guard, consistent with drag); the new guards catch *sign* errors only.
- `evaporation_step_components` runs `cold_shed_velocity_components` on all
  ions, so an ion with `m ≤ m_He` fails the step loudly even if suppressed —
  physically unreachable (m ≥ bare I⁺ ≈ 126.9 amu) and fail-loud is the right
  default; noted that the failure is step-global, not per-fire.

### Follow-ups for Slices X/G (Phase C) — for a future run

Phase C was already implemented when this pass ran, so these were **not**
applied to the X/G surfaces; a future Phase-C touch should:

1. **Treat the ν refuse as part of the Slice-G guard surface.** It lives in
   `check_biphasic_config` (docstring item 3 now covers λ₀ *and* ν;
   `TestNegativeNuRejected` in `test_biphasic_config.py`). Any Slice-G doc that
   enumerates the guard's checks should count both rates.
2. **Draw-order lock dependency:** the new `TestRNGConsumptionParity` (pickup)
   is now load-bearing for the Slice-X frozen shed-then-pickup stream contract.
   Any driver change to how channel draws are consumed must consciously update
   *both* channels' parity tests (evaporation's is
   `TestReviewExtensions::test_rng_parity_with_a_suppressed_ion` + companions).
3. **Override path now raises on bad `n`:** `_gate_threshold_eV` validates
   `n >= 0` + integrality even when the driver threads
   `cfg.gate_onset_override_eV`. The driver's own `n` arrays are non-negative
   ints, so no behaviour change on the production path — but a diagnostic
   harness feeding raw floats through the override now fails loudly.
4. **Config-to-module threading test (Q-review recommendation, still open):**
   add one test that `cfg.gate_onset_override_eV → gate_onset_eV` and
   `cfg.evap_rate_prefactor_per_ps → nu` reach the channels through
   `biphasic_step` (the field-vs-kwarg name split invites silent mis-wiring).
5. **Tier-3 `thermal` arm shape contract:** `capture_velocity_components`
   takes a scalar `u_he` (He vector `(u, u, u)` per ion). The thermal arm needs
   per-ion He velocity components — extend the signature (`ux, uy, uz`) and the
   driver's `_resolve_u_he` threading; do **not** thread a scalar thermal speed.

**Files.** Production: `physics/evaporation.py` (ν guard; override-path n
validation; `effective_dof` integrality; bound docstrings),
`physics/pickup.py` (ρ-ratio guard; p=0 docstring), `physics/helium_density.py`
(finiteness guard; annotation), `physics/_gates.py` (docstring/annotation),
`physics/mass_jump.py` (u_he docstrings), `config.py` (ν refuse in
`check_biphasic_config` + docstrings), `CLAUDE.md` (rule-2 pointer +
guard-live/physics-live convention). Tests: `test_evaporation.py` (+8/−1
fixture), `test_pickup.py` (+6, tolerance fix), `test_helium_density.py` (+4),
`test_biphasic_config.py` (+3), `test_mass_jump.py` (+1).

**Tests:** 14 new guard tests watched RED (`DID NOT RAISE`) → GREEN; the 6 new
lock/characterization tests passed on arrival by design (regression locks on
verified-correct behaviour). Touched suites 228 green; full suite **1813
passed, 0 failed** (1793 → +20 new tests), 27 warnings (26 pre-existing §6.5
mass-pairing / v6→v7 migration + 1 more §6.5 pairing warning from the new
`TestNegativeNuRejected` `validate()` path).

## Phase C — Slice X post-delivery external review pass (2026-07-02)

An extensive post-delivery code review of the Slice-X surface (commit `6988f1e`,
base `49cf08b`) by **three parallel review subagents** — (1) checkpoint schema /
migration cascade, (2) 5-term energy closure, (3) step-state threading + RNG
freeze + doc alignment — each reviewing **both the historical diff and the file
state at HEAD** (`ae68ba5`), since Slice G (`d5a2028`) and the Phase A/B review
pass landed *after* X. Every Important finding was independently re-verified in
code before fixes were applied under the `[PROCEED TO IMPLEMENTATION]` trigger.

**Verdict: zero Critical findings — the Slice-X physics/schema surface is
correct at HEAD.** Independently confirmed: the v5→v6→v7 cascade (v5→v6
sub-transform byte-identical to base; unknown versions still fail loudly;
genuine-v7-missing-`E_int_eV` raises, protected twice), the 5-term summation
(axis/sign/`/n`/units; the miswire tests are jointly discriminating), the
threading contract (every post-X `IonStepState` construction site carries
`E_int_eV`; seams copy-safe; `fixed`/`anchored_discrete` still emit pure-zero
`E_int_eV` because G's t0 seeding is gated on `mass_scenario == "biphasic"`),
and the shed-then-pickup freeze (recorded §2.9, implemented in order, pinned by
real stream-consumption locks; the `ae68ba5` ladder refusal raises before any
draw, so the frozen stream contract is untouched). R4 double-count: not present
— pickup `D_0` enters once, split `f_ret`/`(1−f_ret)`, offset once by the
`e_bind_pair` fold delta.

**Fixes applied (6 Important findings):**

1. **Energy-balance figures now draw `E_int`.** Both figure builders totalled
   the 5-term `E_system` but plotted only four components — on a biphasic run
   the figure would read as a phantom non-closure (`E_int` starts at the
   eV-scale S2 onset). `plot_run_summary._section_ion_energy` plots
   `totals.E_int_eV` behind the same `is not None` guard as `E_mass_transfer`;
   `plot_ion_energy_balance._build_figure` plots it unconditionally (ion-only
   script; module docstring notes MATLAB has no such term). The env-gated ion
   smoke now asserts all six legend labels (a dropped trace fails the smoke).
2. **Plan §3 Slice-X step-state bullet annotated (NB, annotation-not-rewrite):**
   it still said `E_int_eV` "defaults to zeros"; the delivered field is
   required-no-default (design call 2) — the NB points at this log's record.
3. **A/B-review follow-up 4 closed** — `TestConfigThreadingThroughDriver` in
   `test_biphasic_step.py`: (a) `cfg.gate_onset_override_eV → gate_onset_eV`
   lock (baseline forced-shed fires under the default Σ(n) gate — non-vacuity —
   then an override below the in-band `E_int` closes the gate and no ion sheds);
   (b) `cfg.evap_rate_prefactor_per_ps → nu` lock (ν=0 ⇒ the forced-shed state
   never fires). Load-bearing proven by mutation: hard-coding
   `gate_onset_eV=None` in the driver fails (a) with every ion shed (watched
   FAIL → restore → GREEN); previously no test threaded a non-None override
   through `biphasic_step`, so a dropped kwarg would have passed the suite.
4. **`checkpoint.py` version-history register extended to v7** (the comment
   block future bumps consult stopped at v6 while `_ION_SCHEMA_VERSION = 7`).
5. **Checkpoint docstrings aligned with the loader's actual contract:** the
   module docstring claimed "extra fields can be added; the loader uses …
   explicit defaults" / "additions are backward-compatible" — the loader raises
   on any missing field (v6 and v7 were bumped for additions precisely because
   of that); `load_ion_checkpoint` now documents the v5→v6→v7 cascade, the
   synthesized zeros, and the `UserWarning` callers must expect; the
   `_load_checkpoint` shim parenthetical updated; `n_shell` /
   `mass_history_kg` field bullets gained their biphasic semantics (G-era
   drift).
6. **Scalar-load latent defect fixed (pre-existing, on the load path X
   extended):** under `from __future__ import annotations` the loader's
   `f.type is int/float/str` branches were permanently dead and only
   `schema_version`/`mass_scenario` were name-matched, so `num_molecules`
   loaded as a 0-d ndarray, violating the documented "scalar int" contract.
   All scalar fields are now name-matched (`schema_version`, `num_molecules` →
   `int`; `mass_scenario` → `str`), dead branches removed. Watched RED
   (`type(array(5)) is int` fails) → GREEN; regression lock
   `test_scalar_fields_load_as_python_scalars` covers both dataclasses.

**Deliberate non-fixes (Minor, recorded, below the bar):** the `EnergyTotals`
docstring's MATLAB provenance overclaim (MATLAB has no `E_int`; scope the match
to the first three terms); "MASS §6 eq. line 935" line-number citations
(brittle vs NB-annotation drift; cite section names); the `E_int_eV` field
comment omitting the `/num_molecules` normalization; the
`test_E_int_deposit_offset_by_dissip_closes` docstring mis-describing the K2
direction; `IonStepState.n_shell` / read-seam docstring drift after G; the
"≤1 mass event" gloss vs the documented both-fire shed-then-pickup case
("≤1 per channel" is the accurate phrase); the stub-based evaporation parity
lock vs pickup's real-PCG64 post-state check; test tightenings (assert exactly
one migration warning, dtype-preservation assert, non-square axis-lock matrix);
the dead duplicate branch in `_save_checkpoint`; no save-side shape validation
(consistent with the pre-X contract). Follow-up 5 (Tier-3 `thermal` shape
contract) stays open by design.

**A/B follow-up register after this pass:** 1 satisfied, 2 honored
(standing rule), 3 closed (`ae68ba5`), **4 closed (this pass)**, 5 open by
design (Tier 3).

**Files.** Production: `simulation/checkpoint.py` (scalar name-match + v7
register + docstrings), `scripts/post_processing/plot_run_summary.py`,
`scripts/post_processing/plot_ion_energy_balance.py`. Docs:
`TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` (§3 NB), this log. Tests:
`test_checkpoint.py` (+1), `test_biphasic_step.py` (+2),
`test_plot_legacy_debug_smoke.py` (legend assertion).

**Tests:** the scalar-load lock watched RED → GREEN; the gate-onset threading
lock proven load-bearing by driver mutation (FAIL → restore → GREEN); touched
suites green (`test_checkpoint` 23, `test_biphasic_step` 33,
`test_plot_legacy_debug_smoke` 6). Full suite **1816 passed, 0 failed**
(1813 → +3 new tests; env-dependent data present on this machine, so the 9
data-gated tests ran), 27 warnings (all pre-existing §6.5 mass-pairing /
v6→v7 migration).

## Phase C — Slice G post-upstream re-review (2026-07-02)

After the Phase A/B post-delivery external review pass (`ae68ba5`) changed the
upstream physics modules `biphasic_step` composes, and the Slice-X re-review
(`17b21b0`) touched `checkpoint.py` + `test_biphasic_step.py`, Slice G was
re-reviewed at HEAD by **two parallel review subagents**: (1) interface /
behavior consistency of the full `d5a2028..HEAD` diff against the current
Slice-G call sites (incl. RNG draw-order and 5-term-ledger verification +
test runs), (2) the five recorded X/G follow-ups + doc alignment + rule-2
carry audit.

**Verdict: Slice G still sound — no Critical or Important code findings.**
The upstream pass was purely guard-additive + doc-corrective: not one numeric
path a valid biphasic configuration exercises changed. Verified specifically:

- **No new guard is reachable from the production biphasic path.** The ν<0
  refuses (config + `rrk_rate`), the ρ<0 refuse in `lambda_attach`, the
  `_gate_threshold_eV` override-path integer-`n` validation, the
  `effective_dof` fractional-n refuse, and the `helium_density` finiteness
  guard all fire only on invalid config or diagnostic inputs (driver values:
  default ν=2.42; erf-complement ρ∈[0,1]; `np.rint`-seeded ±1-evolved
  integer-valued `n`).
- **Frozen shed-then-pickup RNG contract intact:** each channel still consumes
  exactly one unconditional `rng.random(size=M)` per step; every new guard
  raises (killing the step) rather than skipping a draw. The `ae68ba5`
  `dissociation_ladder="tabulated"` refusal in `biphasic_step` raises *before*
  any draw — and closed a genuine Slice-G hole (that enum previously validated
  at load then silently ran Form-U energetics).
- **5-term ledger sources byte-untouched** (K2 drain, K1 drain, S1 split,
  capture/shed reduced-mass defects, eV conversion); the m↔n consistency
  assert intact; `TestIsolatedEventClosure` still closes at 1e-12.
- `ion.py` / `ion_initial_state.py` byte-identical since `d5a2028`;
  `mass_jump.py` / `_gates.py` / `constants.py` changes are docstring-only.
- **A/B follow-up register re-audited:** 2/3/4/5 exactly as the X-review
  register states; **1 was optimistic** — the guard, its docstring, and its
  tests are complete (four checks: f_int/f_ret, bundle require, λ₀ *and* ν
  sign refuse, λ₀==0 advisory), but the Phase-C plan's item-(4) enumeration
  still described the two-check pre-build guard. Fixed this pass (see below);
  follow-up 1 now **fully satisfied**.
- **Rule-2 carry ledger coherent:** `noise_*` (×4), `validation_histogram_metric`
  (Phase-E activation scheduled), `drag_low_v_floor`, the unbuilt enum arms
  (`sweeping`/`dwell_time`, `thermal`, tabulated density/ladder — all
  point-of-use refusals), and the two non-field carries all recorded; neither
  review commit added, removed, or silently activated a config field; all
  Slice-G retirements verified live in code.

**Doc fixes applied (this pass, docs-only — no code or test changes):**

1. `docs/simulation/ion_module.md` `_check_drag_scope` table de-drifted (was
   stale since Tier-1a and contradicted the pseudocode in the same file):
   `mass_scenario ∉ {fixed, anchored_discrete, biphasic}`,
   `drag_form ∉ REALIZED_FORMS`, mass trip-wire marked `fixed`-only (§6.6
   skip), + a note on the `ion._check_scope_ion_driver` biphasic-bundle
   re-check (G-review fix).
2. `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` — as-built NB after the Slice-G
   decision box counting the guard's **four** checks (closes follow-up 1).
3. `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` — as-built NBs (§Slice-G knobs +
   §4 config contract): the three integrator flags landed as structural
   policy hard-coded in `biphasic_step`, **not** config fields (rule 2); only
   `he_capture_velocity` is a real field.

**Deliberate non-fixes (recorded, below the bar):** the ν=0 / λ₀=0 advisory
asymmetry in `check_biphasic_config` (λ₀=0 warns, ν=0 silent — defensible:
ν defaults to 2.42 so 0 must be deliberate, λ₀ defaults to 0.0); the private
name `_gate_threshold_eV` leaking into the `_as_integer_occupancy` ValueError
context string reached from public entry points; log entry :1974-1975 citing
already-retired carries as guard-live examples (append-only history, accurate
for the Phase-B era it describes).

**A/B follow-up register after this pass:** **1 satisfied (doc gap closed
this pass)**, 2 honored (standing rule), 3 closed (`ae68ba5`), 4 closed
(`17b21b0`), 5 open by design (Tier 3).

**Tests (re-run by the review, no code changed):** targeted slice-G +
upstream suites **361 passed, 0 failed** (2 expected §6.5 pairing warnings);
full suite **1816 passed, 0 failed**, 27 warnings (all pre-existing §6.5
mass-pairing / v6→v7 migration).

### Re-review nits implemented (2026-07-02, user `[PROCEED TO IMPLEMENTATION]`)

The two "deliberate non-fixes" above were promoted to fixes on user request
(test-first; both watched RED → GREEN):

1. **ν=0 inert advisory (symmetry fix).** `check_biphasic_config` now warns
   (does not raise) on `evap_rate_prefactor_per_ps == 0.0`, mirroring the
   λ₀==0 pickup-inert advisory: the RRK channel is structurally inert (k=0,
   P_shed=0), the run is pickup-only growth from n₀ — legitimate as a
   diagnostic, but ν defaults to the Sourced 2.42 so an explicit 0 can no
   longer pass silently. Docstring item 4 renamed to "Rate-zero inert
   warnings" covering both rates. No default-config behavior change (λ₀
   defaults 0.0 → still exactly one advisory; ν default 2.42 → silent). The
   driver-level `test_nu_prefactor_threads_to_evaporation` is unaffected (it
   constructs `SimConfig` directly, no `validate()`).
2. **Private-name leak in the override-path error.** `_gate_threshold_eV` now
   calls `_as_integer_occupancy(n, context="gate threshold")` (was the private
   helper name, which leaked into the `ValueError` a caller of the public
   `rrk_rate`/`is_self_bound`/`gate_margin_eV` sees). Wording locked by
   tightening the two fractional-n matches in `test_evaporation.py` to
   `"gate threshold requires integer occupancy"` (consistent with the
   neighboring `"gate threshold requires n >= 0"`).

**Files.** Production: `config.py` (ν=0 advisory + docstring item 4),
`physics/evaporation.py` (context string). Tests: `test_biphasic_config.py`
(+2, `TestEvapInertWarns`; module docstring), `test_evaporation.py` (2
matches tightened, RED→GREEN with fix 2).

**Tests:** touched suites (`test_biphasic_config`, `test_evaporation`,
`test_biphasic_step`, `test_dissociation_ladder`) **211 passed**; full suite
**1818 passed, 0 failed** (1816 → +2 new tests), 27 warnings (all
pre-existing §6.5 mass-pairing / v6→v7 migration — the new ν=0 advisory
fires in no existing test's `validate()` path).

---

## Phase D — Slice Z pre-build decisions (2026-07-02)

A pre-build refinement session on `TIER2_PHASE_D_IMPLEMENTATION_PLAN.md` (plan-only; no
code — the `[PROCEED TO IMPLEMENTATION]` trigger has not been given for Slice Z). The
plan's interface assumptions were verified against the delivered surfaces (Phase-C driver,
v7 `E_int_eV`/`n_shell` fields, `ladder_cumsum` / `lambda_attach` / `compare_*` /
`ion_ledger_closure` signatures, `tier1a_common` scaffolding) — all sound. Both CE
fragments are ions (`ion_initial_state` seeds all 2N rows at the n=21 complex mass), so
"mean over ions" = mean over all 2N rows with no row-selection convention. Seven decisions
locked by the user; two analytical findings folded into the plan:

1. **Π oracle corrected (plan §2.2/§3 — the headline finding).** The plan's
   "Π > 1 during dense traversal" acceptance oracle is wrong at the 9 Å / 0.80 eV priored
   central point: `Π(n) = λ₀·ρ_ratio·(1−n/n*)₊·f_ret·τ` is **0 at gate-open** (n = 21 = n*
   zeroes the Langmuir cap) and **≤ ~0.2 for the whole 21→14 decline** at central values
   (λ₀ 0.9/ps, f_ret 0.1, τ 6.5 ps; even n = 14 gives factor 1/3). Π > 1 needs the upper
   corner of every band simultaneously (f_ret ≳ 0.4, τ → 16.5, λ₀ → 1.1) — the *strip* end
   of the MASS §6.11 regime axis, exactly what the gentle 9 Å condition is not. Corrected
   expectation: **freeze side** — Π = 0 at gate-open, Π < 1 throughout (the 21→14 decline
   is the initial-`E_int` cascade, not pickup-sustained shedding), Π → 0 at exit; a
   computed Π > 1 is itself a flag. Condition-specific: not to be carried to 2.70 eV.
2. **`t×` has a closed form — the check upgraded from editorial to sharp.** Pre-crossing
   dynamics are fully deterministic (evaporation gate shut; pickup dead at n = n*), so
   `t× = τ·ln(f_int·E_avail/Σ(21))` exactly, identical for every ion. The helper returns
   the **per-ion array** (all-agree = wiring check); the bridge checks the reconstruction
   against the analytic value; and `f_int` ↔ GAH25-prior are recognized as coupled choices.
3. **Staircase target = the deterministic schedule family**, `build_shell_schedule(t★)` for
   t★ ∈ {0.5, 5.0, 9.0} as an overlay band, **t★ = 5.0 primary** — numerically identical to
   the delivered anchored artifacts' `n_shell` but with zero dependence on
   gitignored/machine-local run dirs. "Loaded, not regenerated" honored in spirit.
4. **Representative priored point pinned:** λ₀ = 0.9/ps, τ = 6.5 ps (geometric mid of
   [2.6, 16.5], coinciding with GAH25 shell-1 t₀), **f_int = 0.5** (lands
   t× ≈ 4.8–5.6 ps over the mixture Σ(21) band 0.17–0.19 eV; f_int just above the
   0.21–0.24 floor would give ~1.9 ps, factor ~3 under the prior), f_ret = 0.1,
   picture = `statistical_mixture`, **κ = 1.0** (config default / Slice-L oracle center);
   ν = 2.42/ps, s = 3n−3 Sourced/Derived. Demonstration choice, not calibration (Phase F).
5. **Π helper contract:** `regime_parameter(ckpt, cfg)` **re-derives ρ_ratio(t)** from
   checkpoint positions + `droplet_radii_angstrom` + the cfg density-gate steepness (the
   plan previously took `rho_ratio` as an unspecified input), then `lambda_attach·f_ret·τ`
   per ion; reported mean ± envelope.
6. **Run parameters inherit Tier-1a exactly** (N50 → 100 ions, 30 ps, Tier-0 dt/seed/drag
   bundle, `allow_inconsistent_mass_pairing=True`); only `mass_scenario="biphasic"` + knob
   values differ. Deliverable shape settled: helpers in new `postprocess/bridge_diagnostics.py`
   (Phase-E D2 generalizes it), orchestration mirroring `gen_tier1a_runs`, report script under
   `scripts/post_processing/`, editorial verdict in `TIER2_PHASE_D_BRIDGE_FINDINGS.md` +
   this log.
7. **Phase-C follow-up #4 retired (record correction).** The "still open" config→module
   threading test in the Phase-A/B review-pass follow-ups list is **already delivered** —
   `test_gate_onset_override_threads_to_evaporation` and
   `test_nu_prefactor_threads_to_evaporation` landed in `tests/test_biphasic_step.py`
   during the Slice-G re-review pass. Phase D adds no duplicate.

Plan doc updated in place (`TIER2_PHASE_D_IMPLEMENTATION_PLAN.md` §0/§2/§3/§4/§5/§7/§8).
Slice Z remains behind the `[PROCEED TO IMPLEMENTATION]` trigger.

### Phase D — MASS / CALIBRATION_MAP cross-check of the refined plan (2026-07-02)

A final internal-consistency pass of the refined Phase-D plan against
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` and `CALIBRATION_MAP.md`
(precedent: the Phase-C cross-check above). **Verified consistent:**

- **`t×` chain:** the `t× ≡ GAH25 t₀` identity, the [1, 15] ps sanity band, the
  ±factor-2 prior / factor-10 flag are verbatim MASS §6.11. Every assumption
  behind the closed form holds in the delivered code: gate suppressed while
  `E_int > Σ(n)` (§6.11 Lyapunov — the gate *always* opens, so the helper's
  NaN arm is an edge guard only); pickup dead at `n(0) = 21 = n*`
  (`ANCHOR_N_START = 21 = N_STAR`; `pickup_occupancy_cap` defaults `langmuir`,
  `p = 1.0`, driver threads `cap=cfg.pickup_occupancy_cap`); K2 per-step drain
  `E_int·(1−e^{−dt/τ})` compounds to an exact exponential; cooling precedes the
  draw; noise inert. Σ(21) mixture 0.17–0.19 eV (row 21; ~11% κ-independent
  over [0.3, 5]) → t× ≈ 4.8–5.6 ps at the pinned point.
- **Π regime reading:** the boxed §6.11 `Π(n) = λ(n)·f_ret·τ` and the
  freeze/shed criterion match; the corrected freeze-side oracle is exactly
  MASS's gentle-9 Å regime ("self-binds early → retains a shell"), with Π = 0
  at gate-open being the §6.11 *filling-driven* stabilizing route and Π → 0 at
  exit the *exit-driven* one; the superseded "Π > 1 persists" reading is the
  [Calvo24] strip end.
- **Pinned point vs the map:** λ₀ = 0.9 ∈ central 0.7–1.1 (row 7, "not 2.0
  Na⁺"); τ = 6.5 ∈ [2.6, 16.5] (row 11), geometric mid, coincides with GAH25
  shell-1 t₀ 6.53; f_ret = 0.1 = "prior small" (row 13); ν = 2.42 pinned
  (row 9); s = 3n−3 Derived (row 10); E_avail 0.80 = ½·14.40/9, scenario-keyed
  and guard-stamped (row 15 / MASS A7). **κ = 1.0 is clean:** row 19's "7.5×
  radial cliff" prior is the *physical* shell-1→shell-2 binding drop (R3,
  encoded in Form U's `D₀(1)/D_floor`), not a κ value; κ only shapes the cliff
  sharpness, and 1.0 is the config default / log-center of the oracle sweep.
- **Schema/config/scope:** v7 `E_int_eV` (MASS §7 NB), the `biphasic` config
  literal (§11 NB), no new fields, and the §6.11 mean-field ODE staying
  excluded from D all line up.

**Two findings:**

1. **Plan §8 f_int note completed (fixed in place).** The note cited only the
   "~0.2 soft upper"; MASS A7 (2026-06-21 rework) states the operative
   scenario-invariant **velocity-consistency ceiling ~0.6**
   (`E_trans/E_avail ≈ 40%` at both budgets), so the derived 0.80 eV window is
   [0.21–0.24, ~0.6] and the pinned f_int = 0.5 sits *inside* it (the ~0.2
   advisory is vacuous at 0.80 eV — the floor exceeds it). Strengthens the
   pinned choice; note reworded.
2. **Cross-doc wording drift (reported, not fixed):** CALIBRATION_MAP row 14
   condenses the f_int upper edge as "soft upper ~0.2 is advisory" and omits
   the ~0.6 ceiling MASS A7 states alongside it; MASS A7's own paragraph
   carries both numbers with drafting ambiguity, flagged "provisional pending
   OQ2". Not a Phase-D blocker (D fits nothing); a one-line row-14 touch is
   suggested whenever that row is next edited.

---

## Phase D — Slice Z DELIVERED (2026-07-03) — Phase D COMPLETE

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger. TDD for the
helpers (test→RED→GREEN); the slice composes only delivered modules (Phase-C
`biphasic` driver, Tier-0/1a scaffolding, `compare_*`, 5-term
`ion_ledger_closure`) — **no new physics, no new config fields** (plan §4
honored: the two helpers read config, add nothing). Two small open points were
resolved before the build (user silent on the AskUser round; recommendations
applied as announced): **τ = the config default 6.55 ps** (the delivered
Slice-K geometric mid; the plan's "6.5" is rounded shorthand — closed-form t×
moves 5.19→5.23 ps, immaterial) and **the biphasic cfg builder seeded in
`scripts/tier2_common.py`** (Phase-F F1 extends it rather than migrating a
bridge-local copy).

**Build (4 files + tests):**
- `i2_helium_md/postprocess/bridge_diagnostics.py` — `mean_shell_count`
  (2N-row average), `crossing_time_ps` (per-ion `min{t: E_int < Σ(n(t))}`,
  strict inequality, threshold follows the stored `n(t)`, NaN edge arm),
  `regime_parameter` (Π = λ·f_ret·τ with the §2.2 reconstruction contract:
  ρ_ratio re-derived from checkpoint positions + `droplet_radii_angstrom`
  through the *shared* `simulation.ion._drag_gate_steepness` resolver — no
  formula copy; fail-loud on `tabulated` density profile / unset f_ret).
  Phase-E D2 generalizes this module.
- `scripts/tier2_common.py` — `build_biphasic_cfg` wrapping
  `tier0_common.build_drag_cfg` (the `tier1a_common.build_anchored_cfg`
  mirror): scenario swap, scenario-stamped 0.80 eV, pinned knobs as defaulted
  kwargs (λ₀ 0.9, f_int 0.5, f_ret 0.1; Phase F sweeps by argument), §6.6
  override, n=21 metadata mass; `tier2_bridge_run_dir_name` /
  `TIER2_BRIDGE_TAG`.
- `scripts/gen_tier2_bridge_run.py` — single-run orchestration mirroring
  `gen_tier1a_runs.py` (N=50, 30 ps, dt 0.01, seed 20260604 — Tier-1a
  inherited exactly).
- `scripts/post_processing/tier2_bridge_report.py` — loads the run dir, emits
  the three overlay figures into `<run>/figures/` (mean-n vs
  `build_shell_schedule` family t★∈{0.5,5.0,9.0} with 5.0 primary;
  R/|v| vs the 9 Å CSV on the `compare_*` overlap grid; Π mean±envelope with
  the t× marker) + prints/writes the numeric summary
  (`bridge_summary.txt`).
- `tests/test_bridge_diagnostics.py` (19): hand-average; exact first-crossing
  per ion + NaN arm + strict-inequality + n-driven-threshold cases (one test
  premise was inverted during RED — the gate opens when E_int falls *below*
  Σ(n), so occupancy *growth* raises Σ past a constant E; fixed test-side,
  code unchanged); the §2.2 **analytic Newton-cooling oracle** (reconstruction
  within one dt above `τ·ln(E₀/Σ(21))`); Π formula/ρ-parity on hand-built
  positions via a duck-typed stub + full-shell zero + outside-droplet →0 +
  freeze-side magnitude + the two fail-loud arms; the **few-step driver
  smoke** (fast-τ override pulls t× in-window): 5-term closure < 1e-2,
  helpers consume the real v7 checkpoint, per-ion t× finite and all-agree.

**Bridge run + report (out-of-pytest, plan §6).** Run dir
`9A_drag_shared_pure_cubic_N50_tier2_bridge_biphasic` (machine-local).
Headline numbers (full detail + verdict: **`TIER2_PHASE_D_BRIDGE_FINDINGS.md`**):

- **Sharp oracles all pass:** t× reconstructed 4.960 ps vs closed form
  4.951 ps (one stored dt; all 100 ions agree exactly); Π = 0 at gate-open,
  max mean Π = 0.0054 (deep freeze side), Π → 0 at exit; 5-term residual
  0.0011 % of E_system(0). R/|v| vs TDDFT land at the fixed-null baseline
  (R ratio 1.300 vs 1.299 fixed / 1.297 anchored-t5.0) — no bridge
  regression; the mass scenario is kinematically invisible at this scale
  (confirms the Tier-2 premise that the size distribution is the only
  discriminating observable).
- **The staircase is missed:** emergent mean n(t) = 21 → 20.33 (frozen by
  ~8 ps; per-ion terminal [18, 21]) vs the anchored 21→19→14. Quantified in
  the findings doc: the RRK integral under K2 cooling at the pinned point is
  ≈ 0.7 expected sheds (k ≈ 0.36/ps at gate-open, s−1 = 59 suppression,
  τ = 6.55 drain) — the observed Δn̄ = 0.67 is the mechanism's genuine
  prediction, not a miswire. Cascade budget is kinetically, not
  energetically, limited. **Levers attached:** κ + picture (the free co-fit
  pair, via D₀(n)/Σ(n)), τ within [2.6, 16.5]; λ₀ is *not* the lever (pickup
  dead, Π ≈ 0); ν/s are Sourced/Derived — if Phase F cannot land the
  staircase inside the free bands, the RRK dof convention becomes an OQ-class
  finding, not a silent retune.
- **Stale-artifact note validating decision #3:** the Tier-1a anchored run
  dirs on this machine carry pre-Phase-B cfg fields (`mass_rate_*`) and no
  longer load via `RunDirectory`; the schedule-family overlay needed no run
  artifact (the anchored comparator row was scored by direct `ion.npz` load
  through the v6→v7 shim).

**Out-of-scope guard honored:** no sweep/calibration, no abundance/VMI/
Wasserstein read, no Phase-E machinery, no mean-field ODE, no noise, no
drag-law/RNG/anchored-artifact touch.

**Rule-2 table.** Slice Z adds no config field and retires no carry (the
helpers read live fields only). The standing carries are unchanged from the
Slice-G re-review audit.

**Tests:** `test_bridge_diagnostics.py` 19 green (RED watched at module-absent
collection + one genuine expectation fix during GREEN); full suite
**1837 passed, 0 failed** (1818 → +19), 27 pre-existing warnings (§6.5
mass-pairing / v6→v7 migration). **Phase D (Slice Z) is complete.** Next:
Phase E (Slices E1–E5, `PHASE_E_IMPLEMENTATION_PLAN.md`), with the bridge
flag standing as the first Phase-F calibration target.

### Phase F plan de-drifted after the Slice-Z build (2026-07-03, docs-only)

`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` was amended (as-built NB annotations,
not a redesign) where Slice Z changed its ground truth:

- **§1 prerequisites, row D:** Slice Z marked delivered, with the staircase
  flag spelled out as the concrete meaning of "the 0.80 eV validation lands"
  in the F5 gate.
- **F1 is now an *extension*, not a creation:** `scripts/tier2_common.py`
  exists (bridge-scoped subset — pinned-point kwargs `lambda0_per_ps` /
  `f_int` / `f_ret`; 0.80 eV stamped in-builder; τ/κ/picture on config
  defaults; fixed bridge tag). F1 adds the campaign kwargs (picture, κ,
  tau_ps, budget, relaxation enable) + the knob-encoding tag helpers. Two
  signature corrections recorded: **λ₀ joins the F1 signature** (the planned
  one omitted the Bounded row-7 knob; the delivered kwarg is the precedent),
  and a **back-compat constraint** — the Phase-D scripts import this module,
  so extension is kwargs-with-defaults only (the bridge run must stay
  reproducible). F1's unit tests remain open (Phase D landed none for
  `tier2_common`).
- **F2 reuse row:** `gen_tier2_bridge_run.py` added as the delivered
  single-run template; plus a four-point **bridge-priors NB** on the stage
  design: (1) Stage 1 has a quantitative target (the ~10× under-shed; κ grid
  must span the sharp-cliff end; an in-band failure escalates to the RRK-dof
  OQ, not a ν/s retune); (2) Stage 2's τ-*insensitivity* premise may invert
  at 0.80 eV (bridge estimate: τ→16.5 ps ≈ ×3 shed integral — expect the flag
  arm); (3) the Stage-1 `f_int = floor` pin is a *timing* choice (t× ≈ 0 at
  the floor vs ≈ 5 ps at the bridge point — the cascade budget is f_int-
  independent, pinned at Σ(21)); (4) the 9 Å / 0.80 eV leg carries ~no λ₀
  sensitivity (Π ≤ 0.005), making the 18 Å contrast leg load-bearing. Plus an
  operational note: cfg.json version-skew invalidates campaign run dirs
  (`RunDirectory.load_cfg` refuses unknown fields — the Tier-1a artifacts
  already demonstrate it); score promptly, regenerate rather than migrate.

No code changed; F3–F6 contracts untouched. The `[PROCEED TO IMPLEMENTATION]`
boundary for Phase E/F holds.

## Phase D — Slice Z post-delivery code review + fixes (2026-07-03)

Four parallel review agents over the Phase-D commit (`c0ca2c6..606ebf3`):
helpers physics, scripts, tests (with a live suite run), and a cross-cutting
scope/docs audit. **No Critical findings from any reviewer.** Independently
verified clean: the full-diff sweep (exactly the seven planned files, nothing
outside), the §4 no-new-fields contract, the §7 leak-guard grep (only prose
hits), the pinned-point values (τ = 6.55 confirmed as the plan's own geometric
mid, a rounding non-deviation), Tier-1a run-parameter inheritance
line-for-line, and every checkable findings-doc number recomputed from the
delivered ladder (Σ(21) = 0.18784 eV, closed-form t× = 4.951 ps, Π ceiling
0.196, RRK k ≈ 0.36/ps, test tallies 19 / 1837). Three Important findings +
hygiene; fixes applied same-day (user request, no `[PROCEED]` needed — review
follow-up on delivered code, precedent: the Slice-G post-delivery pass):

- **`crossing_time_ps` gate unified onto the driver's own surface.** The
  helper compared `E < ladder_cumsum(n)` directly — a semi-copy of the gate
  that would *silently misreconstruct* a `gate_onset_override_eV` diagnostic
  run (the driver's `_gate_threshold_eV` swaps Σ(n) for the fixed override;
  such runs exist via the `allow_gate_onset_override=True` warn-through arm).
  Now consumes the public `evaporation.is_self_bound` with a
  `gate_onset_eV=None` passthrough kwarg (plan-§3 signature extended
  backward-compatibly; rule-1 upgrade — the "matching the driver's gate"
  docstring claim is now literal). `tier2_bridge_report.py` threads
  `cfg.gate_onset_override_eV` through and prints a NOTE that the Σ(21)
  closed form does not apply when an override forced the gate.
- **`regime_parameter` τ guard + docstring-contract repair.** The `Raises`
  docstring promised rejection of "unset/non-positive `f_ret` or `τ`" but the
  code checked only `f_ret is None` (a hand-built cfg bypassing `validate()`
  hit a bare `TypeError` on `tau=None`, and τ ≤ 0 returned a sign-flipped Π).
  Fix went *both* directions: τ now fail-louds on `None`/non-positive
  (mirroring the `check_solvation_cooling_config` load guard), while the
  f_ret wording was **narrowed to unset-only** — the promised non-positive
  rejection was wrong, since f_ret = 0 is in-band ([0, 1] per
  `check_internal_energy_budget_config`) and Π ≡ 0 is its *correct*
  reconstruction, so enforcing the old docstring would have refused a
  legitimate run.
- **Π steepness-resolution blind spot closed (test-only).** The ρ-parity test
  ran only on the default gate, where `_drag_gate_steepness(cfg)` collapses
  to `potential_steepness` — a wrong hard-coding would have passed. New
  `erf_independent` variant with a divergent pair (drag_gate_steepness 5.0 vs
  potential_steepness 14.2, near-surface positions, plus a self-check that
  the two ρ fields actually differ) makes the §2.2 "same resolved steepness
  the run used" contract falsifiable.
- **Report fail-loud + docstring.** `tier2_bridge_report.main()` now rejects
  a non-`biphasic` run dir with a clear `ValueError` (was: a bare `TypeError`
  formatting `f_int=None` mid-report) and carries a docstring.
- **Rule-2 hygiene:** dead `from dataclasses import replace` import removed
  from `tests/test_bridge_diagnostics.py`.

**Deliberate non-fixes (recorded, deferred to their owners):** promotion of
the private `_drag_gate_steepness` import to a public helper + the
`postprocess/__init__` export (Phase-E D2, the module's named generalization
point, together with shape/monotonic-`time_ps` guards in the thin helpers);
the third `_run_one` orchestration copy + the override→run-dir-tag coupling
(Phase-F F1/F2 touch `tier2_common` anyway); the `E_system(0)` relative-
residual division (physically far from zero on the ion stage); the loose
`0 < Π < 1` freeze-side assert (the plan's Π oracle is deliberately
qualitative); the findings-doc ad-hoc comparator-baseline provenance (already
disclosed in §3 there; promote to a scripted path only if Phase E/F needs the
row again); the "Pi at gate-open (t=0)" summary label (Π = 0 holds on all of
[0, t×], so the number is right).

**Tests:** `test_bridge_diagnostics.py` **22 green** (19 → +3: override
threading, divergent-steepness ρ-parity, τ rejection); full suite
**1840 passed, 0 failed** (1837 → +3), same 27 pre-existing warnings. No
config field, no schema, no RNG, no drag-law touch; the §7 guard holds.

---

## Phase E — pre-build refinement decisions (2026-07-03)

A pre-build refinement pass over `TIER2_PHASE_E_IMPLEMENTATION_PLAN.md` (plan-only; no
code — the `[PROCEED TO IMPLEMENTATION]` trigger has not been given for E1–E5),
cross-referenced against MASS §2.1–2.2/§4/§6.11/§R5/§7/§11 + CALIBRATION_MAP (rows
7/9/10/11/13/14/19/20/21, anchor coverage), `TIER2_PHASE_D_BRIDGE_FINDINGS.md`, the
Phase-F plan's Phase-E interface bindings, and a full code-surface audit of every
composed surface (`biphasic_step` + its config reads + RNG consumption,
`ion_state_from_checkpoint_column`, v7 `IonCheckpoint`, `ion_ledger_closure`,
`make_ion_accel_fn`/`make_ion_baoab_step`, `bridge_diagnostics.py`, `hedft_loader`,
`lambda_attach`/`rho_he_ratio`, `complex_mass_amu`, `RunDirectory`, the real abundance
CSV, and the test-fixture precedents). Plan rewritten in place. Three decisions locked
(R1/R2 user; R3 recommendation applied with the user idle on the ask — the Phase-D
"applied as announced" precedent):

1. **R1 — E2 composes the delivered `biphasic_step` verbatim (λ₀=0 relaxation view),
   not a shed-only reduced loop.** The plan's "draws only shed Bernoullis" wording
   contradicted the delivered unconditional-draw channel contract: at
   `pickup_rate_coefficient=0` pickup is structurally inert but its draw is still
   consumed, so the Slice-X frozen two-draw stream is preserved byte-identically
   (`check_biphasic_config` already names λ₀=0 the relaxation stage's "close cousin").
   Drag drops out via a zero-γ BAOAB closure (decay=1, `dE_dissip≡0`), keeping the
   delivered bookkeeping. The stage runs on its own generator
   (`SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))`) — draw order extended by a new
   stage, never re-ordered. Audit fact recorded: with λ₀=0 and γ=0 the mass subsystem
   decouples from translation entirely (positions enter only the dead density gate),
   so the `free_flight` arm is exact for the observable; `coulomb` stays the default
   for ledger/asymptotic fidelity, and the two arms must give the identical shed
   sequence under one seed (new oracle).
2. **R2 — E5 = `git mv postprocess/bridge_diagnostics.py → derived_diagnostics.py`.**
   The delivered module's own docstring names Phase-E D2 as its generalization point;
   the parent + Phase-F plans bind to `derived_diagnostics.py`; a parallel new module
   would duplicate t×/Π (rule 1). Helpers move unchanged; the two Phase-D importers
   (`tier2_bridge_report.py`, `test_bridge_diagnostics.py` →
   `test_derived_diagnostics.py`) update mechanically. E5 absorbs its recorded
   Slice-Z-review deferrals: the `_drag_gate_steepness` promotion to a public helper,
   the `postprocess/__init__` exports, and the shape/monotonic-`time_ps` guards.
3. **R3 — abundance-CSV provenance gap → README stub + open item.**
   `integrated_i_he_abundance.csv` has no documented producer (no exporter under
   `data/reference/scripts/`, no README entry; the only repo reference is the consumer
   `plotting_histogram.py`). E3 adds a `data/reference/README.md` data-contract entry
   (columns/units/normalization verified 2026-07-03: 21 rows, n=0–20 contiguous,
   `ionPercent` sums to 100.0000, mass grid step 4.0026 u/e) with provenance marked
   "to be completed". **OPEN ITEM:** user to supply measurement IDs / the producing
   script.

**Corrections + refinements folded into the plan (audit-settled, no fork):**

- **E2 "E_dissip constant" oracle was wrong** — the delivered Slice-G booking sends
  the K2 cooling drain to `E_dissip` every step. Corrected sharp oracle: with γ=0,
  λ₀=0, `ΔE_dissip(t)` ≡ the cumulative K2 drain exactly (drag + pickup-bath
  contributions identically zero). The relaxation driver also reproduces the ion
  arm's `e_bind_pair(n)` E_pot fold (rule 1 — same call), so the 5-term closure holds
  in both force arms.
- **`relaxation_time_ps` = required-when-enabled, no default** — MASS §R5 gives no
  sourced t_exp (only "hundreds of ps"; no µs flight time anywhere in MASS); the value
  is a Phase-F campaign choice; freeze early-exit (`n=0` or `E_int<D_0(n)` ∀ ions —
  k≡0 forever, E_int monotone non-increasing) bounds the cost. `relaxation_dt_ps`
  defaults to `dt_ion` with a load guard ν·dt ≤ 0.1 (one-event-per-step bias
  ≲ (k·dt)²/2 ≈ 0.5%; MASS authorizes no larger step; default ν·dt ≈ 0.024).
  `check_relaxation_config` (the `check_biphasic_config` pattern) additionally
  requires `mass_scenario=="biphasic"` when enabled.
- **E2 artifact = a bona fide v7 `IonCheckpoint` for the relaxation window**, saved
  as `relaxation.npz` via the existing `save_ion_checkpoint` (auto-stride convention
  reused; no schema bump, zero new I/O code, `ion_ledger_closure` applies unchanged),
  wrapped in a thin `RelaxationResult` (checkpoint + `freeze_flags` +
  `time_relaxed_ps`). Seeding via `ion_state_from_checkpoint_column(ion, -1)`.
- **Phase-F interface names pinned into the E plan:** `load_he_abundance_reference`
  (E3 — F3's reuse row already binds the name), matched-time + sim-end W₁ per run
  (E4), `t_cross_ps`/`Pi_t`/`regime_label`/`total_strip_reachable`/`sanity_flags`
  (E5), `relaxation_stage_enabled` (E2).
- **E4:** `chi2`/`ks` Literal arms get point-of-use `NotImplementedError` refusals
  (unbuilt-enum-arm convention); "pure numpy" restated as a style choice (scipy is
  already a hard dependency via `erf`/`curve_fit`); the config dispatch
  (`compare_size_distributions`) is the field's first physics-live reader —
  `validation_histogram_metric` leaves the rule-2 table at the E4 build.
- **E1:** `n_shell` is int-valued float — validate + cast, fail-loud on fractional;
  simulated support is 0–21 (n₀=n*=21, Langmuir-capped; bridge terminal n ∈ [18,21])
  vs the reference's 0–20 — handled by E4's zero-filled union support, not clipped
  in E1.
- **Bridge-findings oracle corrections:** at the 9 Å/0.80 eV pinned point the cascade
  is a ~0.7-shed burst frozen by ~8 ps with max mean Π ≈ 0.005 — E2's shed-to-freeze
  dynamic tests use synthetic hot inputs (the pinned point adds ≈0 relaxation sheds:
  matched-time ≈ sim-end there, recorded as a characterization); E5's regime oracle
  expects `shell_retaining` at 9 Å with Π>1 itself a wiring flag, condition-specific,
  not carried to 2.70 eV. E5's regime label got a concrete documented reporting rule
  (majority-NaN t× or median terminal n ≤ 1 → `total_strip`), a convention, not a
  config knob; `total_strip_reachable` = `e_infinity_eV(0)` under the run's
  picture/κ, annotated as a split-consistency confirmation (K/U tautology
  precedent), not an independent anchor.
- **Prerequisite phasing obsoleted:** Phases A–D are all delivered at HEAD (1840
  green), so every E slice is buildable now; build order E5a (rename first) →
  {E3, E4, E1 sim-end} parallel → E2 → E1 relaxed-input admission + E5b.
- **Relaxation-cost risk downgraded to LOW** (translation decoupled; free-flight arm
  near-free; freeze early-exit; auto-stride storage budget).

**Drift fixed (rename fallout of `061c92b`, live docs only — log history stays
append-only):** `CLAUDE.md` doc list, `CALIBRATION_MAP.md` Phase-E NB,
`TIER2_IMPLEMENTATION_PLAN.md` §4 (both detail-plan pointers), and
`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (entry-docs + cross-links) now cite
`TIER2_PHASE_E/F_IMPLEMENTATION_PLAN.md`.

**Cross-reference verdict:** no contradictions with MASS / CALIBRATION_MAP. The
Wasserstein-on-integer-support choice is MASS §10/§6.9 verbatim; the [1,15] ps t×
sanity band, ±factor-2 GAH25 prior (factor-10 = the genuine flag), Π freeze/shed
criterion, and the Calvo24/OQ6 reachability statement are §6.11 verbatim; the
velocity-observable cut is already recorded in CALIBRATION_MAP's Phase-E NB; the R5
matched-time + upper-bound double-report is §R5's option (a)+(b) combined. Slices
E1–E5 remain behind the `[PROCEED TO IMPLEMENTATION]` trigger.

---

## Phase E — Slice E1 DELIVERED (2026-07-03)

`[PROCEED TO IMPLEMENTATION]` given for E1 after a pre-build discussion that
locked five interface defaults (all user-approved):

1. **Which atoms count → all `2N`, no masking.** The I₂ Coulomb explosion yields
   two I⁺, each carrying its own He shell and each an independent I⁺Heₙ
   detection; the plan oracle (`n=[0,0,1,2,2,2]`, 6 entries) is consistent.
2. **`source` field = mode tag only**, `str ∈ {"sim_end", "relaxed"}` (the R5
   sim-end upper bound vs the E2 matched-time terminal `n`); it is the
   discriminator E4/F3 pair per run. No path stored.
3. **Out-of-range `n` fails loud** (not silently clipped/widened): `n > n_max`
   violates the Langmuir cap (corruption), `n < 0` unphysical; `rint` residual
   checked **exactly** (`== 0`), integer counts stored as float have no
   arithmetic drift — a fractional entry is upstream corruption.
4. **`RelaxationResult` admission = attribute-dispatch, not `isinstance`.** A
   private `_terminal_n_and_source` checks `n_shell` (→ `[:, -1]`, `"sim_end"`)
   then `terminal_n` (→ `"relaxed"`), so E1 ships sim-end now without importing
   the not-yet-built E2 type; the relaxed mode is a one-branch admission.
5. **Degenerate ensemble = empty** (`2N == 0` or `T == 0`) raises loudly.

**New surface.** `postprocess/size_distribution.py`:
- `ShellDistribution` (frozen: `n_values (Nn,) int`, `counts (Nn,) int`,
  `fraction (Nn,) float`, `source str`) + a `mass_amu` property
  (`complex_mass_amu(n_values)`, on-demand rung labelling — kept off the frozen
  tuple so the F3-pinned field set stays exactly the four).
- `compute_terminal_shell_distribution(source, *, n_max=N_STAR=21)` — histogram
  on integer support `0..n_max` via `np.bincount(minlength=n_max+1)`;
  `fraction = counts/counts.sum()`. `n = 21` is a legal rung (support includes
  it; the reference's `0..20` mismatch is E4's zero-filled-union job, not
  clipped here). Reuses `physics/shell_schedule.complex_mass_amu`.
- Exported from `postprocess/__init__` now (E1's two names; E5b adds E3/E4/E5).
  This lands part of the E5.2 deferral early — incremental, no conflict.

**Tests.** `tests/test_size_distribution.py` — **21 green**. Hand-count oracle
(+ terminal-column-not-sentinel guard); integer support `0..21` incl. the legal
`n=21`; custom `n_max`; monotone-falling envelope; both input modes give the
identical histogram (differing only in `source`); a real v7 `IonCheckpoint` via
`test_checkpoint._make_ion_checkpoint`; fail-loud arms (fractional, out-of-range
high/negative, non-finite, empty sim-end, zero-steps, empty relaxed, 1-D
`n_shell`, wrong source type); `mass_amu` matches `complex_mass_amu`; frozen.
Full suite **1861 passed** (1840 → +21), same 27 pre-existing warnings. No
config field, no schema, no RNG, no drag-law touch; the §5 out-of-scope guard
holds. E2–E5 stay behind the trigger.

---

## Phase E — Slice E5a DELIVERED (2026-07-03)

The pure-rename half of E5 (decision R2), kept deliberately minimal so the E5b
`reconstruct_diagnostics` additions land in the final module name with a clean
diff.

**Rename (git-tracked as `R`, history preserved).**
- `git mv i2_helium_md/postprocess/bridge_diagnostics.py →
  postprocess/derived_diagnostics.py` — module body **byte-identical** (helpers
  `mean_shell_count`, `crossing_time_ps`, `regime_parameter` move unchanged; the
  module docstring already names Phase-E D2/E5 as its generalization point, so no
  edit needed).
- `git mv tests/test_bridge_diagnostics.py → tests/test_derived_diagnostics.py`.

**Mechanical importer/reference updates (the only content edits).**
- `tests/test_derived_diagnostics.py`: docstring `bridge_diagnostics.py →
  derived_diagnostics.py`; import `postprocess.bridge_diagnostics →
  .derived_diagnostics`.
- `scripts/post_processing/tier2_bridge_report.py`: same import swap.
- `TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (F1 note): dangling
  `test_bridge_diagnostics.py → test_derived_diagnostics.py` (live-plan drift
  fix, Quality Principle 8). Phase-D historical docs (D plan, D findings) keep
  their `bridge_diagnostics` mentions — append-only history.

**Scope boundary — E5a is rename-only.** The three absorbed deferrals (item 1
`_drag_gate_steepness`→public promotion in `simulation/ion.py`; item 2 the
E3/E4/E5 `postprocess/__init__` exports; item 3 shape/monotonic-`time_ps` guards
on the thin helpers) are **deferred to E5b**, where the new `reconstruct_diagnostics`
code consumes them — item 1 mutates a `simulation/` module (not a rename), and
item 2's E3/E4/E5 surfaces do not exist yet. (E1's own exports already landed
with E1.) Recorded so the E5 acceptance "items 1–3 closed" is met by E5b.

**Verification.** `grep bridge_diagnostics **/*.py` → **no matches**;
`from i2_helium_md.postprocess.derived_diagnostics import ...` resolves;
`test_derived_diagnostics.py` **22 green** (unchanged from the delivered 22);
full suite **1861 passed** (no test lost to the rename), same 27 warnings. No
config/schema/RNG/drag-law touch; the §5 guard holds. E5b + E2–E4 stay behind
the trigger.

---

## Phase E — Slices E3 + E4 DELIVERED (2026-07-03)

`[PROCEED TO IMPLEMENTATION]` given for E3 and E4 together (the parallelizable
pure-function pair — E3 loader on the real CSV, E4 W₁ on synthetic
distributions). Both ship against the real reference + synthetic inputs; no
production runs (Phase F wires those).

### E3 — experimental abundance reference loader

**New surface.** `postprocess/abundance_loader.py :: load_he_abundance_reference(path)
-> HeAbundanceReference` (frozen: `n (Nn,) int`, `label (Nn,) str`,
`mass_center_u (Nn,)`, `ion_counts (Nn,)`, `ion_fraction (Nn,)`, `source_path`).
Mirrors the `hedft_loader` variant contract: order-independent membership header
match (missing **and** extra columns raise), loud `FileNotFoundError`/`ValueError`,
`source_path = p.resolve()`. `ion_fraction = ionPercent / ionPercent.sum()`
(sums to exactly 1.0).

- **CSV read:** `np.genfromtxt(dtype=None, encoding="utf-8", names=True)` — the
  string `label` column coexists with the numeric columns in one structured
  array (the plain `dtype=float` hedft path can't carry the string).
- **Data contract verified against the file (2026-07-03):** 21 rows, `n`
  contiguous 0–20, `ionPercent` sums to **exactly 100.0**, `massCenter` step
  4.0026 u/e, labels `I^+`…`I^+He_20`, spot-checks n=0 → 43.52 %, n=20 →
  0.2645 %. Guard tolerance `|Σ ionPercent − 100| ≤ 1e-3` (the file is exact; a
  modest band absorbs future rounding without masking a broken export).
- **Reference vs sim mass convention (recorded, no conflict):** the reference's
  `massCenter` uses the integer I mass (127 u at n=0); the sim side's
  `complex_mass_amu` uses 126.90. E3 loads the reference column verbatim and
  never recomputes it, so the conventions don't collide (E1/E4 score on integer
  `n`, not mass).
- **R3 provenance (applied).** `data/reference/README.md` gains a data-contract
  entry (columns/units/normalization, consumer `plotting_histogram.py`) with
  **provenance marked "to be completed"**. OPEN ITEM stands: no in-repo exporter
  or measurement-ID record; user to supply.

**Tests.** `tests/test_abundance_loader.py` — **13 green**. Real CSV
(n=arange(21), labels, mass centres, `ion_fraction` sum = 1 + the two plan
spot-checks, `source_path` resolved, int dtype); synthetic happy path +
**column-order-independence**; fail-loud arms (file-not-found, missing column,
extra column, non-contiguous n, n-not-from-0, negative counts, bad percent sum);
frozen.

### E4 — integer-support Wasserstein comparison + config dispatch

**New surface.** `postprocess/distribution_compare.py`:
- `wasserstein_integer_support(sim: ShellDistribution, ref: HeAbundanceReference)
  -> float` — `W₁ = Σ_n |F_sim(n) − F_ref(n)|` over the **contiguous union
  support** (`arange(min, max+1)` guarantees unit spacing, so the CDF-gap sum is
  exact; fractions zero-filled where a side is absent — this is where the sim's
  legal n=21 meets the ref's 0–20). Raises on empty or **disjoint** support (no
  shared `n` ⇒ a units/labelling bug). Pure numpy.
- `compare_size_distributions(sim, ref, *, metric: str) -> float` — the
  `validation_histogram_metric` dispatch; the caller passes
  `cfg.validation_histogram_metric`. `chi2`/`ks` get point-of-use
  `NotImplementedError`; any other value → `ValueError`.

**Rule-2 carry RETIRED.** `SimConfig.validation_histogram_metric`
(`Literal["wasserstein","chi2","ks"]`, default `"wasserstein"`,
`config.py:58/323`) was declared-but-unread (only referenced by its own default
in `test_drag_config.py`). E4 is its **first physics-live consumer**: the value
now flows into a live three-arm dispatch that computes (wasserstein) or refuses
(chi2/ks/unknown). The carry leaves the active rule-2 set at this build — no
config/schema change was needed (the field already existed).

**Tests.** `tests/test_distribution_compare.py` — **16 green**. Hand W₁ oracles
(identical → 0; shift-by-one → 1; shift-by-two → 2; the F=[.5,1,1]/[0,.5,1]
micro-case → 1; symmetry); the **n=21-vs-zero-filled-ref** union case; support
guards (disjoint raises, empty either side raises); real E3 ref × synthetic E1
sim → finite W₁; dispatch (`wasserstein` matches the direct call, `chi2`/`ks`
refuse, unknown raises, and `cfg.validation_histogram_metric` from
`single_pulse_N2000()` flows through); **scipy cross-check** (import-guarded,
tests only — ran, not skipped).

### Shared

- **Exports (progressive E5.2 deferral, per-slice like E1):**
  `postprocess/__init__` now exports `HeAbundanceReference`,
  `load_he_abundance_reference`, `wasserstein_integer_support`,
  `compare_size_distributions` (+ the E1 pair from before). E5b adds the E5
  surfaces.
- **Verification.** E3+E4 suites **29 passed** (13+16, no skips → scipy ran);
  public imports resolve; full suite **1890 passed** (1861 → +29), same 27
  pre-existing warnings, 0 failures.
- **Scope.** No schema, no RNG draw-order, no drag-law, no neutral/ion
  propagation touch. The only config interaction is *activating* an existing
  declared field (no new field, no `validate()` change). §5 out-of-scope guard
  holds. Remaining Phase-E work: **E2** (relaxation stage — the one stochastic/
  driver slice), **E5b** (`reconstruct_diagnostics` + the three absorbed
  deferrals), and E1's relaxed-input admission (gated on E2's `RelaxationResult`).

---

## Phase E — Slice E2 DELIVERED (2026-07-03)

`[PROCEED TO IMPLEMENTATION]` given for E2, the post-ejection relaxation stage
(the R5 mitigation). The driver slice — highest-risk in Phase E — so grounded
against the full delivered seam first (`biphasic_step`, `baoab_propagation_step`,
`make_ion_baoab_step`/`_o_step`, `make_ion_accel_fn`/`AccelFn`,
`ion_state_from_checkpoint_column`/`write_ion_state_to_checkpoint_column`,
`rrk_rate`, `d0_of_n`, `e_bind_pair_eV`, the ion driver's biphasic branch, and the
config check pattern) before writing a line.

**New surface.** `simulation/relaxation_stage.py`:
- `RELAXATION_STREAM_KEY = 0xE2_2026` — the fixed key for the stage's own PCG64
  stream (`SeedSequence((cfg.seed, KEY))`); the ion-stage stream is never touched.
- `RelaxationResult` (frozen: `checkpoint` (v7 `IonCheckpoint`), `freeze_flags
  (2N,) bool`, `time_relaxed_ps`, `terminal_n (2N,)`).
- `run_relaxation_stage(ion, cfg, *, rng=None, save_path=None) -> RelaxationResult`.

**Composition (decision R1 — verbatim reuse).** Per step: `biphasic_step` under a
**relaxation view** `replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)`
(constructed internally, **not** re-`validate()`-d — the λ₀=0 advisory in
`check_biphasic_config` is the sanctioned evaporation-only limit), then a
**zero-γ** BAOAB closure (`_zero_gamma` → decay=1, `dE_dissip≡0`), then the
`e_bind_pair(n)` E_pot fold. Two arms:
- `coulomb` (default): real `make_ion_accel_fn` conservative field; **add** the
  fold to baoab's recomputed E_pot — byte-for-byte the ion driver's biphasic
  branch minus drag.
- `free_flight`: zero-accel closure (ballistic positions, E_kin from post-jump
  (m,v) via the tested baoab path), then **replace** E_pot with `held_MD +
  e_bind_pair(n)` where `held_MD = seed.E_pot − e_bind_pair(seed_n)` (no
  conservative work under free flight → MD potential constant; only the discrete
  fold moves, by `+D_0(n)` per shed). The decoupling fact makes this exact for
  the observable, not an approximation.

**Freeze condition (verified, not guessed).** Read `rrk_rate`: `k=0` when
`E_int ≤ D_0(n)` (via `base=max(0,1−d0/E)`) for n≥2, and `k_direct=0` when
`E_int ≤ D_0(1)` for n=1 (the last-rung barrier). So `freeze = (n==0) | (E_int <
D_0(n))` uniformly captures "k≡0 forever" given E_int monotone non-increasing
(K2 cools toward 0, pickup off) — the plan's exact condition, and it correctly
does **not** freeze a self-bound `E_int > Σ(n)` ion (that un-freezes as it cools).

**Artifact.** A bona fide v7 `IonCheckpoint` for the relaxation window (col 0 =
seed = ion's final column), assembled with the delivered
`write_ion_state_to_checkpoint_column` (no new I/O); `droplet_radii`/`b_ion_outside`
pass-through, `number_of_collisions` zero, `mass_scenario="biphasic"`. Variable
early-freeze length handled by collecting stored states in a list then sizing the
checkpoint to `len(stored)` (no trailing-zero corruption). `ion_ledger_closure`
and `load_ion_checkpoint` apply unchanged.

**Config (new; opt-in; default off → default scope unchanged).**
`relaxation_stage_enabled` (bool, False), `relaxation_time_ps`
(`Optional[float]`, required-when-enabled — no sourced t_exp), `relaxation_dt_ps`
(`Optional[float]` → `dt_ion`), `relaxation_forces` (`Literal["coulomb",
"free_flight"]`, "coulomb"). `check_relaxation_config` wired into `validate()`
after `check_biphasic_config`; no-op when disabled, else requires biphasic +
`relaxation_time_ps>0` + `nu·dt_relax ≤ 0.1` + the forces enum.

**Decisions recorded (defensible calls made inline, per plan latitude):**
1. **`save_path` param, not `RunDirectory` integration.** The stage takes an
   optional `save_path` and writes via `save_ion_checkpoint`; Phase F passes
   `<run_dir>/relaxation.npz`. More decoupled/testable than wiring RunDirectory
   into E2; the "relaxation.npz beside ion.npz" convention is a Phase-F concern.
2. **Private cross-module imports** `_decide_stride_ion` / `_drag_gate_steepness`
   from `simulation.ion` (+ public `DEFAULT_MAX_CHECKPOINT_BYTES_ION`) — reuse the
   auto-stride + gate-steepness single sources (rule 1). Precedent:
   `derived_diagnostics.py` already imports `_drag_gate_steepness` privately.
   E5b's deferral #1 (promote to public) will update both importers together.
3. **E1 relaxed-input admission is already satisfied — no E1 change.** E1's
   `_terminal_n_and_source` duck-types `n_shell` then `terminal_n`;
   `RelaxationResult` exposes `terminal_n` (and no top-level `n_shell`), so
   `compute_terminal_shell_distribution(result)` returns `source="relaxed"`. The
   forward-compatible E1 dispatch (built at the E1 slice) closes this plan item.

**Tests.** `tests/test_relaxation_stage.py` — **15 green**. 5-term invariant
closes on **both** arms with sheds actually firing (asserted `terminal_n < seed
n` so the K1/e_bind/cold-shed bookings are covered, not only K2); coulomb vs
free_flight give the **identical shed sequence** (`n_shell`, `E_int` byte-equal;
the decoupling fact); hot input sheds monotonically (pickup off → `diff(n)≤0`), no
avalanche, `terminal_n ≤ n*`; freeze termination reached (`freeze_flags.all()`,
`time_relaxed < relaxation_time`); cold self-bound input sheds nothing while the
frozen two-draw stream is consumed (real-PCG64 post-state parity — the pickup
`TestRNGConsumptionParity` precedent); the artifact round-trips through
`load_ion_checkpoint` (v7, ledger applies); E1 admits the `RelaxationResult`
(`source="relaxed"`); every `check_relaxation_config` arm fires; the stream key is
guarded against silent change. Fixtures: the delivered `_biphasic_driver_cfg` /
`_tiny_neutral` for the real seed; a hand-built consistent-mass v7 seed
(`m == m_I+ + n·m_He`, paired atoms z-separated so the per-pair Coulomb is finite)
for the dynamics/RNG tests. The native `s=3n−3=60` at n=21 suppresses the RRK rate
to ~0 sheds/window, so the dynamics tests use the delivered `evap_rrk_dof=2`
override (guarded s≥1) to get a non-trivial cascade — a test knob, not a physics
claim.

**Verification.** E2 suite **15 passed**; full suite **1905 passed** (1890 → +15),
same 27 pre-existing warnings, 0 failures. No schema bump (reuses v7), no RNG
draw-order change (the ion stream is untouched; the stage extends it with a new
seeded stream), no drag-law / neutral / ion-propagation touch, no change to
default scope (relaxation defaults off). §5 out-of-scope guard holds.

**Remaining Phase-E work:** **E5b** — `reconstruct_diagnostics` + the three
absorbed deferrals (`_drag_gate_steepness`→public promotion, the E3/E4/E5
`postprocess/__init__` exports batch, the thin-helper shape/monotonic guards).
E1's relaxed-input admission is closed (decision 3 above). After E5b, Phase E is
complete and the program moves to Phase F (the calibration campaign).

---

## Phase E — Slice E5b DELIVERED (2026-07-03) — Phase E COMPLETE

`[PROCEED TO IMPLEMENTATION]` given for E5b, the additive half of E5: the composed
`reconstruct_diagnostics` plus the three absorbed deferrals. Closes Phase E.

**Deferral 1 — `_drag_gate_steepness` → public `drag_gate_steepness`** (single
source, imported without the underscore contract). Renamed the def in
`simulation/ion.py`; updated **all** call sites in one mechanical pass:
`simulation/ion.py` (internal call), `simulation/relaxation_stage.py` (E2 import +
call), `postprocess/derived_diagnostics.py` (import + call + docstring),
`scripts/post_processing/tier0_drag_comparison.py` (import + call),
`tests/test_baoab_propagation_step.py` (import + 4 calls), plus doc/comment
references in `config.py`, `physics/helium_density.py`, `docs/simulation/ion_module.md`,
and `test_derived_diagnostics.py`. No underscore alias kept (clean promotion);
`grep _drag_gate_steepness **/*.py` → **no matches**. The Tier-2 plan docs +
this log keep their historical `_drag_gate_steepness` mentions (append-only).

**Deferral 2 — `postprocess/__init__` Phase-E exports batch.** Added the E5
surface (`Diagnostics`, `TCrossSummary`, `reconstruct_diagnostics`). Combined with
the E1 (E1 slice) and E3/E4 (E3/E4 slice) exports already landed, the full
E1/E3/E4/E5 public surface is now exported — deferral closed.

**Deferral 3 — shape / monotonic-`time_ps` guards on the thin helpers.**
`crossing_time_ps` now rejects a non-strictly-increasing `time_ps`;
`regime_parameter` now rejects positions/`n_shell` that don't share one (2N, T)
shape and a `droplet_radii_angstrom` that isn't (2N,). (`mean_shell_count`
already carried its 2-D guard.)

**New surface (E5b main).** `postprocess/derived_diagnostics.py`:
- `TCrossSummary` (frozen: `median_ps`, `spread_ps`, `all_agree`).
- `Diagnostics` (frozen: `t_cross_ps (2N,)`, `t_cross_summary`, `Pi_t (2N,T)`,
  `regime_label`, `total_strip_reachable: bool`, `sanity_flags: tuple[str,...]`)
  — **field names pinned** for the Phase-F F3/F4 scoreboard.
- `reconstruct_diagnostics(ckpt, cfg) -> Diagnostics` — composes the three
  delivered helpers (re-derives nothing; zero schema cost):
  - **t× + ensemble summary**: `crossing_time_ps` (threading
    `cfg.gate_onset_override_eV`) + median/ptp-spread/`all_agree`
    (`all_agree` = every ion crossed and all within one stored dt — the bridge
    wiring check).
  - **Π(t)**: `regime_parameter` verbatim.
  - **regime label** (a documented reporting rule, **not** a config knob):
    `total_strip` iff the majority of ions have no finite t× **or** median terminal
    n ≤ 1; else `shell_retaining`.
  - **total-strip reachability**: `e_infinity_eV(0)` under the run's picture/κ →
    `isclose(·, 0)`. Structurally true for the Slice-K form (OQ6) — a consistency
    confirmation, not an independent anchor.
  - **sanity flags** (advisory, not raises): t× outside the [1, 15] ps band
    (§6.11; factor-10 miss = the genuine flag, not the GAH25 ±factor-2 prior); Π > 1
    where a freeze-side point is expected (a units/wiring smell).

**Tests.** `tests/test_derived_diagnostics.py` — **35 green** (22 delivered move
unchanged with the E5a rename + regression-lock; **+13** new): the two new guards
fire (non-monotonic time; regime-parameter shape/droplet-radii); `reconstruct_diagnostics`
composes all fields; `t_cross_summary` all-agree vs a constructed disagreement;
regime label flips across shell-retaining / total-strip-by-terminal-n /
total-strip-by-no-crossing constructions; reachability true under all three
pictures; both sanity flags fire on constructed absurd inputs (t× < 1 ps; Π > 1);
the E5 surface is importable from `i2_helium_md.postprocess`.
`test_baoab_propagation_step.py` **7 green** under the renamed helper.

**Verification.** `grep _drag_gate_steepness **/*.py` → none; E5 public imports
resolve; `test_derived_diagnostics.py` **35 passed**; full suite **1918 passed**
(1905 → +13), same 27 pre-existing warnings, 0 failures. No schema, no RNG
draw-order, no drag-law, no propagation touch (the rename is behavior-preserving);
§5 out-of-scope guard holds.

**Phase E is COMPLETE.** All five slices delivered and reviewed: E1
(`size_distribution`), E2 (`relaxation_stage`), E3 (`abundance_loader`), E4
(`distribution_compare`), E5 (`derived_diagnostics` = E5a rename + E5b compose).
The full comparison layer — terminal I⁺Heₙ size-distribution extraction (both
sim-end and relaxed modes), the R5 relaxation stage, the experimental abundance
loader, the integer-support Wasserstein arbiter, and the regime-determination
diagnostics — is built and tested on synthetic checkpoints + the real reference
CSV. The Phase-F pinned interfaces (`load_he_abundance_reference`,
matched+upper-bound W₁ per run, `t_cross_ps`/`Pi_t`/`regime_label`/
`total_strip_reachable`/`sanity_flags`, `relaxation_stage_enabled`) are all
satisfied. Next: **Phase F** — the calibration campaign (F1–F6) that wires
production runs through E1–E5, composing only accepted modules.

---

## Phase E — code review E1/E3/E4/E5 + fixes applied (2026-07-03)

Multi-agent adversarially-verified review of slices E1, E3, E4, E5 (E2 excluded
— dedicated review pending). 14 candidates, 14 verified, 0 refuted; 10 reported.
No physics/units/convention error found; the findings were one provenance bug
and a cluster of validate-early (Quality Principle 4) gaps. All fixes applied
TDD (every guard test watched to fail first):

1. **E1 provenance (the one correctness bug).** A relaxation checkpoint
   reloaded from `relaxation.npz` is a v7 `IonCheckpoint`, so the duck-typed
   dispatch inferred `source="sim_end"` for matched-time relaxed data — F3
   would mis-pair `W1_matched`/`W1_simend_upper`. Fix:
   `compute_terminal_shell_distribution(..., source_tag=...)` explicit
   override (`{"sim_end","relaxed"}`, invalid fails loud); inference default
   unchanged; the reload pitfall documented in the module + function docstrings.
2. **E3 NaN/sign-blind validation.** Blank cells (genfromtxt → NaN) passed the
   sum guard (NaN compares False); negative `ionPercent` summing to 100 loaded.
   Fix: finiteness guards on all numeric columns + `ionPercent ≥ 0`.
3. **E4 unvalidated inputs.** Fractions were never checked, so `.counts` passed
   by mistake (or an E3 NaN) scored a plausible-but-wrong W₁; duplicate support
   silently dropped mass in the scatter. Fix: per-side validation (finite,
   ≥ 0, sums to 1 within 1e-6 — both producers normalize exactly, the band
   only catches factor-level misuse; unique support).
4. **E5 fabricated regime labels.** Empty (2N=0) ensembles and fractional
   `n_shell` (E1's "upstream corruption" condition) produced confident labels.
   Fix: `reconstruct_diagnostics` rejects both, early, mirroring E1's
   conventions on the same v7 field.
5. **E5 Π > 1 flag condition-gated.** The flag fired unconditionally; per the
   plan the smell is freeze-side-specific and "must not be carried to 2.70 eV".
   Fix: `freeze_side_expected: bool = True` param — default keeps pinned-point
   behavior; Phase F passes `False` at shedding-persists conditions.
6. **Cleanup.** `N_STAR` now imported from `physics/constants.py` (shadow copy
   removed, rule 1); the byte-identical header-membership block factored into
   `postprocess/csv_contract.py :: validate_columns` (used by `hedft_loader` +
   `abundance_loader`, messages unchanged); dead `total == 0` branch removed
   (rule 2); `ionCounts` contract corrected in `data/reference/README.md` +
   loader docstring (fractional normalized intensity, **not** raw detector
   counts — the absolute normalization joins the R3 provenance open item);
   `TIER2_PHASE_D_BRIDGE_FINDINGS.md` §4 pointers annotated with the E5a
   rename (append-only convention preserved — annotation, not rewrite).

**Tests.** +19 focused (5 source-tag, 1 N_STAR single-source, 4 abundance
guards, 5 compare guards, 2 reconstruct fail-loud, 2 Π-flag gating); full suite
**1937 passed** (1918 → +19), same 27 pre-existing warnings. No schema, RNG,
drag-law, or propagation touch; no public-surface rename (the Phase-F pinned
names are unchanged; `source_tag` and `freeze_side_expected` are additive
keyword-only params).

---

## Phase E — code review E2 + fixes applied (2026-07-03)

Multi-agent adversarially-verified review of slice E2 (the dedicated review the
E1/E3/E4/E5 pass deferred). Four parallel reviewer dimensions
(physics/energy-ledger, RNG/determinism, architecture/checkpoint-contract, test
adequacy) + an independent probe-driven verification pass on the four
highest-impact candidates. **No physics bug found.** The strongest refutation
attempts all failed: the coulomb arm is verifiably the ion driver's biphasic
branch minus drag; the zero-γ closure is exactly inert (decay = 1.0,
dE_dissip = 0.0, no draw); the `e_bind_pair` fold shifts E_pot by exactly
+D_0(n) per shed; `dt_relax` reaches every dt reader; the two-draw contract,
keyed stream, and coulomb/free-flight decoupling hold. A suspected
freeze-condition bug was **refuted with algebra**: `biphasic_step` passes
`e_inf + E_int` into `newton_cool_step` and subtracts `e_inf` after, so the
asymptote cancels and K2 cools `E_int` toward 0 — E_∞(N) > 0 can never
un-freeze an ion after the early exit. Measured 5-term closure on the delivered
fixtures: ~2e-5 relative.

**Confirmed findings and fixes (all TDD — the four new-behavior tests watched
to fail first; the coverage tests passed immediately against the delivered
code, confirming those behaviors were already correct):**

1. **`time_relaxed_ps` contract (Important, latent).** Returned *absolute*
   time (seed ≈ 20 ps in production + window) while documented window-relative;
   probe showed the natural Phase-F early-freeze check
   `time_relaxed_ps < relaxation_time_ps` answers wrongly on any real seed
   (invisible to the delivered tests only because the fixture seeded
   `time_ps = 0`). Fix: window-relative return
   (`state.time_ps − seed.time_ps`); the checkpoint's `time_ps` axis stays
   absolute (documented). New nonzero-seed-time test.
2. **Seed-column staleness guard (Important, latent).** Exhaustive probe
   (2388 step/stride combos) proved the ion driver's last-column back-fill
   branch (`ion.py:293`) is dead code — the allocation is always exactly
   consumed — so under stride > 1 with `S ≢ 1 (mod stride)` column −1 is a
   stale pre-stride snapshot, and `E_int_eV`/`n_shell` have no `*_final_*`
   fields to recover from. Fix: `run_relaxation_stage` now refuses a seed with
   `mass_history_kg[:, -1] != mass_final_kg` (mass_final is written from the
   true final state; a partial guard — mass events are the observable-relevant
   proxy). Latent today: production ≈ stride 1.
3. **Stale v7 byte-count constant (adjacent live defect found by the
   verifier).** `ion.py::_NUM_2N_T_ARRAYS_ION` was 14 — it omitted the v7
   `E_int_eV` (2N,T) array — undercounting the checkpoint estimate by ~6.7%
   (the stride-1 grant flips at ≈ 2090 molecules on the true size). Fixed to
   15; new schema-count test in `test_ion.py` derives the count from the
   dataclass fields so the constant tracks future schema changes. Focused-bug
   justification for the ion-driver touch: bookkeeping constant only, no
   propagation physics.
4. **Input-provenance guard (Moderate; scenario partially refuted).**
   `ion.mass_scenario` was never checked. Verifier showed realistic
   off-scenario seeds crash *loudly* in `biphasic_step`'s m↔n lockstep
   assertion — but with a "channel/reset bug" message that misdiagnoses the
   cause — while a mass-on-lockstep-grid seed passes silently, and the
   free-flight `held_md_pot` subtraction is mis-based for any seed whose E_pot
   never contained the fold. Fix: explicit `mass_scenario == "biphasic"`
   rejection at entry.

**Test gaps closed (+15 tests, `test_relaxation_stage.py` 15 → 29 plus 1 in
`test_ion.py`):** the plan's per-channel sharp acceptance oracle
**ΔE_dissip ≡ cumulative K2 drain** (reconstructed exactly from the stored
`E_int`/`n_shell` columns; catches channel mis-attribution the 5-term sum
conserves — the one *Critical*-rated coverage gap, found independently by two
reviewers); free-flight `E_pot ≡ held_MD + e_bind_pair(n)` asserted numerically
(a constant-offset bug cancels in the closure); no-avalanche
(`diff(n) ≥ −1`); self-bound gate suppression + not-frozen-at-entry
(un-freeze-into-band as K2 cools); n = 1 direct dissociation to the bare ion
with n = 0 freeze; `relaxation_dt_ps` actually honored (also covers the
time-exhaustion, non-freeze termination arm); the default RNG derivation
`SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))` pinned end-to-end (previously
the default path was never executed — a `default_rng(cfg.seed)` regression
correlating with the ion stream would have stayed green); invariant tolerance
tightened 1e-2 → 1e-4 with the CLAUDE.md-required justification (measured
2e-5; the old band hid per-rung fold errors); `check_relaxation_config` arms
completed (`time ≤ 0`, `dt ≤ 0`, and the ν·dt = 0.1 accept-boundary).

**Corrections to the E2 delivery record above (annotation, append-only):** its
claim "every `check_relaxation_config` arm fires" was overstated (the
`time ≤ 0` / `dt ≤ 0` sub-arms and the boundary-accept were unfired — now
fired), and its test list silently omitted the plan's ΔE_dissip ≡ K2-drain
acceptance oracle (now covered).

**Recorded, not fixed (deferred with rationale):** the negative-ν hole in
`check_relaxation_config` (loud downstream via `rrk_rate`'s 2026-07-02
module-level defense; a config-arm would only move the error earlier); the
strict `<` freeze boundary vs `rrk_rate`'s k = 0 at equality (delays the early
exit by ≥ 1 step at a measure-zero boundary; never a wrong shed); the
`_decide_stride_ion` private-import carry (**re-opened**: E2's delivery
decision 2 promised its public promotion to E5b, but E5b promoted only
`drag_gate_steepness` — carried as an open Phase-F item); the stride > 1
storage branch and native-s (s = 60) nonzero-rate path remain unexercised
(production-scale concerns; the stride-1 assumption is now asserted loudly
inside the new oracles via the time-axis spacing check); `held_md_pot` /
`dt_relax` naming-suffix nits; the pre-loop `freeze_flags` initialization
(defensive init, kept); `b_ion_outside` stays a documented pass-through;
`cfg.seed = None` → entropy seeding stays the documented ion-stage convention.

**Verification.** Fail-first confirmed for all four new-behavior tests
(absolute-time 20.01, two guards DID-NOT-RAISE, 15 ≠ 14);
`test_relaxation_stage.py` 29 green, `test_ion.py` 26 green; full suite
**1952 passed** (1937 → +15), same 27 pre-existing warnings, 0 failures. No
schema bump, no RNG draw-order change (the default-path derivation is now
*pinned*, not changed), no drag-law / neutral / propagation-physics touch
(`_NUM_2N_T_ARRAYS_ION` is a storage-budget constant), no public-surface
rename (`time_relaxed_ps` keeps its pinned name; its semantics now match its
documented contract). §5 out-of-scope guard holds. **Phase E review is now
complete for all five slices.**

---

## Phase F — Slices F1 + F2 DELIVERED (2026-07-03)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (open questions
discussed and resolved by the user first, below). Phase F composes only accepted
A–E modules — **no new physics, no new `SimConfig` fields**; the new surface is the
campaign harness (F1) and the run-matrix generator (F2). Both are orchestration
scripts mirroring the delivered Tier-1a scaffolding.

**Scope decision — 9 A only (user, 2026-07-03).** The campaign runs the **9 A case
only**: there is no reference shell evolution for the 18 A droplet, so the plan's
9/18 A density-contrast route to `f_ret` (F2 "both cases per grid point"; F4 `f_ret`
discrimination) is **unavailable**. Consequence carried forward: `f_ret` is **not
identifiable** from the 9 A size distribution alone — F2 holds it at a prior (0.1),
not swept; F4's `f_ret` leg is a documented non-deliverable, not a computed result.
This narrows the plan's F2/F4 contract (recorded here; the plan text still describes
the two-case design).

**Open questions resolved (user, 2026-07-03, before coding):**
1. **F1 new-kwarg default convention → None-sentinel pass-through.** `picture`,
   `kappa`, `tau_ps` default to `None` and are only injected into `replace(...)`
   when set (ride the config default otherwise). Keeps the Phase-D bridge
   byte-reproducible and avoids duplicating the config defaults (`1.0`, `6.55`,
   `"statistical_mixture"`) in the script layer (Quality Principle 1).
2. **`relaxation_time_ps` = 8530 ns (experimental I⁺ relaxation time in this
   setup).** Stored as `EXPERIMENTAL_RELAXATION_TIME_PS = 8530·10³ = 8.53e6 ps`
   in `tier2_common` (config field is in ps). The E2 stage exits early on
   all-fragments-frozen, so the large cap is a ceiling, not the executed step
   count. (⚠ operational: a full-N ion run whose checkpoint stride drops the true
   final column trips E2's seed-coherence guard — re-run with larger `max_bytes`.)
3. **Relaxation translation arm → `coulomb`** (config default; the two I⁺
   fragments still repel post-dissociation). The E2 plan's "free-flight" wording
   referred to *pickup/drag* being off (which the stage already enforces:
   `λ₀=0`, drag frozen), not to killing I⁺–I⁺ Coulomb.
4. **F2 staging → manual USER-SETTINGS per stage.** F2 is a parameterized
   generator, not an auto-adjudicator: Stage 1 = `KAPPA_GRID × PICTURE_LIST` at
   `F_INT=None` (→ per-point self-unbound floor, the timing pin); the operator
   scores with F3/F4, hand-pins the winner, and re-runs enumerating `TAU_PS` for
   Stage 2. Matches `gen_tier1a_runs` + the "reported, not auto-adjudicated" gate.
5. **κ grid = {0.5, 1, 2, 4, 8}** (gradual → sharp cliff). Deliberately spans the
   cliff end per the Phase-D bridge finding (the under-shed is *kinetic*, the
   s−1 = 59 RRK exponent); if no κ/picture/τ in-band lands the 21→14 staircase,
   that surfaces as the RRK-dof mechanism-level OQ — not a silent ν/s retune.

**Build (F1 — `scripts/tier2_common.py` extended):** `build_biphasic_cfg` gains the
None-sentinel campaign kwargs (`picture`, `kappa`, `tau_ps`, `coulomb_available_eV`
default 0.80, `relaxation_time_ps` → enables E2, `relaxation_forces`); new module
constants `EXPERIMENTAL_RELAXATION_TIME_PS`, `VALIDATION_BUDGET_EV=0.80`,
`PRODUCTION_BUDGET_EV=2.70`; new knob-encoding helpers `tier2_run_tag` /
`tier2_run_dir_name` (encode budget `bNNN`, picture `mix`/`x2`/`cool`, κ/f_int/f_ret
to 2 dp, τ to 1 dp, + `_totalstrip` variant suffix; picture reject-arm). The bridge
tag helper `tier2_bridge_run_dir_name` / `TIER2_BRIDGE_TAG` stays reserved as-is; the
delivered bridge kwargs/defaults (`lambda0_per_ps`/`f_int`/`f_ret`) are unchanged.

**Build (F2 — `scripts/gen_tier2_runs.py`, new):** USER-SETTINGS staged generator,
9 A only, `N=500` single fixed seed; `campaign_grid_points()` (picture × κ, floor
`f_int` when `F_INT=None` via `internal_energy_budget.f_int_floor`),
`build_campaign()` (→ `(label, cfg, run_dir)`, no propagation),
`_run_one()` pipeline **neutral → ion (`biphasic_step`) → E2 relaxation**
(`run_relaxation_stage(..., save_path=<run_dir>/relaxation.npz)`), `OVERWRITE_
EXISTING_RUN` guard. `BUDGET_EV=0.80` first (F5 flips to 2.70).

**Tests:** `tests/test_tier2_common.py` (12 — bridge-reproducibility/ride-default,
campaign overrides, budget stamp, relaxation enable/forces, tag/dir round-trips +
picture reject + grid-uniqueness) and `tests/test_gen_tier2_runs.py` (4 — floor
`f_int` applied, `main()` schedules the expected grid, overwrite refusal, tiny-N
end-to-end writes all four artifacts). Full suite **1968 passed, 0 failed** (34
warnings = the documented §6.5 `anchored_discrete`/`biphasic` mass-pairing path).

**Rule-2 / scope:** no new declared-but-unread config fields (F1/F2 only *read* the
A–E knobs). No schema bump, no RNG draw-order change, no drag-law / neutral /
propagation touch. **Not built (per plan §2, deferred):** F3 scoreboard, F4
identifiability report, F5 production switch + total-strip runs, F6 figures; and the
R6 / p↔κ conditional triggers stay document-only. `_decide_stride_ion` public
promotion (the re-opened E2 carry) is untouched — F2 does not import it.

### Slices F1 + F2 — review + test-hardening pass (2026-07-05)

A workflow-backed high-effort code review (4 finder angles, 7 verifier agents,
9 candidates → 7 verified findings, 0 refuted) plus an independent domain review.
Five findings fixed, two declined with rationale; +net tests (full suite **2007
passed, 0 failed**). All fixes are in the F1/F2 script layer — no physics, no config
contract change.

**Fixed (run-dir naming was the load-bearing risk — F3 attributes each size
distribution *by dir name*, so a tag collision aborts a run or mis-scores it):**
1. **τ tag precision `.1f` → `.2f`** (CONFIRMED). The Stage-1 default τ = 6.55
   formatted to `tau6.5` and aliased a Stage-2 sweep value of 6.5 → identical
   run-dir → `FileExistsError` (or silent mis-attribution under overwrite). The
   Stage-2 τ sweep is a first-class plan step, so this was a real hazard.
2. **λ₀ now encoded in the tag** (CONFIRMED). `tier2_run_tag`/`tier2_run_dir_name`
   gained `lambda0_per_ps` (`_l{λ₀:.2f}`); two campaigns differing only in λ₀ no
   longer collide. (λ₀ is held fixed at 9 A, but the provenance hole was real.)
3. **Generator-level budget guard** (PLAUSIBLE). The config field
   `coulomb_available_eV` keeps its frozen "NO hard refuse" contract (plan §8); the
   *generator* now refuses a `BUDGET_EV` outside the sanctioned {0.80, 2.70} set,
   catching an operator typo (e.g. 8.0) before it silently poisons every RRK gate.
4. **Resume mode** (CONFIRMED). New `SKIP_COMPLETED_RUNS` (default True): a run
   dir with all four artifacts is skipped, so re-running a crashed 15-point grid
   recovers only the failed points instead of redoing every N≈500 pipeline. A
   *partial* dir is not "complete" → still trips the overwrite guard (never
   silently kept).
5. **`_PICTURE_TAGS` lockstep with the config enum** (PLAUSIBLE). An import-time
   assertion checks the abbreviation map covers exactly
   `config.LadderElectronicPicture` (`typing.get_args`), so a 4th picture added to
   the enum fails loudly here instead of silently blocking run-dir naming.

**Declined (with rationale):**
- **`_run_one` duplication vs the bridge script** (CONFIRMED cleanup). The bridge's
  `_run_one` is **frozen / byte-reproducible** (plan F1 NB: "do not repurpose");
  extracting a shared helper would force a touch to that frozen file. The Tier-1a
  multi-generator precedent already tolerates this thin-glue duplication, and the
  bridge won't change — drift risk is low. Left as-is.
- **Double `cfg.validate()` in `_run_one`** (CONFIRMED cleanup). Deliberate
  defense-in-depth: `_run_one` re-validates to catch a post-build `replace()`
  (exactly what the F2 smoke test does when it shrinks the neutral stage). Every
  Tier-1a `_run_one` validates the same way. Kept.

**New tests:** tag τ/λ₀ separation locks, `_PICTURE_TAGS`↔enum lockstep, the
grid-wide floor-`f_int` build guarantee (30 points × both budgets, floor ∈ (0,1]),
bridge-tag reserved regression, `resolve_f_int` literal override, relaxation-cap
default/override, the unsanctioned-budget refusal, and the skip-completed /
refuse-partial resume semantics. **Verification:** the production ion-writer
force-stores the final column (`ion.py:292-305`), so E2's seed-coherence guard is a
fail-loud *corner* (exact stride divisibility), not a systematic abort — the
`max_bytes` caveat stands as documented, no code change.

---

## Phase F — Slice F3 DELIVERED (2026-07-05)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (three open
questions discussed and confirmed by the user first, below). F3 is the deferred
Phase-E **scoreboard assembler**: a *pure scorer* composing only accepted E1/E4/E5
modules + the 5-term `ion_ledger_closure` into **one row per finished campaign run
dir**. It runs nothing (F2 already wrote `cfg.json`/`ion.npz`/`relaxation.npz`).
**No new physics, no new `SimConfig` fields.** TDD throughout (test→RED→GREEN: the
`ModuleNotFoundError` RED was watched before the module existed).

**Three open questions resolved (user, 2026-07-05, before coding):**
1. **Run-matrix enumeration → glob-discovery + read `cfg.json`.** Over a
   settings-mirror or a `build_campaign()` import: the staged campaign accumulates
   run dirs across stages (Stage-1 grid, then the Stage-2 τ sweep), so scoring
   *what exists on disk* keeps the scoreboard in sync without a settings block to
   hand-resync. Every row knob (budget/picture/κ/f_int/f_ret/τ/N) is read from the
   authoritative `cfg.json`, **never parsed from the tag**; only the `case` token is
   recovered from the dir-name prefix (`<case>_drag_…`, the F1/`tier0_common`
   convention). `budget_eV` filters by the run's own scenario stamp post-load.
2. **E5 `freeze_side_expected` → auto from budget.** `_freeze_side_expected(cfg) =
   isclose(cfg.coulomb_available_eV, VALIDATION_BUDGET_EV)` — `True` at 0.80 eV
   (Π>1 is the wiring smell), `False` at 2.70 eV (Π>1 is legitimate physics; the
   E5-review's "must not be carried there"). Compared against the sanctioned
   constant, not a raw float literal.
3. **Extra columns → add the three advisories.** Beyond the 16 pinned plan-F3
   columns: `t_cross_spread_ps`, `t_cross_all_agree`, `sanity_flags`
   (semicolon-joined). Strictly additive; feeds F4's identifiability report.

**Build (`scripts/post_processing/tier2_size_distribution_table.py`, new)**, mirroring
`tier1a_rmse_table.py` (`collect_*_records`/`collect_*_rows`/`score_*_run`/
`format_table`/`write_rows_csv`, USER-SETTINGS block, `main()`):
- **Per run** (`_score_run`): E1 `compute_terminal_shell_distribution` **twice** —
  matched-time from `relaxation.npz` (`source_tag="relaxed"`, the E1 reload-pitfall
  override — a reloaded relaxation checkpoint is a v7 `IonCheckpoint` and would
  else infer `sim_end`) and sim-end upper bound from `ion.npz`
  (`source_tag="sim_end"`); E4 `compare_size_distributions(..., metric=
  cfg.validation_histogram_metric)` for both → `W1_matched`/`W1_simend_upper`; E5
  `reconstruct_diagnostics(ion, cfg, freeze_side_expected=…)` on the **ion** stage
  (t×/Π need the ion trajectory) → `t_cross_ps`(=median)/`regime_label`/
  `total_strip_reachable` + advisories; `ion_ledger_closure(ion).max_abs_residual_eV`;
  and `n_terminal_mean`/`n_terminal_spread` as moments of the **matched (relaxed)**
  distribution.
- **Semantic note (documented in the module header):** `regime_label` reads the
  *ion-end* (R5 upper-bound) n by construction (E5 runs on the ion stage), distinct
  from the relaxed `n_terminal_*` columns — two sources, no conflict.
- `discover_run_dirs` globs `*_tier2_*` dirs with all three read-artifacts;
  `collect_tier2_records`/`_rows` load the reference once, apply the budget filter,
  score each. **Reported, not auto-adjudicated** — no code asserts a fidelity verdict.
- `load_ion_checkpoint(run_dir/"relaxation.npz")` for the matched checkpoint
  (`RunDirectory` only knows `ion.npz`).

**Tests:** `tests/test_tier2_size_distribution_table.py` (**13**) — one genuine
tiny-N (N=2) biphasic run dir built once (module-scoped `neutral→ion→relaxation`),
reused by `copytree`: row carries exactly the 16+3 columns; W₁ both finite/≥0;
ledger residual finite; regime label ∈ the two outcomes; terminal-n moments in
`[0, n*]`; knob columns echo cfg; `RunDirectory`-or-path input; discovery finds
complete runs / skips a run missing `relaxation.npz` / ignores non-`tier2` dirs;
collect scores all discovered + budget filter selects/drops; freeze-side True@0.80
/ False@2.70; CSV round-trips; empty-rows table. `main()` is the manual operator
entry (thin glue over tested units; exercised non-test, the tier1a precedent).

**Rule-2 / scope:** no new declared-but-unread config fields; no schema bump, no RNG
draw-order change, no drag-law / neutral / propagation touch. **Not built (per plan
§2, deferred):** F4 identifiability report, F5 production switch + total-strip runs,
F6 overlay figures (the plan hosts F6 in this module behind `SHOW_FIGURE` — no figure
code yet); R6 / p↔κ conditional triggers stay document-only.

**Verification.** F3 suite **13 passed**; full suite **2020 passed, 0 failed** (2007
→ +13), 69 warnings = the documented §6.5 `biphasic`/`anchored_discrete` mass-pairing
path.

### Slice F3 — review + test-hardening pass (2026-07-05)

An adversarial self-review of the F3 diff (correctness/semantics, fail-loud
coverage, edge cases, reuse/dead-code, test adequacy) surfaced **two real findings**;
both fixed test-first (behavioral RED watched before the fix), plus seven
coverage-lock tests on already-correct behavior.

**Fixed:**
1. **`freeze_side_expected` ignored the `total_strip` variant (correctness).** It
   keyed E5's Π>1 flag purely off budget, but a **0.80 eV `total_strip` run sheds
   by construction** (gate never self-binds), so Π>1 is legitimate there and the
   freeze-side flag was spurious. `_freeze_side_expected(cfg, *, total_strip=False)`
   now returns `False` for a total-strip run at any budget; `_is_total_strip_run`
   reads the F1 `_totalstrip` dir-name suffix (the variant has no `SimConfig` field
   — the marker lives only in the run tag, symmetric to the `case` prefix parse).
2. **`discover_run_dirs` could sweep the reserved bridge run (robustness).** The
   `*_tier2_*` glob also matches the Phase-D `tier2_bridge_biphasic` dir; today it
   is excluded only incidentally (the bridge predates E2 → no `relaxation.npz`), so
   a future relaxation-enabled bridge run would silently enter the campaign
   scoreboard. Now excluded explicitly by `TIER2_BRIDGE_TAG` (single source,
   imported from `tier2_common`).

**Coverage locks added (+9 tests, F3 suite 13 → 22):** the two fixes; the
**source-tag wiring** (a `compare_size_distributions` spy proves matched W₁ scores
the *relaxed* distribution and W1_simend the *sim_end* one, in order — the exact
mis-pairing the E-review flagged as F3's hazard); `_distribution_moments` hand
oracle (two-spike → mean 1 / std 1); `_case_from_run_dir` recover + malformed-name
raise; `_is_total_strip_run` suffix read; a real scored row (t× may be NaN)
surviving the CSV round-trip; the budget filter separating **two** runs at
different cfg stamps; and `main()` returning 0 on an empty root.

**Deliberate non-changes (reviewed, left as-is):** `_score_run` propagates a
per-run E5/load raise rather than swallowing it — fail-loud (Principle 4) over a
resilient batch mode that could mask a corrupt run; the R5 semantic split (E5
regime on ion-end n vs relaxed `n_terminal_*`) is intended and documented; `None`
knobs render as empty CSV cells / `-` in the text table (a diagnostic table, not a
data contract).

**Verification.** F3 suite **22 passed**; full suite **2029 passed, 0 failed**
(2020 → +9), 70 warnings = the documented §6.5 mass-pairing path. Next: **F4**
(identifiability + regime-determination report over the F3 rows).

---

## Phase F — Slice F4 DELIVERED (2026-07-05)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (open questions
discussed and resolved by the user first, below). F4 is the campaign's **central
deliverable**: a *pure assembler* that turns the F3 scoreboard rows into a
**per-quantity identifiability statement** — which of the CALIBRATION_MAP
"Tier-2 load-bearing" quantities the single size-distribution observable can
actually *separate*. It runs nothing and scores nothing new; it reads F3's typed
records and classifies. **No new physics, no new `SimConfig` fields.** TDD
throughout (test→RED `ModuleNotFoundError` watched before the module existed).

**Open questions resolved (user, 2026-07-05, before coding):**
1. **Input contract → import F3's typed collector, not the CSV.** `collect_tier2_
   records` (typed `Tier2RunRecord`) over round-tripping the lossy scoreboard CSV
   (`None`→empty, NaN). One scoring path; the assessment core operates on the F3
   row dicts so it is testable on constructed synthetic rows.
2. **τ-insensitivity → concrete threshold + graceful single-τ degradation.**
   `TAU_SENSITIVITY_THRESHOLD_HE = 0.5` (terminal-n moves > 0.5 He across the τ
   band → `sensitive`, a live dimension; else `insensitive`, confirmation). The
   Phase-D bridge expects the *sensitive* arm at 0.80 eV, so sensitivity is a
   first-class reported outcome, not an anomaly. At F4's first run only Stage-1
   exists (single τ per cell) → `insufficient_data` (`Stage-2 τ sweep not yet
   run`), never an error.
3. **"Reported, not auto-adjudicated" boundary.** Factual minima are allowed
   (argmin κ, best picture, best-W1 regime); a fidelity *verdict* is not. F4
   reports the W₁-vs-κ landscape **and** its argmin, plus whether that minimum is
   well-separated (→ identifiable) or flat/degenerate (→ not); it never declares
   calibration success. The flatness/bracketing/sensitivity thresholds are
   documented as reader aids, not acceptance gates.
4. **Boundary argmin → `not_bracketed`.** If the W₁-vs-κ minimum sits at a grid
   edge (κ=0.5 or 8 — the Phase-D bridge predicts a push toward the sharp-cliff
   end), F4 reports `not_bracketed` (`optimum not enclosed → extend grid / RRK-dof
   mechanism-level OQ`) rather than a clean identifiable κ. Flatness is checked
   *before* bracketing so a degenerate landscape is never mislabelled an edge
   optimum.
5. **9 Å-only scope carried from F1/F2 (2026-07-03).** `f_ret` and `λ_attach` are
   rendered as **deterministic non-deliverables**, not computed: `f_ret`
   `not_identifiable` (the only route is the 9/18 Å density contrast; no 18 Å
   shell reference exists — held at prior 0.1); `λ_attach` `not_identifiable`
   (pickup structurally dead in-window at 9 Å, Π ≤ 0.005, Phase-D finding #4). The
   plan's F4 test wording ("9/18 Å contrast computes") is corrected to "F4 reports
   it as a non-deliverable."
6. **f_int / ν / s deterministic verdicts.** `f_int` `not_identifiable`
   (timing-only knob: moves t×, not the cascade budget → terminal n insensitive by
   construction); `ν`/`s` `held_fixed`; `s` carries the **s↔κ coupling** caveat as
   prose (the κ fit is conditional on the fixed RRK dof), echoed on the κ entry.
7. **Per-budget sectioning.** 0.80 eV (validation, freeze-side) and 2.70 eV
   (production, legitimate Π>1) carry different identifiability meaning; the
   checklist and text report are sectioned per budget. At F4's first run only
   0.80 exists; 2.70 arrives with F5 and degrades to `insufficient_data`.

**Build (`scripts/post_processing/tier2_identifiability_report.py`, new)**, mirroring
the Tier-1a/F3 assembler idiom (`build_*`/`format_*`/`write_*_csv`/`collect_*`/
USER-SETTINGS block/`main()`):
- **Status vocabulary** — `identifiable` / `not_identifiable` / `not_bracketed` /
  `insufficient_data` / `held_fixed` / `sensitive` / `insensitive`; `identifiable`
  + `sensitive` count as "separated" for the headline.
- **Assessment primitives** (pure, synthetic-row-testable): `assess_kappa_
  landscape` (insufficient → flat → not_bracketed → identifiable), `assess_picture_
  separation`, `assess_tau_sensitivity`, `assess_regime` (label counts +
  lowest-W₁ regime, NaN-guarded).
- **Grouping builders** — `build_kappa_landscapes` (cells keyed by budget/picture/
  f_int/f_ret/τ with ≥2 κ → the Stage-1 co-fit cell), `build_tau_sweeps` (cells
  with ≥2 τ). κ is assessed under the **globally best picture** (the co-fit's
  chosen picture); τ under the best cell if a sweep exists there.
- **Checklist** — `build_identifiability_rows` emits one row per (budget,
  load-bearing quantity) over `LOAD_BEARING_QUANTITIES` (κ, picture, τ, regime,
  f_int, f_ret, λ_attach, ν, s); computed first, deterministic after. CSV columns
  `budget_eV, quantity, status, value_or_range, evidence, notes`.
- **Text report** — `format_report` sections per budget with a "K of M separated"
  headline; **reported, not auto-adjudicated** banner.

**Tests:** `tests/test_tier2_identifiability_report.py` (**25**) — all on constructed
synthetic scoreboard rows (the plan's stated F4 surface; no run pipeline, no
figures): κ-landscape interior-min/boundary/flat/single-point; picture
separation/flat/single; τ sensitive/insensitive/single; grouping builders;
regime best-W₁; checklist covers all 9 load-bearing quantities once per budget;
κ+picture identifiable on a clean co-fit grid; static non-deliverable verdicts
(f_ret/λ_attach 9 Å rationale, f_int timing-only, ν/s held + s↔κ note); τ sensitive
when swept / insufficient single-τ; per-budget sectioning; headline/section text;
CSV round-trip; `main()` returns 0 on empty root. One RED→GREEN adjustment (moved
"timing" into the f_int `notes`).

**Rule-2 / scope:** no new declared-but-unread config fields (F4 only *reads* F3
rows + the A–E knobs via cfg); no schema bump, no RNG draw-order change, no
drag-law / neutral / propagation touch. **Not built (per plan §2, deferred):** F5
production switch (2.70 eV) + total-strip secondary runs, F6 overlay figures
(the plan hosts F6 in this module + the F3 module behind `SHOW_FIGURE` — no figure
code yet); R6 / p↔κ conditional triggers stay document-only.

**Verification.** F4 suite **25 passed**; full suite **2054 passed, 0 failed**
(2029 → +25), 70 warnings = the documented §6.5 mass-pairing path. Next: **F5**
(production switch to 2.70 eV + total-strip secondary runs), gated on the 0.80 eV
validation landing (Phase-D bridge cross-check).

### Slice F4 — review + test-hardening pass (2026-07-05)

A workflow-backed high-effort code review (4 finder angles → 21 candidates, 12
independent verifier agents, 21 verified → 10 kept after consolidation, 3
refuted) plus an independent domain pass. **Seven findings fixed, one declined
with rationale; +16 tests** (F4 suite 25 → 41). All fixes are in the F4 script
layer — no physics, no config-contract change, no schema/RNG/drag-law touch.

**Fixed (correctness — the load-bearing risk was mis-attributing the campaign's
central κ/picture/τ verdicts):**
1. **κ verdict no longer falls back to a non-best picture's landscape**
   (CONFIRMED). `_kappa_checklist_entry` used `min(chosen or landscapes, …)`; when
   the W₁-best picture had only single-κ probe runs (no landscape), `chosen` was
   empty and κ* was silently read off a *different* picture's landscape while the
   headline named the best picture. **Root fix:** the picture comparison
   (`_best_w1_by_picture`) is now derived from the **κ landscapes** (min over each
   picture's co-fit cells), not raw rows — so the globally best picture always owns
   a landscape (no cross-picture fallback) and κ/picture share one source. This
   *also* closes an independent domain finding: the old raw-row picture min was
   **contaminated** by engineered `total_strip` runs and Stage-2 τ-sweep rows
   (distinct f_int/τ → single-κ cells → not landscapes), which could flip the
   reported best picture; the landscape-derived comparison structurally excludes
   them.
2. **`assess_kappa_landscape([])` no longer crashes** (CONFIRMED). The argmin ran
   before the `len < 2` guard, so an empty landscape raised `ValueError: min() arg
   is an empty sequence` instead of the documented `insufficient_data`. Guard moved
   ahead of the argmin; empty → `insufficient` with a NaN κ/W₁ sentinel.
3. **τ is read at the co-fit optimum κ cell** (PLAUSIBLE→fixed). `_tau_checklist_
   entry` picked `(preferred or sweeps)[0]` — the first τ sweep under the best
   picture in *insertion order*, ignoring κ. With sweeps at two κ cells (an
   insensitive off-optimum and the sensitive optimum), it could report the wrong
   arm. `best_kappa` (the κ argmin) is now threaded through and the sweep at
   (best picture, best κ) is preferred, with a note when no sweep sits at the
   optimum cell.

**Fixed (robustness / cleanup):**
4. **Duplicate knob-cell collision now warns** (PLAUSIBLE). Two complete run dirs
   at an identical (budget,picture,f_int,f_ret,τ,κ) cell collapsed last-write-wins
   silently. `_accumulate_cell` now emits a `RuntimeWarning` when the incoming
   W₁/terminal-n differs from what is already stored (identical values = a harmless
   deterministic replay, no warn), surfacing a stale/replicate dir instead of
   dropping data silently.
5. **NaN guards on both grouping builders** (CONFIRMED, defensive). `build_kappa_
   landscapes` / `build_tau_sweeps` read `W1_matched` / `n_terminal_mean` with a
   raw `float()` while the sibling consumers (`_best_w1_by_picture`,
   `assess_regime`) filter via `_is_number`. Added the same guard for internal
   consistency (a NaN terminal-n would else make the τ spread NaN → silently
   `insensitive`). *Note:* the reviewers separately **refuted** three variants of
   this as reachable — E1/E4 fail loud on non-finite upstream, so the guard is
   inert in practice; kept purely as defense-in-depth/consistency.
6. **`main()` no longer builds the checklist twice** (CONFIRMED). It called
   `collect_identifiability_rows` (which builds it) *and* `format_report(rows)`
   (which rebuilt it). `format_report` now accepts a pre-built `checklist=`; `main`
   passes the one it already has.

**Declined (with rationale):**
- **`_format_cell` / `write_report_csv` duplicate F3's `_format_value` /
  `write_rows_csv`** (CONFIRMED cleanup). Display-only helpers, not physics —
  Quality Principle 1 targets shared *physics/formulas/constants*, not a 6-line
  formatter; the two CSV writers differ in their column set
  (`IDENTIFIABILITY_COLUMNS` vs `TIER2_TABLE_COLUMNS`), so sharing them needs a new
  column-parameterized util. The Tier-1a/F3 sibling-script family already tolerates
  this thin-glue duplication (the F1/F2 review declined the analogous `_run_one`
  extraction), and each post-processing script is deliberately self-contained. Left
  as-is; a divergence in display precision would be a per-script choice, not a bug.

**New tests (+16, F4 suite 25 → 41):** the empty-landscape guard; low-edge
`not_bracketed` + all-zero-W₁ / zero-min-interior κ branches; the two
`assess_picture_separation` zero-best-W₁ branches; NaN-skip on `assess_regime` /
both builders; the duplicate-cell warn + identical-replay-silent pair; picture
best **not** contaminated by a non-co-fit run + κ attributed to the actual best
picture; τ read at the best-κ cell (not insertion order); `format_report`
pre-built-checklist equivalence; and the `collect_identifiability_rows` F3-record
wrapper wiring (monkeypatched, no real run).

**Verification.** F4 suite **41 passed**; full suite **2070 passed, 0 failed**
(2054 → +16), 70 warnings = the documented §6.5 mass-pairing path (the new
duplicate-cell warning is captured in-test, does not leak). Next: **F5**
(production switch to 2.70 eV + total-strip secondary runs), gated on the 0.80 eV
validation landing (Phase-D bridge cross-check).

---

## Pre-F5 — Staircase capability probe DELIVERED (2026-07-05)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger, after a design
discussion (user question: before production, can the biphasic mechanism
reproduce the 9 Å reference shell evolution — and is it flexible enough to
express various solvation shells after relaxation?). Design doc:
`TIER2_STAIRCASE_PROBE_PLAN.md` (inserted between F4 and F5; modifies no
delivered F-slice). TDD throughout (the `ImportError` RED on both new test
modules watched before any module existed).

**What it is.** A 45-point small-N existence probe extending the Phase-D
bridge's single-point staircase comparison across the full in-band lever space
(`TIER2_PHASE_D_BRIDGE_FINDINGS.md` §2 levers): κ ∈ {0.5, 1, 2, 4, 8} ×
picture (3) × τ ∈ {2.6, 6.55, 16.5} ps, at 9 Å / 0.80 eV / N=50 / bridge
seed+window, pipeline neutral → ion (`biphasic_step`) → E2 relaxation. Scored
against `build_shell_schedule` (zero artifact dependence) on staircase fidelity
(`dn_mean_window` vs 7, `n_ion_end_mean` vs 14, `t_first_shed_ps` vs the
anchored first event, `frac_ions_shed`, trajectory MAD) plus the relaxed
terminal-n moments (flexibility map) and the 5-term ledger residual.
**Reported, not auto-adjudicated** (factual argmin named, no verdict); a
grid-wide miss surfaces the RRK-dof mechanism-level OQ, a landing region points
the N=500 Stage-1/2 campaign — neither discharges the F5 gate by itself.

**Design decisions (user, 2026-07-05, before coding):**
1. **Goals: both** — staircase reproduction *and* the terminal-shell
   flexibility map, in one sweep.
2. **Small-N probe first** (N=50, ≈3/10 of one 15-point N=500 stage) before
   any campaign/production spend.
3. **Full κ×picture×τ grid from the start** (not staged): the bridge predicts
   τ-*sensitivity* at 0.80 eV, so a capability verdict with τ pinned mid-band
   would be unsound.
4. **Standalone probe scripts** in a disjoint run-dir namespace — no touch to
   the reviewed F2/F3 surfaces.
5. **`f_int = 0.5` (bridge pin), NOT the F2 Stage-1 floor** — f_int is
   timing-only (moves t×, not the Σ(21) cascade budget); 0.5 puts t× ≈ 5 ps on
   the anchored first shed, the floor would put gate-open at t ≈ 0 and corrupt
   every timing comparison. *(Flagged in the design presentation; the user
   answered with the implementation trigger while the 0.5 pin stood — flipping
   to the floor is a one-line USER-SETTINGS change if wanted.)*

**Build:**
- `scripts/tier2_common.py` — additive `tier2_probe_run_tag` /
  `tier2_probe_run_dir_name` (prefix `tier2probe_`, same two-decimal knob
  encoding as the campaign tag, picture reject-arm; no `total_strip` variant).
  **Namespace lock:** a probe dir contains no `_tier2_` substring, so the F3
  campaign glob can never sweep it, and vice versa — locked by test in both
  directions. Bridge/campaign helpers and defaults untouched.
- `scripts/gen_tier2_staircase_probe.py` (new) — USER-SETTINGS generator
  mirroring `gen_tier2_runs.py` (`probe_grid_points` / `build_probe` /
  `_run_one` with skip-completed + overwrite guards). Generator-level guard:
  **0.80 eV only** (the anchored staircase does not apply at any other budget,
  including the sanctioned 2.70). The (κ=1, mixture, τ=6.55) grid point is the
  Phase-D bridge configuration — its ion stage must reproduce the bridge
  numbers (Δn̄ ≈ 0.67, t× ≈ 4.96 ps), a free end-to-end wiring oracle (manual
  check at probe-run time, not a pytest).
- `scripts/post_processing/tier2_staircase_probe_report.py` (new) — pure
  scorer mirroring the F3 idiom (`discover_probe_run_dirs` /
  `score_probe_run` / `collect_probe_rows` / `format_table` /
  `format_headline` / `write_rows_csv`, 16 columns). Staircase metrics read
  the **ion** stage (R5: the anchored staircase is in-window), flexibility
  moments read the **relaxed** checkpoint (E1 `source_tag="relaxed"` reload
  override). Figures (mean-n(t) small multiples + relaxed-n heat map) are
  gated behind `SAVE_FIGURES`/`SHOW_FIGURES` with a lazy matplotlib import —
  never touched by pytest.

**Tests:** `tests/test_tier2_staircase_probe.py` (**11** — tag encoding/reject,
two-way namespace disjointness incl. F3-discovery-ignores-probe-dir, 45-point
unique grid, bridge-knob pins + budget/relaxation stamps, wiring-oracle-point
config equality vs the pure bridge build, non-validation-budget refusal,
partial-refuse/complete-skip guards, `main()` schedule, tiny-N end-to-end four
artifacts) and `tests/test_tier2_staircase_probe_report.py` (**9** — the
staircase-metric hand oracles: an ion walking the anchored staircase scores
exactly (7 sheds, terminal 14, MAD 0, first shed = first event time), mixed
and no-shed ensembles (NaN first-shed never a crash), a genuine tiny-N scored
row with knob echo + finite metrics, discovery namespace/incomplete/campaign
exclusion, CSV round-trip, empty table, empty-root `main()`).

**Rule-2 / scope:** no new `SimConfig` fields (the probe only *reads* A–E
knobs through the delivered F1 builder), no schema bump, no RNG draw-order
change, no drag-law / neutral / propagation touch; F2/F3/F4 surfaces
untouched. The 45-run probe itself is an **operator action** (run
`gen_tier2_staircase_probe.py`, then the report script) — not executed as part
of this delivery.

**Verification.** New suites **11 + 9 passed**; full suite **2090 passed,
0 failed** (2070 → +20), 75 warnings = the documented §6.5 mass-pairing path.
Next: run the probe, adjudicate (landing region → point the N=500 Stage-1/2
campaign; grid-wide miss → surface the RRK-dof OQ before F5).

### Staircase probe EXECUTED — outcome (b), RRK-dof OQ fired (2026-07-05)

The 45-run probe was executed and scored (operator run; scoreboard via
`tier2_staircase_probe_report.py`). **Wiring oracle passed:** the (κ=1,
mixture, τ=6.55) point reproduces the Phase-D bridge exactly (Δn̄ = 0.67,
n_end = 20.33); ledger residual uniform 2.23·10⁻⁵ eV across all 45 runs. The
result is therefore the mechanism's genuine prediction.

**Result — outcome (b) of the probe plan §6: no in-band point lands.** Max
reach Δn̄ = 1.7 sheds (at κ=0.5, τ=16.5, `cooling_relaxed`) vs the anchored 7;
relaxed terminal n confined to **[19.3, 20.95]** across the whole grid
(flexibility verdict: the pinned mechanism expresses essentially one shell).

**Three structural findings from the table + ladder/evaporation code:**
1. **The κ lever is inverted — corrects bridge findings §2 lever 1.** Δn̄
   *decreases* monotonically with κ everywhere. The Form-U cliff centre is
   pinned at n*+½ = 21.5 and the stripping range 21→14 lies *below* it:
   σ(21) = sigmoid(−κ/2) ≤ ½ always, so large κ drives the in-shell rungs
   *up* toward D₀(1) (it is the irrelevant n ≥ 22 rungs that shrink). The
   κ→0 direction is also normalisation-capped: **Form-U structurally floors
   D₀(21) at ≈ 0.53·D₀(1)** (shallow minimum near κ ≈ 0.25). No κ can make
   the near-edge rungs cheap. (A de-drift NB was added to the findings doc.)
2. **The picture knob is near-degenerate for this observable.** D₀(n) and
   Σ(n) both scale with (D₀(1) − D_floor), so the RRK argument
   x = D₀(21)/Σ(21) ≈ 0.032 is nearly picture-invariant — the three picture
   blocks of the scoreboard are almost identical. Identifiability note for
   F4: the κ×picture co-fit landscape is flat in picture.
3. **τ is the strongest but capped lever.** Δn̄ grows ≈ linearly over the
   band (0.2 → 1.7), but t× = τ·ln(f_int·E/Σ(21)) grows with τ; the ≈ 80 ps
   that 7 sheds would extrapolate to puts gate-open past the 30 ps window.
   Dead end in-band and out.

**Diagnosis (kinetic, quantitative):** k = ν·(1−x)^(s−1) with x ≈ 0.032 and
s−1 = 59 → k ≈ 0.15·ν ≈ 0.36/ps — the freeze is entirely the RRK exponent on
a small x. The swept levers move x at most ~2× (linear); s acts
exponentially (s_eff−1 = 11 → k ≈ 1.7/ps). Energetics are self-sustaining by
construction (per shed, E_int and Σ(n) drop by the same D₀(n): gate margin
shed-invariant; 7 sheds cost 0.06 of the 0.188 eV budget). **This fires the
pre-registered RRK-dof mechanism-level OQ** (bridge findings §2 lever 3):
s = 3n−3 treats the floppy quantum He₂₁ shell as 60 fully-coupled classical
oscillators.

**Tier-2-blocking consequence:** the crossing construction caps the cascade
budget at Σ(21) at *any* `coulomb_available_eV` (budget moves only t×), and
the experimental abundance reference holds 43 % bare I⁺ / <1 % weight at
n ≥ 19 — the freeze at n ≈ 20 blocks the arbitration observable at 0.80
**and** 2.70 eV. The probe saved the N=500 campaign from a dead band.

**Decision (user, 2026-07-05):** proceed with the **`s_eff` mini-probe**
design addendum — `TIER2_STAIRCASE_PROBE_PLAN.md` Addendum A (12 runs:
s_eff ∈ {None, 30, 20, 12, 8, 5} × τ ∈ {6.55, 16.5} at the bridge pins,
via the already-plumbed `cfg.evap_rrk_dof` override; no new physics, no new
`SimConfig` fields). Falsification/localisation only — a landing s_eff does
NOT auto-promote the convention; promotion (dof-convention enum arm +
Bounded s_eff knob + CALIBRATION_MAP reclassification + F2 re-scope to an
s_eff co-fit) is a subsequent user-level mechanism decision. Implementation
stays behind `[PROCEED TO IMPLEMENTATION]`.

### s_eff mini-probe harness DELIVERED (2026-07-06)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (user: "Please
implement the 12-run miniprobe"). TDD throughout (9 RED failures watched — the
missing `evap_rrk_dof` kwargs, `S_EFF_GRID`, 4-tuple grid, and `s_eff` column
— before any source change). **Additive only: no new physics, no new
`SimConfig` field** — `cfg.evap_rrk_dof` was already declared and physics-live
at the driver's evaporation step (`ion_propagation_step.py:630`), so this is
pure orchestration/report wiring plus the config pass-through.

**Build (all additive to the delivered probe scripts, Addendum A.3):**
- `scripts/tier2_common.py` — `build_biphasic_cfg` gains a None-sentinel
  `evap_rrk_dof` kwarg (injected only when set; `None` rides the config
  default — the biphasic bridge/probe builds stay byte-identical);
  `tier2_probe_run_tag` / `tier2_probe_run_dir_name` gain `evap_rrk_dof`,
  appending `_sNN.NN` (two decimals) when set and **nothing** when `None`
  (delivered 45-run tags byte-identical).
- `scripts/gen_tier2_staircase_probe.py` — new `S_EFF_GRID` dimension;
  `probe_grid_points()` now returns 4-tuples `(picture, κ, τ, s_eff)`;
  `build_probe` threads `evap_rrk_dof` into the cfg + dir name + label +
  operator print. **Active USER SETTINGS flipped to the mini-probe** (κ=1,
  mixture pinned; τ ∈ {6.55, 16.5} × `S_EFF_GRID = [None, 30, 20, 12, 8, 5]`
  → **12 runs**), since the 45-run probe already executed; the full 45-run
  grid is preserved in a "FULL PROBE" comment block (one-swap + `S_EFF_GRID=
  [None]` restore). The wiring-oracle point is now (κ=1, mixture, τ=6.55,
  **s_eff=None**) — the per-n control that must still equal the bridge.
- `scripts/post_processing/tier2_staircase_probe_report.py` — `s_eff` column
  (read from `cfg.json`, `None` → per-n); figure NB that the relaxed-n heat
  map collapses an s_eff sweep (table/CSV authoritative), s_eff added to the
  mean-n(t) small-multiple curve labels.

**Tests (+5 net; probe suites 20 → 25):** `test_tier2_staircase_probe.py`
gained the s_eff tag suffix/omit + namespace check, `build_biphasic_cfg`
stamping + None-default, the 12-point mini-probe default grid, s_eff dir
non-collision, and the **`test_full_grid_back_compat_none_s_eff`** guard (full
grid + `S_EFF_GRID=[None]` → 45 points and byte-identical tag); the delivered
grid tests were updated to the mini-probe config + 4-tuple arity + the
`evap_rrk_dof is None` filter on the wiring oracle.
`test_tier2_staircase_probe_report.py` gained the `s_eff`-column read (per-n
None on the existing tiny run + a set `s_eff=8` tiny run).

**Rule-2 / scope:** no new declared-but-unread config field (the kwarg feeds an
already-live one), no schema bump, no RNG draw-order change, no drag-law /
neutral / propagation / driver touch (the driver already reads
`cfg.evap_rrk_dof`). F2/F3/F4 untouched.

**Verification.** Probe suites **25 passed**; full suite **2095 passed,
0 failed** (2090 → +5), 77 warnings = the documented §6.5 mass-pairing path.
Next (operator action): run the 12-point mini-probe
(`python scripts/gen_tier2_staircase_probe.py`, then the report script),
confirm the s_eff=None point reproduces the bridge freeze (wiring oracle), and
adjudicate the RRK-dof OQ (a landing s_eff → the promotion decision;
no constant s_eff lands → the ladder/gate structure becomes suspect).

#### Review + fix pass (2026-07-06, before the operator run)

A workflow-backed high-effort code review (4 finder angles, 9 candidates → 9
verified, 1 refuted) over the uncommitted probe diff. **Three fixed, one
declined-with-guard; +net tests; full suite 2095 → 2097 passed, 0 failed.**

**Fixed (fidelity — the mini-probe's whole point is the s_eff sweep):**
1. **Headline argmin dropped s_eff (CONFIRMED).** `format_headline`'s
   "lowest trajectory MAD" line printed only picture/κ/τ — identical across all
   12 mini-probe runs (which share those three), so the operator could not read
   *which* s_eff minimised the MAD. Now prints `s_eff=<n|per-n>` via a new
   `_s_eff_label` helper; two headline tests lock it (a numeric winner and the
   per-n control).
2. **Heat-map collapsed the s_eff sweep (CONFIRMED).** `emit_figures` panel (b)
   keyed cells on `(kappa, tau)` only, so the mini-probe's 6 s_eff values at a
   fixed `(κ, τ)` overwrote one cell (last-write-wins, 5 of 6 hidden). The
   column axis is now the **swept dof knob** — κ when it varies (the delivered
   45-run probe, byte-identical), else s_eff. Figure code stays non-test/gated.
3. **`_distribution_moments` duplicated (CONFIRMED, Principle 1).** The probe
   report had copied the terminal-n mean/population-std formula verbatim from
   the F3 scoreboard — the *same reported quantity from the same checkpoint*, so
   a future spread-convention change could silently diverge them. Hoisted to
   **`ShellDistribution.moments()`** (`postprocess/size_distribution.py`, the
   owning module); both report scripts now call `dist.moments()`, the F3 hand
   oracle moved onto the method. (Distinct from the F4-review-declined *display*
   helper duplications: this is a shared physics/statistics formula, squarely
   Principle 1.)

**Declined-with-guard:**
4. **Probe vs campaign tag-encoder duplication (PLAUSIBLE).** The reviewer's
   deeper fix (a shared `_encode_knob_tag(prefix, suffix)`) would touch the
   **frozen, delivered campaign encoder** (`tier2_run_tag`) right before an
   operator run — against the F1/F2 "do not repurpose the delivered encoder"
   stance, and tag-format naming is a convention, not the physics/formula/
   constant Principle 1 governs (the F4 review declined the analogous script
   duplications on the same ground). Instead added a **parity-lock test**:
   `tier2_probe_run_tag(knobs) == tier2_run_tag(knobs).replace("tier2_",
   "tier2probe_", 1)` — so a future change to either encoder's knob body fails
   loudly, catching the exact divergence the finding names, with zero change to
   the frozen encoder.

**Refuted (agreed):** the `staircase_metrics` shed-below-start detection missing
a pre-evaporation pickup — inert at 9 Å / 0.80 eV (pickup structurally dead,
Π ≤ 0.005), confirmed by the verifier and matching the pre-build reasoning.

**Scope:** fixes touch only the two report scripts + the `ShellDistribution`
method (owning module); no new physics, no schema/RNG/drag-law/driver change.
Verification: full suite **2097 passed, 0 failed** (2095 → +2: the per-n
headline + tag-parity lock tests).

### s_eff mini-probe EXECUTED — outcome (a): a constant s_eff lands the staircase (2026-07-06)

The 12-point mini-probe was executed and scored (operator run;
`tier2_staircase_probe_report.py` over 55 probe dirs = the 45-run grid + the
10 new set-`s_eff` points; the two `s_eff=None` controls reuse the delivered
bridge-point dirs). **Wiring oracle passed:** the (κ=1, mixture, τ=6.55,
s_eff=None) point reproduces the Phase-D bridge exactly (Δn̄ = 0.67,
n_end = 20.33); ledger residual uniform 2.23·10⁻⁵ eV across all 55 runs.

**Result — outcome (a) of `TIER2_STAIRCASE_PROBE_PLAN.md` Addendum A.4: a
constant s_eff lands the staircase in magnitude AND timing.** The τ = 6.55 ps
arm (targets: 7 sheds → n_end 14, first shed at the anchored t★ = 5 ps):

| s_eff | Δn̄ | n_ion_end_mean | t_first_shed_ps | n_traj_MAD |
|---|---|---|---|---|
| per-n (60 at n=21) | 0.67 | 20.33 | 6.26 | 3.99 |
| 30 | 2.08 | 18.92 | 5.97 | 3.02 |
| 20 | 3.51 | 17.49 | 5.70 | 2.23 |
| 12 | 5.61 | 15.39 | 5.51 | 1.25 |
| **8** | **7.41** | **13.59** | **5.42** | **1.00** (grid argmin) |
| 5 | 9.35 | 11.65 | 5.37 | 2.41 |

The τ = 16.5 ps arm lands *magnitude only*: s_eff=20 reaches 5.93 sheds →
n_end 15.07 but with t_first ≈ 13.2 ps (timing killed by t× ∝ τ); s_eff=5
strips to n_end 7.38. `frac_ions_shed = 1` at every set-s_eff point.

**Findings:**
1. **The kinetic diagnosis is quantitatively confirmed.** Dropping s from 60
   to 8 moves Δn̄ from 0.67 to 7.41 with *nothing else moving* — the freeze
   was entirely the RRK exponent (s−1 = 59) acting on x ≈ 0.032, exactly as
   the 45-run-probe diagnosis stated.
2. **s↔τ separability worked as designed** (Addendum A.2): first-shed timing
   stays gate-open-governed (≈ 5.4–6.0 ps at τ=6.55; ≈ 12.9–13.5 ps at
   τ=16.5, tracking t×), magnitude is s-governed. The landing region is
   **τ near the GAH25 mid-band × s_eff ∈ ≈ [8, 12]** — the anchored
   staircase alone nearly pins both, the best identifiability result of the
   program so far.
3. **The flexibility question is answered.** Relaxed terminal n now spans
   [7.3, 20.9] across the sweep (monotone in s and τ; spreads 1.3–1.5 He —
   a real distribution, not a delta). The mechanism was never inexpressive;
   the classical dof convention froze it. NB the report headline's
   per-picture ranges ([7.3, 20.9] mixture vs [19.3, 20.9] others) are a
   **design artifact** — only `statistical_mixture` received s_eff overrides
   (the mini-probe pins), not a picture-dependence claim.
4. **Constant s suffices at mid-band τ** (no in-window overshoot past 14 at
   s_eff=8): the outcome-(c) contingency (n-dependent scaled s = α·(3n−3),
   here α ≈ 8/60 ≈ 0.13) is *not demanded* by this data — keep the scaled
   shape as a documented alternative arm, not a build.

**Interpretation (physics).** `s = 3n−3` treats the He₂₁ shell as 60
fully-coupled classical oscillators; the landing s_eff ≈ 8–12 says the
effective heat bath is ~an order of magnitude smaller — the *expected*
direction for a cold, floppy, quantum He shell (most modes quantum-frozen or
uncoupled on the sub-ps shed timescale; cf. evaporative-ensemble treatments
of quantum clusters). The RRK form, energy gate, ladder, and Σ(21) crossing
budget all survive — outcome (b)'s "gate/ladder structure suspect" arm did
**not** fire. The pre-registered mechanism-level OQ resolved to exactly the
lever it named.

**Boundaries (unchanged):** N=50 single seed reads means, not distribution
tails; TDDFT is not ground truth (the 9 Å non-radial flag stands); the
landing does **not** discharge the F5 gate — it points the N=500 Stage-1/2
campaign, which resolves the gate.

**Decision (user, 2026-07-06):** document the finding now (this entry + the
Addendum status flip + a lever-3 resolution NB in
`TIER2_PHASE_D_BRIDGE_FINDINGS.md`); **before any CALIBRATION_MAP change**,
run a picture cross-check (`x2_only` + `cooling_relaxed` with reduced s_eff
at the landing arm) to confirm empirically that the 45-run probe's
picture-near-degeneracy prediction holds under reduced s; promotion
(dof-convention enum arm + Bounded s_eff knob + CALIBRATION_MAP
reclassification + F2 re-scope to an s_eff×τ co-fit) is **rediscussed after**
— no auto-promotion (the Addendum A.4 stance).

### Picture cross-check EXECUTED — near-degeneracy under reduced s_eff (2026-07-06)

10 additional probe runs: picture ∈ {`x2_only`, `cooling_relaxed`} ×
s_eff ∈ {30, 20, 12, 8, 5} at the landing arm (κ=1, τ=6.55 ps, probe pins;
the per-n controls for both pictures already exist from the 45-run grid).
Driven through the delivered `build_probe`/`_run_one` pipeline with an
overridden module-level grid (scratchpad driver) — **no repo-code change**;
the generator's active USER SETTINGS remain the test-locked 12-point
mini-probe. Run dirs land in the `tier2probe` namespace and are scored by the
delivered report script alongside the rest.

**Result — magnitude near-picture-invariance CONFIRMED.** At κ=1, τ=6.55,
Δn̄ per picture (mixture / cooling_relaxed / x2_only):

| s_eff | Δn̄ (mix / cool / x2) | n_ion_end_mean (mix / cool / x2) |
|---|---|---|
| 30 | 2.08 / 2.03 / 2.02 | 18.92 / 18.97 / 18.98 |
| 20 | 3.51 / 3.34 / 3.22 | 17.49 / 17.66 / 17.78 |
| 12 | 5.61 / 5.59 / 5.33 | 15.39 / 15.41 / 15.67 |
| 8 | 7.41 / 7.46 / 7.13 | 13.59 / 13.54 / 13.87 |
| 5 | 9.35 / 9.11 / 9.10 | 11.65 / 11.89 / 11.90 |

Shed-magnitude spread across pictures ≤ ~0.3 sheds (≤ 4 %) at every s_eff —
the 45-run probe's prediction (x = D₀(21)/Σ(21) nearly picture-invariant)
holds under reduced s; **the landing band s_eff ∈ ≈ [8, 12] is
picture-robust**. The one picture-sensitive read is *timing*: t_first at
s_eff=8 is 5.42 / 4.22 / 3.07 ps (mixture / cool / x2; anchored t★ = 5 ps) —
Σ(21) differs per picture, shifting t× = τ·ln(f_int·E/Σ(21)); mixture sits on
the anchor, x2_only opens ~2 ps early. This is weak discrimination, not a
picture selection: t× is degenerate with the Bounded `f_int` (a modest f_int
shift re-aligns any picture's gate-open), so the picture knob stays a
near-flat co-fit dimension for F4, entangled with f_int only through timing.
Trajectory-MAD argmin over all 65 scored runs remains (mixture, κ=1, τ=6.55,
s_eff=8) at 1.00 He; cooling_relaxed and x2_only both reach MAD ≈ 1.36 at
s_eff=8. Ledger residual stays uniform 2.23·10⁻⁵ eV.

**Adjudication:** deferred to the promotion rediscussion (user-level).

### s_eff PROMOTED to a Bounded calibration knob (user decision, 2026-07-06)

Following the mini-probe landing (outcome (a)) and the picture cross-check
(magnitude near-invariance confirmed), the user **confirmed the promotion
direction**. This is the CALIBRATION_MAP reclassification the Addendum A.4
stance reserved for a user-level mechanism decision — executed as
documentation only; **no code changed**.

**Documentation delivered (this entry's companion edits):**
- `CALIBRATION_MAP.md` — new 2026-07-06 update block; **row 10 reclassified
  Derived → Bounded** (constant s_eff, band ≈ [5, 20], staircase landing
  [8, 12] at mid-band τ; classical `s = 3n−3` demoted to the classical-limit
  arm of the dof-convention selection; s ≥ 1 guard + n=1 direct dissociation
  unchanged); rows 19/20 annotated (κ inverted + normalisation-capped for
  21→14; picture magnitude-degenerate, timing-only via Σ(21)→t×, degenerate
  there with f_int); tally Derived 6→5 / Bounded +1; anchor-coverage note —
  the 9 Å staircase (TDDFT prior) co-anchors s_eff×τ via first-shed timing,
  the first load shed from the one-observable Tier-2 stack.
- `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` — §4 dof-bullet NB +
  **A11 resolution NB** (the classical band [3n−3, 3n] survives as the
  rigid-classical limit, not the sweep band; effective reservoir ~8–12 modes
  = quantum mode-freezing / weak coupling, the RRKM-direction difference s
  was declared to absorb; the s↔κ joint coupling resolved *weak*); §9
  register row updated. Mechanism, invariant, and guards unchanged — this
  moves a parameter's *classification*, not the locked mechanism.
- `TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` — §1 row D mechanism-level
  resolution note (the F5 gate still resolves only via the 0.80 eV
  campaign); **F2 Stage-1 re-scope NB**: the κ×picture co-fit is superseded
  by an **s_eff×τ co-fit** (κ pinned with a sensitivity spot-check;
  picture + f_int ride as a timing-degenerate pair, reported not fit);
  §0 postponed-items row annotated.
- `TIER2_STAIRCASE_PROBE_PLAN.md` Addendum status + the
  `TIER2_PHASE_D_BRIDGE_FINDINGS.md` lever-3 NB flipped from
  "pending user decision" to confirmed.

**Deliberately NOT done (stays behind `[PROCEED TO IMPLEMENTATION]`):**
the dof-convention selection surface (NB: `cfg.evap_rrk_dof` already provides
the constant-s override the campaign needs — whether a named enum is added or
the existing Optional field *is* the selection surface is an implementation
decision for that slice); the `S_EFF` grid dimension + tag encoding in
`gen_tier2_runs.py` / `tier2_run_tag` (the probe-tag parity-lock test names
the exact encoding); any F3/F4 scoreboard/report column additions. Rule-2:
no declared-but-unread fields were added by this decision.

**Next:** implement the F2 campaign re-scope (on trigger), then run the
re-scoped N=500 Stage-1/2 validation campaign at 0.80 eV (9 Å + 18 Å)
concentrated on the landing region — that campaign, not the probe, resolves
the F5 gate. *(Correction, same day: the campaign is 9 Å-only per the standing
2026-07-03 scope decision — no 18 Å reference exists; restated by the user
with the implementation trigger.)*

### F2 Stage-1 re-scope IMPLEMENTED — s_eff×τ co-fit campaign (2026-07-06)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger (user:
"PROCEED TO IMPLEMENTATION (just stage 1 … only for 9 Å, no 18 Å
reference)"), scoped to Stage 1 at 9 Å / 0.80 eV — consistent with the
standing 2026-07-03 9 Å-only campaign decision; the old Stage 2
(τ-insensitivity) is absorbed by the in-grid τ dimension and read
report-side by F4. TDD throughout (11 RED failures + 1 collection
`ImportError` watched across five test modules before any source change).

**Selection-surface decision (flagged in the promotion entry, now resolved):**
no new enum — the existing `cfg.evap_rrk_dof` Optional field **is** the
dof-convention selection surface (`None` = per-n `s = 3n−3` classical arm;
a value = the promoted constant s_eff). No new `SimConfig` field, no schema
bump, no driver/physics change (the field was already physics-live).

**Build:**
- `scripts/tier2_common.py` — `tier2_run_tag` / `tier2_run_dir_name` gain a
  None-sentinel `evap_rrk_dof` kwarg appending `_sNN.NN` (two decimals,
  before any `_totalstrip` suffix); `None` appends nothing — pre-promotion
  campaign tags byte-identical. The probe↔campaign **parity-lock test now
  covers the s_eff suffix** in both encoders.
- `scripts/gen_tier2_runs.py` — Stage-1 grid re-scoped: κ=1 / mixture pinned;
  `TAU_GRID_PS = [2.6, 6.55, 16.5]` × `S_EFF_GRID = [None, 5, 8, 12, 16, 20]`
  → **18 runs** (per-n None = in-grid classical-limit control; band [5, 20]
  brackets the [8, 12] landing). `F_INT = 0.5` (the landing-prior pin —
  the s_eff ≈ 8 landing was established there; `None` → floor still
  supported). `campaign_grid_points()` now returns 5-tuples
  `(picture, κ, τ, s_eff, f_int)`. The superseded κ×picture Stage-1 grid is
  preserved in an "ORIGINAL STAGE 1" comment block (one-swap restore).
- `scripts/post_processing/tier2_size_distribution_table.py` (F3) — `s_eff`
  column (from `cfg.json`, `None` → per-n), inserted after `tau_ps` —
  without it the 6 s_eff values at a fixed (κ, τ) would be scoreboard-
  indistinguishable (the probe-report review's bug class).
- `scripts/post_processing/tier2_identifiability_report.py` (F4) — the
  κ- and τ-cell keys gain `s_eff` (via `row.get`, so pre-promotion rows
  keep working) — two runs differing only in s_eff can no longer pool into
  one landscape/sweep cell; new `build_s_eff_landscapes` +
  `SEffLandscape` (per-n control rows excluded from the numeric landscape);
  the generic landscape assessor gained an `axis` parameter (default
  `"kappa"`, back-compatible); the **`s` checklist entry is now computed**
  (W1-vs-s_eff landscape, `axis="s_eff"`) with a `held_fixed` static
  fallback when no sweep is on disk (text updated to the promoted-Bounded
  wording); module docstring bullets updated.

**Tests (+10 net; full suite 2097 → 2107 passed, 0 failed):** campaign-tag
s_eff suffix/omit + dir non-collision (`test_tier2_common`), parity-lock
s_eff extension (`test_tier2_staircase_probe`), re-scoped 18-point grid +
landing-prior f_int pin + floor-under-monkeypatch + 5-tuple schedule with
s_eff/τ stamps (`test_gen_tier2_runs`), F3 s_eff column presence + per-n
None read + set-s_eff=8 tiny-run read, F4 s_eff cell splitting (τ sweeps and
κ landscapes not pooled across s_eff), per-n-control exclusion, and the
assessed `s` entry (identifiable interior / not_bracketed edge; the
held-fixed fallback stays covered by the existing static test).

**Rule-2 / scope:** no new `SimConfig` field, no schema/RNG/drag-law/driver/
propagation change; probe scripts untouched (their delivered tags remain
byte-identical); the F5 production switch, total-strip variant, and 18 Å leg
remain out of scope. The 18-run campaign itself is an **operator action**
(`python scripts/gen_tier2_runs.py`, then the F3/F4 report scripts) — not
executed as part of this delivery.

---

## Cooling spatial gate — `cooling_spatial_gate` model arm + total-strip A/B probe DELIVERED (2026-07-06)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger, after a design
discussion (user question: before the production switch, *can our parameters
produce total stripping* — near-bare I⁺, the experimental 43 % bare-I⁺ peak?).
Design doc: `C:\Users\user\.claude\plans\1-yes-it-deserves-nested-allen.md`.
TDD throughout (RED watched on every new physics/config surface before the fix).

**The finding that drove it (structural trace of the locked mechanism).** No
fixed-asymptote wall exists (`E_inf(N)→0` as N→0 was built for OQ6 reachability).
The limiter is **K2 Newton cooling**: it drains `E_int` toward 0, opening the
self-bound gate but then quenching shedding once `E_int < D_0(n)` — terminal n
floors at n~few (shell-retaining), reached in a few τ (the 8.5 µs flight is
irrelevant; cooling zeroes `E_int` long before). The user then flagged the real
asymmetry: **drag (`γ∝ρ_He`) and pickup (`λ∝ρ_He`) both gate off outside the
bubble via the shared erf-complement surface, but K2 cooling is applied ungated
everywhere** — including after the complex is ejected into vacuum, where there is
no droplet bath to radiate into (the GAH25 τ was fit for a near-droplet *growing*
Na⁺ shell; §6 K2 growing-vs-shrinking caveat). Gating K2 off on ejection lets an
already-self-bound complex keep shedding (shed-invariant margin `G=E_int−Σ(n)`) →
deeper stripping → the total-vaporization end.

**Design decisions (user, 2026-07-06, before coding):**
1. **First-class interchangeable `SimConfig` enum arm**, not a hard-wire (the
   drag-port working method — every model choice behind an enum, arbitrated by the
   size distribution).
2. **Clean erf-complement** (cooling → 0 outside the bubble), reusing the drag gate
   surface; a residual out-of-bubble floor `rho_min` (dragged He cloud) is a
   deferred OQ, not built.
3. **Reuse the drag/pickup gate boundary** `drag_gate_steepness(cfg)` (14.2 Å) — no
   new steepness knob; one bubble boundary for all three He-mediated channels.
4. **Probe-scoped for now** — wired into the staircase-probe path only; the campaign
   encoder, generator, and F3/F4 reports stay untouched until a probe result
   justifies promotion.

**Physics contract.** `none` (default): `Ė_int|K2 = −E_int/τ` (byte-identical to
the pre-arm code). `density_scaled`: `Ė_int|K2 = −(ρ_He/ρ_bulk)·E_int/τ`, i.e.
`τ_eff = τ/rho_ratio`, closed form `E_new = E_inf + (E−E_inf)·exp(−dt·rho_ratio/τ)`.
`[ρ_He/ρ_bulk]=1` → eV/ps unchanged.

**Build (additive; `none` byte-identical everywhere):**
- `physics/solvation_cooling.py` — `newton_cool_step` gains keyword-only
  `rho_ratio=1.0` (`decay = exp(−dt·rho_ratio/τ)`), a `rho_ratio<0` fail-loud guard
  (the anti-cooling class the dt guard already names), and an array-`rho_ratio` →
  ndarray return. `rho_ratio=1.0` is byte-identical (`dt·1.0` exact).
- `config.py` — `CoolingSpatialGate = Literal["none","density_scaled"]`, field
  `cooling_spatial_gate: CoolingSpatialGate = "none"` (inert default), the
  `_KNOWN_COOLING_SPATIAL_GATES` tuple + `check_cooling_spatial_gate_config` reject
  guard wired into `validate()` (mirrors `check_helium_density_config`). Physics-live
  (read by `biphasic_step`), not just guard-live. Old cfg.json without the field
  loads fine (missing key → default).
- `simulation/ion_propagation_step.py` — `biphasic_step` **hoists** the RNG-free
  `_depth`/`rho_he_ratio` gate computation above the K2 block (byte-identical draw
  stream) and threads `cool_rho = rho_ratio if cfg.cooling_spatial_gate ==
  "density_scaled" else 1.0` into `newton_cool_step`. `E_dissip` booking unchanged
  (books the actual before/after drain → the 5-term invariant stays closed).
- `simulation/relaxation_stage.py` — **no code change**; the verbatim `biphasic_step`
  reuse carries the arm through (positions advanced → `depth` live). Behavioral
  change is intended: a gated-off ejected complex sheds without the cooling quench;
  `E_int` stays monotone non-increasing so freeze-termination holds. **[Corrected by
  the 2026-07-06 code-review fixes below — monotone non-increasing does NOT imply
  freeze-termination under `density_scaled`: an ejected fragment stops cooling and can
  sit above `D_0(n)` indefinitely. The probe's finite 1000 ps cap is what bounds it; a
  config guard now rejects a runaway cap for that arm.]**
- `scripts/tier2_common.py` — `build_biphasic_cfg` gains a None-sentinel
  `cooling_spatial_gate` kwarg (injected only when set); `tier2_probe_run_tag` /
  `tier2_probe_run_dir_name` append `_cgds` for `"density_scaled"` (nothing for
  `none`/`None` → delivered probe tags byte-identical). **`tier2_run_tag` (campaign
  encoder) NOT modified** — the parity-lock test exercises only base knobs +
  `evap_rrk_dof`, so it stays green while the probe carries this extra axis.
- `scripts/gen_tier2_staircase_probe.py` — active USER SETTINGS flipped to the
  **total-strip capability A/B grid**: `COOLING_GATE_GRID = ["none","density_scaled"]`
  × `S_EFF_GRID = [1,2,3,5]` (s_eff=1 = max-kinetics bound) × `TAU_GRID_PS =
  [6.55,16.5,30]` (κ=1/mixture pinned) → **24 runs**, 9 Å / 0.80 eV / N=50.
  `RELAXATION_TIME_PS = 1000.0` (finite cap required: a gated-off self-unbound
  complex never freezes, so the 8.53e6 ps default would not terminate). The 12-run
  mini-probe and 45-run full probe preserved as restore comments. `probe_grid_points`
  now yields 5-tuples.
- `scripts/post_processing/tier2_staircase_probe_report.py` — `PROBE_TABLE_COLUMNS`
  gains `cooling_gate`, `n_relaxed_min` (deepest reachable shell = smallest occupied
  n in the relaxed distribution — the total-strip headline), and `frac_frozen` (E2
  freeze completeness, reusing `relaxation_stage._freeze_mask` on the final relaxed
  column so a non-terminating gated cascade is surfaced, not silently read as a
  converged terminal n).

**Tests (+~20; full suite 2126 → 2130 passed, 0 failed):**
`test_solvation_cooling.py` (rho_ratio=1 back-compat / rho_ratio=0 no-op / scaled
decay / negative-rho raise / array broadcast); new
`test_cooling_spatial_gate_config.py` (default-inert, both arms valid, unknown
rejected by guard + validate, Literal↔`_KNOWN` parity); `test_biphasic_step.py`
(`TestCoolingSpatialGate`: deep-inside == ungated exactly, far-outside cooling off /
E_dissip zero, intermediate-depth drain strictly between); `test_relaxation_stage.py`
(E_dissip==cumulative-K2-drain oracle gate-invariant under `density_scaled`; first-step
gated drain < ungated); `test_tier2_staircase_probe*.py` (total-strip 24-point grid,
back-compat 45-run with gate="none", `_cgds` suffix append/omit + namespace lock,
A/B distinct dirs, cfg stamp, report `cooling_gate`/`n_relaxed_min`/`frac_frozen`
columns). Existing relaxation E_dissip oracle stays green unchanged (magnitude-agnostic
reconstruction; default `none`).

**Rule-2 / scope:** one new `SimConfig` field (`cooling_spatial_gate`, physics-live
— not a rule-2 carry); no schema bump, no RNG draw-order change, no drag-law / neutral
/ ion-propagation change; the relaxation stage is untouched (arm flows through). The
`rho_min` floor, campaign wiring (`tier2_run_tag`/`gen_tier2_runs.py`/F3/F4 columns),
and the F5 production switch remain out of scope. Docs updated: MASS §6 K2 + §6.11,
`CALIBRATION_MAP.md` (update block + row 5a), `DRAG_PORT_DESIGN_DECISIONS.md` §5.8.

**Verification.** Full suite **2130 passed, 0 failed**, 80 warnings = the documented
§6.5 mass-pairing path. **Next (operator action):** run the 24-point total-strip A/B
probe (`python scripts/gen_tier2_staircase_probe.py`, then
`tier2_staircase_probe_report.py`), confirm the `none` arm reproduces the
shell-retaining floor and read whether `density_scaled` drives `n_relaxed_min` toward
bare across the s_eff/τ grid (with `frac_frozen` confirming completeness) — the
capability check that gates whether/how the production campaign can express the
experimental 43 % bare-I⁺ peak. Reported, not auto-adjudicated.

---

## Cooling spatial gate — code-review fixes (high-effort review, 2026-07-06)

A high-effort workflow code review (4 finders + 10 verifiers) of the cooling-gate
delivery confirmed the `none`-arm core (byte-identity, frozen two-draw stream after
the RNG-free depth/`rho_ratio` hoist, 5-term closure, `newton_cool_step` guards,
`Literal↔_KNOWN` parity) and surfaced 7 findings, all fixed under the user's explicit
"fix all findings even 6" trigger. TDD for the behavioural fixes.

1. **`density_scaled` breaks the guaranteed relaxation freeze-termination
   (correctness).** The E2 early-exit relies on K2 *actively* draining `E_int` to 0;
   under `density_scaled` an ejected fragment (`rho_He→0`) stops cooling, `E_int` is
   merely held, and evaporation's RRK rate vanishes near threshold, so `E_int` can
   sit above `D_0(n)` forever → the loop runs the full cap (the 8.53e6 ps flight cap
   ⇒ ~8.5e8 steps ⇒ effective hang). The probe already mitigates via `RELAXATION_TIME
   _PS=1000.0`; **fix** adds a `check_relaxation_config` step-budget guard
   (`relaxation_time_ps/dt_relax ≤ _DENSITY_SCALED_MAX_RELAX_STEPS = 1e7`) that turns a
   runaway cap into a load-time error, and corrects the false "freeze-termination
   holds" claim (delivery entry above + the `relaxation_stage` Termination docstring).
2. **`_cfg_matches_row` omitted `cooling_gate` (correctness).**
   `tier2_staircase_probe_report.py` — best-case overlay re-find could load the wrong
   A/B arm; **fix** adds the `cooling_gate` match (`row.get(..., "none")` back-compat).
3. **Report figures collided the A/B arms (correctness).** Small-multiples (curve
   label/linestyle) and the relaxed-n heatmap (subplot-row facet) now thread
   `cooling_gate`; single-gate layouts stay byte-identical, so delivered non-A/B
   figures are unchanged (regenerate the A/B PNGs by re-running the report — operator
   action, figures are not in pytest).
4. **Campaign encoder could collide gated run dirs (latent).** `tier2_run_tag` /
   `tier2_run_dir_name` now append `_cgds` for `density_scaled` (after `_s`, before
   `_totalstrip`); the probe↔campaign parity-lock test is extended to cover it.
   Campaigns not sweeping the gate stay byte-identical.
5. **`newton_cool_step` (cleanup):** `rho = np.asarray(rho_ratio)` bound once (was
   coerced twice); behaviour-neutral.
6. **Enum-reject guards consolidated (cleanup, user-approved).** The 7 near-verbatim
   `if x not in _KNOWN_*: raise ValueError("unknown …")` blocks (drag_form,
   dissociation_ladder, helium_density_profile, cooling_spatial_gate, the 3 pickup
   selectors) collapse into one `config._reject_unknown_enum(value, known, field=…)`
   helper; messages are byte-identical, so the field-name `match=` guard tests still
   pass.

Two candidates were **correctly refuted** by the verify pass (kept as-is): the
`relaxation_stage._freeze_mask` cross-module import (a hard top-level import — fails
loud on rename, not silent drift) and the `== "density_scaled"` suffix special-case
(closed enum).

**Tests (+~10):** F1 density_scaled runaway-cap reject + modest-cap/ungated accept
(`test_relaxation_stage`); F2 matcher gate-distinguish + back-compat
(`test_tier2_staircase_probe_report`); F4 campaign `_cgds` append/omit/order +
parity-lock extension (`test_tier2_staircase_probe`). F5/F6 covered by existing
behaviour-preserving tests. **Rule-2 / scope:** no new `SimConfig` field, no schema/
RNG/drag-law/driver/propagation change; the F6 consolidation is the only pre-existing-
code touch and is message-preserving.

---

## Total-strip A/B probe EXECUTED + probe findings consolidated (2026-07-07)

The 24-point cooling-gate A/B grid (gate ∈ {`none`, `density_scaled`} ×
s_eff ∈ {1, 2, 3, 5} × τ ∈ {6.55, 16.5, 30} ps, κ=1/mixture/probe pins,
relaxation cap 1000 ps) was executed by the operator and scored on 2026-07-07
with the delivered report script (full re-score of all **87** probe dirs;
wiring oracle still exact — the per-n bridge control reproduces Δn̄ = 0.67 /
n_end = 20.33 under the post-gate code, confirming the `none`-arm
byte-identity end-to-end; ledger residual uniform 2.23·10⁻⁵ eV).

**Headline results (full tables + trajectory traces in the findings doc):**
- **Gated + mid-band τ = 6.55: near-total strip.** `density_scaled` strips
  the whole ensemble to n̄ ≈ 2.7–3.5 (relaxed min n = 2, spreads 0.35–0.50)
  vs the ungated n̄ ≈ 8.5–11.9 — the post-ejection K2 quench was the
  shell-retention limiter, as the structural trace predicted.
- **Gated + τ ≥ 16.5: zero sheds, permanently closed gate.** Cooling shuts
  off at droplet ejection (~6 ps) with E_int held at 0.280 eV > Σ(21) ≈
  0.188 eV — the self-bound crossing never happens (frac_frozen = 0.0; the
  finite cap, not freeze, terminates — the code-review correction in
  action). The gated arm is **all-or-nothing in the t×↔ejection race**.
- **Ungated floor confirmed:** even at the max-kinetics bound s_eff = 1,
  `none` floors at n̄ ≈ 8.5 (τ=6.55) / 4.7 (τ=16.5) / relaxed 3.1 (τ=30,
  cascade continuing into the relaxation stage), with first sheds 13–23 ps
  at long τ — deep stripping and staircase timing irreconcilable ungated.
- **Bare I⁺ (n = 0) is NOT reached anywhere** — floor n = 2, energetic
  under the Σ(21)-crossing budget: expressing the experimental 43 % bare
  peak remains open (F5/production discussion item, not a retune).

**Findings doc created:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` — consolidates
all four probe waves (45-run capability probe, s_eff mini-probe, picture
cross-check, total-strip A/B) with every table, the mechanism traces, the
insight register I1–I9, boundaries, open questions (rho_min softening,
bare-peak channel, n-dependent s, standing conditional triggers), and the
87-dir run inventory. Reported, not auto-adjudicated; nothing in it
discharges the F5 gate (that stays with the 0.80 eV N=500 campaign).
Documentation-only entry — no code changed.

### Post-Wave-4 decisions + Wave-5 plan recorded (user, 2026-07-07)

Discussion outcome, recorded as `TIER2_STAIRCASE_PROBE_PLAN.md` **Addendum
B** (documentation only; no code changed):

1. **`density_scaled` is the intended primary production arm** (physical
   self-consistency: one ρ_He bubble boundary for drag/pickup/cooling;
   capability: the only arm that approaches the 43 % bare-I⁺ peak). Stays
   an interchangeable enum with a `none` sensitivity leg; formal
   campaign/production promotion waits for Wave 5.
2. **f_int promoted to a genuine sweeping parameter** — accepted
   consequence of the gated t×↔ejection cliff. Derived F5 prediction
   recorded: at 2.70 eV / f_int = 0.5 the gated arm sheds nothing
   (t× ≈ 12.9 ps > t_eject ≈ 6 ps); production needs roughly
   f_int ≲ 0.17 — the scenario-keyed floor is load-bearing.
3. **Clean erf-complement adopted as correct** (zero out-of-bubble
   cooling); a residual `rho_min` floor judged not physical — kept as a
   documented open question for domain experts, explicitly **no
   implementation**. The τ/f_int cliff is accepted as a real model
   feature.

**Wave 5 planned (Addendum B.2–B.4):** gated landing re-location
mini-probe — s_eff ∈ {8, 12, 16, 20, 30} × τ = 6.55 ps ×
`density_scaled`, probe pins, 5 runs (Wave-4 gated s_eff ≤ 5 + Wave-2
ungated sweep bracket it from disk). Rationale: the [8, 12] landing prior
is gate-conditional (established ungated; the gate ~doubles in-window
shedding at s_eff = 5), so running the F2 campaign on the ungated prior
under a gated arm would repeat the mistake Wave 1 existed to prevent. No
new code/config surface needed; execution via the Wave-3 scratchpad-driver
precedent or a USER-SETTINGS flip behind `[PROCEED TO IMPLEMENTATION]`.
The F2-plan update (gate × s_eff × τ × f_int Stage-1 shape, B.5) is
deferred until the Wave-5 result exists.

### Wave-5 gated landing re-location EXECUTED (2026-07-07)

Under the `[PROCEED TO IMPLEMENTATION]` trigger (user: "PROCEED TO
IMPLEMENTATION of the 5-wave and lets discuss results"). Executed via the
Wave-3 **scratchpad-driver route** (Addendum B.4): the delivered
`gen_tier2_staircase_probe.py` `build_probe`/`_run_one` pipeline driven with
an overridden module-level grid — **zero repo-code change**; the generator's
test-locked active USER SETTINGS (the 24-point A/B grid) are untouched. 5 new
gated dirs (`s_eff ∈ {8,12,16,20,30}` at τ=6.55, `density_scaled`) in the
`tier2probe` namespace, scored by the delivered report script (full re-score,
now **92 probe dirs**). Wiring: the gated s_eff=5 point re-scored unchanged
from Wave 4 (Δn̄ = 17.52); ledger residual uniform 2.23·10⁻⁵ eV.

**Two results (full tables + traces: `TIER2_STAIRCASE_PROBE_FINDINGS.md`
§4b, insights I10–I12):**

1. **Landing shifts to s_eff ≈ 30** (gated Δn̄ 7.20, n_ion_end 13.80,
   MAD 0.96 — grid best), vs the ungated s_eff ≈ 8: the gate multiplies the
   effective in-window shedding by ≈ 3.75× by removing the post-ejection
   quench. Prediction confirmed; the ungated [8,12] prior does **not** carry
   to a gated arm (gated s_eff=8 → n_ion_end 4.84, MAD 4.70). Timing lands
   ~2.3 ps late (t_first 7.28 vs t★ 5), f_int-realignable.
2. **New structural finding — the gated in-window staircase and the terminal
   read decouple.** Ungated: frac_frozen = 1 everywhere, ion-end n = terminal
   n (one converged read). Gated at s_eff=30: frac_frozen = 0, the relaxation
   ran the full 1000 ps cap without any ion freezing (n 13.80 → 6.27, E_int
   0.104 → 0.034 eV, slope −0.0012 n/ps, decelerating). s_eff sets the
   post-ejection cascade *rate*, not the endpoint (energetic floor n ≈ 2 for
   all s_eff). Confirmed by the stored relaxation lengths: gated s_eff=1 =
   2 timesteps (froze at 30 ps), gated s_eff=30 = full 50001. Consequence:
   the gated *terminal* size distribution (the Tier-2 arbitration observable)
   is **flight-time dependent** at the staircase-landing s_eff — the physical
   ~8.5 µs flight is 8500× the tractable cap.

**Interpretation carried forward (I12):** ungated and gated are
terminal-regime opposites (shell-retaining/converged vs
deep-stripping/rate-limited); neither single knob point yields the broad
bimodal experimental distribution — its shape must come from ensemble
heterogeneity across the t×↔ejection race, which the gate makes two-sided.
Documentation-only entry; the campaign-shape decision (B.5) stays open for the
discussion. No code changed; Addendum B status flipped to EXECUTED.

## Detection-time continuation stage — design DELIVERED (2026-07-07)

Post-Wave-5 discussion thread 1 (the I11 flight-time-dependence blocker)
resolved into a full physics-definition design:
`TIER2_DETECTION_STAGE_DESIGN.md`. **Documentation only — no code; the build
stays behind `[PROCEED TO IMPLEMENTATION]`.**

**Problem:** the gated terminal read (the Tier-2 arbitration observable) is
cap-artifactual — fixed-dt integration cannot reach the experimental flight
time structurally (`ν·dt ≤ 0.1` with ν = 2.42/ps ⇒ dt ≤ 0.041 ps ⇒ ≥ 2·10⁸
steps to ~8.5 µs).

**Design (four forks adjudicated by the user, 2026-07-07):**

1. **t_detect = Sourced TOF flight time** (detector-arrival read; apparatus
   provenance to be recorded as a CALIBRATION_MAP row — open input §3.4).
2. **Stochastic event-driven realization** (exponential waiting times;
   integer per-ion n, E1-compatible; the hypoexponential closed form is the
   pytest oracle, not a production arm).
3. **Additive stage after E2** (`simulation/detection_stage.py`; the
   reviewed relaxation stage is untouched; exactness guard at handover —
   the per-ion drain bound `(ρ/ρ_bulk)·(t_detect−t_h)/τ ≤ ε_drain`).
4. **Full per-event ledger** (K1 drain / e_bind fold / cold-shed defect per
   event; 5-term closure extends across the new stage; stored shed times
   make any intermediate-time read a free report-side reconstruction).

**Key exactness fact (grounded in the delivered code):** under the handover
preconditions (λ₀ = 0, γ = 0, gated ρ→0) E_int is constant between sheds, so
the locked mechanism reduces *exactly* to a Markov jump chain — the stage is
a change of integration scheme (≤ ~21 events/ion to any t_detect), not new
physics; the fire path composes the delivered `rrk_rate` / `cold_shed` /
`dE_int_shed_eV` / `e_bind_pair_eV` primitives verbatim (rule-1 reuse).
Permanent-state taxonomy: `frozen` / `suppressed` (delivered gate semantics
faithfully extended — OQ-B untouched) / `time_exhausted` (the live-cascade
detector snapshot Wave 5 could not compute). Translation is exactly out of
scope (cold shed leaves v unchanged; the E2 decoupling fact). Artifact:
`detection.npz` (new small file; **schema v7 untouched**). New OQ-E recorded
(no radiative/electronic channels over the µs flight — domain-expert item).

**Consequences recorded:** the B.5/F2 campaign re-scope is unblocked (scorer
reads `detected`; `relaxed` + `frac_frozen` demoted to convergence
diagnostics); the E2 cap becomes a campaign-economy config decision. Open
before the trigger: the `detection_time_ps` value + provenance (user), and
the stage's slice placement in the Phase-F program.

**NB (2026-07-07, same day):** the Sourced input is resolved — the user
confirmed **t_detect = 8.53 µs** as accurate, from the experimental-setup
publication. Recorded as **CALIBRATION_MAP row 24** (Sourced tally 8→9);
design doc §1/§3.4/§7 updated. The only remaining pre-trigger item is slice
placement in the Phase-F program.

### Waves 6+7 PLANNED — f_int probe + detected-read re-score (user, 2026-07-07)

Probe-program continuation agreed in discussion (explicitly **not**
campaign/production — the B.5 campaign shape stays parked; the goal remains
model understanding: parameter sensibility + flexibility reach). Plan
recorded as `TIER2_STAIRCASE_PROBE_PLAN.md` **Addendum C**; execution awaits
`[PROCEED TO IMPLEMENTATION]`. Documentation-only entry — no code changed.

- **Wave 6 — f_int probe (10 new runs, gated arm):** the one knob every
  wave pinned at 0.5, now swept with four pre-registered predictions:
  P1 gated staircase timing re-aligns at f_int ≈ 0.35 (t★ = 5 ps);
  P2 f_int = 0.65 / τ = 6.55 locates the closed side of the t×↔ejection
  cliff (sheds nothing); P3 f_int = 0.25 re-opens the dead τ = 16.5 arm
  (deep strip); P4 floor-n(f_int) maps the in-bubble cooling leak over
  [t×, t_eject] — the in-model bound on total stripping (the quantified
  OQ-B gap). Grid: τ=6.55: s_eff ∈ {8,30} × f_int ∈ {0.24, 0.30, 0.35,
  0.65} + τ=16.5: s_eff ∈ {8,30} × f_int = 0.25; probe pins; scratchpad
  route; the two on-disk f_int = 0.50 points are the wiring oracle.
- **Wave 7 — detected-read re-score (zero new MD runs):** after Slice DS
  (the detection stage per `TIER2_DETECTION_STAGE_DESIGN.md`), re-read all
  ~102 probe dirs at t_detect = 8.53 µs → the flexibility map *at the
  detector*, resolving Wave 5's open end (gated s_eff = 30 at the detector:
  mid-shell or floor?). DS acceptance criterion recorded: existing probe
  dirs must stay loadable after the config-surface addition (back-compat
  test), else the stale-artifact policy fires and Wave 7 needs
  regeneration — surfaced, never silently migrated.

### Wave-6 f_int probe EXECUTED (2026-07-07)

Under the `[PROCEED TO IMPLEMENTATION]` trigger (user: "PROCEED TO
IMPLEMENTATION of Wave 6 and lets discuss results"). Wave-3/5
scratchpad-driver route through the delivered `build_probe`/`_run_one`
pipeline — **zero repo-code change**; the f_int dimension rides the existing
`_fiX.XX` tag (no namespace work needed). 10 new gated dirs → **102 probe
dirs**; full re-score with the delivered report script; wiring oracle exact
on all three gated f_int = 0.50 companions (s_eff = 5/8/30); ledger residual
uniform ≈ 2.23·10⁻⁵ eV.

**Verdicts on the pre-registered predictions (tables + detail:
`TIER2_STAIRCASE_PROBE_FINDINGS.md` §4d, insights I15–I17):**

1. **P1 qualitatively confirmed, quantitatively refined:** the gated first
   shed is monotone-f_int-tunable across 0.7–7.3 ps (spans t★ = 5), but the
   anchored timing sits at **f_int ≈ 0.40–0.42**, not the predicted 0.35 —
   the gated crossing lag *grows* with t× (≈1.3 ps at 0.35 → ≈2.3 ps at
   0.50); the constant-lag assumption is the localized miss.
2. **P2 confirmed exactly:** f_int = 0.65 / τ = 6.55 sheds nothing at either
   s_eff — the closed side of the t×↔ejection cliff located in-band
   (0.50 < f_int* < 0.65).
3. **P3 confirmed:** f_int = 0.25 re-opens the dead τ = 16.5 arm; s_eff = 30
   lands near-staircase magnitude there (n_end 14.65, MAD 1.46) — a second
   staircase-magnitude region; τ↔f_int trade off along the race.
4. **P4 confirmed:** the relaxed floor is monotone in f_int — earlier
   gate-open is *shallower* (in-bubble leak eats the Σ(21) budget); deepest
   sampled strip at f_int = 0.50 (n̄ ≈ 3.1); **bare I⁺ unreached at any
   f_int** — the OQ-B gap quantified from the f_int side.

**Structural finding (supersedes the "timing-only" note in the gated
regime):** under `density_scaled`, f_int is a **joint timing +
effective-budget knob** — it cannot re-align staircase timing without moving
magnitude; landing both needs a joint (s_eff, f_int) co-fit (interpolated
prior: s_eff ≈ 25, f_int ≈ 0.42 at τ = 6.55). Findings doc §0 wave table,
§4d, §5 (I15–I17), §8 inventory updated; Addendum C status flipped
(Wave 6 EXECUTED; Wave 7 stays planned behind Slice DS). Documentation +
scratchpad execution only — no repo code changed.

### Post-Wave-6 decisions (user, 2026-07-07)

Discussion outcome on the Wave-6 results, recorded as findings **I18** + a
B.5 NB in the probe plan (documentation only; no code changed):

1. **τ and f_int are connected — independent grids are rejected.** Gated
   observables organize along the race margin Δ× = t_eject − t× and "make
   no sense otherwise" (user). Any future campaign grid samples **Δ×
   directly** (f_int derived per τ via f_int = (Σ/E)·e^{t×/τ}); τ stays a
   separate dimension only for its independent role (the in-bubble
   leak/quench strength). Supersedes the B.5 "× τ × f_int" product shape.
2. **The bare-peak channel (OQ-B) is deferred to the F5/production
   discussion** — the probe program has exhausted the kinetic and timing
   levers (I8 + I17: the n ≈ 2 floor is energetic); no further probe work
   on it.
3. **Wave 7 (Slice DS + detected-read re-score) is confirmed as the next
   step, not yet triggered.** The campaign stays parked; probe-only scope
   holds.

### Slice-DS pre-trigger adjudications — E2 retained + re-scoped, DS skippable (user, 2026-07-07)

Before triggering Slice DS the user challenged the stage architecture
directly: *why keep the relaxation stage at all — is E2 dead code once DS
exists?* Resolved by an investigation of `relaxation_stage.py` + discussion.
**Documentation only — no code changed; the build stays behind
`[PROCEED TO IMPLEMENTATION]`.**

**Investigation conclusion (the regime map):** the three stages correspond
to which continuous channels are alive — ion stage (drag+pickup+cooling,
fixed-dt MD), E2 (cooling only, fixed-dt), DS (nothing continuous, exact
jump chain). DS's exactness requires E_int constant between sheds; K2
cooling breaks exactly that, so DS structurally cannot absorb any
live-cooling period. E2 is therefore **not dead code**: it is (a) the only
possible terminal solver for the ungated `none` arm (a
working-method-mandated sensitivity leg — ungated ions arrive at ion-end
neither frozen nor P3-compliant, failing the drain bound by ~12 orders of
magnitude), (b) the P3 guard's remedy path for any in-/near-bubble ion at
handover, and (c) the producer of the 102 existing `relaxation.npz`
artifacts Wave 7 re-reads. What *is* obsolete is E2's original R5 mission
("reach the terminal by matched-time integration" — unachievable on the
gated arm per Wave 5/I11); a time-inhomogeneous single-stage Gillespie
replacing E2+DS was considered and rejected (duplicates the K2 drain law in
hazard integrals; loses the exact-solver/analytic-oracle status; rewrites
the reviewed E2 contract).

**Decisions (user):** keep E2, **make it shorter and skippable**, with one
connected code review:

1. **E2 retained, contract re-framed** — ungated: terminal solver to
   freeze-out (DS no-ops); gated: a short handover **bridge** to a P3-clean
   state. The relaxed read is demoted to a convergence diagnostic; the
   module docstring's mission statement is updated in the slice.
2. **Gated caps shortened by config** (future-facing USER SETTINGS, ~10 ps
   class instead of 1000 ps; the §2.1 guard fails loud if shortened past
   ejection/decoupling). Existing probe dirs untouched.
3. **DS is skippable-E2 capable:** with `relaxation_stage_enabled=False` it
   seeds directly from `ion.npz` (both artifacts are v7 `IonCheckpoint`s —
   one seeding contract, `run_detection_stage(seed_ckpt, ...)`); the P1–P3
   guard is the sole defense on that path, and the `frac_frozen` diagnostic
   is replaced by the DS `state_reason` fractions. Amends design fork 3 to
   "additive after E2 *or* the ion stage".
4. **Slice placement resolved (closes design §7):** standalone pre-F5
   **Slice DS** — DS build + probe-report `n_detect_*` extension + the E2
   re-frame/skip path, one connected code review; the B.5/F2 plan update
   stays parked with the campaign.

Design doc updated in place (`TIER2_DETECTION_STAGE_DESIGN.md`: status
block, §1 item 5, §2.1, §3.1, §3.4, §5, §6, §7 — no open items remain).
Next: `[PROCEED TO IMPLEMENTATION]` for Slice DS, then the manual
production-wiring oracle (gated s_eff = 30 → t_detect) and Wave 7.

### Slice DS DELIVERED — detection-time continuation stage (2026-07-07)

Built under the `[PROCEED TO IMPLEMENTATION]` trigger, per
`TIER2_DETECTION_STAGE_DESIGN.md` + the §1-item-5 amendments (E2 re-frame,
skip path, connected review scope). **Full suite green: 2173 passed.**

**Delivered:**

1. **`i2_helium_md/simulation/detection_stage.py`** (new) —
   `run_detection_stage(seed_ckpt, cfg, *, rng, save_path)` →
   `DetectionResult`. Per-ion Gillespie event loop composing the delivered
   primitives verbatim (`rrk_rate` all-gating, `cold_shed`,
   `dE_int_shed_eV`, `e_bind_pair_eV` fold, `_amu_ang2_ps2_to_eV`); the
   three permanent lanes (`frozen`/`suppressed`/`time_exhausted`, defensive
   raise on an impossible k=0 in-band state); P1–P3 handover guard with the
   **gate-conditional density factor** (erfc ρ under `density_scaled`,
   exactly 1 under `none`), `EPS_DRAIN = 1e-6`, per-ion violator list +
   remedy; E2's two seed-coherence guards; realized-`t_h` axis check;
   stage-private stream `DETECTION_STREAM_KEY = 0xD5_2026` (zero draws on
   permanent seeds); m↔n lockstep assert; full per-event ledger (times,
   pre-shed rungs, K1/fold/defect) as ragged flat arrays + offsets;
   `save_detection_result`/`load_detection_result` for **`detection.npz`**
   (own schema counter v1, `allow_pickle=False`, fail-loud shape/offsets/
   reason validation; **IonCheckpoint v7 untouched**).
2. **Config surface** — `detection_stage_enabled` (default False),
   `detection_time_ps` (Optional, no default; Sourced row 24 supplied per
   run), `check_detection_config` wired into `validate()`: biphasic; time
   set > 0 and beyond the nominal seed-window end (ion window +
   relaxation cap when enabled); **ν > 0** (ν = 0 warns at the biphasic
   guard but is refused here — the permanent-state taxonomy is unsound
   without a live channel). Relaxation stage **not** required (skip path).
3. **E1** — `"detected"` added to the legal source tags; a
   `DetectionResult` is inferred `"detected"` via its `state_reason`
   attribute (duck-typing hook, `terminal_n` property unchanged).
4. **Probe report** — six new columns (`n_detect_mean/spread/min`,
   `frac_det_{frozen,suppressed,time_exhausted}`) read from an **optional**
   per-dir `detection.npz` (`-` when absent; `_REQUIRED_ARTIFACTS`
   unchanged, so every pre-DS dir stays scorable).
5. **E2 contract re-frame** (docstrings only, zero behavior change):
   ungated = terminal solver to freeze-out; gated = short handover bridge;
   relaxed read demoted to convergence diagnostic.
6. **Tests** — `tests/test_detection_stage.py`, 37 tests along the design
   §5 seven-step plan: hand-oracled waiting times on the n=1 direct chain;
   pinned default-stream derivation + zero-draw permanent lanes + stream
   key; the three k=0 lanes incl. the ungated no-op; `time_exhausted`;
   guard arms (non-biphasic, stale stride, t ≤ t_h, P3 in-bubble, P3
   ungated-live); config arms + the **Wave-7 back-compat criterion**
   (pre-DS cfg.json loads with defaults — met, no regeneration);
   **hypoexponential analytic oracle** (400 iid ions on a ≥3-stage
   distinct-rate chain built from the delivered `rrk_rate`; occupancy
   within 3σ + mean within 4σ); cold-shed kick composition
   |v_f| = |v₀|·m₀/m_f + m↔n lockstep; **5-term closure across the stage
   boundary** (atol 1e-9 — pure arithmetic chain, no integrator drift;
   per-event sums == field deltas; E_dissip strictly carried);
   **fixed-dt E2 equivalence** (400 ions, 1 ps horizon, |Δn̄| < 0.45 ≈ 4σ
   + Bernoulli-bias allowance); skip-path lanes (frozen ion-seed ==
   relax-seed chain; real ion.npz-shaped seed accepted with relaxation
   disabled); artifact round-trip + three fail-loud load lanes; E1
   admission; probe-report smoke (columns `-` before, populated after).

**As-built deviations (recorded in the design doc's status block + §2.2/
§2.4/§2.5 in place):** the approved draft's "cold shed leaves v unchanged"
was **wrong** — the delivered `cold_shed` applies the momentum-conserving
`m/(m−m_He)` kick per fire (KE rise cancelled by the booked
`E_mass_transfer` defect; the (E_int, n) jump chain never reads v, so the
exactness claim is untouched). The artifact carries terminal velocities +
per-ion terminal ledger channels beyond the §3.2 list so closure is
checkable from `detection.npz` alone. P3's density factor is
gate-conditional; ν > 0 required.

**Not in this slice (unchanged scope):** Wave-7 execution (behind its own
go-ahead; probe plan C.2 NB), the manual production-wiring oracle (gated
s_eff = 30 → t_detect), campaign/F2 wiring (parked), generator changes
(none — the pipeline hands the stage a cfg view + seed checkpoint).
**Next:** connected code review (DS build + E2 re-frame, per the §1-item-5
review scope), then the production-wiring oracle and Wave 7.

### Slice DS code review — EXECUTED, fixes applied (2026-07-07/08)

The §1-item-5 **connected review** ran (multi-agent workflow, high effort:
4 finder angles → 16 verified candidates → 13 independent verifiers → 10
confirmed findings, **0 refuted**), covering the DS build + the E2 re-frame
together. All 10 findings fixed the same session; **suite after fixes:
2176 passed** (+3 net new tests). Review question ("is the stage boundary
and the demoted E2 contract coherent") answered: the composition (no
duplicate physics) and the 5-term closure held; the **P1–P3 guard did
not** — three confirmed correctness defects, all in the guard/taxonomy or
its consumers:

1. **P3 exemption of suppressed ions was unsound (the substantive physics
   finding).** Suppression is *cooling-reversible* (K2 drains E_int below
   Σ(n) and re-opens the gate — `biphasic_step` cools before the gate for
   exactly this reason), so a suppressed ion with live cooling is NOT
   permanent; the draft guard silently rode such ions to the detector at
   their handover n, biasing the Tier-2 arbitration observable high.
   **Fix:** the cooling bound now applies to every ion not at the energetic
   floor (live + suppressed); only strictly-frozen ions are exempt from it.
   Consequence: ungated non-frozen ions can never hand over (must freeze in
   E2) — now enforced, not just documented.
2. **No bound at all on frozen in-bubble ions (the skip-path hole).**
   Pickup (re-heats E_int by f_ret·D₀ per capture — can *unfreeze* the
   cascade) and drag are density-gated and only *omitted* by the stage.
   **Fix:** a helium-exposure bound `ρ·max(λ₀, 1/τ)·(t_detect−t_h) ≤
   ε_drain` applies to ALL ions, frozen included — the skip path's actual
   defense. The delivered test that showcased in-bubble suppressed
   acceptance was inverted into the rejection oracle it should have been.
3. **RRK-underflow crash lane.** In-band ions with E_int within ~1e-6
   relative of D₀(n) underflow the bracket `(1−D₀/E)^(s−1)` to k = 0
   (s−1 up to 59); `_permanent_reason` fell through and the defensive raise
   aborted the whole run. Reachable (E2's freeze mask is strict `<`).
   **Fix:** classified `frozen` (waiting time exceeds any flight by
   hundreds of orders); the raise now covers only impossible states.
4. **Report dropped skip-path runs.** `relaxation.npz` was still a required
   artifact. **Fix:** optional — skip-path dirs (cfg + ion + detection) are
   discovered and scored, relaxed columns `-`; plus a **stale-detection.npz
   coherence guard** (ensemble-size vs cfg mismatch fails loud — the
   regenerated-dir hazard).

Cleanups applied with them: point-of-use `NotImplementedError` refusals for
the `tabulated` ladder/density arms (the `biphasic_step` precedent);
`fold = −dE_int` (bit-exact K1/fold cancellation; the e_bind-difference
recomputation left ~1e-17 eV per-fire residues + two redundant ladder
lookups); the two copy-pasted stage-boot blocks hoisted to shared helpers
`checkpoint.check_biphasic_seed_checkpoint` / `checkpoint.stage_stream_rng`
(E2 rewired to them — behavior-identical, RNG tests pin it); E1 doc-drift
("Two input modes"/"two legal tags" → three) fixed.

Documentation updated in place: design doc §2.1 (review-hardened guard NB +
new P1/P2 bound), §2.3 (underflow-frozen lane; suppressed-permanence
caveat), status block (review record); probe-report module docstring
(optional artifacts). New tests: ungated-suppressed rejection,
frozen-in-bubble rejection, underflow-frozen classification, skip-path
discovery/scoring, stale-detection coherence. **Slice DS is now delivered
AND reviewed. Next: the manual production-wiring oracle (gated s_eff = 30 →
t_detect) and Wave 7, each behind its own go-ahead.**

### Wave-7 detected-read re-score EXECUTED (2026-07-08)

Under the `[PROCEED TO IMPLEMENTATION]` trigger (user approved the four
discussed open choices: scratchpad-driver route, log-spaced sensitivity
band, all-102 scope, pre-registered predictions W7-P1..P4). Scratchpad
driver through the delivered `run_detection_stage` — **zero repo-code
change, zero new MD runs**: every one of the **102** probe dirs re-read at
the Sourced t_detect = 8.53 µs (CALIBRATION_MAP row 24), seeding from its
`relaxation.npz`; `detection.npz` written per dir; full re-score with the
delivered report (`n_detect_*` + state-reason columns on all 102 rows).
This also discharges the design's manual production-wiring oracle (the
gated s_eff = 30 point continued to t_detect) as part of the sweep.

**Execution note (recorded; findings §4e):** the 63 pre-Wave-4 dirs are
stamped with the legacy relaxation cap 8.53·10⁶ ps (chosen before Slice DS
existed), which the config-load *nominal*-window bound rejects at
t_detect = 8.53·10⁶; the transient cfg view therefore carried the
*realized* relaxation duration (their E2 terminated at all-frozen after
tens of ps) so the nominal and realized bounds state the same fact.
On-disk `cfg.json` untouched; the stage's realized-t_h re-check and the
P1–P3 physics guard ran on every dir — **zero violations**.

**Wiring:** ungated no-op identity exact on all 75 ungated dirs
(`n_detect ≡ n_relaxed`, zero events, all `frozen`); bridge point row
unchanged (0.67 / 20.33 / 20.33 / n_detect 20.33); ledger residual uniform
≈ 2.23·10⁻⁵ eV.

**Verdicts (tables + detail: findings §4e, insights I19–I22):** P1
confirmed; P3 confirmed exactly (all ten dead arms 100 % `suppressed` at
n = 21 — OQ-B's fragmentation question is now literal detector weight); P4
confirmed (early-open leak freezes survive as **converged mid-shell
weight**: n̄_detect 7.4–10.6 at s_eff = 8, f_int 0.24–0.35, 87–93 %
frozen); **P2 refuted in the informative direction** — the gated landing
point (s_eff = 30, τ = 6.55, f_int = 0.50) arrives **100 %
`time_exhausted` at n̄ = 4.09**, not floor-frozen: the descent is
logarithmic (band 6.21 → 4.09 over 10³ → 8.53·10⁶ ps, ≈ 0.4–0.9 He per
time decade), so the n ≈ 2 floor is kinetically asymptotic at high s_eff.

**Headline structural finding (I20): the detector read compresses s_eff.**
At f_int = 0.50 the whole s_eff ∈ [1, 30] range lands n̄_detect ∈
[2.7, 4.1] (the same points span 2.7–13.8 in-window): s_eff decides the
*arrival state* (frozen vs live), the race margin Δ× + in-bubble leak
decide the arrival *n*. Staircase and detector distribution are thereby
near-orthogonal targets (staircase → s_eff; terminal distribution → Δ×) —
the identifiability shape a campaign co-fit wants. The gated arm alone
spans n_detect ≈ 3–13 + the n = 21 suppressed class (I21); bare I⁺ stays
unreached at the detector (min n = 2; I22, OQ-B sharpened).

Documentation: findings §0 wave table / §4e / §5 (I19–I22) / §6 NB
(cap-truncation items historical) / §7 (OQ-B sharpened, OQ-E added) / §8
inventory; probe plan Addendum-C status flipped (Wave 7 EXECUTED, C.2 NB).
Documentation + scratchpad execution only — no repo code changed. The
probe program's Wave-7 deliverable (the detector-level flexibility map) is
complete; discussion resumes from the documented findings (the user
adjudicates; campaign/production scope stays parked).

### Post-Wave-7 physical interpretation recorded (user discussion, 2026-07-08)

The user's expectation ("I expected all to freeze") prompted the physical
synthesis, recorded as findings §4e "Physical interpretation" + **I23** +
two OQ-register updates (documentation only; no code changed):

1. The post-ejection gated cascade is an exactly **closed system**
   (G = E_int − Σ(n) shed-invariant); the detected read sits in the
   **Klots evaporative-ensemble regime** (k·t ≈ 1, logarithmic drift) —
   the never-freezing behavior and the s_eff compression are the textbook
   long-time structure of the mechanism, derived and quantified in §4e.
2. **OQ-F (new):** the eternal cascade is nonetheless a *construction* —
   K1 drains only D₀(n) per shed; the missing per-shed translational
   release ε ~ (E_int−D₀)/s is first-order over a cascade and would make
   cascades self-terminate. Mechanism-convention OQ of the RRK-dof class;
   document-only.
3. **OQ-B hypothesis (sharpened):** the suppressed class
   (E_int > Σ(n), rides inert at n = 21 by construction) physically
   fragments — bare I⁺ may *be* the crossed-after-ejection side of the Δ×
   race, making the experimental bimodality the cliff structure itself.
   Document-level; spec belongs to the OQ-B resolution, behind the
   trigger.

Discussion continued the same day on the second-largest experimental bin
(I⁺He, 17.5 %), recorded as findings §4e "top two experimental bins" +
**I24** + OQ-B/OQ-F coupling notes + **OQ-G (new)**:

4. **Both top experimental bins are structurally outside the cascade
   side's reach** (flat ladder bottom D₀ = 9.22 meV: bare needs |G| = 0
   exactly; n = 1 needs a sub-rung leak and is kinetically strangled —
   consistent with Wave 7's global min n = 2, n = 1 absent in 10 200
   detections).
5. **The never-opened side covers n = 1 only jointly with OQ-F's ε**:
   gateless boil-off without ε goes exactly to bare (G > 0 invariant; the
   n = 1 direct channel has no barrier); with ε, the one-sided G₀
   distribution yields bare (bulk) + a decreasing small-n tail (fringe),
   predicting a **budget-dependent bare peak** (G₀ ≈ 0.09 eV at 0.80 eV
   vs ≈ 0.35 eV at 2.70 eV against 21·ε ≈ 0.1–0.4 eV). OQ-B and OQ-F are
   one coupled discussion.
6. **OQ-G (new):** a deep first rung (real I⁺–He binding vs the flat
   Form-U bottom) is the competing/combinable explanation for n = 1's
   prominence; the small-n abundance tail is the first observable that
   reads the ladder *bottom*. Literature/domain-expert question;
   document-only.

### OQ2 fired (f_int parametrization) + RESEARCH_QUESTIONS.md created (user, 2026-07-09)

Discussion continued on the user's physical objection to the swept f_int
values ("f_int = 0.5 means half the Coulomb explosion in internal energy —
the explosion would be half as violent; actual values ~1 % or lower").
Recorded as findings §4e "OQ2 fires" + **I25** + §7 OQ2 entry
(documentation only; no code changed):

1. **Code fact (verified):** the model implements no partition —
   `E_int(0) = f_int·E_avail` is deposited at t0 with nothing subtracted
   from the Coulomb mechanics (`ion_initial_state.py` → `e_int_onset_eV`).
   At f_int = 0.5 the model books 0.40/1.35 eV with no mechanical source;
   a *true* 0.5 partition would slow fragments by √2 and break the
   Tier-0/1a and VMI velocity physics.
2. **No comfortable literal value exists:** ≥ 0.24 is needed for any
   crossing at 0.80 eV; ≤ ~0.01–0.05 is defensible as literal KER
   coupling but leaves the mechanism inert (8–27 meV vs Σ(21) = 188 meV
   → 1–4 sheds, frozen n ≈ 17–20). The parametrization is broken, not the
   number.
3. **Row-14 floor identity:** the scenario-keyed floors multiply out to
   the same absolute energy (0.235·0.80 ≈ 0.065·2.70 ≈ Σ(21)) — the
   natural variable is an **absolute E_int(0) [eV]**; scenario keying
   dissolves.
4. **Working hypothesis (pending literature validation, NOT
   adjudicated):** E_int(0) ≈ 0.2–0.5 eV budget-independent, from
   ionization reorganization + electronic/spin–orbit relaxation of
   nascent I⁺ (ties the picture knob to E_int(0) provenance). If
   confirmed: f_int reclassifies fraction → absolute Bounded eV; the
   B.1(2) production prediction inverts; the Δ× race coordinate survives.
5. **Phase decision (user):** the next step is **extensive literature
   research / domain cross-validation of the internal-energy handling**.
   Created **`docs/drag_port/Tier2/RESEARCH_QUESTIONS.md`** — the
   consolidated register (RQ1 = OQ2 E_int(0) provenance; RQ2 = OQ-F ε per
   shed; RQ3 = OQ-B suppressed fate/bare peak; RQ4 = OQ-G ladder bottom;
   RQ5 = OQ-E µs-flight channels; RQ6 = s_eff literature cross-check),
   each with model convention, problem statement, literature targets,
   discriminating observables, and couplings; per-RQ findings get
   recorded there with provenance, the user adjudicates, and any model
   change stays behind `[PROCEED TO IMPLEMENTATION]`.

### RQ1 literature research EXECUTED + ADJUDICATED — f_int → absolute E_int(0) (user decision, 2026-07-09)

The RQ1 deep literature run was executed (multi-agent: 5 search angles,
22 primary sources fetched, 90 claims extracted, top 25 through 3-voter
adversarial verification → 20 confirmed / 5 refuted / 0 unverified;
findings synthesized to 11 NBs). Full register with provenance:
`RESEARCH_QUESTIONS.md` → "RQ1 findings — NB register" + "RQ1
adjudication". Document-only; no code, checkpoints, or probe artifacts
touched. Key sourced facts:

1. **Electronic/spin–orbit channel is real and non-mechanical:**
   gas-phase strong-field I₂ CE populates excited I⁺ states *at the
   expense of KER* (Forbes 2022, assignment hedged, no branching
   fractions); levels anchored ³P₀ = 0.799, ³P₁ = 0.879, ¹D₂ = 1.702 eV
   (NIST). Ties the picture knob to E_int(0) provenance as hypothesized.
2. **Solvation reorganization measured:** Na⁺ sudden creation in He
   droplets releases ≈0.22 eV (measured, droplet-size-independent;
   He-DFT total ≈0.38–0.40 eV), dissipated by He ejection following
   Newton-type cooling over ~5 ps (Albrechtsen 2025) — external
   precedent for the K2 term and for the t₀-deposit-then-drain shape.
   Amends the working hypothesis' solvation ceiling upward
   (≤0.19 → ~0.22–0.40 eV).
3. **Kinematic KER→shell leak: exists but unquantified.** Speed–N
   correlation (Braun & Drabbels, *neutral* fragments) contradicts full
   budget-independence, but **both** quantitative coupling estimates
   (~14 % effective-mass; ~1 % TDDFT) were refuted 0-3 — magnitude
   genuinely open. Surface alkali-dimer CE shows no large KER shift.
4. **Applicability caveat (dominant):** no source measures the target
   system (in-droplet I₂ CE → two recoiling I⁺Heₙ); everything transfers
   by analogy (gas-phase CE / surface alkali cations / neutral
   photofragments). Neutral-vs-ion is the sharpest transfer risk. A Ba⁺
   counter-lead (LIF survives in liquid He) tempers the assumption that
   electronic energy must degrade non-radiatively into the shell (→RQ5).

**User adjudication (convention-level; RQ1 closed at that level):**

- Onset convention retained — E_int(0) deposited once at t₀, nothing
  subtracted from the fragment mechanics. The old "no mechanical source"
  objection dissolves: the dominant channels are non-mechanical, so
  no-subtraction is correct bookkeeping, not a naive simplification. The
  kinematic leak stays un-modeled, absorbed by the sweep.
- **Sweep variable reinterpreted fraction → absolute:** the physical
  variable is `E_int(0) [eV] = f_int·E_avail`; sweeps are specified,
  reported, and budget-transferred in absolute eV (0.80 → 2.70 eV holds
  E_int(0) fixed; f_int rescales ×0.80/2.70 ≈ 0.296). Sourced band
  0.2–0.5 eV (floor ≈0.22 eV), discounted tail to ~1 eV.
- The RQ1 "if confirmed" branch is taken: CALIBRATION_MAP **row 16
  Derived → Bounded (absolute), row 14 f_int Bounded → Derived
  coordinate** (one-for-one tally swap, edited 2026-07-09); row-14
  scenario keying dissolves; the B.1(2) production prediction inverts
  (t× budget-independent; production changes via the earlier t_eject).
- Discriminator refined: no smooth VMI KER deficit ∝ f_int expected;
  the electronic channel predicts fine-structure satellites (~0.8 eV
  down-shifted sub-population) — falsifiable both ways.
- Refuted-claims register kept in the RQ1 NBs so dead claims are not
  re-imported (notably "heavy-cation He-tag > bare ion" — do not cite
  toward RQ3).

Rule-2 carry status: unchanged — this is a documentation/convention
change; `e_int_onset_eV` and the config surface are untouched. Any
config rename or absolute-eV input arm is a future slice behind
`[PROCEED TO IMPLEMENTATION]`.

## Wave-8 suppressed-fraction cliff probe — design APPROVED (2026-07-09)

**Decision (user, after brainstorm):** Wave 8 probes OQ-B/RQ3 at the
weight level — does an in-band absolute `E_int(0)` split the fixed-seed
ensemble ≈ 43.5 % suppressed (bare-candidate under the fragmentation
hypothesis) / 56.5 % opened, and how wide is the Wave-6 cliff
(0.40–0.52 eV bracket at τ = 6.55)? Approach selected: **analysis-first**
(Approach B of the discussion) — a zero-MD step 1 computes the per-ion
critical values `E*_i = sup_t[Σ(n_i(t))·e^(I_i(t)/τ)]` from the on-disk
shed-free f_int = 0.65 dirs (per-ion cooling-exposure integrals through
the shared bubble surface), yielding the predicted
`suppressed_fraction(E_int(0))` ECDF; step 2 places only **7 targeted MD
runs** (3 transition points at s_eff = 30, one s-blindness companion at
s_eff = 8, and a 3-run 2.70 eV arm at fixed absolute E_int(0) — the
first-ever probe data at the production budget, testing the RQ1-inverted
"bare grows with budget" prediction). Suppressed stays inert (counted,
not evolved): no fragmentation channel, no new physics, no new config
surface. The optional τ = 16.5 Δ×-universality spot point was **dropped**
at approval (deferred, not rejected).

Full design: `TIER2_STAIRCASE_PROBE_PLAN.md` **Addendum D** (D.1–D.6,
pre-registered W8-P1..P5, outcome shapes a/b/c, boundaries). Sweeps
specified in absolute `E_int(0)` [eV] per the RQ1 adjudication.
Execution stays behind `[PROCEED TO IMPLEMENTATION]`; nothing here
discharges the F5 gate.

## Wave-8 cliff-anatomy probe EXECUTED (2026-07-09)

Executed under the `[PROCEED TO IMPLEMENTATION]` trigger, per Addendum D
(design approved earlier the same day). Scratchpad drivers through the
delivered generator / detection / report pipelines — **zero repo-code
change**; **6 new MD runs** (not the planned 7 — see below) + the zero-MD
step-1 analysis + detection over the new dirs + a full 108-dir re-score.

**Step 1 (zero-MD).** The per-ion critical value from the on-disk
shed-free f_int = 0.65 dirs is **one number, not a distribution**:
E*_i = Σ(21)·E₀/E_int,i(∞) = **0.461219 eV** with ensemble width
4·10⁻¹³ eV. Root cause verified in the artifacts: the probe ensemble is
**100 kinematically congruent replicas** (center-placed deterministic
onsets, one droplet radius → identical cooling exposures
K_tot = 0.898297). Analysis validation: from-positions replication of the
gated exposure exact to 1.9·10⁻¹⁵ eV; pre-open E_int linearity across
f_int to 9.4·10⁻¹⁶ eV with **bit-identical** pre-open positions; the
f0.50 opening time predicted/observed 6.33 ps. The "ECDF crossings" run
placement therefore degenerated: a delta has one crossing → two-point
bracket (the third 0.80 eV point had nothing to measure).

**Step 2 (6 runs, every pre-registered number hit).** b080 f_int 0.57
(E₀ = 0.4560 eV): 100/100 open at exactly the predicted **10.83 ps**,
|G| ∈ [1.79, 2.12] meV (bound 2.12); at s_eff = 8 all 100 ions freeze at
**n = 1** (first n = 1 in the program; n_detect = 1.00, 100 % frozen); at
s_eff = 30 kinetic parking at n_detect = 2.00. b080 f_int 0.58: 100 %
suppressed forever. b270 anchor (E₀ = 0.52 eV absolute): suppressed,
shed-free, and **K₂₇₀ = 0.898297 ≡ K₀₈₀** — `coulomb_available_eV` has
exactly one physics reader (the S2 onset deposit); the mechanics are
budget-blind. b270 f_int 0.17 / 0.18: same absolute step (opens at the
predicted 11.96 ps / never), 2/100 ions at n = 1 even at s_eff = 30.
Re-score: 108 rows, ledger residual uniform ≈ 2.23·10⁻⁵ eV, bridge +
gated fi0.50/fi0.65 wiring oracles exact.

**Verdicts:** W8-P1 **confirmed exactly**; W8-P2 + P3 **refuted** (the
cliff is a step — no E_int(0) splits the ensemble; the 43.5 % bare
fraction is not a knob outcome; heterogeneity must be *injected*:
outcome shape (b), pre-registered as equally decisive); W8-P4 confirmed
sharpened (the sub-rung sliver E₀ ∈ (E* − D₀(1), E*) expresses **n = 1**
— I24's kinetic strangulation holds only at a full-rung leak); W8-P5
**refuted in-model at the mechanism level** → **OQ-H fired** (production
Coulomb kinematics unmodeled; the absolute-E_int(0) gated map transfers
across budgets verbatim; B.1(2) void in both versions; the RQ3-vs-RQ4
bare-vs-budget discriminator currently untestable in-model).

Full record: `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4f, insights I26–I28,
OQ-B post-Wave-8 NB, new OQ-H, §6 congruence boundary NB, §8 inventory
(102 → **108** dirs; first `b270` namespace entries — their staircase
columns are N/A, scored against the 9 Å anchor by the shared scorer).
Addendum D status block updated in `TIER2_STAIRCASE_PROBE_PLAN.md`.
Nothing here discharges the F5 gate; the user adjudicates.

## Wave-9 E0-mixture inversion — design APPROVED (2026-07-09)

**Decisions (user, post-Wave-8 discussion):** (a) the E₀-mixture probe
is worth running before any droplet-distribution slice; (b) the fitted
p(E₀) is treated as a **physics claim** on E_int(0) provenance (RQ1
solvation + electronic branching; a multi-modal fit would be a genuine
fine-structure-branching prediction); (c) the s_eff-from-small-n-bins
identifiability is recorded as **possible but not yet fixed** — no
campaign re-scope.

**Design (Addendum E of `TIER2_STAIRCASE_PROBE_PLAN.md`):** ions are
non-interacting, so a weighted mixture of single-point runs is exactly
an ensemble with p(E₀) supported on the grid — the experimental
histogram inverts into an implied p(E₀) by non-negative least squares
over **empirical per-run detected distributions** (never the analytic
window map). Basis: 13 columns spanning E₀ = 0.08–0.52 eV at s_eff = 8
(floors reached), of which 5 + the suppressed class are already on disk
(Waves 4–8) — **8 new runs** (7 × s8 + 1 × s30 companion at the n = 2
window for the kinetic-wall read). Pre-registered: W9-P1 (top bins fit
with in-band mass), W9-P2 (the deep tail demands below-floor E₀ — the
quantified hand-off to the droplet-size axis), W9-P3 (smooth unimodal
p(E₀), mode ≈ 0.45 eV — multi-modality = RQ1 branching structure),
W9-P4 (s30 re-fit shifts n = 1 → n = 2; identifiability observation
only, per decision (c)). Conditional on RQ3 spec (b), the 9 Å exposure
(RQ7), and the single-radius preset — the droplet axis is deliberately
absent; its demand is W9-P2's output.

The post-Wave-8 heterogeneity discussion (histogram-as-readout map, the
two injection axes E₀ vs droplet-R exposure, the future
droplet-distribution slice and its constraints) is recorded in
`TIER2_STAIRCASE_PROBE_FINDINGS.md` §4f discussion block. Execution
stays behind `[PROCEED TO IMPLEMENTATION]`; nothing here discharges the
F5 gate.

## Wave-9 E0-mixture inversion EXECUTED (2026-07-09)

Executed under the `[PROCEED TO IMPLEMENTATION]` trigger, per Addendum E.
Scratchpad drivers through the delivered generator/detection pipelines +
a pure post-processing simplex fit — **zero repo-code change**; **8 new MD
runs** (7 × `s_eff = 8` at f_int ∈ {0.53, 0.48, 0.45, 0.42, 0.20, 0.15,
0.10} + 1 × `s_eff = 30` companion at f_int = 0.53), each with
`detection.npz` at generation → **108 → 116 dirs**. Wiring oracles exact:
ledger residual uniform 2.231·10⁻⁵ eV on all 8 new dirs; on-disk fi0.50
s8/s30 + fi0.65 re-score unchanged (n_detect 3.05 / 4.09 / 21.00). New
columns land on the flat-bottom window map (E₀ → n̄_detect: 0.424 → 2.04,
0.384 → 3.38, 0.360 → 4.31, 0.336 → 5.22, 0.160 → 12.26, 0.120 → 14.48,
0.080 → 16.62, all frozen; s30@0.424 → 3.02, 100 % time_exhausted — the
kinetic wall).

**Execution decisions (user, confirmed at trigger):** (1) **simplex**
LS (non-negative + Σw = 1), not raw NNLS; (2) suppressed→bare per-ion
under RQ3 spec (b); (3) deep-tail shortfall **reported not forced**;
(4) W9-P4 a localized column swap, **no campaign re-scope**.

**Verdicts.** W9-P1 **confirmed** — the histogram inverts cleanly:
L2 = 0.016, Wasserstein-1 = **0.086 bins**, top six bins < 0.006 (bare
fed solely by the suppressed→bare column, weight 0.436 ≈ 0.435; near-cliff
density on the 22.6 meV rung window 0.780 %/meV at n = 1 vs predicted
0.77, ratio to n = 2 = 2.08× vs 2.2×). W9-P2 **confirmed** — deep tail
forces **≥ 9.6 %** below the 0.22 eV solvation floor (E₀ ≤ 0.192 eV
columns) + n = 18–20 beyond the s = 8 basis floor: the droplet-radius
axis's quantified demand (in-band [0.22, E*] mass 0.467; above-E* bare
class 0.436). W9-P3 **refined** — the implied p(E₀) is smooth and
unimodal but **concentrated at/above E\* ≈ 0.46 eV**, an in-band declining
tail below, **not** an interior 0.45 eV bump; **no multi-modal / electronic-
branching signature** (interior roughness is a fi0.45/0.48/0.50
collinearity artifact) — provenance read: E_int(0) at the upper edge of /
above the RQ1 band. W9-P4 **confirmed** — s = 30 refit of the near-cliff
columns kills the n = 1 bin (resid −0.175), L2 → 0.182, W₁ → 1.398 (16×):
s_eff detector-identifiable via n = 1, **recorded as possibility only**.

Full record: `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4g + insights I29–I32,
§0 wave table, §8 inventory (108 → **116**). Addendum E status block
flipped to EXECUTED in `TIER2_STAIRCASE_PROBE_PLAN.md`. The fit builds
**no p(E₀) sampling surface** and is conditional on RQ3 spec (b); nothing
here discharges the F5 gate; the user adjudicates.

## Waves 8 + 9 — code review + fixes applied (2026-07-09)

Extensive review of the eleven Wave-8/9 scratchpad scripts (Wave 8:
`step1_cliff_anatomy.py`, `w8_driver_080/270_anchor/270_pair.py`,
`w8_detection.py`, `w8_rescore.py`, `w8_verify.py`; Wave 9: `w9_generate.py`,
`w9_fit.py`, `w9_inventory.py`, `w9_ledger_check.py`). Both waves ran the
**zero-repo-change** scratchpad route, so the review targets the analysis/driver
scripts, not a repo diff. Load-bearing assumptions were verified against the
actual repo (the gated K2 cooling law in `ion_propagation_step.py`, the
`drag_gate_steepness` resolver in `ion.py`, `detection_stage` fields, the
`integrated_i_he_abundance.csv` contract, and a real probe `cfg.json`).

**Verdict:** the core science is correct and every documented Wave-8/9 number
reproduces exactly — the `E*_i` derivation (shed-free ride → Σ constant →
`sup_t[Σ·e^{I/τ}] = Σ(21)·E₀/E_end`, read at the converged relaxation-end
`E_end`), the from-positions cooling replication (pre-step positions `rho[:, :-1]`
vs post-step `E[:, 1:]`, matching `newton_cool_step` at fixed n), and the
simplex mixture inversion (normalized columns + simplex weights → normalized
mixture; `W₁ = Σ|ΔCDF|`). The issues found are latent-robustness / reporting-
discipline items, **not** result-invalidating bugs.

**Findings.**
- **F1 (latent correctness).** `step1_cliff_anatomy.py` read `cfg["drag_gate_
  steepness"]` raw, but the sim resolves the gated K2/pickup steepness via
  `drag_gate_steepness(cfg)` = `cfg.potential_steepness` for the default
  `density_proportional`/`erf_tied` gate (only `erf_independent` reads the
  `drag_gate_steepness` field). Correct here **only by coincidence**
  (`drag_gate_steepness == potential_steepness == 14.2`, gate
  `density_proportional`); impact was bounded even so, since the steepness feeds
  only the cross-check + `i_open_pred`, never the deliverable `E*_i` (stored
  `E_end`).
- **F2 (reporting discipline).** Validation residuals were `print`-only:
  `step1`'s `drop_relax`, replication `err`, pre-open `scale_err`/`pos_err`, and
  `w9_fit`'s `total_over` leak count ("expect 0"). A steepness/convention
  mismatch or a leaked `n>20` ion (silently `np.clip`-folded into the top bin,
  distorting the fit) would not halt.
- **F3 (latent).** `step1`'s pre-open position-congruence prefix used ion 0's
  opening column as the "all-suppressed" prefix; correct prefix is the first ion
  to open (`min` over ions). Benign under the delta cliff (congruent ensemble),
  wrong if the ensemble ever spreads.
- **F4 / F5 (cleanups).** `w9_fit.column_distribution` annotated `-> np.ndarray`
  but returns a tuple; `simplex_lsq` ignored `res.success`.
- **F6 (cosmetic, not fixed).** `w8_rescore.py` docstring "109 dirs" (108 after
  Wave 8); in-loop numpy import in `w8_detection.py`; `w9_generate.py` skip-policy
  looser than the Wave-8 drivers' clobber guard. No correctness impact.

**Interpretation caveat (recorded, physics not code).** In the Wave-9 fit the
bare class is carried by the single `f_int = 0.65` column and its weight is placed
at `E₀ = 0.52 eV` in the implied-p(E₀) table, but **any** `E₀ > E* = 0.461` yields
the identical fully-suppressed bare delta — so the above-E\* mass *location* is
unidentifiable by the mixture. The findings already report it as one lumped
"mass above E\*" (0.436); the p(E₀) density near 0.52 eV should not be read as a
resolved feature. Given decision (b) treats p(E₀) as a provenance claim, this
degeneracy is worth an explicit line in §4g.

**Fixes applied (scratchpad scripts only; zero repo-code change).**
- `step1_cliff_anatomy.py`: **F1** — steepness now `drag_gate_steepness(
  RunDirectory(...).load_cfg())`, the sim's own resolver (yields the same 14.2,
  result unchanged); **F2** — asserts (`< 1e-9 eV` / `< 1e-9 Å`) on `drop_relax`,
  the replication `err`, and the pre-open `scale_err`/`pos_err`, so a
  steepness/convention mismatch now halts; **F3** — prefix uses `obs_open.min()`.
- `w9_fit.py`: **F2** — `assert total_over == 0` on the leak check; **F5** —
  `simplex_lsq` raises on `not res.success`; **F4** — return annotation
  `tuple[np.ndarray, int]`.

**Verification.** Both fixed scripts pass all asserts and reproduce the record:
`step1` → `E* = 0.461218816 eV`, `K_tot = 0.898297`, replication `1.9·10⁻¹⁵ eV`,
`pos_err = 0`, opening `6.33 ps`; `w9_fit` → leak `0`, **PRIMARY**
`L2 = 0.01609 / W₁ = 0.08647`, **W9-P4 refit** `L2 = 0.18211 / W₁ = 1.39830` —
matching §4g's `0.016 / 0.086` and `0.182 / 1.398`. No reported result changed.

## RQ2/RQ3 deep-research run + adjudication RECORDED (2026-07-10)

Joint RQ2 (per-shed KE release ε) + RQ3 (fate of the net self-unbound complex)
literature/cross-validation pass — the entry step of the research phase toward
the fragmentation channel that must precede a generative production run
(the discussion also fixed the program order: resolve RQ2/RQ3 → build the
fragmentation channel → build droplet-R heterogeneity → run the full-assembly
production run at R0 = 2.666 Å / 2.70 eV, with the 2.666 Å kinematics /
`n*(R)` check as its opening deterministic gate; RQ7 note below).

**Deep-research run.** 99-agent fan-out (5 angles → 17 primary sources → 39
claims → 25 adversarially verified, 20 confirmed / 5 refuted 0-3 → 7 findings),
recorded as the joint **RQ2+RQ3 NB register** in `RESEARCH_QUESTIONS.md`.
**Primary-source grounding.** The model's foundational paper —
Albrechtsen/Stapelfeldt, *Nature* 623, 319 (2023)
(`Stapelfeld_Paper_Ion_Solvation_in_Helium_Droplets.pdf`) — is itself a He
shell around a **cation** (Na⁺Heₙ, Coulomb-ejected, dissociating in flight,
VMI-detected) and directly settles both questions, upgrading the web
literature's transfer-by-analogy to a same-physics-class reference.

**User adjudication (2026-07-10).**
- **RQ2 — ε ≈ 0, small ε neglected (CLOSED).** A shed drains exactly D₀(n); the
  evaporated He leaves translationally cold. Sourced two ways: Stapelfeldt
  explicitly neglects the dissociation-product KER (near-threshold statistical
  dissociation), and Klots–Hansen gives only ε ≈ D/16 ≈ 0.5–0.6 meV
  (finite-heat-bath reduction toward 0 for the few-mode shell); the feared
  5–20 meV is ruled out (µs survival). Resolves findings OQ-F / the KE_shed
  half of CALIBRATION_MAP OQ2 (the partition half was RQ1) → **OQ2 fully
  closed**. Drain-D₀ convention retained, now sourced; a thesis-ready ε
  statement is recorded in the register.
- **RQ3 — spec (a) inert REJECTED; sequential shed-to-self-termination
  ADOPTED.** The over-energized complex sheds sequentially (Stapelfeldt MD:
  tens–hundreds of ps ≪ 8.5 µs), never rides intact. A genuinely
  net-self-unbound complex (E_int > Σ, margin G > 0) **totally sheds to bare**
  (G-invariant under ε ≈ 0 — the bare-peak source); the marginal class
  self-terminates at small n on the ladder bottom. Endpoint deferred to RQ4.
  Reframes findings OQ-B.
- **RQ4 promoted with a sharpened, falsifiable target:** the incremental-binding
  **ratios** D₀(1):D₀(2):D₀(3) ≈ 2.2:1.5:1.3 (decreasing-from-n=1, which the
  Na⁺ flat-first-shell-plateau analogy does *not* supply), not merely "deep."
  Flat-ladder + structured p(E₀) vs deep-ladder + smooth p(E₀) are
  histogram-degenerate; RQ4 (structure) + the production budget test (RQ7) are
  the two independent handles.

**Doc edits.** `RESEARCH_QUESTIONS.md` — RQ2+RQ3 NB register (7 confirmed NBs,
5 refuted, open questions), primary-source grounding block, RQ2/RQ3 adjudication
block, and the RQ4 sharpened-target subsection; `CALIBRATION_MAP.md` — Update
2026-07-10 blockquote + OQ2 register line marked resolved. **No code** — the
implied model changes (remove the spec-(a) suppression in favour of a
sequential-shed continuation; reshape the ladder bottom per RQ4) are deferred
behind `[PROCEED TO IMPLEMENTATION]`; nothing here discharges the F5 gate. The
next research target is the RQ4 deep-research pass (I⁺Heₙ incremental-binding
ratios).

## Wave-10 Addendum F AMENDED pre-execution — K-reframing + Step 1.5 + W10-P5 (user discussion, 2026-07-10)

Post-RQ4-register discussion ("what is K; does less cooling mean more bare
I⁺?") re-derived the gated fate structure in closed form and fed back into
the planned-not-executed Wave 10 (`TIER2_STAIRCASE_PROBE_PLAN.md`
Addendum F). **User adjudication:** fold into Wave 10 by amending Addendum F
in place (no results depend on it yet) rather than opening a new wave.

**The physics recorded (new F.2b).** The Wave-10 test variable is the
dimensionless in-bubble cooling exposure `K = (1/τ)·∫(ρ_He/ρ_bulk)dt`
(`E_int(t_eject) = E₀·e^{−K}`; measured 0.898297 at the pinned 9 Å
condition, Wave 8). Droplet radius, τ, and the production kinematics (RQ7)
are three handles on the one variable. At fixed sharp E₀ the fate map is
closed-form: **bare ⇔ K < K\* = ln(E₀/Σ(21))** — the *under-cooled* class
(cooling is the only process that can close the gate before ejection; "less
cooling → more bare" holds, as a cliff, not a dial); just-above-K\* → deep
strip n ≈ 1–3 (the n = 1 sliver); K ≫ K\* → shallow strip mid/high-n. The
no-shed leak formula `Σ(21)·(1 − e^{−(K−K*)})` reproduces the Wave-6
relaxed floors (f_int = 0.24 → n ≈ 11 vs measured 10.96; 0.50 → n ≈ 2–3 vs
2.7–3.1) — MD-validated at the pinned K. Consequence: one K-distribution
straddling K\* serves **both** experimental tails (bare below; retained
mid/high-n above) — the width Wave 9 put in p(E₀) may be a K spread.

**Amendments to Addendum F (all pre-execution, document-only):**

1. **F.2b** — the K reframing + closed-form fate map (above).
2. **F.3b / Step 1.5** — a zero-MD semi-analytic K forward model: calibrate
   K(R_droplet) from the existing probe trajectories (anchor K = 0.898297 at
   ~27.9 Å + transit scaling), propagate a stated droplet prior to a bare
   fraction P(K < K\*) + leak→terminal-n histogram + Wasserstein vs
   `integrated_i_he_abundance.csv`, scanned over sharp E₀ ∈ ≈ [0.24, 0.32] eV
   × two K-scales. De-risks Step 2: the droplet slice is built only against
   a live (median, width) target, otherwise demoted to tail-shape duty.
3. **W10-P3 direction corrected** — bare at narrow E₀ comes from the
   **small-droplet / early-ejection / low-K** end (K\*(0.28 eV) ≈ 0.398 vs
   the pinned 0.898), not the original "implausibly deep/slow-cooling
   droplets" (deep droplets *raise* K and serve the shallow-strip mid/high-n
   weight); P3's falsification now requires failure at both K-scales.
4. **W10-P5 (new, pre-registered)** — the K-scale / RQ7 leg: at the
   production-kinematics estimate (fragment speed ×√(2.70/0.80) ≈ 1.84 →
   K₀ ≈ 0.49, E\* ≈ 0.31 eV) the cliff lands near the ensemble median and
   narrow E₀ ≈ 0.28 eV plausibly splits ≈ 43/57 with modest droplet width.
   Prediction: P3 fails at the probe K-scale, succeeds at the production
   scale → resolution "narrow E₀ + RQ7 kinematics", the RQ7 arm a
   **co-requisite** of Step-2 adjudication. F.8 outcome (b) gains the
   corresponding second exit (only both-K-scale failure routes to the RQ1
   super-solvation source). The 0.49 estimate is a speed-scaling bracket,
   not model output — RQ7 stays unmodeled.

**No code, no runs** — Wave 10 remains planned-not-executed; Steps 1/1.5 are
zero-MD scratchpad work under the trigger, Step 2 / F.5 gated as before;
nothing here discharges the F5 gate. Cross-refs: `RESEARCH_QUESTIONS.md` RQ1
(the narrow-E₀ argument), RQ2/RQ3 adjudication (suppressed → total shed makes
the under-cooled class the bare source), RQ4 register (flat ladder → the tail
reverts to reservoir/race), RQ7 (the K-scale conditionality).

## Wave-10 Steps 1 + 1.5 EXECUTED — zero-MD width decomposition + K forward model (2026-07-10)

Executed under the `[PROCEED TO IMPLEMENTATION]` trigger, same day as the
Addendum-F amendment: scratchpad route (`wave10_step1_step15.py` + three CSV
scoreboards in the session scratchpad), **zero new MD, zero repo-code change,
no run artifacts modified** — a pure re-read of the 41 on-disk gated
(`_cgds`) probe dirs against the F.2b closed form and
`data/reference/integrated_i_he_abundance.csv`. Full results: findings
**§4h**, insights **I33–I36**. Headlines:

- **Step 1 (W10-P1 CONFIRMED).** Per-ion exposure integrals off the stored
  trajectories reproduce the Wave-8 anchor (K = 0.897541 vs 0.898297,
  0.08 %; per-ion E\*_i spread 4.15·10⁻¹³ eV = the Wave-8 delta re-derived);
  every gated column's n_detect at fixed (E₀, droplet) is ≤ 1 bin wide
  (MAD 0.04–0.88 He). The histogram width lives on the E₀/K axes only.
  Bonus: the b270 columns re-confirm absolute-E₀ budget-invariance in the
  detected read.
- **Step 1.5 validation.** The (E₀, K) fate map (bare ⇔ K < K\*; no-shed
  leak + exact ε = 0 descent) matches the 13 measured s_eff = 8 columns to
  0.84 He mean, ≈ exact at the cliff; pinned-trajectory
  K(R) has d ln K/d ln R ≈ 1.23 (E\*: 0.276 eV at N = 250 → 0.462 at 2000 →
  1.75 at 16000).
- **W10-P3 CONFIRMED at the probe K-scale — by ~19×, not a margin.** At
  K₀ = 0.898, narrow E₀ ≈ 0.28 eV yields ≤ 2 % bare under every stated
  droplet prior; 43.5 % bare needs log-normal δ ≈ 12 vs Kornilov 0.625.
  The droplet axis also fails the tail *shape* at this scale (W₁ ≥ 2.7
  bins) — the W10-P2 analytic pre-verdict is negative at probe K.
- **W10-P5 CONFIRMED — the split verdict.** At the production speed-scaling
  bracket (K₀ ≈ 0.49) the cliff sits mid-ensemble: (E₀ = 0.28 eV,
  δ = 0.80) lands **bare 42.6 %, W₁ = 0.993 bins, untuned**, monotone
  small-n envelope; n = 1 deficit residual 1.29× ≤ the ~1.4× F.5 taper cap
  (W10-P4 pre-read: coverable); mid-tail n = 2–4 overshoot open.
- **Consequence.** The amended F.8 outcome-(a)-via-RQ7 exit is the live
  resolution: **"narrow E₀ + RQ7 production kinematics"** — the RQ1
  super-solvation burden does not fire on this evidence, and **the RQ7
  kinematics arm is co-requisite for any Step-2 adjudication**; the droplet
  slice's job shrinks to the width/shape of the K-distribution.

**Boundaries (findings §4h):** production K-scale is a ballistic bracket
(RQ7 unmodeled, no 2.70 eV MD); K(R) is the pinned-trajectory approximation;
n_eject = 21 for all R in-model (physical incomplete dressing favours
narrow-E₀ further); conditional on RQ3 sequential-shed + ε ≈ 0. Analytic
pre-verdicts — **Step 2 (droplet slice, gated) arbitrates**; reported, not
auto-adjudicated; nothing here discharges the F5 gate. Doc edits:
findings header/§0 table/§4h/I33–I36; Addendum F status flipped to
"Steps 1+1.5 executed"; this entry. The user adjudicates the follow-through
(Step-2 build order vs the RQ7 kinematics arm).

## Post-Wave-10 program order SET — K₂.₇₀ first; narrow-E₀ reading endorsed conditional on RQ7 (user decision, 2026-07-11)

Discussion of the Wave-10 result (what did we get; can the experimental
signal be reached *without* a broad E₀ sweep). **User adjudications:**

1. **The narrow-E₀ reading is endorsed as the working picture, conditional
   on RQ7.** The E₀ ensemble-*width* demand is dissolved (one sharp E₀ in
   the solvation band + the physical droplet distribution); what survives is
   the scalar E₀ calibration sweep. Recorded as the Wave-10 NB under RQ1.
2. **A physical tens-of-meV solvation smear on E₀ is endorsed as useful**
   (sourced by the NB-RQ1-3 reorganization physics; expected to help the
   n = 1/n = 2 bins) — to be carried as a sensitivity leg of the droplet
   slice, NOT as a return to a broad p(E₀).
3. **RQ7 promoted to the active next target** (status block added under RQ7
   in `RESEARCH_QUESTIONS.md`): K₂.₇₀ is the single number that decides the
   branch — near ~0.5 → narrow E₀ stands; near the probe 0.9 → the
   super-solvation route reopens and the Wave-10 landing inverts.

**Program order (supersedes/refines the 2026-07-10 RQ2/RQ3-entry order
"fragmentation channel → droplet-R → production run" by putting the
kinematics measurement first):**

- **Move 1 — RQ7 kinematics probe (next wave; smallest, most decisive).**
  Deterministic Tier-0-style run at the production separation
  R₀ = 2.666 Å (2.70 eV/fragment; the 9 Å condition is the same Coulomb
  formula at 0.80 eV) with the locked drag bundle; **measure K₂.₇₀** via the
  Wave-10-validated exposure integral, replacing the ballistic ×0.545
  bracket (pre-registered expectation: K₂.₇₀ ∈ ≈ [0.5, 0.7] — drag eats part
  of the extra speed). Independent discriminator: final fragment speeds vs
  `vmi_iplus_he.csv` (pins K₀ with no evaporation physics). Design doc must
  argue drag-calibration validity at ~1.8× speed and include the n\*(R₀)
  opening check. Behind `[PROCEED TO IMPLEMENTATION]`.
- **Move 2 — sequential-shed detection arm (RQ3 follow-through).** Replace
  the spec-(a) "suppressed rides at n = 21" detection semantics with the
  adjudicated sequential shed-to-self-termination continuation (design doc +
  interchangeable enum arm) so suppressed → bare is literal in the MD read,
  not an analytic re-map. Behind the trigger.
- **Move 3 — Step-2 droplet slice, sharpened targets.** Built only after
  Moves 1–2, against the Step-1.5 quantified targets: (i) verify the
  (E₀, K) forward model where valid; (ii) arbitrate the n = 2–4 overshoot
  (no-shed bias vs prior shape vs n_eject(R)); (iii) carry the
  tens-of-meV E₀-smear sensitivity leg. **RNG-draw-order discipline:** both
  new sampling surfaces (droplet radii, E₀ smear) must draw from a separate
  generator or after all existing draws (forbidden-list item) — to be
  stated explicitly in the slice design.
- **Parallel track (non-blocking) — RQ4 external many-body I⁺Heₙ
  calculation** (the n = 1 residual, F.5's non-tuned gate); RQ5 stays
  parked (the detector read behaved textbook-Klots).
- **Before any campaign spend — F5 re-scope adjudication.** The original F2
  Stage-1 s_eff×τ co-fit predates the probe program's findings (detector
  nearly s-blind; staircase = showcase; the histogram's information lives on
  (E₀, K)). The N = 500 campaign is to be re-scoped around
  (E₀, droplet prior, K validation) with s_eff/τ as reported secondaries —
  a standing adjudication item, not yet designed.

**Doc edits this entry covers:** `RESEARCH_QUESTIONS.md` — Wave-10 NB under
RQ1 (narrow-E₀ + smear + measurability consequence) and the RQ7
promoted-status block; this log entry. No code, no runs; every Move stays
behind `[PROCEED TO IMPLEMENTATION]` with its own design document.

---

## Wave-11 (Move 1, RQ7 kinematics probe) DESIGNED — Addendum G added (user discussion, 2026-07-11)

Move 1 of the post-Wave-10 program order detailed into its design document:
**`TIER2_STAIRCASE_PROBE_PLAN.md` Addendum G — Wave 11: measure K₂.₇₀.**
Design decisions recorded there; the load-bearing ones:

- **Primary leg = the anchor-congruent suppressed ride**, not a literal
  fixed-mass Tier-0 run: `R0_GS_angstrom = 2.666` + E_int(0) = 0.52 eV (the
  Wave-8 2.70 eV anchor value; suppressed for any K < K\*(0.52) = 1.02 →
  shed-free, mass-constant m = I⁺He₂₁, no evaporation/pickup draws) — so
  S_K = K₂.₇₀/K₀.₈₀ compares like with like against the Wave-8 K = 0.898297.
  A literal `mass_scenario=fixed` m_eff leg is retained as an *optional*
  convention-sensitivity check only.
- **The two Move-1-mandated validity items are in-wave measurements, not
  prose:** (i) drag-band validity at ~1.84× speed → report the beyond-band
  exposure share of K; (ii) the n\*(R₀) opening check → report the
  overlap-segment share (pair separation < 2×4.67 Å) + the geometry
  statement. Both identified biases push the measured K *up* → K₂.₇₀ is an
  upper bound in those respects, which makes a ≤ 0.7 landing conservative
  and forces any inversion verdict to survive the decompositions first.
- **Step 3 closes the loop zero-MD:** re-run the F.3b scan with the measured
  K₂.₇₀ replacing the ballistic ×0.545 bracket; W11-P5 pre-registers the
  Wave-10 landing surviving the replacement. Decision structure and outcome
  shapes (a)–(d) pre-registered in G.1/G.7, including the honest exit (d):
  validity failure → drag-recalibration research item, no silent
  extrapolation.
- New opening oracles: config-load acceptance at 2.666 Å (a refusing guard
  is a finding), drag-off ballistic integrator check at the ~11× steeper
  Coulomb onset, tag non-collision (first dirs off 9 Å geometry — tag must
  carry budget *and* geometry).

**Doc edits this entry covers:** `TIER2_STAIRCASE_PROBE_PLAN.md` (new
Addendum G, G.1–G.8, predictions W11-P1–P5); `RESEARCH_QUESTIONS.md` — one
pointer line under the RQ7 promoted-status block. **No code, no runs** —
the MD leg (1 primary dir + optional legs) stays behind
`[PROCEED TO IMPLEMENTATION]`; the control re-read and Step 3 are zero-MD
scratchpad re-analysis under the standing convention. Nothing here
discharges the F5 gate.

---

## Wave-11 TRIGGERED — opened companion promoted to default (user adjudication, 2026-07-11)

Two user decisions, same message: (1) the **opened companion
(E_int(0) = 0.28 eV at the 2.666 Å geometry) is promoted from optional to
default** — Wave 11 runs two dirs (suppressed primary at 0.52 eV + opened
companion), leaving only the fixed-mass convention-sensitivity leg
optional/boundary-conditional; (2) **`[PROCEED TO IMPLEMENTATION]` given**
for the Addendum-G MD leg. Addendum G status blockquote, G.3 run matrix,
G.5 Step 4, and G.8 gating updated accordingly. Execution record follows
in this log when delivered.

---

## Wave-11 DELIVERED — K₂.₇₀ = 0.746 measured; W11-P2/P5 refuted as registered; landing re-calibrates in-band; 42 % beyond-band caveat (2026-07-11)

Executed same day as the trigger, per Addendum G as amended (companion
default). Scratchpad route (`wave11_gen_runs.py`, `wave11_measure_K.py`,
`wave11_step3_scan.py` + two scoreboard CSVs, scratchpad); **two new MD
dirs** under `data/runs/` — the first off 9 Å geometry, tag =
probe tag + driver-level `_R2.67` suffix:
`9A_drag_shared_pure_cubic_N50_tier2probe_b270_mix_k1.00_l0.90_fi0.19_fr0.10_tau6.55_s8.00_cgds_R2.67`
(suppressed primary, E₀ = 0.52 eV) and `..._fi0.10_..._cgds_R2.67`
(opened companion, E₀ = 0.28 eV); both carry the full
neutral/ion/relaxation/detection artifact set. **Zero repo-code change.**
Full numerical record: `TIER2_STAIRCASE_PROBE_FINDINGS.md` **§4i** +
insights **I37–I40**. Headlines:

- **Step 0 all green:** config-load accepts 2.666 Å; tags free; the
  pure-cubic b = 0 refusal fired correctly (oracle ran drag-free at
  b = 10⁻¹²); ballistic dt oracle passes at 8.8·10⁻⁵ once the bundle's
  0.1168 eV effective-binding well is included in the analytic reference.
- **W11-P1 exact pass** (control re-read K = 0.897541, spread 10⁻¹³);
  **K₂.₇₀ = 0.74603** — **W11-P2 refuted** (above [0.5, 0.7]; S_K = 0.831,
  the ballistic ×0.545 bracket is wrong — drag eats the extra speed);
  **W11-P3 confirmed** (congruence 10⁻¹³; opened-companion m(t)-feedback
  −2.2 %); **W11-P4 split** (overlap share 8.0 % in-band; **beyond-band
  share 41.7 % — far above pre-registration, the load-bearing caveat**;
  both biases up → 0.746 = upper edge, ballistic 0.49 = hard lower
  bracket); **W11-P5 refuted as registered** (pre-registered
  E₀ ∈ [0.24, 0.32] reaches ≤ 18 % bare at measured K) — but the
  fate-map scan lands the bare peak at **E₀ ≈ 0.37–0.43 eV for every
  stated prior** (inside RQ1's [0.2, 0.5]), untuned best **W₁ = 0.496**
  (E₀ = 0.41, δ = 0.80), better than Wave 10's 0.993; the bare-crossing
  E₀ tracks E\*(K₀) over the whole K bracket (0.30–0.41 eV) so the
  caveat moves the calibration, not the verdict. **Closed form validated
  at production kinematics** (companion MD n_detect 5.76 vs closed-form 6).
  **First VMI-side constraint:** model production speeds undershoot the
  experimental I⁺He peak ~2–2.5× (over-drag direction, consistent with
  the beyond-band caveat).
- Incidental: detection-stage seeding requires the *saved* relaxation
  checkpoint (an `IonCheckpoint`), not the in-memory `RelaxationResult` —
  the design-§1-item-5 convention, hit and honored in the driver.

**Standing adjudication items for the user (reported, not
auto-adjudicated):** (1) whether the re-calibrated landing
(E₀ ≈ 0.38–0.41 eV at measured K) supersedes the Wave-10 narrow-E₀
endorsement's 0.28 eV working value — the RQ1 NB scalar sweep would move
up-band; (2) whether the **drag-recalibration research item** fires (the
42 % beyond-band share + the ~2× VMI speed undershoot both say the locked
pure-cubic law is exercised outside its validity at production speeds);
(3) how the worsened n = 1 deficit (1.5–2.9×, δ-dependent) re-weights RQ4
vs the E₀-smear leg for Move 3. **Doc edits this entry covers:**
`TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4i, I37–I40, status header);
`TIER2_STAIRCASE_PROBE_PLAN.md` (Addendum G status → EXECUTED);
`RESEARCH_QUESTIONS.md` (RQ7 measured-status note). Nothing here
discharges the F5 gate; Moves 2–3 stay behind their own design docs +
trigger.

---

## Two-channel bare-I⁺ reading ENDORSED — arbitration re-targeted to the solvated branch; RQ8 opened; RQ4 promoted to critical path; Waves 12/13 designed (user adjudication, 2026-07-11)

**Experimental input (user, 2026-07-11, post-Wave-11 discussion):** the
detected mean kinetic energy per fragment falls monotonically from
**~3 eV at n = 0 to < 0.1 eV at n = 17**, and the bare-I⁺ KE is
substantially larger than every snowball's — the experimental reading is
that **bare I⁺ is produced by a different ionization channel** than the
I⁺Heₙ ejection. (The delivered `vmi_iplus_he.csv` is the mass-131 = n = 1
gate — its 10.1 Å/ps peak ↔ 0.69 eV; exporter README `sqrt(127/131)`.)

**Discussion findings recorded (scratchpad pre-reads
`wave11_v_of_R.py` / `wave11_solvated_rescore.py`, zero MD):**

- **Class-conditioning check:** under the current law the speed–size
  correlation is compressed to ~1.5× across the whole prior (terminal
  choking; v(R) = 6.0/5.2/3.8/3.2 Å/ps at N = 250/500/2000/8000, detected
  est. 7.7/6.8/5.1/4.2) — no class reaches the n = 1-gated 10.1 Å/ps, so
  the Wave-11 over-drag conclusion survives correct conditioning.
- **Solvated-branch re-score (the two-channel target):** scoring n ≥ 1
  renormalized (experimental solvated n₁ = 31.0 %, n₂₋₄ = 30.5 %), the
  delivered construction **fails structurally at every (E₀, prior)**:
  best W₁ = 1.06, solvated n₁ ≤ 17 %, and only at high E₀ where the
  model over-produces now-surplus bare (54–68 %). Flat rungs feed flat
  small-n bins (n₁ ≈ n₂ ≈ n₃); the experimental steep decline
  (n₁/n₂ = 2.18) reads the **ladder bottom** — converging on RQ4's
  pre-registered target ratios (D₀(1):D₀(2):D₀(3) ≈ 2.2:1.5:1.3).
- **Vindication + cost:** the mechanism's structural refusal to make
  bare (I8/I17/I22) was never a bug — bare was another channel's job;
  conversely the Wave-9/10/11 "bare lands untuned" results become
  channel-branching coincidences, not model observables.

**User adjudications (endorsed as a package):**

1. **Arbitration observable re-targeted to the solvated branch:** the
   n ≥ 1 histogram renormalized + the solvated (n, mean-KE) curve
   (0.69 eV at n = 1 → < 0.1 eV at n = 17). The bare bin becomes a
   channel-branching quantity (upper bound / mixed), **not** a model
   target, pending RQ8.
2. **RQ8 opened** — bare-channel provenance & branching (candidates:
   inner-turning-point vertical ionization ~2.8–3.0 eV singly charged;
   I⁺+I²⁺ born at 5.4 eV partially dissipated; discriminators: I²⁺ at
   m/z ≈ 63.5, bare-KE distribution shape/slow shoulder, size
   correlation).
3. **RQ4 promoted to the critical path** (from "parallel non-blocking"):
   the solvated small-n shape is its direct observable; the flat Form-U
   bottom is falsified at ~2× on that branch, drag-law-robustly.
4. **RQ3 partially un-anchored:** "suppressed → bare" loses its
   experimental anchor (an evaporative bare class would arrive slow,
   ~0.2–0.3 eV — check for a slow shoulder under RQ8); the *marginal*
   self-termination-at-small-n sub-case gains weight.
5. **Program order: Wave 12 (position heterogeneity) → Wave 13 (drag
   saturation)**, designed in `TIER2_STAIRCASE_PROBE_PLAN.md`
   **Addendum H**; the position axis and the RQ4 ladder depth are the
   two candidate sources of the small-n steepness, separated by the
   (n, KE) curve (ladder moves bins, not speeds; position moves both).
6. **Prerequisite export:** the experimental (n, mean-KE) table becomes
   a small provenance-documented reference CSV under `data/reference/`
   (data-contract convention; MATLAB source of truth + exporter) before
   it anchors any calibration.

**Doc edits this entry covers:** `TIER2_STAIRCASE_PROBE_PLAN.md`
(new Addendum H — Waves 12/13 design); `RESEARCH_QUESTIONS.md` (new RQ8;
RQ4 critical-path promotion note; RQ3 un-anchoring note; RQ7
law-conditionality cross-note; coupling-map row). `CALIBRATION_MAP.md`
propagation (re-targeted Tier-2 observable; the proposed Bounded v_c
knob) is **pending** — to be folded in when Wave 13's arm is actually
built. **No code, no runs** — both waves and the new drag arm stay
behind `[PROCEED TO IMPLEMENTATION]`; nothing here discharges the F5
gate.

---

## Addendum H AMENDED — analytic feasibility pass (H.2b) + W12b dressing arm + F.5 graded-floor lever admitted (user endorsement, 2026-07-11)

Discussion while the RQ4 external many-body calculation runs (several
weeks): can the **physically-bounded lever set** reproduce the solvated
targets without waiting? User endorsement of the resulting design,
`TIER2_STAIRCASE_PROBE_PLAN.md` Addendum H amended same day:

1. **H.2b — the zero-MD analytic feasibility pass, now the gate for all
   MD legs.** The validated F.2b fate map extended on paper by the
   bounded levers (L1 position geometry — volume weight + 3–6 Å margin;
   L2 `n_eject(depth)` tied to the existing ρ̂ profile, zero free
   parameters; L3 the F.5 graded-ladder **floor variant**; L4 speed
   scale bracketed, v_c stays W13's). Pre-registered reachability
   question against (solvated 31.0/14.2 % + the 0.69 → < 0.1 eV KE
   curve); outcome shapes: (a) reachable at the floor → MD legs verify,
   (b) reachable only with the slid transition → recorded as a
   *prediction for* the RQ4 calculation, (c) unreachable → the F.5
   escape clause fires (ε/RQ2 or missing mechanism) at zero MD cost.
2. **H.3b — Wave 12b, the `n_eject(depth) = round(21·ρ̂(d_birth))`
   dressing arm** (new physics surface behind its own trigger;
   per-atom arrays already exist → no schema change; `full` vs
   `density_tied` enum with byte-identity + inert-at-center wiring
   oracles; A/B/C chain vs Wave-11 center-pinned and W12 position-only).
   Fingerprint prediction W12b-P1: its n = 1–3 fragments are the
   *fastest* in their bins — the anti-RQ4 discriminator on the (n, KE)
   curve. Retires the §4h-boundary-3 full-dressing strain by modeling
   it.
3. **F.5 floor variant admitted as knob-free** (usage note appended to
   F.5): D₀(1) = 13.3 meV (X₂/³Π) with transition-at-n=1 is fully
   [IHe05]-sourced — the external calculation arbitrates only the
   transition *location* (slid variants = the knobbed family, reported
   separately). The 2.2× solvated residual exceeds the ~1.4× cap alone;
   the taper now stacks with the geometric levers instead of carrying
   the residual.
4. **Sequencing (H.5 updated):** H.2b → W12 → W12b → W13; RQ4 returns
   as arbiter of L3; RQ8 + the (n, KE) export stay arbitration entry
   gates.

**No code, no runs** — H.2b is scratchpad zero-MD under the standing
convention but still awaits the trigger; W12/W12b/W13 each behind their
own `[PROCEED TO IMPLEMENTATION]`. Nothing here discharges the F5 gate.

---

## Addendum H.2b DESIGN FROZEN — user decisions D1–D7 + solvated mean-KE table (user discussion, 2026-07-11)

Pre-execution design discussion of the H.2b analytic feasibility pass
(three rounds), all open points settled and folded into
`TIER2_STAIRCASE_PROBE_PLAN.md` §H.2b ("design freeze" block, D1–D7):

1. **D1** — E₀ stays a sharp scalar scanned over the full RQ1 band
   [0.2, 0.5] eV (outer dimension; not pinned at the Wave-11 landing —
   bare is un-targeted, so the scan re-lands freely).
2. **D2** — the constant-E₀ convention under L2 was flagged as naive
   (user: the surface class would only feed bare). Adopted: the
   parameter-free bracket family E₀ ∝ (Σ(n_eject)/Σ(21))^p, p ∈ {0, 1},
   both reported. **Pre-derived insight recorded:** under BOTH brackets
   the surface class is suppressed (under-cooled either way — p = 0
   raises K\* while the chord K falls; p = 1 leaves K\* uniform and
   short chords still undershoot it), so the W12b-P1 fast n = 1–3
   surface feeder exists only for unsourced super-proportional E₀
   decay; a (depth, p) gate-open-at-birth diagnostic map and a per-bin
   feeder map (birth-depth/chord decomposition) are first-class
   outputs.
3. **D3** — the exposure-law bracket extends to the histogram side: the
   full fate map runs at both the current-law and ballistic speed
   profiles; outcome (c) fires only if unreachable at both (protects
   the verdict from the Wave-11 over-drag caveat).
4. **D4** — droplet priors: the stated Wave-10/11 family (Kornilov
   δ = 0.625 about ⟨N⟩ = 2000 primary; δ = 0.40/0.80; pickup-weighted).
5. **D5** — the margin band 3–6 Å is firm; the lever-interaction map
   (margin ↔ reachable n_eject / chord-K, from the code's actual
   profile) is computed first and heads the scoreboard; a near-inert L2
   inside the firm band is a finding, not grounds to loosen bounds.
6. **D6** — pre-registered reachability bar H2b-T1..T5 (solvated
   W₁ ≤ 0.5 bins; n₁/n₂ ∈ [1.75, 2.6]; n₁ ∈ [26, 36] %; the
   experimental KE curve inside the [current-law, ballistic] envelope
   with a single speed scale to ≤ ×1.5 per bin — the W13 v_c preview;
   bare ≤ 43.5 % hard / flagged > ~15 %). Reachable = T1–T4
   simultaneously at one in-bounds lever point.
7. **D7** — the user-supplied solvated mean-KE table (**mean** kinetic
   energies; provenance-pending, the H.2 export prerequisite stands):
   bare 2.9 / n₁ 0.974 / n₂ 0.545 / n₃ 0.390 / n₄ 0.341 / n₅ 0.289 /
   n₆ 0.248 / n₇ 0.201 / n₈ 0.176 / n₉ 0.154 / n₁₀ 0.138 / n₁₁ 0.122 /
   n₁₂ 0.105 eV, n ≥ 13 < 0.1 eV one-sided. **Supersedes the 0.69 eV
   n = 1 anchor** in H.2 and the 2026-07-11 two-channel entry above.
   Two reads recorded: the Wave-11 companion already matches the
   experimental n = 6 speed within ~6 % (the anomaly is a missing
   *fast small-n class*, not a uniform rescale); bare 2.9 eV sits above
   the ≈ 2.71 eV per-fragment ballistic ceiling of the 2.70 eV channel
   → channel-distinct kinematically, independent of drag (RQ8
   datapoint).

**Doc edits this entry covers:** `TIER2_STAIRCASE_PROBE_PLAN.md`
(Addendum H status blockquote second-amendment note; H.2 item-2 KE
values superseded; new §H.2b design-freeze block D1–D7 + KE table);
`RESEARCH_QUESTIONS.md` (RQ8 experimental-input NB: mean-KE table +
the kinematic-ceiling datapoint).

**No code, no runs** — H.2b execution (scratchpad, zero-MD) still
awaits `[PROCEED TO IMPLEMENTATION]`; W12/W12b/W13 unchanged behind
their own triggers. Nothing here discharges the F5 gate.

---

## H.2b analytic feasibility pass EXECUTED — outcome (c) at the frozen bar; the miss localizes to the RQ4 taper + v_c; trapped class found (2026-07-11)

Executed same day as the design freeze, under the
`[PROCEED TO IMPLEMENTATION]` trigger. Scratchpad route
`h2b_feasibility.py` — **zero MD, zero repo-code change, no artifacts
touched**: a 1D two-body chord forward model (Coulomb + solvation well +
gated pure-cubic drag; ballistic bracket = drag off) importing the
repo's ladder/density/bundle sources; 20 000 molecules → 40 000 fragment
chords per bracket, importance-reweighted over the D4 priors × D5
margins; the unified F.2b fate map (E_ej = E₀ᵢ·e^(−K); suppressed iff
E_ej > Σ(n_eject); exact ε = 0 descent; detected ≈ floor).

- **Wiring oracles all pass:** Σ(21) exact; K = 0.74460 vs the Wave-11
  MD 0.74603 (0.19 %) at production and 0.89767 vs 0.89754 (0.014 %) at
  9 Å; t_exit/v_peak/v_detect 4.60 ps/10.53/4.13 vs 4.62/10.51/4.11;
  ballistic vs analytic 4.3·10⁻⁴; E\* = 0.3955 vs 0.396; companion
  descent n = 6 exact.
- **Verdict — outcome (c) at the frozen bar:** zero of 7 308 cells pass
  T1∧T2∧T3 (hence T1–T4), at **both** exposure brackets (ballistic
  worse everywhere, best W₁ = 1.94 — the D3 protection: not
  law-conditional). Bounded-set best: the knob-free X₂ floor at
  (d080, m3, p1, E₀ = 0.22): W₁ = 0.813, n₁ = 0.223, n₁/n₂ = 1.83.
- **The miss decomposes into the two named items:** (i) the **RQ4
  ladder taper** — sliding the X₂ transition up *worsens* the ratio
  (slid2/3: 1.25/1.29; the n = 1 bin is a one-rung window, equal deep
  rungs widen n = 2/3 equally), while the beyond-bounds **RQ4-graded
  diagnostic (2.2:1.5:1.3) nearly lands** (W₁ = 0.575, n₁ = 0.267,
  ratio 1.76, bins 2–8 to ≲ 0.02; T2+T3 pass, T1 missed by 0.075) —
  the solvated histogram independently demands the RQ4 target ratios;
  the external many-body calculation is confirmed as the **blocking
  arbiter** (plateau ⇒ the F.5 escape clause fires with the geometric
  alternatives exhausted); (ii) the **W13 v_c speed scale** — the KE
  envelope holds at every top cell, the needed lift is speed-selective
  ×1.3–1.5 (energy ×2.2 → ×1.7 over n = 1 → 12) — W13-P1's premise
  pre-confirmed with numbers.
- **New structural finding — the trapped droplet-retained class
  (I44):** 6–11 % of fragments (inward partners of off-center births,
  chord K up to ≈ 18) are dissipated by the current law below the
  0.117 eV solvation barrier and never eject; excluded from all
  detected reads; ballistic-absent; t_end-conditional; MD-arbitrable.
- **Other reads:** E₀ re-lands at the solvation scale 0.22–0.27 eV
  under full geometry (I45); no gate-open-at-birth class in-bounds
  (margin floors n_eject at 13/14/15 for m = 3/4.67/6 Å) → **W12b-P1's
  fast surface feeder predicted absent** under both E₀-law brackets;
  p = 1 wins everywhere (p = 0 dies by over-suppression, as
  pre-derived); T5 satisfied at every top cell (bare 4–21 % ≤ 43.5 %).

**Docs:** findings §4j + I41–I45 + §0 table rows (11, H.2b) + header;
plan Addendum H status blockquote (EXECUTED note). Scoreboards
(`h2b_scan.csv`, `h2b_lever_interaction_map.csv`, `h2b_feeder_map.csv`)
in the session scratchpad.

**Sequencing (user's call, recorded not adjudicated):** per H.5 the
analytic pass gates the MD legs — the geometric MD legs' role shrinks
to verifying the load-bearing geometric components (near-cliff chord-K
density, the trapped class, the KE shape); the RQ4 external calculation
and the W13 v_c arm carry the two quantified residuals. Nothing here
discharges the F5 gate.

---

## Post-H.2b sequencing ENDORSED — H.5 amended; W12 re-scoped + predictions quantified; W12b parked; W13 mandate extended (user decision, 2026-07-11)

The user endorsed the post-H.2b program order; `TIER2_STAIRCASE_PROBE_
PLAN.md` amended same day (H.5 amendment block + H.3/H.3b/H.4 NBs):

1. **(n, mean-KE) export promoted** from arbitration entry gate to
   **W13 calibration prerequisite** (the D7 table → provenance-
   documented CSV under `data/reference/`, data-contract convention;
   error bars would also firm the placeholder ±20 % tolerances).
2. **W12 runs, re-scoped to verification** of the forward model's
   geometric inputs. Pre-registered predictions amended to the H.2b
   forward-model numbers (computed at the delivered physics — full
   dressing, flat ladder, current law, N = 2000, margin band {0, 4.67}
   pending the sampler audit): **W12-P5 (new, load-bearing)** trapped
   fraction 5.9–8.1 % at the 150 ps read (relaxation stage arbitrates
   permanence); K quantiles 0.15/0.49/15.2 (5/50/95 %); detected speed
   range 1.8–8.8 Å/ps (~5×, beyond the original ≥ 2×); speed
   anti-correlation −0.85…−0.92 (the K anti-correlation is weak,
   −0.12 — the fingerprint lives in the speeds); fate reads at
   E₀ = 0.28/0.38: suppressed 43–55 % / 70–77 %, flat small-n bins
   (W12-P3 retired by §4j — no steepening expected at the flat
   ladder). The margin-4.67/E₀ = 0.28 suppressed fraction 0.435
   coincides numerically with the experimental bare bin — recorded as
   a curiosity, not a target (RQ8). MD leg stays behind its own
   trigger.
3. **W12b PARKED** — its W12b-P1 fingerprint is predicted absent
   in-bounds (margin floors n_eject ≥ 13; no gate-open-at-birth class
   under p ∈ {0, 1}); revive only if W12 falsifies the geometry.
4. **W13 after W12, mandate extended** — calibrate v_c on the
   n = 1-class KE (envelope + ×1.3–1.5 speed-selective lift
   pre-confirmed, I43), then **re-run the H.2b solvated scan under the
   calibrated arm** (closes the drag-bracket-interior coverage gap;
   fixes the taper size RQ4 must deliver).
5. **RQ4 external calculation returns as the decisive arbiter** into
   the v_c-updated scan: steep bottom → histogram closes inside
   physics; plateau → the F.5 escape clause (ε/RQ2 or missing
   mechanism) fires with the geometric alternatives exhausted.

The original H.5 "outcome (c) cancels the MD legs" clause is
**superseded** (it assumed an unstructured miss). RQ8 stays the entry
gate for the bare bin. **No code, no runs** — the W12 MD leg and the
W13 arm each remain behind their own `[PROCEED TO IMPLEMENTATION]`;
nothing here discharges the F5 gate.

---

## IHe KED run-summary integration DELIVERED — Tasks 1–5 (2026-07-15)

Implementation under the `[PROCEED TO IMPLEMENTATION]` trigger, delivering the
design spec `docs/superpowers/specs/2026-07-14-ihe-ked-run-summary-design.md`
(commits `67970ae`, `c7329df`, `9697e98`, `111f8a5`, `1947fb9`). Upgrades the
run-summary comparison layer from the `vmi_summary` overlay to the `ihe_ked`
reference (a different, higher-count campaign; see the vmi_summary README's
"different-campaign" note, ⟨E⟩ ~16 % lower there).

- **Task 1 (`67970ae`)** — extracted the shared `select_final_mass_gate`
  mass-gate helper and added projected in-plane speed histograms to the
  run-summary pipeline (the like-for-like counterpart of the Abel-inverted
  2-D experimental slice).
- **Task 2 (`c7329df`)** — new `postprocess/ihe_ked.py`: reference-table
  loader for the mean-KE table (bare / n = 1…12+) with the full error model
  (systematic + statistical bands, correlated across n).
- **Task 3 (`9697e98`)** — `ihe_ked.py` trusted-curve loader (n = 0…4),
  dual-axis 2-D/3-D KE-vs-speed curve construction.
- **Task 4 (`111f8a5`)** — sim-side `fragment_mean_kinetic_energy` /
  `fragment_gate_counts` helpers, built on the shared `select_final_mass_gate`
  gate from Task 1 (one gate, not a re-derived selection).
- **Task 5 (`1947fb9`)** — `plot_run_summary.py`: new `ihe_ked_mean_energy`
  section (mean-to-mean comparison, correlated error bands, gold markers),
  `ihe_ked_curves_3d` / `ihe_ked_curves_2d` (5-panel n = 0…4 overlays with
  `v(⟨E⟩)` markers; the sim 2-D side is the projected in-plane speed from
  Task 1), and a mass-spectrum-vs-abundance grouped-bar panel against
  `integrated_i_he_abundance.csv`. Retired the `radial_velocity_with_vmi`
  overlay panel and its `VMI_REF_*` USER-SETTINGS block (superseded by the
  ihe_ked comparison layer; the `vmi_summary` loader/CSVs stay live for
  `plot_experimental_comparison.py`, which keeps its own workflow per the
  CLAUDE.md HeDFT/experimental separation rule).

**Tests:** `tests/test_ihe_ked.py` (new — reference-table loader, trusted-curve
loader, sim-side `fragment_mean_kinetic_energy` / `fragment_gate_counts`).

**Task 6 — final validation (2026-07-15).** Full suite: **2202 passed, 0
failed, 83 warnings** (0 new warning classes — all trace to the documented
§6.5 mass<->coefficient-pairing `RuntimeWarning` path or the v6→v7
checkpoint-migration `UserWarning`). Two pre-existing failures surfaced in
`tests/test_velocity_distribution.py::TestRealVmi`
(`test_real_vmi_he_loads`, `test_real_vmi_gas_loads`, asserting
`velocity_Aps.size > 100`) were diagnosed as **unrelated to Tasks 1–5** — no
diff in Tasks 1–5 touches `load_vmi_reference`, these tests, or the
`vmi_summary` CSVs. Root cause: commit `c588e77` ("New reference data as well
as updated plans for drag") intentionally re-exported
`data/reference/vmi_summary/{vmi_iplus_he,vmi_iplus_gas,vmi_iplus_he_high_snr}.csv`
together with `export_vmi_reference_data.m` and the `vmi_summary/README.md`,
switching the export to the pyabel `rIbeta()` Abel-inverted 3-D speed
distribution (README "Convention fix (2026-07)" note) — a coarser,
intentionally different velocity grid than the prior raw-pixel radial
binning: row counts confirmed at 72 / 71 / 72 data rows (73/72/73 lines incl.
header), down from >400 pre-re-export (`git show --stat c588e77` diffstat).
The stale `> 100` bound was updated to `> 50` (both new counts clear it with
margin), with an inline comment recording the c588e77/README provenance so
the bound doesn't silently drift again. No file under `data/reference/`
touched.

**Rule-2 / scope:** no schema bump, no RNG draw-order change, no
drag-law/neutral/propagation touch — a post-processing comparison-layer
delivery under the CLAUDE.md drag-model scoped exception. Next: none
scheduled — this closes the ihe_ked run-summary-integration deliverable named
in the 2026-07-14 design spec.

## Wave-13 re-scope + W13-first reorder ENDORSED — Addendum I added; W12 demoted to post-calibration verification (user decision, 2026-07-15)

Post-export direction discussion. The corrected (n, mean-KE) reference
(`data/reference/ihe_ked/`) shows the H.2b D7 target table was the *legacy
moment convention* (coherent +25–35 % correction: n = 0 2.9 → 3.706 eV,
n = 1 0.974 → 1.302 eV — far outside the 4 % ⊕ 6 % correlated bands), so the
T4/I43 premise numbers ("envelope holds; lift ×1.3–1.5") are stale, and a
back-of-envelope places the corrected n = 1 mean at the ballistic edge
(~1.2–1.31 eV zero-drag expectation vs 1.302 eV measured — the n = 1 class is
essentially undecelerated). The user re-directed: **whether any drag form can
reach the experimental mean energies is now the immediate question**; the
histogram side is blocked on the external RQ4 calculation regardless.

Decisions (user, 2026-07-15 — full detail in plan **Addendum I**):

1. **W13 is the immediate priority** (supersedes the H.5-amendment
   W12-before-W13 order). W12 is demoted to *post-calibration verification*
   under the calibrated law (its verification targets — trapped class,
   chord-K, detected speed range — are drag-law-conditional). W12b stays
   parked; RQ4 stays parallel and returns against the v_c-updated scan.
2. **H.4 scope broadened** from the single cubic→quadratic v_c arm to a
   tail-exponent family γ(v > v_c) = b·v_c²·(v/v_c)^p, p ∈ {1, 0, −1} plus
   the hard-cutoff bracket (all exactly `b·v²` in-band; Tier-0 lock and the
   in-band invariance oracle preserved; dimensional analysis in §I.3).
3. **Sweep engine (I-D2):** the extensive form × v_c scan runs in the
   validated H.2b 1D chord forward model; full MD confirms the shortlist
   only. Steps: **0** envelope re-read of the corrected reference against the
   [current-law, ballistic] brackets (zero MD; retires D7, corrects the
   stale T4/I43 numbers, generates the Step-1 pre-registered predictions) →
   **1** tail-family sweep (histogram W₁ carried as a report column) →
   **2** MD confirmation of shortlisted arm(s) behind the drag-form enum
   (own trigger; scored via the delivered ihe_ked run-summary layer) →
   **3** the retained H.4/H.5 loop (K₂.₇₀ re-measure, solvated re-scan,
   n = 1 rung check).
4. **Calibration upgrade (I-D4):** full-curve mean-to-mean fit n = 1…17
   (n = 0 excluded as RQ8 channel-distinct), per-point errors + the two
   correlated bands profiled as coherent scale nuisances; the n = 0…4 P(E)
   trusted curves are validation-only.
5. **Acceptance bar (I-D5):** ×1.25 per bin for n = 1…12, n ≥ 13 within the
   sub-0.1 eV band — a **ceiling of strictness** (loosening with documented
   model-bias justification allowed, tightening not).
6. **Predictions:** W13-P1 carries over; **I-P1** pre-registers the n = 1
   ballistic-edge read; outcome **(c′)** (experimental above ballistic beyond
   the coherent bands → no drag form suffices, W13 halts pre-build) is
   registered as a first-class outcome shape.

No code, no runs, no repo-surface change in this entry — design freeze only.
Steps 0–1 are scratchpad analytics behind `[PROCEED TO IMPLEMENTATION]`
(H.2b precedent); Step 2 is a physics-surface change behind its own trigger.
Nothing discharges F5.

## Addendum I Steps 0–1 EXECUTED — (c′) does not fire; outcome (b): the tail family delivers the scale but not the slope; KE↔histogram anti-correlation through K found; no MD build recommended (2026-07-15)

Executed same day under the trigger. Zero MD, zero repo-code change; the
recovered H.2b forward model (previous session's scratchpad, all five
artifacts intact) revalidated first — **all six wiring oracles exact**; the
new tail integrator reproduces the cached current law at v_c = ∞ to
machine precision (max |ΔK| = 4.6·10⁻¹⁴). Driver `addendum_i_steps01.py` +
`addendum_i_detail.py`, scoreboards `addendum_i_step0_envelope.csv` /
`addendum_i_sweep.csv` (10 416 cells) — scratchpad. Full read-out:
findings **§4k**; insights **I46–I48**; new **OQ-I**.

- **Step 0:** D7 formally retired (the correction is ×1.27–1.34 at n ≤ 3
  shrinking to ×1.02–1.08 at n ≥ 10 — *not* coherent). Envelope verdict:
  **(c′) does not fire** — every bin n = 1…17 inside [current-law,
  ballistic] at all five §4j cells already at f = 1 (I-P1 CONFIRMED).
  Corrected required lift: ×2.6–4.0 (n = 1), near-uniform ×1.6–2.2
  (n = 2…17); current law 0/12 vs the I-D5 bar (worst ×2.34–3.61). I43's
  "speed-selective ×1.3–1.5" superseded (note added in the register).
- **Step 1:** 28 (p, v_c) integrations × 372 fate cells: **zero I-D5 bar
  passes.** The ke_pass ≥ 8 cells are degenerate (empty n = 1 bin, E₀
  pinned at the 0.20 grid edge, f at 0.902, W₁ ≥ 1.31 — the weight-floor
  bin-evasion loophole is exposed in §4k). The genuine compromise
  (p = 0, v_c = 6, floor1, pickup, m6, E₀ 0.21) lands n = 1…6 inside
  ×1.25 incl. n = 1, but n ≥ 7 undershoots progressively to ×0.49 at
  n = 12 — **the miss is the deep-bin slope**, governed by the
  TDDFT-locked in-band law + K-selection (I46). **KE and histogram
  anti-correlate through the single exposure integral K** under the
  frozen fate map + pinned τ (I47); the hard cutoff is excluded as
  physics (I48). Verdicts: I-P2 CONFIRMED, I-P3 first branch confirmed
  with the residual re-localized (n = 1 re-feeds; deep bins fail —
  "n = 1 alone fails" refuted), I-P4 vacuous.
- **Verdict: Addendum-I outcome (b).** Step 2's entry condition (a
  bar-passing cell) is **not met — no drag-arm MD build recommended from
  this sweep.** Conditional shortlist if OQ-I resolves toward joint
  recalibration: p = 0, v_c ≈ 6 (alternate p = −1, v_c ≈ 7–8.5). The
  decision surface is OQ-I: (a) joint (v_c, τ) calibration (τ is a
  Bounded probe pin), (b) fate-map-level change (ε/RQ2, E_ej form — the
  F.5 clause), (c) the deep-bin side may implicate the in-band law
  (Tier-0-locked; would need new TDDFT input). **User adjudicates;
  nothing discharges F5; W12 remains demoted per I-D1; RQ4 unaffected
  and still parallel.**

## Addendum I Step 1b EXECUTED — joint (v_c, τ) closure found at τ = 4.1 ps: first full two-observable landing (T1∧T2∧T3 ∧ 12/12 KE bar) on the rq4graded ladder; Step-2 build re-opens pending three user adjudications (2026-07-15)

User approved OQ-I arm (a) same day; design + predictions J-P1–J-P3
registered in plan **§I.8** pre-execution. Zero MD: τ enters the dynamics
nowhere, so the τ axis is an exact fate-map rescale
(K(τ) = K(6.55)·6.55/τ) over 13 cached tail integrations — 28 899 cells
(`addendum_i_joint.csv`, scratchpad, + `addendum_i_joint.py` /
`addendum_i_joint_detail.py`). Full read-out: findings §4k Step-1b;
insights **I49–I50**; OQ-I arm (a) marked answered.

- **18 contiguous joint closures** (12/12 KE bins within ×1.25 AND
  W₁ ≤ 0.85), all at **τ = 4.1 ps**, all **rq4graded**, spanning
  (p = 0, v_c = 6–7) + (p = −1, v_c = 7–8.5), both priors, all margins,
  E₀ = 0.22–0.26 (interior — the Step-1 edge-pinning resolves).
- **Best cell passes everything:** (p = −1, v_c = 7, τ = 4.1, pickup,
  m3, E₀ = 0.25): W₁ = 0.272, n₁ = 0.280 (T3 ✓), n₁/n₂ = 1.861 (T2 ✓),
  T1 ✓, worst KE ×1.243, deep bins ≥ ×0.80⁻¹, bare 0.149,
  trapped 0.054 — where H.2b found in-bounds unreachability (best
  W₁ 0.575), opening the τ pin closes both observables on the
  RQ4-anticipated taper.
- **Both knobs load-bearing** (v_c = 15 control never joint-closes at
  any τ; Step 1 showed v_c alone fails) — the I47 anti-correlation
  resolves by the pair. The Step-1 deep-bin slope was K-selection, not
  in-band physics (J-P2 refuted in its conservative branch — the
  informative direction). J-P1 and J-P3 CONFIRMED (matched-clock
  scaling; control clean).
- **Form discrimination (I50):** p = −1 wins T2/T3, sags at n ≥ 13;
  p = 0 gives the flattest KE curve of the program (×0.85–1.04 over
  n = 1…17), tops at ratio 1.66. Discriminators: the n ≥ 13 KE tail +
  n₁/n₂.
- **Verdict: Step-1b outcome (a) — the W13 Step-2 MD build re-opens**
  with a two-knob (v_c, τ) arm and target region (τ ≈ 4.1;
  p ∈ {0, −1}; v_c ≈ 6–8.5; E₀ ≈ 0.23–0.26; rq4graded). Three user
  adjudications before any build: (i) the **τ re-classification event**
  (pre-registered in §I.8): 4.1 ps = 0.63 × the GAH25-sourced pin —
  re-argue or reclassify (CALIBRATION_MAP propagation at build);
  (ii) **RQ4-conditionality** — the closure lives on the anticipated
  taper; stakes sharpened (a plateau now breaks both observables);
  (iii) **form choice** (implement both tails for discrimination, or
  one). Caveats: deep-bin passes hug the bar edge (×0.80–0.91); T2
  margin 6 %; τ-grid quantized at 4.1; 1D-model boundaries carry.
  Nothing discharges F5.

## Addendum I Step 1c EXECUTED — closure basin is a plateau (225 joint / 24 full-house cells); refined targets (p=−1: v_c 7.5, τ 3.8 | p=0: v_c 6.0–6.5, τ 3.8–4.0); K-P3 refuted: floor1 joint-closes, the taper buys only the full histogram bar; RQ9 opened & parked (2026-07-15)

User approved the basin scan (#2) and parked the τ-sourcing research as
**RQ9** (`RESEARCH_QUESTIONS.md`, new entry + coupling-map row). Design +
K-P1–K-P3 pre-registered in plan **§I.9**. Zero MD: 8 new tail
integrations (fine v_c grid) + free τ-rescale = 39 312 cells
(`addendum_i_basin.csv`, scratchpad, `addendum_i_basin.py`). Read-out:
findings §4k Step-1c; insight **I51**.

- **K-P1 CONFIRMED — plateau:** joint basins p = 0:
  v_c [6.0, 7.5] × τ [3.2, 4.8]; p = −1: v_c [7.0, 8.5] × τ [3.0, 4.8]
  (≈ 1.5 Å/ps × 1.6 ps each). **24 full-house cells** (joint KE ∧
  T1∧T2∧T3), τ ∈ [3.4, 4.2], all rq4graded; best
  (p = −1, v_c = 7.5, τ = 3.8, d080, m3, E₀ 0.25): W₁ = 0.196,
  n₁ = 0.289, ratio 1.861, 12/12 KE worst ×1.225. Refined targets
  supersede the Step-1b grid point (τ 4.1 → 3.8–4.0 = ×0.58–0.61 of
  the GAH25 pin; RQ9 numbers updated).
- **K-P2 CONFIRMED:** p = 1 zero closures — form family truncates to
  {0, −1}.
- **K-P3 REFUTED (finding):** floor1 holds 16 joint closures (12/12 KE,
  W₁ ≥ 0.575, n₁ ≤ 0.23, zero full houses) — **bounded physics +
  (v_c, τ) reaches the full KE curve + an H.2b-best-level histogram**;
  the RQ4 taper specifically buys n₁ ≥ 0.26 / ratio ≥ 1.75 /
  W₁ 0.58 → 0.20. The Step-1b "a plateau breaks both observables"
  stake is walked back: an RQ4 plateau breaks T2/T3 only (as in H.2b);
  the KE side stands either way (I51).
- Step-2 build decision surface unchanged otherwise; refined targets
  ready. Nothing discharges F5.

## Addendum I Step 2 (MD confirmation build) DESIGN FROZEN — §I.10 added; decisions S2-D1–S2-D4; two build surfaces identified (capped_cubic arm + tabulated-ladder config); 4-config pilot matrix; awaits its own trigger (2026-07-16)

Doc state through Step 1c committed (`2ac509e`). Surface audit via
repo exploration: τ (`internal_energy_cooling_tau_ps`), f_int, s_eff,
cooling gate, biphasic mechanism, relaxation/detection stages,
off-center birth sampling, and the ihe_ked scoring layer all exist;
missing exactly (1) a velocity-capped drag arm (no saturation form is
realized in `physics/drag.py`; `THRESHOLD` is reserved-unimplemented)
and (2) a config path to the existing `TabulatedLadder` (no field wires
`rungs_eV`; neither floor1 nor rq4graded is expressible today).

User decisions (2026-07-16): **S2-D1** build the tabulated-ladder
config surface (receives RQ4 verbatim later); **S2-D2** N = 50 pilot →
N = 500 winner; **S2-D3** focused 4-config matrix (2 targets +
current-law control + floor1 bounded-physics leg); **S2-D4** fixed
N = 2000 droplets, distribution axis deferred.

Design: plan **§I.10** — Slice T1 `capped_cubic` {b, v_c, p_tail}
(dimensional analysis in-plan; v_c = ∞ byte-identity with
`linear_cubic`; in-band invariance oracle; p_tail restricted to
{0, −1} at config load), Slice T2 `tabulated_ladder_rungs_eV` (rung
tables built generator-side), Slice T3 pilot matrix C1–C4 at
production kinematics (E_int(0) = E₀ via partition fraction; s_eff = 8;
density_scaled cooling; position sampling on; margin-0 caveat
recorded), Slice T4 ihe_ked scoring + winner N = 500. Pre-registered
predictions S2-P1–S2-P4 (control anchor; pilot KE ×1.4; the I51
floor1/rq4graded histogram split; N = 500 bar verdict).
CALIBRATION_MAP propagation lands with Slice T1 (v_c/p_tail
Bounded→Derived; τ per RQ9). **Execution awaits a fresh
`[PROCEED TO IMPLEMENTATION]`.** Nothing discharges F5.

## Slice T1 DELIVERED — `capped_cubic` drag arm; both §I.10 oracles pass incl. byte-identity of a delivered probe dir at v_c = 5.3; CALIBRATION_MAP rows 4/4b/11 propagated (2026-07-16)

Executed under a fresh `[PROCEED TO IMPLEMENTATION]` (user: Slice T1
only; the six pre-build clarifications — T2 scope through all three
stage guards, v_c/p_tail as *coefficients* in the bundle, a new
config-tag scheme for C1–C4, the 1000 ps relaxation cap, one shared
pilot seed, slice-by-slice TDD — all confirmed as proposed). Tests
written first and watched fail (ImportError on `CAPPED_CUBIC`, the
missing-feature failure), then the arm implemented.

**Code (2 files):**

- `i2_helium_md/physics/drag.py` — `CAPPED_CUBIC = "capped_cubic"`,
  coefficient keys `("b", "v_c", "p_tail")` [amu·ps/Å², Å/ps,
  dimensionless], added to `REALIZED_FORMS` (the ion driver's scope
  guard and the loader accept it automatically via the single-source
  tables). Branches in `drag_force` / `drag_gamma`: in-band
  `g·b·v³` / `g·b·v²` — **the same arithmetic expression as
  `linear_cubic(a=0)`**, so in-band byte-identity holds bitwise, not
  approximately; tail `g·b·v_c²·v·(v/v_c)^p` / `g·b·v_c²·(v/v_c)^p`,
  continuous at v_c. Implementation detail with physics weight: the
  no-tail case (incl. v_c = ∞) returns the pure-cubic expression
  early, and the tail path substitutes `max(v, v_c)` in the discarded
  branch of `np.where`, so `p_tail = −1` never evaluates `0**−1` at
  rest and v_c = ∞ never produces `inf/inf` — verified warning-free
  under `np.errstate(all="raise")`. Module stays mass-agnostic (no
  `m` anywhere); the FDT amplitude reads the same closed-form γ
  (Tier 3 stays inert).
- `i2_helium_md/config.py` — `capped_cubic` added to the `DragForm`
  Literal + `_KNOWN_DRAG_FORMS`; `check_drag_config` dissipativity arm:
  `b > 0`, `v_c > 0` (v_c = ∞ admissible: the byte-identity limit),
  and **`p_tail ∈ {0, −1}` enforced at config load** (the
  Step-1c-surviving set; any other exponent is refused with a pointer
  to §I.10 — a new adjudication, not a config value).

**Tests (3 files, RED→GREEN):** `tests/test_drag.py` new
`TestCappedCubic` (closed-form pins per branch and per tail;
continuity at v_c; γ ≡ F/v identity; dissipativity + shared gate
factor; γ(0) = 0; warning-free mixed-speed/∞-cap evaluation; **exact
`==` byte-identity** vs `linear_cubic(a=0)` at v_c = ∞ and at finite
v_c ≥ v_max, plus in-band-only byte-identity below a finite in-range
cap; missing-key rejection). `tests/test_drag_config.py` new
`TestCappedCubicGuard` (both tails pass; ∞ cap passes; b ≤ 0, v_c ≤ 0,
p_tail ∉ {0, −1} refused), a `capped_cubic` loader round-trip in
`TestLoaderFormGeneric`, and the `DragForm` enum-completeness update.
`tests/test_ion_drag_smoke.py`: both tails added to
`_FORM_PHASE_PARAMS` with the cap placed *inside* the smoke speed
range (tail branch exercised end-to-end through the BAOAB driver), and
`test_capped_cubic_above_vmax_run_is_byte_identical` — a tiny-N driver
A/B asserting every trajectory/energy array equals the preset
pure-cubic run bitwise. Suites: the three touched files 174 passed;
`test_tier0_drag_comparison.py` 4 passed (Tier-0 18 Å gate untouched);
**full suite 2250 passed** (warnings = the known §6.6 override
RuntimeWarnings).

**§I.10 oracles:**

1. **v_c = ∞ / v_c ≥ v_max byte-identity with `linear_cubic(a=0)`** —
   pytest, exact `==` (see above).
2. **In-band invariance on a delivered 0.80 eV probe dir** — scratchpad
   driver (`slice_t1_inband_oracle.py`, zero repo change, run dir
   untouched): the Phase-D bridge dir
   (`…_tier2probe_b080_mix_k1.00_l0.90_fi0.50_fr0.10_tau6.55`) re-run
   ion-stage from its own `neutral.npz` under `capped_cubic` at
   **v_c = 5.3** for **both** tails. Measured in-window max speed
   **5.2306 Å/ps < 5.3** (the §I.10 premise, confirmed on the
   artifact); every ion-checkpoint array **byte-identical** (the sole
   nominal mismatch was the all-NaN `temperature_diagnostic` fill
   compared without `equal_nan` — a comparison artifact, fixed in the
   oracle, not a physics difference).

**CALIBRATION_MAP propagated:** row 4 rewritten — $v_c$ realizes the
former contingent $v_\text{ceiling}$ (R10-(b)), class **Bounded**
(band 5.3–15, Step-1c targets 7.5 / 6.0–6.5) → Derived at the T4
winner; new **row 4b** $p_\text{tail} \in \{0, −1\}$ (Free choice,
2 arms; p = 1 excluded by K-P2, hard cutoff excluded as physics I48);
row 11 τ carries the **RQ9 reclassification-pending** note (joint
closure at τ ≈ 3.0–4.8 = ×0.6 of the GAH25 pin; pilot at 3.8–4.4);
tally updated.

**Boundaries:** no preset or production config selects `capped_cubic`
(the pilot generator is Slice T3's job); Slices T2–T4 not started;
nothing discharges F5.

## Slice T2 DELIVERED — tabulated-ladder config surface wired through all three stages via a single resolve_ladder bridge; tabulated≡form_u bit-identity oracles pass at step/relaxation/detection level (2026-07-16)

Executed under the user's `[PROCEED TO IMPLEMENTATION]` (T1 committed
first, `eb8398e`; the three pre-T2 discussion points all confirmed as
proposed: (1) resolver-injection wiring shape, (2) config floor ≥ n\*
with loud runtime out-of-table failure + N_STAR + 11 generator
convention, (3) form_u+table refused loudly / picture-κ inert-but-
reportable under tabulated). TDD: the full RED set watched fail
(ImportError on `resolve_ladder`, the missing-feature failure; unknown
field TypeError; DID-NOT-RAISE guards) before any implementation.

**Code (10 files):**

- `physics/dissociation_ladder.py` — `resolve_ladder(dissociation_
  ladder, rungs_eV)`: the single cfg→ladder bridge (form_u → None;
  tabulated → validated `TabulatedLadder`; unknown → reject). Validation
  single-sourced here: tabulated-without-table refused (no silent Form-U
  fallback), table-under-form_u refused (stale-intent hazard), floor
  ≥ N_STAR positive finite rungs (Σ(n\*) must be table-covered).
  `d0_of_n` / `ladder_cumsum` / `gate_threshold` gained a byte-inert
  `ladder=None` kwarg: when injected, the table overrides the Form-U
  parametrisation entirely (picture/κ ignored — documented).
- `physics/solvation_cooling.py` — `ladder=` threaded through the K
  split (`s_collective_eV`, `e_infinity_eV`, `e_bind_pair_eV`,
  `e_electrostriction_eV`) and `newton_cool_step`.
- `physics/internal_energy_budget.py` — threaded through
  `dE_int_pickup_eV`, `pickup_bath_release_eV`, `dE_int_shed_eV`,
  `f_int_floor`, `reconstruct_e_int_eV`.
- `physics/evaporation.py` — threaded through `_gate_threshold_eV`,
  `gate_margin_eV`, `is_self_bound`, `rrk_rate`, `evaporation_step`,
  `evaporation_step_components`.
- `physics/pickup.py` — threaded through `pickup_step` /
  `pickup_step_components` (the S1 heat off the injected table).
- `config.py` — new field `tabulated_ladder_rungs_eV:
  Optional[tuple[float, ...]] = None` (inert default; physics-live at
  this slice — no rule-2 carry). `check_ladder_config` arm 5 delegates
  the pairing/table validation to `resolve_ladder` (rule-1 single
  source with the stages' point-of-use calls).
- `simulation/ion_propagation_step.py` — `biphasic_step` resolves once
  and threads all five energetics calls (K2 cooling asymptote+step,
  evaporation, pickup, bath release). The 2026-07-02 NotImplementedError
  refusal is superseded by the missing-table ValueError.
- `simulation/ion.py` + `simulation/ion_initial_state.py` — the
  per-step and t0 `e_bind_pair` E_pot folds ride the same resolved
  ladder (5-term invariant stays closed under a tabulated run).
- `simulation/relaxation_stage.py` — **the silent-Form-U hazard closed**:
  the stage had no ladder guard pre-T2; it now resolves once and
  threads the freeze mask + both translation-arm folds.
- `simulation/detection_stage.py` — resolves once; threads `rrk_rate`,
  the P1–P3 guard's `d0_of_n`, `_permanent_reason` (d0 + gate margin),
  and the per-fire K1 drain. DS refusal superseded like biphasic_step.
- `simulation/run_directory.py` — `load_cfg` coerces the JSON array
  back to the declared tuple (cfg == loaded holds); invalid payloads
  raise; a pre-T2 `cfg.json` without the key loads with the default
  (back-compat, the Slice-DS acceptance-criterion precedent).

**Tests (RED→GREEN; +43, full suite 2250 → 2293 passed):**
`test_dissociation_ladder.py` (resolver semantics; injection dispatch;
form_u-fed-table bitwise identity on d0/Σ), `test_ladder_config.py`
(guard surface incl. off-diagonal + floor + positivity + inert
default), `test_evaporation.py` / `test_solvation_cooling.py` /
`test_internal_energy_budget.py` (per-module fed-table bit-identity +
distinct-table liveness; three pre-existing monkeypatch stubs updated
to the new signature), `test_biphasic_step.py` (step-level byte
identity with cooling+evaporation+pickup live; short-table loud
pickup-lookup failure — the runtime range convention made concrete;
guard test flipped to the data-path contract),
`test_relaxation_stage.py` (stage byte identity + fail-loud),
`test_detection_stage.py` (detected-read byte identity + fail-loud),
`test_run_directory.py` (tuple round-trip; pre-T2 back-compat; invalid
payload).

**Equivalence argument (why bitwise, not approx):** the Form-U
`ladder_cumsum` prefix table and a `TabulatedLadder` fed
`d0_of_n(1..32)` build the same float array through the same
`np.cumsum`, and the tabulated `d0_of_n` is an exact table lookup of
the same values — so every downstream expression is the identical
arithmetic, and the biphasic two-draw RNG stream is consumed
identically (draws are unconditional). Byte-identity therefore holds
end-to-end, not statistically.

**Boundaries:** no preset/generator selects `tabulated` (the rq4graded
/ floor1 rung tables are Slice T3's generator-side job); the
postprocess `derived_diagnostics` `e_infinity_eV(0, ...)` call is
ladder-independent (Σ(0)/Σ(n\*) = 0 for any table) and stays unwired
*(superseded by the post-delivery review below: the whole
reconstruction layer is now ladder-threaded — the `is_self_bound` call
in `crossing_time_ps` was ladder-dependent and not covered by this
exemption)*; CALIBRATION_MAP unchanged (T2 adds no calibration knob —
rung tables are Derived, generator-side); Slices T3–T4 not started;
nothing discharges F5.


## Slices T1+T2 post-delivery code review EXECUTED — fixes applied; the n=0 bare-ion crash on the tabulated path closed (2026-07-16)

High-effort multi-agent review of `ab723e1..597434e` (Slice T1
`capped_cubic` + Slice T2 tabulated-ladder surface): 4 finders /
10 independent verifiers, 15 candidates verified, 0 refuted, merged to
9 findings. **T1 came out clean** (one cleanup finding); every
substantive finding was in the T2 wiring. All 9 findings fixed under
the user's `[PROCEED TO IMPLEMENTATION]`, test-first (RED watched
before each behavioral fix). Full suite **2304 passed** (2293 + 11 new
tests).

**F1+F2 (CONFIRMED, blocking): bare-ion `n = 0` crash under
`tabulated`.** `rrk_rate` and the `evaporation_step_components`
`dE_int_shed_eV` precompute evaluated `d0_of_n` on unclamped occupancy
arrays; the `n = 0` lanes are discarded by the existing masks (`k = 0`
for `n < 1`; a bare ion can never fire), but `TabulatedLadder.d0_of_n`
raises where the Form-U sigmoid evaluated harmlessly. Since `n = 0` is
reached via the `n = 1` direct fire (the bare-I⁺ peak path!), every
production tabulated run would have died mid-run — independently
reproduced before fixing (scalar, vectorized ensemble, and the
detection cascade all raise). Fix: `np.maximum(n_arr, 1)` at the two
lookup sites, mirroring the adjacent
`effective_dof(np.maximum(n_arr, 2))` precedent; Form-U output is
bitwise unchanged (clamped lanes are mask-discarded). New tests:
scalar/vector `rrk_rate` at `n = 0`, the components bare-ion lane, and
a detection cascade-to-bare equivalence oracle (`n = 1` seed,
`E_int > D_0(1)`, terminal `n = 0` = "frozen").

**F3 (CONFIRMED): reconstruction layer scored tabulated runs against
the dead Form-U ladder.** The pre-T2 `NotImplementedError` refusals had
protected the postprocess layer only vacuously; T2 wired the stages but
not `crossing_time_ps` / `reconstruct_diagnostics` /
`tier2_bridge_report.py`. Fix: `crossing_time_ps` gained a
`ladder=None` kwarg (threaded to `is_self_bound`);
`reconstruct_diagnostics` and the bridge report resolve via
`resolve_ladder(cfg…)` and thread it (incl. `e_infinity_eV(0)` and the
report's `Σ(21)` closed form, for run-exact consistency). Tests:
doubled-table liveness + form_u-fed-table bit-identity at both the
helper and `reconstruct_diagnostics` level.

**F4 (CONFIRMED): `cfg.electrostriction_binding_eV` reported the
Form-U value under a tabulated config.** The derived property now
resolves the config's ladder. Tests: fed-table identity +
distinct-table liveness.

**F8 (PLAUSIBLE, fixed): positive-finite rung validation moved into
`tabulated_ladder`** (a `D_0 ≤ 0` rung opening the RRK bracket at zero
cost is a property of the ladder object, not the config path);
`resolve_ladder` keeps only the pairing + ≥ n\* floor. Direct
construction (tests, the T3 generator) now fails as loudly as the
config path.

**F6 (hot-path): per-step resolve + uncached table math.**
`resolve_ladder` is now memoised on the normalised
`(selector, rungs-tuple)` key (`lru_cache`; exceptions are not cached,
so fail-loud is unchanged; `biphasic_step`'s per-step point-of-use
resolve returns the same immutable instance), and `TabulatedLadder`
caches its rungs array + Σ prefix as read-only `cached_property`
arrays — the Form-U `_sigma_prefix_table` twin. Values bit-identical.

**F5 (T1 cleanup): the duplicated `capped_cubic` tail scaffold**
(mask / fast path / `max(v, v_c)` 0**−1 guard) now lives once in
`drag._capped_cubic_tail_factor`, consumed by both `drag_force` and
`drag_gamma`; each branch keeps its own final multiplication order so
the delivered bitwise tail pins and the in-band byte-identity oracle
are untouched (all 178 drag tests green, Tier-0 gate included).

**F7 (test dedup, the Slice-G precedent):** the 9 per-file Form-U
rung-table builders collapsed into `tests/ladder_feeds.py`
(`form_u_rungs` / `form_u_rungs_for` / `form_u_ladder`); per-file
helpers are now one-line delegates or direct imports.

**F9 (rule 5):** the `ladder` kwarg + the picture/κ-are-ignored
override semantics documented in every public consumer docstring
(solvation_cooling ×5, internal_energy_budget ×5, pickup ×2,
evaporation `gate_margin_eV` / `is_self_bound` / `evaporation_step`;
the `*_components` twins inherit via their "remaining kwargs" clauses).

**Not changed:** no preset selects `tabulated` (T3's job); checkpoint
schema, RNG draw order, constants untouched; the runtime
out-of-table-above lookup stays the accepted fail-loud convention
(only the physically-reachable `n = 0` *below*-table lane was
clamped — it is mask-discarded, not physics).


## Slice T3 DELIVERED — MD-confirmation generator built + C1–C4 pilots EXECUTED at production kinematics (first in-repo off-9 Å MD); dt-halving oracle clean (2026-07-16)

Executed under the user's `[PROCEED TO IMPLEMENTATION]` ("as proposed
for all" on the five pre-T3 discussion points: (1) C-label tag scheme
`tier2probe_conf270_cN` in the probe namespace; (2) ladder-table
construction pins — base mixture κ=1 Form-U via `d0_of_n(1..32)`,
rq4graded rungs 1–3 × (2.2, 1.5, 1.3), floor1 rung 1 := 13.3 meV
absolute; (3) probe pins carried (λ₀ = 0.9/ps, f_ret = 0.1, mixture/κ=1
stamped, 30 ps ion window, Tier-0 dt, exact-quotient f_int) + a
one-config dt-halving spot check; (4) detection stage inline in
`_run_one`; (5) one trigger = build + tests + the 4 × N = 50 pilots).
TDD: both RED sets watched fail first (ImportError on the missing
names, the missing-feature failure).

**Code (2 files):**

- `scripts/tier2_common.py` — `build_biphasic_cfg` gained eight
  None-sentinel byte-inert kwargs (explicit-None call proven equal to
  the bridge call): `drag_form` + `drag_coefficient_overrides` (form
  swap; the new form's coefficient set assembles override-first /
  bundle-second, so `capped_cubic` **structurally inherits the locked
  pure-cubic b** — the T1 in-band byte-identity premise; provenance +
  stamped binding ride unchanged; overrides-without-form, unknown keys,
  and keys available from neither source are refused loudly),
  `dissociation_ladder` + `tabulated_ladder_rungs_eV` (tuple-coerced;
  pairing validated by `resolve_ladder` at `validate()`),
  `R0_GS_angstrom`, `E_coulomb_scale` (stamped *explicitly* — the
  droplet-distribution preset precedent is 0.8, never inherit
  silently), `single_initial_position`, and `detection_time_ps`
  (enables the stage, mirroring the `relaxation_time_ps` pattern). New
  `tier2_confirmation_run_tag` / `tier2_confirmation_run_dir_name`: the
  C-label **is** the run identity; all knobs are read from the
  authoritative `cfg.json` (F3 convention) — at production f_int is an
  exact quotient (0.25/2.70 = 0.0926…, 0.24/2.70 = 0.0889…) and the
  probe tag's two-decimal f_int would alias C1 with C2. The `conf270`
  digits satisfy the Wave-8 "2.70 eV dirs must be tag-distinct" note;
  namespace-locked in both directions (probe glob sweeps it, campaign
  glob and delivered `tier2probe_b…` tags never match).
- `scripts/gen_tier2_md_confirmation.py` (new) — the frozen §I.10
  C1–C4 `ConfirmationSpec` matrix; rq4graded / floor1 rung tables
  **constructed generator-side** from `d0_of_n` × the taper multipliers
  (no taper physics in the package — the T2 boundary; `form_u` key =
  ride-the-default, byte-inert for C3); production-budget guard
  (refuses ≠ 2.70 eV, mirroring the probe generator's 0.80-only guard
  in the opposite direction); `_run_one` runs **all four stages**
  (neutral → ion → E2 relaxation → detection — the first generator to
  run detection inline); five-artifact resume guard
  (`detection.npz` required for skip-complete).

**Tests (+24, RED→GREEN; full suite 2293+11 review = 2304 → 2328
passed):** `test_tier2_common.py` — byte-inert default lock, capped
form-swap with bundle-b inheritance + stamped-binding identity, the
three loud-refusal arms, ladder/production-kinematics/detection stamps,
conf-tag encoding + label validation + namespace lock.
`test_gen_tier2_md_confirmation.py` (new) — Form-U base table ==
`d0_of_n(1..32)` with the **Σ(21) = 0.18783720 eV** wiring pin,
rq4graded/floor1 construction oracles, ladder-key dispatch, the frozen
C1–C4 values, common-pin + per-config stamps (f_int exact, never
rounded), budget guard, unique conf dirs, F3-discovery ignores a
complete conf dir, tiny-N end-to-end through all four stages writing
all five artifacts, partial/skip guards, `main()` schedule.

**Pilots EXECUTED** (N = 50, shared bridge seed 20260604, ion 3000
steps, relaxation full 1000 ps cap, detection at 8.53 µs; dirs
`data/runs/9A_drag_shared_pure_cubic_N50_tier2probe_conf270_c{1..4}`):

| config | drag | τ [ps] | ladder | E₀ [eV] | n_relax_mean | n_detect_mean |
|---|---|---|---|---|---|---|
| C1 | p = −1, v_c = 7.5 | 3.8 | rq4graded | 0.25 | 8.94 | 8.63 |
| C2 | p = 0, v_c = 6.5 | 4.0 | rq4graded | 0.24 | 9.11 | 8.77 |
| C3 | linear_cubic (control) | 6.55 | form_u | 0.25 | 7.38 | 6.97 |
| C4 | p = −1, v_c = 7.5 | 4.4 | floor1 | 0.23 | 8.78 | 8.39 |

**dt-halving spot check** (scratchpad driver, C1, dt 0.01 → 0.005,
ensemble-mean read — the two arms are independent RNG realizations
since halving dt doubles the pickup draw count): n̄_detect
8.630 ± 0.095 → 8.600 ± 0.099 (drift 0.30× SEM); mean detected KE
0.2966 ± 0.0037 → 0.2971 ± 0.0039 eV (drift 0.14× SEM); C1 arrival
state 91 % frozen / 0 % suppressed / 9 % time_exhausted. The Tier-0
dt = 0.01 ps is adequate at production kinematics at the observable
level — the §I.10 "never exercised in-repo above 5.3 Å/ps" flag is
discharged.

**Boundaries / notes for T4:** scoring and the S2-P1–P4 verdicts are
**Slice T4's job** — the means above are generator stdout, not an
adjudication. The delivered staircase-probe report's `*_tier2probe_*`
glob **would sweep the conf dirs** (they are artifact-complete for it)
and would score them against the 0.80 eV anchored comparator — T4
should exclude `*_tier2probe_conf*` from that report or ignore its rows
for these dirs. The margin-0 off-center-birth caveat stands (1D used
3–6 Å). The stale-artifact policy applies to the four pilot dirs.
Nothing discharges F5.


## T3 pilot first-look FINDING + Step 2c geometry-closure slices DESIGNED — §I.11 (T5–T9) added; §I.10 T4 superseded/absorbed (2026-07-16)

First-look read of the four conf detection artifacts (findings **§4l**,
insights **I52–I54**): every config parks in a narrow n ≈ 7–9 cluster,
zero weight below n = 5, bare/n₁ identically zero — **S2-P3 is
structurally silent**. Diagnosis: the Step-1c closure cells were scored
on the full H.2b ensemble (L1 margin, **L2 dressing
n_eject = round(21·ρ̂)**, D2 E₀-law p = 1, D4 priors), while the repo
biphasic seed starts every ion at the full n₀ = 21 shell
(`ANCHOR_N_START` — the Tier-1a validation convention, I54). The §I.10
margin-0 caveat is measured to be load-bearing; the H.3b park premise
is overturned via its own revival clause. What *does* transfer (I53):
the F.2b τ/E₀ race ordering, the capped-tail KE lift (×1.3–2.0 vs
current-law ×0.85–1.15 in populated bins), clean four-stage execution,
dt adequacy.

Consequence (user direction 2026-07-16: "the next slices should be
reproducing what was achieved in the scratchpad"): plan **§I.11 Step
2c** added — slices **T5** (initial-shell dressing arm,
`initial_shell_model ∈ {full, density_tied}` — H.3b revived verbatim),
**T6** (E₀–dressing coupling,
`internal_energy_partition_law ∈ {constant, sigma_proportional}`),
**T7** (hard birth-margin knob — *conditional* on V0-1), **T8**
(droplet-prior family, Kornilov δ + pickup-weighted; supersedes the
S2-D4 deferral), **T9** (staged A/B/C/D oracle chain — each leg A/B'd
against the 1D twin re-scored at exactly that MD configuration — then
the re-centered confirmation re-pilot + the scoring machinery + the
N = 500 winner gate absorbed from §I.10 T4, which is superseded as
written). Pre-slice verification **V0**: (1) the repo birth sampler is
Boltzmann `r²·e^(−U/kT)`, not the 1D bare-r²+margin — reconciliation
decision; (2) the droplet well exerts **no force** on the drag-path ion
(E_pot bookkeeping only) → the §4j trapped class cannot exist in MD —
normalization convention to fix **[REFUTED 2026-07-16 — the force IS
applied; see the V0-2 correction entry below]**; (3) `h2b_feasibility.py` survives only
in expired scratchpads — recover + commit as the reproducible twin.
Pre-registered S2c-P1..P4 (numerics filled from twin re-scores per
leg). **All slices behind their own `[PROCEED TO IMPLEMENTATION]`;
V0-1/V0-3/T8-scope/T9-re-centering are open user decisions.** Nothing
discharges F5.

## V0-2 premise REFUTED — the droplet-well force IS applied on the drag path (Tier-0/1a untouched); T3 trapped-class read EXECUTED: 0/100 bound in all four pilots (2026-07-16)

User challenge ("this would invalidate the whole Tier-0/Tier-1a
trajectory matching — did we really miss this?") prompted a code-level
verification of the V0-2 premise recorded in the entry above. The
premise is **false**; the entry above is flagged in place and plan
§I.11.V0 item 2 is rewritten (this entry is the delivery record).

**The trace.** `physics/leapfrog.py :: _ion_accel_fn` sums the droplet
solvation force (`_droplet_acceleration(use_ion_binding=True)`) with
the ion partner Coulomb; `make_ion_accel_fn` wraps it and is by its
own docstring "the single source of the ion `acc_fn` (Coulomb +
droplet)", consumed by the BAOAB B/A kicks. `simulation/ion.py`
rebuilds exactly that `acc_fn` every step on the drag path for **all**
mass scenarios (`fixed`, `anchored_discrete`, `biphasic`); the E2
relaxation `"coulomb"` mode binds the same `acc_fn`. Numeric check: a
surface ion with its Coulomb partner at 500 Å feels
a_x = −0.5694 Å/ps² vs the droplet-only analytic −0.5667 (residual
exactly the 1.44/500² Coulomb tail).

**The barrier is the Tier-0 one, on both sides.** The well the ion
feels is the Method-B jointly-fitted stamp
`effective_binding_energy_I_ion_eV = 0.11675778 eV`
(`shared_pure_cubic` bundle; §6.5.1 exact-pairing guard; presets copy
it into `binding_energy_I_ion_eV`; confirmed in the T3 `cfg.json`s).
The config-default 0.3 eV is inert in real runs. The 1D twin's §4j
trapped class uses the same ≈ 0.117 eV — twin and MD already agree
(user adjudication 2026-07-16: 0.117 eV is the Tier-0 value; the
barrier-reconciliation rule floated in discussion is **dropped**).
Consequence: **Tier-0/1a trajectory matching is untouched** — the drag
coefficients were extracted and validated *with* the well force on,
which is precisely what the §6.5.1 stamp exists to protect.

**T3 trapped-class read (zero MD, executed under the user's go).**
Escape criterion KE + U(r) ≥ E_b per ion on the stored final states
(complex mass, ion binding, steepness 14.2 Å): **0/100 ions bound in
every config C1–C4**, at both the 30 ps ion-end and 1000 ps
relaxation-end reads. All ions are outside the droplet already at
30 ps (min radius 120–147 Å); minimum escape margin +0.167 eV (the
slowest ion, current-law C3), all others ≥ +0.212 eV. The trapped
droplet-retained class is therefore *expressible* in the MD (force
live, same barrier) but *unpopulated* at the undressed T3
configuration — a geometry-conditional read, to be re-taken at each
T9 leg. Scoring convention fixed: bound ions, if they appear, count
as droplet-retained weight, excluded from the IHe_n histogram and
reported alongside it (matching the twin's bookkeeping).

V0-2 is hereby **closed as corrected + read**; V0-1 (birth law), V0-3
(twin recovery), T8 scope, and T9 re-centering remain the open V0/I.11
decisions. Nothing here discharges F5.

## I.11 immediate goal SET — reproduce the 1D-chord scratchpad closure in 3D MD; V0-1 ADJUDICATED (center-pinned sampler; twin-parity r²+margin arm; Boltzmann stays default); T7 re-scoped unconditional; V0-3 sequenced first (user decision, 2026-07-16)

Discussion outcome (user + assistant, 2026-07-16). The I.11 goal is
sharpened: **reproduce the scratchpad closure — found by the H.2b 1D
chord model (straight chords through the 3D droplet geometry,
closed-form drag/exposure, fate-map evaporation) — in full 3D MD**,
including the twin's ensemble geometry.

**V0-1 adjudicated.** Preliminary quantification of the repo birth
sampler (`sampling/radial_positions.py`, Boltzmann `r²·e^(−U/kT)` at
T = 0.4 K, molecule steepness 14.3324 Å, E_b = 49.4 meV, R = 27.94 Å):
the realized law is **center-pinned** — median birth radius 1.37 Å,
99 % < 3.41 Å (U ≈ 4.2·k_BT at the droplet center; the inward
gradient beats the r² volume factor everywhere) — vs the 1D L1 law's
median ≈ 17–20 Å. The plan's proposed mapping ("adopt Boltzmann into
the twin") is **rejected**: it would leave the closure's position
axis, and with it T5/T6, structurally inert — the T3 lesson repeated.
Decision: Boltzmann stays the physical byte-inert default; the MD
gains the twin's uniform-in-volume r² + hard-margin law as an
explicit arm — recorded as a **twin-parity/capability lever, not a
physical claim** (center-pinning of a heliophilic dopant is
physically defensible; whether physical heterogeneity is instead
carried by the droplet prior or E₀ axis is a later arbitration).

**Consequences (plan §I.11 amended in place):** T7 re-scoped from
conditional margin knob to unconditional **birth-position-law slice**
(`birth_position_law ∈ {boltzmann (default), uniform_volume}` +
margin — exactly the scratchpad's L1 law, hard margin, no smoothing);
the T9 oracle chain gains leg **A′** (= A +
`uniform_volume`) so the position axis flips before the T5 dressing
leg; the "twin adapts to the MD" rule is unchanged (each leg
re-scores at the realized MD law — center-pinned delta at A,
r²+margin from A′ on); the committed twin (V0-3) gains a birth-law
input, and V0-3 is **sequenced first**. New order: V0-3 → T7 → A′ →
T5 → T6 → B/C → T8 → D + re-pilot + scoring. All builds stay behind
their own `[PROCEED TO IMPLEMENTATION]`; nothing discharges F5.

## V0-3 DELIVERED — the 1D twin recovered and committed as `scripts/tier2_h2b_forward_model.py`; all recorded wiring oracles reproduce exactly; birth-law input + `birthlaw` stage added; V0-1 quantification re-issued reproducibly (2026-07-16)

Executed under its own `[PROCEED TO IMPLEMENTATION]`.

**Recovery.** The two surviving scratchpad copies of
`h2b_feasibility.py` are **byte-identical** (cmp) — no version
ambiguity. Committed as `scripts/tier2_h2b_forward_model.py`, code
verbatim except the V0-3 scope additions: (1) repo-relative path
bootstrap (was a hardcoded absolute path); (2) outputs +
per-birth-law fragment cache under the gitignored
`data/runs/h2b_forward_model/` (was the scratchpad dir); (3) the
**birth-law input** `BIRTH_LAW ∈ {"uniform_volume" (native L1 r² +
margin), "boltzmann" (the realized repo sampler via
`sampling/radial_positions.py`)}` threaded through `draw_master`,
with a `margin_weight` guard (margins are uniform_volume-only —
loud `ValueError` under `boltzmann`); (4) the new `birthlaw` stage —
the V0-1 quantification re-issued from the committed script.

**Verification.** All recorded wiring oracles reproduce exactly:
Σ(21) = 0.18783720 eV; production center-pin K = **0.74460**
(exposure·τ 4.8771 ps, t_exit 4.60 ps, v_peak 10.53, v_inf 4.13 —
the Wave-11 landmark row); 9 Å center-pin K = **0.89767**; ballistic
vs analytic 4.3·10⁻⁴; O4 fate map (E₀ = 0.28 → n_det = 6;
0.52 → suppressed/bare); dt-halving ΔK = 3.4·10⁻⁴. The `birthlaw`
stage reproduces the V0-1 adjudication numbers from the committed
code: Boltzmann analytic median 1.36 Å / q95 2.75 / q99 3.40 / mode
1.16 (adjudication's 1.37 = same integral, slightly different grid
convention); the empirical draw through the actual repo sampler
agrees (median 1.35); uniform_volume comparators at margins
{0, 3, 4.67, 6} written to `h2b_birth_law_quantiles.csv`.

**Tests.** `tests/test_tier2_h2b_forward_model.py` — 10 focused
tests (oracle pins with recording-precision tolerances, fate-map O4
hand oracle, ladder variants incl. rq4graded multipliers, r² law
quantiles, margin weights, center-pinned Boltzmann draw, margin
guard, unknown-law guard, `birthlaw` smoke into tmp_path); 10/10
pass in ~3 s; full suite run recorded in this entry's delivery
state.

The twin is now the reproducible oracle for every T5–T9 leg. Next
per the I.11 sequence: Slice T7 (birth-position-law arm), behind its
own trigger. Nothing here discharges F5.

## Slice T7 DELIVERED — `birth_position_law` arm (`boltzmann` default / `uniform_volume` + hard margin, the twin's L1 law); default proven byte-identical to the pre-T7 sampler; CALIBRATION_MAP row 25 (2026-07-16)

Executed under its own `[PROCEED TO IMPLEMENTATION]`; TDD (RED watched
first: ImportError on the missing guard, then 16/16 GREEN).

**As built.**
- `config.py`: `BirthPositionLaw = Literal["boltzmann", "uniform_volume"]`;
  fields `birth_position_law = "boltzmann"` and
  `initial_position_margin_angstrom = 0.0` (declared next to
  `single_initial_position`); `check_birth_position_config` guard wired
  into `SimConfig.validate` — typo guard via the shared
  `_reject_unknown_enum`, finite/≥ 0 margin, and the no-silent-carry
  refusal of a non-zero margin under `boltzmann`.
- `sampling/radial_positions.py`: `_sample_uniform_volume` — the 1D
  twin's L1 law p(r) ∝ r² on [0, R − margin], hard margin, no smoothing
  (user decision); exact inverse-CDF `r = (R − m)·U^(1/3)` per molecule
  (one uniform draw each; no rejection loop); loud `ValueError` when the
  margin consumes any droplet (R − m ≤ 0). The dispatch sits *before*
  the Boltzmann code: the default path is the pre-T7 code, unmodified.
- No preset, generator, or production config selects the arm — the T9
  leg-A′ re-pilot is where it goes live.

**Verification.**
- **Default byte-identity, exact:** default-path draws equal the pre-T7
  sampler's draws elementwise (pre-T7 module reconstructed from
  `git show HEAD`, same cfg/seed/radii, `np.array_equal` True) — the
  RNG draw-order rule holds by construction (the branch reads one cfg
  field and consumes no RNG before dispatch).
- Pre-T7 `cfg.json` loads with the inert defaults (Slice-DS-precedent
  back-compat test); round-trip preserves the new fields.
- `uniform_volume` law: quantiles of (r/cap)³ uniform to < 0.02 at
  n = 20000; margin 0 reaches the surface; hard edge at R − m; the
  committed twin's comparator median 18.47 Å at m = 4.67 reproduced;
  per-droplet support respected under mixed radii; liveness vs the
  center-pinned Boltzmann law (median 1.4 Å vs > 20 Å).
- Neutral `build_initial_state` integration: molecule-centre radii
  bounded by R − m, not center-pinned.
- Tests: `tests/test_birth_position_law.py`, 16 focused tests; narrow
  neighboring suites (radial positions, initial states, config guards,
  twin) 214 passed; **full suite 2354 passed** (2338 + 16).

**CALIBRATION_MAP:** row 25 added — the law is an **arm, not a knob**
(twin-parity/capability lever; the physical default stays `boltzmann`,
center-pinned per V0-1); the margin is **Bounded**, firm band
{3, 4.67, 6} Å (H.2b D5); tally updated (Bounded 5–6 → 6–7).

Next per the I.11 sequence: **T9 leg A′** (delivered T3 configuration +
`uniform_volume` at the closure cell's margin, twin re-scored at that
exact configuration) — behind its own trigger. Nothing here
discharges F5.

## T9 leg A′ TRIGGERED — build delivered + twin predictions PRE-REGISTERED before any MD (2026-07-16)

User trigger: leg A′ with **all four configs** ("I agree this makes us
understand our model better"). Margin fixed by the Step-1c record: the
full-house cells sit at **3 Å** (best cell p = −1, v_c = 7.5, τ = 3.8,
d080, m3, E₀ = 0.25 — findings §4k), not the mid-band 4.67 Å guessed in
discussion.

**Build (delivered under this trigger; TDD on the new surfaces):**
1. Twin (`scripts/tier2_h2b_forward_model.py`): `integrate_pairs` gains
   the Slice-T1 `capped_cubic` tail (`v_c`/`p_tail`;
   `drag_gamma_tail_amu_per_ps`; `v_c=None` keeps the pure-cubic
   arithmetic verbatim — byte-inert, oracles unchanged); new
   `legaprime` stage — the twin re-scored at the leg-A (center-pinned
   delta) and leg-A′ (uniform_volume m = 3 Å) configurations exactly:
   fixed N = 2000, **undressed n_eject = 21** (ANCHOR_N_START; T5 not
   flipped), E₀-coupling p = 0 (T6 not flipped), ladder tables imported
   from the MD generator (the MD's exact tables, incl. floor1's
   13.3 meV — not the twin's 13.25 meV convention), τ as the exact
   K-rescale. Tests +5 (tail hand oracle, in-band parity, tail lifts
   v_inf/lowers K, guard, stage smoke) — twin suite 15/15.
2. `scripts/tier2_common.py :: build_biphasic_cfg` gains two
   None-sentinel byte-inert kwargs (`birth_position_law`,
   `initial_position_margin_angstrom`).
3. Generator (`gen_tier2_md_confirmation.py`): `LEG ∈ {"a", "aprime"}`
   USER SETTING (active: `aprime`; `a` reproduces the four T3 dirs
   byte-identically — locked by tests incl. an exactly-one-lever
   field-diff assertion); A′ dirs are `…_tier2probe_conf270_apc{1..4}`
   (tag charset forbids underscores). Generator + birth-law suites
   32/32.

**Pre-registered twin predictions (m = 20000 chords; CSVs:
`data/runs/h2b_forward_model/h2b_leg_aprime_predictions.csv` + `…_ke.csv`):**

| config | leg-A twin n̄ (MD T3 measured) | leg-A′: trapped | suppressed→bare | n̄_det | top solvated bins |
|---|---|---|---|---|---|
| c1 | 10.0 (8.63) | 0.033 | 0.089 | 6.22 | n1 0.172, n2 0.107, n3 0.082 |
| c2 | 10.0 (8.77) | 0.035 | 0.035 | 6.55 | n1 0.170, n2 0.113, n3 0.091 |
| c3 | 8.0 (6.97) | 0.067 | **0.336** | 4.48 | n1 0.069, n2 0.063 (flat) |
| c4 | 10.0 (8.39) | 0.033 | 0.118 | 6.05 | n1 0.138, n2 0.086, n3 0.079 |

Registered reads (reported-not-adjudicated, per program convention):
- **AP-P1 (S2c-P1 re-anchored):** the MD A′ histograms broaden from the
  T3 ±2-bin cluster to the twin's near-full-range shape, with n₁ the
  top solvated bin on the rq4graded configs and a nonzero bare class.
- **AP-P2 (offset carry):** the leg-A twin−MD offset (twin reads
  +1.2–1.4 He high — real cascade/m(t)/pickup vs fate map) carries to
  A′ at similar magnitude; the *shape* (two-sided race weight) is the
  claim, not the absolute n̄.
- **AP-P3 (trapped class):** first nonzero droplet-retained weight in
  MD (twin: 3–7 %, largest on the current-law c3) — counted per the
  V0-2 convention (excluded from IHe_n, reported alongside).
- **AP-P4 (control split):** c3 (current law, flat ladder) shows the
  largest suppressed/bare fraction (τ = 6.55 → smallest K after
  rescale… the race margin, not the tail, drives it) — if the MD
  instead shows c3 *least* bare-heavy, the twin's K-rescale picture is
  wrong at production geometry.

MD pilots (N = 50, four stages, `apc1..apc4`) launch next; scored
against these numbers. Nothing here discharges F5.

> **NB (2026-07-16, mid-leg): the trapped class arrived in MD on the
> first execution attempt — and collided with the detection stage's
> decoupling contract.** The apc1 pilot failed at the P1–P3 handover
> guard: 5/100 ions still inside helium at the 1000 ps relaxation end
> (worst ρ = 0.67). The V0-2 escape criterion identifies **exactly
> those 5 as energetically bound** (margins −2.9 to −35 meV, radii
> 24–44 Å; every escaped ion sits ≥ 881 Å) — the first droplet-retained
> ions ever expressed in MD, at 5 % vs the twin's pre-registered 3.3 %
> (AP-P3, within N = 50 statistics). A bound ion never decouples, so
> the guard's extend-relaxation remedy cannot apply. **Build (this
> trigger):** `detection_droplet_retained_policy ∈ {refuse (default,
> byte-inert = the delivered loud guard), exclude}` — under `exclude`,
> energetically bound violators (KE + U < E_bind at handover, via the
> shared `droplet_potential`) are classified with the new additive
> state reason `droplet_retained`, carry their handover state verbatim
> (zero events), and are excluded from the IHe_n read; *unbound*
> violators still refuse loudly under either policy (an undecoupled
> escaper is a configuration error, not physics). Threaded as a
> None-sentinel `build_biphasic_cfg` kwarg; the A′ leg stamps
> `exclude` (a bookkeeping convention forced by the leg's own physics,
> not a second lever — the one-lever test asserts exactly
> {birth law, margin, retained policy}). Tests: +4 detection-stage
> (refuse-default parity, bound-exclude with handover-state carry,
> unbound-still-refuses, config guard); detection + generator +
> birth-law suites 79 passed. apc1 regenerated (stale-artifact
> recovery; its partial dir predated the policy).
>
> **Second wrinkle (same day): a quasi-bound capture-in-progress.**
> Under the exclude policy apc1/apc2 completed (n̄_detect 8.24 / 7.25),
> but apc3 refused on one *unbound* violator — ion 9: escape margin
> **+1.14 meV** (barely unbound), at 41.7 Å moving **inward** at
> −0.31 Å/ps (|v| 1.07 — mostly tangential). It is mid-capture: it
> re-enters the droplet, drag eats the ~1 meV margin, and it joins the
> bound class — the 1000 ps cap froze the film before the class was
> decided. The exclude policy is deliberately not widened (excluding
> unbound ions would silently corrupt the histogram); instead the A′
> leg's relaxation cap is extended 1000 → **3000 ps**
> (`APRIME_RELAXATION_TIME_PS` — numerical adequacy for the position
> axis's near-barrier transients, not a physics knob; detection
> re-reads at 8.53 µs regardless; leg A keeps 1000 ps byte-identically)
> and all four apc dirs regenerate at the shared cap. The one-lever
> test now asserts exactly {birth law, margin, retained policy,
> relaxation cap}.

## T9 leg A′ EXECUTED — AP-P1..P4 all CONFIRMED (W₁ 0.55–0.69 bins; suppressed ordering transfers exactly); KE composition diverges through the real cascade (I57); barrier-corrected retained criterion delivered; E2-drag OQ opened (2026-07-16)

Full record: findings **§4m + I55–I57**. Four dirs on disk
(`…_tier2probe_conf270_apc{1..4}`, all five artifacts, cap 8000 ps,
uniform).

**Execution history (honest record).** The relaxation cap chased the
marginal class 1000 → 3000 → 5000 → 8000 ps before the diagnosis
landed: apc3's "ion 9" is a **centrifugal resonance** (+1.14 meV above
the radial escape threshold, below its effective-potential barrier,
~5 Å outward drift per 3000 ps under the **zero-gamma** E2 "coulomb"
mode — conservative dynamics, no capture possible). Consequence
delivered under the leg trigger: the `exclude` policy's bound
criterion upgraded from radial (KE + U < E_b) to the **exact
conservative condition** E_tot < max over the outward path of
[U(r′) + L²/(2 m r′²)] (`_conservatively_bound`; +1 resonance-orbit
test, detection suite 48 passed). apc1–apc3 detections re-issued
standalone from their `relaxation.npz` under the corrected criterion
(the Wave-7 re-read precedent); apc4 via the generator. Retained
counts: 5/6/11/5.

**Verdicts (vs the pre-registered twin numbers, same-day entry
above):** AP-P1 confirmed (every histogram broadens to the two-sided
race shape; first MD weight ever at n = 0 and n = 1); AP-P2 confirmed
sharpened (twin−MD −0.3..−0.6 He; W₁ 0.55–0.69 bins); AP-P3 confirmed
(trapped class real, ~1.6× the twin's 150 ps read, c3 largest both
sides); AP-P4 confirmed (suppressed/bare ordering transfers exactly:
MD 0.315/0.116/0.095/0.032 vs twin 0.336/0.118/0.089/0.035). **The
divergence (I57):** per-bin mean KE at small n is ×2.5 the twin
(n₁ 2.42–2.48 eV vs 0.97) — the fate-map↔real-cascade channel
decouples bin weight from bin occupants; the (n, mean-KE) read is
composition-sensitive and waits for the dressed geometry.

**Open item raised to the user:** whether E2 should carry drag (a
drag-live relaxation arm) — physically an orbiting in-helium ion is
captured, and the conservative convention instead produces long-lived
resonances handled today by classification + an 8000 ps leg cap;
N = 500 (T9 endgame) makes marginal ions a certainty, so this wants
adjudication before then. Next per the I.11 sequence: **Slice T5**
(initial-shell dressing arm), behind its own trigger — leg B's twin
re-score must first re-register predictions at the dressed
configuration. Nothing here discharges F5.

> **NB (2026-07-17): plan §I.11.2 added** — the post-leg-A′ decision
> register (next-session entry point): (1) the n₁-composition re-read
> (zero MD; the v = 19.1 Å/ps > 15.4 symmetric-ballistic argument
> identifying the mass-asymmetric Coulomb split as I57's mechanism,
> with near-edge births as the correlate); (2) the E2 dissipation
> adjudication (zero-gamma / drag-live / Landau-gated arms — the slow
> orbits straddle the Landau cutoff, so naive drag-on may overdamp);
> (3) the leg-B pre-registration amendment (per-bin mean KE mandatory
> alongside histograms). Suggested order (1) → (2) → T5 trigger.

## I.11.2 item 1 (n1-composition re-read) EXECUTED — pre-registered Coulomb-share hypothesis REFUTED as dominant; the driver is the cold-shed momentum convention (OQ-J/RQ10 fired); item 2 ADJUDICATED (arm (c) Landau-gated); item 3 ADOPTED + sharpened (2026-07-17)

All three §I.11.2 items settled in one session (user decisions: "1
PROCEED / 2 go with (c) / 3 agree and amend"). Full record: findings
**§4n + I58–I60 + OQ-J**; `RESEARCH_QUESTIONS.md` **RQ10** + coupling
row; plan §I.11.2 status block.

**Item 1 — executed (zero new MD; scratch read over the four on-disk
apc dirs, probe convention, zero repo-code change).** Method: exact
per-fragment Coulomb-work integral from the stored trajectories
(pair-sum renormalized against U_c(0) − U_c(end)), per-ion birth
radius / chord cosine / first-shed / exit time / mass-at-exit, and the
empirically verified 5-term ledger
(E_kin + E_pot + E_dissip + E_mass_transfer + E_int conserved to
−0.0008 eV mean end-drift) → per-class decomposition
KE_det = 2.70 + coul_x − drag + boost + resid, closing to ≤ 0.06 eV on
every live class.

Verdict — **the §I.11.2 hypothesis ("v = 19.1 > 15.4 Å/ps requires the
mass-asymmetric Coulomb split") is refuted as dominant**:

- measured Coulomb-share excess on n1: **+0.08 eV (≈ 3 %)** — n1 ions
  exit full-shell at 0.7–0.9 ps, mass-symmetric through the
  acceleration (m̄ 210/209 amu self/partner over the first 2 ps), and
  shed only after ejection; the hypothesized light-early compounding
  route does not operate;
- the dominant term is the **cold-shed momentum convention** of the
  delivered evaporation channel
  (`physics/evaporation.py` composes
  `mass_jump.cold_shed_velocity_components`): He left at lab rest,
  v → v·m/m′ per shed (verified event-by-event: v_ratio ≡ m_ratio),
  injecting **+0.91 eV** over the post-exit 21→1 strip
  (×1.611 = m21/m1 in KE); E_mass_transfer books exactly the injected
  energy;
- the **co-moving counterfactual**
  KE_cf = ½·m_det·(v_det·m_det/m_exit)² reproduces the twin per config
  to ≤ 2 % (n1: 0.955/0.952/0.256/0.933 eV vs twin 0.97/—/0.261/—) —
  the twin is the co-moving convention in disguise; **I57 is entirely
  the shed-frame convention**, with the composition axis (near-edge
  outward births, drag deficit 1.19 vs 2.66 eV) deciding occupancy and
  the convention deciding speed.

Consequences recorded: **OQ-J opened** (findings §7) = **RQ10**
(`RESEARCH_QUESTIONS.md`, coupled to RQ2 as one "what leaves with the
He" resolution and to RQ3 via the bare-bin KE fragmentation
discriminator, 1.18 vs 3.3 eV); the experimental n1 mean 1.302 eV is
convention-spanned (cold 2.42–2.48 / co-moving 0.93–0.96 on the
capped-tail configs); **calibration-transfer warning** — the §4k
(v_c, τ) KE closure is co-moving-based and does not transfer to a
cold-shed MD at small n. A8's [Nat23] grounding is an at-rest result
(cold ≡ co-moving there); Tier-1a had adjudicated continuous-velocity
as the physical path. Candidate resolution: an interchangeable
shed-convention enum (both operators already in `mass_jump.py`) behind
its own trigger — **user adjudication required before the T9 endgame**.

**Item 2 — ADJUDICATED (user, 2026-07-17): arm (c) Landau-gated drag**
for the E2 relaxation dissipation question (γ = 0 below the config
Landau cutoff `v_limit`, drag above — sub-Landau superfluid motion is
dissipationless; captures marginal/orbiting ions physically without
overdamping slow motion). Build stays behind its own
`[PROCEED TO IMPLEMENTATION]`, required before any N = 500 run; until
built, the `_conservatively_bound` classification + leg-level caps
remain the handling.

**Item 3 — ADOPTED + sharpened by item 1:** the leg-B (T5) twin
re-score pre-registers per-bin mean detected KE alongside the
histograms/classes, and every KE claim must **state its
shed-convention basis** (twin = co-moving by construction; delivered
MD = cold-shed) — the two KE axes are not comparable without it (I59).

Next per the I.11 sequence: the **OQ-J shed-convention adjudication**
(new, from item 1), then **Slice T5** (initial-shell dressing arm)
behind its own trigger with the amended pre-registration. Nothing here
discharges F5.

## OQ-J working convention ADJUDICATED — co-moving shed for the twin-parity legs; cold retained as the interchangeable bound arm; RQ10 stays open (user decision, 2026-07-17)

User decision, same session as the item-1 re-read ("as we want to
reproduce twin results ... switch to co-moving convention and discuss
why for now; this remains an open research question"). Assistant
concurred on physics grounds independent of the twin-parity motive;
the recorded rationale:

1. **Frame physics.** Evaporation is thermal in the complex rest frame
   (u_thermal ~ 0.5–1 Å/ps ≪ v_lab 10–19 Å/ps) — co-moving is the
   correct zeroth order. Cold shed requires a directed ~75 meV/He
   backward launch (0.5·m_He·v², ~30× thermal, several D₀ rungs) with
   no energy source — the ledger books an unsourced +0.91 eV injection
   on a full strip. Defensible only as a bound, not as physics.
2. **No contradiction with A8/[Nat23]:** the Na⁺ complex is at rest,
   where cold ≡ co-moving; co-moving extends the source correctly to
   moving complexes.
3. **Intra-model consistency:** Tier-1a adjudicated
   continuous-velocity as the physical path; the Tier-2 channel's cold
   composition was a silent divergence, now removed.
4. **Calibration transfer:** the §4k (v_c, τ) closure and the leg-B KE
   bar are co-moving-based; under co-moving MD the I59 gap closes by
   construction.

**Decision shape (build behind its own `[PROCEED TO IMPLEMENTATION]`):**
interchangeable `evaporation_shed_convention ∈ {cold, co_moving}` —
**`cold` stays the byte-inert config default** (every pre-existing
`cfg.json` load remains an honest record of what ran), **the T5+/leg-B
generators stamp `co_moving` explicitly** (the T7 birth-law
precedent); production default revisited at config assembly. Both
operators exist in `mass_jump.py`; the co-moving one needs the
per-ion-mass generalization the cold one got (fires diverge), threaded
as the single mass source through the ion/relaxation/detection stages
(the resolve_ladder single-surface precedent). Ledger: E_mass_transfer
flips to the continuous form (+0.5·m_He·|v|² carried away per shed);
the 5-term closure holds in the same form.

**Consequences recorded:** histograms/classes are essentially
convention-blind (mass events never enter the RRK rate, gate, or
cooling; n₁/bare occupants shed post-exit; deep bins only a small
in-bubble drag feedback via m(t)) — the leg-A′ verdicts stand; only
the KE axis moves (n₁ 2.42–2.48 → ≈ 0.93–0.96 eV, onto the twin). All
existing probe/T3/apc dirs stay valid as records of the cold arm
(their KE reads were never arbitration inputs; the Wave-11 K
measurement was a shed-free suppressed ride — convention-free). New
legs regenerate under the stamped arm per the stale-artifact policy.

**RQ10 stays open** as the physical-resolution question: the true
convention is co-moving + isotropic thermal recoil — the momentum side
of RQ2's ε (how much recoil, correlated how); the two delivered
operators bracket it (ε → 0 vs maximal backward kick). One coupled
RQ2+RQ10 literature/domain-expert discussion.

Next per the amended I.11 sequence: **the shed-convention enum build +
Slice T5** (initial-shell dressing arm) behind their trigger(s), with
the leg-B pre-registration carrying per-bin mean KE on the stated
co-moving basis. Nothing here discharges F5.

## Shed-convention enum DELIVERED — `evaporation_shed_convention ∈ {cold (byte-inert default), co_moving}` threaded through all three stages from one config field; CALIBRATION_MAP row 26 (2026-07-17)

Executed under its own `[PROCEED TO IMPLEMENTATION]` (the OQ-J
adjudication entry above is the design record). TDD: the 18-test suite
watched RED (AttributeError on the missing field) before any
implementation.

**As built.**
- `config.py`: `EvaporationShedConvention = Literal["cold", "co_moving"]`;
  field `evaporation_shed_convention = "cold"` (byte-inert default,
  declared beside the Slice-Q evaporation surface); typo guard via
  `_reject_unknown_enum` added to `check_evaporation_config` (check 3).
- `physics/mass_jump.py`: `continuous_velocity_shed_components`
  generalized to **per-ion mass** (the exact mirror of the cold
  generalization — fires diverge; `_check_masses_shed` + array-aware
  `_check_removed_count`; scalar path behaviour-identical). The
  co-moving ledger term `+0.5·m_He·|v|²` is exact (no reduced-mass
  correction — the He leaves at exactly `v`).
- `physics/evaporation.py`: `evaporation_step` /
  `evaporation_step_components` gain `shed_convention="cold"` and
  dispatch cold vs co-moving operators; `_check_shed_convention`
  point-of-use twin guard for direct physics-module callers.
- `simulation/ion_propagation_step.py`: `biphasic_step` passes
  `cfg.evaporation_shed_convention` (the relaxation stage rides the same
  call through its cfg view — no E2-side change needed).
- `simulation/detection_stage.py`: the Gillespie event loop dispatches
  the same config field (`cold_shed` vs `continuous_velocity_shed`).
- `scripts/tier2_common.py :: build_biphasic_cfg`: None-sentinel
  `evaporation_shed_convention` kwarg (T7 pattern). No preset, generator,
  or production config selects the arm — the T5/leg-B legs stamp it.

**Verification.**
- `tests/test_evaporation_shed_convention.py` — 18 focused tests: config
  default/guard/back-compat (pre-enum `cfg.json` loads with `cold`)/
  round-trip; per-ion co-moving components hand oracle + single-atom
  parity + underweight fail-loud; channel-level RNG contract (same seed →
  identical fire pattern and (E_int, n) chain across arms; exact ledger
  relation `dE_cold = −(m/(m−m_He))·dE_co_moving`); `biphasic_step`
  threading (co-moving keeps v elementwise, books `+0.5·m_He·|v|²`);
  detection stage (co-moving keeps v through all events, positive event
  ledgers, **machine-precision 5-term closure across the stage**; cold
  kicks up, negative ledgers); `build_biphasic_cfg` pass-through.
- Byte-identity of the cold default: the pre-slice suites lock the
  delivered behaviour and pass unchanged — neighboring suites
  (evaporation, mass_jump, biphasic_step, detection, relaxation,
  tier2_common, gen_tier2_md_confirmation, birth law) **373 passed**;
  **full suite 2385 passed**.

**CALIBRATION_MAP:** row 26 — the convention is an **arm, not a knob**;
working convention `co_moving` for the twin-parity legs (per the
adjudication), `cold` retained as the diagnostic bound arm; RQ10 open.

Next per the amended I.11 sequence: **Slice T5** (initial-shell dressing
arm) behind its own trigger, with the leg-B twin re-score pre-registering
per-bin mean detected KE on the stated co-moving basis. Nothing here
discharges F5.

## T9 leg A″ TRIGGERED — A′ + co-moving shed convention; build delivered + predictions PRE-REGISTERED before any MD (2026-07-17)

User trigger ("[PROCEED TO IMPLEMENTATION] leg A″"), following the
discussion: leg B (T5) will run under `co_moving`, so jumping from the
cold A′ baseline to a dressed co-moving leg B would flip two levers at
once — the convention flip needs its own leg first. Leg A″ validates
the delivered enum end-to-end in the real three-stage pipeline at
production kinematics and re-baselines the KE axis on the twin's
co-moving basis.

**Build (delivered under this trigger; TDD on the new lock test):**
`gen_tier2_md_confirmation.py` gains `LEG = "aprime_cm"` (active) —
identical to `aprime` except `evaporation_shed_convention="co_moving"`
stamped through the delivered `build_biphasic_cfg` kwarg; dirs carry
the `apcm` prefix (`…_tier2probe_conf270_apcmc{1..4}`). Lock tests:
the A′ one-lever test re-pinned via monkeypatch (A′ additionally
asserts `cold`); new A″ test asserts the field-diff vs A′ is exactly
`{evaporation_shed_convention}` and the dirs are distinct. Generator +
common + convention suites 94 passed.

**Pre-registered predictions (from the §4n decomposition — no new
analysis; same four configs, same N = 50 bridge seed, margin 3 Å,
exclude policy, 8000 ps cap):**

- **A″-P1 (histogram invariance).** Histograms/classes reproduce A′
  within N = 50 noise: the (E_int, n) jump chain is convention-blind;
  n₁/bare occupants shed entirely post-ejection (bin-identical), deep
  bins may drift slightly via the in-bubble m(t)/v drag feedback. The
  suppressed/bare ordering c3 ≫ c4 > c1 > c2 transfers unchanged. A
  material histogram move falsifies the I59 "convention-blind
  histogram" claim — itself a finding.
- **A″-P2 (the KE landing — the headline).** Per-bin mean detected KE
  at n = 1 lands on the §4n co-moving counterfactuals ≈ the twin:
  **0.955 / 0.952 / 0.256 / 0.933 eV** (c1/c2/c3/c4) vs twin
  0.97 / — / 0.261 / —. I57 closes in the real pipeline; the residual
  twin−MD KE gap is composition-level only.
- **A″-P3 (ledger signature).** `E_mass_transfer` ≥ 0 for every shed
  ion (positive carried-KE booking); 5-term `ion_ledger_closure`
  residual unchanged; suppressed/bare and retained classes carry
  *identical* KE to A′ (no sheds → convention-free) — the free wiring
  oracle.
- **A″-P4 (speed re-scale).** n₁ detected speeds drop from ~19.1 to
  **≈ 11.9 Å/ps** (the momentum base at m_exit; v_det ≈ v_exit under
  co-moving sheds).

MD pilots (N = 50, four stages, `apcmc1..4`) launch next; scored
against A′, the §4n counterfactuals, and the twin. Nothing here
discharges F5.

## T9 leg A″ EXECUTED — A″-P1..P4 all CONFIRMED (P3 sub-claim refined): the co-moving convention lands the twin KE axis in real MD; twin↔MD parity is now two-axis; leg-B baseline re-set to the apcm dirs (2026-07-17)

Full record: findings **§4o + I61**. Four dirs on disk
(`…_tier2probe_conf270_apcmc{1..4}`, all five artifacts; apcmc1 from the
first launch, apcmc2 regenerated after a mid-relaxation kill — the
partial dir removed per the generator's documented recovery path,
apcmc3/4 fresh; resume guard exercised as designed).

**Verdicts vs the pre-registered numbers (same-day entry above):**
- **A″-P1 CONFIRMED** — W₁(A″, A′) = 0.14–0.20 bins over 22; the
  suppressed/bare ordering c3 ≫ c4 > c1 > c2 exact; class fractions
  within ≤ 3 ions. The histogram is convention-blind in real MD (I59
  survives).
- **A″-P2 CONFIRMED** — n₁ mean detected KE **1.005 / 1.025 / 0.312 /
  0.991 eV** (c1..c4) vs twin 0.972/0.963/0.261/0.949 and §4n
  counterfactuals 0.955/0.952/0.256/0.933; the cold ×2.5 divergence
  (I57) collapses to ×1.03–1.20, and the **whole solvated per-bin KE
  curve matches the twin** (c1 n2 0.74/0.74, n3 0.59/0.57, n5
  0.36/0.36). Small uniform +0.04–0.06 eV residual above the
  counterfactual: composition-level (in-bubble sheds), not
  convention-level.
- **A″-P3 CONFIRMED in signature, strict sub-claim refined** — every
  shed ion books E_mass_transfer > 0; 5-term closure drift unchanged
  (max transient 0.11–0.14 eV, end-drift −8·10⁻⁴ eV — the cold arm's
  values). The "suppressed/retained KE identical to A′" oracle holds to
  ≤ 10–20 meV only: **ion pairs are Coulomb-coupled** (an opened
  partner's convention-shifted trajectory perturbs its never-shed twin),
  and **1–3 marginal near-barrier ions per config flip retained**
  (c3: 11 → 14) — meV-scale fragility of the marginal class, sharpening
  the E2-dissipation (Landau-gated) adjudication's stakes before N = 500.
- **A″-P4 CONFIRMED** — n₁ speeds 12.1–12.3 Å/ps (pre-reg ≈ 11.9; cold
  19.1).

**Consequences:** the OQ-J working-convention adjudication is validated
end-to-end at production kinematics; the leg-B (T5) A/B baseline is the
**apcm** dirs — MD and twin now share the co-moving KE basis by
construction, satisfying the item-3 stated-basis amendment. The apc
(cold) dirs remain on disk as the OQ-J bound-arm record. Next per the
I.11 sequence: **Slice T5** (initial-shell dressing arm) behind its own
trigger, with the leg-B twin re-score pre-registering histograms *and*
per-bin mean KE at the dressed configuration. Nothing here discharges F5.

## I11 post-A″ code review EXECUTED — fixes applied; the exclude-policy retained class is now actually excluded from every IHe_n read and its artifacts load; one rule-1 carry recorded (2026-07-18)

High-effort multi-agent review of `c175e2f..069a166` (Slice T3 generator
— delivered alongside the I11 plan, never reviewed — plus the I11 work:
V0-3 twin, Slice T7, T9 leg A′ detection-stage changes, the OQ-J
shed-convention enum, leg A″): 4 finders / 16 independent verifiers,
19 candidates verified, 0 refuted, merged to 10 distinct findings. All
fixed under the user's `[PROCEED TO IMPLEMENTATION]` **except finding 9,
which the user scoped to this log entry only** (the carry below). Full
suite green after the fixes (count in the commit record).

**F1+F2 (CONFIRMED, the severe cluster): the
`detection_droplet_retained_policy="exclude"` class was excluded from
the event loop only.** `STATE_REASONS` was never extended with
`"droplet_retained"`, so (a) `load_detection_result` refused the very
`detection.npz` the stage saves (reproduced round-trip) and (b)
`reason_fractions` silently dropped the retained class (fractions
summed < 1). Worse, retained rows carry the **verbatim in-droplet
handover state** in `n_detected` / `E_kin_detected_eV` /
`mass_detected_kg`, and no numeric consumer masked them — the
generator's `n_detect_mean` and any histogram built from the arrays
folded the trapped subpopulation into the terminal IHe_n read. Fixes:
`STATE_REASONS` gains the fourth member; new
`DetectionResult.detected_mask` property (documented mask contract in
the data-contract docstring); `_terminal_n_and_source`
(`postprocess/size_distribution.py` — the single histogram source)
excludes retained rows on the detected branch (all-retained → loud
empty-ensemble error); the generator print masks and now reports
`droplet_retained=k/2N`; the staircase report gains the
`frac_det_droplet_retained` column. Tests: real save→load round trip
(the old test checked the in-memory array only and missed the load
refusal), fractions-sum-to-1, mask contract, histogram exclusion.
**Residual risk flagged:** the §4m/§4o leg reads were made with scratch
conventions; the delivered apc/apcm `detection.npz` artifacts (which
were unloadable through the delivered loader until this fix) should be
re-read once through the now-masked pipeline to confirm the recorded
histogram/class numbers are unchanged.

**F3 (CONFIRMED): `uniform_volume` under `single_initial_position=True`
was silently inert.** `build_initial_state` zeroes `r0` after sampling,
so the advertised T7 lever would not move positions while still
*shifting the RNG draw stream* (different draw count than the Boltzmann
rejection sampler) — a run reproducing neither leg A nor leg A′.
`check_birth_position_config` now refuses the combination loudly (rule
4 in its docstring); the test fixture selects the position-live arm.

**F4 (CONFIRMED): `_conservatively_bound` omitted the residual pair
Coulomb on an unverified "partner far away" assumption.** The omitted
term (14.4 eV·Å / r_sep ≈ 14 meV at 1000 Å, ~0.3 eV in-droplet)
exceeds the guard's actual meV operating margins (apc3 "ion 9":
+1.1 meV), so a near-threshold escaper with a close partner could be
silently booked `droplet_retained`. The criterion now credits the
**full** pair energy (via the delivered `ion_interaction_potential`,
`[:N]/[N:]` pairing, `E_coulomb_scale` respected) to each fragment —
an over-estimate that keeps `retained` verdicts certain and hands
Coulomb-marginal ions back to the loud violator guard. Directly
relevant to the A″-P3 refinement (1–3 marginal ions flipping retained
per config): the marginal class is now classified with the pair term
in. New test: a close partner (120 Å) flips a well-trapped ion to the
loud refusal; existing retained tests re-pinned with partners at 10⁵ Å.

**F5 (CONFIRMED, leaked debug state):** `tier0_drag_comparison.py`
`ENERGY_FIGURE`/`FORCE_FIGURE` were committed flipped `False→True`
inside `1e64b37` (unrelated to its scope; the post-processing layer is
bug-fix-only). Reverted to the delivered defaults.

**F6 (CONFIRMED):** the twin's `legaprime` stage was missing from both
the module-docstring Usage line and the unknown-mode `SystemExit` list —
the leg-A′/A″ re-score stage was undiscoverable from the script's own
help. Both lists now include it.

**F7 (CONFIRMED, latent):** `margin_weight` branched on the module
global `BIRTH_LAW` while `draw_master` takes a per-call `birth_law` —
the two halves of the importance scheme could silently disagree.
`margin_weight` now takes the same optional `birth_law` parameter with
the same resolution rule (call sites unchanged; numerics identical).

**F8 (CONFIRMED, test hygiene):**
`test_leg_aprime_cm_flips_exactly_one_lever_vs_aprime` asserted the
mutable USER SETTING `script.LEG == "aprime_cm"`; it now monkeypatches
like every sibling, so routine leg rotation no longer turns the suite
red.

**F10 (CONFIRMED, stale data contract):** `DetectionResult` docstrings
updated — `event_dE_mass_transfer_eV` sign is convention-dependent
(cold < 0 / co_moving > 0), `state_reason` documents all four legal
values, and the module-header event-loop paragraph names both shed
operators. Plus the trivial `np.repeat(a, 1)` no-op removed from the
twin's `stage_levers` (bit-identical).

**F9 — rule-1 carry (user decision 2026-07-18: log entry only, no code
change).** The committed twin `scripts/tier2_h2b_forward_model.py`
re-implements physics the package exports: `complex_mass_amu`
(`physics/shell_schedule`), the Coulomb constant `14.39964548`, the
eV↔amu·Å²/ps² conversion, and an integer-support Wasserstein. This is
the deliberate V0-3 verbatim-recovery trade-off (the scratchpad is
committed byte-faithfully so the recorded wiring oracles pin the exact
code that produced the Step-1b/1c closure). **Carry condition:** if any
package constant or the mass table is ever corrected, the twin's copies
drift silently and every T5–T9 oracle comparison scores against stale
physics — the carry is retired by re-syncing (and re-pinning the
oracles) at the first twin edit that touches those numbers, or at the
T9 endgame, whichever comes first.

**Not changed:** checkpoint schema (v7 untouched — `detection.npz` is
not an IonCheckpoint and its own version stays 1: the reason vocabulary
widened but no field changed), RNG draw order (the F3 guard *prevents*
an accidental stream shift), physical constants, presets.

## F1/F2 residual-risk re-read EXECUTED — all eight apc/apcm `detection.npz` load through the fixed pipeline and reproduce the recorded §4m/§4o numbers exactly; the leg-B baseline is certified (2026-07-18)

The review entry above flagged that the delivered leg-A′/A″
`detection.npz` artifacts (unloadable through the delivered loader until
the F1 fix; the §4m/§4o reads were made with scratch conventions) should
be re-read once through the now-masked pipeline. Executed as a zero-MD,
read-only scratch check (probe convention, zero repo-code change):
each of the eight dirs (`…_conf270_apc{1..4}` cold /
`…_conf270_apcmc{1..4}` co_moving) loaded via the delivered
`load_detection_result`, aggregated under `detected_mask` +
`reason_fractions`, and histogrammed through the fixed
`compute_terminal_shell_distribution` (§4m conventions: retained
excluded, suppressed at bin 0, fractions over the detected ensemble).

**Result: ALL PASS — every recorded number reproduces.** Per config:
retained counts 5/6/11/5 (apc) and 6/7/14/6 (apcm; c3 = 14 as the
A″-P3 refinement recorded — the apcm c1/c2/c4 counts were previously
implicit and are now on record); suppressed fractions
0.095/0.032/0.315/0.116 (A′) and 0.096/0.032/0.326/0.117 (A″); n₁
fractions, n̄_det, n₁ mean KE (2.479/2.473/0.666/2.423 cold;
1.005/1.025/0.312/0.991 co-moving) and n₁ speeds (19.10 cold c1;
12.16/12.28/6.77/12.08 co-moving) all match the §4m/§4n/§4o tables to
print precision. `reason_fractions` sums to exactly 1.0 with the
`droplet_retained` member populated on every dir; the histogram source
tags `detected` and drops exactly the retained rows.

**Consequence:** the F1/F2 residual risk is discharged — the scratch
reads behind §4m/§4o and the delivered masked pipeline agree, so the
apcm dirs stand as the certified leg-B (T5) A/B baseline with no
regeneration needed. Script: session scratchpad
(`reread_apc_apcm_detection.py`), not committed (read-only check, no
repo surface). Next per the I.11 sequence: **Slice T5** (initial-shell
dressing arm) behind its own trigger, with the leg-B twin re-score
pre-registering histograms *and* per-bin mean detected KE on the stated
co-moving basis. Nothing here discharges F5.

## Slice T5 DELIVERED — `initial_shell_model ∈ {full (byte-inert default), density_tied}`: the H.3b depth-dependent dressing law in the biphasic seed, per-ion n₀/mass/E_pot fold from the shared erf-complement surface (2026-07-18)

Executed under its own `[PROCEED TO IMPLEMENTATION]` (plan §I.11 Slice
T5, reviving §H.3b verbatim; entry preconditions met — the apcm leg-B
baseline certified by the re-read entry above). TDD: the 16-test suite
watched RED (ImportError on the missing guard) before any
implementation.

**As built.**
- `config.py`: `InitialShellModel = Literal["full", "density_tied"]`;
  field `initial_shell_model = "full"` (byte-inert default, declared in
  its own T5 block beside the Slice-Q evaporation surface);
  `check_initial_shell_config` wired into `SimConfig.validate` — typo
  guard via `_reject_unknown_enum`, and `density_tied` **refused outside
  `mass_scenario='biphasic'`** (the dressing is read only by the
  biphasic column-0 seed; silent inertness would repeat the F3 class —
  refused loudly instead). Center-pinned births + `density_tied` stay
  deliberately legal: physically inert (ρ̂(center) ≈ 1 → n₀ = 21
  exactly, the H.3b oracle) with no RNG-stream shift.
- `simulation/ion_initial_state.py`: under `biphasic` + `density_tied`
  the seed computes `d_birth,i = |r_i(0)| − R_i` and
  `n₀ᵢ = rint(n*·ρ̂(d_birth,i))` through
  `rho_he_ratio(depth, steepness=drag_gate_steepness(cfg))` — the
  **same** erf-complement surface the drag/pickup/cooling gates share
  (function-local import of the steepness resolver: `ion.py` imports the
  seed module). Per-ion consequences exactly per H.3b:
  `mass_kg(0) = complex_mass_amu(n₀ᵢ)·U` (mass history/final ride it),
  the t0 `e_bind_pair(n₀ᵢ)` E_pot fold (the fold call takes
  `n0_initial` — scalar `ANCHOR_N_START` under `full`, byte-identical to
  the delivered call; vectorised array under `density_tied`), and
  `n_shell(0)` via the unchanged mass→n rint rule. **No checkpoint
  schema change** (v7 untouched — all arrays were already per-ion).
  **E_int(0) stays the constant `f_int·E_avail`** — the Σ-ratio coupling
  is Slice T6's p-law, not this arm.
- `scripts/tier2_common.py :: build_biphasic_cfg`: None-sentinel
  `initial_shell_model` kwarg (T7 pattern). No preset, generator, or
  production config selects the arm — the T9 leg-B generator stamps it.

**Verification.** `tests/test_initial_shell_model.py` — 16 focused
tests: config default/typo-guard/biphasic-only refusal/validating
density-tied cfg/back-compat (pre-T5 `cfg.json` loads with `full` —
the Slice-DS precedent)/round-trip; the dressed seed against an
**independent** hand oracle (`math.erf` formula, not `rho_he_ratio`);
per-ion fold and E_kin↔mass consistency by difference-oracle against
the `full` arm (E_pot is mass-independent, so the arms differ exactly
by the fold); the H.3b wiring oracles (center-pinned ⇒ n₀ = 21 exactly
and elementwise-identical seed; n₀ monotone non-increasing in birth
radius); the T5/T6 boundary (E_int(0) uncoupled); `build_biphasic_cfg`
pass-through. Byte-identity of the `full` default: neighboring suites
(ion_initial_state, biphasic_step, tier2_common, shed convention, birth
law, gen_tier2_md_confirmation) **178 passed** unchanged; **full suite
2404 passed**.

**CALIBRATION_MAP:** row 27 — the dressing law is an **arm, not a
knob** (parameter-free tied law; the shell-averaged-ρ̂ variant stays a
stated sensitivity, never a fit).

**Boundaries carried (H.3b):** first-order occupancy statement, no
shell-restructuring dynamics; the dressing↔pickup interplay
(under-dressed ions re-filling via the live Langmuir channel — a
channel the 1D twin does not have) is read from the leg-B runs and
reported as a twin-divergence candidate (S2c-P4 listed channel).

Next per the I.11 sequence: **T9 leg B** behind its own leg trigger —
the leg-B twin re-score at the dressed configuration (pre-registering
histograms *and* per-bin mean detected KE on the stated co-moving
basis, per the item-3 amendment), then the four dressed A/B dirs
against the certified apcm baseline. Nothing here discharges F5.

## T9 leg B TRIGGERED — A″ + T5 `density_tied`; twin re-scored at the dressed configuration; predictions PRE-REGISTERED before any MD (2026-07-18)

User trigger ("Go ahead with leg B run and compare to twin"), following
the T5 delivery entry above. Leg B flips exactly **one lever vs the
certified apcm baseline**: `initial_shell_model = "density_tied"`.

**Builds delivered under this trigger (both TDD, RED watched first):**
- Twin: `stage_legb` in `scripts/tier2_h2b_forward_model.py` — the
  leg-A′ stage with per-fragment
  `n_eject = clip(rint(n*·ρ̂(r₀ − R₂₀₀₀)), 0, n*)` at the shared
  `STEEP_A = 14.2` surface and chords integrated at the dressed birth
  mass `complex_mass_amu(n_eject)` (the H.2b `build_fragment_table`
  convention); same SEED/draws as `stage_legaprime` so the in-stage
  `b_undressed` anchor rows must equal the A′ twin rows — **verified
  exactly at m = 20000** (c1 0.033/0.089/6.22, c3 0.067/0.336/4.48 …
  the §4m twin column verbatim) and locked by the smoke test at
  m = 200. Usage + SystemExit mode lists updated (the F6 lesson).
  Wiring oracles re-verified after the edit (Σ(21) = 0.18783720,
  K = 0.74460 / 0.89767). Twin suite 16 passed.
- Generator: `LEG = "b"` in `gen_tier2_md_confirmation.py` — identical
  to `aprime_cm` except `initial_shell_model="density_tied"` stamped
  through the delivered kwarg; dirs `…_tier2probe_conf270_bc{1..4}`.
  Lock test asserts the field-diff vs A″ is exactly
  `{initial_shell_model}` (monkeypatched — the F8 lesson);
  `test_unknown_leg_rejected` rotated its example to `"z"`. Generator +
  common + T5 suites 93 passed.

**Twin re-score at the dressed configuration (m = 20000; outputs
`h2b_leg_b_predictions.csv` / `h2b_leg_b_ke.csv`).** Dressed
`n_eject`: mean 16.71 (margin-3 floor 13 — the I45 band). The dressing
is a **large lever on the suppressed/bare side** (smaller Σ(n₀) for
near-surface births ⇒ suppression easier at fixed E₀):

| config | supp (dressed / undressed) | n̄_det (dr/undr) | n₁ (dr/undr) | trapped (dr/undr) |
|---|---|---|---|---|
| c1 | **0.389** / 0.089 | 3.94 / 6.22 | 0.101 / 0.172 | 0.055 / 0.033 |
| c2 | **0.349** / 0.035 | 4.17 / 6.55 | 0.105 / 0.170 | 0.060 / 0.035 |
| c3 | **0.533** / 0.336 | 2.78 / 4.48 | 0.049 / 0.069 | 0.093 / 0.067 |
| c4 | **0.436** / 0.118 | 3.67 / 6.05 | 0.077 / 0.138 | 0.055 / 0.033 |

**Pre-registered predictions (all KE claims on the co-moving basis —
twin-native and the MD's A″/leg-B convention; same N = 50 bridge seed,
margin 3 Å, exclude policy, 8000 ps cap):**

- **BP-P1 (suppressed/bare jump — the headline).** MD suppressed
  fractions rise from the A″ 0.096/0.032/0.326/0.117 to ≈ the twin's
  0.389/0.349/0.533/0.436; the ordering stays c3 > c4 > c1 > c2 with
  compressed spacing. Directional caveat (listed divergence): MD may
  land **below** the twin — in-bubble pickup re-fills under-dressed
  shells (Σ grows before gate-open), a channel the twin does not have.
- **BP-P2 (histogram shape).** n̄_det drops to ≈ 3.9/4.2/2.8/3.7; the
  deep-shell tail truncates near the dressed n_eject band (essentially
  no detected weight above n ≈ 14–15, vs A″ weight to n ≈ 18+);
  W₁(MD, dressed twin) at the leg-A′/A″ scale (≲ 0.7 bins over 22).
- **BP-P3 (the KE axis).** Per-bin mean detected KE at n₁ lands ≈
  **0.549 / 0.512 / 0.176 / 0.499 eV** (c1..c4; vs A″ measured
  1.005/1.025/0.312/0.991 — the dressed n₁ occupants are deeper-born,
  slower); bare-class (suppressed, intact-complex) KE ≈
  **0.999 / 0.892 / 0.352 / 0.945 eV**; the whole small-n curve shifts
  down vs A″ per the twin table.
- **BP-P4 (trapped class).** Rises to twin 0.055/0.060/0.093/0.055;
  the MD read has run ~1.6× the twin's 150 ps chord read (A′
  precedent), so ≈ 6–14 ions/100 with c3 the largest.

**Listed twin-divergence channels (S2c-P4 — an unlisted one is a
model-structure finding):** (a) pickup re-filling after under-dressed
birth (direction: MD supp ≤ twin supp, MD deep bins ≥ twin); (b) the
MD dresses per **atom** at ion-t0 positions (±R0/2 chord offset +
neutral drift) vs the twin's molecule-center birth dressing (±1 rung
near the surface); (c) trapped-class dynamics beyond the twin's chord
read. MD pilots (`bc1..4`) launch next; scored against the dressed
twin, the A″ baseline, and these pre-registrations. Nothing here
discharges F5.

## T9 leg B EXECUTED — BP-P1..P4 all CONFIRMED: the T5 dressing transfers quantitatively (W₁ 0.41–0.51 bins, the chain's best); pickup re-filling measured at ≈ 7 ions/100; the bare-bin KE gap resolves into RQ3 bookkeeping; no unlisted divergence (2026-07-18)

Full record: findings **§4p + I62–I64**. Four dirs on disk
(`…_tier2probe_conf270_bc{1..4}`, all five artifacts; bc1 from the
first launch, bc2 regenerated after two background-task kills — the
partial dir removed per the generator's documented recovery path each
time, bc3/4 fresh; final run detached with artifact monitoring; resume
guard exercised as designed).

**Verdicts vs the pre-registered twin numbers (same-day entry above;
all KE on the co-moving basis on both sides):**
- **BP-P1 CONFIRMED, caveat fired as stated** — supp→bare
  0.322/0.267/0.464/0.378 (A″ 0.096/0.032/0.326/0.117; twin
  0.389/0.349/0.533/0.436): the class triples, ordering
  c3 > c4 > c1 > c2 exact, and the uniform 0.06–0.08 deficit vs the
  twin is the listed pickup-re-filling channel — **first direct
  measurement of the dressing↔pickup interplay** (the H.3b boundary;
  I63).
- **BP-P2 CONFIRMED** — n̄_det 4.13/4.40/3.08/3.97; tail truncates at
  the dressed band (weight ends n ≈ 12–14 vs A″ n ≈ 18+);
  **W₁(B, twin_b) = 0.41–0.51 bins, the best twin↔MD agreement of the
  oracle chain**, with the lever's own move 2–4× larger
  (W₁(B, A″) = 1.0–1.9).
- **BP-P3 CONFIRMED (solvated); bare bin resolved** — n₁ KE
  0.661/0.626/0.263/0.601 eV onto the twin's 0.549/0.512/0.176/0.499
  (×1.2; c3's ×1.5 is a 6-ion bin); the apparent ×1.6 bare-bin excess
  collapses to **×1.08–1.12** under the co-moving break-up rescale
  m_I/m_complex (MD reports the intact dressed complex at
  m̄ ≈ 187–191 amu, the twin books the RQ3 spec-b bare fragment) —
  the §4n convention-decided fragmentation read, now with numbers
  (I64). The experimental bare-bin mean KE is a direct fragmentation-
  convention discriminator.
- **BP-P4 CONFIRMED** — trapped 10/10/16/10 per 100 (twin
  5.5/6.0/9.3/5.5 → the A′-precedent ×1.6–1.8 chord-to-MD factor; c3
  largest both sides; dressing raises the class vs A″ 6/7/14/6 in the
  predicted direction).

Seed wiring: MD dressed n₀ q05 13.9 / mean 17.18 vs twin 13.0 / 16.71
(the per-atom ±R0/2 offset ≈ +0.5 He, one atom at n₀ = 12 below the
molecule-center margin floor — listed channel (b), measured small).
**S2c-P4 holds: no unlisted divergence.**

**Consequences:** the Step-2c geometry-closure chain now has legs
A → A′ → A″ → B all landing their pre-registrations — position axis,
shed convention, and shell dressing each validated as one-lever flips
with the twin quantitative at every step. The ensemble-side suppressed/
bare weight problem (OQ-B's Tier-2 face) is now expressible in MD at
production kinematics: the dressed geometry puts 27–46 % in the
bare-candidate class at these knob values, bracketing the experimental
43.5 % — **reported, not adjudicated** (the C-matrix knobs were located
under the 1D ensemble; the re-centered re-pilot after T6/T8 scores the
experimental targets). Next per the I.11 sequence: **Slice T6**
(E₀–dressing coupling p-law; needs T5 live — now is) behind its own
trigger, then T9 legs B/C at fixed N, T8, leg D. The E2 Landau-gated
build remains adjudicated-but-unbuilt, required before any N = 500 run.
Nothing here discharges F5.

## Post-leg-B documentation reconciliation — bare-bin reframe reaffirmed as the standing scoring convention; the (n, mean-KE) data contract is SATISFIED; endgame scoring adopts the committed error model (2026-07-18)

Discussion outcome (post-leg-B; **no MD, no code, no schema change** —
records decisions and reconciles stale doc language so the reframe need
not be re-derived next session):

1. **The bare bin (n = 0) is renormalized out of the arbitration target
   for every remaining leg (C, D) and the re-pilot.** RQ8's standing
   verdict (2026-07-11: bare is a channel-branching quantity, not a model
   target) is reaffirmed by leg B (findings §4p / I62–I64): the dressed
   suppressed/bare class brackets 43.5 % *coincidentally*, and its detected
   KE is a fragmentation-convention bookkeeping read (co-moving vs
   momentum-conserving break-up), not the experimental bare channel. The
   biphasic arbitration surface is the **solvated branch**: the n ≥ 1
   renormalized histogram + the per-n mean-KE curve. Leg-C/T6 and leg-D
   pre-registrations state this explicitly (plan §I.11 T9 + §I.11.2 NBs).

2. **The (n, mean-KE) data-contract prerequisite is SATISFIED, not
   pending.** The committed reference is
   `data/reference/ihe_ked/IHe_KED_reference.csv` (per-n `meanKE_eV` =
   first moment of the 3-D P(E); `COLUMNS.md` +
   `IHe_KED_reference.provenance.json`; full n0/n1 spectra + trusted 3-D
   curves n0–n4). The RQ8 / §H.2b-D7 language calling the export "pending"
   or "the export prerequisite stands" was stale — corrected in
   `RESEARCH_QUESTIONS.md` RQ8.

3. **Committed bare meanKE = 3.706 eV** (mode 4.758), superseding the stale
   ~2.9 eV quote in RQ8/§H.2b-D7. That is ≈ 1.4× the 2.71 eV per-fragment
   ballistic ceiling of the 2.70 eV solvated channel — the foreign-channel
   verdict is on *firmer* committed ground than the stale number implied.
   Committed n = 1 meanKE = 1.302 eV (matches the value used throughout;
   the D7 "0.974" was stale), → 0.066 eV at n = 17.

4. **Endgame (T9) scoring adopts the committed error model as the
   tolerance** (recommended; confirm at the T9 scoring build): per-point
   √(statErr² + sysErr²) for point-to-point scatter, plus the calib (4 %)
   and condition (6 %) fractional bands as two *correlated* whole-curve
   shifts — superseding the ad-hoc flat ×1.25 bar (retained as a coarse
   fallback). `noiseLimited = 0` on all 18 fragments (every value a genuine
   measurement). Plan §I.11 T9 scoring NB added.

Also noted (not a change): the pickup-re-filling channel (leg B I63,
≈ 7 ions/100) makes the 1D twin a **basin-locator, not a predictor**, at
the endgame — the final knobs are MD-located in the re-pilot.

Next per the I.11 sequence: **Slice T6** behind its own
`[PROCEED TO IMPLEMENTATION]` trigger, then T9 leg C. Nothing here
discharges F5.

## Slice T6 DELIVERED — `internal_energy_partition_law ∈ {constant (byte-inert default, p = 0), sigma_proportional (p = 1)}`: the D2 E_int(0)-dressing coupling; onset × (Σ(n₀)/Σ(n*))^p in the biphasic seed; CALIBRATION_MAP row 28 (2026-07-18)

Delivered under its own `[PROCEED TO IMPLEMENTATION]` trigger (TDD, RED
watched first — ImportError on the missing
`check_internal_energy_partition_config` guard). Exactly one new physics
lever, parameter-free, composing on the certified leg-B baseline.

**As built.**

- **Physics** (`physics/internal_energy_budget.py`): new pure helper
  `sigma_partition_factor(n0, *, law, picture, kappa, ladder)` returning the
  dimensionless p-law factor `(Σ(n0)/Σ(n*))^p`. `constant` returns the
  multiplicative identity `1.0` exactly (ignores n0/ladder — the byte-inert
  p = 0 path); `sigma_proportional` returns `ladder_cumsum(n0)/ladder_cumsum(
  N_STAR)` (Slice L's Σ at the cfg picture/kappa/ladder — the *same* Σ the
  §U self-unbound floor and the K-cooling asymptote consume, no duplicate
  physics). Vectorised scalar→float / array→ndarray; fail-loud on an
  unrecognised `law` (CLAUDE.md principle 4). Inert at n0 = n* (ratio 1)
  under either law.
- **Config** (`config.py`): `InternalEnergyPartitionLaw` Literal + field
  `internal_energy_partition_law = "constant"` (byte-inert default, beside
  the Slice-T5 `initial_shell_model` block); guard
  `check_internal_energy_partition_config` wired into `SimConfig.validate`
  (typo reject via the shared `_reject_unknown_enum`; `sigma_proportional`
  refused outside `mass_scenario='biphasic'` — the T5/T7 no-silent-inert
  convention, since the S2 onset is deposited only by the biphasic column-0
  seed).
- **Seed** (`simulation/ion_initial_state.py`): the T5/T6-boundary comment is
  retired — the S2 onset now multiplies `e_int_onset_eV(...)` by
  `sigma_partition_factor(n0_initial, law=..., picture=..., kappa=...,
  ladder=resolve_ladder(...))`, riding the **same** per-ion t0 shell
  `n0_initial` (scalar 21 under `full`, per-ion dressed under `density_tied`)
  and resolved ladder as the E_pot binding fold. E_int(0) is the **only**
  seed quantity T6 touches — mass, n_shell, and the E_pot fold are unchanged
  (Slice T5's). Byte-inert default: under `constant` the factor is a scalar
  `1.0`, so `E_int_eV[:, 0]` is byte-identical to the delivered onset.
- **Pass-through** (`scripts/tier2_common.py`): `build_biphasic_cfg` gains the
  None-sentinel `internal_energy_partition_law` kwarg (the T5/T7 pattern); no
  preset or generator selects the arm — the T9 legs stamp it (leg C flips
  p = 0 → 1 vs the certified leg-B baseline).
- **Tests** (`tests/test_internal_energy_partition_law.py`, 20): the pure
  factor (constant-unity scalar/array + n0-agnostic; `sigma_proportional`
  scalar/array hand oracles against an independent `ladder_cumsum`; inert at
  n*; strictly increasing in n0 with under-dressed < 1; unknown-law reject);
  config default/guards/back-compat (pre-T6 `cfg.json` loads `constant`)/
  round-trip; the coupled seed (constant arm = delivered onset even under a
  live T5 dressing; `sigma_proportional`+`full` inert; `sigma_proportional`+
  `density_tied` per-ion hand oracle with under-dressed ions getting less
  onset; the law touches E_int **only**, not mass/n_shell/E_pot fold); the
  pass-through pair.

Full suite **2426 passed** (byte-inert default locked by the delivered
suites; the T5 boundary test `test_E_int_onset_not_coupled_to_dressing`
stays green — the default `constant` arm keeps that boundary). CALIBRATION_MAP
row 28 (arm, not knob). The p = 1 prior (§4j finding 1: p = 0 over-suppresses)
is **not** hard-wired — the arm is swept in the T9 A/B chain.

Next per the I.11 sequence: **T9 leg C** (flip `internal_energy_partition_law
= sigma_proportional` vs the certified leg-B `bc` baseline, pre-registered
against the 1D twin re-scored on the solvated branch), behind its own leg
trigger. The E2 Landau-gated drag arm (§I.11.2 item 2, arm (c)) stays
adjudicated-but-unbuilt — required before any N = 500 run, not before the
N = 50 legs C/D. Nothing here discharges F5.

## T9 leg C TRIGGERED — B + T6 `sigma_proportional`; twin re-scored at p = 1; predictions PRE-REGISTERED before any MD (2026-07-18)

User trigger ("Go ahead with the leg-C pre-registration"), following the
T6 delivery entry above. Leg C flips exactly **one lever vs the certified
`bc` baseline**: `internal_energy_partition_law = "sigma_proportional"`
(p = 1) — cumulative on leg B (the `density_tied` dressing and the
co-moving shed convention carry over).

**Builds delivered under this trigger (both TDD, RED watched first):**
- Twin: `stage_legc` in `scripts/tier2_h2b_forward_model.py` — the leg-B
  dressed stage with the E₀-coupling exponent flipped `0 → 1`
  (`fate_map` applies `E0_i = E0·(Σ(n_eject)/Σ(21))^p` — **the same closed
  form as the MD's `internal_energy_budget.sigma_partition_factor`**, so
  the A/B is valid). Same SEED/draws as `stage_legb`, so the in-stage
  `c_p0` rows (dressed, p = 0) must equal the leg-B `b` rows — **verified
  exactly at m = 20000** (supp 0.389/0.349/0.533/0.436, n̄ 3.94/4.17/2.78/
  3.67, n₁ KE 0.549/0.512/0.176/0.499 — the §4p leg-B twin verbatim) and
  locked by the smoke/oracle test at m = 200. Usage + SystemExit mode
  lists updated. Twin suite 17 passed.
- Generator: `LEG = "c"` in `gen_tier2_md_confirmation.py` — cumulative on
  leg B (`density_tied` + co-moving carry over) plus
  `internal_energy_partition_law="sigma_proportional"` stamped through the
  delivered T6 kwarg; dirs `…_tier2probe_conf270_cc{1..4}` (distinct from
  the T3 `c{1..4}`). Lock test asserts the field-diff vs leg B is exactly
  `{internal_energy_partition_law}` (monkeypatched). Generator suite 19
  passed.

**Twin re-score at p = 1 (m = 20000; outputs `h2b_leg_c_predictions.csv` /
`h2b_leg_c_ke.csv` under the gitignored `data/runs/h2b_forward_model/`;
the pre-registered numbers are recorded here).** The p = 1 coupling lowers
the onset for under-dressed births (E0_i = E0·Σ(n₀)/Σ(21) ≤ E0), so it
**de-suppresses the over-suppressed side** — the §4j "p = 0 over-suppresses"
direction, now quantified against leg B (= c_p0):

| config | supp (p0 → p1) | n̄_det (p0 → p1) | n₁ frac (p0 → p1) | n₁ KE (p0 → p1) | trapped |
|---|---|---|---|---|---|
| c1 | 0.389 → **0.111** | 3.94 → **4.98** | 0.101 → **0.212** | 0.549 → **0.998** | 0.055 |
| c2 | 0.349 → **0.051** | 4.17 → **5.31** | 0.105 → **0.217** | 0.512 → **0.964** | 0.060 |
| c3 | 0.533 → **0.343** | 2.78 → **3.77** | 0.049 → **0.084** | 0.176 → **0.253** | 0.093 |
| c4 | 0.436 → **0.143** | 3.67 → **4.88** | 0.077 → **0.172** | 0.499 → **0.973** | 0.055 |

The suppressed ordering **c3 > c4 > c1 > c2 transfers** (0.343 > 0.143 >
0.111 > 0.051); trapped is p-invariant (the p-law never touches the chord
dynamics, only the onset).

**Pre-registered predictions (all KE claims on the co-moving basis —
twin-native and the MD's leg-B/C convention; same N = 50 bridge seed,
margin 3 Å, exclude policy, 8000 ps cap; scoring is the solvated branch,
bare renormalised out — post-leg-B decision):**

- **CP-P1 (de-suppression — the headline).** MD `cc` suppressed fractions
  drop from the leg-B `bc` measured 0.322/0.267/0.464/0.378 toward the
  twin p = 1 **0.111/0.051/0.343/0.143**, ordering c3 > c4 > c1 > c2.
  Directional caveat (listed divergence, as leg B): in-bubble pickup
  re-filling means the MD may land **above** the twin (Σ grows before
  gate-open → fewer stay self-unbound; leg B measured ≈ +7 ions/100 of this).
- **CP-P2 (solvated histogram).** n̄_det rises to ≈ **4.98/5.31/3.77/4.88**
  (from bc ≈ 4.13/4.40/3.08/3.97); weight shifts off bin 0 into the shallow
  solvated bins (n₁ ≈ doubles); the deep tail still truncates near the
  dressed n_eject band (top weight to n ≈ 14–15). W₁(MD, twin) expected at
  the chain scale (≲ 0.5 bins on the solvated branch, cf. leg B 0.41–0.51).
- **CP-P3 (KE axis).** Per-bin n₁ mean KE **rises** to ≈
  **0.998/0.964/0.253/0.973 eV** (from bc 0.661/0.626/0.263/0.601) — the
  de-suppressed n₁ occupants ride shallower descents; this lands **near the
  A″ undressed** (1.005/1.025/0.312/0.991) and **closer to the committed
  experimental n₁ = 1.302 eV** than either leg B or the undressed A″. The
  whole solvated KE curve shifts up modestly. Bare-bin KE stays RQ3/RQ8
  bookkeeping (renormalised out).
- **CP-P4 (trapped class).** p-invariant: ≈ 0.055/0.060/0.093/0.055 twin /
  ~1.7× the 150 ps chord read in MD (the A′/leg-B precedent), c3 largest.

**Listed twin-divergence channels (S2c-P4 — an unlisted one is a
model-structure finding):** (a) pickup re-filling after under-dressed birth
(MD supp ≥ twin supp; MD deep bins possibly below twin); (b) per-**atom**
ion-t0 dressing vs the twin's molecule-center birth dressing (±1 rung near
the surface); (c) trapped-class dynamics beyond the twin's 150 ps chord
read. MD pilots (`cc1..4`) launch next under this trigger; scored against
the p = 1 twin, the certified `bc` baseline, and these pre-registrations on
the solvated branch. Nothing here discharges F5.

## T9 leg C EXECUTED — CP-P1..P4 all CONFIRMED: the T6 p-law de-suppresses to the twin (W₁ 0.31–0.47, the chain's best); n̄ + solvated KE transfer; n₁ KE moves toward experiment; no unlisted divergence (2026-07-19)

The four `cc{1..4}` pilots ran to completion (all five artifacts each) and
were scored against the pre-registered p = 1 twin and the certified `bc`
baseline (findings §4q + I65–I67). One lever vs `bc`
(`internal_energy_partition_law = "sigma_proportional"`).

**Execution accommodation (session-specific, physics-neutral — user-approved
2026-07-19).** Each config's E2 relaxation compute (~15–20 min) plus the
~380 MB compressed `relaxation.npz` save exceeded this session's ~25-min
process window (two full-pipeline background runs were terminated mid-save,
leaving a truncated npz). Resolution: a scratchpad per-stage-resume driver
(`run_legc_resume.py`) that (a) resumes each config from its last completed
artifact and (b) shrinks the relaxation checkpoint byte budget so the save
is instant. Verified scored-read-neutral **before** running: detection seeds
from the terminal column only (`detection_stage.py:396`) and the stage always
stores the true final state last (`relaxation_stage.py:439–440`), so the
coarser stride leaves every CP-P1..P4 read bit-exact and touches only the
unused intermediate trajectory (I67). The `cc` `relaxation.npz` are ~2 MB vs
the `bc` dirs' ~380 MB; zero repo-code change (the driver composes delivered
stage functions). All four completed in one pass with the small checkpoint.

**Scored A/B/twin (droplet_retained excluded; suppressed→bin 0; W₁ over 0–21;
KE co-moving):**

| config | supp (C / twin_c / B) | n̄_det (C / twin / B) | n₁ (C / twin) | W₁(C, twin) | n₁ KE (C / twin / B) | trapped (C / twin) |
|---|---|---|---|---|---|---|
| c1 | **0.101** / 0.111 / 0.322 | 4.87 / 4.98 / 4.13 | 0.191 / 0.212 | **0.44** | 1.077 / 0.998 / 0.661 | 0.110 / 0.055 |
| c2 | **0.022** / 0.051 / 0.267 | 5.27 / 5.31 / 4.40 | 0.202 / 0.217 | **0.47** | 1.104 / 0.964 / 0.626 | 0.110 / 0.060 |
| c3 | **0.301** / 0.343 / 0.464 | 3.76 / 3.77 / 3.08 | 0.108 / 0.084 | **0.31** | 0.317 / 0.253 / 0.263 | 0.170 / 0.093 |
| c4 | **0.124** / 0.143 / 0.378 | 4.87 / 4.88 / 3.97 | 0.146 / 0.172 | **0.45** | 1.089 / 0.973 / 0.601 | 0.110 / 0.055 |

**Verdicts (pre-registered CP-P1..P4):**
- **CP-P1 CONFIRMED (de-suppression headline).** p = 1 drops supp to ≈ a
  third of `bc`, landing on the twin (ordering c3 > c4 > c1 > c2 exact). The
  re-filling channel fires in the **corrected** direction: MD ≈ 1–3 ions/100
  **below** twin ("fewer stay self-unbound" — the TRIGGERED "above"/"≥" was a
  wording slip; the parenthetical mechanism was right, and leg B measured the
  same MD-below-twin sign at ≈ 7/100).
- **CP-P2 CONFIRMED.** n̄_det matches the twin to ≤ 0.11 He; **W₁(C, twin) =
  0.31–0.47 — the best twin↔MD histogram agreement of the oracle chain**
  (A′ 0.55–0.69 / A″ 0.51–0.68 / B 0.41–0.51).
- **CP-P3 CONFIRMED.** n₁ KE rises to 1.08/1.10/0.32/1.09 eV (from `bc`
  0.66/0.63/0.26/0.60) — ≈ A″ undressed and the closest of any leg to the
  committed experimental n₁ = 1.302 eV; MD ×1.08–1.15 above twin (small
  composition residual). Bare-class KE (1.99/1.95/0.69/1.94, intact dressed
  complex) stays RQ3/RQ8 bookkeeping, renormalised out.
- **CP-P4 CONFIRMED.** Trapped p-invariant (0.11/0.11/0.17/0.11 ≈ `bc`), c3
  largest, ≈ 1.8× the twin's 150 ps chord read.
- **S2c-P4 holds:** every deviation is a listed channel (re-filling; per-atom
  dressing; trapped dynamics). No unlisted divergence.

Next per the I.11 sequence: **Slice T8** (the droplet-prior D4 family; audit
the legacy `use_single_droplet_size=False` machinery first, then wire the
Kornilov log-normal + pickup-weighted variant behind config), then **T9 leg D**
(+ droplet prior, fixed-N legs A–C carried) → the re-pilot at the MD-located
knobs → the T4/ihe_ked solvated-branch scoring. The E2 Landau-gated drag arm
(§I.11.2 item 2, arm (c)) stays adjudicated-but-unbuilt — required before any
N = 500 run. Nothing here discharges F5.

## E2 Landau-gated drag arm DELIVERED — §I.11.2 item 2, arm (c); byte-inert `relaxation_dissipation` enum; full suite 2443 passed (2026-07-20)

`[PROCEED TO IMPLEMENTATION]` given for the adjudicated arm (c). Built ahead of the
pre-registered I.11 next-step (Slice T8) at the user's call — legitimate because the
arm is **byte-inert** (default preserves every delivered leg) and only needs to exist
before the first `N = 500` run. Two user sign-offs before the trigger: (i) gate on
**speed** with the threshold pinned at the existing legacy `v_limit = 40 m/s`
(0.4 Å/ps) — re-pinning is a separate domain-expert calibration; (ii) `b` drawn from
the config's production drag bundle (single source, no hard-wired constant). TDD, RED
watched first.

**Physics (frozen).** Replaces the E2 coulomb-translate's `_zero_gamma` with a
speed-gated friction coefficient:

```
γ(v, depth) = 0                      for v ≤ v_L      (sub-Landau superfluid, frictionless)
γ(v, depth) = g(depth)·b·v²          for v >  v_L      (the locked pure cubic)
```

`v_L = cfg.v_limit_angstrom_per_ps` (0.4 Å/ps), `b = 2.5153509 amu·ps/Å²` and the erf
spatial gate `g(depth)` both from `cfg.drag_coefficients` via `physics.drag.drag_gamma`
(rule-1 single source). Units balance both branches (`b·v² = amu/ps`). The gate is on
**speed**, not KE: the Landau critical velocity is mass-independent, which is also the
choice that preserves `drag.py`'s mass-agnostic contract (a KE gate would force `m`
into the drag evaluation). Heaviside is strict `>` — frictionless exactly at `v_L`.

**Surface (byte-inert).**
- `config.py`: `RelaxationDissipation = Literal["zero_gamma", "landau_gated_drag"]`;
  field `relaxation_dissipation = "zero_gamma"` (default) beside `relaxation_forces`.
  `check_relaxation_config` gains the typo-reject (`_reject_unknown_enum` against
  `_KNOWN_RELAXATION_DISSIPATIONS`) **and** the no-silent-inert pairing guard —
  `landau_gated_drag` is refused under `relaxation_forces != "coulomb"` (free_flight is
  ballistic, no drag O-step; the T5/T7 `sigma_proportional`-refused-outside-biphasic
  precedent). Back-compat: an old `cfg.json` missing the key loads `zero_gamma`
  (dataclass default, `RunDirectory.load_cfg`'s `SimConfig(**payload)` path).
- `simulation/relaxation_stage.py`: `_make_relaxation_gamma_fn(cfg, gate_steepness)`
  returns the delivered `_zero_gamma` **object** for the default (translate path
  literally unchanged) or the Landau closure `np.where(v > v_L, drag_gamma(...), 0.0)`.
  `_coulomb_translate` takes the `gamma_fn`; `run_relaxation_stage` builds it once from
  `relax_cfg` next to `gate_steepness`. `_free_flight_translate` unchanged (guard
  forbids the landau+free_flight pairing). No change to `drag.py`, `drag_form`, the main
  ion driver, the checkpoint schema, or the RNG draw order.

**Invariant.** The BAOAB step already threads its drag `dE_dissip` into `E_dissip`
(`ion_propagation_step.py:395`) on top of the K2 cooling drain that `biphasic_step`
books (`:630`) — the exact machinery that closes the ion-stage ledger — so the
**5-term invariant closes with drag on**, verified on the real-driver seed
(`test_five_term_invariant_closes_under_landau_arm`). Scope: the marginal
droplet-retained class only; the ejected read is arm-invariant (`g → 0` outside).

**Verification (7-step order).** New suites `tests/test_relaxation_dissipation_config.py`
(6) + `tests/test_relaxation_landau_gamma.py` (9): closed-form γ both sides of `v_L`,
boundary-belongs-to-below, spatial-gate-off-outside, array elementwise, above-cutoff ==
`drag_gamma`; byte-inert identity (`fn is _zero_gamma`); below-cutoff run == zero_gamma
bit-for-bit; above-cutoff isolates the drag channel (n/E_int identical, strictly more
`E_dissip`, colder); real-seed invariant closure. The delivered relaxation suite (36,
incl. `test_e_dissip_equals_cumulative_k2_drain`) stays green — the byte-inert oracle.
**Full suite: 2443 passed** (233 s). CALIBRATION_MAP row 29 (arm, not knob); `v_limit`
gains a physics-live reader beyond the superseded collision path (retires that rule-2
carry's guard-only status for the relaxation stage).

Next per the I.11 sequence is still **Slice T8** (droplet prior) → T9 leg D → re-pilot →
T4/ihe_ked scoring; the arm now stands ready for the first `N = 500` run. `v_L`
re-pinning (domain-expert Landau critical velocity) is the one deferred calibration.
Nothing here discharges F5.

## Post-E2 code review EXECUTED — high-effort review of `1e4b91d..HEAD` (T5 + leg B, T6/T7 + leg C, E2 Landau arm); 4 verified findings + 2 cosmetics fixed, byte-inert defaults preserved, full suite 2445 passed (2026-07-20)

High-effort multi-agent review (4 finder angles + a per-location adversarial
verify pass) of everything delivered after the `1e4b''` post-A″ review: Slice
T5 + leg B, Slice T6 + leg C, and the E2 Landau-gated drag arm. 6 verified
candidates, 1 refuted (the "legb/legc/legaprime are copy-paste" dedup — refused:
they are deliberately-explicit pre-registration oracle stages). All fixed; no
new physics, no schema/RNG/constants touch, every delivered leg stays
byte-identical (the changed paths are the *default* arms only).

**F1 (correctness — the headline). The relaxation freeze-all early-exit
truncated the Landau drag arm for exactly the class it targets.** The loop broke
at `freeze_flags.all()` (`_freeze_mask` = `n ≤ 0 | E_int < D_0`, a pure
*evaporation* freeze). Drag touches neither `n` nor `E_int`, so the freeze
instant is drag-independent — but `landau_gated_drag` exists to damp the
droplet-retained ions' residual super-Landau oscillation, which continues *past*
the evaporation freeze. A retained ion frozen while still oscillating above `v_L`
kept the KE the arm is meant to remove (latent: the arm is byte-inert default and
has not run production, so no committed result moved — but it defeats the arm
ahead of the N = 500 run it gates). Fix: extracted
`_relaxation_converged(state, freeze_flags, *, dissipation, v_limit)` —
`zero_gamma`/`free_flight` still exit on the freeze alone (byte-identical), the
drag arm additionally requires all speeds ≤ `v_limit`. Regression:
`test_relaxation_converged_requires_sub_landau_under_the_drag_arm` (unit) +
`test_landau_arm_damps_past_the_evaporation_freeze` (end-to-end; `time_relaxed`
≫ one step vs the `zero_gamma` arm's one-step stop). Confirmed
`test_above_cutoff_dissipates_and_isolates_the_drag_channel` is undisturbed (its
window hits no preempting freeze, so the arms stay same-length there).

**F2 (correctness — script). Unguarded `w_tot == 0` in the leg-B/leg-C twin
stages wrote NaN into the pre-registered prediction CSVs.** `stage_legb` /
`stage_legc` normalise the histogram + fractions by `w_tot = w.sum()` (non-trapped
weight); an all-trapped config gives `0/0 → NaN`, silently corrupting the
pre-registration record rather than raising. Fix: a `w_tot == 0` guard that raises
in both new stages. `stage_legaprime`'s identical pre-existing copy left untouched
(out of the review range).

**F3 (cleanup/perf). `base_gamma` was evaluated for every ion via `np.where`**,
including the sub-Landau majority (the retained population the arm is about) whose
coefficient is discarded. Fix: masked assignment — the erf-gated cubic is computed
only where `speed > v_limit`; numerically identical to the old `np.where` on the
array contract.

**F4 (cleanup/DRY). The "typo-reject + no-silent-inert pairing" guard pattern was
triplicated** across `check_initial_shell_config`, `check_internal_energy_partition_config`,
and the dissipation block of `check_relaxation_config`. Fix: a shared
`_require_pairing(*, field, value, trigger, dep_field, dep_value, actual, reason)`
helper (mirrors the existing `_reject_unknown_enum` idiom); each guard passes its
own physics `reason`, so the three cannot drift apart. Pinned message substrings
(field names, `biphasic`/`landau_gated_drag`) preserved — the guard suites stay
green.

**Cosmetics.** `_KNOWN_RELAXATION_DISSIPATIONS` moved *before* `check_relaxation_config`
(the other `_KNOWN_*` tuples precede their guards). `sigma_partition_factor`
docstring: "picture/kappa ignored" made conditional on table injection (only true
when an explicit ladder table is injected; the analytic ladder consumes them).

Also cleared on independent inspection (no fix needed): `n0_initial` is always
defined on the biphasic seed path (both sub-branches of the wrapping `elif` bind
it); the MD seed's missing `np.clip(n0, 0, N_STAR)` is harmless because
`rho_he_ratio ∈ [0, 1]` saturates at exactly `1.0`, so `rint(21·ρ̂) ≤ 21` and never
indexes past the 21-rung ladder (the twin's clip is defensive only).

**Verification.** Files changed: `config.py`, `physics/internal_energy_budget.py`,
`simulation/relaxation_stage.py`, `scripts/tier2_h2b_forward_model.py`,
`tests/test_relaxation_landau_gamma.py` (+2 tests). Narrow→broad: Landau (11),
relaxation, the three guard suites, forward-model/partition/generator (56), then
**full suite 2445 passed** (2443 baseline + 2 new; 230 s). No F-register item is
discharged or opened; F5 stands. Next per the I.11 sequence is unchanged:
**Slice T8** (droplet prior) → T9 leg D → re-pilot → T4/ihe_ked scoring.

## Slice T8 pre-build decision block OPENED — rule-1 audit executed; the slice re-scopes downward (plumbing exists end-to-end); five decisions T8-D1..T8-D5 await adjudication (2026-07-20)

The rule-1 audit demanded by the T8 spec (§I.10 / §I.11) ran read-only, zero
MD. Headline findings (full detail + the decision register: plan **§I.11.3**):

1. **The legacy `use_single_droplet_size=False` machinery is already ported
   and live** — `sampling/droplet_sizes.py` is the literal
   `generate_droplet_sizes.m` port (log-normal source + full pickup-cell MC),
   consumed at `initial_state.py:72–78`. δ = 0.625 is a module constant, and
   the source correlation at the droplet-distribution preset's conditions
   gives ⟨N⟩ ≈ 1.86 k — close to but **not** the D4 pin of 2000.
2. **The pipeline is per-ion droplet-radius-aware end to end** (both birth
   laws; every gate depth in the ion step; the T5 dressing seed; relaxation +
   detection via the pass-through `droplet_radii_angstrom` field). F.4's
   feared "new wiring surface" does not exist and **no checkpoint schema
   change is needed** — T8's build surface is the prior selection only.
3. **The twin already implements the D4 family** (H.2b: uniform-in-lnN
   proposal on [250, 16000]; µ_ln = ln 2000 − δ²/2; pickup variant
   = × N^(2/3) un-re-centered, realized mean ≈ 2.6 k) — `stage_legd` is a
   near-free reuse.
4. **The analytic D4 family ≠ the legacy pickup MC** (single-pickup selection
   + evaporation distort the log-normal); twin↔MD parity requires the MD arm
   to sample the analytic family; the MC stays as the legacy arm.

Open decisions (recommendations in §I.11.3): **T8-D1** enum surface
(`droplet_size_prior ∈ {legacy (byte-inert default), kornilov_lognormal,
pickup_weighted_lognormal}` + ⟨N⟩/δ fields, `_require_pairing` against
`use_single_droplet_size=True`); **T8-D2** family pinned to the twin's
convention verbatim incl. the [250, 16000] truncation (logged mass);
**T8-D3** exact draw on the MD side (pickup family is itself log-normal with
µ_ln′ = µ_ln + (2/3)δ²); **T8-D4** leg-D = one lever vs leg C at the δ = 0.625
primary, N = 50, `zero_gamma`, sensitivity family deferred to the re-pilot;
**T8-D5** `stage_legd` twin stage before any MD. Build stays behind its own
`[PROCEED TO IMPLEMENTATION]`; nothing here discharges F5.

## Slice T8 DELIVERED — `droplet_size_prior ∈ {legacy (byte-inert default), kornilov_lognormal, pickup_weighted_lognormal}`: the analytic D4 prior family, exact-draw on the twin's window; twin `stage_legd` with a bit-exact leg-C oracle; CALIBRATION_MAP row 30; full suite 2464 passed (2026-07-20)

T8-D1..T8-D5 adjudicated by the user **as recommended** (same day as the
decision block) and `[PROCEED TO IMPLEMENTATION]` given. TDD, RED watched
first (ImportError/AttributeError on every new surface, then the two
convention fixes: guards fire at config-load via the check function — the E2
test precedent, not at construction — and the checkpoint stores per-**atom**
(2N) radii).

**Physics (frozen — T8-D2/D3, the twin's convention verbatim).**

```
kornilov_lognormal:        lnN ~ Normal(µ_ln, δ),  µ_ln = ln⟨N⟩ − δ²/2
pickup_weighted_lognormal: lnN ~ Normal(µ_ln + (2/3)δ², δ)
```

both truncated to the twin's proposal window **[250, 16000] He** by exact
inverse-CDF (`z = Φ⁻¹(Φ(a) + u·(Φ(b)−Φ(a)))` — one uniform per molecule, no
importance weighting, no rejection). The pickup arm is the exact ln-normal
identity `N^(2/3)·LogNormal(µ_ln) ∝ LogNormal(µ_ln + (2/3)δ²)`, deliberately
**not** re-centered (realized mean ⟨N⟩·e^(2δ²/3) ≈ 2.6 k at δ = 0.625 — the
shift *is* the pickup physics). Defaults ⟨N⟩ = 2000.0 (the D4 pin), δ = 0.625
(Kornilov 2009). R(N) stays the single-source `droplet_radius_bulk_angstrom`.

**Surface (byte-inert).**
- `config.py`: `DropletSizePrior` Literal + window constants
  `DROPLET_PRIOR_N_LO/HI` (rule-1 pinned to the twin's `N_LO/N_HI` by test);
  fields `droplet_size_prior="legacy"`, `droplet_prior_mean_N=2000.0`,
  `droplet_prior_delta=0.625`; `check_droplet_prior_config` wired into the
  validate chain — typo reject + the **both-direction** no-silent-inert
  guard: analytic arms refused under `use_single_droplet_size=True` (shared
  `_require_pairing`), off-default family params refused under `legacy` (the
  T7 margin-under-boltzmann precedent), δ > 0 and ⟨N⟩-inside-window under the
  analytic arms (no-silent-caps). Back-compat: old `cfg.json` without the
  keys loads the defaults (dataclass-default path, field-default test).
- `sampling/droplet_sizes.py`: `sample_droplet_sizes_analytic` (exact draw,
  RuntimeWarning when the truncation discards > 5 % of the family's mass) +
  `analytic_prior_truncated_mass` (closed form: ≈ 0.14 % at δ = 0.625,
  ≈ 1.5 % at δ = 0.80 — the provenance record). Legacy MC untouched.
- `simulation/initial_state.py`: three-way dispatch — non-`legacy` arms call
  the analytic sampler; under the default the branch pair is literally the
  pre-T8 code (monkeypatch test proves the analytic sampler is never touched).
- `scripts/tier2_h2b_forward_model.py`: `stage_legd` (mode `legd`) — same
  SEED and u/µ draw order as `stage_legc`, prior drawn *after* the legc
  stream via the **package sampler** (twin and MD share the prior by
  construction); `d_delta` (delta N = 2000) is the in-stage wiring oracle,
  `d` scales the same u/µ draws to the per-molecule radii (leg D differs from
  leg C only through the droplet axis, never sampling noise); rows add
  N_q05/N_q50/N_q95; the F2 `w_tot == 0` guard carried.

**Verification (7-step order).** New suites `tests/test_droplet_prior_config.py`
(10) + `tests/test_droplet_prior_sampler.py` (8) + the twin oracle test (1):
defaults/back-compat field defaults; typo + all four guard directions;
inverse-CDF uniformity (3σ sample bands, n = 20000); window containment;
pickup-shift = (2/3)δ² within stated tolerance; truncated-mass hand-computed
pins (0.00142 / 0.01528); heavy-truncation warning; initial-state dispatch
(per-atom radii vary, window-bounded); legacy path never touches the analytic
sampler; `stage_legd` `d_delta` rows ≡ `stage_legc` `c` rows **exactly**
(string-equal CSV cells at m = 200) + `d` rows spread the N axis. Narrow →
broad: T8 suites, initial-state/foundations/birth-law/shell neighbors (105),
then **full suite 2464 passed** (2445 baseline + 19 new; 218 s).

The audit finding stands delivered: no checkpoint schema change, no RNG
change on any default path, the legacy pickup MC untouched. CALIBRATION_MAP
row 30 (arm + Sourced δ family, not knobs). Next per the I.11 sequence:
**T9 leg D** behind its own trigger (pre-registration = `stage_legd` at
m = 20000 before any MD; the four `dc` pilots flip exactly
`droplet_size_prior="kornilov_lognormal"` on the certified `cc` baseline,
N = 50, `zero_gamma` retained per T8-D4) → re-pilot → T4/ihe_ked
solvated-branch scoring. The E2 Landau arm stands ready for the N = 500 run.
Nothing here discharges F5.

## T9 leg D TRIGGERED — C + T8 `kornilov_lognormal` droplet prior; twin re-scored at the D4 primary; predictions PRE-REGISTERED before any MD (2026-07-20)

User trigger ("PROCEED TO IMPLEMENTATION of T9 leg D"), following the T8
delivery entry above. Leg D flips exactly **one semantic lever vs the
certified `cc` baseline**: `droplet_size_prior = "kornilov_lognormal"` at
the config-default D4 primary (⟨N⟩ = 2000, δ = 0.625; the guard-paired
`use_single_droplet_size=False` is the same lever) — cumulative on leg C
(`density_tied`, co-moving shed, and the p = 1 onset coupling carry over).
N = 50, `zero_gamma` retained (T8-D4: the Landau arm remains the N = 500
gate, never a leg lever).

**Builds delivered under this trigger (TDD, RED watched first):**
- Generator: `LEG = "d"` in `gen_tier2_md_confirmation.py` — cumulative on
  leg C plus the T8 prior stamped through a new `build_biphasic_cfg`
  passthrough kwarg (`droplet_size_prior`; an analytic arm also stamps the
  paired boolean — one lever). Dirs `…_tier2probe_conf270_dc{1..4}`. Lock
  test asserts the field-diff vs leg C is exactly `{droplet_size_prior,
  use_single_droplet_size}`. Generator suite 20 passed.
- Scorer wiring oracle (session tool, zero repo change): the §4q scoring
  conventions (droplet_retained excluded; suppressed→bin 0; W₁ over 0–21;
  KE co-moving) re-implemented and **verified to reproduce the recorded
  leg-C `cc` numbers exactly** (all four configs, every column) before any
  leg-D data existed.
- Twin: `stage_legd` was delivered with Slice T8 (entry above); run here at
  m = 20000. **The `d_delta` wiring oracle reproduces the leg-C twin
  verbatim at m = 20000** (supp 0.111/0.051/0.343/0.143, n̄ 4.98/5.31/3.77/
  4.88, trapped 0.055/0.060/0.093/0.055 — the TRIGGERED-entry table above),
  so the `d` rows differ from leg C only through the droplet axis.

**Twin re-score at the kornilov prior (m = 20000; outputs
`h2b_leg_d_predictions.csv` / `h2b_leg_d_ke.csv` under
`data/runs/h2b_forward_model/`, committed; the pre-registered numbers are
recorded here).** Realized prior quantiles N_q05/q50/q95 = 596/1638/4636
(the truncated ln-normal about ⟨N⟩ = 2000); dressing softens (n_eject mean
16.71 → 16.53) and the chord-K spread widens hugely (c1: K_q95 → 15.1 —
deep-strip exposure on the large-droplet side; K_q05 → 0.15 on the small):

| config | supp (δ → kor) | n̄_det (δ → kor) | n₁ frac (δ → kor) | n₁ KE (δ → kor) | trapped (δ → kor) |
|---|---|---|---|---|---|
| c1 | 0.111 → **0.127** | 4.98 → **4.42** | 0.212 → **0.236** | 0.998 → **0.992** | 0.055 → **0.051** |
| c2 | 0.051 → **0.058** | 5.31 → **4.77** | 0.217 → **0.246** | 0.964 → **0.958** | 0.060 → **0.056** |
| c3 | 0.343 → **0.383** | 3.77 → **3.35** | 0.084 → **0.087** | 0.253 → **0.245** | 0.093 → **0.081** |
| c4 | 0.143 → **0.163** | 4.88 → **4.36** | 0.172 → **0.191** | 0.973 → **0.967** | 0.055 → **0.051** |

The droplet axis works both ends of the K spectrum at once: the
small-droplet side (short chords, low K) adds suppressed/bare weight and
n₁ occupants (the F.2b low-K direction), the large-droplet side strips
deeper, and the net histogram **softens** (n̄ down 0.42–0.55 He) while the
**KE axis is prior-quiet** (≤ 0.008 eV on every n₁ entry — the descent
speeds are set by the drag law and the onset, not the droplet size).
Suppressed ordering c3 > c4 > c1 > c2 preserved; trapped *drops* slightly
(marginal ions in small droplets clear the shallower well).

**Pre-registered predictions DP-P1..P4 (all KE claims co-moving; same
N = 50 bridge seed, margin 3 Å, exclude policy, 8000 ps cap, `zero_gamma`;
scoring on the solvated branch, bare renormalised out):**

- **DP-P1 (droplet-axis direction).** MD `dc` suppressed fractions rise
  from the `cc` measured 0.101/0.022/0.301/0.124 toward the twin kornilov
  **0.127/0.058/0.383/0.163**, ordering c3 > c4 > c1 > c2. Re-filling
  caveat carries (legs B/C measured MD *below* twin by ≈ 1–7 ions/100 —
  in-bubble pickup closes self-unbound gaps before gate-open).
- **DP-P2 (solvated histogram).** n̄_det **falls** ≈ 0.4–0.55 He vs `cc`
  (to ≈ **4.42/4.77/3.35/4.36**); n₁ rises to ≈
  **0.236/0.246/0.087/0.191**; the deep tail extends (large-droplet
  deep-strip side). W₁(MD, twin) expected at the chain scale (≲ 0.5 bins,
  cf. leg C 0.31–0.47).
- **DP-P3 (KE axis — the twin claims prior-invariance).** Per-bin n₁ mean
  KE ≈ **unchanged vs leg C** (twin 0.992/0.958/0.245/0.967 eV, ≤ 0.008 eV
  off the delta-prior values); MD expected ×1.08–1.15 above twin (the
  leg-C composition residual). A large MD KE move under the prior flip
  would be a **model-structure finding** (the twin says the KE axis cannot
  see the droplet axis).
- **DP-P4 (trapped class).** Twin ≈ 0.051/0.056/0.081/0.051 (down vs
  delta), c3 largest; MD ≈ 1.7–1.9× the twin's 150 ps chord read (the
  A′/B/C precedent).

**Listed twin-divergence channels (S2c-P4 — an unlisted one is a
model-structure finding):** (a) pickup re-filling after under-dressed birth;
(b) per-**atom** ion-t0 dressing vs the twin's molecule-center dressing;
(c) trapped-class dynamics beyond the 150 ps chord read; **(d, new — the
droplet-axis pair):** the MD's well/gate depths follow each ion's R_i
dynamically through the cascade (the twin's chord is frozen at birth
geometry), and small-droplet ejections reach exposure underflow earlier;
(e) the prior draw shifts the neutral-stage RNG stream, so leg D is an
ensemble-level A/B like every T9 leg (no trajectory matching claimed).

MD pilots (`dc1..4`) launch next under this trigger via the leg-C
per-stage-resume accommodation (shrunk relaxation checkpoint,
scored-read-neutral per I67); scored against the kornilov twin, the
certified `cc` baseline, and these pre-registrations on the solvated
branch. Nothing here discharges F5.

## T9 leg D EXECUTED — DP-P3/P4 CONFIRMED (KE prior-quiet; trapped drops), DP-P1/P2 SPLIT: ordering exact + twin band within the listed channels, but the W₁ chain-best breaks (0.54–0.73) on a uniform ≈ 0.5 He softening — the twin's frozen geometry measured (I68–I70) (2026-07-20)

The four `dc{1..4}` pilots ran to completion (all five artifacts each; four
per-config parallel resume drivers under the I67 shrunk-checkpoint
accommodation, ≈ 14 min wall-clock, ~2 MB `relaxation.npz` each) and were
scored against the pre-registered kornilov twin and the certified `cc`
baseline (findings **§4r** + I68–I70). One semantic lever vs `cc`
(`droplet_size_prior="kornilov_lognormal"` + its guard-paired boolean). The
scorer was oracle-locked first: it reproduces the recorded §4q `cc` numbers
exactly on every column.

**Scored A/B/twin (droplet_retained excluded; suppressed→bin 0; W₁ over
0–21; KE co-moving):**

| config | supp (D / twin_d / C) | n̄_det (D / twin / C) | n₁ (D / twin) | W₁(D, twin) | n₁ KE (D / twin / C) | trapped (D / twin) |
|---|---|---|---|---|---|---|
| c1 | 0.087 / 0.127 / 0.101 | 3.89 / 4.42 / 4.87 | 0.228 / 0.236 | **0.71** | 1.037 / 0.992 / 1.077 | 0.080 / 0.051 |
| c2 | 0.077 / 0.058 / 0.022 | 4.19 / 4.77 / 5.27 | 0.154 / 0.246 | **0.73** | 0.982 / 0.958 / 1.104 | 0.090 / 0.056 |
| c3 | 0.393 / 0.383 / 0.301 | 2.82 / 3.35 / 3.76 | 0.067 / 0.087 | **0.54** | 0.317 / 0.245 / 0.317 | 0.110 / 0.081 |
| c4 | 0.098 / 0.163 / 0.124 | 3.93 / 4.36 / 4.87 | 0.196 / 0.191 | **0.68** | 1.038 / 0.967 / 1.089 | 0.080 / 0.051 |

**Verdicts (pre-registered DP-P1..P4; detail in §4r):**
- **DP-P1 SPLIT.** Suppressed ordering exact (c3 > c4 > c1 > c2); c2/c3
  rise as predicted; c1/c4 *fall* vs `cc`, landing 4.0/6.5 ions/100 below
  the twin — inside the pre-listed re-filling band (channel (a)). The
  net-rise-vs-`cc` arm is refuted on c1/c4; the twin-relative claims hold.
- **DP-P2 SPLIT.** n̄ falls on every config but ≈ 2× the predicted drop
  (0.94–1.08 vs 0.42–0.55 He) — a uniform 0.43–0.58 He MD-below-twin
  softening (≈ 1 SE per config, same-signed ×4). **W₁(D, twin) = 0.54–0.73
  — the chain-best trend breaks** (C was 0.31–0.47). Attributed to the
  listed channel (d): the twin's chord is frozen at birth geometry; the
  MD's well/gate/pickup exposure follows R_i through the live cascade.
- **DP-P3 CONFIRMED (the structural claim).** The KE axis is prior-quiet
  in MD as the twin claimed: n₁ KE moves ≤ 0.12 eV under the prior flip
  (c3 unchanged at the printed precision); MD ×1.03–1.07 above twin. The
  solvated KE curve stays owned by the drag law + onset (I68).
- **DP-P4 CONFIRMED (ratio refined).** Trapped drops to 0.080/0.090/0.110/
  0.080 (twin direction exact), c3 largest, ×1.4–1.6 the chord read
  (predicted 1.7–1.9).
- **S2c-P4 holds:** every deviation is a listed channel ((a)/(d)/N = 50
  statistics). No unlisted divergence.

**The chain is complete.** A′→D: every §I.11 physics arm transfers
directionally; the droplet axis is the first where the twin's quantitative
authority ends — the basin-locator stance is now a measurement (I69).

Next per the I.11 sequence: **the confirmation-matrix re-pilot at
MD-located knobs** (re-centering MD-driven; the twin's droplet-axis output
is ordering/direction only) → the T4/ihe_ked **solvated-branch** scoring →
winner at N = 500 with the Landau arm stamped on (`v_L` re-pinning is the
one open domain-expert calibration before that run). Nothing here
discharges F5.
