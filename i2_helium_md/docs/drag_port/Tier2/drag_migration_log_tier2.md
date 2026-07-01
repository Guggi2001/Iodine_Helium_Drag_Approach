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
