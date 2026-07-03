# Tier 2 — Phase F Implementation Plan (Calibration Campaign + Production Switch)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Module descriptions are *interface
> contracts*; run matrices are *specifications*.
>
> **Parent:** `TIER2_IMPLEMENTATION_PLAN.md` §4 (Phase F = Slice R2).
> **Entry docs:** `CALIBRATION_MAP.md` (parameter classes + Tier-2 anchor coverage —
> the postponed-items index); MASS doc §6.11 (regime / `Π` / `t×` / staged
> calibration), §R5 (matched-time), §R6 (structural pairing), §10 (calibration
> targets); `TIER2_PHASE_E_IMPLEMENTATION_PLAN.md` (the comparison layer Phase F drives).

---

## 0. Status and intent

Phase F is where the generative mechanism (Phases A–C), the bridge (Phase D), and
the comparison layer (Phase E) are finally **exercised against experiment**, and
where everything the earlier phases deliberately postponed comes due. It is a
**calibration campaign**: sweep the genuinely-free knobs, score each run's terminal
I⁺Heₙ size distribution against `data/reference/integrated_i_he_abundance.csv` with
the Wasserstein metric, and emit an **identifiability + regime-determination
report** — *reported, not auto-adjudicated* (the Tier-1a reporting-gate stance).

Phase F composes **only accepted A–E modules**; it introduces **no new physics**
and re-litigates **no mechanism** — only calibration *values* move, via config
knobs. New surfaces are **orchestration scripts + report assemblers + figures**
(plus one total-strip run-tag variant). **No new `SimConfig` physics fields.**

### Postponed-items inventory (now in scope)

| Item | Source | Class |
|---|---|---|
| κ (ladder steepness) + electronic picture co-fit | CALIBRATION_MAP rows 19/20 | **Free** (the only true knobs) |
| `f_int` (floor scenario-keyed), `f_ret`, `τ` band [2.6,16.5] ps | rows 14/13/11 | **Bounded** |
| Production onset switch 0.80 → **2.70 eV** | row 15 | Sourced, scenario-keyed |
| Identifiability reporting (`f_ret` via 9/18 Å, `τ` via tail, `s`↔κ) | §"Anchor coverage" | **the central deliverable** |
| Regime determination (shell-retaining ↔ total-strip) | §6.11 | reported output |
| Calvo24 total-strip secondary runs | cross-cutting table | secondary/sensitivity |
| All CLI / orchestration / scoreboard / figure scripts | deferred from Phase E | — |
| R6 9-Å re-extraction; p↔κ cap split | R6 / row 7p | conditional triggers |

### Confirmed scope decisions (2026-06-29)

1. **Staged (clean-fallback) calibration** (MASS §6.11): co-fit only κ+picture; pin
   `f_int` at a fixed `τ`; then a `τ`-insensitivity sweep; `f_ret` via the 9/18 Å
   contrast. The full κ×picture×f_int×f_ret×τ grid is a **documented escalation**
   only if staging shows entanglement.
2. **Larger-N, single-seed runs** (N≈500 → ~1000 fragments/run) for a stable
   size-distribution histogram / Wasserstein score.
3. **Build the total-strip secondary-run path now**; keep R6 re-extraction and the
   p↔κ split as **documented conditional triggers** (built only on demand).

---

## 1. Prerequisites

Phase F runs **last** — it composes accepted modules from every earlier phase:

- **A** — `physics/dissociation_ladder.py` (L), `solvation_cooling.py` (K),
  `internal_energy_budget.py` (U).
- **B** — `physics/helium_density.py` (ρ), `pickup.py` (P), `evaporation.py` (Q),
  `mass_jump.capture`.
- **C** — Slice X (checkpoint v7 `E_int` + 5-term invariant); Slice G
  (`ion_propagation_step.biphasic_step` driver).
- **D** — Slice Z bridge: the **0.80 eV validation cross-check** that gates the
  production switch. **Delivered 2026-07-03** — wiring-clean (all sharp oracles
  pass: t×, Π, 5-term closure, kinematics at baseline), but the pinned priored
  point **misses the anchored staircase** (~0.7 sheds vs 7; evaporation-side
  levers κ/picture/τ). That standing flag is the concrete meaning of "the
  0.80 eV validation lands" in the F5 gate: the campaign's Stage 1/2 must
  resolve it (or surface it as a mechanism-level OQ) before the 2.70 eV
  switch. See `TIER2_PHASE_D_BRIDGE_FINDINGS.md` §2.
- **E** — E1 size-dist extractor, **E2 relaxation stage**, E3 abundance loader, E4
  Wasserstein, E5 diagnostics.

No new `SimConfig` physics fields — Phase F uses the knobs added in A–E.

---

## 2. Slices F1–F6

Phase F is a campaign, not a physics build: most slices are scripts mirroring the
delivered Tier-1a scaffolding (`scripts/tier1a_common.py`, `gen_tier1a_runs.py`,
`scripts/post_processing/tier1a_rmse_table.py`). Tests are **harness smoke tests on
tiny N** — no production-sized runs and no figures in pytest.

### F1 — Campaign harness: `biphasic` config builder + run-dir convention

- **Purpose.** One config factory + naming convention for every campaign run.
- **New surface.** `scripts/tier2_common.py` →
  `build_biphasic_cfg(case, variant, *, picture, kappa, f_int, f_ret, tau_ps,
  num_molecules, ion_time_ps, dt_ion_ps, seed, coulomb_available_eV,
  relaxation_time_ps, coeff_overrides=None, e_bind_override=None)` — mirrors
  `tier1a_common.build_anchored_cfg`. Sets `mass_scenario="biphasic"`,
  scenario-stamped `coulomb_available_eV` (0.80/2.70),
  `allow_inconsistent_mass_pairing=True` (R6 structural), `relaxation_stage_enabled
  =True`, `mass_initial_amu=complex_mass_amu(21)`, and the ladder/pickup/evap/E_int/τ
  knobs. Plus `tier2_run_tag(...)` / `tier2_run_dir_name(...)` mirroring the Tier-1a
  tag helpers (encode picture, κ, f_int, f_ret, τ, budget, case; + a `total_strip`
  variant tag).
- **Reuse.** `scripts/tier0_common.py :: build_drag_cfg`, `run_dir_name`;
  `physics/shell_schedule.py :: complex_mass_amu`; the `replace(...)` +
  `cfg.validate()` pattern from `build_anchored_cfg`.
- **Tests.** `build_biphasic_cfg` validates; the §6.5 pairing guard is tripped and
  cleared by `allow_inconsistent_mass_pairing`; tag/dir round-trips; budget stamp
  matches the scenario.

> **NB (as-built, 2026-07-03 — Slice Z seeded this module; F1 is an *extension*,
> not a creation).** `scripts/tier2_common.py` was delivered at Phase D with the
> bridge-scoped subset of this surface:
> `build_biphasic_cfg(case, variant, *, num_molecules, ion_time_ps, dt_ion_ps,
> seed, lambda0_per_ps=0.9, f_int=0.5, f_ret=0.1, coeff_overrides=None,
> e_bind_override=None)` — pinned-point defaults; `coulomb_available_eV` stamped
> 0.80 inside the builder; τ/κ/picture ride the config defaults (not kwargs) —
> plus the fixed `TIER2_BRIDGE_TAG` / `tier2_bridge_run_dir_name` (single run, no
> knob encoding). F1 therefore **extends** the delivered builder with the
> campaign kwargs this contract lists (`picture`, `kappa`, `tau_ps`,
> `coulomb_available_eV`, `relaxation_time_ps` / the E2 enable) and adds the
> knob-encoding `tier2_run_tag` / `tier2_run_dir_name`; the bridge tag stays
> reserved as-is. Two corrections to the planned signature: (1) **λ₀ joins it**
> (`lambda0_per_ps`, the delivered kwarg — the planned signature omitted the
> Bounded row-7 knob entirely, and the `f_ret` 9/18 Å contrast runs need pickup
> live); (2) **back-compat constraint:** the delivered Phase-D scripts
> (`gen_tier2_bridge_run.py`, `tier2_bridge_report.py`) import this module —
> extend by kwargs-with-defaults so the bridge run stays byte-reproducible; do
> not repurpose the delivered defaults. F1's planned tests remain open (Phase D
> landed no `tier2_common` unit tests; the module is exercised through the
> bridge run and the `test_derived_diagnostics.py` driver smoke only).

### F2 — Run-matrix generator (staged campaign, 0.80 eV, 9 Å + 18 Å, larger N)

- **Purpose.** Produce the self-describing campaign run dirs for the **staged**
  validation campaign.
- **New surface.** `scripts/gen_tier2_runs.py` mirroring `gen_tier1a_runs.py` (USER
  SETTINGS block). Pipeline per grid point: neutral → ion (`biphasic_step` driver,
  C) → relaxation stage (E2). Stages:
  - **Stage 1 (κ×picture co-fit):** κ grid (≈5 points, gradual→cliff) × picture ∈
    {`statistical_mixture`, `x2_only`, `cooling_relaxed`}, at fixed `f_int=`floor,
    `τ=τ_mid`, `f_ret=`prior.
  - **Stage 2 (τ-insensitivity):** at the Stage-1 best κ/picture, sweep
    `τ∈[2.6,16.5]` ps; confirm terminal-n insensitivity (else flag).
  - **`f_ret` contrast:** every grid point runs **both 9 Å and 18 Å** cases so the
    9/18 Å density contrast can discriminate `f_ret`.
  - Budget **0.80 eV first** (validation; cross-checked by the Phase-D bridge).
  - **N≈500, single seed** per grid point.
- **Reuse.** `gen_tier1a_runs.py` `_run_one` pipeline; `RunDirectory`,
  `run_neutral_propagation`, `run_ion_propagation`; `tier2_common` builders;
  **`scripts/gen_tier2_bridge_run.py`** (delivered at Phase D) as the direct
  single-run precedent — its `_run_one` + USER-SETTINGS shape is the template.
- **Tests.** Smoke: tiny-N single-grid-point generation writes a valid run dir (no
  production scale, no figures).

> **NB (bridge priors for the stage design, 2026-07-03 —
> `TIER2_PHASE_D_BRIDGE_FINDINGS.md` §2; annotations, not a redesign).**
> The Phase-D run at the pinned central point updates four premises of the
> staged campaign:
> 1. **Stage 1 (κ×picture) now has a quantitative target.** The cascade
>    under-sheds ~10× at κ = 1 / mixture (Δn̄ = 0.67 vs 7; RRK integral ≈ 0.7
>    expected sheds — kinetically, not energetically, limited). The κ grid must
>    deliberately span the **sharp-cliff end**: κ sets `D₀(n)/Σ(n)` near n*,
>    which enters the RRK factor at the (s−1) = 59 power, so k rises steeply
>    with κ. If no κ/picture/τ combination inside the bands lands the
>    staircase, that is the mechanism-level OQ (RRK dof convention) the
>    findings doc names — surface it, do not silently retune ν or s.
> 2. **Stage 2's premise may invert.** The plan expects a τ-*insensitivity*
>    confirmation; the bridge mean-rate estimate says terminal n **is**
>    τ-sensitive at 0.80 eV (τ → 16.5 ps roughly triples the shed integral).
>    Treat the flag arm as the likely outcome and τ as a live calibration
>    dimension at this budget, not a formality.
> 3. **The Stage-1 `f_int = floor` pin is a timing choice, not a neutral
>    default.** The cascade *budget* is pinned at Σ(21) by the crossing
>    definition — f_int moves only t× (at the floor, t× ≈ 0: the gate is open
>    from onset; at the bridge's 0.5, t× ≈ 5 ps, on the GAH25 prior). Fixing
>    f_int at the floor forfeits the staircase-*timing* comparison; make the
>    choice consciously per stage.
> 4. **The 9 Å / 0.80 eV runs carry ~no λ₀ sensitivity** (pickup structurally
>    dead through the decline: Π ≤ 0.005) — λ₀/`f_ret` discrimination rests on
>    the 18 Å contrast leg, which is therefore load-bearing, not optional.
>
> Operational note (findings §3): `RunDirectory.load_cfg` refuses unknown
> `cfg.json` fields, so campaign run dirs are invalidated by any later config-
> surface change (the delivered Tier-1a artifacts already are). Score campaign
> runs promptly with the code version that generated them; regeneration is the
> recovery path, not artifact migration.

### F3 — Scoreboard + report assembler (the deferred Phase-E assembler)

- **Purpose.** Compose E1+E2+E4+E5 into **one scoreboard row per run**.
- **New surface.** `scripts/post_processing/tier2_size_distribution_table.py`
  mirroring `tier1a_rmse_table.py` (`collect_tier2_records` / `score_tier2_run` /
  `format_table` / `write_rows_csv`). Per run: E1 terminal size distribution
  (matched-time via E2 **and** sim-end upper bound), E4 Wasserstein vs the abundance
  reference (both times), E5 (`t×`, `Π`, regime, `total_strip_reachable`), and the
  5-term `ion_ledger_closure` residual.
- **Columns.** `case, budget_eV, picture, kappa, f_int, f_ret, tau_ps, N,
  W1_matched, W1_simend_upper, t_cross_ps, regime_label, total_strip_reachable,
  ledger_max_resid_eV, n_terminal_mean, n_terminal_spread`.
- **Reuse.** E3 `load_he_abundance_reference`; E1/E4/E5 APIs; `RunDirectory`;
  `energy_balance.ion_ledger_closure` (5-term, X). The Tier-1a CSV/table machinery
  verbatim.
- **Stance.** Reported, **not** auto-adjudicated — no code asserts a fidelity verdict.
- **Tests.** Score one tiny synthetic run → a row with expected keys; CSV round-trips.

### F4 — Identifiability + regime-determination report

- **Purpose.** The campaign's central deliverable: turn the scoreboard into an
  **identifiability** statement, not a point fit.
- **New surface.** `scripts/post_processing/tier2_identifiability_report.py`
  (headless text/CSV). Emits: the κ/picture co-fit landscape (W₁ vs κ per picture);
  `f_ret` discrimination via the **9/18 Å W₁ contrast**; `τ` via the **post-crossing
  tail** (terminal-n insensitivity — confirm or flag); the `s↔κ` coupling note; and
  the **regime determination** (shell-retaining ↔ total-strip, §6.11, with `t×`/`Π`
  as the spine). Headline reported: does **one** observable separate the 8+
  quantities (CALIBRATION_MAP "Anchor coverage")?
- **Reuse.** The F3 scoreboard rows; E5 diagnostics; CALIBRATION_MAP identifiability
  notes as the report's checklist.
- **Tests.** Report assembles from a small synthetic scoreboard; the 9/18 Å contrast
  and τ-insensitivity flags compute correctly on constructed inputs.

### F5 — Production switch (2.70 eV) + total-strip secondary runs

- **Purpose.** Promote to the production budget and exercise the regime-axis far end.
- **Production switch.** Re-run the **staged** campaign at `coulomb_available_eV
  =2.70` (scenario-stamped → §6.5 guard → `allow_inconsistent_mass_pairing=True`),
  **only after** the 0.80 eV validation lands (Phase-D bridge cross-check). Same F2
  generator, budget flag flipped.
- **Total-strip secondary runs (built now).** A `total_strip` config variant pushing
  the regime to total vaporization (e.g. `f_int` high / `τ` long so the gate never
  self-binds in-window → terminal n→small), evaluated against the size distribution
  as a **sensitivity target, not default** (§6.11 / OQ6 reachability via E5).
- **Documented conditional triggers (NOT built).** Record firing criteria only:
  **R6** 9-Å option-3 drag re-extraction — fire only if the 2.70 eV scoreboard shows
  mass↔coefficient sensitivity; **p↔κ** occupancy-cap split — fire only if the
  first-shell cutoff in the size distribution demands it.
- **Reuse.** F1/F2/F3/F4 unchanged; only the budget stamp + the `total_strip` tag.
- **Tests.** Budget stamp guards the pairing as expected; the `total_strip` variant
  builds and validates; conditional triggers are docstring-recorded, not active code.

### F6 — Overlay figures (deferred from Phase E)

- **Purpose.** The visual comparison artifacts (non-test).
- **New surface.** Figure builders in `tier2_size_distribution_table.py` /
  `tier2_identifiability_report.py` behind a `SHOW_FIGURE` gate, mirroring
  `tier1a_rmse_table.py`'s matplotlib pattern: (a) size-distribution overlay
  (simulated integer-n stem vs experimental abundance, matched-time + upper-bound);
  (b) the regime-axis `Π`/`t×` diagnostic; (c) the W₁-vs-κ co-fit landscape.
- **Constraint.** **No figures in pytest** — figure code is import-guarded and only
  runs under `__main__` / `SHOW_FIGURE`.

---

## 3. Dependency graph and build order

```
A (L,K,U) ─ B (ρ,P,Q) ─ C (X,G) ─ D (Z) ─ E (E1–E5)      [all accepted]
                                              │
F1 (harness) ─► F2 (gen runs, 0.80 eV, 9Å+18Å) ─► F3 (scoreboard) ─► F4 (identifiability)
                                                                          │
                                          F5 (2.70 eV switch + total-strip) ◄─ gated on 0.80 eV landing
F6 (figures) ── alongside F3/F4
```

**Build order:** F1 → F2 → F3 → F4; **F5 after F4** (validation-first: 0.80 eV
lands, cross-checked by the Phase-D bridge, before the 2.70 eV switch); F6 alongside
F3/F4.

---

## 4. Cross-cutting gates

- **Validation-first:** the 0.80 eV campaign (Phase-D-bridge-cross-checked)
  **before** the 2.70 eV production switch.
- **Reported, not auto-adjudicated:** no code asserts a fidelity pass/fail (mirrors
  the Tier-1a reporting-gate stance); adjudication is left to the user.
- **Scenario-stamped `coulomb_available_eV`** guards the budget↔drag/shell pairing;
  `allow_inconsistent_mass_pairing=True` is the *normal* production path (R6
  structural, §6.6 mid-window defence).
- **5-term invariant residual** carried per run as a wiring diagnostic.
- **Composes only accepted modules; no new physics, no mechanism re-litigation** —
  only calibration *values* move via config knobs.
- **Out of scope (fails review if it leaks in):** the noise model / Tier 3 (inert);
  any change to the Tier-0 drag law, neutral propagation, RNG draw order, or the
  locked mechanism; the R6 re-extraction and p↔κ split stay **document-only**
  triggers.

---

## 5. Verification

- **Harness smoke (pytest, tiny N, no figures):** `build_biphasic_cfg` validates +
  trips/clears the pairing guard; `gen_tier2_runs` writes one valid tiny-N run dir;
  `score_tier2_run` returns a row on a tiny synthetic run; the identifiability
  report assembles from a small synthetic scoreboard.
- **End-to-end (manual, non-test):** generate a small κ×picture grid at 0.80 eV
  (9 Å + 18 Å) → scoreboard CSV → identifiability report → overlays; confirm W₁ is
  finite, the 5-term invariant closes, and the regime label + total-strip
  reachability are emitted.
- **Suite:** baseline stays green; any new red is a real regression. Interpreter:
  `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q`.

---

## 6. Risks / notes

- **Calibration load / identifiability is the central dependency (HIGH).** One
  observable (the size distribution) carries 8+ quantities (CALIBRATION_MAP "Anchor
  coverage"). Whether it *separates* them is the key open empirical question — F4
  **reports** identifiability, it does not assume it. The staged strategy is chosen
  precisely to keep this tractable.
- **R5 matched-time (carried).** Terminal n is read post-E2-relaxation **and** at the
  sim-end upper bound; F3 reports both; never over-read absolute n.
- **R6 structural pairing.** Every production run trips the §6.5 guard and runs under
  `allow_inconsistent_mass_pairing=True`; the §6.6 mid-window argument is the
  load-bearing defence. The clean fix (9-Å option-3 re-extraction) stays a documented
  conditional trigger.
- **Compute scale.** Larger-N (≈500) × staged grid × {9 Å, 18 Å} × {0.80, 2.70 eV}
  is the heaviest run set in the port; F2 is a script (not a test) and the staged
  structure bounds the run count.
- **Conditional triggers documented, not built:** R6 re-extraction, p↔κ split — each
  with its firing criterion recorded in F5.

**Cross-links:** `TIER2_IMPLEMENTATION_PLAN.md` §4 (Phase F headline);
`CALIBRATION_MAP.md` (parameter classes + anchor coverage);
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §6.11 / §R5 / §R6 / §10;
`TIER2_PHASE_E_IMPLEMENTATION_PLAN.md` (the comparison layer Phase F drives);
`drag_migration_log_tier2.md` (record the Phase F decision + delivery here).
