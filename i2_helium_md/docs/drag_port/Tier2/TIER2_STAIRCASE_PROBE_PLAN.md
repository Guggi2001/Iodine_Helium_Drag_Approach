# Tier 2 — Staircase Capability Probe (pre-F5 design)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Module descriptions are *interface
> contracts*; the run matrix is a *specification*.
>
> **Status:** DELIVERED 2026-07-05 under the `[PROCEED TO IMPLEMENTATION]`
> trigger (delivery record: `drag_migration_log_tier2.md`, "Pre-F5 —
> Staircase capability probe"). The `f_int = 0.5` pin (§2) stood when the
> trigger was given; flipping to the Stage-1 floor is a one-line
> USER-SETTINGS change. The 45-run probe **executed 2026-07-05** — outcome
> (b) of §6: no in-band (κ, picture, τ) point lands (freeze at n ≈ 20,
> kinetic; log entry "Staircase probe EXECUTED"), firing the RRK-dof
> mechanism-level OQ and leading to the Addendum-A `s_eff` mini-probe below.
>
> **Parent:** `TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (inserted between F4 and
> F5; does not modify any delivered F-slice). **Entry docs:**
> `TIER2_PHASE_D_BRIDGE_FINDINGS.md` (the standing staircase miss + levers);
> `drag_migration_log_tier2.md` (F1–F4 as-built state).

---

## 0. Purpose and stance

Before the F5 production switch (2.70 eV) and before the N=500 Stage-1/2
campaign spend, run a **small-N, 45-point knob sweep** answering two questions:

1. **Staircase reproduction.** Does *any* in-band (κ, picture, τ) combination
   reproduce the anchored 9 Å shell evolution **21→19→14** in magnitude *and*
   timing? The Phase-D bridge showed the single pinned priored point sheds
   ~0.7 He vs 7 (kinetically limited, s−1 = 59 RRK exponent) — this probe
   extends that existence test across the full in-band lever space.
2. **Flexibility map.** What range of **post-relaxation terminal shells** can
   the mechanism express across the grid — is it expressive enough to produce
   various solvation shells, not just the frozen n ≈ 20 corner?

**Stance (mirrors the bridge):** an **existence probe, reported not
auto-adjudicated**. TDDFT is not ground truth (the rejected "1b" stance
holds — experiment arbitrates at Tier 2 via the size distribution):

- A clean miss everywhere in-band is a *finding* — the RRK-dof
  mechanism-level OQ named in the bridge findings §2 — never a silent ν/s
  retune.
- A landing region does **not** by itself discharge the F5 gate ("the 0.80 eV
  validation lands"); it tells the N=500 Stage-1/2 campaign where to
  concentrate, and that campaign is what formally resolves the gate.
- The probe composes only accepted A–E modules: **no new physics, no new
  `SimConfig` fields, no mechanism re-litigation** — only calibration values
  move, via existing config knobs.

## 1. Prerequisites and reuse

All delivered; the probe adds orchestration + report scripts only:

- **F1** `scripts/tier2_common.py :: build_biphasic_cfg` — every probe knob
  (`picture`, `kappa`, `tau_ps`, `f_int`, `f_ret`, `lambda0_per_ps`,
  `relaxation_time_ps`, `coulomb_available_eV`) is already a kwarg.
- **F2** `scripts/gen_tier2_runs.py` — the `_run_one` pipeline
  (neutral → ion `biphasic_step` → E2 relaxation), overwrite guard,
  skip-completed resume: the direct template.
- **Phase D** `scripts/post_processing/tier2_bridge_report.py` +
  `postprocess/derived_diagnostics.py` — the mean-n(t) machinery and the
  `build_shell_schedule` anchored comparator (zero artifact dependence — the
  stale-artifact lesson from the bridge findings §3).
- **E1** `compute_terminal_shell_distribution` (matched-time on
  `relaxation.npz` with the `source_tag="relaxed"` reload override — the F3
  pitfall), **E5** diagnostics, and the 5-term `ion_ledger_closure`.
- **F3** table idiom (`format_table` / `write_rows_csv`, USER-SETTINGS,
  `main()`).

## 2. Run matrix

- **Grid (45 runs):** κ ∈ {0.5, 1, 2, 4, 8} × picture ∈ {`statistical_mixture`,
  `x2_only`, `cooling_relaxed`} × τ ∈ {2.6, 6.55, 16.5} ps. Full-band τ from
  the start (not staged): the bridge predicts τ-*sensitivity* at 0.80 eV, and
  a capability probe must show the mechanism's full in-band reach — declaring
  "cannot land" while τ sat pinned at mid-band would be unsound.
- **Pinned:** 9 Å case only; budget **0.80 eV** (the staircase anchor is the
  9 Å validation condition); `f_ret = 0.1`; `λ₀ = 0.9/ps`; **N = 50** single
  fixed seed (bridge-validated mean-n(t) stability); 30 ps ion window; Tier-0
  dt/drag bundle; E2 relaxation enabled (the "after relaxation" leg).
- **`f_int = 0.5` (bridge value) — ⚠ flagged for user confirmation.**
  Rationale: `f_int` is timing-only (it moves t×, not the cascade budget,
  which is Σ(21) by the crossing definition). At 0.5, t× ≈ 5 ps — aligned
  with the anchored first shed; at the F2 Stage-1 floor pin, t× ≈ 0 and the
  gate is open from onset, corrupting the staircase-*timing* comparison at
  every grid point. The alternative (floor, matching Stage-1 exactly) was
  offered and remains selectable before implementation.
- **Built-in wiring oracle:** the (κ=1, `statistical_mixture`, τ=6.55) point
  *is* the Phase-D bridge configuration. Its ion stage must reproduce the
  bridge numbers (Δn̄ ≈ 0.67, t× reconstructed ≈ 4.96 ps, all-ions-agree) —
  a free end-to-end regression check on the whole probe pipeline.

## 3. Generator — `scripts/gen_tier2_staircase_probe.py` (new)

Mirrors `gen_tier2_runs.py`: USER-SETTINGS block, grid enumerator returning
`(label, cfg, run_dir)` with no propagation, `_run_one` pipeline
neutral → ion → E2 relaxation (`relaxation.npz`), `OVERWRITE_EXISTING_RUN`
guard, `SKIP_COMPLETED_RUNS` resume, generator-level budget guard.

- Configs built exclusively through `build_biphasic_cfg` — no new
  `SimConfig` fields, no physics.
- **Distinct namespace:** probe run dirs use a probe tag (working name:
  `_tier2probe_` in place of `_tier2_`), emitted by an additive
  `tier2_probe_run_tag` / `tier2_probe_run_dir_name` helper pair in
  `tier2_common` (kwargs-with-defaults; the bridge tag and campaign tag
  helpers are untouched). The F3 campaign glob `*_tier2_*` does not match
  `_tier2probe_` (no literal `_tier2_` substring) — verified at
  implementation time with a lock test; if the boundary proves leaky, F3's
  explicit-exclusion precedent (the bridge-tag exclusion) is the fallback.

## 4. Report — `scripts/post_processing/tier2_staircase_probe_report.py` (new)

A *pure scorer* (runs nothing) over finished probe dirs
(`cfg.json` + `ion.npz` + `relaxation.npz`), one row per run:

- **Staircase metrics** vs the anchored comparator (`build_shell_schedule`,
  Tier-1a family; zero artifact dependence):
  - `dn_mean_window` — mean total sheds in-window (target ≈ 7),
  - `n_ion_end_mean` — ion-end mean n (target 14; R5: this is the sim-end
    read, distinct from the relaxed read below),
  - `t_first_shed_ps` — mean first-shed time vs the anchored t★,
  - `n_traj_mad` — time-averaged |n̄(t) − n_anchor(t)| over the window.
- **Flexibility metrics** (E1 on `relaxation.npz`, matched-time,
  `source_tag="relaxed"`): `n_relaxed_mean` / `n_relaxed_spread` per point;
  report headline = the reachable relaxed terminal-n range across the grid,
  per picture.
- **Wiring column:** 5-term `ion_ledger_closure` max residual per run.
- **Knob columns** read from the authoritative `cfg.json`, never parsed from
  the tag (the F3 convention).
- **Output:** text table + CSV (F3 idiom). **Figures behind `SHOW_FIGURE`
  only, never in pytest:** (a) small-multiple mean-n(t) overlays — one panel
  per picture×τ cell, one curve per κ, anchored staircase underlaid;
  (b) a relaxed-terminal-n heat map over the grid.
- **Reported, not auto-adjudicated** — no code asserts a capability verdict.

## 5. Tests (pytest, tiny N, no figures)

Per the F-slice precedent:

- grid enumerator yields 45 unique dirs/tags, pinned knobs stamped correctly;
- probe tag never collides with the campaign or bridge namespaces, and F3
  `discover_run_dirs` ignores a probe dir (the namespace lock);
- tiny-N end-to-end `_run_one` writes all artifacts;
- staircase metrics hand-oracled on a constructed synthetic n(t) (known step
  function → known `dn_mean_window`, `t_first_shed_ps`, `n_traj_mad`);
- report assembles from a scored tiny run; CSV round-trips;
- the bridge-point wiring oracle is a *manual* (non-pytest) check — it needs
  the N=50 run, which is production-scale by test standards.

## 6. Outcomes and F5 relation

The user adjudicates (reporting-gate stance). The three outcome shapes:

- **(a) A landing region exists** → point the N=500 Stage-1/2 campaign there;
  proceed toward F5 with the region as the prior.
- **(b) Nothing in-band lands** → surface the RRK-dof mechanism-level OQ
  (bridge findings §2 lever 3) *before* any production spend; the ν/s
  conventions become the suspects — an OQ-class item, not a retune.
- **(c) Partial landing** (magnitude without timing, or vice versa) → the
  per-metric split (`dn_mean_window` vs `t_first_shed_ps`) localizes which
  lever family is short; likely feeds a targeted τ/f_int timing discussion.

## 7. Risks / notes

- **Small-N read.** N=50 stabilizes the *mean* n(t) (bridge-validated), not
  the terminal-n *distribution* tails — the flexibility map reads means and
  spreads, never distribution shape (that is the N=500 campaign's job).
- **R5 carried.** Two terminal reads per run (ion-end vs relaxed) are both
  reported and never conflated; the staircase target (14) is an *in-window*
  ion-end quantity, the flexibility map is the relaxed read.
- **Namespace hygiene is load-bearing** (the F1/F2 review lesson): the
  namespace lock test is first-class, not optional.
- **Stale-artifact policy** (bridge findings §3): probe dirs are invalidated
  by any later config-surface change — score promptly with the generating
  code version; regeneration is the recovery path.
- **Compute:** 45 × N=50 ≈ 3/10 of one 15-point N=500 stage — bounded, and
  it can prevent the full campaign being spent on an infeasible band.

**Cross-links:** `TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (F5 gate);
`TIER2_PHASE_D_BRIDGE_FINDINGS.md` (miss + levers);
`CALIBRATION_MAP.md` (parameter classes);
`drag_migration_log_tier2.md` (record the probe decision + delivery there).

---

# Addendum (2026-07-05) — `s_eff` mini-probe (RRK-dof falsification)

> **Status: DELIVERED 2026-07-06** (harness + tests, under the
> `[PROCEED TO IMPLEMENTATION]` trigger; delivery record:
> `drag_migration_log_tier2.md`, "s_eff mini-probe harness DELIVERED").
> **EXECUTED 2026-07-06 — outcome (a) of §A.4:** a constant s_eff lands the
> staircase in magnitude *and* timing (s_eff = 8 at τ = 6.55 ps: 7.41 sheds →
> n_end 13.59, t_first 5.42 ps vs t★ = 5, trajectory MAD 1.0 He; landing
> region s_eff ∈ ≈ [8, 12] × mid-band τ; the τ = 16.5 arm lands magnitude
> only — timing separates, as designed). The picture cross-check (`x2_only`
> + `cooling_relaxed` under reduced s_eff) executed the same day: magnitude
> near-invariance **confirmed** (≤4 % spread). **Promotion CONFIRMED (user,
> 2026-07-06):** s_eff reclassified Derived → Bounded (CALIBRATION_MAP
> row 10 + update block; MASS doc A11 resolution NB; Phase F plan F2
> Stage-1 re-scope NB to an s_eff×τ co-fit). The selection-surface /
> campaign-generator implementation stays behind
> `[PROCEED TO IMPLEMENTATION]`. Log entries: "s_eff mini-probe EXECUTED",
> "Picture cross-check EXECUTED", "s_eff PROMOTED".

## A.1 Why (probe outcome b)

The 45-point probe **executed** (2026-07-05; results in the migration log):
no in-band (κ, picture, τ) lands the staircase — max Δn̄ = 1.7 sheds vs 7,
relaxed terminal n confined to [19.3, 20.95]. Analysis against the delivered
ladder/evaporation code localises the freeze **kinetically**:

- `k = ν·(1 − x)^(s−1)` with `x = D₀(21)/Σ(21) ≈ 0.032` (κ=1, mixture) and
  `s − 1 = 59` gives `k ≈ 0.15·ν ≈ 0.36/ps` — the suppression is entirely the
  exponent acting on a small `x`.
- The swept levers only move `x` (bounded ≲2×: the Form-U cliff centre is
  pinned at n*+½ = 21.5, so σ(21) ≤ ½ floors D₀(21) at ~0.53·D₀(1); the κ
  direction is *inverted* vs the bridge-findings lever-1 claim) or the window
  (τ, capped by t× = τ·ln(f_int·E/Σ(21)) outrunning the 30 ps window). The
  picture knob is near-degenerate (`x` is nearly picture-invariant since
  D₀ and Σ both scale with D₀(1) − D_floor).
- Energetics are self-sustaining by construction (per shed, E_int and Σ(n)
  drop by the same D₀(n): gate margin shed-invariant), so `s` is the only
  exponential lever. This fires the **pre-registered RRK-dof mechanism-level
  OQ** (bridge findings §2 lever 3): `s = 3n−3` (60 modes at n=21) treats a
  floppy quantum He shell as fully-coupled classical oscillators.
- **Tier-2-blocking, not staircase-only:** the crossing construction caps the
  cascade budget at Σ(21) at *any* `coulomb_available_eV` (budget moves only
  t×), and the experimental abundance holds 43 % bare I⁺ with <1 % weight at
  n ≥ 19 — a mechanism frozen at n ≈ 20 cannot approach the arbitration
  observable at 0.80 **or** 2.70 eV.

## A.2 The experiment

Sweep the **already-plumbed** diagnostic override `cfg.evap_rrk_dof`
(Optional, guarded s ≥ 1 at config-load, physics-live at
`ion_propagation_step.py:630`) at the bridge point. **No new physics, no new
`SimConfig` fields** — one existing config value moves.

- **Grid (12 runs):** `s_eff ∈ {None, 30, 20, 12, 8, 5}` ×
  `τ ∈ {6.55, 16.5} ps` (`None` = the per-n `3n−3` status quo — the in-grid
  control; the τ arm tests s↔τ separability: first-shed *timing* stays
  gate-open-governed ≈ t×, magnitude is s-governed).
- **Pinned:** κ = 1, `statistical_mixture` (both shown flat/capped), f_int =
  0.5, f_ret = 0.1, λ₀ = 0.9/ps, 9 Å / 0.80 eV / N = 50 / bridge seed —
  i.e. the probe pins with only `s_eff` and τ moving. Cost ≈ ¼ of the
  45-run probe.
- **Note on constancy:** the override is a *constant* s (the per-n convention
  varies 60→39 over 21→14). Constant-s is the correct *diagnostic*; the final
  convention shape (constant vs scaled `α·(3n−3)`) belongs to the OQ
  resolution, not this probe.

## A.3 Harness (additive to the delivered probe scripts)

- `build_biphasic_cfg` gains a None-sentinel `evap_rrk_dof` kwarg (the
  established pass-through pattern; `None` = ride the config default —
  back-compat byte-identical).
- `gen_tier2_staircase_probe.py` gains `S_EFF_GRID`; a set value appends
  `_sNN.NN` to the probe tag (collision-safe, stays in the `tier2probe`
  namespace); the per-n `None` appends nothing (tag byte-identity).
- `tier2_staircase_probe_report.py` gains an `s_eff` column read from
  `cfg.json` (`None` → per-n).

> **As-built (2026-07-06).** The generator's **active USER SETTINGS are the
> mini-probe** (κ=1, mixture pinned; τ ∈ {6.55, 16.5} × `S_EFF_GRID =
> [None, 30, 20, 12, 8, 5]` → 12 runs), since the 45-run probe already ran.
> The full 45-run grid is preserved in a "FULL PROBE" comment block for
> one-swap restoration, and the "delivered 45-run behaviour/tags unchanged"
> promise is kept where it is real — **tag byte-identity** — proven by
> `test_full_grid_back_compat_none_s_eff` (full grid + `S_EFF_GRID=[None]` →
> 45 points, `s_eff`-unset tag equals the pre-addendum string). `probe_grid_
> points()` now returns 4-tuples `(picture, κ, τ, s_eff)`. The grid tuple grew
> a dimension but the emitted run dirs/cfgs are byte-identical when `s_eff` is
> `None`. No new `SimConfig` field (`cfg.evap_rrk_dof` was already declared and
> physics-live at `ion_propagation_step.py:630`).

## A.4 Read-out and interpretation boundaries

Same metrics, same **reported-not-adjudicated** stance. Outcome shapes:

- **(a)** some `s_eff` lands magnitude *and* timing → the promotion decision
  goes to the user: demote `s = 3n−3` to one arm of a dof-convention enum and
  promote `s_eff` to a **Bounded** calibration knob (literature: evaporative
  ensembles of quantum clusters), re-scoping the F2 campaign to an s_eff
  co-fit (CALIBRATION_MAP reclassification at that point, not now).
- **(b)** no constant `s_eff` lands both → the gate/ladder structure itself
  becomes suspect (the `TabulatedLadder` declared fallback, or the
  Σ(21)-crossing budget construction) — a deeper mechanism discussion.
- **(c)** landing only with overshoot (moderate s strips past 14) → an
  n-dependent `s` (scaled convention) is indicated.

A landing `s_eff` does **not** auto-promote the convention; this probe
falsifies/localises, the user adjudicates.

---

# Addendum B (2026-07-07) — post-Wave-4 decisions + gated landing re-location mini-probe (Wave 5)

> **Status: EXECUTED 2026-07-07** (under the `[PROCEED TO IMPLEMENTATION]`
> trigger; scratchpad-driver route per B.4, zero repo-code change). Outcome
> shape **(a)+(b) mixed**: a gated s_eff lands the in-window staircase
> (s_eff ≈ 30, ~4× the ungated landing) but timing is ~2.3 ps late (f_int
> re-alignment territory) **and** a new structural finding — the gated
> in-window staircase and the *terminal* read decouple (s_eff sets the
> cascade rate, not the endpoint; terminal is flight-time dependent above the
> freeze-out s_eff). Full results + tables: `TIER2_STAIRCASE_PROBE_FINDINGS.md`
> §4b (insights I10–I12); delivery/decision record:
> `drag_migration_log_tier2.md`, "Wave-5 gated landing re-location EXECUTED".
>
> Follows the Wave-4 cooling-gate total-strip A/B results and their discussion
> (`TIER2_STAIRCASE_PROBE_FINDINGS.md` §4; log entry "Total-strip A/B probe
> EXECUTED", 2026-07-07).

## B.1 Decisions (user, 2026-07-07)

1. **`cooling_spatial_gate = density_scaled` is the intended primary
   production arm.** Grounds: physical self-consistency (drag, pickup, and
   cooling share one ρ_He bubble boundary; ungated K2 applied a
   droplet-bath cooling law in vacuum after ejection — the GAH25 τ was fit
   for a near-droplet growing shell), and capability (the ungated arm
   cannot pass n̄ ≈ 8.5 at mid-band τ even at the max-kinetics bound
   s_eff = 1, structurally blocking the experimental 43 % bare-I⁺ peak;
   the gated arm reaches n̄ ≈ 2.7). This is a *direction*, not a
   hard-wire: the arm stays an interchangeable `SimConfig` enum with a
   `none` sensitivity leg in any campaign (working method), and formal
   promotion into the campaign/production wiring waits for the Wave-5
   result below.
2. **`f_int` is promoted to a genuine sweeping parameter — accepted
   consequence of the gate.** Under `density_scaled` the t×↔ejection race
   is binary (Wave 4: deep strip if t× ≲ t_eject ≈ 6 ps, permanently
   closed gate otherwise), and t× = τ·ln(f_int·E/Σ(21)) makes f_int (with
   τ and the picture-degenerate Σ(21)) the cliff-side selector rather than
   a timing nicety. **Derived prediction recorded for F5:** at 2.70 eV
   with f_int = 0.5, t× ≈ 12.9 ps > t_eject → the gated arm sheds
   *nothing*; opening the gate before ejection at the production budget
   needs roughly f_int ≲ 0.17 — the scenario-keyed f_int floor
   (CALIBRATION_MAP row 14) is load-bearing at production, and ejection
   only gets earlier at higher budget (faster ion).
3. **The clean erf-complement (exact zero out-of-bubble cooling) is
   adopted as correct.** A small residual floor (`rho_min`, dragged He
   cloud) is judged *not physical* by the user; it is retained purely as a
   **documented open question for domain experts** — explicitly **no
   implementation**, now or as an escalation default. (Consequence
   accepted with it: the Wave-4 all-or-nothing τ/f_int cliff is a real
   feature of the model, not an artifact to be softened.)

Picture and κ stay as adjudicated at the s_eff promotion: κ pinned at 1
(near-dead: inverted + normalisation-capped), picture pinned at
`statistical_mixture` and reported — picture enters the dynamics only
through the ratio f_int·E/Σ(21), i.e. it is one effective knob with f_int
(Waves 1+3), now swept via f_int alone.

## B.2 Why Wave 5 — the s_eff landing prior is gate-conditional

The staircase landing **s_eff ∈ ≈ [8, 12] × mid-band τ** (Addendum A) was
established **ungated**. Wave 4 probed the gated arm only at s_eff ≤ 5,
where the gate roughly doubles in-window shedding at τ = 6.55 (s_eff = 5:
Δn̄ 9.35 → 17.52) by removing the post-ejection K2 quench. Extrapolating,
gated s_eff = 8 likely overshoots the anchored 7-shed/n_end-14 staircase
badly — the gated co-fit optimum should sit at **higher s_eff** (slower
kinetics compensating the lost quench). Running the re-scoped F2 campaign
on the ungated prior under a gated production arm would repeat the exact
mistake the Wave-1 probe was built to prevent. Wave 5 re-locates the
landing under the gate before any campaign spend.

## B.3 Run matrix (5 runs)

- **Grid:** `s_eff ∈ {8, 12, 16, 20, 30}` at τ = 6.55 ps,
  `cooling_spatial_gate = density_scaled`.
- **Pinned (probe pins throughout):** κ = 1, `statistical_mixture`,
  f_int = 0.5 (t× ≈ 5 ps on the anchored first shed — the correct pin for
  a staircase-*timing* comparison; the f_int sweep belongs to the
  campaign, B.1 item 2), f_ret = 0.1, λ₀ = 0.9/ps, 9 Å / 0.80 eV / N = 50
  / bridge seed / 30 ps ion window / 1000 ps relaxation cap.
- **Brackets already on disk (no re-run needed):** gated s_eff ∈ {1, 2, 3,
  5} at τ = 6.55 (Wave 4, all overshooting: n_ion_end 2.7–3.5) below, and
  the ungated s_eff sweep (Wave 2) as the quench-on reference.
- **Wiring oracle:** the gated s_eff = 5 / τ = 6.55 point already exists
  from Wave 4 (Δn̄ = 17.52, n_end = 3.48) — a fresh scoring pass must
  reproduce it unchanged.
- **Optional second arm (flagged, not default):** τ = 2.6 ps × the same
  s_eff grid, if the τ = 6.55 arm shows the gated first shed drifting late
  (Wave 4 measured 6.68 ps vs the anchored 5 — cooling slows during exit,
  so the gated crossing lags the ungated t×). Timing re-alignment via
  f_int is the campaign's job, not this probe's.

## B.4 Implementation surface and read-out

**No new code and no new config surface** — every knob is already plumbed
(`build_biphasic_cfg`: `evap_rrk_dof`, `cooling_spatial_gate`; generator
grids `S_EFF_GRID` / `COOLING_GATE_GRID` / `TAU_GRID_PS`; report columns
`s_eff` / `cooling_gate` / `n_relaxed_min` / `frac_frozen`). Two execution
routes, either acceptable:

- **Scratchpad driver** through the delivered `build_probe` / `_run_one`
  pipeline with an overridden module-level grid (the Wave-3 precedent —
  zero repo-code change; run dirs land in the `tier2probe` namespace and
  are scored by the delivered report script), or
- **USER-SETTINGS flip** of the generator's active grid (+ the locked
  grid tests) — a code change, behind `[PROCEED TO IMPLEMENTATION]`.

Same metrics, same reported-not-adjudicated stance. Outcome shapes:

- **(a) A gated landing exists** (some s_eff lands magnitude *and* timing)
  → it becomes the gated campaign prior; feeds the F2 Stage-1 grid
  re-centering and the formal gate-arm promotion decision.
- **(b) Magnitude lands but timing is systematically late** (the ~1.3 ps
  gated crossing lag) → the landing is conditional on a modest f_int
  re-alignment — quantify the lag, hand it to the campaign's f_int
  dimension (B.1 item 2).
- **(c) No constant s_eff lands under the gate** (e.g. every s_eff that
  fixes the magnitude breaks in-window timing or overshoots past 14) →
  tension between the gated arm and the TDDFT staircase prior; since
  TDDFT is not ground truth, this becomes a documented conflict for the
  size-distribution arbitration, not an automatic rejection of either.

## B.5 Campaign-shape implication (recorded for the F-plan, not built)

If Wave 5 lands, the re-scoped Stage 1 becomes **gate = density_scaled
(primary) × s_eff (gated landing band) × τ × f_int**, with a small
`gate = none` sensitivity leg and κ/picture/λ₀ pinned-and-reported. The
identifiability outlook improves: bare fraction ↔ the t×↔ejection race
(f_int·E/Σ), shell-region location ↔ s_eff, timing/tail ↔ τ. Updating
`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (F2 NB) is deferred until the
Wave-5 result exists.

> **NB (post-Wave-6, user conclusion 2026-07-07 — supersedes the "× τ ×
> f_int" grid shape above; findings I18).** τ and f_int must be
> **connected** in any campaign grid: gated observables organize along the
> race margin Δ× = t_eject − t×, so Stage 1 should grid **Δ× directly**
> (f_int derived per τ via f_int = (Σ/E)·e^{t×/τ}), keeping τ as a separate
> dimension only for its independent role (in-bubble leak/quench strength).
> An independent τ × f_int product samples the race incidentally and is
> rejected. The campaign itself stays parked (probe-only scope); OQ-B (the
> bare peak) is explicitly deferred to the F5/production discussion.

---

# Addendum C (2026-07-07) — Wave 6 (f_int flexibility probe) + Wave 7 (detected-read re-score)

> **Status: Wave 6 EXECUTED 2026-07-07** (under the `[PROCEED TO
> IMPLEMENTATION]` trigger; scratchpad-driver route, zero repo-code change;
> 10 new runs → 102 probe dirs; wiring oracle exact on all three f_int = 0.50
> companions). Verdicts: **P2, P3, P4 confirmed; P1 qualitatively confirmed,
> quantitatively refined** (timing re-aligns at f_int ≈ 0.40–0.42, not 0.35 —
> the gated crossing lag grows with t×). New structural finding: **gated,
> f_int is a joint timing + effective-budget knob** (the in-bubble leak over
> [t×, t_eject]); bare I⁺ unreached at any f_int. Full tables + insights
> I15–I17: `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4d; delivery record:
> `drag_migration_log_tier2.md` "Wave-6 f_int probe EXECUTED".
> **Wave 7 EXECUTED 2026-07-08** (under the `[PROCEED TO IMPLEMENTATION]`
> trigger; scratchpad-driver route through the delivered Slice-DS stage,
> zero repo-code change, zero new MD runs): `detection.npz` in all 102 dirs
> at the Sourced t_detect = 8.53 µs, zero P1–P3 guard failures, full
> re-score. Verdicts: **W7-P1, P3, P4 confirmed; W7-P2 refuted in the
> informative direction** (the gated landing point arrives 100 %
> `time_exhausted` at n̄ = 4.09, not floor-frozen — the descent is
> logarithmic). Headline: **the detector read compresses s_eff**; the race
> margin + leak own the arrival n. Full tables + insights I19–I22:
> `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4e; delivery record:
> `drag_migration_log_tier2.md` "Wave-7 detected-read re-score EXECUTED".
>
> Probe-program continuation, explicitly **not** campaign/production work
> (the B.5 campaign shape stays parked). Goal unchanged: understand the
> model — parameter sensibility, flexibility reach (total stripping;
> staircase), now with the detection-time layer. Same stance as all waves:
> existence/capability probes, **reported, not auto-adjudicated**.

## C.1 Wave 6 — the f_int probe (the one knob every wave pinned)

All 92 probe dirs sit at f_int = 0.5. Post-Wave-4 decision B.1(2) promoted
f_int to a first-class lever; Wave 6 actually sweeps it, on the
production-intended gated arm, with three **pre-registered quantitative
predictions** (t× = τ·ln(f_int·E/Σ(21)); E = 0.80 eV, Σ_mix(21) ≈ 0.188 eV;
t_eject ≈ 6 ps; gated crossing lag ≈ +1.8–2.3 ps over t×, from Waves 4/5):

- **P1 (timing re-alignment).** Ungated-equivalent t× at f_int = 0.35 is
  ≈ 2.6 ps; adding the ~2 ps exit lag predicts the gated first shed lands
  ≈ 5 ps — the anchored t★. Tests whether the standing "the ~2.3 ps late
  gated timing is f_int-realignable" claim (Wave-5 finding 4) is
  quantitatively real, and whether the exit lag is itself f_int-independent
  (assumed, now measured).
- **P2 (cliff location, closed side).** f_int = 0.65 at τ = 6.55 puts
  t× ≈ 6.7 ps ≳ t_eject before the lag — predicted to shed **nothing**
  (permanently closed gate, the Wave-4 τ ≥ 16.5 class). Borderline by
  construction: it locates the in-band cliff edge.
- **P3 (re-opening a dead τ-arm).** f_int = 0.25 at τ = 16.5 gives
  t× ≈ 1.0 ps ≪ t_eject — an arm that shed *zero* at f_int = 0.5 (Wave 4)
  is predicted to deep-strip. The race, demonstrated from the f_int side.
- **P4 (floor depth vs the in-bubble leak — the total-strip bound).** The
  n ≈ 2 floor exists because cooling drains E_int *between sheds inside the
  bubble*, i.e. over [t×, t_eject]. f_int moves t× relative to ejection and
  hence the leak: f_int = 0.24 (the row-14 scenario floor; gate open from
  ≈ onset, t× ≈ 0.1 ps) is the max-leak end; f_int with t× ≈ t_eject is the
  min-leak end. The sweep maps floor-n(f_int) and bounds, in-model, how
  close to bare I⁺ the mechanism can get at 0.80 eV — putting a number on
  the structural OQ-B gap.

**Grid (10 new runs; 12 scored points):** gate = `density_scaled`
throughout; τ = 6.55: s_eff ∈ {8, 30} × f_int ∈ {0.24, 0.30, 0.35, 0.65}
(8 new; the f_int = 0.50 pair is on disk from Waves 4/5); τ = 16.5:
s_eff ∈ {8, 30} × f_int = 0.25 (2 new). All other pins unchanged (κ = 1,
mixture, f_ret = 0.1, λ₀ = 0.9/ps, 9 Å / 0.80 eV / N = 50 / bridge seed /
30 ps window / 1000 ps relaxation cap). Execution: the Wave-3/5
scratchpad-driver route (zero repo-code change; `tier2probe` namespace;
scored by the delivered report script). Wiring oracle: the two on-disk
f_int = 0.50 points must re-score unchanged.

s_eff ∈ {8, 30} spans the fast-kinetics and staircase-landing reads so P1
and P4 are read at both; f_int is already a `build_biphasic_cfg` kwarg and
already a tag/report dimension — no new surface anywhere.

## C.2 Wave 7 — the detected-read re-score (needs Slice DS)

**Prerequisite:** the detection stage build per
`TIER2_DETECTION_STAGE_DESIGN.md` (Slice DS — new module + config fields +
tests; behind `[PROCEED TO IMPLEMENTATION]`). t_detect = 8.53 µs is Sourced
(CALIBRATION_MAP row 24).

Then: **re-read every probe dir on disk** (92 + Wave 6's ~10) at the
detector — the detection stage seeds from each dir's existing
`relaxation.npz`, so Wave 7 costs **zero new MD runs**. Deliverable: the
flexibility map *at the detector* — the physically meaningful terminal
reach per (gate, s_eff, τ, f_int) — plus the resolution of Wave 5's open
end: does gated s_eff = 30 arrive at 8.53 µs still mid-shell (~5–6 He) or
at the n ≈ 2 floor? That distinction decides whether the gated arm
contributes mid-shell weight or only deep-strip weight to any eventual
ensemble picture. Report: `n_detect_*` columns + state-reason fractions
(`frozen` / `suppressed` / `time_exhausted`) per design §3.5, and the free
detection-time sensitivity band (§3.3).

**DS acceptance criterion for Wave 7 (recorded now):** adding the two
config fields must leave existing probe dirs loadable (missing new keys →
defaults; verified by an explicit back-compat test at build time). If the
loader rejects them instead, the stale-artifact policy applies and Wave 7
requires regeneration — surface that before executing, do not migrate
artifacts silently.

> **NB (2026-07-07): Slice DS DELIVERED** (log entry "Slice DS DELIVERED";
> design doc as-built block). The acceptance criterion above is **met**:
> `test_pre_ds_cfg_json_loads_with_defaults` proves a pre-DS `cfg.json`
> loads with the fields defaulted — existing probe dirs stay valid, no
> regeneration. The report now carries the `n_detect_*` + state-reason
> columns (populated from an optional per-dir `detection.npz`, `-` when
> absent). **Wave 7 execution remains a separate step** behind its own
> go-ahead: run the detection stage over the ~102 dirs (seeding each from
> its `relaxation.npz` with a detection-enabled cfg view at the Sourced
> t_detect), then re-score.
>
> **NB (2026-07-08): Wave 7 EXECUTED** exactly along this route (all open
> choices adjudicated by the user: scratchpad driver, log-spaced sensitivity
> band 10³–8.53·10⁶ ps, all-102 scope, pre-registered W7-P1..P4). One
> execution wrinkle, recorded in findings §4e: the 63 pre-Wave-4 dirs carry
> the legacy 8.53·10⁶ ps relaxation cap, so the transient cfg view also
> carried the *realized* relaxation duration to satisfy the config-load
> nominal-window bound (the stage's realized-t_h check and the P1–P3 guard
> ran everywhere; on-disk `cfg.json` untouched). Results: findings §4e +
> I19–I22.

## C.3 Sequencing and outcome handling

Wave 6 needs no build and can run first; Slice DS builds in parallel or
after; Wave 7 re-reads everything including Wave 6's dirs. Results go to
`TIER2_STAIRCASE_PROBE_FINDINGS.md` (new sections + insight register), the
decision/delivery record to `drag_migration_log_tier2.md`, and the
discussion resumes from the documented findings — the user adjudicates;
nothing here discharges the F5 gate or touches campaign/production scope.

---

# Addendum D (2026-07-09) — Wave 8: suppressed-fraction cliff anatomy (OQ-B/RQ3 weight-level probe)

> **Status: EXECUTED 2026-07-09** (under the `[PROCEED TO IMPLEMENTATION]`
> trigger; scratchpad-driver route, zero repo-code change; **6** new runs,
> not the planned 7 — step 1 found a **delta-function cliff**
> (E* = 0.461219 eV, ensemble width 4·10⁻¹³ eV; the probe ensemble is
> kinematically congruent), so the three "ECDF-crossing" points collapsed
> into a two-point bracket). Verdicts: **W8-P1 confirmed exactly** (opening
> times 10.83 / 11.96 ps, fractions, |G| windows); **W8-P2, P3 refuted**
> (no width, no 43.5 % crossing — outcome shape (b): heterogeneity must be
> injected); **W8-P4 confirmed sharpened** (the sub-rung sliver puts the
> opened side at n = 1 — the first n = 1 ever expressed); **W8-P5 refuted
> in-model at the mechanism level** (`coulomb_available_eV` is
> bookkeeping-only; K₂₇₀ ≡ K₀₈₀; the absolute-E_int(0) map is
> budget-invariant → new **OQ-H**). Full record:
> `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4f + I26–I28; delivery entry:
> `drag_migration_log_tier2.md` "Wave-8 cliff-anatomy probe EXECUTED".
>
> **Design approved 2026-07-09** (user decision; the optional
> τ = 16.5 Δ×-universality spot point was **dropped** at approval — the
> cross-τ collapse check `E* = Σ·e^(I/τ)` is deferred, not rejected).
> Entry documents:
> findings §4e (I22–I25, the post-Wave-7 fragmentation hypothesis) and
> `RESEARCH_QUESTIONS.md` (RQ1 adjudicated 2026-07-09 — sweeps specified,
> reported, and budget-transferred in **absolute `E_int(0)` [eV]**; RQ3 =
> findings OQ-B). Same stance as all waves: existence/capability probe,
> **reported, not auto-adjudicated**; nothing here discharges the F5 gate.

## D.1 Question and stance

The cascade side is structurally excluded from the experiment's two top
bins (I24: bare needs |G| = 0 exactly; n = 1 sits behind the collapsing
RRK bracket — global min n_detect = 2 in 10 200 detections). If the model
expresses the 43.5 % bare peak at all, it is the **never-opened
(suppressed) side of the Δ× race** under the RQ3 fragmentation
hypothesis. Wave 6 brackets the cliff only coarsely at τ = 6.55:
E_int(0) = 0.40 eV (f_int 0.50) → 0 % suppressed; 0.52 eV (0.65) → 100 %.
Wave 8 asks: **does an in-band E_int(0) split the ensemble ≈ 43.5 %
suppressed / 56.5 % opened — and how wide is the transition?** A
finite-width cliff means a single knob point can carry the bare/shell
*weight ratio*; a step function means bimodality requires ensemble
heterogeneity beyond the fixed-seed spread. Both outcomes are decisive.

**Weight-level only:** `suppressed` stays inert (delivered semantics) and
is *counted* as the bare-candidate class, not evolved. No fragmentation
channel, no new physics, no new `SimConfig` fields. The fragmentation
spec (boil-off, ε-coupling) remains RQ3/RQ2's file.

## D.2 Step 1 — zero-MD cliff anatomy (the E*_i order statistics)

Pre-gate-open no sheds occur, so the gated K2 drain is linear in E_int and
the fragment mechanics are E_int-blind (RQ1 no-partition convention): each
ion's suppression criterion separates from the swept knob. Along its
(fixed) trajectory each ion accrues a **cooling-exposure integral**
`I_i(t) = ∫₀ᵗ ρ_He(r_i(t'))/ρ_bulk dt'` [ps] (converging as the ion exits),
and E_int(t) = E_int(0)·e^(−I_i(t)/τ) until first opening. With the
realized pickup history n_i(t) (Σ non-decreasing under pickup while
suppressed), the per-ion critical value is

**`E*_i = sup_t [ Σ(n_i(t)) · e^(I_i(t)/τ) ]`** — suppressed iff
E_int(0) > E*_i (the sup is effectively the terminal-value product
`Σ(n_i^end)·e^(I_i^tot/τ)`, since e^(I/τ) grows and Σ steps upward).

The pre-open trajectory is what classification needs, and it is shed-free
by construction — so the on-disk **f_int = 0.65 dirs** (100 % suppressed,
zero sheds; the s_eff = 8 and 30 dirs must agree exactly) are the correct
exposure source. Extraction: per-ion r_i(t) and n_i(t) from the stored
artifacts, ρ_He through the same normalized bubble surface the
drag/pickup/cooling gates share. Deliverables (scratchpad analysis, zero
new MD, zero repo-code change):

1. the **E*_i order statistics** (N = 50) → the predicted
   `suppressed_fraction(E_int(0))` ECDF: cliff center, cliff width;
2. **width decomposition**: Poisson-pickup spread in Σ(n_i) vs trajectory
   spread in I_i — which heterogeneity source owns the split;
3. the **43.5 % crossing**: fractions quantize in 2 % steps at N = 50, so
   the natural target is 22/50 = 44 % — E_int(0) placed between the 22nd
   and 23rd order statistics;
4. an **analytic 2.70 eV prescale** (fragment speed ~×√(2.70/0.80) ≈ 1.84
   → shorter exposure; predicted cliff center ≈ 0.3 eV) — rough, pending
   the arm's anchor run (no 2.70 eV trajectory exists yet).

## D.3 Step 2 — targeted MD runs (7 new dirs)

Pins throughout: gated `density_scaled`, τ = 6.55 ps, κ = 1, mixture,
f_ret = 0.1, λ₀ = 0.9/ps, 9 Å, N = 50 bridge seed, 30 ps window, 1000 ps
relaxation cap; `f_int = E_int(0)/E_avail` is the config coordinate.

| Arm | Runs | Spec | Purpose |
|---|---|---|---|
| 0.80 eV transition | 3 | s_eff = 30; E_int(0) at the step-1 ECDF ≈ 20 % / 44 % / 70 % crossings | verify the ECDF where it matters |
| s-blindness companion | 1 | s_eff = 8 at the ≈ 44 % point | classification is decided pre-first-shed → the suppressed **set** (per-ion, not just the fraction) must be identical — wiring-level check |
| 2.70 eV arm | 1 + 2 | s_eff = 30; anchor at E_int(0) = 0.52 eV (f_int ≈ 0.193; predicted fully suppressed → shed-free exposure source), then 2 points at that arm's re-derived transition | first probe data at the production budget; W8-P5 |

Execution route: scratchpad driver through the delivered pipeline (Waves
3/5/6 precedent). Wiring oracles: the on-disk f_int = 0.50 / 0.65 points
re-score unchanged; ledger residual uniform ≈ 2.23·10⁻⁵ eV. Then the
detection stage over the 7 new dirs (Sourced t_detect = 8.53 µs) + full
re-score. **Tag note for implementation:** the probe tag has never carried
the budget (all 102 dirs are 0.80 eV) — the 2.70 eV dirs must be
tag-distinct (budget suffix or equivalent) before generation; verify
non-collision against the existing `_fiX.XX` values.

## D.4 Pre-registered predictions

- **W8-P1 (load-bearing).** The MD suppressed fractions land exactly on
  the step-1 ECDF (same 50 ions, deterministic seed). Any deviation
  falsifies the pre-open separability assumption — candidate culprits:
  pickup statistics differing across f_int (RNG path), trajectory
  feedback the analysis missed.
- **W8-P2.** The cliff has finite width set by the (Σ(n_i), I_i) spread;
  step 1 supplies the number before any run is placed.
- **W8-P3.** A 43.5 % crossing exists inside the sourced band
  [0.2, 0.5] eV (the Wave-6 bracket 0.40–0.52 eV already sits in-band).
- **W8-P4 (shape caveat).** At the split point the opened side arrives
  deep-strip-adjacent (n̄_detect ≈ 3–5 — minimal leak near the cliff):
  the bare/shell **weight ratio** is expressible, the broad experimental
  tail **shape** is not (I12/I21 stand). Report opened-side n_detect vs
  the experimental n ≥ 1 tail (renormalized; Wasserstein per convention)
  as a documented gap, not a target.
- **W8-P5.** The 2.70 eV cliff center sits at lower E_int(0) (~0.3 eV);
  at fixed E_int(0) the suppressed fraction rises with budget → under
  the fragmentation hypothesis **bare grows with budget** — the RQ3-vs-
  RQ4 discriminator (budget-dependent race vs budget-robust ladder
  bottom), now with numbers.

## D.5 Outcome shapes

- **(a) Finite-width cliff, in-band 43.5 % crossing** → OQ-B's weight
  question answered affirmatively at the existence level; feeds the RQ3
  adjudication and the campaign's Δ×-ensemble design.
- **(b) Step-function cliff** (width ≪ the sourced band) → no
  single-point split exists; bimodality *requires* ensemble heterogeneity
  (droplet-size / E_int(0) distribution) — redirects campaign design;
  equally decisive.
- **(c) MD ≠ ECDF** → the separability assumption is wrong; localizes a
  mechanism coupling currently believed absent — a finding in its own
  right (the Wave-1 → Wave-2 pattern).

## D.6 Boundaries

- N = 50, single fixed seed: the ECDF is *that seed's*; fractions
  quantize in 2 % steps. The claim is **existence and location of the
  crossing**, never a fit of the experimental 43.5 %.
- Suppressed = bare remains a **hypothesis** (RQ3 unbuilt); Wave 8 is
  weight-level only and cannot distinguish the fragmentation specs
  (RQ2's ε decides the endpoint distribution, including n = 1).
- Detected reads inherit OQ-E (RRK-evaporation-only flight) and a single
  detection-stage RNG realization per dir.
- Stale-artifact policy (bridge findings §3) applies to the new dirs.

Results go to `TIER2_STAIRCASE_PROBE_FINDINGS.md` (new §4f + insight
register), the decision/delivery record to `drag_migration_log_tier2.md`;
the user adjudicates.

---

# Addendum E (2026-07-09) — Wave 9: the E₀-mixture inversion (implied p(E_int(0)) from the experimental histogram)

> **Status: EXECUTED 2026-07-09** (under the `[PROCEED TO IMPLEMENTATION]`
> trigger; scratchpad drivers through the delivered generator/detection
> pipelines + a pure post-processing fit — zero repo-code change; **8 new MD
> runs** → 116 probe dirs; wiring oracles exact — ledger residual uniform
> 2.231·10⁻⁵ eV on all 8 new dirs, on-disk fi0.50 s8/s30 + fi0.65 re-score
> unchanged). Fit executed under the four confirmed user decisions:
> (1) **simplex-constrained** LS (non-negative + Σw = 1); (2) suppressed→bare
> per-ion under RQ3 spec (b); (3) deep-tail shortfall **reported not forced**;
> (4) W9-P4 a localized column-swap, **no campaign re-scope**. Verdicts:
> **W9-P1 CONFIRMED** (histogram reproduced to L2 = 0.016 / W₁ = 0.086 bins;
> near-cliff density 0.780 %/meV at n = 1 vs predicted 0.77); **W9-P2
> CONFIRMED** (≥ 9.6 % below-floor demand + n = 18–20 beyond the basis floor —
> the droplet-R axis's quantified job); **W9-P3 refined** (p(E₀) smooth and
> unimodal but concentrated **at/above E\* ≈ 0.46 eV**, not an interior 0.45 eV
> bump — **no** multi-modal / electronic-branching signature; interior
> roughness is basis collinearity); **W9-P4 CONFIRMED** (s_eff = 30 refit
> kills the n = 1 bin, L2 → 0.182 / W₁ → 1.398; s_eff detector-identifiable
> via n = 1, recorded as possibility only). Full record:
> `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4g + I29–I32; delivery entry:
> `drag_migration_log_tier2.md` "Wave-9 E₀-mixture inversion EXECUTED".
>
> Design approved 2026-07-09 (user decisions, post-Wave-8 discussion:
> (a) the E₀-mixture probe is worth running *before* any droplet-distribution
> slice; (b) the fitted p(E₀) is treated as a **physics claim** on E_int(0)
> provenance; (c) the s_eff-from-small-n-bins identifiability is **possible but
> not yet fixed** — no campaign re-scope). Entry documents: findings §4f + the
> post-Wave-8 heterogeneity discussion (I26–I28); `RESEARCH_QUESTIONS.md`
> RQ1/RQ3/RQ7. Same stance as all waves: reported, not auto-adjudicated;
> nothing here discharges the F5 gate.

## E.1 Question and idea

Wave 8 proved the arbitration observable cannot come from a single knob
point (delta cliff, I26) — but ions never interact, so **a weighted
mixture of single-point runs is exactly an ensemble with p(E₀) supported
on the grid**. Wave 9 inverts the experimental histogram
(`data/reference/integrated_i_he_abundance.csv`) into the implied
ensemble distribution p(E₀): fit non-negative mixture weights over
empirical per-run detected distributions, read the result as a density
in E₀, and check it against the RQ1-sourced band [0.2, 0.5] eV. The
window map (post-Wave-8 discussion): terminal bins are 22.6 meV-wide E₀
windows below E* = 0.4612 eV (D₀(1)·e^K per rung, flat-bottom regime);
everything above E* is one `suppressed` class (bare under the RQ3
spec-(b) reading).

## E.2 Basis grid (13 columns; 8 new runs)

Primary arm s_eff = 8 (floors *reached*, Wave-8 I27), gated
`density_scaled`, τ = 6.55, probe pins, N = 50 bridge seed, 1000 ps
relaxation cap, detection at t_detect = 8.53 µs. Exact two-decimal f_int
throughout (the Wave-8 tag lesson); E₀ = f_int·0.80 eV:

| f_int | E₀ [eV] | covers | status |
|---|---|---|---|
| 0.58 / 0.65 | 0.464 / 0.52 | suppressed class (the bare column) | on disk (W8/W6) |
| 0.57 | 0.456 | n = 1 window | on disk (W8) |
| **0.53** | 0.424 | n = 2 window | **new** |
| 0.50 | 0.400 | n = 3 window | on disk (W4/5) |
| **0.48** | 0.384 | n = 4 window | **new** |
| **0.45** | 0.360 | n = 5 window | **new** |
| **0.42** | 0.336 | n ≈ 6 | **new** |
| 0.35 | 0.280 | n ≈ 7–8 | on disk (W6) |
| 0.30 | 0.240 | n ≈ 9 | on disk (W6) |
| 0.24 | 0.192 | n ≈ 10–11 | on disk (W6) |
| **0.20** | 0.160 | n ≈ 12–13 | **new** |
| **0.15** | 0.120 | n ≈ 14–15 | **new** |
| **0.10** | 0.080 | n ≈ 16+ (gate open from onset) | **new** |

Plus **one s_eff = 30 companion at f_int = 0.53** (the n = 2 window) —
with the on-disk fi0.57 s8/s30 pair (n_detect 1.00 vs 2.00, the kinetic
wall measured) this gives the clean window-resolved s-sensitivity read
for E.4 P4. Total: **8 new runs**, detection over the new dirs, full
re-score. Deep-window "covers" values are extrapolations — the fit uses
the **measured** per-run n_detect distributions as basis columns, never
the analytic window map.

## E.3 The fit (pure post-processing; no new physics, no new config)

- Basis matrix: per-run empirical n_detect distribution (100 ions,
  1 %-granular) per column; the suppressed columns collapse to one
  bare-class column under RQ3 spec (b).
- Solve non-negative least squares for weights w_j against the
  experimental bins (bare, n = 1 … tail); report per-bin residuals and
  the Wasserstein distance (project convention).
- Read out p(E₀): w_j spread over each column's window width (piecewise
  density); overlay the RQ1 sourced band; report the mass below the
  0.22 eV solvation floor and above E* separately.
- Scratchpad-analysis route (the probe convention); results and the
  fitted weights table go to the findings (§4g).

## E.4 Pre-registered predictions

- **W9-P1:** the top bins (43.5 bare / 17.5 / 8.0 / 5.2 / 4.0 / 3.3 %)
  are reproducible with in-band E₀ mass alone (0.35–0.55 eV region);
  near-cliff implied density ≈ 0.77 %/meV declining ≈ 2.2× over the
  first window.
- **W9-P2:** the deep-shell tail (n ≳ 11) demands E₀ mass **below** the
  0.22 eV sourced floor — the below-floor weight quantifies, in advance,
  how much the droplet-size axis (findings §4f discussion) must supply.
- **W9-P3 (the physics claim, decision (b)):** the implied p(E₀) is
  smooth and unimodal (mode ≈ 0.45 eV, scale ~50 meV). If instead the
  fit demands multi-modal structure, that reads as electronic-branching
  structure (RQ1 fine-structure channels) — a genuine, falsifiable
  provenance prediction.
- **W9-P4 (decision (c) — possible, not fixed):** re-fitting with the
  s_eff = 30 columns shifts the n = 1 window mass into n = 2 (the
  kinetic wall) and degrades the top-bin fit — evidence that the
  detector's small-n bins select the kinetics band. Recorded as an
  identifiability observation only; **no campaign re-scope**.

## E.5 Boundaries

- **Conditional on:** RQ3 spec (b) (suppressed → cleanly bare; an ε-type
  channel would add suppressed-side small-n weight and change the fit),
  the 9 Å exposure K = 0.898297 (RQ7 — all window edges scale with e^K),
  and the fixed-K single-droplet preset (the droplet axis is deliberately
  *absent*; its demand is an output, W9-P2, not an input).
- N = 50 single seed per column (1 % granularity per column; 2 %
  fraction quantization); the flat-bottom window map degrades at high n
  (rung shrink toward the cliff) — hence empirical columns.
- The fit does not discharge the F5 gate and does not build a p(E₀)
  sampling surface — the mixture lives entirely in post-processing.

Results → findings §4g + insight register; decision/delivery record →
`drag_migration_log_tier2.md`; the user adjudicates.

---

# Addendum F (2026-07-10) — Wave 10: E₀-width decomposition (narrow-E₀ + droplet-axis vs broad p(E₀)) + n-graded ladder-bottom candidate

> **Status: Steps 1 + 1.5 EXECUTED 2026-07-10** (zero-MD, scratchpad route,
> under the trigger; results in `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4h +
> insights I33–I36: W10-P1 confirmed, W10-P3 confirmed at the probe K-scale,
> W10-P5 split verdict confirmed → RQ7 co-requisite). **Step 2 (droplet
> slice) and F.5 (n-graded ladder arm) remain NOT built and gated.**
> Designed 2026-07-10 out of the RQ4
> deep-research adjudication discussion (`RESEARCH_QUESTIONS.md` RQ4 findings
> NB register + the flat-vs-tapered/narrow-E₀ follow-on). Document-only:
> Step 1 (F.3) is a zero-MD re-analysis runnable under the trigger; Step 2
> (F.4) **requires the not-yet-built droplet-size slice** and the n-graded
> ladder arm (F.5) requires new config surface — both behind
> `[PROCEED TO IMPLEMENTATION]`, and F.5 additionally behind an external
> many-body I⁺Heₙ calculation to be non-tuned. Entry documents:
> `RESEARCH_QUESTIONS.md` RQ4 (the killed tapered-ladder hypothesis) + RQ1
> (E₀ provenance) + RQ2/RQ3 (ε ≈ 0 adjudication); `TIER2_STAIRCASE_PROBE_
> FINDINGS.md` §4f (Wave-8 delta cliff, I26) + §4g (Wave-9 implied p(E₀),
> I29–I32). Same stance as all waves: existence/capability probe, **reported,
> not auto-adjudicated**; nothing here discharges the F5 gate or touches
> campaign/production scope.
>
> **Amended 2026-07-10 (same day, pre-execution; user-adjudicated cooling/K
> discussion):** the test axis is reframed from droplet-R to the **cooling
> exposure K** (new F.2b — closed-form fate map; both experimental tails from
> one axis); a zero-MD **Step 1.5** semi-analytic K forward model is inserted
> (new F.3b, de-risks the Step-2 slice); the **W10-P3 direction is corrected**
> (bare at narrow E₀ comes from the *small*-droplet / early-ejection / low-K
> end, not deep droplets); and a pre-registered **K-scale / RQ7 leg (W10-P5)**
> is added, with the corresponding second exit in F.8. Decision entry:
> `drag_migration_log_tier2.md` (2026-07-10).

## F.1 Why — RQ4 reverts the tail to the reservoir/race, and the required p(E₀) may be unphysically broad

The RQ4 deep-research pass (2026-07-10) found **no analog supports a
tapered-from-n=1 first-shell ladder** (every studied cation plateaus; the
downslope sits at shell closure, not the bottom rungs). So the small-n
abundance tail is **not** a budget-robust ladder-bottom readout — it reverts
to the reservoir-shape + Δ×-race explanation (Wave-9's implied p(E₀)). That
raises the question this wave targets:

**Is a broad, smooth p(E₀) spanning ≈ 0.16 → > 0.46 eV — the width Wave-9's
mixture fit demanded — physically justified, or is a *narrow* E₀ more
accurate with the width carried by other ensemble axes?**

The two are observationally degenerate on a single-budget histogram (Wave-9
L2 = 0.016 is a mixture fit; interior weights are basis-collinear). Wave 10
breaks the degeneracy by asking which *axis* carries the heterogeneity.

## F.2 The physics motivation (from the RQ1 register, sharpened)

Two quantities were being conflated. The RQ1 "0.2–0.5 eV band" is the
**uncertainty on the *central* value** of E₀, not the **ensemble width** that
a multi-bin tail needs. Term by term, the *spread* is small:

- **Solvation reorganization** — the best-sourced term (NB-RQ1-3) is
  **measured droplet-size-*independent*** (Na⁺: E_disp ≈ 224 meV flat across
  3600–9000 He), hence a **sharp** local first-shell quantity; for I⁺,
  |S_I⁺| ≈ 0.308 eV / Σ(21) ≈ 0.188 eV put it near ~0.25–0.3 eV with narrow
  spread (tens of meV).
- **Electronic / spin–orbit** — **discrete** (³P₀ 0.799, ³P₁ 0.879, ¹D₂
  1.702 eV): it can only add *modes*, not smooth width. **Wave-9 W9-P3
  found p(E₀) unimodal, no second mode** — so this term is small, or
  radiates (Ba⁺ counter-lead), or is unresolved; in no case a smooth-width
  source.
- **Kinematic shell-heating** — the only candidate smooth-broad term — has
  magnitude **refuted 0-3 twice** in RQ1 (genuinely unknown).

So E₀ physics justifies a **narrow, sharp E₀ ≈ E_solv ≈ 0.28 eV**, not a
broad 0.16–0.46 eV distribution — the user's instinct.

**The tension this creates (load-bearing).** The 43.5 % bare peak is the
*suppressed / never-opened* side of the Δ× race (RQ3 spec-b), which requires
**E₀ > E\* = 0.4612 eV**. A sharp E₀ at the solvation scale (~0.28 eV) sits
**below E\***, giving **zero** suppressed → **no bare peak**. So a narrow E₀
can only work if a *different* axis pushes ~43 % of ions into permanent
suppression at E₀ ≈ 0.28 eV — i.e. the width Wave-9 assigned to p(E₀) must be
re-carried by that axis. Identifying whether such an axis exists, and whether
it reaches 43 %, is the whole test.

## F.2b The K reframing (amendment 2026-07-10) — one exposure variable, closed-form fate map

The Wave-8/9 machinery makes the test variable explicit: everything the
droplet axis would do is carried by the **dimensionless in-bubble cooling
exposure**

$$K \;=\; \frac{1}{\tau}\int_0^{t_{\rm eject}} \frac{\rho_{\rm He}(r(t))}{\rho_{\rm bulk}}\,dt
\qquad\Rightarrow\qquad E_{\rm int}(t_{\rm eject}) = E_0\,e^{-K},$$

measured **K = 0.898297** at the pinned 9 Å condition (Wave 8; one number for
the whole ensemble by kinematic congruence, I26). Droplet radius (via
t_eject), τ, and the production kinematics (RQ7 — faster fragments, earlier
ejection) are three *handles on the same variable*. At fixed sharp E₀ the
gated fate structure is closed-form in K, with the suppression edge at
**K\*(E₀) = ln(E₀/Σ(21))** and the post-crossing in-bubble leak (no-shed
approximation) `Σ(21)·(1 − e^{−(K−K*)})`:

| K vs K\*(E₀) | in-bubble leak | fate |
|---|---|---|
| K < K\* (under-cooled) | — (never crosses) | suppressed → **bare** (RQ3-adjudicated total shed, G > 0) |
| K ≳ K\* (crosses just before ejection) | ≲ D₀-scale | deep strip, **n ≈ 1–3** (the n = 1 sliver = within ~one leak-rung above K\*) |
| K ≫ K\* (crosses early) | large | shallow strip, **mid/high-n** |

Cooling is the *only* process that can pull E_int below Σ before ejection
(evaporation is suppressed above Σ), so the bare class is exactly the
**under-cooled** class — "less cooling → more bare" holds, as a **cliff, not
a dial** (at fixed (E₀, K) the split is 0 % or 100 %; a smooth fraction needs
heterogeneity in E₀ or K). The no-shed leak formula is already MD-validated
at the pinned K by the Wave-6 relaxed floors (f_int = 0.24: leak ≈ 0.110 eV
→ n ≈ 11 vs measured 10.96; f_int = 0.50: leak ≈ 0.025 eV → n ≈ 2–3 vs
measured 2.7–3.1).

**Consequence (both tails from one axis).** A *distribution* over K
straddling K\* produces the bare peak (below-K\* mass), the small-n bins
(just-above-K\* mass), *and* the retained mid/high-n weight (K ≫ K\*) — the
two demands Wave 9 assigned separately (super-E\* E₀ mass for bare; the
droplet axis for the ≥ 9.6 % below-floor deep-shell weight) are opposite
ends of one K spread. Wave 10's question therefore sharpens to: **does a
physically plausible K-distribution straddle K\*(E₀ ≈ 0.28 eV) ≈ 0.398 with
≈ 43.5 % below it, while its upper side reproduces the retained tail?** At
the probe K-scale (0.898) the bare side needs droplets at roughly *half* the
pinned exposure — see the corrected W10-P3 and the K-scale leg W10-P5.

## F.3 Step 1 — zero-MD spread decomposition (which axis is already delta)

Wave 8 established that at **fixed E₀ and fixed droplet** the suppressed/opened
split is a **delta cliff** (E\* = 0.461219 eV, ensemble width 4·10⁻¹³ eV,
I26) — the Poisson-pickup and trajectory-exposure spreads across the N = 50
ions do **not** smear it. Step 1 makes this rigorous *as a width-budget
statement* and extends it to the opened side, purely by re-analysing on-disk
Wave-8/9 dirs (zero new MD, zero repo-code change, scratchpad analysis):

1. **Suppressed-split width at fixed (E₀, droplet):** re-confirm the Wave-8
   delta from the E\*_i order statistics (already extracted) — restate as
   "Poisson + trajectory carry ≈ 0 tail width at fixed droplet."
2. **Opened-side terminal-n width at fixed (E₀, droplet):** for each Wave-9
   basis column (one E₀, N = 50, detected at 8.53 µs) measure the n_detect
   spread. Expectation from I21/W8-P4: near-delta (deep-strip-adjacent at the
   cliff), i.e. a single E₀ populates ≈ one detected bin.
3. **Conclusion:** if both are near-delta, then at fixed droplet the entire
   Wave-9 tail width came from the **E₀ mixture axis** — so under a *narrow*
   E₀ the only remaining carrier is the **droplet-size/radius (K) axis**,
   which every probe dir holds fixed. Step 1 quantifies exactly how much
   width must be re-sourced (the full Wave-9 span) and thereby *defines the
   Step-2 target*.

## F.3b Step 1.5 — semi-analytic K forward model (zero MD; amendment 2026-07-10)

Between the width bookkeeping (Step 1) and the gated droplet slice (Step 2):
test W10-P2/P3 *on paper* with the F.2b closed form, before any new surface
is built.

1. **Calibrate K(R_droplet).** One MD anchor exists (K = 0.898297 at the
   ~27.9 Å probe droplet); extend by transit scaling t_eject ≈ R/v with the
   density-taper correction read from the stored probe trajectories r(t)
   (scratchpad re-analysis of existing dirs; no new MD).
2. **Propagate a droplet prior.** For stated candidate droplet-size
   distributions (the experimental log-normal prior — stated, not silently
   pinned, per F.7), compute the implied K-distribution and:
   - **bare fraction** = P(K < K\*(E₀)),
   - the **opened-side histogram** via the leak → terminal-n map (F.2b
     table, validated on the Wave-6 anchors), with the kinetic/detected
     correction cross-checked against the Wave-9 `s_eff = 8` basis columns,
   - the **Wasserstein distance** to `integrated_i_he_abundance.csv`.
3. **Scan** sharp E₀ ∈ ≈ [0.24, 0.32] eV × the two K-scales of W10-P5
   (probe 0.898 / production-kinematics estimate ≈ 0.49).

**Output, either way decisive for sequencing:** (a) no plausible droplet
prior reaches ≈ 43 % bare (the analytic W10-P3 verdict) → the Step-2 slice
is demoted to tail-*shape* duty only; or (b) a **(median, width) target**
the Step-2 slice must hit — the slice is then built to *verify* a
quantified prediction, not to explore. Runnable under the trigger
(scratchpad route, zero repo-code change), same stance as Step 1.

**What the analytic pass cannot see (recorded for Step 2):** incomplete
shells in small droplets (n_eject < 21 → smaller Σ(n_eject) → E\* lower
still, compounding toward bare), the pickup/drag-exposure coupling to final
speeds (a VMI side-constraint), and non-ballistic exit. The closed form
brackets; the MD slice arbitrates.

## F.4 Step 2 — the droplet-axis substitution test (needs the droplet slice)

**Prerequisite (new build, gated):** a droplet-size axis — a distribution over
droplet radius / K (leak depth) and the associated ejection-timing/exposure —
wired through the biphasic + gated-cooling + detection pipeline. The existing
`single_pulse_droplet_distribution` preset is a candidate reuse for the
radius sampling, but the biphasic-mass + `density_scaled`-cooling + detection
wiring is new surface (behind `[PROCEED TO IMPLEMENTATION]`, its own slice
plan and back-compat/namespace tests per the Slice-DS precedent). This is the
droplet-R axis W9-P2 already named as "the droplet-size axis's quantified
job" for the ≥ 9.6 % below-floor demand — Wave 10 promotes it from *demand*
to *test*.

**Design:** hold **E₀ fixed and sharp** (anchor E₀ ≈ 0.28 eV ≈ E_solv, plus a
narrow ±spread sensitivity leg), sweep only the **droplet-radius distribution**
(K axis), score the terminal I⁺Heₙ distribution at the detector, and ask:

- Does the spread in K (deeper leak in larger droplets → shallower strip /
  higher n; earlier ejection in smaller droplets → low K → suppression,
  F.2b) reproduce the experimental tail *shape* (bare → n = 1 → … monotone
  envelope) at narrow E₀?
- **Can it reach ≈ 43 % suppressed (bare) at E₀ ≈ 0.28 eV** — i.e. does the
  droplet axis push a large fraction of ions to never open the gate before
  ejection, substituting for the > 0.46 eV E₀ mass Wave-9 required? Per
  F.2b this is the **small-droplet / low-K side** — K < K\*(0.28 eV) ≈ 0.40
  vs the pinned 0.898; Step 1.5 quantifies the required prior before this
  slice is built, and Step 2 is built only against a live Step-1.5 target
  (F.3b output (b)).

Reuse everything downstream: the delivered detection stage (Slice DS), the
`tier2_staircase_probe_report` metrics + state-reason fractions, the
Wasserstein comparison against `integrated_i_he_abundance.csv`.

## F.5 The n-graded ladder-bottom — a *bounded* candidate arm for the residual only

RQ4 leaves one escape hatch open: **open-shell I⁺ (³P) could taper its first
rungs where closed-shell alkalis plateau.** The defensible, *capped* version
is electronic, not a free taper, and lives inside the existing
`ladder_electronic_picture` fork:

The He–I⁺ interaction has a deep ³Π (X₂, D₀ ≈ 13.3 meV / 106.9 cm⁻¹) and a
shallow ³Σ⁻ (≈ 7.9 meV / 63.6 cm⁻¹) ([IHe05]). The **first** He binds on the
deepest lobe (X₂); as the anisotropic ³P shell fills, later He are forced by
geometric/electronic frustration toward the statistical mixture (~9.2 meV) and
the shallow ³Σ⁻. That gives a **mild, ab-initio-capped** first-rung taper:

$$D_0(1):D_0(2):D_0(3)\ \approx\ 13.3:9.2:7.9\ \text{meV}\ \approx\ 1.4:1.0:0.86,$$

a first-rung elevation of at most **~1.4–1.7×** (X₂ / mixture … X₂ / ³Σ⁻) —
**not** the 2.2× the tail wants. So it is a **contribution of a known, bounded
size**, invoked only for a residual n = 1/n = 2 step that Steps 1–2 leave
uncovered, never the whole tail. It is an **n-graded picture** (X₂ for the
first, least-frustrated atom relaxing to the mixture), a refinement of the
existing three-arm fork, not a new mechanism.

**Non-tuned only with the calculation.** Fixing *where* the X₂→mixture
transition sits from the histogram bins would relocate Wave-9's overfitting
from the reservoir to the ladder (the RQ4-doc warning). To stay principled the
transition must come from a **many-body I⁺Heₙ frustration calculation** (how
many He fit on the deep ³Π lobe before the mixture takes over) — RQ4 open
question 1, an external ab-initio / PIMC / ⁴He-DFT item. Without it the graded
picture is *motivated but a knob*, and must be reported as such.

**Sequencing discipline (keeps the ladder honest):** (1) Steps 1–2 with the
**flat 9.2 meV bottom** first — clean baseline; (2) read the **residual**
n = 1/n = 2 gap the droplet axis + flat ladder cannot cover at narrow E₀;
(3) *only if* that residual is ≲ 1.4× (coverable by the capped electronic
taper) does the graded arm earn its place — matched to the ab-initio ceiling,
not fit to the bins. A residual larger than the cap points instead to
non-zero ε (reopening RQ2/RQ3) or a missing mechanism.

> **Floor-variant usage note (2026-07-11, Addendum H.2b amendment):**
> the **transition-at-n=1 floor variant** (D₀(1) = 13.3 meV X₂;
> D₀(n ≥ 2) = mixture) is **knob-free** — both depths are
> [IHe05]-sourced and the picture argument pins the first,
> least-frustrated atom unconditionally — so it is exercisable *now* as
> a bounded lower-bound lever (H.2b L3) without violating this clause.
> What remains calculation-gated is only the *transition location*
> (how many He fit the deep ³Π lobe): slid-transition variants are the
> knobbed family, reported separately, and the RQ4 external many-body
> calculation returns as their **arbiter**. Under the 2026-07-11
> solvated-branch re-targeting the residual (~2.2×) exceeds this cap
> alone — the taper now stacks with the Addendum-H geometric levers
> (position, n_eject(depth)) instead of carrying the residual by itself.

## F.6 Pre-registered predictions

- **W10-P1 (Step 1, load-bearing).** At fixed (E₀, droplet) both the
  suppressed split and each opened column's n_detect are near-delta
  (≤ 1 bin), confirming Poisson + trajectory carry negligible tail width →
  under narrow E₀ the droplet-K axis is the sole remaining carrier.
- **W10-P2 (Step 2, the decomposition verdict).** At sharp E₀ ≈ 0.28 eV the
  droplet-radius distribution *can* reproduce the small-n tail *shape*
  (larger droplet → deeper leak → deeper bin). Failure mode: the tail comes
  out too narrow (K-spread insufficient), falsifying narrow-E₀-carries-tail.
- **W10-P3 (the bare-peak tension, decisive; direction corrected
  2026-07-10).** At E₀ ≈ 0.28 eV < E\* the droplet axis **cannot** reach
  ≈ 43 % suppressed (bare) *at the probe K-scale* — because suppression
  needs the gate to never open before ejection, i.e. **K < K\*(0.28 eV) ≈
  0.398 vs the pinned 0.898**: the bare side is the *small*-droplet /
  early-ejection / slow-cooling **low-K** end of the distribution (the
  original "implausibly deep/slow-cooling droplets" wording had the deep
  half inverted — deep droplets *raise* K and serve the shallow-strip
  mid/high-n weight instead), and putting ≈ 43 % of the experimental
  droplet prior at roughly *half* the pinned exposure is not expected. If
  confirmed **at both K-scales (W10-P5)**, **narrow E₀ is falsified for the
  bare peak**: the ≥ 43 % bare mass genuinely needs E₀ > 0.46 eV, i.e. a
  smooth super-solvation source (kinematic coupling, magnitude open) or a
  larger-than-Na⁺ I⁺ reorganization — back to RQ1. If *refuted* (droplet
  axis does reach 43 % at narrow E₀), narrow E₀ stands and is the more
  physical model.
- **W10-P4 (ladder residual).** After Steps 1–2 with the flat bottom, the
  residual n = 1/n = 2 step is ≲ 1.4× (electronic-taper-coverable) rather
  than the full 2.2×, i.e. the droplet axis carries most of it and the
  graded picture only closes a bounded gap.
- **W10-P5 (the K-scale / RQ7 leg; added 2026-07-10).** The W10-P3 verdict
  is **K-scale-conditional**. At the probe scale (K₀ = 0.898,
  E\* = 0.4612 eV) narrow E₀ ≈ 0.28 eV sits far below the cliff and P3 is
  predicted to **confirm** (fail to reach 43 % bare); at the
  production-kinematics estimate (fragment speed ×√(2.70/0.80) ≈ 1.84 →
  K₀ ≈ 0.49, E\* ≈ 0.31 eV) the cliff lands **near the ensemble median**
  (K\*(0.28 eV) ≈ 0.40 vs K₀ ≈ 0.49) and a *modest* droplet width plausibly
  splits ≈ 43/57 — narrow E₀ **succeeds**. If Step 1.5 confirms this split
  verdict, the resolution is "**narrow E₀ + RQ7 kinematics**": the RQ7 arm
  becomes a **co-requisite** of any Step-2 adjudication (a P3 confirmation
  at the probe K-scale alone is inconclusive), and the droplet slice's job
  shrinks to the *width*, not the location, of the K-distribution. (The
  0.49 estimate is a speed-scaling bracket, not model output — RQ7 stays
  unmodeled; no 2.70 eV MD is implied by this leg.)

## F.7 Boundaries and gating

- **Step 1** is zero-MD re-analysis (runnable under the trigger, scratchpad
  route, `tier2probe` namespace untouched). **Step 2 and F.5 are new surface**
  behind `[PROCEED TO IMPLEMENTATION]`; F.5 additionally behind the external
  many-body I⁺Heₙ calculation.
- **Step 1.5** is zero-MD re-analysis + closed-form propagation (runnable
  under the trigger, scratchpad route). Its K-scale leg is **analytic only**
  — the production K₀ ≈ 0.49 is a kinematic bracket, not a modeled arm; RQ7
  stays open and gated, and no 2.70 eV MD is run.
- **Conditional on** RQ3 spec-(b) (suppressed → bare; an ε channel would add
  suppressed-side small-n weight and change the read — reopening RQ2/RQ3),
  the ε ≈ 0 adjudication, and the 9 Å exposure K = 0.898297 (RQ7 — all window
  edges scale with e^K; the 2.70 eV budget is a separate arm, RQ7-gated).
- The droplet axis is the deliberately-absent input of every prior wave
  (Wave-9 E.5); introducing it is the point of Step 2, but it also imports a
  droplet-size *distribution* assumption (experimental log-normal vs a
  modelled prior) that must be stated, not silently pinned.
- N and seed conventions inherit the probe defaults; the claim is
  **existence and axis-attribution of the tail width**, never a production fit.
- Stale-artifact policy (bridge findings §3) applies to any Step-2 dirs.

## F.8 Outcome shapes

- **(a) Narrow E₀ + droplet axis reproduces the tail *and* reaches ~43 %
  bare** (W10-P2 confirmed, W10-P3 refuted) → option (1) holds in its more
  physical form; the "broad p(E₀)" was a reduction artifact; RQ1 burden
  eases (sharp E₀), droplet-R slice becomes the campaign's width carrier.
- **(b) Droplet axis reproduces the *tail shape* but cannot reach 43 % bare
  at narrow E₀** (W10-P3 confirmed) → the bare peak genuinely needs E₀ mass
  above E\*; narrow E₀ is falsified *for the bare class specifically*, and
  the > 0.46 eV smooth source becomes the open RQ1 item (kinematic coupling
  magnitude, or I⁺-specific super-solvation) — a real, localized gap.
- **(c) Droplet axis cannot carry the tail width at all** → the width must be
  E₀-intrinsic after all (broad p(E₀) stands despite the physics argument),
  or a mechanism is missing (non-zero ε, in-flight channels RQ5); the most
  consequential outcome — it would mean neither the ladder (RQ4-flat) nor a
  physically-narrow E₀ nor the droplet axis explains the tail.

Any of the three is decisive for whether option (1) (accept flat ladder,
tail = reservoir/race) survives, and localizes the next research target
(RQ1 super-solvation, the droplet-size slice, or the RQ2/RQ3 ε reopening).

**Amendment (2026-07-10) — outcome (b) has two exits, separated by W10-P5.**
If the droplet axis fails to reach 43 % bare at the probe K-scale but
succeeds at the production K-scale, the bare peak needs **RQ7 kinematics**
(earlier production ejection → lower K → lower E\*), not super-E\* E₀ mass:
"narrow E₀ + RQ7 kinematics" stands (effectively outcome (a) in its
RQ7-flavoured form) and the RQ1 super-solvation burden never fires. Only
failure at **both** K-scales routes to the original (b) (super-solvation /
kinematic-coupling E₀ source).

**Cross-links:** `RESEARCH_QUESTIONS.md` RQ4 (killed taper) / RQ1 (E₀
provenance, the narrow-vs-broad argument) / RQ2/RQ3 (ε ≈ 0) / RQ7 (budget /
K axis); `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4f (Wave-8 delta) + §4g
(Wave-9 p(E₀)); `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`
(`ladder_electronic_picture` fork, [IHe05] X₂/³Σ⁻ rungs, Form-U bottom);
`CALIBRATION_MAP.md` (E_int(0) Bounded absolute, row 16). Results →
`TIER2_STAIRCASE_PROBE_FINDINGS.md` (new §4h + insight register) +
`drag_migration_log_tier2.md`; the user adjudicates.

---

# Addendum G (2026-07-11) — Wave 11: RQ7 production-kinematics probe (measure K₂.₇₀)

> **Status: EXECUTED 2026-07-11** (triggered and delivered same day; the
> opened companion promoted to default pre-execution, user adjudication
> 2026-07-11). Results: `TIER2_STAIRCASE_PROBE_FINDINGS.md` **§4i** +
> insights I37–I40 — **K₂.₇₀ = 0.746 measured** (W11-P1/P3 confirmed,
> W11-P2/P5 **refuted as registered**, W11-P4 split: overlap 8 % in-band,
> beyond-band share 41.7 % — the load-bearing validity caveat); the
> fate-map landing re-calibrates to E₀ ≈ 0.38–0.41 eV (in-RQ1-band,
> untuned W₁ = 0.496); closed form validated at production kinematics
> (0.24 He). Delivery record: `drag_migration_log_tier2.md` 2026-07-11.
> This is Move 1 of the post-Wave-10
> program order (`drag_migration_log_tier2.md` 2026-07-11). Entry documents:
> `RESEARCH_QUESTIONS.md` RQ7 (promoted status block) +
> `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4h (I34–I36, the split verdict) +
> the F.2b/F.3b machinery this wave re-uses. Same stance as all waves:
> a **measurement**, reported, not auto-adjudicated; nothing here
> discharges the F5 gate or touches campaign/production scope.

## G.1 Why — K₂.₇₀ is the single number that decides the Wave-10 branch

Wave 10 (W10-P5, I35/I36) left the narrow-E₀ reading of the experimental
histogram standing on exactly one unmeasured quantity: the **production
cooling exposure K₂.₇₀**. At the pinned 9 Å exposure (K = 0.898) narrow E₀
is excluded ~19×; at the ballistic speed-scaling bracket (K₀ ≈ 0.49 =
0.898·√(0.80/2.70)) it lands the bare peak untuned. The bracket is
kinematics-free extrapolation, not model output (§4h boundary 1). The
delivered model can *measure* this number: the Coulomb mechanics already
integrate any initial separation, the drag bundle is locked, and the
Wave-10-validated exposure integral reads K off stored trajectories with
0.08 % discretization error against the Wave-8 bisection anchor. One
deterministic dir replaces the bracket with a measurement.

**Decision structure (pre-registered):** K₂.₇₀ near ~0.5 → the Wave-10
landing stands ("narrow E₀ + faster ejection"); near the probe 0.9 → the
landing inverts and the RQ1 super-solvation route reopens (F.8 original
outcome (b)); in between → the F.3b scan re-run at the measured K
adjudicates (zero-MD, Step 3 below).

## G.2 The measurement — definition, channel, dimensional analysis

**Exposure integral (F.2b, unchanged):**

$$K \;=\; \frac{1}{\tau}\int_0^{t_{\rm exit}} \hat\rho\big(r(t)\big)\,dt,
\qquad \hat\rho = \rho_{\rm He}/\rho_{\rm bulk},$$

integrated over the stored per-ion trajectory until droplet-surface exit.
Units: τ [ps], dt [ps], ρ̂ [–] → **K dimensionless** ✓. K\*(E₀) =
ln(E₀ [eV] / Σ(21) [eV]) [–] ✓. K is defined *per the same τ* that the
Newton-cooling term uses (pinned mid-band τ = 6.55 ps), so E_int(t_exit) =
E₀·e^(−K) holds by construction and the ratio S_K = K₂.₇₀/K₀.₈₀ is
τ-independent.

**Channel:** fragments born at the I₂ ground-state separation
R₀ = 2.666 Å; the same Coulomb formula both budgets already share,
`14.4 eV·Å / R₀`: 14.4/2.666 = **5.40 eV shared = 2.70 eV/fragment**
(I⁺+I⁺, `E_coulomb_scale = 1`); the 9 Å condition is 14.4/9 = 1.60 eV
shared = 0.80 eV/fragment. Config coordinates: `R0_GS_angstrom = 2.666`
(the one genuinely new knob value) + `coulomb_available_eV = 2.70`
(provenance stamp, no hard refuse per plan §8).

**Speed scales (pre-registered estimates, 1 eV = 9648.5 amu·Å²/ps²):**
ballistic asymptotes v∞ = √(2E/m): bare I⁺ (126.9 amu) 11.0 → 20.3 Å/ps
(0.80 → 2.70 eV); dressed I⁺He₂₁ (211.0 amu) 8.6 → 15.7 Å/ps. Either mass
convention gives the same ratio **√(2.70/0.80) = 1.84** — the "~1.8×
speed" the validity argument (G.4) must cover.

## G.3 Run matrix (2 default dirs + one optional leg; zero-cost control)

Pins throughout (identical to the Wave-8/9 probe convention): gated
`density_scaled`, τ = 6.55 ps, s_eff = 8, κ = 1, mixture, f_ret = 0.1,
λ₀ = 0.9/ps, probe droplet (~27.9 Å / N = 2000 He), N = 50 bridge seed,
30 ps window, 1000 ps relaxation cap; absolute E_int(0) is the physics
coordinate (`f_int` the config coordinate).

| Arm | Runs | Spec | Purpose |
|---|---|---|---|
| control (zero new MD) | 0 | re-read the on-disk 9 Å suppressed dir (fi = 0.65) through the Wave-10 exposure script | wiring oracle: reproduce K = 0.897541 (trajectory route) |
| **production-kinematics primary** | 1 | `R0_GS_angstrom = 2.666`, budget stamp 2.70 eV, **E_int(0) = 0.52 eV** (f_int ≈ 0.193 — the Wave-8 2.70 eV anchor value) | measure **K₂.₇₀**: suppressed for any K < K\*(0.52) = 1.02 → shed-free, mass-constant ride, kinematically congruent with the Wave-8 anchor convention (m = I⁺He₂₁ throughout) |
| **opened companion (default; promoted 2026-07-11)** | 1 | same geometry, E_int(0) = 0.28 eV (the Wave-10 landing E₀) | bounds the m(t)-feedback on K (opened cascade lightens the ride); first opened data point at production kinematics |
| fixed-mass leg (optional) | 1 | literal Tier-0 convention: `mass_scenario = fixed` at m_eff | convention-sensitivity of K to the mass convention (anchor rode m = 21 He, Tier-0 pins 19) |

The **primary deliverable** is the pair (K₂.₇₀, S_K = K₂.₇₀/K₀.₈₀) under
anchor-congruent conventions, plus t_exit and the exit/detector speeds.
Pre-registered expectation: **K₂.₇₀ ∈ ≈ [0.5, 0.7]** — drag eats part of
the extra speed, so above the ballistic 0.49, below the probe 0.898.

**Tag note (extends D.3):** these are the first dirs off 9 Å geometry —
the tag must carry budget *and* geometry (e.g. an `_R2.67` suffix on top
of the Wave-8 budget suffix); verify non-collision against all existing
probe tags before generation.

**RNG:** no new sampling surfaces. The onset is deterministic and
center-placed; the suppressed ride has no evaporation draws and no pickup
(shell full at n\* = 21, occupancy cap → rate 0). Draw order untouched
(forbidden-list item respected structurally, not by exception).

**Guards:** the `mass_scenario`↔`drag_coefficients` §6.5 guard trips as in
every gated probe run (`allow_inconsistent_mass_pairing = True`, §6.6
mid-window defense unchanged). New check: the config-load path must
*accept* `R0_GS_angstrom = 2.666` (confirm no validation bound or
geometry guard assumes R₀ ≳ several Å) — if a guard refuses, that is a
finding, not something to silently relax.

## G.4 Validity arguments (required by the Move-1 spec)

**(i) Drag-calibration band at ~1.8× speed.** The pure-cubic γ = g·b·v²
(force g·b·v³) was calibrated in the 9/18 Å TDDFT windows, whose in-window
speeds top out at the 0.80 eV scale (~11 Å/ps bare-mass equivalent). The
2.666 Å run spends its late transit at up to ~1.84× that. Direction
argument: a cubic force law is *steeper* than classical form drag
(F ∝ v²), so extrapolation beyond the calibrated band plausibly
**over-drags** → slower transit → **K biased high**. The wave does not
resolve the extrapolation (that would need TDDFT at production
kinematics); it **quantifies the exposure at risk**: report the
speed-band decomposition of K — the fraction of ∫ρ̂ dt accrued while the
ion's speed exceeds the 9 Å-window maximum. Small share → the measured
K₂.₇₀ is band-safe; large share → outcome (d).

**(ii) The n\*(R₀) opening check.** At R₀ = 2.666 Å the two ions sit far
inside one first-shell radius (~4.67 Å, GAH25) — the t = 0 full-dressing
convention (each ion dressed with n_eject = 21, §4h boundary 3) is
geometrically a *shared* bubble, and the model's per-ion drag environment
books ρ̂ ≈ 1 between the fragments where the real density is depleted.
Check: report the overlap-segment share of K (exposure accrued while the
pair separation < 2×4.67 ≈ 9.3 Å). Pre-registered estimate: small
(≈ few %) — the Coulomb ramp at 2.666 Å reaches ~1 eV within ~1.3 Å, so
the overlap segment lasts ~0.1–0.2 ps of a ≈ 4 ps exposure. Direction:
the booked density is an over-count there → **K biased high**.

Both identified biases point *up*: the measured K₂.₇₀ is an **upper
bound** in these two respects. Consequence for the decision structure: if
even the upper bound lands ≤ ~0.7, the narrow-E₀ reading survives
conservatively; an inversion verdict (K ≈ 0.9) must first survive the two
decomposition reports before it kills the Wave-10 landing.

**(iii) Integrator stability at the steeper Coulomb onset.** The 2.666 Å
start multiplies the initial Coulomb force ~11× ((9/2.666)²); the fixed dt
was never exercised there. Opening oracle: the standard ledger residual +
a drag-off control integration against the analytic ballistic limit
(v∞ table in G.2, exact to discretization). A residual blow-up is a dt
finding and blocks the K read.

## G.5 Steps

1. **Step 0 — opening checks (with the primary run):** config-load
   acceptance at 2.666 Å; drag-off ballistic oracle + ledger residual
   (G.4 iii); tag non-collision.
2. **Step 1 — the measurement:** primary dir; K₂.₇₀ via the Wave-10
   exposure script (scratchpad route, Waves 3/5/6/8 precedent); control
   re-read as wiring oracle.
3. **Step 2 — diagnostics from the same trajectories:** speed-band and
   overlap-segment decompositions of K (G.4 i–ii); exit and detector
   speeds vs `data/reference/vmi_iplus_he.csv` (the independent
   discriminator — pins the kinematics with **no evaporation physics**;
   state the mass convention of the compared species explicitly, dressed
   complex vs detected small-n fragments under RQ3 sequential shed with
   ε ≈ 0 ⇒ speed ≈ preserved).
4. **Step 3 — zero-MD verdict re-read:** re-run the F.3b scan
   (E₀ ∈ [0.24, 0.32] eV × stated droplet priors) with the **measured**
   K₂.₇₀ replacing the ×0.545 bracket; read bare fraction + W₁ against
   `integrated_i_he_abundance.csv`. The optional fixed-mass leg only if
   Step 1's decompositions or the scan sit near a boundary.

## G.6 Pre-registered predictions

- **W11-P1 (wiring oracle).** The exposure script reproduces the 9 Å
  control at K = 0.897541 exactly (same code path, same dir); the primary
  run's ledger residual matches the probe-standard ≈ 2.2·10⁻⁵ eV scale.
- **W11-P2 (load-bearing).** **K₂.₇₀ ∈ [0.5, 0.7]** — above ballistic
  (drag eats part of the extra speed), well below the probe 0.898
  (faster transit dominates).
- **W11-P3 (congruence).** The suppressed production-kinematics ensemble
  is kinematically congruent like Wave 8 (per-ion K spread ≲ 10⁻¹²); the
  opened companion's K differs from the primary only at the few-%% level
  (m(t)-feedback small).
- **W11-P4 (validity decompositions).** Beyond-band exposure share and
  overlap-segment share are each ≲ 10–15 % of K; both corrections point
  down (K measured = upper bound).
- **W11-P5 (the verdict re-read).** At the measured K₂.₇₀ the F.3b scan
  still lands the bare peak from a stated prior within
  E₀ ∈ [0.24, 0.32] eV (bare ≈ 43.5 % crossed, W₁ ~ 1 bin) — the Wave-10
  landing survives measurement replacing the bracket.

## G.7 Outcome shapes

- **(a) K₂.₇₀ ≲ 0.7 and the scan lands** (P2+P5 confirmed) → "narrow E₀ +
  RQ7 kinematics" is confirmed at the measurement level; the ballistic
  bracket is retired; the RQ7 kinematics *number* is closed (the modeling
  arm stays unneeded for the K question); Step-2/Move-3 targets are
  re-quantified at the measured K.
- **(b) K₂.₇₀ ≈ 0.8–0.9** (survives the G.4 decompositions) → the Wave-10
  landing inverts; F.8-original outcome (b) fires: super-E\* E₀ mass
  (RQ1 super-solvation / kinematic coupling) becomes the open bare-peak
  source again.
- **(c) intermediate / boundary case** → the scan adjudicates per prior;
  if priors split, the droplet slice (Move 3) inherits the tie-break and
  must carry K-width as a fitted, reported quantity.
- **(d) validity failure** — beyond-band exposure share large, overlap
  share large, or the dt oracle fails → K₂.₇₀ is not measurable inside
  the locked calibration; the honest exit is a drag-recalibration
  research item (TDDFT at production kinematics), *not* a silent
  extrapolation. The bracket stays the working number, flagged.

## G.8 Boundaries and gating

- **The MD leg is new surface behind `[PROCEED TO IMPLEMENTATION]`**
  (two default dirs + the optional fixed-mass leg; the control and Step 3
  are zero-MD re-analysis, scratchpad route, `tier2probe` namespace
  untouched). Trigger given 2026-07-11.
- The drag law stays 0.80-eV-band-calibrated: this wave replaces the
  *ballistic bracket* with an *in-model measurement*; it does not
  recalibrate γ. The G.4 decompositions bound, not remove, the
  extrapolation risk.
- n_eject = 21 for all geometries stays a convention (§4h boundary 3);
  the n\*(R₀) check documents its production-geometry cost, it does not
  fix it. Suppressed detection semantics remain spec-(a) until Move 2.
- No production run, no fragmentation channel, no droplet slice — this
  wave measures one number and re-reads one scan.
- N = 50 / bridge-seed / single-droplet conventions inherit the probe
  defaults; the claim is **the location of K₂.₇₀ relative to the two
  K-scales**, never a production fit.
- Stale-artifact policy (bridge findings §3) applies to the new dirs.

**Cross-links:** `RESEARCH_QUESTIONS.md` RQ7 (promoted status; this wave
is its Resolution-path Move 1) / RQ1 (the branch that reopens on
inversion) / RQ3 (sequential shed, the speed-preservation caveat in
Step 2); `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4h (F.2b fate map, W10-P5
split verdict, boundaries 1–3); Addendum F (F.2b, F.3b — machinery
re-used verbatim); `drag_migration_log_tier2.md` 2026-07-11 (post-Wave-10
program order). Results → `TIER2_STAIRCASE_PROBE_FINDINGS.md` (new §4i +
insight register) + `drag_migration_log_tier2.md`; the user adjudicates.

---

# Addendum H (2026-07-11) — two-channel re-targeting: solvated-branch arbitration + Wave 12 (position heterogeneity) + Wave 13 (drag saturation)

> **Status: DESIGNED, endorsed 2026-07-11** (the two-channel bare-I⁺
> adjudication — log entry same day). **All MD legs and the Wave-13 drag
> arm stay behind `[PROCEED TO IMPLEMENTATION]`.** Entry documents: the
> 2026-07-11 log adjudication (experimental (n, KE) input; solvated
> re-score pre-read), `RESEARCH_QUESTIONS.md` RQ8 (new) + RQ4
> (critical-path promotion) + RQ3 (un-anchoring), findings §4i (Wave-11
> measurement + I37–I40). Same stance as all waves: **reported, not
> auto-adjudicated**; nothing here discharges the F5 gate.
>
> **Amended 2026-07-11 (same day, user endorsement — the "bounded-lever
> playground" discussion while the RQ4 external calculation runs):**
> (1) a **zero-MD analytic feasibility pass** is inserted as the gating
> step for all MD legs (new H.2b) — the validated F.2b fate map extended
> on paper by the three physically-bounded levers; (2) a **Wave 12b
> depth-dependent dressing arm** `n_eject(depth)` is designed (new H.3b,
> new physics surface, gated); (3) the **F.5 graded-ladder floor variant**
> (transition-at-n=1; knob-free, both depths [IHe05]-sourced) is admitted
> as an exercisable lever now, with the external many-body calculation
> re-cast as the *arbiter* of the transition location rather than a
> blocker (usage note appended to F.5). H.5 sequencing updated.
>
> **Second amendment 2026-07-11 (H.2b design freeze, user decisions
> D1–D7 — log entry same day):** the analytic pass is fully specified
> pre-trigger: sharp-E₀ scan over the RQ1 band (D1); the E₀–dressing
> coupling bracket family p ∈ {0, 1} + the gate-open-at-birth diagnostic
> map (D2); the exposure-law bracket extended to the histogram side
> (D3); the Wave-10/11 droplet-prior family (D4); the firm 3–6 Å margin
> band + lever-interaction map first (D5); the pre-registered
> reachability bar H2b-T1..T5 (D6); and the solvated mean-KE target
> table, superseding the 0.69 eV n = 1 anchor (D7). Details: the H.2b
> design-freeze block below.
>
> **H.2b EXECUTED 2026-07-11** (same day, under the trigger; scratchpad
> `h2b_feasibility.py`, zero MD, zero repo-code change; all wiring
> oracles pass — K 0.7446/0.74603 production, 0.89767/0.89754 at 9 Å,
> fate map exact). **Outcome (c) at the frozen bar:** no cell of the
> 7 308-point scan passes T1–T4 at either exposure bracket (ballistic is
> worse — not law-conditional). The miss is *structured*: bins 2–8 and
> the KE envelope land; the residuals decompose into (i) the **RQ4
> ladder taper** — the slid-X₂ family is rejected as the n₁ repair, the
> knob-free floor reaches n₁/n₂ = 1.83, and the RQ4-graded 2.2:1.5:1.3
> diagnostic nearly lands (W₁ 0.575 vs 0.5) → the external calculation
> is the blocking arbiter; and (ii) the **W13 v_c speed scale** (lift
> ×1.3–1.5, speed-selective). New structural finding: a **trapped
> droplet-retained ion class** (6–11 %, current-law only, inward
> partners below the 0.117 eV solvation barrier). W12b-P1's fast surface
> feeder is predicted absent in-bounds (margin floors n_eject ≥ 13).
> Full record: findings §4j + I41–I45; delivery entry:
> `drag_migration_log_tier2.md` "H.2b analytic feasibility pass
> EXECUTED". MD-leg sequencing goes to the user.

## H.1 Why — the two-channel reading re-targets the arbitration

Experiment: mean fragment KE falls ~3 eV (n = 0) → < 0.1 eV (n = 17),
and bare I⁺ is KE-distinct → **a different ionization channel produces
bare I⁺**. Consequences (log 2026-07-11): the biphasic mechanism owns the
**solvated branch only**; its structural refusal to make bare
(I8/I17/I22) is vindicated; the Wave-8→11 E₀ ≈ E\*(K) fine-tuning is
released. The solvated re-score shows the *new* gap: the flat-bottom
ladder feeds flat small-n bins (n₁ ≈ n₂ ≈ n₃ ≈ 13 %) where the
experiment shows a steep 31/14/9 decline — a ~2× n₁ deficit that no
(E₀, prior, drag-law) choice repairs. Two candidate mechanisms carry
that steepness, separated by the (n, KE) curve:

- **RQ4 ladder depth** (D₀ decreasing from n = 1, ratios ≈ 2.2:1.5:1.3)
  widens the n = 1 bin *energetically* — moves the histogram, not the
  speeds;
- **position heterogeneity** (Wave 12) reshapes the per-fragment
  exposure density p(K) *geometrically* — moves histogram **and** speeds
  together.

## H.2 The re-targeted arbitration surface

1. **Solvated histogram:** n ≥ 1 renormalized
   (`integrated_i_he_abundance.csv` with bin 0 dropped: n₁ = 31.0 %,
   n₂ = 14.2 %, n₂₋₄ = 30.5 %). Wasserstein per convention.
2. **Solvated (n, mean-KE) curve:** the user-supplied mean-KE table
   (H.2b D7, 2026-07-11 design freeze — **supersedes** the 0.69 eV
   n = 1 value quoted here previously): 2.9 eV at n = 0
   (channel-distinct, RQ8), 0.974 eV at n = 1 (the mass-131 VMI gate)
   falling to < 0.1 eV for n ≥ 13. **Prerequisite:** export the experimental
   table as a provenance-documented reference CSV under
   `data/reference/` (data-contract convention: MATLAB source of truth,
   exporter under `data/reference/scripts/`, columns + units + measurement
   IDs) before any calibration consumes it.
3. **The bare bin is a channel-branching quantity** (upper bound /
   mixed), interpreted only after RQ8; the model's own bare production
   is reported, not targeted.

## H.2b The analytic feasibility pass (Step 0.5 — zero MD; gates every MD leg; amendment 2026-07-11)

Before any MD or new physics surface: the F.2b fate map is MD-validated
at both kinematics (0.84 He at 9 Å; 0.24 He at production), and all
three bounded levers extend it **on paper**. Build the full forward
model — position × dressing × ladder × droplet prior → fate map →
solvated histogram + (n, KE) curve — and ask the pre-registered
reachability question:

> **Is (solvated n₁ = 31.0 %, n₂ = 14.2 %, n₂₋₄ = 30.5 % + the
> 0.69 → < 0.1 eV KE curve) reachable inside the all-bounded lever set?**

**The levers and their bounds (nothing floats):**

- **L1 — position geometry:** birth radius p(r) ∝ r² on [0, R − margin],
  margin ≈ one shell radius (4.67 Å); isotropic fragment axis → per-
  fragment chord through the sphere (partner chords anti-correlated).
  Bound: the volume-weight form is parameter-free; only the margin is a
  stated sensitivity (3–6 Å).
- **L2 — depth-dependent dressing:** n_eject(d) = round(21 · ρ̂(depth at
  birth)), **tied to the existing erf density profile** (zero new
  parameters; shell-averaged-ρ̂ variant reported as sensitivity, not
  fit). Σ(n_eject) shrinks accordingly (suppression easier for the
  surface class — the §4h-boundary-3 direction, now modeled instead of
  caveated).
- **L3 — F.5 graded-ladder floor:** D₀(1) = 13.3 meV (X₂/³Π), D₀(n ≥ 2)
  = the delivered mixture Form-U (flat 9.2 meV bottom) — the
  **transition-at-n=1 floor variant is knob-free** (both depths
  [IHe05]-sourced; the picture argument pins the first, least-frustrated
  atom unconditionally, and the external calculation can only move the
  transition *up*). Slid-transition variants (2, 3, …) are computed but
  **reported separately as the knobbed family** the RQ4 calculation will
  arbitrate; the ³Σ⁻ 7.9 meV tail is a stated sensitivity.
- **L4 — speed scale (v_c):** enters only the KE-curve read. The
  analytic pass scores the histogram exactly and **brackets** the KE
  curve between the current-law speed profile (stored trajectory) and
  the ballistic limit — v_c calibration stays W13's job.

**Approximation stance (F.3b precedent — brackets, MD arbitrates):**
off-center chords are integrated against the *pinned speed-vs-distance
profile* of the stored production trajectory (speed as a function of
path length, re-based per chord); the exposure per chord uses the same
ρ̂ and τ. Biases are direction-mixed and stated in the read-out; the
W12/W12b MD legs verify the analytically-landed region, not explore.

**Outcome shapes (each decisive for sequencing):**

- **(a) Reachable at the floor** (L3 transition-at-1, L1/L2 at stated
  bounds) → the MD legs verify a quantified prediction; the RQ4
  calculation arrives as arbitration with margin; flat-ladder burial is
  reversed only by that calculation.
- **(b) Reachable only with the transition slid up** → explicitly
  recorded as "waits on the RQ4 calculation"; the slid value is a
  prediction *for* that calculation, not a calibration.
- **(c) Unreachable inside all bounds** → fires the F.5 escape clause:
  the solvated steepness needs non-zero ε (reopen RQ2/RQ3) or a missing
  mechanism — the most consequential outcome, obtained at zero MD cost.

Route: scratchpad, zero new MD, zero repo-code change (the standing
Steps-1/1.5 convention), runnable under the trigger. Scoreboard CSV +
findings subsection; the user adjudicates which MD legs then run.

### H.2b design freeze (2026-07-11, user decisions D1–D7 — pre-trigger)

Settled in the pre-execution discussion; the pass is built exactly to
this specification. Execution still awaits `[PROCEED TO IMPLEMENTATION]`.

- **D1 — E₀ policy.** E₀ is a **sharp scalar scanned over the full RQ1
  band [0.2, 0.5] eV** as an outer dimension — not pinned at the
  Wave-11 landing (0.38–0.41 eV). With bare un-targeted (H.2 item 3)
  the scan is free to re-land wherever the solvated shape + KE brackets
  want it.
- **D2 — E₀–dressing coupling (the fifth lever bound).** The naive
  constant-E₀ convention was flagged in discussion (user: with L2 and
  nothing else changing, the surface class only feeds bare). Adopted:
  the parameter-free bracket family
  **E₀(d) = E₀·(Σ(n_eject(d))/Σ(21))^p with p ∈ {0, 1}**, both reported
  (p = 0 constant = the max-suppression bracket; p = 1 proportional
  occupancy — the same first-order tie as L2 itself). **Pre-derived
  structural insight (recorded):** under *both* brackets the surface
  class is suppressed — under-cooled either way (p = 0: Σ(n_eject)
  shrinks so K\* rises while the chord K falls; p = 1: K\* is uniform
  and short chords still undershoot it) — so the W12b-P1 "fast n = 1–3
  surface feeder" exists **only** for super-proportional E₀ decay
  (E₀ < Σ(n_eject) at birth, gate open from onset), which is unsourced.
  A **(depth, p) gate-open-at-birth diagnostic map** is computed and
  reported as a map, not a lever. Under p ∈ {0, 1} the n = 1–3 weight
  must come from the **near-cliff band of the chord-K distribution**
  (the r²-weighted density just above K\*); the **feeder map** (per
  detected bin, the birth-depth/chord-class decomposition) is a
  first-class scoreboard output — it is what separates the geometric
  feeder from the RQ4 ladder feeder on the (n, KE) curve (H.1).
- **D3 — the exposure-law bracket extends to the histogram side** (L4
  amended). K per chord is drag-law-conditional (Wave-11: 42 %
  beyond-band exposure, both identified biases pushing K up), so the
  **whole fate map runs at both brackets** — the current-law speed
  profile and the ballistic limit. Outcome (c) ("unreachable") fires
  only if it holds at **both** brackets; a single-bracket miss is
  law-conditional and is reported as such.
- **D4 — droplet priors.** The stated Wave-10/11 family: Kornilov
  log-normal δ = 0.625 about ⟨N⟩ = 2000 (primary), δ = 0.40 / 0.80
  sensitivity, pickup-weighted ∝ N^(2/3) variant. Positions per droplet
  on [0, R − margin] with the margin in absolute Å at every R.
- **D5 — the margin band 3–6 Å is firm.** The L1 margin caps L2's reach
  (ρ̂ at margin depth may floor n_eject near ~15–19 — erf-width
  dependent): the **lever-interaction map** (reachable n_eject and
  chord-K ranges vs margin, from the code's actual density profile) is
  computed *first* and heads the scoreboard. L2 turning out near-inert
  inside the firm band is a **finding** (the surface class is small by
  construction; the steepness must come from chord-K geometry and/or
  the ladder), never grounds to loosen the bound.
- **D6 — pre-registered reachability bar (H2b-T1..T5).** "Reachable" =
  at least one lever point inside all stated bounds satisfying T1–T4
  simultaneously; T5 is a consistency check, not a target. Outcome
  mapping unchanged ((a) at the L3 floor variant / (b) only with the
  slid transition / (c) nowhere in bounds).
  - **T1 (shape, primary):** solvated (n ≥ 1 renormalized) W₁ ≤ 0.5
    bins — half a rung, the quality scale of the best prior landing
    (Wave-11 untuned 0.496); the pass is closed-form, so no N = 50
    quantization allowance applies.
  - **T2 (steepness — the failing observable):** n₁/n₂ ∈ [1.75, 2.6]
    (experimental 2.18 ± 20 %; placeholder tolerance, no error bars
    supplied).
  - **T3 (top-bin absolute):** n₁ ∈ [26, 36] % (31 ± 5 abs) — blocks a
    T1 pass by mid-tail smearing.
  - **T4 (KE feasibility, bracketed):** every experimental mean-KE
    point n = 1–12 lies inside the [current-law, ballistic] envelope at
    the landed point, **and** a single speed scale inside the envelope
    reproduces the curve to ≤ ×1.5 per bin — deliberately a preview of
    W13's v_c: T4 passing means W13 has a solution to find; T4 failing
    at every lever point means no velocity cap reconciles the geometry
    (a finding).
  - **T5 (bare consistency, reported):** biphasic bare fraction
    ≤ 43.5 % hard (it stacks under the experimental bare bin together
    with the RQ8 channel); flagged above ~15 % as RQ8 tension (a
    visible slow/fast shoulder in the bare-KE distribution would be
    predicted).
- **D7 — the solvated mean-KE target table** (user-supplied 2026-07-11;
  **mean** kinetic energies; provenance-pending — the H.2 item-2 export
  prerequisite stands for arbitration/calibration; H.2b scores against
  the table as given). **Supersedes the 0.69 eV n = 1 anchor** quoted
  in H.2 and the 2026-07-11 two-channel log entry. n ≥ 13 is a
  one-sided < 0.1 eV band.

  | n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13–21 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | mean KE [eV] | 2.9 | 0.974 | 0.545 | 0.390 | 0.341 | 0.289 | 0.248 | 0.201 | 0.176 | 0.154 | 0.138 | 0.122 | 0.105 | < 0.1 |

  In complex-speed terms (m = 127 + 4n amu): n = 1 ≈ 12.0 Å/ps,
  n = 2 ≈ 8.8, n = 3 ≈ 7.4, n = 12 ≈ 3.4 — a ≈ 3.5× solvated speed
  span. Two immediate reads: (i) the Wave-11 center-birth companion
  (n_detect ≈ 6 at 5.31 Å/ps) already sits within ~6 % of the
  experimental n = 6 speed (5.63 Å/ps) — at face value the per-bin
  curve demands a **fast small-n class** more than a uniform speed
  rescale, exactly the H.1 geometric-feeder question and readable from
  the T4 single-scale test; (ii) **RQ8 kinematic datapoint
  (recorded):** bare at 2.9 eV sits above the ≈ 2.71 eV per-fragment
  ballistic ceiling of the 2.70 eV channel (at face value; no error
  bars) — bare is channel-distinct *kinematically*, independent of any
  drag argument.

## H.3 Wave 12 — position-heterogeneity probe (config flip; MD leg gated)

**Physics.** Sample the molecule's radial offset at ionization from the
volume-weighted interior distribution (4πr² weight, flat confinement to
leading order) with a surface-exclusion margin of order one shell radius
— the physically defensible dopant-delocalization ansatz. Per-fragment
paths then span ~(R − d) to (R + d): near-surface births give K ≈ 0
(under-cooled, fast, ≈ full energy) and their partners cross up to 2R
(over-cooled, slow) — the **partner anti-correlation** is the
mechanism's fingerprint, and the per-fragment K density near K\* is what
can steepen the small-n bins geometrically.

**Surface.** `single_initial_position = False` — the field, the sampler
(`sample_radial_positions`), and the per-atom depth plumbing (drag /
cooling / pickup gates all read per-atom depth vs per-atom radius)
already exist; **RNG draw order is safe by construction** (r0 is always
sampled and only zeroed under the center pin — MATLAB-parity decision).
**Audit before running:** the legacy sampler's distribution shape
(uniform-in-volume?) and surface-margin handling; a margin knob, if
needed, is new config surface and named in the trigger request.

**Run matrix (2 dirs + the Wave-11 A/B baseline):** production
kinematics (2.666 Å) with the Wave-11 pins; E₀ ∈ {0.28, 0.38} eV
(bracketing the re-opened scalar), position sampling ON; the Wave-11
center-pinned dirs are the controlled A/B partners. N = 50 bridge seed.

**Read-out:** per-fragment K distribution (shape of p(K), density near
K\*(E₀)); the (n_detect, v_detect) joint scatter; solvated-histogram
shift vs the center-pinned A/B; detected speed range vs the compressed
~1.5×; per-pair (K_i, K_partner) anti-correlation.

**Pre-registered predictions:**
- **W12-P1.** Per-fragment K spread becomes O(0.1–1) (vs the 10⁻¹³
  congruence), p(K) shaped by the volume weight (mass piled toward
  short exposures).
- **W12-P2.** The detected speed range stretches to ≥ 2× (fast tail
  from near-surface births retaining ~1.2–1.7 eV under the *current*
  law — geometry's share of the speed anomaly, measured).
- **W12-P3.** The solvated small-n bins steepen relative to the A/B
  (direction claim; magnitude is the wave's deliverable).
- **W12-P4.** (K_i, K_partner) anti-correlates within pairs.

**Boundaries:** single droplet radius (position axis deliberately
isolated from the size axis); current drag law (over-drag stands — W12
measures the *geometric* share of the anomalies); n_eject = 21 is now
position-strained (near-surface ions physically dress less → suppression
easier → reinforces the fast/low-n class; convention boundary, direction
favorable — **retired by W12b**, which models it); N = 50, single seed.

> **NB (2026-07-11, post-H.2b — predictions amended to the forward-model
> numbers; H.5 amendment item 2).** W12's role is now *verification* of
> the H.2b geometric inputs at the delivered physics (full dressing,
> flat ladder, current law, N = 2000, position on). Forward-model
> predictions (20 000 chords; margin band {0, 4.67} Å pending the
> sampler audit; 150 ps read):
>
> - **W12-P5 (trapped class — the load-bearing one):** 5.9–8.1 % of
>   fragments never eject (inward partners dissipated below the
>   0.117 eV solvation well); the relaxation stage decides whether the
>   park is permanent.
> - **W12-P1 (amended, quantified):** per-fragment K quantiles
>   5/50/95 % = 0.15/0.49/15.2 (margin 4.67) vs the 10⁻¹³ center-pin
>   congruence — the q95 is the trapped tail.
> - **W12-P2 (amended):** detected speed range (5–95 %) = 1.8–8.8 Å/ps
>   (margin 4.67; 1.9–10.5 at margin 0) vs the center-pin 4.13 — a ~5×
>   span, well beyond the original ≥ 2× claim.
> - **W12-P4 (amended):** strong *speed* anti-correlation within pairs
>   (corr(v_A, v_B) ≈ −0.85…−0.92); the K anti-correlation is weak
>   (≈ −0.12, diluted by the trapped tail) — the fingerprint lives in
>   the speeds, not the exposures.
> - **Fate reads (flat ladder, so no n₁ steepening expected — W12-P3 is
>   retired by §4j):** E₀ = 0.28 → suppressed 43–55 %, n̄_det ≈ 3.0–3.8,
>   near-flat small-n bins ≈ 0.05; E₀ = 0.38 → suppressed 70–77 %,
>   n̄_det ≈ 1.4–1.8. (The margin-4.67/E₀ = 0.28 suppressed fraction
>   0.435 numerically coincides with the experimental bare bin — a
>   curiosity, not a target: bare is channel-branching under RQ8.)

## H.3b Wave 12b — depth-dependent dressing `n_eject(depth)` (new physics surface, gated; amendment 2026-07-11)

**Physics.** An ion born near the surface cannot assemble a full
21-atom shell — the He is locally not there. Replace the t = 0
full-dressing convention with the **tied** law

$$n_{\rm eject}(d) \;=\; \operatorname{round}\!\big(21\cdot\hat\rho(d_{\rm birth})\big),$$

ρ̂ the *existing* erf density gate at the ion's birth depth (zero new
free parameters; the shell-averaged-ρ̂ variant is a stated sensitivity,
never a fit). Consequences: the surface class is capped at small n **by
construction** and exits fast (short path, little dissipation) — a
direct, *fast* n = 1–3 feeder, which is the fingerprint separating it
from the RQ4/L3 ladder feeder (slow, cascade-side) on the (n, KE) curve;
Σ(n_eject) shrinks per-ion (suppression easier for the surface class);
the maximally-strained §4h-boundary-3 convention is retired rather than
caveated.

**Implementation surface (behind its own trigger):** per-ion n(0) and
mass at the biphasic seed (`ion_initial_state.py` S2 block: n_shell(0),
`mass_kg` per atom — arrays are already per-atom, **no checkpoint schema
change**), the E_pot binding fold seeded with e_bind_pair(n_eject(d))
per ion (the 5-term invariant closes per-ion as before), behind an
interchangeable enum (`full` = delivered 21-for-all vs `density_tied`).
Dimensional analysis trivial (counts; depths in Å against the existing
gate). Wiring oracle: `full` arm reproduces a delivered dir
byte-identically; center-pinned + `density_tied` reproduces n_eject = 21
exactly (ρ̂(center) ≈ 1) — the arm is inert without the position axis.

**Run matrix:** 2 dirs = the W12 matrix re-run with `density_tied` —
the A/B/C chain (Wave-11 center-pinned / W12 position-only / W12b
position+dressing) isolates each mechanism's share.

**Pre-registered predictions:**
- **W12b-P1.** The n = 1–3 bins gain weight *beyond* W12's geometric
  share, fed by near-surface births, and those fragments are the
  **fastest** in their bins (the anti-RQ4 fingerprint).
- **W12b-P2.** Model-side bare rises from the surface class (n_eject
  small + K small ⇒ easy total shed) — reported against the RQ8
  branching bound, not targeted.
- **W12b-P3.** Deep-interior physics is untouched (mid/high-n bins and
  the deep-shell tail shift < 1 bin vs W12).

**Boundaries:** the tied law is a first-order occupancy statement (no
shell-restructuring dynamics); the dressing–pickup interplay
(under-dressed slow ions re-filling via the live Langmuir channel) is
read from the runs, not separately modeled; single radius, N = 50,
current drag law.

> **NB (2026-07-11, post-H.2b): W12b is PARKED** (H.5 amendment item 3).
> The H.2b lever-interaction map shows the W12b-P1 fingerprint is
> predicted **absent** in-bounds: at the code's 14.2 Å erf width the L1
> margin floors n_eject at 13/14/15 (margin 3/4.67/6 Å), and no
> gate-open-at-birth class exists under either bounded E₀ law
> (findings §4j, I45). Building this surface would confirm an inertness
> prediction. Revive only if the W12 verification falsifies the
> forward model's geometry (e.g. the trapped class or the chord-K
> distribution comes out materially different).

## H.4 Wave 13 — drag saturation (the v_c arm; repo-code change, gated)

**Form (strict dimensional analysis):** with b [amu·ps/Å²] the locked
pure-cubic coefficient (2.5154) and v_c [Å/ps] the single new knob,

γ(v) = b·v² for v ≤ v_c; γ(v) = b·v_c·v for v > v_c   [amu/ps ✓]
|F| = b·v³ → b·v_c·v² (cubic below, quadratic — the physical
Newton-drag asymptotic — above; continuous at v_c) [amu·Å/ps² ✓].

**Physical motivation:** cubic force extrapolated to 10.5 Å/ps is ~8×
the classical form-drag estimate (½·ρ_He·C_d·A_shell·v²); Wave-11
measured 41.7 % of the production exposure beyond the calibrated band
and a ~2× detected-speed undershoot on the correctly-conditioned n = 1
gate (0.18 vs 0.69 eV). TDDFT cannot arbitrate at production kinematics
(it breaks at 2.666 Å — the program's premise); **the n = 1-gated VMI
peak is the experiment's own measurement of the high-v drag integral**
and becomes the calibration target. Scope decision (endorsed): the VMI
*peak position* moves from validation to calibration; the speed
*distribution shape* and the whole solvated histogram remain validation.

**Implementation shape:** a new interchangeable arm behind the drag-form
enum in `physics/drag.py` (naming per DRAG_PORT_DESIGN_DECISIONS), with
the **in-band invariance oracle**: v_c ≥ 5.3 Å/ps (above the measured
0.80 eV in-window max 5.23) ⇒ every delivered 0.80 eV artifact
byte-identical; v_c = ∞ reproduces a delivered probe dir exactly
(regression test). v_c classification: **Bounded** (5.3 ≲ v_c ≲ 15),
→ Derived once pinned.

**Calibration loop (pre-registered to converge ≤ 2 iterations):**
(1) pin v_c on the n = 1-class detected speed under the droplet prior +
position axis (the W12-informed conditional, not the pinned-droplet
speed); (2) re-measure K₂.₇₀ under the calibrated arm (the Wave-11
exposure machinery verbatim — the 0.746 value is law-conditional);
(3) re-run the solvated-branch scan; (4) confirm the n = 1 class
location moved < 1 rung.

**Pre-registered prediction W13-P1:** v_c alone does **not** fix the
solvated small-n steepness (ladder-geometric, drag-law-robust) — v_c
owns the (n, KE) curve's *scale*, RQ4/position own its *shape*. If v_c
calibration *does* repair the histogram shape, the two-mechanism
decomposition of H.1 is wrong — a finding.

> **NB (2026-07-11, post-H.2b — mandate extended; H.5 amendment items
> 1+4).** (a) The H.2b pass pre-confirms W13-P1's premise with numbers:
> the experimental mean-KE curve sits inside the [current-law,
> ballistic] envelope at every top cell, and the required lift is
> speed-selective ×1.3–1.5 (I43) — a v_c solution exists to find.
> (b) The **(n, mean-KE) reference export is now a calibration
> prerequisite** (the D7 table, provenance-documented, before the
> calibration loop runs). (c) Step 3 of the calibration loop gains a
> concrete deliverable: **re-run the H.2b solvated scan under the
> calibrated arm** — the exposure K per chord is drag-law-conditional,
> H.2b scored only the two extremes, and the calibrated law lands
> between them; the re-scan fixes the size of the taper the RQ4
> calculation must deliver, and the RQ4 ladder is then arbitrated
> against the v_c-updated scan (H.5 amendment item 5).

## H.5 Sequencing, gates, outcome shapes

- **Order (amended 2026-07-11): H.2b first, and it gates everything.**
  The analytic pass decides *which* MD legs are worth running (outcome
  (c) cancels them in favour of the ε/mechanism question). Then
  W12 → W12b (the A/B/C chain needs the position-only baseline) → W13
  (its calibration conditional inherits the W12/W12b geometry); RQ8 +
  the (n, KE) export are entry gates for *arbitration* (the A/B/C
  read-outs are self-contained and can run first); RQ4's external
  many-body input proceeds in parallel and returns as the **arbiter of
  the L3 transition location** (H.2b outcome (b) turns it into a tested
  prediction).
- **Gating:** W12 MD leg (2 dirs + possible margin knob) and the entire
  W13 arm (new physics surface) each behind their own
  `[PROCEED TO IMPLEMENTATION]`. Nothing discharges F5; the F5 re-scope
  adjudication (log 2026-07-11, post-Wave-10 entry) now also inherits
  the solvated-branch re-targeting.
- **Outcome shapes:**
  - **(a)** Position supplies most of the speed range *and* small-n
    steepness → v_c is a mild correction; ladder may stay flat; RQ4
    demoted again.
  - **(b)** Position supplies range but not steepness → RQ4 depth is
    load-bearing; W13 calibrates scale; the ladder-bottom external
    calculation becomes the blocking item.
  - **(c)** Neither suffices → mechanism-level OQ on the solvated shape
    (reopen ε/RQ2, picture, or the E₀-smear) — the most consequential
    outcome.
  - **(d)** RQ8 finds the bare bin partially biphasic → a branching
    nuisance parameter enters the arbitration; the solvated targets
    stay primary.

**Cross-links:** log 2026-07-11 (two-channel adjudication);
`RESEARCH_QUESTIONS.md` RQ8/RQ4/RQ3/RQ7; findings §4i + I37–I40
(Wave-11 measurement); Addendum G (K machinery); `CALIBRATION_MAP.md`
propagation pending at W13 build. Results → findings (new §4j/§4k +
insight register) + the log; the user adjudicates.

### H.5 amendment (2026-07-11, post-H.2b sequencing — user endorsed)

The H.2b execution returned the *structured* outcome (c) (findings §4j:
residuals = the RQ4 taper + the v_c scale; ε contingent on the RQ4
calculation). The original "outcome (c) cancels the MD legs for the ε
question" clause assumed an unstructured miss and is **superseded**; the
endorsed sequence:

1. **The (n, mean-KE) reference export** (H.2 item 2) is promoted from
   arbitration entry gate to **calibration prerequisite** — W13 pins v_c
   on the n = 1-class detected KE, so the D7 table must become the
   provenance-documented CSV (convention, mass gates, error bars → they
   also firm the placeholder ±20 % tolerances) before the W13
   calibration loop runs.
2. **W12 runs, re-scoped from exploration to verification** of the
   forward model's load-bearing geometric inputs — foremost the trapped
   droplet-retained class (I44; the MD relaxation stage arbitrates
   whether the 150 ps "park" is permanent), plus the chord-K
   distribution, the detected speed range, and the pair
   anti-correlation. Pre-registered predictions amended to the H.2b
   numbers (H.3 NB below). MD leg stays behind its own
   `[PROCEED TO IMPLEMENTATION]`.
3. **W12b is PARKED** (H.3b note): its W12b-P1 fingerprint is predicted
   absent in-bounds (the margin floors n_eject ≥ 13; no gate-open-at-
   birth class under p ∈ {0, 1}) — building a new physics surface to
   confirm an inertness prediction is poor spend. Revive only if W12
   falsifies the pass's geometry.
4. **W13 follows W12** (calibration conditional unchanged) with an
   **extended mandate**: after calibrating v_c, re-run the H.2b solvated
   scan under the calibrated law (the pre-registered W13 step 3) — this
   also closes the one H.2b coverage gap (the histogram was scored only
   at the two drag extremes; the calibrated law lands between them) and
   fixes the size of the taper the RQ4 calculation must deliver.
5. **RQ4 external calculation proceeds in parallel and returns as the
   decisive arbiter:** its ladder goes verbatim into the v_c-updated
   scan. Steep bottom (≈ 2 : 1.4 or more) → the histogram closes inside
   physics; plateau → the F.5 escape clause (ε / RQ2, missing
   mechanism) fires with the geometric alternatives exhausted.

RQ8 remains the entry gate for any interpretation of the bare bin.
Nothing in this amendment discharges the F5 gate.

---

# Addendum I (2026-07-15) — Wave 13 re-scoped: drag-form reachability & calibration (W13-first reorder)

User decision 2026-07-15, following the delivery of the provenance-documented
(n, mean-KE) reference (`data/reference/ihe_ked/`, deployed 2026-07-14,
run-summary comparison layer integrated 2026-07-15). **Supersedes** (a) the
H.5-amendment ordering "W12 before W13" and (b) the single-arm scope of H.4:
Wave 13 is promoted to the immediate priority and broadened from a one-knob
v_c calibration to a drag-form *family* investigation against the corrected
experimental mean-KE curve. Design frozen here; every execution step stays
behind `[PROCEED TO IMPLEMENTATION]`.

## I.1 Why — the export invalidates the D7 premise numbers and re-ranks the observables

1. **The D7 table is stale, coherently.** The H.2b T4 verdict was scored
   against the user-supplied D7 means (n = 0: 2.9 eV, n = 1: 0.974 eV,
   < 0.1 eV for n ≥ 13). The corrected export
   (`IHe_KED_reference.csv`) sits ~25–35 % higher across the curve
   (n = 0: 3.706 eV, n = 1: 1.302 eV) — the D7 values match the *legacy
   moment convention* the export's changelog identifies as the central
   pipeline bug (2-D slice density, missing dE Jacobian). The shift is far
   outside the two correlated bands (4 % calib ⊕ 6 % condition), so the
   T4/I43 premise numbers ("envelope holds n = 1–12; required lift
   speed-selective ×1.3–1.5") are unverified against the real reference.
2. **The n = 1 class sits at the ballistic edge (indicative).** Two-body
   2.70 eV KER → ~1.35 eV on the I⁺ fragment; mass-corrected to the
   detected I⁺He complex ~1.31 eV; minus the ~0.117 eV solvation barrier
   ≈ 1.2 eV zero-drag expectation, vs the measured 1.302 eV mean. Whether
   *any* drag form can reach the small-n bins is genuinely open — the
   question must be answered before a calibration loop is built.
3. **Observable re-ranking.** The histogram arbitration is blocked on the
   external RQ4 calculation regardless; the (n, mean-KE) curve is now the
   best-constrained observable in the program (gold points at 1–4 %).
4. **W12's verification targets are drag-law-conditional** (trapped class,
   chord-K, detected speed range), so verifying them under a law about to
   be replaced is half-wasted; verification moves to *after* calibration.

## I.2 Decisions (user, 2026-07-15)

- **I-D1 (ordering):** W13 immediate; W12 demoted to post-calibration
  verification under the calibrated law; W12b stays parked; RQ4 external
  calculation stays parallel and returns against the v_c-updated scan
  (H.5-amendment items 3 and 5 otherwise unchanged).
- **I-D2 (sweep engine):** the extensive form × knob scan runs in the
  validated H.2b 1D chord forward model (oracle-matched to MD landmarks at
  ≤ 0.2 %); full MD confirms the shortlist only. MD remains the arbiter.
- **I-D3 (n = 0 excluded):** the bare bin is RQ8 channel-distinct
  (suspected different process, not yet simulable) and enters neither the
  calibration nor the acceptance bar.
- **I-D4 (calibration observable):** mean-to-mean ⟨E⟩(n), full-curve fit
  over n = 1…17 (upgrade from the single n = 1 pin of H.4), weighted by
  the per-point stat ⊕ sys errors with the 4 % calib and 6 % condition
  bands profiled as coherent scale nuisances. The n = 0…4 P(E) trusted
  curves are *validation-only* (mean-to-mean is the calibration
  convention per the export README; never mean-to-peak).
- **I-D5 (acceptance bar):** per-bin energy ratio within **×1.25** for
  n = 1…12 after profiling the coherent bands; n ≥ 13 need only remain in
  the sub-0.1 eV band. The bar is a **ceiling of strictness**: it may be
  loosened with documented model-bias justification (pinned mass −2.2 %,
  no-shed leak, bin reshuffling), never tightened.
- **I-D6 (form family):** all candidate forms reduce **exactly** to the
  locked `b·v²` in-band (Tier-0 lock and the H.4 in-band invariance
  oracle preserved); the investigation space is the high-v tail only.

## I.3 The tail family (strict dimensional analysis)

One two-knob family covers the physical candidates. With
b = 2.5153509 [amu·ps/Å²] the locked pure-cubic coefficient,
v_c [Å/ps] the cap, and p [dimensionless] the tail exponent:

γ(v) = b·v²                    for v ≤ v_c        [amu/ps ✓]
γ(v) = b·v_c²·(v/v_c)^p        for v > v_c        [amu/ps ✓]

continuous at v_c for every finite p (both branches → b·v_c²). Force
|F| = γ·v = b·v_c²·v·(v/v_c)^p [amu·Å/ps² ✓]:

| p | γ tail | F tail | physical reading |
|---|---|---|---|
| 1 | b·v_c·v | ∝ v² | Newton form-drag asymptotic (= the H.4 arm) |
| 0 | b·v_c² | ∝ v | Stokes-like linear force |
| −1 | b·v_c³/v | constant | saturated (plastic) drag force |
| −∞ (limit) | 0 above v_c | 0 | hard cutoff — extreme no-high-v-drag bracket; γ discontinuous at v_c, admissible as a 1D diagnostic bracket only (needs smoothing if ever promoted to MD) |

The sweep scans p ∈ {1, 0, −1} × v_c (Bounded, 5.3 ≲ v_c ≲ 15 per H.4,
the ≥ 5.3 floor protecting the 0.80 eV in-window byte-identity), plus the
p → −∞ bracket. Continuous-p refinement is admitted if the discrete set
brackets a landing. The drag module stays mass-agnostic (γ in amu/ps; m
enters only at the integrator O-step) per the friction convention.

## I.4 Steps

- **Step 0 — envelope re-read (zero MD, zero repo code).** Score the two
  existing brackets — current pure-cubic law and ballistic — against
  `IHe_KED_reference.csv` under the I-D4 error treatment. Deliverables:
  per-n reachability verdict (inside/outside the [current-law, ballistic]
  band after coherent-band profiling), the required-lift curve, formal
  retirement of the D7 table, and correction of the stale T4/I43 premise
  numbers in the findings doc. Caveat carried explicitly: the envelope is
  a sharp guide, not a per-bin theorem — changing the drag law also
  reshuffles fragments between n bins. Step 0's numbers generate the
  quantitative pre-registered predictions for Step 1 *before* Step 1 runs.
- **Step 1 — tail-family sweep (1D chord model, scratchpad).** §I.3
  family × the H.2b geometry (droplet prior + position axis). Score per
  I-D4/I-D5 over n = 1…17. Every cell also records its abundance-histogram
  W₁ as a **report column** (not a pass bar — the RQ4 taper is pending),
  so the KE fix's histogram side-effect is visible. Output: shortlist of
  1–3 (p, v_c) candidates + the discrimination read (does the curve shape
  separate the tails, or only the scale?).
- **Step 2 — MD confirmation (repo change, own trigger).** Implement only
  the shortlisted arm(s) behind the drag-form enum in `physics/drag.py`
  (naming per DRAG_PORT_DESIGN_DECISIONS; v_c = ∞ byte-identity
  regression; in-band invariance oracle). A small set of 2.70 eV
  production runs under the droplet prior with the biphasic mass
  mechanism + detection stage, scored through the delivered `ihe_ked`
  run-summary comparison layer (mean-KE panel primary; n = 0…4 curve
  overlays as validation).
- **Step 3 — retained from H.4/H.5 (unchanged in substance).** Re-measure
  K₂.₇₀ under the calibrated arm (Wave-11 machinery verbatim); re-run the
  solvated-branch scan under the calibrated law (fixes the taper size RQ4
  must deliver); confirm the n = 1-class location moved < 1 rung. Then
  the re-scoped W12 verification runs under the calibrated law, and RQ4
  arbitrates against the v_c-updated scan.

## I.5 Pre-registered predictions

- **W13-P1 carries over unchanged:** the drag form owns the (n, KE)
  curve's scale, not the solvated histogram's small-n steepness
  (ladder-geometric). If a tail choice *does* repair the histogram shape,
  the H.1 two-mechanism decomposition is wrong — a finding.
- **I-P1:** Step 0 finds the n = 1 bin at or inside the ballistic bracket
  only after coherent-band profiling (i.e. the n = 1 class is essentially
  undecelerated; the required high-v drag is near zero). If n = 1 exceeds
  the ballistic bracket *beyond* the coherent bands, outcome (c′) fires.
- Numeric per-bin predictions for Step 1 are generated by Step 0 and
  registered in the findings doc before Step 1 executes (program
  convention).

## I.6 Outcome shapes

- **(a)** Some (p, v_c) lands n = 1…12 within the I-D5 bar → calibrated
  law fixed; proceed Step 2 → Step 3.
- **(b)** Scale lands but the curve shape fails beyond the bar for every
  tail → drag alone is insufficient; the position/E₀-smear/RQ4 couplings
  re-enter (maps onto H.5 outcomes (b)/(c)).
- **(c′)** Experimental means exceed the ballistic bracket beyond the
  coherent bands at some n ≤ 12 → **no drag form suffices**; W13 halts
  before any build and a mechanism-level OQ opens (kick-energy share,
  detection convention, added-energy channel).
- **(d)** Degenerate landing (multiple tails pass) → v_c pinned per form;
  discrimination deferred to P(E)-shape validation and Tier 3 second
  moments.

## I.7 Boundaries and gating

1. Steps 0–1 are scratchpad-only analytics (H.2b precedent) but still
   execute only under `[PROCEED TO IMPLEMENTATION]`; Step 2 is a repo
   physics-surface change behind its **own** trigger.
2. The 1D model's boundaries carry verbatim from §4j (pinned-mass
   straight chords, no-shed leak baked into the closed form, detected ≈
   energetic floor at s_eff = 8, ε = 0, 150 ps trapped-class read).
3. Nothing here discharges the F5 gate; RQ8 remains the entry gate for
   any bare-bin interpretation.
4. Tier-0 lock untouched: every candidate form is byte-identical physics
   in the TDDFT-calibrated band.

**Cross-links:** findings §4j (H.2b verdict + boundaries); H.4 (form
motivation, in-band oracle, calibration-loop skeleton); H.5 amendment
(items 3/5 retained); `data/reference/ihe_ked/README.md` (error model,
mean-to-mean convention, trust tiers); `CALIBRATION_MAP.md` propagation
pending at the Step-2 build; log entry 2026-07-15 (this endorsement).

## I.8 Step 1b — joint (v_c, τ) mini-sweep (OQ-I arm (a); user-approved 2026-07-15, post-Step-1)

Step 1 returned outcome (b) with the K-coupling finding (I47): the drag
tail and the evaporative descent ride the same exposure integral, so the
KE lift starves the histogram. Arm (a) of OQ-I tests whether a faster
cooling clock decouples them. **Zero MD, scratchpad; same trigger
regime as Steps 0–1.**

**Mechanics.** τ enters the chord dynamics nowhere — only the
bookkeeping K = ∫ρ̂ dt/τ. So K(τ) = K(6.55) · (6.55/τ) per chord, and
the τ axis is a pure fate-map rescale of cached tail trajectories:
13 integrations (p ∈ {1, 0, −1} × v_c ∈ {5.3, 6, 7, 8.5} Å/ps + one
v_c = 15 current-law control), then the scan over
τ ∈ {6.55, 5.2, 4.1, 3.3, 2.6, 2.0} ps × the Step-1 fate-cell subset is
free. Scoring unchanged (I-D4/I-D5 KE bar; histogram W₁/n₁/ratio as
report columns; joint-closure read-out at W₁ ≤ 0.85 — the floor1
current-law level — with the T2 ratio still RQ4-reserved).

**Physicality note (recorded before execution):** τ = 6.55 ps is the
GAH25-sourced probe pin (Bounded). A landing at τ ≪ 6.55 is a
*re-classification event* for the calibration map (the pin becomes a
fitted knob and its GAH25 sourcing must be re-argued or dropped), not a
free pass — the sweep measures *what τ closes*; whether that τ is
physical is a separate adjudication.

**Pre-registered predictions (2026-07-15, before execution):**

- **J-P1 (matched clock):** the histogram recovers where τ restores the
  current-law K scale: τ\*(p, v_c) ≈ 6.55 · K_q50(tail)/K_q50(current)
  — numerically τ\* ≈ 4.1 ps at (p = 0, v_c = 6), ≈ 3.4 ps at
  (p = −1, v_c = 5.3), ≈ 5.3–6.3 ps for the mildest tails. E₀ optima
  return into the interior (0.21–0.25 eV) off the 0.20 edge, and
  W₁ recovers to the §4j current-law level (0.58–0.85) at the matched τ.
- **J-P2 (the decoupling test proper):** at the matched τ the binning
  returns near-current-law while the per-chord speeds keep the tail
  lift, so the KE curve becomes a near-uniform lift of the current-law
  curve — deep bins improve markedly (from ×0.49 toward ≳ 0.7 at
  n = 12). Whether they clear the strict ×1.25 bar is the open
  question: the tail lifts high-K chords *least* (they spend longest
  below v_c where the law is locked), so residual deep-bin failure at
  the bar remains the likelier outcome. The diagnostic is the worst
  deep-bin ratio, not the pass count.
- **J-P3 (control):** the v_c = 15 current-law control shows no KE
  movement on the τ axis (KE is τ-independent by construction); only
  its binning shifts. Any apparent KE-vs-τ trend at v_c = 15 beyond
  binning reshuffles would flag a wiring error.

**Outcome shapes:** (a) joint closure (KE bar + W₁ ≤ 0.85) in a
contiguous (v_c, τ) region → the W13 Step-2 MD build re-opens with a
two-knob arm + the τ re-classification question; (b) histogram recovers
but deep bins stay < ×1.25⁻¹ → the slope is *not* K-coupling — it moves
to OQ-I arms (b)/(c) (fate map / in-band), drag calibration closes as
"scale-only"; (c) histogram does not recover at any τ ≥ 2.0 → the
K-rescale picture itself is wrong (binning is not K-monotone under
tails) — a model-structure finding.

## I.9 Step 1c — closure-basin refinement scan (user-approved 2026-07-15, post-Step-1b; RQ9 parked in parallel)

Step 1b landed on a coarse grid (τ quantized at 4.1; v_c at 7 points;
deep-bin passes hugging the bar edge). Before any Step-2 build, map the
joint-closure basin finely — is it a plateau or a knife-edge, and what
are the calibration targets with uncertainties? **Zero MD, scratchpad,
same trigger regime.** The τ-sourcing question is **RQ9 (parked, user
decision)** — this scan supplies its model-side sensitivity input.

**Grid.** p ∈ {0, −1}: v_c ∈ {5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5} Å/ps
(4 new integrations per form; 3 cached); p = 1 control at the cached
{5.3, 6.0, 7.0, 8.5}; τ ∈ 3.0–5.6 ps step 0.2 (free rescale);
E₀ ∈ 0.20–0.32 step 0.01 (the closure lives at 0.22–0.26); ladders
rq4graded + floor1 (control); both priors; margins {3, 4.67, 6} Å.
Scoring identical to §I.8 (I-D4/I-D5 + W₁ ≤ 0.85 joint criterion).
Deliverables: basin occupancy map over (p, v_c, τ); per-form best-fit
(v_c, τ, E₀) with basin-width uncertainties; sharpened p = 0 vs p = −1
discrimination read (n ≥ 13 tail).

**Pre-registered predictions (2026-07-15, before execution):**

- **K-P1 (plateau):** the joint basin is contiguous per form, centered
  near τ ≈ 4.1 ps, with width ≥ 0.4 ps in τ and ≥ 1.0 Å/ps in v_c —
  a calibratable plateau, not a knife-edge (based on the Step-1b
  bar-pass spread τ ∈ [3.3, 5.2] and closure spread v_c ∈ [6, 8.5]).
- **K-P2 (form exclusion):** p = 1 (Newton tail) remains excluded — no
  joint closure at any (v_c, τ) (its Step-1b W₁ floor among KE-passing
  cells was 1.11, far above 0.85).
- **K-P3 (taper still required):** floor1 yields zero joint closures
  anywhere on the fine grid — the rq4graded conditionality is a ladder
  property, not a grid artifact.

## I.10 Step 2 — MD confirmation build (design frozen 2026-07-16, user decisions S2-D1–S2-D4; execution behind its own `[PROCEED TO IMPLEMENTATION]`)

Surface audit (2026-07-16): τ, f_int, s_eff, cooling gate, biphasic
mechanism, relaxation + detection stages, off-center birth sampling
(`single_initial_position=False`) and the ihe_ked scoring layer all
exist. Missing exactly two repo surfaces: the capped drag tail and a
config path to the existing-but-unwired `TabulatedLadder`.

**Decisions (user, 2026-07-16):**

- **S2-D1:** build the tabulated-ladder config surface (wires the
  existing `tabulated_ladder(rungs_eV)` path; later receives the RQ4
  answer verbatim).
- **S2-D2:** statistics = N = 50 pilot on all configs → N = 500 on the
  winning cell.
- **S2-D3:** focused matrix — 2 targets + 2 controls (4 configs).
- **S2-D4:** fixed N = 2000 droplets (probe convention, §6.6 defense);
  the droplet-distribution axis is deferred.

### Slice T1 — drag arm `capped_cubic` (physics/drag.py)

Strict dimensional analysis. Coefficients
{b [amu·ps/Å²], v_c [Å/ps], p_tail [dimensionless]}:

γ(v) = g(depth)·b·v²                    v ≤ v_c    [amu/ps ✓]
γ(v) = g(depth)·b·v_c²·(v/v_c)^p_tail   v > v_c    [amu/ps ✓]
|F| = γ·v [amu·Å/ps² ✓]; continuous at v_c for all finite p_tail.

Contracts preserved: mass-agnostic (no m argument); gate single-sourced
via `spatial_gate`; γ = F/v identity; FDT amplitude uses γ directly
(Tier-3 stays inert); dissipativity guard (b > 0, v_c > 0); form/coeff
cross-check extended to the new form; p_tail restricted to {0, −1} at
config load (the Step-1c-surviving set; others need a new adjudication).

Oracles/tests: closed-form pins per branch; **v_c = ∞ (or v_c ≥ v_max)
byte-identity with `linear_cubic` (a = 0, same b)**; in-band invariance
(a delivered 0.80 eV probe dir reproduces byte-identically at
v_c ≥ 5.3 — in-window max speed 5.23); existing suites green
(`test_drag.py`, `test_drag_config.py`, `test_ion_drag_smoke.py`,
`test_tier0_drag_comparison.py` — the Tier-0 18 Å gate is untouched
because the in-band law is identical).

> **NB (2026-07-16): Slice T1 DELIVERED** under its own trigger (log
> entry "Slice T1 DELIVERED"). All contracts as specified; both oracles
> pass — the delivered bridge probe dir reproduces **byte-identically**
> at v_c = 5.3 for both tails (measured in-window max speed
> 5.2306 Å/ps), and the v_c = ∞ / v_c ≥ v_max identity is exact (`==`)
> because the in-band branch is literally the pure-cubic arithmetic.
> Full suite 2250 passed. CALIBRATION_MAP rows 4 (v_c realizes the R10
> v_ceiling, Bounded), 4b (p_tail Free-choice {0, −1}), and 11 (τ RQ9
> pending) propagated. No preset/production config selects the arm yet
> (Slice T3's job); T2–T4 await their builds.

### Slice T2 — ladder config surface (config.py + dissociation_ladder)

New optional field `tabulated_ladder_rungs_eV` (active only with
`dissociation_ladder="tabulated"`; validation: ≥ n\* positive entries,
loud failure otherwise; inert default None). The rq4graded / floor1 rung
tables are *constructed in the generator script* (probe-side, from
`d0_of_n` × the taper multipliers) — no taper physics enters the
package. Tests: tabulated == form_u when fed the form_u rungs;
guard behavior; inert-default regression.

> **NB (2026-07-16): Slice T2 DELIVERED** under its own trigger (log
> entry "Slice T2 DELIVERED"; TDD, tests watched fail first). Scope ran
> through **all three stages** as approved pre-trigger:
> `resolve_ladder(dissociation_ladder, rungs_eV)` in
> `physics/dissociation_ladder.py` is the single cfg→ladder bridge
> (pairing + ≥ n\* positive-finite-rung validation; called by
> `check_ladder_config` at config load and by ion `biphasic_step` /
> relaxation / detection at point-of-use), and every ladder consumer
> gained a byte-inert `ladder=None` kwarg — the injected table replaces
> the Form-U parametrisation entirely (picture/κ ignored when set). The
> two pre-T2 `NotImplementedError` stage refusals became the loud
> missing-table `ValueError`; the relaxation stage — previously
> guardless, silently Form-U — now resolves. `cfg.json` round-trips the
> table (list→tuple in `load_cfg`); a pre-T2 `cfg.json` loads with the
> inert `None` default (back-compat test, the Slice-DS precedent).
> Equivalence oracle: tabulated fed the Form-U rungs is **bit-identical**
> through `biphasic_step`, relaxation, and detection; liveness proven by
> distinct-table tests + the short-table loud pickup-lookup failure.
> Range convention: config floor ≥ n\* = 21 entries; generator tables at
> N_STAR + 11 = 32 (the Form-U cache height); out-of-table occupancy
> fails loudly at the lookup. Full suite 2293 passed. No preset selects
> `tabulated` (T3's job); T3–T4 await their builds.

### Slice T3 — generator + pilot runs (N = 50, production kinematics)

Common config: biphasic; R0_GS = 2.666 Å, E_coulomb_scale = 1.0
(2.70 eV/fragment channel); E_int(0) via
`internal_energy_partition_fraction` = E₀ / 2.70; s_eff = 8;
`cooling_spatial_gate="density_scaled"`; `single_initial_position=False`
at fixed N = 2000 (margin 0 — 1D used 3–6 Å; recorded caveat);
relaxation + detection stages on (8.53 µs); fixed seeds; probe dir
namespace.

| config | drag | τ [ps] | ladder | E₀ [eV] | role |
|---|---|---|---|---|---|
| C1 | p = −1, v_c = 7.5 | 3.8 | rq4graded | 0.25 | full-house target |
| C2 | p = 0, v_c = 6.5 | 4.0 | rq4graded | 0.24 | form discrimination |
| C3 | linear_cubic (current) | 6.55 | flat (form_u) | 0.25 | baseline control |
| C4 | p = −1, v_c = 7.5 | 4.4 | floor1 | 0.23 | bounded-physics claim (I51) |

> **NB (2026-07-16): Slice T3 DELIVERED + pilots EXECUTED** under its
> own trigger (log entry "Slice T3 DELIVERED"; TDD, RED watched first;
> full suite 2328 passed). As-built: `build_biphasic_cfg` gained eight
> None-sentinel byte-inert kwargs (drag form-swap inheriting the locked
> bundle `b` structurally; ladder selector + table; `R0_GS_angstrom`;
> `E_coulomb_scale` stamped explicitly; `single_initial_position`;
> `detection_time_ps`); the C-runs are tagged by config label
> (`tier2probe_conf270_c1`…, probe namespace, knobs read from
> `cfg.json` — a knob-encoding tag would alias C1/C2 at two-decimal
> f_int); `scripts/gen_tier2_md_confirmation.py` builds the rq4graded /
> floor1 tables generator-side (Σ(21) Form-U pin 0.18783720 eV
> asserted) and runs neutral → ion → E2 → detection inline. **All four
> N = 50 pilots are on disk** (first in-repo off-9 Å MD; relaxation
> full 1000 ps cap): n̄_detect = 8.63 / 8.77 / 6.97 / 8.39
> (C1/C2/C3/C4). **dt-halving spot check clean** at production
> kinematics (C1, dt 0.01 → 0.005: n̄_detect drift 0.30× SEM, detected
> KE drift 0.14× SEM) — the Tier-0 dt carries. Scoring + S2-P1–P4
> verdicts await Slice T4 (note there: the staircase report's probe
> glob sweeps conf dirs — exclude or ignore).

### Slice T4 — scoring, winner selection, N = 500 confirmation

Pilot scored through the ihe_ked run-summary layer (mean-to-mean,
correlated bands per I-D4) + solvated W₁/n₁/ratio. Winner (best joint
score) re-run at N = 500 for the bar-level verdict. Findings §4l + log.

**Pre-registered predictions (2026-07-16, before any build):**

- **S2-P1 (control anchors):** C3 reproduces the known current-law
  behavior — detected KE curve ≈ the §4k Step-0 current-law bracket
  (0/12 bar bins), solvated read per the Wave-7/8 probe results.
- **S2-P2 (pilot, KE):** C1/C2 land n = 1…6 within ×1.4 of the
  reference at N = 50 (bar ×1.25 relaxed for pilot statistics + the
  known 1D↔MD deltas: real RRK cascade vs frozen fate map, m(t)
  feedback −2.2 % on K, real shed vs no-shed). Deep bins reported, not
  gated, at N = 50.
- **S2-P3 (histogram split):** C1 lands n₁ ∈ [0.24, 0.38] with
  ratio ≥ 1.6; C4 reproduces the KE curve but caps at n₁ ≤ 0.24 —
  the floor1/rq4graded split of I51 survives the real cascade.
- **S2-P4 (bar verdict):** the N = 500 winner passes the I-D5 bar
  (×1.25, n = 1…12) after profiling the coherent bands. If it fails
  while the pilot passed at ×1.4, the discrepancy localizes to the
  frozen-fate-map ↔ RRK mapping — a §4l finding either way.

**Gating:** Slices T1–T4 are repo-code + run-artifact changes — they
execute only under a fresh `[PROCEED TO IMPLEMENTATION]`. CALIBRATION_MAP
propagation (v_c, p_tail enter as Bounded→Derived; τ reclassified per
RQ9) lands with Slice T1. Nothing discharges F5.

> **NB (2026-07-16, post-T3):** the T3 pilots executed and surfaced the
> §4l structural finding (I52–I54): the undressed geometry parks every
> config at n ≈ 7–9 and S2-P3 is structurally silent. **Slice T4 as
> written is superseded** — scoring the undressed pilots would rank
> configs on a distribution the closure never predicted for them. The
> T4 scoring machinery and the winner/N = 500 gate are **absorbed into
> §I.11 Slice T9**; the four T3 dirs stay on disk as the A-leg of the
> §I.11 oracle chain.

## I.11 Step 2c — geometry closure: reproducing the Step-1c closure ensemble in MD (slices T5–T9; designed 2026-07-16)

> **Boundary.** Plan, not code — the strict Physics-Definition /
> Software-Implementation boundary holds; every slice executes only
> under its own `[PROCEED TO IMPLEMENTATION]`. Entry documents:
> findings §4l (I52–I54, the T3 structural finding), §4j (the H.2b
> lever set), §H.3b (the parked dressing design, revived here),
> `drag_migration_log_tier2.md` "Slice T3 DELIVERED".

### I.11.0 Why and what "reproduce" means

The Step-1b/1c joint closure (I49/I51: 24 full-house cells around
(p = −1, v_c ≈ 7–7.5, τ ≈ 3.8–4.1, rq4graded, E₀ ≈ 0.23–0.25)) was
scored on the full H.2b ensemble: **L1** birth position with a 3–6 Å
margin, **L2** depth-dressing n_eject(d) = round(21·ρ̂(d)), **D2**
E₀-coupling E₀(d) = E₀·(Σ(n_eject)/Σ(21))^p with p = 1, **D4** droplet
priors, plus the delivered L3/L4 (ladders, capped drag). Slices T1–T3
put L3/L4 and the position *axis* in the repo; the ensemble levers
(L2, D2, margin realization, priors) exist only in the 1D model. Step
2c builds them as interchangeable, byte-inert-default arms and re-runs
the confirmation **stagewise, each stage A/B'd against the 1D twin
re-scored at exactly that MD configuration** — the twin adapts to the
MD (delta prior, realized birth law), never the reverse, so every
stage is a quantitative oracle rather than a leap.

> **NB (2026-07-16, user decision) — immediate goal sharpened:
> reproduce the scratchpad (1D chord-model) closure in full 3D MD,
> twin-parity birth law included.** The V0-1 quantification found the
> repo Boltzmann sampler **center-pinned** at T = 0.4 K (see the V0-1
> adjudication below), so the position axis the closure ensemble rides
> is structurally inert under the physical sampler — T5's dressing and
> T6's coupling would stay inert even after being built. Decision: the
> Boltzmann law stays the physical default (byte-inert); the MD gains
> the 1D ensemble's uniform-in-volume r² law + hard surface margin as
> an explicit interchangeable arm (Slice T7, re-scoped unconditional
> below) — a **twin-parity/capability lever, not a physical claim**;
> whether physical heterogeneity is instead carried by the droplet
> prior (T8) or the E₀ axis stays a later arbitration. The oracle
> chain gains a leg **A′** (= A + `uniform_volume`) so the position
> axis flips *first* — and the "twin adapts to the MD" rule is
> unchanged: each leg's twin re-score uses the MD's *realized* birth
> law (center-pinned delta at leg A, r² + margin from A′ on).
> Sequencing: **V0-3 first** (commit the twin, with a birth-law
> input), then T7 → A′ → T5 → T6 → B/C → T8 → D.

### I.11.V0 Pre-slice verification (zero-MD investigation; may run with the T5 build)

1. **Birth-law reconciliation (decision).** The repo's off-center
   sampler is Boltzmann-weighted — `r²·exp(−U(r)/k_B T)` at
   T = 0.4 K through `droplet_potential`
   (`sampling/radial_positions.py`) — not the L1 bare-r² + hard
   margin. Quantify the realized radial law and its effective surface
   exclusion vs the firm 3–6 Å band. Proposed mapping (user decides):
   **accept the Boltzmann law as the physical realization of L1**
   (temperature-sourced, parameter-free — more physical than the 1D's
   hard margin) and adopt it into the 1D twin; the hard-margin knob
   (T7) is then built only if the realized exclusion falls outside the
   firm band.

   **ADJUDICATED 2026-07-16 (quantified + decided; proposed mapping
   REJECTED).** Preliminary quantification (scratch integration of
   the repo density at R = 27.94 Å / N = 2000, T = 0.4 K, molecule
   steepness 14.3324 Å, E_b = 49.4 meV; to be re-issued from the
   committed twin at V0-3 for reproducibility): the realized law is
   **center-pinned**, not r²-with-soft-exclusion — median birth
   radius 1.37 Å, 95 % < 2.75 Å, 99 % < 3.41 Å, mode 1.16 Å. The
   wide erf tail leaves U ≈ 4.2·k_BT at the droplet *center*, so the
   inward gradient beats the r² volume factor everywhere. The law
   shape is the *opposite* of L1 (L1 median r ≈ 17–20 Å at margins
   3–6 Å), so adopting Boltzmann into the twin would leave the
   closure's position axis — and with it T5's dressing and T6's
   coupling spreads — structurally inert (the T3 lesson repeated one
   level up). Decision (user): **reproduction-first** — Boltzmann
   stays the physical byte-inert default; the MD gains the 1D law as
   an explicit arm (T7, unconditional); the twin gains a birth-law
   input at the V0-3 commit so every T9 leg re-scores at the realized
   MD law. Physically, center-pinning of a heliophilic dopant is
   defensible — which is exactly why the r² arm is recorded as a
   capability lever, not a physics claim.
2. **Trapped-class convention (premise REFUTED + read EXECUTED
   2026-07-16).** The claim originally recorded here ("the droplet
   well enters the ion stage as E_pot *bookkeeping only* — no force on
   the drag path") is **false** and is withdrawn. The drag-path
   conservative acceleration is Coulomb **+ droplet well**
   (`_ion_accel_fn` via `make_ion_accel_fn`, `physics/leapfrog.py` —
   the documented single source of the ion `acc_fn`, consumed by the
   BAOAB B/A kicks and rebuilt every step in `simulation/ion.py` for
   *all* mass scenarios incl. `biphasic`; the E2 relaxation
   `"coulomb"` mode binds the same `acc_fn`). Verified numerically
   (surface ion, Coulomb partner at 500 Å: a_x = −0.5694 Å/ps² vs
   −0.5667 droplet-only analytic; residual = the Coulomb tail). The
   well the ion feels is the Method-B jointly-fitted
   `effective_binding_energy_I_ion_eV = 0.11676 eV` (bundle stamp,
   §6.5.1 exact-pairing guard; presets copy it into
   `binding_energy_I_ion_eV`; confirmed in the T3 `cfg.json`s) — the
   *same* ≈ 0.117 eV barrier the twin's §4j trapped class uses, so
   twin and MD already share one barrier and no reconciliation rule
   is needed. **Tier-0/1a trajectory matching is untouched** — it was
   calibrated with the well force on (that is what the §6.5.1 stamp
   records). **Trapped-class read of the T3 pilots (zero-MD artifact
   re-read):** 0/100 ions bound in every config C1–C4 at both the
   30 ps and 1000 ps reads (escape margin KE + U − E_b positive
   everywhere; slowest ion C3 at +0.167 eV, all others ≥ +0.212 eV;
   all ions outside the droplet already at 30 ps). The class is
   therefore *expressible* in MD but *unpopulated* at the undressed
   configuration — a geometry-conditional read, re-taken at each T9
   leg. Scoring convention fixed now: if bound ions appear in any
   leg, they are counted as droplet-retained weight, excluded from
   the IHe_n histogram, and reported alongside it (matching the
   twin's trapped-class bookkeeping).
3. **The 1D twin becomes reproducible (decision).** `h2b_feasibility.py`
   survives only in expired session scratchpads (two copies located).
   Proposal: recover, verify its recorded wiring oracles (K 0.74460 /
   0.89767, Σ(21) = 0.18783720), and commit it (e.g.
   `scripts/tier2_h2b_forward_model.py`) so every T5–T9 oracle
   comparison is reproducible in-repo. Without this the twin re-scores
   are unverifiable hand-me-downs. **Scope addition (2026-07-16):** the
   committed twin gains a **birth-law input** (r² + margin — its native
   L1 law — plus the realized Boltzmann/center-pinned law) so each T9
   leg re-scores at the MD's realized configuration; the V0-1
   quantification is re-issued from the committed script as its first
   reproducible output. **Sequenced first** among all V0/I.11 items
   (user decision) — both remaining V0 items and every twin re-score
   depend on it.

   > **NB (2026-07-16): V0-3 DELIVERED** under its own trigger (log
   > entry "V0-3 DELIVERED"). The two surviving scratchpad copies are
   > **byte-identical** (cmp); the recovery is committed as
   > `scripts/tier2_h2b_forward_model.py` — code verbatim except the
   > V0-3 scope: repo-relative bootstrap, outputs to the gitignored
   > `data/runs/h2b_forward_model/`, the `BIRTH_LAW` input
   > (`uniform_volume` native / `boltzmann` realized-repo-law with a
   > margin guard, per-law fragment cache), and the new `birthlaw`
   > stage. **All recorded wiring oracles reproduce exactly:**
   > Σ(21) = 0.18783720 eV, production center-pin K = 0.74460, 9 Å
   > K = 0.89767 (plus O-L2/O2/O3/O-B/O4/O-dt as recorded). The V0-1
   > quantification re-issued reproducibly: Boltzmann analytic median
   > 1.36 Å / q95 2.75 / q99 3.40 / mode 1.16 (the adjudication's 1.37
   > was the same integral on a slightly different grid convention);
   > empirical draw through the actual repo sampler agrees (median
   > 1.35). Tests: `tests/test_tier2_h2b_forward_model.py` (10 —
   > oracle pins, fate-map O4, ladder variants, both birth laws,
   > margin guard, birthlaw smoke).

### Slice T5 — initial-shell dressing arm (revives H.3b verbatim)

The §H.3b design carries over unchanged, now with its revival clause
met (the NB there is superseded by §4l: the MD *did* falsify the
undressed geometry). Config enum `initial_shell_model ∈ {"full"
(default, byte-inert = delivered 21-for-all), "density_tied"}`:

$$n_{0,i} = \operatorname{round}\big(n^{*}\cdot\hat\rho(d_{{\rm birth},i})\big)$$

ρ̂ through the **shared erf-complement surface** (single source with the
drag/pickup/cooling gates; dimensionless count — dimensional analysis
trivial). Per-ion consequences at the biphasic seed
(`ion_initial_state.py`): `n_shell(0)`, `mass_kg(0) =
complex_mass_amu(n₀ᵢ)`, the t0 `e_bind_pair(n₀ᵢ)` E_pot fold — arrays
are already per-ion, **no checkpoint schema change**; the 5-term
invariant closes per-ion as before. Biphasic-only config guard.
Oracles (H.3b): `full` arm byte-identical to a delivered dir;
center-pinned + `density_tied` ≡ n₀ = 21 exactly (inert without the
position axis); n₀ monotone non-increasing in birth radius. Boundary
carried: first-order occupancy statement; the dressing↔pickup interplay
(under-dressed ions re-filling via the live Langmuir channel — a
channel the 1D twin does *not* have) is read from the runs and reported
as a twin-divergence candidate.

### Slice T6 — E₀–dressing coupling (the D2 p-law)

Config enum `internal_energy_partition_law ∈ {"constant" (default,
byte-inert ≡ p = 0), "sigma_proportional" (p = 1)}`:

$$E_{{\rm int},i}(0) = f_{\rm int}\cdot E_{\rm avail}\cdot\big(\Sigma(n_{0,i})/\Sigma(n^{*})\big)^{p}$$

ladder-resolved (the ratio rides the injected table under
`tabulated`; dimensionless — units unchanged). Structurally inert
without T5 (ratio ≡ 1 at n₀ = n*), byte-inert default, hand-oracled
per-ion values. Both arms built; the 1D verdict (p = 1 wins, p = 0
over-suppresses — §4j finding 1) is a prior, and the arm is swept in
the T9 A/B chain, not hard-wired.

### Slice T7 — birth-position-law arm (re-scoped UNCONDITIONAL 2026-07-16; was: conditional hard-margin knob)

The V0-1 conditional fired in its strongest form: the realized
Boltzmann law is center-pinned — a *law-shape* mismatch, not a margin
mismatch — so the margin-only knob is superseded by a law selector.
Config enum `birth_position_law ∈ {"boltzmann" (default, byte-inert =
the delivered `sample_radial_positions` path), "uniform_volume"}` plus
`initial_position_margin_angstrom` (float ≥ 0, default 0.0; read only
under `uniform_volume`, guard-refused if set under `boltzmann` — no
silent carry):

    p(r) ∝ r²  on  [0, R − margin],  0 outside    [dimensionless]

— uniform in droplet volume with a hard surface exclusion, exactly
the 1D twin's L1 law; margins per the H.2b firm band {3, 4.67, 6} Å.
The arm is a **twin-parity/capability lever, not a physical claim**
(the Boltzmann default remains the physical arm; see the I.11.0 NB).
Oracles/tests: `boltzmann` default byte-identical to a delivered dir;
`uniform_volume` realized quantiles match the analytic r² CDF on
[0, R − m] (fixed-seed quantile pins); margin 0 reaches the surface;
guard behavior. The law is exactly the scratchpad's — hard margin, no
smoothing (user decision 2026-07-16).

> **NB (2026-07-16): Slice T7 DELIVERED** under its own trigger (log
> entry "Slice T7 DELIVERED"; TDD, RED watched first — ImportError on
> the missing guard). As built: `BirthPositionLaw` Literal +
> `birth_position_law`/`initial_position_margin_angstrom` fields
> (config.py, next to `single_initial_position`);
> `check_birth_position_config` guard wired into `SimConfig.validate`
> (typo guard, finite/≥ 0 margin, margin-under-boltzmann refusal);
> `_sample_uniform_volume` in `sampling/radial_positions.py` — exact
> inverse-CDF `r = (R − m)·U^(1/3)` per molecule, loud failure when
> the margin consumes any droplet — dispatched *before* the Boltzmann
> code so the default path is untouched. **Byte-identity proven
> exactly**: the default-path draws equal the pre-T7 (`git show HEAD`)
> sampler's draws elementwise at fixed seed (stronger than the
> delivered-dir oracle, which it implies at the sampler level);
> pre-T7 `cfg.json` loads with the inert defaults (Slice-DS
> precedent). Tests: `tests/test_birth_position_law.py` (16 — config
> surface/guards/back-compat/round-trip, r² CDF quantile pins, hard
> margin edge, twin-comparator median 18.47 Å at m = 4.67,
> per-droplet support, boltzmann-vs-uniform liveness, initial-state
> integration). Full suite **2354 passed**. CALIBRATION_MAP row 25 +
> tally propagated (law = arm; margin Bounded {3, 4.67, 6} Å). No
> preset or generator selects the arm yet — that is the T9 leg-A′
> re-pilot's job.

### Slice T8 — droplet prior (the D4 family; supersedes the S2-D4 deferral)

The closure cells are per-(prior, margin); the deep tail carries the
W9-P2 ≥ 9.6 % below-floor demand — full reproduction needs the droplet
axis. Rule-1 first: audit the legacy `use_single_droplet_size=False`
machinery before adding surface; then wire the Wave-10/11 family —
Kornilov log-normal δ ∈ {0.40, 0.625, 0.80} about ⟨N⟩ = 2000 +
the pickup-weighted ∝ N^(2/3) variant — behind config, with
R(N) = 2.2173·N^(1/3) (the existing convention). Scope note: the T9
oracle chain runs **fixed-N first** (the twin re-scores at a delta
prior exactly), so T8 sequences last of the physics slices and its
pilot cost is bounded by the chain's final leg.

### Slice T9 — dressed confirmation re-run + scoring (absorbs §I.10 T4)

The staged A/B/C/D oracle chain at the C1 knob point, one arm flipped
per leg, each leg pre-registered against the committed 1D twin
re-scored at exactly that configuration:

- **A** = the delivered T3 dirs (all arms off) — byte-identity
  regression anchor;
- **A′** = + `uniform_volume` birth law at the closure cell's margin
  (T7) — the position axis flips first (2026-07-16 re-scope: T5's
  dressing is inert under the center-pinned default law; S2c-P1's
  small-n/broadening statement reads at this first position-live leg);
  **EXECUTED 2026-07-16 — AP-P1..P4 all CONFIRMED** (findings §4m +
  I55–I57; W₁(MD, twin) 0.55–0.69 bins, suppressed ordering exact;
  the KE composition diverges through the real cascade — I57; the
  V0-2 retained policy + barrier-corrected bound criterion delivered
  under the leg trigger; E2-drag OQ opened);
- **B** = + `density_tied` (T5);
- **C** = + `sigma_proportional` (T6);
- **D** = + droplet prior (T8; fixed-N legs A–C).

Then the confirmation matrix re-pilot at the final configuration
(C-values re-centered per the twin's re-score — the Step-1c basin
was located under the 1D ensemble and need not sit at the same
(v_c, τ, E₀) under the MD-realized one), N = 50 → the **T4 scoring
machinery built here** (ihe_ked run-summary mean-to-mean with the
I-D4 correlated bands + solvated W₁/n₁/ratio; conf-namespace-aware —
the staircase report's `*_tier2probe_*` glob must exclude
`*_tier2probe_conf*` or its rows be disregarded) → winner at N = 500
for the bar-level verdict (the absorbed §I.10 S2-P4 gate, bar ×1.25
n = 1…12 after band profiling).

**Pre-registered predictions (numeric values filled from the twin
re-scores before each leg executes — program convention):**

- **S2c-P1.** Leg B develops small-n weight fed by the **near-cliff
  chord-K band** (I45: not a fast surface class — none exists
  in-bounds); the n-histogram broadens from the ±2 cluster toward the
  twin's per-cell shape.
- **S2c-P2.** Leg C moves weight between the suppressed/bare side and
  the shallow bins per the p = 1 uniform-K* mechanism; p = 0 remains
  over-suppressed (the §4j direction).
- **S2c-P3 (revived S2-P3).** At the final configuration the
  rq4graded/floor1 split engages: rq4graded reaches n₁ ≥ 0.26 /
  ratio ≥ 1.75 territory; floor1 caps near ratio ≈ 1.8 (I51) — the
  I51 discrimination, now in MD.
- **S2c-P4.** Twin-divergence channels are the *listed* ones only
  (pickup re-filling after under-dressed birth; no trapped class in
  MD; real RRK cascade vs ε = 0 floor map) — an unlisted divergence
  is a model-structure finding.

### I.11.1 Sequencing, gating, boundaries

**V0-3 first** (commit the twin + birth-law input; V0-1 adjudicated,
V0-2 closed as corrected — see the log) → T7 → T9 leg A′ → T5 → T6
(needs T5 to be live) → T9 legs B/C (fixed N) → T8 → T9 leg D +
re-pilot + scoring (re-sequenced 2026-07-16; supersedes the "T7
conditional on V0-1" clause — the conditional fired, T7 is
unconditional). Every slice TDD behind its own trigger; every new field/enum
byte-inert by default with a back-compat cfg.json test (the Slice-DS
precedent). CALIBRATION_MAP: T5/T6 add **arms, not knobs** (the tied
law and the p-law are parameter-free); T8 adds the prior selection
(Bounded, δ family Sourced from Kornilov). Nothing here discharges F5;
the bare bin stays RQ8-gated; RQ4 remains the external arbiter of the
rq4graded taper.

### I.11.2 Post-leg-A′ open decisions (discussion 2026-07-17; entry point for the next session)

Leg A′ is EXECUTED (findings §4m + I55–I57; log entry "T9 leg A′
EXECUTED"). AP-P1..P4 all confirmed; the one divergence is **I57** —
histogram parity without composition parity, visible on the KE axis.
Three items follow, in the suggested order:

> **Status (2026-07-17, user decisions + execution):** **item 1
> EXECUTED** (zero MD; findings §4n + I58–I60; the pre-registered
> Coulomb-share hypothesis is *refuted as dominant* — the driver is the
> **cold-shed momentum convention** of the delivered evaporation
> channel, +0.91 eV on n₁; the co-moving counterfactual matches the
> twin to ≤ 2 % → **OQ-J / RQ10 fired**; **shed convention ADJUDICATED
> same day: co-moving is the working convention for the twin-parity
> legs** — enum `evaporation_shed_convention ∈ {cold (byte-inert
> default), co_moving (stamped by T5+ legs)}`, **build DELIVERED
> 2026-07-17** (log entry "Shed-convention enum DELIVERED"; row 26);
> RQ10 stays open as the physical-resolution question,
> coupled to RQ2's ε). **Item 2 ADJUDICATED:
> arm (c) Landau-gated drag** selected by the user; build stays behind
> its own `[PROCEED TO IMPLEMENTATION]`, required before any N = 500
> run. **Item 3 ADOPTED and sharpened by item 1:** leg-B twin
> re-scores pre-register per-bin mean detected KE alongside the
> histograms, **stating the shed-convention basis of every KE claim**
> (the twin is co-moving by construction; the delivered MD is
> cold-shed) — without the stated basis the two KE axes are not
> comparable (I59).

**1. The n₁-composition re-read (zero MD, on-disk apc dirs; do this
before T5).** The undressed-MD n₁ bin sits at mean KE 2.42–2.48 eV
(experiment 1.302 eV; twin at A′ 0.97 eV). The kinematics already
prove short chords / low drag cannot be the whole story: at the n₁
detected mass (~131 amu), 2.48 eV means v ≈ 19.1 Å/ps — *above* the
symmetric zero-drag limit of 15.4 Å/ps (both fragments at m₂₁ through
the Coulomb acceleration). The excess requires the **mass-asymmetric
Coulomb split**: E_A = E_pair·m_B/(m_A+m_B), so a fragment that sheds
*early, mid-acceleration* is light while the force acts and takes a
larger share (instant-bare vs m₂₁ partner → 3.37 eV, v = 22.6 Å/ps;
the measured 19.1 sits between). The A′ twin is structurally blind to
this route (both chords pinned at m₂₁ — the fate map has no in-flight
m(t)), which is I57's mechanism. The near-edge-birth hypothesis
(short outward chords, early gate-open) is the plausible *correlate*
— early shed ⇒ light early ⇒ less drag AND bigger share, the two
mechanisms compound. The re-read decomposes the n₁ (and bare) bins:
per ion, birth radius, chord cosine, shed-time history (from the
stored n(t)), and the partner's mass history → attribute the KE
excess drag-deficit vs Coulomb-share. Zero new MD; probe-convention
scratch read; sharpens exactly what leg B must move.

**2. The E2 dissipation adjudication (user decision; required before
any N = 500 run).** The E2 relaxation translation is **zero-gamma**
(conservative: Coulomb + well; pickup also off, λ₀ = 0) — correct
by construction for ejected ions (ρ̂ underflows ⇒ all arms
identical), but the A′ position axis produces still-helium-coupled
ions at handover, for which E2 omits real physics (leg-A′ execution
record: the centrifugal-resonance "ion 9", the 1000→8000 ps cap
chase, the barrier-criterion classification). Candidate arms:
  (a) **zero-gamma** (status quo): resonances persist; handled by the
      `_conservatively_bound` classification + a leg-level cap;
  (b) **drag-live E2**: extrapolates the pure-cubic γ = b·v² below
      its Tier-0 calibration band (in-window speeds ~5–15 Å/ps) —
      captures every marginal ion dynamically, but may *overdamp*
      slow motion;
  (c) **Landau-gated drag**: γ = 0 below the config's Landau cutoff
      (`v_limit` — the legacy E_min machinery), drag above — the
      slow orbits in question (|v| ~ 0.3–1 Å/ps) straddle exactly
      that regime, and sub-Landau superfluid motion is physically
      dissipationless, so (a) may be *more* physical than (b) below
      the cutoff. Arguably the most physical arm; needs the
      domain-expert read.
Scope: affects only the droplet-retained/marginal class (5–11 % at
A′) — never the ejected read; but N = 500 makes marginal ions a
certainty, so the arm must be selected (and built behind a byte-inert
enum) before the T9 endgame.

**3. Leg-B (T5) pre-registration amendment.** After I57, histogram
agreement alone no longer counts as twin parity: the leg-B twin
re-score must pre-register **per-bin mean detected KE** alongside the
histograms/classes, and the A/B verdict reads both axes. Note the
1D stake: the Step-1c full-house cells passed the KE bar (12/12
within ×1.25, incl. n₁ vs the experimental 1.302 eV) on the *dressed*
ensemble — so the dressed twin claims n₁ ≈ 1.0–1.6 eV, and T5 also
weakens the early-shed asymmetry route (same-depth birth ⇒ symmetric
birth masses per pair; less shell to shed near the surface). Whether
the real cascade agrees is leg B's question. Bracket worth keeping in
view: undressed MD n₁ KE spans 0.67 eV (c3, current law) to 2.48 eV
(capped tails) — the experimental 1.302 sits inside the form bracket.

**Suggested order: (1) → (2) → T5 trigger** (leg B with the amended
per-bin-KE pre-registration). All three stay behind their own
adjudication/trigger; nothing here discharges F5.
