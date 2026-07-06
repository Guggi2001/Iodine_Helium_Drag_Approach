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
