# THESIS_FINAL_GOALS — agreed changes and tests before thesis completion

Settled in the final discussion round (2026-10-02). This document holds
only agreed goals. Results go to the migration log and D0 as usual.

**Status (2026-10-13)** — entry point for a new chat; details in each
section's own status, D0 §1.G4.x (`TIER2_PARAMETER_INFLUENCE.md`) and the
log (`drag_migration_log_tier2.md`).

1. **G4 / G4.8** (Stage-2 refit at s_ρ = 3.5 Å) executed; adoption
   superseded by G4.9.
2. **G4.9** (18 Å-only refit, 9 Å held out): winner `cubic`
   (`srho3p5A18A_cubic`, b 2.2286, E_bind 0.1214 eV).
3. **G4.10** (force balance): 18 Å b_A 1.618 outside [1.803, 2.654] —
   inconclusive (references not derivative-grade); post-hoc June
   pipeline: cubic b 18 Å 2.449 / 9 Å 2.634, both in band.
4. **G4.11** (cross-case drag-law test): only the cubic consistent
   (s_X 1.182, s_M 1.099). **Contested** by the radial 9 Å runs: the cubic
   over-dissipates (tracks transverse atom 2). **Adoption of
   `srho3p5A18A_cubic` on hold.**
5. **G4.12** (TDDFT-tied ridge F = a·v + c·v², uncapped): steps 1–3 done.
   Every ridge point fits 18 Å; the 9 Å exit bound fails only the pure
   quadratic; all 11 passing points land in the twin (after the frozen
   near-pass refinement), on one KE₁ / midHot dial.
6. **G4.13** (MD, production pipeline, c = 2.0): 7/16 Stage-A cells gate
   on one seed; the escalated cell t64e310 (N = 6000) misses n₁ (0.186),
   ties h405 on KE₁ (0.638 vs 0.640), W₁ 1.10 vs 0.755. **"Dominated"
   holds for that cell only** — the grid could not resolve the landing
   (user-caught correction).
7. **G4.14** (fine τ-line sweep, τ the primary axis): **Steps 0a–1
   EXECUTED 2026-10-06** (74 MD runs, independent seeds). No c / τ / E₀
   beats h405 on the frozen ranking; τ trades the histogram's large-n tail
   against KE₁ at every c. The size-histogram misfit is spread, not
   location (no n ≥ 17; n₁ low): a mass-model question, not a drag / τ / E₀
   one. Amendments A1–A5; A5 = KE-led criterion + one pooled confirmation
   (c 2.0, τ 6.4): **loses on KE** (KE₁ 0.617 vs 0.640 / 0.643; χ² 507
   vs 456 / 636). Nothing adopted.

8. **Thesis discussion (2026-10-06):** lq vs capped compared; n = 0
   excluded on physical grounds; the n ≥ 16 limit explained (escape
   cliff, TDDFT-pinned). See G4.14 "Discussion outcome".

9. **G8 bare-start screen** (2026-10-06): fragments born with no He
   shell (n₀ = 0, E_int(0) = 0), lq A5 cell and h405, one run each —
   what leaves the droplet, at which n and KE. **Executed:** bare ions
   land at n = 2–4 by pickup; KE₁ 0.49 (lq) / 0.87 (h405) vs exp median
   1.13 — the mixture cannot supply fast n = 1 ions. See G8.

   **Decision (user, 2026-10-06):** mixture not pursued; the n = 1 KE
   deficit goes to the thesis outlook (levers exhausted inside the
   model; see G8 status).

10. **G9 cold-onset minority screen** (2026-10-08, frozen; awaits
    `[PROCEED TO IMPLEMENTATION]`): lq A5 cell, a fraction w up to
    0.05 eV below the main E₀ (point / uniform), main E₀ re-placed to
    n̄ 3.89; 6 seed-s1 columns, mixtures by exact pooling.
    **Executed:** no tail (n ≥ 15 ≤ 0.002); the tail needs cold ions at
    E₀ ≲ 0.18 eV. See G9.

11. **Decision (user, 2026-10-08): ladder and high-n tail → outlook.**
    Ladder: no fix (edge width inert; n* beyond 21 ± 1–2 not
    defensible; D0 §1.G9b); outlook = sourced I⁺Heₙ D₀(n) and n*
    (RQ4 / M1) — the bottom taper is data-derived. Tail: n and KE are
    both set by exposure, so no in-model knob gives fast high-n ions;
    outlook = a mechanism setting n independently of exposure (per-ion
    energy scatter, residence-independent cooling).
    **Discussion (2026-10-08):** `THESIS_DISCUSSION.md` §7 — the tail
    needs a higher escaper exposure ceiling *and* a blurred escape
    threshold (the model has a gap at X/τ 1.86–5); outlook = collision
    shot noise (Tier 3), the same missing ingredient as the narrow KE
    width at fixed n.

12. **Shallow well + drag retune (2026-10-08, measured, pooled N = 6000;
    not Tier-0-consistent):** E_bind 0.0482 eV with a 19.5 (c 2.0) at
    midHot ≈ 1 — n ≥ 15 0.0053 (A5 0.0006), max n 18, W₁ 0.905 (A5
    1.097), KE₁ 0.597 (A5 0.617). Thesis framing (user): A5 pool =
    TDDFT-based result, 0.0482 / a 19.5 = effective improvement for higher
    fragments, departing from the 18 Å trace (Tier-0 RMSE 0.557 vs A5
    0.122; no a rescues E_bind 0.0482, best 0.226 > 0.19). Table and
    caveats: G10 below; D0 §1.G4.14.

13. **G5 executed (2026-10-08/09, pooled N = 6000):** `post_pickup`
    (realistic droplets, ⟨N⟩ 16.7k vs 12.9k) on lq A5 raises n̄ by +0.70
    and leaves KE / shape flat (the chord → (n, KE) map is
    sampler-invariant; the n̄ gap is the n = 1 deficit). E₀ or τ restore
    n̄ only at midHot ≈ 0.85; **the drag does it:** ridge c 1.6 (a 18.007,
    E_bind 0.1173, still TDDFT-tied) + E₀ 0.3297 holds A5 on every
    observable and improves W₁ / S. The G10 recipe on top (E_bind 0.0482,
    a 19.8) gives the effective high-tail variant. Comparison table: G5.

14. **G11 (agreed 2026-10-09):** thesis reframed as
    "TDDFT constrains (a ridge), the experiment selects (c, E₀)"; Tier-0
    extraction mass → m_eff = 21 He (18 Å shell mass; supersedes G3);
    Tier 1 kept as the fixed-mass consistency check (not void; G8 is a
    different question).

15. **G11.4 steps 2–3 executed (2026-10-09/10):** Tier-0 ridge refit at
    m = I⁺ + 21 He + 9 Å exit bound — B2 / A5 pass unchanged (step 4
    re-landing dropped); h405 fails the 9 Å bound (−1.21 Å/ps). Tier-0
    figures: `tier0_drag_comparison.py` (thesis trace + R(t), B2 at
    m = 21) and `tier0_thesis_figure.py` (B2 / B3 / h405). The lq bundle
    `srho3p5A18A_linear_quadratic` now carries B2. Records: G11.4 status.

16. **Tier 1 retired (user, 2026-10-13):** no Tier-1a redo (G11.4 step 5
    dropped); the old Tier-1a results are dropped from the thesis. The
    mass question moves from the validation chapter to the discussion
    (mass-influence section, compiled from existing results, zero MD;
    the bare start = outlook). See G11.3 / G11.5.

**Next:** G11.4 step 6 — G6–G7 (G6 deferred by the user) + storyline;
mass-influence discussion section per G11.5 (steps 2–3 done, steps 4–5
dropped; see G11.4 status). G4.14 recap point (a) settled by the decision below; (b)–(c)
parked.
**Thesis production law (user, 2026-10-09; supersedes A5 of 2026-10-08):
B2 — `post_pickup` sampler, linear + quadratic ridge c 1.6** (a 18.007,
E_bind 0.11727, 3.5 Å gate, τ 6.4, E₀ 0.3297; reference
`…N6000_detfix_conf270_g5lqppc160e3297pooled_th500`). **Effective
high-tail variant: B3** — B2 with E_bind 0.0482 eV, a 19.8, E₀ 0.3389
(`…g5lqppc160a198w0482e3389pooled_th500`); not Tier-0-consistent (G10.5).
A5 (raw sampler, c 2.0) = the previous production law; h405
(`capped_cubic` v_c 5.5) = the contrast law for the identifiability
argument. Generator: `gen_tier2_detfix_postpickup.py --cell lq --b2` /
`--b3`; `gen_tier2_detfix_battery.py` still builds h405 (repo default not
switched).

## 0. Why the whole chain is redone

Three agreed changes together invalidate the current drag bundle and the
h405 standing point:

1. a sourced He **density** width, which changes the Tier-0 refit (G2);
2. the `post_pickup` droplet-size sampler, which moves the landing (G1);
3. the bundle change makes h405 void anyway, so the re-arbitration
   (G5) is done once, under the new sampler.

The order below keeps an attribution ledger: each step changes one thing,
so the thesis can report the bundle effect, the sampler effect and the
re-arbitration separately (CLAUDE.md validation hierarchy: one tier at a
time).

## G1. Production sampler → `post_pickup`

- **Change.** `droplet_size_sampler_mode="post_pickup"` in production
  (with `droplet_size_prior="legacy"`, as now).
- **Rationale.** The experiment selects doped droplets. `post_pickup`
  models that selection (pickup cross-section ∝ N^(2/3), exactly one
  pickup, pickup-induced evaporation) and is the legacy MATLAB production
  sampler. The basis for `raw` was a circumstantial quantile match whose
  discriminating margin (≈ 8 % in R) equals the N→R convention spread.
  Since a re-arbitration is needed anyway (G5), there is no cost argument
  for keeping `raw`.
- **Evidence.** D0 §15.8 reads 1–4 (A/B at h405: n̄ +0.72, n₁ −0.033,
  W₁ +0.24, trap +0.057, supp −0.059, KE₁ +0.007 eV; h405 fails both MD
  gate clauses under `post_pickup`); `detfix_postpickup_table.csv`;
  log entry 2026-10-02; `CALIBRATION_MAP.md` row 30.
- **Known risk.** The `post_pickup` r6 run tripped the D4
  marginal-partner safeguard (one slow escaper still coupled at
  t_h = 500 ps). Larger droplets may make this more frequent. Track D4
  firings in every `post_pickup` battery.
- **Cost.** Config decision only; the generator
  `gen_tier2_detfix_postpickup.py` already exists.

## G2. He density width s_ρ = 3.5 Å (density only, not the potential)

- **Change.** `drag_spatial_gate="erf_independent"` +
  `drag_gate_steepness=3.5`. `potential_steepness` stays 14.2 Å.
- **Rationale.** 14.2 Å is the DFT fit of the solvation **potential**
  (ρ ⊗ V_I–He + cavity) and is correct for U. Its reuse as the density
  width gives a 10–90 % interface of 25.7 Å. Literature for ⁴He droplets
  (N = 10³–10⁴) gives 5.7 Å (DFT) / 6–8 Å (experiment), Harms, Toennies,
  Dalfovo, PRB 58, 3341 (1998). For the erf-complement profile the
  10–90 % width is 1.81·s, so s = 3.5 Å ⇒ 6.3 Å, inside that range.
  Used for both the Tier-0 refit and production.
- **What it moves.** Every reader of `drag_gate_steepness(cfg)`; treated
  as one coherent physical correction (all readers together, atlas §6.8
  T3), not decomposed:

  | reader | where | Tier 0 | production |
  |---|---|---|---|
  | drag gate | `ion.py` → `ion_propagation_step` | **yes (only reader)** | yes |
  | pickup density gate | `ion_propagation_step.py` `rho_ratio` | no (fixed mass) | yes |
  | Newton cooling gate | same, `cool_rho` (`density_scaled`) | no | yes |
  | detection-stage cooling / exposure at handover | `detection_stage.py` | no | yes (≈ 0 for ejected ions) |
  | birth dressing n₀ = round(21·ρ̂) | `ion_initial_state.py` (`density_tied`, live in production) | no | **yes** |

  Unchanged: the solvation potential and well force, the neutral stage,
  the Boltzmann birth law (`potential_steepness_molecule`), E2 (off).
  Expected: at 3.5 Å, ρ̂ ≈ 1 deeper than ≈ 7 Å, so almost every
  Boltzmann birth (≈ 33 Å deep) is dressed to n₀ = 21 and the dressing
  arm becomes nearly inert (D0 §17 "inert at the parent geometry"). G5
  step A measures the combined size.
- **Where it bites.** A saturated drag is width-blind on total
  dissipation at production R (0.0 %) but not at the 9 Å calibration
  droplet (+11 %), so the effect enters mainly through the Tier-0 refit
  (b, E_bind; E_bind currently inflated ≈ 2.1×). A new bundle changes the
  chord → (n, KE) map, which is the only thing that can move KE₁ (D0
  §15.8 read 6).
- **Evidence.** D0 §17 rows "`potential_steepness` reused as the He
  density width" and "its measured leverage"; D0 §9.4–§9.5.
- **Cost.** Config only.

## G3. Tier-0 extraction mass stays m_eff = 19 He (≈ 203 amu)

> **Superseded 2026-10-09 (user): m_eff = 21 He — see G11.2.** The text
> below is kept as the record of the earlier decision.

- **Decision.** Keep the constant-mass extraction at m_eff = 202.954 amu.
  No change to 21, no time-resolved extraction.
- **Rationale (thesis statement).** The 9 Å TDDFT shell declines
  21 → 19 (10 ps) → 14 (14 ps). In the Method-B balance
  F_drag = m·a − F_C the constant-mass error tracks |m(t) − m_eff|/m_eff,
  so the window-representative 19 He is the minimal difference to the
  TDDFT run over the window (≈ 0 mid-window, ≈ ⅓ only at the ends).
  21 He is correct only at t = 0 and would bias b upward. A time-resolved
  m(t) extraction is not possible for both systems: no 18 Å shell-loss
  curve is available.
- **Evidence.** `DRAG_PORT_DESIGN_DECISIONS.md` §2.2, §6.6;
  `TIER0_FINDINGS.md`.
- **Cost.** None.

## G4. Tier-0 refit (Stage 2, three families) and Tier-1a re-check

This is atlas §6.8 **T5** (`TIER2_SENSITIVITY_ATLAS_PLAN.md`). Code
requires the `[PROCEED TO IMPLEMENTATION]` trigger.

### G4.1 Scope

- **Method-B Stage 2 only** (joint fit over the 9 Å + 18 Å TDDFT traces,
  m_eff = 19 He, s_ρ = 3.5 Å). No per-case refits, no Stage-1 held-out
  prediction, no T6 per-case E_bind diagnostic. The bundles stamp
  `held_out_axis: none (stage-2 only)`.
- **Three families only** (the physically sensible ones), each fitting
  its coefficients + a shared E_bind:

  | family (new name) | force | free | fixed |
  |---|---|---|---|
  | `shared_linear` | F = a·v | a [amu/ps], E_bind | b = 0 |
  | `shared_cubic` | F = b·v³ (γ = b·v²) | b [amu·ps/Å²], E_bind | a = 0 |
  | `shared_linear_cubic` | F = a·v + b·v³ | a, b, E_bind | — |

  Dimensions: a·v and b·v³ both give amu·Å/ps² (force), γ = F/v in
  amu/ps; mass-agnostic (friction convention). All three run on the
  existing `linear_cubic` force function with a or b pinned; the config
  guard already accepts (a > 0, b = 0) and (a = 0, b > 0).
- `diagnostic_4param` stays on as the automatic follow-up if a family
  fails its bands (localizes drag vs binding).

### G4.2 Pre-registered bands and selection rule (frozen before the run)

- **Bands per family** (unchanged METHOD_B §9.4, first-runs rule):
  RMSE 18 Å ≤ 0.19 Å/ps, RMSE 9 Å ≤ 0.45 Å/ps, escape fraction = 1 in
  both cases. A family failing its bands is out.
- **Nested-model test.** The single-term families are the boundaries of
  `shared_linear_cubic`. With the objective in Å/ps:
  - `T_a0`: obj(cubic) − obj(linear_cubic) ≤ 0.005 → cubic equivalent;
  - `T_b0`: obj(linear) − obj(linear_cubic) ≤ 0.005 → linear equivalent.
- **Winner.** Exactly one boundary equivalent → that single-term
  family. Neither → `shared_linear_cubic`. Both → **cubic** (continuity
  with Tier 0; the production cap is defined on cubic).
- **Predictions on record** (frozen before the run, not gates): from D0
  §9.4 the sharper gate moves T → 1 and roughly halves the binding
  refund, so E_bind is expected to **fall**; b's shift is expected small.
  Atlas T2 measured the ensemble response at h405 as null/small: this
  refit is justified by model correctness, not as a KE₁ fix.

### G4.3 Code changes

1. **Gate pass-through.** `build_case_setup` gains a pass-through for
   `drag_spatial_gate` / `drag_gate_steepness` (or a general
   `cfg_overrides`); the driver gets them as USER SETTINGS. Default:
   identical to today.
2. **Escape criterion on the potential width.** `_escape_fraction` uses
   `potential_steepness` (the well's reach), not the density width. At
   3.5 Å the current code would count ions as escaped 7 Å outside the
   surface, inside the 14.2 Å well. Bit-inert for the old bundle (both
   widths 14.2), which the oracle shows.
3. **`shared_linear` joint variant** (b = 0; bounds, unpack, prescan,
   starts, writer), and the three families under their new names.
   `shared_linear_cubic` keeps b ≥ 10⁻³ (the linear boundary is covered
   by its own variant).
4. **Bundle stamps the gate.** New bundles carry `drag_spatial_gate`,
   `drag_gate_steepness_angstrom` (3.5), `potential_steepness_angstrom`
   (14.2), the source (Harms et al. 1998, flagged *borrowed* per atlas
   §6.8.2) and `held_out_axis`. `DragCoefficients` gets an optional
   `extraction_gate_steepness_angstrom`; legacy bundles without the key
   read as 14.2, gate tied to the potential.
5. **Gate↔bundle guard** at config load (alongside the §6.5 mass and
   §6.5.1 binding guards): the run's resolved `drag_gate_steepness(cfg)`
   must equal the bundle's stamped value. Bundle dataclass + loader
   change only; checkpoint schema untouched.
6. **Bundle selection via config.** One `SimConfig` selector chooses
   the drag bundle (family + density width) and resolves to its
   directory in the presets layer, so the catalog, the production
   generator and the scorers switch through config instead of
   per-script constants. Legacy bundles stay selectable (old runs
   reproduce). The frozen presets stay on the old bundle until adoption.
7. **Downstream:** `tier0_common.CATALOG`, the production generator
   (`gen_tier2_detfix_battery.py` / `build_cell`) and the scorers read
   the selector. Run-dir names include the family, so new runs cannot
   overwrite old ones.

Tests: pass-through (default identical), escape uses the potential
width, `shared_linear` variant, loader accepts the new keys, guard
accepts matching / refuses mismatched / reads legacy as 14.2, config
selector resolves every family.

### G4.4 Storage

Old reference data stays untouched. New sibling folders, named after
the physical change:

```
data/reference/drag/shared/trajectory_matching_srho3p5A/
    shared_linear/fit_parameters.json
    shared_cubic/fit_parameters.json
    shared_linear_cubic/fit_parameters.json
    verdict.json            (bands, T_a0 / T_b0, winner)
data/reference/drag/{9A,18A}/tier0_srho3p5A/md_mean_trajectory*.csv
```

### G4.5 Execution order

1. **Oracle:** after changes 1–2, the default gate reproduces the
   committed `shared_pure_cubic` objective (0.12926494 Å/ps at the
   committed b, E_bind) bit-exactly — anchored on the committed artifact.
2. Unit tests (G4.3).
3. `DRY_RUN` at 3.5 Å: determinism + timing per joint evaluation (sets
   the real runtime).
4. Commit the pre-registration (G4.2), then run the three families →
   `verdict.json`.
5. Tier-0 forward check (`gen_tier0_runs.py` +
   `tier0_drag_comparison.py`) with the winner; export the mean
   trajectories to `tier0_srho3p5A/`.
6. **Tier-1a re-check** with the winner (anchored `fixed` vs
   `anchored_discrete` A/B, §2.9 ledger closure); gate: R(t), |v(t)| and
   ledger closure no worse than the delivered result.

### G4.6 Consequence for production (depends on the winner)

Production runs `capped_cubic` {b, v_c, p_tail}, which has no linear
term.

- **cubic wins:** no form change; the new b and E_bind enter via the
  selector.
- **linear_cubic wins:** the capped law needs a linear term (new form /
  enum arm, guard, twin mirror), ≈ half a day.
- **linear wins:** nothing to cap; whether the high-velocity tail still
  needs a cap becomes its own physics question before G5.

### G4.7 Evidence and cost

- **Evidence.** `METHOD_B_trajectory_matching_extraction.md` §9;
  `TIER0_FINDINGS.md`; atlas plan §6.8.2–§6.8.4; D0 §9.4, §17;
  `drag_migration_log_tier1a.md`.
- **Cost.** Code + tests ≈ 4–5 h (changes 1–5 ≈ 3 h, `shared_linear`
  ≈ 1 h, config selector + downstream ≈ 1 h); oracle + DRY_RUN
  ≈ 15 min; extraction runtime unknown until DRY_RUN (rough guess
  1–3 h for three families); Tier-0 forward check ≈ 15 min; Tier-1a
  re-check ≈ 1 h.

### G4.8 Status — EXECUTED 2026-10-02 (measured; adoption DECISION NEEDED)

Code `ebda074`, artifacts `c1f9bb5`; full record in
`drag_migration_log_tier2.md` (2026-10-02 G4 entry) and D0 §1.G4.

- **Oracle:** bit-exact against the pre-`aacc71f` 9 Å reference; the
  committed reference was re-smoothed on 2026-09-29 (9 Å RMSE −1.6·10⁻⁴).
  Refit on the current reference (user decision) + a 14.2 Å cubic
  control: the smoothing effect is nil (Δb −0.001).
- **Winner `shared_cubic`:** b 2.2300, E_bind 0.1277 eV (obj 0.1661; RMSE
  18 Å 0.125 / 9 Å 0.207). `shared_linear_cubic` collapses to a = 0
  (T_a0 ≈ 0); `shared_linear` fails its bands (T_b0 +0.509). G4.6: no
  production form change.
- **Predictions on record:** E_bind *rose* (+9 %, not fell); b −11 %.
- **Tier-0 forward check:** reproduces the fit.
- **Tier-1a re-check:** R(t) and ledger closure unchanged; the |v(t)|
  clause **fails as worded** (0.207–0.222 vs 0.124–0.136 Å/ps), mostly
  inherited from the Tier-0 fit (fixed row +0.083); the variable-mass
  increment is comparable to legacy.

> **DECISION NEEDED.** (a) Record the Tier-1a |v| failure as
> Tier-0-inherited, accept `srho3p5A_shared_cubic`, proceed to G5; or
> (b) treat it as blocking and open the 9 Å misfit (per-case wells
> 0.132 / 0.193 eV in `diagnostic_4param`) before G5. Production stays on
> `shared_pure_cubic` at 14.2 Å until decided.

### G4.9 18 Å-only refit (agreed 2026-10-05; supersedes the G4.8 decision)

- **Why.** The 9 Å reference is non-radial: the compared ion stalls,
  oscillates and exits sideways. A radial MD can match it only by
  absorbing the missing dimension into the drag/binding (METHOD_B §5).
  The 14.2 Å gate's long soft drag region could imitate the long dwell
  time; the 3.5 Å gate cannot. 9 Å is therefore dropped from the fit.
  The G4 Tier-1a |v| failure (scored on the same 9 Å trace) is
  non-informative and no longer gates.
- **Step 1 — confirm the mechanism (no new fits).** 9 Å |v| residual
  vs time at the committed shared-cubic points (14.2 Å control: b 2.5144,
  E_bind 0.1168; 3.5 Å G4: b 2.2300, E_bind 0.1277), plotted against
  the reference transverse velocity. Confirmed if the extra 3.5 Å
  residual sits in the stall / sideways phase.
- **Step 2 — 18 Å-only Method-B fit at s_ρ = 3.5 Å.** Families
  `shared_linear`, `shared_cubic`, `shared_linear_cubic` (same engine,
  bounds, starts and anchors as G4; single case). Bands: 18 Å RMSE
  ≤ 0.19 Å/ps, escape = 1. Selection: the G4.2 nested tests and winner
  rule unchanged (T_a0, T_b0 ≤ 0.005 Å/ps). Control: `shared_cubic` at
  14.2 Å, same reference. Oracle first: a single-case evaluation at the
  G4 point and the control point reproduces the committed 18 Å RMSEs
  (0.125 / 0.1346).
- **Step 3 — 9 Å escape check** of the winner, no refit: the 9 Å escape
  fraction at the 18 Å-fitted E_bind (trapping would be a real failure).
  Recorded, not gating. *Amended 2026-10-05 before any fit result:* the
  9 Å RMSE is not used as evidence. Step 1 showed the 9 Å ion is radial
  only to ≈ 2 ps (< t* = 2.67 ps), so no radial early window exists, and
  a radial model's full-window RMSE on a non-radial trace neither
  validates nor falsifies the law. It is recorded once, descriptively.
- **Known costs (thesis statements).** No held-out case remains;
  validation rests on Tier 2. 18 Å alone did not identify a power-law
  exponent (§10.8), so the family choice is limited to the three forms
  above by decision. The b ↔ E_bind trade-off on one curve is reported
  via the sensitivity half-widths. m_eff = 19 He is borrowed from the
  9 Å shell-loss curve (G3).
- **Projection on record (not a gate).** E_bind falls toward the old
  18 Å-only value (0.071 eV, §8); via the system slope ≈ −0.53 this
  would raise KE₁ by ≈ +0.02–0.03 eV.
- **Storage.** `data/reference/drag/18A/trajectory_matching_srho3p5A/`
  `{family}/` + `verdict.json`; the step-1 table and the 9 Å prediction
  beside it. Code requires `[PROCEED TO IMPLEMENTATION]`.

### G4.10 Force-balance cross-check of the 18 Å fit (agreed 2026-10-05)

- **Why.** G4.9 calibrates the drag on one curve (18 Å), with no held-out
  case. A local force balance on the smoothed TDDFT data is an independent
  check of b. It does not need a radial trajectory, so 9 Å can contribute
  again and widen the speed range (≈ 2.7–4.5 Å/ps). This revives the old
  Method A (`Drag_function/drag_function/drag_calculation.py`) with one
  correction.
- **Correction: project the Coulomb force onto the motion.** The old
  scalar balance `F_drag = F_C − m·d|v|/dt` assumes F_C ∥ v. Only
  F_C·cos θ changes the speed (θ = angle between v₂ and r₂ − r₁). At 9 Å
  the scalar form books F_C(1 − cos θ) as drag: 15–65 % of the cubic
  drag (b = 2.23) at 2–5 ps, 60 % at 7–8 ps, 12–18 % after 9 ps
  (per-ps bin means; cos θ 0.23–0.83, measured from the TDDFT
  positions and 0.5 ps-smoothed velocities). At 18 Å cos θ ≈ 1, so both
  forms agree there (consistency check).
- **Balance along v̂, per time point** [amu·Å/ps²]:

  $$\gamma(v)\,\hat\rho(r)\,v = F_C\cos\theta - \frac{dU}{dr}\,(\hat r\cdot\hat v) - m_\text{eff}\,\frac{d|v|}{dt}$$

  (*Amended 2026-10-05, before any run:* the well term was written
  `+F_well (r̂·v̂)`. `droplet_force` returns dU/dr > 0 and the force on
  the ion is −dU/dr·r̂, inward (`physics/leapfrog.py`), so the term
  enters with a minus sign. In the 18 Å interior it is ≈ 40 against
  F_C cos θ ≈ 68 amu·Å/ps², so the sign is material.)

  - |v| and d|v|/dt from the committed smoothed reference
    (`velocity_smoothed/cleaned_data_long.csv`, same pipeline as Method B);
    the derivative by Savitzky–Golay (primary window 0.5 ps).
  - F_C = k_e e²/R², with R and θ from the raw TDDFT positions; the
    velocity direction from the vector components smoothed over 0.5 ps.
  - F_well from the package `droplet_force` (no duplicate physics) at the
    G4.9 E_bind = 0.1214 eV, potential width 14.2 Å, the N = 2000 droplet
    radius of the Tier-0 setup, centred at the origin (assumption: the
    droplet centre stays at the I₂ centre of mass).
  - m_eff = 19 He (G3).
  - Points: each case's Method-B window, interior only (ρ̂ ≥ 0.9 at
    s_ρ = 3.5 Å), divided by ρ̂.
- **Fits.** F_drag(v) = a·v, b·v³, a·v + b·v³ (mass-agnostic convention),
  on 18 Å alone, 9 Å alone and both together.
- **Pre-registered agreement (frozen before any run).** The 18 Å cubic
  b_A agrees if it lies within the G4.9 b ± 2 half-widths
  (2.2286 ± 0.426 → [1.80, 2.65] amu·ps/Å²). Agreement → adopt
  `srho3p5A18A_cubic` for G5. Disagreement → no automatic adoption; open
  as a decision.
- **Recorded, not gating.** The 9 Å and joint b and best form; the scalar
  (old) vs projected b for both cases; sensitivity of b to SG window
  (0.25 / 0.5 / 1.0 ps), to E_bind ± its half-width, and for 9 Å to the
  anchored m(t) 21 → 14 He; the apparent drag at points outside the
  surface (ρ̂ < 0.1) as a direct measure of the late slowing seen in
  G4.9 step 1.
- **Amendments (2026-10-05, before any run; user-agreed).**
  1. Well-term sign corrected (see the balance above).
  2. *Derivative resolution.* The committed references are already
     order-1 SG-smoothed over 3.4 ps (18 Å) / 2.4 ps after IMF (9 Å),
     so their derivative is a window-mean slope and the 0.25 / 0.5 /
     1.0 ps SG sensitivity on them is near-inert. Added row (recorded,
     not gating): |v| and d|v|/dt from the **raw** TDDFT |v₂| with SG
     windows 0.25 / 0.5 / 1.0 ps. The committed reference stays primary.
  3. *18 Å form identifiability.* The 18 Å interior (ρ̂ ≥ 0.9) spans
     4.54–6.20 ps with |v| 2.95–3.02 Å/ps (2 % range): the 18 Å cubic
     b_A is effectively F_drag/v³ at v ≈ 2.98, and the 18 Å-alone form
     comparison cannot separate linear from cubic (stated before the
     run). Form information comes from 9 Å (|v| 3.08–4.94) and joint.
  4. *b_A ↔ E_bind.* In the 18 Å interior the drag is the residual of
     two comparable terms, so b_A is linear in E_bind
     (≈ −12.5 amu·ps/Å² per eV); the check tests the (b, E_bind) pair,
     not b alone. Added recorded sensitivity: droplet radius ± 8 % (the
     N → R convention spread, G1 above).
- **Cost.** No MD. A small module under `i2_helium_md/extraction/`, a
  driver and tests, ≈ 2–3 h. Code requires `[PROCEED TO
  IMPLEMENTATION]`.

- **Status — EXECUTED 2026-10-05 (measured; INCONCLUSIVE).** 18 Å
  cubic b_A = 1.618 fails the band as worded; 9 Å 2.013 and joint 1.991
  lie inside. The committed references (order-1 SG over 3.4 / 2.4 ps) are
  not derivative-grade: at 18 Å they sit 0.15–0.2 Å/ps below the
  oscillating raw speed in the 1.66 ps interior, so no force-balance
  coefficient is trusted yet. The first reading (depth-dependent drag)
  is withdrawn. Full record: D0 §1.G4.10, log 2026-10-05 G4.10 entry +
  correction, artifacts `data/reference/drag/force_balance_srho3p5A/`.

  > **NEXT.** Design and validate a derivative-grade smoothing for both
  > cases (discussion first), rerun the three forms on 9 Å, 18 Å and
  > joint, and compare form by form with trajectory matching. Adoption of
  > `srho3p5A18A_cubic` for G5 waits on it.

- **Post-hoc set (2026-10-05, after the run; user-agreed).** The June
  Method-A pipeline (short references, [4.54, 8.0] / [2.67, 6.0], MAD
  spline, 0.5 ps trim) with only the physics corrected (projection, well,
  F/ρ̂ at 3.5 Å): cubic b **18 Å 2.449, 9 Å 2.634 — both in the band**;
  9 Å (held out of G4.9) at the upper edge, leaving it under + 8 % radius
  (2.705) and the anchored mass (2.691); 18 Å robust except − 8 % radius.
  No form information (`linear_cubic` sign-unstable). The pre-registered
  gap came mostly from the long reference's end-effect (0.27 Å/ps at
  4.54 ps). Reading: the size of the cubic b is supported by an
  independent local method on both cases. D0 §1.G4.10-PH.

### G4.11 Cross-case drag-law test across families (PRE-REGISTERED 2026-10-05, user-checked; EXECUTED 2026-10-05 — see status at the end)

**Question.** Is there one friction law F(v) that describes the TDDFT
ion at both speed ranges — 18 Å (≈ 2.7–3.2 Å/ps) and 9 Å (≈ 3.4–4.6 Å/ps)?
Families: `linear` F = a·v, `quadratic` F = c·v², `linear_quadratic`
F = a·v + c·v², with `cubic` F = b·v³ as the control (G4.10-PH).
Units: a [amu/ps], c [amu/Å], b [amu·ps/Å²]; each term gives
amu·Å/ps² (force); γ = F/v [amu/ps] (friction convention).

**Two extraction routes (unchanged physics: Coulomb, well −dU/dr,
density gate s_ρ = 3.5 Å, U width 14.2 Å, m_eff = 19 He, R = 27.94 Å).**

- *Trajectory matching:* simulate the radial explosion with a candidate
  law, fit its coefficients + E_bind so the simulated |v(t)| matches the
  TDDFT curve (18 Å only, as G4.9).
- *Force balance:* no simulation; at each instant the drag is what
  remains of Newton's law along the motion,
  γρ̂v = F_C cos θ − dU/dr (r̂·v̂) − m d|v|/dt, then fit F(v). The
  G4.10-PH pipeline (June data handling: short references,
  [4.54, 8.0] / [2.67, 6.0] ps, MAD spline, 0.5 ps trim).

**Main criterion: cross-case consistency.** A single-case trajectory
fit matches its own curve for almost any law, and on 18 Å alone every
family reproduces the drag size at v ≈ 3 Å/ps. A law is tested only by
carrying it to the other speed range.

- **C-X (main, cross-case):** the family's 18 Å trajectory-matched law
  L₁₈(v), evaluated at the 9 Å force-balance points. Scale
  s_X = Σ F₉·L₁₈(v) / Σ L₁₈(v)² (least-squares factor; s_X = 1 ⇔ the
  18 Å law predicts the 9 Å drag). **Pass: |s_X − 1| ≤ 0.19.**
- **C-M (secondary, cross-method):** the same scale on the 18 Å
  force-balance points, s_M. **Pass: |s_M − 1| ≤ 0.19.**
- **C-A (recorded, single-coefficient families only):** force balance
  alone, 9 Å coefficient / 18 Å coefficient — cross-case with no
  trajectory fit involved.
- **Band 0.19:** the cubic's relative trajectory-fit uncertainty,
  2 · 0.2127 / 2.2286 (the G4.10 band), applied to every family alike.
- **Shape (recorded):** RMS of (F₉ − s_X·L₁₈)/⟨F₉⟩ — a two-coefficient
  law can match the 9 Å scale with the wrong slope.
- **E_bind pairing:** the force balance for family f uses f's own
  18 Å trajectory-fit E_bind (b-given-E_bind, as G4.10).

**Decision rule.** A family is *consistent* if C-X and C-M pass.
Only the cubic consistent → the cubic form is supported by the TDDFT
data (not only chosen). Several consistent → the TDDFT data do not
discriminate those forms; the cubic stays by decision (continuity,
production cap). Cubic not consistent → open as a decision.

**Known before writing this (disclosed, not hidden).** Cubic control:
s_X = 2.634 / 2.2286 = 1.182 (passes, near the edge), s_M = 1.099.
Linear, from the G4.10-PH rows (force balance at the *cubic's* E_bind):
9 Å a 41.97 vs 18 Å trajectory-fit a 23.15 → s_X ≈ 1.81. Quadratic and
linear + quadratic at 3.5 Å: not known.

**Step 1 — 18 Å trajectory-matching fits at s_ρ = 3.5 Å (new for the
quadratic families).**

- Reuse G4.9 `linear` (a 23.15, E_bind 0.1165) and `cubic` (b 2.2286,
  E_bind 0.1214). New: `quadratic` (`lq_shared_pure_quadratic`) and
  `linear_quadratic` (`lq_shared_3param`), 18 Å only, same engine, N = 50,
  seed 20260604, 20 ps at 0.01 ps, §10.4.1 locked anchors (a0 14.5556,
  c0 11.0161), sensitivity half-widths as G4.9.
- Bands as G4.9: 18 Å RMSE ≤ 0.19 Å/ps, escape fraction = 1 (a family
  failing them is reported, still scored).
- **Oracle first:** (i) the G4.9 cubic point reproduces its committed
  18 Å RMSE 0.1092 through the new driver; (ii) the June per-case
  quadratic point (§10.8: c 11.27, E_bind 0.060) reproduces its 18 Å
  RMSE 0.0993 at the 14.2 Å gate.
- **Code:** a G4.11 driver with the gate pass-through, writing to
  `data/reference/drag/18A/trajectory_matching_srho3p5A/{quadratic,linear_quadratic}/`
  (never the 14.2 Å folders); the per-case form writer stamps the gate
  (`drag_spatial_gate`, `drag_gate_steepness_angstrom`,
  `potential_steepness_angstrom`), as G4.3 change 4. Diagnostic bundles,
  not preset-wired.

**Step 2 — force balance + scoring.** `fit_drag_form` gains `quadratic`
and `linear_quadratic`; a G4.11 driver evaluates C-X / C-M / C-A per
family on the G4.10-PH pipeline, with the G4.10-PH sensitivity rows
(ρ̂ cut, E_bind ± hw, radius ± 8 %, spline s, 9 Å anchored m(t)) —
recorded; the verdict is on the primary row. Artifacts
`data/reference/drag/force_balance_srho3p5A/g411_cross_case/`
(`verdict.json`, `scores.csv`).

**Not in scope.** 9 Å trajectory fits (non-radial, G4.9); power law
(18 Å alone rails n, §10.8); production changes.

**Tests.** Writer gate stamp; single-case lq fit under the 3.5 Å gate
override (stub landscape, no MD); `fit_drag_form` exact recovery for the
two new forms; s_X / s_M on synthetic forces (exact law → 1).

**Cost.** Step 1 ≈ two single-case fits (G4.9 scale, < 1 h incl.
sensitivity); step 2 zero MD. Code ≈ 2 h. Requires
`[PROCEED TO IMPLEMENTATION]`. Frozen: criteria, band, decision rule,
E_bind pairing, pipeline (this commit, before any run).

- **Status — EXECUTED 2026-10-05 (measured; rule outcome: cubic
  supported).** Step 1: `quadratic` c 7.224, E_bind 0.1191, 18 Å RMSE
  0.1168 (bands pass); `linear_quadratic` collapses onto it (a 0.004).
  Step 2, primary row, s_X / s_M: linear 1.823 ✗ / 0.981, quadratic and
  linear + quadratic 1.476 ✗ / 1.035, **cubic 1.182 ✓ / 1.099 ✓**. Only
  the cubic is consistent. Marginal: the cubic C-X leaves the band under
  + 8 % radius (1.214) and the 9 Å anchored m(t) (1.207); linear and
  quadratic fail in every row. Oracles bit-exact. Full record: D0
  §1.G4.11, log 2026-10-05 G4.11 entry, artifacts
  `data/reference/drag/force_balance_srho3p5A/g411_cross_case/`.

  > **NEXT (user).** Adopt `srho3p5A18A_cubic` for G5 (gate cleared on
  > the rule), or not.

- **Post-hoc caveat (2026-10-05, user's Tier-0 9 Å runs; visual, not yet
  quantified).** Run radially on 9 Å with the 18 Å-only laws, the
  quadratic tracks the less-transverse HeDFT atom (atom 1) and the cubic
  tracks the transverse, stalling atom 2. A radial ion has no transverse
  loss path, so with the right law it should leave at least as fast as
  atom 1: the cubic then over-dissipates at 9 Å speeds. This contradicts
  the C-X reading ("every law under-predicts the 9 Å drag"), which rests
  on the force balance of **atom 2** alone. The G4.11 rule outcome stands
  as run; its physical reading is contested and adoption is **on hold**.
  Discussion: log 2026-10-05 "G4.11 re-read". Follow-up: G4.12.

### G4.12 TDDFT-tied linear + quadratic ridge vs the experiment, uncapped (user-checked 2026-10-05; open items (a)–(c) settled below; frozen before any run)

**Why.** 18 Å fixes one number, the drag at v ≈ 3 Å/ps (§10.8; G4.10
amendment 3). Linear (a 23.15, RMSE 0.1239) and quadratic (c 7.224,
RMSE 0.1168) both fit it; G4.11's free lq fit fell onto the quadratic
(a 0.004). So a one-parameter ridge of mixed laws F = a·v + c·v² should
fit 18 Å. Along it the fast-class force varies widely (projection, F(3)
≈ 67 amu·Å/ps² held fixed): F(9.7) ≈ 217 + 65·c, i.e. Φ = F(9.7)/418.5
≈ 0.52 (c = 0) → 1.00 (c ≈ 3.1) → 1.62 (c = 7.22). The window
c ≈ 0.8–3 amu/Å (a ≈ 13–20 amu/ps) puts the n = 1 force inside the band
where linear-family systems landed in MD (Φ 0.64–0.99), **without a cap**.
The free-form arm 2 "rejects the ram term" result does not exclude this:
it added c·v² at fixed a = 35 (raising the force everywhere), not along
the 18 Å ridge. The free-form family closed 2026-08-13 with a free a; here
a and c are tied to TDDFT, so this is a different question.

**Form and units.** F = g(ρ̂)·(a·v + c·v²), `linear_quadratic`, no cap;
a [amu/ps], c [amu/Å]; each term amu·Å/ps²; γ = a + c·v [amu/ps]
(friction convention; mass-agnostic). s_ρ = 3.5 Å, U width 14.2 Å,
m_eff = 19 He.

**What is honest about it.** 18 Å permits the ridge but does not select a
point on it; c is chosen by the experiment, as v_c is today. The free
element moves from a cap above 5.5 Å/ps to the curve's shape over all
speeds — kink-free, not more physical. Known risk (projection from the
free-form ring): along the ridge F(5.5) ≈ 120–165 vs h405 418.5, the
mid-band softness that overheated midHot (1.4–2.2) at a = 27.5.

**Step 1 — map the 18 Å ridge (Tier-0 engine, 18 Å only).** For fixed
c ∈ {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.224} fit
(a ≥ 0, E_bind) on 18 Å: same engine, N = 50, seed 20260604, 20 ps at
0.01 ps, long reference, 3.5 Å gate. Record RMSE, a, E_bind, escape,
F(3), F(5.5), F(9.7), Φ per point. Bands (G4.9): RMSE ≤ 0.19 Å/ps,
escape = 1. The grid step (0.4–0.5) is finer than the band that the
landing selects (Φ 0.64–1.0 ⇔ Δc ≈ 2.3).
- *Oracles first:* (i) the lq form at the G4.9 linear point (a 23.154,
  c 0, E_bind 0.11647) reproduces RMSE 0.12389 bit-exactly; (ii) at the
  G4.11 quadratic point (a 0, c 7.2241, E_bind 0.11910) RMSE 0.11680
  bit-exactly; (iii) the fixed-c = 0 fit returns the G4.9 linear point
  within optimizer tolerance (|Δa| ≤ 0.2, |ΔE| ≤ 0.002 eV).

**Step 2 — 9 Å exit bound (Tier-0 radial runs, recorded per point).**
- *Anchor first:* score the user's existing 9 Å runs of the 18 Å-only
  linear / quadratic / cubic bundles and confirm the visual reading in
  numbers (quadratic ≈ atom 1, cubic ≈ atom 2).
- *Atom assignment:* the "less-transverse" atom is the one with the
  smaller mean |v_⊥|/|v| over the 9 Å window, from the raw HeDFT
  velocities (expected atom 1; verified, not assumed).
- *Exit speed:* mean |v| over [13.0, 14.0] ps, MD ensemble mean and each
  HeDFT atom. If either HeDFT atom is still inside the surface there
  (ρ̂ > 0.1), stop and return to the user before scoring.
- *Bound:* pass if v_exit,MD ≥ v_exit,atom 1. Margins are recorded;
  |margin| ≤ 0.05 Å/ps is flagged "at the bound". Also recorded:
  in-window RMSE against |v1|, |v2| and their mean.
- *One-sided (disclosed).* The bound rejects over-dissipation only; a law
  that is too weak also passes. Along the ridge F(4) ≈ 89.3 + 4c, so the
  9 Å drag falls as c drops below 7.22: if the quadratic sits at atom 1,
  every c < 7.22 passes by construction. Step 2 is then a consistency
  check on the 9 Å reading, not a selector inside the window; the
  selection is step 3's.
- *Mass bias (disclosed).* The MD ion stays at m_eff = 19 He; HeDFT
  atom 1 sheds ≈ 21 → 14 He over the window. A lighter ion decelerates
  more slowly under the same drag, so the comparison is biased against
  the MD (makes the bound harder to pass), not in its favour.

**Step 3 — experiment, zero MD (the twin), uncapped.** Points passing
steps 1 and 2 enter `scripts/tier2_h2b_forward_model.py` (the `linq`
branch, no cap) at the h405 pipeline settings. Hard gate = §3.5c verbatim
(n₁_solv ∈ [0.19, 0.30] AND n̄ ∈ [4.4, 7.1]). Reported, not gating:
KE₁ (licensed, Step-0 ρ 0.998), Φ, midHot, deepKE (blind-spot
disclosure), W₁, trap floor. **Settled (user, 2026-10-05):**
(a) gate: ridge points run with `rho_steepness = 3.5` (the width they
were fitted under; potential stays 14.2 Å). Anchor first: h405 through
the twin at both 14.2 and 3.5 Å, so the gate's own shift is measured
before any ridge point is scored. (b) E_bind: each ridge point's own
18 Å E_bind is primary (TDDFT pairing); the twin's well grid is a
recorded sensitivity only. (c) τ / E₀: linscan's full grid, not one
narrowed around h405 (a cubic-centred grid could fake "no point lands");
any near-pass is refined before a fail is recorded. Execution stops
after step 2 for a user recap before step 3.

**Decision rule.** Some c > 0 passes steps 1–3 → an uncapped
TDDFT-tied mixed law is a landing candidate; an MD battery follows under
its own trigger (nothing adopted from the twin). No point passes →
recorded: the TDDFT-tied ridge cannot land uncapped, and the cap stays a
permanent effective element. Step 2 failures are reported per point;
whether the bound gates step 3 is fixed here: **it does** (a law that
over-dissipates the radial 9 Å ion is not carried forward).

**Not in scope.** Power law; capped variants of the ridge; production
changes; any re-reading of G4.11's frozen criteria.

**Code (requires `[PROCEED TO IMPLEMENTATION]`).** A fixed-c single-case
lq fit (a, E_bind free) in the form engine or as a driver objective; a
G4.12 driver (oracle / ridge / bound modes); a 9 Å exit-bound scorer;
twin stage wiring. Tests: fixed-c fit on a stub landscape; the c = 0
seam with the linear form; exit-bound arithmetic on synthetic traces.

**Cost.** Step 1 ≈ 12 single-case fits (≈ 1 h, at most two processes);
step 2 ≈ 12 × N = 50 9 Å runs (minutes each); step 3 zero MD.

- **Status — steps 1–2 EXECUTED 2026-10-05 (measured).** Oracles pass
  (bit-exact; c = 0 seam exact). Step 1: all 12 points pass the bands,
  RMSE 0.117–0.124; window c 0.8–3 → Φ 0.66–0.99. Step 2: atom 1 is the
  less-transverse atom; margin vs atom 1 +1.27 (c 0) → +0.20 (c 5) →
  −0.14 (c 7.224, quadratic): 11 of 12 pass. Caveat: the window passes
  by +0.56 to +1.06 Å/ps, i.e. it under-dissipates relative to atom 1
  (one-sided bound). Full record: D0 §1.G4.12, log 2026-10-05 G4.12 entry.

- **Status — step 3 EXECUTED 2026-10-05 (twin, measured).** Anchors
  pass (h405 at 14.2 Å string-identical). After the frozen near-pass
  refinement all 11 points gate → **rule outcome: landing candidate**
  (MD battery only under its own trigger). Reads: one KE₁ / midHot dial
  (c 0 → 5: KE₁ 0.95 → 0.31, midHot 2.35 → 0.53); soft band at c ≈ 2–2.5
  with KE₁ 0.62–0.55 vs h405 0.725 at midHot 1.06 — h405 dominates at
  matched midHot. Record: D0 §1.G4.12, log G4.12 step-3 entry.

  > **NEXT (user, 2026-10-06).** MD battery at c = 2.0 → G4.13.

### G4.13 MD battery for the G4.12 ridge point c = 2.0, uncapped (PRE-REGISTERED 2026-10-06, user trigger "extensive MD battery for c = 2.0 around a suitable parameter grid"; frozen before any MD)

**Why c = 2.0.** The one ridge point whose twin midHot sits in the soft
band (1.11) on the refined surface (G4.12 step 3). The law:
F = g(ρ̂)·(a·v + c·v²), a 16.717 amu/ps, c 2.0 amu/Å, own 18 Å E_bind
0.11748 eV (`g412_ridge/c2p000.json`, values read at build time, never
retyped), density gate `erf_independent` 3.5 Å, no cap. Units as G4.12.

**Pipeline = production** (`gen_tier2_detfix_battery.py`): h405's frozen
cell (`build_cell(h405)`), Stage I to 500 ps, E2 off, partner-aware
closure, N = 1000 per member, max two concurrent. Only the drag law, the
well and τ / E₀ move. Bundle `srho3p5A18A_linear_quadratic` with
{a, c, E_bind} overridden (its stamped gate 3.5 Å comes with it).

**Stage A — grid ring (seed 20260730 = h405 detfix s1, CRN).** Four τ
arms × an E₀ ladder of four, placed on the twin's gated ridge with the
measured twin→MD transfer (Block-0 n̄ model a 0.581 / b −0.238, resid
sd 0.29; MD n₁ ≈ 0.01–0.015 below twin → ladders lean to higher E₀):

| τ [ps] | E₀ [eV] |
|---|---|
| 4.8 | 0.39, 0.40, 0.41, 0.42 |
| 5.2 | 0.36, 0.37, 0.38, 0.39 |
| 5.6 | 0.335, 0.345, 0.355, 0.365 |
| 6.4 | 0.30, 0.31, 0.32, 0.33 |

Plus **h405g35**: h405's own cell with the 3.5 Å density gate (the only
change), same seed — the MD gate shift on the incumbent. 17 runs. The
14.2 Å partner is the committed `detfix_battery_table.csv` row `s1_new`
(anchored: its run dir rescored must reproduce it, 4 dp).

**MD hard gate (unchanged, §3.5d):** n₁_solv ∈ [0.19, 0.30] AND
n̄ ∈ [3.77, 4.37]. Reported, not gating: KE₁ (mean/median), KE₂, midHot,
deepKE, W₁, χ²_med, trap (bound / marginal), S — full vector (scorer
CSVs carry every column).

**Stage B — battery (selection frozen here).** Among Stage-A gated
cells: midHot ∈ [0.85, 1.15] and maximal KE₁_mean; none in band → the
gated cell with midHot nearest the band; ties → lower W₁. No gated cell
→ stop, record, no Stage B. The selected cell runs at the remaining
production seeds (s2–s5 20260731–34, r6 20260729) → pooled N = 6000
with the Stage-A s1 member, scored against the production reference row
`pooled_new6000` (nbar 3.890, n₁ 0.2093, W₁ 0.7551, midHot 0.948,
deepKE 0.571, χ² 456, KE₁ 0.640).

**Decision rule (Stage B pooled).** *Matches or beats h405* ⇔ pooled
gate passes AND KE₁_mean ≥ 0.640 AND midHot ∈ [0.85, 1.15]. Otherwise
*dominated*: recorded, the cap stays as the effective element; nothing
is adopted in either case (adoption is a separate user decision).
**Twin prediction (disclosed):** at the ridge cells twin KE₁ 0.53–0.60
(twin ~1–2 % cold) → MD KE₁ ≈ 0.54–0.61 → *dominated*.

**Oracles before any MD.** (i) the c = 2.0 coefficients and E_bind
read from the committed ridge JSON; (ii) each cell's cfg diffs against
the h405 detfix member of the same seed in exactly the pre-registered
key set (drag law, well, gate, τ, E₀ — set fixed at the dry run, any
other key aborts); (iii) handover alignment as the production
generator; (iv) scorer anchor: h405 detfix s1 rescored = committed
`s1_new` row (4 dp).

**Not in scope.** Other c; capped variants; 18 Å / Tier-0 changes;
production changes.

**Cost.** ≈ 3 min per N = 1000 member (detfix precedent) → Stage A
≈ 30 min, Stage B ≈ 10 min. Code: generator + scorer (+ tests).

- **Status — EXECUTED 2026-10-06 (measured): decision *dominated* at
  the selected cell t64e310 only.** Stage A: 7/16 cells gate (all
  midHot in band); selected t64e310. Stage B pooled N = 6000: n₁ 0.186
  (gate ✗), KE₁ 0.638 (h405 0.640), midHot 1.032, W₁ 1.100 (h405
  0.755), χ² 558 (456). Record: D0 §1.G4.13, log 2026-10-06 G4.13 entry.
- **Correction (2026-10-06, user-caught).** The verdict does not extend
  to the law. Stage A selected on one seed with E₀ step 0.01 and τ step
  0.4–0.8, while the seed sd of n₁ (0.010, Stage B) equals the whole
  gate margin and the gated E₀ band moves −0.07 eV per ps of τ: the grid
  could not resolve the landing. "The cap stays" is withdrawn as a
  conclusion; the law-level question moves to G4.14.

### G4.14 Fine τ-line sweep for ridge-law landings (FROZEN 2026-10-06, user-agreed; nothing run; code requires `[PROCEED TO IMPLEMENTATION]`)

**Why.** G4.13 resolved neither the gate (n₁ seed noise ≈ margin) nor τ,
the dominant histogram lever (W₁ jumps +0.13–0.17 between τ 5.6 and 6.4
on seed s1, unresolved; h405's own τ 4.4 never sampled). Measured
transfer for the c = 2.0 law (16 Stage-A cells vs twin at the same cell;
seed sd from the 6 Stage-B seeds):

| | MD − twin | seed sd (N = 1000) | placeable from the twin? |
|---|---|---|---|
| KE₁ | +0.020 ± 0.001 | 0.002 | yes |
| midHot | × 0.945 ± 0.014 | 0.003 | yes |
| n̄ | −0.37 ± 0.035 | 0.07 | yes (E₀: −0.25 per 0.01 eV) |
| n₁ | −0.011 ± 0.014 | 0.010 | no — decides the gate |
| W₁ | +0.31 ± 0.06 | 0.036 | no — set by τ, MD only |

(Committed: `data/runs/h2b_forward_model/g413_transfer.csv`, from
`tier2_g413_table.py --transfer`; anchored on `atlas_g412_twin.csv`.)

**Design inputs settled in discussion (2026-10-06; read from committed
CSVs — measured unless marked).**
- *s1 is a high seed for this law.* t64e310: s1 n₁ 0.204, the other five
  0.177–0.190 (mean 0.183, sd 0.005). Every transfer offset above is
  s1-only, and CRN makes the whole Stage-A surface share s1's draw: grid
  averaging under one shared seed does **not** remove it. (h405: s1 only
  +0.0035 above its pool — the bias is law-specific.) *Projection:* if
  the s1 bias is surface-wide, MD − twin n₁ ≈ −0.03; the twin's n₁
  ceiling inside the MD n̄ gate is 0.20–0.22 at every c 0.8–2.5,
  τ 4–6.4 (`atlas_g412_twin/refine.csv`), so the MD ceiling would sit at
  ≈ 0.17–0.19 — possibly no landing with margin at any τ.
- *The n₁ offset is structured:* within each τ, d(n₁) runs ≈ 0 → −0.03
  as E₀ rises (MD n₁ saturates; τ 4.8: 0.193 / 0.193 / 0.181 / 0.188).
- *W₁ is not monotone in τ:* gated MD W₁ 0.87–0.93 over τ 4.8–5.6, 1.10
  at 6.4; the twin minimum is at τ 5.2–5.6 and rises again at 4.4.
- *Like-for-like partner:* the 3.5 Å gate alone moves h405's W₁ 0.80 →
  0.866 (s1, `h405g35`) — the ridge's plateau. G2 moves production to
  3.5 Å, so pooled h405g35 is a required reference.
- *KE₁ is n̄-dependent:* at τ 6.4, dKE₁/dn̄ ≈ 0.11; the t64e310 "tie"
  (0.638) sat at n̄ 4.085 vs h405 3.890. *Projection:* ≈ 0.62 at h405's
  n̄. KE₁ is therefore ranked at matched n̄.
- *midHot bounds c:* twin midHot × 0.945 ≈ 1.23 (c 1.6), 1.07 (c 2.0),
  0.92 (c 2.5).

**Design (τ is the primary axis; E₀ is slaved to target n̄ along it).**

0. **Step 0a — four fine ridge points (Tier-0, zero MD).**
   c ∈ {1.8, 1.9, 2.1, 2.2} (user, 2026-10-06: a finer set around 2.0
   to check whether anything better comes up) through G4.12 steps 1–2
   unchanged: `scripts/extraction/method_b_g412_ridge.py`, fixed-c 18 Å
   fit (a ≥ 0, E_bind free; N = 50, seed 20260604, 20 ps at 0.01 ps,
   3.5 Å gate), bands RMSE ≤ 0.19 Å/ps and escape = 1, then the 9 Å exit
   bound (gates, as in G4.12). Output `g412_ridge/c1p800.json` …
   `c2p200.json` in the committed format. *Oracle first:* the c = 2.0 fit
   rerun reproduces the committed `c2p000.json` (a, E_bind, RMSE)
   bit-exactly.
1. **Step 0 — calibrated twin (zero MD).** c ∈ {1.2, 1.6, 1.8, 1.9, 2.0,
   2.1, 2.2, 2.5} (cached chords + the four new; τ is an exact rescale
   K → K·6.55/τ) × τ 3.6–7.2 step 0.1 × E₀ step 0.0025, uncapped, own
   E_bind, `rho_steepness` 3.5. Output per c: E₀ at predicted MD n̄ 3.85
   and 4.05 along τ (n̄ offset −0.37), and predicted MD KE₁ (+0.020),
   midHot (× 0.945), n₁ (twin − 0.011, flagged *gate-undecided*).
   *Oracles:* twin rows at shared grid points reproduce
   `atlas_g412_twin.csv` / `atlas_g412_refine.csv`; the offsets
   recomputed reproduce `g413_transfer.csv`.
   **midHot filter (frozen):** a c enters MD only if its predicted MD
   midHot is in [0.85, 1.15] at some τ in 4.8–6.6 on its n̄ targets.
   Offsets are measured at c = 2.0 only; Step 1 measures them per c.
   **Stop after Step 0 for a user recap** (filter outcome, MD run list).
2. **Reference — pooled h405g35 (MD, 5 runs).** h405's cell with the
   3.5 Å density gate (the G4.13 `h405g35` definition) at s2–s5
   (20260731–34) and r6 (20260729), pooled with the committed s1 member →
   N = 6000. *Anchor:* the s1 run dir rescored reproduces its committed
   `g413_table.csv` row (4 dp). Runs before Step 1.
3. **Step 1 — MD τ lines (fitted, independent seeds).** Production
   pipeline exactly as G4.13 (`build_cell(h405)`, Stage I 500 ps, E2 off,
   partner-aware closure, N = 1000 per run, ≤ 2 concurrent; only drag
   law, well, gate, τ, E₀, seed move).
   - c = 2.0: τ ∈ {4.0, 4.4} ∪ {4.8, 5.0, …, 6.6} (12 τ) × two E₀ (Step-0
     targets n̄ 3.85 / 4.05) → 24 runs.
   - Every other c passing the filter: τ ∈ {4.8, 5.2, 5.6, 6.0, 6.4} ×
     the same two targets → 10 runs each.
   - **Seeds:** one fresh seed per run, consecutive from 20261101 in
     generator order (disjoint from production 20260729–34), recorded per
     run. No CRN: the noise is independent across cells, so it averages
     in the fit, s1's bias stays out, and the Step-2 production seeds
     stay unused for selection (no winner's curse).
   - The 16 Stage-A s1 cells are an out-of-sample check only (the c = 2.0
     fit's mean residual on them estimates s1's offset); never fit data.
4. **Fit model (frozen).** Per c and observable (n₁, n̄, KE₁, W₁, midHot,
   χ², deepKE, trap), least squares in τ and ΔE = E₀ − E₀,3.85(τ): cubic
   in τ (quadratic for 10-run c) + ΔE + τ·ΔE; n₁ additionally + ΔE²
   (the saturation above). **Residual criterion:** residual sd ≤ 1.5 ×
   per-member seed sd (n₁ 0.015, n̄ 0.105, KE₁ 0.003, W₁ 0.054, midHot
   0.0045). Fail with residuals showing no τ structure → a second
   independent-seed set at the same cells (refit on both); fail with τ
   structure → one more τ order, once; still failing → recorded
   *unresolved* for that c, no Step 2.
5. **Step 2 — confirmation (production seeds).** Per c, on the fitted
   surface, the feasible set is: predicted n₁ ≥ 0.200, n̄ ∈ [3.83, 4.31],
   midHot ∈ [0.85, 1.15]. If it is non-empty, two points: τ_KE (max KE₁
   at n̄ 3.89) and τ_W (min W₁), merged if ≤ 0.2 ps apart; E₀ placed by
   the fit at n̄ 3.89. Each runs at s1–s5 + r6 → pooled N = 6000. Empty
   feasible set → recorded "no landing at c", no Step 2.

**Landing (frozen, pooled N = 6000).** n₁ ≥ 0.200 (gate 0.19 + one
per-member seed sd ≈ 2.5 pooled SE) AND n̄ ∈ [3.83, 4.31].

**Ranking (frozen; separate from landing), against both references
reported separately:** `pooled_new6000` (14.2 Å; n̄ 3.890, KE₁ 0.640,
W₁ 0.755) and pooled h405g35 (3.5 Å). Axes: KE₁ at n̄ 3.89 (pooled KE₁
shifted by the fitted local dKE₁/dn̄ — disclosed), W₁, midHot in band.
*Beats the reference* ⇔ landing AND KE₁ ≥ ref AND W₁ ≤ ref + pooled SE
(0.036/√6 ≈ 0.015) AND midHot ∈ [0.85, 1.15]. Otherwise the fitted KE₁
vs W₁ trade-off curve along τ per c is the recorded result. deepKE, χ²,
trap: reported, not ranked. Nothing adopted in either case (adoption is
a separate user decision); production unchanged.

**Prediction (frozen before MD; projections).** MD n₁'s gate ceiling
rises with τ; if the s1 bias is surface-wide, c = 2.0 has no landing
with margin below τ ≈ 6. KE₁ at n̄ 3.89 ≈ 0.62 at τ 6.4. W₁ is minimal
(≈ 0.87) near τ 5–5.6, not at low τ. Against h405g35 the ridge ties W₁
near τ 5 but loses ≈ 0.05 KE₁; against `pooled_new6000` it is dominated.
The fine c points move along the same KE₁ / midHot dial: c ≤ 1.8 gains
KE₁ but leaves the midHot band at its upper edge (≈ 1.15).

**Oracles before any MD.** (i) Ridge coefficients and E_bind read from
the committed JSONs, never retyped; (ii) each cfg diffs against the h405
detfix member in exactly {drag law, a, c, E_bind, gate, τ, E₀, seed} —
any other key aborts; (iii) handover alignment as the production
generator; (iv) scorer anchor: h405 s1 rescored = committed `s1_new`
(4 dp).

**Execution order.** Step 0a → Step 0 → user recap → h405g35 pool →
Step 1 → fit + residual check (+ second seed set if triggered) → Step 2
→ D0 §1.G4.14 first, then log and status.

**Not in scope.** c outside 1.2–2.5 in MD; capped ridge variants; 18 Å /
Tier-0 changes beyond the four Step-0a fits; production changes.

**Code (requires `[PROCEED TO IMPLEMENTATION]`).** Step-0a reuse of
`method_b_g412_ridge.py` for the new c (+ chords); Step-0 twin sweep
with offsets; the G4.14 generator (seed block, n̄-targeted E₀, h405g35
pool) on `gen_tier2_g413_ridge_battery.py`; the surface fitter +
residual criterion + selection on `tier2_g413_table.py` (scorer CSVs
carry the full vector). Tests: seed-block assignment and disjointness;
n̄-target E₀ inversion on a stub twin; the fitter recovering a synthetic
surface within noise, and the residual criterion firing on a stub with
a missing τ term; the feasible-set selection on a synthetic surface.

**Cost.** Step 0a ≈ 4 fits + 4 9 Å bounds (≈ 30 min); Step 0 zero MD;
h405g35 5 runs (≈ 15 min); Step 1 24 + 10 per further c (≈ 3 min per run,
two concurrent → ≈ 45 min at c = 2.0); Step 2 6 runs per point
(≈ 10 min each).

- **Status — Steps 0a, 0, h405g35 pool, Step 1 EXECUTED 2026-10-06
  (measured).** Step 0a: c 1.8/1.9/2.1/2.2 pass (RMSE 0.1219–0.1222, 9 Å
  margins +0.82 → +0.73); oracle bit-exact. Step 0: anchors pass; midHot
  filter keeps c 1.8–2.5, drops 1.2/1.6 → 74 Step-1 runs. h405g35 pooled
  N = 6000: n₁ 0.207, n̄ 3.857, KE₁ 0.643, W₁ 0.818, χ² 636. Step 1: all
  74 runs; no c beats h405 on the frozen ranking (best W₁ 0.831 at c 1.8
  τ 4.9 with KE₁ 0.577; best KE₁ 0.632 at W₁ 1.05). Record: D0 §1.G4.14.

**Amendments (2026-10-06, user-agreed after Step 1; frozen before the
confirmation MD).**
- *A1 — n₁ fit model.* The extra ΔE² term is dropped: with two E₀ per τ
  it is not identifiable (design condition 118–294 vs 6–9 without it;
  residual sd unchanged) and swung n₁ between the levels (c 1.9: 0.117 →
  0.250; c 2.2: −0.22 → 0.28). All c refitted without it.
- *A2 — c = 2.2.* Trips the second-seed rule (W₁ resid 0.0576 > 0.0545,
  no τ structure); recorded as triggered and **skipped** (dominated on
  both axes).
- *A3 — histogram bands* (reported, not gating): per-band solvated mass
  and W₁ share for n 1 / 2–3 / 4–6 / 7–9 / 10–14 / ≥ 15 in every scorer row.
- *A4 — disk.* Step-1 runs drop `ion.npz` after `detection.npz` verifies
  (T: full; the scorer reads cfg + detection only).
- *A5 — KE-led criterion + one pooled confirmation.* Reason (discussion
  2026-10-06): n₁ / n₁/n₂ are de-emphasised; 93/95 ridge runs are already
  monotone within noise, and the high-τ ridge cells match h405's per-member
  KE. **Criterion (pooled N = 6000):** n̄ ∈ [3.83, 4.31]; solvated
  histogram monotone within 2σ Poisson over n = 1…17; midHot ∈ [0.85, 1.15].
  **Ranked on KE** against each pooled reference separately (h405 14.2 Å
  `pooled_new6000`: χ² 456, KE₁ 0.640; h405g35 3.5 Å: χ² 636, KE₁ 0.643):
  *matches or beats on KE* ⇔ criterion passes AND χ²_med ≤ ref AND
  KE₁ ≥ ref − 0.002. W₁, bands, n₁ gate reported, not ranked. **Cell:**
  c = 2.0, τ 6.4, E₀ 0.3155 eV (fitted n̄ = 3.89 point of the committed
  trade-off curve, read at build time), production seeds s1–s5 + r6.
  **Prediction (fitted, frozen):** KE₁ ≈ 0.616, midHot ≈ 0.99, n₁ ≈ 0.203,
  W₁ ≈ 1.10; single-run-level χ² ≈ 281 → pooled ≈ 420–470 (h405's pooled /
  member-mean ratio 1.5–1.66). Expected verdict: **loses on KE** against
  both references through KE₁ (−0.024 / −0.027), while χ² is expected to
  beat the 3.5 Å reference and roughly tie the 14.2 Å one. Nothing is
  adopted in any case.

- **Status — A5 EXECUTED 2026-10-06 (measured): loses on KE** (as
  predicted). Pooled c 2.0 / τ 6.4: n̄ 3.886, n₁ 0.204, monotone, midHot
  0.995, KE₁ 0.617, χ² 507, W₁ 1.097. Record: D0 §1.G4.14; pooled figures
  container `…_N6000_detfix_conf270_g414c200t64pooled_th500`.

**Recap for the thesis-choice discussion (2026-10-06; entry point for the
next chat).** The question: what the uncapped linear + quadratic (lq)
ridge achieved against the capped production law h405, its strengths and
weaknesses, and how the thesis presents it. Nothing is decided yet.

*KE side (measured, pooled N = 6000).*

| | lq c 2.0 τ 6.4 | h405 3.5 Å | h405 14.2 Å | experiment |
|---|---|---|---|---|
| midHot (n 2–8 KE, sim/exp) | 0.995 | 0.981 | 0.948 | 1 |
| χ²_med (KE curve) | 507 | 636 | 456 | — |
| KE₁ [eV] | 0.617 | 0.643 | 0.640 | 1.302 ± 0.071 |
| KE₂ [eV] | 0.526 | 0.557 | 0.550 | 0.706 |
| W₁ (size histogram) | 1.097 | 0.818 | 0.755 | 0 |

*Strengths (argued in discussion).* (1) Standard drag decomposition:
viscous a·v + ram c·v², units balanced, friction convention unchanged.
(2) TDDFT-tied: a and E_bind fitted on 18 Å for every c (RMSE
0.117–0.124); c ≤ 5 passes the 9 Å exit bound. (3) Retires the effective
v_c cap (D0 §17 class **E**): smooth, kink-free at every speed. (4) The
KE-best cell needs τ ≈ 6.4 ps, essentially the sourced GAH25 6.55 ps
(h405: 4.4 ps). (5) KE: n = 2–8 matched (midHot 0.995); χ² better than
h405 at the same 3.5 Å gate (production moves to 3.5 Å under G2). The
KE₁ gap (−0.026 eV) is ≈ 2 % of the experimental value; every law sits at
about half of experiment — the flight-mass question (M1), not drag.

*Weaknesses / caveats.* (1) c is not fixed by TDDFT (18 Å pins only
F(v ≈ 3 Å/ps)); experiment chooses c via midHot (window ≈ 1.8–2.5): one
free element (v_c) traded for another (c), with a physical shape. (2)
Size histogram worse (W₁ 1.10 vs 0.82): good KE needs high τ, which
empties the n ≥ 10 tail; the common misfit (no n ≥ 17, n₁ low) is a
mass-model issue shared with h405, but h405 compromises better. (3) The
9 Å bound is one-sided: the lq law dissipates less than HeDFT atom 1
(margin ≈ +0.7 Å/ps). (4) Landing marginal (n₁ 0.204 vs 0.200).

*Candidate thesis claim (draft, not agreed).* "An uncapped, TDDFT-tied
linear + quadratic drag law reproduces the measured I⁺(He)ₙ KE
observables at the level of the effective capped law, at the literature
cooling time; it removes the cap's arbitrary kink but does not improve
the size distribution, whose remaining misfit lies in the mass model."

*Open for the discussion.* (a) lq as co-equal physical alternative next
to h405, or the recommended law (→ G5 attribution on it)? (b) Report c
as the midHot window 1.8–2.5 with c = 2.0 representative? (c) Optional
pooled run at the KE-best c 1.8 cell (τ ≈ 5.9, fitted KE₁ 0.632)?
(d) Histogram (from the shape discussion): question n* = 21 and the top
ladder rungs first? Re-include the bare fraction n = 0 (experiment 0.435
vs model ≈ 0.19; the solvated read excludes it by design, RQ8)?
**Bare part settled 2026-10-06 (user): n = 0 stays excluded, now on
physical grounds** (kinematic ceiling, opposite power response,
gas-phase spectral identity — `RESEARCH_QUESTIONS.md` RQ8 NB
2026-10-06); the model's bare class need not be the largest bin.
Session-local analyses (shape scan of all 200 runs, W₁ band split) were
not committed as scripts; their numbers are in D0 §1.G4.14.

**Discussion outcome (2026-10-06).** Numbers, anchors and provenance:
D0 §1.G4.14; log 2026-10-06; RQ8 NB 2026-10-06.

*1. lq vs capped cubic.*
- Same KE₁ lever in both laws: lowering E₀ (raising τ) moves fast
  suppressed ions (intact n = 21) into n = 1–3. The capped law gains
  more KE₁ that way (0.76 at τ 6.4) but mid-n KE falls out of band
  (midHot 0.81); lq keeps midHot ≈ 1. Both lose the n ≥ 10 tail.
- lq advantage: no cap, literature τ, consistent KE curve over n.
  Not a higher KE₁. KE₂/KE₁ ≈ 0.85 for every law vs 0.54 measured:
  the n = 1 deficit is not a drag question.
- n = 0 stays excluded (user decision): bare KE 3.7 eV exceeds the
  model ceiling 2.7 eV, shifts with laser power unlike n ≥ 1, and
  matches the gas-phase I⁺ spectrum. Part of bare is the no-He-contact
  twin of a solvated fragment (covariance), a channel not in the model.

*2. Why no n ≥ 16.*
- Pickup keeps n = 21 inside; E_int at exit = E₀·exp(−X/τ); after exit
  one He is shed per ≈ 9 meV. High n needs a cold exit, i.e. long
  exposure, i.e. a slow ion. High-n ions are marginal escapers (exit
  KE ≈ 0.05 eV); beyond X/τ ≈ 2 ions are retained. This caps n ≈ 15
  and makes the deep bins slow. Same ceiling for h405 and lq.
- Not levers: drag law, birth geometry / chord, `post_pickup` sampler
  (map reweighting only); twins already independent.
- Escape well: lowering E_bind alone (drag fixed) brings deep-bin KE to
  experiment (0.98 at 0.048 eV) but n ≥ 15 only to 0.007 (exp 0.057).
  Tier-0-consistent (18 Å refit) a shallower well forces stronger drag
  (0.080 eV: a +36 %; 0.048 eV fails the band) and the gain vanishes.
  TDDFT pins the escape cost.
- Per-ion E₀ spread would decouple cooling from slowing (projection);
  not adopted — single E₀ kept.
- Measured gap: KE width at fixed n ≈ 8× narrower than experiment.

*Thesis statement.* The model reproduces the solvated histogram to
n ≈ 15 and the n = 2–8 KE curve. The tail ends because cooling and
slowing share one channel (residence in He) and TDDFT fixes the escape
cost. Outlook: per-ion energy scatter or residence-independent cooling.

## G5. `post_pickup` for the thesis production law; h405 as the identifiability contrast (REPLACED 2026-10-08, user; EXECUTED 2026-10-09 — production law B2, effective variant B3)

**Decision (user, 2026-10-08): the linear + quadratic law is production
for the thesis** — the A5 cell: lq ridge c 2.0 amu/Å, a 16.717 amu/ps,
E_bind 0.1175 eV (18 Å fit, 3.5 Å density gate, uncapped), τ 6.4 ps,
E₀ 0.3155 eV; pooled reference `g414c200t64pooled` (N = 6000).
`post_pickup` is the physically motivated sampler (G1). h405
(`capped_cubic` v_c 5.5) is kept as the contrast law, not re-arbitrated.
The previous G5 (attribution A/B on the G4 bundle, then a ≥ 1-day
re-arbitration of h405) is superseded.

**Plan** (production pipeline: Stage I to 500 ps, E2 off, partner-aware
closure; seeds s1–s5 — r6 excluded as in the existing h405 `post_pickup`
pool, where it tripped the D4 marginal-partner safeguard; N = 1000 per
member, ≤ 2 concurrent, ion.npz dropped after detection verifies):

| step | cell | sampler | runs | purpose |
|---|---|---|---|---|
| A | lq A5, unchanged | `post_pickup` | 5 (s1–s5) | sampler effect on lq, directly comparable to h405's (D0 §15.8) |
| B1 | lq A5, E₀ column | `post_pickup` | ≈ 3 (s1) | place E₀ to restore n̄ ≈ 3.89 (projection from the h405 A/B: ≈ +0.03 eV); MD scan as in G9, not the twin (validated on the raw ensemble only) |
| B2 | lq A5 at the chosen E₀ | `post_pickup` | 5 (s1–s5) | the thesis production result under the physical sampler; pooled N = 5000 |
| — | h405 | `raw` / `post_pickup` | 0 | already pooled: `h405pooled` (s1–s5), `h405pppooled` (s1–s5), `h405g35` (3.5 Å gate) |

- **Comparison set.** lq raw (A5 pool) vs lq `post_pickup` (A, B2) vs
  h405 raw / `post_pickup`: full observable vector, realized droplet
  sizes, D4 firings tracked per run (G1 known risk).
- **Gate.** Existing MD clauses (n̄ ∈ [3.83, 4.31], n₁) reported; KE₁
  reported, not fitted.
- **Not in scope.** τ or drag re-tuning under `post_pickup` (only E₀
  re-placed); the G10 effective cell under `post_pickup`. *(Both lifted
  by the user during execution — see the outcome.)*

### G5 Outcome (EXECUTED 2026-10-08/09; pooled N = 6000; D0 §15.9)

**Decisions (user, 2026-10-09).** Thesis production law = **B2**
(`post_pickup`, lq ridge c 1.6); effective high-tail variant = **B3**
(B2 + shallow well). Both replace A5 / G10 in those roles.

What was run (production pipeline, s1–s5 + r6, no D4 firing anywhere):

| step | cell | result |
|---|---|---|
| A | A5, `post_pickup`, unchanged | n̄ +0.70, n₁ −0.04; W₁, midHot, KE₁, S flat |
| B1 | A5 + E₀ 0.3355 / 0.3455 / 0.3555 (s1) | n̄ restored at E₀ ≈ 0.337 but midHot ≈ 0.87, KE₁ 0.55 |
| B1′ | A5 + τ 6.8 / 7.2 / 7.6 (s1) | n̄ restored at τ ≈ 6.96 but midHot ≈ 0.85 |
| B1c | ridge c 1.6 / 1.8 + E₀ re-placed (s1) | c 1.6: n̄ 3.94, midHot 0.994 — holds both |
| **B2** | c 1.6, E₀ 0.3297, N = 6000 | **production** (table below) |
| B1w | B2 + E_bind 0.0482, a 19.8 / 20.8 / 21.8 (s1) | a 19.8: midHot 1.004, n ≥ 15 ×7 |
| **B3** | B2 + E_bind 0.0482, a 19.8, E₀ 0.3389, N = 6000 | **effective high-tail variant** |

**Why the drag and not E₀ / τ.** The chord → (n, KE) map is identical for
both samplers; the larger droplets only move ions from short chords (the
n = 1 window, 33–37 Å) to long ones. E₀ and τ change the map by extra
evaporation, which drops slow, heavy clusters into the low-n bins (KE per
n falls). The drag changes the map in KE directly: lower c (the ridge
freedom TDDFT leaves open) restores midHot at re-placed E₀. **c 1.6 is a
committed TDDFT ridge point** (`g412_ridge/c1p600.json`: 18 Å RMSE
0.1224 Å/ps vs band 0.19, escape 1, 9 Å exit margin 0.866 Å/ps; same as
A5 at F(3 Å/ps), 7 % lower at 9.7 Å/ps); all six B2 cfgs carry its a, c,
E_bind and the 3.5 Å gate bit-for-bit.

**Comparison (pooled; one read path, `tier2_g5_winners_table.py`; every
row anchored on its committed pooled row, 4 dp).**

| quantity | A5 (raw, c 2.0) | h405 (raw) | h405 (`post_pickup`, N = 5000) | G10 eff. (raw, 0.0482 / a 19.5) | A5 `post_pickup` unchanged | **B2 production** | **B3 effective** | experiment |
|---|---|---|---|---|---|---|---|---|
| sampler | raw | raw | post_pickup | raw | post_pickup | post_pickup | post_pickup | — |
| drag | lq c 2.0, a 16.72 | capped cubic b 2.515, v_c 5.5 | capped cubic b 2.515, v_c 5.5 | lq c 2.0, a 19.5 | lq c 2.0, a 16.72 | lq c 1.6, a 18.01 | lq c 1.6, a 19.8 | — |
| E_bind [eV] | 0.1175 | 0.1168 | 0.1168 | 0.0482 | 0.1175 | 0.1173 | 0.0482 | — |
| E₀ [eV] / τ [ps] | 0.3155 / 6.4 | 0.405 / 4.4 | 0.405 / 4.4 | 0.3216 / 6.4 | 0.3155 / 6.4 | 0.3297 / 6.4 | 0.3389 / 6.4 | — |
| density gate [Å] | 3.5 | 14.2 | 14.2 | 3.5 | 3.5 | 3.5 | 3.5 | — |
| Tier-0 (18 Å) | ✓ | old gate | old gate | ✗ (G10.5) | ✓ | ✓ | ✗ (G10.5) | — |
| n̄ (MD gate 3.83–4.31) | 3.886 | 3.890 | 4.594 | 4.009 | 4.581 | 3.903 | 3.956 | — |
| n₁ (solvated) | 0.204 | 0.209 | 0.177 | 0.207 | 0.163 | **0.211** | **0.211** | 0.310 |
| W₁ (size histogram) | 1.097 | **0.755** | 1.005 | 0.905 | 1.090 | 0.967 | 0.795 | 0 |
| midHot (n 2–8 KE) | 0.995 | 0.948 | 0.968 | 0.997 | 1.003 | **1.002** | 1.016 | 1 |
| KE₁ mean [eV] | 0.617 | 0.640 | **0.647** | 0.597 | 0.618 | 0.609 | 0.590 | 1.302 |
| KE₂ mean [eV] | 0.526 | 0.550 | **0.555** | 0.511 | 0.528 | 0.521 | 0.509 | 0.706 |
| deepKE (n 10–17) | 0.492 | 0.571 | 0.535 | 0.577 | 0.535 | 0.477 | **0.601** | 1 |
| KE n 10–14 / exp | 0.598 | 0.586 | 0.611 | 0.733 | 0.624 | 0.642 | **0.772** | 1 |
| n ≥ 15 (solvated) | 0.0006 | 0.0086 | **0.0095** | 0.0053 | 0.0012 | 0.0008 | 0.0075 | 0.057 |
| max n | 16 | 17 | 18 | 18 | 16 | 17 | 18 | ≥ 20 |
| trapped | 0.035 | 0.080 | 0.135 | 0.023 | 0.067 | 0.049 | 0.028 | — |
| suppressed | 0.067 | 0.192 | 0.135 | 0.090 | 0.038 | 0.111 | 0.151 | — |
| S | 1.958 | 1.705 | 1.991 | 1.605 | 1.928 | 1.705 | **1.505** | — |
| χ²_med (KE curve) | 507 | **456** | 572 | 581 | 511 | 527 | 674 | — |

Bold = best model value in the row (closest to experiment; S and χ²
lowest); n̄, trapped and suppressed carry no target. Experiment row from the references
(`ihe_ked`, abundance; solvated renormalization).

**Reading.**
- **B2 vs A5:** realistic droplets at no cost — n̄, n₁, midHot held; W₁
  −0.13 and S −0.25 better; KE₁ −0.008 eV. Both TDDFT-tied.
- **B3 vs B2:** n ≥ 15 ×9, max n 18, deep bins faster (0.77 of
  experiment), trapping halved, W₁ −0.17; costs KE₁ −0.019 eV, χ²_med
  +28 %, suppressed +0.04. Effective model (same framing as G10.4).
- **vs h405:** h405 still has the best size histogram (W₁ 0.755) and KE
  curve (χ² 456), at midHot 0.948 and with twice the trapping; B3
  matches its tail (n ≥ 15 0.0075 vs 0.0086) at midHot ≈ 1.
- **Still open (outlook):** n₁ ≈ 0.21 vs 0.31 and KE₁ ≈ half the
  experiment for every law (G8 / status 9–11); n ≥ 15 ≈ 8× below
  experiment even for B3 (status 11, `THESIS_DISCUSSION.md` §7).

Records: `data/runs/h2b_forward_model/g5_postpickup_table.csv` (all steps,
per member), `g5_winners_table.csv` (this table, full vector),
`g5_postpickup_chord.csv`; generator `gen_tier2_detfix_postpickup.py
--cell lq` (`--b1`, `--b1tau`, `--b1c`, `--b1w`, `--b2`, `--b3`); pool
targets `g5_lqpp6000`, `g5_lqppb2`, `g5_lqppb3`.

**h405 for the thesis: identifiability argument.** TDDFT (18 Å) fixes
F(v ≈ 3 Å/ps) and, below ≈ 0.08 eV, the well / drag split (G10.5) — not
the high-speed shape: a whole ridge F = a·v + c·v² passes the trace
(G4.12). The experiment selects c on the ridge (midHot window
c ≈ 1.8–2.5), and a differently shaped law (capped cubic, v_c 5.5)
lands comparably (KE₁ 0.640 vs 0.617, midHot 0.948 vs 0.995, W₁ 0.755
vs 1.097). KE₁ sits at ≈ half the experimental value for every law — a
mass-model question, not drag. Conclusion: neither TDDFT nor the detected
I⁺(He)ₙ data identify a unique drag function in He nanodroplets.
G5 adds a second instance: the c the experiment selects on the ridge
depends on the droplet description (raw → c 2.0, `post_pickup` → c 1.6).
Caveat: h405 is not a controlled A/B partner — it also differs in density
gate (14.2 vs 3.5 Å), τ (4.4 vs 6.4 ps) and Tier-0 status (bundle fitted
at the old gate); `h405g35` (3.5 Å gate: W₁ 0.818, KE₁ 0.643) is the
closer like-for-like comparison.

**Code (requires `[PROCEED TO IMPLEMENTATION]`).** An lq A5 cell option in
`gen_tier2_detfix_postpickup.py` (+ E₀ column), pool targets and scorer
rows (full vector, realized sizes, D4 firings).

**Cost.** ≈ 13 runs, < 1 h MD.

## G6. Chord-map figure at the new standing point (after G5 B2)

- **Change.** Re-run the `emission_chords` read
  (`tier2_detfix_postpickup_table.py`) on the new N = 6000 battery and
  plot the chord → (n, KE) map, overlaid on the old-bundle map
  (`data/runs/h2b_forward_model/detfix_postpickup_chord.csv`). Panels:
  n vs chord, KE vs chord with the n = 1 window marked, chord
  distribution per sampler.
- **Rationale.** The figure carries the thesis argument: geometry
  reweights the chord distribution; only the map (drag bundle, budget,
  flight mass) moves KE₁. The old-vs-new overlay shows the G2/G4 bundle
  change acting on the map directly.
- **Evidence.** D0 §15.8 reads 5–6; §14.3 refinement note.
- **Cost.** Zero MD; scorer re-run + plot script, < 1 h.

## G7. KE₁ geometry-invariance caveat (after G6)

- **Change.** Thesis text only. Re-check on the G6 map that the n = 1
  window is still narrow and the map near-deterministic; then state:
  - "KE₁ is geometry-invariant" holds within the deterministic model;
  - the experimental KE width at fixed n shows the real map scatters
    (model KE₁ SD 0.038 eV at the old bundle), and with scatter a
    chord-window selection could bias KE₁ slightly;
  - map scatter (Tier-3 noise) is the outlook, not modelled here.
- **Rationale.** The near-determinism is partly a consequence of the
  missing fluctuation channel; "identical map in both arms" is expected
  by construction (same physics parameters). The claim needs its scope.
- **Evidence.** D0 §15.8 read 6 (caveats); CLAUDE.md Tier 3.
- **Cost.** None beyond G6.

## G8. Bare-start channel screen (agreed 2026-10-06; frozen before any run)

**Question.** If some fragments start with no He shell, which of them
leave the droplet, at which n, and with what KE? This is a screen for a
later P1/P2 mixture (bare class + intact class) that could supply the
missing fast low-n ions. It is a modelling assumption here; the physical
origin is decided later.

**Bare class (Case 1, user decision).** n₀ = 0 for every fragment, both
fragments of every pair (P1). E_int(0) = 0: the T6 `sigma_proportional`
law gives (Σ(0)/Σ(21)) = 0 (no shell, no reservoir). Pickup adds
f_ret·D₀ ≈ 1–2 meV per He, below D₀(1) = 20.3 meV, so final n ≈ number of
pickups (sourced λ₀ = 0.9/ps). Case 2 (bare ion keeps E₀) is not run.

**Runs (2, MD).** Production pipeline unchanged (Stage I 500 ps, E2 off,
partner-aware closure), N = 1000, seed s1 (the parent's neutral stage,
so births are identical to the parent run). Parents:
- lq: the A5 cell (c 2.0, τ 6.4, E₀ 0.3155, 3.5 Å gate), s1 member;
- h405: production detfix s1 (14.2 Å gate).
Each cfg diffs from its parent in exactly `initial_shell_model`
(`density_tied` → `bare`). The parents already run T5 `density_tied` +
T6 `sigma_proportional`, both inert for them (n₀ = 21 at production
depths); found at the dry-run oracle, before any MD.

**Code (`[PROCEED TO IMPLEMENTATION]` given 2026-10-06).**
`initial_shell_model = "bare"` (biphasic-only; requires
`sigma_proportional`, so the deposit cannot silently stay on a bare
ion); generator `scripts/gen_tier2_bare_screen.py`; scorer
`scripts/post_processing/tier2_bare_screen_table.py`. Before `ion.npz`
is dropped, a small per-ion summary is saved (n at handover, pickup
exposure ∫ρ̂ dt, exit time, max E_int).

**Readout (reported; no weight scan, no pass/fail rule).** Per law:
fates (scored / retained / suppressed); detected histogram (n = 0 share,
solvated n ≥ 1, share at n ≥ 4, bands); KE per n (mean, median, IQR,
count) beside the experiment (mean, median) and the same-law intact
reference (lq A5 pool; `pooled_new6000`); n vs exposure. The user judges
from the tables whether the mixture is worth a design (P1 fit, then P2).

**Anchors.** The scorer rescoring the two parent runs reproduces their
committed rows (h405 s1 = `s1_new`; lq A5 s1 member = its
`g414_table.csv` row, 4 dp).

- **Status — EXECUTED 2026-10-06 (measured; anchors pass).** Bare ions
  never evaporate (E_int ≤ 7 meV); n = pickup count, peaking at n = 2–3
  (n₁ share 0.15 lq / 0.20 h405; 35–45 % at n ≥ 4; max n 11–12).
  KE₁: lq 0.491 (below intact 0.621), h405 0.865 (above intact 0.644,
  below the experimental median 1.128); h405 bare overshoots n = 2–4;
  KE₂/KE₁ 0.87–0.90 (exp 0.54). Only the KE width at fixed n moves
  toward experiment. The KE gain is a property of the cap (lq's ram drag
  takes it back). A bare/intact mixture cannot supply the fast n = 1
  ions; mixture design: user decision. Record: D0 §1.G8.
- **Decision (user, 2026-10-06): not pursued; the n = 1 KE deficit goes
  to the outlook.** Tested inside the model and ruled out as the n = 1
  lever: drag-law shape (capped, lq; G4.12–G4.14), the in-surface τ / E₀
  freedom (§3.5g, τ lines), a uniform flight-mass change (midHot,
  `TIER2_MASS_SCENARIOS` §5), exit stripping and the (C) channel mixture,
  velocity-dependent sticking (§13.4), birth geometry (reweights the map
  only), the escape well (Tier-0-pinned, G4.15), and a bare start (G8).
  The graded partial strip (Route A, §12) rests on the cap-only §4
  identity and has no early-shedding mechanism. Remaining candidates are
  outside this model: the source channel (the experimental n = 1 KE is a
  two-population mean, a fast tail plus a core; KE₂/KE₁ 0.54 vs ≈ 0.85
  for every model variant; RQ8 / RQ1), map scatter (KE width ≈ 8× too
  narrow, Tier-3 noise), and a sourced I⁺–He binding (M1).

## G9. Cold-onset minority screen at the lq A5 cell (agreed 2026-10-08; frozen before any run)

**Question.** Does a small onset-energy spread below E₀ (most ions at
the main E₀, a fraction w colder by ≤ 0.05 eV) extend the high-n tail
without breaking n̄ / KE₁ / midHot? The earlier two-E₀ projection (D0
§1.G4.14, "Onset-energy spread": 10–20 % at ½ E₀ ≈ −0.16 eV) reached
n ≥ 15 but broke n̄ (5.3–6.4) and overshot deep-bin KE, because each
cold ion jumped ≈ 6 rungs and the main E₀ was not re-placed. Modelling
assumption; physical origin decided later.

**Projection (frozen, rule of thumb).** Δn ≈ ΔE₀·e^(−X/τ)/9.2 meV:
−0.05 eV → +3 rungs on short chords, +2 typical, +0.7 deep. n̄ shift
≈ +1.25·w (point mass at −0.05) / ≈ +0.6·w (uniform), so the main E₀
moves up by ≈ 0.05·w to hold n̄. Expected: n ≥ 15 rises, deep-bin KE
rises, KE width at fixed n widens, KE₁ ≈ unchanged; no mass above
n ≈ 17–18 (weak top rungs, n* = 21).

**Columns (MD; seed s1 = the A5 s1 births, CRN).** The A5 law (c 2.0,
a / E_bind from the committed ridge JSON, 3.5 Å gate, τ 6.4, uncapped),
production pipeline, N = 1000, E₀ grid step 0.0125 eV:
0.2655, 0.2780, 0.2905, 0.3030, **0.3155 (= committed A5 s1)**, 0.3280,
0.3405 → **6 new runs**. Each cfg diffs from the A5 s1 cfg in exactly
`internal_energy_partition_fraction` (= E₀ / 2.70). The committed G4.13
s1 cells at E₀ 0.300 / 0.310 / 0.320 / 0.330 are a cross-check only
(n̄ monotone across the merged grid within seed noise), never mixture
columns.

**Mixture (zero MD, exact pooling).** Ions are independent apart from
their Coulomb partner, and the columns share births. A molecule's two
fragments are assigned to one column (P1). The assignment uses a
dedicated `np.random.default_rng(20261008)` permutation, not the
simulation RNG. The mixed ensemble is the masked reads
(`include_mask`) of the columns pooled.
- **Main peak:** a blend of two adjacent columns (width ≤ 0.0125 eV),
  with the blend fraction set so that the expected n̄ (weighted column
  histograms) = **3.89** (user: re-place the main E₀).
- **Shape P (point):** fraction w at main − 0.05 eV (the same blend
  shifted down 4 grid steps).
- **Shape U (uniform):** fraction w uniform over [main − 0.05, main]:
  trapezoid over the 5 shifted nodes (1/8, 1/4, 1/4, 1/4, 1/8).
- w ∈ {0, 0.05, …, 0.50}; w = 0 is the single-E₀ reference at n̄ 3.89.

**Readout per (shape, w)** (reported; best rows marked, nothing
adopted): realised n̄, n₁, n ≥ 15 share, max n, histogram bands, W₁,
deepKE, KE n 10–14 / experiment, KE₁, KE₂/KE₁, midHot, KE IQR/median
at n = 1–4, trap. Marked: min W₁ and max n ≥ 15 among rows with midHot
∈ [0.85, 1.15] and KE₁ ≥ (w = 0 row) − 0.01. Against the A5 s1 member,
the A5 pool and experiment.

**Anchors.** The A5 s1 column rescored = committed `c200t64confs1`
(4 dp). The w = 0 mixture at main = 0.3155 alone reproduces that row
exactly.

**Code (requires `[PROCEED TO IMPLEMENTATION]`).** Generator
`scripts/gen_tier2_g9_cold_onset.py` (column cfgs from the A5 s1 build,
diff oracle, `ion.npz` dropped after detection verifies); scorer
`scripts/post_processing/tier2_g9_cold_onset_table.py` (the mixture
builder + readout; CSVs force-added). Tests: assignment partitions the
molecules with both fragments in one column; trapezoid / blend weights
sum to 1; a stub mixture's expected n̄ = the weighted column mean.

**Cost.** 6 runs, ≤ 2 concurrent; scorer zero MD. Then a user
discussion on how to continue (pooled seeds, shape, physical origin).

- **Status — EXECUTED 2026-10-08 (measured; anchors pass).** No tail:
  n ≥ 15 stays 0.0011–0.0018 for every shape and w (exp 0.057); W₁ and
  deepKE move within seed noise; KE₁ falls (P 0.623 → 0.585, U → 0.605)
  because the main E₀ must rise. Even the whole ensemble at −0.05 eV
  (E₀ 0.2655) gives n ≥ 15 only 0.003 (n̄ 5.86). Reading: an escaper's
  exit E_int ≥ 0.17·E₀ (lq exposure ceiling), so n ≥ 13 needs
  E₀ ≲ 0.18 eV — a distinct cold population, not a small spread.
  Record: D0 §1.G9.

## G10. Shallow escape well + drag retune: higher fragments at matched midHot (MEASURED 2026-10-08; pooled N = 6000)

Sources: D0 §1.G4.14 (well-only bullet); `data/runs/h2b_forward_model/
g415_table.csv` (pooled rows), `g414_table.csv` (A5 / h405 pools),
`atlas_g415_summary.csv` (twin); log entries 2026-10-08.

### G10.1 What was run

- **A5 pool (TDDFT-tied reference).** lq ridge law c 2.0 amu/Å,
  a 16.717 amu/ps, E_bind 0.1175 eV — all three from the 18 Å fit at the
  3.5 Å density gate; τ 6.4 ps, E₀ 0.3155 eV.
- **Shallow well, A5 law (0.0482, c 2.0).** Only E_bind lowered to
  0.0482 eV: fewer trapped ions, faster deep bins, but midHot 1.13.
- **Shallow well + drag retune.** The twin placed midHot ≈ 1 at c 2.3
  (A5 a) or a 19.5 (c 2.0), with E₀ re-placed to hold n̄
  (0.3235 / 0.3216 eV). MD hit the twin's predictions (midHot 0.988 /
  0.996 vs 0.988 / 0.998 on seed s1).
- All arms: production pipeline (Stage I 500 ps, E2 skipped, partner-aware
  closure), pooled N = 6000 (s1–s5 + r6, production seeds). Figure
  containers `…_N6000_detfix_conf270_g414c200t64pooled_th500` and
  `…_g415wo048{,c230,a195}pooled_th500`.

### G10.2 Pooled result (MEASURED, N = 6000 each)

| | A5 pool (TDDFT-tied) | 0.0482, c 2.0 | 0.0482, c 2.3 | **0.0482, a 19.5** | h405 pool | experiment |
|---|---|---|---|---|---|---|
| midHot (n 2–8 KE) | 0.995 | 1.130 | 0.990 | **0.997** | 0.948 | 1 |
| χ²_med (KE curve) | 507 | 614 | 557 | 581 | 456 | — |
| n ≥ 15 (solvated) | 0.0006 | 0.0062 | 0.0059 | **0.0053** | 0.0086 | 0.057 |
| max n | 16–17 | 17 | 18 | **18** | — | ≥ 20 (window) |
| KE n 10–14 / exp | ≈ 0.6 (s1) | 0.93 | 0.77 | **0.73** | — | 1 |
| deepKE (n 10–17) | 0.49 | 0.74 | 0.61 | **0.58** | 0.57 | 1 |
| W₁ (size histogram) | 1.097 | 0.950 | 0.879 | **0.905** | 0.755 | 0 |
| trapped | 0.035 | 0.013 | 0.019 | **0.023** | 0.080 | — |
| n₁ (solvated) | 0.204 | 0.201 | 0.209 | **0.207** | 0.209 | 0.310 |
| n̄ | 3.89 | 4.05 | 4.01 | **4.01** | 3.89 | — |
| KE₁ [eV] | 0.617 | 0.653 | 0.584 | **0.597** | 0.640 | 1.302 |
| suppressed | 0.067 | 0.069 | 0.095 | **0.090** | — | — |

Member spread of n ≥ 15 for the a 19.5 cell: 0.0045–0.0079 (N = 1000
each). "—" = not read for that pool.

### G10.3 Reading

- **At the same midHot as A5, the shallow well + retuned drag gives:**
  ≈ 9× more n ≥ 15, max n 18 instead of 16–17, faster deep bins
  (KE n 10–14 0.73 vs ≈ 0.6 of experiment), fewer trapped ions, and a
  better size histogram (W₁ 0.905 vs 1.097).
- **Cost:** KE₁ −0.02 eV (0.597 vs 0.617), χ²_med 581 vs 507, suppressed
  share +0.02.
- **c vs a:** tie on the tail; c 2.3 is slightly ahead on W₁ / deep KE /
  trapping, a 19.5 keeps 0.013 eV more KE₁. The thesis uses a 19.5.
- **Mechanism (`THESIS_DISCUSSION.md` §7):** the lower escape cost lets more exposed, colder
  ions out, which raises the escaper exposure ceiling (`THESIS_DISCUSSION.md` §7.2). The drag
  retune then restores the mid-n KE without giving that back.
- **What it does not do.** The n ≥ 15 abundance is still ≈ 10× below
  experiment, and the deficit is spread over n = 15–18, not only n = 19–20:
  the model's n 15–18 is ≈ 0.005, experiment 0.046. Whether the top bins
  (n = 19–20, 0.011 of solvated) carry experimental systematics
  (background, mass-window edge) is a question for the reference data;
  the smooth-tail argument of `THESIS_DISCUSSION.md` §7.3 (missing threshold fluctuations)
  remains the model-side explanation.

### G10.4 Status for the thesis (framing)

- **A5 pool = TDDFT-based result.** Every drag and well parameter comes
  from the 18 Å fit; τ ≈ the sourced 6.55 ps.
- **0.0482 / a 19.5 = improvement for higher fragments, an effective
  model.** It is **not Tier-0-consistent**: a is raised and the well is
  lowered independently of the 18 Å trace, and the Tier-0 refit at
  E_bind 0.0482 eV fails the 18 Å band (a 28.0, RMSE 0.226 > 0.19;
  D0 §1.G4.14). **Thesis statement (user, 2026-10-08):** the 0.0482 /
  a 19.5 cell is an effective model that improves the higher fragments,
  at the price of departing from the 18 Å TDDFT trace (G10.5). Do not
  claim that the well / drag split is free at this depth: G10.5 shows the
  trace excludes it.

### G10.5 Tier-0 check: the high-tail fix cannot match the 18 Å trace (MEASURED 2026-10-08, zero MD)

Each law scored on the 18 Å Tier-0 trace with the Tier-0 engine
(`evaluate_form_objective`; N = 50, seed 20260604, 20 ps at 0.01 ps,
3.5 Å gate; RMSE of |v₂| over the 4.5–14.8 ps window; band RMSE ≤
0.19 Å/ps and escape = 1). Anchors: the committed A5 point (`c2p000.json`)
and the committed 0.0482 refit (`c2p000_e0p0482.json`) reproduce their
stored RMSE exactly. Session-local evaluator (scratchpad), no files
written.

| law | a [amu/ps] | c [amu/Å] | E_bind [eV] | RMSE [Å/ps] | band | bias MD − TDDFT [Å/ps], 3–6 / 6–10 / 10–14 / 14–14.8 ps |
|---|---|---|---|---|---|---|
| A5 (TDDFT-tied) | 16.72 | 2.0 | 0.1175 | **0.122** | pass | +0.21 / −0.09 / 0.00 / +0.16 |
| well-only, A5 law | 16.72 | 2.0 | 0.0482 | 0.719 | fail | +0.42 / +0.45 / +0.89 / +1.13 |
| G10 c 2.3 | 16.72 | 2.3 | 0.0482 | **0.661** | fail | +0.38 / +0.38 / +0.82 / +1.06 |
| **G10 a 19.5** | 19.5 | 2.0 | 0.0482 | **0.557** | fail | +0.30 / +0.27 / +0.70 / +0.95 |
| best a at c 2.0 (committed refit) | 28.02 | 2.0 | 0.0482 | 0.226 | fail | −0.04 / −0.24 / +0.14 / +0.39 |
| best a at c 2.3 (1-D refit) | 27.14 | 2.3 | 0.0482 | 0.227 | fail | −0.04 / −0.24 / +0.14 / +0.39 |
| refit at 0.080 (c 2.0, for reference) | 22.76 | 2.0 | 0.080 | 0.159 | pass | — |

- **The G10 cells fail the band by ≈ 3×** (RMSE 0.56–0.66 vs 0.19; A5
  0.12). The ion is too fast throughout the window, and the excess grows
  to ≈ +1 Å/ps by 14 ps: without the deep well the TDDFT ion is not slowed
  enough near the surface.
- **No a rescues E_bind 0.0482.** The best a (≈ 27–28, c 2.0 or 2.3)
  still fails (0.226). Its residual has a fixed shape: too slow at
  6–10 ps (stronger drag), too fast at 10–15 ps (missing well toll). The
  18 Å trace therefore sees the well itself in its late window, not only
  the total escape cost.
- **Scope of the well / drag degeneracy.** Stronger drag compensates a
  shallower well within the band down to ≈ 0.080 eV (refit RMSE 0.159);
  at 0.0482 eV the split is excluded by the trace. The band edge lies
  between 0.048 and 0.080 eV (not searched).
- **Consequence.** The higher-fragment improvement of G10.2 cannot be
  obtained from a TDDFT-consistent escape well at this depth. It is
  reported as an effective model with a measured Tier-0 departure
  (RMSE 0.557 for a 19.5 vs 0.122 for A5).

## G11. Thesis reframing, Tier-0 mass m_eff = 21 He, and the remaining plan (agreed 2026-10-09; nothing run)

### G11.1 Reframed storyline — "TDDFT constrains, the experiment selects"

| tier | old framing | new framing |
|---|---|---|
| Tier 0 | fit a drag law to the 9 and 18 Å TDDFT trajectories | the 18 Å trace pins F near v ≈ 3 Å/ps and the well / drag split at low E (G10.5), not the high-speed shape: the output is a **one-parameter ridge** F = a·v + c·v² (G4.12), not a point. 9 Å = validity check only (exit bound; non-radial, transverse motion). The cross-case cubic preference (G4.11) is contested by the radial 9 Å runs → caveat, not a foundation |
| Tier 1 | does variable mass matter? (it did not) | consistency check of the Tier-0 fixed-mass extraction against production's variable mass (G11.3) |
| Tier 2 | sweep v_c, E₀, τ until the experiment matches | two knobs only: **c on the ridge** and **E₀**; τ 6.4 ≈ the sourced 6.55 ps, droplet ensemble physical (`post_pickup`, G1/G5). The experiment selects c, and the selected c depends on the droplet description (raw → 2.0, `post_pickup` → 1.6). With h405 (differently shaped law, comparable landing) this is the identifiability result: neither TDDFT nor the detected I⁺(He)ₙ data pin a unique drag function |
| effective | — | B3: a shallow well buys the high tail at the price of departing from the TDDFT trace (G10.5); kept as an explicit departure |
| outlook | — | n₁ share and KE₁, the n ≥ 15 deficit, the KE width at fixed n — all outside the deterministic model (source channel, Tier-3 collision noise, sourced I⁺–He D₀ / M1) |

Fewer tuned knobs than before; the non-uniqueness is a result, not a
weakness.

### G11.2 Tier-0 extraction mass → m_eff = 21 He (user decision 2026-10-09; supersedes G3)

- **Decision.** Tier-0 extraction at m_eff = 21 He (I + 21 He ≈ 210.96
  amu) instead of 19 He (202.954 amu).
- **Rationale.** G3's argument was a 9 Å argument (window-representative
  mass over the 9 Å shell decline 21 → 19 → 14; "no 18 Å loss curve").
  Tier 0 is now 18 Å-only (G4.9), and in the 18 Å TDDFT run the shell
  stays at 21 throughout — 21 is the actual mass there; 19 was an
  effective choice. Expected effect small but not zero.
- **Consequences.**
  1. **Ridge refit.** In F_drag = F_C − m·dv/dt a heavier ion with the
     same observed deceleration needs more drag: a and E_bind per c rise
     by a few %. Refit the full ridge grid (c 0–7.2 as G4.12) at 18 Å,
     3.5 Å gate, plus the 9 Å exit bound. Tier-0 engine only, no battery
     MD. The m = 19 fits stay committed as the record.
  2. **B2 / B3 re-landing.** Both were built on the m = 19 ridge; ≈ 4 %
     more drag lowers midHot slightly. One s1 run at the refit c 1.6
     point (E₀ unchanged); a short c re-scan only if midHot moves (the
     landing may slide toward c ≈ 1.5); re-pool N = 6000. B3 likewise (a
     re-scanned at the 0.0482 well).
  3. **Mass-pairing defence simplifies** (CLAUDE.md §6.5–6.7: constant-mass
     coefficients under variable-mass production). The current defence
     is "m ≈ 19 He = m_eff mid-window"; at m_eff = 21 He the fit mass
     equals production's n₀ = 21 start exactly at birth.
  4. **9 Å exit bound at m = 21 is slightly less strict** (a heavier ion
     decelerates less under the same drag and escapes more easily). Run
     at 21 for consistency; state the direction in the thesis.
  5. **Code.** `M_EFF_AMU` (202.954) is stamped into every bundle and
     checked exactly by the loader (`tier0_common.py`,
     `presets._DRAG_M_EFF_AMU`). Make the extraction mass a parameter of
     the refit and write the m = 21 bundles to a new directory — nothing
     committed is overwritten. Requires `[PROCEED TO IMPLEMENTATION]`.

### G11.3 Tier 1 is not void (relation to the bare-start result G8)

- **Tier-1a** asks whether the TDDFT-window kinematics (first ~15–30 ps)
  change when the ion sheds 21 → 19 → 14 He instead of carrying a fixed
  mass: a ~25 % change *inside* the range TDDFT shows. With
  continuous-velocity shedding the answer was no
  (`drag_migration_log_tier1a.md`, Slice C).
- **G8 (bare start)** changes the initial condition far *outside* that
  range — n₀ = 0 (127 vs ~211 amu), no E_int reservoir, no evaporation —
  and acts at the detector (~500 ps) through pickup statistics and a
  lighter ion.
- **No contradiction; both belong in the thesis.** Within the
  TDDFT-supported shedding the fixed-mass extraction is adequate
  (validates Tier 0 → Tier 2); the initial shell state, which TDDFT does
  not constrain, is a strong lever on the observables (model scope; feeds
  the n = 1 outlook).
- **Redo required.** The old Tier-1a used the old law and the 14.2 Å
  gate. At m_eff = 21 the question sharpens: at 18 Å the shell is
  constant (nothing to test); at 9 Å it is exactly "does ignoring the
  21 → 14 shedding change the trajectory" — fixed 21 vs anchored
  continuous-velocity shedding, at the B2 law (open: also A5 as a ridge
  contrast). Deterministic, cheap.
- **Superseded (user, 2026-10-13): Tier 1 retired.** Tier-1a existed to
  show that the assumed mass does not matter for the 9 Å extraction
  (shedding 21 → 14 at known times). Since G4.9 the extraction is 18 Å
  only, where the shell stays at 21 He, so m = 21 is the actual mass, not
  an assumption; 9 Å is a one-sided exit bound (margin +0.82 Å/ps), not a
  fit target. Nothing is left for Tier 1 to validate. The old Tier-1a
  results (old law, 14.2 Å gate) are dropped from the thesis (the later
  ridge / gate changes leave them without meaning). The mass question
  moves to the discussion: G11.5.

### G11.4 Plan (order)

| # | step | MD | depends on |
|---|---|---|---|
| 1 | G3 rewritten to m_eff = 21 He (done here: G11.2) | — | — |
| 2 | Tier-0 ridge refit at m = 21 + 9 Å exit bound — **done 2026-10-09** | none (Tier-0 engine) | 1 |
| 3 | Tier-0 figure: 18 Å trace with the m = 21 ridge points c 1.6 / 2.0 + h405; 9 Å check — **done 2026-10-09/10** (B2 / B3 / h405) | none | 2 |
| 4 | ~~B2 re-landing at the m = 21 c 1.6 point (s1; c re-scan only if needed), re-pool N = 6000; B3 likewise; update the G5 comparison table~~ — **dropped** (B2 / A5 pass at m = 21 unchanged) | ~1–2 h | 2 |
| 5 | ~~Tier-1a at 9 Å: fixed 21 vs anchored continuous shedding, at the B2 law~~ — **dropped** (Tier 1 retired, G11.3; mass question → discussion, G11.5) | — | 2 |
| 6 | G6 chord map, G7 caveat, one-page storyline section (`THESIS_DISCUSSION.md`) per G11.1; mass-influence discussion section per G11.5 | none | 4 |

Steps 3 and 5 can run in parallel with 4. Code for 2, 4, 5 requires
`[PROCEED TO IMPLEMENTATION]`.

- **Status — step 2 EXECUTED 2026-10-09** (`method_b_g412_ridge.py`
  modes `*_m21`; m = I⁺ + 21 He = 210.9546 amu, production m(n);
  `g11_ridge_m21/ridge_m21.csv`, `9A/g11_exit_bound_m21/`; m = 19
  artifacts untouched; oracle: mass path at m = 19 bit-exact). All 12
  ridge points pass 18 Å (RMSE 0.113–0.119); 9 Å cut unchanged (c ≤ 5
  pass, 7.224 fails, −0.130). a per c **falls** 1.4–2.3 % (c 1.6
  17.757, c 2.0 16.480), E_bind +0.0004 eV — G11.2 consequence 1 and 4
  had the sign wrong: the heavier ion also gains less Coulomb speed, so
  less drag fits, and the 9 Å margin shrinks 0.00–0.04 Å/ps.
  **B2 / A5 re-scored unchanged at m = 21 pass** (RMSE 0.1182 / 0.1178,
  escape 1, 9 Å margin +0.823 / +0.739): both stay Tier-0-consistent at
  m = 21 → **step 4 (re-landing / re-pool) dropped** (user request: keep
  the pools); B3 unaffected (already not Tier-0-consistent). D0 §1.G4.12.
- **Status — step 3 EXECUTED 2026-10-09** (`scripts/post_processing/
  tier0_thesis_figure.py` → `data/figures/thesis_tier0/tier0_B2_A5_h405
  .{png,pdf}` + series / scores CSV; figure not committed). B2 / A5
  (production coefficients, m = 21, 3.5 Å gate) and h405 (capped cubic,
  bundle m = 19, 14.2 Å gate); anchors bit-exact (B2 / A5 check JSONs,
  h405 bundle per_case RMSE). 18 Å RMSE 0.118 / 0.118 / 0.134 (all
  pass). **9 Å exit margin +0.82 / +0.74 / −1.21 Å/ps: h405 fails** — it
  over-dissipates and tracks the transverse atom 2 (same reading as the
  cubic in G4.11); the lq laws under-dissipate vs atom 1 (one-sided bound).
  **Revised (user, 2026-10-09):** A5 replaced by B3 in the figure
  (`tier0_B2_B3_h405.*`; B3 at m = 21: 18 Å RMSE 0.580 (fails, G10.5),
  9 Å margin +1.23); raw TDDFT |v1|, |v2| behind panel (a). Thesis Tier-0
  trace + R(t) figure stays `tier0_drag_comparison.py`, now set to B2 at
  m = 21 (`gen_tier0_runs.py` RUN_TAG `B2m21`; gate RMSE 0.1182, = check).
- **Where the results live (committed):** ridge + B2 / A5 checks
  `data/reference/drag/18A/trajectory_matching_srho3p5A/g11_ridge_m21/`
  (`ridge_m21.csv` = summary table; `c*.json`, `check_B2/A5.json`); 9 Å
  bound `data/reference/drag/9A/g11_exit_bound_m21/`; m = 19 record
  unchanged in `g412_ridge/` + `9A/g412_exit_bound/`. Figures
  (`data/figures/thesis_tier0/`) are not committed — regenerate with
  `tier0_thesis_figure.py` / `gen_tier0_runs.py` + `tier0_drag_comparison.py`.

**Open question (user).** Storyline section: drafted here or written
directly into the thesis? (Tier-1a scope: moot, step 5 dropped.)

### G11.5 Mass influence — discussion section, not validation (agreed 2026-10-13; zero MD)

Placed in the thesis discussion (after the Tier-2 results), not in the
validation-strategy chapter: with an 18 Å-only extraction at the true
shell mass there is nothing mass-related left to validate; what remains
is model scope. Compiled from existing results only (option (a)).

| level | question | result (existing) | source |
|---|---|---|---|
| extraction | sensitivity of the fitted law to the fit mass | 19 → 21 He: a −1.4 to −2.3 % per c, E_bind +0.0004 eV, 9 Å margin −0.00 to −0.04 Å/ps; B2 / A5 pass unchanged | G11.4 step 2, D0 §1.G4.12 |
| transfer | fixed-mass law under production's variable mass | exact at birth (n₀ = 21 = fit mass); afterwards rests on the mass-agnostic force law (γ never takes m) — **stated assumption**: the drag force does not depend on the shell size n, untestable with the TDDFT references (no small-n trace) | G11.2 (3), CLAUDE.md friction convention |
| initial shell (outlook) | ions born bare (n₀ = 0) | see below | G8, D0 §1.G8 |

**Bare start — outlook, defined.** The drag law is extracted for one
complex only, I⁺ + 21 He (the 18 Å shell). G8 applies that law unchanged
to ions born bare (n₀ = 0, E_int(0) = 0). Measured (N = 1000, s1, lq A5
cell = the ridge point next to B2; h405 as contrast):
- bare ions do not evaporate; n = pickup count (peak n = 2–3);
- KE₁: lq 0.491 eV, below the intact 0.621 and far below the
  experimental median 1.128; the lighter ion's larger Coulomb speed is
  taken back by the ridge law's v² (ram) term;
- h405 (capped high-speed drag) keeps the gain (KE₁ 0.865) — so the
  bare-start KE₁ is set by the high-speed drag, but even the capped law
  stays below the experiment.

Reading for the thesis: a bare start, run with a drag law fitted for the
21-He complex, does not produce the fast n = 1 ions. This is
**consistent with** the 21-He law over-dissipating on an ion it was not
extracted for (a bare I⁺ has a different solvation structure / bubble),
**not shown**: G8 cannot separate the drag from pickup and the lighter
flight mass, and no TDDFT reference for a bare / small-n ion exists.
Outlook: a drag law for small n (TDDFT for I⁺Heₙ, n ≪ 21) before a
bare / partly-solvated source channel can be tested; ties to the n = 1
outlook (G8 decision, RQ1 / RQ8, M1).

The old Tier-1a results (`drag_migration_log_tier1a.md`) are not used.

## Thesis framing implied

- Thesis table: bundle effect (A − reference), sampler effect (B − A),
  re-arbitration (C − B).
- The sampler is a stated modelling choice with its measured
  sensitivity (D0 §15.8), not a validated property.
- m_eff = 21 He (the 18 Å shell mass) is the Tier-0 extraction mass
  (G11.2; supersedes the G3 window-representative 19 He).
