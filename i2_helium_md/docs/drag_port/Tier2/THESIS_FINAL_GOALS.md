# THESIS_FINAL_GOALS — agreed changes and tests before thesis completion

Settled in the final discussion round (2026-10-02). This document holds
only agreed goals. Results go to the migration log and D0 as usual.

**Status (2026-10-05):** G4 **executed** (G4.8). The adoption decision
is superseded by **G4.9** (18 Å-only refit; 9 Å held out as
non-radial), **executed 2026-10-05: winner `cubic` (`srho3p5A18A_cubic`,
b 2.2286, E_bind 0.1214 eV)**; adoption for G5 is gated by **G4.10**
(force-balance cross-check), **executed 2026-10-05: 18 Å b_A = 1.618
outside the band [1.803, 2.654] → DECISION NEEDED** (G4.10 status).
G1, G2 (production switch), G5–G7 not started.

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

- **Status — EXECUTED 2026-10-05 (measured; DECISION NEEDED).** 18 Å
  cubic b_A = 1.618 → **disagrees** with the band; 9 Å 2.013 and joint
  1.991 lie inside. In the 18 Å interior the drag falls with depth at
  constant v (F/v³ 2.16 → 1.01; the deepest bin matches 2.23). Full record:
  D0 §1.G4.10, log 2026-10-05 G4.10 entry, artifacts
  `data/reference/drag/force_balance_srho3p5A/`.

  > **DECISION NEEDED.** (a) Adopt `srho3p5A18A_cubic` for G5 and record
  > the depth-dependent near-surface drag; or (b) open the near-surface
  > well (14.2 Å) vs density (3.5 Å) width inconsistency before G5.

## G5. Attribution batteries, then re-arbitration under `post_pickup`

Run in this order on the production pipeline (Stage I to 500 ps, E2
off, partner-aware closure; `exclude_all_coupled`):

| step | bundle | sampler | point | purpose | cost |
|---|---|---|---|---|---|
| A | new (G4) | `raw` | h405 | bundle effect alone (incl. on KE₁), vs the N = 6000 reference | ≈ 15 min MD (s1–s5) |
| B | new (G4) | `post_pickup` | h405 | sampler effect under the new bundle, vs D0 §15.8 | ≈ 15 min MD |
| C | new (G4) | `post_pickup` | re-arbitrated | new standing point | ≥ 1 day |

Step C:

1. **Twin ring re-validation** at the new surface first. The twin's
   ranking authority (D0 §14.4: W₁ and midHot licensed, deepKE not) was
   measured at the old surface.
2. (v_c, τ, E₀) twin scan, MD finals, then a pooled N = 6000 battery as
   the new production reference.
3. Gate: the existing MD gate clauses (n̄, n₁) and the W₁/midHot
   ranking; KE₁ reported, not fitted (anti-circularity, as in
   `TIER2_MASS_SCENARIOS.md`).

Each scorer CSV carries the full observable vector (incl. KE columns)
and realized sizes.

## G6. Chord-map figure at the new standing point (after G5 C)

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

## Thesis framing implied

- Thesis table: bundle effect (A − reference), sampler effect (B − A),
  re-arbitration (C − B).
- The sampler is a stated modelling choice with its measured
  sensitivity (D0 §15.8), not a validated property.
- m_eff = 19 He is stated as the window-representative extraction mass
  (G3).
