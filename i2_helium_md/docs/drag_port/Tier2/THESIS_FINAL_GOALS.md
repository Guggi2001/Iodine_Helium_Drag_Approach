# THESIS_FINAL_GOALS — agreed changes and tests before thesis completion

Settled in the final discussion round (2026-10-02). This document holds
only agreed goals. Results go to the migration log and D0 as usual.

**Status (2026-10-05):** G4 **executed** (G4.8). The adoption decision
is superseded by **G4.9** (18 Å-only refit; 9 Å held out as
non-radial), **executed 2026-10-05: winner `cubic` (`srho3p5A18A_cubic`,
b 2.2286, E_bind 0.1214 eV)**; adoption for G5 is gated by **G4.10**
(force-balance cross-check), **executed 2026-10-05: 18 Å b_A = 1.618
outside the band [1.803, 2.654], read as INCONCLUSIVE** (references not
derivative-grade). **Post-hoc (June Method-A pipeline + corrections):
cubic b 18 Å 2.449, 9 Å 2.634 — both inside the band**, 9 Å at the edge
(G4.10 status). **G4.11** (cross-case drag-law test), **executed
2026-10-05: only the cubic is consistent (s_X 1.182, s_M 1.099) → cubic
supported by the TDDFT data** (marginal: C-X at the band edge); linear and
quadratic fail C-X. The adoption gate for `srho3p5A18A_cubic` in G5 is
cleared on the rule, but the reading is **contested** by the user's
Tier-0 9 Å runs (cubic over-dissipates the radial ion; G4.11 post-hoc
caveat) — **adoption on hold. G4.12** (TDDFT-tied linear + quadratic
ridge vs the experiment, uncapped): **steps 1–2 executed 2026-10-05** —
the whole ridge fits 18 Å; the 9 Å exit bound passes every c ≤ 5 and
fails only the pure quadratic. Step 3 (twin): all 11 land after the frozen refinement (landing candidate), but h405 dominates at matched midHot. **G4.13** (MD battery, c = 2.0, executed 2026-10-06): dominated **at the selected cell t64e310 only** (pooled n₁ 0.186 misses the gate, KE₁ ties h405, W₁ 1.10 vs 0.755); the grid could not resolve the landing (correction, user-caught). **Next: G4.14** fine τ-line sweep, draft under discussion.
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

### G4.14 Fine τ-line sweep for ridge-law landings (DRAFT 2026-10-06 — agreed in outline, under discussion; nothing run)

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

**Design (τ is the primary axis; E₀ is slaved to the gate along it).**

1. **Step 0 — calibrated twin, zero MD.** c ∈ {1.2, 1.6, 2.0, 2.5}
   (cached chords; τ is an exact rescale K → K·6.55/τ, so free) ×
   τ 3.6–7.2 step 0.1 × E₀ step 0.0025, with the offsets above. Output
   per c: the gate ridge E₀*(τ) and predicted MD KE₁ / midHot / n̄ along
   it. The offsets are measured at c = 2.0 only; other c need one MD
   anchor each in step 1.
2. **Step 1 — MD fine τ line (fitted, not replicated).** c = 2.0:
   τ 3.6–6.4 step 0.2 (15 τ) × two E₀ bracketing E₀*(τ), one shared seed
   (CRN) → 30 runs (≈ 70 min). n₁, n̄, W₁, KE₁ are fitted as smooth
   surfaces in (τ, E₀) across all runs, so the seed noise averages over
   the grid instead of per-cell replicates; a second seed (s2) follows
   if the fit residuals exceed the gate margin. Other c: τ step 0.4,
   only where step 0 puts their front (plus the anchor cell).
3. **Step 2 — confirmation.** The fitted best τ on the ridge, per c
   that lands, at 6 seeds → pooled N = 6000 vs `pooled_new6000`.

**Landing / ranking (to freeze).** Landing ⇔ pooled hard gate passes
with an n₁ margin ≥ one pooled seed sd. Ranking vs h405 on KE₁, W₁,
midHot — separate from landing.

**Prediction (to freeze before MD).** W₁ falls as τ falls; gated
landings exist along the ridge; W₁ ≤ 0.80 needs τ ≲ 4.4, where KE₁ at
c = 2.0 drops to ≈ 0.53, so lower c is needed to recover KE₁ — the
(c, τ) trade the line measures.

**Open for discussion (user).** c set; the selection / landing margin;
the fit model (and its residual criterion); seeds in step 1; whether
deepKE / χ² enter the ranking.

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
