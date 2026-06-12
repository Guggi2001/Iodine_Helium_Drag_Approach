# Drag-Model Port — Migration Log

**Purpose.** Detailed record of the drag-model port's implementation slices and
the Tier-0 investigation. This is the **reference detail** that `CLAUDE.md`
summarises: CLAUDE.md carries the compact current status and the live rules; this
log carries the per-slice delivery records and the full Tier-0 diagnosis history
(including the withdrawn readings, kept for the audit trail). Consult this when
the summary in CLAUDE.md is insufficient.

**Companion docs:** `DRAG_PORT_DESIGN_DECISIONS.md` (the frozen design),
`METHOD_B_trajectory_matching_extraction.md` (the active extraction method),
`TIER0_FINDINGS.md` (the Tier-0 verdict + diagnosis), `TIER0_SCRIPTS.md`,
`tier0_comparison_tasks_left.md`, and the per-slice specs
(`SLICE{1,2,3,4}_GOALS_*.md`, `slice_4.md`, `SLICE4_FIX_*.md`).

---

## Implementation slices — delivery records

### Slice 1 — complete

Delivered `physics/drag.py`: three pure, config-free, mass-free functions —
`drag_force(v, depth, …)` [amu·Å/ps²], `drag_gamma(v, depth, …)` [amu/ps,
closed form `g·(a+b·v²)`], `spatial_gate(depth, steepness)` [dimensionless
erf-complement] — plus the `DragCoefficients` bundle type (form tag,
coefficients, `extraction_mass_model` + `extraction_mass_amu` provenance).
`linear_cubic` realised; `linear_quadratic` / `threshold` / `power_law`
reserved and raise `NotImplementedError`. `physics/collisions.py` left intact
and importable (additive, parallel — not a deletion). Specs:
`SLICE1_GOALS_gated_drag_module.md`, `drag_module.md`.

Two upstream verifications were flagged and remain the user's to confirm at the
extraction source (a self-consistent refit cannot detect either):

- the stamped `extraction_mass_amu ≈ 202.954` is the mass the force balance
  *actually ran under*, not a relabelled value (an earlier literal was
  `179.912`);
- `drag_data.csv`'s `F_drag` column stores the **positive magnitude** (signed
  force balance is `F_drag = m_eff·a − F_C`, negative during the explosion).

Post-extraction empirical finding folded into the decisions doc: the
`power_law` exponent is `n ≈ +2` (not the anticipated `−2`), so the form is
regular at `v→0` and `drag_low_v_floor` is inert for the real coefficients
(retained only for a hypothetical `n<0` re-extraction).

### Slice 2 — complete

Delivered `physics/baoab.py`: the BAOAB operator-split ion-stage stepper
(decision §4.6, `make_ion_baoab_step`), replacing `velocity_verlet_step` for
the ion stage only. Key properties locked in:

- **Scheme B–A–O–A–B:** B/A are the baseline kick/drift (conservative force =
  Coulomb + droplet via `_ion_accel_fn`); O is the new physics — drag as
  multiplicative velocity damping `v ↦ e^(−γ·dt/m)·v` plus a dormant
  Langevin-noise site. Mass enters **only** here (in amu), in the O-step
  exponent and the energy bookkeeping.
- **Asymmetric γ-freeze (intentional, in the module docstring):** γ's velocity
  frozen at the O-step input velocity; γ's depth/gate at the current O-step
  position. Keeps the never-adds-energy dissipativity exact.
- **`_kick`/`_drift` extracted from `leapfrog.py`** to avoid duplicate physics —
  the one knowing touch to the frozen integrator, behaviour-preserving and
  self-verified by the anchor test.
- **Energy:** O-step returns `ΔE_dissip` in **amu·Å²/ps²** (4-tuple
  `(pos, vel, E_pot, ΔE_dissip)`); Slice 4 converts to eV. The noise-injection
  energy channel is a deliberate future signature bump, not a reserved slot.
- **`dt` is per-call** (`step(pos, vel, dt)`), mirroring `make_ion_step`; the
  per-step closure rebuild is driven by mass, not `dt`.
- **Anchor (killer) test:** at `γ=0`, noise off, BAOAB ≡ baseline
  `velocity_verlet_step` to round-off, with `acc_fn` call count asserted
  (`calls == 2·n == verlet_calls`) — proving integrator correctness and that
  the kick/drift extraction left the baseline unchanged.

Force-eval caching is deferred to Slice 4 (profiling-contingent) with one
recorded trap: cache the conservative *force* `F_cons`, never the
*acceleration* — under mass dynamics `a = F/m` goes stale at fixed position,
a bug invisible to the fixed-mass anchor test. Specs:
`SLICE2_GOALS_baoab_ion_stepper.md`, `baoab.md`.

### Slice 3 — complete

Delivered the declarative + validation layer on `config.py` / `presets.py`,
additively and with **no behavioral change** (the collision path still runs;
Slice 4 swaps it). Key properties locked in:

- **~18 drag fields** as named module-scope `Literal` aliases (house style,
  matching `CollisionMode = Literal[1,2,3]`), all with **inert** defaults
  (`mass_scenario=fixed`, `noise_form=none`, `drag_coefficients=None`) — *not*
  the design's "primary." A config left untouched runs the hard-sphere path
  unchanged. `mass_initial_amu` defaults to `m_eff_amu` (avoids a `None`-resolve
  branch). `m_eff_amu = 202.953908` (full precision, for the loader's exact
  match).
- **`check_drag_config(cfg)`** — separate function, called from `validate()`,
  no-ops when `drag_coefficients is None`. Runs (0) an **unconditional**
  `drag_form` typo-reject (recovers the runtime safety `Literal` gives up vs.
  `enum.Enum`) + `drag_form`↔`coeffs.form` agreement; (1) §3.3 per-form
  dissipativity (`linear_cubic`: `a>0`, turnover assert-and-skip while `b>0`,
  max-speed ceiling unsourced/recorded); (2) §6.5 mass↔coefficient consistency
  (`fixed`↔constant within ~8 amu; non-`fixed`↔time-resolved; inconsistent →
  refuse unless `allow_inconsistent_mass_pairing`). The non-`fixed` branch is
  exercised by synthetic construction in tests, not left untested-green.
- **`load_drag_coefficients(coeff_dir, *, expected_m_eff_amu)`** in the
  presets/config layer (keeps `physics/` I/O-free): content-validating, stamps
  `extraction_mass_amu` **from the JSON** (single source of truth), refuses on
  provenance mismatch (exact 1e-6) / missing / malformed. Case lives in the
  **presets** via `REFERENCE_DRAG_ROOT` — no `DragCase` enum, no helper.
- **Two distinct mass tolerances, kept separate:** loader's exact-match
  provenance identity vs. the guard's ~8 amu physics band
  (`_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU = 8.0`). Distinct constants, distinct
  tests.
- **Two drag-enabled presets** (`single_pulse_N2000_drag`,
  `single_pulse_N2000_18Angst_drag`) wire the real coefficients; existing
  presets unmutated; hard-sphere path bit-identical.

Specs: `SLICE3_GOALS_config_and_guard.md`, `config_and_preset.md`.

### Slice 4 — complete

Specified in `SLICE4_GOALS_ion_driver_rewiring.md`; delivered in `slice_4.md`
(+ `SLICE4_FIX_initial_mass_consistency.md`). Wired the BAOAB stepper into the
ion stage so a drag preset runs end-to-end at Tier 0 (deterministic, fixed
mass) — the first runnable drag trajectory, and the slice that **exercises** the
scoped collision-physics exception. Key properties locked in:

- **Parallel per-step function `baoab_propagation_step`** added to
  `ion_propagation_step.py` as a third sibling (alongside
  `neutral_propagation_step` / `ion_propagation_step`). `ion.py` dispatches
  **once per run** on `drag_coefficients is not None` → BAOAB; else → the
  existing collision step, **bit-identical and uncalled-for-drag**. The
  collision path (and its mass attachment) is not deleted or made dormant — it
  simply stays in the branch drag configs don't enter (the hard-sphere path is
  calibration data §6 depends on).
- **`ion.py` owns closure construction:** builds the `gamma_fn` + spatial gate
  and the BAOAB closure via `make_ion_baoab_step`, **rebuilt every step**
  (matching the `make_ion_step` pattern; Tier-1-ready though mass is fixed now).
  The gate assembly implements the §5.5 collapse — `density_proportional`
  (default) and `erf_tied` both build the erf-complement gate until a density
  profile exists. `baoab_propagation_step` receives the ready `step` and does
  thin per-step accounting.
- **Energy:** `ΔE_dissip` (amu·Å²/ps² from the stepper) → eV via the **baseline
  idiom** (amu→kg via `U`, ×100², ÷`EV`), accumulated into `E_dissip_eV`.
  Closure `E_kin+E_pot+E_dissip` is **tight** (Verlet-drift level — dissipation
  is exact, not approximate).
- **Checkpoint v5 retained** (v6 is Tier-1 mass dynamics):
  `E_mass_attach_defect_eV=0`, `temperature_diagnostic=NaN`,
  `number_of_collisions=0` under the drag branch.
- **`_check_drag_scope`** — drag-branch analog of `_check_scope`, rejects the
  out-of-Tier-0 envelope (`T_eff>0`, `mass_scenario≠fixed`,
  `drag_form≠linear_cubic`).
- **Deterministic smoke-run harness** — the "thorough debug" artifact: tiny
  fixed-seed run asserting finite trajectories, tight energy closure, monotone
  `E_dissip>0`, v5 checkpoint round-trip, scope-guard rejection. Explicitly
  **not** a TDDFT check.
- **Shared scaffolding (rule 1):** only behavior-preserving, physics-free lifts
  out of `ion_propagation_step` (depth, eV `E_kin`, possibly `E_pot`), each
  gated on `test_ion_propagation_step.py` staying green; write-your-own if a
  lift isn't clean.

Slice 4 scope fence — does **not** touch: the checkpoint schema (v5 retained);
mass dynamics; active noise; TDDFT validation/threshold-setting; new `SimConfig`
fields; `physics/` (`drag.py`/`baoab.py` consumed unchanged); the neutral stage.

---

## Tier-0 outcome — full record

### Tier-0 investigation and verdict (detail)

The Tier-0 comparison has **run** (`TIER0_COMPARISON_spec.md` →
`TIER0_FINDINGS.md`; scripts `TIER0_SCRIPTS.md`; status
`tier0_comparison_tasks_left.md`). It converged through several **withdrawn
readings** (full audit trail in the findings: "different regimes" → "frame
systematic" → "windowing + bubble-mode") onto a settled picture once the
reference was re-exported with real 3D per-atom velocities **and** positions:

- **18 Å — clean pass.** `linear_cubic` reproduces the trace (t\*-seeded distance
  0.30 Å, mean |v| RMSE ~0.09–0.16 Å/ps); its window stays in the clean radial
  regime throughout. Committed regression floor
  (`tests/test_tier0_drag_comparison.py`).
- **9 Å — the reference is genuinely NON-RADIAL (first-class finding).** Real 3D
  data shows atom 2 (the clean extraction atom) carries a sustained ~4 Å/ps
  transverse drift in a **radial↔transverse oscillation** (a y-peak follows each
  radial peak, ~twice across the window; the 2-D visualization had hidden this).
  At t\* the split is ~99.6% transverse. The extraction takes `|v|` and the MD
  projects it **radially** — a **defined modelling convention**, not a
  purely-radial drag law. The ~0.39 Å/ps same-smoothed residual is the
  central-force MD being structurally unable to carry the real transverse
  co-translation: a **model-dimensionality** statement, **not** a drag-form error,
  **not** a harness artifact (that was already removed), **not** the withdrawn
  frame story.

**Two method changes follow (2026-06-08):**

1. **Extraction A → B (`METHOD_B_trajectory_matching_extraction.md`).** Adopt
   **trajectory-matching calibration**: fit `{a,b}` by minimizing the
   forward-integrated in-window trajectory RMSE against the same-smoothed
   reference, rather than the Method-A direct `F_drag`-vs-`v` regression.
   Formalizes the hand-tuning observation; optimizes the observable that matters.
2. **Tier 0 repurposed: consistency → held-out generalization.** Under B the
   trajectory match *is* the fit objective, so the consistency-check framing is
   circular and **retired**. Tier 0's infrastructure (the `window=` parameter, the
   harnesses, the gate) survives but now scores **held-out** data: a held-out
   sub-window, the cross-case shared-form check (§3.6, the transport-physics
   signal), and the downstream observables (Tier 2 VMI, Tier 3). This role is
   *more* necessary under B (which can overfit) than the consistency check was.

**Active task: the Method-B extraction with held-out validation.** Mandatory
guard for **9 Å specifically**: trajectory-matching a radial-projected MD onto the
non-radial 9 Å reference risks the coefficients **absorbing the transverse
discrepancy** as a dimensionality fudge that fits perfectly and generalizes badly.
Held-out validation (held-out case / VMI) is **non-optional** for 9 Å; 18 Å
(genuinely radial) calibrates safely.

**Tier 1 is gated on the 9 Å Method-B fit surviving held-out validation** (not on
a frame re-extraction). Do not start Tier 1 (the `mass_rate_*` fields, the
`IonCheckpoint` v6 bump, mass-scenario A/B/biphasic) until then.
`EXTRACTION_FRAME_FIX_milestone.md` is **further demoted** — the non-radial
reality is handled by the radial-projection convention flag + held-out validation;
the He-field relative-velocity route is a contingency only if held-out validation
shows the radial-projection convention cannot generalize.

The Tier-0 *infrastructure* is reusable as-is. `linear_cubic` stands; the
coefficients become **Method-B-extracted** (a coefficient swap behind the same
interchangeable surface — `drag.py`/`baoab.py` consume them unchanged).

---

## Method-B extraction — delivery record (2026-06-11)

Implemented per the approved plan (slices B0–B8); user decisions locked at
planning: **3-parameter joint fit** `{a, b, E_bind}` on from-onset trajectory
RMSE (VMI purely held-out, never in the objective), **from-onset** integration
per evaluation (not t\*-seeded), and the 9 Å Method-A artifact **left as-is,
flagged provenance-broken** (below).

### Delivered

- **B0** — `18A/velocity_smoothed/cleaned_data_18_long.csv` renamed
  `cleaned_data_long.csv` (symmetric with 9 Å; was staged-uncommitted).
- **B1** — `DragCoefficients` gains optional `extraction_method`
  (`force_balance` default | `trajectory_matching`) and
  `effective_binding_energy_I_ion_eV` (None default); `load_drag_coefficients`
  requires the Method-B provenance keys (`extraction_mass_model`, binding,
  window, `reference_file`) when a JSON declares
  `extraction_method="trajectory_matching"`. Legacy Method-A files load
  byte-identically and are **not** re-stamped.
- **B2** — §6.5.1 guard: third arm in `check_drag_config` — stamped binding
  must equal `cfg.binding_energy_I_ion_eV` exactly (1e-9 eV wiring identity);
  an unstamped (Method-A) bundle is refused as "pairing not jointly
  validated"; `SimConfig.allow_unvalidated_binding_pairing=False` escape hatch
  downgrades to a loud `RuntimeWarning`. The two drag presets carry the hatch
  **transitionally** (documented in-code) until a validated Method-B bundle is
  wired in.
- **B3** — `i2_helium_md/extraction/trajectory_matching.py`: `CaseSetup`
  (one cached neutral checkpoint per case/seed; the drag ion path is RNG-free
  so the objective is bit-deterministic — verified by the DRY_RUN repeat
  assert), `evaluate_objective` (ensemble-mean |v2| vs `cleaned_SG` via
  `compare_speed_to_reference`, plus escape penalty
  `2.0·(1−escape_fraction)` A/ps), `fit_trajectory_matching` (normalized
  bounded Nelder-Mead, E_bind pre-scan, 3 starts probing the degeneracy
  ridge), `sensitivity_halfwidths`, `write_fit_parameters` (refuses trapped
  fits). Tests: `tests/test_extraction_trajectory_matching.py` (16, stub-MD;
  one real N=2/0.2 ps wiring eval).
- **B4** — `scripts/extraction/method_b_extraction.py` (USER SETTINGS;
  DRY_RUN determinism mode). Production fits written to
  `data/reference/drag/<case>/trajectory_matching/fit_parameters.json`
  (sibling of `linear_and_cubic/`; Method-A artifacts untouched).
- **B5/B7** — `scripts/extraction/method_b_cross_case_check.py` →
  `held_out_validation.json` per case, `"vmi": "pending"` stamped (Tier 1
  stays gated).

### Fit results (calibrated — NOT validated)

| case | a [amu/ps] | b [amu·ps/Å²] | E_bind [eV] | in-window \|v2\| RMSE [Å/ps] | escape |
|---|---|---|---|---|---|
| 18 Å | 1.456 (**pinned at the 0.1·a₀ lower bound**) | 3.316 | 0.0709 | 0.0968 over [4.54, 14.77] | 1.0 |
| 9 Å  | 7.50 (**ridge: starts span 3.68–7.50 at Δobjective ≈ 2e-4**) | 1.751 | 0.1543 | 0.0409 over [2.67, 14.08] | 1.0 |

Both fits escape fully with effective bindings far below the static 0.308 eV —
the §6.5.1 trap resolved as designed, and the E_bind values are consistent
with the TDDFT "escape with less than the static barrier" cross-check (c).
Seed sweep (5 neutral seeds): parameter std ≈ 0 — the HeDFT-comparison preset
(single initial position, fixed droplet radius, 0.4 K thermal velocities vs an
eV-scale explosion) is insensitive to the neutral seed, so the band is carried
by the RMSE-sensitivity half-widths instead (labeled
`seed_sweep_std_plus_rmse_sensitivity`; a sensitivity band, not a CI).

### Cross-case held-out verdicts — BOTH CASES FAIL the provisional bands

- **(i) 18 Å A↔B agreement: FAIL** — `a_B/a_A = 0.10` (band [0.5, 2.0]);
  `b_B/b_A = 1.62` passes.
- **(ii) residual asymmetry: PASS** — 9 Å/18 Å RMSE ratio 0.42 (≤ 3.0). Note
  the 9 Å fit is *better* than the clean-radial 18 Å — itself the §5 warning
  signature (a radial MD matching a non-radial reference too well).
- **(iii) cross-case E_bind: FAIL** — |0.154 − 0.071| = 0.083 eV > 0.05 eV.

**Honest reading (recorded, not threshold-tuned post hoc):** the dominant
finding is that **the linear coefficient `a` is weakly identified by the
full-window trajectory objective** — over both windows the speeds sit where
`b·v² ≫ a`, so the cubic term dominates: 18 Å's `a` ran to the imposed bound
with a sensitivity half-width equal to the parameter itself, and 9 Å's
multi-start ridge spans a factor ~2 in `a` at indistinguishable objectives.
The (i) failure is therefore at least partly an **identifiability** artifact
rather than pure overfit — but per the anti-laundering discipline the
provisional thresholds were not adjusted after seeing the data. The (iii)
E_bind disagreement is entangled with the same ridge (a↔E trade-off family).

### Consequences and open decisions (for the user)

- **B6 preset re-wiring NOT executed** (gated on the verdicts): the presets
  remain on legacy Method-A bundles behind the transitional §6.5.1 warning.
  Both `trajectory_matching/` bundles stand flagged **not-yet-usable for the
  radial MD** (METHOD_B §8 honest-failure path).
- Candidate next moves (decision needed): (a) **constrain `a`** (e.g. fix to
  the Method-A value and refit `{b, E_bind}`), (b) **cross-case shared-form
  joint refit** (shared `E_bind` and/or shared `a` across 9 Å + 18 Å — the
  strongest answer to the ridge), (c) accept a **pure-cubic reduced form**
  (drop `a`; new form behind the existing enum surface), or (d) proceed
  directly to the **VMI held-out tier** to let the production observable
  arbitrate among the degenerate pairs.
  **→ RESOLVED 2026-06-11 (same day, second session): (b) chosen, (c) folded
  in as a variant — see "Shared-form joint refit — decision record" below.**
- **Tier 1 stays gated**; `"vmi": "pending"` recorded in both
  `held_out_validation.json` files.

### 9 Å Method-A artifact: provenance-broken (user decision 2026-06-11)

The 9 Å `linear_and_cubic/fit_parameters.json` was hand-edited across recent
commits (`a`: 10.88 → 36.80 → 3.80; `b`/`a_err`/`b_err` stale) during the
hand-tuning experiments that motivated Method B. An lstsq refit on the
committed `drag_data.csv` scatter recovers `a = 36.7996` exactly (the a7fc41c
value) — so the authentic Method-A output is recoverable from git if ever
wanted. Decision: **leave as-is, flag broken**; the 9 Å A↔B agreement check is
dropped from the cross-case validation; the two 9 Å force-balance tests in
`tests/test_drag.py` (`TestForceBalanceReproduction[9A]`) are the **expected
known-red baseline** until the artifact is restored or superseded.

### Deferred

- The **Tier-0-docs rewrite** to the held-out framing (METHOD_B §8's separate
  doc pass) remains deferred — recorded here, not dropped.
- The **VMI held-out tier** (drag-enabled production droplet-distribution
  ensemble + Wasserstein vs `vmi_iplus_he.csv`) is the explicit follow-up
  phase and the ultimate arbiter. *(2026-06-11 update: recognized as
  unreachable until Tier 1 — see the decision record below.)*

---

## Shared-form joint refit — decision record (2026-06-11, second session)

Resolves the open decision above. Full specification with pre-registered
thresholds: `METHOD_B_trajectory_matching_extraction.md` §9. Implementation
awaits `[PROCEED TO IMPLEMENTATION]`.

### User decisions (locked)

- **Direction: (b) cross-case shared-form joint refit**, chosen because a
  successful shared fit speaks for a *general working principle* of the drag
  law (one size-independent transport law). Option (c) pure-cubic is folded in
  as an explicit variant rather than a separate path.
- **Sharing scope — both, primary + diagnostic:** primary hypothesis =
  fully-shared `{a, b, E_bind}` (3 parameters / 2 trajectories,
  over-constrained ⇒ falsifiable); diagnostic (only on primary fail) =
  shared `{a, b}` + per-case `E_bind`, localizing drag-vs-binding failure.
- **Stage 1 first:** fit 18 Å only (clean-radial), score the **untouched 9 Å
  prediction** (zero refit) — the only strictly held-out result obtainable
  before the joint fit consumes the cross-case axis. Recorded before Stage 2.
- **`a` handling:** free with lower bound **0**, plus an explicit `a ≡ 0`
  pure-cubic variant; indistinguishable objectives ⇒ pure-cubic recorded as
  the empirical conclusion on `a`'s identifiability.
- **Gate policy:** shared-form pass under the pre-registered §9.4 bands
  **ungates Tier 1**; the VMI tier moves to **post-Tier-1 final arbiter**.

### New structural finding — the VMI gate was circular

User clarification: **the VMI comparison cannot run in this phase** — the VMI
channels are mass-selected and final velocity depends on the mass history, so
a fixed-`m_eff` ensemble cannot be honestly scored against `vmi_iplus_he.csv`.
VMI therefore requires Tier 1 mass evolution — which was itself gated on
held-out survival, whose arbiter was VMI. **Cross-case is the only
falsification axis reachable now**; the two-stage design preserves held-out
content (Stage 1) before spending the axis (Stage 2), and the gate policy
above resolves the circularity on the record.

### Design parameters (detail in METHOD_B §9)

- Joint objective: per-case (RMSE + escape penalty) first, then the
  equal-weight mean across cases (prevents long-window dominance).
- Normalization anchors `{a0, b0}` from the **18 Å Method-A bundle only**
  (the 9 Å Method-A artifact is provenance-broken).
- Pre-registered bands (anchors recorded in §9.4): `S1_PRED_RMSE_MAX`
  0.45 Å/ps, `S2_RMSE_18A_MAX` 0.19 Å/ps, `S2_RMSE_9A_MAX` 0.45 Å/ps,
  escape = 1.0 everywhere, `T_a0` 0.005 Å/ps.
- §3.3 guard relaxation decided: `a ≥ 0` when `b > 0` for `linear_cubic`
  (still strictly dissipative; realizes the pure-cubic variant without a new
  form tag).
- Planned implementation surface (Phase B): extend
  `i2_helium_md/extraction/trajectory_matching.py` (`evaluate_joint_objective`,
  `fit_shared_trajectory_matching` with the three variants, Stage-1 prediction
  via the existing `evaluate_objective`); new
  `scripts/extraction/method_b_shared_refit.py` (USER SETTINGS + DRY_RUN
  pattern); artifacts under `data/reference/drag/shared/trajectory_matching/`;
  guard change in `config.py`; test extensions in
  `tests/test_extraction_trajectory_matching.py` and
  `tests/test_drag_config.py`. **Presets NOT re-wired** (gated on the verdict
  + a separate user decision).
- On fail: diagnostic variant localizes; per-case bundles stay flagged;
  `EXTRACTION_FRAME_FIX_milestone.md` remains the recorded contingency.

### §9 pre-run clarifications (2026-06-11, third session — locked before any run)

Six ambiguities surfaced in spec review and resolved by the user; folded into
METHOD_B §9 (first-runs rule: fixed before the first run, not after):

1. **Gating: only the Stage-2 bands carry verdict power.** The Stage-1 rows
   (`S1_PRED_RMSE_MAX`, `S1_PRED_ESCAPE`) are recorded for the held-out audit
   trail but are non-gating — a Stage-1 fail does not block the Tier-1 ungate
   if Stage 2 passes.
2. **Stage 1 reuses the delivered §8 18 Å single-case bundle** (`a = 1.456`,
   `b = 3.316`, `E_bind = 0.0709`; old `0.1·a₀` bound) — no fresh refit under
   the new `a ≥ 0` bound. Stage 1 is pure prediction (one 9 Å forward
   integration + scoring), so it has no optimizer and no variants.
3. **If the two Stage-2 variants differ beyond `T_a0`, record BOTH bundles**
   (`variant` provenance field disambiguates); the production-candidate
   choice is deferred to the preset-rewiring decision. If equivalent, it
   won't matter (pure-cubic conclusion recorded per §9.2).
4. **`T_a0` sign convention corrected:** Δ = objective(`a≡0`) −
   objective(`a`-free) ≤ 0.005 Å/ps (nested models ⇒ Δ ≥ 0 up to optimizer
   noise; the original wording had the difference inverted and would have
   been vacuously true).
5. **The Stage-1↔Stage-2 parameter-shift check is qualitative,
   recorded-only** — no pre-registered band, no gate power.
6. **The §9.4 bands supersede the §8 provisional cross-case bands.** The
   18 Å A↔B `a`-ratio check ([0.5, 2.0]) is not re-applied to the shared fit
   (it would auto-fail any `a → 0` fit; dead under the weak-`a` finding);
   the per-case `held_out_validation.json` updates score against §9.4 only.

Implementation still awaits `[PROCEED TO IMPLEMENTATION]`.
*(Given later the same session — delivery record below.)*

---

## Shared-form joint refit — delivery record (2026-06-11, third session): **PASS, Tier 1 ungated**

`[PROCEED TO IMPLEMENTATION]` given after the pre-run clarifications were
locked; implemented, tested, and **run to verdict the same day**. Summary
verdict: **the shared form generalizes — every §9.4 band passed — and the
law is empirically pure-cubic.** Full numbers in METHOD_B §9.7;
machine-readable record in
`data/reference/drag/shared/trajectory_matching/verdict.json`.

### Delivered (code)

- **§9.5 guard relaxation** (`config.py`): §3.3 `linear_cubic` arm now
  accepts `a >= 0` (refuses `a < 0` always, and the non-dissipative
  `a == 0, b <= 0` corner) — the pure-cubic variant stays `linear_cubic`
  with `a = 0`, no new form tag, no loader/enum churn.
- **`i2_helium_md/extraction/trajectory_matching.py`** extensions:
  `evaluate_joint_objective` (per-case RMSE + escape penalty first, then
  the equal-weight mean — prevents long-window dominance),
  `fit_shared_trajectory_matching` with the three §9.2 variants
  (`shared_3param` with `a` lower bound **0**, `shared_pure_cubic` with
  `a ≡ 0` fixed, `diagnostic_4param` with per-case `E_bind`; explicit
  `anchors=` keyword so the 9 Å setup's provenance-broken Method-A anchors
  can never condition the fit), `joint_sensitivity_halfwidths` (refuses
  per-case-binding results; `a`-half-width 0 by convention under
  pure-cubic), `write_shared_fit_parameters` (loader-compatible §9.6
  superset schema; refuses the diagnostic variant and trapped fits). The
  half-width factor scan was extracted into a shared
  `_objective_rise_halfwidths` core (rule 1; per-case
  `sensitivity_halfwidths` delegates, behavior preserved).
- **`scripts/extraction/method_b_shared_refit.py`** (USER SETTINGS +
  DRY_RUN pattern): Stage 1 → Stage 2 in order (held-out content recorded
  before the joint fit spends the axis), pre-registered §9.4 bands as
  named constants, Stage-2-only gating, corrected-`T_a0` classification,
  qualitative parameter-shift record, diagnostic only on primary failure,
  both variant bundles written, per-case `held_out_validation.json`
  updated in place (provisional-§8-band records left intact, superseded).
- **Tests:** `test_drag_config.py` dissipativity class extended (negative-`a`
  refused, `a = 0, b > 0` passes, `a = 0, b <= 0` refused);
  `test_extraction_trajectory_matching.py` +19 (joint objective equal-weight
  mean / key mismatch; shared fit recovers a shared stub law; `a`-bound-0
  reachability; pure-cubic fixes `a = 0`; diagnostic localizes per-case
  bindings; provenance/variant/single-case refusals; joint sensitivity;
  shared writer loads through the Slice-3 loader + refusal paths; committed
  shared-bundle artifact checks). Suite: **641 passed, 2 failed** — the two
  failures are the recorded known-red 9 Å force-balance baseline, unchanged.

### Run record (production, 2026-06-11)

N=50, seed 20260604, ion 20 ps @ 0.01 ps; anchors `a0=14.5556, b0=2.0534`
(18 Å Method-A only); DRY_RUN determinism verified (bitwise-identical joint
evals, 1.57 s each = 2 ion runs).

| check | value | band | verdict |
|---|---|---|---|
| Stage-1 9 Å prediction RMSE | **0.2685 Å/ps** | ≤ 0.45 | PASS (recorded, non-gating) |
| Stage-1 9 Å escape | 1.0 | = 1.0 | PASS (recorded, non-gating) |
| Stage-2 `shared_3param` 18 Å RMSE | **0.1345 Å/ps** | ≤ 0.19 | PASS (gating) |
| Stage-2 `shared_3param` 9 Å RMSE | **0.1240 Å/ps** | ≤ 0.45 | PASS (gating) |
| Stage-2 escape, both cases | 1.0 / 1.0 | = 1.0 | PASS (gating) |
| `T_a0` = obj(`a≡0`) − obj(`a`-free) | −0.000000 Å/ps | ≤ 0.005 | EQUIVALENT |

Best shared fit: `a = 0.0002` (ran to the 0 bound; all three multi-starts
agree — the §8 `a`↔`E_bind` ridge resolved at `a → 0`), `b = 2.516
amu·ps/Å²`, `E_bind = 0.1168 eV` (between the per-case 0.071/0.154);
pure-cubic variant indistinguishable (`b = 2.515`, same `E_bind`).
**Empirical conclusion on the §8 weak-`a` finding: the linear term carries
no information over these windows — the transport law is effectively
`γ = g·b·v²`.** Stage-1→Stage-2 shift recorded (Δa −1.456, Δb −0.800,
ΔE +0.046 eV; qualitative, non-gating). Diagnostic variant not needed.

### Consequences

- **Tier 1 UNGATES** (gate policy §9.1); the **VMI tier moves to
  post-Tier-1 final arbiter** (still stamped `pending`).
- **Presets NOT re-wired** — the production-candidate choice between the
  two equivalent recorded variants (`shared_3param` at the bound vs.
  `shared_pure_cubic`) plus the re-wiring itself await a separate user
  decision; the Method-A bundles and the transitional §6.5.1 hatch remain.
  *(→ both RESOLVED 2026-06-12 — decision record below.)*
- The per-case §8 bundles stay flagged as before (the shared bundle is the
  candidate, not them); the 9 Å Method-A artifact remains
  provenance-broken/known-red.
- Still deferred: the Tier-0-docs rewrite to the held-out framing.

---

## Production candidate + alternative-form phase — decision record (2026-06-12)

Resolves the two open §9 frontier decisions and defines the next phase.
Full specification: `METHOD_B_trajectory_matching_extraction.md` §10.
Implementation awaits `[PROCEED TO IMPLEMENTATION]`.

### User decisions (locked)

- **Production candidate = `shared_pure_cubic`.** The §9.7 variants are
  `T_a0`-equivalent; `a = 0` exactly is the honest encoding of the
  identifiability conclusion (`shared_3param`'s `a = 0.0002` is bound
  noise, not a measured linear coefficient).
- **Preset re-wiring APPROVED** (execution pending, after the doc pass +
  milestone commit): the two drag presets move to the `shared_pure_cubic`
  bundle; its jointly-calibrated stamped binding means the transitional
  §6.5.1 `allow_unvalidated_binding_pairing` hatch comes off the presets
  and the guard runs at full strength.
- **Milestone commit APPROVED**: the §9 implementation + shared bundles +
  this documentation pass go in as one coherent commit on
  `drag_implementation` (the per-case §8 state is no longer separable).
- **Next implementation phase — alternative drag-form discrimination
  (METHOD_B §10):** realize the reserved `power_law` and `linear_quadratic`
  forms (incl. the pure-quadratic `a ≡ 0` variant) and run each family
  through the same shared-form joint-refit machinery, compared against the
  pure-cubic incumbent under the reused §9.4 bands plus a pre-registered
  form-equivalence threshold `T_form` (and `power_law` `n` bounds) to be
  locked before the first run. `threshold` stays reserved (out of this
  phase's scope). Motivation: the Method-A power-law export's `n ≈ +2`
  vs §9.7's pure-cubic `n = 3` — can the trajectory objective discriminate
  the exponent at all? The free-`n` `power_law` fit nests both incumbents
  and its `n̂` is the direct identifiability measurement.
- **Ordering:** form phase precedes Tier 1 (still ungated-but-not-started);
  the form comparisons are model selection on seen data (the cross-case
  axis was spent by §9 Stage 2) — recorded honestly as ranking, not fresh
  held-out validation; VMI stays the post-Tier-1 final arbiter.

### Milestone commit + preset re-wiring — delivery record (2026-06-12)

- **Milestone commit `696bfa7`** (29 files): the full Method-B arc (per-case
  §8 + shared §9 implementation, bundles, tests) + the §10 documentation
  pass + the upstream `Drag_function` 18 Å export. **Pre-commit finding:**
  a 17:14 post-verdict experiment (2026-06-11) had overwritten both
  Method-A `linear_and_cubic/fit_parameters.json` files with a shared-fit
  multi-start row (the *third* start, not the chosen best) and regenerated
  both Tier-0 `md_mean_trajectory.csv` regression references — contradicting
  the recorded state and turning the 18 Å force-balance tests red. **User
  decision: restored all four reference artifacts to HEAD before
  committing** (the experiment values remain available in the shared
  bundles / `verdict.json`); suite verified back at 641 passed /
  2 known-red.
- **Preset re-wiring EXECUTED** (`[PROCEED TO IMPLEMENTATION]` given
  2026-06-12): both drag presets (`single_pulse_N2000_drag`,
  `single_pulse_N2000_18Angst_drag`) now load the one shared
  `shared_pure_cubic` bundle (`_SHARED_DRAG_BUNDLE_DIR` in `presets.py`)
  and set `binding_energy_I_ion_eV` **from the loaded bundle's stamp**
  (single source of truth — the §6.5.1 exact-identity guard passes by
  construction); the transitional `allow_unvalidated_binding_pairing=True`
  hatch is **removed** from the presets (B2's transitional state retired).
  Tests: `test_transitional_drag_presets_warn_not_refuse` flipped to
  `test_rewired_drag_presets_validate_silently` (warnings escalate to
  errors); `TestDragPresets` pins the re-wired bundle identity
  (`extraction_method`, `a == 0.0`, binding identity); the redundant
  hatch override removed from `test_run_directory.py`. Full suite:
  **641 passed, 2 known-red** (unchanged baseline), transitional
  RuntimeWarnings gone. The legacy Method-A bundles stay on disk as the
  frozen independent cross-reference (loader contract unchanged).

