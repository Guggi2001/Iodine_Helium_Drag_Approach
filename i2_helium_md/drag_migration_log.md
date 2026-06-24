# Drag-Model Port — Migration Log

**Purpose.** This log is the **single home for the drag-model port's decision
history** — the per-slice delivery records, the dated decision records, and the
full Tier-0 diagnosis history including every withdrawn reading. The companion
docs (`CLAUDE.md`, `METHOD_B_…`, `TIER0_FINDINGS.md`,
`DRAG_PORT_DESIGN_DECISIONS.md`) carry only the **present state** plus the live
rules/specs and point here for history. Consult this when you need to know *how*
a decision was reached or *what was tried and withdrawn*.

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

### Withdrawn Tier-0 readings (audit trail)

The 9 Å diagnosis converged through several readings that were each adopted then
withdrawn. Recorded here (the one home for them) so a later reader knows what was
tried and why it was wrong; `TIER0_FINDINGS.md` keeps only the settled picture.

1. **"Frame systematic" (withdrawn).** The 9 Å mismatch was first read as a
   lab-vs-relative velocity-frame / COM-drift error and spawned
   `EXTRACTION_FRAME_FIX_milestone.md`. **Wrong:** atom 1's sideways motion is a
   TDDFT artifact *discarded at extraction* (γ was fit to the single clean atom 2,
   not the `½(v₁+v₂)` COM the story rested on), and the break is local/directional,
   not the global bias a frame error would produce.
2. **"Windowing + bubble-mode" (partly right, superseded).** Next read as three
   non-law effects: the artifact atom in the mean comparison, the window running
   past atom 2's ~6 ps drift onset, and the 1.2 ps bubble oscillation the
   extraction denoised away. The same-smoothed comparison confirmed bubble-mode as
   the dominant *raw* contributor (residual 0.88 → ~0.40) — but it did **not**
   collapse to the 18 Å-class ~0.09, so a real residual survived.
3. **"Different-regime / over-damping" (refined, not final).** The surviving
   ~0.40 was briefly read as a near-constant ~10% over-damping (prime suspect: the
   9 Å clean-window `a` nearly doubling). Superseded once real 3D positions showed
   the residual is dimensional, not a magnitude miscalibration.
4. **Settled — 9 Å reference is genuinely non-radial (current).** Re-exporting the
   reference with real 3D per-atom velocities *and* positions showed atom 2 carries
   a sustained ~4 Å/ps transverse co-translation; the central-force MD cannot
   represent it, so the ~0.39 Å/ps residual is a **model-dimensionality** statement,
   not a drag-form/frame/magnitude error. This is the first-class finding in
   `TIER0_FINDINGS.md`.

### Tier-0 investigation and verdict (detail)

The Tier-0 comparison has **run** (`TIER0_COMPARISON_spec.md` →
`TIER0_FINDINGS.md`; scripts `TIER0_SCRIPTS.md`; status
`tier0_comparison_tasks_left.md`). It converged through the withdrawn readings
above onto a settled picture once the reference was re-exported with real 3D
per-atom velocities **and** positions:

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
*(→ RESTORED 2026-06-12 — decision record below; the known-red baseline is
cleared.)*

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

---

## Method-A future role + 9 Å artifact restoration — decision record (2026-06-12)

Resolves the "is a Method-A fit still needed?" discussion. User decisions
locked via review questions; the artifact restoration was executed the same
session.

### The Method-A future-role verdict

**No new Method-A (force-balance) fits are needed.** The §10 form phase
fits the new drag forms by trajectory matching (the design-doc §3.7
"outstanding `linear_quadratic` fit pass" route is superseded), and Tier 1
reads the time-resolved `m(t)` shell trajectory directly off the TDDFT
reference (design doc §2/§6.6) — the time-resolved re-extraction stays an
optional refinement. The existing Method-A artifacts retain three roles:

1. **Frozen regression-test fixtures** (never expires):
   `tests/test_drag.py::TestForceBalanceReproduction` pins
   `physics/drag.py` against the committed `drag_data.csv` scatter — the
   guarantee that the module IS the extracted law.
2. **Anchor *values*** for trajectory-matching optimizers — numerical
   conditioning only, now decoupled (below).
3. **Recorded evidence**: the `power/` exports' `n ≈ +2` is the §10
   motivation and revives once as the A↔B cross-reference **on the
   exponent** (trajectory-matched `n̂` vs Method-A's `n ≈ +2`). After §10
   the METHOD_B §7 A↔B cross-reference role is fully retired (for
   `linear_cubic` it is already played out: the `a`-ratio check is dead
   under weak-`a`; the surviving content was `b` agreeing within ~20%).

### Anchors → pre-registered constants (user decision)

The §10 form phase sources its optimizer anchors as **named constants with
provenance in the script's USER SETTINGS**, read off the Method-A bundles
once and locked at the §10 pre-run session (alongside `T_form` and the `n`
bounds). This removes the last *runtime* Method-A dependency and resolves
the post-re-wiring anchor trap (`build_case_setup`'s `setup.a0/b0` reads
`a0 = 0` from the shared bundle since the re-wiring — documented in
`docs/extraction/trajectory_matching_module.md` and
`scripts/extraction/scripts_explained.md`). Docs-only now; no code change
until the form phase.

### 9 Å Method-A artifact RESTORED (supersedes the 2026-06-11 leave-as-is decision)

With Method-A's role reduced to frozen-fixture + anchor-values + evidence,
a broken fixture was the one live inconsistency. Restoration executed:

- `git show a7fc41c` recovered the candidate file; **cross-verified** by an
  independent lstsq refit on the committed `drag_data.csv` trusted interior
  (the test's exact recipe): refit `a = 36.79962483`, `b = 1.38728042` —
  both match the a7fc41c values to ~1e-9 (test tolerance 1e-4), residual
  band 3.3% (test band 7%). The a7fc41c file is the authentic,
  self-consistent Method-A output; the broken state differed **only in
  `a`** (3.7996 — a dropped leading digit of 36.7996; `b`/errors/window
  were already authentic).
- `data/reference/drag/9A/linear_and_cubic/fit_parameters.json` updated:
  `a` 3.79962487165992 → 36.79962487165992. Non-fit fields unchanged
  (`t_start`/`t_end` feed `test_tier0_drag_comparison._read_window`).
- **The known-red baseline is CLEARED: full suite 643 passed, 0 failed**
  (was 641/2 since the Method-B delivery). The 9 Å A↔B-check drop and the
  §9.3 18 Å-only anchoring remain recorded history (correct decisions at
  the time — the artifact really was broken when they were made).

---

## §10 form-phase pre-run constants — decision record (2026-06-12, pre-run session)

The §10.6 pre-run lock executed: all open §10.4 design questions resolved
by the user and the pre-registered constants fixed **before any form-phase
run** (first-runs rule). Authoritative table: METHOD_B §10.4.1. Docs-only —
implementation still awaits `[PROCEED TO IMPLEMENTATION]`.

### User decisions (locked)

- **`T_form` = two-threshold scheme** (supersedes the single-0.005
  candidate): Δ = (family **best-variant** objective) − (incumbent
  `shared_pure_cubic` 0.1293 Å/ps), classified as genuinely-better /
  marginally-better / equivalent / worse via `T_FORM_EQUIV_APS = 0.005`
  (optimizer-noise/`T_a0` scale) and `T_FORM_BETTER_APS = 0.013` (the
  10%-sensitivity scale of the incumbent objective). Escalation to the
  user only beyond −0.013; the −0.013 < Δ < −0.005 zone is recorded as
  "marginally better, not escalation-worthy". Rationale: keeps a
  within-flatness improvement from entering the verdict record stamped
  "genuinely better".
- **`power_law` optimizer parameterization = pivot `(γ_ref, n)`**,
  `γ_ref = C·v_ref^(n−1)` at the locked pivot `V_REF_APS = 3.0 Å/ps`
  (inside both in-window speed ranges) — axis-aligns the
  `log C ≈ const − n·log v̄` matching ridge so `n̂`'s sensitivity
  half-width is meaningful. Optimizer-internal only; bundles stamp raw
  `{C, n}` (loader contract unchanged).
- **`n` bounds locked `[1, 4]`** (nesting points 2 and 3 interior;
  `n ≥ 1` keeps `γ` finite at `v = 0`, `drag_low_v_floor` stays inert).
- **Anchors 18 Å-only, as pre-registered constants:**
  `a0 = 14.555626399148123`, `b0 = 2.0534044239692157`
  (`18A/linear_and_cubic/`), `c0 = 11.016050300970692`
  (`18A/quadratic/`, the `d1ca029` Method-A quadratic export),
  `C0 = 10.36139380949775`, `n0 = 2.0557526931077588` (`18A/power/`).
  18 Å-only is now a convention choice (the 9 Å artifact is restored),
  kept for consistency with the §9.3 record.
- **Conventions:** Stage-1 analog per family runs the family's full
  variant (`lq_shared_3param` / `pl_shared_3param`), scored against the
  same 0.45 Å/ps number for comparability with pure-cubic's 0.2685
  (recorded, non-gating); the `T_form` comparison uses each family's
  best variant.

### Interpretive context (recorded pre-run)

The in-window smoothed-reference speed ranges are narrow — 18 Å
2.54–3.02 Å/ps (max/min 1.19), 9 Å 2.83–4.95 Å/ps (1.75) — so the
exponent leverage comes mainly from the from-onset transient and the
**cross-case speed-scale difference** (one shared law must serve windows
at ~2.75 vs ~3.9 Å/ps characteristic speed), not from in-window shape. A
single-case fit could barely discriminate `v²` from `v³`; the shared
machinery is what makes the question answerable. A large `n̂` half-width
(the degeneracy outcome) is therefore a live expectation and is already
a recorded-finding path in §10.2, not a failure.

## §10 form-phase — delivery record (2026-06-14): **incumbent `shared_pure_cubic` CONFIRMED**

The four-slice arc's §10 alternative-form discrimination phase implemented and
run. The exponent is settled empirically: the full-window trajectory objective
**does** discriminate the drag exponent, and it picks `n ≈ 3` (pure-cubic),
**not** the Method-A `n ≈ 2.06`.

### Delivered (code)

- **`physics/drag.py`** — `linear_quadratic` (`F = g·(a·v + c·|v|·v)`,
  `γ = g·(a + c·|v|)`) and `power_law` (`F = g·C·|v|ⁿ`, `γ = g·C·|v|ⁿ⁻¹`,
  numpy `0.0**0.0 == 1.0` realizes the `n=1` limit) realized behind the dispatch
  branch (Slice-1 promise kept — no signature change). `_raise_unrealised_form`
  is now `threshold`-only; `REALIZED_FORMS = (linear_cubic, linear_quadratic,
  power_law)` is the single source for loader + scope guard; module docstring
  carries the per-form governing equations, §10.3 dimensional analysis, and the
  nesting identities. **`power_law` coefficient key is `"C"`** (resolves the
  spec-vs-code discrepancy; §10.3/§10.5 stamp raw `{C, n}`; the historical
  Method-A `18A/power/fit_parameters.json` keeps its legacy `"gamma"` key —
  frozen evidence read once for `C0`, never loaded through the loader).
- **`config.py`** — §3.3 guard arms (§10.3): `linear_quadratic` `a ≥ 0`,
  `c ≥ 0`, `a + c > 0` (pure-quadratic corner `a=0, c>0`); `power_law` `C > 0`,
  `n ≥ 1` (n<1 refused → `drag_low_v_floor` stays inert).
- **`presets.py`** — form-generic `load_drag_coefficients` (optional JSON
  `"form"` key; absent ⇒ legacy `linear_cubic`, byte-identical; per-form
  required keys from `_REQUIRED_COEFF_KEYS`); `_FIT_PARAM_REQUIRED_KEYS` now a
  per-form dict.
- **`simulation/ion_propagation_step.py`** — `_check_drag_scope` relaxed to the
  realized form set `{linear_cubic, linear_quadratic, power_law}` (`threshold`
  still refused; T_eff>0 / mass_scenario≠fixed arms unchanged). `ion.py` needed
  no change — its `gamma_fn = partial(drag_gamma, …)` was already form-generic.
- **`extraction/trajectory_matching.py`** — the §10 layer beside the frozen §9
  API: `_run_and_score` core (extracted from `evaluate_objective`, which now
  delegates; bitwise-equal), the variant tags + `FORM_VARIANTS`, `_C_to_pivot`/
  `_pivot_to_C`, `FormObjectiveResult`/`JointFormObjectiveResult`/
  `SharedFormTrajectoryMatchingFit`, `evaluate_form_objective` /
  `evaluate_joint_form_objective`, `fit_shared_form_trajectory_matching` (≥2
  cases) + `fit_form_trajectory_matching` (single-case Stage-1 analog),
  `form_sensitivity_halfwidths` (lq scans `(a,c,E)`; pl scans the PIVOT
  `(γ_ref,n,E)`), and `write_shared_form_fit_parameters` (refuses single-case /
  trapped / unfilled-uncertainty; pl bundles add a `pivot` block).
- **Drivers** — `scripts/extraction/form_phase_common.py` (LOCKED §10.4.1
  constants, reused §9.4 band check, two-threshold `T_form` classifier,
  combined-verdict read-modify-write) + two siblings
  `method_b_form_refit_linear_quadratic.py` / `method_b_form_refit_power_law.py`
  (USER SETTINGS + DRY_RUN pattern; anchors passed as pre-registered constants,
  **never** `setup.a0/b0`).
- **Tests** — `TestFormDispatch` (threshold-only NIE), `TestFormPhaseFamilies`,
  `TestNestingIdentities`, `TestFormPhaseDissipativityGuard`,
  `TestLoaderFormGeneric`, per-family BAOAB smoke runs, and the extraction
  `TestPivotTransform`/`TestFormObjective`/`TestFormFit`/`TestFormSensitivity`/
  `TestWriteSharedFormFitParameters`. Full suite **705 passed / 0 failed**.

### Session decisions (user, 2026-06-12 — bound the remaining work)

1. Driver structure = sibling script per family (the §9 `method_b_shared_refit.py`
   stays frozen).
2. Multi-start = implementation freedom (lq mirrors §9's 3 starts; pl uses 4
   starts probing the exponent at n0, 2.0, 3.0).
3. Stage-1-analog 18 Å-only fits = record-only JSON, no loadable bundle (the
   writer refuses single-case fits).
4. `_check_drag_scope` relaxed to the realized form set; `threshold` still
   refused.
5. `power_law` coefficient key renamed `"gamma"` → `"C"` (spec-vs-code gap).

### Run record (production, 2026-06-14; first-runs rule honoured)

N=50, seed 20260604, 20 ps @ 0.01 ps; anchors the §10.4.1 LOCKED 18 Å-only
constants (`a0=14.5556`, `c0=11.0161`, `C0=10.3614`, `n0=2.0558`, `V_REF=3.0`);
joint objective bitwise-deterministic (DRY_RUN verified both siblings).
Incumbent objective `0.1292649398514104 Å/ps`.

- **`power_law` (`pl_shared_3param`, free `n`) — EQUIVALENT → pure-cubic
  confirmed.** All 4 starts → `n̂ = 2.927` (C ≈ 2.835, E_bind ≈ 0.113 eV),
  objective **0.129171**, Δ = **−0.0001** (within ±0.005); bands PASS
  (18 Å 0.1305, 9 Å 0.1278, escape 1.0/1.0); `n_err` half-width **0.279**.
  Stage-1 analog 9 Å predict 0.411 ≤ 0.45 PASS. The free-exponent fit
  **independently recovers `n ≈ 3`**, not Method-A's `n ≈ 2.06`.
- **`linear_quadratic` (forced `n = 2`) — WORSE → rejected-by-objective.** Both
  variants pass the §9.4 bands and collapse to the pure-quadratic corner
  (`a → 0`, `T_a0`-analog −1e-6 EQUIVALENT, `c ≈ 12.8`, E_bind ≈ 0.048 eV), but
  objective **0.16315**, Δ = **+0.0339** (past +0.005) → recorded
  rejected-by-trajectory-objective. Stage-1 analog 9 Å predict **0.699 > 0.45
  FAIL** (non-gating) — does not generalize the way pure-cubic did (0.2685).
- **Verdict: no family beats `shared_pure_cubic`; no escalation.** The §10.2
  exponent tension is resolved in favour of `n = 3`; the trajectory objective
  discriminates the exponent (unlike the linear term, §9). VMI (post-Tier-1)
  is the final external check. **Presets stay on `shared_pure_cubic` (NOT
  re-wired by this phase). Tier-1 start is now on the table.**
- **Artifacts** under `data/reference/drag/shared/trajectory_matching/`:
  `{lq_shared_3param,lq_shared_pure_quadratic,pl_shared_3param}/fit_parameters.json`
  (all load through `load_drag_coefficients`), `stage1_analog_*.json`,
  `verdict_{linear_quadratic,power_law}.json`, `form_comparison_verdict.json`.

## §10.8 per-case alternative-form fits — delivery record (2026-06-14): **DIAGNOSTIC, presets unchanged**

Filled the missing diagonal: per-case (9 Å-only, 18 Å-only) fits for
`linear_quadratic` / `power_law` (§10.7 had only the shared joint refits;
`linear_cubic` already had §8 per-case bundles). **Strictly documentary** (user
decision 2026-06-14): no outcome reopens the form question; the incumbent stays
`shared_pure_cubic`; these are entirely before any Tier-1 start.

### Decisions (user, 2026-06-14)

1. Artifact = loadable bundles (one per `(case, variant)`); both `lq` variants +
   the free-`n` `pl` ⇒ 6 bundles.
2. Strictly documentary — never preset-wired; honesty flags load-bearing.
3. Sibling driver reusing `form_phase_common`, **not** a flag on the
   `linear_cubic`-wired §8 `method_b_extraction.py`.

### Delivered (code)

- **`extraction/trajectory_matching.py`** — new `write_per_case_form_fit_parameters`
  (form-generic, single-case; mirrors both existing writers' refusals — trapped /
  unfilled uncertainty — and adds the §10.8 honesty flags
  `per_case_calibrated_not_validated`, `cross_case_axis_not_applied`; `pl` stamps
  the pivot block + a note that `C_err` is the fixed-`n` partial band). The shared
  writer's single-case refusal is **left intact** (the two are distinct entry
  points). Reuses `fit_form_trajectory_matching` (single-case) and
  `form_sensitivity_halfwidths`.
- **`scripts/extraction/method_b_per_case_form.py`** — sibling driver; builds each
  case setup once, loops `(case, variant)`, abort-and-skip on trapped (recorded,
  not a failure), sensitivity-only band (seed sweep omitted per §8), writes to
  `data/reference/drag/<case>/trajectory_matching/<variant>/` (extends the §8
  per-case location with a variant subdir; no collision with the flat
  `linear_cubic` per-case file).
- **Tests** — `TestPerCaseFormBundleWriter` (round-trip via `load_drag_coefficients`;
  honesty + transverse flags; `pl` pivot block; refuses multi-case / trapped /
  unfilled). Full suite **712 passed / 0 failed**.

### Run record (production, 2026-06-14)

N=50, seed 20260604, 20 ps @ 0.01 ps; §10.4.1 LOCKED anchors (conditioning-only);
DRY_RUN determinism verified both cases. **6/6 bundles written, escape 1.0, 0
trapped skips.**

- **18 Å** — `lq_3param` `a=16.63±1.59, c=6.13±0.59` (RMSE 0.0994); `lq_pure_quad`
  `c=11.27±0.55` (RMSE 0.0993 — equal, so the `a=16.6` is the §10.4.1 weak-`a`
  ridge, globally degenerate with the `a=0` corner); `pl_3param` **`n̂=4.000±1.333`
  RAILED to the n=4 bound** (exponent unidentified within the 18 Å curve).
- **9 Å** (transverse-contaminated) — `lq_3param` **collapses `a→0`**, `c=10.32±0.10`
  (RMSE 0.0446); `pl_3param` **`n̂=2.651±0.026`** (exponent identified *tightly*
  within the single curve, `C=3.59`).
- **Interpretation.** Per-case exponent leverage is **case-asymmetric** — it
  *refines* the §10.4.1 "unidentified within a single curve" expectation: 9 Å's
  broader in-window speed range pins `n` tightly; 18 Å's narrow range rails to the
  bound. The shared fit's `n̂=2.927≈3` arises from the cross-case speed-scale
  combination (9 Å→2.65, 18 Å→unidentified) — neither single case lands at 3.
  **No reopening of the incumbent; `shared_pure_cubic` stands; presets unchanged;
  Tier-1 not started.**
- **Artifacts** — the 6 `<case>/trajectory_matching/<variant>/fit_parameters.json`
  + `shared/trajectory_matching/per_case_form_summary.json`.

## Tier-0 comparison gate → same-smoothed `|v2|` — delivery record (2026-06-15)

Completed a mid-implementation that had been left half-wired: switching the
`tier0_drag_comparison.py` **scored gate** from the raw-HeDFT distance + mean
`|v|` RMSEs to the **same-smoothed `|v2|` RMSE against `cleaned_data_long.csv`**
(column `cleaned_SG`). Rationale: under Method B the coefficients were fit by
minimizing exactly that same-smoothed trajectory RMSE
(`METHOD_B_trajectory_matching_extraction.md` §6), so the Tier-0 report must
score *what the parameters were fit to*, not an unrelated raw-HeDFT metric.

### State recovered

The prior session had already rewritten `score()` (new `smoothed:
SmoothedSpeedReference` parameter; computes `compare_speed_to_reference(ion,
atom="I2", …)` against the smoothed reference; returns `v_I2_smoothed_rmse_Aps`
/ `v_I2_smoothed_mean_ratio` / `n_scored_smoothed`) and added the imports
(`SmoothedSpeedReference`, `load_smoothed_speed_reference`,
`compare_speed_to_reference`), but stopped before wiring the call site —
`main()` still called `score(ion, hedft, window)` (a `TypeError`: the smoothed
arg was missing and never loaded) and `print_summary()` still labelled the raw
distance / mean `|v|` as the GATE.

### Delivered (code)

- **`scripts/post_processing/tier0_drag_comparison.py`** — two edits completing
  the switch (no API/signature change beyond the already-present `score()`):
  - `main()` now loads `smoothed = load_smoothed_speed_reference(
    CLEANED_VELOCITIES_PATH_2)` (the per-case `cleaned_data_long.csv`) and passes
    it to `score(ion, hedft, window, smoothed)`.
  - `print_summary()` promotes `v_I2_smoothed_rmse_Aps` to the GATE line
    (labelled same-smoothed Method-B objective, with sample count) and demotes
    raw distance + mean `|v|` into the diagnostics block.
- The smoothed loader accepts both layouts — 18 Å 2-column (`time,cleaned_SG`)
  and 9 Å 3-column (`+IMF_cleaned`) — so a single load path covers both cases.

### Verification

- `py_compile` OK; `pytest tests/test_tier0_drag_comparison.py -q` → **4 passed**
  (the regression test uses its own `_score` on the package `compare_*` APIs, so
  it is independent of the script's `score()` signature — the change is safe).
- `load_smoothed_speed_reference` confirmed on both files: 9 Å 11412 pts
  `[2.670, 14.081]`, 18 Å 10229 pts `[4.540, 14.768]`.

The gate is now I2-velocity-only by construction (no smoothed distance or
smoothed I1 reference exists); raw distance + per-atom `|v|` survive as reported
diagnostics. No physics, presets, bundles, or thresholds changed.

### Doc — `TIER0_SCRIPTS.md` rewrite (2026-06-15)

Rewrote `scripts/TIER0_SCRIPTS.md` to the **three remaining** Tier-0 scripts at
current state, dropping the stale sections for the two scripts removed in
`d93e3a3` (`tier0_tstar_seeded_comparison.py`,
`tier0_same_smoothed_comparison.py`) and the dead `TIER0_COMPARISON_spec.md`
reference. The doc now describes: `tier0_common.py` (the shared `CATALOG` +
`build_drag_cfg` + naming + per-case `window_source_dir`), `gen_tier0_runs.py`
(catalog-driven run generation with optional inline hand-tuning), and
`tier0_drag_comparison.py` (scores the **same-smoothed `|v2|`** gate, raw
metrics demoted to diagnostics). Intro updated for the retired
consistency-check → held-out-generalization role and the `shared_pure_cubic`
production law. No code referenced by the doc changed beyond the gate switch
above.

## Tier-1a plan refinement — decision record (2026-06-23)

Plan-doc refinement of `TIER1A_IMPLEMENTATION_PLAN.md` after a code + cross-doc
audit. No code written (the `[PROCEED TO IMPLEMENTATION]` boundary holds); these are
documentation + scoping decisions for the eventual build.

- **SQ1 re-scoped as reused, not new.** The plan framed SQ1–SQ3 as the integrator
  "work." Audit confirmed against MASS A13 ("SQ1 built, accepted as-is") and
  `physics/baoab.py`: the O-step already applies `e^(−γ·dt/m)` with γ frozen at
  `v_in`, books exact dissipation, and `make_ion_baoab_step` is rebuilt every step so
  `m` can vary. Genuine Tier-1a integrator work is **SQ2 (mass-jump operator) + SQ3
  (post-jump `m⁺`) + the `m(t)` plumbing** into the existing per-step rebuild. The
  `fixed`-mode bit-for-bit regression is repositioned as the SQ1-untouched guard.

- **Mass-scenario enum: `scenario_A_accretion` and `scenario_B_stripping` RETIRED;
  `anchored_discrete` ADDED.** User decision (2026-06-23). New literal:
  `MassScenario = {fixed, biphasic, anchored_discrete}`. A and B are superseded by the
  locked `biphasic` mechanism (DESIGN §2.5/§2.8) — they were inert enum members read
  only by the `check_drag_config` guard. **Build touch-points:** the `MassScenario`
  literal in `config.py`; the guard's non-`fixed` branch set (`config.py` ~417–448) →
  `{biphasic, anchored_discrete}`; any preset/test referencing the A/B names. DESIGN
  §2.8 (which still says "A, B retained as baselines") to be updated at build time.
  `anchored_discrete` is a non-`fixed` scenario → trips the §6.5 `time_resolved`
  pairing guard structurally against the constant-`m_eff` Tier-0 coefficients → runs
  under `allow_inconsistent_mass_pairing=True` on the §6.6 mid-window defence (R6).

- **Checkpoint v5 → v6 delta clarified.** The four-term ledger arrays already exist
  (`E_kin_eV`/`E_pot_eV`/`E_dissip_eV`, `(2N,T)`). v6 = rename
  `E_mass_attach_defect_eV → E_mass_transfer_eV` (DESIGN §2.9), add per-atom
  `n_shell (2N,T)`, drop the `mass_history_kg` non-decreasing assumption, add a
  `mass_scenario` metadata field. No `E_int` at 1a.

- **Three open items resolved (user, 2026-06-23):** (1) a numeric `t_star` is **not
  required for Tier 1a** — the class verdict and segment-2 timings are `t*`-independent;
  leave `t_star_ps` a config parameter with a placeholder. (2) **No new run-size
  decision** — reuse the existing standard run preset, invoked with
  `mass_scenario=anchored_discrete` (and `fixed` for the null). (3)
  `coulomb_available_eV=0.80` is stamped **for provenance only — no hard refuse**.

## Tier-1a plan refinement — decision record (2026-06-24)

Second planning pass (study + discussion of CLAUDE.md / DESIGN / TIER1A plan). Still
code-free; the `[PROCEED TO IMPLEMENTATION]` boundary holds. Three decisions, two of
which **supersede** parts of the 2026-06-23 record above.

- **Deliverable reframed — the telescoping boost is NOT the verdict (supersedes
  2026-06-23 item).** The 2026-06-23 record left `t_star_ps` a placeholder on the
  ground that "the class verdict (telescoping boost 1.153) … [is] `t*`-independent."
  That conflates the *force-free* kinematic ceiling with the *full-force* trajectory:
  1.153 telescopes timing-independently **only** when drag and Coulomb are off. Under
  full drag+Coulomb each shed's firing *time* feeds back through the $\propto v^3$ drag
  work, so $R(t)$, $|v(t)|$ and the endpoint are genuinely `t*`-sensitive. The 1.153
  boost is demoted to a **force-free integration-test sanity ceiling**; the substantive
  Tier-1a result is the trajectory's response to shed timing. (User correction.)

- **`t*` sweep pinned to `{0.5, 5, 9}` ps (wide span; supersedes the 2026-06-23
  placeholder).** Run matrix: one `fixed` null (t*-independent, no sheds) + three
  `anchored_discrete` runs, one per swept `t*`. Values chosen for maximal contrast in
  onset timing; physically-admissible window confirmed `t* ∈ (0,10)` ps (user).

- **Checkpoint v6 = back-compat load shim (new detail).** v6 writer as before
  (rename `E_mass_attach_defect_eV → E_mass_transfer_eV`, add `n_shell (2N,T)`, drop
  the `mass_history_kg` non-decreasing assumption, add `mass_scenario` metadata). The
  v6 **loader accepts legacy v5** — maps the renamed field and synthesizes an absent
  `n_shell` (constant at the run's fixed shell count) — so the 14 existing v5 `ion.npz`
  run dirs (incl. the Tier-0 `shared_pure_cubic` runs) still load instead of failing
  the version check. (User decision.)

- **Docs updated this pass:** `TIER1A_IMPLEMENTATION_PLAN.md` §1/§2/§8/§9/§11
  (deliverable reframing, `t*` sweep + run matrix, back-compat shim note);
  `DRAG_PORT_DESIGN_DECISIONS.md` §2.8 (A/B recorded as retired, not "retained").

