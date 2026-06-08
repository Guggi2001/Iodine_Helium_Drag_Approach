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

