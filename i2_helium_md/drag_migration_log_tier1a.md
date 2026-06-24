# Drag-Model Port — Migration Log: Tier 1a

**Purpose.** This log is the **decision + delivery history for Tier 1a** (anchored
kinematic mass-dynamics validation). Tier-0 and extraction history (Slices 1–4, the
Tier-0 diagnosis, Method-B / shared-form / §10 form-discrimination) lives in the
companion `drag_migration_log_tier0.md`. The other companion docs (`CLAUDE.md`,
`TIER1A_IMPLEMENTATION_PLAN.md`, `DRAG_PORT_DESIGN_DECISIONS.md`,
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`) carry the present state plus the
live plan/specs and point here for history.

**Companion docs:** `drag_migration_log_tier0.md` (Tier-0 + extraction history),
`TIER1A_IMPLEMENTATION_PLAN.md` (the active build plan),
`DRAG_PORT_DESIGN_DECISIONS.md` (the frozen design),
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (the live mass-model detail).

---

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

---

## Tier-1a Slice S — delivery record (2026-06-24): **first Tier-1a code, the schedule generator**

The first Tier-1a build unit, behind the `[PROCEED TO IMPLEMENTATION]` trigger.
Slice S only — Slices M / I⋆ / B, the `SimConfig` anchor/`t_star`/`anchored_discrete`
fields, and the v5→v6 schema bump remain deferred. Plan: `TIER1A_IMPLEMENTATION_PLAN.md`
§4 (Slice S), §10 (oracle).

### Delivered (code)

- **`physics/shell_schedule.py`** (new, pure/stateless): `build_shell_schedule(
  t_star_ps, crossing_fraction=0.5)` → frozen `ShellSchedule`. Exposes **two distinct
  shell quantities** (the clarification that drove the post-build fix):
  - `n_bar(t)` — the **continuous** anchored loss curve (piecewise-linear interpolant
    of the TDDFT waypoints `(t*,21),(10,19),(14,14)`). Fractional **by construction**;
    its only role is to locate the half-integer downward crossings `n_bar = n−½` (the
    S4 rule). It is *not* the physical count.
  - `n_of_t(t)` — the **physical integer** shell count: a piecewise-constant staircase,
    21 until the first shed, −1 at each of the 7 shed times, flat 14 past 14 ps. This —
    never `n_bar` — is what the complex mass `m(t)=MASS_I_ION_AMU+n(t)·MASS_HE_AMU`
    consumes downstream (mass jump, integrator).
  - `events` — 7 ordered `ShedEvent`s (`time_ps`, `n_before→n_after`, pre/post mass,
    `kick_factor = m/(m−m_He)`), solved analytically (no root-finding). Fail-loud guards
    on `t_star_ps∈[0,10)`, `crossing_fraction∈(0,1)`, count==7, strict-increasing times
    in `(t*,14]`. Mass-agnostic to the integrator (takes no `m` from outside).
- **`physics/constants.py`** (additive): `MASS_HE_AMU = 4.0026`, `MASS_I_ION_AMU =
  126.90` — the Tier-1a shell-schedule iodine-ion reference, **intentionally distinct**
  from the rounded MD `MASS_I_AMU = 127.0` (commented; do not unify). Both exported via
  `physics/__init__.py`.
- **`scripts/post_processing/plot_shell_schedule.py`** (new, standalone — no run dir /
  `SimConfig`): foregrounds the integer `n(t)` staircase for the `t*∈{0.5,5,9}` sweep,
  `n_bar(t)` demoted to a thin dashed guide, sheds marked.
- **`tests/test_shell_schedule.py`** (new, 50 tests): structure, segment-1/2 timing
  (analytic, tight), masses/kicks vs §10, telescoping invariant, the integer-count
  staircase, `n_bar` evaluator, fail-loud.

### Oracle note (recorded — not a bug)

The §10 **kick factors** (1.0193…1.0219) and **telescoping product** (1.1532) are
ratios `m/(m−m_He)` and reproduce the table to its 4 printed decimals **exactly**
(asserted tight) — they are the real oracle. The §10 **absolute pre-shed masses** are
internally rounded: the `n=19` row is pinned to config `m_eff=202.953908` (precise
iodine 126.9045) while the formula labels I⁺ as `126.90`, so the interior rows (n=16–19)
sit ~0.0045 amu off the pure `126.90+n·4.0026` formula this module uses. Asserted
against the formula tight and against the §10 table loose (~5e-3 amu); documented in the
module/test docstrings.

### Verification

- `tests/test_shell_schedule.py`: 50 passed. Full suite: **787 passed, 0 failed** (no
  regression; prior baseline 643/0, suite has since grown). Plot script runs headless
  (Agg); generated PNG not committed (CLAUDE.md figure rule).

### Post-build correction (same session)

User flagged that a *fractional* shell count is unphysical. Resolved: `n_bar` is correct
as the continuous loss-curve scaffolding, but the **integer** `n_of_t(t)` staircase was
promoted to a first-class output and the plot reworked to foreground it. The two
quantities are now explicitly separated in the module docstring (`|n − n_bar| ≤ ½` tie
asserted in tests).

---

## Tier-1a Slice M — delivery record (2026-06-24): **the cold-shed mass-jump operator (SQ2)**

The second Tier-1a build unit, behind the `[PROCEED TO IMPLEMENTATION]` trigger.
Selected as the next slice because it is the only remaining *independent* unit on the
critical path to I⋆ (S, M, B are independent per plan §6; B is the other independent
leaf but is entangled with the v5→v6 schema bump). Plan: `TIER1A_IMPLEMENTATION_PLAN.md`
§4 (Slice M), §2 (locked physics), §10 (oracle).

### Delivered (code)

- **`physics/mass_jump.py`** (new, pure/stateless — mirrors `shell_schedule.py`):
  - `cold_shed(v_minus, m_minus_amu, *, m_he_amu=MASS_HE_AMU) → ShedResult` — the SQ2
    momentum-conserving reset (He shed at rest, `m·v` invariant): `v⁺ =
    m/(m−m_He)·v⁻`, `m⁺ = m−m_He`, and `dE_mass_transfer = −½·(m·m_He)/(m−m_He)·‖v⁻‖²`.
    The energy term is the **exact reduced-mass form**, booked negative so the §2.9
    four-term ledger closes by construction — **not** the heavy-ion `½·m_He·v²`
    approximation (the two agree only as `m→∞`; ~3 % apart at `n=1`).
  - `kick_factor(m_minus_amu, m_he_amu)` — the `m/(m−m_He)` primitive (≡ Slice S
    `ShedEvent.kick_factor`).
  - `apply_shed(..., mode=...)` — the Tier-1a A/B selector: `fixed` null (no reset,
    mass held, zero defect) vs `anchored_discrete` (delegates to `cold_shed`). **Mode
    is a function argument, deliberately NOT `SimConfig.mass_scenario`** (config-surface
    work deferred to I⋆).
  - `ShedResult` frozen dataclass (`v_plus`, `m_plus_amu`, `dE_mass_transfer`).
    Mechanical-amu units (`amu·Å²/ps²`), mass-agnostic to the drag law. Fail-loud
    guards (`m_He>0`, `m⁻>m_He` so `m⁺>0`, finite inputs).
- **`tests/test_mass_jump.py`** (new, 32 tests): momentum conservation to machine
  precision; post-jump mass; direction preservation; kick vs §10 oracle (4 figures)
  **and** cross-checked against live `shell_schedule.ShedEvent.kick_factor` (ties M to
  the delivered Slice S without touching its internals); telescoping product
  `m(21)/m(14)=1.1532`; exact-vs-heavy-ion ~3 % flag at `n=1`; heavy-ion limit recovery
  at large `m`; `fixed`/`anchored_discrete` modes; fail-loud guards.

No `physics/__init__.py` export (intentional, matching `shell_schedule.py`'s isolation;
I⋆ imports directly). No constants added (Slice S already landed `MASS_HE_AMU`,
`MASS_I_ION_AMU`).

### Deferred (explicitly out of Slice M scope — unchanged by this delivery)

- **SimConfig `MassScenario` enum surgery** — retiring `scenario_A_accretion` /
  `scenario_B_stripping`, adding `anchored_discrete`, and editing the `config.py`
  ~418–448 consistency guard — stays with **I⋆** (user decision, 2026-06-24). Slice M
  models the two A/B modes by a *function argument*, so the config surface is untouched
  and M stays bit-isolated.
- **`IonCheckpoint` v5→v6 schema bump** (rename `E_mass_attach_defect_eV →
  E_mass_transfer_eV`, add `n_shell`, drop the `mass_history_kg` non-decreasing
  assumption, add `mass_scenario` metadata, back-compat v5 load shim) — **Slice B /
  wiring**, not touched here.
- **Integrator wiring, drag evaluation, real trajectory** — **I⋆** (now unblocked: it
  composes S + M into the existing per-step `make_ion_baoab_step` rebuild).
- **`E_int` internal-energy reservoir** — out of scope for all of Tier 1a; the shed
  boost is energetically unsourced in the ledger by agreement (plan §1).

### Verification

- `tests/test_mass_jump.py`: 32 passed. Full suite: **819 passed, 0 failed** (no
  regression; prior baseline 787/0 after Slice S). No figures/checkpoints/config
  generated (CLAUDE.md test rules).

### Follow-up (next slice)

- **I⋆** is the natural next build (composes S + M; carries the SimConfig enum surgery
  + the integrator-rebuild `m(t)` plumbing + SQ3 post-jump `m⁺`). **B** remains
  independently buildable but is better sequenced once M emits real
  `dE_mass_transfer` into a ledger, and it owns the v5→v6 schema bump.

---

## Tier-1a Slice I⋆ — delivery record (2026-06-24): **variable-mass integrator wiring (SQ2–SQ3 + `m(t)` plumbing)**

The third Tier-1a build unit, behind the `[PROCEED TO IMPLEMENTATION]` trigger.
Composes the finished S + M into the existing per-step BAOAB rebuild and carries the
config-surface enum surgery deferred from M. Plan: `TIER1A_IMPLEMENTATION_PLAN.md`
§4 (Slice I⋆), §2 (SQ upgrade), §8 (config contract).

### Design point resolved — jump at the **step seam** (not inside B-A-O-A-B)

The plan's §2 "B/A → jump → O" ideal is realized as **jump-then-(BAOAB)** at step
granularity: the cold shed is applied to `(v, m)` **before** the per-step
`make_ion_baoab_step` rebuild, so the rebuilt closure reads `m⁺` for *both* the
conservative kicks and the drag O-step (SQ3). `baoab.py` is **untouched** (SQ1
reused verbatim). The order reduction (first half-kick uses `m⁺` rather than `m⁻`)
is `O(dt)` in mass and benign as `dt→0` (jumps are `dt`-independent in count; plan
§2). This is the recommended locus from the plan's open-question.

### Delivered (code)

- **`physics/mass_jump.py`** (additive): `cold_shed_velocity_components(vx, vy, vz,
  m_minus_amu, *, m_he_amu)` — the per-atom `(2N,)` vectorization of `cold_shed`
  (uniform scalar mass/kick across atoms, per-atom `|v_i|²` defect), and a shared
  `_reduced_mass_defect_coeff` so the exact reduced-mass form lives in **one** place
  (scalar `cold_shed` refactored to use it; numerically identical, M's 32 tests
  unchanged). Mass-agnostic to the drag law.
- **`simulation/ion_propagation_step.py`** — `shed_step(state, schedule,
  next_shed_idx, dt)` (the SQ2/SQ3 pre-step): fires **≤1** scheduled shed per step
  (next pending event with `time ≤ t+dt`), applies the momentum-conserving reset,
  emits `m⁺` and books the per-atom reduced-mass defect (eV) into the **existing**
  `E_mass_attach_defect_eV` field (semantic generalization — **no schema change**;
  the rename is Slice B). `_check_drag_scope` relaxed to admit `anchored_discrete`
  (still refuses `biphasic`); the realized-mass `m_eff` trip-wire scoped to `fixed`
  only (anchored mass legitimately runs n=21→14, defended by §6.6, not the band).
- **`simulation/ion.py`** — driver builds the schedule once
  (`build_shell_schedule(cfg.t_star_ps)`, `anchored_discrete` only) and calls
  `shed_step` before each rebuild; the `fixed` path is byte-for-byte the Tier-0 path
  (schedule never built → regression guard).
- **`simulation/ion_initial_state.py`** — `anchored_discrete` initial mass override
  to the n=21 complex mass (210.955 amu), distinct from the `fixed` m_eff override
  and the inherited neutral mass.
- **`config.py`** — `MassScenario = {fixed, biphasic, anchored_discrete}`
  (`scenario_A_accretion`/`scenario_B_stripping` **retired**); `_EVOLVING_MASS_SCENARIOS
  = (biphasic, anchored_discrete)`; new `AnchorMode` alias; new fields `t_star_ps`
  (default 5.0; sweep {0.5,5,9}), `anchor_mode="time"`, `coulomb_available_eV=0.80`
  (provenance stamp, **no hard refuse**). The guard routes `anchored_discrete` into
  the `time_resolved`-requiring arm → trips against constant m_eff coeffs → runs
  under `allow_inconsistent_mass_pairing=True` (R6, §6.6).

### Tests

- **`tests/test_ion_variable_mass.py`** (new, 11): vectorized cold shed vs the
  scalar M oracle; `shed_step` SQ2 reset / SQ3 `m⁺` / ≤1-per-step under a dense
  window / no-fire identity / exhausted-schedule no-op; jump-step **measure-zero**
  (7 sheds invariant to `dt∈{0.001,0.01,0.1}`); **telescoping** speed boost ×1.1532.
- **`tests/test_ion_drag_smoke.py`** (+3): `fixed` run independent of `t_star_ps`
  (SQ1-untouched regression); a full `anchored_discrete` run sheds 7 He n=21→14;
  the **four-term ledger closes** (`E_kin+E_pot+E_dissip+E_mass_transfer`, <5e-2) with
  the defect channel carrying energy — the M reset is the exact form, not a relabel.
- **`tests/test_ion_initial_state.py`** (+1): `anchored_discrete` starts at the
  n=21 mass. Retired A/B names updated in `test_drag_config.py`,
  `test_ion_drag_smoke.py`, `test_ion_initial_state.py`, `test_mass_jump.py`.

### Deferred (explicitly out of Slice I⋆ scope — to Slice B)

- The `E_mass_attach_defect_eV → E_mass_transfer_eV` rename, the per-atom `n_shell`
  array, the `IonCheckpoint` v5→v6 bump + back-compat v5 load shim, the standalone
  four-term ledger module, and the full `t*`-sweep RMSE table (the §7/§9 integration
  phase). I⋆ stays schema-neutral (still writes/loads v5).

### Verification

- `tests/test_ion_variable_mass.py`: 11 passed; targeted config/init/smoke: 117
  passed. Full suite: **834 passed, 0 failed** (no regression; prior baseline
  819/0 after Slice M; +15). No figures/checkpoints/config generated.

### Follow-up (next slice)

- **B** is now the natural next build: M+I⋆ emit real `dE_mass_transfer`; B owns the
  v5→v6 schema bump (field rename, `n_shell`, back-compat shim) and the four-term
  closure module, then the integration-phase `t*`-sweep RMSE table.
