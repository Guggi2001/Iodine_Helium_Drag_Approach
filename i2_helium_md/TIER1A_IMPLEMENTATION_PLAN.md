# Tier 1a — Implementation Plan (Anchored Kinematic Validation)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no LAMMPS, no pseudo-code until
> the explicit `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked*
> formulations from the project docs (LaTeX + dimensional checks); module
> descriptions are *interface contracts*, not implementations.
> 
> **Entry doc:** `DRAG_PORT_DESIGN_DECISIONS.md` §6.4 (Tier 1 → 1a). Mass model and
> SQ preconditions: `MASS_DYNAMICS_LOCKED_…md` (§4, §6, A13, 2026-06-21 cont. 3).
> Parameter classes: `CALIBRATION_MAP.md`.

---

## 0. Status and intent

Tier 1a asks **one** question: *does mass dynamics change the ion's translational
trajectory relative to a constant-mass ion, and if so, how does that change appear
in the anchored 9 Å comparison?*

The He shell schedule $n(t)$ is **anchored** to the 9 Å TDDFT loss curve (read in,
not generated), so the kinetics — electronic picture, ladder $D_0$, $\kappa$, $\nu$,
$s$ — **never enter**. With the schedule fixed, the only thing that varies between
the two runs is whether the ion's mass is held constant or follows $m(t)$. Tier 1a
is therefore a **controlled A/B comparison**:

- **`fixed`** — constant $m_\text{eff}$ (the null).
- **`anchored_discrete`** — anchored variable $m(t)$, losing one He per shed event with a
  **continuous velocity** (`v^+=v^-`). **Status after Slice C:** this is now the
  physical Tier-1a anchored-mass comparison. The older cold-shed reset is retained
  only as an explicitly labelled diagnostic upper-bound / stress-test because it
  needs the full `E_int` energy-gated tier to be physically meaningful.

The comparison isolates the *influence of mass dynamics* on $R(t)$, $|v(t)|$, and on
whether the kinematic energy/momentum bookkeeping closes. It does **not** test the
generative mechanism that *produces* the schedule — anchored away here, arbitrated
by experiment at Tier 2.

**Central engineering fact.** The BAOAB ion integrator's drag O-step is **already
built and accepted** (MASS A13: SQ1 "built, accepted as-is"). `physics/baoab.py`
applies `e^{-γ·dt/m}` with `γ` frozen at the step-entry velocity, books exact
dissipation, and — critically — `make_ion_baoab_step` is **rebuilt every step** by
the driver expressly so the mass can change (docstring, `baoab.py:94–97`). So the
substance of this build is **not** a from-scratch O-step. The genuine new work is
**SQ2** (the scheduled mass-shed operator; continuous-velocity in the current
Tier-1a driver), **SQ3** (the post-jump O-step reading `m⁺`), and the **`m(t)` plumbing** that feeds the schedule-driven mass into
the existing per-step rebuild. SQ1's mechanism is **reused**; the constant-mass
bit-for-bit regression (§9) is its untouched-guard, not a re-derivation.

---

## 1. Scope lock — what 1a does and does not test

**Tests:**

- Whether an anchored variable mass $m(t)$ changes the translational trajectory of
  the 9 Å TDDFT ion relative to constant mass, in the scalar observables **$R(t)$
  (radial distance) and $|v(t)|$ (speed)** — the native coordinates, since the 9 Å
  signal is not clean-radial (transverse wobble, the known 9 Å model-dimensionality
  flag). Matches the Tier-0 distance-/$|v|$-RMSE convention.
- Whether the **four-term** kinematic ledger closes:
  $E_\text{kin}+E_\text{pot}+E_\text{dissip}+E_\text{mass\_transfer}=\text{const}$.
- The **trajectory response to shed timing**: a `t*` sweep ($t^*\in\{0.5,5,9\}$ ps)
  under full drag+Coulomb, comparing how $R(t)$ and $|v(t)|$ of `anchored_discrete`
  move relative to the `fixed` null. After Slice C the shed event itself does **not**
  inject a speed kick; any trajectory change comes from the post-shed mass feeding
  subsequent conservative and drag dynamics. The old force-free $|v|\propto1/m$
  telescoping boost (1.153) is now a cold-shed-bound diagnostic only.

**Does NOT test (out of scope by agreement — do not let these leak in):**

- The biphasic **generative mechanism** (Poisson pickup + RRK + self-bound gate).
  Anchoring the schedule bypasses it; unfalsified until Tier 2.
- The **internal-energy reservoir** $E_\text{int}$. Dropped for 1a. Consequence:
  true cold-shed evaporation is not modeled as the Tier-1a physical path. The
  current production path removes a co-moving He atom without a velocity boost;
  cold-shed waits for the full energy-gated tier or remains a labelled bound.
- Any **picture/ladder/$\kappa$/$\nu$/$s$** value. Bypassed.
- **Shell-timing prediction ("1b").** Rejected: TDDFT is not ground truth
  (experiment arbitrates at Tier 2), and a timing match would calibrate the
  early-window flow to the 0.80 eV validation regime, not the 2.70 eV production one.

---

## 2. Locked physics — constants, equations, SQ upgrade

All dimensional checks below balance; they are the test oracles for the new modules.

**Constants.**

| Symbol                      | Value                      | Units         | Source                                        |
| --------------------------- | -------------------------- | ------------- | --------------------------------------------- |
| $m_{\mathrm I^+}$           | 126.90                     | amu           | atomic                                        |
| $m_\text{He}$               | 4.0026                     | amu           | atomic                                        |
| $b$                         | 2.5154                     | amu·ps/Å²     | Tier-0 locked (pure cubic)                    |
| $g$                         | $\in[0,1]$                 | dimensionless | density-proportional gate (erf complement)    |
| $E_\text{avail}^\text{ion}$ | **0.80**                   | eV            | validation, $d{=}9$ Å, $\tfrac12\cdot14.40/9$ |
| anchors                     | $(t^*,21),(10,19),(14,14)$ | ps, count     | 9 Å TDDFT (author-confirmed)                  |
| $m_\text{eff}$              | 202.95                     | amu           | 19 He (drag-extraction reference)             |

**Drag force law (Tier-0, plug-in).** Pure cubic, applied with no $M$-scaling
(`MASS` §2.2):
$$F_\text{drag}=g\,b\,v^3,\quad \gamma(v)=g\,b\,v^2,\quad a_\text{drag}=-\frac{\gamma(v)\,v}{m(t)}$$
$[\gamma]=\text{amu/ps}$, $[F_\text{drag}]=\text{amu·Å/ps}^2$ (force) ✓,
$[a_\text{drag}]=\text{Å/ps}^2$ ✓. The force *value* is Tier-0; how it enters the
O-step under a changing $m(t)$ is SQ1 (below). The lighter post-shed complex
decelerates *harder* per unit drag force.

**Production Tier-1a shed (continuous velocity; co-moving He).**
$$v^+=v^- ,\qquad m^+=m-m_\text{He}$$
The removed He is booked as leaving co-moving with the complex at the instant of
the mass update. Total momentum and kinetic energy of "remaining complex + removed
He" are conserved, while the tracked complex alone loses the kinetic energy carried
away by the removed atom.

**Production mass-transfer ledger term.**
$$\Delta E_\text{mass\_transfer}=+\tfrac12\,m_\text{He}\,\lVert v^-\rVert^2$$
amu·Å²/ps² = energy ✓. The sign is positive because the tracked complex's
post-shed `E_kin` drops by exactly this amount when its mass changes at unchanged
velocity, and the ledger term compensates that drop. The downstream eV conversion is
the same mechanical-units path used for drag dissipation.

**Cold-shed diagnostic bound (not the production Tier-1a path).**
$$v^+=\frac{m}{m-m_\text{He}}\,v^- ,\qquad
\Delta E_\text{mass\_transfer}=-\tfrac12\,\frac{m\,m_\text{He}}{m-m_\text{He}}\,\lVert v^-\rVert^2$$
This reset corresponds to He leaving at rest and gives the force-free telescoping
speed ceiling
$$\prod_{k}\frac{m_k}{m_k-m_\text{He}}=\frac{210.955}{182.936}=\mathbf{1.153}.$$
Slice R showed that applying this bound as the Tier-1a production model creates
unphysical discontinuous `|v2|` jumps because `E_int` is absent. It is retained as a
diagnostic upper bound and deferred to the later full energy-gated evaporation tier.

### The SQ upgrade — fixed-mass BAOAB → variable-mass-correct

SQ1 is **built and reused** (MASS A13); SQ2 and SQ3 are the **new work**, plus the
`m(t)` plumbing. The integrator's O-step already divides by a per-step `m` and the
factory already rebuilds every step — so making `m` *follow the schedule* is a
driver-side change, not an O-step rewrite.

- **SQ1 — drag O-step under varying $m(t)$ (BUILT, reused).** The O-step already
  freezes $\gamma$ at the step-entry velocity and applies $e^{-\gamma\,dt/m}$ at the
  per-step mass (`baoab.py`), giving exact dissipation at any $dt$,
  $\Delta E_\text{dissip}=\tfrac12 m(\lVert v_\text{in}\rVert^2-\lVert v_\text{out}\rVert^2)$,
  and unconditional dissipativity ($\lVert v_\text{out}\rVert\le\lVert v_\text{in}\rVert$).
  $O(dt)$ on the cubic (second-order claim retired); residual one-signed over-braking
  bias (R10). **The only Tier-1a change is that the `m` passed to each per-step
  rebuild now follows the schedule** — no `_o_step` edit. The constant-mass regression
  (§9) proves this path is byte-identical under `fixed`.
- **SQ2 — scheduled mass shed, new.** The continuous-velocity mass update above —
  the **invariant-closure precondition** for the current Tier-1a production path.
  The cold-shed reset remains available only as a diagnostic bound. Ordering inside the step:
  **jump-then-O**, **at most one mass event per step** (sheds-only at 1a). Jump-step
  order reduction is benign (jumps are $dt$-independent in count → measure-zero as
  $dt\to0$).
- **SQ3 — post-jump O-step uses $m^+$, new.** After a jump, the drag O-step (SQ1)
  must use the **post-jump mass** $m^+$ in the friction term (noise off at 1a, so
  friction only). Forced by SQ2; a definiteness requirement.

**Discarded — continuous (Meshchersky) realization.** A smooth
$m(t)\dot v=F-v\dot m$ mass-loss term is **unphysical** here: He is shed as discrete
cold atoms, not a continuous jet, and the smooth form misrepresents the per-event
momentum bookkeeping. Removed; recorded so it is not re-litigated.

---

## 3. Component inventory — existing vs new

**Existing (plug-in, no build) — §5 for contracts:**

- **Drag force law** (Tier-0): the value $b$ and the cubic $\gamma(v)$, $a_\text{drag}$
  *evaluation* (`physics/drag.py`, mass-agnostic). Tested.
- **Drag O-step (SQ1)**: the frozen-$\gamma$ damping $e^{-\gamma\,dt/m}$ + exact
  dissipation in `physics/baoab.py`. Built and accepted (A13); **reused**, not rebuilt.
- **Validation scripts**: RMSE of $R(t)$, $|v(t)|$ vs the smoothed reference, per
  mode, relative ranking. Existing; plugged in, no work.

**New (build + test) — §4 slices:**

- **Slice S** — schedule generator (anchored $\bar n(t)$ → 7 shed events).
- **Slice M** — mass-jump operator and $m(t)$ state (the SQ2 reset + defect; emits
  $m^+$).
- **Slice I⋆** — **variable-mass integrator wiring (SQ2–SQ3 + `m(t)` plumbing)**:
  feed the schedule-driven $m(t)$ into the existing per-step `make_ion_baoab_step`
  rebuild; insert the mass jump jump-then-O at ≤1/step (SQ2, calling M); the post-jump
  rebuild reads $m^+$ (SQ3). The reused SQ1 O-step is unchanged.
- **Slice B** — four-term energy/momentum ledger and closure gate.
- **Slice R** — Tier-1a run orchestration, RMSE table, and trajectory diagnostics.
- **Slice C** — continuous-velocity Tier-1a shedding; supersedes cold-shed as the
  production `anchored_discrete` driver path and retags anchored runs to avoid stale
  cold-shed artifacts.
- **Slice V** — onset-violent stripping stress diagnostics; separate
  `onset_strip` schedule family for sensitivity/stress runs, not a replacement for
  the physical anchored-continuous Tier-1a RMSE table.

**Testing philosophy (every new slice):**

- Pure functions (S) tested against closed-form oracle values (§10).
- Stateful modules (M, B) tested on analytic micro-cases (single shed, synthetic
  stream) with neighbours **mocked**.
- I⋆ tested on analytic limits (constant-mass reduction, single-jump) before any
  real run.
- R tested on temporary run dirs / synthetic checkpoints; production `data/runs`
  generation stays out of pytest.
- No new slice's suite depends on another's implementation; shared fixtures in §10.
- Integration (§7) composes only **accepted** modules.

---

## 4. New-work slices

### Slice S — Schedule generator *(pure; fully independent)*

> **IMPLEMENTED (2026-06-24).** `physics/shell_schedule.py` —
> `build_shell_schedule(t_star_ps, crossing_fraction=0.5)` → frozen `ShellSchedule`
> with the continuous `n_bar(t)` loss-curve evaluator, the **integer** shell-count
> staircase `n_of_t(t)` (the physical count the mass consumes), and the 7 ordered
> `ShedEvent`s. Constants `MASS_HE_AMU`, `MASS_I_ION_AMU=126.90` added to
> `constants.py`. Oracle pytest `tests/test_shell_schedule.py` (50 tests) green; plot
> `scripts/post_processing/plot_shell_schedule.py`. Full suite 787/0. Delivery detail:
> `drag_migration_log_tier1a.md` (Slice S record, 2026-06-24).

**Purpose.** Turn the three anchors into the deterministic, monotone, sheds-only
integer schedule and its 7 fire events. No physics, no state.

**Interface.** Consumes anchors + $t^*$; emits (a) the piecewise-linear $\bar n(t)$
evaluator and (b) the ordered event list (fire-time, $n\to n-1$, pre-/post-shed
mass, kick factor).

**Encoded form.**


$$

\bar n(t)=
\begin{cases}
21, & t\le t^* \\
21-\dfrac{2(t-t^*)}{10-t^*}, & t^* < t\le10 \\
19-1.25\,(t-10), & 10 < t\le14 \\
14, & t>14
\end{cases}

$$

$
Segment loss rates  $2/(10-t^*)$ and $1.25  ps⁻¹ — a real ~5× acceleration kept as a
slope kink at 10 ps (flagged, not smoothed). **S4 discretization:** half-integer
downward crossing ($\bar n=n-\tfrac12$).

**Independence.** Depends only on constants + $t^*$. Built and tested first.

**Test spec.**

- Event count = 7; monotone; all fire-times in $(t^*,14]$.
- Segment-2 events at $\{10.4,11.2,12.0,12.8,13.6\}$ ps, uniform 0.8 ps, $t^*$-indep.
- Segment-1 events at $t^*+0.25(10-t^*)$, $t^*+0.75(10-t^*)$; sweep $t^*$, both $<10$.
- Kick factors match the §10 oracle to 4 figures.
- **Telescoping invariant:** product of kick factors $=m(t^*)/m(\text{end})=1.153$
  for any admissible $t^*$ and any S4 tie-break — verdict-robustness as a unit test.

**Acceptance.** All pass; tail flat at $n=14$ past 14 ps.

---

### Slice M — Mass-jump operator & state (SQ2) *(stateful; mockable schedule)*

> **IMPLEMENTED (2026-06-24).** `physics/mass_jump.py` — `cold_shed(v_minus,
> m_minus_amu, *, m_he_amu=MASS_HE_AMU) → ShedResult` (the SQ2 reset, $m^+$, and the
> exact reduced-mass `dE_mass_transfer`), `kick_factor(...)`, and `apply_shed(...,
> mode=...)` (`fixed` null vs `anchored_discrete`; mode a function arg, **not** the
> SimConfig enum). Pure/stateless, mechanical-amu, mass-agnostic to the drag law.
> Oracle pytest `tests/test_mass_jump.py` (32 tests) green; full suite 819/0. The
> SimConfig `mass_scenario` enum surgery, the v5→v6 schema bump, and the integrator
> wiring stay deferred (see below / §8). Delivery detail: `drag_migration_log_tier1a.md`
> (Slice M record, 2026-06-24). **Superseded for the production driver by Slice C:**
> cold-shed remains as a diagnostic bound; `anchored_discrete` now uses
> continuous-velocity shedding.

**Purpose.** Originally held $m(t)$ and performed the SQ2 cold-shed reset; after
Slice C the same module also owns the production continuous-velocity shed. It emits
the increment the ledger needs and the $m^+$ the integrator needs (SQ3).

**Interface.** Consumes a (mockable) fire event + current state; emits
$m^+=m-m_\text{He}$, $v^+=\tfrac{m}{m-m_\text{He}}v^-$, and
$\Delta E_\text{mass\_transfer}$. Two run modes:

- `fixed` — $m\equiv m_\text{eff}$, no reset, zero defect (null).
- `anchored_discrete` — reset + defect on each fire event.

**Independence.** Reset + defect testable on a **mock 1–2 event schedule**; no
integrator, no drag.

**Test spec.**

- Single-shed momentum conservation: $(m-m_\text{He})v^+=mv^-$ to machine precision.
- KE rises by exactly $\tfrac12\tfrac{m\,m_\text{He}}{m-m_\text{He}}\lVert v^-\rVert^2$;
  $\Delta E_\text{mass\_transfer}$ is its negative.
- Reduced-mass vs heavy-ion: assert the $(m-m_\text{He})$ form; flag the ~3 % $n{=}1$
  divergence from $\tfrac12 m_\text{He}v^2$.
- `fixed` applies no reset, holds $m_\text{eff}$, zero defect.

**Acceptance.** Momentum exact; defect exact; $m^+$ emitted for SQ3.

---

### Slice I⋆ — Variable-mass integrator wiring (SQ2–SQ3 + `m(t)` plumbing) *(composes S, M)*

> **IMPLEMENTED (2026-06-24).** Jump applied at the **step seam** (jump-then-BAOAB;
> `baoab.py` untouched, SQ1 reused). `physics/mass_jump.py` gains the per-atom
> `cold_shed_velocity_components`; `simulation/ion_propagation_step.py` gains
> `shed_step` (SQ2/SQ3, ≤1/step, books the defect into the existing
> `E_mass_attach_defect_eV` — schema-neutral); `simulation/ion.py` builds the
> schedule and calls it per step; `simulation/ion_initial_state.py` starts
> `anchored_discrete` at the n=21 mass; `config.py` retires A/B, adds
> `anchored_discrete` + `t_star_ps`/`anchor_mode`/`coulomb_available_eV`. Tests:
> `tests/test_ion_variable_mass.py` (11) + run-level additions in
> `tests/test_ion_drag_smoke.py`. Full suite 834/0. The v5→v6 schema bump, the
> field rename, and the `t*`-sweep RMSE table stay with **Slice B**. Delivery
> detail: `drag_migration_log_tier1a.md` (Slice I⋆ record, 2026-06-24).
> **Superseded for the production driver by Slice C:** the seam/timing/mass plumbing
> remains, but `shed_step` now applies the continuous-velocity primitive instead of
> `cold_shed_velocity_components`.

**Purpose.** Wire the schedule-driven $m(t)$ and the mass jump into the existing
per-step BAOAB rebuild. The drag O-step itself (SQ1) is reused unchanged; the work is
the jump branch, the post-jump mass, and feeding $m(t)$ through the seam that
`baoab.py:94–97` already exposes.

**Interface.** Per step: conservative B/A kicks (existing `acc_fn`: Coulomb +
droplet) → **query schedule** → if a fire event, apply the **M** jump
(**jump-then-O**, ≤1/step, **SQ2**) → drag O-step with frozen $\gamma$ at $v_\text{in}$
on the **current/post-jump mass** (**SQ1 + SQ3**). Emits next state, per-step
$\Delta E_\text{dissip}$ (SQ1), and $\Delta E_\text{mass\_transfer}$ (from M).

**What changes from the fixed-mass integrator.**

- The driver passes the schedule-derived $m(t)$ into each per-step
  `make_ion_baoab_step` rebuild (the seam already exists; SQ1 O-step untouched).
- A jump branch is inserted before the O-step, calling M (SQ2 ordering).
- The post-jump rebuild reads $m^+$ (SQ3 — automatic once the jump updates the mass
  state before the rebuild).

**Independence.** Composes S (events) and M (reset), both mockable; the existing
conservative kicks and the drag-force evaluation are reused unchanged.

**Test spec.**

- **Constant-mass reduction:** with `fixed` mode (no jumps), the upgraded step
  reproduces the existing fixed-mass integrator bit-for-bit (regression guard).
- **SQ1 (reused, exercised under running $m$):** drag-only decay (mock M no-op,
  Coulomb 0) follows the exact frozen-$v$ exponential at the running $m$;
  unconditional stability at large $dt$. Confirms the `m(t)` plumbing, not new O-step
  code.
- **Exact dissipation identity** at $dt\in\{0.001,0.01,0.05\}$ ps.
- **SQ2 ordering:** a fire event applies the reset **before** the O-step; ≤1 event
  per step even under a forced-dense schedule.
- **SQ3:** the post-jump O-step damping uses $m^+$, not $m^-$.
- **Jump-step measure-zero:** refining $dt$ leaves the event *count* fixed.

**Acceptance.** Constant-mass regression exact; SQ1–SQ3 all demonstrated on the
analytic limits above.

---

### Slice B — Four-term bookkeeping & closure *(stateful; synthetic-stream testable)*

> **IMPLEMENTED (2026-06-24, core only).** Checkpoint **v5→v6**
> (`simulation/checkpoint.py`): `E_mass_attach_defect_eV → E_mass_transfer_eV`,
> new `n_shell (2N,T)` + `mass_scenario` metadata, the non-decreasing-mass
> assumption dropped, and a **back-compat v5 load shim** (synthesizes `n_shell`
> from `mass_history_kg`, defaults `mass_scenario=fixed`) so existing v5 run dirs
> still load. Driver/state write the new fields (`ion_initial_state.py`,
> `ion_propagation_step.py`, `ion.py`). The closure gate is
> `postprocess/energy_balance.py::ion_ledger_closure` (+ `LedgerClosure`), reusing
> `ion_energy_totals`. Tests: `tests/test_checkpoint.py` (`TestIonSchemaV6` + v5
> shim) and `tests/test_energy_balance.py` (`TestLedgerClosure` incl. relabel-fault
> injection). Full suite 841/0. Delivery detail:
> `drag_migration_log_tier1a.md` (Slice B record, 2026-06-24).

**Purpose.** Accumulate the four-term ledger and assert closure.

**Interface.** Consumes per-step $\Delta E_\text{dissip}$ (SQ1) and per-event
$\Delta E_\text{mass\_transfer}$ (SQ2), plus $E_\text{kin}$/$E_\text{pot}$ from state;
emits the running invariant + closure residual.

**Encoded ledger.**
$$E_\text{kin}+E_\text{pot}+E_\text{dissip}+E_\text{mass\_transfer}=\text{const}$$
**$E_\text{int}$ explicitly absent.**

**What it certifies (and what it does not).**

- *Certifies:* the SQ2 reset is the reduced-mass form, not a relabel (label-only $v$
  → $E_\text{kin}$ jumps with no matching $E_\text{mass\_transfer}$ → residual blows
  up); and SQ1 drag work integrates to $\Delta E_\text{dissip}$ exactly.
- *Does not certify:* the energetics of shedding (unsourced at 1a). Closure is
  "by construction" given the reset — a **wiring correctness gate**, not a physics
  result. State this in the run report so closure is not over-read.

**Independence.** Tested on a **synthetic event stream** with known increments — no
real trajectory.

**Test spec.** Residual within Verlet-drift tolerance on the synthetic stream; a
deliberate relabel-instead-of-reset fault is **caught** (residual diverges).

**Acceptance.** Residual ≤ drift bound for correct streams; fault injection detected.

---

### Slice R — RMSE table + Tier-1a trajectory diagnostics *(reporting; consumes accepted S/M/I⋆/B)*

> **IMPLEMENTED (2026-06-24).** The deferred §5/§9 reporting layer is delivered:
> `scripts/tier1a_common.py`, `scripts/gen_tier1a_runs.py`, and
> `scripts/post_processing/tier1a_rmse_table.py`, with coverage in
> `tests/test_tier1a_scripts.py`. Full suite after the final plotting change:
> **846 passed, 1 expected warning** (the intentional `anchored_discrete`
> constant-coefficient pairing warning). Delivery detail:
> `drag_migration_log_tier1a.md` (Slice R record, 2026-06-24).

**Purpose.** Generate and inspect the four-run Tier-1a matrix: one `fixed` null plus
`anchored_discrete` at `t*={0.5,5.0,9.0}` ps, all on the 9 Å
`shared_pure_cubic` production law. The code reports diagnostics only; it does not
rank, threshold, or adjudicate the physics.

**Implemented surface.**

- `scripts/tier1a_common.py` reuses `scripts.tier0_common.build_drag_cfg` and the
  Tier-0 run-name convention, adding only the anchored-discrete config mutation and
  `tier1a_run_tag` / `tier1a_run_dir_name`.
- `scripts/gen_tier1a_runs.py` writes the four self-describing run dirs under
  `data/runs/`, defaulting to `CASE="9A"`, `VARIANT="shared_pure_cubic"`, `N=50`,
  `ION_TIME_PS=30.0`, `DT_ION_PS=0.01`.
- `scripts/post_processing/tier1a_rmse_table.py` emits a rich table with traceability
  columns (`case`, `variant`, `N`, `run_tag`, `scenario`, `t_star_ps`, scored-window
  endpoints), raw `R` RMSE, same-smoothed I2 `|v2|` RMSE, `v2` mean ratio,
  `ledger_max_resid_eV`, and `n_shell` start/end/shed count. It can optionally save
  CSV, export Tier-0-format mean-series CSVs, and build diagnostic figures.

**Visualization behavior.** The trajectory figure is intentionally narrower than the
original all-run overlay: it plots **only `|v2|`**, comparing the `fixed` run against
one user-selected anchored case (`PLOT_T_STAR_PS`, default `5.0`). This avoids hiding
pre-window anchored data under overlapping curves and directly shows the effect of
the selected shed schedule. Optional positions, energy, and radial-force figures
remain available for deeper inspection.

**Result and conclusion.** The generated Tier-1a checkpoints are full-length
trajectories (`0.0→29.99` ps for the default 30 ps ion run), and all anchored `|v2|`
traces are finite before the 2.67 ps scoring window. The clarified `|v2|` plots show
large discontinuous jumps at the scheduled shed events. Those jumps are the expected
mathematical consequence of the cold-shed reset `v⁺=m/(m−m_He)v⁻`, but they are **not
physically acceptable as the main Tier-1a anchored comparison** because Tier 1a
omits `E_int`. The cold-shed operator is therefore reclassified for Tier 1a as a
diagnostic upper-bound / stress-test; the physical anchored-mass Tier-1a comparison
should use a continuous-velocity mass update and defer true cold-shed evaporation to
the full energy-gated tier. **Slice C implements that correction; Slice R's cold-shed
run artifacts are now stale and must be regenerated under the `continuous` run tags.**

**Tests.** `tests/test_tier1a_scripts.py` covers anchored cfg validation, shared run
naming, rich scorer rows and CSV output, selected-`t*` `|v2|` plotting, and
mean-series export using temporary/synthetic data.

---

### Slice C — Continuous-velocity Tier-1a shedding *(physics correction; completed)*

> **IMPLEMENTED (2026-06-24).** The production `anchored_discrete` shed primitive
> now keeps velocity continuous at each scheduled mass event. The anchored schedule
> remains unchanged (`n=21→14` at the same TDDFT-anchored event times), BAOAB is
> still rebuilt after the mass update, and cold-shed helpers remain available only
> as an explicit diagnostic bound.

**Purpose.** Supersede the unphysical Tier-1a cold-shed production path exposed by
Slice R, without changing the shell schedule, Tier-0 drag law, checkpoints, or
integrator structure.

**Implemented behavior.**

- `physics/mass_jump.py` adds `continuous_velocity_shed(...)` and
  `continuous_velocity_shed_components(...)`.
- `apply_shed(..., mode="anchored_discrete")` delegates to the continuous-velocity
  primitive.
- `simulation/ion_propagation_step.py::shed_step` calls the continuous vectorized
  primitive, copies `vx/vy/vz` unchanged across the shed event, drops mass by one He,
  and books `+0.5*m_He*|v|^2` into `E_mass_transfer_eV`.
- `cold_shed(...)`, `cold_shed_velocity_components(...)`, and `kick_factor(...)`
  remain tested as the cold-shed diagnostic bound; they are no longer the Tier-1a
  driver path.
- `scripts/tier1a_common.py` retags anchored runs as
  `tier1a_anchored_continuous_t{...}`. Old cold-shed run dirs are left untouched but
  are no longer selected by the generator/scorer naming convention.

**Tests.**

- `tests/test_mass_jump.py`: continuous-velocity primitive invariants (velocity
  unchanged, one-He mass drop, total momentum/KE of remaining complex + removed
  co-moving He conserved, tracked-complex KE drop exactly booked positive).
- `tests/test_ion_variable_mass.py`: driver-facing `shed_step` leaves velocities
  unchanged, preserves the <=1 shed-per-step rule, keeps event count invariant under
  `dt`, and no longer has a telescoping speed boost.
- `tests/test_ion_drag_smoke.py`: anchored run still sheds 7 He (`n=21→14`), ledger
  closes with positive mass-transfer bookkeeping, and the old cold-shed velocity
  kick is absent from shed transitions.
- `tests/test_tier1a_scripts.py`: scorer/generator naming uses the new continuous
  anchored run tags, and selected-`t*` `|v2|` plotting remains fixed plus one
  anchored case.

**Acceptance.** Targeted Slice C suite green: 76 passed, 1 expected warning. Full
suite after this slice: 855 passed, 1 expected warning.

---

### Slice V — Onset-violent stripping stress diagnostics *(diagnostic; completed)*

> **IMPLEMENTED (2026-06-29).** A separate `onset_strip` diagnostic schedule family
> now stress-tests abrupt early stripping at `t_strip=0.5 ps` with endpoint
> `n_final in {14,2,1,0}`. This does **not** supersede the physical Tier-1a
> anchored-continuous comparison or its RMSE table; it is a stress/sensitivity
> family for probing onset-violent stripping behavior.

**Purpose.** Exercise extreme stripping endpoints under the existing drag/Coulomb
driver and ledger checks, while keeping the physical Tier-1a TDDFT-anchored schedule
and continuous-velocity interpretation separate.

**Implemented behavior.**

- Config accepts `anchor_mode="onset_strip"` plus `anchor_n_final`; the driver routes
  onset-strip schedules independently from the TDDFT-anchored `time` schedule.
- The onset event keeps velocity continuous, so there is no cold-shed boost at
  `t_strip`; batch mass transfer is positive and accounts for the co-moving removed
  He atoms.
- `scripts/gen_tier1a_stress_runs.py` generates stress run configs/tags such as
  `tier1a_stress_onset_strip_n0_t0.5`.
- `scripts/post_processing/tier1a_stress_table.py` scores stress rows only, validates
  cfg/ion metadata fail-closed, keeps the fixed null only for comparison
  plot/export, and provides the selected-`n_final` `|v2|` diagnostic plot.

**Tests.**

- Run-level smoke verifies the 21→0 two-level staircase, positive mass transfer, and
  finite velocities.
- Stress script coverage verifies cfg validation, metadata validation, stress-only
  table rows, fixed-null comparison handling, and selected-`n_final` plotting.

**Acceptance.** Focused verification after the final fix:
`pytest tests/test_tier1a_scripts.py tests/test_ion_drag_smoke.py -q` passed
(`64 passed, 17 warnings`). `py_compile` for the stress generator and scorer passed.
The full suite remains pending Task 9, so no full-suite claim is made here.

---

## 5. Plug-in components (existing — no build)

- **Drag force law (Tier-0).** $\gamma(v)=g\,b\,v^2$, $a_\text{drag}=-\gamma v/m$
  *evaluation* (`physics/drag.py`). Reused by Slice I⋆; not re-implemented.
- **Drag O-step (SQ1).** The frozen-$\gamma$ damping + exact dissipation in
  `physics/baoab.py`. Built and accepted (A13); reused unchanged. Slice I⋆ only feeds
  it the schedule-driven $m(t)$.
- **Tier-1a validation/reporting scripts.** Slice R now owns the run matrix,
  same-smoothed I2 `|v2|` metric, raw `R` metric, ledger diagnostic, table output,
  mean-series export, and selected-`t*` trajectory visualization. The output remains
  **reported, not adjudicated** by code; Tier-0's absolute floor does not transfer.
- **Tier-1a stress diagnostics.** Slice V owns the separate onset-strip stress
  generator/scorer surface. Its table is intentionally separate from the physical
  Tier-1a anchored-continuous RMSE table.

> **Note.** The BAOAB O-step (SQ1) *is* a plug-in above. What Slice I⋆ adds is the
> jump branch (SQ2), the post-jump mass (SQ3), and the `m(t)` plumbing into the
> existing per-step rebuild — not an O-step rewrite.

---

## 6. Dependency graph and build order

```
contract layer (data shapes + §2 constants)            ← shared, behaviourless
        │
   ┌────┼────────────┬───────────────┐
   S    M            B                (S, M, B independent — mock the rest)
   │    │
   └────┴──► I⋆ (SQ2–SQ3 + m(t) plumbing; reuses existing kicks + drag eval + SQ1 O-step)
              │
        integration (§7): S+M+I⋆+B  →  R (run matrix + RMSE table + diagnostics)
```

**Independent (parallelizable):** S, M, B — each behind its interface, others
mocked.
**Core composed build:** I⋆ (SQ2–SQ3 + `m(t)` plumbing) over S + M + the reused
fixed-mass skeleton (incl. the SQ1 O-step). 
**Plug-in:** drag force evaluation. 
**Delivered reporting:** Slice R consumes accepted checkpoints and postprocess APIs; it
does not change physics.
**Deferred:** radial depth-anchored cross-check (fire sheds at TDDFT $R$-values via
the numeric $R(t)$; build only after the time-anchored null is green).

---

## 7. Integration tests (compositions only)

- **Constant-mass regression** (`fixed`): I⋆ reproduces the legacy fixed-mass
  trajectory exactly — proves the upgrade is non-destructive.
- **Force-free anchored run** (Coulomb = drag = 0, `anchored_discrete`): velocity is
  unchanged through all 7 sheds, mass follows the 21→14 staircase, and the ledger
  closes with positive co-moving-He bookkeeping. Isolates the mass channel end-to-end.
- **Drag-only anchored run** (Coulomb = 0): trajectory changes only through the
  post-shed mass entering the drag acceleration, not through event-local speed kicks.
- **Full 1a A/B run** (Coulomb + drag + schedule), `fixed` vs `anchored_discrete`:
  produce $R(t)$, `|v2|`, ledger, RMSE table, and selected-`t*` trajectory figure.
  **Delivered result:** the cold-shed realization runs and closes its four-term
  wiring ledger, but the selected-`t*` `|v2|` plot reveals unphysical discontinuous
  velocity jumps; cold-shed is therefore only a Tier-1a diagnostic bound.

---

## 8. Config / scenario contract

> **All `SimConfig` fields named here that do not already exist are NEW additions**
> (verified against `config.py`): `coulomb_available_eV`, `anchor_mode`, `t_star_ps`,
> the shed anchors, and the `anchored_discrete` enum value. They are config-surface
> work for the build, not existing knobs.

- **`mass_scenario`** literal becomes **`{fixed, biphasic, anchored_discrete}`** —
  add `anchored_discrete`; **retire `scenario_A_accretion` and `scenario_B_stripping`**
  (superseded by `biphasic`; DESIGN §2.5/§2.8). Touch-points: the `MassScenario`
  literal, the `check_drag_config` non-`fixed` branch set (`config.py` ~417–448), and
  any preset/test referencing A/B. Retirement recorded in `drag_migration_log_tier1a.md`
  (2026-06-23/24) and DESIGN §2.8 (updated 2026-06-24: A/B recorded as retired).
- `anchor_mode` = `time` (NEW field; radial depth-anchored cross-check deferred, §6).
- Diagnostic `anchor_mode="onset_strip"` is separate from the physical
  TDDFT-anchored `time` mode. It uses `anchor_n_final` and the fixed
  `t_strip=0.5 ps` onset event for stress endpoints `n_final∈{14,2,1,0}`; it is not
  used to populate the physical Tier-1a RMSE table.
- `coulomb_available_eV` = **0.80** (validation, $d{=}9$ Å), NEW field stamped to the
  scenario tag **for provenance only**. **No hard refuse** — the value is recorded,
  not enforced as a load-time gate (DESIGN §6.5/§6.5.1 context).
- **R6 pairing:** the `anchored_discrete` run uses the constant-$m_\text{eff}$ drag
  coefficients (clean time-resolved re-extraction blocked by the 9 Å transverse flag).
  As a non-`fixed` scenario it enters the §6.5 `time_resolved`-requiring arm, trips the
  guard structurally → runs under `allow_inconsistent_mass_pairing = True`, on the §6.6
  mid-window defence: anchored $m\approx19$ He ($=m_\text{eff}$) mid-window → near-zero
  error; the $n{=}21$/$n{=}14$ ends ($\sim\tfrac13$ extrapolation) sit in the §6.7
  free-zone.
- `t_star_ps` ($t^*$): NEW field; schedule onset, $n=21$ held for $t\le t^*$. Tier 1a
  **sweeps** $t^*\in\{0.5,5,9\}$ ps (wide span: extremes + midpoint) to expose the
  trajectory's sensitivity to shed timing under full forces. The segment-2 event times
  and the force-free telescoping boost (1.153) are $t^*$-independent (§10); the two
  segment-1 event *times* and the full-force trajectory are not. **Run matrix:** one
  `fixed` null (t*-independent, no sheds) + three `anchored_discrete` runs, one per
  swept $t^*$.
- **Tail:** flat $n=14$ past 14 ps (no data → no anchored extrapolation). One-signed:
  if the true cascade continues (R5, terminal $n$ an upper bound), late sheds are
  omitted, so **1.153 is a floor on the boost, not a point estimate.**
- **`IonCheckpoint` schema bump v5 → v6** (current is v5, `_ION_SCHEMA_VERSION`). The
  four-term ledger arrays **already exist** — `E_kin_eV`, `E_pot_eV`, `E_dissip_eV`,
  all `(2N,T)`. v6 delta: (i) **rename** `E_mass_attach_defect_eV → E_mass_transfer_eV`
  (same `(2N,T)`, sign now covers shedding; DESIGN §2.9); (ii) **add** `n_shell`
  shell-count `(2N, num_steps)` (per-atom, matching the existing convention and load
  shape-validation); (iii) **drop** the `mass_history_kg` non-decreasing assumption
  (comment + any guard); (iv) add a `mass_scenario` metadata field. **No
  $E_\text{int}$ field at 1a.** Slice B's ledger reuses
  `postprocess/energy_balance.py:ion_energy_totals`, swapping the defect term for
  `E_mass_transfer_eV`. **Back-compat load shim:** the v6 loader accepts legacy v5
  checkpoints — maps `E_mass_attach_defect_eV → E_mass_transfer_eV` and synthesizes an
  absent `n_shell` (constant at the run's fixed shell count) — so the existing v5
  `ion.npz` run dirs (incl. the Tier-0 `shared_pure_cubic` runs) still load rather
  than failing the version check.
- **Run size:** no new decision — reuse the existing standard run preset, invoked with
  `mass_scenario=anchored_discrete` (and `fixed` for the null). The deterministic shed
  schedule applies identically to all atoms.

---

## 9. Acceptance criteria

> **Scientific deliverable (not a code-side gate).** The end product of a Tier-1a
> run is a **reported diagnostic package** — `fixed` null plus the three
> `anchored_discrete` $t^*$ runs, a rich RMSE table, mean-series exports, and a
> selected-`t*` `|v2|` figure — **emitted and left to the user**. The code asserts
> **no** verdict on which mode is better. The criteria below are
> **wiring-correctness / reporting-completeness gates only**, not the scientific
> conclusion.

**Per slice:** the slice's own suite (§4) green with all other slices mocked.

**SQ upgrade (Slice I⋆):**

- Constant-mass regression exact — the **SQ1-untouched guard** (`fixed` reproduces
  the current fixed-mass trajectory bit-for-bit; proves the reused O-step is unchanged).
- New code demonstrated on analytic limits: SQ2 (jump-then-O, ≤1/step,
  continuous-velocity shed and positive co-moving-He transfer), SQ3 (post-jump $m^+$);
  SQ1 (reused) re-confirmed under a running $m$. Cold-shed reset tests are retained
  as diagnostic-bound coverage only.

**Integration:**

- Four-term ledger closes to the Verlet-drift bound in both modes.
- Relabel-instead-of-reset fault caught by Slice B.
- The $t^*$ sweep ($\{0.5,5,9\}$ ps) **runs to completion and emits the RMSE table**
  plus the selected-`t*` `|v2|` figure (Slice R). This is a reporting gate, **not**
  an automated comparison verdict.
- **Post-Slice-C scientific conclusion:** the force-free telescoping boost (1.153)
  and the cold-shed velocity reset are diagnostic-bound only. The physical
  anchored-mass Tier-1a comparison is the continuous-velocity mass update; cold-shed
  waits for the full `E_int` / energy-gated evaporation tier.
- **Stress diagnostics:** onset-strip runs are reported in their own stress table
  and plots. They are sensitivity checks, not physical Tier-1a verdict rows.

**Out-of-scope guard:** any code path that reads a $D_0$ rung, $\kappa$, $\nu$, $s$,
the electronic picture, an $E_\text{int}$ value, or a continuous mass-loss term
**fails review** — bypassed/removed at 1a by construction.

---

## 10. Test-oracle fixtures (golden values)

**Masses (amu):** $m@n=126.90+n\cdot4.0026$. $m@21=210.955$,
$m@19=202.954\;(\equiv m_\text{eff})$, $m@14=182.936$.

**Shed event table** (pre-shed $m$; kick $=m/(m-m_\text{He})$):

| event | crossing | time (ps)          | pre-shed $m$ | kick   |
| ----- | -------- | ------------------ | ------------ | ------ |
| 21→20 | 20.5     | $t^*+0.25(10-t^*)$ | 210.955      | 1.0193 |
| 20→19 | 19.5     | $t^*+0.75(10-t^*)$ | 206.952      | 1.0197 |
| 19→18 | 18.5     | 10.4               | 202.954      | 1.0201 |
| 18→17 | 17.5     | 11.2               | 198.951      | 1.0205 |
| 17→16 | 16.5     | 12.0               | 194.949      | 1.0210 |
| 16→15 | 15.5     | 12.8               | 190.946      | 1.0214 |
| 15→14 | 14.5     | 13.6               | 186.939      | 1.0219 |

**Cumulative boost:** $\prod=210.955/182.936=1.1532$ (telescoping; timing- and
tie-break-independent).

---

## 11. Interpretation guardrails (carry into every run report)

- **No event-local speed kicks in production Tier 1a.** The current
  `anchored_discrete` driver path keeps `v^+=v^-` at each shed; any trajectory change
  comes from the lower post-shed mass in subsequent Coulomb/drag dynamics. The
  cold-shed reset raises $|v|$ discontinuously ($\propto1/m$) and remains a
  **kick-only upper-bound diagnostic**, not the physical Tier-1a mass-dynamics model.
- **Transverse contamination is common-mode.** The non-radial 9 Å signal inflates
  RMSE($R$) and RMSE($|v|$) for both modes roughly equally (same Coulomb + drag, only
  mass differs), so it cancels in the **ranking**. A large absolute RMSE is not a
  failed run; the verdict is the relative order.
- **One-signed tail.** Flat $n=14$ past 14 ps omits any continuing cascade. This
  remains a limitation of the anchored schedule.
- **Closure is wiring, not physics.** Four-term closure is by construction once the
  SQ2 bookkeeping matches the selected shed primitive; it certifies plumbing. The
  full cold-shed evaporation interpretation still requires the `E_int` reservoir
  and five-term invariant of the later energy-gated tier.
- **Verdict scope.** There is **no automated pass/fail on fidelity.** The run emits
  the §5/§9 RMSE table (`fixed` vs `anchored_discrete` across the $t^*$ sweep) and a
  selected-`t*` `|v2|` figure for user evaluation. What the delivered table can now
  speak to is the continuous-velocity anchored-mass response, while the cold-shed
  diagnostic bound remains separate and the biphasic generative mechanism remains
  unfalsified until Tier 2.
- **Onset-strip stress scope.** `onset_strip` is an intentionally separate diagnostic
  family: stress rows stay out of the physical anchored-continuous RMSE table, the
  fixed null appears only as a comparison baseline for plots/exports, and production
  stress artifacts must be generated before interpreting any stress result.
