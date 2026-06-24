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
trajectory relative to a constant-mass ion, and if so, in the direction and
magnitude the cold-shed momentum signature predicts?*

The He shell schedule $n(t)$ is **anchored** to the 9 Å TDDFT loss curve (read in,
not generated), so the kinetics — electronic picture, ladder $D_0$, $\kappa$, $\nu$,
$s$ — **never enter**. With the schedule fixed, the only thing that varies between
the two runs is whether the ion's mass is held constant or follows $m(t)$. Tier 1a
is therefore a **controlled A/B comparison**:

- **`fixed`** — constant $m_\text{eff}$ (the null).
- **`anchored_discrete`** — anchored variable $m(t)$, losing one He per shed event with a
  momentum-conserving cold-shed reset.

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
**SQ2** (the momentum-conserving mass-jump operator), **SQ3** (the post-jump O-step
reading `m⁺`), and the **`m(t)` plumbing** that feeds the schedule-driven mass into
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
  move relative to the `fixed` null. The force-free $|v|\propto 1/m$ telescoping boost
  (1.153) is a **sanity ceiling** the force-free integration test must hit — **not**
  the verdict (it is timing-insensitive only when drag and Coulomb are off; see §2).

**Does NOT test (out of scope by agreement — do not let these leak in):**

- The biphasic **generative mechanism** (Poisson pickup + RRK + self-bound gate).
  Anchoring the schedule bypasses it; unfalsified until Tier 2.
- The **internal-energy reservoir** $E_\text{int}$. Dropped for 1a. Consequence: the
  shed velocity boost is energetically *unsourced* in the ledger — acceptable here
  precisely because the generative mechanism is out of scope.
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

**Cold-shed reset (momentum-conserving, $u_\text{He}=0$).**
$$v^+=\frac{m}{m-m_\text{He}}\,v^-$$
Shed He leaves at rest, so complex momentum $mv$ is **invariant across the
instantaneous jump**; the velocity vector is scaled, direction preserved. Ratio ×
Å/ps = Å/ps ✓.

**Force-free telescoping ceiling (sanity check, NOT the verdict).** Because $mv$ is
conserved per shed, the cumulative *free-flight* **speed** boost collapses to the
endpoint mass ratio:
$$\prod_{k}\frac{m_k}{m_k-m_\text{He}}=\frac{m(t^*)}{m(\text{end})}=\frac{210.955}{182.936}=\mathbf{1.153}$$
i.e. $|v|\propto 1/m$, with $R(t)$ its integral. This **~15.3 %** boost is
independent of the shed count and placement (S4) **only in the force-free limit** —
it is the ceiling the force-free integration test (§7) must hit. Under full
drag+Coulomb the boost is *eroded* by the $1/m(t)$ drag-acceleration and the firing
*times* feed back through the $\propto v^3$ drag work, so the trajectory $R(t)$,
$|v(t)|$ **and** the endpoint become genuinely $t^*$-sensitive. The Tier-1a result is
therefore the **$t^*$ sweep** ($\{0.5,5,9\}$ ps, §8) showing that timing response —
not the single telescoping number.

**Mass-transfer ledger term (reduced-mass defect).** Per shed the complex KE rises
by the reduced-mass defect; with $E_\text{int}$ dropped it is booked as a negative
increment so the ledger closes by construction:
$$\Delta E_\text{mass\_transfer}=-\tfrac12\,\frac{m\,m_\text{He}}{m-m_\text{He}}\,\lVert v^-\rVert^2$$
amu·Å²/ps² = energy ✓. **Exact reduced-mass form, not the heavy-ion approximation
$\tfrac12 m_\text{He}v^2$** (~3 % closure error at $n{=}1$; `MASS` SQ2). A label-only
$v$ — relabelling mass without the reset — voids invariant closure.

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
- **SQ2 — momentum-conserving jump, new.** The reset + reduced-mass defect above —
  the **invariant-closure precondition** (A13). Ordering inside the step:
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

**Testing philosophy (every new slice):**

- Pure functions (S) tested against closed-form oracle values (§10).
- Stateful modules (M, B) tested on analytic micro-cases (single shed, synthetic
  stream) with neighbours **mocked**.
- I⋆ tested on analytic limits (constant-mass reduction, single-jump) before any
  real run.
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
> (Slice M record, 2026-06-24).

**Purpose.** Hold $m(t)$ and perform the SQ2 cold-shed reset; emit the increment the
ledger needs and the $m^+$ the integrator needs (SQ3).

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
> injection). Full suite 841/0. **Deferred to a later plan:** the §5/§9
> `t*∈{0.5,5,9}` RMSE-table *run* deliverable. Delivery detail:
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

## 5. Plug-in components (existing — no build)

- **Drag force law (Tier-0).** $\gamma(v)=g\,b\,v^2$, $a_\text{drag}=-\gamma v/m$
  *evaluation* (`physics/drag.py`). Reused by Slice I⋆; not re-implemented.
- **Drag O-step (SQ1).** The frozen-$\gamma$ damping + exact dissipation in
  `physics/baoab.py`. Built and accepted (A13); reused unchanged. Slice I⋆ only feeds
  it the schedule-driven $m(t)$.
- **Validation scripts.** RMSE($R$), RMSE($|v|$) vs the smoothed TDDFT curve, per
  mode. The output is an **RMSE table emitted for user evaluation — no automated
  verdict** (no pass/fail threshold, and the relative ranking is *reported, not
  adjudicated* by the code; Tier-0's absolute floor does not transfer). Report RMSE
  over the full window *and* the post-free-zone window (excluding the first several
  ps, §6.7) as a **reporting split**, not an acceptance gate. Plugged in as-is.

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
        integration (§7): S+M+I⋆+B  →  validation scripts (plug-in)
```

**Independent (parallelizable):** S, M, B — each behind its interface, others
mocked.
**Core composed build:** I⋆ (SQ2–SQ3 + `m(t)` plumbing) over S + M + the reused
fixed-mass skeleton (incl. the SQ1 O-step). 
**Plug-in:** drag force evaluation, validation. 
**Deferred:** radial depth-anchored cross-check (fire sheds at TDDFT $R$-values via
the numeric $R(t)$; build only after the time-anchored null is green).

---

## 7. Integration tests (compositions only)

- **Constant-mass regression** (`fixed`): I⋆ reproduces the legacy fixed-mass
  trajectory exactly — proves the upgrade is non-destructive.
- **Force-free anchored run** (Coulomb = drag = 0, `anchored_discrete`): $|v(t)|$ a pure
  7-step staircase; endpoint speed = initial × 1.153; ledger closes. Isolates the
  mass channel end-to-end.
- **Drag-only anchored run** (Coulomb = 0): the ~15.3 % kick ceiling is *eroded* by
  the $1/m(t)$ drag-acceleration; net boost $<15.3\%$ and monotone in drag strength.
- **Full 1a A/B run** (Coulomb + drag + schedule), `fixed` vs `anchored_discrete`: produce
  $R(t)$, $|v(t)|$, ledger; hand to validation. Modes separated by the boost; closure
  holds in both.

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
> run is a **single RMSE evaluation table** — `fixed` null plus the three
> `anchored_discrete` $t^*$ runs, RMSE($R$) and RMSE($|v|$) vs smoothed 9 Å TDDFT,
> full-window *and* post-free-zone columns (§5, §6.7) — **emitted and left to the
> user**. The code asserts **no** verdict on which mode is better. The criteria
> below are **wiring-correctness gates only** (build acceptance), not the scientific
> conclusion.

**Per slice:** the slice's own suite (§4) green with all other slices mocked.

**SQ upgrade (Slice I⋆):**

- Constant-mass regression exact — the **SQ1-untouched guard** (`fixed` reproduces
  the current fixed-mass trajectory bit-for-bit; proves the reused O-step is unchanged).
- New code demonstrated on analytic limits: SQ2 (jump-then-O, ≤1/step, reduced-mass
  reset/defect), SQ3 (post-jump $m^+$); SQ1 (reused) re-confirmed under a running $m$.

**Integration:**

- Four-term ledger closes to the Verlet-drift bound in both modes.
- Relabel-instead-of-reset fault caught by Slice B.
- The $t^*$ sweep ($\{0.5,5,9\}$ ps) **runs to completion and emits the RMSE table**
  (§5 deliverable) — `anchored_discrete` vs the `fixed` null in $R(t)$/$|v(t)|$ across
  shed timings. This is a *production gate* (the table is produced), **not** an
  automated comparison verdict — the ranking is left to the user. The force-free
  telescoping boost (1.153) is a *ceiling* — eroded by drag and a *floor* given the
  flat-$n=14$ tail (opposing effects) — reported as a sanity bound, not the verdict.

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

- **Opposing effects.** Shed kicks raise $|v|$ ($\propto 1/m$); the $1/m(t)$
  drag-acceleration lowers it. "Reproduces $R,|v|$" tests the **superposition**, not
  the kicks alone. The 1.153 boost is a *kick-only ceiling*.
- **Transverse contamination is common-mode.** The non-radial 9 Å signal inflates
  RMSE($R$) and RMSE($|v|$) for both modes roughly equally (same Coulomb + drag, only
  mass differs), so it cancels in the **ranking**. A large absolute RMSE is not a
  failed run; the verdict is the relative order.
- **One-signed tail.** Flat $n=14$ past 14 ps omits any continuing cascade ⇒ boost is
  a **floor**.
- **Closure is wiring, not physics.** Four-term closure is by construction once the
  SQ2 reset is correct; it certifies plumbing, not the energetics of shedding
  (unsourced at 1a).
- **Verdict scope.** There is **no automated pass/fail on fidelity.** The run emits
  the §5/§9 RMSE table (`fixed` vs `anchored_discrete` across the $t^*$ sweep) and the
  **user evaluates it** — the code's only gates are wiring-correctness (§9). What the
  table can speak to is the **influence of mass dynamics on the kinematics** (does
  $R(t)$/$|v(t)|$ move off the null, and how with shed timing), *not* the biphasic
  generative mechanism (unfalsified until Tier 2). The telescoping boost is a
  force-free sanity ceiling, not the verdict.
