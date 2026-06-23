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
- **`discrete`** — anchored variable $m(t)$, losing one He per shed event with a
  momentum-conserving cold-shed reset.

The comparison isolates the *influence of mass dynamics* on $R(t)$, $|v(t)|$, and on
whether the kinematic energy/momentum bookkeeping closes. It does **not** test the
generative mechanism that *produces* the schedule — anchored away here, arbitrated
by experiment at Tier 2.

**Central engineering fact.** The existing BAOAB ion integrator was built for
**fixed mass**. Making it correct under a time-varying, jump-discontinuous mass is
the substance of this build, and that is exactly what SQ1–SQ3 specify. SQ1–SQ3 are
**not pre-done** — they are the variable-mass upgrade.

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
- The **mass-model class verdict**: variable mass (`discrete`) vs constant (`fixed`),
  via the $|v|\propto 1/m$ cold-shed signature.

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

| Symbol | Value | Units | Source |
|---|---|---|---|
| $m_{\mathrm I^+}$ | 126.90 | amu | atomic |
| $m_\text{He}$ | 4.0026 | amu | atomic |
| $b$ | 2.5154 | amu·ps/Å² | Tier-0 locked (pure cubic) |
| $g$ | $\in[0,1]$ | dimensionless | density-proportional gate (erf complement) |
| $E_\text{avail}^\text{ion}$ | **0.80** | eV | validation, $d{=}9$ Å, $\tfrac12\cdot14.40/9$ |
| anchors | $(t^*,21),(10,19),(14,14)$ | ps, count | 9 Å TDDFT (author-confirmed) |
| $m_\text{eff}$ | 202.95 | amu | 19 He (drag-extraction reference) |

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

**Telescoping class discriminator (timing-insensitive).** Because $mv$ is conserved
per shed, the cumulative free-flight **speed** boost collapses to the endpoint mass
ratio:
$$\prod_{k}\frac{m_k}{m_k-m_\text{He}}=\frac{m(t^*)}{m(\text{end})}=\frac{210.955}{182.936}=\mathbf{1.153}$$
i.e. $|v|\propto 1/m$, with $R(t)$ its integral. This **~15.3 %** boost is
independent of the number of intermediate steps and of their placement (S4), so the
*class verdict* is robust to the schedule; only the *trajectory-fit phase* is
timing-sensitive.

**Mass-transfer ledger term (reduced-mass defect).** Per shed the complex KE rises
by the reduced-mass defect; with $E_\text{int}$ dropped it is booked as a negative
increment so the ledger closes by construction:
$$\Delta E_\text{mass\_transfer}=-\tfrac12\,\frac{m\,m_\text{He}}{m-m_\text{He}}\,\lVert v^-\rVert^2$$
amu·Å²/ps² = energy ✓. **Exact reduced-mass form, not the heavy-ion approximation
$\tfrac12 m_\text{He}v^2$** (~3 % closure error at $n{=}1$; `MASS` SQ2). A label-only
$v$ — relabelling mass without the reset — voids invariant closure.

### The SQ upgrade — fixed-mass BAOAB → variable-mass-correct

These three are the **work**, not a precondition already met. The integrator
currently assumes a constant $m$; all three must hold once $m(t)$ varies and jumps.

- **SQ1 — drag O-step under varying $m(t)$.** Freeze $\gamma$ at the step-entry
  velocity and apply $e^{-\gamma\,dt/m}$ using the **instantaneous** $m(t)$ (the
  Tier-0 frozen-$v_\text{in}$ result, now carried under a changing mass). Buys exact
  dissipation at any $dt$,
  $\Delta E_\text{dissip}=\tfrac12 m(\lVert v_\text{in}\rVert^2-\lVert v_\text{out}\rVert^2)$,
  and unconditional dissipativity ($\lVert v_\text{out}\rVert\le\lVert v_\text{in}\rVert$).
  $O(dt)$ on the cubic (second-order claim retired). Residual: a one-signed
  over-braking bias (R10).
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
  *evaluation*. Tested. (Its use inside the O-step under variable mass is SQ1 — new.)
- **Validation scripts**: RMSE of $R(t)$, $|v(t)|$ vs the smoothed reference, per
  mode, relative ranking. Existing; plugged in, no work.

**New (build + test) — §4 slices:**
- **Slice S** — schedule generator (anchored $\bar n(t)$ → 7 shed events).
- **Slice M** — mass-jump operator and $m(t)$ state (the SQ2 reset + defect; emits
  $m^+$).
- **Slice I⋆** — **variable-mass integrator upgrade (SQ1–SQ3)**: modify the
  fixed-mass BAOAB so the drag O-step uses $m(t)$ with exact dissipation (SQ1), the
  mass jump enters jump-then-O at ≤1/step (SQ2 ordering, calling M), and the
  post-jump O-step uses $m^+$ (SQ3).
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

**Purpose.** Turn the three anchors into the deterministic, monotone, sheds-only
integer schedule and its 7 fire events. No physics, no state.

**Interface.** Consumes anchors + $t^*$; emits (a) the piecewise-linear $\bar n(t)$
evaluator and (b) the ordered event list (fire-time, $n\to n-1$, pre-/post-shed
mass, kick factor).

**Encoded form.**
$$\bar n(t)=\begin{cases}21 & t\le t^*\\ 21-\dfrac{2(t-t^*)}{10-t^*} & t^*<t\le10\\ 19-1.25\,(t-10) & 10<t\le14\\ 14 & t>14\end{cases}$$
Segment loss rates $2/(10-t^*)$ and $1.25$ ps⁻¹ — a real ~5× acceleration kept as a
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

**Purpose.** Hold $m(t)$ and perform the SQ2 cold-shed reset; emit the increment the
ledger needs and the $m^+$ the integrator needs (SQ3).

**Interface.** Consumes a (mockable) fire event + current state; emits
$m^+=m-m_\text{He}$, $v^+=\tfrac{m}{m-m_\text{He}}v^-$, and
$\Delta E_\text{mass\_transfer}$. Two run modes:
- `fixed` — $m\equiv m_\text{eff}$, no reset, zero defect (null).
- `discrete` — reset + defect on each fire event.

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

### Slice I⋆ — Variable-mass integrator upgrade (SQ1–SQ3) *(core build; composes S, M)*

**Purpose.** Convert the fixed-mass BAOAB ion step into a variable-mass-correct step.
This is the central modification, not a seam on untouched code.

**Interface.** Per step: conservative B/A kicks (existing `acc_fn`: Coulomb +
droplet) → **query schedule** → if a fire event, apply the **M** jump
(**jump-then-O**, ≤1/step, **SQ2**) → drag O-step with frozen $\gamma$ at $v_\text{in}$
on the **current/post-jump mass** (**SQ1 + SQ3**). Emits next state, per-step
$\Delta E_\text{dissip}$ (SQ1), and $\Delta E_\text{mass\_transfer}$ (from M).

**What changes from the fixed-mass integrator.**
- The $m$ in the O-step damping $e^{-\gamma\,dt/m}$ and in $\Delta E_\text{dissip}$
  becomes the instantaneous $m(t)$ (SQ1).
- A jump branch is inserted before the O-step, calling M (SQ2 ordering).
- The post-jump O-step reads $m^+$ (SQ3).

**Independence.** Composes S (events) and M (reset), both mockable; the existing
conservative kicks and the drag-force evaluation are reused unchanged.

**Test spec.**
- **Constant-mass reduction:** with `fixed` mode (no jumps), the upgraded step
  reproduces the existing fixed-mass integrator bit-for-bit (regression guard).
- **SQ1:** drag-only decay (mock M no-op, Coulomb 0) follows the exact frozen-$v$
  exponential at the running $m$; unconditional stability at large $dt$.
- **Exact dissipation identity** at $dt\in\{0.001,0.01,0.05\}$ ps.
- **SQ2 ordering:** a fire event applies the reset **before** the O-step; ≤1 event
  per step even under a forced-dense schedule.
- **SQ3:** the post-jump O-step damping uses $m^+$, not $m^-$.
- **Jump-step measure-zero:** refining $dt$ leaves the event *count* fixed.

**Acceptance.** Constant-mass regression exact; SQ1–SQ3 all demonstrated on the
analytic limits above.

---

### Slice B — Four-term bookkeeping & closure *(stateful; synthetic-stream testable)*

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
  *evaluation*. Reused by Slice I⋆; not re-implemented. (Its O-step application under
  $m(t)$ is SQ1 — new, in I⋆.)
- **Validation scripts.** RMSE($R$), RMSE($|v|$) vs the smoothed TDDFT curve, per
  mode, with the relative ranking as the verdict (no pass/fail threshold — Tier-0's
  absolute floor does not transfer). Report RMSE over the full window *and* the
  post-free-zone window (excluding the first several ps, §6.7) as a **reporting
  split**, not an acceptance gate. Plugged in as-is.

> **Note.** The BAOAB integrator is *not* in this list. It exists for fixed mass but
> is **modified** by Slice I⋆; SQ1–SQ3 are that modification.

---

## 6. Dependency graph and build order

```
contract layer (data shapes + §2 constants)            ← shared, behaviourless
        │
   ┌────┼────────────┬───────────────┐
   S    M            B                (S, M, B independent — mock the rest)
   │    │
   └────┴──► I⋆ (variable-mass upgrade, SQ1–SQ3; reuses existing kicks + drag eval)
              │
        integration (§7): S+M+I⋆+B  →  validation scripts (plug-in)
```

**Independent (parallelizable):** S, M, B — each behind its interface, others
mocked.
**Core composed build:** I⋆ (SQ1–SQ3) over S + M + the reused fixed-mass skeleton. 
**Plug-in:** drag force evaluation, validation. 
**Deferred:** radial depth-anchored cross-check (fire sheds at TDDFT $R$-values via
the numeric $R(t)$; build only after the time-anchored null is green).

---

## 7. Integration tests (compositions only)

- **Constant-mass regression** (`fixed`): I⋆ reproduces the legacy fixed-mass
  trajectory exactly — proves the upgrade is non-destructive.
- **Force-free anchored run** (Coulomb = drag = 0, `discrete`): $|v(t)|$ a pure
  7-step staircase; endpoint speed = initial × 1.153; ledger closes. Isolates the
  mass channel end-to-end.
- **Drag-only anchored run** (Coulomb = 0): the ~15.3 % kick ceiling is *eroded* by
  the $1/m(t)$ drag-acceleration; net boost $<15.3\%$ and monotone in drag strength.
- **Full 1a A/B run** (Coulomb + drag + schedule), `fixed` vs `discrete`: produce
  $R(t)$, $|v(t)|$, ledger; hand to validation. Modes separated by the boost; closure
  holds in both.

---

## 8. Config / scenario contract

- `coulomb_available_eV` = **0.80** (validation, $d{=}9$ Å), stamped to the scenario
  tag; a 2.70 eV onset must **not** run against these references (DESIGN §6.5/§6.5.1).
- `mass_scenario` ∈ {`fixed`, `discrete`}; `anchor_mode` = `time` (radial
  depth-anchored cross-check deferred).
- **R6 pairing:** the `discrete` run uses the constant-$m_\text{eff}$ drag
  coefficients (clean time-resolved re-extraction blocked by the 9 Å transverse
  flag). It trips the §6.5 guard structurally → `allow_inconsistent_mass_pairing =
  True`, on the §6.6 mid-window defence: anchored $m\approx19$ He ($=m_\text{eff}$)
  mid-window → near-zero error; the $n{=}21$/$n{=}14$ ends ($\sim\tfrac13$
  extrapolation) sit in the §6.7 free-zone.
- `t_star` ($t^*$): schedule onset = the surface-crossing transient cut; $n=21$ held
  for $t\le t^*$.
- **Tail:** flat $n=14$ past 14 ps (no data → no anchored extrapolation). One-signed:
  if the true cascade continues (R5, terminal $n$ an upper bound), late sheds are
  omitted, so **1.153 is a floor on the boost, not a point estimate.**
- `IonCheckpoint` schema: add $n(t)$ + the four-term ledger fields; **no
  $E_\text{int}$ field at 1a**. Drop the mass-monotonicity guarantee.

---

## 9. Acceptance criteria

**Per slice:** the slice's own suite (§4) green with all other slices mocked.

**SQ upgrade (Slice I⋆):**
- Constant-mass regression exact (non-destructive upgrade).
- SQ1 (frozen-$v$ O-step under $m(t)$ + exact dissipation), SQ2 (jump-then-O, ≤1/step,
  reduced-mass reset/defect), SQ3 (post-jump $m^+$) all demonstrated on analytic
  limits.

**Integration:**
- Four-term ledger closes to the Verlet-drift bound in both modes.
- Relabel-instead-of-reset fault caught by Slice B.
- Class verdict separates `fixed` from `discrete` on the boost; the boost is a
  *floor* (tail) and is *eroded by drag* (opposing effects) — both stated in report.

**Out-of-scope guard:** any code path that reads a $D_0$ rung, $\kappa$, $\nu$, $s$,
the electronic picture, an $E_\text{int}$ value, or a continuous mass-loss term
**fails review** — bypassed/removed at 1a by construction.

---

## 10. Test-oracle fixtures (golden values)

**Masses (amu):** $m@n=126.90+n\cdot4.0026$. $m@21=210.955$,
$m@19=202.954\;(\equiv m_\text{eff})$, $m@14=182.936$.

**Shed event table** (pre-shed $m$; kick $=m/(m-m_\text{He})$):

| event | crossing | time (ps) | pre-shed $m$ | kick |
|---|---|---|---|---|
| 21→20 | 20.5 | $t^*+0.25(10-t^*)$ | 210.955 | 1.0193 |
| 20→19 | 19.5 | $t^*+0.75(10-t^*)$ | 206.952 | 1.0197 |
| 19→18 | 18.5 | 10.4 | 202.954 | 1.0201 |
| 18→17 | 17.5 | 11.2 | 198.951 | 1.0205 |
| 17→16 | 16.5 | 12.0 | 194.949 | 1.0210 |
| 16→15 | 15.5 | 12.8 | 190.946 | 1.0214 |
| 15→14 | 14.5 | 13.6 | 186.939 | 1.0219 |

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
- **Verdict scope.** A pass validates the **influence of mass dynamics on the
  kinematics** — *not* the biphasic generative mechanism, which remains unfalsified
  until Tier 2.
