# Implementation — The Ion-Stage Pipeline

> **What this document is.** How the ion-stage physics is actually built in
> code: the module map, the unit and mass contracts, what happens in one
> timestep, where the seams are, and what each design decision costs. It is the
> methods-chapter companion to `PROJECT_OVERVIEW.md` (the narrative) and
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (the mechanism physics).
>
> **Why it exists.** The two design documents are implementation-agnostic *by
> project rule* — they carry LaTeX and decision registers, never code. And
> module-level documentation under `docs/` all but stops at the pre-drag-port
> baseline: `drag.py` and `baoab.py` have entries, but the thirteen modules
> that implement the mass mechanism and the post-MD stages have **none**. This
> document covers that gap.
>
> **Scope.** The ion stage only. The neutral stage (velocity-Verlet, unchanged
> from the MATLAB port) and the post-processing layer are out of scope.

---

## 1. Module map

Everything below the line was built for the drag/mass port. The two starred
entries have `docs/physics/` pages; the rest are documented only here and
in their own docstrings.

| module | role | lines |
|---|---|---|
| `physics/leapfrog.py` | shared `_kick` / `_drift` primitives, acceleration closures | 423 |
| `physics/potentials.py`, `interactions.py`, `constants.py` | conservative forces, unit constants | — |
| — | — | — |
| `physics/drag.py` ★ | the drag laws: `γ(v)`, `F_drag(v)`, the coefficient bundle | 512 |
| `physics/_gates.py` | the shared erf-complement gate primitive | 61 |
| `physics/baoab.py` ★ | the B-A-O-A-B ion stepper; the O-step | 207 |
| `physics/mass_jump.py` | shed / capture velocity-and-mass operators | 587 |
| `physics/pickup.py` | Poisson He-capture channel | 408 |
| `physics/evaporation.py` | energy-gated RRK evaporation channel | 517 |
| `physics/dissociation_ladder.py` | Form-U `D₀(n)` and `Σ(n)`; tabulated ladders | 408 |
| `physics/internal_energy_budget.py` | the per-event `E_int` arithmetic (S1/S2/K1) | 207 |
| `physics/solvation_cooling.py` | Newton cooling K2; the binding split | 163 |
| `physics/helium_density.py` | the `ρ_He/ρ_bulk` occupancy gate | 161 |
| `physics/shell_schedule.py` | the Tier-1a anchored `n(t)` schedule | 326 |
| `physics/state_coupling.py` | the `s(n)` drag-state coupling probe | 213 |
| `physics/exit_strip.py` | the depth-graded exit-strip knock probabilities | 120 |
| `simulation/ion_propagation_step.py` | **the spine** — one-step orchestration | 1150 |
| `simulation/ion.py` | the driver: dispatch, loop, storage stride | 567 |
| `simulation/relaxation_stage.py` | E2, fixed-dt post-ejection continuation | 635 |
| `simulation/detection_stage.py` | the exact event-driven µs cascade | 996 |
| `simulation/checkpoint.py` | schema v8, load shims | 711 |

Test coverage for these lives in ~30 focused modules
(`test_drag.py`, `test_baoab.py`, `test_mass_jump.py`, `test_pickup.py`,
`test_evaporation.py`, `test_biphasic_step.py`, `test_detection_stage.py`, …)
out of 115 in the suite.

---

## 2. Two contracts that shape everything

### 2.1 The friction convention: γ is a force coefficient

`γ(v)` has units **amu/ps** and is defined `γ(v) = |F_drag(v)| / v`, so the
friction force is `γ(v)·v` with **no leading mass**. This is unusual — most
Langevin codes carry a friction *rate* — and it was chosen deliberately.

The consequence is that **`physics/drag.py` never takes a mass**. It is fully
mass-agnostic: 512 lines of drag physics with no `m` anywhere in the signature.
Mass enters the model at exactly **one** site, the BAOAB O-step, as one
division inside the damping exponent `e^(−γ·dt/m)`.

This is what made the variable-mass work tractable. When Tier 1a introduced
`m(t)` and Tier 2 made it stochastic and per-ion, the drag module needed no
change at all — only the integrator's mass argument moved.

### 2.2 Units: mechanical amu inside, eV at the boundary

`physics/` runs in **mechanical amu**: `m` in amu, `v` in Å/ps, energy in
`amu·Å²/ps²`. There is no kg and no eV below the simulation layer.
`simulation/ion_propagation_step.py` owns the conversion, through a single
shared helper:

```python
def _amu_ang2_ps2_to_eV(dE):
    return dE * U * (100.0 ** 2) / EV
```

Every drag-path energy booking routes through it — the O-step's dissipation,
the shed defect, the channel mass-transfer terms. The operation order is pinned
to the pre-refactor inline sequence so that Tier-0 and Tier-1a outputs stay
byte-identical across the refactor that introduced it.

**Why this matters for the thesis:** the ledger closes to Verlet drift rather
than to float noise, and that is only true because the conversion happens once,
in one place, in one order.

---

## 3. One timestep, end to end

The driver (`simulation/ion.py`) dispatches once per run on whether a validated
drag-coefficient bundle exists:

```python
use_drag = cfg.drag_coefficients is not None
```

If not, the legacy hard-sphere collision path runs unchanged. If so, the ion
stage runs BAOAB. Under the production `biphasic` scenario one timestep is:

```
 1. biphasic_step(state, rng, cfg, ...)          ← the mass/energy seam
      a. K2 Newton cooling of E_int
      b. evaporation draw  (frozen FIRST)
      c. pickup draw       (frozen SECOND, on the post-shed state)
 2. rebuild acc_fn and the BAOAB closure at the post-event mass m⁺
 3. baoab_propagation_step(...)   →  B A O A B, positions and time advance
 4. fold e_bind_pair(n) into E_pot
 5. exit_strip_step(...)          ← only at outbound droplet crossings
 6. every stride-th step: write the state to a checkpoint column
```

Three structural decisions are visible in that ordering and each is load-bearing.

### 3.1 The mass event happens at the step seam, not inside the step

`biphasic_step` performs **no position or time advance**. It mutates
`(v, m, n, E_int, E_dissip, E_mass_transfer)` at the current time and returns.
The driver then rebuilds the integrator closure, so BAOAB reads the **post-jump
mass `m⁺`** for both the conservative kicks and the O-step damping exponent.

This is the SQ2/SQ3 pair from the Tier-1a plan. The alternative placement —
`B/A → jump → O` inside the step — was rejected because it puts the first
half-kick at the pre-jump mass. The cost of the chosen placement is an `O(dt)`
order reduction in that same half-kick, which is benign because **the number of
jumps is `dt`-independent**: as `dt → 0` the jump steps become measure-zero.

The BAOAB closure is therefore rebuilt **every step**, which looks wasteful and
is the single largest cost in the inner loop. It exists so that mass can change.

### 3.2 The draw order is frozen: shed first, then pickup

Both channels draw exactly one vectorised uniform per step:

```python
n_e, m_e, ... = evaporation_step_components(rng=rng, ...)   # draw 1
n_p, m_p, ... = pickup_step_components(rng=rng, n=n_e, ...) # draw 2, post-shed
```

The evaporation draw runs on the **post-cooling** `E_int`, so the self-bound
gate opens in the same step in which cooling drains `E_int` below `Σ(n)` — not
one step later. The pickup draw runs on the **post-shed** `(n, m, v)`, so an ion
whose both draws fire applies shed-then-pickup: net `n` unchanged, but *both*
sets of energy bookings applied.

Both channels draw **unconditionally**, even when the rate is structurally zero.
That is what keeps the scalar single-ion functions exact oracles for the
vectorised ensemble forms, and it is why the relaxation stage can set
`λ₀ = 0` and still preserve the byte-identical RNG stream.

### 3.3 The binding fold is applied by the driver, not the channel

A shed drains `E_int` by `D₀(n)`. For the invariant to close, `E_pot` must rise
by the same `D₀(n)`. But `E_pot` is recomputed from scratch by the BAOAB step,
so the fold has to be re-applied *after* it:

```python
new_state = replace(new_state, E_pot_eV=new_state.E_pot_eV + e_bind_pair_eV(
    new_state.n_shell, picture=..., kappa=..., ladder=ladder))
```

`e_bind_pair_eV(n) = −Σ(n)` is a pure function of the post-event occupancy, so
this is an absolute recompute, not an increment — which means it cannot drift.

---

## 4. The drag law

### 4.1 Forms and dispatch

`physics/drag.py` implements six realised forms behind one tag, each exposing
both `F_drag` and `γ` in **closed form**:

| tag | γ(v) | coefficients |
|---|---|---|
| `linear_cubic` | `g·(a + b·v²)` | a [amu/ps], b [amu·ps/Å²] |
| `linear_quadratic` | `g·(a + c·v)` | a, c [amu/Å] |
| `power_law` | `g·C·v^(n−1)` | C, n |
| `pure_linear` | `g·a` | a |
| `capped_cubic` | `g·b·v²` below `v_c`; `g·b·v_c²·(v/v_c)^p_tail` above | b, v_c, p_tail |
| `capped_linear_quadratic` | analogous | a, c, v_c, p_tail |

`threshold` is declared but deliberately unrealised and raises
`NotImplementedError`.

Two implementation details are worth the thesis's attention.

**γ is never computed as `|F|/v`.** The two are analytically equal, but the
division manufactures a `0/0` singularity at rest that the closed forms do not
have. Each form's γ is written out separately. The cost is that the two
functions could drift apart under a one-sided edit — which matters because the
Tier-3 FDT noise amplitude `√(2γk_BT)` and the O-step both ride on the identity.
That risk is contained by putting the only tricky part, the cap's tail factor,
in a single shared helper `_capped_cubic_tail_factor` used by both.

**The cap's tail is guarded against `np.where` evaluating both branches.**
`np.where` is not lazy: the tail expression is evaluated on the whole array
including the in-band entries. At `p_tail = −1` a literal `(v/v_c)**-1` would
hit `0**-1` at rest. The helper therefore evaluates the tail on
`v_t = max(v, v_c)`, which equals `v` in the tail and is a safe discarded
stand-in below it. It also returns `None` when no sample is above the cap, so
an entirely in-band step takes the pure-cubic fast path — which is what makes
`capped_cubic` **byte-identical to the Tier-0 locked law in-band**.

### 4.2 The gate

One formula, one place (`physics/_gates.py`):

```python
def _erf_complement(depth, steepness):
    return 0.5 * (1.0 - erf(depth / steepness))
```

Two callers with genuinely different physical roles share it: the drag gate
`g(depth)` in `drag.py`, and the `ρ_He/ρ_bulk` occupancy gate in
`helium_density.py` that governs pickup and (on the `density_scaled` arm)
cooling. The helper is a neutral leaf precisely so that `helium_density` never
imports the drag law — it is a density quantity, independent of γ.

### 4.3 The coefficient bundle carries its own provenance

`DragCoefficients` is frozen and validates itself on construction: form-tagged
required keys, extraction mass model, extraction method, and the effective
binding energy it was **jointly calibrated with**. One guard is worth quoting
because it encodes a scientific claim, not a type check:

```python
if self.extraction_method == "free_form" and \
        self.effective_binding_energy_I_ion_eV is not None:
    raise ValueError("extraction_method='free_form' forbids an "
                     "effective_binding_energy_I_ion_eV stamp: ...")
```

A free-form coefficient (the linear counterfactual family) was never jointly
calibrated with any binding energy, so carrying the stamp would claim a
validation that does not exist. The type system is being used to prevent a
provenance lie.

---

## 5. BAOAB and the O-step

`physics/baoab.py` is 207 lines, of which the actual new physics is about
fifteen. B and A reuse `leapfrog._kick` and `leapfrog._drift` — the same
primitives velocity-Verlet uses — so there is no duplicated integrator physics.
The conservative acceleration is evaluated twice per step, fresh, with no
cross-step force caching.

**Why not velocity-Verlet.** Verlet assumes `F = F(x)`. Drag is `F(v)`, so the
final half-kick would need the velocity it is computing. Operator splitting
sidesteps this by giving the velocity-dependent part its own exactly-solvable
sub-step.

**The O-step, and its deliberate asymmetry:**

```python
gamma = gamma_fn(speed_in, depth)      # amu/ps
decay = np.exp(-gamma * dt / m)        # in (0, 1]
new_vel = (decay*vx, decay*vy, decay*vz)
dE_dissip = 0.5 * m * speed_in_sq * (1.0 - decay**2)
```

γ is frozen at a **mixed** point, and this is documented in-code as intentional:

- the **velocity** is frozen at the O-step input velocity — the only choice
  available without an implicit solve;
- the **depth** is taken at the *current* post-half-drift position, which in
  B-A-O-A-B is already legitimately updated.

Three properties follow, and they are the reason the ledger works:

1. `dE_dissip` is **exact for the frozen γ** at any `dt` — it is not a
   first-order estimate of the dissipated energy but the closed-form energy
   difference of the analytic damping.
2. Because `γ ≥ 0` for every guard-validated coefficient set, `decay ∈ (0, 1]`
   and the O-step **can never add kinetic energy**. Dissipativity is exact, not
   approximate.
3. The scheme is `O(dt)` on the cubic law. An earlier second-order claim was
   retired; a one-signed over-braking bias remains and is recorded as a known
   limitation.

The Langevin noise term has a marked, empty site inside `_o_step`. At `T_eff = 0`
**no random number is drawn at all**, and `T_eff > 0` raises
`NotImplementedError` rather than silently running an unvalidated path. That is
Tier 3.

---

## 6. The mass operators

`physics/mass_jump.py` supplies four primitives, all pure and stateless: they
take a velocity and a mass and return a velocity, a mass, and a ledger term.
They perform no scheduling, no integration, and never see the drag law.

**Capture** (pickup fires):

```
v⁺ = (m·v⁻ + m_He·u_He)/(m + m_He)
m⁺ = m + m_He
ΔE_mass_transfer = +½·(m·m_He)/(m + m_He)·|v⁻ − u_He|²      > 0
```

**Cold shed** (He leaves at rest in the lab frame):

```
v⁺ = m/(m − m_He)·v⁻       ← a discontinuous speed kick
m⁺ = m − m_He
ΔE_mass_transfer = −½·(m·m_He)/(m − m_He)·|v⁻|²             < 0
```

**Continuous-velocity shed** (He leaves co-moving) — the Tier-1a production path:

```
v⁺ = v⁻                    ← no kick at all
m⁺ = m − n·m_He
ΔE_mass_transfer = +½·n·m_He·|v⁻|²                          > 0
```

Two things to carry into a write-up.

**The reduced-mass coefficient is computed once**, in
`_reduced_mass_defect_coeff`, as an *unsigned magnitude* — the sign is applied
at each call site (negative for a cold shed, positive for a capture). The same
formula `½·(m·m_He)/m_plus` serves both channels because `m_plus` is `m ± m_He`
respectively. Using the exact reduced-mass form rather than the heavy-ion
`½·m_He·v²` limit is what makes the invariant close *exactly*; the naive form
is off by ~3 % at `n = 1`.

**The co-moving shed needs no reduced-mass correction** — the He leaves at
exactly `v`, so the carried kinetic energy is literally `½·m_He·v²`. The
asymmetry between the two shed operators is physical, not an inconsistency.

The cold shed remains in the codebase, selectable by
`cfg.evaporation_shed_convention`, as an explicitly-labelled diagnostic upper
bound. The `(E_int, n)` jump chain reads neither `v` nor `m`, so **the fire
pattern under a fixed seed is identical across the two conventions** — which
makes the pair a clean controlled comparison.

---

## 7. The two mass channels

Both are "stateful-per-event but storage-free" primitives: they take
`(n, v, m)` plus an injected RNG and return the post-event tuple. Neither
persists `E_int`, performs integration, nor closes the invariant — those belong
to the driver.

### 7.1 Pickup (`physics/pickup.py`)

```python
λ_attach = λ₀ · (ρ_He/ρ_bulk) · max(0, 1 − n/n*)**p
P_attach = 1 − exp(−λ_attach · dt)
```

On a fire: capture reset, `n → n+1`, and the S1 split — `+f_ret·D₀(n+1)` into
`E_int`, the remaining `(1−f_ret)·D₀(n+1)` into `E_dissip`. The two injections
are distinct and there is no double count.

The guards are aggressive and each names the failure it prevents. `λ₀ < 0`,
`ρ_ratio < 0`, `p < 0`, and `n* ≤ 0` all raise, because each produces a channel
that is **silently dead or divergent** rather than obviously broken: a negative
rate gives `P_attach < 0`, and uniform draws live in `[0,1)`, so the channel
would simply never fire — indistinguishable from "pickup on" without the
refusal.

Only `pickup_rate_form="density_only"` and `he_capture_velocity="at_rest"` are
built. `sweeping`, `dwell_time` and `thermal` are valid config members that
round-trip through `cfg.json` but raise `NotImplementedError` **lazily at
point-of-use**, never at config load. This is the project's "rule-2 carry"
convention: an enum arm may be declared before it is built, but it must fail
loudly when reached.

### 7.2 Evaporation (`physics/evaporation.py`)

The rate function encodes *all* gating, so a downstream draw can fire on
`rng.random() < P_shed` alone:

```python
k(E_int, n) = ν · (1 − D₀(n)/E_int)**(s−1)     for D₀(n) < E_int < Σ(n),  n ≥ 2
            = ν                                 for n = 1 and E_int > D₀(1)
            = 0                                 otherwise
```

Three branches, three distinct pieces of physics:

- **`n ≥ 2`** — the RRK bracket inside the self-bound band. The upper bound
  `E_int < Σ(n)` is the parameter-free self-bound gate; the lower bound is
  bracket positivity. Bounded to `k ≤ ν`, so there is no gate-open avalanche.
- **`n = 1`** — the degenerate boundary. The band collapses (`Σ(1) = D₀(1)`) and
  `s = 3n−3 = 0`, so the diatomic is modelled as direct dissociation `k = ν`
  under the *inverted* gate `E_int > D₀(1)`. This is flagged in-code as the
  boundary of validity of the statistical picture.
- **`n ≤ 0`** — no rung, `k = 0`.

The numerical guarding is worth reading as a template. No divide-by-zero, no
`0**0`, and no `0**negative` is ever *formed*, not merely masked:

```python
E_safe   = np.where(E_pos, E, 1.0)
base     = np.where(E_pos, np.maximum(0.0, 1.0 - d0/E_safe), 0.0)
base_safe = np.where(base > 0.0, base, 1.0)
k_rrk    = nu * np.where(base > 0.0, base_safe ** (s - 1.0), 0.0)
```

Because `np.where` evaluates both branches, sanitising the *inputs* rather than
the outputs is the only correct pattern here. The same reasoning appears twice
more: the `D₀` lookup is clamped to `n ≥ 1` and the dof lookup to `n ≥ 2`,
because a tabulated ladder raises outside its range where the analytic Form-U
sigmoid merely evaluated harmlessly. Both clamps are recorded as review fixes.

`effective_dof(n)` returns `4` at `n = 2` (the linear three-atom complex) and
`3n − 3` above, and **raises** for `n < 2` rather than returning a degenerate
value. The empirically-preferred constant `s_eff` arrives through
`cfg.evap_rrk_dof`, which overrides the `n ≥ 2` bracket only and is itself
guarded at `s ≥ 1`.

---

## 8. The internal-energy budget

`E_int` is the fifth invariant term and the only genuinely new state variable in
the ion stage. Its arithmetic lives in `physics/internal_energy_budget.py`
(pure, stateless, no reservoir state) and its cooling in
`physics/solvation_cooling.py`.

The index convention is a place where an off-by-one would be invisible, so it is
stated explicitly in both modules: **a pickup forms the bond at rung `n+1` and
releases `D₀(n+1)`; a shed breaks the current rung `n` and consumes `D₀(n)`**;
in both cases `n` is the pre-event occupancy.

**Cooling (K2)** relaxes not `E_int` but the composite GAH25 variable
`E_solv.struct = E_bind(N) + E_int`, using a closed-form exponential step:

```python
E_new = E_inf + (E_solv_struct - E_inf) * exp(-dt * rho_ratio / tau)
```

Three properties matter. It is **exact**, so one large step equals many small
ones and no integrator state is needed. The asymptote `E_∞(N) = −|S(N)|` is
**occupancy-resolved** — a fixed full-shell asymptote would mechanically halt
shedding, so `E_∞ → 0` as `N → 0` keeps the total strip dynamically reachable.
And `rho_ratio` implements the spatial gate as `τ_eff = τ/ρ̂`, with
`rho_ratio = 1.0` **byte-identical** to the ungated form.

The driver's use of this is a small trick worth noting: it cools
`e_inf + E_int` and then subtracts `e_inf` back off. At fixed `n` the asymptote
cancels exactly, leaving `E_int·exp(−dt·ρ̂/τ)` — so one function serves both the
composite-variable physics and the simple `E_int` decay, with no second
implementation.

---

## 9. What closes the invariant, and what checks it

The five-term invariant is

```
E_kin + E_pot + E_dissip + E_mass_transfer + E_int ≈ const
```

Per channel, the closure works like this:

| event | E_int | E_pot | E_dissip | E_mass_transfer | E_kin |
|---|---|---|---|---|---|
| drag (O-step) | — | — | `+ΔE_dissip` | — | `−ΔE_dissip` |
| K2 cooling | `−drain` | — | `+drain` | — | — |
| shed (K1) | `−D₀(n)` | `+D₀(n)` via fold | — | `+½m_He v²` | `−` same |
| pickup (S1) | `+f_ret·D₀` | `−D₀` via fold | `+(1−f_ret)·D₀` | `+`defect | `−` defect |

Two runtime checks defend it.

**The `m ↔ n` consistency assert**, inside `biphasic_step`:

```python
if not np.allclose(m_p_amu, MASS_I_ION_AMU + n_p * MASS_HE_AMU,
                   rtol=0.0, atol=1e-6):
    raise AssertionError("biphasic_step m<->n consistency drift: ...")
```

The mass and the integer occupancy counter are advanced by *different* code
paths — the mass by the reset operators, `n` by the channels — so they can in
principle disagree. Since the resets produce exact `± m_He`, the residual is
float rounding only, and a `1e-6` amu tolerance is a genuine structural catch
rather than a fudge. This is why `n_shell` is carried as **genuine state** under
`biphasic` rather than re-derived from the mass, as it is on the deterministic
paths.

**The checkpoint writer refuses the silent fallback.** `mass_scenario` is a
required keyword, and a `biphasic` caller with `n_shell is None` raises — because
falling back to `rint((m/U − m_I⁺)/m_He)` would mask exactly the channel/reset
drift the genuine state exists to expose.

**What ledger closure does and does not certify.** It proves the shed is the
reduced-mass form and not a relabel — verified by deliberate fault injection, a
label-only velocity change blows the residual up — and that drag work integrates
to `ΔE_dissip` exactly. It does **not** certify that the energetics of shedding
are physically right. The Tier-1a formulation is worth keeping: *closure is
wiring, not physics.*

---

## 10. Three stages, three integrators

The scored ensemble is not the MD output. The pipeline runs three stages with
genuinely different numerics, and choosing differently in each is a deliberate
decision rather than an inconsistency.

**Stage 1 — MD (`simulation/ion.py`), ~30 ps, fixed `dt = 0.01 ps`.**
Drag live, both mass channels live, cooling live. Full BAOAB.

**Stage 2 — relaxation (`simulation/relaxation_stage.py`), to ~8 ns, fixed dt.**
Fixed-dt is *required* here: K2 cooling makes `E_int` decay continuously between
sheds, which an event-driven solver structurally cannot represent. The stage is
explicitly **not new physics** — it calls the delivered `biphasic_step`
verbatim under a relaxation *view* of the config
(`replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)`), so the locked
K2 → gate → RRK → K1 sequence runs unchanged. Pickup is structurally inert at
`λ₀ = 0` while **its RNG draw is still consumed**, preserving the frozen
two-draw stream byte-identically. Translation is drag-off via a zero-γ BAOAB
closure rather than a separate code path.

A property established by audit and worth stating: with `λ₀ = 0` and `γ = 0`,
the mass subsystem is **completely independent of translation**, so the
`coulomb` and `free_flight` arms produce the identical shed sequence under the
same seed. Translation is retained only for ledger and asymptotic-state
fidelity.

**Stage 3 — detection (`simulation/detection_stage.py`), to 8.53 µs,
event-driven.** An exact Gillespie solver: exponential waiting times per shed,
at most ~n events per ion, *any* detection time reachable. **There is no
timestep, hence no discretization bias** — the fixed-dt Bernoulli chain is the
biased approximation of *this* stage, not the reverse.

Its exactness is conditional and the conditions are checked, not assumed. Three
guards (P1–P3) verify at handover that pickup, drag, and cooling are physically
dead *by position* for every ion — not merely that the stage does not call them.
The distinction is made explicitly in-code: *"omission is not proof of
absence."* Under the verified guards, `E_int` is constant between sheds, the
mechanism reduces exactly to a Markov jump chain at the delivered `rrk_rate`,
and the self-bound margin `G = E_int − Σ(n)` is shed-invariant — so a suppressed
or frozen ion is suppressed or frozen *permanently*. That permanence claim is
true **because** the guard verified the cooling channel is dead, which is why
the guard cannot be skipped.

---

## 11. Checkpoints, RNG, and reproducibility

**Schema.** `IonCheckpoint` is at **v8**. The evolution tracks the physics:
v5 the collision-era baseline; v6 renamed `E_mass_attach_defect_eV` →
`E_mass_transfer_eV` (it can now be positive), added `n_shell`, and dropped the
non-decreasing-mass assumption; v7 added `E_int_eV`; v8 added the CE-channel
fields. Load shims upgrade older files in place, and an unknown version fails at
load rather than being coerced. Loading uses `allow_pickle=False` throughout.

The v8 fields are three per-ion (C)-design surfaces: `ce_channel` (the sampled
CE channel code), `ce_E_m_eV` (its kinetic-energy release, which the driver
turns into the per-molecule Coulomb pair scale), and `ce_strip_count` (the
cumulative exit-strip knock counter). The v7→v8 shim is **silent** rather than
warning, because its synthesized sentinels — no channel sampled, no strip run —
are exact for a pre-v8 file, unlike the v6→v7 arm which has to approximate a
missing reservoir with zeros and says so.

**RNG.** Every stochastic entry point takes an **injected**
`numpy.random.Generator`; there is no module-level global and no implicit
seeding. Stage-private streams are derived by fixed key —
`stage_stream_rng(cfg.seed, EXIT_STRIP_ION_STREAM_KEY)` — so that a stage which
is off consumes nothing and the pre-existing streams are untouched. Draw order
is frozen and treated as part of the physics contract; the project's forbidden
list names changing it explicitly.

**Storage stride.** The driver estimates checkpoint size from a hard-coded array
count (`_NUM_2N_T_ARRAYS_ION = 15`) and downsamples *storage* when a full-
resolution checkpoint would exceed 1 GB. Internal steps still run at `dt_ion`.
The array count is guarded by a schema-count test, because an undercount would
silently grant stride 1 past the byte budget.

---

## 12. Challenges, and what a reviewer would attack

Stated plainly, because they are real.

**The mass↔coefficient pairing guard is run under an explicit exception.** The
drag law was extracted at a constant effective mass `m_eff`, and a config-load
guard enforces that constant-mass coefficients pair only with
`mass_scenario=fixed`. Both the Tier-1a and Tier-2 production runs trip this
guard *structurally* and run under `allow_inconsistent_mass_pairing=True`. The
defence is a mid-window argument — `m ≈ 19 He ≈ m_eff` mid-window, with the
`n = 21` and `n = 14` ends in the free-extrapolation zone. It is a defensible
argument, not a proof, and it is the single largest standing assumption in the
integrator wiring.

**`f_ret` has never been swept.** The S1 retained fraction sits at 0.1 in every
production run. It is a genuine free parameter and is documented as a gap.

**`s_eff` is empirical.** The classical mode count `3n − 3` is falsified in-model
(the cascade freezes at `n ≈ 20`); a constant `s_eff ≈ 8` reproduces the
anchored staircase. The literature justification for an effective degrees-of-
freedom count an order of magnitude below the mode count is an open question
(RQ6).

**The twin forward model cannot represent the mass channels.** Most of the
parameter scanning was done with a fast twin that fixes `n` at birth, integrates
at constant mass, and has **no pickup channel at all**. Every twin number in the
mass thread therefore omits both the mass growth and the ~3.6 amu/ps of capture
friction — an asymmetry outside the measured twin↔MD authority box, and one to
state explicitly wherever twin results are quoted.

**The pickup channel's own drag is never booked.** `γ_pickup = λ₀·ρ̂·m_He`
evaluates to about 3.6 amu/ps at the production pin — of order 5–10 % of the
drag law — and does not appear in the ledger. Open as M8.

**The `O(dt)` over-braking bias is one-signed.** The O-step freezes γ at the
input velocity, and since γ increases with `v` on the cubic, the frozen value
over-estimates the friction the ion actually experiences as it decelerates. The
bias therefore does not average out.

**The three-stage timescale separation was calibrated at the wrong geometry.**
It assumed ejection is effectively instantaneous, which holds at `R ≈ 27 Å`. At
the corrected `R ≥ 49 Å` escape takes ~0.1–1 µs with helium still present, and a
third to a half of ions never leave. This is logged as scaffolding for a
geometry that has since been corrected — the first ledger row the correction
created rather than retired.

**And the largest one, which is not a code defect at all:** the observable is
scored at bare mass while the trajectory is integrated dressed, so every
pre-evaporation energy is attenuated by `m(1)/m(21) = 0.6205`. The bookkeeping
is self-consistent; whether 21 helium atoms can co-move through a 2.7 eV
Coulomb explosion is the open physics question
(`TIER2_MASS_SCENARIOS.md`, gated by M1).

---

## 13. Where to look

| to see | read |
|---|---|
| the whole per-step sequence | `simulation/ion_propagation_step.py` — `biphasic_step`, then `baoab_propagation_step` |
| dispatch, the loop, closure rebuilds | `simulation/ion.py::run_ion_propagation` |
| the drag laws and the tail guard | `physics/drag.py` — `drag_gamma`, `_capped_cubic_tail_factor` |
| the only place mass enters the drag | `physics/baoab.py::_o_step` |
| the shed/capture arithmetic | `physics/mass_jump.py::_reduced_mass_defect_coeff` |
| all the gating logic in one function | `physics/evaporation.py::rrk_rate` |
| the invariant's structural guard | the `m ↔ n` assert in `biphasic_step` |
| the exactness conditions of the µs stage | `simulation/detection_stage.py` module docstring, P1–P3 |
| the mechanism physics (not the code) | `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` |
| the architecture decisions (not the code) | `DRAG_PORT_DESIGN_DECISIONS.md` |
| the narrative and the results | `PROJECT_OVERVIEW.md` |
