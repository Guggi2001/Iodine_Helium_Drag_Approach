# Slice 2 — Goals: BAOAB Ion-Stage Stepper

**Status:** Specification. Physics definition + integrator scheme + interface
contract only. No implementation code until `[PROCEED TO IMPLEMENTATION]`.

**Scope:** Slice 2 *only* — the ion-stage time-stepper that replaces the
baseline `velocity_verlet_step` for the ion stage with a BAOAB operator-split
integrator (decision §4.6). It consumes Slice 1's `drag_gamma` as its only new
physics input. The neutral stage is untouched. No `SimConfig` wiring (Slice 3),
no driver rewiring (Slice 4), no mass dynamics, no active noise.

**Why this is next (not Slice 3).** Same discipline as Slice 1: build the leaf
with closed-form killer tests in isolation before any integration.

- Slice 2 needs nothing from Slice 3 — it takes its parameters as plain
  arguments (a conservative `acc_fn`, a `gamma_fn`, `m`, `dt`, `droplet_radii`,
  and a dormant `T_eff`). Config wiring is later.
- Its failure modes are integrator-correctness against analytic solutions —
  orthogonal to Slice 4's integration failures (real run, energy closes,
  checkpoint round-trips). Folding them would couple integrator bugs with
  wiring bugs.

Order stays **2 → 3 → 4**.

---

## 1. The Tier-0 reduction (the key simplification)

At Tier 0 the entire stochastic apparatus is dormant and the stepper is
deterministic, yet it is still the *production* Tier-0 integrator:

- **mass fixed** at $m_\text{eff}$ (`mass_scenario = fixed`),
- **noise off** ($T_\text{eff} = 0$),
- **drag** = `linear_cubic` (Slice 1, the only realised form).

Under these, the O-step collapses from the full Ornstein–Uhlenbeck update to
**pure multiplicative velocity damping**:

$$v \;\mapsto\; e^{-\gamma(v,\text{depth})\,dt/m}\;v.$$

So Slice 2 is built and tested **fully deterministically**. Noise is added
later by un-pinning $T_\text{eff}$ and adding one RNG draw inside the O-step —
a localised change to that step alone, not a restructuring.

---

## 2. Friction convention (carried over; mass enters *here*)

The unified convention from Slice 1 / the decisions doc holds:

- $\gamma(v)$ is a **force coefficient**, units **amu/ps**,
  $\gamma(v) = |F_\text{drag}(v)|/v$. Slice 1 exposes it in closed form.
- The friction **rate** is $\gamma/m$ [1/ps] and appears **only** in the BAOAB
  damping exponent $e^{-\gamma\,dt/m}$ — which lives in *this* module.

**This is the one place mass enters.** Slice 1 is mass-free; the stepper
divides by $m$ exactly once, in the O-step exponent. Kinetic-energy
bookkeeping (the dissipated-energy return, §6) also uses $m$. Inertia in the
B-step uses $m$ via the conservative `acc_fn` (which already divides force by
mass in the baseline).

---

## 3. The BAOAB scheme

Symmetric Strang composition **B–A–O–A–B** over one step $dt$, ion stage only.
Let `acc_fn(x) → (a_cons, E_pot)` be the *existing* conservative ion
acceleration (Coulomb + droplet, `leapfrog._ion_accel_fn`), unchanged.

```
state (x0, v0)
  B   v ← v0 + (dt/2)·a_cons(x0)          # half-kick, conservative force only
  A   x ← x0 + (dt/2)·v                   # half-drift
  O   v ← e^(−γ·dt/m)·v  (+ noise, off)   # full OU: drag (+ dormant noise)
  A   x ← x  + (dt/2)·v                   # half-drift
  B   v ← v  + (dt/2)·a_cons(x)           # half-kick, conservative force only
return (x1, v1, E_pot, ΔE_dissip)
```

- **B (kick)** and **A (drift)** are exactly velocity-Verlet's kick and drift —
  *the* fact the anchor test (§7) exploits.
- **O** is the only genuinely new physics: drag as damping, plus the dormant
  noise site. Conservative forces never enter O; drag never enters B.
- Two `acc_fn` evaluations per step — at `x0` (first B) and at `x1` (last B) —
  matching the two evaluations `velocity_verlet_step` makes. The canonical
  efficient form caches the final-B force as the next step's first-B force;
  caching vs. re-evaluating must not change the result (guarded by the anchor
  test), so the choice is left to implementation.

### Governing O-step (Tier-0, noise off)

$$v_\text{out}^{O} = e^{-\gamma(v_\text{in}^{O},\,\text{depth})\,dt/m}\;
v_\text{in}^{O},\qquad
\Delta E_\text{dissip} = \tfrac12\,m\,\big(\lVert v_\text{in}^{O}\rVert^2 -
\lVert v_\text{out}^{O}\rVert^2\big)\;\ge 0.$$

### Governing O-step (full form, for reference — noise dormant at Tier 0)

$$v_\text{out}^{O} = e^{-\gamma\,dt/m}\,v_\text{in}^{O}
+ \sqrt{\frac{k_B T_\text{eff}}{m}\Big(1 - e^{-2\gamma\,dt/m}\Big)}\;\boldsymbol\xi,
\qquad \boldsymbol\xi \sim \mathcal N(0,\mathbb 1).$$

$\boldsymbol\xi$ is a **dimensionless standard-normal** draw — *not* the
continuous-time $\xi(t)$ (units ps$^{-1/2}$) of decisions-doc §1.2; they are
different objects. At $T_\text{eff}=0$ the second term vanishes and this
reduces to the Tier-0 form above. The gate $g(\text{depth})$ lives inside
$\gamma$ (Slice 1: $\gamma = g\,(a+bv^2)$), so it scales **both** the damping
and the noise amplitude — the FDT coupling of §5.2, automatic and free.

---

## 4. The two freeze decisions (as discussed — stated so they are not "fixed"
later as bugs)

The nonlinear O-step is frozen to an $O(dt)$ approximation (§4.4). Its two
arguments are frozen **asymmetrically and deliberately**:

- **Velocity → frozen at the O-step input velocity** $v_\text{in}^{O}$ (the
  post-first-B velocity; A does not change velocity). Compute
  $\gamma(v_\text{in}^{O},\cdot)$ once, apply the exponential. This is the only
  choice available without an implicit solve, and it makes the dissipated
  energy *exact for the frozen $\gamma$*: the energy removed is computed with
  the same $\gamma$ that damped the velocity.
- **Depth (gate) → evaluated at the current O-step position** (after the first
  half-drift A). In B–A–O–A–B the position is *already* updated when O runs, so
  `depth = r_atom − r_droplet` is freshly and legitimately known — no reason to
  stale it to step entry.

**Consequence (intentional):** $\gamma$ is evaluated at a *mixed*
$(v_\text{old}, \text{depth}_\text{new})$ point. This is correct — they are
independent arguments to a pure function — and must be stated in the module
docstring so it is not "corrected" later.

**Why it matters beyond reproducibility.** Slice 1 guarantees
$\gamma(v,\text{depth}) \ge 0$ for `linear_cubic` with $a,b>0$, so
$e^{-\gamma dt/m} \in (0,1]$ and $\lVert v_\text{out}^{O}\rVert \le
\lVert v_\text{in}^{O}\rVert$ — **the O-step can never add kinetic energy,
exactly.** Freezing depth at step-entry instead would let a fast
surface-crossing gate with a stale wrong-side $g$ and forfeit that exactness
for no gain.

---

## 5. Dimensional analysis (gate before adoption)

| expression | units | balance |
|---|---|---|
| damping exponent $\gamma\,dt/m$ | $(\text{amu/ps})(\text{ps})/\text{amu}$ | $=1$, dimensionless ✓ |
| damped velocity $e^{-\gamma dt/m}\,v$ | (dimensionless) × Å/ps | Å/ps ✓ |
| half-kick $(dt/2)\,a_\text{cons}$ | (ps)(Å/ps²) | Å/ps ✓ |
| half-drift $(dt/2)\,v$ | (ps)(Å/ps) | Å ✓ |
| $\Delta E_\text{dissip} = \tfrac12 m(v_\text{in}^2 - v_\text{out}^2)$ | (amu)(Å/ps)² | amu·Å²/ps², a mechanical energy; $\ge0$ ✓ |
| noise term $\sqrt{(k_BT_\text{eff}/m)(1-e^{-2\gamma dt/m})}\,\xi$ | $\sqrt{\text{Å}^2/\text{ps}^2}$ × (dimensionless) | Å/ps ✓ |

Noise-amplitude inner check: $k_B T_\text{eff}$ in amu·Å²/ps² (an energy),
$\div m$ in amu → Å²/ps², root → Å/ps. Balances. **No formulation rejected.**

---

## 6. Energy bookkeeping at the interface

- **Dissipated energy is returned in mechanical units `amu·Å²/ps²`**, computed
  analytically from the O-step velocity damping
  $\tfrac12 m(\lVert v_\text{in}^{O}\rVert^2 - \lVert v_\text{out}^{O}\rVert^2)$,
  per atom, $\ge 0$. **FLAG FOR SLICE 4:** Slice 4 converts to eV for the
  checkpoint `E_dissip_eV` accumulator (the baseline conversion uses mass in kg
  and $v\times100$ for m/s, then $\div\,$`EV`; see
  `ion_propagation_step.py` defect term). The stepper stays pure-mechanical and
  does **not** import the eV conversion.
- **Return it now, even though nothing consumes it until Slice 4.** The O-step
  produces $\Delta E_\text{dissip}$ as a byproduct; returning it from the step
  signature now keeps the interface stable. **FLAG: consumed by Slice 4** (the
  §2.9 energy invariant
  $E_\text{kin}+E_\text{pot}+E_\text{dissip}+E_\text{mass\_transfer}$); inert
  in Slice 2 beyond being asserted correct by the §7 tests.
- When noise turns on (later), the noise *injects* energy that feeds the
  thermal floor and must be tracked **separately** from $E_\text{dissip}$ (a
  sibling channel, §4.5 of the decisions doc). Out of Slice 2 scope; the return
  signature should leave room for a second energy term without breaking.

---

## 7. Placement and the `leapfrog.py` helper extraction

- **New module `physics/baoab.py`.** Rationale: the baseline already houses its
  integrator (`leapfrog.py`) under `physics/`, so an integrator belongs there
  by precedent; and §4.6's deliberate neutral/ion fork reads more clearly as a
  separate file than as a second factory inside the file whose docstring says
  it is model-agnostic.
- **Kick/drift helpers extracted from `leapfrog.py`.** BAOAB's B and A *are*
  velocity-Verlet's kick and drift. To avoid duplicate physics
  (`CLAUDE.md` rule 1), extract `_kick` / `_drift` helpers from the monolithic
  `velocity_verlet_step` and have both it and `baoab.py` call them.
- **This is the one knowing touch to the frozen baseline integrator.** It
  crosses the baseline's "the drag port should not need to touch `leapfrog.py`"
  note (PHYSICS_BASELINE §5) — accepted because §4.6 (the newer, more specific
  decision) added an integrator at all, and because the extraction is a pure
  refactor that the anchor test (below) makes **self-verifying**: if extracting
  kick/drift changed `velocity_verlet_step`'s behaviour at all, the anchor test
  fails loudly. Record the touch explicitly in the commit and the module
  docstring.

---

## 8. Module interface (sketch — contract, not code)

- **`make_ion_baoab_step(m, droplet_radii, acc_fn, gamma_fn, dt, *,
  T_eff=0.0, rng=None)`** → closure `step(pos, vel) → (pos', vel', E_pot,
  ΔE_dissip)`. Mirrors `leapfrog.make_ion_step`. `gamma_fn(v, depth)` is the
  Slice 1 `drag_gamma` closed over `coeffs` + `steepness`; `acc_fn` is the
  conservative ion acceleration. `T_eff=0` and `rng=None` keep noise dormant.
- The closure computes `depth = r_atom − droplet_radii` internally (mirroring
  the baseline `r1 − droplet_radii`) and passes it to `gamma_fn`, so the gate
  is evaluated at the current O-step position per §4.
- **Per-step closure rebuild (confirmed).** Slice 4 will rebuild this closure
  every step (mass changes under future mass scenarios), exactly as
  `ion_propagation_step.py:184-193` rebuilds `make_ion_step`. The factory is
  built to support that even though Tier-0 mass is fixed — so mass dynamics drop
  in later with no restructuring.
- Returns dissipated energy in **amu·Å²/ps²** (§6).

---

## 9. Noise — dormant at Tier 0 (recorded so the off-switch is unambiguous)

- Noise is **off** when `T_eff = 0` (the Tier-0 default); the second O-step term
  is identically zero and no RNG is drawn.
- When activated (Slice ≥3): `T_eff > 0`, one **dimensionless standard-normal**
  draw per atom per active spatial channel (T1 longitudinal-only first, §1.4),
  amplitude $\sqrt{(k_BT_\text{eff}/m)(1-e^{-2\gamma dt/m})}$ with the gate
  already inside $\gamma$. The **RNG draw-order pinning** (to preserve
  determinism / match any reference stream) is a Slice ≥3 concern, not settled
  here.
- Slice 2 only ensures the O-step is *structured* so the noise term slots in
  at one site without touching B, A, or the energy-return shape.

---

## 10. Acceptance criteria (test surface)

| Quantity | Expected | Tolerance |
|---|---|---|
| **Anchor test** — $\gamma=0$, noise off: `baoab_ion_step` ≡ baseline `velocity_verlet_step` on identical state + `acc_fn` | exact agreement (turning drag off recovers the frozen baseline) | round-off (rtol ~1e-12) |
| Analytic decay — constant linear $\gamma$, no conservative force: $v(t)=v_0\,e^{-\gamma t/m}$ over many steps | matches closed form | tight (rtol ~1e-6, $dt$-limited) |
| Dissipated-energy identity — returned $\Delta E_\text{dissip}$ equals $\tfrac12 m(\lVert v_\text{in}^O\rVert^2-\lVert v_\text{out}^O\rVert^2)$ | exact by construction | round-off |
| Dissipativity — O-step never increases $\lVert v\rVert$ for $\gamma\ge0$ | $\lVert v_\text{out}^O\rVert \le \lVert v_\text{in}^O\rVert$ | exact |
| Gate-off limit — deep vacuum ($g\to0\Rightarrow\gamma\to0$): O-step is identity, $\Delta E_\text{dissip}=0$ | no drag, no dissipation outside droplet | round-off |
| Helper-extraction safety — `velocity_verlet_step` output unchanged vs. pre-refactor (covered by the anchor test + existing `test_*` for the neutral/ion Verlet path) | byte-identical behaviour | round-off |
| Energy-return units — $\Delta E_\text{dissip}$ in amu·Å²/ps², not eV | enforced (no eV conversion in module) | exact |
| Noise dormant — `T_eff=0` path draws no RNG and equals the deterministic O-step | identical with/without `rng` passed | exact |
| Mass enters only in O-step + energy — `gamma_fn` and `acc_fn` carry no extra mass coupling | enforced | exact |

The anchor test is the killer test: it proves the swap is safe by recovering
the frozen baseline in the no-drag limit.

---

## 11. Open items surfaced (recorded, not blockers)

- **Force-evaluation caching** (cache final-B force as next first-B vs.
  re-evaluate). Performance/cleanliness choice; the anchor test guards
  correctness either way. Decide at implementation.
- **Second energy channel for noise.** The return signature should anticipate a
  separate noise-injection energy term (§6) so adding noise later does not
  reshape the interface. Shape it now or accept a signature bump at Slice ≥3 —
  flag, decide at implementation.
- **RNG draw-order** for the future noise draw (§9) — deferred to Slice ≥3.

---

## 12. Scope fence — what Slice 2 does NOT touch

- **No `SimConfig` fields.** Parameters arrive as function arguments; the enum
  surfaces and the §6.5 guard are Slice 3.
- **No driver wiring.** Replacing the collision call sites in
  `ion_propagation_step.py:211-256` and rebuilding the closure per step in the
  *driver* is Slice 4. Slice 2 provides the factory; Slice 4 calls it.
- **No checkpoint changes / energy rename.** `IonCheckpoint` v6 and
  `E_mass_attach_defect_eV → E_mass_transfer_eV` come with mass dynamics.
- **No active noise.** Dormant at `T_eff=0`; structured for later (§9).
- **No mass dynamics.** Mass fixed at Tier 0; the per-step rebuild pattern is in
  place so dynamics drop in later.
- **No eV conversion.** Stepper stays mechanical (amu·Å²/ps²); Slice 4 converts.
- **Neutral stage untouched.** Keeps `velocity_verlet_step` (§4.6 accepted
  asymmetry). The only `leapfrog.py` change is the kick/drift helper extraction,
  which leaves the neutral path behaviourally identical (§7).

---

## 13. Definition of done

- `physics/baoab.py` with the B–A–O–A–B stepper and `make_ion_baoab_step`
  factory, consuming Slice 1's `drag_gamma` and the baseline conservative
  `acc_fn`.
- Kick/drift helpers extracted from `leapfrog.py`; `velocity_verlet_step`
  behaviourally unchanged (anchor test + existing Verlet tests green).
- Asymmetric freeze (velocity at O-input, depth at O-position) implemented and
  documented in the module docstring as intentional.
- $\Delta E_\text{dissip}$ returned in amu·Å²/ps², flagged for Slice 4 eV
  conversion.
- Noise dormant at `T_eff=0`, O-step structured for later activation.
- All §10 acceptance criteria green, anchor test first.
- The two §11 implementation-time open items recorded in-repo.
- No `SimConfig`, driver, checkpoint, or neutral-stage changes.
