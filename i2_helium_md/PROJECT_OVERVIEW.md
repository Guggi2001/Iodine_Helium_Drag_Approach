# Project Overview — Coulomb Explosion of Iodine in Helium Nanodroplets

> **What this document is.** The single reader-facing account of the whole
> project: what is being simulated, what the predecessor model did, why it was
> replaced, how the replacement was calibrated, and where it stands today.
> Written as thesis background, so every number carries an inline pointer to
> the document that owns it.
>
> **What it is not.** It is not a status file. Status, decisions and history
> live in the detail docs (§15). Where this document and
> `drag_migration_log_tier2.md` disagree, the log wins.
>
> **Standing status at time of writing (2026-08-13).** Production point
> `finc1v725` stands. Successor candidate `h405` is measured but **not
> adopted** — the G4 adjudication is open. Tier-0's in-band form verdict
> stands. The free-form linear counterfactual family is **closed**. The
> current open arm is the flight mass (`TIER2_MASS_SCENARIOS.md`), gated by a
> literature question (M1).

---

## 1. What is being simulated

An iodine molecule I₂ is embedded in a superfluid helium nanodroplet of a few
thousand to a few tens of thousands of He atoms. A laser pulse doubly ionizes
it. The two resulting I⁺ ions repel each other under the bare Coulomb force —
a **Coulomb explosion** — and fly apart through the surrounding helium, each
one dressed in a solvation shell of roughly 21 He atoms. As they travel they
dissipate energy into the droplet, shed helium, and eventually escape into
vacuum still carrying a finite shell. A velocity-map imaging (VMI) detector
records them **8.53 µs** later, resolving both their kinetic energy and, in
the companion mass spectrum, how many helium atoms each ion kept.

The simulation reproduces this experiment. Each iodine atom feels three
distinct forces simultaneously (`docs/physics/physics_background.md` §1):

- **droplet solvation** — a smooth erf profile
  `V(r) = ½(erf((r − offset)/steepness) + 1)·E_bind`, because the He surface
  has finite thickness rather than being a step;
- **the partner atom** — a Morse bond in the neutral stage, bare Coulomb
  repulsion in the ion stage;
- **the helium bath** — the dissipation channel, which is the subject of this
  whole project.

The run is two-staged (`PHYSICS_BASELINE.md` §1): a **neutral stage**
thermalizes the embedded molecule, and at `t_neutral_end` the charges and
potentials switch to the **ion stage**, which propagates the Coulomb-driven
dissociation. Internal units throughout are Å, ps, Å/ps, eV and amu.

The two questions the experiment poses, and the model must answer, are:

1. **How much kinetic energy does a fragment carry when it arrives?**
2. **How many helium atoms is it still wearing?**

Everything below is about getting both right at once.

---

## 2. The predecessor: a hard-sphere collision model

The project began as a Python port of a legacy MATLAB molecular-dynamics code.
Its dissipation model is documented in `PHYSICS_BASELINE.md` §6, a section
titled **"THE BASELINE BEING REPLACED"**.

Energy loss to the helium was modelled as a **stochastic sequence of discrete
elastic binary hard-sphere collisions** with individual He atoms — a dilute-gas
kinetic picture. The cross-section was velocity-dependent,

```
σ(v) = σ₀ · v^exponent   [Å²],   σ₀ = 2500 Å²,   exponent = −2.0
```

with collisions sampled by a Bernoulli gate `p_scatter = d·σ·ρ_droplet` over
the previous step's path length `d`, subject to a Landau cutoff at
`v_limit = 40 m/s` (which is what keeps σ → ∞ as v → 0 from mattering). Each
event drew an impact parameter by inverse CDF (`b/R = √u`), applied the
hard-sphere centre-of-mass scattering relations, and re-assembled the velocity
on a random orthonormal basis.

**The whole-shell-stripped-at-start assumption.** The ion stage inherits the
*bare* iodine mass from the neutral checkpoint — 127 amu
(`ion_initial_state.py`). Physically this means the entire ~21-atom solvation
shell is assumed to be removed at the instant of ionization, after which the
ion **re-accretes** helium one atom at a time through a second Bernoulli gate
riding on the collision gate:

```
b_mass_attach = (uniform < mass_attach_probability) & b_collision
mass_attach_probability = 0.09 (9 Å preset) / 0.005 (18 Å preset)
```

Mass history was strictly non-decreasing, and because the He attached at the
ion's own velocity the recomputed `E_kin` gained a fictitious `½Δm·v²`, which
the code subtracted back through a cumulative accumulator so that
`E_kin + E_pot + E_dissip + E_mass_attach_defect ≈ const` closed.

**It worked.** This model kept the whole validation surface green
(`PHYSICS_BASELINE.md` §6.4, §13): HeDFT trajectory comparison inside
tolerance at both R₀ = 9 Å and 18 Å, VMI final-velocity histogram agreement,
and a closed energy ledger. The predecessor's results are the reason the
project continued, and any account that treats the collision model as simply
wrong is misreading the history.

---

## 3. Why it was replaced

Four reasons, in increasing order of how much they matter.

**1. It is the wrong physics for a superfluid.** Discrete binary scattering off
individual atoms is a dilute-gas picture. The ion actually travels inside a
*bubble* in a superfluid, where dissipation is continuous.

**2. Parameter laundering.** The model reproduced the references by co-fitting
three knobs — the ion binding energy, the geometric cross-section, and the
velocity exponent (`PHYSICS_BASELINE.md` §3.5.1). The evidence that these were
fit parameters and not constants is direct: the 18 Å preset drops the *same
species'* I⁺ binding from **0.30 eV to 0.05 eV**, a factor of six, and the
attach probability drops 0.09 → 0.005 with droplet size. A cross-section of
2500 Å² is not a measured quantity.

**3. The escape was brute force.** `TIER0_FINDINGS.md` Finding 3 puts it
bluntly: *"The old hard-sphere model only escaped because its collisions
over-accelerated the ions past the barrier."*

**4. Discrete events cannot generate the observable.** The production target is
the per-fragment I⁺Heₙ distribution. A fixed-mass run yields exactly one
species and cannot generate fragment channels
(`DRAG_PORT_DESIGN_DECISIONS.md` §2.1). The Bernoulli gate also fused three
physically distinct behaviours — the energy/angle update, sticky-vs-elastic
character, and helium pickup — into one coin flip.

**The newer literature that made a replacement possible.** Three results
underpin the new model:

- **Albrechtsen, Stapelfeldt et al., *Nature* 623, 319 (2023)** — "Ion
  solvation in helium droplets". Na⁺Heₙ Coulomb-ejected, dissociating in
  flight, VMI-detected. This is the model's foundational paper: it supplies the
  D₀(N) ladder picture, Poisson pickup, the "shed one He when E_int > D₀(N)"
  rule, and Newton-type dissipation.
- **Time-resolved alkali TDDFT, *J. Chem. Phys.* 160, 164308 (2024)** —
  dissipation follows **Newton's law of cooling** for the first ~5 ps, roughly
  ion-independent across the alkali series. The direct precedent for the
  cooling term.
- **González group, *PCCP* 2016 (10.1039/C6CP04315A)** — TDDFT of halogen-dimer
  photodissociation in helium; a two-regime mechanism ending in **viscous-flow
  friction** of "probably general character". The closest microphysical
  justification for a drag term.

**And one standing caveat, carried on every page of the design docs**
(`DRAG_PORT_DESIGN_DECISIONS.md` §3.3):

> **No literature law exists for the velocity dependence of drag on an ion
> travelling through superfluid helium**, so every analytic form here is a
> *hypothesis* to be cross-checked against the TDDFT curves, not a fitted
> truth.

Together with its companion — **TDDFT is not ground truth; experiment
arbitrates** — this defines the epistemic stance of the entire program.

---

## 4. The drag model

**Friction convention (unified, and deliberately unusual).** `γ(v)` is a
**force coefficient** in **amu/ps**, defined `γ(v) = |F_drag(v)|/v`. The
friction force is `γ(v)·v` — with **no leading mass**. The friction *rate*
`γ/m` [1/ps] appears only inside the BAOAB damping exponent `e^(−γ·dt/m)`.

The consequence is structural: the drag module is **mass-agnostic**. It never
takes a mass. Mass enters at exactly one place, the integrator's O-step, as one
explicit division by `m(t)`. This matters later, when the mass itself becomes
the open question (§14).

**Spatial gating.** Drag acts only inside the droplet, through an
erf-complement gate on the signed depth `d = r_atom − r_droplet`:

```
g(d) = ½ · (1 − erf(d / steepness)),    g ∈ [0, 1], dimensionless
```

**The candidate forms.** Because no literature law exists, several hypotheses
were built behind a config enum, chosen deliberately to disagree most in
extrapolation:

| form | F_drag (magnitude) | γ(v) | coefficient units |
|---|---|---|---|
| `linear_cubic` | `g·(a·v + b·v³)` | `g·(a + b·v²)` | a [amu/ps], b [amu·ps/Å²] |
| **pure cubic** (a ≡ 0) | `g·b·v³` | `g·b·v²` | b [amu·ps/Å²] |
| `linear_quadratic` | `g·(a·v + c·v²)` | `g·(a + c·\|v\|)` | c [amu/Å] |
| `power_law` | `g·C·\|v\|ⁿ` | `g·C·\|v\|^(n−1)` | C [amu·Å^(1−n)·ps^(n−2)] |
| `capped_cubic` | `g·b·v³` for v ≤ v_c; `g·b·v_c²·v·(v/v_c)^p_tail` above | — | v_c [Å/ps], p_tail dimensionless |
| `pure_linear` | `g·a·v` | `g·a` | a [amu/ps] |

Dimensional check, one worked example: `b·v³ = (amu·ps·Å⁻²)(Å/ps)³ =
amu·Å/ps²` — a force in the internal unit system. ✓

Two **exact nesting identities** are unit-tested: `power_law(n=2, C=c)` is the
pure quadratic and `power_law(n=3, C=b)` is the pure cubic. The free-exponent
fit therefore nests both incumbents, which makes n̂ a *direct measurement* of
exponent identifiability rather than a model comparison.

**Integrator.** Velocity-Verlet breaks when the force depends on velocity, so
the ion stage uses **BAOAB** operator splitting
(`DRAG_PORT_DESIGN_DECISIONS.md` §4), with γ frozen at the step-entry velocity
inside the O-step. That gives exact dissipation bookkeeping at any timestep,
`ΔE_dissip = ½m(‖v_in‖² − ‖v_out‖²)`, and unconditional dissipativity.

---

## 5. Tier 0 — calibrating the drag against TDDFT trajectories

**Why validation is staged.** The parameters are entangled: drag strength,
binding energy, mass dynamics and noise all push the same trajectory around.
Fitting them simultaneously fits nothing. The program therefore validates in
**tier order**, freezing each tier's winner before introducing the next
unknown: Tier 0 the drag *form* (deterministic, fixed mass, TDDFT traces);
Tier 1a the *mass kinematics*; Tier 2 the *size distribution* against
experiment; Tier 3 the *second moments* (noise).

**Method B — trajectory matching.** The original extraction (Method A)
computed `F_drag(t) = m_eff·a(t) − F_C(t)` point by point and regressed the
result against speed. It never saw a forward-integrated trajectory. Method B
instead chooses the coefficients that minimize the forward-integrated,
radial-projected, in-window trajectory RMSE, using the **actual BAOAB ion
driver** (`METHOD_B_trajectory_matching_extraction.md`):

```
{a,b}★ = argmin RMSE_[t*, t_end]( |v₂|^MD(a,b) , |v₂|^smoothed-ref )
```

The fit objective *is* the trajectory match, so it absorbs the whole pipeline —
O-step, spatial gate, mass treatment — and optimizes the quantity that matters.
Its structural cost is that the trajectory match can no longer also serve as
the validation (that would be training on the test set), so Tier 0's
consistency-check role was retired and replaced by explicit **held-out**
validation: a held-out *case* (fit 18 Å, predict 9 Å) and a held-out
*observable* (the experimental distributions, which is Tier 2).

**The two references and the non-radial finding.** The traces are TD-He-DFT
Coulomb-explosion trajectories of I₂²⁺ in droplets of radius 9 Å and 18 Å
(`data/reference/9A_All_Data.csv`, `18A_All_Data.csv`). They differ physically,
and that difference is the entire leverage:

| case | window | in-window speed range | character |
|---|---|---|---|
| 18 Å | [4.54, 8.0] ps | 2.54–3.02 Å/ps | clean radial throughout |
| 9 Å | [2.67, 14] ps | 2.83–4.95 Å/ps | genuinely non-radial |

At the window start, the 9 Å atom is **~99.6 % transverse** (|v| 4.90, radial
−0.29, transverse 4.89) while the 18 Å atom is ~99.9 % radial. The honest
reading, recorded as a first-class finding, is that the surviving ~0.39 Å/ps
residual at 9 Å is a **model-dimensionality statement**: a central-force MD
(Coulomb + radial droplet + radial drag) structurally cannot carry a sustained
transverse co-translation. Three earlier readings of that residual — a frame
systematic, a windowing/bubble-mode effect, and a ~10 % over-damping — were
adopted and then **withdrawn** as better data arrived. The withdrawal history
is preserved deliberately in `drag_migration_log_tier0.md`.

**The verdict.** The shared-form joint refit fits *one* size-independent law
plus *one* effective binding across *both* cases — three parameters against two
full trajectories, over-constrained and therefore falsifiable. Results:

- **`shared_pure_cubic` wins:** `γ = g·b·v²` with **b = 2.5154 amu·ps/Å²** and
  **E_bind = 0.1168 eV**. 18 Å RMSE 0.1345, 9 Å 0.1240, escape 1.0/1.0.
- **The linear term carries no information** — `a` ran to the zero bound
  (0.00025 amu/ps) and the equivalence test returned T_a0 = 0.000000.
- **The free exponent independently recovers cubic:** all four `power_law`
  starts converged to **n̂ = 2.927 ± 0.279**, objective within 0.0001 of the
  incumbent. This resolved a standing 3-vs-2 tension in favour of 3 and
  retired the old Method-A n ≈ 2.06 as an artifact of the hard-sphere
  σ ∝ v⁻² cross-section bleeding into the extraction.
- **Forced quadratic is rejected:** `linear_quadratic` collapses to the pure
  quadratic corner (a → 0, c ≈ 12.8), objective worse by **Δ +0.0339**, and its
  held-out 9 Å prediction fails at **0.699 > 0.45**. Pure cubic's held-out
  prediction was 0.2685.

**The velocity window that owns all of this: 2.54–4.95 Å/ps.** Outside it the
traces say nothing. That single fact drives §12 and §13, so it is worth
stating early and plainly.

**Finding 3, and the trap it names.** With the in-window-correct drag, ions
arrive at the droplet boundary with radial kinetic energy *below* the static
solvation barrier, reverse, and never escape. This is not a bug in the drag: the
predecessor's own TD-HeDFT analysis notes that real ions escape despite having
less than the solvation energy, because the helium *reorganizes* — "ion ejection
cannot be predicted with static solvation potential values alone". The MD has
no dynamical bubble, so it imposes the full static barrier that the real
dynamics bypasses. The resolution is that `binding_energy_I_ion_eV` becomes an
**effective parameter calibrated jointly with the drag**, stamped in the same
bundle and guarded at config load against unvalidated pairings. The canonical
trap this names — *do not reduce the drag to force escape* — is the reason the
held-out validation is mandatory: **in-window match ≠ correct production
behaviour**.

---

## 6. Tier 1a — anchored mass dynamics

Tier 1a asks one question: does a *variable* mass change the ion's trajectory
relative to a constant-mass ion, and does the energy ledger still close?

**"Anchored" means the shell schedule is read, not generated.** The He count
n(t) is anchored to the 9 Å TDDFT loss curve through three author-confirmed
points — **(t*, 21), (10 ps, 19), (14 ps, 14)** — interpolated piecewise
linearly and discretized on half-integer downward crossings, giving exactly
**7 shed events**. Because the schedule is fixed externally, the entire
*generative* mechanism (pickup rates, RRK evaporation, the D₀ ladder, the
electronic picture, the internal-energy reservoir) is bypassed. The run is a
clean controlled A/B: `fixed` mass vs `anchored_discrete`.

**The integrator upgrade (SQ1–SQ3).** The BAOAB step is rebuilt every step so
the mass can change; the mass jump is applied at the step seam, at most one
event per step; and the post-jump O-step uses the *new* mass `m⁺`. The
continuous (Meshchersky) alternative `m(t)v̇ = F − v·ṁ` was explicitly
discarded and recorded as discarded: helium leaves as discrete cold atoms, not
as a continuous jet.

**The physics correction that mattered.** The original *cold shed* (He leaves at
rest in the lab frame) implies `v⁺ = m/(m − m_He)·v⁻` — a discontinuous speed
kick at every event, with a telescoping ceiling of ×1.153 over the full 21→14
schedule. Mathematically expected, but **not physically acceptable at Tier 1a**,
where there is no internal-energy reservoir to source the boost. The production
path is therefore **continuous-velocity shedding**: the He leaves co-moving,
`v⁺ = v⁻`, with `ΔE_mass_transfer = +½·m_He·‖v⁻‖²`. Cold shed was demoted to an
explicitly labelled diagnostic upper bound. With this correction there is no
event-local speed kick at all — every trajectory change comes from the lower
post-shed mass feeding subsequent dynamics, and since drag is applied with no
mass scaling (`a_drag = −γ(v)v/m(t)`), the lighter complex decelerates *harder*
per unit force.

**The ledger.** Four terms, with the internal-energy reservoir explicitly
absent:

```
E_kin + E_pot + E_dissip + E_mass_transfer = const
```

Closure was verified by deliberate fault injection: a label-only velocity
change with no matching `ΔE_mass_transfer` blows the residual up. But the
scope of that certificate is stated carefully and should be quoted in any
write-up: it proves the shed is the reduced-mass form and not a relabel, and
that drag work integrates to `ΔE_dissip` exactly. It does **not** certify the
energetics of shedding, which are unsourced at this tier. **"Closure is wiring,
not physics."**

A predictive variant that would have *generated* the shell timing ("Tier 1b")
was **rejected**: TDDFT is not ground truth, and a timing match would have
calibrated the model to the validation regime rather than the production one.

---

## 7. Tier 2 — the generative mass mechanism

Tier 2 replaces the anchor with a mechanism. The `biphasic_energy_gated` model
(`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §1, §4, §6) evolves n(t) by
**two independent discrete stochastic channels per integrator sub-step**. The
equilibrium shell size is *emergent* — there is no target size and no loss-rate
parameter.

**Channel 1 — Poisson pickup**, density-gated and Langmuir-capped:

```
λ_attach = λ₀ · (ρ_He(depth)/ρ_bulk) · (1 − n/n*)₊^p     [ps⁻¹]
P_attach(dt) = 1 − e^(−λ_attach·dt)
```

Each event takes n → n+1 with a **mandatory momentum-conserving velocity
reset**, `v⁺ = (m/(m + m_He))·v⁻` for helium captured at rest. That reset is
what *defines* `E_mass_transfer` in the invariant; a label-only mass change
voids closure. The density-only form was selected because a resting ion still
accretes, which rules out both pure sweeping (Ṁ ∝ v·ρ) and unregularized
dwell-time (Ṁ ∝ ρ/v).

**Channel 2 — energy-gated, RRK rate-limited evaporation.** Two conditions on
the top rung of the I⁺Heₙ dissociation ladder:

- a **parameter-free self-bound gate**: all shedding is suppressed while
  `E_int > Σᵢ₌₁ⁿ D₀(i)`, i.e. while the complex is net self-unbound and is held
  together only by the droplet. Newton cooling drains E_int below the
  integrated binding after a few ps, and the gate opens at the crossing time.
  The threshold is the integrated ladder, not a free onset parameter.
- once self-bound, a classical **RRK rate**:

```
k(E_int, n) = ν · (1 − D₀(n)/E_int)^(s−1)    [ps⁻¹]
```

with **ν = 2.42 ps⁻¹ pinned** from the [IHe05] EPAPS fit
(V″(R_e) = 748.1 cm⁻¹/Å², reduced mass 3.880 amu). Because k saturates at ν,
per-step shedding is bounded and there is **no gate-open avalanche** — the
cascade timescale becomes a prediction rather than an assumption.

The mode count `s` is where the model was **falsified and repaired in public**:
the classical value s = 3n − 3 (= 60 at n = 21) freezes the cascade at n ≈ 20,
while a constant **s_eff ≈ 8** reproduces the anchored 21→19→14 staircase.
Why an effective degrees-of-freedom count an order of magnitude below the mode
count is defensible remains an open literature question (RQ6).

**The internal-energy reservoir E_int.** A per-ion state variable in eV, with
four terms:

- **S2 onset** — `E_int(0) = E₀`, deposited once at ionization;
- **S1 pickup heating** — each attachment releases D₀(n+1), a fraction `f_ret`
  retained;
- **K1 evaporation sink** — each shed drains exactly D₀(n), and the He leaves
  cold;
- **K2 Newton cooling** — `dE/dt = −(E − E_∞)/τ`, applied to the structural
  solvation energy with an **occupancy-resolved** E_∞ (a fixed full-shell value
  would mechanically halt shedding). A no-double-count finding makes this
  clean: a cold shed is energy-neutral for the cooled variable, so K1 and K2 do
  not overlap.

**The five-term invariant**, the cross-cutting correctness gate for every slice:

```
E_kin + E_pot + E_dissip + E_mass_transfer + E_int ≈ const
```

It closes exactly only with the **reduced-mass** capture defect
`ΔE_cap = ½·(m·m_He/(m + m_He))·‖v⁻ − u_He‖²`; the naive `½m_He v²` limit
overstates by ~3 % at n = 1.

**Three stages, not one.** The scored ensemble is not the MD output. A ~30 ps
drag-active MD window hands over to a conservative Coulomb relaxation stage
(to ~8 ns), which hands over to an event-driven RRK cascade running to the
detection time **t_detect = 8.53 µs**. The difference is decisive: at the old
standing point, handover was n̄ = 7.77 with *zero* bare ions, while detection
was n̄ = 4.00 with 16.1 % suppressed, and per-n mean KE differed by up to 26×.
Comparing the MD window directly to experiment would be a category error.

### A mechanism correction worth reading twice

On 2026-07-30 a user challenge caught a wiring error in the project's own
description of itself (`TIER2_PARAMETER_INFLUENCE.md` §2). The prior statement
was "E_int is fed by drag work." It is not: `biphasic_step` books drag loss to
`E_dissip`, and E_int's only sources are the E₀ seed and pickup heat.

The real coupling is **temporal**:

```
E_int^exit ≈ E₀ · exp(−t_res / τ_eff)
```

The drag does not *heat* the complex — it controls **how long the complex stays
inside the droplet**, over which time the internal energy decays exponentially.
Change the drag, change the residence time; change the residence time, change
how much E_int survives to the exit; and how much survives sets how far the
evaporative cascade descends, hence the detected size.

The consequence reframes the model's central tension. The observed
anti-correlation between kinetic energy and the size histogram is **not** two
observables fighting over one energy budget. Reaching n = 1 requires exiting
*hot*, and the ions that exit hot are exactly those that exit *fast* (short
residence ⇒ little drag loss *and* little cooling). It is a **co-selection on
residence time: exit fast ⇔ exit hot.** That single insight unifies four
previously separate facts — why v_c and τ close only as a pair, why softening
the drag tail "un-damps" the cascade, why the suppressed fraction responds to
E₀, and why the source kinetic-energy release is the unique lever that moves
both observables the right way at once.

The episode is worth recording in the thesis as a methodological case study: a
plausible energy-flow story survived for weeks because it predicted the right
correlation, and was caught only by reading the actual wiring.

---

## 8. The experimental comparison layer

**Reference data.**

- `data/reference/integrated_i_he_abundance.csv` — the experimental I⁺Heₙ
  **size distribution**, n = 0…20. The Tier-2 arbiter. Bare I⁺ 43.5 %, I⁺He
  17.5 %; on the solvated branch n₁ = 0.310 and n₁/n₂ = 2.18.
- `data/reference/ihe_ked/IHe_KED_reference.csv` — per-n **mean fragment
  kinetic energy**, falling 3.706 eV at n = 0 to 0.066 eV at n = 17.
- `data/reference/vmi_summary/` — VMI velocity references, retained as
  diagnostics after the 2026-06-29 scope decision restricted Tier-2 scoring to
  the size distribution (VMI is aggregate, not per-fragment).

**The metric.** Histogram agreement uses the 1-D Wasserstein distance, which on
unit-spaced integer support reduces to a CDF-gap sum:

```
W₁ = Σₙ |F_sim(n) − F_ref(n)|        units: helium atoms
```

so a W₁ of 0.5 means the simulated distribution sits about half a rung away
from the experimental one. The headline quantity is **W₁_solv**, computed on
the *solvated* (n ≥ 1) branch, because bare I⁺ is a separate channel (RQ8).

**The scored observable vector.** Fate split (suppressed / trapped / bare),
n̄_det, **n₁_solv**, **midHot** (mean of sim/ref mean-KE over n = 2–8; want 1),
**W₁_solv**, **deepKE** (the same ratio over n = 10–17), χ²_med, and — added
mid-campaign and pinned hard afterwards — **KE₁**, the mean kinetic energy of
detected n = 1 fragments.

KE₁ deserves its own line because it became the sharpest constraint in the
program. The experimental n = 1 kinetic-energy distribution has mean 1.302 eV,
median 1.128, mode 0.891, σ 0.697, with **48.7 % of its mass above 1.15 eV**.
The adjudicated comparison anchor is the peak, **1.00 eV**.

---

## 9. The first landing — and why it is now historical

The first configuration to land the experimental histogram was **`finc1v725`**:
`capped_cubic` with v_c = 7.25 Å/ps and p_tail = −1, τ = 3.2 ps, E₀ = 0.27 eV,
the `rq4graded` ladder, Landau cutoff 0.58 Å/ps, E_bind 0.1168 eV. Its pooled
5 × N = 1000 battery (`TIER2_STAIRCASE_PROBE_FINDINGS.md` §4cc) gave:

| observable | pooled value | reference |
|---|---|---|
| W₁_solv | 0.571 ± 0.095 | 0 |
| n₁_solv | 0.243 | 0.310 |
| midHot | 1.014 | 1.0 |
| n̄_det | 4.37 | 4.889 |
| KE₁ | 1.034 ± 0.128 eV | 1.00 (peak anchor) |
| χ²_med | 242 | gate was ≤ 30 |

What landed: the mid-solvated kinetic-energy curve was correct in **absolute
eV** to 1.4 %, per-bin within ~10 % over n = 2–5, and the solvated histogram
sat about half a rung from experiment — for a mechanism with only three fitted
knobs (v_c, τ, E₀), a genuinely good result. What did not land: χ²_med = 242
against its own ≤ 30 gate, caused by deeply-solvated survivors arriving 30–60 %
too cold as a **monotone slope in n**, not a uniform scale (RQ11).

**Then the ensemble it was scored on turned out to be wrong.** The droplet
geometry is the subject of the next section, and it supersedes these numbers.
One consequence is worth stating here, because it is counter-intuitive:
`finc1v725`'s apparently-excellent KE₁ = 1.034 eV was an **artifact of the
wrong geometry**. The model holds two disjoint n = 1 channels
(`TIER2_SENSITIVITY_ATLAS_FINDINGS.md`:1905-1935) — a shallow-birth channel
whose kinetic energy sits right at the experimental peak (1.04–1.09 eV across
R = 26.6–68.3 Å), and a deep-cascade channel at 0.64 eV. The old geometry
sampled almost exclusively the shallow one. Corrected, there are no shallow
births to offer.

`finc1v725` remains the *standing* production point by decision, because
nothing has been formally adopted to replace it. But as a physical result it is
a milestone, not the answer.

---

## 10. Realistic droplet sizes: the geometry correction

This section is the core of the recent work.

### How droplets are modelled

Two independent laws (`TIER2_PARAMETER_INFLUENCE.md` §15.1):

- **Size.** A nozzle correlation `⟨N⟩ = k₁·p^0.97·T^(−3.88)·d²` (Lackner) sets
  the mean from the experiment's own source conditions, and the distribution
  about it is **log-normal with δ = 0.625** (Kornilov):
  `ln N ~ Normal(ln⟨N⟩ − δ²/2, δ)`. A pickup-and-evaporation chain then
  optionally modifies it.
- **Birth position.** Either a thermal-equilibrium ansatz
  `p(r) ∝ r²·exp(−U_drop(r − R)/k_BT)`, or a `uniform_volume` draw with a hard
  exclusion margin at the surface.

Radius follows from size by the bulk convention
`R = (3N/4πn_He)^(1/3) = 2.2173·N^(1/3)` Å at n_He = 0.0219 Å⁻³.

### What was wrong

The ⟨N⟩ = 2000 pin was **inherited from the 9 Å TDDFT droplet and never
re-derived from the source conditions** (§15.2). The scale of the resulting
error is large:

| production | ⟨N⟩ | R̄ | mean birth depth |
|---|---|---|---|
| legacy MATLAB / port | ~16.4k | ~54 Å | — |
| Tier-2 standing point | 2.0k | 26.6 Å | 9.0 Å |
| **experiment's own source conditions** | **~12 794** | **~49.4 Å** | **29–40 Å** |

A factor of **8.2 in N** — larger than the entire geometry grid that had been
proposed to explore the question. The external anchor is fixed independently
and cleanly: the parent experiment's quoted smallest and largest droplet radii
(34 Å and 68.3 Å) sit at quantiles 0.043 and 0.949 of the **raw** source-
condition log-normal at ⟨N⟩ = 12794 (§15.5), which also settles that the parent
ensemble is the raw distribution rather than the pickup-weighted one.

Crucially, **correcting this introduces no new fitted number.** The nozzle
correlation supplies ⟨N⟩ = 12794 from the preset's own 40 mbar / 14 K. The mean
was never a knob; it had simply been overridden by an inherited pin.

### Why E₀ tracks the geometry — a coherence check

The onset internal energy E₀ enters through the fate cliff `E_ej = E₀·e^(−K)`,
where `K = ∫ρ dt` is the cooling exposure along the ion's path. Longer chords
mean more exposure and demand a larger E₀ to keep ions on the un-suppressed
side. So E₀ is not a free physical constant being fitted — **it is tracking the
geometry**. The measured chain has four points (`TIER2_PARAMETER_INFLUENCE.md`
§4):

```
centre-pinned droplet        E₀ ≈ 0.38–0.41 eV
standing mixture (R̄ 26.6 Å)  E₀ ≈ 0.22–0.27
corrected geometry           E₀ ≈ 0.30–0.37
h405 (corrected, final)      E₀ = 0.405   ← back at the pinned scale
```

The centre-pinned convention traverses the full radius and so accumulates
maximum exposure; the corrected geometry, with its much larger droplets, does
much the same. That the correction restores E₀ to the pinned-droplet scale is a
coherence check on the mechanism, not a coincidence — and it is one of the
stronger physical arguments the program has produced.

### The staged correction, and what it costs

The correction was executed in five stages, G0–G4
(`TIER2_SENSITIVITY_ATLAS_PLAN.md` §3.5). G0 froze the specification —
notably **rejecting** the analytic Kornilov prior arm, because at ⟨N⟩ = 12794
roughly a quarter of the probability mass falls outside its [250, 16000]
window. G2 **adopted** the corrected geometry as the target on 2026-07-27,
with `finc1v725` standing until a successor was delivered. The realized
ensemble is R_q50 ≈ 49.5 Å with mean birth depth 34.6 Å.

**What it buys.** Five entries retire from the model's calibration ledger with
no new knob introduced: the ⟨N⟩ = 2000 pin, the `uniform_volume` birth law, the
3 Å birth margin (which had been measured as load-bearing at 3.4σ on n₁_solv),
and two scaling conventions.

**What it costs.** The standing point does not survive it. Measured at the
corrected geometry (`TIER2_SENSITIVITY_ATLAS_FINDINGS.md`:1125-1145),
`finc1v725` goes:

```
trap  0.0416 → 0.3080      supp  0.2019 → 0
n̄     4.387  → 16.99       n₁    0.2429 → 0
W₁    0.678  → 12.11
```

That is not a degradation, it is a collapse: at realistic droplet sizes the ions
no longer escape. The whole (v_c, τ, E₀) surface had to be re-arbitrated.

The correction also **creates** one new problem, honestly logged as the first
ledger row it added rather than retired (§17, §14.2): the three-stage timescale
separation (30 ps drag-active MD → 8 ns conservative relaxation → µs free
flight) was calibrated at R ≈ 27 Å where ejection is effectively instantaneous.
At R ≥ 49 Å, escape takes ~0.1–1 µs with helium still present, and one third to
one half of ions never leave at all.

### What was actually run

At realistic droplet sizes, roughly:

- a controlled fixed-R grid — 11 cells × N = 500 at R = 26.6 / 34.0 / 49.4 /
  68.3 Å across three birth laws;
- the corrected-geometry re-arbitration — a 14-cell confirmation ring at
  N = 500, six finalists and a nine-cell E₀ ladder arm at N = 1000;
- the successor's verification battery — 5 × N = 1000 across independent seeds;
- four mechanism probes at N = 1000 each (drag tail exponent, drag-state
  coupling, source-budget slope, and the CE-channel/exit-strip programme);
- the free-form linear campaign of §13 — a further 12 cells of MD at N = 500.

Order **47 500 ion trajectories**, each carrying 30 ps of drag-active MD plus
8 ns of relaxation plus the µs cascade — tens of hours of wall-clock at
concurrency 3, supported by twin forward-model scans of order 10⁴–10⁵ cells at
zero MD cost.

**One caveat to state explicitly.** The ensemble-level verdict that the
corrected *size distribution* breaks the landing (trap 0.31–0.42, W₁ 9.0–9.8)
came from a zero-MD grid re-weighting whose own pre-registered oracle was ruled
**inadmissible** — 53.8 % of the standing density lies below the grid's support
edge (§15.7). The MD campaigns above, not the re-weighting, are the
load-bearing evidence.

---

## 11. h405 — the current best result

**h405 is not new physics.** It is the same `capped_cubic` model relocated on
the (v_c, τ, E₀) surface at the corrected geometry, with every other knob at
the standing pins:

| knob | `finc1v725` | **h405** |
|---|---|---|
| v_c | 7.25 Å/ps | **5.50 Å/ps** |
| τ | 3.2 ps | **4.4 ps** |
| E₀ | 0.27 eV | **0.405 eV** |
| geometry | old (R̄ 26.6 Å) | **corrected (R_q50 49.5 Å)** |
| form, p_tail, ladder, Landau, E_bind, mechanism | — | unchanged |

Pooled 5 × N = 1000 battery, seeds 20260730–34
(`TIER2_SENSITIVITY_ATLAS_FINDINGS.md`:1749-1800):

| observable | h405 pooled | reference |
|---|---|---|
| n̄_det | 3.877 | 4.889 |
| n₁_solv | 0.209 | 0.310 |
| **W₁_solv** | **0.765** (per-seed SD 0.034) | 0 |
| midHot | 0.946 | 1.0 |
| deepKE | 0.504 | 1.0 |
| suppressed / trapped | 0.194 / 0.078 | — |
| **KE₁** | **0.641 ± 0.003 eV** | 1.00 (peak anchor) |
| KE₂ | 0.549 | 0.706 |
| χ²_med | 461 | — |

*(Three vintages of these numbers circulate and are not interchangeable: the
pooled battery values above, a single-seed committed row at W₁ 0.709 /
KE₁ 0.637, and the twin forward-model value KE₁ 0.6236. This document uses the
pooled battery throughout.)*

**On the comparison with `finc1v725`.** W₁ 0.765 against the old 0.571 is
**not** a regression — the two are not the same physical ensemble, and on this
ensemble the incumbent scores 12.11. The seed scatter also tightened by 2.8×
(SD 0.034 vs 0.095), which is what one expects when the sampled geometry stops
straddling a cliff.

**Why it is a candidate and not production.** The adoption trigger was
pre-registered as three gates. Two passed — every battery member gates
individually, and the pooled W₁ landed inside its band. The third **failed**:
the pooled joint score 1.734 came out above the comparison cell's 1.683, the
single-seed 1.559 having been a favourable draw. So nothing fired
automatically, and the decision reverted to a user call that has not been made.
Every session since has closed with the same sentence: *"Nothing adopted;
`finc1v725` stands; h405 candidacy + G4 adjudications open."*

### The KE₁ deficit — the sharpest open constraint

h405's 0.641 eV against the 1.00 eV anchor is a **0.36 eV shortfall**,
contributing +5.84 to the joint score — the largest single defect on the
surface. Stated without reference to mass: the model's n = 1 ions exit at
9.58 Å/ps where the anchor wants 12.14 Å/ps, a 27 % velocity shortfall.

And it is **unreachable within the calibrated surface**. The registered
reachability verdict bounds the single-knob contributions at E₀ ≤ +0.038,
τ ≤ +0.094, v_c ≤ +0.054 — an additive bound of **+0.186 eV against a required
+0.309**. Every escape route beyond the surface has since been measured and
closed:

- **the drag tail exponent** — softening p_tail does reach KE₁ 0.917–1.428, but
  every cell that clears 0.95 carries midHot 3.0–4.4, far outside band. Kill
  fired.
- **drag-state coupling** — refuted. The exit toll is paid *dressed*: no ion
  reaches n ≤ 8 while still inside the droplet, so the coupling is gate-clipped.
- **exit stripping** — a counterfactual on the committed h405 checkpoints,
  converting the *entire* suppressed population to n = 1 at existing velocity,
  reaches only **KE₁ 0.708**. Necessary, not sufficient.
- **the CE channel mixture** — a source-side kinetic-energy spread reaches
  KE₁ 0.920 and reproduces the high-energy tail (above-1.15 fraction 0.418 vs
  the experimental 0.487), genuinely breaking the needle — but traps the slow
  channel and wrecks the histogram (W₁ → 2.97). The line was shelved
  2026-07-30 as a measured boundary.

The positive result that survives all of this: **KE₁ position is source
physics, not drag or mass-mechanism physics.**

---

## 12. Is the velocity cap physically defensible?

The production drag law is `capped_cubic`: the Tier-0 cubic below v_c, and a
saturated tail above it. The obvious objection — recorded verbatim in the docs
as the user's own — is *"can we reach the h405 result without capping? This is
so unphysical."* It launched the campaign in §13. The honest answer has two
halves.

**The defence** (`TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md` §6.4). v_c = 5.5
Å/ps is **Mach 2.3** against the helium sound speed of ~2.4 Å/ps. Above it, a
constant force means a constant energy loss per ångström — the signature of
wave-drag and vortex-shedding regimes. Saturating drag on a supersonic bubble
is respectable physics, not a fudge. What *is* idealized is the sharp corner —
and even that is forced rather than chosen: any smooth saturation family, e.g.
`F = b·v³/(1 + (v/v_c)^m)`, departs from the Tier-0 cubic near the top of the
calibrated band unless m is large, and **m → ∞ is exactly the cap**. The one
closed-form alternative tried, a Padé saturating cubic, is excluded by
arithmetic: in-band cubicity forces its saturation scale above ~15 Å/ps, which
forfeits saturation at the production velocities entirely.

**The critique, in its strong form.** `TIER2_PARAMETER_INFLUENCE.md` §17 is a
*physical-sensibility ledger* that classifies every element of the model as
**P** (physics-constrained), **C** (convention), **E** (effective/fitted),
**S** (scaffolding) or **M** (missing). It classes the in-band form and b as
**P, band-limited** — and the capped tail (v_c, p_tail) as **E: effective, no
external constraint**. The column headed "what would retire it" names *TDDFT at
production kinematics*, and a domain input recorded 2026-07-30 states exactly
that as **infeasible** — TDDFT itself breaks down there.

So the cap is not unconstrained-for-now, the way τ and E₀ are. It is
**unconstrainable in principle**: a permanent effective element parked in a
velocity range no instrument can ever reach. Two supporting tensions point the
same way: three different observables prefer three different v_c (KE χ² wants
≥ 8.5, the n₁ anchor ~6.2–6.3, W₁ ~7.25), which is what absorbing missing
physics looks like; and p_tail is not a usable second lever, since a more
negative exponent re-heats the kinetic-energy curve and re-strips the shells
together rather than orthogonally.

**The verdict to carry into the thesis.** The criticism that lands is not
"the cap is unphysical" — it has a regime story. It is that the cap is a fitted
element hidden where nothing can check it. The measurement that would convert
its idealized corner from an aesthetic judgement into a number is a Hill-type
smooth-saturation twin scan (b locked at the Tier-0 value, v_c and m free),
measuring *how sharp the corner has to be*. It is designed, posted, and not
taken.

---

## 13. The free-form linear counterfactual — and why it fell short

**The premise, which is measured rather than assumed.** Detected n = 1 ions
exit at ~12.3 Å/ps and n = 2 at ~10.4 Å/ps; at the corrected geometry h405's
KE₁ corresponds to ~9.7 Å/ps. All of them live far above v_c = 5.5. So for the
fast class — the ions that produce the observables — the operative drag curve
is *already entirely experiment-fitted*, and there is no Tier-0 authority in
play. Sweeping a free-form linear coefficient is therefore "not worse than
fitting v_c": it trades the in-band anchor for one fewer structural assumption,
with no cap and no kink.

**The form.** `γ(v) = ρ̂·a`, a in amu/ps — a single parameter setting the force
at every velocity. The maximal-rigidity limit, and deliberately so: it has **no
second lever**, so the fast class cannot be softened without committing the
slow class to the consequences.

**The campaign.**

| stage | scale | verdict |
|---|---|---|
| Step 0 — twin authority | zero MD, 26-cell replay | KE₁ ranking **licensed**, Spearman ρ = 0.9979 |
| Arm 1 — `linscan` | 12 312 cells, zero MD | basin exists: 35 cells gate; best twin KE₁ 0.892 vs h405's 0.624 |
| Arm 2 — `linqscan` | 5 184 cells, zero MD | **the ram term is rejected**; c = 0 carries every best cell; pure linear selected within the family |
| MD ring | 8 × N = 500 | **it lands** — 6/7 gate under both policy arms; best MD KE₁ **0.903** vs 0.637, the first drag-side KE₁ above the exit-strip ceiling |
| h405-clone battery | 3 × N = 500 | **KE-equivalent, not landing-equivalent**: ΔKE₁ −0.006 eV, midHot 0.891 (cooler than h405's 0.946) — but W₁ 1.079 vs 0.765 |
| (a, τ) joint ring | 6 × N = 500 | the joint-improvement region is **empty** |
| τ refinement | 14 364 cells, zero MD | basin triples; still no cell better on both axes |
| decisive pair | 2 × N = 500 | pre-committed **kill fires**; family closed 2026-08-13 |

**How close it got.** The nearest approach, cell `t1` (a = 35, τ = 4.0,
E₀ = 0.45), lands **W₁ 0.7805 against h405's 0.7675** — a gap of +0.0130
against a per-seed SD of 0.031, i.e. **~0.4σ, statistically indistinguishable**
— while genuinely beating h405 on KE₁ by +0.0205 eV (≈9σ on that observable's
seed SD). It fails only because the pre-registered test required a strict beat
on W₁, and because another cell dominates it outright. **The family never
produced a single cell better than h405 on both W₁ and KE₁.** The Pareto front
is unchanged: h405 at (0.767, 0.637) and `lr6` at (0.785, 0.715).

Two costs are worth recording. The mid-band overheating that the MD ring showed
(midHot 1.39–1.98) turned out **not** to be a property of constant γ — the
clone at higher a came back at 0.891, cooler than h405 — so it is a coordinate
effect, not a form effect. But midHot in-band and KE₁ above h405 remain mutually
exclusive across the entire family.

**The structural diagnosis**, which is the durable product of the arm: the
failure is neither the twin transfers nor the extrapolation. The licensed legs
extrapolated essentially perfectly outside their training range (KE₁ predicted
0.659 / 0.506 against measured 0.6576 / 0.5071). It is the functional form of
the predictor — **W₁ is not a function of (n₁, tail)**. Cell `t1` carries a tail
17 % better than h405's and still loses, because the n = 2–9 interior, which
neither input constrains, holds the remaining probability mass. This is the
same lesson the clone battery gave, recurring a third time: **the moments
match, the shape does not.**

### Why linear is not TDDFT-based at all

Two layers, and the second is the sharper one.

**First, Tier 0 rejects it in-band.** Given a free exponent, the TDDFT traces
choose n̂ = 2.927, and the linear coefficient runs to zero whenever a higher
power is allowed to compete. That is a direct measurement, not a preference.

**Second — the deviation-accounting argument.** The linear line crosses the
Tier-0 cubic at `√(a/b) = 4.11 Å/ps` (a = 42.5, b = 2.5153) — which is
**inside** the TDDFT-calibrated window of 2.54–4.95 Å/ps. Measured against the
Tier-0 law, that line runs

- **×4.2** the cubic at v = 2 Å/ps,
- **×1.9** at v = 3 Å/ps,
- **31 % soft** at the band top, 4.95 Å/ps.

So removing the cap does not remove the model's free element — it **relocates**
it, from above 4.95 Å/ps where nothing can ever check it, into the one velocity
window where the program has ab-initio authority and has already ruled.

The program's own adjudication states it better than a paraphrase can:

> **Both forms carry a fiction; the cap hides it where nothing can see it, the
> smooth line hides it where the traces can. By the validation hierarchy that
> is a worse trade, bought with aesthetics.**

---

## 14. The current open arm: what mass actually flies?

`TIER2_MASS_SCENARIOS.md`, opened 2026-08-11. Status: **open investigation,
nothing adopted, no implementation authorised.**

### The governing identity

Trajectories are integrated at the **dressed** mass — and at the corrected
geometry every fragment is born with the full n = 21 shell, because the pickup
gate is saturated at those depths. But the scored observable uses the **bare**
mass at detection. Therefore any energy gained or lost while dressed is
attenuated in the observable by

```
∂KE₁/∂E_dressed = m(1)/m(n_flight) = 130.903/210.955 = 0.6205
```

This is not a bug. It follows from the Tier-1a velocity-preserving shed (an
evaporating atom leaves at the complex velocity carrying its own ½m_He v²)
combined with the measured fact that ions exit the droplet fully dressed. Pay
the toll dressed, shed outside — the bookkeeping is self-consistent. **The
question is whether the premise is physically defensible: that 21 helium atoms
co-move through a 2.7 eV Coulomb explosion.**

### It also solves the "only 55 % of the binding energy shows up" puzzle

A long-standing oddity: the measured transfer function is
`dKE₁/dE_bind = −0.542` (capped) and `−0.549` (linear), MD-confirmed at
−0.5308. An ion that pays a 0.1168 eV solvation well appears to lose only ~55 %
of it from its kinetic energy. The earlier reading called the missing 45 % a
"refund of un-incurred drag." Three measured invariances killed that story
(`TIER2_PARAMETER_INFLUENCE.md` §9.6): the quantity is invariant to drag
*form*, to birth *depth*, and to density *width*. **A quantity indifferent to
everything about the drag is not a drag effect.**

The answer is the **mass-frame partition**. The toll is paid in full — 98.1 %
in the dressed frame, and the energy ledger reads 0.1167 eV against E_bind
0.1168, i.e. 99.9 %. The observable simply only ever sees the iodine's share of
it, 0.6205. The rest was paid by helium that evaporates afterwards and is no
longer part of the scored particle. The full decomposition of the apparent
45 % shortfall is 38.0 % mass partition + 6.7 % cascade/bin membership + 1.9 %
genuine drag refund.

The consequence is a hard structural ceiling:

```
∂KE₁/∂(any pre-evaporation energy) ≤ m(1)/m(21) = 0.6205
```

### The energy ledger, and the shedding-cost law

For h405, over 6375 scored ions:

```
source budget                       2.7006 eV
− solvation toll                    0.1167   (99.9 % of E_bind)
− drag dissipation                  1.5790   (58 % of budget)
= exit kinetic energy, dressed      1.0050
× m(1)/m(21) = 0.6205
= scored KE₁                        0.6236   ← matches the committed value
```

Re-integrating the same chord at bare mass gives exit energy invariant to
0.09 % and KE₁ = **1.0059 eV**. The invariance is not luck: under p_tail = −1
the drag force above v_c is constant, so the work done depends on the *path*,
not on how fast it is walked. Generalizing:

```
KE₁ = E_exit − ½ · Δm · v_shed²
```

The entire KE₁ question collapses to **one number — the speed at which the
helium departs.** Because the cost scales as v², late shedding is expensive and
the velocity peak is the worst possible moment; reaching KE₁ = 1.00 eV requires
v_shed ≲ 1.2 Å/ps, i.e. effectively never having been dressed. Both existing
measurements fall on the curve, including the twin-vs-MD KE₁ gap that had
previously gone unexplained.

Worth noting for anti-circularity: the drag surface was arbitrated against the
*size distribution*, never against KE₁, and it independently leaves
E_exit = 1.006 eV — matching the kinetic-energy anchor to 0.6 % once the mass
frame is corrected. That is one coincidence. Suggestive, not a fit.

### The obstruction

The identity applies to *every* bin, so a uniform flight-mass change lifts every
bin equally — precisely the failure mode that killed the p_tail route. What the
data requires is not a boost but a **U**: ×1.60 at n = 1, ×1.31 at n = 2, ×1.05
at n = 3, then *cooling* by 8–16 % across n = 4–8, rising again to ×1.31 at
n = 15 and ×2.45 at n = 17.

Inverted into flight masses, that U is a near-**threshold**: strip hard above
~8.1–8.8 Å/ps and not at all below. An independent estimate from co-moving
energetics gives the same crossing (½m_He v² = 15.0 meV at 8.5 Å/ps against a
charge-induced-dipole binding of ~18.2 meV at 3 Å) — from physics that knows
nothing about the kinetic-energy spectrum. **But the model's own ladder
disagrees**: D₀ ≈ 9.2 meV for the core rungs implies a threshold at 6.7 Å/ps,
which would strip into the mid band and reproduce the p_tail failure.

### Hence M1 gates everything

Once D₀ is known the threshold is not free: `v_threshold = √(2D₀/m_He)`. The
candidate range spans a factor of four, and **three of the four outcomes are
failures**:

| D₀ source | value | threshold | outcome |
|---|---|---|---|
| the model's own ladder | 9.2 meV | 6.7 Å/ps | strips into the mid band — the p_tail failure |
| what the KE data needs | ~15 meV | 8.5 Å/ps | works |
| charge-induced dipole at 3 Å | 18 meV | 9.3 Å/ps | marginal |
| charge-induced dipole at 2.5 Å | 38 meV | 13.6 Å/ps | nothing strips; route dead |

The document's own summary: *"This is not a robust mechanism; it is a narrow
window in a quantity nobody has sourced."* **M1 — obtaining a sourced I⁺–He
dissociation energy — is a literature question, not a compute one, and it
precedes everything.**

The remaining open items, briefly: **M2** re-score the cooling exposure, since a
faster ion crosses in less time and the size distribution must be re-checked;
**M3** test whether exit-energy mass-invariance survives outside the capped
form; **M4** measure the evaporation-recoil channel (estimated ~20–40 meV over a
full cascade, biased toward low n); **M5** the pivot — *is the deep n = 10–17
tail load-bearing?*; **M6** design a speed-gated strip, which cannot be
specified before M1 and M5.

**M5 deserves a supervisor conversation, not a modelling decision.** Preserve
the deep tail and the ion must fly heavy, which fixes KE₁ near 0.62 and leaves
the deficit standing. Release it and a graded partial strip reaches ~0.83 with
the mechanism intact, but caps the detected size at n ≈ 8. The evidence to bring:
the n = 17 bin holds 487 of 20 000 scored ions; deepKE is already an
acknowledged 16 % cold; the n = 20–21 spike is scoring scaffolding rather than
measured structure; and the experimental n = 19–20 population is 0.59 % against
the model's 17.6 %, a factor of ~30. As the document says, this is *an
experimental judgement, not a modelling one.*

---

## 15. Open questions, standing cautions, and where things live

### Research questions

| # | question | status |
|---|---|---|
| RQ1 | Provenance and magnitude of E_int(0) | partially relieved — the histogram inverts into a p(E₀) coinciding with the sourced band untuned |
| RQ2 | Per-shed kinetic-energy release ε | open; bounded small, estimated 20–40 meV over a cascade (see M4) |
| RQ3 | Fate of a net self-unbound complex | partially un-anchored — the bare bin now reads as channel branching |
| RQ4 | Ladder-bottom depth and the n₁/n₂ = 2.18 ratio | open; the flat Form-U bottom is falsified at ~2× |
| RQ5 | µs-flight channels beyond RRK evaporation | open |
| RQ6 | Why s_eff ≈ 8–12 rather than the classical 3n−3 | open literature question |
| RQ7 | Production Coulomb kinematics | superseded in part by the mass-frame ceiling |
| RQ8 | Bare-I⁺ channel provenance and branching | open |
| RQ9 | Physicality of the fitted cooling time τ | parked by decision |
| RQ10 | Shed-frame / momentum convention | resolved to co-moving in production |
| RQ11 | Deep-bin KE undershoot (monotone slope over n = 10–17) | open; the mass thread is its current owner |

Three facts localize what remains (`RESEARCH_QUESTIONS.md` §2): the two largest
experimental bins are structurally unreachable from the cascade side; the gated
cascade never self-terminates, because a shed drains only D₀(n); and the
initial reservoir has no defensible provenance at the values that make the
mechanism work.

### Two standing cautions

**Form-blindness, in its final precise form.** The solvated-histogram
Wasserstein distance **cannot see the drag form**, even at matched binding well
(Δ +0.009, not significant). But the fate split and the kinetic-energy/size
observables separate cubic from quadratic at 8–35σ. So the forms *are*
physically distinguishable — just not by the arbiter the tier chose. Two rules
follow, and both are load-bearing:

1. **W₁ must never be used as a form discriminator.**
2. **The experimental landing can never be cited as evidence for any drag
   form, cubic included.** Form authority lives solely in the Tier-0 held-out
   traces.

**Tier 3 is measured but not started.** The model ensemble is under-dispersed
against the VMI panels — the second moments are wrong. The noise model remains
stubbed behind its enum, awaiting its own plan.

### Reading order

| to understand | read |
|---|---|
| the pre-drag physics baseline | `docs/matlab_port/PHYSICS_BASELINE.md` |
| the drag architecture and its alternatives | `DRAG_PORT_DESIGN_DECISIONS.md` |
| which parameters are sourced, derived, bounded or free | `CALIBRATION_MAP.md` |
| the mass mechanism in full | `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` |
| the drag-form extraction and its verdict | `docs/drag_port/Tier0/METHOD_B_trajectory_matching_extraction.md`, `TIER0_FINDINGS.md` |
| anchored mass validation | `docs/drag_port/Tier1/TIER1A_IMPLEMENTATION_PLAN.md` |
| open physics questions | `docs/drag_port/Tier2/RESEARCH_QUESTIONS.md` |
| per-knob influence, and the physical-sensibility ledger (§17) | `docs/drag_port/Tier2/TIER2_PARAMETER_INFLUENCE.md` |
| the current open arm | `docs/drag_port/Tier2/TIER2_MASS_SCENARIOS.md` |
| chronology and every verdict | `docs/drag_port/Tier2/drag_migration_log_tier2.md` |

Four of these are very large — the Tier-2 migration log, the sensitivity-atlas
findings and plan, and the staircase-probe findings run to thousands of lines
each. Open them by section heading, not from the top.

### External references cited across the program

- Albrechtsen, Stapelfeldt et al., *Nature* **623**, 319 (2023) — ion solvation
  in helium droplets; the biphasic mechanism's foundational paper.
- Albrechtsen et al. (2025), arXiv:2502.11783 — dissipated solvation energy,
  droplet-size independent.
- Time-resolved alkali TDDFT, *J. Chem. Phys.* **160**, 164308 (2024) —
  Newton's law of cooling for the first ~5 ps.
- González group, *PCCP* (2016), 10.1039/C6CP04315A — two-regime ICVF
  mechanism; viscous-flow friction.
- Braun & Drabbels, *J. Chem. Phys.* **127**, 114303 & 114304 (2007) — fragments
  escape dressed in finite Heₙ; size correlates with speed.
- Buchachenko et al., *J. Chem. Phys.* **122**, 194311 (2005) — the sole
  ab-initio I⁺–He datum; anchors D₀(1).
- Bartl, Scheier & Echt, *J. Phys. Chem. A* (2013), 10.1021/jp406540p, and the
  *Int. Rev. Phys. Chem.* (2020) review — snowball distributions as evaporative-
  ensemble products; first shell 20 for I₂⁺.
- Kornilov et al. — log-normal droplet-size distribution, δ = 0.625.
- Lackner et al. — ⟨N⟩ nozzle scaling.
- *Phys. Rev. A* **107**, 023104 (2023) and arXiv:2401.09211 — surface Coulomb
  explosion as the limiting case.
- García-Alfonso (2024) — single-ion energy lowering.
- NIST ASD — I⁺ spin–orbit levels.
