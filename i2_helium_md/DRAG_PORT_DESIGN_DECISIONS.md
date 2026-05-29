# Drag-Model Port — Design Decisions

**Status:** Open planning document. Each section captures one architectural
choice in the migration from the hard-sphere collision model
(`i2_helium_md/physics/collisions.py`) to a TDDFT-calibrated drag-force
model for I⁺ in a helium bubble.

**Format.** Each section has the same structure:

- *Physical question* — what is being decided and why it matters.
- *Primary* — the approach chosen for the first implementation.
- *Secondary* — viable fallbacks to keep pluggable behind the same
  interchangeability surface.
- *Discarded* — options ruled out, with the reason recorded so they do
  not get re-litigated.
- *Interchangeability surface* — the module boundary or `SimConfig`
  field behind which the choice lives, so swap-in design stays visible.

The document is intentionally implementation-agnostic. No code, no
pseudo-code. Mathematical formulation in LaTeX where needed.

---

## 1. Noise model

### 1.1 Physical question

The hard-sphere model is not "deterministic dissipation plus noise" —
it *is* a stochastic process whose mean happens to look like drag.
Replacing it with a deterministic $F_\text{drag}(v)$ alone discards the
entire fluctuation channel: trajectory ensemble spread, transverse
randomisation, rare-event tails, and the thermal floor at the He bath
temperature $T_\text{particles\_K} = 0.4\,\text{K}$.

Adding a Langevin noise term restores stochasticity, but three sub-
decisions follow immediately and they couple:

- **Functional form of the noise term** (additive vs. multiplicative
  vs. empirically-anchored).
- **Calibration purpose** (what physical thing the amplitude is meant
  to reproduce).
- **Geometric structure** (longitudinal-only, isotropic 3D, or
  anisotropic with separate transverse component).

These three are decided independently below.

---

### 1.2 Functional form

The Langevin equation under consideration:

$$m\,\dot v = F_\text{ext}(t) - \gamma(v)\,v + \mathcal{N}(v, t)$$

**Friction-coefficient convention (used document-wide).** $\gamma(v)$ is
a *force coefficient* with units $[\text{amu/ps}]$, defined directly
from the extracted drag law as $\gamma(v) = |F_\text{drag}(v)|/v$ — so
the friction force is $\gamma(v)\,v$ (units $\text{amu·Å/ps}^2$, a
force, with **no** leading $m$). The corresponding friction *rate* is
$\gamma(v)/m$ $[\text{1/ps}]$; it appears only inside the BAOAB damping
exponent $e^{-\gamma\,dt/m}$ (§4.3), never as a multiplier on the force.
This is the single source of the factor-of-$m$ that distinguishes the
two appearances.

$\gamma(v)$ is $v$-dependent under either drag form (§3): for the
primary linear+cubic drag, $\gamma(v) = a + b\,v^2$ (finite at $v=0$,
units amu/ps since $[a]=\text{amu/ps}$, $[b\,v^2]=\text{amu·ps/Å}^2
\cdot \text{Å}^2/\text{ps}^2 = \text{amu/ps}$); for the secondary
power-law drag with $n<0$, $\gamma(v) = \gamma_\text{PL}\,v^{\,n-1}$
(where $\gamma_\text{PL}$ is the power-law coefficient from §3.8's
$\{\gamma,n\}$ bundle; divergent at $v=0$). This $v$-dependence is what
makes the choice non-trivial.

#### Primary — N2: multiplicative noise with local FDT

$$\mathcal{N}(v,t) = \sqrt{2\,\gamma(v)\,k_B T_\text{eff}}\;\xi(t)$$

with $\langle\xi(t)\rangle = 0$, $\langle\xi(t)\xi(t')\rangle =
\delta(t-t')$ (so $\xi$ is white noise, $[\xi(t)]=\text{ps}^{-1/2}$).
There is **no** leading $m$ inside the root: the textbook form
$\sqrt{2\,m\,\gamma_\text{rate}\,k_BT}$ uses the friction *rate*
$\gamma_\text{rate}=\gamma/m$, so $m\,\gamma_\text{rate}=\gamma$ and the
$m$ cancels in this convention. Dimensional check:
$[\,2\gamma k_BT\,] = (\text{amu/ps})(\text{amu·Å}^2/\text{ps}^2) =
\text{amu}^2\text{Å}^2/\text{ps}^3$; its root times
$[\xi]=\text{ps}^{-1/2}$ gives $\text{amu·Å/ps}^2$, a force. Balances.

The noise amplitude tracks the *local* friction coefficient, so the
fluctuation–dissipation relation is honoured instantaneously rather
than in some averaged sense. The price is that the SDE is now
non-trivial:

- **Itô/Stratonovich ambiguity.** Multiplicative noise SDEs give
  different equilibrium distributions under different stochastic
  calculi. Stratonovich is the standard physics choice (limit of
  smooth-noise OU processes) and is what the BAOAB-family integrators
  assume.
- **Low-$v$ regularisation — only under the secondary power-law drag.**
  Under the primary linear+cubic drag, $\gamma(v) = a + b\,v^2$ is
  finite at $v=0$, so the noise amplitude is well-defined down to rest
  and no regulariser is needed. Under the secondary power-law drag with
  $n<0$, $\gamma(v)\to\infty$ as $v\to 0$ and the noise amplitude
  diverges along with the drag; the same regulariser used for the
  deterministic drag (see the §3.8 regularisation note) must then apply to the noise.
- **$T_\text{eff}$ is not necessarily $T_\text{particles\_K}$.** See
  §1.3 for the calibration discussion — the effective temperature in
  the FDT relation is a fit parameter under noise purpose (b).

*Why primary:* the extracted drag law is nonlinear, so the friction is
$v$-dependent in any case. Anchoring the noise to that same $v$-
dependence is the most internally consistent choice. The Itô/
Stratonovich choice is settle-once and the integrator handles it.

#### Secondary — N3: empirical noise from TDDFT residual variance

$$\sigma_\xi^2 = \text{Var}\big(v(t) - v_\text{smoothed}(t)\big)\Big|_{\text{post bubble-mode removal}}$$

The amplitude is read directly off the residual of the CEEMDAN+SG
smoothing pipeline already in `Drag_extraction_code.md`. No FDT, no
temperature, no model assumption beyond "what the data wasn't
smoothed away into the drag fit is, by definition, noise."

*Why secondary:* model-free and anchored to the same dataset as the
drag itself — same provenance, same systematic biases. But:

- Only one TDDFT trajectory per case (9 Å, 18 Å). Cannot separate
  trajectory-to-trajectory variance from within-trajectory variance.
- Conflates numerical noise, mean-field error, and any residual
  bubble dynamics that escaped the CEEMDAN IMF drop. Not the same
  thing as a He bath kicking the ion.
- Gives one number per case, not a function of $v$. So it would be
  applied as an additive (constant amplitude) noise — losing the
  $v$-dependence that motivated N2 in the first place.

Keep as a trial alternative when comparing against the TDDFT curves —
it's the natural "atheoretical" baseline.

#### Discarded — N1: additive noise FDT-anchored at $\gamma_0 = \gamma(v\to 0)$

$$\mathcal{N}(v,t) = \sqrt{2\,\gamma_0\,k_B T}\;\xi(t),\qquad \gamma_0 = \lim_{v\to 0}\gamma(v)$$

Rejected on the dynamical-scale ground, with a note on form-dependence:

- Under the primary linear+cubic drag, $\gamma_0 = a$ is finite, so
  N1 is mathematically well-defined. But anchoring the noise scale at
  the low-velocity limit while the ion spends most of its time at
  eV-scale KE sets the noise by the wrong dynamical scale — the
  fluctuations the ion actually experiences are governed by $\gamma$ at
  its working speed, not at rest.
- Under the secondary power-law drag with $n = -2$, $\gamma_0\to\infty$
  and the amplitude is undefined outright. N1 is then not merely
  ill-scaled but mathematically incompatible with the drag form.

In both cases N2 (local FDT) dominates: it reduces to the *correct*
limit of N1 at low $v$ for the linear+cubic form while tracking the
working-speed friction everywhere else.

#### Interchangeability surface

`SimConfig.noise_form ∈ {multiplicative_local_fdt, empirical_residual, none}`.
The integrator dispatches on this. Adding new variants later requires
adding a branch but not changing the field signature.

---

### 1.3 Calibration purpose

Even with the functional form fixed (N2), the noise amplitude depends
on a temperature-like parameter $T_\text{eff}$ (or, for N3, on the
residual variance). What that parameter is *calibrated to* is a
separate choice.

#### Primary — (b): calibrate to hard-sphere trajectory variance

Run the existing hard-sphere model on the reference HeDFT cases (9 Å,
18 Å) and on the `single_pulse_droplet_distribution` preset. Extract
the empirical variance of:

- final ion speed across the ensemble,
- longitudinal velocity at intermediate checkpoints,
- (for transverse noise, if T3 active — see §1.4) angular spread of
  final velocity vectors.

Fit $T_\text{eff}$ such that the new Langevin model reproduces those
variances on the same presets.

*Why primary:* this is the only choice that operationally preserves
what the hard-sphere model produced. The drag law itself is calibrated
to reproduce the *mean* TDDFT trajectory; calibrating the noise to
hard-sphere variance fills in the second moment using the only
ensemble-resolved data we have. Imperfect but pragmatic.

*Caveat to flag explicitly:* this is mildly circular — we are
enshrining the discarded model's statistical fingerprint as the target
for the new model. Acceptable because hard-sphere was itself tuned
against VMI data, so its variance is at least indirectly anchored to
experiment. But the relationship is not tight.

#### Secondary — (c): calibrate from TDDFT residual

Use the same residual-variance machinery as the secondary functional
form (N3) to set $T_\text{eff}$ even when the *form* is N2. I.e., back
out $T_\text{eff}$ such that

$$\sqrt{2\,\gamma(v_\text{typ})\,k_B T_\text{eff}} = \sigma_\text{residual}$$

at some representative velocity $v_\text{typ}$ (median or RMS over the
extraction window).

*Why secondary:* lets purpose-(c) and form-N3 share a calibration
pathway. If primary calibration (b) gives suspicious results, this is
the natural cross-check using independent data.

#### Secondary — (a): strict FDT at $T_\text{particles\_K} = 0.4\,\text{K}$

Take $T_\text{eff} = 0.4\,\text{K}$ as a hard physical constraint and
let the noise amplitude be whatever it is.

Order-of-magnitude estimate: with $m \approx 203\,\text{amu}$,
$k_B T \approx 3.4\times 10^{-5}\,\text{eV} \approx 0.33\,
\text{amu·Å}^2/\text{ps}^2$, $\gamma_0 \sim \text{few amu/ps}$,
$dt = 0.01\,\text{ps}$, the per-step OU kick (coefficient convention,
$\Delta v \sim \sqrt{2\gamma_0 k_B T\,dt}\,/\,m$) is

$$\Delta v_\text{noise} \sim \frac{\sqrt{2\gamma_0 k_B T\,dt}}{m}
\sim 7\times10^{-4}\,\text{\AA/ps}$$

versus typical ion speeds $\sim 10\,\text{\AA/ps}$ — four to five orders
of magnitude smaller, **dynamically null on the simulation timescale**.

*Why kept as secondary not discarded:* this is the only choice that is
formally physically principled (a real bath at a real temperature
satisfying real FDT). If the question "what does *thermal* noise do
here" ever needs answering, this is the only option that answers it.
Acknowledge upfront that it is expected to be invisible in the
output — it's a *correctness* anchor, not a *behaviour* anchor.

#### Discarded — none.

All three calibration purposes are physically meaningful. Hardware
cost of keeping (a) and (c) pluggable is small (different value of
one parameter); no reason to remove them.

#### Interchangeability surface

`SimConfig.noise_calibration ∈ {hard_sphere_variance, tddft_residual,
strict_fdt_bath}` with the resulting $T_\text{eff}$ (or equivalent
amplitude) stored alongside. The calibration itself is a
pre-processing step that produces a number; the simulation does not
need to know which procedure produced it.

---

### 1.4 Geometric structure (transverse noise)

The drag force acts along $-\hat v$. The noise can act along $\hat v$
only, isotropically in 3D, or anisotropically with separate
amplitudes parallel and perpendicular to $\hat v$.

The TDDFT extraction provides longitudinal information only —
trajectories are effectively 1D along the Coulomb dissociation axis.
Transverse noise is *unanchored* by the extraction pipeline; whatever
is chosen needs an external calibration source.

Strategy: **start minimal, escalate if needed.**

#### Primary — T1: longitudinal noise only

$$\vec{\mathcal{N}}(v,t) = \mathcal{N}_\parallel(v,t)\,\hat v$$

Only the speed fluctuates; the direction is determined entirely by
deterministic forces (drag + Coulomb + droplet).

*Why primary:*

- Matches the dimensionality of the extracted drag law (1D in →
  1D out).
- Zero free parameters beyond the longitudinal amplitude already
  decided in §1.2-§1.3.
- Minimum viable noise model. Easiest to validate.
- The hard-sphere ion–He mass ratio is $\rho = m_I/m_\text{He}
  \approx 32$, so per-collision lateral deflection is small. Whether
  *cumulative* transverse spread matters is an empirical question
  that T1 will answer by failing or not failing the VMI angular
  distribution validation.

*Validation criterion for escalation:* compare the final-velocity
angular distribution against the experimental VMI references
(`vmi_iplus_he.csv`, `vmi_iplus_gas.csv`). If T1 fits, T2/T3 are
unnecessary. If T1 underestimates angular spread, escalate to T3.

#### Secondary — T3: anisotropic with separate transverse amplitude

$$\vec{\mathcal{N}}(v,t) = \mathcal{N}_\parallel\,\hat v + \mathcal{N}_{\perp,1}\,\hat e_1 + \mathcal{N}_{\perp,2}\,\hat e_2$$

with $\hat e_1, \hat e_2$ spanning the plane perpendicular to $\hat v$,
and $\sigma_\perp \neq \sigma_\parallel$ in general.

*Calibration target for $\sigma_\perp$:* **VMI experimental angular
distribution**, not hard-sphere runs. Reasoning: $\sigma_\parallel$ is
already FDT-tied to the drag (§1.2), and the drag is anchored to
TDDFT. So $\sigma_\parallel$ has TDDFT provenance. $\sigma_\perp$ has
no internal anchor and should be tied to the most authoritative
external dataset available, which is the experimental VMI. Calibrating
$\sigma_\perp$ against hard-sphere runs would compound model error.

*Open issue — basis at $v\to 0$:* the basis $(\hat v, \hat e_1,
\hat e_2)$ is ill-defined at $v = 0$.

- **Primary low-$v$ behaviour:** $\sigma_\perp \to 0$ smoothly as
  $v \to 0$. Concretely, gate $\sigma_\perp$ by some smooth function
  of $v$ that vanishes at the origin. Loses the transverse channel at
  rest but is unambiguous and continuous.
- **Alternative to test:** smoothly blend to isotropic
  ($\sigma_\perp \to \sigma_\parallel$ as $v \to 0$, with the basis
  becoming irrelevant in the isotropic limit). Avoids losing
  stochasticity at rest but introduces a velocity-dependent
  isotropisation knob.

Both should be testable behind the same flag; pick "vanish" as the
first try.

*Why secondary not primary:* introduces a second amplitude parameter
that requires its own calibration source. Worth the cost only if T1
demonstrably fails the VMI angular validation.

#### Secondary — T2: isotropic 3D noise

$$\vec{\mathcal{N}}(v,t) = (\mathcal{N}_x, \mathcal{N}_y, \mathcal{N}_z), \quad \langle\mathcal{N}_i\mathcal{N}_j\rangle \propto \delta_{ij}$$

Same amplitude in all three Cartesian directions, independent
components. The drag still acts along $\hat v$, but the noise has no
privileged direction.

*Why secondary not primary:* doesn't match the hard-sphere structure
(which is privileged along $\hat v_\text{in}$ per collision). At
$v = 0$ it gives transverse kicks where hard-sphere gives nothing
defined. It also breaks the natural FDT pairing between drag (along
$\hat v$) and noise (also along $\hat v$).

*Why kept as backup:* trivially compatible with BAOAB integrators
(additive isotropic noise is the textbook case). If T3's basis
regularisation at $v\to 0$ turns out to be a numerical headache, T2 is
the easy escape valve. Conceptually clean; one parameter.

#### Discarded — none.

T1, T2, T3 are all physically defensible. Keep all three behind a
single enum.

#### Interchangeability surface

`SimConfig.noise_geometry ∈ {longitudinal, isotropic, anisotropic}`
with amplitude-related fields downstream of that choice. For
`anisotropic`, an additional field
`SimConfig.noise_low_v_behavior ∈ {vanish, blend_to_isotropic}`
controls the $v\to 0$ regulariser.

---

### 1.5 Summary of Section 1 choices

| Sub-decision | Primary | Secondary | Discarded |
|---|---|---|---|
| Functional form | N2: multiplicative, local FDT | N3: empirical residual | N1: additive at $\gamma_0$ |
| Calibration purpose | (b): hard-sphere variance | (c): TDDFT residual; (a): strict FDT at 0.4 K | — |
| Geometric structure | T1: longitudinal only | T3: anisotropic, then T2: isotropic | — |
| Low-$v$ behaviour (T3 only) | vanish | blend to isotropic | — |

### 1.6 Validation criterion specific to noise

The noise model is the only piece of the port whose effect lives
entirely in **second moments** of the output. The mean trajectory is
set by the drag; the noise only changes the spread around it.
Validation must therefore exercise *ensemble* statistics:

- Final-velocity histogram width against `vmi_iplus_he.csv`.
- Angular distribution of final velocities (if T3) against the same.
- Cross-trajectory variance at fixed time within an ensemble run on
  `single_pulse_droplet_distribution` (8000 atoms, variable droplet).

The existing `compare_distance` / `compare_velocity_magnitude` checks
look at single-trajectory matching and will **not** discriminate
noise-model variants. A new validation hook will be required.

---

## 2. Mass attachment

### 2.1 Physical question

What happens to the helium solvation shell during and after the Coulomb
explosion of I₂⁺ inside the droplet? This is genuinely open physics, not
just a modelling choice. Three things are known experimentally:

- **Pre-ionisation:** neutral I₂ in the droplet ground state carries
  ~42 He in a structured solvation shell — i.e. **~21 He per iodine
  atom**. This is a measured configuration of the *molecule*, not the
  ion-stage initial condition.
- **Mid-flight, 9 Å case (time-resolved, TDDFT-derived):** the per-I⁺
  shell *declines* monotonically along the trajectory — **~21 He
  pre-explosion → ~19 He at 10 ps → ~14 He at 14 ps**. The
  drag-extraction reference mass is taken at **~19 He ≈ 203 amu**, the
  window-representative value of this trajectory; this is the constant
  $m_\text{eff}$ the drag law is fit against (§2.2). It sits near the
  mean of the declining shell rather than at either extreme, so the
  drag law is calibrated at a mass the ion genuinely carries through the
  middle of the fit window.
- **Detector, post-ejection:** the I⁺(He)ₙ size distribution falls off
  monotonically — most ions are bare, then $n=1$, then $n=2$, etc.
  This is the *terminal* mass constraint.

The unknown is everything between explosion and mid-flight. Two
hypotheses bracket the possibilities, with a third occupying the middle:

- **Scenario A — full strip:** the Coulomb explosion at $R = R_e$
  blows the entire shell off. Bare I⁺ emerges into the droplet and
  rebuilds the shell during traversal.
- **Scenario B — coherent shell, surface stripping:** the shell
  largely survives the explosion (the dressed ion accelerates as a
  unit), with shedding concentrated at the droplet surface.
- **Biphasic — partial impulse strip, then equilibration:** part of
  the shell is lost in the initial impulse, then mass relaxes toward
  a conditions-dependent equilibrium with both gain and loss channels.

The Coulomb energy released per atom (~2.7 eV) is roughly 2× the total
shell binding (~1.3 eV for 42 He at ~30 meV each), placing the
situation in the borderline between impulsive and adiabatic regimes.
Neither A nor B can be ruled out a priori. The architecture must
support all three so they can be compared against the experimental
size distribution.

### 2.2 The $m_\text{eff}$ framing

A clarifying observation that makes the rest of this section
coherent: **$m_\text{eff} \approx 203\,\text{amu}$ (~19 He) is not "the
true mass of the ion." It is the mass that was *assumed* during drag
extraction.** It is a parameter of the *drag law*, not of the ion.

Consequences:

- The simulation's instantaneous $m(t)$ is the *physical* mass — what
  the dressed ion currently weighs.
- $m_\text{eff}$ is a property of the fitted force, like a reference
  state in a transport coefficient.
- When $m(t) \neq m_\text{eff}$, the drag law is being *extrapolated*
  outside its calibration domain. This is unrigorous but unavoidable
  (the extraction provides no $M$-scaling for $\gamma$).

This framing dictates the architectural decisions:

- **Integrator inertia:** uses $m(t)$ always. $F = m(t)\,a$. No
  ambiguity.
- **Drag force:** applied as $F_\text{drag} = -\gamma(v)\,v$
  regardless of $m(t)$. The extracted $\gamma(v)$ is taken at face
  value, with no invented $M$-dependence.
- **Kinetic energy bookkeeping:** uses $m(t)$ always.
- **`SimConfig.m_eff_amu` exists** as the drag-law reference mass,
  exposed and named separately from the physical mass trajectory.
  This makes the distinction visible and supports future re-
  extractions at different reference masses.

**Provenance and optional refinement (see §6.5–§6.6).** The extraction
mass is ~19 He ≈ 203 amu, the window-representative value of the
declining ~21→~14 He shell (§2.1) — i.e. a sensible constant near the
middle of the fit window, not an off-target guess. The only residual
unrigour is that it is a *constant* standing in for a quantity that
varies by ~⅓ across the trajectory; the $m(t)\neq m_\text{eff}$
extrapolation above is therefore mild near mid-window and grows toward
the ends. §6.6 records the optional time-resolved-$m(t)$ re-extraction
that removes even this residual; §6.5 makes the mass-scenario↔coefficient
pairing an enforced config-load guard and requires the extraction to
stamp each coefficient bundle with its `extraction_mass_model` metadata.

### 2.3 Scenario A — full strip with continuous accretion

- **Initial mass:** ~127 amu (bare I⁺).
- **Mass evolution:** monotone increase via continuous accretion ODE.
- **Functional form (primary velocity scaling):**
  $$\dot M(t) = \kappa_0 \cdot \rho_\text{He}(\text{depth}(t))$$
  density-driven, velocity-independent. He attaches at a rate
  proportional to local availability.
- **Alternative velocity scalings (secondary, interchangeable):**
  - $\dot M \propto v\,\rho_\text{He}$ — sweeping geometry; fast ion
    encounters more He per unit time.
  - $\dot M \propto \rho_\text{He}/v$ — dwell-time picture; slow ion
    lingers per He. This is what the hard-sphere model *accidentally*
    implements via the $\sigma \propto v^{-2}$ scaling.
- **Calibration target:** $\kappa_0$ fit jointly against (i) the
  time-resolved per-I⁺ shell trajectory at 9 Å (~21→~19→~14 He, §2.1),
  and (ii) the experimental I⁺(He)ₙ size distribution stored in the
  local Python history.
- **Energy bookkeeping:** kinetic-energy defect becomes a continuous
  integral $E_\text{defect}(t) = -\tfrac{1}{2}\int_0^t \dot M\,v^2\,dt'$.
  Always ≤ 0 (cold He arrives at $v$, fake KE injected, subtracted
  out).

**Known limitation, recorded explicitly.** Under Scenario A, the
simulation propagates an ion of mass ~127 amu through the bubble-exit
transient (~0.5 ps), applying a drag law calibrated at $m_\text{eff}
\approx 203$ amu. The drag law is being *extrapolated outside its
calibration domain at $t = 0$*. This is precisely the window the
drag extraction excluded via $t^*$, so no calibration data exists
inside it. The transient nevertheless determines the initial
conditions for the well-calibrated mid-flight phase. **Validation
tolerances on $v(t)$ during the first ~0.5 ps should not be tight
for Scenario A specifically.** Mid- and late-trajectory observables
remain the validation anchor.

### 2.4 Scenario B — coherent shell, surface stripping

- **Initial mass:** the measured ion-stage-onset shell, **~21 He per
  I⁺ ≈ 211 amu** (§2.1), *not* $m_\text{eff}$. The 42-He number is the
  whole-molecule pre-ionisation shell (~21 per atom); under B's
  coherent-shell premise the shell largely survives the explosion, so
  the ion stage begins near the measured ~21-He value and *declines*
  from there (~19 at 10 ps, ~14 at 14 ps). This initial physical mass
  (~211 amu) is distinct from the drag-law reference $m_\text{eff}
  \approx 203$ amu (~19 He, §2.2); per §2.2 the integrator uses the
  physical mass while the drag is applied at face value.
- **Mass evolution:** monotone decrease via continuous loss ODE.
  $\dot M < 0$ throughout.
- **Functional form (primary):**
  $$\dot M(t) = -\eta_0 \cdot \rho_\text{He}(\text{depth}(t))$$
  loss rate driven by *gradient* or *position* relative to droplet
  surface — concretely, mass loss concentrated where He density is
  changing rapidly (the bubble wall). Possible alternative forms
  include $\dot M \propto -|a|$ (loss driven by acceleration
  magnitude) — kept as a secondary trial.
- **Calibration target:** $\eta_0$ fit against the experimental
  size distribution.
- **Drag-law extrapolation:** B starts at ~211 amu, just *above* the
  $m_\text{eff} \approx 203$ amu the drag law was extracted at — only
  ~4% off, since both sit near the top of the measured shell range.
  Far less than Scenario A's bare-127-amu start (~37% light). As the
  shell strips, $m(t)$ falls *through* $m_\text{eff}$ toward lighter
  values; the extrapolation is near-zero in the first part of the
  trajectory and grows only as the shell drops well below 19 He late.
  The most benign drag-law-extrapolation profile of the evolving
  scenarios.
- **Energy bookkeeping:** the kinetic-energy "defect" term changes
  sign and meaning. With $\dot M < 0$:
  $$E_\text{transfer}(t) = -\tfrac{1}{2}\int_0^t \dot M\,v^2\,dt' \geq 0$$
  represents *real energy carried away* by shed He, not a correction
  for fake injection. The name `E_mass_attach_defect_eV` no longer
  fits; rename to `E_mass_transfer_eV` to cover both signs cleanly.
  Schema bump to `IonCheckpoint` v6.

### 2.5 Biphasic — partial strip with equilibration

- **Initial mass:** the measured ion-stage-onset shell, **~21 He per
  I⁺ ≈ 211 amu** (§2.1), same as Scenario B — the 42-He figure is the
  whole-molecule pre-ionisation shell (~21 per atom), and the ion stage
  begins with the measured surviving shell. The biphasic model takes
  this as the starting point and allows both gain and loss channels.
  As in §2.4, this physical initial mass is distinct from the drag-law
  reference $m_\text{eff} \approx 203$ amu (§2.2).
- **Mass evolution:** relaxation toward a conditions-dependent
  equilibrium:
  $$\dot M(t) = \frac{M_\text{eq}(v, \text{depth}) - M(t)}{\tau}$$
  with $M_\text{eq}(\cdot)$ a function and $\tau$ a relaxation time.
- **Functional form choices for $M_\text{eq}$:** undetermined
  without further data. Simplest first try: $M_\text{eq}$ depends on
  local He density only, with the equilibrium shell being heavier in
  denser surroundings. Velocity dependence would represent shell
  stripping by ion speed (fast ion can't hold its shell).
- **Calibration:** at least two parameters ($\tau$ and the
  scaling of $M_\text{eq}$). Cannot be pinned down by mid-flight
  mass alone; needs the experimental size distribution *and* the
  time-resolved shell trajectory (~21→~19→~14 He, §2.1) as separate
  constraints.
- **Drag-law extrapolation:** least problematic of the three — the
  trajectory $m(t)$ is closest to $m_\text{eff}$ on average if the
  equilibrium is set near 203 amu.
- **Mass is not monotone.** Implications for downstream tooling that
  may assume `mass_history_kg` is non-decreasing. Documentation flag.

**Why biphasic is interesting but not primary.** Most parameters,
most degrees of freedom, hardest to calibrate from available data. It
is the most flexible model and likely the most physically realistic,
but starts as a secondary option until A and B have been compared.

### 2.6 Sub-decision selection

**Primary — Scenario A, density-driven accretion.**

The simplest hypothesis that matches the "violent explosion strips
everything" intuition. Single parameter $\kappa_0$ to calibrate
against two constraints (mid-flight mass, terminal size distribution).
The known drag-law extrapolation issue during the bubble-exit
transient is explicitly recorded as a tolerance relaxation, not a
fatal flaw.

**Secondary — Scenario B, density-driven stripping.**

Symmetric to A in structure (one parameter, same calibration targets)
but with opposite sign on $\dot M$. Worth implementing in parallel so
A and B can be compared on the same validation surface. Weak
empirical hint in its favour: the existing model's factor-of-18 drop
in attachment rate between 9 Å and 18 Å droplets (0.09 → 0.005) is
more naturally explained by *net loss being faster in smaller
droplets* than by *gross accretion being faster in smaller droplets*.

**Secondary — Biphasic relaxation.**

Held in reserve for the case where neither A nor B individually fits
the experimental size distribution adequately. More parameters, more
flexibility, harder to constrain.

**Tertiary / null hypothesis — Fixed mass at $m_\text{eff}$.**

Pin $m(t) = m_\text{eff} \approx 203$ amu throughout the ion stage. No
mass dynamics. Useful as a *baseline*: if M1 matches the validation
surface adequately, mass dynamics don't matter and A/B/biphasic are
over-modelling. If M1 fails, the comparison tells us *how much* mass
dynamics matters. Also: M1 is the only model where the drag law is
never extrapolated outside its calibration domain — it is the
internal-consistency reference. §6.5 sharpens this: M1 is the *only*
scenario for which the in-hand constant-mass coefficients form a
self-consistent pair, and the §6.5 config-load guard enforces that
pairing.

**Discarded — collision-gated Bernoulli attachment (current model).**

Cannot survive without collisions to gate on. The physical content
that was hidden in its $\sigma(v)$ coupling (slow ions pick up faster
via the implicit $1/v$ scaling) gets re-expressed cleanly in the
Scenario A velocity-scaling alternatives.

**Empirical lean toward Scenario B (recorded; see §6.8).** The TDDFT
shell counts (~21 He pre-explosion → ~19 at 10 ps → ~14 at 14 ps,
9 Å case) show monotone *loss*, pointing at Scenario B over the
primary Scenario A. B is additionally the only scenario whose
time-resolved $m(t)$ can be read directly off the reference rather than
modelled. This strengthens the 0.09→0.005 attach-rate hint above with a
second independent argument. **Not auto-promoting:** A remains the
recorded primary and simplest hypothesis until the §6.4 hierarchy
(Tier 1/Tier 2) adjudicates A against B.

### 2.7 Sub-decision: velocity scaling within Scenario A

Independent of the A/B/biphasic choice, *within* Scenario A the
velocity scaling of $\dot M$ encodes a real physical claim.

- **Primary:** $\dot M \propto \rho_\text{He}(\text{depth})$, no $v$
  dependence. Simplest. He attaches based on local availability;
  speed irrelevant.
- **Secondary:** $\dot M \propto v\,\rho_\text{He}$ — sweeping
  geometry. Fast ions encounter more He per unit time.
- **Secondary:** $\dot M \propto \rho_\text{He}/v$ — dwell-time. Slow
  ions linger and stick. This is what hard-sphere accidentally
  implements; preserving it is a hedge against losing some emergent
  behaviour from the old model.
- **Discarded:** none. All three are physically defensible and the
  data does not yet discriminate.

### 2.8 Interchangeability surfaces

The architecture exposes mass dynamics through the following surfaces:

- `SimConfig.mass_scenario ∈ {fixed, scenario_A_accretion,
  scenario_B_stripping, biphasic}` — top-level scenario selector.
- `SimConfig.mass_initial_amu` — initial ion-stage mass; defaults to
  ~127 amu (bare I⁺) under `scenario_A_accretion`, to the measured
  onset shell ~211 amu (~21 He, §2.1) under `scenario_B_stripping` and
  `biphasic`, and to $m_\text{eff}$ under `fixed` (the only scenario
  whose initial mass is the drag-law reference, per the §6.5
  consistency guard).
- `SimConfig.m_eff_amu` — drag-law reference mass, default ~203 (19 He).
  Used by the drag extraction and exposed for transparency; **not**
  the same as `mass_initial_amu`.
- `SimConfig.mass_rate_form ∈ {density_only, sweeping, dwell_time}` —
  velocity scaling for $\dot M$ under Scenarios A and B.
- `SimConfig.mass_rate_coefficient` — $\kappa_0$ or $\eta_0$, sign
  determined by scenario.
- For biphasic: `SimConfig.mass_relaxation_tau_ps` and a
  functional-form selector for $M_\text{eq}(\cdot)$.

These are decoupled from the noise-model surfaces (§1) and from the
drag-functional-form surfaces (§3, future).

### 2.9 Schema and energy-bookkeeping changes

**`IonCheckpoint` schema bump to v6** is required under any scenario
other than `fixed`. Changes:

- Rename `E_mass_attach_defect_eV` → `E_mass_transfer_eV`. Same array
  shape `(2N, T)`, same cumulative semantics, but the sign convention
  now covers both accretion (negative cumulative) and stripping
  (positive cumulative).
- `mass_history_kg` retains its shape but its monotonicity
  guarantee is dropped (biphasic and Scenario B both violate
  monotonicity).
- Add a small metadata field recording which scenario produced the
  run, so downstream tools can interpret `E_mass_transfer_eV` and
  `mass_history_kg` correctly.

**Energy invariant** under continuous mass dynamics:
$$E_\text{kin}(t) + E_\text{pot}(t) + E_\text{dissip}(t) + E_\text{mass\_transfer}(t) \approx \text{const}$$
modulo Verlet drift. The continuous form replaces the per-step
discrete defect from the existing model. No new dissipation channels;
just a smooth version of the existing accounting.

### 2.10 Validation criterion specific to mass attachment

Unlike the noise model (§1.6), mass attachment affects both **first
and second moments**:

- First moment: the trajectory $v(t)$ is modulated by $m(t)$ via
  $F = m\,a$. Comparison against the TDDFT reference distance and
  velocity traces (post-transient phase) is sensitive to mass
  evolution.
- Distribution: the final I⁺(He)ₙ size distribution at the detector
  is *the* discriminating experimental observable between scenarios.
  This is the new validation target the existing code does not yet
  exercise.

The local Python history of experimental size-distribution data needs
to be plumbed into the validation surface as a new comparison routine.

---

## 3. Drag functional form — analytic vs. tabulated

### 3.1 Physical question

How is $F_\text{drag}(v)$ represented in the integrator? The extraction
pipeline (`drag_calculation.py`) produces both a closed-form analytic
fit and the raw force-balance scatter $(v, F_\text{drag})$ on the
trusted interior, so three representations are available at zero extra
extraction cost: the power-law fit (`fit_variant=1`), the
linear+cubic fit (`fit_variant=2`), and a tabulated/interpolated form
sampled directly from the scatter. Both analytic variants have already
been extracted for both the 9 Å and 18 Å cases.

The choice matters because the drag form is consumed in three distinct
places, each with different requirements:

- the integrator acceleration term (§4) — wants cheap evaluation;
- the N2 noise amplitude (§1.2), which needs an analytic
  $\gamma(v) = |F_\text{drag}(v)| / v$ — wants a closed form,
  ideally finite at $v=0$;
- the low-velocity regulariser (§3.8) — only triggered by a form that is
  singular at $v=0$.

### 3.2 Dimensional analysis

All candidate forms must produce a force in
$\text{amu}\cdot\text{Å/ps}^2$ from a speed in $\text{Å/ps}$. Each
coefficient must absorb whatever power of velocity it multiplies.

**Power law** $\;|F_\text{drag}| = \gamma\,v^{\,n}$:
$$[\gamma]\,(\text{Å/ps})^{n} = \text{amu}\cdot\text{Å/ps}^2
\;\Rightarrow\; [\gamma] = \text{amu}\cdot\text{Å}^{\,1-n}\cdot\text{ps}^{\,n-2}.$$
For the extracted $n = -2$: $[\gamma] = \text{amu}\cdot\text{Å}^{3}
\cdot\text{ps}^{-4}$. Balances. Singular at $v\to 0$ ($n<0
\Rightarrow F\to\infty$).

**Linear + cubic** $\;F_\text{drag} = a\,v + b\,v^3$:
$$[a]\,\text{Å/ps} = \text{amu}\cdot\text{Å/ps}^2
\;\Rightarrow\; [a] = \text{amu/ps},$$
$$[b]\,(\text{Å/ps})^3 = \text{amu}\cdot\text{Å/ps}^2
\;\Rightarrow\; [b] = \text{amu}\cdot\text{ps}\cdot\text{Å}^{-2}.$$
Matches the extraction doc's stated units ($a$ in amu/ps, $b$ in
amu·ps/Å²). Balances. Regular at $v=0$ ($F\to 0$).

**Linear + quadratic** $\;F_\text{drag} = a\,v + c\,v\,|v|$ (the
$v|v|$ keeps the term odd, so it is sign-correct/dissipative for both
signs of $v$):
$$[c]\,(\text{Å/ps})^2 = \text{amu}\cdot\text{Å/ps}^2
\;\Rightarrow\; [c] = \text{amu}\cdot\text{Å}^{-1}.$$
Balances. Regular at $v=0$ ($F\to 0$). High-$v$ wing $\sim v^2$
(inertial / form drag) is physically gentler than the cubic.

**Threshold / saturating** $\;F_\text{drag} = F_\text{sat}\,
\tanh(v/v_0)$:
$$[F_\text{sat}] = \text{amu}\cdot\text{Å/ps}^2,\qquad
[v_0] = \text{Å/ps}.$$
$\tanh$ is dimensionless. Balances. Regular at $v=0$
($F\to F_\text{sat}\,v/v_0$, i.e. linear with effective
$a_\text{eff} = F_\text{sat}/v_0$ in amu/ps). Drag *saturates* at
$F_\text{sat}$ for $v\gg v_0$ — the only candidate whose force is
*bounded* at high $v$.

**Tabulated** $\;F_\text{drag} = \text{interp}(v;\{v_i, F_i\})$:
dimensionless interpolation over dimensioned samples; result carries
the units of the stored $F_i$, i.e. $\text{amu}\cdot\text{Å/ps}^2$.
Balances trivially. Behaviour at $v\to 0$ and beyond $v_\text{max}$ is
an extrapolation-policy choice, not fixed by the form.

No formulation is rejected on dimensional grounds.

### 3.3 Physical trade-offs

**Analytic vs. tabulated — degrees of freedom.** Going analytic →
tabulated *gains* fidelity to the extracted shape inside the
extraction window but *loses* a clean analytic $\gamma(v)$ for the FDT
noise and a clean parametric uncertainty band. The extraction window
is bounded in velocity ($[t^*, t_\text{end}]$ maps to a finite speed
range), so a table must extrapolate at both ends — relocating the
arbitrariness from "functional form" to "extrapolation rule" rather
than removing it. A table also forces the noise amplitude to
re-differentiate $\gamma(v)$ numerically, reintroducing precisely the
noise-amplification problem the extraction spline exists to suppress.

**Power-law vs. linear+cubic — degrees of freedom.** Power-law → l+c
*gains* $v=0$ regularity and a finite $\gamma_0 = a$ (which
re-enables a finite-amplitude thermal floor and resurrects the
physical content of the discarded N1 noise as a low-$v$ limit), and
*gains* an analytic derivative $dF/dv = a + 3bv^2$ for §4 operator
splitting and any implicit-step Jacobian. It *loses* the single-power
simplicity and the direct match to the extraction's historical default
(`fit_variant=1`).

**The form is an extrapolation, not an interpolation — and no
ion-in-superfluid drag law is known.** The extraction window is bounded
in velocity, but the ion operates (and especially the bubble-exit
transient, §2.3) *outside* that window. No literature law exists for the
velocity dependence of drag on an ion travelling through superfluid
helium, so every analytic form here is a *hypothesis* to be
cross-checked against the TDDFT curves, not a fitted truth. The four
analytic forms are deliberately chosen to *disagree most in the
extrapolation regime* — their high-$v$ wings span saturating
($\tanh$), quadratic ($v|v|$), cubic ($v^3$), and decaying
($v^{-2}$ power law) behaviour. Comparing their simulation outputs is
the point of keeping them interchangeable; the form that best matches
the references *is* a simulation result, not an a-priori choice. This
reframes the per-form discussion below: "primary" means "first
hypothesis to run," not "believed correct."

**Asymptotic signatures (for orientation, not selection):** linear
Stokes-like drag at low $v$; an inertial/form-drag $\sim v^2$ wing at
high $v$ is the generic expectation when an object sheds fluid; a
critical-velocity (Landau) threshold is the superfluid-specific
possibility hinted at by the baseline's existing `v_limit`/`E_min`
cutoff. The power-law $n\approx-2$ is most likely a fitting artifact of
the old hard-sphere $\sigma\propto v^{-2}$ cross-section bleeding into
the extraction window rather than a transport law — retained as a
hypothesis precisely so the cross-check can confirm or reject that.

**Dissipativity guard (physical validity).** For $F_\text{drag}$ to
oppose motion at all relevant speeds, each form has a validity
condition checked at config-load (the extraction's `curve_fit` carries
no bounds):

- **linear+cubic:** requires $a + b\,v^2 \ge 0$ across the operating
  range. $a>0, b\ge0$ is monotone and clean. If $b<0$ there is a
  turnover speed $v_\dagger = \sqrt{-a/b}$ beyond which the cubic flips
  the force sign (anti-dissipative); assert $a>0$ and $v_\dagger$ above
  the maximum trajectory speed.
- **linear+quadratic:** requires $a + c\,|v| \ge 0$; with $a,c>0$ this
  holds for all $v$ with no turnover. Assert $a>0, c\ge0$.
- **threshold:** $F_\text{sat}/v_0 > 0$ (i.e. both $F_\text{sat}>0$ and
  $v_0>0$) guarantees dissipative and monotone everywhere; bounded
  force means no high-$v$ stiffness.
- **power law:** $\gamma>0$; the $v=0$ singularity for $n<0$ is deferred
  to the §3.8 regularisation note rather than guarded here.

The need for a "maximum trajectory speed" to test $v_\dagger$ against
is itself an open item — it has no home in the current baseline config
and must be sourced from the TDDFT trajectories or set as a generous
ceiling (flagged for §6 validation/tolerances).

### 3.4 The analytic form set — four swappable hypotheses

All four analytic forms are first-class and interchangeable behind the
`drag_form` enum (§3.8). Because no ion-in-superfluid drag law is known
(§3.3), "primary" denotes the *first hypothesis to run and the default*,
not a belief about correctness; the empirical cross-check against the
TDDFT references selects among them.

**Primary — linear+cubic** $\;F_\text{drag} = a\,v + b\,v^3$. Regular
at $v=0$, finite $\gamma_0 = a$ for the N2 noise, analytically
differentiable ($dF/dv = a + 3bv^2$). Coefficients $\{a,b\}$ (with
$1\sigma$ errors and $R^2$) already extracted for both 9 Å and 18 Å,
so it is the form with data in hand and is kept as the default first
run. Caveats recorded honestly: the $v^3$ wing extrapolates more
steeply than the physically-generic quadratic, and it carries the
$b<0$ turnover risk (guarded, §3.3).

*Note — linear+quadratic is the stronger physical default and may
supersede this on cross-check.* On both integrator-neutrality and
asymptotics grounds, linear+quadratic is arguably the better primary:
its high-$v$ wing is the generic inertial $\sim v^2$ rather than the
unmotivated $v^3$, it has no turnover risk, and its $\gamma(v)$ grows
only linearly (milder stiffness). It is kept *secondary* here only
because linear+cubic is the form already extracted and chosen; if the
empirical cross-check does not favour linear+cubic, linear+quadratic is
the first alternative to promote. Flagged so the choice is revisited
with data rather than left implicit.

**Secondary — linear+quadratic** $\;F_\text{drag} = a\,v + c\,v\,|v|$.
Regular and dissipative at $v=0$ with no turnover (monotone for
$a,c>0$), inertial $\sim v^2$ high-$v$ wing, $\gamma(v) = a +
c\,|v|$ growing only linearly. Finite $\gamma_0 = a$. Needs a fit
pass (not yet extracted). See the promotion note above.

**Secondary — threshold/saturating** $\;F_\text{drag} =
F_\text{sat}\tanh(v/v_0)$. The superfluid-physics-anchored hypothesis:
linear at low $v$ (effective $\gamma_0 = F_\text{sat}/v_0$ in amu/ps),
*bounded* force at high $v$, and a built-in velocity scale $v_0$ that
can represent a critical-velocity-like onset. $\gamma(v)$ *decays* at
high $v$ — essentially zero stiffness in the operating regime, the
friendliest of all forms for any integrator. Cost: an extra parameter
and harder calibration from two trajectories. Not yet extracted; needs
a fit pass.

**Secondary — power law** $\;|F_\text{drag}| = \gamma\,v^{\,n}$. The
extraction's historical default form and the "drag falls with speed"
hypothesis. Most likely a hard-sphere $\sigma\propto v^{-2}$ artifact
(§3.3); retained so the cross-check can confirm or reject that.
Singular at $v=0$ for $n<0$, which re-enables the §3.8 low-velocity floor and
the matching noise regulariser (§1.2) — the only form that does so.
Coefficients $\{\gamma,n\}$ already extracted for both cases.

Since all four are closed-form, swapping among them is a change behind
the form enum with no structural impact on the integrator (the chosen
BAOAB integrator is form-agnostic and unconditionally stable for all
four — see §4).

*Consequence for regularisation.* Only `power_law` ($n<0$) needs a
low-velocity floor; the other three forms are finite at $v=0$ and need
none. This obligation lives in the §3.8 regularisation note
(`drag_low_v_floor`), not a standalone section.

*Consequence for §1.* The N1-discard reasoning in §1.2 is rescoped: the
"$\gamma_0\to\infty$" objection applies only under `power_law`; under
the other three, $\gamma_0$ is finite and N2 reduces to the correct N1
limit at low $v$.

### 3.5 Discarded — tabulated $F(v)$

Discarded for three reasons, recorded so the option is not
re-litigated:

1. The N2 noise amplitude needs an analytic $\gamma(v)$; a table forces
   numerical re-differentiation and reintroduces the noise-amplification
   the extraction spline was built to avoid.
2. The table is faithful only inside the bounded extraction velocity
   window; outside it the representation must extrapolate, merely
   relocating the arbitrariness rather than removing it.
3. No clean parametric uncertainty band (the 10-seed extraction sweep
   maps directly onto $\{\gamma,n\}$ or $\{a,b\}$, not onto a table).

### 3.6 Per-case form policy — shared-form preferred

The two bubble sizes *may* sit in different drag regimes (the 18 Å
droplet has a longer dense-He traversal and, in the old model, a
different cross-section and binding). The architecture therefore
allows the drag *form*, not just its coefficients, to differ per case:
`drag_form` is a per-preset value.

The *default policy*, however, is shared-form-preferred:

- Use the primary form (`linear_cubic`) for both cases. Report the
  shared-form fit as the headline result.
- A shared functional form fitting both cases with only the
  coefficients moving is weak evidence the law captures real transport
  physics rather than per-case curve-fitting — this cross-case
  consistency signal is worth preserving as the anchor result.
- Diverge to a different form for one case *only* if it trips the
  §3.3 dissipativity or an $R^2$ acceptance guard. Such a divergence is
  itself a finding (the cases are in different regimes) and is recorded
  with its reason rather than chosen by default.

Fully independent per-case selection is *not* adopted as the default:
it would bake in the "different regimes" assumption before testing it
and forfeit the consistency signal, while gaining nothing the
guard-triggered fallback does not already provide. The architecture
supports it (the field is per-preset either way); only the default
policy is constrained.

*Trade-off flagged:* per-case forms maximise per-case fit quality at
the cost of the one global consistency check. Shared-form-preferred
keeps that check as the headline and treats divergence as a documented
exception.

### 3.7 Extraction-side consumption note

The primary form is `linear_cubic` (extraction `fit_variant=2`), but
`Drag_extraction_code.md` documents `fit_variant=1` (power law) as the
extraction's *default* output. Both variants have now been extracted
for both cases, so no pipeline change is required — but the MD side must
**consume the `fit_variant=2` (linear+cubic) coefficients as the
default**, and the extraction's default-output setting should be
aligned (or the linear+cubic product explicitly carried forward) so
that the primary drag form is the one actually plumbed into the
simulation. The `linear_quadratic` and `threshold` forms (§3.4) are
*not yet extracted* and require their own fit passes before they can be
cross-checked. Recorded as a cross-document consumption dependency plus
two outstanding fit passes, not a re-extraction of existing forms.

### 3.8 Interchangeability surface

- `SimConfig.drag_form ∈ {linear_cubic, linear_quadratic, power_law,
  threshold}` — default/primary `linear_cubic` (data in hand);
  `linear_quadratic` flagged as the likely physical supersedor;
  `threshold` the superfluid-anchored hypothesis; `power_law` is the
  only form requiring a low-velocity floor (see the regularisation note
  below). Per-preset, enabling per-case forms.
- `SimConfig.drag_coefficients` — form-tagged coefficient bundle of
  variable arity by form, per-case:
  - `linear_cubic`: $\{a, b\}$, units amu/ps and amu·ps/Å².
  - `linear_quadratic`: $\{a, c\}$, units amu/ps and amu/Å.
  - `power_law`: $\{\gamma, n\}$, $\gamma$ in
    amu·Å$^{1-n}$·ps$^{n-2}$, $n$ dimensionless.
  - `threshold`: $\{F_\text{sat}, v_0\}$, units amu·Å/ps² and Å/ps.
- **Drag-validity guard (config-load):** per-form dissipativity check
  from §3.3 — `linear_cubic` ($a>0$, $v_\dagger$ above max speed);
  `linear_quadratic` ($a>0, c\ge0$); `threshold` ($F_\text{sat}>0,
  v_0>0$); `power_law` ($\gamma>0$). The "maximum trajectory speed"
  needed for the `linear_cubic` turnover check is an open config item
  (flagged for §6 validation/tolerances).
- **Low-velocity regularisation (only `power_law`).** The three
  finite-at-zero forms (`linear_cubic`, `linear_quadratic`,
  `threshold`) need no regulariser — $F_\text{drag}\to0$ and
  $\gamma(v)$ stays finite as $v\to0$. Only `power_law` with $n<0$
  diverges; it then requires a floor on $\gamma(v)$ (equivalently a
  small $v_\text{floor}$ below which $\gamma$ is held constant) to keep
  both the drag and the FDT noise amplitude finite. The chosen BAOAB
  integrator (§4) partially defuses the divergence on its own — the
  damping factor $e^{-\gamma\,dt/m}\in[0,1]$ stays bounded even as
  $\gamma\to\infty$ — but the noise amplitude
  $\sqrt{2\gamma k_BT_\text{eff}}$ still needs the floor. Exposed as
  `SimConfig.drag_low_v_floor` (Å/ps), active only when
  `drag_form = power_law`; ignored otherwise. This absorbs what was
  previously a standalone low-velocity-regularisation section.
- Cross-references: $\gamma(v) = |F_\text{drag}|/v$ (units amu/ps)
  feeds the §1.2 N2 noise amplitude; $m$ here is the physical $m(t)$ per §2.2, and the
  drag is applied as $-\gamma(v)v$ independent of $m(t)$.

---

## 4. Integrator coupling — additive vs. operator-split

### 4.1 Physical question

How do the drag force (§3) and, once active, the Langevin noise (§1)
enter the time-stepping relative to the baseline velocity-Verlet
integrator? The baseline `velocity_verlet_step` (`leapfrog.py:74`) is
kick-drift-kick, second-order, and model-agnostic: it consumes any
`acc_fn(pos) → (acc, E_pot)`. The deeper issue is that this clean
structure assumes a *position-dependent* force $F = F(x)$, whereas drag
is $F_\text{drag}(v)$ and noise makes the system a stochastic
differential equation. Both break the assumptions velocity-Verlet
rests on.

### 4.2 Why drag breaks the velocity-Verlet assumption

Velocity-Verlet's second-order accuracy and clean error structure rely
on the force being a function of position only. A velocity-dependent
force has three consequences:

- The second kick $v_1 = v_0 + \tfrac12(a_0 + a_1)\,dt$ becomes
  *implicit*: $a_1 = a_1(x_1, v_1)$ depends on the $v_1$ being solved
  for. $v_1$ appears on both sides.
- Explicit (naïve) treatment of a dissipative $F(v)$ term is only
  *conditionally* stable: stable when $\gamma\,dt/m \ll 1$, divergent
  when drag is stiff. Stiffness is form-dependent (§3): for
  `linear_cubic`/`linear_quadratic`/`threshold` over the working range,
  $\gamma\,dt/m \sim 10^{-4}$ (very non-stiff); for `power_law` with
  $n=-2$, $\gamma\to\infty$ as $v\to0$ and the term becomes arbitrarily
  stiff at low speed — exactly where the ion ends up at late times.
- With §1 noise, the system is an SDE; standard Verlet has no noise
  term and cannot inject one with correct fluctuation–dissipation
  statistics at finite $dt$.

Time-reversibility is lost under *any* coupling — dissipation is
irreversible by construction. That is physics, not a defect of the
scheme.

### 4.3 Dimensional check

Drag enters the equation of motion as an acceleration
$a_\text{drag} = F_\text{drag}(v)/m(t)$:
$$[a_\text{drag}] = \frac{\text{amu}\cdot\text{Å/ps}^2}{\text{amu}}
= \text{Å/ps}^2,$$
matching the existing Coulomb and droplet accelerations. With §2.2,
$m(t)$ is the physical mass and the drag is applied as
$-\gamma(v)\,v/m(t)$ with $\gamma(v)$ in amu/ps. The OU damping factor
below, $e^{-\gamma\,dt/m}$, has a dimensionless exponent:
$[\gamma\,dt/m] = (\text{amu/ps})(\text{ps})/\text{amu} = 1$. Balances.

### 4.4 The option space

**Option I — naïve additive inside `_ion_accel_fn`.** Add
$-\gamma(v)v/m$ to the acceleration using whatever velocity is
available at each Verlet sub-evaluation; integrator structure
untouched. Minimal diff. But it treats the velocity-dependent term
explicitly, so stability is conditional and form-dependent (safe for
the finite-$\gamma$ forms over their working range, unstable for
`power_law` near $v=0$), and it loses second-order accuracy in the drag
part. Noise can only be bolted on as an explicit kick, which gets the
finite-$dt$ FDT statistics wrong.

**Option II — semi-implicit additive.** Same entry point, but solve the
implicit second kick. For the linear part this is closed-form; for the
cubic/quadratic part it is a cheap scalar Newton step per atom.
Unconditionally stable for dissipative drag, but the integrator must
now *know* drag is velocity-dependent and treat it specially (the clean
"drag is just another acceleration" abstraction is broken). Still no
natural exact-OU site for the noise.

**Option III — operator splitting (BAOAB family).** Split the evolution
operator into a conservative Verlet kick/drift (B/A, handling
Coulomb+droplet via the existing `acc_fn`) and an
Ornstein–Uhlenbeck step (O, handling drag+noise), composed symmetrically
(half-B, half-A, full-O, half-A, half-B). For *linear* drag the O-step
is exact:
$$v \mapsto e^{-\gamma\,dt/m}\,v
+ \sqrt{\frac{k_B T_\text{eff}}{m}\left(1 - e^{-2\gamma\,dt/m}\right)}\;\xi,$$
unconditionally stable, exact in the drag part for linear $\gamma$, and
the canonical site for injecting the §1 Langevin noise with correct FDT
statistics at finite $dt$. For the nonlinear forms, $\gamma(v)$ varies
within the step; freeze $\gamma$ at step-entry velocity (an $O(dt)$
error in the nonlinear correction, acceptable since it is a correction
on top of the linear floor) or do a short implicit sub-solve inside O.
The bounded damping factor $e^{-\gamma\,dt/m}\in[0,1]$ stays
well-behaved even as $\gamma\to\infty$, so the O-step *partially defuses*
the `power_law` stiffness (§3.8).

### 4.5 Trade-offs

*Structure.* I keeps the model-agnostic integrator and minimal diff but
loses unconditional stability, a principled noise site, and second-order
drag accuracy. II gains unconditional drag stability and keeps the
conservative forces in the clean callable but breaks the abstraction and
still has no exact-OU noise site. III gains unconditional stability,
exact linear-drag propagation, correct finite-$dt$ FDT noise, and
partial `power_law` stiffness relief, at the cost of a *new ion-stage
integrator path* and the Stratonovich/Itô bookkeeping (which BAOAB
resolves cleanly — see §1.2).

*Conservation.* The §2.9 energy invariant
$E_\text{kin}+E_\text{pot}+E_\text{dissip}+E_\text{mass\_transfer}$
needs per-step dissipated energy. Under I/II this is the drag work
$\int F_\text{drag}\cdot v\,dt \approx F_\text{drag}\cdot v\,dt$. Under
III the O-step yields the dissipated energy analytically from the
velocity damping, and separates it cleanly from the noise energy
injection (which feeds the thermal floor, not $E_\text{dissip}$) —
cleaner accounting.

*Coupling.* §1 noise: III is the only option with a native,
statistically-correct noise site; I/II require a separate stochastic
bolt-on with wrong finite-$dt$ statistics. §3 form: I's stability is
form-dependent; II/III are form-robust. §2 mass: $m(t)$ changes per
step, so $\gamma/m$ in the O-step uses the current $m(t)$ — consistent
with the per-step closure rebuild (`ion_propagation_step.py:184-193`),
no conflict. §3.8 floor: III reduces but does not eliminate the `power_law`
regulariser (a floor is still needed to keep the noise amplitude finite
and avoid $e^{-\infty}$ pathologies).

### 4.6 Sub-decision selection

**Primary — III, BAOAB operator splitting.** Chosen on the principle
that the integrator must be **form-agnostic and unconditionally stable**
so that no drag-form choice can destabilise it and corrupt the
empirical form cross-check (§3.3): "I changed the form and the
trajectory changed" must be a clean physics statement, not an artifact
of a shifting stability margin. III is unconditionally stable for all
four forms including `power_law`'s singularity (bounded exponential),
its O-step is exact for the linear floor every form shares, and it is
the only option with a native FDT-correct noise site — which the §1.2
noise decision already implicitly committed to (its repeated references
to BAOAB-family integrators). The cost (a separate ion-stage integrator
path, more implementation work) is accepted as the price of protecting
the comparison.

*Accepted architectural asymmetry.* The neutral stage keeps the
baseline `velocity_verlet_step`; only the ion stage gets the BAOAB
path. This is a deliberate fork — the neutral stage has no drag and no
Langevin noise, so there is nothing for BAOAB to buy there, and §14 of
the baseline keeps the neutral integrator off-limits regardless.

**Discarded — I, naïve additive.** Its only virtue is minimal diff,
which evaporates the moment noise is switched on (the default intent).
Form-dependent stability would let `power_law` runs diverge at low $v$,
contaminating the form cross-check with integration artifacts. Recorded
discarded so it is not revisited.

**Discarded — II, semi-implicit additive.** Unconditionally stable for
deterministic drag and simpler than standing up BAOAB, so it was a
candidate secondary for a *deterministic-drag-only* mode. Discarded
because the strict-FDT noise (§1.3 option (a)) is kept as a permanent
correctness anchor and is always at least nominally active, so a
truly noise-free mode is not a first-class use case; and because
maintaining two integrator paths for the ion stage (II and III) doubles
the validation surface for no behavioural gain that III does not already
cover. If a deterministic-drag debugging mode is ever needed, it is
recovered by setting the noise amplitude to zero within III, not by
re-introducing II.

### 4.7 Interchangeability surface

- The integrator choice is **not** exposed as a routine `SimConfig`
  swap. Unlike the force-form and noise enums, the integrator is a
  deeper structural object; III is hard-wired for the ion stage and
  I/II are recorded as discarded rather than kept pluggable. A
  convergence-testing hook may temporarily expose alternatives, but the
  production path is single.
- `SimConfig.integrator_dt` reuses the existing `dt_ion` (0.01 ps); the
  BAOAB step consumes the same timestep.
- The O-step reads `drag_form`/`drag_coefficients` (§3.8) to build
  $\gamma(v)$ and the noise surfaces (`noise_form`,
  `noise_calibration`, `noise_geometry`, §1) to build the fluctuation
  term. All physics choices remain behind their own enums; III is the
  fixed *mechanism* that consumes them.
- Per-step dissipated and noise-injected energies are written to the
  §2.9 accumulators (`E_dissip_eV`, and the thermal-floor channel),
  computed analytically from the O-step rather than by finite-difference
  drag work.

---

## 5. Spatial gating — does drag turn off at the droplet surface?

### 5.1 Physical question

The drag law (§3) was extracted from TDDFT trajectories with the ion
*inside* the droplet. But the ion's trajectory ends with ejection into
vacuum on its way to the detector, and in vacuum there is no helium to
drag against — drag (and the FDT noise that accompanies it) must go to
zero. The question is *how* drag transitions from full (deep inside,
dense He) to off (outside, vacuum), not *whether* it does.

**This decision has no viable null.** Unlike the noise, mass, and form
decisions — each of which had a defensible "do nothing" baseline —
leaving drag on outside the droplet is *physically wrong*: the ion
would keep decelerating in vacuum, corrupting precisely the final
detector velocity the model validates against (`vmi_iplus_he.csv`).
Some gating is mandatory; only its shape is open.

### 5.2 The gating factor and its dimensional status

Introduce a dimensionless gating factor $g(\text{depth}) \in [0,1]$
multiplying the drag, with depth $= r_\text{atom} - r_\text{droplet}$
(negative inside, positive outside) as in the baseline droplet
potential (`potentials.py:73`):
$$F_\text{drag,gated}(v,\text{depth}) = g(\text{depth})\cdot
F_\text{drag}(v).$$
Because $g$ is dimensionless, the gated force keeps units
$\text{amu}\cdot\text{Å/ps}^2$ for any choice of $g$ — no candidate is
dimensionally excluded; the choice is purely physical.

**The gate applies to drag *and* noise (hard FDT coupling).** The
fluctuation–dissipation relation ties the noise amplitude to the same
friction the gate modulates. The gated noise amplitude is therefore
$$\sqrt{2\,\gamma(v)\,g(\text{depth})\,k_B T_\text{eff}},$$
carrying the *same* $g(\text{depth})$. If the gate applied to drag only,
the ion would receive thermal kicks in vacuum where it feels no
friction — an FDT violation producing spurious heating outside the
droplet. Gating both is mandatory, not optional.

### 5.3 Existing geometry in the baseline

Two notions of "inside" already exist and bound the choices:

- a **sharp boolean** `depth < 0`, used as a hard cutoff by the
  hard-sphere collision sampler (baseline §6.2);
- a **smooth erf transition** over `potential_steepness = 14.2 Å`, used
  by the droplet confining potential
  $V_\text{drop}(r) = \tfrac12(\text{erf}(\text{depth}/\text{steepness})
  +1)\,E_b$.

The gating decision largely reduces to which of these drag follows, or
whether it gets its own profile.

### 5.4 Candidate gating forms

**G1 — sharp boolean** $g = \mathbb{1}[\text{depth} < 0]$. Full inside,
zero outside, step at the surface. Maximal consistency with the
collision model being replaced. But a step discontinuity in force
breaks the smoothness the BAOAB O-step (§4) assumes and injects a
spurious impulse as the ion crosses the surface — the same
integrator-artifact concern that drove the §4 choice. Physically crude:
He density does not vanish at a mathematical surface.

**G2 — erf gate tied to the confining potential**
$g(\text{depth}) = \tfrac12\big(1 - \text{erf}(\text{depth}/
\text{steepness})\big)$, reusing `potential_steepness`. Smooth,
differentiable, 1 deep inside → 0.5 at the nominal surface → 0 outside,
over the *same* ~14 Å shell where confinement turns on. Zero new
parameters. Physically coherent: drag and confinement are both proxies
for local He, so tying them to one profile is the parsimonious claim.

**G3 — erf gate with independent drag steepness** same shape as G2 but
with a separate `drag_gate_steepness`. Admits that the He shell relevant
to *momentum transfer* may have a different effective range than the
one relevant to *binding*. One new parameter — but the extraction window
deliberately excludes the surface-crossing transient ($t^*$ cut, §2.3),
so there is likely *no clean TDDFT data on the gating profile itself*,
making this extra parameter effectively uncalibratable for now.

**G4 — density-proportional** $g(\text{depth}) =
\rho_\text{He}(\text{depth})/\rho_\text{bulk}$. The most physically
honest: drag *is* momentum transfer to He, so it scales with local He
number density. Its decisive advantage is internal consistency — the
*same* $\rho_\text{He}(\text{depth})$ would gate drag, the FDT noise,
*and* the §2.3 mass-accretion rate $\dot M \propto
\rho_\text{He}(\text{depth})$. One density profile would then drive
every helium-coupling channel in the model.

### 5.5 Sub-decision selection

**Primary — G4, density-proportional, collapsing to G2 for now.**
Chosen as the principled framing: drag, noise, and mass-rate all become
functions of one local He density, removing gating as an independent
modelling choice (it falls out of the density profile). **No
$\rho_\text{He}(\text{depth})$ profile is currently defined** beyond the
flat interior `DENSITY_DROPLET = 0.8\,\rho_\text{bulk}` and the erf
surface, so today G4 is *implemented as* the erf complement (G2),
reusing `potential_steepness`. The two are operationally identical until
a measured or computed density profile exists. When one does, G4 becomes
genuinely distinct without an architecture change — a different
$\rho_\text{He}(\text{depth})$ feeds the same gate. **Flagged as a
future improvement up for testing.**

**Secondary — G2 as the explicit erf-tied gate.** The same function,
framed as "tied to confinement geometry" rather than "tied to density."
Identical today; kept as the named fallback if the density framing's
profile dependence proves awkward.

**Secondary — G3, independent drag steepness.** Held in reserve for the
case where the cross-check shows the surface falloff matters *and* a
calibration source for the separate steepness becomes available
(e.g. a re-extraction that retains the surface-crossing region). Not
calibratable from the current extraction.

**Discarded — G1, sharp boolean.** Discontinuous force breaks BAOAB
smoothness and injects a surface-crossing impulse that would contaminate
the final-velocity observable. Recorded discarded with the note that it
is what the old collision model did (`depth < 0`), so the discontinuity
is a known property of the model being replaced, not a regression.

### 5.6 Scope boundary — no bubble / near-field term

The TDDFT calibration was run with the ion *inside* the droplet, so the
helium-bubble physics (the near-field density depletion around the ion)
is already baked into the extracted $\gamma(v)$. Spatial gating
therefore handles **only** the macroscopic droplet-surface cutoff
(am-I-still-in-the-droplet), **not** a near-field bubble density hole.
Recorded explicitly so a future contributor does not add a near-field
density term and double-count the bubble physics already inside the
drag law.

### 5.7 Interchangeability surface

- `SimConfig.drag_spatial_gate ∈ {density_proportional, erf_tied,
  erf_independent, sharp}` — primary `density_proportional` (G4),
  implemented as the erf complement until a density profile exists;
  `erf_tied` (G2) the explicit-erf secondary; `erf_independent` (G3)
  the separate-steepness secondary; `sharp` (G1) discarded, retained in
  the enum only for regression comparison against the old model.
- `SimConfig.drag_gate_steepness` — used only by `erf_independent`;
  defaults to `potential_steepness` (14.2 Å) so it reduces to G2 when
  unset.
- `SimConfig.helium_density_profile` — placeholder for the future
  $\rho_\text{He}(\text{depth})$ that promotes G4 beyond the erf
  collapse; until defined, `density_proportional` uses the erf
  complement internally.
- The gate multiplies both the deterministic drag and the FDT noise
  amplitude (§5.2); the BAOAB O-step (§4) reads $g(\text{depth})$ at
  each step alongside $\gamma(v)$.

---

## 6. Validation surface and tolerances

### 6.1 Physical question

Every "primary" chosen in §1–§5 was deferred to empirical cross-check
against the TDDFT and experimental references. §6 specifies the
instrument that resolves those deferrals. It is not a postscript: if the
validation surface is weak, every upstream primary stays unfalsified and
the interchangeability apparatus produces nothing. The section has three
jobs — extend the baseline validation surface to observables it does not
currently exercise, order the comparisons so the large undetermined
parameter space stays separable, and enforce the consistency constraint
between the mass scenario (§2) and the extracted drag coefficients (§3).

### 6.2 What the baseline validation surface cannot do

The baseline (PHYSICS_BASELINE §13) ships trajectory comparison
(`compare_distance`, `compare_velocity_magnitude`), final-velocity
histograms against the VMI references, and energy-balance plots. These
were built for a *deterministic* model matched to *single* reference
trajectories. The drag port breaks three of their assumptions:

- **Noise lives in second moments (§1.6).** `compare_distance` matches
  one trajectory to one reference and is blind to ensemble spread; it
  cannot discriminate any noise variant.
- **Mass scenarios need a new observable (§2.10).** The discriminating
  measurement between Scenarios A/B/biphasic is the terminal I⁺(He)ₙ
  size distribution, which the code does not currently exercise.
- **The form cross-check needs the comparison to be clean (§3, §4).**
  The whole reason BAOAB was chosen (§4.6) is so a trajectory difference
  is attributable to the drag form rather than to an integration
  artifact. That payoff is only realised if the validation can isolate
  the form from the other unknowns.

### 6.3 The attribution problem

The undetermined parameters interact: drag form (4), noise
form/calibration/geometry (§1), mass scenario (4+) and its velocity
scaling (3), spatial gate (4). The observables are coupled — final
velocity depends on form *and* mass *and* gate, while its *spread*
depends on noise *and* ensemble size. A flat "run everything, compare to
VMI" cannot separate these. The resolution is a **validation hierarchy
ordered by separability**: each tier isolates as few unknowns as
possible and fixes its winner before the next tier's unknown is
introduced. Validation is therefore **sequential, not simultaneous** —
itself a recorded design decision.

### 6.4 The sequential validation hierarchy

**Tier 0 — drag form, deterministic, fixed mass.** Run
`mass_scenario = fixed`, noise amplitude zero, against the TDDFT
`9A`/`18A` distance and velocity traces *inside the extraction window
only*. Isolates the **drag form** (§3) with no mass evolution, no noise,
no transient. The cleanest possible test; must come first because every
later tier is contaminated until the form is pinned.

*Cleanliness condition (see §6.6).* Tier 0 is cleanest when the
simulation mass treatment matches the one the reference trajectory was
generated under. The TDDFT shell mass is *not* constant (it declines
~21→14 He across the window, §6.6), so a fixed-$m_\text{eff} \approx
203$ amu run matches the reference closely mid-window (where the shell
is ~19 He) but drifts at the ends. For a *strictly* clean form
isolation, impose the reference $m(t)$ (prescribed, not dynamically
evolved) with coefficients extracted under that same $m(t)$; *then* a
mismatch is purely the form. For a first pass the fixed-$m_\text{eff}$
run is adequate, since the residual mass mismatch is small over most of
the window.

**Tier 1 — mass scenario, deterministic.** With the Tier-0 form fixed
and noise still off, turn on each mass scenario (A/B/biphasic) and
compare the full post-transient trajectory plus the time-resolved shell
trajectory (~21→~19→~14 He, §2.1). Isolates the **mass scenario** (§2). The transient-tolerance
relaxation (§6.7) lives here. Each non-`fixed` scenario requires its own
coefficient extraction (§6.5) before its trajectory comparison carries
meaning.

**Tier 2 — terminal mass distribution.** Compare the simulated I⁺(He)ₙ
size distribution against the experimental detector data. This is the
*only* observable that sharply separates A/B/biphasic and is a
distribution-vs-distribution comparison the code cannot currently do.
The reference is a histogram of n-counts per integer snowball mass
(I⁺, I⁺He₁, I⁺He₂, …), i.e. a *discrete* distribution over integer n.

**Tier 3 — ensemble second moments.** With noise on, compare the
final-velocity histogram *width* and (if T3 geometry, §1.4) angular
spread against the VMI references, plus cross-trajectory variance on the
8000-atom `single_pulse_droplet_distribution`. Isolates the **noise
model** (§1) given form + mass + gate from Tiers 0–2.

### 6.5 The mass-scenario / coefficient consistency constraint

The force balance underlying extraction,
$$F_\text{drag}(t) = m(t)\,a(t) - F_C(R(t)),$$
constrains only the *combination* $m(t)\,a(t)$ against the measured
$a(t), R(t)$. What is attributed to drag depends entirely on the assumed
$m(t)$: two mass assumptions applied to the *same* TDDFT trajectory yield
two *different* drag laws, each self-consistent only when re-simulated
under its own mass assumption. Consequences:

- A drag law extracted at constant $m_\text{eff}$ is self-consistent
  **only** when re-applied at constant $m_\text{eff}$. Running it under
  an evolving-mass scenario applies the law to a mass trajectory
  different from its extraction — the trajectory will miss the reference
  for a reason that is not physics.
- `mass_scenario = fixed` at $m_\text{eff}$ is therefore the **only**
  scenario for which the in-hand constant-mass coefficients are
  self-consistent — a stronger statement than §2.6's "no-extrapolation
  baseline": it is the only scenario matching the extraction's own mass
  assumption.
- `mass_scenario ∈ {A, B, biphasic}` requires coefficients
  **re-extracted under that scenario's time-varying $m(t)$**. Using
  constant-mass coefficients here is the §2.2 extrapolation, now
  diagnosed precisely: the law was solved under a *different* $m(t)$ and
  is no longer the unique drag consistent with the data.

`mass_scenario` and `drag_coefficients` are thus a **coupled pair**, not
freely-mixable enums. **Enforced as a config-load consistency guard:**

- Each coefficient bundle carries
  `drag_coefficients.extraction_mass_model ∈ {constant, time_resolved}`
  plus the constant value (or the $m(t)$ reference) it was extracted
  under. This requires the extraction pipeline to **stamp** each bundle
  with its mass metadata — an extraction-side action item (§6.8), since
  the current `fit_result` records only $\{a,b\}$ / $\{\gamma,n\}$.
- `fixed` at value $M$ requires `extraction_mass_model = constant` at
  the same $M$ within a tolerance of **~1–2 He (~4–8 amu)** — the regime
  where the curve is genuinely insensitive to mass. This is where the
  "±1–2 He is negligible" argument lives legitimately: as the *width* of
  the consistency check, not as licence to ignore a larger gap.
- non-`fixed` scenarios require `extraction_mass_model = time_resolved`
  under a matching $m(t)$.
- Pairing a non-`fixed` scenario with constant-mass coefficients is the
  inconsistent case → **refuse to run** (hard error). One escape hatch:
  `SimConfig.allow_inconsistent_mass_pairing = True` (default `False`)
  downgrades the refusal to a loud warning, for deliberate exploratory
  runs only.

### 6.6 Effective mass is a varying quantity collapsed to a constant

TDDFT shell counts (per iodine, 9 Å case): ~21 He pre-explosion, ~19 at
10 ps, ~14 at 14 ps — a monotone *decline* of ~⅓ across the trajectory.
The drag extraction collapses this to a constant $m_\text{eff} \approx
203$ amu (~19 He), the window-representative value. Two facts follow:

- **The constant is well-chosen but still a constant.** 19 He sits near
  the mean of the declining shell, so the drag law is calibrated at a
  mass the ion genuinely carries through the middle of the fit window —
  there is no large level error to correct. The only residual unrigour
  is the collapse of a ~⅓-varying $m(t)$ to a single number.
- **The residual error is mild and bounded.** Because $F_\text{drag} =
  m\,a - F_C$ is a *difference* of comparable terms, the fractional
  error on $F_\text{drag}$ from the constant-mass approximation tracks
  $|m(t) - m_\text{eff}|/m_\text{eff}$, which is near-zero mid-window
  and reaches ~⅓ only at the trajectory ends (~14 He late, ~21 He
  early). This is within the regime where the curve is relatively
  insensitive to mass, so the in-hand coefficients are usable as-is for
  the constant-mass scenarios; time-resolved re-extraction is an
  *optional refinement*, not a correction of a flaw.

**Constant mass is not the fully conservative choice it appears to be.**
The original rationale for a constant — distrust that the gentle
shell evolution (9 Å bubble) transfers to the violent
equilibrium-distance onset ($R_e \approx 2.6$ Å, Coulomb ~12× stronger)
— is legitimate, but a constant is *also* a mass law ($\dot M = 0$),
measured nowhere near the violent onset and not matching the declining
trend even in the gentle case. The genuinely conservative statement is
that the **first ~0.5 ps violent transient is uncalibrated for both drag
and mass under every scenario** — generalising §2.3's Scenario-A-only
concession.

**Optional time-resolved re-extraction (extraction-side, §6.10).** The
one refinement that removes even the residual constant-mass unrigour:
interpolate $m(t)$ from the shell counts and use it directly in the
force balance $F_\text{drag}(t) = m(t)\,a(t) - F_C(R(t))$, rather than
the constant $m_\text{eff} = 203$ amu. This is *not required* for the
constant-mass scenarios (`fixed`, and A/B/biphasic validated against the
constant-mass coefficients within the §6.5 tolerance), but it is the
most defensible extraction and is the natural pairing for an
evolving-mass scenario whose $m(t)$ matches the reference trajectory
(notably Scenario B, whose measured loss *is* the $m(t)$ to use). The
architecture is unaffected — the sim still applies $-\gamma(v)v$ at face
value per §2.2; only the coefficients change.

### 6.7 Transient tolerance — generalised free-extrapolation zone

No TDDFT data exists at the violent bubble-exit onset (the $t^*$ cut,
§2.3, excludes it), so neither the drag law nor the mass law is anchored
in the first ~0.5 ps — under *any* scenario, not only Scenario A.
Validation tolerances on $v(t)$ and $R(t)$ during this window must be
**loose for all scenarios**; the mid- and late-trajectory observables
remain the anchor. This generalises the §2.3 Scenario-A-specific note.

### 6.8 Empirical lean toward Scenario B (recorded, not auto-promoting)

The TDDFT shell data (21→19→14, monotone loss) empirically points to
Scenario B (stripping), the §2 *secondary*, over Scenario A (accretion),
the §2 *primary* — a pre-simulation evidence point. B is additionally
privileged: it is the only scenario whose time-resolved $m(t)$ can be
read **directly off the reference** (the measured loss *is* the $m(t)$
to use in option-3 extraction), whereas A requires *modelling* an
accretion $m(t)$ the gentle TDDFT does not provide. This strengthens the
§2.6 hint (the 0.09→0.005 attach-rate drop) with a second, independent
argument. **Not auto-flipped:** the hierarchy (§6.4) exists precisely to
test A against B; A remains the simplest hypothesis and the recorded
primary until Tier 1/Tier 2 adjudicate. Recorded as a note against §2.6.

### 6.9 Histogram comparison metric

Tiers 2 and 3 compare distributions, requiring a divergence metric. The
choice has consequences and is tiered like the rest.

**Primary — Wasserstein (earth-mover) distance.** Reports the cost of
transporting one distribution onto the other in the *physical units of
the observable* — Å/ps for final velocity, integer He-count for the size
distribution. Binning-free, sensitive to overall shape and location, and
directly interpretable ("the distributions differ by ⟨Δ⟩ He on
average"). For the discrete integer-n size distribution (§6.4) the
1-D Wasserstein on integer support is exact and natural. Best default
for both distribution tiers.

**Secondary — binned $\chi^2$.** Respects the existing VMI bin
convention directly (0.04 Å/ps internal bins, 15-bin moving mean,
display to 2800 m/s; PHYSICS_BASELINE §13), so it slots into the current
histogram tooling with least change. Sensitive to per-bin disagreement
and gives a familiar goodness-of-fit number, but depends on binning and
needs care with low-count bins in the size-distribution tail (most ions
bare, falling off monotonically). Kept for continuity with the baseline
VMI comparison.

**Secondary — Kolmogorov–Smirnov (max CDF gap).** Distribution-free,
parameter-free, good as a quick scalar screen. But it is most sensitive
near the distribution median and relatively insensitive in the tails —
a poor fit for the size distribution, whose discriminating content is in
the *tail* (how many ions retain large shells). Kept as a cheap
first-pass screen, not the adjudicating metric.

### 6.10 Open items and interchangeability surface

- **Numeric acceptance thresholds are deferred.** §6 fixes the metrics
  and the hierarchy; the pass/fail numbers are set once the first runs
  are seen, since everything upstream is being cross-checked empirically.
- **New validation hooks required** (the baseline does not provide
  them): an ensemble-variance comparison (Tier 3) the single-trajectory
  `compare_distance` / `compare_velocity_magnitude` cannot give (§1.6);
  a discrete-n size-distribution comparison routine plumbing the
  experimental I⁺(He)ₙ history into the surface (§2.10, Tier 2).
- **Maximum trajectory speed** for the `linear_cubic` turnover guard
  (§3.8) is sourced here — from the TDDFT trajectories or a generous
  ceiling — as it has no home in the baseline config.
- **Extraction-side action items** consolidated: stamp coefficient
  bundles with `extraction_mass_model` metadata (§6.5); the *optional*
  time-resolved-$m(t)$ re-extraction (§6.6) — a refinement, not a
  correction, since the constant $m_\text{eff} \approx 203$ amu is
  already window-representative; the outstanding `linear_quadratic` /
  `threshold` fit passes (§3.7).
- `SimConfig.allow_inconsistent_mass_pairing` (default `False`) — the
  only new field §6 introduces; the consistency guard (§6.5) is
  otherwise a config-load check over existing fields plus the coefficient
  metadata.
- `SimConfig.validation_histogram_metric ∈ {wasserstein, chi2, ks}` —
  primary `wasserstein`; selects the Tier 2/3 divergence metric.
