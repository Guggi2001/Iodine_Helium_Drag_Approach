# Drag-Model Port — Design Decisions

**Status:** Open planning document. Each section captures one architectural
choice in the migration from the hard-sphere collision model
(`i2\_helium\_md/physics/collisions.py`) to a TDDFT-calibrated drag-force
model for I⁺ in a helium bubble.

**Format.** Each section has the same structure:

* *Physical question* — what is being decided and why it matters.
* *Primary* — the approach chosen for the first implementation.
* *Secondary* — viable fallbacks to keep pluggable behind the same
interchangeability surface.
* *Discarded* — options ruled out, with the reason recorded so they do
not get re-litigated.
* *Interchangeability surface* — the module boundary or `SimConfig`
field behind which the choice lives, so swap-in design stays visible.

The document is intentionally implementation-agnostic. No code, no
pseudo-code. Mathematical formulation in LaTeX where needed.

\---

## 0\. Document map (entry point)

**This document is the single entry point** for the I⁺/He drag-port +
mass-dynamics planning. Three artifacts, with a strict ownership rule:

* **`DRAG\_PORT\_DESIGN\_DECISIONS.md` (this doc)** — entry point and the **drag
spine** (noise, drag form, integrator, spatial gate, validation tiers). Stable:
Tier 0 is complete and the drag law is locked.
* **`MASS\_DYNAMICS\_LOCKED\_energy\_gated\_evaporation.md`** — the **live mass-model
detail**: two-channel mechanism, energy budget, ladder, early window, all
risks (R1–R12), assumptions (A1–A11), and open questions (OQ1–OQ6). Mechanism
locked, calibration open.
* **`CALIBRATION\_MAP.md`** — the **cross-doc parameter index**: every parameter →
class (locked/sourced/derived/bounded/free) → Tier anchor → cross-check.

**Navigation rule.** Anything carrying a *revision date*, an *OQ flag*, or a
*calibration value that can still move* lives in **exactly one place — the MASS
doc** — and is reached from here by reference. Sections in this doc describe
*structure* (mechanism shape, which doc owns what), not drifting values. §2
below is a structural summary of the mass model; its detail, risks, and open
calibration are owned by the MASS doc.

**Consolidation status (2026-06-15).** This is the *light* entry-point
consolidation: the superseded three-scenario framing (former §2) is replaced by
the summary below; the MASS doc remains the live detail doc. A **full merge**
(absorbing MASS wholesale, retiring it) is deliberately **deferred to the
OQ-freeze**, which triggers when:

1. **OQ1** resolves (drag electronic-state provenance — can still flip A10);
2. the **ladder shape** and **electronic picture** (the only two genuinely-free
knobs) are Tier-2-arbitrated; and
3. **OQ6** is checked ($E\_\\infty$ stripping reachability).
Until then MASS is still moving and is kept separate by design.

\---

## 1\. Noise model

### 1.1 Physical question

The hard-sphere model is not "deterministic dissipation plus noise" —
it *is* a stochastic process whose mean happens to look like drag.
Replacing it with a deterministic $F\_\\text{drag}(v)$ alone discards the
entire fluctuation channel: trajectory ensemble spread, transverse
randomisation, rare-event tails, and the thermal floor at the He bath
temperature $T\_\\text{particles\_K} = 0.4,\\text{K}$.

Adding a Langevin noise term restores stochasticity, but three sub-
decisions follow immediately and they couple:

* **Functional form of the noise term** (additive vs. multiplicative
vs. empirically-anchored).
* **Calibration purpose** (what physical thing the amplitude is meant
to reproduce).
* **Geometric structure** (longitudinal-only, isotropic 3D, or
anisotropic with separate transverse component).

These three are decided independently below.

\---

### 1.2 Functional form

The Langevin equation under consideration:

$$m,\\dot v = F\_\\text{ext}(t) - \\gamma(v),v + \\mathcal{N}(v, t)$$

**Friction-coefficient convention (used document-wide).** $\\gamma(v)$ is
a *force coefficient* with units $\[\\text{amu/ps}]$, defined directly
from the extracted drag law as $\\gamma(v) = |F\_\\text{drag}(v)|/v$ — so
the friction force is $\\gamma(v),v$ (units $\\text{amu·Å/ps}^2$, a
force, with **no** leading $m$). The corresponding friction *rate* is
$\\gamma(v)/m$ $\[\\text{1/ps}]$; it appears only inside the BAOAB damping
exponent $e^{-\\gamma,dt/m}$ (§4.3), never as a multiplier on the force.
This is the single source of the factor-of-$m$ that distinguishes the
two appearances.

$\\gamma(v)$ is $v$-dependent under either drag form (§3): for the
primary linear+cubic drag, $\\gamma(v) = a + b,v^2$ (finite at $v=0$,
units amu/ps since $\[a]=\\text{amu/ps}$, $\[b,v^2]=\\text{amu·ps/Å}^2
\\cdot \\text{Å}^2/\\text{ps}^2 = \\text{amu/ps}$); for the secondary
power-law drag with $n<0$, $\\gamma(v) = \\gamma\_\\text{PL},v^{,n-1}$
(where $\\gamma\_\\text{PL}$ is the power-law coefficient from §3.8's
${\\gamma,n}$ bundle; divergent at $v=0$). This $v$-dependence is what
makes the choice non-trivial.

#### Primary — N2: multiplicative noise with local FDT

$$\\mathcal{N}(v,t) = \\sqrt{2,\\gamma(v),k\_B T\_\\text{eff}};\\xi(t)$$

with $\\langle\\xi(t)\\rangle = 0$, $\\langle\\xi(t)\\xi(t')\\rangle =
\\delta(t-t')$ (so $\\xi$ is white noise, $\[\\xi(t)]=\\text{ps}^{-1/2}$).
There is **no** leading $m$ inside the root: the textbook form
$\\sqrt{2,m,\\gamma\_\\text{rate},k\_BT}$ uses the friction *rate*
$\\gamma\_\\text{rate}=\\gamma/m$, so $m,\\gamma\_\\text{rate}=\\gamma$ and the
$m$ cancels in this convention. Dimensional check:
$\[,2\\gamma k\_BT,] = (\\text{amu/ps})(\\text{amu·Å}^2/\\text{ps}^2) =
\\text{amu}^2\\text{Å}^2/\\text{ps}^3$; its root times
$\[\\xi]=\\text{ps}^{-1/2}$ gives $\\text{amu·Å/ps}^2$, a force. Balances.

The noise amplitude tracks the *local* friction coefficient, so the
fluctuation–dissipation relation is honoured instantaneously rather
than in some averaged sense. The price is that the SDE is now
non-trivial:

* **Itô/Stratonovich ambiguity.** Multiplicative noise SDEs give
different equilibrium distributions under different stochastic
calculi. Stratonovich is the standard physics choice (limit of
smooth-noise OU processes) and is what the BAOAB-family integrators
assume.
* **Low-$v$ regularisation — only under the secondary power-law drag.**
Under the primary linear+cubic drag, $\\gamma(v) = a + b,v^2$ is
finite at $v=0$, so the noise amplitude is well-defined down to rest
and no regulariser is needed. Under the secondary power-law drag with
$n<0$, $\\gamma(v)\\to\\infty$ as $v\\to 0$ and the noise amplitude
diverges along with the drag; the same regulariser used for the
deterministic drag (see the §3.8 regularisation note) must then apply to the noise.
* **$T\_\\text{eff}$ is not necessarily $T\_\\text{particles\_K}$.** See
§1.3 for the calibration discussion — the effective temperature in
the FDT relation is a fit parameter under noise purpose (b).

*Why primary:* the extracted drag law is nonlinear, so the friction is
$v$-dependent in any case. Anchoring the noise to that same $v$-
dependence is the most internally consistent choice. The Itô/
Stratonovich choice is settle-once and the integrator handles it.

**Correction under the locked pure-cubic form (2026-06-15).** Tier-0
locked the $a=0$ pure-cubic drag (§3.4), so
$\\gamma\_0 \\equiv \\lim\_{v\\to0}\\gamma(v) = 0$ — **not** the finite $a$ this
N2 argument assumed above. The FDT noise amplitude
$\\sqrt{2,\\gamma(v),g,k\_B T\_\\text{eff}}$ therefore $\\to 0$ as $v\\to0$:
the ion at rest feels neither friction nor thermal kick. This is the
*opposite* failure mode from the discarded N1 / $n<0$ case (divergence
at rest); here the noise vanishes. **Accepted**, because (i) the
strict-FDT bath kick was already shown dynamically null
($\\sim 7\\times10^{-4}\\ \\text{Å/ps}$ vs $\\sim10\\ \\text{Å/ps}$ working
speed, §1.3(a)) — a floor that vanishes at rest removes nothing
observable; (ii) no finite rest-friction was ever data-anchored;
(iii) $v=0$ is a regular point ($F\\to0$ smoothly), so no regulariser is
needed and `drag\_low\_v\_floor` stays inert (already inert for the
$n\\approx+2$ export; now also for pure cubic). The N2 *expression* is
unchanged; only "$\\gamma\_0=a$ finite" is corrected to "$\\gamma\_0=0$,"
and N2 still reduces correctly (trivially) at low $v$. *Blanket:* read
every remaining §1–§3 reference to "$\\gamma\_0=a$" or "finite $\\gamma$
at rest" for `linear\_cubic` with $a=0$ substituted.

#### Secondary — N3: empirical noise from TDDFT residual variance

$$\\sigma\_\\xi^2 = \\text{Var}\\big(v(t) - v\_\\text{smoothed}(t)\\big)\\Big|\_{\\text{post bubble-mode removal}}$$

The amplitude is read directly off the residual of the CEEMDAN+SG
smoothing pipeline already in `Drag\_extraction\_code.md`. No FDT, no
temperature, no model assumption beyond "what the data wasn't
smoothed away into the drag fit is, by definition, noise."

*Why secondary:* model-free and anchored to the same dataset as the
drag itself — same provenance, same systematic biases. But:

* Only one TDDFT trajectory per case (9 Å, 18 Å). Cannot separate
trajectory-to-trajectory variance from within-trajectory variance.
* Conflates numerical noise, mean-field error, and any residual
bubble dynamics that escaped the CEEMDAN IMF drop. Not the same
thing as a He bath kicking the ion.
* Gives one number per case, not a function of $v$. So it would be
applied as an additive (constant amplitude) noise — losing the
$v$-dependence that motivated N2 in the first place.

Keep as a trial alternative when comparing against the TDDFT curves —
it's the natural "atheoretical" baseline.

#### Discarded — N1: additive noise FDT-anchored at $\\gamma\_0 = \\gamma(v\\to 0)$

$$\\mathcal{N}(v,t) = \\sqrt{2,\\gamma\_0,k\_B T};\\xi(t),\\qquad \\gamma\_0 = \\lim\_{v\\to 0}\\gamma(v)$$

Rejected on the dynamical-scale ground, with a note on form-dependence:

* Under the primary linear+cubic drag, $\\gamma\_0 = a$ is finite, so
N1 is mathematically well-defined. But anchoring the noise scale at
the low-velocity limit while the ion spends most of its time at
eV-scale KE sets the noise by the wrong dynamical scale — the
fluctuations the ion actually experiences are governed by $\\gamma$ at
its working speed, not at rest.
* Under the secondary power-law drag with $n = -2$, $\\gamma\_0\\to\\infty$
and the amplitude is undefined outright. N1 is then not merely
ill-scaled but mathematically incompatible with the drag form.

In both cases N2 (local FDT) dominates: it reduces to the *correct*
limit of N1 at low $v$ for the linear+cubic form while tracking the
working-speed friction everywhere else.

#### Interchangeability surface

`SimConfig.noise\_form ∈ {multiplicative\_local\_fdt, empirical\_residual, none}`.
The integrator dispatches on this. Adding new variants later requires
adding a branch but not changing the field signature.

\---

### 1.3 Calibration purpose

Even with the functional form fixed (N2), the noise amplitude depends
on a temperature-like parameter $T\_\\text{eff}$ (or, for N3, on the
residual variance). What that parameter is *calibrated to* is a
separate choice.

#### Primary — (b): calibrate to hard-sphere trajectory variance

Run the existing hard-sphere model on the reference HeDFT cases (9 Å,
18 Å) and on the `single\_pulse\_droplet\_distribution` preset. Extract
the empirical variance of:

* final ion speed across the ensemble,
* longitudinal velocity at intermediate checkpoints,
* (for transverse noise, if T3 active — see §1.4) angular spread of
final velocity vectors.

Fit $T\_\\text{eff}$ such that the new Langevin model reproduces those
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
form (N3) to set $T\_\\text{eff}$ even when the *form* is N2. I.e., back
out $T\_\\text{eff}$ such that

$$\\sqrt{2,\\gamma(v\_\\text{typ}),k\_B T\_\\text{eff}} = \\sigma\_\\text{residual}$$

at some representative velocity $v\_\\text{typ}$ (median or RMS over the
extraction window).

*Why secondary:* lets purpose-(c) and form-N3 share a calibration
pathway. If primary calibration (b) gives suspicious results, this is
the natural cross-check using independent data.

#### Secondary — (a): strict FDT at $T\_\\text{particles\_K} = 0.4,\\text{K}$

Take $T\_\\text{eff} = 0.4,\\text{K}$ as a hard physical constraint and
let the noise amplitude be whatever it is.

Order-of-magnitude estimate: with $m \\approx 203,\\text{amu}$,
$k\_B T \\approx 3.4\\times 10^{-5},\\text{eV} \\approx 0.33,
\\text{amu·Å}^2/\\text{ps}^2$, $\\gamma\_0 \\sim \\text{few amu/ps}$,
$dt = 0.01,\\text{ps}$, the per-step OU kick (coefficient convention,
$\\Delta v \\sim \\sqrt{2\\gamma\_0 k\_B T,dt},/,m$) is

$$\\Delta v\_\\text{noise} \\sim \\frac{\\sqrt{2\\gamma\_0 k\_B T,dt}}{m}
\\sim 7\\times10^{-4},\\text{\\AA/ps}$$

versus typical ion speeds $\\sim 10,\\text{\\AA/ps}$ — four to five orders
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

`SimConfig.noise\_calibration ∈ {hard\_sphere\_variance, tddft\_residual, strict\_fdt\_bath}` with the resulting $T\_\\text{eff}$ (or equivalent
amplitude) stored alongside. The calibration itself is a
pre-processing step that produces a number; the simulation does not
need to know which procedure produced it.

\---

### 1.4 Geometric structure (transverse noise)

The drag force acts along $-\\hat v$. The noise can act along $\\hat v$
only, isotropically in 3D, or anisotropically with separate
amplitudes parallel and perpendicular to $\\hat v$.

The TDDFT extraction provides longitudinal information only —
trajectories are effectively 1D along the Coulomb dissociation axis.
Transverse noise is *unanchored* by the extraction pipeline; whatever
is chosen needs an external calibration source.

Strategy: **start minimal, escalate if needed.**

#### Primary — T1: longitudinal noise only

$$\\vec{\\mathcal{N}}(v,t) = \\mathcal{N}\_\\parallel(v,t),\\hat v$$

Only the speed fluctuates; the direction is determined entirely by
deterministic forces (drag + Coulomb + droplet).

*Why primary:*

* Matches the dimensionality of the extracted drag law (1D in →
1D out).
* Zero free parameters beyond the longitudinal amplitude already
decided in §1.2-§1.3.
* Minimum viable noise model. Easiest to validate.
* The hard-sphere ion–He mass ratio is $\\rho = m\_I/m\_\\text{He}
\\approx 32$, so per-collision lateral deflection is small. Whether
*cumulative* transverse spread matters is an empirical question
that T1 will answer by failing or not failing the VMI angular
distribution validation.

*Validation criterion for escalation:* compare the final-velocity
angular distribution against the experimental VMI references
(`vmi\_iplus\_he.csv`, `vmi\_iplus\_gas.csv`). If T1 fits, T2/T3 are
unnecessary. If T1 underestimates angular spread, escalate to T3.

#### Secondary — T3: anisotropic with separate transverse amplitude

$$\\vec{\\mathcal{N}}(v,t) = \\mathcal{N}*\\parallel,\\hat v + \\mathcal{N}*{\\perp,1},\\hat e\_1 + \\mathcal{N}\_{\\perp,2},\\hat e\_2$$

with $\\hat e\_1, \\hat e\_2$ spanning the plane perpendicular to $\\hat v$,
and $\\sigma\_\\perp \\neq \\sigma\_\\parallel$ in general.

*Calibration target for $\\sigma\_\\perp$:* **VMI experimental angular
distribution**, not hard-sphere runs. Reasoning: $\\sigma\_\\parallel$ is
already FDT-tied to the drag (§1.2), and the drag is anchored to
TDDFT. So $\\sigma\_\\parallel$ has TDDFT provenance. $\\sigma\_\\perp$ has
no internal anchor and should be tied to the most authoritative
external dataset available, which is the experimental VMI. Calibrating
$\\sigma\_\\perp$ against hard-sphere runs would compound model error.

*Open issue — basis at $v\\to 0$:* the basis $(\\hat v, \\hat e\_1,
\\hat e\_2)$ is ill-defined at $v = 0$.

* **Primary low-$v$ behaviour:** $\\sigma\_\\perp \\to 0$ smoothly as
$v \\to 0$. Concretely, gate $\\sigma\_\\perp$ by some smooth function
of $v$ that vanishes at the origin. Loses the transverse channel at
rest but is unambiguous and continuous.
* **Alternative to test:** smoothly blend to isotropic
($\\sigma\_\\perp \\to \\sigma\_\\parallel$ as $v \\to 0$, with the basis
becoming irrelevant in the isotropic limit). Avoids losing
stochasticity at rest but introduces a velocity-dependent
isotropisation knob.

Both should be testable behind the same flag; pick "vanish" as the
first try.

*Why secondary not primary:* introduces a second amplitude parameter
that requires its own calibration source. Worth the cost only if T1
demonstrably fails the VMI angular validation.

#### Secondary — T2: isotropic 3D noise

$$\\vec{\\mathcal{N}}(v,t) = (\\mathcal{N}\_x, \\mathcal{N}\_y, \\mathcal{N}\_z), \\quad \\langle\\mathcal{N}\_i\\mathcal{N}*j\\rangle \\propto \\delta*{ij}$$

Same amplitude in all three Cartesian directions, independent
components. The drag still acts along $\\hat v$, but the noise has no
privileged direction.

*Why secondary not primary:* doesn't match the hard-sphere structure
(which is privileged along $\\hat v\_\\text{in}$ per collision). At
$v = 0$ it gives transverse kicks where hard-sphere gives nothing
defined. It also breaks the natural FDT pairing between drag (along
$\\hat v$) and noise (also along $\\hat v$).

*Why kept as backup:* trivially compatible with BAOAB integrators
(additive isotropic noise is the textbook case). If T3's basis
regularisation at $v\\to 0$ turns out to be a numerical headache, T2 is
the easy escape valve. Conceptually clean; one parameter.

#### Discarded — none.

T1, T2, T3 are all physically defensible. Keep all three behind a
single enum.

#### Interchangeability surface

`SimConfig.noise\_geometry ∈ {longitudinal, isotropic, anisotropic}`
with amplitude-related fields downstream of that choice. For
`anisotropic`, an additional field
`SimConfig.noise\_low\_v\_behavior ∈ {vanish, blend\_to\_isotropic}`
controls the $v\\to 0$ regulariser.

\---

### 1.5 Summary of Section 1 choices

|Sub-decision|Primary|Secondary|Discarded|
|-|-|-|-|
|Functional form|N2: multiplicative, local FDT|N3: empirical residual|N1: additive at $\\gamma\_0$|
|Calibration purpose|(b): hard-sphere variance|(c): TDDFT residual; (a): strict FDT at 0.4 K|—|
|Geometric structure|T1: longitudinal only|T3: anisotropic, then T2: isotropic|—|
|Low-$v$ behaviour (T3 only)|vanish|blend to isotropic|—|

### 1.6 Validation criterion specific to noise

The noise model is the only piece of the port whose effect lives
entirely in **second moments** of the output. The mean trajectory is
set by the drag; the noise only changes the spread around it.
Validation must therefore exercise *ensemble* statistics:

* Final-velocity histogram width against `vmi\_iplus\_he.csv`.
* Angular distribution of final velocities (if T3) against the same.
* Cross-trajectory variance at fixed time within an ensemble run on
`single\_pulse\_droplet\_distribution` (8000 atoms, variable droplet).

The existing `compare\_distance` / `compare\_velocity\_magnitude` checks
look at single-trajectory matching and will **not** discriminate
noise-model variants. A new validation hook will be required.

\---

## 2\. Mass model (summary — full detail in MASS doc)

> \*\*Supersedes the former §2 three-scenario framing (2026-06-15).\*\* The original
> Scenario A / B / biphasic trade-off analysis is superseded by the locked mass
> model in `MASS\_DYNAMICS\_LOCKED\_energy\_gated\_evaporation.md` (retained in version
> history). This section keeps only the \*\*stable spine\*\* and the \*\*shared
> infrastructure the drag spine references\*\* (the $m\_\\text{eff}$ framing §2.2 and
> the energy invariant/schema §2.9); the superseded scenario subsections (§2.3–§2.7)
> are stubbed and point to the MASS doc. Mechanism detail, risks (R1–R12),
> assumptions (A1–A11), and open calibration (OQ1–OQ6) are owned by the MASS doc;
> every parameter's provenance by `CALIBRATION\_MAP.md`. \*\*Subsection numbers are
> preserved as stable cross-reference anchors\*\* even where content is condensed.

**Locked mechanism — `biphasic\_energy\_gated` (the stable spine).** Two
discrete-stochastic channels on an integer, non-monotone shell count $n(t)$:
**pickup** (Poisson, rate $\\lambda\_\\text{attach}(\\rho\_\\text{He}(\\text{depth}))$,
$n\\to n+1$) and **evaporation** (energy-gated + RRK rate-limited: shedding
suppressed while net self-unbound $E\_\\text{int}>\\sum\_i D\_0^{,\\mathrm{I^+}}(i)$,
then saturating rate $k=\\nu(1-D\_0(n)/E\_\\text{int})^{,s-1}$, cold shed); plus an
**internal-energy reservoir** $E\_\\text{int}$ cooled (Newton's law) via the
GAH25-matched variable $E\_\\text{solv.struct}=E\_\\text{bind}+E\_\\text{int}$. Terminal
$n$ at 20 ps spans a **regime axis** (${f\_\\text{int},\\tau\_\\text{dissip},$ ladder
depth$}$): shell-retaining (default) ↔ total stripping (Calvo24 limit —
reachable, evaluated in secondary runs, not excluded; OQ6). **Live
cross-dependency:** the drag effective binding (§6.5.1 / Tier 0) and the mass-side
electronic picture / ladder depth (MASS A10) are coupled through **OQ1**, owned by
the MASS doc; this doc's Tier-0 provenance flag (§6.4) points there.

### 2.1 Physical question

What happens to the He solvation shell during/after the Coulomb explosion of I₂⁺
inside the droplet — retain, strip, or equilibrate? Genuinely open physics. Known:
pre-ionisation \~21 He/iodine; mid-flight 9 Å (TDDFT) shell **declines \~21→\~19 (10
ps)→\~14 (14 ps)** with drag-extraction reference mass \~19 He ≈ 203 amu (§2.2);
detector size distribution falls off monotonically. The original A/B/biphasic
bracket is superseded by the locked biphasic mechanism (above; MASS doc). This is
the **generative process for the experimental observable** — the per-fragment
I⁺Heₙ velocity histograms — so mass evolution is load-bearing: a fixed-mass run
yields one species and cannot generate the fragment channels.

### 2.2 The $m\_\\text{eff}$ framing  *(retained — shared infrastructure)*

**$m\_\\text{eff}\\approx203$ amu (\~19 He) is not the ion's true mass; it is the mass
*assumed during drag extraction*** — a parameter of the *drag law*, not the ion.
Consequences: the simulation's instantaneous $m(t)$ is the physical mass and is
used for integrator inertia ($F=m(t)a$) and KE bookkeeping; the drag force is
applied as $F\_\\text{drag}=-\\gamma(v)v$ at face value with no invented
$M$-dependence; when $m(t)\\neq m\_\\text{eff}$ the drag law is **extrapolated**
outside its calibration domain (unrigorous but unavoidable — the extraction
provides no $M$-scaling for $\\gamma$), mild near mid-window and growing toward the
ends. `SimConfig.m\_eff\_amu` exists as the named drag-law reference mass, separate
from the physical $m(t)$. §6.5 makes the mass-scenario↔coefficient pairing a
config-load guard (`extraction\_mass\_model` metadata); §6.6 records the optional
time-resolved-$m(t)$ re-extraction. *(This is the R6 consistency anchor — MASS R6.)*

### 2.3 Scenario A (full strip + accretion)  *(superseded → MASS doc)*

Superseded by the locked biphasic mechanism. **Retained references:** the $t^\*$
surface-crossing **transient cut** (the transient is excluded from drag
extraction; no literature law exists there) and the **density-driven accretion
rate** $\\dot M\\propto\\rho\_\\text{He}$ — both now live in the MASS pickup channel
(MASS §4) and the §6.7 transient free-zone (widened to several ps there, MASS
A7/R9).

### 2.4 Scenario B (coherent shell, surface stripping)  *(superseded → MASS doc)*

Superseded. Retained as a comparison/regression baseline only.

### 2.5 Biphasic (partial strip + equilibration)  *(promoted → MASS doc)*

The biphasic structure is the **production mechanism**, promoted and made precise
(energy-gated, RRK rate-limited) in the MASS doc. Detail there.

### 2.6 Sub-decision selection  *(superseded → MASS doc)*

Production default is `biphasic\_energy\_gated` (MASS doc); `fixed`/`A`/`B` retained
as comparison baselines. The earlier attach-rate-drop hint (0.09→0.005) is recorded
against this subsection and folded into the §6 validation notes.

### 2.7 Velocity scaling within pickup  *(superseded → MASS doc)*

Density-only $\\lambda\_\\text{attach}\\propto\\rho\_\\text{He}$ confirmed as primary
(resting-ion velocity scaling not required); pickup $v$-dependence kept pluggable
(MASS R7). Detail: MASS §4 / §10.

### 2.8 Interchangeability surface

`SimConfig.mass\_scenario` — production `biphasic\_energy\_gated`; `fixed` and
`anchored\_discrete` (Tier 1a) as comparison/regression baselines. **`scenario\_A\_accretion`
and `scenario\_B\_stripping` are retired** (2026-06-23/24; superseded by the locked
`biphasic` mechanism, §2.5) — the literal is now
`{fixed, biphasic, anchored\_discrete}`. Field-level config (pickup rate,
$f\_\\text{ret}$, $f\_\\text{int}$, $\\nu$, $s$, ladder shape, electronic picture, …):
MASS doc §11.

### 2.9 Schema and energy-bookkeeping changes  *(retained — shared infrastructure)*

**`IonCheckpoint` schema bump to v6** under any non-`fixed` scenario: rename
`E\_mass\_attach\_defect\_eV`→`E\_mass\_transfer\_eV` (same `(2N,T)` shape, sign now
covers accretion and stripping); drop the `mass\_history\_kg` monotonicity
guarantee; add a scenario-metadata field. **Energy invariant** under continuous
mass dynamics:
$$E\_\\text{kin}+E\_\\text{pot}+E\_\\text{dissip}+E\_\\text{mass\_transfer}\\approx\\text{const}$$
(modulo Verlet drift). **The MASS model extends this to a five-term invariant**
by the $E\_\\text{int}$ reservoir:
$E\_\\text{kin}+E\_\\text{pot}+E\_\\text{dissip}+E\_\\text{mass\_transfer}+E\_\\text{int}=\\text{const}$
(MASS §6); cold shedding is energy-neutral for $E\_\\text{solv.struct}$ (no K1/K2
double-count).

### 2.10 Validation criterion specific to mass

Mass affects **both moments**: the first moment via $F=m(t)a$ (TDDFT distance/
velocity traces, **Tier 1**, shell trajectory \~21→19→14 He), and the
**distribution** — the terminal I⁺Heₙ size distribution at the detector is *the*
discriminating observable (**Tier 2**). See §6.4.



## 3\. Drag functional form — analytic vs. tabulated

> \*\*Empirical finding (2026, post-extraction) — `power\_law` exponent is
> `n ≈ +2`, not `n ≈ −2`.\*\* Throughout §1–§3 the power-law form was
> reasoned about under the anticipated $n \\approx -2$ (a hard-sphere
> $\\sigma \\propto v^{-2}$ artifact, singular at $v\\to 0$). The actual
> exported extraction gives $n \\approx +2$ for both cases (9 Å and 18 Å
> similar). Consequences, propagated to the affected passages below:
> - \*\*The drag is regular at $v=0$, not singular.\*\* $n>0 \\Rightarrow
>   F\_\\text{drag}\\to 0$ and $\\gamma(v) = \\gamma\_\\text{PL} v^{\\,n-1}\\to 0$
>   as $v\\to0$. The $\\gamma\_0\\to\\infty$ objection that justified
>   discarding N1 (§1.2) and the low-$v$ regularisation requirement
>   (§3.8 `drag\_low\_v\_floor`) \*\*do not apply to the real coefficients\*\*.
>   `drag\_low\_v\_floor` is retained \*architecturally\* for a hypothetical
>   $n<0$ re-extraction but is \*\*inert for the in-hand export\*\*.
> - \*\*An $n\\approx+2$ wing is the inertial/form-drag $\\sim v^2$ behaviour
>   §3.3 named as the generic physical expectation\*\* — so the power-law
>   form, far from being the discardable artifact, now \*coincides\* with
>   the `linear\_quadratic` high-$v$ wing. This is a Tier-0 cross-check
>   result, not a prior.
> - \*\*The "$n\\approx-2$ is a hard-sphere artifact" hypothesis (§3.3, §3.4)
>   is partly resolved:\*\* the extraction does \*not\* reproduce the old
>   $\\sigma\\propto v^{-2}$ scaling, weakening the artifact concern. Left
>   recorded rather than deleted, since the surrounding hypothesis-framing
>   was the reasoning that motivated keeping the form interchangeable.
> - \*\*No Slice 1 impact:\*\* `power\_law` is deferred regardless of sign, and
>   Slice 1's low-$v$ "contrast against divergence" test is asserted
>   against the \*hypothetical\* $n<0$ form, not the real export.
>
> The original $n\\approx-2$ reasoning is left in place below as the
> recorded prior; read every "$n<0$ / singular at rest / needs a floor"
> claim about `power\_law` as \*\*conditional on a sign the real export does
> not have.\*\*

> \*\*Status note (2026-06-12) — the form cross-check is now scheduled, via
> trajectory matching, not Method-A fit passes.\*\* The shared-form joint
> refit (METHOD\_B §9.7) passed with an effectively \*\*pure-cubic\*\* law
> (`a → 0`, `γ = g·b·v²`), and the next phase (METHOD\_B §10) realizes
> `linear\_quadratic` (incl. pure-quadratic `a ≡ 0`) and `power\_law`
> (free `n`, bounded `n ≥ 1`) and fits each family with the same
> shared trajectory-matching machinery against the pure-cubic incumbent —
> superseding the §3.7 "outstanding `linear\_quadratic` fit pass" route.
> `threshold` stays reserved (not in the §10 scope). The motivating
> tension: the Method-A power-law export's $n\\approx+2$ vs the
> trajectory-matched pure-cubic $n=3$.

### 3.1 Physical question

How is $F\_\\text{drag}(v)$ represented in the integrator? The extraction
pipeline (`drag\_calculation.py`) produces both a closed-form analytic
fit and the raw force-balance scatter $(v, F\_\\text{drag})$ on the
trusted interior, so three representations are available at zero extra
extraction cost: the power-law fit (`fit\_variant=1`), the
linear+cubic fit (`fit\_variant=2`), and a tabulated/interpolated form
sampled directly from the scatter. Both analytic variants have already
been extracted for both the 9 Å and 18 Å cases.

The choice matters because the drag form is consumed in three distinct
places, each with different requirements:

* the integrator acceleration term (§4) — wants cheap evaluation;
* the N2 noise amplitude (§1.2), which needs an analytic
$\\gamma(v) = |F\_\\text{drag}(v)| / v$ — wants a closed form,
ideally finite at $v=0$;
* the low-velocity regulariser (§3.8) — only triggered by a form that is
singular at $v=0$.

### 3.2 Dimensional analysis

All candidate forms must produce a force in
$\\text{amu}\\cdot\\text{Å/ps}^2$ from a speed in $\\text{Å/ps}$. Each
coefficient must absorb whatever power of velocity it multiplies.

**Power law** $;|F\_\\text{drag}| = \\gamma,v^{,n}$:
$$\[\\gamma],(\\text{Å/ps})^{n} = \\text{amu}\\cdot\\text{Å/ps}^2
;\\Rightarrow; \[\\gamma] = \\text{amu}\\cdot\\text{Å}^{,1-n}\\cdot\\text{ps}^{,n-2}.$$
For the extracted $n = -2$: $\[\\gamma] = \\text{amu}\\cdot\\text{Å}^{3}
\\cdot\\text{ps}^{-4}$. Balances. Singular at $v\\to 0$ ($n<0
\\Rightarrow F\\to\\infty$).

**Linear + cubic** $;F\_\\text{drag} = a,v + b,v^3$:
$$\[a],\\text{Å/ps} = \\text{amu}\\cdot\\text{Å/ps}^2
;\\Rightarrow; \[a] = \\text{amu/ps},$$
$$\[b],(\\text{Å/ps})^3 = \\text{amu}\\cdot\\text{Å/ps}^2
;\\Rightarrow; \[b] = \\text{amu}\\cdot\\text{ps}\\cdot\\text{Å}^{-2}.$$
Matches the extraction doc's stated units ($a$ in amu/ps, $b$ in
amu·ps/Å²). Balances. Regular at $v=0$ ($F\\to 0$).

**Linear + quadratic** $;F\_\\text{drag} = a,v + c,v,|v|$ (the
$v|v|$ keeps the term odd, so it is sign-correct/dissipative for both
signs of $v$):
$$\[c],(\\text{Å/ps})^2 = \\text{amu}\\cdot\\text{Å/ps}^2
;\\Rightarrow; \[c] = \\text{amu}\\cdot\\text{Å}^{-1}.$$
Balances. Regular at $v=0$ ($F\\to 0$). High-$v$ wing $\\sim v^2$
(inertial / form drag) is physically gentler than the cubic.

**Threshold / saturating** $;F\_\\text{drag} = F\_\\text{sat},
\\tanh(v/v\_0)$:
$$\[F\_\\text{sat}] = \\text{amu}\\cdot\\text{Å/ps}^2,\\qquad
\[v\_0] = \\text{Å/ps}.$$
$\\tanh$ is dimensionless. Balances. Regular at $v=0$
($F\\to F\_\\text{sat},v/v\_0$, i.e. linear with effective
$a\_\\text{eff} = F\_\\text{sat}/v\_0$ in amu/ps). Drag *saturates* at
$F\_\\text{sat}$ for $v\\gg v\_0$ — the only candidate whose force is
*bounded* at high $v$.

**Tabulated** $;F\_\\text{drag} = \\text{interp}(v;{v\_i, F\_i})$:
dimensionless interpolation over dimensioned samples; result carries
the units of the stored $F\_i$, i.e. $\\text{amu}\\cdot\\text{Å/ps}^2$.
Balances trivially. Behaviour at $v\\to 0$ and beyond $v\_\\text{max}$ is
an extrapolation-policy choice, not fixed by the form.

No formulation is rejected on dimensional grounds.

### 3.3 Physical trade-offs

**Analytic vs. tabulated — degrees of freedom.** Going analytic →
tabulated *gains* fidelity to the extracted shape inside the
extraction window but *loses* a clean analytic $\\gamma(v)$ for the FDT
noise and a clean parametric uncertainty band. The extraction window
is bounded in velocity ($\[t^\*, t\_\\text{end}]$ maps to a finite speed
range), so a table must extrapolate at both ends — relocating the
arbitrariness from "functional form" to "extrapolation rule" rather
than removing it. A table also forces the noise amplitude to
re-differentiate $\\gamma(v)$ numerically, reintroducing precisely the
noise-amplification problem the extraction spline exists to suppress.

**Power-law vs. linear+cubic — degrees of freedom.** Power-law → l+c
*gains* $v=0$ regularity and a finite $\\gamma\_0 = a$ (which
re-enables a finite-amplitude thermal floor and resurrects the
physical content of the discarded N1 noise as a low-$v$ limit), and
*gains* an analytic derivative $dF/dv = a + 3bv^2$ for §4 operator
splitting and any implicit-step Jacobian. It *loses* the single-power
simplicity and the direct match to the extraction's historical default
(`fit\_variant=1`).

**The form is an extrapolation, not an interpolation — and no
ion-in-superfluid drag law is known.** The extraction window is bounded
in velocity, but the ion operates (and especially the bubble-exit
transient, §2.3) *outside* that window. No literature law exists for the
velocity dependence of drag on an ion travelling through superfluid
helium, so every analytic form here is a *hypothesis* to be
cross-checked against the TDDFT curves, not a fitted truth. The four
analytic forms are deliberately chosen to *disagree most in the
extrapolation regime* — their high-$v$ wings span saturating
($\\tanh$), quadratic ($v|v|$), cubic ($v^3$), and decaying
($v^{-2}$ power law) behaviour. Comparing their simulation outputs is
the point of keeping them interchangeable; the form that best matches
the references *is* a simulation result, not an a-priori choice. This
reframes the per-form discussion below: "primary" means "first
hypothesis to run," not "believed correct."

> \*\*Tier-0 finding (2026-06-04) — the dominant in-window error was \*frame\*,
> not \*form\*.\*\* The Tier-0 cross-check (`TIER0\_FINDINGS.md`) found that the
> 9 Å lab-frame mismatch is \*\*not\*\* discriminated by the form set above:
> `linear\_quadratic`'s $v^2$ wing would not fix it, because the error is an
> \*over-large effective friction from using lab speed\* on a COM-drifting
> trajectory, not a high-$v$ wing-shape problem. So before the form
> cross-check can mean anything, the \*\*extraction frame\*\* (lab vs. ion–He
> relative velocity) must be corrected (`EXTRACTION\_FRAME\_FIX\_milestone.md`).
> Until then the in-hand coefficients are provisional and the "which form
> best matches" question is premature — a relative-velocity re-extraction may
> make a single shared `linear\_cubic` fit both cases (§3.6 headline), in which
> case no form switch is warranted. The form set stays interchangeable as
> designed; the point is only that the \*first\* discriminating error found was
> upstream of the form.
>
> \*\*CORRECTION (2026-06-04, same day, superseding the note above).\*\* The
> "frame, not form" reading was itself \*\*substantially wrong\*\* and is
> withdrawn. Two facts overturned it: (1) the COM that the relative-velocity
> story rested on was `½(v₁+v₂)` — but atom 1's sideways motion is a \*TDDFT
> artifact\* (discarded at extraction; the law was fit to the single clean
> atom 2, not an average), so a molecular COM was never in the extraction and
> one of its inputs is distrusted; (2) the real 9 Å mismatch is a \*\*windowing
> + bubble-mode\*\* effect, not a frame rotation — the extraction window runs
> past \~6 ps where atom 2 develops a directional drift a central-force MD
> cannot represent, and the residual in-window RMSE (\~0.86) is dominated by
> the 1.2 ps bubble oscillation present in the \*raw\* reference but removed by
> the extraction's denoising. Re-extracting on the clean short window barely
> moved the RMSE (0.86→0.88), confirming the \*form/fit\* is fine. So: the form
> set conclusion stands (no form switch warranted), the coefficients are \*\*not\*\*
> frame-provisional in the way claimed, and `EXTRACTION\_FRAME\_FIX\_milestone.md`
> is \*\*demoted to a contingency\*\* (see `TIER0\_FINDINGS.md` and
> `tier0\_comparison\_tasks\_left.md`). The genuine residual question is settled by
> comparing against the \*same-smoothed\* reference (an internal-consistency
> check), not by a frame re-extraction.

**Asymptotic signatures (for orientation, not selection):** linear
Stokes-like drag at low $v$; an inertial/form-drag $\\sim v^2$ wing at
high $v$ is the generic expectation when an object sheds fluid; a
critical-velocity (Landau) threshold is the superfluid-specific
possibility hinted at by the baseline's existing `v\_limit`/`E\_min`
cutoff. The power-law $n\\approx-2$ is most likely a fitting artifact of
the old hard-sphere $\\sigma\\propto v^{-2}$ cross-section bleeding into
the extraction window rather than a transport law — retained as a
hypothesis precisely so the cross-check can confirm or reject that.

**Dissipativity guard (physical validity).** For $F\_\\text{drag}$ to
oppose motion at all relevant speeds, each form has a validity
condition checked at config-load (the extraction's `curve\_fit` carries
no bounds):

* **linear+cubic:** requires $a + b,v^2 \\ge 0$ across the operating
range. $a\\ge0, b>0$ is monotone and clean. **The locked production law
has $a=0$** (pure cubic, §3.4): $F=bv^3>0$ and $dF/dv=3bv^2\\ge0$ for
all $v\\ge0$ — strictly dissipative and monotone, with the turnover
$v\_\\dagger=\\sqrt{-a/b}$ **undefined** (guard trivially satisfied).
Assert $a\\ge0$ and $b>0$. (A future $b<0$ re-extraction would
reintroduce a turnover $v\_\\dagger$ to assert above the maximum
trajectory speed.)
* **linear+quadratic:** requires $a + c,|v| \\ge 0$; with $a,c>0$ this
holds for all $v$ with no turnover. Assert $a>0, c\\ge0$.
* **threshold:** $F\_\\text{sat}/v\_0 > 0$ (i.e. both $F\_\\text{sat}>0$ and
$v\_0>0$) guarantees dissipative and monotone everywhere; bounded
force means no high-$v$ stiffness.
* **power law:** $\\gamma>0$; the $v=0$ singularity for $n<0$ is deferred
to the §3.8 regularisation note rather than guarded here.

The need for a "maximum trajectory speed" to test $v\_\\dagger$ against
is itself an open item — it has no home in the current baseline config
and must be sourced from the TDDFT trajectories or set as a generous
ceiling (flagged for §6 validation/tolerances).

### 3.4 The analytic form set — four swappable hypotheses

All four analytic forms are first-class and interchangeable behind the
`drag\_form` enum (§3.8). Because no ion-in-superfluid drag law is known
(§3.3), "primary" denotes the *first hypothesis to run and the default*,
not a belief about correctness; the empirical cross-check against the
TDDFT references selects among them.

**Primary (LOCKED 2026-06-15) — pure cubic**, the $a=0$ instance of
linear+cubic: $;F\_\\text{drag} = g(\\text{depth}),b,v^3$, with
$b = 2.5154\\ \\text{amu·ps/Å}^2$ and $g$ the §5 gate. Tier-0
trajectory-matching (Method B, $n=3$) was selected over Method A
($n\\approx2$) and drove the linear coefficient to zero, so the linear
term is **absent** in the production law, not merely small. Regular at
$v=0$ ($F\\to0$ smoothly), analytically differentiable
($dF/dv = 3bv^2$), shared across both cases (§3.6). Coefficients are in
hand for 9 Å and 18 Å. **Two consequences of $a=0$, recorded:**
(i) $\\gamma\_0 \\equiv \\lim\_{v\\to0}\\gamma(v) = 0$ — *not* the finite $a$
the §1.2 N2 argument assumed; the FDT noise amplitude vanishes at rest
(addressed in §1.2's correction note, and accepted). (ii) The "$v^3$
wing extrapolates more steeply than the physically-generic quadratic"
caveat is now **load-bearing**, not cosmetic — beyond the fit window the
cubic over-brakes; see R10 (high-velocity ceiling) in
`MASS\_DYNAMICS\_LOCKED\_energy\_gated\_evaporation.md` §8. The $b<0$
turnover risk does **not** apply ($b>0$, no real $v\_\\dagger$).

*Note — linear+quadratic is the stronger physical default and may
supersede this on cross-check.* On both integrator-neutrality and
asymptotics grounds, linear+quadratic is arguably the better primary:
its high-$v$ wing is the generic inertial $\\sim v^2$ rather than the
unmotivated $v^3$, it has no turnover risk, and its $\\gamma(v)$ grows
only linearly (milder stiffness). It is kept *secondary* here only
because linear+cubic is the form already extracted and chosen; if the
empirical cross-check does not favour linear+cubic, linear+quadratic is
the first alternative to promote. Flagged so the choice is revisited
with data rather than left implicit.

**Secondary — linear+quadratic** $;F\_\\text{drag} = a,v + c,v,|v|$.
Regular and dissipative at $v=0$ with no turnover (monotone for
$a,c>0$), inertial $\\sim v^2$ high-$v$ wing, $\\gamma(v) = a +
c,|v|$ growing only linearly. Finite $\\gamma\_0 = a$. Needs a fit
pass (not yet extracted). See the promotion note above.

**Secondary — threshold/saturating** $;F\_\\text{drag} =
F\_\\text{sat}\\tanh(v/v\_0)$. The superfluid-physics-anchored hypothesis:
linear at low $v$ (effective $\\gamma\_0 = F\_\\text{sat}/v\_0$ in amu/ps),
*bounded* force at high $v$, and a built-in velocity scale $v\_0$ that
can represent a critical-velocity-like onset. $\\gamma(v)$ *decays* at
high $v$ — essentially zero stiffness in the operating regime, the
friendliest of all forms for any integrator. Cost: an extra parameter
and harder calibration from two trajectories. Not yet extracted; needs
a fit pass.

**Secondary — power law** $;|F\_\\text{drag}| = \\gamma,v^{,n}$. The
extraction's historical default form and the "drag falls with speed"
hypothesis. Was expected to be a hard-sphere $\\sigma\\propto v^{-2}$
artifact (§3.3); retained so the cross-check can confirm or reject that.
**The cross-check has run (§3 finding note): the real export is
$n\\approx+2$, not $-2$ — drag *rises* as $\\sim v^2$, the artifact
hypothesis is not borne out, and the form is regular at $v=0$ with no
floor needed.** The "singular at $v=0$ for $n<0$" property and the floor
it re-enabled are therefore conditional on a sign the in-hand
coefficients do not have. Coefficients ${\\gamma,n}$ already extracted
for both cases.

Since all four are closed-form, swapping among them is a change behind
the form enum with no structural impact on the integrator (the chosen
BAOAB integrator is form-agnostic and unconditionally stable for all
four — see §4).

*Consequence for regularisation.* Only `power\_law` ($n<0$) needs a
low-velocity floor; the other three forms are finite at $v=0$ and need
none. This obligation lives in the §3.8 regularisation note
(`drag\_low\_v\_floor`), not a standalone section.

*Consequence for §1.* The N1-discard reasoning in §1.2 is rescoped: the
"$\\gamma\_0\\to\\infty$" objection applies only under `power\_law`; under
the other three, $\\gamma\_0$ is finite and N2 reduces to the correct N1
limit at low $v$.

### 3.5 Discarded — tabulated $F(v)$

Discarded for three reasons, recorded so the option is not
re-litigated:

1. The N2 noise amplitude needs an analytic $\\gamma(v)$; a table forces
numerical re-differentiation and reintroduces the noise-amplification
the extraction spline was built to avoid.
2. The table is faithful only inside the bounded extraction velocity
window; outside it the representation must extrapolate, merely
relocating the arbitrariness rather than removing it.
3. No clean parametric uncertainty band (the 10-seed extraction sweep
maps directly onto ${\\gamma,n}$ or ${a,b}$, not onto a table).

### 3.6 Per-case form policy — shared-form preferred

The two bubble sizes *may* sit in different drag regimes (the 18 Å
droplet has a longer dense-He traversal and, in the old model, a
different cross-section and binding). The architecture therefore
allows the drag *form*, not just its coefficients, to differ per case:
`drag\_form` is a per-preset value.

The *default policy*, however, is shared-form-preferred:

* Use the primary form (`linear\_cubic`) for both cases. Report the
shared-form fit as the headline result.
* A shared functional form fitting both cases with only the
coefficients moving is weak evidence the law captures real transport
physics rather than per-case curve-fitting — this cross-case
consistency signal is worth preserving as the anchor result.
* Diverge to a different form for one case *only* if it trips the
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

The primary form is `linear\_cubic` (extraction `fit\_variant=2`), but
`Drag\_extraction\_code.md` documents `fit\_variant=1` (power law) as the
extraction's *default* output. Both variants have now been extracted
for both cases, so no pipeline change is required — but the MD side must
**consume the `fit\_variant=2` (linear+cubic) coefficients as the
default**, and the extraction's default-output setting should be
aligned (or the linear+cubic product explicitly carried forward) so
that the primary drag form is the one actually plumbed into the
simulation. The `linear\_quadratic` and `threshold` forms (§3.4) are
*not yet extracted* and require their own fit passes before they can be
cross-checked. Recorded as a cross-document consumption dependency plus
two outstanding fit passes, not a re-extraction of existing forms.

### 3.8 Interchangeability surface

* `SimConfig.drag\_form ∈ {linear\_cubic, linear\_quadratic, power\_law, threshold}` — default/primary `linear\_cubic` (data in hand);
`linear\_quadratic` flagged as the likely physical supersedor;
`threshold` the superfluid-anchored hypothesis; `power\_law` is the
only form requiring a low-velocity floor (see the regularisation note
below). Per-preset, enabling per-case forms.
* `SimConfig.drag\_coefficients` — form-tagged coefficient bundle of
variable arity by form, per-case:

  * `linear\_cubic` (LOCKED, pure cubic): ${a=0,\\ b=2.5154}$, units
amu/ps and amu·ps/Å².
  * `linear\_quadratic`: ${a, c}$, units amu/ps and amu/Å.
  * `power\_law`: ${\\gamma, n}$, $\\gamma$ in
amu·Å$^{1-n}$·ps$^{n-2}$, $n$ dimensionless.
  * `threshold`: ${F\_\\text{sat}, v\_0}$, units amu·Å/ps² and Å/ps.
* **Drag-validity guard (config-load):** per-form dissipativity check
from §3.3 — `linear\_cubic` ($a\\ge0$, $b>0$; with the locked $a=0$ the
turnover $v\_\\dagger$ is undefined and the guard is trivially met);
`linear\_quadratic` ($a>0, c\\ge0$); `threshold` ($F\_\\text{sat}>0,
v\_0>0$); `power\_law` ($\\gamma>0$). The "maximum trajectory speed"
needed for the `linear\_cubic` turnover check is an open config item
(flagged for §6 validation/tolerances).
* **Low-velocity regularisation (only `power\_law`, and only if $n<0$).**
The three finite-at-zero forms (`linear\_cubic`, `linear\_quadratic`,
`threshold`) need no regulariser — $F\_\\text{drag}\\to0$ and
$\\gamma(v)$ stays finite as $v\\to0$. Only `power\_law` with $n<0$
diverges; it then requires a floor on $\\gamma(v)$ (equivalently a
small $v\_\\text{floor}$ below which $\\gamma$ is held constant) to keep
both the drag and the FDT noise amplitude finite. **The real extracted
exponent is $n\\approx+2$ (§3 finding note), which is regular at $v=0$
and needs no floor — so `drag\_low\_v\_floor` is inert for the in-hand
coefficients and is retained only for a hypothetical $n<0$
re-extraction.** The chosen BAOAB integrator (§4) partially defuses any
divergence on its own — the damping factor $e^{-\\gamma,dt/m}\\in\[0,1]$
stays bounded even as $\\gamma\\to\\infty$ — but a divergent noise
amplitude $\\sqrt{2\\gamma k\_BT\_\\text{eff}}$ would still need the floor.
Exposed as `SimConfig.drag\_low\_v\_floor` (Å/ps), active only when
`drag\_form = power\_law` *and* $n<0$; ignored otherwise. This absorbs
what was previously a standalone low-velocity-regularisation section.
* Cross-references: $\\gamma(v) = |F\_\\text{drag}|/v$ (units amu/ps)
feeds the §1.2 N2 noise amplitude; $m$ here is the physical $m(t)$ per §2.2, and the
drag is applied as $-\\gamma(v)v$ independent of $m(t)$.

\---

## 4\. Integrator coupling — additive vs. operator-split

### 4.1 Physical question

How do the drag force (§3) and, once active, the Langevin noise (§1)
enter the time-stepping relative to the baseline velocity-Verlet
integrator? The baseline `velocity\_verlet\_step` (`leapfrog.py:74`) is
kick-drift-kick, second-order, and model-agnostic: it consumes any
`acc\_fn(pos) → (acc, E\_pot)`. The deeper issue is that this clean
structure assumes a *position-dependent* force $F = F(x)$, whereas drag
is $F\_\\text{drag}(v)$ and noise makes the system a stochastic
differential equation. Both break the assumptions velocity-Verlet
rests on.

### 4.2 Why drag breaks the velocity-Verlet assumption

Velocity-Verlet's second-order accuracy and clean error structure rely
on the force being a function of position only. A velocity-dependent
force has three consequences:

* The second kick $v\_1 = v\_0 + \\tfrac12(a\_0 + a\_1),dt$ becomes
*implicit*: $a\_1 = a\_1(x\_1, v\_1)$ depends on the $v\_1$ being solved
for. $v\_1$ appears on both sides.
* Explicit (naïve) treatment of a dissipative $F(v)$ term is only
*conditionally* stable: stable when $\\gamma,dt/m \\ll 1$, divergent
when drag is stiff. Stiffness is form-dependent (§3): for
`linear\_cubic`/`linear\_quadratic`/`threshold` over the working range,
$\\gamma,dt/m \\sim 10^{-4}$ (very non-stiff); for `power\_law` with
$n=-2$, $\\gamma\\to\\infty$ as $v\\to0$ and the term becomes arbitrarily
stiff at low speed — exactly where the ion ends up at late times.
* With §1 noise, the system is an SDE; standard Verlet has no noise
term and cannot inject one with correct fluctuation–dissipation
statistics at finite $dt$.

Time-reversibility is lost under *any* coupling — dissipation is
irreversible by construction. That is physics, not a defect of the
scheme.

### 4.3 Dimensional check

Drag enters the equation of motion as an acceleration
$a\_\\text{drag} = F\_\\text{drag}(v)/m(t)$:
$$\[a\_\\text{drag}] = \\frac{\\text{amu}\\cdot\\text{Å/ps}^2}{\\text{amu}}
= \\text{Å/ps}^2,$$
matching the existing Coulomb and droplet accelerations. With §2.2,
$m(t)$ is the physical mass and the drag is applied as
$-\\gamma(v),v/m(t)$ with $\\gamma(v)$ in amu/ps. The OU damping factor
below, $e^{-\\gamma,dt/m}$, has a dimensionless exponent:
$\[\\gamma,dt/m] = (\\text{amu/ps})(\\text{ps})/\\text{amu} = 1$. Balances.

### 4.4 The option space

**Option I — naïve additive inside `\_ion\_accel\_fn`.** Add
$-\\gamma(v)v/m$ to the acceleration using whatever velocity is
available at each Verlet sub-evaluation; integrator structure
untouched. Minimal diff. But it treats the velocity-dependent term
explicitly, so stability is conditional and form-dependent (safe for
the finite-$\\gamma$ forms over their working range, unstable for
`power\_law` near $v=0$), and it loses second-order accuracy in the drag
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
Coulomb+droplet via the existing `acc\_fn`) and an
Ornstein–Uhlenbeck step (O, handling drag+noise), composed symmetrically
(half-B, half-A, full-O, half-A, half-B). For *linear* drag the O-step
is exact:
$$v \\mapsto e^{-\\gamma,dt/m},v

* \\sqrt{\\frac{k\_B T\_\\text{eff}}{m}\\left(1 - e^{-2\\gamma,dt/m}\\right)};\\xi,$$
unconditionally stable, exact in the drag part for linear $\\gamma$, and
the canonical site for injecting the §1 Langevin noise with correct FDT
statistics at finite $dt$. For the nonlinear forms, $\\gamma(v)$ varies
within the step; freeze $\\gamma$ at step-entry velocity (an $O(dt)$
error in the nonlinear correction, acceptable since it is a correction
on top of the linear floor) or do a short implicit sub-solve inside O.
The bounded damping factor $e^{-\\gamma,dt/m}\\in\[0,1]$ stays
well-behaved even as $\\gamma\\to\\infty$, so the O-step *partially defuses*
the `power\_law` stiffness (§3.8).

### 4.5 Trade-offs

*Structure.* I keeps the model-agnostic integrator and minimal diff but
loses unconditional stability, a principled noise site, and second-order
drag accuracy. II gains unconditional drag stability and keeps the
conservative forces in the clean callable but breaks the abstraction and
still has no exact-OU noise site. III gains unconditional stability,
exact linear-drag propagation, correct finite-$dt$ FDT noise, and
partial `power\_law` stiffness relief, at the cost of a *new ion-stage
integrator path* and the Stratonovich/Itô bookkeeping (which BAOAB
resolves cleanly — see §1.2).

*Conservation.* The §2.9 energy invariant
$E\_\\text{kin}+E\_\\text{pot}+E\_\\text{dissip}+E\_\\text{mass\_transfer}$
needs per-step dissipated energy. Under I/II this is the drag work
$\\int F\_\\text{drag}\\cdot v,dt \\approx F\_\\text{drag}\\cdot v,dt$. Under
III the O-step yields the dissipated energy analytically from the
velocity damping, and separates it cleanly from the noise energy
injection (which feeds the thermal floor, not $E\_\\text{dissip}$) —
cleaner accounting.

*Coupling.* §1 noise: III is the only option with a native,
statistically-correct noise site; I/II require a separate stochastic
bolt-on with wrong finite-$dt$ statistics. §3 form: I's stability is
form-dependent; II/III are form-robust. §2 mass: $m(t)$ changes per
step, so $\\gamma/m$ in the O-step uses the current $m(t)$ — consistent
with the per-step closure rebuild (`ion\_propagation\_step.py:184-193`),
no conflict. §3.8 floor: III reduces but does not eliminate the `power\_law`
regulariser (a floor is still needed to keep the noise amplitude finite
and avoid $e^{-\\infty}$ pathologies).

### 4.6 Sub-decision selection

**Primary — III, BAOAB operator splitting.** Chosen on the principle
that the integrator must be **form-agnostic and unconditionally stable**
so that no drag-form choice can destabilise it and corrupt the
empirical form cross-check (§3.3): "I changed the form and the
trajectory changed" must be a clean physics statement, not an artifact
of a shifting stability margin. III is unconditionally stable for all
four forms including `power\_law`'s singularity (bounded exponential),
its O-step is exact for the linear floor every form shares, and it is
the only option with a native FDT-correct noise site — which the §1.2
noise decision already implicitly committed to (its repeated references
to BAOAB-family integrators). The cost (a separate ion-stage integrator
path, more implementation work) is accepted as the price of protecting
the comparison.

*Accepted architectural asymmetry.* The neutral stage keeps the
baseline `velocity\_verlet\_step`; only the ion stage gets the BAOAB
path. This is a deliberate fork — the neutral stage has no drag and no
Langevin noise, so there is nothing for BAOAB to buy there, and §14 of
the baseline keeps the neutral integrator off-limits regardless.

**Discarded — I, naïve additive.** Its only virtue is minimal diff,
which evaporates the moment noise is switched on (the default intent).
Form-dependent stability would let `power\_law` runs diverge at low $v$,
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

* The integrator choice is **not** exposed as a routine `SimConfig`
swap. Unlike the force-form and noise enums, the integrator is a
deeper structural object; III is hard-wired for the ion stage and
I/II are recorded as discarded rather than kept pluggable. A
convergence-testing hook may temporarily expose alternatives, but the
production path is single.
* `SimConfig.integrator\_dt` reuses the existing `dt\_ion` (0.01 ps); the
BAOAB step consumes the same timestep.
* The O-step reads `drag\_form`/`drag\_coefficients` (§3.8) to build
$\\gamma(v)$ and the noise surfaces (`noise\_form`,
`noise\_calibration`, `noise\_geometry`, §1) to build the fluctuation
term. All physics choices remain behind their own enums; III is the
fixed *mechanism* that consumes them.
* Per-step dissipated and noise-injected energies are written to the
§2.9 accumulators (`E\_dissip\_eV`, and the thermal-floor channel),
computed analytically from the O-step rather than by finite-difference
drag work.

\---

## 5\. Spatial gating — does drag turn off at the droplet surface?

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
detector velocity the model validates against (`vmi\_iplus\_he.csv`).
Some gating is mandatory; only its shape is open.

### 5.2 The gating factor and its dimensional status

Introduce a dimensionless gating factor $g(\\text{depth}) \\in \[0,1]$
multiplying the drag, with depth $= r\_\\text{atom} - r\_\\text{droplet}$
(negative inside, positive outside) as in the baseline droplet
potential (`potentials.py:73`):
$$F\_\\text{drag,gated}(v,\\text{depth}) = g(\\text{depth})\\cdot
F\_\\text{drag}(v).$$
Because $g$ is dimensionless, the gated force keeps units
$\\text{amu}\\cdot\\text{Å/ps}^2$ for any choice of $g$ — no candidate is
dimensionally excluded; the choice is purely physical.

**The gate applies to drag *and* noise (hard FDT coupling).** The
fluctuation–dissipation relation ties the noise amplitude to the same
friction the gate modulates. The gated noise amplitude is therefore
$$\\sqrt{2,\\gamma(v),g(\\text{depth}),k\_B T\_\\text{eff}},$$
carrying the *same* $g(\\text{depth})$. If the gate applied to drag only,
the ion would receive thermal kicks in vacuum where it feels no
friction — an FDT violation producing spurious heating outside the
droplet. Gating both is mandatory, not optional.

### 5.3 Existing geometry in the baseline

Two notions of "inside" already exist and bound the choices:

* a **sharp boolean** `depth < 0`, used as a hard cutoff by the
hard-sphere collision sampler (baseline §6.2);
* a **smooth erf transition** over `potential\_steepness = 14.2 Å`, used
by the droplet confining potential
$V\_\\text{drop}(r) = \\tfrac12(\\text{erf}(\\text{depth}/\\text{steepness})
+1),E\_b$.

The gating decision largely reduces to which of these drag follows, or
whether it gets its own profile.

### 5.4 Candidate gating forms

**G1 — sharp boolean** $g = \\mathbb{1}\[\\text{depth} < 0]$. Full inside,
zero outside, step at the surface. Maximal consistency with the
collision model being replaced. But a step discontinuity in force
breaks the smoothness the BAOAB O-step (§4) assumes and injects a
spurious impulse as the ion crosses the surface — the same
integrator-artifact concern that drove the §4 choice. Physically crude:
He density does not vanish at a mathematical surface.

**G2 — erf gate tied to the confining potential**
$g(\\text{depth}) = \\tfrac12\\big(1 - \\text{erf}(\\text{depth}/
\\text{steepness})\\big)$, reusing `potential\_steepness`. Smooth,
differentiable, 1 deep inside → 0.5 at the nominal surface → 0 outside,
over the *same* \~14 Å shell where confinement turns on. Zero new
parameters. Physically coherent: drag and confinement are both proxies
for local He, so tying them to one profile is the parsimonious claim.

**G3 — erf gate with independent drag steepness** same shape as G2 but
with a separate `drag\_gate\_steepness`. Admits that the He shell relevant
to *momentum transfer* may have a different effective range than the
one relevant to *binding*. One new parameter — but the extraction window
deliberately excludes the surface-crossing transient ($t^\*$ cut, §2.3),
so there is likely *no clean TDDFT data on the gating profile itself*,
making this extra parameter effectively uncalibratable for now.

**G4 — density-proportional** $g(\\text{depth}) =
\\rho\_\\text{He}(\\text{depth})/\\rho\_\\text{bulk}$. The most physically
honest: drag *is* momentum transfer to He, so it scales with local He
number density. Its decisive advantage is internal consistency — the
*same* $\\rho\_\\text{He}(\\text{depth})$ would gate drag, the FDT noise,
*and* the §2.3 mass-accretion rate $\\dot M \\propto
\\rho\_\\text{He}(\\text{depth})$. One density profile would then drive
every helium-coupling channel in the model.

### 5.5 Sub-decision selection

**Primary — G4, density-proportional, collapsing to G2 for now.**
Chosen as the principled framing: drag, noise, and mass-rate all become
functions of one local He density, removing gating as an independent
modelling choice (it falls out of the density profile). **No
$\\rho\_\\text{He}(\\text{depth})$ profile is currently defined** beyond the
flat interior `DENSITY\_DROPLET = 0.8\\,\\rho\_\\text{bulk}` and the erf
surface, so today G4 is *implemented as* the erf complement (G2),
reusing `potential\_steepness`. The two are operationally identical until
a measured or computed density profile exists. When one does, G4 becomes
genuinely distinct without an architecture change — a different
$\\rho\_\\text{He}(\\text{depth})$ feeds the same gate. **Flagged as a
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
is already baked into the extracted $\\gamma(v)$. Spatial gating
therefore handles **only** the macroscopic droplet-surface cutoff
(am-I-still-in-the-droplet), **not** a near-field bubble density hole.
Recorded explicitly so a future contributor does not add a near-field
density term and double-count the bubble physics already inside the
drag law.

> \*\*Tier-0 finding (2026-06-04) — partial contradiction, recorded.\*\* §5.6
> assumed the bubble/relative physics was \*correctly\* absorbed into
> $\\gamma(v)$. The Tier-0 comparison (`TIER0\_FINDINGS.md`) shows the
> extraction absorbed it \*\*incorrectly when COM drift is present\*\*: $\\gamma$
> was fit against \*\*lab\*\* speed, so on a trajectory where the ion co-drifts
> with the He (9 Å: \~4 Å/ps COM drift), $\\gamma(|v\_\\text{lab}|)$ over-damps
> because lab speed $\\gg$ relative speed. The "it's baked into $\\gamma$"
> claim holds \*\*only for a relative-velocity-correct extraction\*\*, which the
> in-hand lab-speed coefficients are not. The fix is a relative-velocity
> re-extraction (`EXTRACTION\_FRAME\_FIX\_milestone.md`), \*not\* adding a
> near-field density term — so the §5.6 "do not double-count" guidance still
> stands; what changes is that the current $\\gamma$ is provisional until
> re-extracted on $|v\_\\text{rel}|$. Whether the \*simulation\* must then track a
> local $v\_\\text{He}$ to form $v\_\\text{rel}$ at run time (vs. a purely
> extraction-side correction) is the open question deferred in the milestone.
>
> \*\*CORRECTION (2026-06-04, same day).\*\* The frame reading above was
> superseded the same day — see the correction under §3.3 and
> `TIER0\_FINDINGS.md`. The 9 Å mismatch was traced to \*\*windowing + bubble-mode
> oscillation\*\*, not a lab-vs-relative frame error built on the (artifactual)
> atom-1 COM. §5.6's original "no near-field term, it's baked into γ" guidance
> therefore \*\*stands unmodified\*\* — there is no demonstrated frame contamination
> requiring a near-field term, and `EXTRACTION\_FRAME\_FIX\_milestone.md` is
> demoted to a contingency. This super-annotation is retained for the audit
> trail (we hypothesised frame, then found windowing).

### 5.7 Interchangeability surface

* `SimConfig.drag\_spatial\_gate ∈ {density\_proportional, erf\_tied, erf\_independent, sharp}` — primary `density\_proportional` (G4),
implemented as the erf complement until a density profile exists;
`erf\_tied` (G2) the explicit-erf secondary; `erf\_independent` (G3)
the separate-steepness secondary; `sharp` (G1) discarded, retained in
the enum only for regression comparison against the old model.
* `SimConfig.drag\_gate\_steepness` — used only by `erf\_independent`;
defaults to `potential\_steepness` (14.2 Å) so it reduces to G2 when
unset.
* `SimConfig.helium\_density\_profile` — placeholder for the future
$\\rho\_\\text{He}(\\text{depth})$ that promotes G4 beyond the erf
collapse; until defined, `density\_proportional` uses the erf
complement internally.
* The gate multiplies both the deterministic drag and the FDT noise
amplitude (§5.2); the BAOAB O-step (§4) reads $g(\\text{depth})$ at
each step alongside $\\gamma(v)$.

\---

## 6\. Validation surface and tolerances

### 6.1 Why this section exists

Every "primary" in §1–§5 is deferred to empirical cross-check against the TDDFT
and experimental references; §6 is the instrument that resolves those deferrals
and the **order** in which it does so. Three jobs: extend the baseline surface to
observables it does not exercise (§6.2), order the comparisons so the entangled
parameter space stays separable (§6.3–§6.4), and enforce the two consistency
constraints coupling drag to {mass, binding} (§6.5–§6.5.1). The per-parameter map
of *which tier anchors what* lives in `CALIBRATION\_MAP.md` ("anchor coverage by
tier"); this section is the rationale and current status.

### 6.2 What the baseline surface cannot do

The baseline (PHYSICS\_BASELINE §13) ships single-trajectory comparison, VMI
final-velocity histograms, and energy-balance plots — built for a deterministic
model matched to single references. The drag port breaks three assumptions, each
needing a new hook: noise lives in **second moments** (§1.6, single-trajectory
compare is blind to spread); the mass model is discriminated only by the terminal
**I⁺Heₙ size distribution** (a discrete-$n$ distribution the code does not yet
compare); and the form cross-check needs the comparison **clean** (the BAOAB
choice §4.6 only pays off if validation isolates the form). New hooks required: an
ensemble-variance comparison (Tier 3) and a discrete-$n$ size-distribution routine
(Tier 2).

### 6.3 The attribution problem → sequential validation

The unknowns interact (drag form, noise, mass model, gate) and the observables are
coupled — final velocity depends on form *and* mass *and* gate; its spread on noise
*and* ensemble. A flat "run everything, compare to VMI" cannot separate them.
Resolution: a **hierarchy ordered by separability** — each tier isolates as few
unknowns as possible and fixes its winner before the next is introduced.
Validation is **sequential, not simultaneous**.

### 6.4 The sequential validation hierarchy — status and goals

Current state (2026-06-15): **Tier 0 complete · Tier 1 ungated and next · Tiers 2–3
downstream.** Full records: `TIER0\_FINDINGS.md`,
`METHOD\_B\_trajectory\_matching\_extraction.md`. Per-parameter tier anchors:
`CALIBRATION\_MAP.md`.

**Tier 0 — drag form, deterministic, fixed mass (COMPLETE).** Method B
trajectory-matching settled the drag law: `shared\_pure\_cubic` ($a=0$,
$b=2.5154$ amu·ps/Å², effective binding $E\_\\text{bind}=0.1168$ eV) — production;
`power\_law` with free exponent independently recovered $n\\approx3$;
`linear\_quadratic` ($n{=}2$) rejected by objective and held-out 9 Å. **18 Å is
the clean-radial regression floor** (distance RMSE ≤ 3.0 Å, mean $|v|$ RMSE ≤ 0.25
Å/ps); **9 Å carries the transverse-contamination flag** (a model-dimensionality
limit, not a drag error). Two findings carry forward: that 9 Å flag, and "correct
drag traps the ions" → the binding becomes an effective jointly-calibrated
parameter (§6.5.1, with the OQ1 provenance flag there).

**Tier 1 — mass scenario, deterministic (UNGATED — next). Goal: *which mass
model*, with no noise and no Tier-2 commitment.** With the Tier-0 form fixed and
noise off, run the locked `biphasic\_energy\_gated` mechanism (and `fixed` as
baselines) and compare the full post-transient trajectory **and the time-resolved
shell trajectory** (\~21→\~19→\~14 He, §2.1) against the 9 Å TDDFT reference.
Isolates the **mass model** (§2 → MASS doc). **This tier is OQ-independent** — the
shell trajectory does not depend on the electronic picture or ladder *values* — so
it can run now. The transient free-zone (§6.7) lives here; each non-`fixed`
scenario needs its own coefficient extraction (§6.5) before its comparison carries
meaning. The deterministic invariant/consistency checks (five-term energy
conservation, cold-shed neutrality, no gate-open avalanche; MASS §6) are the
build's correctness gate and hold **regardless of OQ values** — the cheap insurance
to build first.

**Tier 2 — terminal size distribution. Goal: *arbitrate the genuinely-free knobs*.**
Compare the simulated I⁺Heₙ size distribution (discrete integer-$n$) against the
experimental detector histogram. This is the *only* observable that sharply
separates the mass scenarios and the **two genuinely-free knobs — the ladder shape
and the electronic picture** (`CALIBRATION\_MAP` tally; MASS R3/A10) — plus the
bounded early-window scalars ($f\_\\text{int}, f\_\\text{ret}, \\tau\_\\text{dissip}, \\nu$).
It also evaluates the **total-stripping limit** (Calvo24, secondary; MASS §6.11,
OQ6). *Load-bearing caveat:* Tier 2 carries **8+ quantities on one observable**;
whether the size distribution actually separates them is the central calibration
risk (documented per-parameter in `CALIBRATION\_MAP`). **OQ-gated:** *committing*
the ladder/binding calibration here depends on **OQ1** (the drag $E\_\\text{bind}$ vs
mixture-ladder double-count, §6.5.1 / MASS A10/OQ6) — build and run Tier 2 freely,
but do **not lock** a Tier-2 ladder fit while OQ1 is open.

**Tier 3 — ensemble second moments.** With noise on, compare the final-velocity
histogram **width** and (if T3 geometry, §1.4) angular spread against VMI, plus
cross-trajectory variance on the 8000-atom `single\_pulse\_droplet\_distribution`.
Isolates the **noise model** (§1) given form + mass + gate from Tiers 0–2.

*Numeric acceptance thresholds are deferred:* §6 fixes the metrics (§6.9) and the
order; pass/fail numbers are set once first runs are seen, since everything
upstream is cross-checked empirically. The one committed floor is the 18 Å Tier-0
regression above.

### 6.5 Consistency constraint I — mass scenario ↔ drag coefficients (R6)

The force balance $F\_\\text{drag}(t)=m(t),a(t)-F\_C(R(t))$ pins only the
combination $m(t)a(t)$, so what is attributed to drag depends on the assumed
$m(t)$: a law extracted at constant $m\_\\text{eff}$ is self-consistent **only**
re-applied at $m\_\\text{eff}$. `mass\_scenario` and `drag\_coefficients` are a
**coupled pair**, enforced as a config-load guard — coefficient bundles stamped
with `extraction\_mass\_model ∈ {constant, time\_resolved}` (an extraction-side
stamping action, §6.6); `fixed`@$M$ needs `constant`@$M$ within \~1–2 He;
non-`fixed` needs `time\_resolved` under a matching $m(t)$; a mismatched pairing
**refuses to run**, with `SimConfig.allow\_inconsistent\_mass\_pairing` (default
`False`) downgrading the refusal to a loud warning for exploration. **Production
status:** `biphasic\_energy\_gated` runs a non-monotone $m(t)$ and therefore trips
this guard structurally — see MASS **R6** for the production handling (the §6.6
mid-window defence makes the inconsistent pairing defensible; the clean
re-extraction is blocked by the 9 Å transverse flag).

### 6.5.1 Consistency constraint II — drag ↔ droplet binding (OQ1)

"Correct drag traps the ions": the in-window-correct drag delivers sub-barrier
surface KE, so the **static** 0.308 eV solvation well prevents the ejection that
the TD-HeDFT ions achieve **dynamically** (the He reorganises; the static barrier
is bypassed). This is not a drag error and not a re-measured solvation energy —
the binding depth becomes an **effective parameter jointly calibrated with the
drag** against the VMI final-velocity distribution (not the static 0.308 eV; not a
hand-picked "just-escapes" threshold; **do not weaken the drag to force escape** —
that detunes the validated quantity). Enforced like the mass pair: the effective
binding is stamped alongside the coefficients (`effective\_binding\_energy\_I\_ion\_eV`

* provenance) and the §6.5 guard refuses an un-jointly-validated drag↔binding
pairing. *Deferred principled alternative:* a **dynamical** (velocity-dependent)
barrier — the honest long-term fix, decided by the VMI distribution.

> > \*\*Provenance flag (2026-06-15, OQ1 — see mass-doc §10A).\*\* The static
> `0.308 eV` solvation barrier here is the \*\*$X\_2$/³Π deep-snowball\*\* value
> (\[IHe05] He–I⁺(³Π), $S\_{\\mathrm{I^+}}=-3578$ K in the \[I2-notes] He-DFT); the
> jointly-fit effective `0.1168 eV` is a factor \~2.6 shallower. \*\*If\*\* that
> lowering is purely \*dynamical\* (the in-window-correct drag delivers
> sub-barrier surface KE, as argued just above) — i.e. the trajectory was run on
> $X\_2$ and the ion simply never equilibrates to the deep well — then `0.1168 eV`
> is an $X\_2$-input number reduced by dynamics, and the mass-doc A10
> statistical-SO-mixture ladder (also shallower than $X\_2$) must \*\*not\*\*
> double-count the same reduction. This couples the drag-binding provenance to
> the mass-side electronic-picture choice; both are provisional pending author
> confirmation (mass-doc OQ1). User is in contact with the authors.

### 6.6 The $m\_\\text{eff}$ constant and optional time-resolved re-extraction

TDDFT shell (9 Å): \~21→19→14 He; extraction collapses this to the constant
$m\_\\text{eff}\\approx203$ amu (\~19 He), the window-representative value. Because
$F\_\\text{drag}=m,a-F\_C$ is a difference of comparable terms, the fractional error
from the constant-mass approximation tracks $|m(t)-m\_\\text{eff}|/m\_\\text{eff}$ —
**near-zero mid-window, \~⅓ only at the trajectory ends** (which coincide with the
§6.7 free-zone). So the in-hand constant-mass coefficients are usable as-is
mid-window; this is the **load-bearing defence for the R6 inconsistent pairing**
(MASS R6). **Option-3 (time-resolved) re-extraction** uses $m(t)$ directly in the
force balance $F\_\\text{drag}(t)=m(t)a(t)-F\_C(R(t))$ — an optional refinement, the
natural pairing for an evolving-mass scenario (notably the measured-loss $m(t)$),
not required for constant-mass scenarios. *Extraction-side action item:* verify the
stamped `extraction\_mass\_amu` is the mass the balance actually ran under, not a
relabelled value (a self-consistent refit cannot detect a wrong-but-consistent
mass; an earlier literal 179.912 vs the correct ≈202.954).

### 6.7 Transient free-extrapolation zone

No TDDFT data exists at the violent bubble-exit onset (the $t^\*$ cut, §2.3), so
neither drag nor mass law is anchored in the first \~0.5 ps under any scenario —
**widened to several ps on the mass side** by the early self-instability (MASS
A7/R9). Validation tolerances on $v(t)$, $R(t)$ here are **loose for all
scenarios**; mid/late observables remain the anchor. The pure-cubic high-$v$
over-braking in this same window is MASS **R10**.

### 6.8 Shell-loss evidence (informs Tier 1)

The TDDFT shell data (21→19→14, monotone loss) is the **Tier-1 target** and was the
original empirical lean toward stripping over accretion. The three-scenario
framing is superseded by the locked biphasic mechanism (§2.5 → MASS doc), but the
observation stands and is privileged for extraction: the directly-readable $m(t)$
(the measured loss) is the natural pairing for §6.6 option-3. The locked biphasic
mechanism subsumes both gain and loss; Tier 1/2 adjudicate where on the regime
axis the ion lands (MASS §6.11).

### 6.9 Histogram comparison metric (Tier 2/3)

**Primary — Wasserstein (earth-mover)** in the observable's physical units (Å/ps
for velocity, integer He-count for size): binning-free, exact on integer-$n$
support, directly interpretable. **Secondary — binned $\\chi^2$** (respects the VMI
bin convention, PHYSICS\_BASELINE §13; care with the low-count size-distribution
tail). **Secondary — KS** (cheap scalar screen; tail-insensitive, so not the
size-distribution adjudicator). `SimConfig.validation\_histogram\_metric ∈ {wasserstein, chi2, ks}`, primary `wasserstein`.

