# TIER2_MASS_SCENARIOS — what mass actually flies, and what that does to KE

**Status: OPEN INVESTIGATION, opened 2026-08-11. Nothing adopted, nothing
built. This document is the fresh-session entry point for the flight-mass
thread.** All results here are twin-side (`scripts/tier2_h2b_forward_model.py`)
and were produced by anchored scratchpad probes; no committed stage was
modified and no committed artifact was rewritten. Every number below is
reproducible from cached chord families with zero new MD.

> **Why this is a separate document.** It grew out of the RQ12 /
> helium-density-width thread (`TIER2_PARAMETER_INFLUENCE.md` §9.6) but it
> is not about the density width, not about E_bind, and not about the drag
> law. It is about the **mass of the object that flies through the Coulomb
> explosion** — a modelling commitment inherited from Tier-1a that has
> never been examined and that turns out to control the single largest
> unclosed residual in Tier 2, the KE₁ deficit.

> **Where the thread stands (end of 2026-08-11).** §1–§8 establish the
> kinematics: the deficit *is* the flight mass, arithmetically and
> measurably. §10–§11 dimension the obvious cure (born light, re-dressed
> in transit) and find it **half-works**: the kinematic half is fine, but
> pickup is linear in exposure and cannot generate the size distribution's
> width **at any λ₀**, while a stripped complex loses the E_int that fuels
> the evaporation cascade.
>
> **The governing result is §4.1, the shedding-cost law:**
> `KE₁ = 1.006 eV − ½·Δm·v_shed²`. The whole KE₁ question is **one number —
> the speed at which the helium departs** — and shedding cost scales as v²,
> so late shedding is expensive and the velocity peak (measured 13.44 Å/ps,
> identical for every bin) is the worst possible moment.
>
> Two routes were considered. **Route B (§13)** — velocity-dependent
> sticking — is **CLOSED (§13.4)**: a threshold sheds at v_c by
> construction, an n = 1 bin requires v_c > 9.584 Å/ps, and at that speed
> the cost already equals the status quo. No gain at any threshold, for any
> D₀. **Route A (§12)** — a graded partial strip of the outer shell,
> **KE₁ ≈ 0.83 eV** with the mechanism intact — survives, but now carries a
> sharper difficulty than its deep-tail price: the shell must leave at
> **v ≲ 3 Å/ps**, and no mechanism for that is identified. Open decisions:
> **M5 (§12.1)**, the deep tail, user + supervisor; **M1**, a sourced
> I⁺–He binding, which still gates Route A.
>
> **§14 records three withdrawn claims from this session — read it before
> reusing any number here.**

---

## 1. The governing identity

The twin integrates every trajectory at the **dressed** mass
`complex_mass_amu(ne_mol)`, and at the corrected geometry every fragment
carries `n₀ = 21` (measured 20000/20000; the shallowest birth in the
committed master is at depth −19.95 Å, so the pickup gate is saturated).
The scored observable uses the **bare** mass at detection:

```python
ke = kinetic_energy_eV(complex_mass_amu(n_det[mask]), v_inf[mask])
```

Hence, for any energy the complex gains or loses while dressed,

$$\frac{\partial KE_1}{\partial E_\text{dressed}} \;=\; \frac{m(1)}{m(n_\text{flight})}
\;=\; \frac{130.903}{210.955} \;=\; \mathbf{0.6205}$$

This factor is not specific to E_bind. It applies to the well toll, to
drag losses, to the source budget — to **everything pre-evaporation**.
Its discovery is recorded in D0 §9.6; the consequences are here.

**It is not an error.** It follows from the Tier-1a velocity-preserving
shed (an evaporating atom leaves at the complex velocity, carrying its own
½m_He v² away), combined with the measured fact — from the s(n) coupling
probe — that ions exit the droplet **fully dressed** (n̄ 19; no ion reaches
n ≤ 8 inside). Pay the toll dressed, shed outside. The bookkeeping is
self-consistent. The question this document opens is whether the *premise*
— that 21 helium atoms co-move through a 2.7 eV Coulomb explosion — is
physically defensible.

---

## 2. The §3.5g cross-check — the attenuation does NOT explain "exhausted"

**Result: NEGATIVE. Hypothesis refuted, recorded so it is not re-proposed.**

The first thing tried was whether the 0.6205 attenuation explains §3.5g's
verdict that the in-surface (v_c, τ, E₀) freedom is exhausted on the KE₁
axis. It does not. §3.5g's verdict is a **ratio** — in-gate bound
+0.186 eV against a required +0.309 eV — and the attenuation scales
numerator and denominator identically:

| frame | available | required | ratio |
|---|---|---|---|
| bare (scored) | 0.186 | 0.309 | **0.6019** |
| dressed | 0.300 | 0.498 | **0.6019** |

Frame-invariant. The in-surface freedom is exhausted for reasons that have
nothing to do with the mass frame, and §3.5g's conclusion stands
unmodified.

**What the frame conversion does buy** is a restatement of the deficit in
physical rather than observable units — see §3.

---

## 3. The energy ledger of the scored n = 1 ions

h405 (`capped_cubic` v_c 5.5 / τ 4.4 / E₀ 0.405), corrected geometry,
6375 scored ions, m = 20000:

| term | eV | note |
|---|---|---|
| source budget `KE_COUL / r0_sep / 2` | 2.7006 | per fragment, r0_sep 2.666 Å |
| − solvation toll actually paid | **0.1167** | 99.9 % of E_bind = 0.1168 |
| − drag dissipation | 1.5790 | **58 % of the budget** |
| = exit KE, dressed frame | 1.0050 | |
| × m(1)/m(21) | 0.6205 | |
| = scored KE₁ | 0.6236 | matches the committed value |

Deficit against the adjudicated 1.00 eV anchor: **0.3764 eV bare =
0.6066 eV dressed = 38.4 % of the drag dissipation = 22.5 % of the source
budget.** §3.5g's in-gate bound is 19.0 % of the drag dissipation. So the
frame-free statement of §3.5g is: **the gates permit recovering 19 % of
the drag loss; the anchor requires 38 %.** Factor 2.0 — the same 0.60
ratio, in units that mean something.

The cleanest statement of all, with no mass bookkeeping: **the model's
n = 1 ions exit at 9.584 Å/ps; the anchor wants 12.14 Å/ps. A 27 %
velocity shortfall.**

---

## 4. The flight-mass measurement — exit ENERGY is mass-invariant

**Analytic expectation.** Under h405's cap with `p_tail = −1`, the drag
force above v_c is `γ·v = b·v_c³` — **constant, velocity-independent**.
The well toll and the path through ρ are both geometric. So the work done
on an escaping fragment depends only on the path it takes, not on how fast
it walks it: the exit *energy* should not depend on the flight mass at all.
Only the velocity changes, and therefore only the scored KE.

**Measured** (h405 chord re-integrated at two flight masses; membership
frozen at the production leg; the m(21) leg verified **bit-identical** to
the committed chord npz before reading the m(1) leg):

| flight mass | E_exit [eV] | v_inf [Å/ps] | scored KE₁ [eV] | implied drag loss [eV] |
|---|---|---|---|---|
| m(21) = 210.955 (production) | 1.0050 | 9.584 | 0.6236 | 1.5789 |
| m(1) = 130.903 (bare) | **1.0059** | **12.172** | **1.0059** | 1.5780 |

**Exit energy invariant to 0.09 %; drag loss invariant to 0.06 %.**
Confirmed as analytic, not coincidental.

### The law

$$\boxed{\;KE_1 \;=\; E_\text{exit}\cdot\frac{m(1)}{m(n_\text{flight})}
\;=\; 1.006\ \text{eV}\times\frac{130.903}{m(n_\text{flight})}\;}$$

| n_flight | 1 | 2 | 3 | 4 | 8 | 12 | 19 | 21 |
|---|---|---|---|---|---|---|---|---|
| predicted KE₁ [eV] | 1.006 | 0.976 | 0.948 | 0.921 | 0.829 | 0.753 | **0.649** | **0.624** |

### 4.1 The shedding-cost law — the same identity, generalised (2026-08-12)

The law above assumes a mass held constant through the flight. The useful
generalisation lets the helium leave **anywhere** along the trajectory.
Because the drag work is mass-invariant (measured) and the Coulomb budget
is fixed, the only thing that changes with *when* helium departs is the
kinetic energy it carries off:

$$\boxed{\;KE_1 \;=\; E_\text{exit} \;-\; \tfrac12\,\Delta m\,v_\text{shed}^2
\;=\; 1.006\ \text{eV} \;-\; \tfrac12\,\Delta m\,v_\text{shed}^2\;}$$

It is algebraically the *same* identity — substituting
`v = √(2E_exit/m_flight)` turns one into the other — but expressed in the
variable that actually matters physically. **The entire KE₁ question is one
number: the speed at which the helium departs.**

| when the Δm = 80.05 amu leaves | v_shed [Å/ps] | cost [eV] | KE₁ [eV] | status |
|---|---|---|---|---|
| v = 0 (never dressed / bare flight) | 0 | 0 | 1.006 | **measured 1.0059** |
| at exit — the current model | 9.584 | 0.381 | 0.625 | **measured 0.6236** |
| at the velocity **peak** | 13.440 | 0.749 | 0.257 | law-predicted, not measured |

Validated at both endpoints; the third row is an extrapolation *within* the
validated form, not an independent measurement.

**Consequences, and they order the whole landscape:**

- **Shedding cost scales as v², so late shedding is expensive and the
  velocity peak is the worst possible moment.** Jettisoning mass when the
  ion is fast is not a saving — the departing atoms take ½Δm v² with them.
- **The optimum is to shed as early as possible**, i.e. at v ≈ 0 — which
  means never having been dressed. KE₁ = 1.00 requires
  v_shed ≲ **1.2 Å/ps**.
- **It retro-explains (A) exit stripping's measured 0.708 ceiling**:
  stripping at the exit sheds at essentially the same speed the model
  already sheds at, so it buys nothing.
- **It gives Route A its real criterion** (§12): the outer shell must go
  **early and slow**, not at the surface and not at the peak.

### Both existing measurements fall on this curve

- The **twin** flies at n = 21 → predicted 0.624, measured **0.6236**.
- **MD** sheds ~2 He inside the droplet and exits at n̄ ≈ 19 (s(n) probe)
  → predicted 0.649, measured **0.6406** (committed
  `atlas_ke_lown_scan.csv`, h405 pooled). 1.2 %.

The law therefore reproduces the twin↔MD KE₁ gap, which had never been
explained.

### What this means

The 0.36 eV KE₁ deficit — which survived the §3.5g retro-scan, the §3.5h
p_tail ring (kill fired), the s(n) coupling (all predictions refuted) and
the entire (C) channel-mixture programme (registration failed) — is
**arithmetically identical to one assumption**: that the ion drags 21
helium atoms through the Coulomb explosion. The required flight mass is
not knife-edge: any **n_flight ≈ 1–3** lands at 0.95–1.006 eV. What is
excluded is n_flight ≈ 21.

**Anti-circularity note.** The drag surface (v_c, τ, E₀, b) was arbitrated
against the **size distribution** — n̄, n₁, W₁ — and never against KE₁,
which has been the open residual throughout. That arbitration
independently leaves an exit energy of 1.006 eV, matching the KE anchor to
0.6 % once the mass frame is corrected. This is not a fit. It is,
however, a *single* coincidence and should be treated as suggestive until
the shape problem in §5 is resolved.

---

## 5. THE SHAPE PROBLEM — a uniform flight-mass change is the p_tail failure mode

**This is the decisive constraint and it kills the naive version.** The
identity applies to *every* bin, so lowering the flight mass uniformly
boosts every bin by the same factor — precisely what killed the §3.5h
p_tail lever (uniform tail softening → midHot 1.9–4.4).

Measured KE-vs-n at h405 under three flight-mass conventions (membership
frozen at the production leg; reference = `g3_ref_mean_ke()`):

| n | count | ⟨v_inf⟩ [Å/ps] | ref KE | **A** fly@21 | **B** fly@m(n) | **C** fly@1 |
|---|---|---|---|---|---|---|
| 0 | 7254 | 10.665 | 3.706 | 0.749 | 1.246 | 1.207 |
| 1 | 6375 | 9.584 | 1.302 | 0.624 | 1.005 | 1.005 |
| 2 | 4304 | 8.792 | 0.706 | 0.541 | 0.846 | 0.871 |
| 3 | 3369 | 8.116 | 0.496 | 0.474 | 0.721 | 0.765 |
| 4 | 2334 | 7.528 | 0.390 | 0.420 | 0.620 | 0.677 |
| 6 | 1901 | 6.383 | 0.267 | 0.319 | 0.446 | 0.514 |
| 8 | 1322 | 5.181 | 0.192 | 0.221 | 0.294 | 0.357 |
| 10 | 1006 | 4.218 | 0.142 | 0.154 | 0.195 | 0.249 |
| 12 | 832 | 3.466 | 0.107 | 0.109 | 0.132 | 0.176 |
| 14 | 644 | 2.733 | 0.081 | 0.071 | 0.082 | 0.115 |
| 17 | 487 | 1.587 | 0.066 | 0.027 | 0.029 | 0.043 |

**A** = production (fly dressed at 21, score at m(n)).
**B** = shed early to the final composition, then fly.
**C** = fully bare flight for everyone (the uniform strawman).

| convention | KE₁ | midHot (n 2–8) | deepKE (n 10–17) |
|---|---|---|---|
| **A** fly@21 (production) | 0.624 | **1.068** | 0.844 |
| **B** fly@m(n) | **1.005** | 1.529 | **0.999** |
| **C** fly@1 uniform | **1.005** | **1.722** | 1.360 |

C is a disaster (midHot 1.72 — the p_tail signature almost exactly).
B is the physically-motivated version and still lands midHot 1.53.
Notably **B nails both ends** (KE₁ 1.005, deepKE 0.999) and wrecks only
the middle — but wrecking the middle is disqualifying.

### The required correction is a U, not a boost

`needed / model A` per bin (n = 1 against the adjudicated 1.00 anchor,
all others against the reference mean):

| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 13 | 15 | 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| required factor | 4.95 | **1.60** | **1.31** | 1.05 | 0.93 | 0.86 | **0.84** | 0.87 | 0.92 | 0.98 | 1.06 | 1.31 | 2.45 |

Boost is required **only at n ≤ 3 and n ≥ 13**. The middle needs
*cooling* by 8–16 %. Inverting into flight masses (`m_flight = m(21)/b`):

| n | required b | implied m_flight | interpretation |
|---|---|---|---|
| 1 | 1.60 | 132 ≈ m(1) | strip essentially everything |
| 2 | 1.31 | 161 ≈ m(8) | strip ~13 He |
| 3 | 1.05 | 201 ≈ m(18) | strip ~3 He |
| ≥ 4 | < 1 | > m(21) | **impossible** — needs a different mechanism |

So the required profile is a near-**threshold**: strip hard above some
speed, nothing below, with the threshold between the n = 2 mean
(8.79 Å/ps) and the n = 3 mean (8.12 Å/ps).

---

## 6. The co-moving argument — an independent estimate of the same threshold

To co-move with the ion at speed v, a helium atom needs kinetic energy
½·m_He·v² = 2.074×10⁻⁴·v² eV (v in Å/ps):

| v [Å/ps] | 12.17 (anchor) | 9.58 (n = 1) | 8.79 (n = 2) | 8.5 | 8.12 (n = 3) |
|---|---|---|---|---|---|
| co-moving KE [meV] | 30.7 | 19.0 | 16.0 | **15.0** | 13.7 |

Against the charge-induced-dipole binding of He to I⁺,
`U = −α q²/(2r⁴) = −1.476/r⁴ eV` with α(He) = 0.2050 Å³:

| r [Å] | 2.5 | 3.0 | 3.5 | 5.0 |
|---|---|---|---|---|
| binding [meV] | 37.8 | **18.2** | 9.8 | 2.4 |

The first shell is **marginally** able to follow; everything outside it is
bound an order of magnitude too weakly. The crossing sits at ~15–18 meV,
i.e. at **v ≈ 8–9 Å/ps** — the same threshold the KE data demands, from
physics that knows nothing about the KE spectrum.

**Caveat, explicit: these polarizability numbers are a back-of-envelope
estimate, not sourced, and a snowball's binding is not a single pair
potential.** Sourcing a real I⁺–He interaction and a snowball binding
profile is a prerequisite for any build.

> **And the model disagrees with the envelope.** The committed ladder
> (`d0_of_n`) gives D₀ ≈ **9.2 meV** for the core rungs and 5.5–6.0 meV at
> n = 21 — roughly **half** the 18 meV used above, putting the threshold at
> **6.7 Å/ps** rather than 9.3. That is below the n = 3–5 population and
> would strip into the mid band. The full four-way comparison and its
> consequences are §13.3 (R1/R2); the numbers in this section should be
> read against that table, not on their own.

### Why a speed gate is not the p_tail lever

p_tail softened the drag at *all* speeds, so it lifted every bin. A
co-moving-mass threshold acts only above ~8.5 Å/ps, and n ≥ 3 never
reaches that (8.12 Å/ps and falling monotonically to 1.59 at n = 17). The
mid and deep bands are structurally untouched. **The selectivity comes
from gating on velocity, and velocity is tightly tied to n in this model.**

---

## 7. Why the shed convention gives no boost, and what the recoil channel is worth

Under the velocity-preserving shed, shedding produces **zero** kinetic
boost by construction: the atom leaves at the complex velocity, v is
unchanged, and `KE = ½m(n)v²` merely gets *smaller* for fragments that
shed more. The current implementation therefore **penalises low n**, which
is the opposite of what the reference demands. The flight-mass finding is
the same issue viewed from the other side.

The alternative — evaporation **recoil** — is real but small. In a
two-body separation releasing ε, the complex retains
`m_He/(m_He + m_c) ≈ 2 %` of it; directions are random, so the kicks
cancel in momentum but **add coherently in energy**. Over a full cascade
this is an estimated **~20–40 meV**, biased toward low n (they shed most).
A genuine differential low-n lever worth ~10 % of the 376 meV deficit —
not the answer, but it belongs in the ledger, and it is also a candidate
contributor to the Tier-3 second-moment under-dispersion.

*(Σε over the ladder is estimated, not measured. Cheap to measure.)*

---

## 8. Three discrepancies, not one (and a withdrawn fourth)

| # | band | status | addressed by flight mass? |
|---|---|---|---|
| 1 | n ≤ 2 too cold (1.60×, 1.31×) | explained exactly | **yes**, if speed-gated |
| 2 | n = 4–12 too hot by 8–16 % | untouched | no |
| 3 | n ≥ 13 too cold, up to 2.45× | untouched | no |
| — | n = 0 needs 4.95× | known slow-bare problem (§3.5l) | no |

> **CORRECTION (same session).** An earlier draft of this section read:
> *"the model's KE spans a factor 28, the reference 56 — the dressing↔velocity
> correlation is a factor two too weak."* **That was an artifact of including
> n = 0.** Over the bins that enter the histogram:
>
> | span n = 1 → 17 | model | reference |
> |---|---|---|
> | KE ratio | 0.624 / 0.027 = **23.1** | 1.302 / 0.066 = **19.7** (15.2 vs the 1.00 anchor) |
>
> The model's dynamic range is **not too small — if anything slightly too
> large.** The "factor two too weak" claim is withdrawn. It was the single
> statement making simultaneous landing look structurally hopeless, and it
> was wrong.

**n = 0 is separately and definitively unreachable, and it is not a mass
problem.** The reference bare-ion KE is **3.706 eV** against a per-fragment
source budget of **2.7006 eV**. No mass scenario, drag law or cascade can
give a fragment more kinetic energy than the source deposits. The bare peak
requires a higher-budget channel — the paper's Q3 at 4.32 eV per I⁺ — so it
is a **source-channel observable** (consistent with §3.5l "bare KED is
fast-fed") and belongs on its own axis, excluded from the
histogram-plus-KE question.

What is left is therefore not a curve needing global re-shaping. It is a
boost that must be **confined to n ≤ 2**, with the middle left alone. That
is a sharp, testable requirement rather than a structural impossibility —
and it makes the *sharpness of the velocity gate* the whole question
(§6, and M6).

---

## 9. The central tension

**The histogram wants the helium present during flight; the kinematics
wants it gone.**

The candidate that could serve both is **born light + re-dressed during
transit + evaporation**, rather than the current birth-fully-dressed +
evaporation-only. §10 dimensions it. §11 records what it does to the
mass mechanism. The short version: the kinematic half works, the
width-generating half does not, and §12 is the variant that survives both.

---

## 10. The re-dressing question — dimensioned (2026-08-11)

### 10.1 "Instantaneous stripping" is the wrong framing

At t = 0 every ion is identical: same Coulomb impulse, same
`r0_sep = 2.666 Å`, born deep where ρ̂ ≈ 1. Nothing distinguishes them. So
**any rule applied at t = 0 strips the whole ensemble equally**, which is
convention C — midHot 1.72, the §3.5h p_tail failure mode. A selective
strip cannot be imposed at birth.

### 10.2 The three-phase picture

The differentiation is not in the stripping. It is in the **re-capture**.

1. **Everyone accelerates, everyone strips (~first 1.5 ps).** Coulomb
   dominates drag by >10× while the fragments are close (at r_sep = 5 Å:
   ≈ 5600 vs ≈ 420 amu·Å/ps²). Every ion is driven to ~10–15 Å/ps, where
   no He can co-move. Uniform, by construction.

   > **MEASURED 2026-08-12.** `v_peak = 13.44 Å/ps` for **every detected
   > bin, identical to four digits** (n = 0 through n = 17), reached by
   > **t ≈ 0.75 ps** (11.81 at 0.25 ps, 13.28 at 0.50, plateaued by 0.75).
   > Phase 1 is confirmed: early in the flight the ensemble really is one
   > object, and all of the eventual spread in v_inf (10.665 → 1.587)
   > develops *afterwards*, in Phase 2. Headroom above the exit speed runs
   > from +26 % (n = 0) to +747 % (n = 17).
2. **Coulomb dies off past r_sep ≈ 20 Å; drag takes over.** The
   acceleration consumed only the first ~9 Å of a ~34 Å median birth
   depth, so ~25 Å of droplet remains. *Here* the ions diverge, set by
   birth position and direction.
3. **The slow re-dress; the fast do not.** Short remaining path ⇒ exits
   fast ⇒ nothing sticks ⇒ **n = 0–2, high KE**. Long remaining path ⇒
   decelerated while still buried ⇒ He sticks again ⇒ **n = 10–17, low KE**.

So it is a **universal boost with a selective clawback** — the clawback
being momentum-sharing on capture (`v → v·m/(m+m_He)`, the code's
`he_capture_velocity="at_rest"` reset) plus the longer drag path. Both hit
the middle and deep bins; neither touches n = 1, which by definition is an
ion that never re-dressed. No mixture is imposed; birth geometry sorts the
ensemble.

### 10.3 MEASURED: pickup cannot generate the size distribution, at any λ₀

Expected re-dressing computed from the model's own committed trajectories
(`n_pickup = λ₀·∫ρ̂ dt`, with `∫ρ̂ dt = K_raw·τ_ref` per fragment; h405,
m = 20000, landmark oracle green):

| n | count | ⟨t_exit⟩ [ps] | ⟨∫ρ̂ dt⟩ [ps] | pickup at λ₀ = 0.9 | Langmuir-capped |
|---|---|---|---|---|---|
| 0 | 7254 | 2.62 | 2.65 | 2.4 | 2.3 |
| 1 | 6375 | 3.15 | 3.19 | 2.9 | 2.7 |
| 2 | 4304 | 3.54 | 3.59 | 3.2 | 3.0 |
| 4 | 2334 | 4.13 | 4.21 | 3.8 | 3.5 |
| 8 | 1322 | 5.22 | 5.40 | 4.9 | 4.3 |
| 12 | 832 | 6.68 | 7.04 | 6.3 | 5.5 |
| 17 | 487 | 9.81 | 10.68 | 9.6 | 7.7 |

Exposure rises with n as the picture requires — but only from 2.65 to
10.68 ps, **a factor 4.0**. The deep ions reach 7.7 atoms, not 17.

> **The obstruction is λ₀-independent.** Pickup is **linear** in exposure
> (`n ≈ λ₀·∫ρ̂ dt`; the Langmuir cap only compresses further), so the
> spread in n it can generate *equals the spread in exposure* — measured
> 4.0×. The size distribution spans ~17×. Raising λ₀ multiplies every bin
> alike: at λ₀ = 3.3/ps the deep ions reach 17 but the fast ions reach 7,
> not 0. **Pickup cannot generate the width at any rate constant.**

**Why the existing architecture can.** Evaporation runs off
`E_ej = E₀·e^(−K)` — **exponential** in the same exposure. A factor-4
spread in K becomes an enormous spread in ejection energy and hence in
final n. That is how the model produces the width, and presumably why it
was designed as *dress fully, then evaporate down* rather than *start bare,
pick up*.

**Consequence: the ion must be dressed at some point, or the exponential
has nothing to act on.**

### 10.4 What remains open

The two requirements sit at **different times** — light during the ~1.5 ps
Coulomb impulse, dressed during the remaining 3–10 ps of transit. Those do
not have to conflict. Subtracting the impulse window from the exposures
above leaves the deep ions ~9 ps, needing λ₀ ≈ 3–4× the pin to reach
n ≈ 21 — at a cost of 12–16 amu/ps of extra friction (§11.4), ~30 % on top
of the drag law. A three-way trade: KE gain from the light window, width
from the dressed window, drag cost from the higher λ₀. **Not settleable by
hand** (§14).

---

## 11. Does the mass mechanism survive a strip-and-rebuild? (2026-08-11)

Checked against the actual bookkeeping —
`physics/internal_energy_budget.py`, `physics/pickup.py`, `fate_map`:

```
S2 onset:   E_int(0) = f_int · E_avail       (0.5 × 2.70 ≈ 1.35 eV, once, at ionization)
S1 pickup:  dE_int   = +f_ret · D_0(n+1)     (f_ret = 0.1)
K1 shed:    dE_int   = −D_0(n)               (energy-GATED: needs E_int ≥ D_0)
```

### 11.1 Survives: the ladder

`D₀(n)` / `Σ(n)` is a **state function of occupancy**, history-independent.
A complex at n = 8 is the same whether it evaporated down from 21 or built
up from 0. The ladder picture is untouched.

### 11.2 Survives as physics: post-ejection evaporation

Once ρ̂ → 0, pickup and drag stop and a hot cluster relaxes in vacuum by
RRK over µs. Unchanged. What changes is only the **initial condition**
handed to it — and per §11.3 that is much colder and smaller, which would
make the µs stage nearly inert. The detected distribution would then *be*
the pickup distribution, i.e. §10.3's 4× spread. The obstruction arrives
twice, by independent routes.

### 11.3 Does NOT survive: the fuel

- **S2's deposit leaves with the shell.** `E_int(0) = f_int·E_avail` is
  deposited into the *complex*; E_int is the vibrational/rotational
  reservoir of I⁺Heₙ, and a bare I⁺ has only electronic states. Strip the
  shell at ~1 ps and that energy departs with the atoms. The cascade's fuel
  tank is emptied. **This is exactly RQ1** ("provenance and magnitude of
  E_int(0)", working hypothesis 0.2–0.5 eV) — a strip scenario does not
  perturb RQ1, it makes the model convention untenable.
- **The rebuild cannot refill it.** Condensation does release latent heat
  and S1 books it, but at the pinned **f_ret = 0.1 a capture banks a tenth
  of what re-evaporating that atom costs** (`+f_ret·D₀` in, `−D₀` out). A
  rebuilt complex is **sub-critically heated and freezes at what it picked
  up.** No cascade, no width.

### 11.4 Would require a channel the model does not have

K1 is *energy-gated* — thermal evaporation. Mechanical stripping is
**athermal**: a cold, fast ion strips regardless of E_int, and the bond
energy comes out of **kinetic** energy. That is a new channel, ungated by
E_int, drawing on E_kin, and it would force a re-derivation of the 5-term
invariant.

**And it has a measured price.** The (C) programme already built one
(`physics/exit_strip.py`). The probe measured it **over-tolls**: ε × knock
counts at p90 ≈ 20 gave 0.5–0.7 eV against the 0.2 eV anchor, and it *fed*
the suppressed gate (supp 0.18 → 0.55 in A-only). Stripping ~21 atoms costs
roughly the KE the strip was meant to buy. The mitigating difference: (C)
stripped at the **exit**, paying the toll with no mass benefit, whereas
stripping during acceleration would at least buy the lighter flight mass.

Also relevant: the pickup channel's own capture reset contributes a
`ρ̂`-gated **linear** friction `γ_pickup = λ₀·ρ̂·m_He` = **3.6 amu/ps** at
the pinned λ₀ — 5–10 % of the explicit drag law, apparently never
accounted anywhere. At λ₀ = 3–4× it becomes 12–16 amu/ps.

---

## 12. Route A — the graded partial strip

The ladder is **graded**: inner atoms bound far more tightly than outer
ones. So the physically natural strip is not total — it is the
loosely-bound outer shell:

- outer atoms strip cheaply (low D₀, small toll, little E_int carried off);
- the tightly-bound **core survives**, holding most of Σ and most of E_int
  — **the cascade keeps its fuel**;
- flight mass drops 21 → ~8, i.e. `m(21)/m(8) = 1.33`, giving
  **KE₁ ≈ 0.624 × 1.33 ≈ 0.83 eV**.

This keeps S2, keeps the ladder, keeps the µs-flight physics, and adds only
the athermal outer-shell channel. It lands in the same **0.83–0.86 eV**
window as the two unrelated estimates in §5 and §10.4.

> **The criterion this route must satisfy, from §4.1 (2026-08-12).** The
> gain is `−½·Δm·v_shed²`, so the outer shell must leave **early and slow**
> — not at the surface (that is (A), measured ceiling 0.708) and not at the
> velocity peak (predicted 0.257, worse than doing nothing). Getting the
> full ~0.83 needs the ~13 outer atoms gone at **v ≲ 3 Å/ps**, i.e. inside
> the first few tenths of a picosecond, while the ion has barely started
> moving. **No physical mechanism for that has been identified** — at
> v ≈ 0 the helium has no reason to leave. This is now Route A's central
> difficulty, and it is sharper than the deep-tail question.

**But it caps n_det at ~8** — nothing can evaporate down from a core it
never had, so n = 10–17 disappears.

### 12.1 The pivot: is the deep tail load-bearing?

Every version of this thread converges on one decision:

| | deep tail preserved | deep tail released |
|---|---|---|
| ion must fly | heavy (n ≈ 21) | light-ish (core ≈ 8) |
| KE₁ | ~0.62 (the standing deficit) | **~0.83** |
| n = 10–17 | reproduced | lost |
| mechanism | untouched | S2 + ladder intact, athermal channel added |

Evidence bearing on it, for the supervisor conversation:

- the n = 17 bin holds **487 of 20000** scored ions;
- **deepKE is already an acknowledged 16 % cold** (0.844);
- §3.5j identified the n = 20–21 spike as the **`suppressed` fate class** —
  scoring scaffolding, not measured structure;
- the experimental n = 19–20 population is **0.59 %** against the model's
  17.6 % (≈ 30×, §3.5i.2).

**This is an experimental judgement, not a modelling one, and it is
explicitly the user's + supervisor's call (posed 2026-08-11).** It is
recorded here because it gates M5/M6 and the build decision.

---

## 13. Route B — velocity-dependent sticking (the `dwell_time` arm)

**Status: POSSIBLE ROUTE, NOT RECOMMENDED FOR BUILD. Recorded because it
is the only proposal in this thread that attacks the §10.3 width
obstruction rather than working around it — and because, on the model's
own numbers, it currently fails. Read §13.3 before §13.2.**

Proposed by the user (2026-08-11): make the pickup rate fall with speed,
so slow ions dress far more than fast ones.

### 13.1 What it is, and what it is not

$$\lambda_\text{attach} = \lambda_0\,\hat\rho\,\Big(1-\tfrac{n}{n^*}\Big)_+^{\,p}\cdot S(v)$$

with `S(v)` a **sticking probability** → 1 at rest, → 0 above a threshold
set by `½·m_He·v² = D₀`.

- **This is the `dwell_time` enum arm** (`PickupRateForm ∈ {density_only,
  sweeping, dwell_time}`) — already declared, rule-2 unbuilt.
- **It is NOT the gas-phase law.** Langevin ion–neutral capture gives
  σ ∝ v⁻¹, hence a *constant* rate coefficient `k = σv` — which is exactly
  `density_only`. (Generally σ ∝ v^(−4/s) for an r^(−s) potential; s = 4
  for ion-induced-dipole.) Getting λ ∝ 1/v requires σ ∝ v⁻², which needs
  s = 2, a Coulomb potential — no standard basis. **So this route is an
  accommodation/sticking claim, a different physical assertion from
  capture-cross-section scaling, and it must be argued on its own terms.**
- **Formulate it as sticking, not as ρ/v.** MASS §5 rules out
  *unregularized* dwell-time because ρ_He/v **diverges at rest**, against
  [Nat23]'s finite 2.0/ps resting-ion datum. A sticking factor → 1 at
  v → 0 satisfies that datum automatically. The 1/v spelling does not.

### 13.2 Why it is attractive

- **It breaks the λ₀-independence of the §10.3 obstruction.** With
  `n_pickup = λ₀∫ρ̂ S(v) dt`, fast ions are multiplied by ≈ 0 and slow ions
  by ≈ 1, so raising λ₀ raises **only the slow bins**. The proportionality
  the obstruction rested on is gone. Nothing else raised in this thread
  does that.
- **One function, three effects.** A single statement — *helium cannot
  stay attached above v_c* — produces the universal strip (§10.2 phase 1),
  the differential re-dressing (phase 3), and the n ≤ 2 KE boost.
- **The threshold is not a free parameter** once D₀ is known:
  `v_c = √(2D₀/m_He)`.
- It is a **rate-form change to an existing enum**, not a new mechanism.

### 13.3 RISKS — and the first one is currently disqualifying

**R1 — On the model's own ladder, the threshold is in the wrong place and
this route destroys the low-n bins.** Reading `d0_of_n`:

| n | D₀ [meV] | threshold √(2D₀/m_He) [Å/ps] |
|---|---|---|
| 1–12 (core) | 9.2 | **6.67** |
| 18 | 8.0–9.0 | 6.2–6.6 |
| 21 (outermost) | 5.5–6.0 | **5.1–5.4** |

Bin velocities are n = 1 → 9.58, n = 2 → 8.79, n = 3 → 8.12, n = 4 → 7.53,
n = 6 → 6.38, n = 8 → 5.18. A threshold at 6.67 Å/ps strips **everything
above n ≈ 5** — straight into the mid band — emptying n = 1–5 into n = 0.
That is the p_tail failure mode, not a fix. **Evaluated with committed
model quantities, Route B fails.**

**R2 — The scenario is a knife-edge in an uncertain parameter.** The
outcome is a strong function of D₀, and the candidate values span a factor
four:

| source of D₀ | value | threshold | consequence |
|---|---|---|---|
| model's ladder (`d0_of_n`) | 9.2 meV | 6.7 Å/ps | strips into the middle — **low-n bins destroyed** |
| **what the KE data needs** | **≈ 15 meV** | **8.5 Å/ps** | strips only n ≤ 2 — **works** |
| charge-induced-dipole envelope, 3 Å | 18 meV | 9.3 Å/ps | strips n ≤ 1 — marginal |
| CID at 2.5 Å | 38 meV | 13.6 Å/ps | **nothing strips — route dead** |

Three of the four outcomes are failures, and the model's own value is one
of them. This is not a robust mechanism; it is a narrow window in a
quantity nobody has sourced. **M1 is therefore not a gate on this route —
it is the route's entire content.**

**R3 — The friction lands where it is least wanted.** `γ_pickup =
λ₀·ρ̂·S(v)·m_He` is concentrated where `S ≈ 1`, i.e. on the **slow** ions —
exactly the deep bins, which need *heating* (required factor 1.06–2.45),
not extra drag. At λ₀ ≈ 2/ps that is ~8 amu/ps applied preferentially to
the population already 16 % cold.

**R4 — The E_int fuel problem (§11.3) is untouched.** A stripped complex
still loses S2's `f_int·E_avail`, and rebuild heat at `f_ret = 0.1` is
still an order of magnitude short of re-evaporation. Route B changes *when*
helium attaches, not *where the cascade's energy comes from*.

**R5 — The middle and n = 0 are untouched.** Still 8–16 % hot and still
energetically unreachable respectively.

**R6 — It requires overturning a doubly-argued locked choice.** MASS §5
supports `density_only` with **two independent arguments**: the [Nat23]
resting-ion datum, and (as established in this discussion) the gas-phase
Langevin cancellation. A rate-form change has to defeat both, and MASS
§1920 admits a v-dependence *"only if density-only fails Tier 1/2"* — a
condition that has not been demonstrated.

**R7 — The surface-retention constraint, unsatisfied by any version.** The
experimental n = 1 peak at 1.00 eV is v = 12.14 Å/ps, at which a co-moving
He carries **30.6 meV** against the model's 9.2 meV rung. There is a
legitimate resolution — the detachment criterion applies only *inside* the
droplet, where relative motion exists; in vacuum a co-moving atom is
stable — but it forces the surviving atom to be acquired or retained
**right at the surface, at nearly full speed**. No variant in this document
yet satisfies that.

**R8 — All sizing here is envelope arithmetic**, and §14 records this
thread going 0-for-3 on envelope estimates. Treat every number in §13.2
as indicative only.

### 13.4 R9 — the structural kill (2026-08-12): Route B cannot gain KE at ANY threshold

R1 said Route B fails at the *current* ladder value. The shedding-cost law
(§4.1) makes the failure independent of D₀ altogether.

A sticking threshold sheds at `v = v_c` **by construction**, and every ion
crosses it — measured, since all reach the same v_peak = 13.44 Å/ps. So
every ion pays `½·Δm·v_c²`. Now impose the two constraints together:

- **For an n = 1 bin to exist at all**, the threshold must sit *above* the
  n = 1 exit speed (9.584 Å/ps) — otherwise those ions cannot hold their
  atom on the way out. So **v_c > 9.584**.
- **At v_c = 9.584 the strip costs 0.381 eV — exactly what the current
  model's exit-shed already costs.** Every higher threshold costs more.

$$v_c < 9.584 \Rightarrow \text{no } n=1 \text{ bin};\qquad
v_c \ge 9.584 \Rightarrow KE_1 \le 0.625\ \text{eV}$$

**There is no corner. Route B yields no KE₁ gain at any threshold, for any
value of D₀.** The user's specific variant — put the threshold just below
the *peak* so everyone strips while fastest — is the worst case of all: it
sheds at 13.44 Å/ps for a predicted KE₁ ≈ **0.257 eV**, *worse than doing
nothing*. The intuition is exactly inverted; fast is when jettisoning mass
costs the most.

### 13.5 Verdict

Route B was the **only attack on the width obstruction** in this thread,
and it is now closed on kinematics rather than on binding energies — a
cleaner kill than R1, because it does not depend on M1's outcome. It
should not be built.

M1 remains the decisive measurement for **Route A**, whose viability still
turns on whether an outer shell is loosely enough bound to leave early. But
it no longer decides Route B, which fails regardless.

---

## 14. Corrections issued during this thread

Recorded because two were reported before being caught, and the pattern
matters more than either error.

1. **Langevin capture rate — WITHDRAWN.** An imported ion–neutral rate
   (k_L ≈ 538 Å³/ps ⇒ λ ≈ 11.7/ps) was used to argue that a stripped ion
   re-dresses inside the Coulomb impulse, and that accretion derives the
   `pure_linear` drag law with γ ≈ 47 amu/ps. It is **13× above the
   project's pinned λ₀ = 0.9/ps** and 6× above [Nat23]'s loose 2.0/ps
   upper bound. It also violates **assumption A2** (MASS doc), which
   explicitly forbids importing a pickup *rate*: *"no alkali rate can be
   lifted directly — even the ordering argument forbids it … calibrate λ₀
   entirely against this work's I⁺ TDDFT shell build-up and the size
   distribution."* Everything hung on it is withdrawn — including the claim
   that the linear drag family gains a physical derivation. **Nothing here
   supports the §6.3/§6.4 adoption gate.**
2. **"Slow ions have far longer exposure, so pickup delivers plenty" —
   WITHDRAWN.** Right in direction, wrong in magnitude: the exposure spread
   is only 4.0× and the deep bins reach 7.7 atoms, not 17 (§10.3).
3. **"Dynamic range a factor two too weak" — WITHDRAWN** (§8): an artifact
   of including the n = 0 bin.

What survives all three: the governing identity (§1), the mass-invariance
measurement (§4), the KE-vs-n table (§5), and the §3.5g negative (§2) —
none of which came from an imported number.

**Standing lesson:** every one of these came from reasoning by plausible
scaling instead of reading the model's own numbers. Envelope calculations
in this thread have a poor track record; prefer a read against a committed
artifact.

---

## 15. Open questions, in the order they gate

| # | question | cost | gate |
|---|---|---|---|
| **M1** | **Source D₀ for I⁺–He.** Upgraded 2026-08-11 from a gating question to **the decisive measurement**: §13.3 shows the outcome is a strong function of D₀ across a factor-four candidate range (model ladder 9.2 meV / KE-data requirement ≈ 15 / CID envelope 18–38), and **three of the four outcomes are failures**. Also answers whether a rigid first shell can co-move at 10 Å/ps at all. Needs a sourced I⁺–He / snowball binding, not §6's envelope. | literature | **user / physics discussion — precedes everything** |
| **M2** | **The K re-score.** A faster flight crosses the droplet in less time ⇒ `K = ∫ρ dt` falls ⇒ `E_ej = E₀e^(−K)` rises ⇒ more evaporation ⇒ n̄ moves. The size distribution is currently matched at n̄ ≈ 4.5. Parked from the 2026-08-11 session. | zero MD | trigger |
| **M3** | Does exit-energy mass-invariance hold on the **linear** arm? It follows from the *cap's* constant force; under `γ = ρ·a` a faster ion loses more, so invariance should be partial. Bears on the §6.3/§6.4 adoption gate. | zero MD | trigger |
| **M4** | Measure Σε over the ladder and turn §7's 20–40 meV recoil estimate into a number. | zero MD | trigger |
| **M5** | **Is the deep tail load-bearing?** (§12.1) The pivot every version converges on: preserve n = 10–17 and the ion must fly heavy (KE₁ ~0.62); release it and the graded partial strip gives ~0.83 with the mechanism intact. | judgement | **user + supervisor — posed 2026-08-11** |
| **M6** | A speed-gated strip is pre-registrable (threshold ≈ 8.5 Å/ps, or a binding cut ≈ 15–18 meV) with hard predictions on all three bands. **It must leave midHot within a few percent of 1.068 or it is dead on arrival, exactly as p_tail was.** | design + MD | **not before M1 + M5** |
| **M7** | **Where does E_int(0) live if the shell is stripped?** (§11.3) RQ1 becomes structural, not just uncertain: S2's `f_int·E_avail` departs with the shell, and rebuild heat at `f_ret = 0.1` is 10× short of re-evaporation. Any strip variant needs a new answer before it can be dimensioned. | design | with M1 |
| **M8** | **Account the pickup channel's own drag.** `γ_pickup = λ₀·ρ̂·m_He` = 3.6 amu/ps at the pin (§11.4) — 5–10 % of the drag law, apparently never booked. Also: does it double-count with Tier-0's `b`, which was fitted at `mass_scenario=fixed`? (Not asserted — Tier-1a says the 9 Å traces *shed* in-window, which cuts the other way.) | zero MD | trigger |

---

## 16. Standing rules for this thread

- **Nothing here is adopted, and no build is authorised.** New drag-program
  code stays behind `[PROCEED TO IMPLEMENTATION]` (CLAUDE.md).
- **Do not import a pickup rate.** MASS assumption **A2** forbids it: only
  the *mechanism* (Poisson, ∝ ρ̂) transfers, never the rate. λ₀ is
  calibrated against this work's I⁺ TDDFT shell build-up and the size
  distribution. Violated once in this thread (§13.1) — do not repeat.
- **The twin cannot test any strip/re-dress scenario.** It fixes `ne_mol`
  at birth and integrates at constant mass, with no pickup channel at all.
  Every twin number in this document therefore omits both the mass growth
  and the 3.6 amu/ps of capture friction — a twin↔MD asymmetry not in the
  authority box.
- **The shed convention is now load-bearing on a primary observable.** The
  0.6205 follows from Tier-1a's velocity-preserving shed, chosen on
  Tier-1a grounds with cold-shed retained as a diagnostic bound. It should
  be re-affirmed as an explicit model commitment rather than carried as an
  inherited default.
- **One E_bind is charged for the whole I⁺He₂₁ complex.** The 20 bound He
  leave the droplet too and are charged nothing (~0.4 meV/atom ⇒ ~8 meV,
  ~7 % of E_bind). Small, real, and previously unnamed.
- **Any probe anchors on a committed artifact before a new number is
  read** (the rule adopted after the two RQ12 probe defects). Every result
  in this document did so; the anchors are named in §4 and D0 §9.6.
- **Atlas stance unchanged:** nothing here moves `finc1v725`, h405
  candidacy, or the free-form linear adoption gate.

---

## 17. Cross-links

- `TIER2_PARAMETER_INFLUENCE.md` **§9.6** — the mass-partition finding,
  the corrected E_bind decomposition, the density-width null; **§9.2/§9.4/
  §9.5** carry the superseding notes.
- `TIER2_SENSITIVITY_ATLAS_PLAN.md` §6.8 — the RQ12 programme this grew
  out of (T1/T2/T3 closed by the same session).
- `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` — the §6.8 results section.
- `RESEARCH_QUESTIONS.md` RQ12 — density-width provenance.
- `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` — the locked mass
  mechanism this thread questions the *kinematic* premise of; **A2** (no
  imported pickup rate), **S1/S2/K1** (§11), the Langmuir cap.
- `RESEARCH_QUESTIONS.md` **RQ1** — E_int(0) provenance, which §11.3
  upgrades from uncertain to structural for any strip variant.
- `physics/pickup.py`, `physics/internal_energy_budget.py`,
  `physics/exit_strip.py` — the three modules any strip/re-dress build
  would touch; `exit_strip.py` already carries the (C) programme's measured
  over-toll.
- `docs/drag_port/Tier1/TIER1A_IMPLEMENTATION_PLAN.md` — the SQ1–SQ3
  variable-mass integrator and the velocity-preserving shed decision.
- `drag_migration_log_tier2.md` — the 2026-08-11 delivery/decision record.
