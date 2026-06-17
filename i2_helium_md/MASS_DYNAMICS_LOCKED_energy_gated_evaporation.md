# Mass Dynamics — Locked Decision: Discrete-Stochastic Pickup + Energy-Gated (RRK Rate-Limited) Evaporation

> **Entry point:** `DRAG_PORT_DESIGN_DECISIONS.md` (§0 document map, §2 mass summary). This is the **live mass-model detail doc**; parameter provenance index is `CALIBRATION_MAP.md`.

**Status:** Decision locked 2026-06-15. This document fixes the *mechanism* of
the ion-stage mass process for the drag-model port and records why it was
chosen, what it assumes, and what risks it carries. It **supersedes the
"primary = Scenario A, density-driven accretion" framing** of
`DRAG_PORT_DESIGN_DECISIONS.md` §2.6: the production mass model is now a
**two-channel discrete-stochastic process** (Poisson pickup + energy-gated
evaporation), which is the §2.5 *biphasic* structure with a physically grounded
loss channel rather than an abstract relaxation time.

**Scope.** Locks the mechanism and the energy-accounting structure (including the
Newton's-law-of-cooling *form* of internal-energy dissipation, §6 K2, fixed by
[GAH25]). Does **not** lock the numerical rate/ladder values (calibration
targets, §10), the partition fraction $f_\text{ret}$, or the cooling time
$\tau_\text{dissip}$. Maintains the strict Physics-Definition /
Software-Implementation boundary: no code here.

**Revision 2026-06-15 (folded in).** Two refinements to the evaporation
channel, both leaving the two-channel structure intact:
1. **Self-bound gate, parameter-free (R9 resolved).** Evaporation is
   suppressed entirely while the complex is net self-unbound,
   $E_\text{int} > \sum_{i=1}^{n} D_0^{\,\mathrm{I^+}}(i)$ — the regime an
   instantaneous gate would strip wholesale. The threshold is the
   *integrated ladder* already constrained by §6.5.1, **not** a free
   `evap_gate_onset_eV` knob. The several-ps self-unbound window ([GAH25])
   becomes a prediction, not a fit.
2. **RRK rate-limited shedding (R5 reframed, gate-open burst removed).**
   Once self-bound, the top rung sheds as a stochastic unimolecular
   (classical-RRK) process with a *saturating* rate, replacing the
   instantaneous `while`-cascade. This bounds per-step shedding (no
   avalanche when the gate opens) and gives the cascade a physical
   timescale, so completion within 20 ps is a prediction (R5), not an
   assumption. Introduces two kinetics parameters $\{\nu, s\}$ — both
   **sourced/derived rather than free** (see bullet 6, A11).
   See §4, §6 (K1), §8 (R5/R9/R10), §10, §11.
3. **Early-window parameterization (S2 reparameterized; bound-and-check, not
   lock).** $E_\text{int}(0)=f_\text{int}E_\text{avail}$ with $f_\text{int}$
   bounded by the self-unbound floor and 1; $\tau_\text{dissip}$ carried as a
   sweep band $[2.6,16.5]$ ps. The self-bound crossing time $t_\times$ is a
   **derived diagnostic** cross-checked against GAH25's ~5–6.5 ps at
   ±factor-2 (a foreign-system prior, not ground truth) and arbitrated by the
   Tier-2 size distribution. The two scalars are reported as a **regime
   determination** (shell-retaining vs total-vaporization, R1), not fitted
   point values. See §6.11, §10, §11.
4. **K2 cools the GAH25 variable $E_\text{solv.struct}=E_\text{bind}+E_\text{int}$,
   not $E_\text{int}$ alone (correctness fix on reading GAH25 §IV).** Their
   Newton fit (Eq. 12) is to the solvation-structure energy; matching it removes
   a quantity-mismatch error and lets $\tau_\text{GAH25}$ transplant directly.
   Cold shedding is **energy-neutral** for $E_\text{solv.struct}$, so K1 and K2
   do not double-count — but this neutrality is *contingent on the $\sum D_0$
   gate* that suppresses hot ejection (A8). The gate crossing is *identically*
   GAH25's self-bound onset $t_0$. $f_\text{ret}$ is thereby disentangled from
   K2 and is now purely the pickup partition. New risks R11 (≈0.5 ps fit/gate
   overlap, **accepted**) and R12 (hot-structure binding ≠ ladder sum, **gated**
   by A9); see §6 K2, §8, §9 (A8/A9).
The drag form referenced throughout is the locked pure cubic
$F_\text{drag}=g\,b\,v^3$ ($a=0$, $b=2.5154$ amu·ps/Å²); its high-$v$
ceiling caveat is R10.
5. **Binding-ladder electronic picture: statistical SO mixture is the production
   default (A10).** The He–I⁺ ground interaction has a deep ³Π ($\equiv X_2$,
   $D_e=143.9$ cm⁻¹) and a shallow ³Σ⁻ ($63.6$ cm⁻¹) curve ([IHe05]); production
   binds via the equal-weight $X_2{+}I_1{+}I_0$ statistical mixture (shallower),
   motivated by the violent Coulomb-explosion birth, with deep $X_2$-only kept as
   the comparison case. This sets $\sum_i D_0$, the gate, and $t_\times$. The
   supporting $E_\text{bind}$ evidence and the default itself are **provisional
   pending OQ1** (was the drag run on $X_2$ with purely dynamical lowering?) —
   §10A. Label correction: $X_2\equiv³Π\equiv V_\Pi$ (molecular ground state);
   ³P₂ is the atomic sublevel it correlates to. See R3, A10, §10A, §11.
6. **RRK $\{\nu,s\}$ sourced, not free (A11).** $\nu$ gets a stretch-frequency
   prior $\sim\mathcal{O}(1)$ ps⁻¹ from the [IHe05] $X_2$ curvature; $s$ is
   mode-counted $3n-6$ by default (effective-scalar override available),
   calibrated *jointly* with `ladder_shape` since both probe shell rigidity.
   Classical RRK is a stated simplification with $s$ absorbing the RRKM/quantum
   difference. Cascade-timing cross-check from [I2-notes] flagged as OQ5. See
   §4, A11, §10/§10A, §11.

**External references.**

- **[Nat23]** — Albrechtsen, Schouder, Viñas Muñoz et al., *"Observing the
  primary steps of ion solvation in helium droplets,"* Nature **623**, 319
  (2023), doi:10.1038/s41586-023-06593-5 (Stapelfeldt group). Experiment.
  Poissonian He-binding rate 2.0 atoms/ps for Na⁺ at rest; energy-gated
  post-ejection dissociation with a tabulated Na⁺Heₙ ladder.
- **[Calvo24]** — F. Calvo, *"Concurrent processes in the time-resolved solvation
  and Coulomb ejection of sodium ions in helium nanodroplets,"* J. Chem. Phys.
  **161**, 121101 (2024), doi:10.1063/5.0230829. Atomistic ring-polymer MD of the
  *full pump–probe sequence incl. Coulomb ejection*. Supplies fragment-size and
  fragment-temperature distributions and the violent-ejection (total-vaporization)
  limit.
- **[GAH25]** — García-Alfonso, Barranco, Pi, Halberstadt, *"Time-resolved
  solvation of alkali ions in superfluid helium nanodroplets: Theoretical
  simulation of a pump–probe study,"* J. Chem. Phys. **163**, 144309 (2025),
  doi:10.1063/5.0291643. ⁴He-TDDFT of *both* pump and probe steps. Supplies the
  internal-energy budget, the Newton's-law-of-cooling dissipation form, the
  non-monotone shell evolution during ejection, and the early self-instability of
  the solvation structure.
- **[IHe05]** — Buchachenko, Tscherbul, Kłos, Szczęśniak, Chałasiński, Webb &
  Viehland, *"Interaction potentials of the RG–I anions, neutrals, and cations
  (RG = He, Ne, Ar),"* J. Chem. Phys. **122**, 194311 (2005),
  doi:10.1063/1.1900085 (+ erratum J. Chem. Phys. **161**, 149901 (2024)).
  UCCSD(T) ab initio He–I⁺ pair potentials with relativistic small-core
  pseudopotential; ~3% well-depth accuracy from ZEKE/mobility cross-checks
  (firmer than the [I2-notes] working draft). Two ground-correlating NR curves
  (Table IV): **³Π ≡ $X_2$**, $D_e=143.9$ cm⁻¹ at $R_e=3.25$ Å (deep), and
  **³Σ⁻**, $D_e=63.6$ cm⁻¹ at $R_e=3.76$ Å (shallow). SO-coupled ground-sublevel
  states $X_2,I_1,I_0$ (Eq. 9); the **I⁺ transport** treatment (§III.B) uses an
  equal-weight statistical mixture of these three — the basis for the A10 default.
  Source for the first ladder rung $D_0^{\,\mathrm{I^+}}(1)$, the electronic
  picture (A10), and the open-shell / SO structure of the I⁺–He interaction.
- **[I2-notes]** — García-Alfonso, Barranco, Halberstadt, Hauser & Pi,
  *"I₂ molecules in superfluid He nanodroplets,"* **working document** (dated
  Feb 2025, broken eq/fig refs, internal interpolation discrepancy). I⁺-specific
  He-DFT/TDDFT: I⁺(³Π/$X_2$)@He₂₀₀₀ first-shell solvation $S_{\mathrm{I^+}}=-3578$
  K, $n^*=21$ (the $X_2$-only deep-snowball anchor, A10 picture 1); and a direct
  I⁺–I⁺ Coulomb-ejection trajectory from 9 Å (Figs. 21–23: $z(t)$, $d(t)$,
  $v(t)$, peak $v\approx5$–6 Å/ps) — an I⁺-specific kinematic check for the 9 Å
  drag case and the R10 ceiling. Treat numbers as provisional working-draft
  anchors (OQ4).

[Calvo24] and [GAH25] are **closer to this work's regime than [Nat23]** — they
model the *probe/ejection* stage where the ion is moving and the shell is hot and
non-monotone, not just at-rest accretion. They are still alkali, not I⁺ (see §3
and the §13 external-validation table for what transfers vs. warns).

---

## 1. The decision (locked)

The ion-stage helium count $n(t)$ evolves by **two independent discrete
stochastic channels** operating per BAOAB sub-step:

- **Pickup** — a **Poisson process**, rate $\lambda_\text{attach}(\rho_\text{He}(\text{depth}))$, each event $n \to n+1$, $m \to m + m_\text{He}$.
- **Evaporation** — **energy-gated and RRK rate-limited**, governed by two
  conditions on the tracked complex internal energy $E_\text{int}$ and an
  I⁺Heₙ dissociation-energy ladder $D_0^{\,\mathrm{I^+}}(n)$:
  - **Self-bound gate (parameter-free).** Shedding is suppressed entirely
    while $E_\text{int} > \sum_{i=1}^{n} D_0^{\,\mathrm{I^+}}(i)$ (the
    complex is net self-unbound; the droplet holds it together).
  - **RRK shed rate.** Once self-bound and above the top rung
    ($D_0^{\,\mathrm{I^+}}(n) < E_\text{int} < \sum_i D_0^{\,\mathrm{I^+}}(i)$),
    the top rung sheds with a *saturating* unimolecular rate (§4); each
    event $n \to n-1$, $m \to m - m_\text{He}$, $E_\text{int} \mathrel{-}=
    D_0^{\,\mathrm{I^+}}(n)$, shed He leaving **cold** (≈ zero kinetic
    energy). The bounded rate caps per-step shedding — no instantaneous
    cascade.

Mass is therefore **integer-valued in He count**, **non-monotone**, and the
equilibrium shell size is an **emergent balance** between Poisson gain and
energy-gated loss — not a free parameter. Tracking the gate requires a **new
per-ion state variable** $E_\text{int}(t)$ (the complex internal/vibrational
energy), which rides in `IonStepState` and is persisted in a bumped
`IonCheckpoint` v6.

---

## 2. Why this approach — the chain of reasoning

**2.1 The observable demands a generative integer-$n$ process.** The final
experimental comparison is a set of *per-fragment* velocity histograms — one
distribution each for I⁺, I⁺He₁, I⁺He₂, … with a monotone-falling envelope
(bare highest, then $n=1$, then $n=2$, …). A fixed-mass run produces a single
ion species and **cannot generate the fragment channels at all** — it is not
"adequate vs inadequate," it is structurally incapable of producing the
observable. Mass evolution is therefore the *generative process that creates the
histogram channels*, not a correction to a trajectory. This is the load-bearing
reason the null/`fixed` hypothesis cannot be the production model (it remains
valid only as the §6.5 self-consistency reference and a Tier-1 sensitivity
anchor).

**2.2 Discrete-stochastic over continuous-then-bin.** A continuous $\dot M$ ODE
yields a real-valued $m(t_\text{end})$; the histogram needs integer $n$. Binning
a continuous endpoint makes the *spread in $n$* a second-moment quantity driven
by trajectory-to-trajectory variance, risking a near-delta in $n$ if the
deterministic ensemble is tight. **[Nat23] settles this:** the measured $S_t(n)$
distributions (their Fig. 4a) are *broad* at every fixed time, and that breadth
is the **Poisson spread of independent binding events**, not ensemble scatter.
Discrete-stochastic pickup therefore *natively* produces the envelope width the
experiment shows. Continuous-then-bin is rejected.

**2.3 Energy-gated over rate-gated evaporation.** This is the choice made this
session. [Nat23]'s dissociation is **energy-gated, not rate-gated**: a complex
leaving with internal energy $E_\text{int}$ sheds He until $E_\text{int}$ drops
below the dissociation threshold $D_0(n)$; the number of shed atoms $K$ is
whatever the energy budget allows, then "further dissociation is energetically
forbidden." This is the actual physics of the "hot complex evaporates on the way
out" picture that motivated reaching for the paper. A rate-gated alternative
(both channels Poisson) is cheaper and needs no $E_\text{int}$ state, but it
*discards the energy-gating that is the paper's distinctive result* and
reintroduces a free loss-rate parameter. We accept the cost of tracking
$E_\text{int}$ to keep the mechanism faithful.

**Refinement (2026-06-15): RRK rate-limiting is not a regression to
rate-gating.** The §1 locked mechanism originally shed *instantaneously*
(a per-step `while`-cascade down to the first sub-threshold rung). That is
unphysical on two counts: (i) [Nat23]/[Calvo24] show the cascade takes
tens-to-hundreds of ps — evaporation is **rate-limited**, not instant; and
(ii) after the early self-bound gate opens onto a shell grown by pickup, an
instantaneous cascade produces a one-step **avalanche** (a modelling
artifact). Both are fixed by shedding the top rung at a classical-RRK rate
(§4). This is **still energy-gated** — the rate vanishes below $D_0(n)$ and
rises with $E_\text{int}$ — so it does *not* revert to the discarded
rate-gated option (a free, energy-blind Poisson loss rate). The two new
constants $\{\nu, s\}$ are standard unimolecular-dissociation kinetics
parameters (prefactor and effective vibrational DOF), physically meaningful
and calibratable, not an arbitrary loss rate. The energy gate that §2.3
chose over rate-gating is retained; only the *kinetics* of that gated
channel sharpen from instantaneous to finite-rate.

**2.4 Why this is biphasic, and why that is now acceptable.** Gain (pickup) +
loss (evaporation) is the §2.5 biphasic structure. §2.5 held biphasic as
*secondary* because it had "the most degrees of freedom, hardest to calibrate"
and an *undetermined* $M_\text{eq}(\cdot)$ and $\tau$. The energy-gated loss
channel **removes that objection**: the loss side is set by the $D_0^{\,\mathrm
{I^+}}(n)$ ladder and the internal-energy budget, not by a free relaxation time;
$M_\text{eq}$ is emergent, not fitted. Biphasic is promoted to production *on the
strength of the loss channel being physically grounded rather than parametric.*

---

## 3. Provenance from [Nat23] — what transfers, and the regime caveat

**This is a different physical regime, and the distinction governs how the
numbers may be used.**

| | **[Nat23] — Na⁺** | **This work — I⁺** |
|---|---|---|
| Ion creation | at surface dimple, **at rest** | post-Coulomb-explosion, **fast (~10 Å/ps)** |
| Process measured | **accretion** onto stationary ion | net **shell loss** during high-speed traversal (21→19→14 He, 9 Å) |
| Ion mass | 23 amu | 127 amu (bare); $m_\text{eff}\approx203$ amu |
| Ion–He interaction | deep well (~40× He–He) | electrostriction trap, different magnitude |
| Evaporation regime | **post-ejection**, gentle complex | **in-flight + post-ejection**, violent onset |

**Consequence:** [Nat23] supplies **mechanism and parameter scale**, never a
drop-in calibration. Every number imported below is an **order-of-magnitude
anchor (treat as ±factor-of-2)** for the I⁺ case, to be tightened against this
project's own TDDFT shell counts and the experimental I⁺Heₙ size distribution.

**What [Nat23] genuinely transfers:**

1. **Poisson is the right stochastic structure** for independent, constant-rate
   binding (their two stated Poisson criteria: independent binding, constant
   rate).
2. **A rate scale:** 2.0 He/ps at $v=0$ for Na⁺ — a low-velocity anchor.
3. **The evaporation mechanism:** energy-gated cascade shedding until
   $E_\text{int} < D_0(n)$, with a tabulated $D_0(N)$ *structure* (their
   Extended Data Table 1).
4. **Two evaporation facts:**
   - **Timescale:** dissociation completes within "tens of picoseconds," the
     vast majority finished by 400 ps. For a 20 ps ion stage this means
     evaporation **overlaps the simulated window** — it is *not* a separable
     post-process applied after droplet exit.
   - **Shed-He energetics:** shed He leave with very low KE (statistical
     redistribution localizes just enough energy to break one bond). So a loss
     event removes $\approx D_0(n)$ from the internal budget and carries away
     **negligible translational KE** — a clean energy-bookkeeping closure.

**What does NOT transfer (must be re-sourced or calibrated):**

- The numeric $D_0^\text{Na^+}(N)$ values (Na⁺–He specific). The I⁺ ladder
  $D_0^{\,\mathrm{I^+}}(n)$ must come from literature (I⁺Heₙ binding energies) or
  be treated as a calibrated quantity, consistent with §6.5.1's treatment of
  ion binding as an *effective* parameter.
- The 2.0 He/ps magnitude directly (resting Na⁺, light, deep well). Usable only
  as a low-$v$ OOM anchor for $\lambda_\text{attach}$.

---

## 4. The two-channel mechanism (physics definition)

Per BAOAB sub-step of size $dt = dt_\text{ion} = 0.01$ ps:

**Pickup (Poisson → Bernoulli per step):**
$$
P_\text{attach}(dt) = 1 - e^{-\lambda_\text{attach}\,dt} \approx \lambda_\text{attach}\,dt,
\qquad
\lambda_\text{attach} = \lambda_0 \,\frac{\rho_\text{He}(\text{depth})}{\rho_\text{bulk}}.
$$
A single independent Bernoulli draw per ion per step; **no collision gate** (this
is *not* the discarded $b_\text{collision}$-coupled model of baseline §7.1).

**Evaporation (energy-gated + RRK rate-limited):** two conditions on the
top rung, evaluated per step.

*Self-bound gate (R9, parameter-free).* Suppress all shedding while
$$
E_\text{int} > \sum_{i=1}^{n} D_0^{\,\mathrm{I^+}}(i)
\qquad(\text{complex net self-unbound; held only by the droplet}).
$$
This is exactly the regime in which an instantaneous gate would strip the
whole shell; Newton-cooling (§6 K2) drains $E_\text{int}$ below the
integrated binding after several ps, at which point the structure is
self-bound and shedding is allowed. The threshold is the integrated ladder
constrained by §6.5.1 — no free onset parameter.

*RRK shed rate (once self-bound).* For
$D_0^{\,\mathrm{I^+}}(n) < E_\text{int} < \sum_i D_0^{\,\mathrm{I^+}}(i)$,
the top rung sheds as a stochastic unimolecular process with the classical
(Rice–Ramsperger–Kassel) rate
$$
k(E_\text{int}, n) =
\nu\left(1 - \frac{D_0^{\,\mathrm{I^+}}(n)}{E_\text{int}}\right)^{\,s-1},
\qquad
P_\text{shed}(dt) = 1 - e^{-k\,dt},
$$
with $\nu$ a prefactor (ps⁻¹) and $s$ the effective number of vibrational
degrees of freedom of the I⁺Heₙ complex. **Parameter sourcing (2026-06-15):**
- $\nu$ is the **I⁺–He stretch attempt frequency**, with a physical prior from
  the [IHe05] $X_2$ curve curvature,
  $\nu=\tfrac{1}{2\pi}\sqrt{V''(R_e)/\mu}$, $\mu\approx m_\text{He}$ (He-dominated
  reduced mass): for a ~144 cm⁻¹ well at $R_e=3.25$ Å this lands at
  **$\nu\sim\mathcal{O}(1)$ ps⁻¹** (low end of the earlier 1–10 band).
- $s$ is **mode-counted, not free:** $s=3n-6$ (the vibrational DOF of the
  $n$-He shell; $3n-5$ if treated linear), so $s$ is $n$-dependent
  ($\sim57$ at $n{\sim}21$, $\sim9$ at $n{\sim}5$) and parameter-free by
  default, with a single effective-scalar override available for sensitivity
  (mirroring the gate-onset: derive-by-default, override-for-test). The
  mode-count and the ladder shape are **not independent** — for a blurred
  /gradual shell (A10/R3) the effective $s$ may fall below the naïve $3n-6$, so
  $s$ and `ladder_shape` are calibrated *jointly* against the size distribution,
  not separately (both probe the same shell-rigidity question). $s$ also absorbs
  the classical-RRK-vs-RRKM simplification (A11).

*Dimensional check on $\nu$:* $[\sqrt{V''/\mu}]=\sqrt{(\text{amu·Å}^2
\text{ps}^{-2}/\text{Å}^2)/\text{amu}}=\text{ps}^{-1}$ ✓.

Below threshold ($E_\text{int}\le D_0(n)$) the rate is zero. On a shed event:
$n \to n-1$, $E_\text{int} \mathrel{-}= D_0^{\,\mathrm{I^+}}(n)$, shed He
cold. At most one rung is drawn per step; the rate **saturates**
($k \to \nu$ as $E_\text{int} \to \infty$), so per-step shedding is bounded
by $\sim\nu\,dt$ — **no gate-open avalanche**, and the cascade carries a
physical timescale set by $\{\nu, s\}$ (R5 becomes a prediction, not an
assumption).

**Emergent equilibrium:** steady state where the expected Poisson gain rate
equals the RRK-rate-limited shed rate; no $M_\text{eq}$ or $\tau$ parameter.

### Dimensional analysis (mandatory)

| Quantity | Expression | Units | Check |
|---|---|---|---|
| Pickup rate | $\lambda_\text{attach}$ | $\text{ps}^{-1}$ | $[\lambda\,dt]=1$ ✓ |
| Gating density ratio | $\rho_\text{He}/\rho_\text{bulk}$ | dimensionless | ✓ |
| Pickup KE injection | $\tfrac12 m_\text{He} v^2$ | $\text{amu·Å}^2/\text{ps}^2$ = energy | ✓ |
| Dissociation rung | $D_0^{\,\mathrm{I^+}}(n)$ | eV (energy) | ✓ |
| Integrated binding (gate) | $\sum_i D_0^{\,\mathrm{I^+}}(i)$ | eV | compared to $E_\text{int}$ (eV) ✓ |
| RRK prefactor | $\nu$ | $\text{ps}^{-1}$ | $[\nu\,dt]=1$ ✓ |
| RRK bracket | $1 - D_0(n)/E_\text{int}$ | dimensionless | eV/eV ✓ |
| RRK effective DOF | $s=3n-6$ | dimensionless | mode-counted (A11); exponent $s-1$ dimensionless ✓ |
| RRK rate | $k=\nu(\cdot)^{s-1}$ | $\text{ps}^{-1}$ | $[k\,dt]=1$ ✓ |
| Shed-He KE | $\approx 0$ | — | no $\tfrac12 m v^2$ on loss ✓ |
| Internal energy | $E_\text{int}$ | eV | ✓ |

No term fails dimensional balance.

---

## 5. What the references resolve in the existing design docs

- **§2.7 velocity scaling (pickup) → density-only confirmed.** A *resting* ion
  still accretes at a substantial finite rate (2.0/ps at $v=0$, [Nat23]). This
  **rules out pure-sweeping** $\dot M \propto v\,\rho_\text{He}$ (→ 0 at rest,
  contradicts the datum) and is **inconsistent with unregularized dwell-time**
  $\dot M \propto \rho_\text{He}/v$ (diverges at rest). It is **consistent with
  the §2.7 primary, density-only** $\lambda_\text{attach}\propto\rho_\text{He}$.
  Recorded as a second independent argument for density-only.
- **Integer-$n$ spread (earlier open question) → Poisson spread, native.**
  Resolves the continuous-then-bin risk in favour of the discrete channel
  (§2.2 above).
- **The two-channel non-monotone mechanism → independently confirmed at TDDFT
  level ([GAH25]).** During the probe/ejection step the Na⁺He₅ shell count is
  explicitly non-monotone: it *grows* (5→7 in the first ~0.4 ps as the still-
  inward-moving ion sweeps up He) then *oscillates and slowly declines* (to ~4–5
  by end of simulation). [GAH25] attribute this to exactly the two channels
  locked here — collisional/sweeping gain while the ion travels, and
  internal-energy-driven loss ("it can cool down by dissociating He atoms; this
  can occur inside as well as outside the droplet on the way to the detector").
  This is direct theoretical support for (a) discrete gain+loss, (b)
  non-monotonicity (validating the §7 dropped-monotonicity schema change), and
  (c) *energy-gated* rather than rate-gated loss.
- **The $E_\text{int}$ reservoir → measured, with known budget structure
  ([GAH25]).** [GAH25] Eqs. 10–11 write $E_\text{init}=E_\text{dissip}+
  E_\text{bind}(N)+E_\text{int}(N)$, i.e. *this document's §6 budget*, validated
  against TDDFT. Critically, the Newton cooling they fit (Eq. 12, Table III) is
  applied to $E_\text{solv.struct}=E_\text{bind}+E_\text{int}$, **not** to
  $E_\text{int}$ alone — which fixes both the K2 *form and variable* (revision
  folded into §6 K2) and means their self-instability onset $t_0$ is exactly our
  gate crossing $t_\times$ (§6.11). It also supplies the early self-instability
  constraint for S2. *Caveat carried forward:* the energy relaxation
  ($t_0$, $\tau$, $E_\infty$) was computed for **Na⁺ only** — the least
  I⁺-like alkali (§6.11, R8).
- **§6.8 empirical lean toward loss → reinforced, and widened at the violent end.**
  The 9 Å TDDFT 21→19→14 decline is *net* loss; the new model reproduces net loss
  as gain−(energy-gated loss). [Calvo24] adds the violent-ejection bracket: under
  strong Coulomb ejection from a small droplet the cation is ejected so hard that
  *any attached solvent atom is vaporized* (total stripping). Since the production
  I⁺ condition is a Coulomb explosion at $R_e$ (far more violent than the gentle
  9 Å shell trajectory), the loss channel may dominate more than the 9 Å decline
  implies — see R1 and §13.

---

## 6. The internal-energy budget $E_\text{int}(t)$ — structure locked, partition open

Evaporation gating requires a tracked internal-energy reservoir, which a
point-mass MD does not natively have. Its **structure is locked**; the
**partition coefficients are an open calibration item.**

**Sources (heat $E_\text{int}$):**
- **S1 — pickup binding release.** Each attachment releases the bond energy
  $D_0^{\,\mathrm{I^+}}(n{+}1)$ inward; a fraction $f_\text{ret}$ is retained as
  internal heat, the remainder dissipates to the droplet bath. $f_\text{ret}$ is
  **open** (§10).
- **S2 — Coulomb-explosion onset deposit (parameterized as a partition
  fraction).** The initial internal energy is written
  $$
  E_\text{int}(0) = f_\text{int}\,E_\text{avail},
  \qquad f_\text{int}\in[0,1],
  $$
  with $E_\text{avail}$ the energy liberated at the dication onset (the I–I
  Coulomb release at $R_e$, ~5.5 eV/pair, of which almost all goes to
  *translational* dissociation) and $f_\text{int}$ the small fraction that
  couples into shell-*internal* modes. The fraction is preferred over a bare
  $E_\text{int}(0)$ because it is physically interpretable and intrinsically
  bounded. **Two-sided physics bounds (not a hardcoded value):**
  - *Lower:* the one GAH25 fact that survives the alkali→I⁺ /
    gentle→violent regime transfer is the **qualitative** self-unbound onset
    — so $E_\text{int}(0) > \sum_{i=1}^{n_0} D_0^{\,\mathrm{I^+}}(i)$
    (order 0.1–0.3 eV for a ~21-He shell, ladder-shape dependent). The
    mechanism is trusted even where the timescale is not.
  - *Upper:* $f_\text{int}\le 1$ caps it at $E_\text{avail}$; physically
    $f_\text{int}$ is small, giving an illustrative band ~0.1–1 eV.
  The earlier "static-$D_0$ gate strips the shell at $t=0$" hazard is
  **resolved** by the self-bound suppression gate
  ($E_\text{int}>\sum_i D_0$, §4 / R9), not by special-casing S2. What
  remains is that the *crossing time* $t_\times$ (when cooling brings
  $E_\text{int}$ below $\sum_i D_0$ and the gate opens) depends jointly on
  $f_\text{int}$ and $\tau_\text{dissip}$ — handled as a **bounded,
  cross-checked prediction**, not a lock (§6.11).

**Sinks (drain $E_\text{int}$):**
- **K1 — evaporation (RRK rate-limited, self-bound gated).** Each shed
  event consumes one rung $D_0^{\,\mathrm{I^+}}(n)$ from $E_\text{int}$
  (energy spent breaking the bond); shed He leaves cold (≈ 0 KE, per
  [Nat23]). Shedding is suppressed while the complex is net self-unbound
  ($E_\text{int} > \sum_i D_0^{\,\mathrm{I^+}}(i)$) and otherwise proceeds
  at the saturating RRK rate $k(E_\text{int},n)$ (§4) — bounded per step,
  so $E_\text{int}$ drains through this channel gradually rather than in a
  single-step cascade.
- **K2 — bath dissipation, Newton's law of cooling, applied to the GAH25
  variable $E_\text{solv.struct}$ (form *and* variable now fixed by [GAH25]).**
  **Key correction (2026-06-15, on reading GAH25 §IV):** GAH25 fit Newton's
  law not to internal energy but to the **solvation-structure energy**
  $$
  E_\text{solv.struct}(t) = E_\text{bind}(N) + E_\text{int}(N)
  $$
  (their Eqs. 11–12), computed as $\langle\Psi|H|\Psi\rangle$ over a fixed
  sphere; the asymptote $E_\infty$ they report is a *binding* energy
  ($-3424$ K shell-1, $-4144$ K shell-2), not an internal-energy floor. So
  the cooled master variable is $E_\text{solv.struct}$, and
  $$
  \left.\frac{dE_\text{solv.struct}}{dt}\right|_\text{K2}
  = -\frac{E_\text{solv.struct}-E_\infty}{\tau_\text{dissip}},
  \qquad \tau_\text{dissip}=\tau_\text{GAH25},\;\; E_\infty=E_\text{bind}^\text{eq}.
  $$
  Internal energy is recovered as
  $E_\text{int} = E_\text{solv.struct}-E_\text{bind}(N)$ with
  $E_\text{bind}(N) = -\sum_{i=1}^{N} D_0^{\,\mathrm{I^+}}(i)$ — **valid only
  after the self-bound crossing $t_\times$** (locked rule, A9; see R12). The
  single relaxation time $\tau_\text{dissip}$ is now applied to the *correct*
  quantity, removing the quantity-mismatch error of cooling $E_\text{int}$
  alone. $\tau$ is **transplanted directly** ($\tau_\text{K2}=\tau_\text{GAH25}$,
  no inflation) — see the no-double-count finding below. It remains bracketed
  but not pinned (§10, R8): experiment gives $\tau\approx2.6\pm0.4$ ps (Na⁺,
  $\langle N\rangle{=}3600$); TDDFT gives $\tau\approx7.3$ (shell-1) – $16.5$
  (shell-2) ps — treat as a ±factor-3 sweep band, with the residual that these
  are **Na⁺** numbers transplanted to a **Rb⁺-like** I⁺ (§6.11, R8).

  **No-double-count finding (resolves the earlier $\tau_\text{K2}\!\ge\!\tau_\text{GAH25}$
  worry).** With cold shedding (A8), an evaporation event is *energy-neutral*
  for $E_\text{solv.struct}$: it removes $D_0$ from $E_\text{int}$ and adds
  $D_0$ to $E_\text{bind}$ (one fewer bond), the atom carrying ≈0:
  $\Delta E_\text{solv.struct}=(+D_0)+(-D_0)=0$. Therefore K1 (discrete
  evaporative drain) does **not** overlap K2 (continuous bath relaxation) in
  the cooled variable, and GAH25's $\tau$ can be used directly with no
  inflation. This neutrality is **contingent on the cold-shed assumption and
  on the gate that protects it** (A8, R11) — not free-standing.

**Dimensional check (K2):**
$[(E_\text{solv.struct}-E_\infty)/\tau] = \text{eV/ps}$ = power ✓.

**Jump-consistency check (the variable swap leaves S1/K1 unchanged).** The
$E_\text{solv.struct}$ formulation must reproduce the locked S1/K1 rules on
$E_\text{int}$. It does, exactly:
- *Cold-shed (K1):* $E_\text{solv.struct}$ unchanged (neutral, above);
  $E_\text{bind}\to E_\text{bind}+D_0 \Rightarrow E_\text{int}\to E_\text{int}-D_0$
  — the K1 rule. ✓
- *Pickup (S1):* the bath radiates $(1-f_\text{ret})D_0$, so
  $\Delta E_\text{solv.struct} = -(1-f_\text{ret})D_0$; with
  $E_\text{bind}\to E_\text{bind}-D_0$ this gives
  $E_\text{int}\to E_\text{int}+f_\text{ret}D_0$ — the S1 rule. ✓
So only K2's variable changes; S1, K1, and $f_\text{ret}$ are untouched.
$f_\text{ret}$ is now **purely the pickup partition** — no longer entangled
with the cooling channel (which was the §6-internal double-count concern).

**Relationship to the translational drag channel — the double-counting guard.**
The drag force (§3 of the design doc) removes the ion's *translational* KE and
that energy is booked in `E_dissip_eV`. $E_\text{int}$ is the *internal/relative*
energy of the complex and is a **separate reservoir**. They must not draw from
the same energy twice: drag work → `E_dissip`; binding release → $E_\text{int}$
(via $f_\text{ret}$) and `E_dissip` (the remainder). This separation is a
**locked accounting rule**; getting it wrong is the chief energy-balance risk
(§8 R4).

**Global invariant (locked, must close to Verlet drift):**
$$
E_\text{kin} + E_\text{pot} + E_\text{dissip} + E_\text{mass\_transfer} + E_\text{int} \approx \text{const}.
$$
This extends the §2.9 invariant by the new $E_\text{int}$ reservoir.
`E_mass_transfer_eV` (the renamed §2.9 field) absorbs the pickup KE-injection
defect (cold He swept up at ion speed, $-\tfrac12 m_\text{He}v^2$, as in baseline
§7.1.3) and the cold-shed bookkeeping. $E_\text{int}$ remains the tracked
reservoir in the invariant; K2 now cools $E_\text{solv.struct}$, but at fixed
$N$ (between shed events) $E_\text{bind}(N)$ is constant so
$dE_\text{solv.struct}=dE_\text{int}$, and the drained energy is booked to
`E_dissip_eV` (bath) exactly as before — the variable swap changes only the
*driving force* (now the gap to equilibrium *binding*, not to an internal
floor), not the invariant's structure.

### 6.11 The two early scalars, the self-bound crossing, and regime discrimination

The early window is governed by exactly two uncalibrated scalars —
$f_\text{int}$ (S2, sets $E_\text{int}(0)$) and $\tau_\text{dissip}$ (K2, sets
how fast it cools). Neither is locked; both are **bounded by physics and
cross-checked against data**, in keeping with R1/R8 (imported magnitudes are
priors with wide bars, never ground truth).

**Bounds (inputs).**
- $\tau_\text{dissip}\in[2.6,\,16.5]$ ps — the existing R8 bracket, treated as
  a **sweep band**, not a fit target. Externally anchored, so it need not be
  pinned from this work's size distribution.
- $f_\text{int}$ — bounded below by the self-unbound floor and above by 1 (S2).

**The crossing time $t_\times$ is a derived diagnostic, not an input.** Define
$t_\times$ as the first time $E_\text{int}(t)$ falls below
$\sum_i D_0^{\,\mathrm{I^+}}(i)$ — the instant the suppression gate opens.
**This is literally GAH25's $t_0$.** Their self-bound criterion is
$E_\text{solv.struct}=E_\text{bind}+E_\text{int}<0$; our gate suppresses while
$E_\text{int}>\sum_i D_0 \Leftrightarrow E_\text{bind}+E_\text{int}>0$ — the
*same inequality*, crossing zero at the *same instant*. So $t_\times$ is not a
loose analog of $t_0$; it is the same quantity, and GAH25 **measured** it
(Na⁺: $t_0 = 5.0\pm0.1$ ps shell-2, $6.53\pm0.06$ ps shell-1, Table III).
$t_\times$ is reconstructable post-hoc from the `E_int_eV` $(2N,T)$ array
already in the v6 schema (§7), so it carries **zero schema cost**. Cross-check
it, in increasing authority:
1. a broad sanity band — flag only if absurd ($t_\times<1$ ps or $>15$ ps);
2. the GAH25 $t_0$ = 5.0–6.5 ps as a **±factor-2 prior** — but note this is a
   **Na⁺** number and I⁺ is Rb⁺-like (table below), so even a factor-2 offset
   is unsurprising; a factor-10 miss is the genuine flag;
3. the **Tier-2 size distribution** as the actual arbiter — longer suppression
   → larger shell at gate-open → different terminal-$n$ envelope.

GAH25's 5–6.5 ps is thus a Na⁺-anchored prior the model should land *near for
the right reasons*, never a hardcoded threshold.

**Why I⁺ is not Na⁺-like — and why the energy timescale is the least
transferable number.** Converting GAH25's ion–He well depths (their Table I,
in K) and dropping I⁺ in from [IHe05] ($D_e=143.9\ \text{cm}^{-1}=207$ K at
$R_e=3.25$ Å):

| ion | $D_e$ (K) | $R_e$ (Å) | sinking energy (K) |
|---|---|---|---|
| Na⁺ | 410.3 | 2.41 | 4461 |
| K⁺ | 236.5 | 2.90 | 3329 |
| **I⁺** | **207** | **3.25** | (≈ Rb⁺/Cs⁺ range, not computed) |
| Rb⁺ | 204.1 | 3.10 | 3120 |
| Cs⁺ | 168.9 | 3.37 | 2981 |

I⁺ sits **on Rb⁺** in well depth and between Rb⁺/Cs⁺ in size — roughly *half*
Na⁺'s well depth and ~30% less sinking energy. But GAH25 ran the energy
relaxation (the $t_0$, $\tau$ numbers) **only for Na⁺** — the most strongly
bound, hottest, fastest-relaxing case. So both early scalars are imported from
the *worst-matched* alkali. This is the central reason the bound-and-check
stance (vs hardcoding) is correct, and why R8's ±factor-3 $\tau$ band and the
$t_\times$ cross-check are load-bearing, not ceremonial.

**Identifiability and the clean fallback.** $t_\times$ depends only on the
*combination* (how fast $E_\text{int}$ falls from $E_\text{int}(0)$ to
$\sum_i D_0$), so $(f_\text{int},\tau)$ pairs that reach the crossing together
are degenerate *if the size distribution constrained only $t_\times$*. It does
not: $\tau$ keeps acting **after** the gate opens, setting post-crossing
$E_\text{int}(t)$, which drives the RRK bracket $(1-D_0/E_\text{int})$ and hence
the **cascade speed** — a signature in the envelope's spread/tail that
$f_\text{int}$ (which mainly sets *where* the envelope peaks, via accreted shell
size at gate-open) does not carry. Whether the data separates them is an
empirical identifiability check; the safe fallback is to **pin $f_\text{int}$
from the size distribution at a fixed $\tau$, then sweep
$\tau\in[2.6,16.5]$ and confirm terminal $n$ is insensitive.** If insensitive,
the degeneracy is harmless; if not, lean on the tail-shape separation.

**Why this is more than two-scalar bookkeeping — regime discrimination.**
During suppression there is no evaporation drain, so $E_\text{int}$ evolves
under cooling (down) vs pickup heating (up, via $f_\text{ret}D_0$). Two limits
fall out:
- self-binds early → retains a shell → moderate terminal $n$ (the gentle,
  9 Å-like regime);
- pickup heating outpaces cooling → never self-binds within 20 ps → strips at
  exit → the [Calvo24] total-vaporization limit.

So $f_\text{int}$ and $\tau_\text{dissip}$ are the knobs that place I⁺ on the
spectrum between "deep electrostriction trap holds the shell" and "violent
onset strips everything" — **exactly the open physics question flagged in R1**
— and the size distribution tells you which. The calibration is therefore
reported as a **regime determination** (with $t_\times$ and terminal-$n$
envelope as the observables), not a point estimate of two fitted scalars.

**Total stripping is a reachable limit, retained for secondary evaluation (not
the default).** The shell-retaining end is the production-default expectation,
but the [Calvo24] total-vaporization limit is **the far end of this same
biphasic axis, not a separate model** — it is *not ruled out*. It is kept as a
**secondary / sensitivity-run evaluation target** (terminal $n\to$ small, a
narrow histogram near bare I⁺), checked against the Tier-2 size distribution; we
do not pre-judge it. *Caveat to flag (OQ6, §10A):* whether secondary runs can
actually *reach* full strip is contingent on the K2 asymptote $E_\infty=
E_\text{bind}^\text{eq}$ — as locked, the bound-equilibrium target may
mechanically leave a residual shell even when the dynamics drives toward
stripping. This is a thing to **check** when the secondary runs are set up, not
a change made now; if confirmed limiting, $E_\infty$ would need revisiting (no
mechanism change is made here).

---

## 7. Schema and config implications (locked)

- **`IonCheckpoint` → v6** (already required by §2.9 for any non-`fixed`
  scenario). Additional to the §2.9 changes:
  - new per-step field `E_int_eV (2N, T)` — internal-energy trajectory;
  - `mass_history_kg` monotonicity guarantee **dropped** (non-monotone by
    construction);
  - `E_mass_attach_defect_eV` → `E_mass_transfer_eV` (sign covers gain and loss);
  - scenario-metadata field records `mass_scenario = biphasic_energy_gated` so
    downstream tools interpret the arrays correctly.
- Do **not** silently overload existing fields (baseline §12).

---

## 8. Risks (enumerated, with severity)

**R1 — Regime transfer (HIGH).** Every imported number is from alkali solvation,
not I⁺-after-Coulomb-explosion. [Nat23] is at-rest accretion; [GAH25]/[Calvo24]
reach the moving/ejection regime but remain alkali. The pickup rate and the $D_0$
*structure* are borrowed across a regime boundary, and the well-depth axis is
demonstrably *not* transferable (see R-A2 / A2). The one strong-ejection alkali
case ([Calvo24], He₁₀₀₀) *fully vaporizes* the shell, so the production I⁺
regime may sit closer to strip-dominated than the gentle 9 Å decline suggests —
an open physics question (deeper I⁺ electrostriction trap holding the shell vs.
more violent $R_e$ onset stripping it) that only Tier 2 resolves. *Mitigation:*
use the references for mechanism + OOM scale only; calibrate all numbers against
this project's own TDDFT shell counts (Tier 1) and the experimental I⁺Heₙ size
distribution (Tier 2). Treat imported magnitudes as ±factor-of-2 priors.

**R2 — $E_\text{int}$ sourcing (HIGH).** A point-mass MD has no native internal
reservoir; $E_\text{int}(t)$ is *constructed*, and its source partition
($f_\text{ret}$, S2) and drain (K2) are not directly measured. If the
construction is wrong, evaporation fires at the wrong times and terminal $n$ is
miscalibrated. *Mitigation:* the cascade is bounded by the $D_0$ ladder
regardless of partition (energetics cap the shed count), so errors are bounded,
not divergent; partition pinned in Tier 2 against the size distribution.

**R3 — $D_0^{\,\mathrm{I^+}}(n)$ ladder: first rung now sourced, $n>1$ shape
open, electronic picture forked (MEDIUM, partly resolved).** The dimer rung is
in hand from [IHe05]: the ground-state He–I⁺ well depth is **$D_e = 143.9$
cm⁻¹ ≈ 17.8 meV ≈ 207 K at $R_e = 3.25$ Å**. *Electronic-state labelling
(corrected, do not conflate):* this is the **He–I⁺(³Π) NR curve**, which by
[IHe05] Eq. (9) is identical to the molecular SO-coupled ground state
$X_2$ ($V_{X_2}=V_\Pi$); $X_2$ correlates to the *atomic* ground sublevel
I⁺(³P₂). So "³Π," "$X_2$," and "$V_\Pi$" are the **same curve** — the one the
production snowball ([I2-notes]) and (per current understanding) the drag
extraction were built on. The other ground-correlating NR curve is the
**³Σ⁻**, much shallower ($D_e=63.6$ cm⁻¹ at $R_e=3.76$ Å), and it mixes into
the two further SO states $I_1, I_0$. **Whether production binds via $X_2$
alone (deep) or a statistical $X_2{+}I_1{+}I_0$ mixture (shallower) is a real
fork — see A10**, defaulted to the mixture. Two caveats on the number, two on
the ladder:

- *$D_e$ vs $D_0$.* 143.9 cm⁻¹ is the potential well depth $D_e$, not the
  zero-point-corrected dissociation energy $D_0(1)$ the cascade gate needs. He is
  light, so the ZPE is non-negligible; by analogy with the Na⁺–He case (where
  $D_e\approx285$ cm⁻¹ vs the [Nat23] $D_0(1)=270$ cm⁻¹, a ~15 cm⁻¹ gap),
  $D_0^{\,\mathrm{I^+}}(1)\approx125$–$135$ cm⁻¹ ≈ 15–17 meV. Use $D_e$ as the
  upper bound and apply a ZPE correction (or compute $D_0$ from the published
  analytic fit) before wiring the first rung.
- *Open-shell, spin-orbit-split.* Unlike the closed-shell alkali cations, I⁺ is
  ³P₂ and [IHe05] resolves *six* SO-coupled He–I⁺ curves; three correlate with
  the ground ³P₂ sublevel. The single scalar $D_0(1)$ used here is the lowest
  ($X_2$) of these — an approximation that collapses the multiplet (A5).
- *Only the dimer is computed.* [IHe05] is a pair potential; no I⁺Heₙ cluster
  ladder ($n>1$) exists in the literature. The rungs for $n>1$ must be modelled
  (template, below) or folded into the §6.5.1 effective-binding calibration.

**Where I⁺ sits, and the consequent ladder-shape decision.** Converting [IHe05]
to the [GAH25] Table I scale, the He–I⁺ pair well (207 K) sits **between Rb⁺
(204 K) and K⁺ (237 K)** — i.e. squarely in the *heavier-alkali* regime, far
shallower than Na⁺ (410 K) or Li⁺ (852 K). This says: **do not template the
ladder on Na⁺** (the [Nat23] table); template on the K⁺/Rb⁺ end, where [GAH25]
already observed the shell structure is *less* clean (oscillations around linear
binding, "bigger size allowing several He atoms to bind at the same time").

**Open-shell → likely gradual structural evolution without pronounced shells.**
The decisive physical caveat: I⁺ is open-shell and SO-split, like the heavy
noble-gas cations and Pb⁺. For Pb⁺Heₙ the literature finding is that **spin-orbit
coupling causes a gradual structural evolution without pronounced solvation
shells** — the closed-shell "fill shell, sharp drop at closure" picture (Na⁺) is
replaced by a smoother, more monotone binding-energy ladder with no sharp magic
number. If I⁺Heₙ behaves the same way, $D_0^{\,\mathrm{I^+}}(n)$ is a **smoothly
decaying** sequence rather than the flat-then-cliff Na⁺ shape. This would
*simplify* the cascade (no special-cased shell-closure rung) but means borrowing
the Na⁺ ladder shape would be actively misleading.

*Resolution — test and compare, don't commit.* We do not pre-judge shell-
structured vs gradual. Both ladder shapes are kept pluggable behind
`dissociation_ladder` (§11): (i) a **gradual/smooth** decay anchored at the
[IHe05] first rung (the open-shell-consistent first hypothesis), and (ii) a
**shell-structured** K⁺/Rb⁺-templated shape (the closed-shell-analogue fallback).
The Tier-2 size distribution discriminates them: a smooth vs a shell-structured
ladder produces a measurably different terminal-$n$ envelope (a sharp ladder
piles probability at the closure $n$; a smooth ladder spreads it). The §6.5.1
effective-binding calibration constrains the *integrated* ladder
$\sum_i D_0(i)$ regardless of which shape is used. *Mitigation status:* first rung
sourced; shape is a deliberate two-way test, not an unresolved gap.

**R4 — Energy double-counting between drag and $E_\text{int}$ (MEDIUM).** Drag
dissipation and internal-energy heating could erroneously draw the same energy.
*Mitigation:* the §6 locked accounting rule (drag work → `E_dissip`; binding
release → split $E_\text{int}$/`E_dissip`); the global invariant (§6) is the
detector — if it drifts beyond Verlet tolerance, the partition is leaking.

**R5 — Evaporation overlaps the simulated window, and is NOT complete within it
(MEDIUM→HIGH).** Unlike [Nat23]'s post-ejection treatment, evaporation here runs
*during* the 20 ps flight, so the per-step gate must run every step (handled by
construction). But [Calvo24] and [Nat23] agree the cascade is **not finished at
20 ps**: fragments remain hot and unequilibrated at 20 ps, and statistical
evaporation of loosely-bound outer He continues out to *hundreds of ps* (smaller
fragments expected at the experimental timescale). **Therefore terminal $n$ at
the end of the 20 ps ion stage is an *upper bound*, not the detected $n$.**
**Reframed (2026-06-15):** with the RRK rate (§4) the cascade now has a
*physical* timescale set by $\{\nu, s\}$ rather than completing
instantaneously, so "how far the cascade has progressed at 20 ps" is a
genuine model **prediction**, not an artifact of an instantaneous gate.
This does not remove the upper-bound caveat — the cascade still continues
past 20 ps — but it makes the truncation physically meaningful and
calibratable ($\nu, s$ against the [Nat23]/[Calvo24] cascade timescale and
the Tier-2 size distribution).
*Mitigation:* Tier 2 must either (a) append a post-ejection relaxation stage that
runs the energy-gated cascade forward to the experimental timescale, or (b)
compare simulation and experiment at *matched* time. The cascade-order assumption
(shed the loosest/top rung first) is consistent with [Calvo24]'s "outer atoms
preferentially evaporate." Flagged for the Tier 2 comparison routine.

**R6 — Mass↔coefficient inconsistency is STRUCTURAL for production
(MEDIUM→MEDIUM-HIGH).** The locked drag $b=2.5154$ was extracted at constant
$m_\text{eff}\approx203$ amu. The production model `biphasic_energy_gated`
runs a non-monotone $m(t)$ starting at ~211 amu (~21 He) and falling through
$m_\text{eff}$. Therefore **every production run trips the §6.5 consistency
guard** and must set `allow_inconsistent_mass_pairing=True`: the loud-warning
path is the *normal* production path, not an exploratory exception.
*The clean fix is blocked.* Re-extracting drag under the biphasic $m(t)$ via
the design §6.6 option-3 force balance
$F_\text{drag}(t)=m(t)\,a(t)-F_C(R(t))$ needs a time-resolved shell
trajectory, which exists only for the **9 Å case — the one carrying the
transverse-contamination flag**. So a clean extraction and the production
mass model are mutually blocked: the only data that could de-trip the guard
is the flagged case. *Options, none free:* (i) **accept the inconsistent
pairing as production reality**, justified by the design §6.6 argument that
the constant-mass error is mild mid-window (fractional error
$\sim|m(t)-m_\text{eff}|/m_\text{eff}$: near-zero mid-window, $\sim\tfrac13$
at the ends, where tolerances are already loose per §6.7/A7/R10); (ii)
**attempt the 9 Å option-3 re-extraction despite the transverse flag**,
treating transverse contamination as a separate bounded error; (iii)
**defer** — hold both pairings and let the Tier-2 size distribution reveal
whether the inconsistency materially moves terminal $n$. Do **not** flip the
guard default. *Recommended interim:* (i)+(iii); escalate to (ii) only if
Tier-2 shows sensitivity. The §6.6 mid-window argument is the load-bearing
defence — the mass error is smallest exactly where the law is calibrated and
largest only where tolerances are already loose, so the pairing is
*defensible* even though it is structural.

**R10 — Pure-cubic over-braking above the fit window (MEDIUM).** The locked
drag $bv^3$ has no turnover ($a=0,\ b>0$), so the design §3.3 $b<0$ guard is
moot — but the cubic wing extrapolates **super-linearly** beyond the
extraction velocity window $v_\text{max,fit}$: steeper than the power-law
export ($n\approx+2$) and steeper than the generic inertial $\sim v^2$
expectation. The violent Coulomb onset at $R_e\approx2.6$ Å can drive
$v>v_\text{max,fit}$ during the first several ps — exactly the $t^*$-excluded
window with no TDDFT anchor — where the law is simultaneously **uncalibrated
and steepest**, risking transient over-deceleration that mis-sets the initial
conditions for the calibrated mid-flight phase. *Relation to R5:* R5 bounds
the **late** end (cascade incomplete at 20 ps → terminal $n$ upper bound);
R10 bounds the **early** end (drag too steep at onset). Together they bracket
the reliable interval to mid-flight, consistent with the §6.7/A7 several-ps
free-zone. *Mitigation, two tiers:* **(a) tolerance-only (default)** — the
early window is already a loose-tolerance free-zone; accept the over-brake
and rely on the transient washing out before mid-flight (zero new
parameters); **(b) ceiling cap (guard, only if the transient proves
sensitive)** — define $v_\text{ceiling}$ (sourced from the TDDFT trajectories'
peak speed or a generous ceiling; owned here by R10) and hold $\gamma$ constant above it,
$$
\gamma(v)=\begin{cases} b\,v^2 & v\le v_\text{ceiling}\\[2pt]
b\,v_\text{ceiling}^2 & v> v_\text{ceiling},\end{cases}
$$
so the force above ceiling becomes linear, $F=b\,v_\text{ceiling}^2\,v$,
converting the runaway cubic to a linear tail. *Dimensional + continuity
check:* $[b\,v_\text{ceiling}^2]=\text{amu/ps}$ ✓;
$F(v_\text{ceiling}^+)=b\,v_\text{ceiling}^2\cdot v_\text{ceiling}
=b\,v_\text{ceiling}^3=F(v_\text{ceiling}^-)$ ✓ ($C^0$ continuous; $C^1$
kink at the join — acceptable for a guard, not a physical claim).
*Promotion trigger:* Tier-1 9 Å trajectory showing transient $v$ excursions
above $v_\text{max,fit}$.

**R7 — Pickup-rate velocity dependence unknown (LOW–MEDIUM).** The 2.0/ps anchor
is at $v=0$; the in-flight rate at ~10 Å/ps is unmeasured. Density-only assumes
no explicit $v$-dependence beyond what depth/density supply. *Mitigation:*
density-only is the [Nat23]-consistent default (§5); sweeping/dwell-time remain
pluggable (design §2.7) if Tier 1/2 demand $v$-dependence.

**R8 — Newton-cooling $\tau_\text{dissip}$ uncertain at factor-3, and a Na⁺
number on a Rb⁺-like ion (MEDIUM).** The K2 *form and variable* are now fixed
([GAH25] fit $E_\text{solv.struct}$, matched in §6 K2), so the earlier
quantity-mismatch error (cooling $E_\text{int}$ alone) is gone and $\tau$ is
transplanted directly with no inflation. What remains uncertain is the *rate
value* and its *transferability*: GAH25 report $\tau\approx7.3$ ps (shell-1) /
$16.5$ ps (shell-2) for TDDFT vs $2.6$ ps for experiment, and the
shell-definition spread is itself a factor-2; more importantly, all of these
are **Na⁺** — the most strongly bound alkali — while I⁺ is Rb⁺-like (§6.11
table, roughly half Na⁺'s well depth). *Mitigation:* carry $\tau_\text{dissip}$
as a ±factor-3 sweep band $[2.6,16.5]$ ps and confirm terminal $n$ is
insensitive across it (§6.11 fallback); pin against this work's size
distribution if sensitive. Do not adopt any literature value as ground truth.
*Weak reassurance:* the shell-2 fit overlaps the unstable phase *more* than
shell-1 yet has a *longer* $\tau$ (16.5 vs 7.3 ps) — the wrong sign for
hot-ejection contamination inflating the cooling rate (see R11).

**R9 — Static-$D_0$ gate invalid during early self-instability (RESOLVED,
2026-06-15).** [GAH25] shows the shell is net self-unbound for the first
several ps; a naïve per-rung $E_\text{int}>D_0(n)$ gate would strip the
entire shell at onset. *Resolution (parameter-free):* suppress all shedding
while $E_\text{int} > \sum_{i=1}^{n} D_0^{\,\mathrm{I^+}}(i)$ — i.e. exactly
while the complex is net self-unbound (the would-strip-everything regime).
This is the physical self-bound criterion, keyed to the *integrated* ladder
already constrained by §6.5.1, with **no free onset parameter** (the earlier
`evap_gate_onset_eV` knob is retired). It auto-scales: early pickup grows the
shell → larger $\sum_i D_0$ → longer suppression, and Newton-cooling (§6 K2)
brings $E_\text{int}$ below the integrated binding after several ps, at which
point shedding turns on at the RRK rate (§4). The several-ps GAH25 window is
thereby a **prediction** (does the crossing land at ~5–6.5 ps?), and the gate
condition $E_\text{solv.struct}>0$ is *identical* to GAH25's self-bound
criterion, so $t_\times\equiv t_0$ (§6.11) — directly measured, not analogized.
The RRK rate additionally prevents a one-step avalanche at the moment the gate
opens onto the grown shell. *Residual:* the crossing time still depends on
$f_\text{int}$ (S2) and $\tau_\text{dissip}$ (K2) — see §6.11; the *mechanism*
is resolved, the *two early scalars* remain calibration targets (§10).

**R11 — Hot-ejection / Newton-fit boundary overlap (LOW–MEDIUM, ACCEPTED).**
Cold-shed neutrality (A8) underpins the no-double-count $\tau$ transplant, but
it fails for *ballistic* (hot) ejection during the self-unbound phase. In our
dynamics the $\sum D_0$ gate suppresses all shedding through that phase, so the
hazard is quarantined — **except** at the boundary: GAH25's shell-1 Newton fit
starts at 6.0 ps while self-binding completes at $t_0=6.53$ ps, a ~0.5 ps
window where their fit and the unstable phase overlap, so a little hot ejection
may sit inside the data we imported $\tau$ from. *Decision: accepted, not
guarded.* Rationale: (i) the overlap is ~0.5 ps of an 11 ps fit interval;
(ii) the shell-2 vs shell-1 $\tau$ ordering is the wrong sign for significant
hot-ejection inflation (R8); (iii) we deliberately do **not** add a hard guard
decoupling the fit window from the gate, to avoid complexity for a bounded
boundary effect. *Open item (cannot resolve from the paper):* the actual
$KE_\text{shed}$ distribution over 5–6.5 ps, which sets the leak magnitude —
needs GAH25 movies/supplementary or a value from the authors.

**R12 — Hot-structure binding $\neq$ static ladder sum (LOW–MEDIUM, GATED).**
The $E_\text{int}$ reconstruction $E_\text{int}=E_\text{solv.struct}+\sum_i
D_0(i)$ assumes equilibrium binding; the hot early structure's true
$E_\text{bind}(N)$ deviates from $-\sum_i D_0(i)$ (multi-peak, distorted,
GAH25 Fig. 3). *Mitigation (locked rule, A9):* never reconstruct $E_\text{int}$
from the ladder inside the gate window; trust only $E_\text{solv.struct}$
there. The gate enforces this automatically. *Open item:* quantify the
$E_\text{bind}(N)$ vs $-\sum_i D_0$ deviation right at $t_\times$, where the
structure is self-bound but not yet equilibrated ($E_\infty$ reached only at
~11+ ps in GAH25) — bounds the reconstruction error at the one instant it is
first used.

---

## 9. Assumptions (enumerated, each with how to falsify/tighten)

**A1 — Pickup is Poisson (independent, rate set by local density).** *Tighten:*
check that simulated per-step attachment statistics reproduce a Poissonian
$P_n(t)$ shape against the [Nat23] form; falsified if binding shows memory/
correlation.

**A2 — Only the pickup *mechanism* (Poisson, density-driven) transfers; the
*rate* does not (DOWNGRADED).** The earlier assumption that accretion rate is
approximately ion-mass-independent is **explicitly contradicted on the
well-depth axis** by [GAH25]: their Table I shows the Ak⁺–He well depth falling
Na⁺ 410 K → Cs⁺ 169 K, and the measured solvation rate tracks it (Na⁺ fastest,
Cs⁺ slowest). The rate is also model/grid-dependent within TDDFT (Na⁺ 1.33/ps at
fine grid vs [0.74, 0.79] coarse) and differs from RPMD (~2× larger
coordination). Since I⁺–He sits in a different (much deeper, electrostriction)
regime than any alkali here, **no alkali rate can be lifted directly** — even the
*ordering* argument forbids it. *Consequence:* import only the *form* (Poisson,
$\propto\rho_\text{He}$); calibrate $\lambda_0$ *entirely* against this work's I⁺
TDDFT shell build-up and the size distribution. The 2.0/ps figure is a loose
upper-bound-ish OOM sanity check, nothing more.

**A3 — Shed He leave cold (≈ 0 KE).** Directly from [Nat23]'s statistical-
dissociation argument. *Tighten:* sanity-check that relaxing this (giving shed He
thermal KE at 0.4 K) changes terminal $v$ negligibly — expected null given the
tiny bath energy.

**A4 — $E_\text{int}$ is a single scalar reservoir per ion.** Collapses all
internal/vibrational/rotational modes to one number. *Tighten:* adequate if
terminal-$n$ statistics match Tier 2; falsified if the size distribution
requires mode-resolved structure.

**A5 — $D_0^{\,\mathrm{I^+}}(n)$ exists as a well-defined scalar ladder
(QUALIFIED).** Assumes a single ground-state binding ladder analogous to
[Nat23] Eq. (1), $E_\text{bind}(n)=\sum_i D_0(i)$, collapsing the I⁺(³P₂)
spin-orbit multiplet ([IHe05] resolves six He–I⁺ curves) to one scalar rung per
$n$. *Two things being tested rather than assumed:* (i) the *shape* — smooth
gradual decay (open-shell/Pb⁺-like, the first hypothesis) vs flat-then-cliff
shell structure (closed-shell/alkali-like) — discriminated by the Tier-2 size
distribution (R3); (ii) the scalar-multiplet collapse — which **electronic
picture** sets the rung depth (statistical mixture vs $X_2$-only) is itself a
decision (A10), and the single-scalar collapse is adequate if terminal-$n$
statistics match Tier 2, falsified if the size distribution shows structure only
an SO-resolved ladder reproduces. *Tighten:* source/compute I⁺Heₙ rungs for
$n>1$; cross-check the integrated ladder against the §6.5.1 effective binding;
resolve OQ1 (§10A).

**A6 — Pickup and evaporation are independent channels.** Both draw per step
without cross-gating. *Tighten:* physical; the only coupling is via $E_\text{int}$
(pickup heats it, evaporation drains it), which is explicit in §6.

**A7 — The early transient is uncalibrated AND self-unbound for *several* ps, not
~0.5 ps (EXTENDED).** No TDDFT anchors the violent onset for *either* drag or
mass, and [GAH25] shows the solvation structure is net self-unbound for ~5–6.5 ps
(held together only by the droplet). *Consequences:* (i) loose validation
tolerances on $v(t)$, $R(t)$, $E_\text{int}(t)$ across this window under all
scenarios; (ii) the evaporation gate is **suppressed automatically** here by
the self-bound criterion $E_\text{int}>\sum_i D_0$ (R9 resolved), so it no
longer strips the shell at onset; the crossing out of this window
($t_\times$) is a cross-checked prediction (§6.11). This widens the §6.7
Scenario-A-only ~0.5 ps concession to a several-ps, all-scenario free zone
for the *drag/tolerance* side (R10), while the *mass* side is now governed by
the parameter-free gate rather than a loose tolerance.

**A8 — Cold-shed energy-neutrality for $E_\text{solv.struct}$, contingent on
the gate (NEW, load-bearing).** Two sub-claims: (1) a shed atom carries ≈0 KE
(A3, [Nat23] statistical evaporation), and (2) bond breaking moves exactly
$D_0$ from $E_\text{int}$ to $E_\text{bind}$ with nothing else moving.
Together they give $\Delta E_\text{solv.struct}=0$ per shed event, which is
what licenses transplanting GAH25's $\tau$ onto K2 with no double-count (§6
K2). **This is not free-standing — it is protected by the $\sum D_0$ gate.**
The neutrality fails during the violent early phase, where atoms can leave
*ballistically* with real KE ($\Delta E_\text{solv.struct}=-KE_\text{shed}\neq0$,
an un-modelled extra cooling); but the gate suppresses *all* shedding while
$E_\text{solv.struct}>0$, i.e. through exactly that phase, so no hot-ejection
event fires in our dynamics. Thus the gate and cold-shed are a **coupled pair
of assumptions**, each load-bearing for the other: the gate makes cold-shed
valid, and cold-shed makes the $\tau$-transplant valid. *Residual soft spot:*
the GAH25 shell-1 Newton fit starts at 6.0 ps while self-binding completes at
$t_0=6.53$ ps — a ~0.5 ps window where their fit and the unstable phase
overlap (R11, accepted). *Tighten:* obtain the $KE_\text{shed}$ distribution
during 5–6.5 ps (movies/supplementary or from Halberstadt) to bound the leak.

**A9 — Static-ladder $E_\text{int}$ reconstruction is valid only after
$t_\times$ (NEW, locked rule).** Recovering $E_\text{int}=E_\text{solv.struct}
-E_\text{bind}(N)$ uses $E_\text{bind}(N)=-\sum_i D_0(i)$, the *equilibrium*
ladder. For a hot, distorted, multi-peak early structure (GAH25 Fig. 3) the
true binding $E_\text{bind}(N)\neq-\sum_i D_0(i)$, so the reconstruction errs
exactly in the early window. **Rule:** never use the per-rung $E_\text{int}$
vs $D_0(n)$ comparison (or the ladder reconstruction) inside the gate window
($E_\text{solv.struct}>0$); only $E_\text{solv.struct}$ itself is trusted
there. Post-$t_\times$, as the structure relaxes toward equilibrium, the
ladder sum is valid. The gate enforces this automatically (R12). *Tighten:*
check how far $E_\text{bind}(N)$ for the hot structure deviates from
$-\sum_i D_0$ right at $t_\times$, where the structure is self-bound but still
far from equilibrium ($E_\infty$ not reached until ~11+ ps in GAH25).

**A10 — Binding-ladder electronic picture: statistical SO mixture (production
DEFAULT), $X_2$-only (documented alternative) (NEW).** I⁺ is open-shell, atomic
ground term ³P, ground SO sublevel ³P₂. The He–I⁺(³P) NR interaction gives two
curves ([IHe05] Table IV): **³Π ≡ $X_2$** ($D_e=143.9$ cm⁻¹, $R_e=3.25$ Å,
deep) and **³Σ⁻** ($D_e=63.6$ cm⁻¹, $R_e=3.76$ Å, shallow). SO coupling yields
six molecular states; three — $X_2, I_1, I_0$ — correlate to the ground ³P₂
sublevel, with $V_{X_2}=V_\Pi$ and $V_{I_1},V_{I_0}$ mixing in the shallower
$V_\Sigma$ ([IHe05] Eq. 9). Two pictures for the production binding ladder:

- **(1) $X_2$-only (adiabatic ground, deep "snowball").** Bind via the single
  deepest curve. This is what the [I2-notes] He-DFT used → first-shell solvation
  $S_{\mathrm{I^+}}=-3578$ K $=-0.308$ eV, $n^*=21$. Predicts a *deep* integrated
  ladder.
- **(2) Statistical mixture $X_2{+}I_1{+}I_0$, equal weight.** This is
  Buchachenko's own treatment for I⁺ **transport** ([IHe05] §III.B): ions
  created with enough energy populate a statistical mix of the three
  ground-sublevel SO states. Mixing in $V_\Sigma$ makes the *effective* binding
  **shallower** than $X_2$ alone.

**DEFAULT = statistical mixture (2).** Rationale: production I⁺ is born from I₂
double-ionization / Coulomb explosion — a violent, high-energy creation, exactly
the regime [IHe05] invokes to justify the statistical SO population for
transport. The $X_2$-only picture (1) is retained as the **comparison case**
(it matches the [I2-notes] snowball structure and is the natural deep-binding
bound).

*Supporting evidence — and an honest contingency (see OQ1, §10A).* The
trajectory-matched effective binding $E_\text{bind}=0.1168$ eV sits a factor
~2.6 **below** the $X_2$-only first-shell solvation ($0.308$ eV). The direction
(production binding ≪ $X_2$ snowball) is *consistent* with the shallower mixture
picture — so the electronic-picture question and the $S_{\mathrm{I^+}}$-vs-
$E_\text{bind}$ reconciliation are partly the **same** question. **However**, if
the drag trajectory was run on $X_2$ and the lowering to 0.1168 eV is purely a
*dynamical* effect (hot, fast, non-equilibrated passage — the current working
hypothesis, OQ1), then 0.1168 eV is **not independent evidence** for the mixture,
and adopting a shallower static mixture ladder *on top of* a dynamically-lowered
binding would risk **double-counting** the same reduction. The mixture default
therefore rests on the birth-violence physics; the $E_\text{bind}$ corroboration
is **provisional pending OQ1**.

*Consequence.* The picture sets the integrated ladder $\sum_i D_0(i)$, hence the
self-bound gate threshold and $t_\times$ (§4, §6.11): the mixture gives a
shallower $\sum_i D_0$ → lower gate, earlier crossing. *Falsify/tighten:*
integrated ladder under each picture vs the §6.5.1 effective binding; the
terminal-$n$ envelope (Tier 2); and confirmation from the authors (OQ1).
*Config:* `ladder_electronic_picture ∈ {statistical_mixture (default),
x2_only}` (§11).

**A11 — Evaporation kinetics are classical RRK, with $s$ absorbing the
simplification (NEW).** The shed rate (§4) uses classical Rice–Ramsperger–Kassel
form, $k=\nu(1-D_0/E_\text{int})^{s-1}$, the simplest unimolecular kinetics.
This deliberately ignores the *quantum* mode structure and zero-point energy of
the I⁺Heₙ complex — which matters here because He modes are soft and the system
is cold, so the more defensible form would be RRKM (explicit density/sum of
states). RRKM is not adopted because it needs the I⁺Heₙ vibrational spectrum,
which does not exist ($n>1$, A10/R3). *Position:* classical RRK is a stated
simplification; the **effective DOF $s$ absorbs the RRKM/quantum/ZPE
difference** — it is mode-counted at $3n-6$ by default but calibratable as an
effective scalar, so the classical form is a parametrized stand-in, not a
first-principles claim. *Coupling (flagged):* $s$ and the ladder shape (A10/R3)
both probe shell rigidity/separability and are **not independent** — a
blurred/gradual shell lowers the effective $s$ below $3n-6$; calibrate them
**jointly** against the size distribution. *Cross-check / tighten:* the I⁺
cascade timescale from [I2-notes] (OQ5, §10A) pins the $\{\nu,s\}$ combination
the same way GAH25's window pinned $\{f_\text{int},\tau\}$; falsified if the
size distribution requires an RRKM-shaped, non-power-law switch-on.

---

## 10. Open items / calibration targets

| Item | Symbol | Source / target | Tier |
|---|---|---|---|
| Pickup rate coefficient | $\lambda_0$ | OOM prior only ([GAH25] well-depth dep.); pin from I⁺ TDDFT + size dist. | 1 / 2 |
| First ladder rung | $D_0^{\,\mathrm{I^+}}(1)$ | **sourced:** [IHe05] $D_e=143.9$ cm⁻¹ at $R_e=3.25$ Å, minus ZPE → ~125–135 cm⁻¹ | — |
| Ladder shape ($n>1$) | gradual vs shell-structured | open-shell→gradual (1st try) vs K⁺/Rb⁺ template; discriminated by size dist. | 2 |
| Integrated ladder | $\sum_i D_0(i)$ | cross-checked against §6.5.1 effective binding | 2 |
| Binding-release retained fraction | $f_\text{ret}$ | size distribution | 2 |
| Newton-cooling relaxation time | $\tau_\text{dissip}$ | **sweep band $[2.6,16.5]$ ps** (R8), externally anchored; not pinned from this work's size dist. | sweep |
| Newton-cooling asymptote (binding) | $E_\infty=E_\text{bind}^\text{eq}$ | equilibrium-shell binding ([GAH25] Table III: $-3424$/$-4144$ K for Na⁺ shell-1/2); compute for I⁺ shell | 2 |
| Onset partition fraction | $f_\text{int}$ | $E_\text{int}(0)=f_\text{int}E_\text{avail}$ (§6.11/S2); bounded by self-unbound floor and 1; pin from size dist. | 2 |
| Self-bound crossing time | $t_\times$ | **derived diagnostic, not fitted**; cross-check vs GAH25 ~5–6.5 ps (±factor-2) and size dist. (§6.11) | — |
| Early-instability gate | **derived, not fitted** | self-bound criterion $E_\text{int}<\sum_i D_0(i)$ (R9); the ~several-ps onset is now a *prediction* vs [GAH25], cross-checked by size dist. | — |
| RRK prefactor | $\nu$ | **stretch-frequency prior** $\sim\mathcal{O}(1)$ ps⁻¹ from [IHe05] $X_2$ curvature (§4); cross-check vs [I2-notes] cascade timing (OQ5) + size dist. | 2 |
| RRK effective DOF | $s$ | **mode-counted $s=3n-6$, not free** (A11); effective-scalar override for sensitivity; calibrated *jointly* with `ladder_shape` (A11 coupling) | derived |
| Pickup $v$-dependence (if needed) | sweeping/dwell | only if density-only fails Tier 1/2 | 1 / 2 |
| Total-stripping limit (Calvo24) | terminal $n\to0$ | reachable far end of the biphasic regime axis (§6.11); **evaluated in secondary/sensitivity runs, not excluded, not default**; check vs size dist. (and OQ6) | secondary |

**Adjudicating observable:** the experimental I⁺Heₙ velocity-histogram set
(per-fragment, falling envelope) is the discriminating Tier 2 target, **compared
at matched time or after a post-ejection relaxation stage** (R5). The Wasserstein
metric on integer-$n$ support (design §6.9 primary) is the natural divergence
measure for the size distribution.

### 10A. Open questions / provenance flags (awaiting external confirmation)

These are unresolved at the level of *provenance*, not modelling — they depend
on information from the source authors (García-Alfonso / Halberstadt / Barranco
/ Pi) rather than on calibration. **User is in direct contact with the authors;
update on confirmation.**

- **OQ1 — Electronic-state provenance of the drag extraction (HIGH priority for
  A10).** The locked drag ($b=2.5154$, effective binding $E_\text{bind}=0.1168$
  eV) was, on current understanding, extracted from a trajectory run on the
  **$X_2$ (³Π) curve**, with the low effective binding arising from *dynamical*
  effects (the ion never equilibrates to the deep snowball during its hot, fast
  passage) — **not** from a shallower input curve. *If confirmed:* (i) 0.1168 eV
  is an $X_2$-input number reduced by dynamics, so it is **not** independent
  evidence for the A10 statistical-mixture picture; (ii) there is a
  **double-count hazard** — a shallower static mixture ladder must not be stacked
  on top of an already dynamically-lowered binding; the A10 default would then
  rest *only* on the birth-violence argument, and the mixture-vs-$X_2$ choice
  should be re-examined for whether the dynamics already delivers the mixture-like
  effective depth. *If instead the extraction already SO-averaged:* 0.1168 eV is
  mixture-like and self-consistent with the A10 default. **Until resolved:** treat
  the A10 mixture default as physically motivated but provisional; do not cite
  $E_\text{bind}$ as independent support for it; carry both possibilities.
- **OQ2 — $KE_\text{shed}$ distribution over 5–6.5 ps** (bounds R11 hot-ejection
  leak; from GAH25 movies/supplementary).
- **OQ3 — $E_\text{bind}(N)$ vs $-\sum_i D_0$ deviation at $t_\times$** (bounds
  R12 reconstruction error at first use).
- **OQ4 — $S_{\mathrm{I^+}}=-3578$ K and $n^*=21$ provenance** ([I2-notes] is an
  undated working draft with broken refs and an internal interpolation
  discrepancy; trust the numbers as $X_2$-only structural anchors, but confirm
  against a published version when available).
- **OQ5 — I⁺ He-ejection timing for the $\{\nu,s\}$ cross-check.** [I2-notes]
  Figs. 20–23 give the I⁺ *ion* kinematics ($z,d,v$ vs $t$) but not obviously a
  He-count-vs-time trace. The cascade-timing anchor that would pin the
  $\{\nu,s\}$ combination (how many He leave over the ~5 ps window, per-event
  spacing) may live only in the supplementary **movies**. *If extractable:*
  in-hand cross-check for $\nu\sim\mathcal{O}(1)$ ps⁻¹ and the effective $s$.
  *If movie-only:* author-contact item. Confirm whether $N_\text{He}(t)$ is
  readable.
- **OQ6 — Does $E_\infty$ permit full stripping in secondary runs?** The
  total-vaporization limit (Calvo24) is retained as a secondary-run evaluation
  target (§6.11, §10 table), but the K2 asymptote $E_\infty=E_\text{bind}^\text{eq}$
  is a *bound-equilibrium* target that may mechanically leave a residual shell
  even when the dynamics drives toward bare I⁺. **Check at secondary-run
  setup** whether the locked $E_\infty$ caps terminal $n$ above 0; if so, a
  regime-dependent or shell-tracking $E_\infty$ would be needed to let the
  stripping end be genuinely reachable. Flagged only — no mechanism change made
  now.

---

## 11. Interchangeability surface (config fields)

Extends `DRAG_PORT_DESIGN_DECISIONS.md` §2.8:

- `SimConfig.mass_scenario` gains value `biphasic_energy_gated` (production).
  Existing `fixed`, `scenario_A_accretion`, `scenario_B_stripping`, `biphasic`
  retained for comparison/regression.
- `SimConfig.mass_initial_amu` — initial physical mass; for
  `biphasic_energy_gated` defaults to the measured ion-stage-onset shell
  (~21 He ≈ 211 amu, design §2.1), distinct from $m_\text{eff}\approx203$ amu.
- `SimConfig.pickup_rate_coefficient` — $\lambda_0$ ($\text{ps}^{-1}$).
- `SimConfig.pickup_rate_form ∈ {density_only, sweeping, dwell_time}` — default
  `density_only` ([Nat23]-supported, §5).
- `SimConfig.dissociation_ladder` — $D_0^{\,\mathrm{I^+}}(n)$ bundle (eV per
  rung), or pointer to the effective-binding calibration. First rung anchored at
  [IHe05] (ZPE-corrected); for $n>1$ the shape is selected by:
- `SimConfig.ladder_electronic_picture ∈ {statistical_mixture, x2_only}` —
  `statistical_mixture` (equal-weight $X_2{+}I_1{+}I_0$, shallower; the
  production default per A10, motivated by violent Coulomb-explosion birth) vs
  `x2_only` (deep $X_2/³Π$ snowball, the [I2-notes]-consistent comparison). Sets
  the rung depths and hence $\sum_i D_0$, the gate, and $t_\times$. Default
  provisional pending OQ1 (§10A).
- `SimConfig.ladder_shape ∈ {gradual, shell_structured}` — `gradual` (smooth
  decay, open-shell/Pb⁺-consistent) is the first hypothesis; `shell_structured`
  (K⁺/Rb⁺-templated, flat-then-cliff) is the closed-shell-analogue alternative.
  Discriminated by the Tier-2 size distribution (R3, A5). Not pre-judged.
- `SimConfig.internal_energy_retained_fraction` — $f_\text{ret}\in[0,1]$.
- `SimConfig.internal_energy_cooling_tau_ps` — $\tau_\text{dissip}$, the
  Newton's-law-of-cooling relaxation time for $E_\text{solv.struct}$ (§6 K2;
  replaces the earlier generic "bath rate"). Carried as the $[2.6,16.5]$ ps
  sweep band (R8).
- `SimConfig.solv_struct_asymptote_eV` — $E_\infty=E_\text{bind}^\text{eq}$,
  the equilibrium-shell *binding* energy that $E_\text{solv.struct}$ relaxes
  toward (§6 K2). NOT an internal-energy floor (corrects the earlier
  `internal_energy_equilibrium_eV`, which is retired).
- `SimConfig.internal_energy_partition_fraction` — $f_\text{int}\in[0,1]$,
  the fraction of onset energy coupled into shell-internal modes;
  $E_\text{int}(0)=f_\text{int}\cdot E_\text{avail}$ (§6.11/S2).
  Calibration target (§10), bounded below by the self-unbound floor.
- `SimConfig.coulomb_available_eV` — $E_\text{avail}$, the onset energy
  liberated at $R_e$ (I–I Coulomb, ~5.5 eV/pair scale). Sets the upper
  reference for $E_\text{int}(0)$; sits in the uncalibrated $t^*$ window so
  it is a fixed reference, not a fit target.
- `SimConfig.internal_energy_initial_eV` — **retired** as a free field;
  derived as $f_\text{int}\cdot E_\text{avail}$. Retained only as an optional
  manual override (default: unset → use the derived value) for sensitivity
  exploration.
- `SimConfig.evap_rate_prefactor_per_ps` — RRK prefactor $\nu$ (ps⁻¹), the
  I⁺–He stretch attempt frequency; physical prior $\sim\mathcal{O}(1)$ ps⁻¹ from
  the [IHe05] $X_2$ curvature (§4), cross-checked vs [I2-notes] cascade timing
  (OQ5). Governs the shed rate once self-bound.
- `SimConfig.evap_rrk_dof` — RRK effective DOF $s$, the exponent $s-1$ in the
  rate (§4). **Derived by default as $s=3n-6$** (mode-count, $n$-dependent,
  A11); set to a fixed scalar only to override for sensitivity. Calibrated
  jointly with `ladder_shape` (A11 coupling), not independently.
- `SimConfig.evap_gate_onset_eV` — **retired** as a free knob (R9 resolved).
  The suppression threshold is now the *derived* integrated ladder
  $\sum_i D_0^{\,\mathrm{I^+}}(i)$. Retained only as an optional manual
  override (default: unset → use the derived $\sum_i D_0$ criterion); set it
  only for deliberate sensitivity exploration.
- `SimConfig.allow_inconsistent_mass_pairing` (default `False`, inherited §6.5) —
  required `True` to run `biphasic_energy_gated` against constant-$m_\text{eff}$
  drag coefficients (the R6 exploratory path).

---

## 13. External validation — what each reference confirms vs. warns

Traceability for every imported claim, so provenance of the mechanism, the
Newton-cooling form, $\tau$, and the downgraded A2 is auditable.

| Reference | Regime | **Confirms (supports the lock)** | **Warns (constrains calibration)** |
|---|---|---|---|
| **[Nat23]** Albrechtsen, Nature 623 (2023) | Na⁺, at-rest accretion (pump only) | Poisson pickup structure; energy-gated dissociation cascade; cold-shed He (≈0 KE); broad $S_t(n)$ → native integer-$n$ spread; $D_0$ *ladder structure* | Rate 2.0/ps is at-rest, light, deep-well Na⁺ — OOM only; $D_0$ *values* are Na-specific (not transferable) |
| **[Calvo24]** Calvo, JCP 161 (2024) | Na⁺, full pump+probe incl. ejection (RPMD) | Non-monotone shell during ejection; fragment-size *and* fragment-temperature distributions as the discriminating observable; outer (loose) He evaporate preferentially → cascade order; violent-ejection → total vaporization bracket | Cascade not complete at 20 ps → terminal $n$ is upper bound (R5); strong-ejection alkali strips *everything* → I⁺ loss may dominate more than 9 Å decline implies (R1); RPMD lacks superfluidity/exchange |
| **[GAH25]** García-Alfonso/Halberstadt, JCP 163, 144309 (2025) | Na⁺,K⁺,Rb⁺,Cs⁺ solvation + Na⁺ pump–probe (⁴He-TDDFT, He₂₀₀₀) | The §6 budget (Eqs. 10–11); **Newton's-law-of-cooling applied to $E_\text{solv.struct}=E_\text{bind}+E_\text{int}$** (Eq. 12, Table III: $t_0$=6.53/5.0 ps, $\tau$=7.3/16.5 ps, $E_\infty$=−3424/−4144 K for shell-1/2) → fixes K2 *form and variable*; self-bound onset $t_0\equiv$ our gate crossing $t_\times$; cold-shed → no K1/K2 double-count (A8); non-monotone gain-then-loss confirmed | **Energy analysis is Na⁺-only** — the least I⁺-like alkali; I⁺ is Rb⁺-like in well depth (207 vs 204 K), so $t_0$/$\tau$ are imported from the worst-matched case (R8, §6.11); well-depth-dependent rate → A2 (rate not transferable); ~0.5 ps fit/unstable-phase overlap (R11); hot-structure binding ≠ ladder sum early (R12); still alkali, not I⁺ |
| **[IHe05]** Buchachenko/Viehland, JCP 122 (2005) | He–I⁺ ab initio **pair** potential (the only I⁺ source) | First ladder rung $D_0^{\,\mathrm{I^+}}(1)$: $D_e=143.9$ cm⁻¹, $R_e=3.25$ Å; places I⁺ between Rb⁺/K⁺ → ladder template choice; validated vs I⁺ mobility | $D_e$ not $D_0$ (ZPE correction needed); open-shell ³P₂, SO-split (six curves) → scalar collapse (A5) and likely *gradual, shell-less* ladder (R3); **no $n>1$ cluster data** |

**Net effect on the lock:** the *mechanism* (two-channel discrete, Poisson
pickup, energy-gated loss, $E_\text{int}$ reservoir) is **strengthened** — two
independent theory papers reproduce its qualitative structure. The
*calibration-side* claims are **tightened and made more honest:** K2 form fixed,
$\tau$ and A2 explicitly bracketed/downgraded, and the transient free-zone widened
to several ps with a concrete correctness constraint (R9). No mechanism change.

---

## 14. Cross-references

- Baseline mass attachment being replaced: `PHYSICS_BASELINE.md` §7.1
  (collision-gated Bernoulli — discarded; this model is its clean
  re-expression).
- Mass-scenario taxonomy and §2.9 schema/energy invariant:
  `DRAG_PORT_DESIGN_DECISIONS.md` §2 (esp. §2.5 biphasic, §2.7 velocity scaling,
  §2.8 surfaces, §2.9 schema, §2.10 validation).
- Mass↔coefficient and drag↔binding consistency guards:
  `DRAG_PORT_DESIGN_DECISIONS.md` §6.5, §6.5.1.
- Validation hierarchy (Tier 1 mass scenario, Tier 2 size distribution, Tier 3
  ensemble): §6.4. Histogram metric: §6.9.
- Transient free-extrapolation zone: §2.3, §6.7 (widened to several ps here, A7/R9).
- External:
  - [Nat23] — Poisson rate 2.0/ps (Fig. 3b), energy-gated dissociation
    (Methods, Eqs. 1–6), $D_0^\text{Na^+}(N)$ ladder (Extended Data Table 1),
    cold-shed argument (Methods), dissociation timescale (Extended Data Fig. 4).
  - [Calvo24] — fragment-size distributions and violent-ejection vaporization
    (Fig. 4, p.4 text); fragment temperatures and outer-shell evaporation
    (Fig. 4b, p.4–5); non-monotone capture during ejection.
  - [GAH25] — $E_\text{int}$ budget Eqs. 10–11; Newton's-law cooling Eqs. 9, 12
    and Table III ($\tau$, $t_0$, $E_\infty$); early self-instability (Fig. 7,
    §IV); well-depth-dependent rate (Table I, Table II); non-monotone $n_1(t)$
    during probe (Fig. 6).
  - [IHe05] — He–I⁺ pair potential, ground ³P/$X_2$ component $D_e=143.9$ cm⁻¹,
    $R_e=3.25$ Å (Table IV); SO-coupled six-curve structure and ³P₂ multiplet
    splittings (Eqs. 9, §III A); ion-mobility validation (Fig. 1, §III B). Heavy-
    alkali comparison via [GAH25] Table I scale. Open-shell gradual-evolution
    precedent: Pb⁺Heₙ (Slattery/González-Lezana snowball literature).
