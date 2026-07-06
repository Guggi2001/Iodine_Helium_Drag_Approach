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

**Revision 2026-06-21 — consistency-check pass (4 fixes); mechanism unchanged,
lock holds.** Cross-reading the loaded sources against the pinned ansatz surfaced
three scenario-independent bookkeeping fixes and one scenario-keyed budget split.
None touch the two-channel mechanism or any lock; all reversible.
1. **$n^*$ leak corrected: per-atom collective binding now uses the cation count
   21 throughout (was silently 20).** The figures "$|S|/n^*\approx124$ cm⁻¹
   $\approx179$ K" (K2 split, $\partial|S|/\partial n$ marginal; R3
   electrostriction-dominance) divide $|S|$ by $\approx20$ — the **neutral**
   I@He₂₀₀₀ first-shell count ([I2-notes] Fig. 9), not the **cation** $n^*=21$
   (Fig. 10). With $n^*=21$: $|S|/n^*=2484/21=\mathbf{118}$ cm⁻¹ $=170$ K, and the
   collective-vs-pair-at-radius ratio is $118/25=\mathbf{4.7\times}$ (was $5\times$).
   Qualitative conclusions stand (electrostriction still dominant; marginal still
   $\approx D_0(1)=107$). **OQ8 upgraded LOW→LOW–MEDIUM:** the $20/21$ ambiguity is
   no longer "geometry only" — it had leaked into the K2 marginal-release
   bookkeeping (A8) and the $E_\text{elec}$ magnitude. Adopt $n^*=21$ for all
   per-atom energetics; carry $n^*=21^{+0}_{-1}$ only where the GAH25 $R_e$-scaling
   ($\approx20$) is the explicit cross-check. (CALIBRATION_MAP rows 12b, 22 mirror.)
2. **$X_2$ Form U integrated band corrected to 0.25–0.28 eV; the 0.12 lower half
   was an unstated crowding/drag-binding cross-check, not a $\kappa$ extreme.**
   Evaluating Form U analytically with the cliff at $n^*+\tfrac12=21.5$,
   $\sum_{n=1}^{21}D_0=21D_\text{floor}+(D_0(1)-D_\text{floor})\,W(\kappa)/N(\kappa)$,
   $W=\sum_n(1-\sigma(n))$, $N=1-\sigma(1)$, over $\kappa\in[0.3,5]$ gives **$X_2$
   0.25–0.28 eV** and **mixture 0.17–0.19 eV**. The mixture span (~11%) **confirms
   the "nearly $\kappa$-independent" gate-threshold claim** (R9, S2, row 21) for
   Form U as defined. But the previously quoted $X_2$ band "0.12–0.28" imported a
   **crowding reduction** whose lower end coincides with the drag binding
   0.117 eV — i.e. the one effect that breaks $\kappa$-independence sat inside the
   band advertised as $\kappa$-independent. *Fix:* report Form U as **0.25–0.28
   ($X_2$) / 0.17–0.19 (mix)**; list the crowding-reduced value and the §6.5.1 drag
   binding (0.117 eV) as **separate cross-checks**, not band endpoints. The
   decoupling (gate/floor $\kappa$-independent → pinnable ahead of the Tier-2 fit)
   survives, now cleanly. (CALIBRATION_MAP row 21 mirror.)
3. **Third electronic sub-case added: cooling-relaxation, between the mixture and
   $X_2$ rungs.** The A10 default forms the mixture rung as the equal-weight mean
   of dissociation energies ($\tfrac13(106.9{+}62.4{+}54.0)=74.4$ cm⁻¹). But
   [IHe05] §III.B averages **transport cross-sections** over a *frozen* SO
   population; the dissociation-energy mean is this project's own ansatz, not
   IHe05's operation. Physically, an evaporating complex **cools**, so it may relax
   toward the deepest curve ($X_2$, 106.9) rather than stay in the birth mixture —
   the gate-relevant rung could sit **above 74.4, below 106.9**. *Fix:* the
   `ladder_electronic_picture` fork carries a **third option**, `cooling_relaxed`
   (rung between mixture and $X_2$, e.g. a relaxation-weighted blend), in addition
   to `statistical_mixture` (default) and `x2_only`. Strictly between the two
   existing brackets, so it cannot widen the Tier-2 search beyond current bounds;
   it fills the gap the binary fork skips. See A5, A10, §11.
4. **$E_\text{avail}^\text{ion}$ split into a scenario-keyed pair (validation vs
   production); A7 reworded to "robust, narrower margin."** The pinned 2.70 eV is
   the **production** onset (vertical double-ionization at $R_e(\mathrm I_2)=2.666$ Å,
   $\tfrac12\cdot14.40/2.666$). Tier 0 and Tier 1 instead validate against the
   [I2-notes] protocol that ionizes after the I atoms have drifted to **$d\approx9$
   Å** (A-state dissociation then double-ionization, §III.B.2; author-confirmed),
   giving $E_\text{avail}^\text{ion}=\tfrac12\cdot14.40/9=\mathbf{0.80\ eV}$ — a
   $3.375\times$ smaller budget. Since $\sum_i D_0$ is scenario-independent, the S2
   floor $f_\text{int}^\text{floor}=\sum_iD_0/E_\text{avail}^\text{ion}$ scales as
   $1/E_\text{avail}$: mixture **0.065–0.072 @ 2.70 eV → 0.21–0.24 @ 0.80 eV**;
   $X_2$ (Form U) **0.093–0.104 → 0.31–0.35**. **A7 reworded:** the self-unbound
   onset is **robust across both scenarios** (floor stays far below the hard cap
   $f_\text{int}\le1$), but the headroom shrinks from $\sim14\times$ (floor 0.07 @
   2.70 eV) to $\sim4.5\times$ (floor 0.22 @ 0.80 eV) — keep "robust," drop fixed-margin
   language. The **soft upper edge $f_\text{int}\lesssim0.2$ is a velocity-consistency
   plausibility bound, not a constraint** (S2): the floor is the only derived
   early-window edge. The heuristic is **scenario-invariant** —
   $E_\text{trans}/E_\text{avail}\approx40\%$ at both budgets ($v\propto\sqrt{E_\text{avail}}$:
   0.32 eV at $v\approx5.5$ Å/ps for 0.80 eV, [I2-notes] Fig. 23; $\sim1.1$ eV at
   $\sim10$ Å/ps for 2.70 eV), so the soft ceiling reads $\sim0.6$ at **both**;
   only the *floor* moves with scenario, narrowing the window from below. *Config:*
   `coulomb_available_eV` becomes scenario-keyed — **0.80 eV (validation, $d{=}9$ Å) /
   2.70 eV (production, $R_e$)** — stamped to the same scenario tag that guards the
   drag↔mass pairing (DESIGN §6.5/§6.5.1), so a 2.70-eV onset cannot run against
   0.80-eV-calibrated drag/shell references without tripping the guard. Add Fig. 23's
   peak $v\approx5$–6 Å/ps as an explicit R10 / early-window cross-check **at the
   0.80 eV budget**. **This A7/S2 partition is provisional pending OQ2** — the
   $KE_\text{shed}$/energy-partition from the GAH25 movies or Halberstadt would
   convert the soft ceiling into a real budget-derived headroom and may move the
   window's upper edge and the translational fraction. See A7, S2, §6.11, R10, OQ2;
   CALIBRATION_MAP rows 14, 15.

**Revision 2026-06-21 (cont.) — Method-level consistency proof: invariant closes;
RRK mode-count corrected.** Auditing the *mechanism* (not just parameters) as a
closed system gave one pass and one fix.
1. **Five-term invariant verified closed under every channel (PASS, no change).**
   Tracking $\Delta(E_\text{kin},E_\text{pot},E_\text{dissip},E_\text{mass\_transfer},
   E_\text{int})$ through drag (continuous), pickup S1, cold-shed K1, and K2 cooling:
   each channel sums to zero (drag: KE↔dissip; S1: $-D_0+f_\text{ret}D_0+(1{-}f_\text{ret})D_0=0$
   plus capture-KE defect → mass_transfer; K1: $E_\text{int}{-}D_0$, $E_\text{pot}{+}D_0$,
   A8 marginal → bath; K2: $dE_\text{int}=-dE_\text{dissip}$ at fixed $N$). K1∩K2
   non-overlap (cold-shed neutrality) and drag↔$E_\text{int}$ disjointness hold. The
   §6 invariant is a **genuine closed conservation law for the Method** (modulo
   Verlet drift), not approximate. No edit — recorded as the audit result.
2. **RRK effective DOF corrected $3n-6\to3n-3$; last atom is direct dissociation
   (MEDIUM).** Two coupled errors: (i) §4 defined $s$ as the DOF "of the I⁺Heₙ
   **complex**" ($n{+}1$ atoms → $3n-3$) but wrote $3n-6$ (an $n$-atom count),
   disagreeing with its own definition by 3; (ii) $3n-6$ goes $s{=}0$ at $n{=}2$ and
   $s{=}{-}3$ at $n{=}1$, where $k=\nu(1-D_0/E_\text{int})^{s-1}$ **diverges at
   threshold** ($s<1$), breaking the saturating-rate / no-avalanche guarantee (§4)
   and R5's physical-timescale claim — in the *small-$n$ regime that is the only
   place $\{\nu,s\}$ are observable* (at large $n$ the huge exponent makes $k$ a
   near-step function set by cooling, not kinetics; the size-distribution tail
   sensitivity lives at $n\lesssim$ few, which the total-strip OQ6 run also
   traverses). *Fix:* (a) $s=3n-3$ (full complex, $\ge3$ for $n\ge2$); (b) $n{=}1$
   modelled as **direct dissociation $k=\nu$** gated by $E_\text{int}>D_0(1)$
   (single mode → statistical RRK degenerate; also the boundary where cold-shed
   neutrality A8 is weakest — flagged); (c) **config-load guard $s\ge1$** on any
   effective-scalar override (dissipativity-style, §3.3-analog). Bounded
   $k\in[0,\nu)$ and the no-avalanche guarantee now hold end-to-end through the last
   atom. Touches §4, A11, §10/§11; CALIBRATION_MAP row 10. Mechanism and invariant
   unchanged.

**Revision 2026-06-21 (cont. 2) — pickup↔gate loop proven stable; occupancy cap
(Form B) added for completeness.** Continuing the Method-level consistency proof to
the feedback loop between accretion and the self-bound gate.
1. **Pickup↔gate loop is stable and well-posed — "equilibrium is emergent, not a
   parameter" now proven (§6.11 stability note).** On the mean-field flow
   $\langle\dot n\rangle,\langle\dot E_\text{int}\rangle$: (a) the self-unbound
   margin $G\equiv E_\text{int}-\Sigma(n)$ is a **pathwise Lyapunov function** —
   between pickups $\dot G=-E_\text{int}/\tau\le0$, and *at* a pickup
   $\Delta G=-(1-f_\text{ret})D_0(n{+}1)\le0$ because accretion adds a full binding
   rung but only $f_\text{ret}$ of it as heat — so $G$ monotonically falls and the
   gate **always** crosses in finite time (pickup is *stabilizing*; the earlier
   runaway worry had the sign backwards). $G(0)>0$ is exactly the $f_\text{int}$
   floor (S2). (b) Terminal $n$ is a **stable freeze-out attractor** governed by the
   dimensionless $\Pi(n)\equiv\lambda(n)f_\text{ret}\tau$: $\Pi>1$ shedding persists,
   $\Pi<1$ freeze; with $\Pi\to0$ at exit ($\rho_\text{He}\to0$) guaranteeing
   termination on every trajectory. $\Pi$ is the quantitative spine of the R1
   regime axis; no limit cycle (K2 dissipative, source decays). No edit to the
   mechanism — recorded as the audit result, with $\Pi$ and $t_\times$ as the two
   early diagnostics.
2. **Occupancy cap added (Langmuir "Form B"), A12 — closes the only gap the proof
   surfaced.** Density-only $\lambda_\text{attach}$ had no $n$-dependence, so a
   resting/slow ion would accrete unboundedly ($n\to\infty$), contradicting
   [Nat23]'s resting-ion leveling-off. Fix: $\lambda_\text{attach}=\lambda_0
   (\rho_\text{He}/\rho_\text{bulk})(1-n/n^*)_+^{\,p}$ — a site-saturation factor
   →0 at $n^*$. **Caps the rate only**, so the §6 invariant and the $G$-crossing are
   untouched and the freeze-out attractor *gains* a second stabilizing route
   ($n\to n^*$ as well as exit). **Inert for production** (the ion exits before
   saturation; Form B ≡ density-only for ejection) — added for resting-ion
   correctness, reviewer-defensibility, and the Tier-1-tail case. Costs the exponent
   $p$, **default tied to $\kappa$** *(NB 2026-07-01: an **inverse** tie — rigid
   shell = large $\kappa$ = sharper cutoff = **smaller** $p$; **not** a literal
   $p=\kappa$. Phase B holds $p=1$ fixed; see the §11 config NB / A12 / Tier-2 §3.2.)*
   (both encode first-shell abruptness — not
   independent, like $s\!\leftrightarrow\!\kappa$), so **zero net new free
   parameters** unless the size distribution forces the split. Rejected: hard wall
   (discontinuous; makes the $\pm1$–2-uncertain $n^*$ a hard input) and shell-2
   two-reservoir ("Form C", deferred). Touches §4 (rate + dim table), §6.11, A12
   (new), §10/§11; CALIBRATION_MAP rows $\lambda_\text{attach}$, $p$ (new), $\Pi$
   (new).

**Revision 2026-06-21 (cont. 3) — integrator↔mass-jump operator split specified;
Method proof closes (conditional).** The last structural item: where the
continuous BAOAB integrator meets the discrete mass jumps (A13, new).
1. **SQ1 (built, accepted): drag-on path is $O(dt)$, not $O(dt^2)$.** The
   state-dependent $\gamma(v)=gbv^2/m$ is frozen at $v_\text{in}$ and applied as
   $e^{-\gamma dt}$ — locally first-order on the cubic. **Retire the BAOAB
   second-order claim for production** (holds only in the constant-$\gamma$ limit).
   The trade is correct: freezing $v_\text{in}$ buys **exact dissipation
   bookkeeping** ($\Delta E_\text{dissip}=\tfrac12 m(\|v_\text{in}\|^2-\|v_\text{out}\|^2)$,
   any $dt$) and **unconditional dissipativity** ($\|v_\text{out}\|\le\|v_\text{in}\|$,
   no large-$dt$ blow-up on the stiff cubic) — both worth more than an order here.
   Residual: a **one-signed over-braking bias** (folded into R10; expect simulated
   peak $v$ slightly below Fig. 23 at fixed $dt$).
2. **SQ2 (spec, unbuilt): momentum-conserving velocity reset is an invariant
   precondition.** Jumps reset $v^+=(m v^-\pm m_\text{He}u_\text{He})/(m\pm m_\text{He})$
   (He at rest). **The §6 five-term invariant closure presupposes this** —
   $E_\text{mass\_transfer}$ *is* the reduced-mass KE defect
   $\tfrac12\tfrac{m\,m_\text{He}}{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2$ this reset
   produces; a label-only $v$ voids the proof. This also **corrected the invariant
   text** from the heavy-ion $\tfrac12 m_\text{He}v^2$ (a ~3% closure error at
   $n{=}1$) to the exact reduced-mass form. Jump-step order reduction is benign
   (jump-steps are $dt$-independent in number → measure-zero as $dt\to0$); ordering
   fixed jump-then-O; at most one mass event per step (shed before pickup).
3. **SQ3 (spec, unbuilt): post-jump O-step uses $m^+$** in friction and the FDT
   noise amplitude (forced by SQ2; small but a definiteness requirement).
   *Net:* no structural defect. **The Method consistency proof now closes,
   conditional on SQ2(a) being implemented as specified** — recorded as an
   implementation precondition, not an open physics question. Touches §4 (jump
   rules + dim table), §6 (invariant, reduced-mass correction), A13 (new), §11
   (4 config fields); CALIBRATION_MAP $E_\text{mass\_transfer}$ note.

**Revision 2026-06-17 (folded in) — EPAPS fit closes the first rung and the RRK
prefactor; A10 mixture rung now numeric.** The [IHe05] EPAPS analytic fit
parameters (corrected Eq. (3), erratum [IHe05-E]) were obtained, so $V''(R_e)$ for
the He–I⁺ curves is now exact (the fit reproduces Table IV to four figures;
$R_e=3.2527$ Å, $D_e=143.89$ cm⁻¹; the Degli-Esposti–Werner switching function is
inert at the well, $f(R_e)=1$). Three quantities move from order-of-magnitude
prior to **pinned**, and the A10 fork acquires concrete first-rung numbers:
1. **First rung $D_0^{\,\mathrm{I^+}}(1)$, $X_2$/³Π = 106.9 cm⁻¹** (13.3 meV,
   0.01325 eV). Computed from the **exact $J{=}0$ vibrational ground state** of the
   fitted curve (radial Schrödinger, $\mu=3.880$ amu), **not** a harmonic ZPE: the
   well holds 5 bound levels ($-106.9,-52.3,-20.9,-6.2,-1.1$ cm⁻¹), and the true
   ZPE $G(0)=37.0$ cm⁻¹ is **26% of $D_e$** — far larger than the ~10% Na⁺-analogy
   fraction the old $D_0\approx125$–135 cm⁻¹ estimate assumed (R3). Residual
   uncertainty is now the [IHe05] ±3% well-depth accuracy, $\approx\pm3$ cm⁻¹; the
   ZPE-method error is eliminated. See R3.
2. **RRK prefactor $\nu=2.42$ ps⁻¹** ($\omega_e=80.6$ cm⁻¹, well-bottom curvature
   $V''(R_e)=748.1$ cm⁻¹/Å²) — same order as the old $\sim\mathcal{O}(1)$ prior but
   now sourced from the fit, not estimated. *Caveat for the $s$/`ladder_steepness` ($\kappa$)
   joint calibration (A11):* this is the bottom-of-well frequency; the
   near-threshold level spacing collapses (55→31→15→5 cm⁻¹ up the anharmonic
   ladder), so the effective attempt frequency of a near-dissociation complex is
   lower and is absorbed into the effective $s$. See §4, A11.
3. **A10 electronic-picture fork is now numeric.** Via [IHe05] Eq. (9) with the
   atomic ³P$_j$ splittings ($D_0^{at}=6451$, $D_1^{at}=7090$ cm⁻¹), the three
   ground-correlating SO dimer rungs are $D_0$: $X_2=106.9$, $I_1=62.4$,
   $I_0=54.0$ cm⁻¹, so the equal-weight **statistical-mixture first rung $=74.4$
   cm⁻¹** (9.23 meV, 0.00923 eV) — **~70% of the $X_2$-only rung**. This quantifies
   the OQ1 double-count hazard: a ~30% reduction at the dimer level is precisely
   what must *not* be stacked on top of an already dynamically-lowered drag
   binding. See A10, §10A OQ1. *Three objects kept distinct (do not equate):*
   $D_0(1)$ (single dimer bond, above) ≠ $\sum_i D_0(i)$ (integrated pair ladder =
   gate threshold) ≠ $S_{\mathrm{I^+}}=-0.308$ eV (many-body first-shell solvation,
   pair sum + He–He + electrostriction) ≠ $E_\text{bind}=0.1168$ eV (ion↔droplet
   effective binding from drag).
The $n>1$ ladder *shape* (R3/A5) remains open; only the first rung and the
electronic-picture rung scale are pinned here. Added reference [IHe05-E].

**Revision 2026-06-17 (cont.) — Form U ladder adopted; $|S|$ is collective; the
shell cliff is geometric.** The $n>1$ ladder shape is now fixed *as a form* (one
continuous knob), and two findings reshape R3/A5:
1. **Form U adopted.** The discrete `ladder_shape ∈ {gradual, shell_structured}`
   is retired in favor of a single sigmoid family,
   $D_0(n)=D_\text{floor}+(D_0(1)-D_\text{floor})(1-\sigma(n))/(1-\sigma(1))$,
   $\sigma(n)=[1+e^{-\kappa(n-n^*-\tfrac12)}]^{-1}$, anchored at the pinned rung, a sourced
   bulk-He floor $D_\text{floor}=|\mu_\text{He}^\text{bulk}|\approx4.97$ cm⁻¹, and
   sourced $n^*\approx21$. The **single Free shape knob is $\kappa$** (gradual
   $\kappa\!\to\!0$ ↔ cliff $\kappa\!\gg\!1$), Tier-2-arbitrated jointly with the
   electronic picture and $\{\nu,s\}$. New config field `ladder_steepness`;
   tabulated ladder kept as declared fallback. See R3, A5, §11.
2. **$|S_{\mathrm{I^+}}|$ is collective, not $\sum_i D_0$ (corrects the old
   integrated-ladder cross-check).** A monotone pair ladder anchored at the pinned
   $D_0(1)=106.9$ cm⁻¹ **cannot** reach the DFT $|S|=2484$ cm⁻¹ (flat ceiling
   2245; shortfall ≥239 cm⁻¹, ≥922 for the mixture). $|S|$ therefore carries
   collective content (electrostriction + DFT correlation) and is an **upper
   bound**, not a rung-sum target — do not calibrate rungs to $|S|/n^*$. The
   reachable integrated cross-check is the §6.5.1 drag binding (942 cm⁻¹). New
   author-contact item: the $|S|$ energy reference (§10A).
3. **The radial cliff is geometric, surviving open-shell blurring.** Pure
   charge-induced-dipole binding drops **7.5×** from shell-1 to shell-2 (radial,
   not electronic), and the shell sits at the [GAH25] radius $r_1^e\approx4.67$ Å
   (not the pair $R_e$), where He–He is roomy (NN 3.97 Å > 2.97 Å, mildly
   attractive), so the physical prior is **mild in-shell decline then cliff**
   (large $\kappa$), revising the earlier "open-shell ⇒ gradual/shell-less" lean.
   See R3. *(Geometry corrected 2026-06-17 on re-reading GAH25 — an earlier note
   using the pair $R_e$ wrongly found a compressed 2.70 Å shell and a 13.4× cliff.)*
All three are **reversible** if OQ1 or the $|S|$-reference question reopens.

**Revision 2026-06-17 (cont. 2) — early-window scalars pinned; gate/floor
decouple from $\kappa$.** With Form U in place, three coupled early-window
quantities are now numeric:
1. **$E_\text{avail}^\text{ion}=2.70$ eV (pinned, per-ion convention adopted).**
   Pair release $e^2/R_{\mathrm{II}}=14.40/2.666=5.40$ eV at $R_{\mathrm{II}}=
   R_e(\mathrm{I_2})=2.666$ Å, split equally → 2.70 eV/ion. Switched S2 from the
   pair value (where $f_\text{int}$ silently carried the ½) to the per-ion budget;
   *reversible* (pair value doubles the floor). See S2, §11.
2. **Integrated first-shell ladder $\sum_{i=1}^{21}D_0$ now numeric:** $X_2$
   0.12–0.28 eV (flat-shell↔crowding-reduced), mixture 0.17–0.19 eV. **Nearly
   $\kappa$-independent** (~10%) once the cliff is centered at $n^*+\tfrac12$ — so
   $\kappa$ shapes only the terminal-$n$ histogram, **not** the self-bound gate
   threshold (row 21, R3, S2).
3. **$f_\text{int}$ self-unbound floor $\approx0.04$–$0.10$ ($X_2$, per-ion),
   ~0.06–0.07 (mixture).** Small ⇒ **the GAH25 self-unbound onset is robust, not
   fine-tuned** (strengthens A7/A8). Floor is picture-set and pinnable now, ahead
   of the $\kappa$/picture Tier-2 fit.
Also adopted: the **Form U $n^*+\tfrac12$ centering** refinement (cliff between
shell-1 and shell-2). All reversible.

**Revision 2026-06-17 (cont. 3) — GAH25 full paper studied; $E_\infty$ split
locked; ladder geometry corrected.** Re-reading the actual [GAH25] (not the prior
summary) confirmed the K2 machinery verbatim (Eqs. 11–12, Table III) and forced:
1. **$E_\infty$ binding split LOCKED (decision: split).** $E_\text{solv.struct}=
   E_\text{bind}^\text{pair}(N)+E_\text{elec}(N)+E_\text{int}(N)$, cooling to the
   **occupancy-resolved** $E_\infty(N)=-|S(N)|$. GAH25 geometry shows
   electrostriction is the *dominant* binding term (collective $\approx124$ vs
   pair-at-shell-radius $\approx25$ cm⁻¹/atom), so it is carried explicitly. Resolves
   **OQ6** (occupancy-resolution removes the stripping cap — a fixed $E_\infty$
   drives $E_\text{int}^\text{eq}$ to $-0.28$ eV by $N{=}2$ and halts shedding) and
   the equilibrium layer of **R12** ($E_\text{int}^\text{eq}=0$ exactly). A8 amended
   (marginal electrostriction release → bath). See K2, R12, A8, A9, OQ6.
2. **Ladder geometry corrected (was wrong in cont./cont.2).** The shell sits at
   the GAH25 radius $r_1^e\approx4.67$ Å, **not** the pair $R_e=3.25$ Å: radial
   cliff **7.5×** (not 13.4×); shell-1 He–He **roomy/mildly attractive** (NN
   3.97 Å > 2.97), **not** compressed. Form U and the large-$\kappa$ prior survive;
   only the two numbers and the crowding sign change. See R3.
3. **Secondary GAH25 cross-checks (logged):** $n^*\approx20$ by $R_e$-scaling (vs
   [I2-notes] 21; OQ8); $\lambda_\text{attach}$ for I⁺ should center ~0.7–1.1/ps
   (Rb⁺/Cs⁺) not 2.0 (Na⁺); Calvo K⁺ PIMC as a possible $\kappa$ anchor. See OQ8,
   R1, references.
All reversible pending OQ1/OQ6/OQ7/OQ8.

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
   value $\nu=2.42$ ps⁻¹ from the [IHe05] $X_2$ curvature (**pinned 2026-06-17**;
   was $\sim\mathcal{O}(1)$ prior); $s$ is
   mode-counted $3n-6$ by default (effective-scalar override available),
   calibrated *jointly* with `ladder_steepness` ($\kappa$) since both probe shell rigidity.
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
  the solvation structure. **Full paper studied 2026-06-17:** Table I (pair
  $D_e/R_e$ — places I⁺ at Rb⁺), Table II ($r_1^e,r_2^e,n_1^e$ — fixes the shell
  *geometry*: $r_1^e\gg$ pair $R_e$, used to correct the R3 cliff to 7.5× and the
  crowding sign), Table III (Newton fit $t_0,\tau,E_\infty$ — confirms K2 exactly).
  Energy analysis is **Na⁺-only**; cools a *growing* shell (vs our shrinking).
- **[Calvo24/25]** — F. Calvo, J. Chem. Phys. **161**, 121101 (2024) (Na⁺) and
  J. Low Temp. Phys. **51**, 453 (2025) (K⁺), PIMC/RPMD pump–probe. **Lead (not
  yet obtained):** K⁺ is the nearest alkali with a potentially published
  per-rung evaporation-energy ladder $D_0(n)$ — a candidate to anchor the Form U
  steepness $\kappa$ rather than leaving it fully free (Tier 2). Also the
  total-stripping (Calvo) limit referenced for the secondary-run target.
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
  **EPAPS analytic fit obtained (2026-06-17):** doc E-JCPSA6-122-018521 supplies
  the $\{g_l,\alpha,\beta,\delta,D_4,D_6,D_8\}$ for each curve, evaluated with the
  corrected short-range form ([IHe05-E]); this pins $V''(R_e)$ and hence
  $D_0(1)$ and $\nu$ exactly (2026-06-17 revision).
- **[IHe05-E]** — Buchachenko *et al.*, *Erratum*, J. Chem. Phys. **161**, 149901
  (2024), doi:10.1063/5.0237596. Corrects Eq. (3) of [IHe05] to
  $V_\text{SR}(R)=\sum_{l=0}^{8} g_l R^l \exp[-\alpha R-\beta]$ for consistency with
  the EPAPS fit parameters (flagged by N. Halberstadt — the author contact for
  OQ1/OQ4). "Results and conclusions unaffected"; Table IV stands. Required to
  evaluate the EPAPS fit correctly.
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

- **Pickup** — a **Poisson process**, rate $\lambda_\text{attach}(\rho_\text{He}(\text{depth}),n)$ (density-gated, with a Langmuir occupancy cap toward $n^*$; A12), each event $n \to n+1$, $m \to m + m_\text{He}$, with a **momentum-conserving velocity reset** $v^+=\dfrac{m\,v^-+m_\text{He}u_\text{He}}{m+m_\text{He}}$ (captured He taken **at rest** in the droplet frame, $u_\text{He}=0$, so $v^+=\tfrac{m}{m+m_\text{He}}v^-$; thermal $u_\text{He}$ deferred). This reset is **mandatory** — it is the operation whose KE defect *defines* $E_\text{mass\_transfer}$ in the invariant (A13); carrying $v$ through unchanged would break closure.
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
    energy, $u_\text{He}\approx0$), with the **symmetric momentum-conserving
    reset** $v^+=\tfrac{m\,v^- - m_\text{He}u_\text{He}}{m-m_\text{He}}\to
    \tfrac{m}{m-m_\text{He}}v^-$ (the lighter complex retains the momentum;
    A13). The bounded rate caps per-step shedding — no instantaneous
    cascade. At most **one mass event per step** (pickup or shed; if both
    Bernoulli draws fire, shed is applied first, then pickup — a fixed order
    for unambiguous bookkeeping, joint prob. $\sim\lambda\nu\,dt^2\sim4\times10^{-4}$; A13).

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

- The numeric $D_0^{\mathrm{Na^+}}(N)$ values (Na⁺–He specific). The I⁺ ladder
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
\lambda_\text{attach} = \lambda_0 \,\frac{\rho_\text{He}(\text{depth})}{\rho_\text{bulk}}\,
\Big(1-\tfrac{n}{n^*}\Big)_+^{\,p}.
$$
The final factor is a **Langmuir-style occupancy (site-saturation) cap** (added
2026-06-21, A12): a sticking coefficient that falls as the first shell fills,
$(\cdot)_+\equiv\max(\cdot,0)$, exponent $p\ge0$. It reduces to the density-only
form for $n\ll n^*$ (preserves the [Nat23] resting-ion anchor 2.0/ps, §2.7/§5) and
goes continuously to zero at $n{=}n^*$, so a resting/slow-exit ion saturates at
$n^*$ instead of accreting unboundedly. **For the production ejection problem it is
physically inert** — the ion exits ($\rho_\text{He}\to0$) before the shell
saturates, so $\Pi$ crosses 1 via density first (§6.11); it earns its keep only if
the Tier-1 tail lingers near saturation, and for resting-ion correctness (A12).
Units: the factor is dimensionless, $[\lambda_\text{attach}]=\text{ps}^{-1}$ ✓.
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
degrees of freedom of the I⁺Heₙ complex. **Parameter sourcing (updated 2026-06-17):**
- $\nu$ is the **I⁺–He stretch attempt frequency**,
  $\nu=\tfrac{1}{2\pi}\sqrt{V''(R_e)/\mu}$, with the **true** reduced mass
  $\mu=m_\text{He}m_{\mathrm{I^+}}/(m_\text{He}+m_{\mathrm{I^+}})=3.880$ amu (not
  $\approx m_\text{He}$; ~3% lighter, $+1.5\%$ on $\omega$). The [IHe05] EPAPS fit
  gives $V''(R_e)=748.1$ cm⁻¹/Å² exactly, so $\omega_e=80.6$ cm⁻¹ and
  **$\nu=2.42$ ps⁻¹ (pinned)** — same order as the earlier $\sim\mathcal{O}(1)$
  prior, now sourced from the fit rather than estimated. *Caveat:* this is the
  **bottom-of-well** frequency; the anharmonic level spacing collapses toward
  threshold (55→31→15→5 cm⁻¹ up the $X_2$ ladder), so a near-dissociation
  complex's effective attempt frequency is lower — absorbed into the effective
  $s$ (A11), reinforcing the joint $\{\nu,s\}$/`ladder_steepness` ($\kappa$) calibration.
- $s$ is **mode-counted, not free:** $s=3(n{+}1)-6=\mathbf{3n-3}$ — the
  vibrational DOF of the **full I⁺Heₙ complex** ($n$ He $+$ 1 ion $=n{+}1$ atoms,
  nonlinear), **corrected 2026-06-21** from the earlier $3n-6$, which counted an
  $n$-atom object and disagreed with this line's own "complex" definition by 3 and
  went $\le0$ for $n\le2$ (see the 2026-06-21 (cont.) revision and A11). So $s$ is
  $n$-dependent ($\sim60$ at $n{\sim}21$, $\sim12$ at $n{\sim}5$, $=3$ at $n{=}2$)
  and parameter-free by default, with a single effective-scalar override available
  for sensitivity (derive-by-default, override-for-test), **guarded $s\ge1$ at
  config-load** (a dissipativity-style bound, §3.3-analog: $s<1$ makes $k$ diverge
  at threshold). **Physical band $[3n-3,\,3n]$ (A11):** $3n-3$ is the conservative
  free-complex default; hindered shell rotations retained as **librations** push
  the effective $s$ *up* toward $3n$, while a blurred/separable shell (A10/R3)
  pushes it *down* (fewer effectively-coupled modes) — opposite drifts, both the
  same shell-rigidity question. So $s$ and `ladder_steepness` ($\kappa$) are
  calibrated *jointly* against the size distribution, not separately. $s$ also
  absorbs the classical-RRK-vs-RRKM simplification (A11).
  > **NB (2026-07-06 — probe resolution + promotion; supersedes the sweep-band
  > reading above).** The pre-F5 staircase probe fired this empirically: at
  > $s=3n-3$ (=60 at $n{=}21$) the cascade is kinetically frozen at
  > $n\approx20$ everywhere in the κ×picture×τ bands, while a **constant
  > $s_\text{eff}=8$** lands the anchored 21→19→14 staircase in magnitude and
  > timing, picture-robustly. **User decision (2026-07-06): $s$ is now a
  > Bounded constant-$s_\text{eff}$ knob** (band ≈[5, 20], landing [8, 12];
  > CALIBRATION_MAP row 10), with the mode count $s=3n-3$ demoted to the
  > classical-limit arm of the dof-convention selection. The $[3n{-}3,3n]$
  > band survives as the rigid-classical limit, not the sweep band — the
  > effective reservoir is ~an order of magnitude smaller (quantum
  > mode-freezing / weak coupling; the RRKM-direction difference this very
  > bullet declared $s$ absorbs). The $s\ge1$ guard, the $n{=}1$ direct
  > dissociation rule, and the mechanism/invariant are unchanged. Records:
  > `docs/drag_port/Tier2/drag_migration_log_tier2.md` ("s_eff mini-probe
  > EXECUTED" + the promotion decision); A11 resolution NB below.
- **Last-atom step ($n{=}1$) is direct dissociation, not statistical RRK
  (2026-06-21).** I⁺He₁ is a diatomic with a *single* vibrational mode, so the
  statistical phase-space picture is degenerate ($s=3n-3=0$ at $n{=}1$; even
  $3n-6$ failed here). The $n{=}1\to0$ loss is therefore modelled as **direct
  dissociation at the fixed attempt frequency**, $k=\nu$, gated only by
  $E_\text{int}>D_0(1)$ (no RRK bracket). This keeps the per-step rate bounded by
  $\nu\,dt$ (avalanche guarantee intact) through the last atom and is exactly the
  regime where cold-shed neutrality (A8) is weakest, so it is flagged as the
  **boundary of validity of the statistical evaporation picture** (A11).

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
| Occupancy cap factor | $(1-n/n^*)_+^{\,p}$ | dimensionless | $n,n^*$ counts; $p$ pure ✓ |
| Pickup KE defect (reduced-mass) | $\tfrac12\tfrac{m\,m_\text{He}}{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2$ | $\text{amu·Å}^2/\text{ps}^2$ = energy | exact under momentum reset (A13) ✓ |
| Dissociation rung | $D_0^{\,\mathrm{I^+}}(n)$ | eV (energy) | ✓ |
| Integrated binding (gate) | $\sum_i D_0^{\,\mathrm{I^+}}(i)$ | eV | compared to $E_\text{int}$ (eV) ✓ |
| RRK prefactor | $\nu$ | $\text{ps}^{-1}$ | $[\nu\,dt]=1$ ✓ |
| RRK bracket | $1 - D_0(n)/E_\text{int}$ | dimensionless | eV/eV ✓ |
| RRK effective DOF | $s=3n-3$ ($n\ge2$) | dimensionless | mode-counted, full complex (A11); $s{-}1\ge1$; guarded $s\ge1$ ✓ |
| RRK rate | $k=\nu(\cdot)^{s-1}$ ($n\ge2$) | $\text{ps}^{-1}$ | $[k\,dt]=1$; bounded $k\in[0,\nu)$ ✓ |
| Last-atom rate | $k=\nu$ ($n{=}1$, direct) | $\text{ps}^{-1}$ | gated $E_\text{int}>D_0(1)$; bounded by $\nu$ ✓ |
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
  E_\text{int}(0) = f_\text{int}\,E_\text{avail}^\text{ion},
  \qquad f_\text{int}\in[0,1],
  $$
  with **$E_\text{avail}^\text{ion}$ the per-ion share** of the I–I Coulomb release
  (**scenario-keyed, 2026-06-21**) and $f_\text{int}$ the small fraction coupling
  into *this* ion's shell-*internal* modes. *Source:* the pair release is
  $E_\text{avail}^\text{pair}=e^2/R_{\mathrm{II}}=14.40/R_{\mathrm{II}}$, split
  equally by equal-mass dissociation, so $E_\text{avail}^\text{ion}=7.20/R_{\mathrm{II}}$
  eV ($[\text{eV·Å}/\text{Å}]=\text{eV}$ ✓). **Two scenarios at different $R_{\mathrm{II}}$:**
  - *Validation (Tier 0/1):* the [I2-notes] protocol dissociates I₂ on the A state
    and double-ionizes after the atoms drift to $R_{\mathrm{II}}\approx9$ Å →
    $E_\text{avail}^\text{ion}=7.20/9=\mathbf{0.80\ eV}$. **This is the budget the
    drag and the 21→19→14 shell references were generated under**, so it is the one
    the mass model uses when reproducing them (author-confirmed, §III.B.2).
  - *Production:* vertical double-ionization at $R_{\mathrm{II}}=R_e(\mathrm{I_2})
    =2.666$ Å → $E_\text{avail}^\text{ion}=7.20/2.666=\mathbf{2.70\ eV}$ (pair 5.40),
    a $3.375\times$ hotter onset.
  *Convention flag (reversible):* the doc previously used the **pair** value with
  $f_\text{int}$ silently absorbing the ½; the per-ion value is cleaner and is the
  default — switching back doubles the floor below. The fraction is preferred over a
  bare $E_\text{int}(0)$ because it is physically interpretable and intrinsically
  bounded. **Two-sided bounds, scenario-keyed (2026-06-21):**
  - *Lower (self-unbound floor):* $E_\text{int}(0) > \sum_{i=1}^{n_0}
    D_0^{\,\mathrm{I^+}}(i)$ at $n_0=21$, i.e. $f_\text{int} > f_\text{int}^\text{floor}
    =\sum_i D_0/E_\text{avail}^\text{ion}$ ($[\text{eV}/\text{eV}]$ dimensionless ✓).
    The integrated first-shell ladder from Form U (pure $\kappa$-range, corrected
    2026-06-21): **$X_2$ 0.25–0.28 eV, mixture 0.17–0.19 eV** (both ~11% over
    $\kappa\in[0.3,5]$ — confirms near-$\kappa$-independence). *Separate
    cross-checks (not band endpoints):* the crowding-reduced value and the §6.5.1
    drag binding 0.117 eV; the earlier "0.12–0.28" $X_2$ band conflated the
    crowding reduction into the band. **Scenario-keyed floor (2026-06-21):** at
    $E_\text{avail}^\text{ion}=2.70$ eV (production), $f_\text{int}^\text{floor}
    \approx0.09$–$0.10$ ($X_2$), ~0.065 (mix); at $E_\text{avail}^\text{ion}=0.80$
    eV ($d{=}9$ Å, the **Tier-0/1 validation** budget), $\approx0.31$–$0.35$
    ($X_2$), $0.21$–$0.24$ (mix). **Self-unbound onset is robust across both
    scenarios** (floor far below the hard cap $f_\text{int}\le1$), but the headroom
    shrinks from $\sim14\times$ (2.70 eV) to $\sim4.5\times$ (0.80 eV); the
    mechanism does not hinge on a tuned partition (A7, A8), though the Tier-1 margin
    is the narrower one.
  - *Upper:* hard cap $f_\text{int}\le 1$. The ~0.2 edge is a **velocity-consistency
    plausibility bound, not a constraint** — Tier 2 may exceed it with a flag, not a
    rejection. Velocity sanity check, scenario-invariant in *fraction*
    ($E_\text{trans}/E_\text{avail}\approx40\%$, since $v\propto\sqrt{E_\text{avail}}$):
    ejection ~10 Å/ps → ~1.1 eV of the 2.70 eV (production); ~5–6 Å/ps → ~0.32 eV of
    the 0.80 eV ([I2-notes] Fig. 23, validation) — both read a soft ceiling ~0.6.
    Only the *floor* moves with scenario; the ceiling does not. **Provisional pending
    OQ2** (the $KE_\text{shed}$/partition would convert the soft ceiling into a real
    budget-derived headroom). Working window: firm lower edge, soft ceiling.
  - *Decoupling (load-bearing for build order):* with the Form U cliff at the
    shell boundary $n^*+\tfrac12$, the first-shell sum is **nearly $\kappa$-independent**
    (varies ~10% over $\kappa\in[0.3,5]$) — $\kappa$ shapes only the cross-cliff
    region (terminal-$n$ histogram), **not** the gate threshold or this floor. So
    the floor is **picture-set and pinnable now**, without waiting on the $\kappa$
    arbitration (R9, §6.11).
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
  = -\frac{E_\text{solv.struct}-E_\infty(N)}{\tau_\text{dissip}},
  \qquad \tau_\text{dissip}=\tau_\text{GAH25}.
  $$

  **Binding split — pair + electrostriction (LOCKED 2026-06-17, GAH25-vindicated;
  resolves OQ6 and R12).** The equilibrium binding is *not* the pair ladder sum:
  re-reading GAH25 (Table I/II), the first shell sits at $r_1^e\approx4.67$ Å (far
  outside the pair $R_e=3.25$ Å), where pair polarization is only
  $D_4/r_1^{e4}\approx25$ cm⁻¹/atom, yet the DFT per-atom shell binding is
  $|S|/n^*\approx170$ K $\approx118$ cm⁻¹/atom ($n^*=21$, corrected 2026-06-21;
  was 179 K/124 cm⁻¹ off $n^*{=}20$) — so **collective electrostriction
  (snowball compression) is the dominant binding term, ~4.7× the pair-at-radius**,
  not a small correction. We therefore carry it explicitly:
  $$
  E_\text{solv.struct}(N)=\underbrace{-\textstyle\sum_{i\le N}D_0^{\,\mathrm{I^+}}(i)}_{E_\text{bind}^\text{pair}(N)\ \text{[IHe05]}}
  +\underbrace{-\big(|S(N)|-\textstyle\sum_{i\le N}D_0\big)}_{E_\text{elec}(N)\le0\ \text{[I2-notes]/GAH25}}
  +\;E_\text{int}(N),
  \qquad E_\infty(N)=-|S(N)|.
  $$
  *Dim:* every term in eV; $E_\infty(N)$ in eV ✓. Two locked consequences:
  - **$E_\infty$ is occupancy-resolved (OQ6 resolved).** A *fixed* full-shell
    $E_\infty$ would drive the reconstructed $E_\text{int}^\text{eq}=E_\infty-E_\text{bind}(N)$
    increasingly negative as the shell strips ($-0.03$ eV at $N{=}21$ →
    $-0.28$ eV at $N{=}2$), pushing $E_\text{int}$ below every rung and
    **mechanically halting shedding** — a hard cap. With $E_\infty(N)=-|S(N)|$
    tracking the current shell ($\to0$ as $N\to0$), full stripping stays
    reachable. $|S(N)|=|S|\cdot\sum_{i\le N}D_0/\sum_{i\le n^*}D_0$ distributes the
    collective full-shell value over the ladder shape (*flagged assumption*; the
    GAH25 near-constant collective marginal $\partial|S|/\partial n\approx118$ cm⁻¹
    supports a near-linear distribution).
  - **$E_\text{int}$ reconstruction is now exact (R12 systematic eliminated).**
    Recovering $E_\text{int}=E_\text{solv.struct}-E_\text{bind}^\text{pair}(N)-E_\text{elec}(N)$
    uses the *same* collective binding K2 cools toward, so $E_\text{int}^\text{eq}=0$
    with no spurious offset — the $|S|-\sum D_0$ "collective excess" lives in
    $E_\text{elec}$, not dumped into $E_\text{int}$ (cf. the old R12 hazard).
    Still **valid only after $t_\times$** (A9).

  The single relaxation time $\tau_\text{dissip}$ is now applied to the *correct*
  quantity, removing the quantity-mismatch error of cooling $E_\text{int}$
  alone. $\tau$ is **transplanted directly** ($\tau_\text{K2}=\tau_\text{GAH25}$,
  no inflation) — see the no-double-count finding below. It remains bracketed
  but not pinned (§10, R8): experiment gives $\tau\approx2.6\pm0.4$ ps (Na⁺,
  $\langle N\rangle{=}3600$); TDDFT gives $\tau\approx7.3$ (shell-1) – $16.5$
  (shell-2) ps — treat as a ±factor-3 sweep band, with the residual that these
  are **Na⁺** numbers transplanted to a **Rb⁺-like** I⁺ (§6.11, R8), *and* that
  GAH25 cools a shell *growing* toward equilibrium whereas the I⁺ shell *shrinks*
  (a direction-of-evolution transplant — rate assumed symmetric).

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
`E_mass_transfer_eV` (the renamed §2.9 field) absorbs the pickup KE defect and
the cold-shed bookkeeping. **Exact (reduced-mass) form (corrected 2026-06-21,
A13):** the momentum-conserving capture defect is
$$
\Delta E_\text{cap}=\tfrac12\,\frac{m\,m_\text{He}}{m+m_\text{He}}\,\|v^- - u_\text{He}\|^2
\;\xrightarrow{u_\text{He}=0}\;\tfrac12\,\frac{m\,m_\text{He}}{m+m_\text{He}}\,\|v^-\|^2,
$$
$[\text{amu·Å}^2/\text{ps}^2]$ ✓ — **not** the heavy-ion approximation
$\tfrac12 m_\text{He}v^2$ (baseline §7.1.3), which is the $m\gg m_\text{He}$ limit
and overstates the defect by $\sim3\%$ at $n{=}1$ (reduced mass 3.88 vs 4 amu),
shrinking as the complex grows. Using the reduced-mass form makes the invariant
close **exactly** under the §4 momentum reset, not to $\sim3\%$. $E_\text{int}$
remains the tracked reservoir in the invariant; K2 now cools $E_\text{solv.struct}$,
but at fixed $N$ (between shed events) $E_\text{bind}(N)$ is constant so
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
already in the v6 schema (§7; **v6→v7** in the live build — see the §7 NB), so
it carries **zero schema cost** beyond the `E_int_eV` field itself. Cross-check
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
do not pre-judge it. *Reachability (OQ6, RESOLVED 2026-06-17):* whether secondary
runs can *reach* full strip was contingent on the K2 asymptote — and a *fixed*
full-shell $E_\infty$ would indeed mechanically leave a residual shell (it drives
the reconstructed $E_\text{int}^\text{eq}$ negative as the shell shrinks, halting
shedding). The **occupancy-resolved** $E_\infty(N)=-|S(N)|$ now adopted (K2 split)
removes this cap ($E_\infty\to0$ as $N\to0$), so total strip stays dynamically
reachable. The remaining open item is only *energetic favorability* — where the
liberated electrostriction energy goes on full collapse (OQ6/OQ7), not a
mechanical block.

**Stability of the pickup↔gate loop — the equilibrium is provably emergent
(2026-06-21).** The coupled mean-field flow of the two locked channels is

$$
\langle \dot n \rangle
=
\lambda(n)-k\,\mathbf{1}_{\mathrm{open}}.
$$

$$
\langle \dot E_{\mathrm{int}} \rangle
=
-\frac{E_{\mathrm{int}}}{\tau}
+\lambda(n)\,f_{\mathrm{ret}}\,D_0(n+1)
-k\,D_0(n)\,\mathbf{1}_{\mathrm{open}}.
$$

with

$$
\lambda(n)
=
\lambda_0
\left(\frac{\rho_{\mathrm{He}}}{\rho_{\mathrm{bulk}}}\right)
\left(1-\frac{n}{n^*}\right)_+^{\,p},
\qquad
k
=
\nu\left(1-\frac{D_0(n)}{E_{\mathrm{int}}}\right)^{s-1}.
$$

The first equation has units of $\mathrm{ps}^{-1}$ and the second of $\mathrm{eV\,ps}^{-1}$. Two results close the
"emergent, not a parameter" claim:
- **Crossing is guaranteed (pathwise Lyapunov).** Let the self-unbound margin be
  $G\equiv E_\text{int}-\Sigma(n)$, $\Sigma(n)=\sum_{i\le n}D_0(i)$ (eV); the gate
  is suppressed while $G>0$. Between pickups $\dot G=-E_\text{int}/\tau\le0$; at a
  pickup, $\Sigma$ gains a **full** rung while $E_\text{int}$ gains only
  $f_\text{ret}$ of it, so $\Delta G=-(1-f_\text{ret})D_0(n{+}1)\le0$. Hence **$G$
  is monotone non-increasing along every sample path**, strictly decreasing
  whenever $E_\text{int}>0$ — so cooling drives $G\to0$ in finite time and the gate
  **always opens**. Pickup is on the *stabilizing* side: accretion deposits a whole
  rung of binding but only a fraction of heat, the $(1-f_\text{ret})$ remainder
  radiating to the bath (so the earlier "longer suppression → runaway pickup" worry
  had the sign backwards). The $f_\text{int}$ floor (S2) is exactly $G(0)>0$ —
  necessary and sufficient for this suppression-then-crossing structure.
- **Terminal $n$ is a stable freeze-out attractor, set by one dimensionless group.**
  Treating $E_\text{int}$ as fast ($\lambda\tau\lesssim1$), its quasi-steady value
  at fixed $n$ is $E_\text{int}^\text{qs}=\Pi(n)\,D_0(n)$ with
  $$
  \boxed{\;\Pi(n)\equiv\lambda(n)\,f_\text{ret}\,\tau\;}\qquad(\text{dimensionless: }\text{ps}^{-1}\!\cdot1\cdot\text{ps}\ \checkmark).
  $$
  Shedding requires $E_\text{int}^\text{qs}>D_0(n)$, i.e. $\Pi(n)>1$; for $\Pi(n)<1$
  shedding shuts off and $n$ **freezes** (stable: a fluctuation above $D_0$ is
  drained by cooling; an extra pickup cools back; the large-$s$ near-step sharpens
  the shut-off). **Two stabilizing routes to $\Pi<1$:** exit-driven
  ($\rho_\text{He}\to0$ as the ion leaves) and filling-driven ($n\to n^*$ via the
  A12 cap). The first guarantees termination for **every** ejection trajectory
  regardless of the cap; the second makes a resting/slow ion saturate rather than
  diverge. So terminal $n$ is well-posed and stable for all parameters — no
  fine-tuning, no runaway, no limit cycle (K2 is dissipative and the only source
  decays through the exit/fill gates). **$\Pi$ is the quantitative spine of the R1
  regime axis:** shell-retaining vs stripping is just how long $\Pi>1$ persists
  (dense-traversal time vs $\tau$) before it crosses 1. $\Pi(t)$ and $t_\times$ are
  the two derived diagnostics of the early dynamics (CALIBRATION_MAP).

---

## 7. Schema and config implications (locked)

- **`IonCheckpoint` → v6** (already required by §2.9 for any non-`fixed`
  scenario). Additional to the §2.9 changes:
  *(NB 2026-07-01, Tier-2 Phase-C reconcile: the live schema slot **v6** was
  already consumed by the Tier-1a mass-dynamics bump, so `E_int_eV` actually
  lands at **v6→v7** — a v6→v7 back-compat shim synthesizes zeros (+ a load-time
  warning) for pre-reservoir files. The "→ v6" wording here predates the Tier-1a
  v6 use of the slot; **version numbers only, mechanism unchanged**. See
  `docs/drag_port/Tier2/TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §2.2/§8 +
  `docs/drag_port/Tier2/drag_migration_log_tier2.md`.)*
  - new per-step field `E_int_eV (2N, T)` — internal-energy trajectory;
  - `mass_history_kg` monotonicity guarantee **dropped** (non-monotone by
    construction);
  - `E_mass_attach_defect_eV` → `E_mass_transfer_eV` (sign covers gain and loss);
  - scenario-metadata field records `mass_scenario = biphasic_energy_gated` so
    downstream tools interpret the arrays correctly.
    *(NB 2026-07-01, Tier-2 Phase-C: the production **config literal is
    `biphasic`** — `biphasic_energy_gated` is the documentation-only full name;
    the delivered `MassScenario` set is `{fixed, biphasic, anchored_discrete}`.
    See the §11 NB below.)*
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

- *$D_e$ vs $D_0$ — RESOLVED (2026-06-17, EPAPS fit).* 143.9 cm⁻¹ is the well
  depth $D_e$, not the zero-point-corrected dissociation energy $D_0(1)$ the
  cascade gate needs. The old estimate $D_0\approx125$–135 cm⁻¹ borrowed a
  Na⁺-like ~10% ZPE fraction (Na⁺–He: $D_e\approx285$ vs [Nat23] $D_0=270$, a
  ~15 cm⁻¹ gap) and was **too high**: He–I⁺ is a shallow, long-bond, light-$\mu$
  well, so its ZPE is a much larger fraction of $D_e$. Solving the **exact $J{=}0$
  radial Schrödinger equation** on the fitted [IHe05] $X_2$ curve ($V''(R_e)=748.1$
  cm⁻¹/Å², 5 bound levels) gives a true ground-state ZPE $G(0)=37.0$ cm⁻¹
  (**26% of $D_e$**, vs the harmonic 40.3), so
  $$
  D_0^{\,\mathrm{I^+}}(1)\big|_{X_2} = 106.9\ \text{cm}^{-1}
  = 13.3\ \text{meV} = 0.01325\ \text{eV}\quad(\pm3\ \text{cm}^{-1}\ \text{from the}
  \ \text{[IHe05] }\pm3\%\ \text{well-depth accuracy}).
  $$
  This is **pinned** — no harmonic assumption, no Na⁺ analogy. The shape for $n>1$
  remains open (below); only the first rung is fixed. *Electronic-picture
  dependence (A10):* 106.9 cm⁻¹ is the $X_2$-only rung; the statistical-mixture
  rung is $74.4$ cm⁻¹ (~70%; see A10).
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

**Open-shell blurs angular sub-structure, but the radial cliff is geometric and
survives (REVISED 2026-06-17).** The earlier lean was "open-shell ⇒ gradual,
shell-less" (Pb⁺Heₙ: SO coupling smears solvation shells). That argument is now
**only half right**, and the corrected picture *favors* shell structure:
- *Angular (within-shell):* anisotropy/SO mixing does wash out sharp angular
  sub-structure, supporting a **mild, gradual decline within shell 1** (the
  Pb⁺-like part of the argument survives).
- *Radial (between-shell):* the dominant binding is the isotropic
  charge-induced-dipole $-D_4/R^4$ ($D_4=11852$ cm⁻¹·Å⁴, [IHe05] EPAPS, $=\tfrac12
  \alpha_\text{He}e^2$ to 0.4%), which is **not** affected by electronic openness.
  *Geometry from [GAH25] Table II (corrected 2026-06-17):* the first shell sits at
  $r_1^e\approx4.67$ Å (interpolating Rb⁺ 4.5 / Cs⁺ 4.8 Å) — far outside the pair
  $R_e=3.25$ Å — and shell-2 at $r_2^e\approx7.72$ Å. So the pair binding drops by
  $(r_2^e/r_1^e)^4=(7.72/4.67)^4\approx\textbf{7.5×}$ across the shell boundary
  (lower still with inner-shell screening). A **radial cliff at $n^*$ is therefore
  expected even for an open-shell ion** — it is geometric, not electronic.
- *He–He is roomy at the true shell radius (CORRECTED 2026-06-17):* at the
  [GAH25] shell radius $r_1^e\approx4.67$ Å (not the pair $R_e$), the shell-1
  He–He nearest-neighbour spacing for $n\!\sim\!20$ is **3.97 Å $>$ He–He $R_e$
  2.97 Å → roomy, mildly *attractive*** (an earlier note using the pair $R_e$
  wrongly found 2.70 Å / compressed). So He–He gives a small positive
  contribution, not a penalty. The within-shell decline instead comes from the
  He being **pushed outward** to $r_1^e$, where pair polarization is only
  $D_4/r_1^{e4}\approx25$ cm⁻¹ — but the collective snowball lifts the *net*
  per-atom binding back to $\approx118$ cm⁻¹ ($|S|/n^*$, $n^*{=}21$), close to the lone-He
  $D_0(1)=107$. This is exactly why the pair ladder undercounts $|S|$
  (electrostriction-dominated; see K2 split, R12).
Net: the physically-motivated $D_0^{\,\mathrm{I^+}}(n)$ is **mild within-shell
decline, then a ~7.5× radial cliff to a bulk-He floor** — closer to *structured*
than *gradual*, with the in-shell slope the genuinely uncertain part. (The pair
ladder anchored at $D_0(1)$ remains the right *shed-cost* object: the collective
marginal $\partial|S|/\partial n\approx118$ cm⁻¹ ≈ $D_0(1)$, a fortunate
near-cancellation — K2.)

*Resolution — Form U, one continuous steepness knob (ADOPTED 2026-06-17).* The
old discrete `{gradual, shell_structured}` choice is replaced by a **single
sigmoid family** spanning both as limits, anchored exactly at the pinned rung and
a sourced floor:
$$
D_0^{\,\mathrm{I^+}}(n)=D_\text{floor}+\big(D_0(1)-D_\text{floor}\big)\,
\frac{1-\sigma(n)}{1-\sigma(1)},\qquad \sigma(n)=\big[1+e^{-\kappa(n-n^*-\tfrac12)}\big]^{-1},
$$
with $D_0(1)$ pinned (picture-dependent: 106.9 / 74.4 cm⁻¹), $D_\text{floor}=
|\mu_\text{He}^\text{bulk}|\approx4.97$ cm⁻¹ (7.15 K, sourced), $n^*\approx21$
(sourced, [I2-notes]), and **the single Free shape knob $\kappa$**:
$\kappa\!\to\!0$ → gradual quasi-linear decay; $\kappa\!\gg\!1$ → sharp cliff at
$n^*$. *Dimensional check:* $\sigma$ dimensionless ($\kappa$ per-unit-$n$,
$(n{-}n^*{-}\tfrac12)$ dimensionless), bracket dimensionless, $D_0(n)$ in cm⁻¹ ✓.
*Centering at $n^*+\tfrac12$ (refinement 2026-06-17):* the cliff sits **between**
the last in-shell atom ($n^*$) and the first shell-2 atom ($n^*{+}1$), so shell-1
rungs are counted at full depth and the drop lands on shell 2. *Physical prior:*
the 7.5× radial cliff argues **large $\kappa$**. The Tier-2 size
distribution arbitrates $\kappa$ continuously (broad histogram → small $\kappa$;
magic peak at $n^*$ → large $\kappa$), co-fit with the electronic picture (the two
are **not** separable — the mixture shares the *same* floor but a lower top, so it
is a more-compressed ladder at fixed $\kappa$). *Decoupling (2026-06-17):* with the
cliff at $n^*+\tfrac12$, the **integrated first-shell sum $\sum_{i=1}^{21}D_0$ is
nearly $\kappa$-independent** (~10% over $\kappa\in[0.3,5]$) — $\kappa$ shapes the
cross-cliff/shell-2 region (the histogram) but **not** the self-bound gate
threshold or the S2 floor, which are therefore picture-set and pinnable without
the $\kappa$ fit (§6 R9, S2). *Fallback (declared, not default):*
a single sigmoid is monotone and single-cliff; if the histogram shows a peak at a
non-closure $n$ (a geometric magic number, e.g. icosahedral 12/13) or a second
cliff, revert to a tabulated ladder.

**Integrated-ladder cross-check — corrected (2026-06-17): $\sum_i D_0 \neq |S|$.**
The DFT first-shell solvation $|S_{\mathrm{I^+}}|=0.308$ eV $=2484$ cm⁻¹ implies a
mean rung 118 cm⁻¹ — *above* the pinned $D_0(1)=106.9$. A monotone pair ladder
therefore **cannot reach $|S|$** (flat ceiling $21\times106.9=2245$ cm⁻¹, short by
≥239 cm⁻¹; ≥922 for the mixture). So $|S|$ is a **collective** quantity
(electrostriction/snowball compression + DFT correlation) that pair dissociation
energies structurally undercount; it is an **upper bound, not a rung-sum target**,
and rungs must **not** be calibrated to $|S|/n^*$. The reachable integrated
cross-check is instead the §6.5.1 *drag effective binding* $E_\text{bind}=0.1168$
eV $=942$ cm⁻¹, which a moderate-$\kappa$ ladder meets naturally (representative
sums 835–1571 cm⁻¹) — retained as the consistency probe (with the OQ1 double-count
caveat), while $|S|$ bounds the many-body excess (~240+ cm⁻¹). *Mitigation status:*
first rung + floor + closure sourced; one continuous shape knob $\kappa$; all
reversible if OQ1 or the $|S|$-reference question (new OQ, §10A) reopens.

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

**R12 — Hot-structure binding $\neq$ static ladder sum (LOW–MEDIUM, GATED;
equilibrium part RESOLVED 2026-06-17).** Two layers, now separated:
- *Equilibrium layer (resolved by the K2 binding split).* The old reconstruction
  $E_\text{int}=E_\text{solv.struct}+\sum_i D_0$ undercounted the binding by the
  collective excess $|S|-\sum_i D_0$ (electrostriction; ~0.03–0.06 eV at large
  $\kappa$, and per GAH25 the *dominant* shell-binding term). Using the **split**
  binding $E_\text{bind}^\text{pair}+E_\text{elec}$ in both the cooling target and
  the reconstruction makes $E_\text{int}^\text{eq}=0$ exactly — the systematic is
  eliminated, not just gated (K2).
- *Hot-transient layer (still gated, A9).* Even with the correct equilibrium
  split, the hot early structure's *instantaneous* binding deviates from the
  equilibrium value (multi-peak, distorted, GAH25 Fig. 3). *Mitigation (A9):*
  never reconstruct $E_\text{int}$ inside the gate window; trust only
  $E_\text{solv.struct}$ there. *Open item:* quantify the deviation right at
  $t_\times$ ($E_\infty$ reached only ~11+ ps in GAH25) — bounds the error at the
  one instant it is first used.

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
$n$. *Two things being tested rather than assumed:* (i) the *shape* — now a
**single Form-U sigmoid** with one continuous steepness knob $\kappa$ spanning
gradual ($\kappa\!\to\!0$) to sharp-cliff ($\kappa\!\gg\!1$), anchored at the
pinned rung, sourced bulk floor, and sourced $n^*$ (R3, adopted 2026-06-17);
physical prior favors large $\kappa$ (geometric radial cliff) — discriminated by
the Tier-2 size distribution; (ii) the scalar-multiplet collapse — which
**electronic picture** sets the rung depth (statistical mixture vs $X_2$-only) is
itself a decision (A10), co-fit with $\kappa$ (not separable: shared floor, lower
top → more-compressed mixture ladder), and the single-scalar collapse is adequate
if terminal-$n$ statistics match Tier 2, falsified if the size distribution shows
structure only an SO-resolved ladder reproduces. *Cross-check, corrected
(2026-06-17):* the integrated ladder $\sum_i D_0 \neq |S_{\mathrm{I^+}}|$ — a
monotone pair ladder cannot reach the DFT $|S|=2484$ cm⁻¹ (collective excess
≥239 cm⁻¹), so check $\sum_i D_0$ against the reachable §6.5.1 drag binding
(942 cm⁻¹), treating $|S|$ as an upper bound only. *Tighten:* resolve OQ1 and the
$|S|$-reference question (§10A).

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
the parameter-free gate rather than a loose tolerance. *Reworded (2026-06-21,
scenario-keyed):* the self-unbound onset is **robust across both scenarios** but
with **scenario-dependent margin**. The S2 floor
$f_\text{int}^\text{floor}=\sum_i D_0/E_\text{avail}^\text{ion}$ scales as
$1/E_\text{avail}$: $\approx0.065$ (mix) / $0.09$–$0.10$ ($X_2$) at the 2.70 eV
production budget, rising to $\approx0.21$–$0.24$ / $0.31$–$0.35$ at the 0.80 eV
$d{=}9$ Å **validation** budget (Tier 0/1). The floor stays far below the hard cap
$f_\text{int}\le1$ in **both** cases, so the violent explosion clears it and the
qualitative GAH25 onset transfers without a tuned partition — but the headroom
shrinks from $\sim14\times$ (production) to $\sim4.5\times$ (validation), so the
*tight* self-unbound test is the Tier-1 one, not the production one. Keep "robust";
the earlier "at most ~10%, fine-tuned-free for essentially the entire range" holds
only at the 2.70 eV budget (S2).

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
*Amendment under the binding split (2026-06-17):* with $E_\text{bind}$ now
pair + electrostriction (K2), a cold shed changes $E_\text{bind}^\text{pair}$ by
the pair $D_0$ *and* $E_\text{elec}$ by the marginal electrostriction
($\partial|S|/\partial n-D_0$). To keep $E_\text{int}$-neutrality, the **marginal
electrostriction release books to the bath** (`E_dissip`) per shed event — small
($\partial|S|/\partial n\approx118$ cm⁻¹ vs the lone-He $D_0(1)=107$, a ~15%
marginal), so A8 is now "neutral for $E_\text{int}$ up to a bath-booked marginal,"
not exactly neutral. On *full* stripping the accumulated $E_\text{elec}$ (the
snowball-collapse energy) is liberated; whether it radiates to the droplet or
adds to ejected-ion KE is open (OQ6, bundled with OQ7).

**A9 — Static-ladder $E_\text{int}$ reconstruction is valid only after
$t_\times$ (locked rule; reconstruction binding updated 2026-06-17).** Recovering
$E_\text{int}=E_\text{solv.struct}-E_\text{bind}^\text{pair}(N)-E_\text{elec}(N)$
now uses the **split** equilibrium binding (pair + electrostriction, K2), so at
equilibrium $E_\text{int}^\text{eq}=0$ exactly — the old collective-excess offset
is gone (R12 resolved). It nonetheless remains a reconstruction of the
*equilibrium* binding: for a hot, distorted, multi-peak early structure (GAH25
Fig. 3) the true instantaneous binding $\neq$ the equilibrium split, so the
reconstruction still errs in the early window. **Rule unchanged:** never use the
per-rung $E_\text{int}$ vs $D_0(n)$ comparison (or the reconstruction) inside the
gate window ($E_\text{solv.struct}>0$); only $E_\text{solv.struct}$ itself is
trusted there. Post-$t_\times$, the split binding is valid. The gate enforces this
automatically (R12). *Tighten:* check how far the hot-structure binding deviates
from the equilibrium split right at $t_\times$ ($E_\infty$ not reached until
~11+ ps in GAH25).

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

**First-rung numbers, both pictures (2026-06-17, EPAPS + [IHe05] Eq. 9).**
Evaluating the SO-coupled curves from $V_\Pi$ (³Π) and $V_\Sigma$ (³Σ⁻) with the
atomic ³P$_j$ splittings $D_0^{at}=6451$, $D_1^{at}=7090$ cm⁻¹, and ZPE-correcting
each via its exact $J{=}0$ ground state:

| SO state | $D_e$ (cm⁻¹) | $D_0$ (cm⁻¹) | $D_0$ (eV) |
|---|---|---|---|
| $X_2$ (³Π) | 143.9 | **106.9** | 0.01325 |
| $I_1$ | 89.5 | 62.4 | 0.00773 |
| $I_0$ | 79.1 | 54.0 | 0.00670 |
| **equal-weight mixture** | — | **74.4** | **0.00923** |

So the production-default (mixture) first rung is **$D_0(1)=74.4$ cm⁻¹ ≈ 9.23 meV**,
**~70% of the $X_2$-only rung (106.9 cm⁻¹)**. The two pictures are now separated by
a concrete factor, not a qualitative "shallower."

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
x2_only, cooling_relaxed}` (§11; `cooling_relaxed` added 2026-06-21 — rung
between mixture and $X_2$, see the 2026-06-21 revision item 3).

**A11 — Evaporation kinetics are classical RRK, with $s$ absorbing the
simplification; mode-count and small-$n$ validity corrected 2026-06-21.** The shed
rate (§4) uses classical Rice–Ramsperger–Kassel form,
$k=\nu(1-D_0/E_\text{int})^{s-1}$, the simplest unimolecular kinetics. This
deliberately ignores the *quantum* mode structure and zero-point energy of the
I⁺Heₙ complex — which matters here because He modes are soft and the system is
cold, so the more defensible form would be RRKM (explicit density/sum of states).
RRKM is not adopted because it needs the I⁺Heₙ vibrational spectrum, which does not
exist ($n>1$, A10/R3). *Position:* classical RRK is a stated simplification; the
**effective DOF $s$ absorbs the RRKM/quantum/ZPE difference** — mode-counted at
$s=3(n{+}1)-6=\mathbf{3n-3}$ (full $n{+}1$-atom complex) by default, calibratable
as an effective scalar, so the classical form is a parametrized stand-in, not a
first-principles claim. **DOF counting — why $3n-3$, not $3n-6$ (recorded
2026-06-21).** $s$ must count the *internal vibrational* modes among which
$E_\text{int}$ randomizes — i.e. the standard $3N-6$ (nonlinear): total $3N$
Cartesian DOF minus 3 c.o.m. translations minus 3 overall rotations, neither of
which is part of the dissociating reservoir. The complex is **I⁺Heₙ $=n{+}1$
atoms**, so $s_\text{vib}=3(n{+}1)-6=3n-3$ — the full subtraction *is* applied;
the $-6$ is present, the ion is simply included in $N$. The retired $3n-6$ was
$3N-6$ with $N=n$, dropping the ion from the atom count — which then
double-removes symmetry: a He-only cluster moving against a fixed ion anchor has
**no free overall translation or rotation to subtract** (the ion potential breaks
both), so that picture counts $3n$, not $3n-6$. The honest bracket is therefore
$$
3n-3 \;\le\; s_\text{vib} \;\le\; 3n,
$$
with $3n-3$ the **conservative free-complex default** (the 3 shell rotations
treated as free and removed) and $3n$ the stiff limit (those 3 treated as soft
**librations** that *do* hold and exchange energy on dissociation timescales,
since a real snowball is a hindered, not free, rotor). The old $3n-6$ sat *below*
this entire physical band. The libration question — how far up the $[3n-3,3n]$
band the effective $s$ rides — is **exactly the shell-rigidity question that the
effective-$s$ override and $\kappa$ co-calibrate** (coupling below), so the band
is the sensitivity range, not an error bar to eliminate. *Boundary cases:* at
$n{=}2$ (I⁺He₂, 3 atoms) the default $3n-3=3$ assumes a **bent** geometry; a
**linear** He–I⁺–He gives $3N-5=4$, a $+1$ shift that is non-negligible at
single-digit $s$ — flag linear as the alternative. At $n{=}1$ (diatomic) $3N-6$ is
undefined and $3N-5=1$: statistical RRK with $s{=}1$ gives constant $k=\nu$, which
**is** the direct-dissociation rule adopted in §4 — so that special case is not a
patch but the $s{=}1$ linear-diatomic limit of the same formula.
**Boundedness requires $s\ge1$ (corrected 2026-06-21):**
the rate saturates $k\in[0,\nu)$ only for $s>1$; $s=1$ gives constant $k=\nu$; and
$s<1$ makes $k\to\infty$ at threshold ($E_\text{int}\to D_0^+$), reintroducing the
gate-open avalanche the RRK rate was added to remove. The earlier $s=3n-6$ went
$s{=}0$ at $n{=}2$ and $s{=}{-}3$ at $n{=}1$ — divergent — so it broke the
boundedness guarantee in the small-$n$ regime. The corrected $3n-3$ is $\ge3$ for
all $n\ge2$; a config-load guard enforces $s\ge1$ on any override. **The last atom
($n{=}1$) is handled as direct dissociation $k=\nu$** (single mode → statistical
RRK degenerate; §4), the explicit boundary of the statistical picture (and where
cold-shed neutrality A8 is weakest). *Why this regime matters:* at large $n$ the
huge exponent makes $k$ a near-step function, so terminal $n$ there is set by
*cooling crossing the gate*, not by RRK kinetics — the genuine $\{\nu,s\}$
sensitivity of the size-distribution tail lives in the **small-$n$ cascade**
($n\lesssim$ few), exactly the regime the $3n-6$ error corrupted and the
total-stripping secondary run (OQ6) must traverse. *Coupling (flagged):* $s$ and
the ladder shape (A10/R3) both probe shell rigidity/separability and are **not
independent** — a blurred/gradual shell lowers the effective $s$ below $3n-3$;
calibrate them **jointly** against the size distribution. *Cross-check / tighten:*
the I⁺ cascade timescale from [I2-notes] (OQ5, §10A) pins the $\{\nu,s\}$
combination the same way GAH25's window pinned $\{f_\text{int},\tau\}$; falsified if
the size distribution requires an RRKM-shaped, non-power-law switch-on.
**Resolution NB (2026-07-06 — staircase probe + s_eff mini-probe + picture
cross-check).** The flagged down-drift is the empirical outcome, and it is
large: at $s=3n-3$ no in-band (κ, picture, τ) point sheds more than 1.7 of the
anchored 7 He (kinetic freeze — the $s{-}1=59$ exponent on
$x=D_0(21)/\Sigma(21)\approx0.032$), while a **constant
$s_\text{eff}\approx8$–12 — an order of magnitude below the classical band —
lands magnitude and timing simultaneously** (s_eff=8, τ=6.55: 7.41 sheds →
n\_end 13.59, first shed 5.42 ps vs t★=5, trajectory MAD 1.0 He), robust
across all three electronic pictures (≤4 % magnitude spread). The classical
bracket $[3n-3,3n]$ survives as the *rigid-classical limit*, not the sweep
band; the effective reservoir is ~8–12 modes (quantum freezing / weak coupling
on the shed timescale — exactly the RRKM/quantum/ZPE difference $s$ was
declared to absorb). The κ-coupling flagged above resolved **weak**: κ is
inverted and normalisation-capped for the stripping range (Form U floors
$D_0(21)\approx0.53\,D_0(1)$), so it cannot substitute for $s$. **Promotion
(user, 2026-07-06):** $s$ reclassified Derived → **Bounded** constant
$s_\text{eff}$ (band ≈[5, 20]; CALIBRATION_MAP row 10 + 2026-07-06 update
block); $s\ge1$ guard, $n{=}1$ direct-dissociation limit, and the small-$n$
boundary flags all unchanged. Execution records:
`docs/drag_port/Tier2/drag_migration_log_tier2.md`.

**A12 — Pickup is occupancy-capped (Langmuir site saturation), inert for ejection
but required for resting-ion correctness (NEW 2026-06-21).** The pickup rate
carries a blocking factor $\lambda_\text{attach}=\lambda_0(\rho_\text{He}/
\rho_\text{bulk})(1-n/n^*)_+^{\,p}$ (§4): a sticking coefficient that falls
smoothly as the first shell fills, $\to0$ at $n{=}n^*$. *Why needed:* the
density-only form has no $n$-dependence, so termination of accretion is purely
exit-driven — correct for the ejection problem (the ion always leaves), but a
resting/slow-exit ion would accrete unboundedly ($n\to\infty$), which is
unphysical and contradicts [Nat23]'s leveling-off of resting-Na⁺ accretion. The
cap closes that gap. *Form choice (Langmuir, "Form B"):* caps the **rate**, never
the per-event energetics — so it leaves the §6 invariant and the $G$-Lyapunov
crossing (§6.11) untouched (each pickup still deposits a full rung of binding,
$f_\text{ret}$ of it as heat) and only *strengthens* the freeze-out attractor by
adding a second stabilizing route $\Pi(n)<1$ at $n\to n^*$ (§6.11). Rejected
alternatives: a hard wall $\mathbb 1_{n<n^*}$ (makes the $\pm1$–2-uncertain $n^*$
a discontinuous hard input, OQ4/OQ8; bad for the BAOAB/jump split) and a
two-reservoir shell-2 routing ("Form C", faithful but doubles bookkeeping and
reopens the cliff-timing OQ — deferred unless Tier-2 shows shell-2 population).
*Inert for production:* the ion exits before the shell saturates, so $\Pi$ crosses
1 via $\rho_\text{He}\to0$ first; Form B and density-only give identical ejection
trajectories. It earns its keep only if the Tier-1 $21\to19\to14$ tail lingers
near saturation, and for in-principle resting-ion validity. *Coupling (flagged):*
the exponent $p$ and `ladder_steepness` $\kappa$ both encode first-shell abruptness
from different observables (pickup shut-off vs binding cliff), so they are **not
independent** — same situation as $s\leftrightarrow\kappa$ (A11). **Default: $p$
tied to $\kappa$** (single shell-rigidity parameter, zero net new free knobs).
*NB 2026-07-01 (Tier-2 Phase-B §3.2): this tie is **inverse** — a rigid shell =
large $\kappa$ = **smaller** $p$ — so it is **not** a literal $p=\kappa$; if ever
tied, match cutoff **slopes**, not values. Phase B holds $p=1$ fixed.* Split into
an independent bounded/free $p$ only if the size distribution demands
it (Tier-2). *Falsify/tighten:* a size distribution whose first-shell cutoff
sharpness is inconsistent with the $\kappa$-implied $p$ would force the split;
[Nat23] resting-ion saturation level bounds $n^*$ and $p$ jointly.

**A13 — Integrator/mass-jump operator split: accuracy, ordering, and the
momentum-reset invariant precondition (NEW 2026-06-21).** The continuous BAOAB
integrator (DRAG doc) and the discrete mass-jump process (this doc, §4) meet once
per step; their composition has three pieces, settled here as physics-definition
ahead of the (later) jump implementation.

*SQ1 — velocity-dependent drag in the O-step (built, accepted as-is).* The cubic
gives a state-dependent effective friction $\gamma(v)=g\,b\,v^2/m$ ($[\gamma]=
\text{ps}^{-1}$ ✓), so the O-step freezes $\gamma(v_\text{in})$ and applies
$e^{-\gamma dt}$ rather than the exact nonlinear flow $v(t)=v_0/\sqrt{1+2(gb/m)v_0^2t}$.
**Consequence:** the drag-on path is **globally $O(dt)$**, not $O(dt^2)$ — the
BAOAB second-order/configurational-superconvergence claim holds only in the
constant-$\gamma$/no-drag limit and **must not be claimed for production**. *Why
the trade is correct here:* freezing $v_\text{in}$ buys two unconditional
properties worth more than an order — (i) **exact dissipation bookkeeping**
$\Delta E_\text{dissip}=\tfrac12 m(\|v_\text{in}\|^2-\|v_\text{out}\|^2)$ for any
$dt$ (the invariant's drag term closes by construction, not to $O(dt^2)$), and
(ii) **unconditional dissipativity** $\|v_\text{out}\|\le\|v_\text{in}\|$ since
$e^{-\gamma dt}<1$, i.e. no large-$dt$ energy-injection blow-up on a stiff cubic.
*Residual:* a **one-signed over-braking bias** (frozen $\gamma$ uses the largest
$v$ in the step), monotone in $dt$, folded into R10 tolerance — so the $dt$
convergence check should expect simulated peak $v$ slightly **below** TDDFT
(Fig. 23) at fixed $dt$, and "tolerance-only" should note the bias has a sign.

*SQ2 — mass-jump placement and the conservation precondition (spec; unbuilt).*
(a) **Momentum-conserving velocity reset is mandatory, not an accuracy choice**
(§4): $v^+=(m v^-+m_\text{He}u_\text{He})/(m\pm m_\text{He})$, He at rest. **The
five-term invariant (§6) closure *presupposes* this reset** — $E_\text{mass\_transfer}$
is *defined* as the reduced-mass KE defect it produces, $\tfrac12\tfrac{m\,m_\text{He}}
{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2$. Carrying $v$ through the jump as a mere
label injects/removes KE and **voids the closure proof**; the integrator's reset
and the invariant's $E_\text{mass\_transfer}$ term must use the *same* $v^+$.
(b) **Jump-step order reduction is benign — record, don't fix.** A jump makes the
step non-palindromic → locally $O(dt)$; but jump-steps number $\sim(\lambda+\nu)T$,
**independent of $dt$**, so they are a vanishing fraction as $dt\to0$ and the
global order is unchanged from the $O(dt)$ SQ1 already concedes. Precondition
$\lambda dt,\nu dt\ll1$ holds ($\sim0.01$–$0.024$). (c) **Fixed ordering:** jump
*then* O (so the O-step sees $v^+$ and $\gamma(v^+)$); the $[\mathcal L_M,\mathcal
L_O]dt^2$ commutator is within tolerance. (d) **At most one mass event per step**
(shed before pickup if both fire; §4).

*SQ3 — post-jump mass in the O-step (spec; unbuilt, forced by SQ2).* Any O-step
following a jump uses the **post-jump mass** $m^+$ in both friction and the FDT
noise amplitude $\sqrt{(1-e^{-2\gamma dt})k_BT_\text{eff}/m^+}$, consistent with
the $v^+$ it perturbs. Mismatch breaks FDT balance by $m^+/m^-$. Small (N2 bath
noise is weak, §1.3a), but a definiteness requirement; the $1/m$ also correctly
makes a heavier post-pickup complex receive a smaller thermal kick.

*Net:* no structural defect — SQ1 is a documented accuracy-for-conservation trade,
SQ2/SQ3 are now specified. **The Method's consistency proof closes conditional on
the SQ2(a) momentum reset being implemented as specified** (the one item that, if
missed, retroactively breaks §6). Tracked as an implementation precondition, not an
open physics question. *Config:* §11 (`mass_jump_velocity_reset`,
`he_capture_velocity`, `one_mass_event_per_step`, `jump_o_step_ordering`).

---

## 10. Open items / calibration targets

| Item | Symbol | Source / target | Tier |
|---|---|---|---|
| Pickup rate coefficient | $\lambda_0$ | OOM prior only ([GAH25] well-depth dep.); pin from I⁺ TDDFT + size dist. | 1 / 2 |
| Pickup occupancy exponent | $p$ | **NEW 2026-06-21 (A12):** Langmuir cap $(1-n/n^*)_+^{\,p}$; **default tied to $\kappa$** (one shell-rigidity knob, 0 net new free) — **NB 2026-07-01: the tie is *inverse* (large $\kappa$ → sharper cutoff → *smaller* $p$), NOT literal $p=\kappa$; Phase-B holds $p=1$ fixed (§3.2)**; split to bounded/free only if size dist. demands; [Nat23] saturation bounds $n^*,p$ | 2 (conditional) |
| Pickup↔gate order parameter | $\Pi(t)=\lambda(n)f_\text{ret}\tau$ | **NEW 2026-06-21, derived diagnostic (§6.11):** $\Pi{>}1$ shedding persists / $\Pi{<}1$ freeze; regime-axis spine (R1); $\Pi\to0$ at exit guarantees termination; reconstructable post-hoc like $t_\times$ | — |
| First ladder rung | $D_0^{\,\mathrm{I^+}}(1)$ | **pinned (2026-06-17):** [IHe05] EPAPS fit, exact $J{=}0$ ZPE → $X_2$ **106.9 cm⁻¹** (0.01325 eV); mixture **74.4 cm⁻¹** (0.00923 eV); $\pm3$ cm⁻¹ | — (sourced) |
| Ladder shape ($n>1$) | $\kappa$ (Form U sigmoid) | single steepness knob, gradual↔cliff; prior large (7.5× radial cliff, R3); co-fit w/ picture + $\{\nu,s\}$; discriminated by size dist. | 2 (Free) |
| Ladder floor | $D_\text{floor}$ | $\|\mu_\text{He}^\text{bulk}\|\approx4.97$ cm⁻¹ (7.15 K), bulk superfluid; picture-independent | — (sourced) |
| Integrated ladder | $\sum_i D_0(i)$ | **Form U, corrected 2026-06-21:** $X_2$ **0.25–0.28** eV / mix 0.17–0.19 eV (pure $\kappa$-range, ~11%, nearly $\kappa$-independent); drag binding 0.117 eV and the crowding-reduced value are **separate cross-checks**, not band ends; $\neq|S|=0.308$ (collective UB) | derived |
| Binding-release retained fraction | $f_\text{ret}$ | size distribution | 2 |
| Newton-cooling relaxation time | $\tau_\text{dissip}$ | **sweep band $[2.6,16.5]$ ps** (R8), externally anchored; not pinned from this work's size dist. | sweep |
| Newton-cooling asymptote (binding) | $E_\infty(N)=-|S(N)|$ | **occupancy-resolved (2026-06-17, split):** $E_\text{bind}^\text{pair}+E_\text{elec}$; full-shell $|S_{\mathrm{I^+}}|{=}0.308$ eV ([I2-notes]); GAH25 Na⁺ $-3424/-4144$ K calibrates form+$\tau$; OQ6 resolved | 2 (shape via $\kappa$) |
| Onset Coulomb budget | $E_\text{avail}^\text{ion}$ | **scenario-keyed (2026-06-21):** **0.80 eV** validation ($d{=}9$ Å, ½·14.40/9) / **2.70 eV** production ($R_e=2.666$ Å, ½·5.40); fixed reference, stamped to scenario guard | — (sourced) |
| Onset partition fraction | $f_\text{int}$ | $E_\text{int}(0)=f_\text{int}E_\text{avail}^\text{ion}$ (§6.11/S2); **floor scenario-keyed (2026-06-21): 0.065 (mix)/0.09–0.10 ($X_2$) @ 2.70 eV; 0.21–0.24/0.31–0.35 @ 0.80 eV**, picture-set & $\kappa$-indep; soft upper ~0.2 (advisory, not a constraint); pin from size dist. | 2 |
| Self-bound crossing time | $t_\times$ | **derived diagnostic, not fitted**; cross-check vs GAH25 ~5–6.5 ps (±factor-2) and size dist. (§6.11) | — |
| Early-instability gate | **derived, not fitted** | self-bound criterion $E_\text{int}<\sum_i D_0(i)$ (R9); the ~several-ps onset is now a *prediction* vs [GAH25], cross-checked by size dist. | — |
| RRK prefactor | $\nu$ | **pinned $2.42$ ps⁻¹ (2026-06-17):** $\omega_e=80.6$ cm⁻¹ from [IHe05] EPAPS $V''(R_e)=748.1$ cm⁻¹/Å² (§4); cross-check vs [I2-notes] cascade timing (OQ5) + size dist. | 2 |
| RRK effective DOF | $s$ | **Bounded constant $s_\text{eff}$ (promoted 2026-07-06; was mode-counted derived):** band ≈[5, 20], staircase landing [8, 12]; $s{=}3n{-}3$ (corrected 2026-06-21 from $3n-6$) demoted to the classical-limit arm; $n{=}1$ direct dissociation $k=\nu$ and the $s\ge1$ guard unchanged; κ-joint resolved weak (A11 resolution NB) | 2 + 9 Å staircase (s↔τ via timing) |
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
  effective depth. **Magnitude now pinned (2026-06-17):** the mixture/$X_2$ ratio
  at the dimer level is $74.4/106.9=0.70$, so the static-mixture reduction is
  **~30%**; that is the size of the reduction that must not be double-applied if
  the drag's 0.1168 eV already carries a dynamical lowering. *If instead the
  extraction already SO-averaged:* 0.1168 eV is mixture-like and self-consistent
  with the A10 default. **Until resolved:** treat the A10 mixture default as
  physically motivated but provisional; do not cite
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
  in-hand cross-check for the pinned $\nu=2.42$ ps⁻¹ and the effective $s$.
  *If movie-only:* author-contact item. Confirm whether $N_\text{He}(t)$ is
  readable.
- **OQ6 — Full-stripping reach (mechanism RESOLVED 2026-06-17; energy release
  still open).** *Resolved:* the capping risk is real — a *fixed* full-shell
  $E_\infty$ drives the reconstructed $E_\text{int}^\text{eq}$ negative as the
  shell strips ($-0.03$ eV at $N{=}21$ → $-0.28$ eV at $N{=}2$) and halts shedding.
  The **occupancy-resolved $E_\infty(N)=-|S(N)|$** (K2 split) removes the cap
  ($E_\infty\to0$ as $N\to0$), so the Calvo24 total-strip limit stays reachable.
  *Still open:* on full stripping the accumulated $E_\text{elec}$
  (snowball-collapse energy, ~$|S|-\sum D_0$) is liberated — whether it radiates
  to the droplet or adds to ejected-ion KE bears on whether total strip is
  energetically *favored*. **Bundle with OQ7** (both need the $N$-resolved
  structure/reference of $|S|$); ask Halberstadt.
- **OQ7 — Energy reference of the DFT solvation $S_{\mathrm{I^+}}$ (NEW
  2026-06-17, LOW–MEDIUM).** The finding that a monotone pair ladder cannot reach
  $|S|=2484$ cm⁻¹ interprets cleanly as *collective excess* only if $|S|$ is
  referenced to **free I⁺ + free He**. If [I2-notes] references it to bulk He
  (chemical-potential zero), snowball-formation, or a larger structure than the
  first shell, the floor/closure anchors, the $E_\text{elec}$ magnitude, and the
  size of the "collective excess" shift. **Ask Halberstadt** (same contact as
  OQ1/OQ4) for the $S_{\mathrm{I^+}}=-3578$ K reference convention. Does not block
  Form U or the K2 split (both anchored on $D_0(1)$, floor, $n^*$, and the
  *shape* of $|S(N)|$); only sets the collective magnitude.
- **OQ8 — I⁺ shell occupancy $n^*$ and the Tier-1 endpoint (NEW 2026-06-17;
  upgraded LOW→LOW–MEDIUM 2026-06-21).** GAH25 $R_e$-scaling (I⁺ $R_e=3.25$ Å, 56%
  Rb⁺→Cs⁺ whose $n_1^e=18\to21$) gives $n^*\approx20$, corroborating [I2-notes]'s 21
  to ±1–2 but at the high edge. **Stakes raised:** the $20$ vs $21$ ambiguity is
  **no longer geometry-only** — it had leaked into the *energetics*, with the
  per-atom collective binding quoted off $n^*{=}20$ ($124$ cm⁻¹/$179$ K) where the
  cation value is $21$ ($118$ cm⁻¹/$170$ K); corrected throughout 2026-06-21 (K2,
  R3, A8, see the 2026-06-21 revision item 1). Adopt $n^*=21$ for all per-atom
  energetics. The Tier-1 endpoint "$14$" in $21\to19\to14$ (which equals Na⁺'s
  $n_1^e$) is **author-confirmed I⁺-specific via direct exchange with the [I2-notes]
  creators** (2026-06-21) — *not* a Na⁺ template; the absent §III.B.2 figure is not
  a provenance gap. Remaining open: only the $n^*=20$-vs-$21$ scaling, a clean
  GAH25 cross-check to close.

---

## 11. Interchangeability surface (config fields)

Extends `DRAG_PORT_DESIGN_DECISIONS.md` §2.8:

- `SimConfig.mass_scenario` gains value `biphasic_energy_gated` (production).
  Existing `fixed`, `scenario_A_accretion`, `scenario_B_stripping`, `biphasic`
  retained for comparison/regression.
  *(NB 2026-07-01, Tier-2 Phase-C reconcile: the delivered `MassScenario` literal
  set is **`{fixed, biphasic, anchored_discrete}`** — the production value is
  **`biphasic`** (kept as the literal; `biphasic_energy_gated` is the
  documentation-only full name), and `scenario_A_accretion` /
  `scenario_B_stripping` were superseded by the Tier-1a `anchored_discrete`
  scenario and never implemented. **Config literal names only; the §11 knobs and
  mechanism are unchanged.** See
  `docs/drag_port/Tier2/TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` §0/§8 +
  `docs/drag_port/Tier2/drag_migration_log_tier2.md`.)*
- `SimConfig.mass_initial_amu` — initial physical mass; for
  `biphasic_energy_gated` defaults to the measured ion-stage-onset shell
  (~21 He ≈ 211 amu, design §2.1), distinct from $m_\text{eff}\approx203$ amu.
- `SimConfig.pickup_rate_coefficient` — $\lambda_0$ ($\text{ps}^{-1}$).
- `SimConfig.pickup_rate_form ∈ {density_only, sweeping, dwell_time}` — default
  `density_only` ([Nat23]-supported, §5).
- `SimConfig.pickup_occupancy_cap ∈ {langmuir, none}` — default **`langmuir`**
  (Form B, A12): multiplies $\lambda_\text{attach}$ by $(1-n/n^*)_+^{\,p}$, a
  site-saturation factor → 0 at $n{=}n^*$. `none` recovers pure density-only
  (identical for production ejection, which exits before saturation). Caps the
  **rate** only; per-event energetics, the §6 invariant, and the §6.11 $G$-crossing
  are untouched.
- `SimConfig.pickup_occupancy_exponent` — $p\ge0$, the cap sharpness. **Default:
  tied to `ladder_steepness` ($\kappa$)** as a single shell-rigidity parameter
  (zero net new free knobs; A12 $p\!\leftrightarrow\!\kappa$ coupling); set to an
  independent value only to split it as a bounded/free Tier-2 knob when the size
  distribution's first-shell cutoff sharpness demands it.
  *(NB 2026-07-01, Tier-2 Phase-B refinement: the $p\!\leftrightarrow\!\kappa$ tie is
  **inverse** — a rigid shell = large $\kappa$ = sharper cutoff = **smaller** $p$ — so it is
  **NOT** a literal $p=\kappa$; if ever tied, match cutoff **slopes**, not values. The
  Phase-B production default is **$p=1$ held fixed**, freed at Phase F only if the size-dist
  first-shell edge can't be met with $p=1$ + $\kappa$. See
  `docs/drag_port/Tier2/TIER2_PHASE_B_IMPLEMENTATION_PLAN.md` §3.2.)*
- `SimConfig.mass_jump_velocity_reset ∈ {momentum_conserving, label_only}` —
  **must be `momentum_conserving`** (A13): $v^+=(m v^-\pm m_\text{He}u_\text{He})
  /(m\pm m_\text{He})$. `label_only` (carry $v$ unchanged) is **forbidden in
  production** — it voids the §6 invariant closure; retained only as a deliberate
  diagnostic to *demonstrate* non-closure. The reset and the invariant's
  `E_mass_transfer_eV` term share the same $v^+$.
- `SimConfig.he_capture_velocity ∈ {at_rest, thermal}` — default **`at_rest`**
  ($u_\text{He}=0$, droplet frame; matches cold-shed A8). `thermal` (sample
  $u_\text{He}$) deferred — adds a second noise channel, re-touches FDT (A13).
- `SimConfig.one_mass_event_per_step` — default **`true`** (A13): if both pickup
  and shed Bernoulli draws fire in a step, apply shed then pickup (fixed order).
- `SimConfig.jump_o_step_ordering` — fixed **jump-then-O** (A13): the post-jump
  O-step uses $m^+$ and $\gamma(v^+)$ in both friction and the FDT noise amplitude.
- `SimConfig.dissociation_ladder` — $D_0^{\,\mathrm{I^+}}(n)$, **Form U** (sigmoid,
  adopted 2026-06-17): $D_0(n)=D_\text{floor}+(D_0(1)-D_\text{floor})\,(1-\sigma(n))
  /(1-\sigma(1))$, $\sigma(n)=[1+e^{-\kappa(n-n^*-\tfrac12)}]^{-1}$. Anchors: first rung
  **pinned** from [IHe05] EPAPS (exact ZPE), $X_2$ 0.01325 eV / mixture 0.00923 eV
  (picture-dependent, below); floor $D_\text{floor}=|\mu_\text{He}^\text{bulk}|
  \approx4.97$ cm⁻¹ (sourced); closure $n^*\approx21$ (sourced, [I2-notes]).
  Tabulated-ladder override retained as the declared fallback (R3).
- `SimConfig.ladder_steepness` — $\kappa$ (per-unit-$n$), the **single Free shape
  knob**: $\kappa\!\to\!0$ gradual, $\kappa\!\gg\!1$ sharp cliff at $n^*$. Physical
  prior large (7.5× geometric radial cliff, R3). Tier-2-arbitrated **jointly**
  with `ladder_electronic_picture` (not separable) and the RRK $\{\nu,s\}$.
  Replaces the retired discrete `ladder_shape ∈ {gradual, shell_structured}`.
- `SimConfig.ladder_electronic_picture ∈ {statistical_mixture, x2_only,
  cooling_relaxed}` — `statistical_mixture` (equal-weight $X_2{+}I_1{+}I_0$,
  shallower; the production default per A10, motivated by violent
  Coulomb-explosion birth) vs `x2_only` (deep $X_2/³Π$ snowball, the
  [I2-notes]-consistent comparison) vs `cooling_relaxed` (**added 2026-06-21**:
  rung between mixture and $X_2$, e.g. a relaxation-weighted blend — captures a
  complex that cools toward the deepest curve rather than staying frozen in the
  birth mixture; strictly inside the existing bracket, A10/§A5). Sets the rung
  *top* (74.4 / 106.9 cm⁻¹, `cooling_relaxed` in between; shared floor → mixture
  is more compressed at fixed $\kappa$), hence $\sum_i D_0$, the gate, and
  $t_\times$. Default provisional pending OQ1 (§10A).
- `SimConfig.internal_energy_retained_fraction` — $f_\text{ret}\in[0,1]$.
- `SimConfig.internal_energy_cooling_tau_ps` — $\tau_\text{dissip}$, the
  Newton's-law-of-cooling relaxation time for $E_\text{solv.struct}$ (§6 K2;
  replaces the earlier generic "bath rate"). Carried as the $[2.6,16.5]$ ps
  sweep band (R8).
- `SimConfig.solv_struct_asymptote` — $E_\infty(N)=-|S(N)|$, the
  **occupancy-resolved** equilibrium-shell *binding* that $E_\text{solv.struct}$
  relaxes toward (§6 K2, split, 2026-06-17). NOT a fixed full-shell constant
  (that caps stripping, OQ6) and NOT an internal-energy floor. Built as
  $|S(N)|=|S_{\mathrm{I^+}}|\cdot\sum_{i\le N}D_0/\sum_{i\le n^*}D_0$ with
  full-shell $|S_{\mathrm{I^+}}|=0.308$ eV ([I2-notes]); $\to0$ as $N\to0$.
  Retires the fixed `solv_struct_asymptote_eV`.
- `SimConfig.electrostriction_binding` — $E_\text{elec}(N)=-(|S(N)|-\sum_{i\le N}
  D_0)\le0$, the collective snowball-compression binding beyond the pair ladder
  (**dominant term**, ~5× pair-at-shell-radius per GAH25). Part of
  $E_\text{solv.struct}$; its marginal release on shedding books to `E_dissip`
  (A8 amendment). Sourced (not free) from $|S_{\mathrm{I^+}}|$ + the ladder shape.
- `SimConfig.internal_energy_partition_fraction` — $f_\text{int}\in[0,1]$,
  the fraction of *this ion's* onset energy coupled into shell-internal modes;
  $E_\text{int}(0)=f_\text{int}\cdot E_\text{avail}^\text{ion}$ (§6.11/S2).
  Calibration target (§10), **bounded below by the self-unbound floor
  $f_\text{int}^\text{floor}=\sum_iD_0/E_\text{avail}^\text{ion}$** — scenario-keyed
  (2026-06-21): **0.09–0.10 ($X_2$) / ~0.065 (mix) at 2.70 eV (production); 0.31–0.35
  ($X_2$) / 0.21–0.24 (mix) at 0.80 eV ($d{=}9$ Å validation)** — picture-set,
  $\kappa$-independent. Soft upper edge ~0.2 is a **velocity-consistency
  plausibility bound, not a constraint** (Tier-2 may exceed it with a flag; S2).
- `SimConfig.coulomb_available_eV` — $E_\text{avail}^\text{ion}$,
  **scenario-keyed (2026-06-21):** $\mathbf{0.80}$ eV for **validation** (Tier 0/1,
  ionization after drift to $d\approx9$ Å, $\tfrac12\cdot14.40/9$; the budget the
  drag and 21→19→14 shell references were generated under, [I2-notes] §III.B.2,
  author-confirmed) and $\mathbf{2.70}$ eV for **production** (vertical
  double-ionization at $R_{\mathrm{II}}=R_e(\mathrm{I_2})=2.666$ Å,
  $\tfrac12\cdot14.40/2.666=\tfrac12\cdot5.40$). Fixed reference (sits in the
  uncalibrated $t^*$ window), not a fit target. **Stamped to the run's scenario tag
  and tied to the DESIGN §6.5/§6.5.1 drag↔mass guard** — a 2.70-eV onset cannot run
  against 0.80-eV-calibrated drag/shell references without tripping it. *Reversible:*
  use the pair value (1.60 / 5.40 eV) for the old per-pair convention (doubles the
  floor).
- `SimConfig.internal_energy_initial_eV` — **retired** as a free field;
  derived as $f_\text{int}\cdot E_\text{avail}$. Retained only as an optional
  manual override (default: unset → use the derived value) for sensitivity
  exploration.
- `SimConfig.evap_rate_prefactor_per_ps` — RRK prefactor $\nu$ (ps⁻¹), the
  I⁺–He stretch attempt frequency; **default 2.42 ps⁻¹** ($\omega_e=80.6$ cm⁻¹,
  pinned from [IHe05] EPAPS $V''(R_e)$, §4), cross-checked vs [I2-notes] cascade
  timing (OQ5). Governs the shed rate once self-bound.
- `SimConfig.evap_rrk_dof` — RRK effective DOF $s$, the exponent $s-1$ in the
  rate (§4). **Derived by default as $s=3n-3$** (full $n{+}1$-atom complex,
  $n$-dependent; **corrected 2026-06-21** from $3n-6$, which went $\le0$ for
  $n\le2$). **$n{=}1$ uses direct dissociation $k=\nu$** (no RRK bracket; single
  vibrational mode). Override allowed as a fixed scalar but **guarded $s\ge1$ at
  config-load** (a dissipativity-style bound — $s<1$ diverges at threshold and
  reintroduces the avalanche). Calibrated jointly with `ladder_steepness`
  ($\kappa$) (A11 coupling), not independently.
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
| **[IHe05]** Buchachenko/Viehland, JCP 122 (2005) + erratum [IHe05-E] | He–I⁺ ab initio **pair** potential (the only I⁺ source) | First ladder rung $D_0^{\,\mathrm{I^+}}(1)$ **pinned** via EPAPS fit + exact ZPE ($X_2$ 106.9 / mixture 74.4 cm⁻¹); $R_e=3.25$ Å; places I⁺ between Rb⁺/K⁺ → ladder template choice; RRK $\nu=2.42$ ps⁻¹ from $V''(R_e)$; validated vs I⁺ mobility | ~~$D_e$ not $D_0$ (ZPE)~~ **resolved 2026-06-17** ($G(0)=37$ cm⁻¹ exact); open-shell ³P₂, SO-split (six curves) → scalar collapse (A5) and likely *gradual, shell-less* ladder (R3); **no $n>1$ cluster data** |

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
    (Methods, Eqs. 1–6), $D_0^{\mathrm{Na^+}}(N)$ ladder (Extended Data Table 1),
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
