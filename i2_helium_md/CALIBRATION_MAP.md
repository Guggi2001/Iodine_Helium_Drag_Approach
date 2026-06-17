# Calibration Map — I⁺/He Drag-Port + Mass-Dynamics

> **Entry point:** `DRAG_PORT_DESIGN_DECISIONS.md` (§0 document map). This is the **cross-doc parameter index**; mass-model detail is in `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`.

One table mapping every model parameter to its **classification**, its
**Tier-1 vs Tier-2 anchor**, and its **cross-check**. Spans both
`DRAG_PORT_DESIGN_DECISIONS.md` (drag, gate, noise) and
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (pickup, evaporation,
cooling, ladder, early window). Physics-Definition only; no code.

## Classification scheme

- **Locked** — fixed by completed Tier 0; not re-fit downstream.
- **Sourced** — value taken from external literature / ab initio, carried with a
  prior (and its transfer caveat).
- **Derived** — computed from other quantities; not independently fit.
- **Bounded** — physical bounds + sweep band; reported as a range/regime, not a
  point fit.
- **Free** — genuinely determined by fitting to data (the only true knobs).

## Validation tiers (sequential, ordered by separability; DESIGN §6.4)

- **Tier 0 (COMPLETE)** — drag form + effective binding, deterministic, fixed
  mass, scored *in-window* against the `9A`/`18A` TDDFT traces (Method B
  trajectory-matching).
- **Tier 1 (next, ungated)** — mass scenario, deterministic, full post-transient
  trajectory + the time-resolved shell trajectory (~21→19→14 He). The
  sensitivity layer.
- **Tier 2** — terminal I⁺Heₙ **size distribution** (+ per-fragment velocity
  histograms) vs experiment. The sharp arbiter for mass parameters.
- **Tier 3** — ensemble second moments: final-velocity histogram **width** +
  (if applicable) VMI angular spread. Isolates the noise model.

## Flag legend

`R#` risk · `A#` assumption · `OQ#` open question (author-contact/provenance).
All in the MASS doc unless marked **[D]** = DESIGN doc.

---

## The map

| # | Domain | Parameter (units) | Role | Class | Primary anchor (Tier) | Cross-check | Flag |
|---|---|---|---|---|---|---|---|
| 1 | Drag | $b=2.5154$ (amu·ps/Å²) | cubic drag coefficient, $F=g\,b\,v^3$ | **Locked** | Tier 0 (9A/18A trajectory-match) | held-out cross-case shared-form (Tier 1); `power_law` $n̂{=}2.93$ recovers cubic | 9 Å transverse flag **[D]** |
| 2 | Drag | $a=0$ (amu/ps) | linear drag coefficient (pure cubic) | **Locked** | Tier 0 (Method B drove $a\to0$) | $\gamma_0=0$ accepted; noise null | §1.2 **[D]** |
| 3 | Drag | $E_\text{bind}=0.1168$ (eV) | effective ion binding, co-fit with drag | **Locked** (Tier-0), VMI **pending** | Tier 0 trajectory-match → Tier 2 VMI | TDDFT escape-energy sanity; integrated ladder $\sum_i D_0$ (§6.5.1) | OQ1 (electronic provenance) |
| 4 | Drag | $v_\text{ceiling}$ (Å/ps) | high-$v$ cubic ceiling cap (contingent) | **Bounded** (only if R10-(b) invoked) | TDDFT peak speed (Tier 0/1) | Tier-1 transient $v$ excursions above $v_\text{max,fit}$ | R10 |
| 5 | Gate | $g(\text{depth})$ (dimensionless) | drag spatial gate, G4→G2 | **Derived** (erf-tied G2 until $\rho_\text{He}$ profile exists) | confining-potential steepness (14.2 Å) | Tier-1 trajectory; promote to G4 with measured $\rho_\text{He}$ | §5 **[D]** |
| 6 | Noise | $T_\text{eff}$, FDT amplitude | multiplicative local-FDT bath kick (N2) | **Bounded/Derived** (tied to $\gamma(v)$; $\propto\sqrt{\gamma g k_BT_\text{eff}}$) | Tier 3 (ensemble width / VMI) | strict-FDT shown dynamically null (§1.3a) | §1 **[D]** |
| 7 | Pickup | $\lambda_\text{attach}(\rho_\text{He})$ (ps⁻¹) | Poisson He capture rate | **Sourced + Bounded** (Nat23 2.0/ps Na⁺ at rest, ±factor-2) | Tier 2 size dist | Tier-1 shell trajectory (21→19→14); 9/18 Å density contrast | R1 (Na⁺→I⁺), R7 ($v$-dep) |
| 8 | Pickup | $\rho_\text{He}(\text{depth})$ (Å⁻³) | density profile gating capture | **Sourced** (baseline/TDDFT density) | baseline | Tier 1 | — |
| 9 | Evap | $\nu$ (ps⁻¹) | RRK prefactor (I⁺–He stretch freq) | **Sourced** ($\sim\mathcal{O}(1)$, IHe05 $X_2$ curvature) | Tier 2 size dist | I2-notes cascade timing (OQ5); curvature recompute | OQ5, A11 |
| 10 | Evap | $s=3n-6$ (dimensionless) | RRK effective vibrational DOF | **Derived** (mode-count; effective-scalar override) | Tier 2 — **joint with `ladder_shape`** | size-dist tail/spread shape | A11 (classical RRK; coupling) |
| 11 | Cooling | $\tau_\text{dissip}$ (ps) | Newton-cooling time of $E_\text{solv.struct}$ | **Bounded** (sweep $[2.6,16.5]$, externally anchored) | external (GAH25 Table III / exp) + Tier 2 if sensitive | $t_\times$ vs GAH25 $t_0$; terminal-$n$ insensitivity sweep | R8 (Na⁺-only) |
| 12 | Cooling | $E_\infty=E_\text{bind}^\text{eq}$ (eV) | K2 asymptote (equilibrium-shell binding) | **Sourced** (GAH25 −3424/−4144 K Na⁺; compute I⁺ shell) | equilibrium-shell binding | I⁺ shell solvation $S_{\mathrm{I^+}}{=}{-}3578$ K [I2-notes] | OQ6 (stripping reach) |
| 13 | Budget | $f_\text{ret}\in[0,1]$ (dimensionless) | S1 pickup binding-release retained fraction | **Bounded** (prior small) | Tier 2 size dist | 9/18 Å density contrast (feedback gain → density-dependence of terminal $n$) | — |
| 14 | Early | $f_\text{int}\in[0,1]$ (dimensionless) | S2 onset partition, $E_\text{int}(0){=}f_\text{int}E_\text{avail}$ | **Bounded** (self-unbound floor → 1) | Tier 2 size dist | GAH25 $t_0$ via $t_\times$; identifiability fallback (pin at fixed $\tau$) | R1 regime axis |
| 15 | Early | $E_\text{avail}$ (eV) | Coulomb onset energy available | **Sourced/fixed ref** (~5.5 eV/pair scale) | fixed reference ($t^*$-window) | — | — |
| 16 | Early | $E_\text{int}(0)$ (eV) | initial internal energy | **Derived** ($=f_\text{int}E_\text{avail}$) | — | — | — |
| 17 | Early | $t_\times$ (ps) | self-bound crossing (gate-open time) | **Derived diagnostic** | — | GAH25 $t_0$ = 5.0/6.53 ps (±factor-2, Na⁺) | §6.11 |
| 18 | Ladder | $D_0^{\mathrm{I^+}}(1)=143.9$ cm⁻¹ ≈ 17.8 meV | first dissociation rung ($X_2/³Π$) | **Sourced** (IHe05, ~3% accuracy) | external ab initio | mobility/ZEKE-validated curve | OQ1 (which curve) |
| 19 | Ladder | $D_0^{\mathrm{I^+}}(n{>}1)$ / `ladder_shape` | rung profile: gradual vs shell-structured | **Free** (gradual default, templated) | Tier 2 size dist (gap vs smooth) | $n^*{=}21$ [I2-notes]; joint with $s$ | R3, A5, OQ4 |
| 20 | Ladder | `ladder_electronic_picture` | statistical-mixture (default) vs $X_2$-only | **Free choice** | Tier 2 size dist | $\sum_i D_0$ vs $E_\text{bind}{=}0.1168$; $S_{\mathrm{I^+}}$ | A10, OQ1 |
| 21 | Ladder | $\sum_i D_0^{\mathrm{I^+}}(i)$ (eV) | integrated ladder = self-bound gate threshold | **Derived** (from ladder + picture) | — | §6.5.1 effective binding 0.1168 eV | A10 |
| 22 | Ladder | $n^*=21$ | first-shell equilibrium occupancy | **Sourced** ([I2-notes] $X_2$-only) | external | alkali $R_e$-scaling (Rb⁺/Cs⁺ ≈ 18/21) | OQ4 |

### Cross-cutting (not single parameters)

| Item | Nature | Status | Tier / check | Flag |
|---|---|---|---|---|
| Mass↔drag coefficient pairing | consistency, not a parameter | structural: production runs `allow_inconsistent_mass_pairing` | §6.6 mid-window defence; Tier 2 sensitivity | R6 |
| Total-stripping (Calvo24) limit | reachable far end of regime axis | secondary/sensitivity-run evaluation, not default | Tier 2 size dist (terminal $n\to$ small) | OQ6 ($E_\infty$ reach) |
| Terminal regime ($f_\text{int},\tau$, ladder depth) | reported output, not a knob | regime determination (shell-retaining ↔ total strip) | Tier 2 size dist | R1, §6.11 |

---

## Tally — how few true knobs remain

- **Locked (3):** $b$, $a$, $E_\text{bind}$ (Tier-0 complete).
- **Sourced (7):** $\lambda_\text{attach}$, $\rho_\text{He}$, $\nu$, $E_\infty$,
  $E_\text{avail}$, $D_0(1)$, $n^*$.
- **Derived (5):** $g(\text{depth})$, $s$, $E_\text{int}(0)$, $t_\times$,
  $\sum_i D_0$.
- **Bounded (4–5):** $\tau_\text{dissip}$, $f_\text{ret}$, $f_\text{int}$,
  $v_\text{ceiling}$, ($T_\text{eff}$/noise).
- **Free (2):** the ladder **shape** ($D_0(n{>}1)$) and the **electronic
  picture** — both choices the Tier-2 size distribution arbiters.

**Headline:** after sourcing/derivation/bounding, the genuinely free
calibration collapses to **the ladder shape and the electronic picture**, both
discriminated by one observable (the Tier-2 size distribution), with the
early-window scalars bounded and cross-checked rather than fit. Everything else
is locked, imported with a prior, or computed.

## Anchor coverage by tier

- **Tier 0:** drag form + effective binding (done).
- **Tier 1:** drag/mass-scenario sensitivity via the shell trajectory (21→19→14)
  — touches $\lambda_\text{attach}$, $g$, $v_\text{ceiling}$, the early-window
  tolerance.
- **Tier 2 (load-bearing):** $\lambda_\text{attach}$, $f_\text{ret}$,
  $f_\text{int}$, $\nu$, $s$, ladder shape, electronic picture, terminal regime,
  $\tau$ (if sensitive). **Single observable carrying most of the calibration**
  — its discriminating power is the project's central dependency.
- **Tier 3:** the noise model.
- **External / author-contact:** OQ1 (drag electronic state), OQ2
  ($KE_\text{shed}$), OQ3 ($E_\text{bind}(N)$ vs ladder), OQ4 ($S_{\mathrm{I^+}}$/$n^*$
  provenance), OQ5 (cascade timing), OQ6 ($E_\infty$ stripping reach).

*Note:* Tier 2 carries a heavy load (8+ quantities on one observable). The
identifiability arguments are documented per-parameter (e.g. $f_\text{ret}$ via
9/18 Å contrast, $\tau$ via post-crossing tail, $s$ jointly with ladder shape) —
but whether the size distribution actually separates all of them is itself the
key open empirical question for the calibration campaign.
