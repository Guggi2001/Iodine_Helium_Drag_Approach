# Calibration Map — I⁺/He Drag-Port + Mass-Dynamics

> **Entry point:** `DRAG_PORT_DESIGN_DECISIONS.md` (§0 document map). This is the **cross-doc parameter index**; mass-model detail is in `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`.

One table mapping every model parameter to its **classification**, its
**Tier-1 vs Tier-2 anchor**, and its **cross-check**. Spans both
`DRAG_PORT_DESIGN_DECISIONS.md` (drag, gate, noise) and
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (pickup, evaporation,
cooling, ladder, early window). Physics-Definition only; no code.

> **Update 2026-06-17.** [IHe05] EPAPS analytic fit obtained → rows 9, 18, 20
> tightened from order-of-magnitude prior to **pinned**: $D_0^{\mathrm{I^+}}(1)$
> ($X_2$ 106.9 / mixture 74.4 cm⁻¹, exact $J{=}0$ ZPE) and RRK $\nu=2.42$ ps⁻¹
> ($V''(R_e){=}748.1$ cm⁻¹/Å²). Detail + the SO-coupled $X_2/I_1/I_0$ depths in the
> MASS doc 2026-06-17 revision. The two **Free** knobs (ladder shape, electronic
> picture) are unchanged in count; the electronic-picture *rung scale* is now
> numeric (mixture/$X_2=0.70$).
>
> **Update 2026-06-17 (cont.).** Ladder shape **Form U** adopted: the discrete
> `ladder_shape ∈ {gradual, shell_structured}` (old row 19) becomes a single
> sigmoid with one continuous knob $\kappa$ (`ladder_steepness`), anchored at the
> pinned rung, a sourced bulk-He floor (new row 23), and $n^*$ (row 22). The Free
> count stays 2 ($\kappa$ + picture, co-fit). Cross-check corrected (row 21):
> $\sum_i D_0 \neq |S_{\mathrm{I^+}}|$ — $|S|$ is collective (pair ladder cannot
> reach 2484 cm⁻¹), an upper bound only; reachable target is the drag binding
> 0.1168 eV. New OQ7 (the $|S|$ energy reference). Reversible.
>
> **Update 2026-06-17 (cont. 3).** Full [GAH25] studied → confirms K2; **$E_\infty$
> split locked**: occupancy-resolved $E_\infty(N)=-|S(N)|$ with explicit
> electrostriction term (row 12, new 12b — the *dominant* binding per GAH25
> geometry), resolving OQ6 (stripping cap) and the equilibrium layer of R12.
> Ladder geometry corrected from GAH25 Table II: shell radius ~4.67 Å (not pair
> $R_e$) → cliff **7.5×** (row 19, was 13.4×), He–He roomy not compressed.
> Secondary: $n^*\approx20$ cross-check (row 22, OQ8); $\lambda_\text{attach}$
> central ~0.7–1.1/ps for I⁺ (row 7); Calvo K⁺ PIMC as a $\kappa$ lead. Reversible.
>
> **Update 2026-06-17 (cont. 2).** Early-window scalars pinned: $E_\text{avail}^
> \text{ion}=2.70$ eV (row 15, per-ion convention; ½ of 5.40 eV pair, reversible);
> integrated ladder $\sum_i D_0$ now numeric and **nearly $\kappa$-independent**
> (row 21); $f_\text{int}$ self-unbound **floor $\approx0.04$–$0.10$** (row 14),
> small ⇒ the GAH25 self-unbound onset is robust, not fine-tuned. Form U cliff
> recentered to $n^*+\tfrac12$. The gate threshold and floor are **picture-set,
> $\kappa$-independent** — pinnable ahead of the Tier-2 $\kappa$/picture fit.

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
| 7 | Pickup | $\lambda_\text{attach}(\rho_\text{He})$ (ps⁻¹) | Poisson He capture rate | **Sourced + Bounded** — central ~0.7–1.1/ps (GAH25 Rb⁺/Cs⁺; I⁺ is Rb⁺-like), **not** 2.0 (Na⁺ exp); ±factor-2 | Tier 2 size dist | Tier-1 shell trajectory; 9/18 Å density contrast | R1 (Na⁺→I⁺), R7 ($v$-dep) |
| 8 | Pickup | $\rho_\text{He}(\text{depth})$ (Å⁻³) | density profile gating capture | **Sourced** (baseline/TDDFT density) | baseline | Tier 1 | — |
| 9 | Evap | $\nu$ (ps⁻¹) | RRK prefactor (I⁺–He stretch freq) | **Sourced (pinned 2.42, 2026-06-17)** — $\omega_e{=}80.6$ cm⁻¹ from IHe05 EPAPS $V''(R_e){=}748.1$ cm⁻¹/Å² | Tier 2 size dist | I2-notes cascade timing (OQ5); near-threshold spacing → effective $s$ | OQ5, A11 |
| 10 | Evap | $s=3n-6$ (dimensionless) | RRK effective vibrational DOF | **Derived** (mode-count; effective-scalar override) | Tier 2 — **joint with `ladder_steepness` ($\kappa$)** | size-dist tail/spread shape | A11 (classical RRK; coupling) |
| 11 | Cooling | $\tau_\text{dissip}$ (ps) | Newton-cooling time of $E_\text{solv.struct}$ | **Bounded** (sweep $[2.6,16.5]$, externally anchored) | external (GAH25 Table III / exp) + Tier 2 if sensitive | $t_\times$ vs GAH25 $t_0$; terminal-$n$ insensitivity sweep | R8 (Na⁺-only) |
| 12 | Cooling | $E_\infty(N)=-\|S(N)\|$ (eV) | K2 asymptote, **occupancy-resolved** | **Sourced (split, 2026-06-17)** — full-shell $\|S_{\mathrm{I^+}}\|{=}0.308$ eV [I2-notes]; GAH25 −3424/−4144 K Na⁺ fixes form+$\tau$ | equilibrium-shell binding | OQ6 resolved (fixed $E_\infty$ caps strip); shape via $\kappa$ | OQ6, OQ7 |
| 12b | Cooling | $E_\text{elec}(N)$ (eV) | electrostriction binding (collective excess) | **Sourced/Derived** ($-(\|S(N)\|-\sum_i D_0)$; **dominant** per GAH25 geometry) | — | marginal release → bath on shed (A8); $\partial\|S\|/\partial n{\approx}124$ cm⁻¹ ≈ $D_0(1)$ | K2, A8, OQ7 |
| 13 | Budget | $f_\text{ret}\in[0,1]$ (dimensionless) | S1 pickup binding-release retained fraction | **Bounded** (prior small) | Tier 2 size dist | 9/18 Å density contrast (feedback gain → density-dependence of terminal $n$) | — |
| 14 | Early | $f_\text{int}\in[0,1]$ (dimensionless) | S2 onset partition, $E_\text{int}(0){=}f_\text{int}E_\text{avail}^\text{ion}$ | **Bounded** — floor $\approx0.04$–$0.10$ ($X_2$), ~0.06–0.07 (mix), picture-set & $\kappa$-indep (2026-06-17); soft upper ~0.2 | Tier 2 size dist | GAH25 $t_0$ via $t_\times$; floor robust → onset not fine-tuned (A7) | R1 regime axis |
| 15 | Early | $E_\text{avail}^\text{ion}$ (eV) | per-ion Coulomb onset budget | **Sourced/fixed (pinned 2.70, 2026-06-17)** — ½·$e^2/R_e(\mathrm{I_2})$, $R_e{=}2.666$ Å | fixed reference ($t^*$-window) | per-ion convention (reversible to 5.40 pair) | — |
| 16 | Early | $E_\text{int}(0)$ (eV) | initial internal energy | **Derived** ($=f_\text{int}E_\text{avail}^\text{ion}$) | — | — | — |
| 17 | Early | $t_\times$ (ps) | self-bound crossing (gate-open time) | **Derived diagnostic** | — | GAH25 $t_0$ = 5.0/6.53 ps (±factor-2, Na⁺) | §6.11 |
| 18 | Ladder | $D_0^{\mathrm{I^+}}(1)$: $X_2$ 106.9 / mix 74.4 cm⁻¹ | first dissociation rung (picture-dependent) | **Sourced (pinned 2026-06-17)** — IHe05 EPAPS fit, exact $J{=}0$ ZPE $G(0){=}37$ cm⁻¹ ($D_e{=}143.9$ is **not** $D_0$); $\pm3$ cm⁻¹ | external ab initio | mobility/ZEKE-validated curve; mixture/$X_2$=0.70 | OQ1 (which curve) |
| 19 | Ladder | $\kappa$ (Form U) / `ladder_steepness` | sigmoid steepness, gradual↔cliff (1 knob) | **Free** (prior large: 7.5× radial cliff) | Tier 2 size dist (broad vs magic-peak) | $n^*{=}21$ [I2-notes]; joint w/ picture + $s$; tabulated fallback | R3, A5, OQ4 |
| 20 | Ladder | `ladder_electronic_picture` | statistical-mixture (default) vs $X_2$-only | **Free choice** | Tier 2 size dist | rung scale pinned: mix 74.4 vs $X_2$ 106.9 cm⁻¹ (=0.70); $\sum_i D_0$ vs $E_\text{bind}{=}0.1168$; $S_{\mathrm{I^+}}$ | A10, OQ1 |
| 21 | Ladder | $\sum_i D_0^{\mathrm{I^+}}(i)$ (eV) | integrated ladder = self-bound gate threshold | **Derived** — $X_2$ 0.12–0.28 / mix 0.17–0.19 eV (2026-06-17); **nearly $\kappa$-indep** (cliff at $n^*{+}\tfrac12$) | — | drag binding 0.1168 eV (reachable); $|S|{=}0.308$ is **collective UB, NOT a sum target** | A10, OQ7 |
| 22 | Ladder | $n^*=21$ | first-shell equilibrium occupancy | **Sourced** ([I2-notes] $X_2$-only) | external | **GAH25 $R_e$-scaling → ~20** (Rb⁺ 18/Cs⁺ 21; corroborates 21 to ±1–2); OQ8 | OQ4, OQ8 |
| 23 | Ladder | $D_\text{floor}$ (cm⁻¹) | outer-shell rung floor (Form U) | **Sourced** ($\|\mu_\text{He}^\text{bulk}\|\approx4.97$ cm⁻¹, 7.15 K) | external (bulk superfluid) | picture-independent; sets sigmoid lower anchor | R3 |

### Cross-cutting (not single parameters)

| Item | Nature | Status | Tier / check | Flag |
|---|---|---|---|---|
| Mass↔drag coefficient pairing | consistency, not a parameter | structural: production runs `allow_inconsistent_mass_pairing` | §6.6 mid-window defence; Tier 2 sensitivity | R6 |
| Total-stripping (Calvo24) limit | reachable far end of regime axis | secondary/sensitivity-run evaluation, not default | Tier 2 size dist (terminal $n\to$ small) | OQ6 ($E_\infty$ reach) |
| Terminal regime ($f_\text{int},\tau$, ladder depth) | reported output, not a knob | regime determination (shell-retaining ↔ total strip) | Tier 2 size dist | R1, §6.11 |

---

## Tally — how few true knobs remain

- **Locked (3):** $b$, $a$, $E_\text{bind}$ (Tier-0 complete).
- **Sourced (8):** $\lambda_\text{attach}$, $\rho_\text{He}$, $\nu$, $E_\infty$,
  $E_\text{avail}^\text{ion}$, $D_0(1)$, $n^*$, $D_\text{floor}$.
- **Derived (5):** $g(\text{depth})$, $s$, $E_\text{int}(0)$, $t_\times$,
  $\sum_i D_0$.
- **Bounded (4–5):** $\tau_\text{dissip}$, $f_\text{ret}$, $f_\text{int}$,
  $v_\text{ceiling}$, ($T_\text{eff}$/noise).
- **Free (2):** the ladder **steepness** $\kappa$ (Form U sigmoid, $D_0(n{>}1)$)
  and the **electronic picture** — both choices the Tier-2 size distribution
  arbiters, co-fit (not separable).

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
  provenance), OQ5 (cascade timing), OQ6 ($E_\infty$ stripping reach), OQ7
  ($S_{\mathrm{I^+}}$ energy reference — interprets the $|S|$-vs-$\sum D_0$ gap).

*Note:* Tier 2 carries a heavy load (8+ quantities on one observable). The
identifiability arguments are documented per-parameter (e.g. $f_\text{ret}$ via
9/18 Å contrast, $\tau$ via post-crossing tail, $s$ jointly with ladder shape) —
but whether the size distribution actually separates all of them is itself the
key open empirical question for the calibration campaign.
