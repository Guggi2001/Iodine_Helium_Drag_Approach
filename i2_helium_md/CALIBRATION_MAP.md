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
>
> **Update 2026-06-21 (consistency-check pass).** Four fixes, no lock/mechanism
> change (MASS doc 2026-06-21 revision): (1) **$n^*$ leak corrected** — per-atom
> collective binding adopts the cation $n^*{=}21$ throughout ($|S|/n^*{=}118$ cm⁻¹
> / 170 K, was 124/179 off $n^*{=}20$); ratio-to-pair $4.7\times$ (rows 12b, 22);
> **OQ8 raised LOW→LOW–MEDIUM** (the ambiguity had leaked into energetics, not
> geometry only). (2) **Integrated $X_2$ band corrected to 0.25–0.28 eV** (Form U
> pure-$\kappa$ range); the crowding-reduced value and the drag binding 0.117 eV are
> **separate cross-checks, not band ends** (row 21). (3) **`ladder_electronic_picture`
> gains `cooling_relaxed`** (between mix & $X_2$; row 20) — still one selection, Free
> count unchanged. (4) **$E_\text{avail}^\text{ion}$ scenario-keyed** — 0.80 eV
> validation ($d{=}9$ Å) / 2.70 eV production ($R_e$); the $f_\text{int}$ floor
> triples at the validation budget (0.065→0.22 mix), self-unbound onset still robust
> but margin $14\times{\to}4.5\times$ (rows 14, 15; A7 reworded). The soft ceiling
> ~0.2 is **advisory, not a constraint**. Tier-1 endpoint 14 author-confirmed
> I⁺-specific. **Item 4 provisional pending OQ2** ($KE_\text{shed}$/partition).
> Reversible.
>
> **Update 2026-06-21 (cont.) — Method-level consistency proof.** Audited the
> mechanism as a closed system (MASS doc 2026-06-21 cont. revision): (1) the
> **five-term energy invariant verified closed** under drag, pickup S1, cold-shed
> K1, and K2 cooling — a genuine conservation law, no change. (2) **RRK $s$
> corrected $3n-6\to3n-3$** (full $n{+}1$-atom complex; row 10): the old count went
> $\le0$ for $n\le2$, diverging the rate at threshold and breaking the
> no-avalanche guarantee in the small-$n$ tail — the only regime where $\{\nu,s\}$
> are observable. $n{=}1$ is now **direct dissociation $k=\nu$**; override guarded
> $s\ge1$. Mechanism/invariant unchanged.
>
> **Update 2026-06-21 (cont. 2) — pickup↔gate stability + occupancy cap.**
> Extended the Method proof to the accretion↔gate feedback loop (MASS doc cont. 2
> revision, §6.11 stability note): (1) **loop proven stable** — the self-unbound
> margin $G=E_\text{int}-\Sigma(n)$ is a pathwise Lyapunov function (gate always
> crosses; pickup is stabilizing), and terminal $n$ is a stable freeze-out
> attractor set by $\Pi=\lambda f_\text{ret}\tau$ ($>1$ shed / $<1$ freeze;
> $\Pi\to0$ at exit ⇒ termination universal). "Equilibrium is emergent, not a
> parameter" is now **proven**. $\Pi$ added as a derived diagnostic (row 17b;
> Derived 5→6). (2) **Occupancy cap added** (Langmuir Form B, A12):
> $\lambda_\text{attach}{\times}(1-n/n^*)_+^{\,p}$ closes the resting-ion
> $n\to\infty$ gap; **inert for production** (ion exits first), rate-only so the
> invariant and $G$-crossing are untouched. New exponent $p$ (row 7p) **default
> tied to $\kappa$** → 0 net new free knobs. Rows 7, 17b, 7p, tally updated.
> Reversible.
>
> **Update 2026-06-21 (cont. 3) — integrator↔mass-jump operator split (A13).**
> Final structural item (MASS doc cont. 3 revision): **SQ1** the velocity-dependent
> drag freezes $\gamma(v_\text{in})$ → drag-on path is **$O(dt)$, not $O(dt^2)$**
> (retire the BAOAB 2nd-order claim for production); accepted because it buys exact
> dissipation bookkeeping + unconditional dissipativity, at the cost of a one-signed
> over-braking bias (R10). **SQ2** mass jumps **must** use the momentum-conserving
> reset $v^+=(m v^-\pm m_\text{He}u_\text{He})/(m\pm m_\text{He})$ — an **invariant
> precondition**: $E_\text{mass\_transfer}$ is the reduced-mass defect
> $\tfrac12\tfrac{m m_\text{He}}{m+m_\text{He}}\|v^-{-}u_\text{He}\|^2$ it produces
> (corrected from the heavy-ion $\tfrac12 m_\text{He}v^2$, ~3% closure error at
> $n{=}1$); He at rest. **SQ3** post-jump O-step uses $m^+$. SQ2/SQ3 unbuilt, now
> specified. **Method proof closes conditional on the SQ2 reset.** Cross-cutting row
> added. Reversible.
>
> **Update 2026-07-06 — RRK-dof promotion (staircase-probe outcome).** The
> pre-F5 staircase capability probe + s_eff mini-probe
> (`TIER2_STAIRCASE_PROBE_PLAN.md` Addendum A; execution records in
> `docs/drag_port/Tier2/drag_migration_log_tier2.md`) fired and resolved the
> RRK-dof mechanism-level OQ: at the classical mode count $s=3n-3$ (=60 at
> $n^*{=}21$) the cascade is kinetically frozen at $n\approx20$ everywhere in
> the κ×picture×τ bands (max 1.7 sheds vs the anchored 7), while a **constant
> $s_\text{eff}=8$ at mid-band τ lands the anchored 21→19→14 staircase in
> magnitude and timing** (7.41 sheds → n\_end 13.59, first shed 5.42 ps vs
> t★=5 ps, trajectory MAD 1.0 He), **picture-robust** (≤4 % magnitude spread
> across all three pictures, cross-check 2026-07-06). **User decision
> (2026-07-06): $s$ promoted** — row 10 reclassified **Derived → Bounded**
> (constant effective $s_\text{eff}$, band ≈ [5, 20], landing [8, 12] at
> mid-band τ); the classical $s=3n-3$ is demoted to the classical-limit arm
> of a dof-convention selection (retained as diagnostic/limit; the $s\ge1$
> guard and $n{=}1$ direct dissociation unchanged). Consequences: the A11
> classical band $[3n{-}3,3n]$ is **superseded as the sweep band** (the
> landing sits far *below* it — the quantum/blurred-shell down-drift A11
> itself flagged); κ (row 19) is empirically **inverted and
> normalisation-capped** for the stripping range (Form U floors $D_0(21)$ at
> ≈0.53·$D_0(1)$); the picture (row 20) is **near-degenerate in shed
> magnitude** ($x=D_0/\Sigma$ picture-invariant), sensitive only through
> $\Sigma(21)\to t_\times$ timing, degenerate there with $f_\text{int}$.
> Free count unchanged (κ + picture — both now known near-flat on the
> staircase observable; the size distribution remains their arbiter);
> Derived 6→5, Bounded +1. s↔τ are separated by first-shed timing, so the
> 9 Å staircase (a TDDFT prior, not ground truth) co-anchors the pair —
> the first load shed from the one-observable Tier-2 stack. Implementation
> (selection surface + F2 campaign re-scope to an $s_\text{eff}$×τ co-fit)
> stays behind the implementation trigger. Reversible.
>
> **Update 2026-07-06 — cooling spatial gate arm (`cooling_spatial_gate`).** K2
> Newton cooling ("bath dissipation") was applied **ungated everywhere**, unlike
> drag ($\gamma\propto\rho_\text{He}$) and pickup ($\lambda\propto\rho_\text{He}$),
> which gate off outside the bubble via the shared erf-complement surface. Added a
> first-class interchangeable enum arm (**new row 5a**): `none` (default,
> **byte-identical** to the locked ungated form) / `density_scaled`
> ($\tau_\text{eff}=\tau/\rho_\text{ratio}$, reusing the **row-5** drag gate
> steepness — one bubble boundary for all three He-mediated channels; cooling
> switches off in vacuum). Physically motivated: the GAH25 $\tau$ was fit near the
> droplet, but an ejected I⁺Heₙ cluster has no external bath. This is the
> **ejection-regime total-vaporization lever** — ungated cooling self-quenches the
> cascade at $n\sim$ few (shell-retaining), the gate opens the near-bare end the
> experimental 43 % bare-I⁺ peak requires. Probed by the total-strip A/B capability
> grid (`gen_tier2_staircase_probe.py`, **probe-scoped**; not yet a campaign knob).
> Clean erf-cutoff is the first build; a residual out-of-bubble floor
> $\rho_\text{min}$ is a **deferred OQ**. See MASS §6 K2 / §6.11,
> DRAG_PORT_DESIGN_DECISIONS §"Cooling spatial gate",
> `drag_migration_log_tier2.md`. Reversible.
>
> **Update 2026-07-07 — detection time sourced (terminal-read design).** Wave 5
> (I11) showed the gated terminal read is **flight-time dependent** (the
> post-ejection cascade is live at any tractable fixed-dt cap), so the Tier-2
> observable needs an explicit detection time. Design delivered
> (`docs/drag_port/Tier2/TIER2_DETECTION_STAGE_DESIGN.md`): an event-driven
> post-ejection continuation reads the per-ion size distribution at
> $t_\text{detect}$ exactly (the mechanism reduces to a Markov jump chain
> after ejection — a scheme change, not new physics). **New row 24:**
> $t_\text{detect}=8.53$ µs $=8.53{\cdot}10^6$ ps, **Sourced** from the
> experimental-setup publication (TOF flight to detector; user-confirmed
> 2026-07-07). Sourced tally 8→9. Detection-time *sensitivity* is a free
> report-side read of the stored shed events (design §3.3). Build stays
> behind the implementation trigger. Reversible.
>
> **NB (2026-07-07, Slice DS delivered):** the carrying config field
> `detection_time_ps` is now live (`check_detection_config`;
> required-when-enabled, no baked-in default — the constant stays
> calibration data supplied per run). Row-24 value unchanged.

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
  *(NB 2026-07-01, Tier-2 Phase-E scope: the delivered observable is the
  **integer-n size distribution only** — per-fragment velocity histograms were
  **cut** (the `vmi_summary/*.csv` references are aggregate I⁺He / I⁺-gas
  distributions, not per-n). See
  `docs/drag_port/Tier2/TIER2_PHASE_E_IMPLEMENTATION_PLAN.md`.)*
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
| 5a | Gate | `cooling_spatial_gate` (arm) | K2 cooling spatial gate: `none` (locked, ungated $-E_\text{int}/\tau$) / `density_scaled` ($\tau_\text{eff}{=}\tau/\rho_\text{ratio}$, off outside bubble) | **Free choice** (2 arms, 1 selection; default `none` **byte-identical**; added 2026-07-06) | Tier 2 size dist (total-strip A/B capability probe) | reuses **row-5** drag gate surface (`drag_gate_steepness`, one bubble boundary); ejection-regime total-vaporization lever (shell-retaining ↔ near-bare); $\rho_\text{min}$ floor **deferred OQ**; probe-scoped, not yet a campaign knob | MASS §6 K2/§6.11 |
| 6 | Noise | $T_\text{eff}$, FDT amplitude | multiplicative local-FDT bath kick (N2) | **Bounded/Derived** (tied to $\gamma(v)$; $\propto\sqrt{\gamma g k_BT_\text{eff}}$) | Tier 3 (ensemble width / VMI) | strict-FDT shown dynamically null (§1.3a) | §1 **[D]** |
| 7 | Pickup | $\lambda_\text{attach}(\rho_\text{He},n)$ (ps⁻¹) | Poisson He capture rate (density-gated, occupancy-capped) | **Sourced + Bounded** — central ~0.7–1.1/ps (GAH25 Rb⁺/Cs⁺; I⁺ is Rb⁺-like), **not** 2.0 (Na⁺ exp); ±factor-2; ×Langmuir cap $(1-n/n^*)_+^{\,p}$ (A12, 2026-06-21) | Tier 2 size dist | Tier-1 shell trajectory; 9/18 Å density contrast; cap inert for ejection | R1 (Na⁺→I⁺), R7 ($v$-dep), A12 |
| 7p | Pickup | $p\ge0$ | occupancy-cap sharpness (Langmuir exponent, A12) | **Free (conditional)** — Phase-B default **$p=1$ fixed** (NB 2026-07-01: the $p\!\leftrightarrow\!\kappa$ tie is **inverse** — rigid shell = large $\kappa$ = smaller $p$; **NOT** $p=\kappa$; match cutoff slopes if tied); split only if size dist. demands | Tier 2 size dist | [Nat23] resting-ion saturation bounds $n^*,p$; $p\!\leftrightarrow\!\kappa$ coupling (inverse); Phase-B §3.2 | A12, OQ4/OQ8 |
| 8 | Pickup | $\rho_\text{He}(\text{depth})$ (Å⁻³) | density profile gating capture | **Sourced** (baseline/TDDFT density) | baseline | Tier 1 | — |
| 9 | Evap | $\nu$ (ps⁻¹) | RRK prefactor (I⁺–He stretch freq) | **Sourced (pinned 2.42, 2026-06-17)** — $\omega_e{=}80.6$ cm⁻¹ from IHe05 EPAPS $V''(R_e){=}748.1$ cm⁻¹/Å² | Tier 2 size dist | I2-notes cascade timing (OQ5); near-threshold spacing → effective $s$ | OQ5, A11 |
| 10 | Evap | $s_\text{eff}$ (dimensionless) | RRK effective vibrational DOF | **Bounded (promoted 2026-07-06; was Derived)** — constant $s_\text{eff}$, band **≈[5, 20]** (staircase landing [8, 12] at mid-band τ; picture-robust ≤4 %); classical mode-count $s{=}3n{-}3$ (corr. 2026-06-21 from $3n-6$) demoted to the classical-limit arm of the dof-convention selection; $n{=}1$ direct dissociation $k=\nu$ and the $s\ge1$ guard unchanged | 9 Å anchored staircase (magnitude; s↔τ separated by first-shed timing; TDDFT prior) + Tier 2 size dist | classical band $[3n{-}3,3n]$ superseded (landing far below: quantum mode-freezing / blurred shell — the A11 down-drift); κ-joint resolved **weak** (κ near-flat on the staircase); small-$n$ tail still the size-dist regime | A11 (resolution NB 2026-07-06) |
| 11 | Cooling | $\tau_\text{dissip}$ (ps) | Newton-cooling time of $E_\text{solv.struct}$ | **Bounded** (sweep $[2.6,16.5]$, externally anchored) | external (GAH25 Table III / exp) + Tier 2 if sensitive | $t_\times$ vs GAH25 $t_0$; terminal-$n$ insensitivity sweep | R8 (Na⁺-only) |
| 12 | Cooling | $E_\infty(N)=-\|S(N)\|$ (eV) | K2 asymptote, **occupancy-resolved** | **Sourced (split, 2026-06-17)** — full-shell $\|S_{\mathrm{I^+}}\|{=}0.308$ eV [I2-notes]; GAH25 −3424/−4144 K Na⁺ fixes form+$\tau$ | equilibrium-shell binding | OQ6 resolved (fixed $E_\infty$ caps strip); shape via $\kappa$ | OQ6, OQ7 |
| 12b | Cooling | $E_\text{elec}(N)$ (eV) | electrostriction binding (collective excess) | **Sourced/Derived** ($-(\|S(N)\|-\sum_i D_0)$; **dominant** per GAH25 geometry) | — | marginal release → bath on shed (A8); $\partial\|S\|/\partial n{\approx}118$ cm⁻¹ ($n^*{=}21$, corr. 2026-06-21; was 124 off $n^*{=}20$) ≈ $D_0(1)$ | K2, A8, OQ7, OQ8 |
| 13 | Budget | $f_\text{ret}\in[0,1]$ (dimensionless) | S1 pickup binding-release retained fraction | **Bounded** (prior small) | Tier 2 size dist | 9/18 Å density contrast (feedback gain → density-dependence of terminal $n$) | — |
| 14 | Early | $f_\text{int}\in[0,1]$ (dimensionless) | S2 onset partition, $E_\text{int}(0){=}f_\text{int}E_\text{avail}^\text{ion}$ | **Bounded** — floor scenario-keyed (2026-06-21): **0.065 (mix)/0.09–0.10 ($X_2$) @ 2.70 eV; 0.21–0.24/0.31–0.35 @ 0.80 eV** ($d{=}9$ Å), picture-set & $\kappa$-indep; soft upper ~0.2 is **advisory, not a constraint** | Tier 2 size dist | GAH25 $t_0$ via $t_\times$; floor robust both scenarios, margin $14\times{\to}4.5\times$ (A7) | R1 regime axis, OQ2 |
| 15 | Early | $E_\text{avail}^\text{ion}$ (eV) | per-ion Coulomb onset budget | **Sourced/fixed, scenario-keyed (2026-06-21)** — **0.80 eV** validation (½·14.40/9, $d{=}9$ Å) / **2.70 eV** production (½·14.40/2.666, $R_e$) | fixed reference ($t^*$-window) | stamped to scenario guard (§6.5); reversible to pair (1.60/5.40) | OQ-none |
| 16 | Early | $E_\text{int}(0)$ (eV) | initial internal energy | **Derived** ($=f_\text{int}E_\text{avail}^\text{ion}$) | — | — | — |
| 17 | Early | $t_\times$ (ps) | self-bound crossing (gate-open time) | **Derived diagnostic** | — | GAH25 $t_0$ = 5.0/6.53 ps (±factor-2, Na⁺) | §6.11 |
| 17b | Early | $\Pi(t)=\lambda(n)f_\text{ret}\tau$ (dimensionless) | pickup↔gate order parameter (NEW 2026-06-21) | **Derived diagnostic** (§6.11) | — | $\Pi{>}1$ shed / $\Pi{<}1$ freeze; regime-axis spine; $\Pi\to0$ at exit ⇒ termination guaranteed; reconstructable post-hoc | §6.11, R1 |
| 18 | Ladder | $D_0^{\mathrm{I^+}}(1)$: $X_2$ 106.9 / mix 74.4 cm⁻¹ | first dissociation rung (picture-dependent) | **Sourced (pinned 2026-06-17)** — IHe05 EPAPS fit, exact $J{=}0$ ZPE $G(0){=}37$ cm⁻¹ ($D_e{=}143.9$ is **not** $D_0$); $\pm3$ cm⁻¹ | external ab initio | mobility/ZEKE-validated curve; mixture/$X_2$=0.70 | OQ1 (which curve) |
| 19 | Ladder | $\kappa$ (Form U) / `ladder_steepness` | sigmoid steepness, gradual↔cliff (1 knob) | **Free** (prior large: 7.5× radial cliff) | Tier 2 size dist (broad vs magic-peak) | $n^*{=}21$ [I2-notes]; joint w/ picture + $s$; tabulated fallback; **staircase probe 2026-07-06: inverted + normalisation-capped for 21→14** (cliff centre $n^*{+}\tfrac12$ above the stripping range; Form U floors $D_0(21)$ ≈ 0.53·$D_0(1)$) — near-flat on the staircase, size dist remains its arbiter | R3, A5, OQ4 |
| 20 | Ladder | `ladder_electronic_picture` | statistical-mixture (default) / $X_2$-only / cooling_relaxed | **Free choice** (3 options, 1 selection; `cooling_relaxed` added 2026-06-21, between mix & $X_2$) | Tier 2 size dist | rung scale: mix 74.4 / $X_2$ 106.9 cm⁻¹ (=0.70), `cooling_relaxed` in between; $\sum_i D_0$ vs $E_\text{bind}{=}0.1168$; $S_{\mathrm{I^+}}$; **cross-check 2026-07-06: near-degenerate in shed magnitude** ($x=D_0/\Sigma$ picture-invariant, ≤4 % under reduced $s_\text{eff}$); timing-only sensitivity via $\Sigma(21)\to t_\times$, degenerate there with $f_\text{int}$ | A10, OQ1 |
| 21 | Ladder | $\sum_i D_0^{\mathrm{I^+}}(i)$ (eV) | integrated ladder = self-bound gate threshold | **Derived** — $X_2$ **0.25–0.28** / mix 0.17–0.19 eV (Form U pure-$\kappa$ range, corr. 2026-06-21; ~11%, **nearly $\kappa$-indep**, cliff at $n^*{+}\tfrac12$) | — | drag binding 0.1168 eV + crowding-reduced value are **separate cross-checks, not band ends**; $\|S\|{=}0.308$ collective UB | A10, OQ7 |
| 22 | Ladder | $n^*=21$ | first-shell equilibrium occupancy | **Sourced** ([I2-notes] $X_2$-only; **adopt 21 for all per-atom energetics**, 2026-06-21) | external | GAH25 $R_e$-scaling → ~20 (corroborates to ±1–2); Tier-1 endpoint 14 **author-confirmed I⁺-specific** (2026-06-21); OQ8 raised LOW→LOW–MEDIUM (had leaked into energetics) | OQ4, OQ8 |
| 23 | Ladder | $D_\text{floor}$ (cm⁻¹) | outer-shell rung floor (Form U) | **Sourced** ($\|\mu_\text{He}^\text{bulk}\|\approx4.97$ cm⁻¹, 7.15 K) | external (bulk superfluid) | picture-independent; sets sigmoid lower anchor | R3 |
| 24 | Detection | $t_\text{detect}=8.53$ µs ($8.53{\cdot}10^6$ ps) | terminal-read detection time (TOF flight to detector) — where the Tier-2 size distribution is read (`detection_time_ps`) | **Sourced (2026-07-07)** — experimental-setup publication, user-confirmed | Tier 2 size dist (the *detected* read is the arbitration observable) | log-sensitivity band = free report-side re-read of stored shed events (detection-stage design §3.3); gated reads at $s_\text{eff}\gtrsim12$ are undefined without it (I11) | — |

### Cross-cutting (not single parameters)

| Item | Nature | Status | Tier / check | Flag |
|---|---|---|---|---|
| Mass↔drag coefficient pairing | consistency, not a parameter | structural: production runs `allow_inconsistent_mass_pairing` | §6.6 mid-window defence; Tier 2 sensitivity | R6 |
| Total-stripping (Calvo24) limit | reachable far end of regime axis | secondary/sensitivity-run evaluation, not default | Tier 2 size dist (terminal $n\to$ small) | OQ6 ($E_\infty$ reach) |
| Terminal regime ($f_\text{int},\tau$, ladder depth) | reported output, not a knob | regime determination (shell-retaining ↔ total strip) | Tier 2 size dist | R1, §6.11 |
| Integrator↔mass-jump operator split | accuracy + conservation, not a parameter | **SQ1** drag-on path $O(dt)$ (frozen $\gamma$, exact bookkeeping + dissipativity); **SQ2** momentum reset = invariant precondition (unbuilt); **SQ3** post-jump $m^+$ in O-step (unbuilt) | implementation precondition; $dt$-convergence (over-brake bias sign) | A13 |

---

## Tally — how few true knobs remain

- **Locked (3):** $b$, $a$, $E_\text{bind}$ (Tier-0 complete).
- **Sourced (9):** $\lambda_\text{attach}$, $\rho_\text{He}$, $\nu$, $E_\infty$,
  $E_\text{avail}^\text{ion}$, $D_0(1)$, $n^*$, $D_\text{floor}$,
  $t_\text{detect}$ (8.53 µs, added 2026-07-07).
- **Derived (5):** $g(\text{depth})$, $E_\text{int}(0)$, $t_\times$,
  $\sum_i D_0$, $\Pi$ (pickup↔gate order parameter, §6.11). *(Was 6 — $s$
  promoted to Bounded 2026-07-06.)*
- **Bounded (5–6):** $\tau_\text{dissip}$, $f_\text{ret}$, $f_\text{int}$,
  $s_\text{eff}$ (promoted 2026-07-06, band ≈[5, 20], staircase landing
  [8, 12]), $v_\text{ceiling}$, ($T_\text{eff}$/noise).
- **Free (2):** the ladder **steepness** $\kappa$ (Form U sigmoid, $D_0(n{>}1)$)
  and the **electronic picture** (now 3 options — mix/$X_2$/cooling_relaxed — but
  one selection) — both choices the Tier-2 size distribution
  arbiters, co-fit (not separable). *(2026-07-06: both shown **near-flat on the
  9 Å staircase** — κ inverted/capped, picture magnitude-degenerate — so the
  size distribution is now their only realistic anchor; the Stage-1 campaign
  co-fit re-scopes to $s_\text{eff}$×τ, Phase F plan F2 NB.)* *Conditional:*
  the occupancy exponent $p$
  (A12) is **tied to $\kappa$ by default** (0 net new free knobs); it becomes a
  3rd free knob only if the size-distribution first-shell cutoff forces the
  $p\!\leftrightarrow\!\kappa$ split.

**Headline:** after sourcing/derivation/bounding, the genuinely free
calibration collapses to **the ladder shape and the electronic picture**, both
discriminated by one observable (the Tier-2 size distribution), with the
early-window scalars bounded and cross-checked rather than fit. Everything else
is locked, imported with a prior, or computed. *(2026-07-06 amendment: the
dominant cascade-magnitude lever is now the Bounded $s_\text{eff}$, co-anchored
with τ by the 9 Å staircase via first-shed timing — the size distribution
arbitrates the remaining shape/picture questions plus the $s_\text{eff}$
small-$n$ tail.)*

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
key open empirical question for the calibration campaign. *(2026-07-06: first
load shed — the 9 Å anchored staircase (TDDFT prior) co-anchors
$s_\text{eff}$×τ, separated by first-shed timing; and the probe showed κ and
picture near-flat on the staircase, so their weight rests entirely on the size
distribution. The $s$↔ladder-shape joint fit resolved weak: κ cannot move the
near-edge rungs enough to matter.)*
