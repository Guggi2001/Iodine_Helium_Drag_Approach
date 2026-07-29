# Tier-2 Drag State Coupling — Design Draft (s(n): the drag learns the ion has stripped)

**Status: DRAFT FOR DISCUSSION (2026-07-29). Physics definition only — no
code exists, none is written until the open questions below are adjudicated
and the user gives `[PROCEED TO IMPLEMENTATION]`. Nothing here moves
`finc1v725` or any G4 adjudication.**

---

## 1. Motivation — the evidence chain (pointers, not re-derivation)

1. **The corrected geometry halves n = 1 KE** (0.641 ± 0.003 vs the
   experimental peak ≈ 1.0 eV; findings "G4 Step 2" Block V, D0 §14.3) and
   the user's data-trust axiom fixes the target: n = 1 at ~1.0 eV must be
   reachable under central births.
2. **The (v_c, τ, E₀) surface cannot reach it** — §3.5g retro-scan R5:
   additive in-gate bound +0.19 eV vs required +0.31 eV.
3. **No re-selection knob can reach it** — the composition ceiling
   (§3.5h registration): the h405 n = 1 KE population is capped at
   0.734 eV; a mean of 1.0 cannot be re-weighted into existence.
4. **Velocity-only softening is measured-killed** — §3.5h ring: uniform
   tail relief un-damps the entire cascade (midHot 1.9–4.4, deepKE up to
   12.5×, PT-P3 kill fired); the surviving γ(v) shape (band-limited
   v_c2 ≈ 9.5) is rejected as overfitting (no physical scale at 9.5;
   break placed between two observable bins; user adjudication
   2026-07-29).
5. **The needle degeneracy** — under any γ(v)-only law, terminal n and
   final KE are two readouts of one path integral, so per-bin KE
   collapses to a needle (measured: 0.64 ± 0.03 over 1556 ions across
   R ≈ 30–70 Å droplets). The experiment shows n = 1 both *fast* and
   *wide* — n and KE are not locked in reality. Breaking the lock
   requires a second, per-ion state variable in the force law.
6. **The physical candidate is already in the mechanism**: the shell
   count n(t). Our γ(v) is the drag of the *dressed* snowball (TDDFT
   calibration; extraction mass I + 19 He stamped in the bundle) applied
   unchanged to ions that have stripped to n ≤ 2. Drag scales with
   effective cross-section; a near-bare I⁺ is geometrically ~40 % of the
   dressed area — the same order as the required toll reduction
   (0.71 → ~0.3 eV). The incumbent's shallow-birth artifact was hiding
   exactly this missing dependence (short paths made the error small).

## 2. Physical picture

The moving object is not "the iodine ion" but "the iodine snowball
currently carrying n He." As the RRK cascade strips the shell, the
object shrinks; the ram/displacement coupling to the surrounding liquid
shrinks with its cross-section. Ions that strip early therefore spend
most of their path as a small object → lower toll → exit fast → land at
low n. Ions that keep their shell brake hard → thermalize → land deep
and slow. The coupling variable *is* the bin variable — selectivity by
construction, and stochastic strip-timing decorrelates n from KE (the
needle breaks).

## 3. Governing form and strict dimensional analysis

The friction convention is unchanged (γ is a force coefficient,
`F_drag = γ·v`, no mass in the drag module). One dimensionless state
factor joins the existing dimensionless spatial gate:

```
γ(v, d, n) = g(d) · s(n) · γ_form(v)
F_drag(v, d, n) = γ(v, d, n) · v
```

| symbol | meaning | units |
|---|---|---|
| γ_form(v) | the locked form (capped_cubic: b·v² in-band, b·v_c²·(v/v_c)^p_tail tail) | amu/ps |
| g(d) | spatial (density) gate, existing | dimensionless, ∈ [0, 1] |
| s(n) | **new** state factor, shell count n | dimensionless, > 0 |
| F_drag | force magnitude | amu·Å/ps² |

FDT noise amplitude `√(2·γ·k_B·T_eff)` uses the full γ — s enters the
(Tier-3, stubbed) noise channel automatically; no separate convention
needed. Drag power `P = γ·v²` [amu·Å²/ps³] feeds the E_int partition
exactly as today; the 5-term invariant is structurally unchanged (only
the magnitude of the drag-work term is modulated).

**Geometric closure (the proposed s):**

```
s(n) = ( R_eff(n) / R_eff(n_ref) )²
R_eff(n) = ( R_core³ + 3n / (4π·ρ_shell) )^(1/3)
```

| parameter | meaning | units | class (CALIBRATION_MAP) | prior |
|---|---|---|---|---|
| R_core | effective collision radius of the bare I⁺ snowball core | Å | **Bounded** (I⁺–He potential minimum / snowball literature) | ≈ 3.0–3.6 |
| ρ_shell | He number density in the attached shell | Å⁻³ | **Bounded** (bulk 0.0218 … snowball-enhanced ≈ 2×) | ≈ 0.03 |
| n_ref | normalization shell count | – | **Derived** (the bundle's extraction mass 202.95 = I + 19 He) | 19 |
| exponent 2 | area scaling of a geometric cross-section | – | **Fixed by geometry** (a free q only as a diagnostic arm) | 2 |

Units check: `R_core³` [Å³] + `n/ρ_shell` [Å³] → R_eff [Å] → s
dimensionless ✓. Illustrative magnitudes at the priors (R_core 3.2,
ρ_shell 0.03): R_eff(19) ≈ 5.68 Å, and

| n | 21 | 19 | 14 | 8 | 5 | 2 | 1 | 0 |
|---|---|---|---|---|---|---|---|---|
| s(n) | 1.06 | 1.00 | 0.85 | 0.65 | 0.54 | 0.40 | 0.37 | 0.32 |

The bare-end value s ≈ 0.32–0.37 is the **prior**, not a fit — and it
sits at the order the KE₁ deficit requires (toll ratio ≈ 0.45 for a
full-path-bare idealization; progressive stripping needs somewhat less).

## 4. Tier-0 / Tier-1a legitimacy

- **Tier-0 is byte-identical by construction:** Tier-0 runs are
  `mass_scenario=fixed` with no shell state; the driver passes s ≡ 1
  (the `off` behavior). The form lock and the calibrated b are untouched.
- **The TDDFT calibration is the dressed state:** normalizing at
  n_ref = 19 (the extraction-mass shell) makes s a *relative* correction
  to exactly what was calibrated. In the Tier-1a/2 window the shell runs
  ≈ 21 → 14, i.e. s ∈ [0.85, 1.06] — a ±15 % in-window modulation, NOT
  byte-identical for Tier-2 runs. Consequence: the corrected-geometry
  basin (v_c 5.5 / τ 4.4 / E₀ ≈ 0.4) may shift mildly; the probe runs at
  the h405 pins first and measures the shift instead of assuming none.
- The large-|effect| regime (s < 0.6) activates only after heavy
  stripping — dominantly beyond the TDDFT window, i.e. in the same
  extrapolated regime the tail occupies. Unlike the tail, the
  extrapolation variable (n) is mechanism-owned and independently
  constrained by the size distribution itself.

## 5. Architecture placement (minimal violence)

- The drag module stays **mass-agnostic and state-blind**: it continues
  to expose γ_form(v) with the gate; the **driver** computes s(n(t))
  from the per-ion shell state it already tracks (the biphasic m(t) uses
  the same variable — no new state is introduced) and multiplies, exactly
  as it does the spatial gate. Seam: the `gamma_fn(speed, depth)`
  closure built in `simulation/ion.py` gains the per-ion factor.
- **Enum surface:** `drag_state_coupling ∈ {"off", "shell_area"}`,
  default `"off"` = bit-identical current behavior (every existing run
  reproduces). Coefficients `{R_core_angstrom, rho_shell_per_A3}` live
  beside it with a config-load guard (positive; s(0) ∈ (0, 1);
  n_ref from the bundle stamp — consistency-checked, not free).
- The rule-2 carry convention applies (fields declared before the
  implementing slice reads them → recorded in the tier log).

## 6. Trade-offs and lost physics (working-method disclosure)

1. **Separability assumption:** γ = g(d)·s(n)·γ_form(v) assumes state
   and velocity dependence factorize. TDDFT cannot check this above
   5.58 Å/ps; the probe's KE(n) grading is the empirical check.
2. **Pickup asymmetry:** λ_attach keeps its own occupancy/capture model
   (Langevin-type physics, not geometric ram) — the two cross-sections
   are deliberately NOT tied (OQ-D). Documented asymmetry.
3. **Cooling contact not scaled (v1):** Newton cooling τ is He-contact
   mediated; physically a stripped ion also cools slower. Deliberately
   out of v1 to keep one new mechanism at a time (OQ-F); revisit if the
   probe under-delivers or over-strips.
4. **No shape/wake anisotropy:** s is scalar; leading-edge/wake
   asymmetries stay unmodeled.
5. **Noise inherits s** via the FDT convention (consistent, but noted
   for Tier-3).

## 7. Sign-level predictions (magnitudes pre-registered at probe design, PT-P4 lesson applied)

- **P-A (target):** KE₁ rises toward ~1.0 at the geometric priors —
  without fitting them to it.
- **P-B (the discriminator):** the n = 1 KE needle **breaks** — per-bin
  SD grows from 0.03 to ≳ 0.1 eV. No γ(v) form can do this; it is the
  coupling's signature prediction.
- **P-C (grading):** KE gain ordered n = 1 > n = 2 > n = 3; deep bins
  (n ≥ 10, s ≥ 0.7 throughout) small.
- **P-D (the risk, stated honestly):** mid-n enders also strip
  substantially en route (s(5) ≈ 0.54), so midHot WILL rise; the claim
  is only that the differential is stronger than p_tail's (relief tracks
  the cascade state instead of a shared velocity band). The probe
  carries a midHot kill criterion exactly like PT-P3.
- **P-E (no confident sign):** the n̄/n₁/supp back-reaction couples the
  residence, pickup and heating channels with opposite signs — after the
  PT-P4 inversion, these are declared **exploratory reads**, not
  predictions.

## 8. Probe sketch (registered in detail only after the OQs are adjudicated)

2–3 cells × N = 1000 at the h405 pins, seed 20260729 (CRN-paired):
`shell_area` at the geometric priors; one weaker-coupling arm (e.g.
ρ_shell at the bulk bound, or the q = 1 radius-scaling diagnostic); the
committed h405 row as the s = 1 baseline. Oracles: cfg-diff vs committed
h405 in exactly the new enum + coefficients; scorer O1/O2 as in §3.5h.
MD-first (KE axes are not twin-scannable; Block 0). Kill criterion:
midHot, PT-P3 form.

## 9. Open questions for adjudication (blocking; discuss before any code)

- **OQ-A — E2 stage:** does s apply inside the relaxation stage's
  `landau_gated_drag` too? *Recommendation: yes — one force law
  everywhere; a stage-split coupling would be a new convention to defend.*
- **OQ-B — form behind the enum:** geometric two-parameter R(n) closure
  (recommended: physical priors, exponent fixed) vs a one-parameter
  power law s = (n_eff/n_ref)^q (fewer symbols, but q is a free fit
  knob — the overfitting shape again).
- **OQ-C — normalization:** n_ref = 19 from the bundle's extraction-mass
  stamp (recommended: s is then a correction *relative to the calibrated
  state*) vs n₀ = 21 (the initial dressing).
- **OQ-D — pickup:** leave λ_attach's cross-section untied (recommended)
  — accept the documented asymmetry?
- **OQ-E — noise:** confirm s enters the FDT amplitude via γ (automatic;
  Tier-3 consequence only).
- **OQ-F — cooling contact:** keep τ un-scaled in v1 (recommended), with
  the physical caveat on record?
- **OQ-G — the design un-freeze:** the drag surface gains a per-ion
  state input (the v-only freeze dates from Tier-0 fixed-mass). Formal
  user adjudication required to un-freeze.
- **OQ-H — §6.5 guard:** the mass↔coefficient pairing guard reads the
  bundle's extraction stamps; s is a separate surface and should not
  touch it — confirm no new guard coupling is wanted beyond the
  n_ref-consistency check.
- **OQ-I — basin re-finding:** if the probe lands KE₁ but shifts the
  gate observables, is the follow-up a local (v_c, τ, E₀) re-tune around
  h405 (the §3.5g slopes are pre-s and must be re-measured), or a fresh
  G3-style twin+MD pass? (Cost question; can wait for the probe.)
