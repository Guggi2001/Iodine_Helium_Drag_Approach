# The twin — 1D chord forward model (zero MD)

Thesis-writing reference (2026-10-21) for the results chapter §5.5: how the
twin obtains size-resolved energies without full MD, what dynamics it keeps,
its approximations, its terminal-shell treatment, and its measured authority.
No new analysis — every statement is read from the code or quoted from a
committed record.

**Code:** `scripts/tier2_h2b_forward_model.py` (the G4.12–G4.14 stages use the
"corrected master" path). Line numbers refer to HEAD 5f877e6.

## 1. Principle

The twin never integrates a shell. Per fragment it integrates **one
fixed-mass 1D trajectory along the explosion axis**, which yields three
numbers: the cooling exposure K, the asymptotic speed v_inf and a trapped
flag. A **closed-form fate map** turns K into a detected size n_det. The
size-resolved energy is ½·m(n_det)·v_inf², averaged per n.

## 2. Retained dynamics

- **Ensemble** (`_g3_corrected_master`, :2336): `raw` droplet sizes
  (⟨N⟩ 12 794, nozzle 40 bar / 14 K), Boltzmann birth radii in the 313.2 K
  parent well, uniform axis cosine μ; m = 20 000 molecules, fixed seed,
  cached once.
- **Pair kinematics** (`integrate_pairs`, :215): two fragments start back to
  back at separation 2.666 Å and move along ±axis on straight lines. The
  radius r(s) along the chord is computed exactly (:308), so density and
  well are sampled along the true chord through the droplet.
- **Forces** (`accel`, :314–336):
  - pair Coulomb repulsion k_C / r_sep²;
  - erf solvation well (depth E_bind, width 14.2 Å), projected onto the
    direction of motion;
  - gated drag −ρ̂(depth; s_ρ)·γ(v)·v, forms `cubic` / capped cubic / `lin` /
    `linq` (γ = a + c·|v|, the ridge). The density width s_ρ (3.5 Å for
    G4.12) is decoupled from the well width, mirroring production's
    `erf_independent` gate.
- **Integrator:** Heun / RK2, dt 0.01 ps, t_end 150 ps. The residual Coulomb
  energy at t_end is added analytically, split equally between the two
  fragments (:358–362).
- **Cooling exposure:** K = ∫ρ̂ dt / τ (:342). K scales exactly as 1/τ, so
  one integration per drag law serves every τ by rescaling (:3827); E₀ costs
  nothing extra. This makes (τ × E₀) surface scans almost free.
- **Trapped:** a fragment within 2 × 14.2 Å of the surface at t_end (:367);
  excluded from detection.

## 3. Approximations relative to MD

| Twin | MD (production, B2) |
|---|---|
| **Fixed mass** m = I⁺ + n_eject·He over the whole trajectory (:3745–3747) | variable mass, SQ1–SQ3 mass jumps |
| **No pickup** during flight | Poisson pickup λ₀ 0.9; re-pickup beyond chord ≈ 40 Å |
| **Birth dressing** n_eject = rint(21·ρ̂(r₀ − R)) on the 14.2 Å width | n₀ = 21 for every ion |
| **No E_int reservoir, no RRK kinetics**: closed-form Newton cooling + instantaneous ladder descent (ε = 0) | E_int ODE, energy-gated RRK (s_eff 8, ν 2.42), 5-term ledger |
| Shed He carries no momentum: velocity unchanged, KE at the final mass | `co_moving` shed convention (equivalent) |
| Deterministic, symmetric pair (same r₀, opposite cosines) | Stage I to 500 ps, partner-aware Coulomb closure, t_detect 8.53 µs |
| Straight-line 1D motion, no transverse deflection | full 3D |

Instantaneous descent to the energetic floor agrees with MD only while the
RRK cascade finishes before detection; per `THESIS_DISCUSSION.md` §10 this
holds for s_eff ≲ 8 (production: 8).

## 4. Terminal shell occupancy (`fate_map`, :380)

1. Energy at ejection: E_ej = E₀ · (Σ(n_eject)/Σ(21))^p · e^(−K); p = 1 on
   the corrected path.
2. **Suppressed** (n = 0) if E_ej > Σ(n_eject): the shell is self-unbound at
   ejection.
3. Otherwise **exact descent**: n_det = smallest k with
   Σ(k) > Σ(n_eject) − E_ej, capped at n_eject. Σ is the cumulative rung sum
   from the MD generator's own rung tables (`rq4graded` / `floor1`,
   `_g3_md_rung_tables`, :2116), so the ladder is identical in twin and MD.
4. **Scoring** (`g3_score`, :2057) follows the MD scorer's conventions:
   KE = ½·m(n_det)·v_inf² (:2086); per-n means over bins with ≥ 20
   fragments; midHot geometric over n 2–8; deepKE arithmetic over n 10–17;
   W₁ on the solvated histogram.

## 5. Measured authority (twin vs MD)

Transfer over the 14-cell MD family
(`TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md:1023`):

| Observable | Agreement | Bias (twin − MD) | Verdict |
|---|---|---|---|
| KE₁ | Pearson 0.9998 | −0.014 eV | licensed as a ranking; separately licensed at Spearman 0.998 over 26 cells (`ke1auth`) |
| tail n ≥ 10 | Pearson 0.964 | ×1.7–1.9 too high | licensed |
| n₁ | Pearson 0.943 | +0.0008 | usable |
| n̄ | Pearson 0.921 | +0.59 | usable after offset |
| midHot | ratio MD/twin 0.906 ± 0.004 | — | usable after offset |
| **W₁** | Pearson 0.319 | — | **refuted**; the sign of the τ effect comes out wrong |

## 6. Role in the path to B2

- **G4.12 step 3:** the twin ranked the 11 ridge points that pass the 9 Å
  bound (one KE₁ / midHot dial).
- **G4.14 Step 0:** twin values moved to predicted MD values with offsets
  measured at c 2.0 (`g413_transfer.csv`, `G414_OFFSETS` :4432) and used to
  filter which c entered MD.
- **G5:** c and E₀ were placed by direct MD, not by the twin, because the
  twin was validated on the `raw` ensemble only (G5 plan, row B1).

Every reported B2 number is MD. The twin's role was **screening and ranking**
(KE₁, midHot, tail), never landing or the size histogram.

## 7. Caveats for the thesis text

- The 14-cell transfer was measured on the free-form linear family; the
  ridge-specific offsets come from the 16 G4.13 cells at c 2.0 only.
- The twin's birth dressing uses the 14.2 Å width; production starts every
  ion at n₀ = 21.

**Suggested framing:** *"a deterministic fixed-mass chord model with
closed-form cooling and instantaneous ladder descent, licensed for ranking
KE₁, midHot and the n ≥ 10 tail, but not for the size-histogram distance;
all reported results are MD."*
