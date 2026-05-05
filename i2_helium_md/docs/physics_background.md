# Physics Background

This document explains the physical model behind the simulation so that the
code reads as clear physics rather than unmotivated formulas.

> The simulation reproduces a **pump–probe ionization experiment** on iodine
> molecules embedded inside superfluid helium nanodroplets. A laser pulse
> dissociates and ionizes I₂ inside the droplet; the simulation tracks the
> resulting atoms/ions as they fly through the helium and eventually escape.
> Measured velocity distributions of ejected species are compared to the
> simulation output.

---

## 1. The three regimes of interaction

Each iodine atom in the simulation feels **three distinct forces simultaneously**.
Each force has its own potential function — none of them are redundant.

### Regime A: droplet solvation  (`droplet_potential`)

Question: *"Am I inside or outside the helium droplet?"*

Governs the cost of escaping the droplet. Modelled as a smooth erf profile:

```
V(r) = ((erf((r - offset) / steepness) + 1) / 2) * binding_energy
```

- `r < 0`  (atom inside droplet)  →  V ≈ 0
- `r ≈ 0`  (atom at the surface)  →  V = binding_energy / 2
- `r ≫ 0`  (atom outside)         →  V ≈ binding_energy

**Parameters**

| Parameter | Physical meaning | Typical value |
|---|---|---|
| `steepness`  | width of the fuzzy droplet surface  | ≈ 14 Å (a few atomic layers) |
| `binding_energy`  | asymptotic cost to fully leave the droplet  | 27 meV (I atom), 49 meV (I₂) |
| `offset`  | shift of the transition point along r  | 0 Å |

**Why erf and not a step function?** Because the droplet surface has finite
thickness — superfluid helium's density profile is smooth, not discontinuous.

---

### Regime B: chemistry with the partner atom  (`morse_X`, `morse_I2plus_state_select`)

Question: *"What's my other-half iodine doing to me?"*

The two atoms of an I₂ molecule feel each other through the chemical bond. This
is **short-range** physics (electron exchange, Pauli repulsion, bond breaking).

#### Neutral I–I: `morse_X(r, cfg)`

Standard Morse potential with spectroscopic parameters for the I₂ X ground
state:

```
U(r) = D_e * (1 - exp(-a·(r - R_e)))²
```

| Parameter | Physical meaning | Value (I₂ X state) |
|---|---|---|
| `D_e`  | dissociation energy  | 1.556 eV |
| `R_e`  | equilibrium bond length  | 2.666 Å |
| `ω_e` → `a`  | vibrational stiffness via `_morse_a()`  | derived |

**What this potential contains — and what it does NOT:**

- ✅ Covalent bonding (attractive well near R_e)
- ✅ Pauli repulsion at very short r (steep exponential wall)
- ✅ Van der Waals at long r (asymptote to D_e)
- ❌ **No Coulomb** — neutral atoms have zero net charge

#### The `Xdip` correction

When `cfg.Xdip_active = True` (default), a narrow Gaussian dip is subtracted
from the X-state potential at r = 9 Å:

```
U(r) ← U(r) - 0.9 · exp(-(r - 9)² / (2 · 0.3²))
```

This is an **empirical correction** modelling an avoided crossing between
electronic states. Without it, the neutral MD dynamics do not match the
He-DFT reference calculations at R₀ = 9 Å. The thesis discusses this fix in
the context of the "9 Å problem."

#### Ionized I₂⁺: four electronic states

The `morse_I2plus_state_select()` function handles a subtle case: after
ionization, I₂⁺ can end up in any of four low-lying electronic states. Each
state has its own Morse curve with different D_e, R_e, ω_e. The function
randomly assigns one state per molecule and uses the right curve.

| State ID | D_e (eV) | R_e (Å) | IP_rel (eV) |
|---|---|---|---|
| 0 (ground)  | 2.70  | 2.61  | 0.00 |
| 1  | 2.03  | 2.61  | 0.63 |
| 2  | 1.26  | 2.95  | 1.68 |
| 3  | 0.56  | 2.95  | 2.44 |

Baseline IP₀ = 9.36 eV (lowest ionization potential of I₂).

> **Important — when these Morse I₂⁺ curves are actually used:** *not* for the
> standard I⁺–I⁺ case. See §2 below.

---

### Regime C: Coulomb repulsion  (added inline in `add_partner_interaction_ion`)

Question: *"Are we both charged? If so, how hard do we repel?"*

For two singly-charged iodine ions:

```
U(r) = E_coulomb_scale · q1 · q2 · 14.39964548 / r    [eV, r in Å]
```

The constant `14.39964548 eV·Å` is just `e² / (4πε₀)` converted to these units.
`E_coulomb_scale` is an empirical knob (default 1.0) the author uses to
"dial down" the Coulomb repulsion to study screening effects from the
helium environment.

---

## 2. Morse vs. Coulomb for I⁺–I⁺: the design decision

This is a subtle point that confused us during the port.

**The standard I⁺–I⁺ ion–ion interaction does NOT use the Morse I₂⁺ curves.**
It uses pure Coulomb:

```matlab
U = E_coulomb_scale * q1 * q2 * 14.4 / r;   % ion_interaction_potential.m
```

The Morse I₂⁺ state-select code path is only triggered when
`single_charge_ionization_allowed = True` **and** exactly one atom of the pair
is ionized (i.e. the I–I⁺ asymmetric case). In the standard double-ionization
scenario used for the HeDFT comparison, the Morse I₂⁺ curves remain inactive.

### Why pure Coulomb is defensible for I⁺–I⁺

1. **Length-scale:** ions starting at R ≈ 2.666 Å fly past R ≈ 5 Å within
   ~0.1 ps of ionization. Beyond 5 Å, all four Morse curves have
   asymptoted to flat. The Morse↔Coulomb difference lives in a narrow
   early-time window.

2. **Energy-scale:** at the starting distance R = 2.666 Å, the Coulomb term
   is 14.4 / 2.666 ≈ 5.4 eV, about 2× the Morse well depth (2.7 eV for the
   lowest I₂⁺ state). Coulomb dominates the final ejection velocity.

3. **Modeling philosophy:** `E_coulomb_scale` is a free parameter the author
   fits to data. Keeping the ion-ion term simple (one parameter) and
   absorbing unknowns into that scale factor is cleaner than layering Morse
   curves on top with poorly known branching ratios.

### When pure Coulomb starts breaking down

From the thesis summary:

> "once the I⁺ are further apart than about 9 Å, the MD ion velocities are
> much lower than the expected ion velocities at 200 ps. Furthermore, there
> is no good way to model the recombination of the molecule at 9 Å."

This is an open research problem. If you ever want to:

- Study early-time dynamics (first 0.1 ps after ionization)
- Model vibrationally-excited bound ion states
- Investigate recombination near R ≈ 9 Å

…then the Morse I₂⁺ state-select code **is** the right tool. It's been ported
and tested in Step 4; enabling it is a one-line config change.

---

## 3. Summary table of interactions

| Scenario  | Code path  | Physics |
|---|---|---|
| Neutral I–I (short r)  | `morse_X(r, cfg)`  | Chemistry only, no Coulomb |
| Neutral I–I (long r)  | `morse_X` asymptote  | V → D_e, no long-range force |
| I⁺–I⁺ (standard)  | pure Coulomb `14.4/r`  | No chemistry, scaled by `E_coulomb_scale` |
| I–I⁺ (if `single_charge_ionization_allowed`)  | Coulomb + `morse_I2plus_state_select`  | Chemistry + Coulomb mix |
| Atom in droplet field  | `droplet_potential`  | Solvation, independent of partner |

---

## 4. The story of one molecule, start to finish

Following a single I atom through the simulation makes the layering explicit:

```
t = 0, before laser
  - Position: ~30 Å from droplet center (inside a ~40 Å droplet)
  - Partner: other I atom at R = 2.666 Å (equilibrium)
  - Forces: droplet_potential ≈ 0 (deep inside)
            morse_X(R) near its minimum
            no Coulomb

t = 0+, laser fires
  - Photon energy hc/λ_pump promotes to antibonding state
  - Kinetic energy ≈ (hc/λ - E_diss) / 2 per atom
  - Atoms start flying apart

t ≈ 1–10 ps, dissociation phase (neutral scope)
  - Forces: droplet_potential still ≈ 0 (still inside)
            morse_X(R) repulsive (R > R_e, climbing the Morse wall)
            hard-sphere collisions with He atoms slow the atom down

t ≈ 50–200 ps, approaching droplet surface
  - Forces: droplet_potential contributes real force at r ≈ droplet radius
            morse_X(R) negligible (R now ~20 Å, past the asymptote)
            fewer He collisions as density drops near surface

--- if simulation includes ionization ---

t = t_ionize, ionization event (handled by Step 11 `ion_propagation`)
  - Neutral becomes I⁺
  - Forces: droplet_potential with different binding energy for ions
            Coulomb repulsion 14.4/R (no more Morse)
            He collisions with different (larger) cross section
            Possible He attachment ("snowball") if enabled

t ≈ 200 ps, escape
  - Atom/ion leaves droplet with some final velocity
  - This final velocity is what the VMI experiment measures
```

---

## 5. Author's open questions (from thesis)

Kept here for context — these are NOT problems we aim to solve in the port,
but help explain why parameters like `E_coulomb_scale`, `sigma_ion_exponent`,
and `binding_energy_I_ion_eV` have the defaults they do.

1. **The σ(v) question.** For I⁺ + He scattering, a cross section of
   σ₀ = 1900 Å² matches the 18 Å HeDFT data but gives too-fast ions from
   2.666 Å. A cross section of σ₀ = 5000 Å² matches 2.666 Å but over-brakes
   at 18 Å. This suggests σ(v) = σ₀ · v⁻ⁿ with n > 2, i.e. slower ions scatter
   more often (current exponent: -2).

2. **The 9 Å problem.** Once ions are past ~9 Å, MD velocities
   fall short of observed values. Candidate explanations include non-Coulomb
   ion-ion forces at intermediate range and solvation effects not captured
   by the erf droplet model.

3. **Recombination.** Current code has no mechanism for I⁺ + I⁺ → I₂⁺
   recombination or I⁺ + e⁻ → I capture. This limits the model's accuracy
   at long times and when charges can neutralize.

---

## References

- J. Chem. Phys. **107**, 9046 (1997). doi: 10.1063/1.475194
  — I₂ spectroscopic constants (D_e, ω_e, R_e) for X and I₂⁺ states.
- Phys. Rev. B **58**, 3341 (1998). doi: 10.1103/PhysRevB.58.3341
  — bulk liquid helium density 0.0219 Å⁻³.
- Stadlhofer thesis (unpublished) — original MATLAB code and experimental
  context. See legacy repo `README.md` for scope notes.
