# Slice 1 — Goals: Pure Gated-Drag Physics Module

**Status:** Specification. Physics definition + data contract only. No
implementation code until `[PROCEED TO IMPLEMENTATION]`.

**Scope:** Slice 1 *only* — the config-free, mass-free drag-physics module
that occupies the `collisions.py` swap-point seam. The BAOAB step (Slice 2),
the `SimConfig` surface and §6.5 guard (Slice 3), and the driver rewiring
(Slice 4) are explicitly out of scope and listed under "Scope fence" below.

**Why this is the first task.** Zero upstream dependencies; it is the literal
swap-point content; it mirrors the baseline's cleanest existing seam
(`collisions.py` is pure physics with config injected as kwargs); and it
front-loads the highest-information physics question — *does deterministic
gated `linear_cubic` reproduce the TDDFT force-balance scatter in-window?* —
into the cheapest testable unit, before any integrator exists to confound it.

---

## 1. Convention (unified, per the corrected decisions doc)

A single friction convention now runs through §1.2 / §3.8 / §4.3 / §4.4:

- **$\gamma(v)$ is a force coefficient**, units **amu/ps**, defined
  $\gamma(v) = |F_\text{drag}(v)| / v$. The friction force is
  $\gamma(v)\,v$ — **no leading $m$**.
- The friction **rate** is $\gamma/m$ [1/ps]; it appears *only* inside the
  BAOAB damping exponent $e^{-\gamma\,dt/m}$ (Slice 2, not here).
- The FDT noise amplitude $\sqrt{2\,\gamma\,k_B T_\text{eff}}$ (Slice 1.x /
  §1.2, not active in Slice 1) uses $\gamma$ [amu/ps] directly.

**Consequence for this module:** because $\gamma$ is a force coefficient, the
module is **fully mass-agnostic** — it never takes $m$ as an argument. The
only site where mass enters is the integrator's O-step, as one explicit
division by $m(t)$. This removes the factor-of-$m_\text{eff}$ hazard from the
module boundary entirely.

---

## 2. Governing equations

Primary form `linear_cubic`, with the erf-complement spatial gate
(§5.5 G4-collapsing-to-G2):

$$F_\text{drag}(v,\text{depth}) = g(\text{depth}) \cdot \big(a\,v + b\,v^3\big)$$

$$\gamma(v,\text{depth}) = g(\text{depth}) \cdot \big(a + b\,v^2\big)$$

$$g(\text{depth}) = \tfrac12\Big(1 - \mathrm{erf}\big(\text{depth}/\text{steepness}\big)\Big)$$

with $\text{depth} = r_\text{atom} - r_\text{droplet}$ (negative inside the
droplet, positive outside), as in the baseline droplet potential
(`potentials.py:73`).

The drag force opposes motion: applied as $-F_\text{drag}$ along $\hat v$ by
the consumer (sign restored at the integrator, per the extraction-side
convention).

---

## 3. Strict dimensional analysis (gate before adoption)

| quantity | units | balance check |
|---|---|---|
| $v$ | Å/ps | — |
| $a$ | amu/ps | $a\,v$: $(\text{amu/ps})(\text{Å/ps}) = \text{amu·Å/ps}^2$ ✓ |
| $b$ | amu·ps/Å² | $b\,v^3$: $(\text{amu·ps·Å}^{-2})(\text{Å/ps})^3 = \text{amu·Å/ps}^2$ ✓ |
| $g$ | dimensionless | erf and argument $(\text{depth}/\text{steepness})$ both dimensionless ✓ |
| $\text{depth}, \text{steepness}$ | Å | ratio dimensionless ✓ |
| $F_\text{drag}$ | amu·Å/ps² | a force; matches Coulomb/droplet accel after $\div m$ ✓ |
| $\gamma$ | amu/ps | $g\,(a + b v^2)$: $a$ already amu/ps, $b v^2 = (\text{amu·ps·Å}^{-2})(\text{Å}^2\text{ps}^{-2}) = \text{amu/ps}$ ✓ |

**Units balance. No formulation rejected on dimensional grounds.**

---

## 4. Module interface (three pure functions, zero mass)

- `drag_force(v, depth, coeffs, steepness)` → amu·Å/ps²
- `drag_gamma(v, depth, coeffs, steepness)` → amu/ps (force coefficient;
  later consumed by *both* the O-step rate $\gamma/m$ and the FDT noise
  amplitude)
- `spatial_gate(depth, steepness)` → dimensionless $\in [0,1]$

No `SimConfig` dependency; coefficients and steepness injected as arguments
(the `collisions.py` kwarg pattern). The module holds no mass and no state.

**Closed-form $\gamma$, not definitional division.** Expose $\gamma$ via its
**closed form** $g\,(a + b v^2)$, *not* via $|F_\text{drag}|/v$. They are
analytically equal, but the division is numerically singular at $v \to 0$
even though the `linear_cubic` limit is finite ($\gamma \to g\,a$). The O-step
and FDT noise both need $\gamma$ near rest, where the physics is well-behaved;
computing by division would manufacture a singularity the form does not have.
This is a physics-definition point, not a coding nicety — and it is exactly
the $v\to0$ boundary where `power_law` genuinely *does* diverge and needs the
§3.8 floor. Making the closed-form-vs-division choice **explicit per form** is
the correct factoring.

---

## 5. Form dispatch (logic tree)

```
drag_form?
├─ linear_cubic     → implement fully.
│                      Coeffs {a, b} in hand for both 9 Å and 18 Å.
│                      Guard: a > 0; turnover v† = √(−a/b) above max speed.
│                      [OPEN: max trajectory speed — see §8]
├─ linear_quadratic → raise NotImplemented (no fit pass yet, §3.7)
├─ threshold        → raise NotImplemented (no fit pass yet, §3.7)
└─ power_law        → explicit deferral, not silent.
                       Coeffs {γ, n} exist, but needs §3.8 low-v floor and
                       the divergent-γ handling; out of Tier-0 scope.
```

Only `linear_cubic` is realised in Slice 1. The other three are reserved
behind the same dispatch so adding them later is a branch, not a signature
change.

---

## 6. Data contract — coefficient-bundle *type*

Slice 1 **owns the bundle type**, not the §6.5 guard (which lives in Slice 3).
Written against the *type*, with values plumbed in at test time.

The frozen record carries:

- `form` — the drag-form tag (`linear_cubic` for Slice 1).
- coefficients, form-tagged and variable-arity: `{a, b}` for `linear_cubic`,
  units amu/ps and amu·ps/Å².
- `extraction_mass_model ∈ {constant, time_resolved}` plus the constant
  value (or the $m(t)$ reference) it was extracted under — the metadata the
  §6.5 consistency guard will later read. Slice 1 *defines and carries* this
  field; it does not *enforce* anything with it.

The module consumes a bundle; it never constructs one from config.

---

## 7. Acceptance criteria (test surface — values plumbed at test time)

Mirrors `tests/physics/test_collisions.py` in spirit: pure-physics, no stepper.

- **Force-balance reproduction (the killer test).** Given the in-hand
  $\{a,b\}$ for 9 Å and 18 Å, `drag_force` re-evaluated on the trusted
  interior reproduces the extraction's own $(v, F_\text{drag})$ scatter
  within fit tolerance. This is an integrator-free physics milestone:
  it re-derives the extraction's fit and confirms the module *is* the
  extracted law.
- **Dissipativity.** $F_\text{drag}\cdot v \ge 0$ in magnitude opposing
  motion across the operating range (i.e. drag never adds energy). If
  $b < 0$, assert $v_\dagger = \sqrt{-a/b}$ exceeds the (provisional)
  max trajectory speed.
- **Gate limits.** $g(-\infty) = 1$, $g(0) = \tfrac12$, $g(+\infty) = 0$;
  monotone-decreasing in depth; $C^1$ (the continuity the discarded sharp
  gate G1 lacked).
- **FDT coupling carrier.** `drag_gamma` carries the *same* $g(\text{depth})$
  that scales `drag_force` (§5.2 hard FDT coupling) — assert identical gate
  factor, so the noise amplitude the future O-step reads is gated
  consistently with the drag.
- **Low-$v$ regularity.** $\gamma(0^+, \text{depth}) \to g\,a$, finite — no
  floor needed for `linear_cubic`. Contrast asserted against the documented
  `power_law` divergence (which is deferred, §5).
- **Mass-agnosticism.** No function takes $m$. $F_\text{drag}$ and $\gamma$
  depend only on $(v, \text{depth}, \text{coeffs}, \text{steepness})$.

---

## 8. Open items surfaced (recorded, not blockers)

- **Max trajectory speed** for the `linear_cubic` turnover guard
  (§3.8 / §6.10). Needed before the dissipativity assertion is *exact*;
  use a generous ceiling for now. Has no home in the baseline config; §6.10
  sources it from the TDDFT trajectories or a ceiling.
- **Guard-not-vacuous check.** Confirm the extracted $\{a,b\}$ actually
  satisfy $b \ge 0$ (or $v_\dagger$ comfortably above any speed the ion
  reaches). A 5-minute numeric check at test time once real values are
  plumbed in — validates that the dissipativity guard is not testing an
  empty condition.

---

## 9. Scope fence — what Slice 1 does NOT touch

- **No BAOAB step.** O-step, damping exponent, OU update → Slice 2.
- **No `SimConfig` field additions.** Only the coefficient-bundle *type* is
  defined here; the enum surfaces and the §6.5 consistency guard → Slice 3.
- **No `ion_propagation_step` wiring.** Call-site drops/replacements
  (`:211–256`) → Slice 4.
- **No checkpoint v6, no energy rename.** `E_mass_attach_defect_eV` →
  `E_mass_transfer_eV` and the schema bump → later, with mass dynamics (§2).
- **No noise.** FDT amplitude is defined-but-inactive; Slice 1 only ensures
  $\gamma$ is exposed in the gated, closed form the noise will consume.
- **No mass dynamics.** Module is mass-free by construction; $m(t)$ lives
  entirely in the integrator.
- **`collisions.py` left intact and importable.** Parallel, not deleted —
  CLAUDE.md forbids removal without explicit approval. Slice 1 is additive.

---

## 10. Definition of done

- The three pure functions specified, `linear_cubic` realised, other forms
  dispatched-but-reserved.
- The coefficient-bundle type defined with `extraction_mass_model` metadata.
- All §7 acceptance criteria green against plumbed-in 9 Å and 18 Å values.
- The two §8 open items recorded in-repo (not necessarily resolved).
- `collisions.py` untouched; no config, integrator, checkpoint, or driver
  changes.
