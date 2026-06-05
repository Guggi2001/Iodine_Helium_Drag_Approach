# Task — t\*-Seed From Real 3D Velocities (retire the reconstruction)

**Status:** Task specification. Replaces the rotation-invariant *reconstruction*
in the t\*-seeded comparison harness with a read of the **actual 3D per-atom
velocity components** (now available). No implementation until
`[PROCEED TO IMPLEMENTATION]`. Data-in-hand, MD-side; **no external dependency**.

---

## 1. Why this task exists

The t\*-seeded 9 Å "residual" (same-smoothed |v2| ≈ 0.40 Å/ps, characterised as a
near-constant ~10% over-damping in `TIER0_FINDINGS.md` / `tier0_comparison_tasks_left.md`)
was traced to a **harness reconstruction artifact**, confirmed by a zero-cost
control:

- The seed builder (`seed_neutral_at_reference`) could not recover the full 3D
  velocity, so it reconstructed the scored rotation-invariants: radial from
  `dR/dt` (an `h=0.1` central difference on **raw** `R(t)`), and **all leftover
  speed dumped into a made-up transverse (x) component**.
- At 9 Å t\*=2.67: `|v2|=4.895`, radial `0.22`, **transverse `4.89`** — i.e. the
  seeded speed is **~99.9% inert transverse reconstruction**.
- The central-force MD (Coulomb + radial droplet + radial drag) has **no force
  along that transverse direction**, so the drag (acting on `|v|`) immediately
  damps the 4.89 Å/ps of inert KE → a **strong startup downcurve**, then the
  radial-Coulomb dynamics take over → a **slope parallel to the reference**, with
  the damped transverse KE left as a **permanent offset**. Exactly the observed
  shape.
- **Control (zero cost, run):** placing the leftover speed **radially** instead
  of transverse **removed the startup damping collapse entirely.** Confirmed: the
  residual is a seed-construction artifact, **not** a drag-coefficient error.

**Consequences already established (do not re-litigate):**

- The doubled 9 Å clean-window `a` (13.86→24.876) is **off the critical path** —
  it is not the cause of the Tier-0 residual. (It may remain a separate
  extraction-side curiosity; not this task.)
- `dR/dt` from raw `R(t)` is **not trustworthy for derivatives** (raw `R` carries
  the 1.2 ps bubble oscillation; the `h=0.1` difference reads the oscillation, not
  the separation rate). This task **retires the finite difference**, it does not
  refine it.
- The "~10% magnitude residual blocking Tier 1" framing is **obsolete** and must
  be corrected in the docs (separate doc-update pass, not this code task).

---

## 2. What changes

Replace the **reconstruction** of per-atom velocity in
`seed_neutral_at_reference` (in `tier0_tstar_seeded_comparison.py`) with a
**read of the real 3D per-atom velocity components** at t\*, from the new 3D
reference data.

Specifically retire:

- `vz1, vz2 = +dRdt/2, -dRdt/2` — the radial-from-finite-difference step.
- `vperp = sqrt(|v|² - vz²)` with the leftover **dumped transverse** — the step
  that injected the inert KE.
- the `h=0.1` central difference on raw `R(t)` entirely.

Replace with: interpolate the actual `(vx, vy, vz)` of each atom at t\* from the
3D reference and seed them directly.

---

## 3. The decision the 3D data forces — radial-projected vs. full-3D seed

With true per-atom 3D velocity, the seed *could* be the faithful full vector. But
the control result already tells us the central-force MD **cannot sustain a
transverse component** — it will damp any non-radial velocity on startup,
artifact or not, because there is no transverse force. So seeding the *real*
transverse component faithfully would still produce a (smaller, but real-magnitude)
startup damp. The choice:

**Primary — radial-projected seed (recommended).** Project the real 3D velocity
onto the I–I axis and seed only the component the MD can faithfully evolve:
$$v_{\text{rad},i} = (\vec v_i \cdot \hat R)\,\hat R, \qquad
\hat R = (\vec r_1 - \vec r_2)/|\vec r_1 - \vec r_2|.$$
Seed `vz_i = v_rad,i` (along the I–I axis = z), transverse components **zero**.
Rationale: the t\*-seeded comparison is a *clean-form* diagnostic isolating the
**drag law's behaviour on the radial dynamics the MD actually models**. Feeding it
a transverse component the model structurally cannot carry only re-introduces the
startup-damp confound (now with a real but still-unmodellable magnitude). The
control confirmed radial-only is the clean choice.

- The radial rate is now **exact from 3D** ((v₁−v₂)·R̂), replacing the discredited
  `dR/dt` finite difference. This is the concrete win: a trustworthy radial seed.

**Secondary — full-3D faithful seed (diagnostic only).** Seed the real
`(vx,vy,vz)`. Use **only** to *quantify* how much transverse KE the real data
carries at t\* (vs. the 4.89 artifact) and how much the MD damps it — i.e. to
measure whether atom 2's large non-radial speed is physically real. This is the
genuine, better-founded version of the co-drift question (now on the *clean* atom,
from *real* 3D data, not the discarded artifact atom). It is a **measurement**, not
the comparison seed. If it shows a large *real* transverse component, that is a
physics finding (atom 2 genuinely co-translates) to record — but it does **not**
change that the central-force MD models only the radial dynamics.

**Discarded — keep the reconstruction.** The `dR/dt` + transverse-dump
reconstruction is the artifact source; remove it, do not retain it behind a flag.

---

## 4. The check this immediately enables — the real radial fraction at t\*

The artifact diagnosis hinged on `|v2|` being 99.9% transverse *by reconstruction*.
With real 3D, compute the **true** radial vs. transverse split at t\*:

- `v_rad,2 = (v₂·R̂)`; `v_perp,2 = |v₂ − v_rad,2 R̂|`.
- **Print both**, as the harness already prints the reconstructed split. This
  answers the question the reconstruction could not: *is atom 2 genuinely
  near-radial at t\* (small real transverse), or does it really carry a large
  non-radial speed?*
  - small real transverse → the 4.89 was purely reconstruction error; radial seed
    is unambiguously right; residual fully resolved.
  - large real transverse → atom 2 really co-translates; record as a physics
    finding (the well-founded co-drift observation), but the radial-projected seed
    is *still* the right comparison seed (the MD models radial dynamics).

---

## 5. Data schema — the 3D reference the seed reads

The seed needs, at minimum, per-atom 3D **velocity** at t\* (and the positions /
`R̂` to project onto). Expected (align to the file actually produced; the existing
CSV uses `V1_z,V2_z,V1_x,V2_x,...` so extend rather than rename if practical):

```
Time_ps,
X1, Y1, Z1,  X2, Y2, Z2,           # per-atom positions (Å) — for R̂ and projection
VX1, VY1, VZ1,  VX2, VY2, VZ2,     # per-atom velocities (Å/ps) — the new y-components
R_distance                          # |r1 − r2| (Å); cross-check |r1−r2| == R_distance
```

- The loader change is additive: a reader for the 3D columns. If the 3D data
  lives in a **new** file (e.g. `9A_3D_All_Data.csv`), add a loader path; if the
  y-components were **added to the existing** file, extend `load_hedft_trajectory`
  to expose `vy`. Either way, keep the existing rotation-invariant columns working
  (the scored comparison still uses `R`, `|v1|`, `|v2|`) — this task changes the
  **seed**, not the scored quantities.
- **Cross-check on load:** assert `|r1 − r2|` from the positions matches the
  stored `R_distance` (provenance sanity, same discipline as elsewhere).

---

## 6. Scope fence

- **Only the seed builder changes.** The scored quantities (`R`, `|v1|`, `|v2|`
  via `compare_*`), the windowing, the from-onset harness, the regression gate,
  and the drag physics are **untouched**.
- **No coefficient change.** `{a,b}` are not refit here; the doubled-`a` question
  is explicitly out of scope.
- **No production-path change.** The real simulation seeds from the explosion
  onset with physical ICs; it never used this reconstruction. This task fixes the
  *diagnostic harness* only.
- **Doc updates are a separate pass** (TIER0_FINDINGS, tier0_comparison_tasks_left,
  CLAUDE.md — to record the artifact resolution and unblock Tier 1). Flagged, not
  done here.

---

## 7. Definition of done

- `seed_neutral_at_reference` reads real per-atom 3D velocity at t\*; the `dR/dt`
  finite difference and the transverse-dump are removed.
- The seed uses the **radial-projected** velocity (primary, §3); transverse seeded
  zero. The radial rate comes from `(v₁−v₂)·R̂`, not from differencing `R(t)`.
- The harness prints the **real** radial/transverse split at t\* (§4).
- Loader reads the 3D columns with the `|r1−r2|==R_distance` cross-check (§5).
- Re-run 9 Å and 18 Å t\*-seeded same-smoothed: confirm the 9 Å startup downcurve
  and the ~0.40 offset are **gone** (the control already showed radial-only
  removes the collapse; this makes it the permanent, real-data-backed seed).
- A one-line verdict: with the artifact removed, is the 9 Å t\*-seeded
  same-smoothed residual now ~18 Å-class? If yes → Tier 0 is a clean
  internal-consistency pass and Tier 1 unblocks (pending the separate doc pass).
