"""Compare matlab_forced.csv vs python_forced.csv at the driver level.

CLAUDE.md validation target 5 (stochastic forced-event driver
plumbing). The collision kernel itself is covered statistically by
``scripts/collision_comparison_test/`` and is NOT retested here.

Three groups of checks, each gated separately:

* **CROSS_MATCH** -- driver bookkeeping that must be deterministic
  on both sides regardless of RNG stream:
  - integer event identifiers (``step``, ``atom``)
  - simulated time ``t_ps``
  - per-step event flags (``b_collision``, ``b_attach``)
  - cumulative ``number_of_collisions``
  - mass over time ``mass_kg`` (within constants-rounding tolerance)
  - per-step ``sigma_used_A2`` at t=0 (analytic, no collision yet)
  - ``depth_A`` at t=0 (deterministic at start)

* **PER_SIDE_INVARIANTS** -- physics invariants that must hold on
  each side independently (no cross-language identity required):
  - ``mass_history`` non-decreasing per atom
  - ``b_attach`` only where ``b_collision`` is True
  - ``E_dissip_eV`` non-decreasing per atom
  - conservation:
    ``E_kin + E_pot + E_dissip + E_mass_attach_defect``
    drifts by less than 1% over the run

* **NOT_COMPARED** -- RNG-dependent quantities printed for
  inspection only: post-collision velocities, post-collision
  energies, depth/sigma after the first collision, and the cumulative
  ``E_mass_attach_defect_eV`` (which depends on RNG-driven ``v_post``).

Exits non-zero if any CROSS_MATCH or PER_SIDE_INVARIANTS check fails.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent
PY_CSV = SCRIPT_DIR / "python_forced.csv"
ML_CSV = SCRIPT_DIR / "matlab_forced.csv"

# Tolerances
ABS_TOL_TIME = 1e-12                    # ps
ABS_TOL_INT = 0                         # integers must match exactly
SCALED_REL_TOL_CONSTANT_ROUNDING = 5e-4 # 500 ppm (constants budget)
ABS_TOL_MASS_KG = 1e-30                 # kg, mass differs only by U rounding
CONSERVATION_DRIFT_TOL_REL = 1e-2       # 1% of |total energy| over the run


def _load(path: Path) -> tuple[list[str], np.ndarray]:
    if not path.exists():
        raise SystemExit(
            f"ERROR: {path.name} not found in {path.parent}. "
            f"Run the corresponding export script first."
        )
    with open(path) as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return header, data


def main() -> int:
    py_hdr, py = _load(PY_CSV)
    ml_hdr, ml = _load(ML_CSV)

    if py_hdr != ml_hdr:
        raise SystemExit(
            f"ERROR: header mismatch.\n  python: {py_hdr}\n  matlab: {ml_hdr}"
        )
    if py.shape != ml.shape:
        raise SystemExit(
            f"ERROR: shape mismatch. python={py.shape}, matlab={ml.shape}"
        )

    cols = {n: i for i, n in enumerate(py_hdr)}

    # Sort both by (step, atom) so row-by-row diff aligns.
    py = py[np.lexsort((py[:, cols["atom"]], py[:, cols["step"]]))]
    ml = ml[np.lexsort((ml[:, cols["atom"]], ml[:, cols["step"]]))]

    for id_col in ("step", "atom"):
        if not np.array_equal(py[:, cols[id_col]], ml[:, cols[id_col]]):
            raise SystemExit(f"ERROR: {id_col} column differs between files")

    n_steps = int(py[:, cols["step"]].max()) + 1
    n_atoms_per_step = py.shape[0] // n_steps

    failures: list[str] = []
    print(f"Comparing {PY_CSV.name} vs {ML_CSV.name}")
    print(f"  rows: {py.shape[0]}  ({n_steps} steps x {n_atoms_per_step} atoms)\n")

    # ------------------------------------------------------------------
    # CROSS_MATCH group
    # ------------------------------------------------------------------
    print("=" * 75)
    print("CROSS_MATCH: must agree across MATLAB and Python")
    print("=" * 75)
    print(f"{'quantity':<22s} {'max |abs diff|':>16s}  {'tolerance':>12s}  {'category':<24s}")
    print("-" * 75)

    def _cross_match_row(name: str, tol: float, category: str) -> None:
        col = cols[name]
        d = np.abs(py[:, col] - ml[:, col])
        ok = d.max() <= tol
        print(
            f"{name:<22s} {d.max():>16.3e}  {tol:>12.3e}  "
            f"{category:<24s}{'  PASS' if ok else '  FAIL'}"
        )
        if not ok:
            failures.append(f"{name}: max abs_diff={d.max():.3e} > tol={tol:.3e}")

    _cross_match_row("t_ps",                 ABS_TOL_TIME, "exact (sim time)")
    _cross_match_row("number_of_collisions", ABS_TOL_INT,  "exact (forced events)")
    _cross_match_row("b_collision",          ABS_TOL_INT,  "exact (forced events)")
    _cross_match_row("b_attach",             ABS_TOL_INT,  "exact (forced events)")
    _cross_match_row("mass_kg",              ABS_TOL_MASS_KG, "constants rounding")

    # sigma_used and depth: only step 0 is RNG-independent on both sides
    # (positions/velocities at t=0 come from inputs.json directly). Later
    # steps depend on RNG-driven velocities. Compare only step 0.
    step0_mask_py = py[:, cols["step"]] == 0
    step0_mask_ml = ml[:, cols["step"]] == 0
    for name in ("sigma_used_A2", "depth_A"):
        col = cols[name]
        d = np.abs(py[step0_mask_py, col] - ml[step0_mask_ml, col])
        scale = max(np.abs(ml[step0_mask_ml, col]).max(), 1e-30)
        rel = d.max() / scale
        ok = rel <= SCALED_REL_TOL_CONSTANT_ROUNDING
        print(
            f"{name + ' (step 0)':<22s} {d.max():>16.3e}  "
            f"{SCALED_REL_TOL_CONSTANT_ROUNDING:>12.3e}  "
            f"{'constants rounding':<24s}{'  PASS' if ok else '  FAIL'}"
        )
        if not ok:
            failures.append(
                f"{name} (step 0): scaled_rel={rel:.3e} > "
                f"{SCALED_REL_TOL_CONSTANT_ROUNDING:.3e}"
            )

    # ------------------------------------------------------------------
    # PER_SIDE_INVARIANTS group (run separately for each side)
    # ------------------------------------------------------------------
    print()
    print("=" * 75)
    print("PER_SIDE_INVARIANTS: must hold independently on each side")
    print("=" * 75)

    def _check_invariants(label: str, data: np.ndarray) -> list[str]:
        side_failures: list[str] = []
        # Reshape to (n_steps, n_atoms_per_step) per quantity.
        def _per_atom(name: str) -> np.ndarray:
            return data[:, cols[name]].reshape(n_steps, n_atoms_per_step)

        mass     = _per_atom("mass_kg")
        b_coll   = _per_atom("b_collision").astype(int)
        b_att    = _per_atom("b_attach").astype(int)
        E_diss   = _per_atom("E_dissip_eV")
        E_kin    = _per_atom("E_kin_eV")
        E_pot    = _per_atom("E_pot_eV")
        E_def    = _per_atom("E_mass_attach_defect_eV")

        # mass non-decreasing per atom
        if (np.diff(mass, axis=0) < -1e-32).any():
            side_failures.append(f"{label}: mass decreased over time on at least one atom")

        # b_attach => b_collision
        if ((b_att > 0) & (b_coll == 0)).any():
            side_failures.append(f"{label}: an attachment fired without a collision")

        # E_dissip non-decreasing per atom
        if (np.diff(E_diss, axis=0) < -1e-12).any():
            side_failures.append(f"{label}: cumulative E_dissip decreased over time")

        # Conservation invariant: E_kin + E_pot + E_dissip + E_defect.
        # Sum across atoms to get system total per step.
        total_per_step = (E_kin + E_pot + E_diss + E_def).sum(axis=1)
        drift = total_per_step[-1] - total_per_step[0]
        rel = abs(drift) / max(abs(total_per_step[0]), 1e-30)
        print(f"  {label}: total[0]={total_per_step[0]:+.6f} eV  total[end]={total_per_step[-1]:+.6f} eV  "
              f"drift={drift:+.3e} eV ({rel*100:+.4f}%)  ", end="")
        if rel > CONSERVATION_DRIFT_TOL_REL:
            print("FAIL")
            side_failures.append(
                f"{label}: conservation drift {rel:.3e} > {CONSERVATION_DRIFT_TOL_REL:.3e}"
            )
        else:
            print("PASS")

        if not side_failures:
            print(f"  {label}: mass-monotonic, attach-implies-collide, E_dissip-monotonic  PASS")
        return side_failures

    failures.extend(_check_invariants("python", py))
    failures.extend(_check_invariants("matlab", ml))

    # ------------------------------------------------------------------
    # NOT_COMPARED group (informational)
    # ------------------------------------------------------------------
    print()
    print("=" * 75)
    print("NOT_COMPARED (RNG-dependent): printed for inspection only")
    print("=" * 75)
    for name in ("vx_Aps", "vy_Aps", "vz_Aps",
                 "E_kin_eV", "E_pot_eV", "E_dissip_eV", "E_mass_attach_defect_eV"):
        col = cols[name]
        # restrict to steps where collisions have occurred at least once
        # (step >= 2 in our scenario): pre-collision values match, post-collision differ
        mask = py[:, cols["step"]] >= 2
        d = np.abs(py[mask, col] - ml[mask, col])
        scale = max(np.abs(ml[mask, col]).max(), 1e-30)
        print(
            f"  {name:<26s} max |abs diff| = {d.max():>10.3e}  "
            f"scaled-rel = {d.max() / scale:>10.3e}  (post-collision; RNG-dependent)"
        )

    print()
    if failures:
        print("VERDICT: FAIL")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("VERDICT: PASS")
    print(
        "  Driver bookkeeping (event counts, mass history, t_ps, sigma at t=0) "
        "agrees on both sides. Per-side conservation invariant holds."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
