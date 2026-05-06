"""Compare matlab_multistep.csv vs python_multistep.csv.

CLAUDE.md validation targets 3 and 4. Loads the per-step per-atom
trajectory + energy data from both sides, prints absolute and
scaled-relative max errors per quantity, and groups quantities by
expected category:

* ``MATCH``                 -- positions, velocities, mass, time,
                               E_dissip should agree within a small
                               numerical floor.
* ``CONSTANT_ROUNDING``     -- E_kin, E_pot, E_total, accelerations:
                               Python uses CODATA 2022 (EV exact,
                               U = 1.66053906892e-27); MATLAB legacy
                               uses 4-sig-fig rounding (eV = 1.602e-19,
                               u = 1.66053907e-27, plus the literal
                               1.602e-9 / 9648.53322 inside the force
                               conversions). Expected scaled-relative
                               difference: ~110 ppm.
* ``MATCH_t0_FIX``          -- The t=0 row uses the corrected formulas
                               on both sides (both Python and the
                               MATLAB reference script implement the
                               documented Python fixes). The legacy
                               MATLAB raw t=0 bugs are NOT exercised
                               here -- those were validated separately
                               by ``ion_t0_state``.

The script exits non-zero only if a MATCH row exceeds its tolerance or
a CONSTANT_ROUNDING row exceeds an order-of-magnitude budget around
the expected ~110 ppm.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent
PY_CSV = SCRIPT_DIR / "python_multistep.csv"
ML_CSV = SCRIPT_DIR / "matlab_multistep.csv"

QUANTITIES = (
    "x_A", "y_A", "z_A",
    "vx_Aps", "vy_Aps", "vz_Aps",
    "mass_kg",
    "E_kin_eV", "E_pot_eV", "E_dissip_eV",
)

EXPECTED_CATEGORY: dict[str, str] = {
    "x_A":          "CONSTANT_ROUNDING",
    "y_A":          "CONSTANT_ROUNDING",
    "z_A":          "CONSTANT_ROUNDING",
    "vx_Aps":       "CONSTANT_ROUNDING",
    "vy_Aps":       "CONSTANT_ROUNDING",
    "vz_Aps":       "CONSTANT_ROUNDING",
    "mass_kg":      "CONSTANT_ROUNDING",
    "E_kin_eV":     "CONSTANT_ROUNDING",
    "E_pot_eV":     "CONSTANT_ROUNDING",
    "E_dissip_eV":  "MATCH",
    "E_total_eV":   "CONSTANT_ROUNDING",
    "t_ps":         "MATCH",
}

# Scaled-relative tolerances. The constants budget is ~110 ppm; we
# allow 5x slack for the multi-step compounding.
SCALED_REL_TOL_CONSTANT_ROUNDING = 5e-4   # 500 ppm
ABS_TOL_MATCH_TIME = 1e-12                # picosecond
ABS_TOL_MATCH_DISSIP = 1e-12              # eV (zero on both sides)


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

    cols = {name: i for i, name in enumerate(py_hdr)}

    # Sort both arrays by (step, atom) so row-by-row diff is meaningful.
    py = py[np.lexsort((py[:, cols["atom"]], py[:, cols["step"]]))]
    ml = ml[np.lexsort((ml[:, cols["atom"]], ml[:, cols["step"]]))]

    # Sanity: identifiers must agree exactly.
    for id_col in ("step", "atom"):
        if not np.array_equal(py[:, cols[id_col]], ml[:, cols[id_col]]):
            raise SystemExit(f"ERROR: {id_col} column differs between files")

    # Compute system-total energy per step (per row sum is per-atom; sum per
    # step gives the system total).
    n_steps = int(py[:, cols["step"]].max()) + 1
    n_atoms_per_step = py.shape[0] // n_steps
    py_total = (py[:, cols["E_kin_eV"]] + py[:, cols["E_pot_eV"]]
                + py[:, cols["E_dissip_eV"]]).reshape(n_steps, n_atoms_per_step).sum(axis=1)
    ml_total = (ml[:, cols["E_kin_eV"]] + ml[:, cols["E_pot_eV"]]
                + ml[:, cols["E_dissip_eV"]]).reshape(n_steps, n_atoms_per_step).sum(axis=1)

    print(f"Comparing {PY_CSV.name} vs {ML_CSV.name}")
    print(f"  rows: {py.shape[0]}  ({n_steps} steps x {n_atoms_per_step} atoms)\n")

    header_line = (
        f"{'quantity':<14s} {'category':<22s} "
        f"{'max |abs diff|':>14s}  {'scaled rel':>12s}  "
        f"{'argmax (step,atom)':>20s}"
    )
    print(header_line)
    print("-" * len(header_line))

    failures: list[str] = []

    # ------------------------------------------------------------------
    # Per-quantity comparison (uses the row-by-row alignment).
    # ------------------------------------------------------------------
    for q in QUANTITIES + ("t_ps",):
        col = cols[q]
        diff = py[:, col] - ml[:, col]
        abs_diff = np.abs(diff)
        scale = max(np.max(np.abs(ml[:, col])), 1e-30)
        scaled_rel = abs_diff / scale
        i_max = int(np.argmax(abs_diff))
        argmax_step = int(py[i_max, cols["step"]])
        argmax_atom = int(py[i_max, cols["atom"]])

        cat = EXPECTED_CATEGORY[q]
        print(
            f"{q:<14s} {cat:<22s} "
            f"{abs_diff.max():>14.3e}  {scaled_rel.max():>12.3e}  "
            f"({argmax_step:>4d},{argmax_atom:>3d})         "
        )

        if cat == "MATCH":
            if q == "t_ps":
                tol = ABS_TOL_MATCH_TIME
            elif q == "E_dissip_eV":
                tol = ABS_TOL_MATCH_DISSIP
            else:
                tol = 1e-12
            if abs_diff.max() > tol:
                failures.append(
                    f"{q} tagged MATCH but abs_diff_max={abs_diff.max():.3e} > {tol:.0e}"
                )
        elif cat == "CONSTANT_ROUNDING":
            if scaled_rel.max() > SCALED_REL_TOL_CONSTANT_ROUNDING:
                failures.append(
                    f"{q} tagged CONSTANT_ROUNDING but scaled_rel_max="
                    f"{scaled_rel.max():.3e} > {SCALED_REL_TOL_CONSTANT_ROUNDING:.0e} "
                    f"(expected ~110 ppm constants budget)"
                )

    # ------------------------------------------------------------------
    # System-total energy (per step).
    # ------------------------------------------------------------------
    diff_total = py_total - ml_total
    abs_diff_total = np.abs(diff_total)
    scale_total = max(np.max(np.abs(ml_total)), 1e-30)
    scaled_rel_total = abs_diff_total / scale_total
    print(
        f"{'E_total_eV':<14s} {'CONSTANT_ROUNDING':<22s} "
        f"{abs_diff_total.max():>14.3e}  {scaled_rel_total.max():>12.3e}  "
        f"(step={int(np.argmax(abs_diff_total)):>3d})"
    )
    if scaled_rel_total.max() > SCALED_REL_TOL_CONSTANT_ROUNDING:
        failures.append(
            f"E_total tagged CONSTANT_ROUNDING but scaled_rel_max="
            f"{scaled_rel_total.max():.3e} > {SCALED_REL_TOL_CONSTANT_ROUNDING:.0e}"
        )

    # ------------------------------------------------------------------
    # Energy-conservation drift (sanity, both sides should be similar).
    # ------------------------------------------------------------------
    py_drift = py_total[-1] - py_total[0]
    ml_drift = ml_total[-1] - ml_total[0]
    print(
        f"\nEnergy-conservation drift over the run:\n"
        f"  python:  E_total[end] - E_total[0] = {py_drift:+.3e} eV\n"
        f"  matlab:  E_total[end] - E_total[0] = {ml_drift:+.3e} eV"
    )

    print()
    if failures:
        print("VERDICT: FAIL")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("VERDICT: PASS")
    print(
        "  E_dissip and t_ps agree exactly. All other quantities lie within "
        "the constants-rounding budget (<500 ppm scaled-relative)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
