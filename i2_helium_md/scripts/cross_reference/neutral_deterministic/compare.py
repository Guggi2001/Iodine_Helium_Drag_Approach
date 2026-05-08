"""Compare python_trajectory.csv vs matlab_trajectory.csv.

Loads both CSV files, prints a side-by-side summary at key timesteps,
and reports max absolute deviations across the whole trajectory.

Run with:
    python compare.py

Expected behavior:
    With identical initial conditions, the two integrators should
    agree to ~ppm precision modulo the small differences in the
    physical constants used:

        Constant        Python              MATLAB legacy        diff
        --------        ------              -------------        ----
        EV [J]          1.602176634e-19     1.602e-19            ~110 ppm
        U  [kg]         1.66053906892e-27   1.66053907e-27       ~0.05 ppm

    The leapfrog converts forces eV/A -> A/ps^2 using EV. So a 110 ppm
    difference in EV produces a ~110 ppm difference in acceleration,
    which compounds over the 100 steps but should remain in the 0.01%
    range for positions and velocities, and similar for energies.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).parent
PYTHON_CSV = SCRIPT_DIR / "python_trajectory.csv"
MATLAB_CSV = SCRIPT_DIR / "matlab_trajectory.csv"


def load_csv(path: Path) -> tuple[list[str], np.ndarray]:
    with open(path) as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return header, data


def main() -> None:
    if not PYTHON_CSV.exists():
        print(f"ERROR: {PYTHON_CSV} not found. Run run_python.py first.")
        return
    if not MATLAB_CSV.exists():
        print(f"ERROR: {MATLAB_CSV} not found. Run run_matlab.m first.")
        return

    py_hdr, py_data = load_csv(PYTHON_CSV)
    ml_hdr, ml_data = load_csv(MATLAB_CSV)

    if py_hdr != ml_hdr:
        print(f"WARNING: header mismatch")
        print(f"  python: {py_hdr}")
        print(f"  matlab: {ml_hdr}")

    if py_data.shape != ml_data.shape:
        print(f"ERROR: shape mismatch: python {py_data.shape} vs matlab {ml_data.shape}")
        return

    print(f"Loaded {py_data.shape[0]} rows × {py_data.shape[1]} columns from each.\n")

    # ------------------------------------------------------------------
    # Side-by-side at selected timesteps
    # ------------------------------------------------------------------
    sample_idx = [0, 1, 10, 25, 50, 75, 100]
    print("=" * 100)
    print("Side-by-side at selected timesteps:")
    print("=" * 100)
    for col in [1, 2, 4, 7, 13, 14, 15]:    # x1, y1, x2, vx1, E_kin, E_pot, E_total
        name = py_hdr[col]
        print(f"\n--- column {col}: {name} ---")
        print(f"{'t [ps]':>8s}  {'python':>22s}  {'matlab':>22s}  {'diff':>14s}  {'rel diff':>12s}")
        for i in sample_idx:
            if i >= py_data.shape[0]:
                continue
            t = py_data[i, 0]
            p = py_data[i, col]
            m = ml_data[i, col]
            d = p - m
            rd = abs(d) / max(abs(m), 1e-30)
            print(f"{t:>8.3f}  {p:>22.12e}  {m:>22.12e}  {d:>+14.3e}  {rd:>12.3e}")

    # ------------------------------------------------------------------
    # Max absolute differences across all timesteps
    #
    # NOTE: pure relative differences blow up when a quantity passes
    # through zero (e.g. E_pot at the Morse equilibrium). We use a
    # "scaled" relative measure that normalizes by the typical magnitude
    # of the column over the whole trajectory, not by the instantaneous
    # value. This gives a meaningful "fraction of typical signal" number.
    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("Max absolute / scaled-relative deviations across all 101 rows:")
    print("(scaled-relative = |diff| / max|column| -- avoids blow-up at zero crossings)")
    print("=" * 100)
    print(f"{'column':<14s} {'max |abs diff|':>16s}  {'scaled-rel':>14s}  "
          f"{'argmax t':>10s}  {'note':<25s}")
    for col in range(1, py_data.shape[1]):
        diff = py_data[:, col] - ml_data[:, col]
        abs_diff = np.abs(diff)
        # Scale by the typical magnitude of the column (max abs over the run);
        # this avoids relative-error blow-up when the quantity passes zero.
        scale = max(np.max(np.abs(ml_data[:, col])), 1e-30)
        rel_scaled = abs_diff / scale
        i_max = int(np.argmax(abs_diff))
        # Flag columns that pass near zero (relative-by-instantaneous-value
        # would have blown up).
        min_abs = np.min(np.abs(ml_data[:, col]))
        note = "(passes near zero)" if min_abs < 0.01 * scale else ""
        print(f"{py_hdr[col]:<14s} {abs_diff.max():>16.6e}  "
              f"{rel_scaled.max():>14.6e}  {py_data[i_max, 0]:>10.3f}  {note:<25s}")

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    pos_max_rel = np.max([
        np.max(np.abs(py_data[1:, c] - ml_data[1:, c]))
        / max(np.max(np.abs(ml_data[1:, c])), 1e-30)
        for c in [1, 2, 3, 4, 5, 6]
    ])
    en_max_rel = np.max([
        np.max(np.abs(py_data[1:, c] - ml_data[1:, c]))
        / max(np.max(np.abs(ml_data[1:, c])), 1e-30)
        for c in [13, 14, 15]
    ])
    print(f"Position scaled-relative max deviation: {pos_max_rel*1e6:.2f} ppm")
    print(f"Energy   scaled-relative max deviation: {en_max_rel*1e6:.2f} ppm")
    if pos_max_rel < 5e-3 and en_max_rel < 5e-3:
        print("\nVERDICT: Python and MATLAB agree within expected tolerance "
              "(<0.5%). Constants differences (~110 ppm in EV, ~0 ppm in U) "
              "account for most of the deviation.")
    else:
        print("\nVERDICT: Python and MATLAB DISAGREE beyond expected tolerance. "
              "Check the integrator, force functions, and unit conversions.")


if __name__ == "__main__":
    main()
