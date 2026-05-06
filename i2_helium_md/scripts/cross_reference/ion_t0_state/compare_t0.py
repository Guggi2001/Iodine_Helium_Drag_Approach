"""Compare matlab_t0.csv vs python_t0.csv for the ion t=0 state.

Loads both CSVs, prints a per-quantity diff table tagged with the
expected category (MATCH / INTENTIONAL_FIX_* / CONSTANT_ROUNDING),
and exits non-zero only if a row tagged ``MATCH`` exceeds tolerance.

Run from anywhere:

    python scripts/cross_reference/ion_t0_state/compare_t0.py

Categories
----------
* ``MATCH``                 -- both sides should agree to ~1e-12.
* ``INTENTIONAL_FIX_KE``    -- Python correct, legacy MATLAB has the
                               t=0 kinetic-energy bug at line 289 of
                               vmi_sim_3d_ion_propa.m.
* ``INTENTIONAL_FIX_PE``    -- Python correct, legacy MATLAB has the
                               t=0 potential-energy bug at line 291
                               (2D radial + missing partner Coulomb).
* ``CONSTANT_ROUNDING``     -- Python uses CODATA 2022; MATLAB legacy
                               uses 4-sig-fig rounding of the same
                               constants, giving sub-ppm differences.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent
PYTHON_CSV = SCRIPT_DIR / "python_t0.csv"
MATLAB_CSV = SCRIPT_DIR / "matlab_t0.csv"

# Per-quantity expected category. Used to label rows in the printed
# table and to decide which rows are gated for pass/fail.
EXPECTED_CATEGORY: dict[str, str] = {
    "x_A":              "MATCH",
    "y_A":              "MATCH",
    "z_A":              "MATCH",
    "vx_Aps":           "MATCH",
    "vy_Aps":           "MATCH",
    "vz_Aps":           "MATCH",
    "mass_kg":          "CONSTANT_ROUNDING",
    "droplet_radius_A": "MATCH",
    "time_ps":          "MATCH",
    "E_kin_eV":         "INTENTIONAL_FIX_KE",
    "E_pot_eV":         "INTENTIONAL_FIX_PE",
    "E_dissip_eV":      "MATCH",
}

# Tolerance for rows tagged MATCH.
ABS_TOL_MATCH = 1e-12

# Sub-ppm tolerance for CONSTANT_ROUNDING. CODATA U is
# 1.66053906892e-27 kg vs MATLAB legacy 1.66053907e-27 kg, a relative
# difference of ~6.5e-8. mass_kg = 127*U so the absolute diff is
# ~1.4e-32 kg. We use a generous 1e-30 kg cap.
ABS_TOL_CONSTANT_ROUNDING_MASS_KG = 1e-30


def _load(path: Path) -> dict[str, tuple[float, float]]:
    if not path.exists():
        raise SystemExit(
            f"ERROR: {path.name} not found in {path.parent}. "
            f"Run the corresponding export script first."
        )
    out: dict[str, tuple[float, float]] = {}
    with open(path) as f:
        header = f.readline().strip().split(",")
        if header != ["quantity", "atom_0", "atom_1"]:
            raise SystemExit(f"ERROR: unexpected header in {path.name}: {header}")
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, v0, v1 = line.split(",")
            out[name] = (float(v0), float(v1))
    return out


def main() -> int:
    py = _load(PYTHON_CSV)
    ml = _load(MATLAB_CSV)

    if set(py) != set(ml):
        only_py = set(py) - set(ml)
        only_ml = set(ml) - set(py)
        raise SystemExit(
            f"ERROR: quantity mismatch. only_python={only_py}, only_matlab={only_ml}"
        )

    # Stable order matches the order in EXPECTED_CATEGORY.
    quantities = list(EXPECTED_CATEGORY.keys())
    missing = [q for q in py if q not in EXPECTED_CATEGORY]
    if missing:
        raise SystemExit(
            f"ERROR: CSV contains quantities without an expected category: {missing}"
        )

    print(f"Comparing  python_t0.csv  vs  matlab_t0.csv\n")
    header = (
        f"{'quantity':<18s} {'atom':>4s}  "
        f"{'matlab':>22s}  {'python':>22s}  "
        f"{'abs diff':>12s}  {'expected':<22s}"
    )
    print(header)
    print("-" * len(header))

    failures: list[str] = []

    for q in quantities:
        cat = EXPECTED_CATEGORY[q]
        m0, m1 = ml[q]
        p0, p1 = py[q]

        for atom_idx, (m, p) in enumerate(((m0, p0), (m1, p1))):
            d = p - m
            ad = abs(d)
            print(
                f"{q:<18s} {atom_idx:>4d}  "
                f"{m:>+22.12e}  {p:>+22.12e}  "
                f"{ad:>12.3e}  {cat:<22s}"
            )

            if cat == "MATCH" and ad > ABS_TOL_MATCH:
                failures.append(
                    f"{q}[atom={atom_idx}] tagged MATCH but abs_diff={ad:.3e} "
                    f"> {ABS_TOL_MATCH:.0e}"
                )
            elif cat == "CONSTANT_ROUNDING" and q == "mass_kg":
                if ad > ABS_TOL_CONSTANT_ROUNDING_MASS_KG:
                    failures.append(
                        f"{q}[atom={atom_idx}] tagged CONSTANT_ROUNDING but "
                        f"abs_diff={ad:.3e} > {ABS_TOL_CONSTANT_ROUNDING_MASS_KG:.0e}"
                    )
            elif cat.startswith("INTENTIONAL_FIX") and ad == 0.0:
                failures.append(
                    f"{q}[atom={atom_idx}] tagged {cat} but Python and MATLAB "
                    f"agree exactly -- the documented Python fix may have regressed."
                )

    print()
    if failures:
        print("VERDICT: FAIL")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("VERDICT: PASS")
    print(
        "  All MATCH rows agree to <1e-12. "
        "Both INTENTIONAL_FIX_* rows show non-zero diffs as expected."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
