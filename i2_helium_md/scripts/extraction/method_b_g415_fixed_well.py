"""Fixed-well 18 A refit of the lq ridge point c = 2.0 (thesis discussion 2026-10-06).

D0 §1.G4.14 "surface-well depth E_bind fixes both": the committed ring runs at
E_bind 0.154 / 0.117 / 0.048 eV show the shallow well raising the n >= 15
abundance and the deep-bin KE toward experiment. E_bind is Tier-0 co-fitted
with the drag (Locked, CALIBRATION_MAP row 3; Tier-0 family spread 0.048-0.154
eV). This driver asks the Tier-0 question only: does the 18 A TDDFT trace still
accept the lq law F = a v + c v^2 at c = 2.0 when E_bind is FIXED at a shallow
value, with a refitted? Same setup, objective and bands as G4.12 step 1
(``method_b_g412_ridge.py``: N = 50, seed 20260604, 20 ps at 0.01 ps, 3.5 A
gate; band RMSE <= 0.19 A/ps and escape = 1).

Units: a [amu/ps], c [amu/A], E_bind [eV]; objective / RMSE [A/ps].

Modes (USER SETTINGS ``MODE``):

* ``"oracle"`` -- 1-D fit of a at the committed c2p000.json E_bind reproduces
  the committed a and RMSE (|da| <= ORACLE_DA, |dRMSE| <= ORACLE_DRMSE: the
  committed point is a 2-D Nelder-Mead optimum, this is a 1-D bounded Brent
  search, so agreement is to optimizer tolerance, not bit-exact). Writes nothing.
* ``"fit"`` -- the 1-D fit at each E_bind in ``E_BIND_GRID_EV``; one JSON per
  point in ``OUT_DIR``. Refuses to overwrite.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_g415_fixed_well.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

from scipy.optimize import minimize_scalar

# =============================================================================
# USER SETTINGS
# =============================================================================
MODE = "oracle"   # "oracle" | "fit"

C_FIXED = 2.0                         # amu/A, the A5 ridge point
E_BIND_GRID_EV = (0.0482, 0.080)      # shallow end of the Tier-0 spread; midway
A_BOUNDS = (0.0, 60.0)                # amu/ps, a >= 0 as in lq_fixed_c
XATOL_A = 0.01                        # amu/ps, Brent tolerance
ORACLE_DA = 0.05                      # amu/ps
ORACLE_DRMSE = 2.0e-4                 # A/ps

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import evaluate_form_objective  # noqa: E402
from i2_helium_md.physics.drag import LINEAR_QUADRATIC  # noqa: E402
from scripts.extraction.form_phase_common import dump_json, load_json, stamp  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    CASE, FORCE_PROBE_SPEEDS_APS, H405_PLATEAU, RIDGE_DIR, RMSE_18A_MAX_APS, TM_ROOT,
    _setup, lq_force, point_label,
)

OUT_DIR = TM_ROOT / "g415_fixed_well"


def well_label(c: float, e_bind_eV: float) -> str:
    """File tag, e.g. (2.0, 0.0482) -> ``c2p000_e0p0482``."""
    return f"{point_label(c)}_e{e_bind_eV:.4f}".replace(".", "p")


def fit_a_at_fixed_well(setup, c: float, e_bind_eV: float) -> dict:
    """Bounded 1-D minimization of the Method-B objective over a [amu/ps] at
    fixed c [amu/A] and E_bind [eV]. Returns the point record (bands included)."""
    evals = []

    def objective(a: float) -> float:
        r = evaluate_form_objective(LINEAR_QUADRATIC, {"a": float(a), "c": c}, e_bind_eV, setup)
        evals.append((float(a), r))
        return float(r.objective)

    t0 = time.perf_counter()
    opt = minimize_scalar(objective, bounds=A_BOUNDS, method="bounded",
                          options={"xatol": XATOL_A})
    a = float(opt.x)
    r = evaluate_form_objective(LINEAR_QUADRATIC, {"a": a, "c": c}, e_bind_eV, setup)
    forces = {f"F_{v:g}".replace(".", "p"): float(lq_force(a, c, v))
              for v in FORCE_PROBE_SPEEDS_APS}
    rec = {
        **stamp(),
        "goal": "D0 §1.G4.14 fixed-well Tier-0 check (thesis discussion 2026-10-06)",
        "form": LINEAR_QUADRATIC, "variant": "lq_fixed_c_fixed_well",
        "c_amu_per_A": c, "e_bind_eV": e_bind_eV, "a_amu_per_ps": a,
        "rmse_18A_Aps": float(r.rmse_Aps), "escape_fraction": float(r.escape_fraction),
        "objective_Aps": float(r.objective), **forces,
        "Phi": float(lq_force(a, c, 9.7)) / H405_PLATEAU,
        "bands_pass": bool(r.rmse_Aps <= RMSE_18A_MAX_APS and r.escape_fraction == 1.0),
        "n_evaluations": len(evals) + 1, "converged": bool(opt.success),
        "a_bounds": list(A_BOUNDS), "xatol_a": XATOL_A,
    }
    print(f"    c={c} E_bind={e_bind_eV}: {time.perf_counter() - t0:.0f} s, a={a:.3f}, "
          f"rmse={rec['rmse_18A_Aps']:.4f}, escape={rec['escape_fraction']:.3f}, "
          f"bands {'PASS' if rec['bands_pass'] else 'FAIL'}", flush=True)
    return rec


def run_oracle() -> None:
    committed = load_json(RIDGE_DIR / f"{point_label(C_FIXED)}.json")
    rec = fit_a_at_fixed_well(_setup(), C_FIXED, float(committed["e_bind_eV"]))
    da = abs(rec["a_amu_per_ps"] - committed["a_amu_per_ps"])
    dr = abs(rec["rmse_18A_Aps"] - committed["rmse_18A_Aps"])
    print(f"  oracle: |da| {da:.4f} (<= {ORACLE_DA}), |dRMSE| {dr:.2e} (<= {ORACLE_DRMSE})")
    if da > ORACLE_DA or dr > ORACLE_DRMSE:
        raise AssertionError("oracle failed: 1-D fit does not reproduce c2p000.json")
    print("  ORACLE PASS")


def run_fit() -> None:
    paths = {e: OUT_DIR / f"{well_label(C_FIXED, e)}.json" for e in E_BIND_GRID_EV}
    existing = [p.name for p in paths.values() if p.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite {existing}")
    setup = _setup()
    for e, path in paths.items():
        dump_json(fit_a_at_fixed_well(setup, C_FIXED, float(e)), path)


def main(mode: str = MODE) -> int:
    {"oracle": run_oracle, "fit": run_fit}[mode]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else MODE))
