"""G4.12 step 3 twin wiring (`scripts/tier2_h2b_forward_model.py`).

- the 3.5 A density gate reaches the linq chord (toy ensemble): a narrower
  gate changes the trajectory, and ``rho_steepness = STEEP_A`` is
  bit-identical to the default (byte-inert pass-through);
- the refinement grid contains the linscan grid (the refine anchor relies
  on it);
- ``_g412_point_summary`` on synthetic rows: best-KE1 pick, empty gate.

Toy ensembles / synthetic rows only — no caches, no corrected master.
"""

from __future__ import annotations

import numpy as np

from scripts.tier2_h2b_forward_model import (
    G3SCAN_E0_GRID,
    G3SCAN_TAU_PS,
    G412_REFINE_E0_GRID,
    G412_REFINE_TAU_PS,
    STEEP_A,
    _g412_point_summary,
    integrate_pairs,
)

TOY = dict(
    r0=np.array([5.0, 15.0, 25.0]),
    mu=np.array([0.3, -0.5, 0.9]),
    R_drop=np.array([40.0, 40.0, 40.0]),
    mass_amu=np.array([203.0, 203.0, 203.0]),
)
LINQ = dict(drag_form="linq", lin_a=16.7, linq_c=2.0, t_end=5.0)


def test_rho_steepness_reaches_linq_and_default_is_inert():
    base = integrate_pairs(**TOY, **LINQ)
    same = integrate_pairs(**TOY, **LINQ, rho_steepness=STEEP_A)
    narrow = integrate_pairs(**TOY, **LINQ, rho_steepness=3.5)
    for k in ("K", "v_inf"):
        assert np.array_equal(base[k], same[k])
    assert not np.array_equal(base["v_inf"], narrow["v_inf"])


def test_refine_grid_contains_linscan_grid():
    assert set(G3SCAN_TAU_PS) <= set(G412_REFINE_TAU_PS)
    assert {float(e) for e in G3SCAN_E0_GRID} <= {float(e) for e in G412_REFINE_E0_GRID}


def _row(gate, ke1, midhot=1.0, w1=0.6, tau=4.8, e0=0.35):
    return {"gate": gate, "n1_ke_eV": ke1, "midhot_geo": midhot,
            "w1_solv": w1, "deepke": 0.7, "tau_flag": 0, "phi": 0.84,
            "trapped_frac": 0.02, "tau_ps": tau, "E0_eV": e0}


def test_point_summary_picks_best_gated_ke1():
    pt = {"c": 2.0, "a": 16.7, "e_bind": 0.1175, "margin": 0.78}
    rows = [_row(1, 0.55, midhot=0.9), _row(0, 0.99), _row(1, 0.62, midhot=1.1, e0=0.31)]
    s = _g412_point_summary(pt, rows)
    assert s["n_gated"] == 2 and s["n_cells"] == 3
    assert s["best_KE1"] == 0.62 and s["best_KE1_midhot"] == 1.1
    assert s["midhot_min"] == 0.9 and s["midhot_max"] == 1.1


def test_point_summary_no_gate():
    pt = {"c": 3.0, "a": 13.5, "e_bind": 0.118, "margin": 0.56}
    s = _g412_point_summary(pt, [_row(0, 0.5)])
    assert s["n_gated"] == 0 and np.isnan(s["best_KE1"])
