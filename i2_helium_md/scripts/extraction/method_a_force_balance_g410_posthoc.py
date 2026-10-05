"""G4.10 post-hoc set: the June Method-A pipeline + the agreed corrections.

POST-HOC (decided 2026-10-05 after the pre-registered G4.10 run, which
replaced the June pipeline's windows, edge trim and derivative). This set
keeps the June Method-A data handling (``Drag_function/drag_function/
drag_calculation.py``) and changes only the physics of the balance:

* data: the short committed references ``velocity_smoothed/cleaned_data.csv``
  (18 A: order-1 SG over 3.4 ps on [4.54, 8.0]; 9 A: IMF + SG on
  [2.67, 6.0]; the June run used an uncommitted local 9 A file whose
  cubic b differs by 0.012 -- the repo file is used here);
* derivative: cubic smoothing spline with the MAD-based ``s``, trusted
  interior = drop 500 samples (0.5 ps) at each end;
* physics (the corrections): Coulomb projected on the motion (cos theta),
  solvation well -dU/dr (r_hat . v_hat) at the G4.9 E_bind, density gate
  s_rho = 3.5 A (F / rho_hat), m_eff.

The G4.10 band (2.2286 +- 2 * 0.2127) is applied descriptively to the 18 A
and 9 A cubic b; it is not a pre-registered gate for this set.

Rows: ``primary`` (all trusted points); rho_hat cut >= 0.5 / >= 0.9;
E_bind +- 1 half-width; radius +- 8 %; spline s x 0.5 / x 2; 9 A anchored
m(t) (t* 2.67 ps); decomposition rows (old: no well, scalar, no gate; + well;
+ well + projection, no gate).

Writes ``data/reference/drag/force_balance_srho3p5A/posthoc_june_pipeline/``:
``fits.csv`` and ``summary.json``.

Usage::

    python scripts/extraction/method_a_force_balance_g410_posthoc.py
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from i2_helium_md.extraction.force_balance import (  # noqa: E402
    FORMS,
    balance_geometry,
    fit_drag_form,
    force_balance_trace,
    spline_speed_and_derivative,
)
from i2_helium_md.physics.constants import droplet_radius_bulk_angstrom  # noqa: E402
from i2_helium_md.physics.shell_schedule import (  # noqa: E402
    build_shell_schedule,
    complex_mass_amu,
)
from i2_helium_md.postprocess.hedft_loader import (  # noqa: E402
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.presets import single_pulse_N2000_18Angst_drag  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
CASES = ("18A", "9A")
REFERENCE_FILENAME = "cleaned_data.csv"   # June short references
TRIM_POINTS = 500                         # June trusted interior: 0.5 ps per end
GATE_STEEPNESS_A = 3.5                    # s_rho (G2)
DIRECTION_WINDOW_PS = 0.5                 # velocity-direction smoothing (G4.10)
RADIUS_SPREAD = 0.08                      # +- 8 % N -> R convention spread (G1)
T_STAR_9A_PS = 2.67                       # anchored m(t) onset
BAND_HALF_WIDTHS = 2.0                    # G4.10 band, applied descriptively
# =============================================================================

DRAG_ROOT = REPO_ROOT / "data" / "reference" / "drag"
G49_BUNDLE = DRAG_ROOT / "18A" / "trajectory_matching_srho3p5A" / "cubic" / "fit_parameters.json"
OUT_DIR = DRAG_ROOT / "force_balance_srho3p5A" / "posthoc_june_pipeline"


def _variants(bundle: dict) -> dict[str, dict]:
    """Every configuration, keyed by name (``primary`` first)."""
    cfg = single_pulse_N2000_18Angst_drag()
    E = float(bundle["effective_binding_energy_I_ion_eV"])
    hw_E = float(bundle["effective_binding_energy_err_eV"])
    R = float(droplet_radius_bulk_angstrom(cfg.single_droplet_size))
    base = {
        "projection": "projected",
        "gate_steepness_A": GATE_STEEPNESS_A,
        "rho_min": 0.0,
        "binding_energy_eV": E,
        "droplet_radius_A": R,
        "potential_steepness_A": float(cfg.potential_steepness),
        "m_eff_amu": float(cfg.m_eff_amu),
        "mass": "m_eff",
        "spline_scale": 1.0,
    }
    rows = {
        "primary": {},
        "rho_cut_0p5": {"rho_min": 0.5},
        "rho_cut_0p9": {"rho_min": 0.9},
        "ebind_minus_hw": {"binding_energy_eV": E - hw_E},
        "ebind_plus_hw": {"binding_energy_eV": E + hw_E},
        "radius_minus8pct": {"droplet_radius_A": R * (1.0 - RADIUS_SPREAD)},
        "radius_plus8pct": {"droplet_radius_A": R * (1.0 + RADIUS_SPREAD)},
        "spline_s_x0p5": {"spline_scale": 0.5},
        "spline_s_x2": {"spline_scale": 2.0},
        "mass_anchored_9A": {"mass": "anchored"},
        "decomp_old_nowell_scalar_nogate": {
            "projection": "legacy_no_well", "gate_steepness_A": None},
        "decomp_well_scalar_nogate": {"projection": "scalar", "gate_steepness_A": None},
        "decomp_well_projected_nogate": {"gate_steepness_A": None},
    }
    return {k: {**base, **d} for k, d in rows.items()}


def _trace(s: dict, geometry, ref):
    """Trusted-interior force-balance trace of one case under settings ``s``."""
    t = np.asarray(ref.time_ps, dtype=float)
    speed, dvdt, s_used = spline_speed_and_derivative(
        t, ref.speed_Aps, smoothing_scale=s["spline_scale"]
    )
    keep = slice(TRIM_POINTS, t.size - TRIM_POINTS)
    if s["mass"] == "anchored":
        mass = complex_mass_amu(build_shell_schedule(T_STAR_9A_PS).n_of_t(t[keep]))
    else:
        mass = s["m_eff_amu"]
    fb = force_balance_trace(
        geometry,
        t[keep],
        speed[keep],
        dvdt[keep],
        mass_amu=mass,
        binding_energy_eV=s["binding_energy_eV"],
        potential_steepness_A=s["potential_steepness_A"],
        droplet_radius_A=s["droplet_radius_A"],
        gate_steepness_A=s["gate_steepness_A"],
        projection=s["projection"],
        interior_rho_min=s["rho_min"],
    )
    return fb, s_used


def main() -> None:
    bundle = json.loads(G49_BUNDLE.read_text())
    b49, hw_b = float(bundle["b"]), float(bundle["b_err"])
    band = (b49 - BAND_HALF_WIDTHS * hw_b, b49 + BAND_HALF_WIDTHS * hw_b)
    data = {}
    for case in CASES:
        traj = load_hedft_trajectory(REPO_ROOT / "data" / "reference" / f"{case}_All_Data.csv")
        ref = load_smoothed_speed_reference(
            DRAG_ROOT / case / "velocity_smoothed" / REFERENCE_FILENAME
        )
        data[case] = (ref, balance_geometry(traj, direction_window_ps=DIRECTION_WINDOW_PS))

    rows, info = [], {}
    for name, s in _variants(bundle).items():
        cases = ("9A",) if name == "mass_anchored_9A" else CASES
        per_case = {}
        for c in cases:
            fb, s_used = _trace(s, data[c][1], data[c][0])
            per_case[c] = fb
            if name == "primary":
                info[c] = {
                    "trusted_window_ps": [float(fb.time_ps[0]), float(fb.time_ps[-1])],
                    "rho_hat_range": [float(fb.rho_hat.min()), float(fb.rho_hat.max())],
                    "spline_s": s_used,
                }
        sets = {c: (c,) for c in cases}
        if len(cases) > 1:
            sets["joint"] = cases
        for set_name, members in sets.items():
            v = np.concatenate([per_case[c].speed_Aps[per_case[c].interior] for c in members])
            F = np.concatenate([per_case[c].F_drag_bulk[per_case[c].interior] for c in members])
            for form in FORMS:
                f = fit_drag_form(v, F, form)
                in_band = (
                    bool(band[0] <= f.b <= band[1])
                    if form == "cubic" and set_name in CASES else None
                )
                rows.append({"variant": name, "case_set": set_name, **asdict(f),
                             "cubic_in_G410_band": in_band})

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "fits.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cubic = {
        f"{r['variant']}/{r['case_set']}": r["b"]
        for r in rows if r["form"] == "cubic"
    }
    summary = {
        "date": "2026-10-05",
        "status": "POST-HOC (after the pre-registered G4.10 run); band applied descriptively",
        "goal": "THESIS_FINAL_GOALS.md G4.10 post-hoc set: June Method-A pipeline + corrections",
        "band_G410": list(band),
        "pipeline": {
            "reference": f"velocity_smoothed/{REFERENCE_FILENAME}",
            "derivative": "cubic smoothing spline, MAD-based s (June)",
            "trim_points_per_end": TRIM_POINTS,
            "projection": "projected", "gate_steepness_A": GATE_STEEPNESS_A,
            "direction_window_ps": DIRECTION_WINDOW_PS,
        },
        "cases": info,
        "cubic_b": cubic,
        "cubic_in_band_all_sensitivity_rows": {
            c: all(
                r["cubic_in_G410_band"] for r in rows
                if r["form"] == "cubic" and r["case_set"] == c
                and not r["variant"].startswith("decomp_")
            )
            for c in CASES
        },
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ("band_G410", "cases",
                                              "cubic_in_band_all_sensitivity_rows")}, indent=2))
    for r in rows:
        print(f"{r['variant']:32s} {r['case_set']:6s} {r['form']:13s} a={r['a']:8.3f} "
              f"b={r['b']:7.4f} rmse={r['rmse_force']:7.3f} n={r['n_points']:5d} "
              f"v=[{r['v_min']:.2f},{r['v_max']:.2f}] band={r['cubic_in_G410_band']}")


if __name__ == "__main__":
    main()
