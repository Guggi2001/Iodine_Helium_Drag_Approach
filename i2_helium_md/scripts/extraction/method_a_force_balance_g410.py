"""Force-balance (Method A) cross-check of the 18 A Method-B fit (G4.10).

THESIS_FINAL_GOALS.md G4.10 (pre-registered cbefdb7, amended 114c1ca before
any run). The drag is extracted point by point from the TDDFT traces with the
Coulomb force projected on the motion (``extraction.force_balance``), then
fit with F = a v, b v^3, a v + b v^3 on 18 A, 9 A and joint.

Pre-registered gate (FIXED; do not re-tune): the primary 18 A cubic b_A must
lie in the G4.9 b +- 2 sensitivity half-widths (read from the committed
``srho3p5A18A_cubic`` bundle: 2.2286 +- 2*0.2127 -> [1.80, 2.65]).
Agreement -> adopt ``srho3p5A18A_cubic`` for G5; disagreement -> decision.

Primary configuration: committed smoothed reference, SG derivative (order 2,
0.5 ps), projected balance, m_eff, G4.9 E_bind, potential width 14.2 A,
N = 2000 bulk radius, droplet centre at the origin, s_rho = 3.5 A, interior
rho_hat >= 0.9, velocity direction from 0.5 ps-smoothed components.

Recorded, not gating: SG 0.25 / 1.0 ps; raw-|v2| derivative 0.25 / 0.5 /
1.0 ps (amendment 2); scalar and legacy projections; E_bind +- 1 half-width;
droplet radius +- 8 % (amendment 4); 9 A anchored m(t) 21 -> 14 He
(t* = 2.67 ps); the apparent drag outside the surface (rho_hat < 0.1).
The 18 A-alone form comparison cannot separate linear from cubic
(amendment 3: 2 % speed range).

Writes ``data/reference/drag/force_balance_srho3p5A/``: ``verdict.json``,
``fits.csv`` (every variant x case set x form), ``bins_1ps.csv`` (primary
per-ps force terms, both cases), ``trace_{9A,18A}.csv`` (primary, every 10th
sample).

Usage::

    python scripts/extraction/method_a_force_balance_g410.py
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
    savgol_speed_and_derivative,
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
GATE_STEEPNESS_A = 3.5            # s_rho (G2)
INTERIOR_RHO_MIN = 0.9            # G4.10 interior cut
EXTERIOR_RHO_MAX = 0.1            # G4.10 "outside the surface" read
DIRECTION_WINDOW_PS = 0.5         # velocity-direction smoothing (G4.10)
SG_PRIMARY_PS = 0.5               # derivative window, primary (amendment 1)
SG_WINDOWS_PS = (0.25, 0.5, 1.0)  # sensitivity windows
SG_POLYORDER = 2
RADIUS_SPREAD = 0.08              # +- 8 % N -> R convention spread (G1)
T_STAR_9A_PS = 2.67               # anchored m(t) onset = 9 A window start
TRACE_STRIDE = 10                 # export every 10th sample (0.01 ps)
# --- Pre-registered gate (FIXED before any run; do not re-tune) -------------
GATE_HALF_WIDTHS = 2.0
# =============================================================================

DRAG_ROOT = REPO_ROOT / "data" / "reference" / "drag"
G49_BUNDLE = DRAG_ROOT / "18A" / "trajectory_matching_srho3p5A" / "cubic" / "fit_parameters.json"
OUT_DIR = DRAG_ROOT / "force_balance_srho3p5A"


def _base_settings() -> dict:
    """Primary settings; physics constants from the Tier-0 preset + G4.9 bundle."""
    cfg = single_pulse_N2000_18Angst_drag()
    bundle = json.loads(G49_BUNDLE.read_text())
    return {
        "speed_source": "reference",
        "sg_window_ps": SG_PRIMARY_PS,
        "projection": "projected",
        "mass": "m_eff",
        "m_eff_amu": float(cfg.m_eff_amu),
        "binding_energy_eV": float(bundle["effective_binding_energy_I_ion_eV"]),
        "potential_steepness_A": float(cfg.potential_steepness),
        "droplet_radius_A": float(droplet_radius_bulk_angstrom(cfg.single_droplet_size)),
        "gate_steepness_A": GATE_STEEPNESS_A,
        "_bundle": bundle,
    }


def _variants(base: dict) -> dict[str, dict]:
    """Every recorded configuration, keyed by name (``primary`` first)."""
    hw_E = float(base["_bundle"]["effective_binding_energy_err_eV"])
    R = base["droplet_radius_A"]
    v = {"primary": {}}
    for w in SG_WINDOWS_PS:
        if w != SG_PRIMARY_PS:
            v[f"ref_sg{w:g}ps"] = {"sg_window_ps": w}
    for w in SG_WINDOWS_PS:
        v[f"raw_sg{w:g}ps"] = {"speed_source": "raw", "sg_window_ps": w}
    v["scalar"] = {"projection": "scalar"}
    v["legacy_no_well"] = {"projection": "legacy_no_well"}
    v["ebind_minus_hw"] = {"binding_energy_eV": base["binding_energy_eV"] - hw_E}
    v["ebind_plus_hw"] = {"binding_energy_eV": base["binding_energy_eV"] + hw_E}
    v["radius_minus8pct"] = {"droplet_radius_A": R * (1.0 - RADIUS_SPREAD)}
    v["radius_plus8pct"] = {"droplet_radius_A": R * (1.0 + RADIUS_SPREAD)}
    v["mass_anchored_9A"] = {"mass": "anchored"}
    return {k: {**{kk: vv for kk, vv in base.items() if kk != "_bundle"}, **d} for k, d in v.items()}


def _case_trace(case: str, s: dict, geometry, traj, ref):
    """Force-balance trace of one case under settings ``s``."""
    t = np.asarray(ref.time_ps, dtype=float)
    if s["speed_source"] == "reference":
        speed, dvdt = savgol_speed_and_derivative(
            t, ref.speed_Aps, s["sg_window_ps"], smooth_speed=False, polyorder=SG_POLYORDER
        )
    else:
        raw = np.interp(t, traj.time_ps, traj.v2_magnitude_Aps)
        speed, dvdt = savgol_speed_and_derivative(
            t, raw, s["sg_window_ps"], smooth_speed=True, polyorder=SG_POLYORDER
        )
    if s["mass"] == "anchored":
        mass = complex_mass_amu(build_shell_schedule(T_STAR_9A_PS).n_of_t(t))
    else:
        mass = s["m_eff_amu"]
    return force_balance_trace(
        geometry,
        t,
        speed,
        dvdt,
        mass_amu=mass,
        binding_energy_eV=s["binding_energy_eV"],
        potential_steepness_A=s["potential_steepness_A"],
        droplet_radius_A=s["droplet_radius_A"],
        gate_steepness_A=s["gate_steepness_A"],
        projection=s["projection"],
        interior_rho_min=INTERIOR_RHO_MIN,
    )


def _bins(case: str, fb, bin_ps: float = 1.0) -> list[dict]:
    """Per-bin means of the force terms over the whole window (primary)."""
    rows = []
    edges = np.arange(np.floor(fb.time_ps[0]), fb.time_ps[-1] + bin_ps, bin_ps)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (fb.time_ps >= lo) & (fb.time_ps < hi)
        if not m.any():
            continue
        rows.append(
            {
                "case": case,
                "t_lo_ps": lo,
                "t_hi_ps": hi,
                "rho_hat": fb.rho_hat[m].mean(),
                "speed_Aps": fb.speed_Aps[m].mean(),
                "cos_theta": fb.cos_theta[m].mean(),
                "F_coulomb": fb.F_coulomb[m].mean(),
                "F_coulomb_along": fb.F_coulomb_along[m].mean(),
                "F_well_along": fb.F_well_along[m].mean(),
                "F_inertia": fb.F_inertia[m].mean(),
                "F_drag_apparent": fb.F_drag_apparent[m].mean(),
                "F_coulomb_offaxis": (fb.F_coulomb[m] * (1.0 - fb.cos_theta[m])).mean(),
            }
        )
    return rows


def main() -> None:
    base = _base_settings()
    bundle = base["_bundle"]
    variants = _variants(base)
    data = {}
    for case in CASES:
        traj = load_hedft_trajectory(REPO_ROOT / "data" / "reference" / f"{case}_All_Data.csv")
        ref = load_smoothed_speed_reference(
            DRAG_ROOT / case / "velocity_smoothed" / "cleaned_data_long.csv"
        )
        data[case] = (traj, ref, balance_geometry(traj, direction_window_ps=DIRECTION_WINDOW_PS))

    fit_rows, traces = [], {}
    for name, s in variants.items():
        cases = ("9A",) if name == "mass_anchored_9A" else CASES
        per_case = {c: _case_trace(c, s, data[c][2], data[c][0], data[c][1]) for c in cases}
        traces[name] = per_case
        sets = {c: (c,) for c in cases}
        if len(cases) > 1:
            sets["joint"] = cases
        for set_name, members in sets.items():
            v = np.concatenate([per_case[c].speed_Aps[per_case[c].interior] for c in members])
            F = np.concatenate([per_case[c].F_drag_bulk[per_case[c].interior] for c in members])
            for form in FORMS:
                f = fit_drag_form(v, F, form)
                fit_rows.append({"variant": name, "case_set": set_name, **asdict(f)})

    # --- Pre-registered gate ------------------------------------------------
    b49 = float(bundle["b"])
    hw_b = float(bundle["b_err"])
    band = (b49 - GATE_HALF_WIDTHS * hw_b, b49 + GATE_HALF_WIDTHS * hw_b)
    prim = {(r["case_set"], r["form"]): r for r in fit_rows if r["variant"] == "primary"}
    b_A = prim[("18A", "cubic")]["b"]
    agree = band[0] <= b_A <= band[1]

    # --- Exterior apparent drag (primary) ------------------------------------
    exterior = {}
    for case in CASES:
        fb = traces["primary"][case]
        m = fb.rho_hat < EXTERIOR_RHO_MAX
        exterior[case] = {
            "n_points": int(m.sum()),
            "t_range_ps": [float(fb.time_ps[m].min()), float(fb.time_ps[m].max())] if m.any() else None,
            "mean_F_drag_apparent": float(fb.F_drag_apparent[m].mean()) if m.any() else None,
            "mean_speed_Aps": float(fb.speed_Aps[m].mean()) if m.any() else None,
            "mean_ratio_to_G49_bulk_cubic": float(
                (fb.F_drag_apparent[m] / (b49 * fb.speed_Aps[m] ** 3)).mean()
            ) if m.any() else None,
        }

    # The nested a v + b v^3 always has the lowest RMSE, so the form read is
    # the per-form RMSE plus the sign/size of the linear_cubic coefficients.
    rmse_by_form = {
        set_name: {r["form"]: r["rmse_force"] for r in fit_rows
                   if r["variant"] == "primary" and r["case_set"] == set_name}
        for set_name in ("18A", "9A", "joint")
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "fits.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fit_rows[0]))
        w.writeheader()
        w.writerows(fit_rows)
    bins = [row for c in CASES for row in _bins(c, traces["primary"][c])]
    with (OUT_DIR / "bins_1ps.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bins[0]))
        w.writeheader()
        w.writerows(bins)
    cols = ("time_ps", "speed_Aps", "dvdt_Aps2", "mass_amu", "depth_A", "rho_hat",
            "cos_theta", "F_coulomb", "F_coulomb_along", "F_well_along", "F_inertia",
            "F_drag_apparent", "interior", "F_drag_bulk")
    for case in CASES:
        fb = traces["primary"][case]
        arr = np.column_stack([np.asarray(getattr(fb, c), dtype=float) for c in cols])[::TRACE_STRIDE]
        np.savetxt(OUT_DIR / f"trace_{case}.csv", arr, delimiter=",", header=",".join(cols),
                   comments="", fmt="%.9g")

    verdict = {
        "date": "2026-10-05",
        "goal": "THESIS_FINAL_GOALS.md G4.10: force-balance cross-check of the 18 A fit",
        "preregistration": "G4.10 committed cbefdb7; amendments 114c1ca (before any run)",
        "primary_settings": {k: v for k, v in variants["primary"].items()},
        "settings_fixed": {
            "interior_rho_min": INTERIOR_RHO_MIN,
            "direction_window_ps": DIRECTION_WINDOW_PS,
            "sg_polyorder": SG_POLYORDER,
            "droplet_centre": "origin (I2 centre of mass)",
            "t_star_9A_anchored_ps": T_STAR_9A_PS,
        },
        "gate": {
            "G49_b": b49,
            "G49_b_halfwidth": hw_b,
            "band": list(band),
            "b_A_18A_cubic_primary": b_A,
            "agree": bool(agree),
            "consequence": (
                "adopt srho3p5A18A_cubic for G5" if agree
                else "no automatic adoption; open as a decision"
            ),
        },
        "primary_fits": {f"{k[0]}/{k[1]}": v for k, v in prim.items()},
        "rmse_by_form_primary": rmse_by_form,
        "form_identifiable_18A": False,
        "exterior_rho_below_0p1": exterior,
        "anchor": (
            "legacy_no_well + old spline derivative reproduces the committed "
            "Method-A linear_and_cubic fits (18 A repo reference; 9 A old local "
            "cleaned_data.csv) to <= 5e-6 relative (9648.5 rounding)"
        ),
    }
    (OUT_DIR / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict["gate"], indent=2))
    for r in fit_rows:
        print(f"{r['variant']:18s} {r['case_set']:6s} {r['form']:13s} a={r['a']:8.3f} "
              f"b={r['b']:7.4f} rmse={r['rmse_force']:7.3f} n={r['n_points']:5d} "
              f"v=[{r['v_min']:.2f},{r['v_max']:.2f}]")


if __name__ == "__main__":
    main()
