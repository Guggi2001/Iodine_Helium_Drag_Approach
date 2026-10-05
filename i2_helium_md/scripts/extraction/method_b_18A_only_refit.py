"""Method-B 18 A-only refit at s_rho = 3.5 A, 9 A held out (G4.9).

THESIS_FINAL_GOALS.md G4.9. The 9 A TDDFT reference is non-radial (the
compared ion stalls, oscillates and exits sideways); a radial MD can match it
only by absorbing the missing dimension into the drag / binding (METHOD_B §5).
9 A therefore leaves the fit: the drag law is calibrated on 18 A alone and 9 A
is run only as an escape check (its RMSE is descriptive, never evidence).

Same engine, bounds, starts, anchors (18 A Method-A bundle), N, seed and
reference files as G4 (``method_b_srho3p5_refit.py``); single case via
``fit_shared_trajectory_matching(..., allow_single_case=True)``.

Pre-registered (G4.9, committed fc4ece9 before any run; do NOT re-tune):

* bands per family: 18 A RMSE <= 0.19 A/ps, escape fraction == 1;
* nested tests and winner rule unchanged from G4.2 (``select_winner``):
  T_a0, T_b0 <= 0.005 A/ps;
* step 3 (9 A escape check of the winner, no refit) is recorded, not gating.

Modes (USER SETTINGS ``MODE``):

* ``"oracle"`` -- single-case evaluations at the committed G4 ``shared_cubic``
  point (3.5 A) and the committed 14.2 A control point must reproduce the
  committed per-case RMSEs bit-exactly, both cases. Writes nothing.
* ``"residual"`` -- step 1: 9 A |v2| residual vs time at both committed points,
  next to the raw HeDFT transverse velocity of atom 2. Writes
  ``step1_9A_residual_vs_transverse.{csv,json}``.
* ``"dry_run"`` -- one 18 A single-case evaluation twice: determinism + timing.
* ``"fit"`` -- step 2: fit each family in ``FAMILIES`` on 18 A at 3.5 A; writes
  ``<family>/fit_summary.json`` and, unless trapped, the loadable
  ``<family>/fit_parameters.json``.
* ``"control"`` -- ``cubic`` on 18 A at the legacy 14.2 A gate; summary only.
* ``"verdict"`` -- bands + winner rule, then step 3: the winner's (and the
  control's) 9 A escape check (full-window RMSE recorded, descriptive only).
  Writes ``verdict.json``.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_18A_only_refit.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np

# =============================================================================
# USER SETTINGS
# =============================================================================
MODE = "oracle"   # "oracle" | "residual" | "dry_run" | "fit" | "control" | "verdict"
FAMILIES = ("cubic", "linear_cubic", "linear")

# --- Pre-registered G4.9 bands (FIXED before any run; do not re-tune) --------
RMSE_18A_MAX_APS = 0.19
# Escape fraction must be exactly 1.0. Nested tests: G4.2 (select_winner).

# Step 1 (descriptive): late-phase split [ps] (METHOD_B §3.5 dated the
# start of the transverse drift to ~6 ps; step 1 found it at ~2 ps).
LATE_SPLIT_PS = 6.0

# Step 1 (descriptive): reference samples with |v_perp| / |v| >= this count as
# the transverse (stall / sideways) phase.
TRANSVERSE_FRACTION_MIN = 0.5

HELD_OUT_AXIS = "9A (escape check only, non-gating; G4.9 step 3)"

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    SHARED_CUBIC,
    SHARED_LINEAR,
    SHARED_LINEAR_CUBIC,
    evaluate_objective,
    fit_shared_trajectory_matching,
    joint_sensitivity_halfwidths,
    write_shared_fit_parameters,
)
from i2_helium_md.postprocess.compare_trajectories import (  # noqa: E402
    compare_speed_to_reference,
)
from i2_helium_md.postprocess.hedft_loader import load_hedft_trajectory  # noqa: E402
from i2_helium_md.presets import DRAG_BUNDLE_CATALOG, REFERENCE_DRAG_ROOT  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from scripts.extraction.method_b_srho3p5_refit import (  # noqa: E402
    CONTROL_DIR as G4_CONTROL_DIR,
    DENSITY_WIDTH_SOURCE,
    GATE_OVERRIDES,
    MAXFEV_PER_START,
    OUT_ROOT as G4_OUT_ROOT,
    UNCERTAINTY_MODEL,
    _dump_json,
    _load_json,
    _stamp,
    anchors,
    build_setups,
    select_winner,
)

BUNDLE_TAG = "srho3p5A18A"
OUT_ROOT = REFERENCE_DRAG_ROOT / "18A" / "trajectory_matching_srho3p5A"
CONTROL_OUT = OUT_ROOT / "control_gate14p2A_cubic"
# Short family name (folder / bundle id) -> engine variant (G4 names).
FAMILY_VARIANT = {
    "linear": SHARED_LINEAR,
    "cubic": SHARED_CUBIC,
    "linear_cubic": SHARED_LINEAR_CUBIC,
}
HEDFT_9A_RAW = REFERENCE_DRAG_ROOT.parent / "9A_All_Data.csv"


def bundle_id(family: str) -> str:
    """Catalog id of a family's bundle (must exist in DRAG_BUNDLE_CATALOG)."""
    bid = f"{BUNDLE_TAG}_{family}"
    if bid not in DRAG_BUNDLE_CATALOG:
        raise ValueError(f"bundle id {bid!r} is not in DRAG_BUNDLE_CATALOG")
    if (REFERENCE_DRAG_ROOT / DRAG_BUNDLE_CATALOG[bid]).resolve() != (
        OUT_ROOT / family
    ).resolve():
        raise ValueError(f"catalog maps {bid!r} elsewhere than {OUT_ROOT / family}")
    return bid


def band_checks_18a(r18: dict) -> dict:
    """Apply the pre-registered G4.9 bands to one family's 18 A score."""
    checks = {
        "RMSE_18A": {
            "value_Aps": r18["rmse_Aps"], "band_max_Aps": RMSE_18A_MAX_APS,
            "pass": bool(r18["rmse_Aps"] <= RMSE_18A_MAX_APS),
        },
        "ESCAPE_18A": {
            "value": r18["escape_fraction"], "band": 1.0,
            "pass": bool(r18["escape_fraction"] == 1.0),
        },
    }
    checks["all_pass"] = all(c["pass"] for c in checks.values())
    return checks


def residual_shares(
    t_ps: np.ndarray,
    residual_14_Aps: np.ndarray,
    residual_35_Aps: np.ndarray,
    transverse_fraction: np.ndarray,
    *,
    fraction_min: float,
    t_split_ps: float,
) -> dict:
    """Where does the extra squared residual of the 3.5 A gate sit? (step 1)

    ``d = r35**2 - r14**2`` per reference sample [(A/ps)^2]. Returns, for the
    transverse mask (``transverse_fraction >= fraction_min``) and the late mask
    (``t >= t_split_ps``), the share of ``sum(d)`` inside the mask next to the
    share of samples inside it (a share well above the sample share means the
    extra residual concentrates there). Non-finite samples are dropped.

    Raises
    ------
    ValueError
        On unequal lengths or a non-positive total ``sum(d)`` (no extra
        residual to locate).
    """
    arrays = [np.asarray(a, dtype=float) for a in
              (t_ps, residual_14_Aps, residual_35_Aps, transverse_fraction)]
    if len({a.shape for a in arrays}) != 1:
        raise ValueError("all inputs must have the same shape")
    t, r14, r35, frac = arrays
    ok = np.isfinite(t) & np.isfinite(r14) & np.isfinite(r35) & np.isfinite(frac)
    t, r14, r35, frac = t[ok], r14[ok], r35[ok], frac[ok]
    d = r35 ** 2 - r14 ** 2
    total = float(d.sum())
    if not total > 0.0:
        raise ValueError(f"no extra residual to locate: sum(d) = {total}")
    out = {"sum_d_Aps2": total, "n_samples": int(t.size)}
    for name, mask in (("transverse", frac >= fraction_min), ("late", t >= t_split_ps)):
        out[name] = {
            "share_of_extra_sq_residual": float(d[mask].sum() / total),
            "share_of_samples": float(mask.mean()),
        }
    return out


def _single(setups: dict, case: str) -> dict:
    return {case: setups[case]}


def _committed_points() -> dict:
    """The committed G4 3.5 A shared_cubic point and the 14.2 A control point."""
    g4 = _load_json(G4_OUT_ROOT / SHARED_CUBIC / "fit_summary.json")
    ctrl = _load_json(G4_CONTROL_DIR / "fit_summary.json")
    return {
        "3p5A": {"theta": (g4["a"], g4["b"], g4["e_bind_eV"]), "per_case": g4["per_case"],
                 "gate": GATE_OVERRIDES},
        "14p2A": {"theta": (ctrl["a"], ctrl["b"], ctrl["e_bind_eV"]),
                  "per_case": ctrl["per_case"], "gate": None},
    }


def run_oracle() -> None:
    for label, pt in _committed_points().items():
        setups = build_setups(gate_overrides=pt["gate"])
        for case, setup in setups.items():
            r = evaluate_objective(pt["theta"], setup)
            committed = pt["per_case"][case]["rmse_Aps"]
            ok = r.rmse_Aps == committed
            print(f"  {label} {case}: rmse {r.rmse_Aps!r} vs committed {committed!r} "
                  f"-> {'BIT-EXACT' if ok else 'MISMATCH'}")
            if not ok:
                raise AssertionError(f"oracle failed at {label} {case}")
    print("  ORACLE PASS: single-case evaluations reproduce the committed RMSEs.")


def _md_speed_on_ref(setup, theta) -> tuple:
    """Re-run one committed point and return the comparison on the 9 A grid."""
    captured = {}

    def run_fn(cfg, neutral):  # capture the ion run the objective scores
        captured["ion"] = run_ion_propagation(cfg, neutral)
        return captured["ion"]

    r = evaluate_objective(theta, setup, run_fn=run_fn)
    comp = compare_speed_to_reference(
        captured["ion"], atom="I2", t_ref_ps=setup.t_ref_ps,
        ref_speed_Aps=setup.ref_speed_Aps, window=setup.window,
    )
    return comp, r


def run_residual() -> None:
    pts = _committed_points()
    comps = {}
    for label, pt in pts.items():
        setup = build_setups(gate_overrides=pt["gate"])["9A"]
        comp, r = _md_speed_on_ref(setup, pt["theta"])
        committed = pt["per_case"]["9A"]["rmse_Aps"]
        if r.rmse_Aps != committed:
            raise AssertionError(f"{label}: rmse {r.rmse_Aps!r} != committed {committed!r}")
        comps[label] = comp
    t = comps["3p5A"].t_overlap_ps
    if not np.array_equal(t, comps["14p2A"].t_overlap_ps):
        raise AssertionError("the two comparisons use different reference grids")
    ref = comps["3p5A"].hedft_on_overlap
    r35 = comps["3p5A"].md_on_hedft_grid - ref
    r14 = comps["14p2A"].md_on_hedft_grid - ref

    raw = load_hedft_trajectory(HEDFT_9A_RAW)
    v_perp_raw = np.hypot(raw.v2_x_Aps, raw.v2_y_Aps)
    v_perp = np.interp(t, raw.time_ps, v_perp_raw)
    v_z = np.interp(t, raw.time_ps, raw.v2_z_Aps)
    v_mag = np.interp(t, raw.time_ps, raw.v2_magnitude_Aps)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(v_mag > 0.0, v_perp / v_mag, np.nan)

    shares = residual_shares(t, r14, r35, frac, fraction_min=TRANSVERSE_FRACTION_MIN,
                             t_split_ps=LATE_SPLIT_PS)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_ROOT / "step1_9A_residual_vs_transverse.csv"
    np.savetxt(
        csv_path,
        np.column_stack([t, ref, comps["14p2A"].md_on_hedft_grid,
                         comps["3p5A"].md_on_hedft_grid, r14, r35,
                         v_perp, v_z, v_mag, frac]),
        delimiter=",", fmt="%.8g", comments="",
        header=("Time_ps,ref_v2_smoothed_Aps,md_v2_gate14p2A_Aps,md_v2_gate3p5A_Aps,"
                "resid_gate14p2A_Aps,resid_gate3p5A_Aps,hedft_v2_perp_Aps,"
                "hedft_v2_z_Aps,hedft_v2_mag_Aps,hedft_transverse_fraction"),
    )
    print(f"  wrote {csv_path}")
    summary = {
        **_stamp(),
        "goal": "THESIS_FINAL_GOALS.md G4.9 step 1 (no new fits)",
        "points": {k: {"a": v["theta"][0], "b": v["theta"][1], "e_bind_eV": v["theta"][2]}
                   for k, v in pts.items()},
        "rmse_9A_Aps": {k: float(c.rmse) for k, c in comps.items()},
        "transverse_source": str(HEDFT_9A_RAW.relative_to(PROJECT_ROOT).as_posix())
        + " (raw HeDFT, atom 2: |v_perp| = hypot(V2_x, V2_y))",
        "transverse_fraction_min": TRANSVERSE_FRACTION_MIN,
        "late_split_ps": LATE_SPLIT_PS,
        "shares": shares,
    }
    _dump_json(summary, OUT_ROOT / "step1_9A_residual_vs_transverse.json")
    for name in ("transverse", "late"):
        s = shares[name]
        print(f"  {name}: {100 * s['share_of_extra_sq_residual']:.1f} % of the extra "
              f"squared residual on {100 * s['share_of_samples']:.1f} % of the samples")


def dry_run(anch) -> None:
    setup = build_setups()["18A"]
    theta = (anch[0], anch[1], 0.2)
    t0 = time.perf_counter()
    r1 = evaluate_objective(theta, setup)
    t1 = time.perf_counter()
    r2 = evaluate_objective(theta, setup)
    print(f"  18A single-case eval: {t1 - t0:.2f} s; rmse {r1.rmse_Aps:.4f}, "
          f"escape {r1.escape_fraction:.3f}")
    if r1.objective != r2.objective:
        raise AssertionError("single-case objective is NOT deterministic")
    print("    repeated evaluation bitwise identical -> deterministic OK")


def _summary(fit, setup, family: str | None) -> dict:
    r = fit.best.per_case["18A"]
    per = {"rmse_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction,
           "objective_Aps": r.objective, "num_overlap_points": r.num_overlap_points}
    return {
        **_stamp(),
        "family": family, "variant": fit.variant,
        "bundle_id": bundle_id(family) if family else None,
        "calibration_cases": fit.cases,
        "a": fit.best.a, "b": fit.best.b, "e_bind_eV": fit.best.e_bind_eV["18A"],
        "objective_Aps": fit.best.objective,
        "per_case": {"18A": per},
        "bands": band_checks_18a(per),
        "n_evaluations": fit.n_evaluations, "converged": fit.converged,
        "start_minima": [{"a": s.a, "b": s.b, "e_bind_eV": s.e_bind_eV["18A"],
                          "objective_Aps": s.objective} for s in fit.starts],
        "anchors": {"a0": fit.anchors[0], "b0": fit.anchors[1]},
        "gate": {"drag_spatial_gate": fit.drag_spatial_gate,
                 "drag_gate_steepness_angstrom": fit.gate_steepness_angstrom,
                 "potential_steepness_angstrom": fit.potential_steepness_angstrom},
        "reference_file": str(setup.reference_path.relative_to(PROJECT_ROOT).as_posix()),
        "n_molecules": fit.n_molecules, "seed": fit.seed, "meff_amu": fit.meff_amu,
    }


def _fit(setups, family, anch):
    t0 = time.perf_counter()
    fit = fit_shared_trajectory_matching(
        setups, variant=FAMILY_VARIANT[family], anchors=anch,
        maxfev_per_start=MAXFEV_PER_START, allow_single_case=True,
    )
    b = fit.best
    print(f"    {family}: {time.perf_counter() - t0:.0f} s, {fit.n_evaluations} evals; "
          f"a={b.a:.4f}, b={b.b:.4f}, E_bind={b.e_bind_eV['18A']:.4f} eV, "
          f"rmse={b.per_case['18A'].rmse_Aps:.4f}")
    return fit


def run_fits(anch) -> None:
    setups = _single(build_setups(), "18A")
    for family in FAMILIES:
        if family not in FAMILY_VARIANT:
            raise ValueError(f"unknown family {family!r}; expected {sorted(FAMILY_VARIANT)}")
        fit = _fit(setups, family, anch)
        summary = _summary(fit, setups["18A"], family)
        trapped = fit.best.per_case["18A"].penalty_Aps != 0.0
        if not trapped:
            hw_a, hw_b, hw_e = joint_sensitivity_halfwidths(setups, fit.best)
            fit.a_err, fit.b_err, fit.e_bind_err_eV = hw_a, hw_b, hw_e
            fit.uncertainty_model = UNCERTAINTY_MODEL
            summary["sensitivity_halfwidths"] = {"a": hw_a, "b": hw_b, "e_bind_eV": hw_e}
            write_shared_fit_parameters(
                fit, OUT_ROOT / family, stage="single_case_18A",
                bundle_id=bundle_id(family),
                density_width_source=DENSITY_WIDTH_SOURCE, held_out_axis=HELD_OUT_AXIS,
            )
        summary["bundle_written"] = not trapped
        _dump_json(summary, OUT_ROOT / family / "fit_summary.json")


def run_control(anch) -> None:
    setups = _single(build_setups(gate_overrides=None), "18A")
    fit = _fit(setups, "cubic", anch)
    summary = _summary(fit, setups["18A"], None)
    summary["purpose"] = ("G4.9 control: cubic on 18 A at the legacy 14.2 A gate, same "
                          "reference; gate effect = 3.5 A cubic - control. Not a bundle.")
    summary["bundle_written"] = False
    _dump_json(summary, CONTROL_OUT / "fit_summary.json")


def check_9a(theta, gate) -> dict:
    """Step 3: run a fixed 18 A-calibrated point on 9 A (no refit).

    The escape fraction is the check; the full-window RMSE is descriptive only
    (a radial model on the non-radial 9 A trace, G4.9 amendment).
    """
    setup = build_setups(gate_overrides=gate)["9A"]
    r = evaluate_objective(theta, setup)
    return {"escape_fraction": r.escape_fraction,
            "rmse_full_window_Aps_descriptive_only": r.rmse_Aps,
            "window_ps": list(setup.window)}


def run_verdict() -> None:
    summaries = {f: _load_json(OUT_ROOT / f / "fit_summary.json") for f in FAMILY_VARIANT}
    variants = {f: FAMILY_VARIANT[f] for f in summaries}
    objectives = {variants[f]: float(s["objective_Aps"]) for f, s in summaries.items()}
    passes = {variants[f]: bool(s["bands"]["all_pass"]) for f, s in summaries.items()}
    selection = select_winner(objectives, passes)
    winner_variant = selection["winner"]
    winner = next((f for f, v in variants.items() if v == winner_variant), None)
    print(f"  T_a0 = {selection['T_a0_Aps']:+.6f}, T_b0 = {selection['T_b0_Aps']:+.6f}; "
          f"bands {passes}; WINNER {winner} ({selection['reason']})")

    check = {}
    if winner is not None:
        w = summaries[winner]
        check["winner"] = check_9a((w["a"], w["b"], w["e_bind_eV"]), GATE_OVERRIDES)
    ctrl = None
    if (CONTROL_OUT / "fit_summary.json").is_file():
        ctrl = _load_json(CONTROL_OUT / "fit_summary.json")
        check["control_14p2A"] = check_9a((ctrl["a"], ctrl["b"], ctrl["e_bind_eV"]), None)
    for k, v in check.items():
        print(f"  9A check ({k}): escape {v['escape_fraction']}, full-window rmse "
              f"{v['rmse_full_window_Aps_descriptive_only']:.4f} A/ps (descriptive)")

    verdict = {
        **_stamp(),
        "goal": "THESIS_FINAL_GOALS.md G4.9: 18 A-only refit at s_rho = 3.5 A",
        "preregistration": "THESIS_FINAL_GOALS.md G4.9, committed fc4ece9 before any run",
        "thresholds": {"RMSE_18A_MAX_Aps": RMSE_18A_MAX_APS, "escape_band": 1.0,
                       "T_a0_Aps": selection["T_a0_threshold_Aps"],
                       "T_b0_Aps": selection["T_b0_threshold_Aps"]},
        "families": {f: {k: s[k] for k in ("bundle_id", "a", "b", "e_bind_eV",
                                             "objective_Aps", "per_case", "bands",
                                             "converged", "bundle_written")}
                     for f, s in summaries.items()},
        "selection": {**selection, "winner_family": winner},
        "control_14p2A": None if ctrl is None else {
            k: ctrl[k] for k in ("a", "b", "e_bind_eV", "objective_Aps")},
        "gate_effect_cubic_minus_control": None if ctrl is None else {
            "db": summaries["cubic"]["b"] - ctrl["b"],
            "de_bind_eV": summaries["cubic"]["e_bind_eV"] - ctrl["e_bind_eV"]},
        "check_9A_non_gating": check,
        "check_9A_note": ("escape fraction is the check; the RMSE is descriptive only "
                          "(9 A radial only to ~2 ps < t*; G4.9 amendment)"),
    }
    _dump_json(verdict, OUT_ROOT / "verdict.json")


def main(mode: str = MODE) -> int:
    print(f"===== Method-B 18A-only refit (G4.9), MODE={mode} =====")
    if mode == "oracle":
        run_oracle()
    elif mode == "residual":
        run_residual()
    elif mode == "verdict":
        run_verdict()
    elif mode in ("dry_run", "fit", "control"):
        anch = anchors()
        print(f"  anchors (18A Method-A only): a0={anch[0]:.4f}, b0={anch[1]:.4f}")
        {"dry_run": dry_run, "fit": run_fits, "control": run_control}[mode](anch)
    else:
        raise ValueError(f"unknown MODE {mode!r}")
    return 0


if __name__ == "__main__":
    # Optional command-line override of MODE: ``python ... fit``.
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else MODE))
