"""G4.12 steps 1-2: the TDDFT-tied linear + quadratic ridge and the 9 A exit bound.

THESIS_FINAL_GOALS.md G4.12 (frozen beb5440; do NOT re-tune). 18 A fixes one
number, the drag at v ~ 3 A/ps, so ``F = a v + c v^2`` along a one-parameter
18 A ridge fits it. This driver maps that ridge at the s_rho = 3.5 A gate and
bounds each point by the radial 9 A exit speed. Step 3 (the twin) is NOT here.

Units: a [amu/ps], c [amu/A]; each term amu*A/ps^2; gamma = a + c|v| [amu/ps].

Modes (USER SETTINGS ``MODE``; argv override ``python ... ridge 0 0.4``):

* ``"oracle"`` -- the lq form reproduces (i) the G4.9 linear point's 18 A RMSE
  (a 23.154, c 0) and (ii) the G4.11 quadratic point's (a 0, c 7.224),
  bit-exactly, at the 3.5 A gate. Writes nothing.
* ``"ridge"`` -- fixed-c 18 A fits (``lq_fixed_c``: a >= 0, E_bind free) for
  each c in ``C_GRID`` (or argv); one JSON per point in ``RIDGE_DIR``.
* ``"assemble"`` -- collect the per-point JSONs into ``ridge.csv``; oracle
  (iii): the c = 0 point is the G4.9 linear point within |da| <= 0.2,
  |dE| <= 0.002 eV.
* ``"anchor"`` -- score the user's existing 9 A Tier-0 runs of the 18 A-only
  linear / quadratic / cubic bundles (no MD) -> ``bound_anchor.csv``.
* ``"bound"`` -- for each assembled ridge point (or argv c values): radial
  9 A Tier-0 run (N = 50, 30 ps, the user's gen_tier0_runs settings), then the
  exit-bound score -> ``bound.csv``.

G4.14 Step 0a (frozen 7726e43): the same steps 1-2, unchanged, for the fine
points ``G414_C_FINE``. The G4.12 artifacts are never rewritten.

* ``"oracle414"`` -- the c = 2.0 fixed-c fit rerun in memory reproduces the
  committed ``c2p000.json`` (a, E_bind, RMSE) bit-exactly. Writes nothing.
* ``"ridge414"`` -- fixed-c fits for ``G414_C_FINE`` (or argv); refuses to
  overwrite an existing point JSON.
* ``"bound414"`` -- the 9 A exit bound for ``G414_C_FINE`` (or argv).
* ``"assemble414"`` -- the four points (bands + bound) -> ``ridge_g414.csv``.

G11.4 step 2 (THESIS_FINAL_GOALS G11.2, agreed 2026-10-09): steps 1-2 again at
the 18 A shell mass m = I+ + 21 He (production ``m(n)`` convention,
``G11_M_EFF_AMU``), written to ``G11_RIDGE_DIR`` / ``G11_BOUND_DIR``; the m = 19
artifacts stay the record and are never rewritten.

* ``"oracle_m21"`` -- the mass-parameter path run at m = 19 reproduces the
  committed c = 1.6 point's 18 A RMSE bit-exactly. Writes nothing.
* ``"check_m21"`` -- the committed production points (``G11_CHECK_POINTS``:
  B2, A5) re-scored unchanged on 18 A at m = 21 (bands; no refit).
* ``"ridge_m21"`` -- fixed-c fits at m = 21 for ``C_GRID`` (or argv).
* ``"bound_m21"`` -- the 9 A exit bound at m = 21 for argv labels (floats =
  refit ridge points, names = check points; default all).
* ``"assemble_m21"`` -- ridge + checks + bound -> ``ridge_m21.csv``.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_g412_ridge.py
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import time

import numpy as np

# =============================================================================
# USER SETTINGS
# =============================================================================
MODE = "oracle"   # "oracle" | "ridge" | "assemble" | "anchor" | "bound"

# Pre-registered ridge grid [amu/A] (G4.12 step 1; FIXED).
C_GRID = (0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.224)
# G4.14 Step 0a fine points around c = 2.0 [amu/A] (frozen 7726e43).
G414_C_FINE = (1.8, 1.9, 2.1, 2.2)
G414_ORACLE_C = 2.0

# --- Pre-registered bands / constants (FIXED; do not re-tune) -----------------
RMSE_18A_MAX_APS = 0.19            # G4.9 band; escape fraction must be 1.0
SEAM_DA_MAX = 0.2                  # oracle (iii) [amu/ps]
SEAM_DE_MAX_EV = 0.002             # oracle (iii) [eV]
EXIT_WINDOW_PS = (13.0, 14.0)      # step-2 exit-speed window
AT_BOUND_APS = 0.05                # |margin| flag
RHO_INSIDE_MAX = 0.1               # HeDFT atom "still inside" stop threshold
H405_PLATEAU = 418.5               # amu*A/ps^2, Phi = F(9.7) / plateau
FORCE_PROBE_SPEEDS_APS = (3.0, 5.5, 9.7)

MAXFEV_PER_START = 400             # G4.9 / G4.11 comparability
REFERENCE_FILENAME = "cleaned_data_long.csv"
A0_FIXED_C = 14.5556               # normalization anchor (= LOCKED §10.4.1 a0)

# Step-2 radial runs: the user's gen_tier0_runs.py settings.
BOUND_N = 50
BOUND_ION_TIME_PS = 30.0
BOUND_DT_ION_PS = 0.01
BOUND_VARIANT = "srho3p5A18A_linear_quadratic"   # lq bundle; {a, c, E} overridden
ANCHOR_VARIANTS = ("srho3p5A18A_linear", "srho3p5A18A_quadratic", "srho3p5A18A_cubic")

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    LQ_FIXED_C,
    build_case_setup,
    evaluate_form_objective,
    fit_form_trajectory_matching,
)
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU  # noqa: E402
from i2_helium_md.physics.drag import LINEAR_QUADRATIC  # noqa: E402
from i2_helium_md.physics.helium_density import rho_he_ratio  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402
from scripts.extraction.form_phase_common import (  # noqa: E402
    DT_ION_PS,
    ION_TIME_PS,
    N,
    SEED,
    dump_json,
    load_json,
    stamp,
)
from scripts.extraction.method_b_srho3p5_refit import GATE_OVERRIDES  # noqa: E402

CASE = "18A"
TM_ROOT = REFERENCE_DRAG_ROOT / CASE / "trajectory_matching_srho3p5A"
RIDGE_DIR = TM_ROOT / "g412_ridge"
BOUND_DIR = REFERENCE_DRAG_ROOT / "9A" / "g412_exit_bound"
G49_LINEAR_SUMMARY = TM_ROOT / "linear" / "fit_summary.json"
G411_QUADRATIC_SUMMARY = TM_ROOT / "quadratic" / "fit_summary.json"
HEDFT_9A_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
GATE_WIDTH_A = float(GATE_OVERRIDES["drag_gate_steepness"])

# G11.2: extraction mass = I+ + 21 He under the production m(n) convention
# (constants.py), so the fit mass equals production's n0 = 21 start exactly.
G11_N_HE = 21
G11_M_EFF_AMU = MASS_I_ION_AMU + G11_N_HE * MASS_HE_AMU   # 210.9546 amu
G11_RIDGE_DIR = TM_ROOT / "g11_ridge_m21"
G11_BOUND_DIR = REFERENCE_DRAG_ROOT / "9A" / "g11_exit_bound_m21"
G11_RIDGE_CSV = G11_RIDGE_DIR / "ridge_m21.csv"
# Committed production points re-scored unchanged at m = 21 (label -> m = 19
# ridge point; THESIS_FINAL_GOALS G5 Outcome: B2 = c 1.6, A5 = c 2.0).
G11_CHECK_POINTS = {"B2": 1.6, "A5": 2.0}


# =============================================================================
# Pure helpers (unit-tested in tests/test_g412_ridge.py)
# =============================================================================
def lq_force(a: float, c: float, v_Aps) -> np.ndarray:
    """Ungated bulk law ``F = a v + c v^2`` [amu*A/ps^2] at speeds ``v`` [A/ps]."""
    v = np.asarray(v_Aps, dtype=float)
    return a * v + c * v * v


def point_label(c: float) -> str:
    """File/run tag of a ridge point, e.g. ``c0p400`` (3 decimals, no dot)."""
    return f"c{c:.3f}".replace(".", "p")


def window_mean(t_ps, y, window) -> float:
    """Mean of ``y`` over samples with ``t`` in the closed ``window`` [ps]."""
    t = np.asarray(t_ps, dtype=float)
    mask = (t >= window[0]) & (t <= window[1])
    if not mask.any():
        raise ValueError(f"no samples in window {window} (t in [{t[0]}, {t[-1]}])")
    return float(np.mean(np.asarray(y, dtype=float)[mask]))


def transverse_fraction(pos, vel) -> np.ndarray:
    """``|v_perp| / |v|`` per sample, perp to the radial direction ``pos/|pos|``.

    ``pos``, ``vel``: shape (T, 3), positions [A] relative to the droplet
    centre and velocities [A/ps]. Returns shape (T,).
    """
    pos = np.asarray(pos, dtype=float)
    vel = np.asarray(vel, dtype=float)
    r_hat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    v_rad = np.sum(vel * r_hat, axis=1, keepdims=True) * r_hat
    return np.linalg.norm(vel - v_rad, axis=1) / np.linalg.norm(vel, axis=1)


def rmse_on_grid(t_md, y_md, t_ref, y_ref, window) -> float:
    """RMSE of ``y_md`` vs ``y_ref`` (linearly interpolated onto ``t_md``) in ``window``."""
    t_md = np.asarray(t_md, dtype=float)
    mask = (t_md >= window[0]) & (t_md <= window[1])
    ref = np.interp(t_md[mask], np.asarray(t_ref, float), np.asarray(y_ref, float))
    return float(np.sqrt(np.mean((np.asarray(y_md, float)[mask] - ref) ** 2)))


def exit_bound(v_exit_md: float, v_exit_atom1: float) -> dict:
    """G4.12 step-2 bound: pass iff ``v_exit_md >= v_exit_atom1``; margin flagged."""
    margin = float(v_exit_md - v_exit_atom1)
    return {"margin_Aps": margin, "pass": bool(margin >= 0.0),
            "at_bound": bool(abs(margin) <= AT_BOUND_APS)}


# =============================================================================
# Step 1: the 18 A ridge
# =============================================================================
def _setup(mass_amu: float | None = None):
    """18 A case setup at the 3.5 A gate; ``mass_amu`` = fixed ion mass (None: preset m = 19)."""
    overrides = dict(GATE_OVERRIDES)
    if mass_amu is not None:
        overrides.update(m_eff_amu=float(mass_amu), mass_initial_amu=float(mass_amu))
    print(f"  building {CASE} case setup (one neutral run) at the 3.5 A gate, "
          f"m = {mass_amu or 'preset'} amu ...")
    setup = build_case_setup(
        CASE, num_molecules=N, seed=SEED, ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS, reference_filename=REFERENCE_FILENAME,
        cfg_overrides=overrides,
    )
    if mass_amu is not None and setup.cfg_base.m_eff_amu != float(mass_amu):
        raise AssertionError(f"setup mass {setup.cfg_base.m_eff_amu} != {mass_amu}")
    return setup


def _check(label: str, got: float, committed: float) -> None:
    ok = got == committed
    print(f"  {label}: rmse {got!r} vs committed {committed!r} -> "
          f"{'BIT-EXACT' if ok else 'MISMATCH'}")
    if not ok:
        raise AssertionError(f"oracle failed: {label}")


def run_oracle() -> None:
    setup = _setup()
    lin = load_json(G49_LINEAR_SUMMARY)
    r = evaluate_form_objective(LINEAR_QUADRATIC, {"a": lin["a"], "c": 0.0},
                                lin["e_bind_eV"], setup)
    _check("(i) G4.9 linear via lq", r.rmse_Aps, lin["per_case"][CASE]["rmse_Aps"])
    quad = load_json(G411_QUADRATIC_SUMMARY)
    r = evaluate_form_objective(LINEAR_QUADRATIC, quad["coefficients"],
                                quad["e_bind_eV"], setup)
    _check("(ii) G4.11 quadratic", r.rmse_Aps, quad["per_case"][CASE]["rmse_Aps"])
    print("  ORACLE PASS: (i) and (ii) bit-exact.")


_G412_STAMP = ("THESIS_FINAL_GOALS.md G4.12 step 1 (18 A ridge, s_rho = 3.5 A)",
               "THESIS_FINAL_GOALS.md G4.12, frozen beb5440 before any run")
_G414_STAMP = ("THESIS_FINAL_GOALS.md G4.14 Step 0a (fine 18 A ridge points, s_rho = 3.5 A)",
               "THESIS_FINAL_GOALS.md G4.14, frozen 7726e43 before any run")


def _point_record(fit, c: float, goal_stamp=_G412_STAMP) -> dict:
    b = fit.best
    a = b.coefficients["a"]
    r = b.per_case[CASE]
    forces = {f"F_{v:g}".replace(".", "p"): float(lq_force(a, c, v))
              for v in FORCE_PROBE_SPEEDS_APS}
    return {
        **stamp(),
        "goal": goal_stamp[0],
        "preregistration": goal_stamp[1],
        "variant": fit.variant, "form": fit.form,
        "c_amu_per_A": c, "a_amu_per_ps": a, "e_bind_eV": b.e_bind_eV[CASE],
        "rmse_18A_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction,
        "objective_Aps": b.objective,
        **forces,
        "Phi": float(lq_force(a, c, 9.7)) / H405_PLATEAU,
        "bands_pass": bool(r.rmse_Aps <= RMSE_18A_MAX_APS and r.escape_fraction == 1.0),
        "n_evaluations": fit.n_evaluations, "converged": fit.converged,
        "start_minima": [{"a": s.coefficients["a"], "e_bind_eV": s.e_bind_eV[CASE],
                          "objective_Aps": s.objective} for s in fit.starts],
        "anchors": dict(fit.anchors),
        "gate": {"drag_spatial_gate": fit.drag_spatial_gate,
                 "drag_gate_steepness_angstrom": fit.gate_steepness_angstrom,
                 "potential_steepness_angstrom": fit.potential_steepness_angstrom},
        "n_molecules": fit.n_molecules, "seed": fit.seed, "meff_amu": fit.meff_amu,
    }


def _fit_point(setup, c: float, goal_stamp=_G412_STAMP) -> dict:
    t0 = time.perf_counter()
    fit = fit_form_trajectory_matching(
        setup, variant=LQ_FIXED_C,
        anchors={"a0": A0_FIXED_C, "c_fixed": float(c)},
        maxfev_per_start=MAXFEV_PER_START,
    )
    if fit.gate_steepness_angstrom != GATE_WIDTH_A:
        raise AssertionError(f"fit ran at gate {fit.gate_steepness_angstrom} A")
    rec = _point_record(fit, float(c), goal_stamp)
    print(f"    c={c:.3f}: {time.perf_counter() - t0:.0f} s, a={rec['a_amu_per_ps']:.3f}, "
          f"E={rec['e_bind_eV']:.4f}, rmse={rec['rmse_18A_Aps']:.4f}, "
          f"Phi={rec['Phi']:.3f}", flush=True)
    return rec


def run_ridge(c_values) -> None:
    setup = _setup()
    for c in c_values:
        dump_json(_fit_point(setup, float(c)), RIDGE_DIR / f"{point_label(float(c))}.json")


def run_oracle414() -> None:
    """G4.14 Step-0a oracle: the c = 2.0 refit equals the committed point."""
    committed = load_json(RIDGE_DIR / f"{point_label(G414_ORACLE_C)}.json")
    rec = _fit_point(_setup(), G414_ORACLE_C)
    for k in ("a_amu_per_ps", "e_bind_eV", "rmse_18A_Aps"):
        _check(f"G4.14 oracle c={G414_ORACLE_C} {k}", rec[k], committed[k])
    print("  ORACLE PASS: the c = 2.0 refit reproduces c2p000.json bit-exactly.")


def run_ridge414(c_values) -> None:
    paths = {c: RIDGE_DIR / f"{point_label(float(c))}.json" for c in c_values}
    existing = [p.name for p in paths.values() if p.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite committed ridge points {existing}")
    setup = _setup()
    for c, path in paths.items():
        dump_json(_fit_point(setup, float(c), _G414_STAMP), path)


RIDGE_COLUMNS = ("c_amu_per_A", "a_amu_per_ps", "e_bind_eV", "rmse_18A_Aps",
                 "escape_fraction", "F_3", "F_5p5", "F_9p7", "Phi", "bands_pass",
                 "converged", "n_evaluations")


def load_ridge(c_values=C_GRID) -> list[dict]:
    recs = [load_json(RIDGE_DIR / f"{point_label(c)}.json") for c in c_values]
    return sorted(recs, key=lambda r: r["c_amu_per_A"])


def run_assemble() -> None:
    missing = [c for c in C_GRID if not (RIDGE_DIR / f"{point_label(c)}.json").exists()]
    if missing:
        raise FileNotFoundError(f"ridge points not yet fitted: {missing}")
    recs = load_ridge()
    with open(RIDGE_DIR / "ridge.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(RIDGE_COLUMNS)
        for r in recs:
            w.writerow([r[k] for k in RIDGE_COLUMNS])
    lin = load_json(G49_LINEAR_SUMMARY)
    seam = recs[0]
    da = abs(seam["a_amu_per_ps"] - lin["a"])
    de = abs(seam["e_bind_eV"] - lin["e_bind_eV"])
    ok = seam["c_amu_per_A"] == 0.0 and da <= SEAM_DA_MAX and de <= SEAM_DE_MAX_EV
    print(f"  oracle (iii) c=0 seam: |da|={da:.4f} (<= {SEAM_DA_MAX}), "
          f"|dE|={de:.5f} (<= {SEAM_DE_MAX_EV}) -> {'PASS' if ok else 'FAIL'}")
    for r in recs:
        print("   " + ", ".join(f"{k}={r[k]:.4g}" if isinstance(r[k], float) else f"{k}={r[k]}"
                                for k in RIDGE_COLUMNS))
    if not ok:
        raise AssertionError("oracle (iii) failed: the c = 0 ridge point is not the G4.9 linear point")


# =============================================================================
# Step 2: the 9 A exit bound
# =============================================================================
def hedft_atom_summary(hedft, droplet_radius_A: float, window_9a) -> dict:
    """Per-atom HeDFT numbers: transverse fraction, exit speed, exit-window rho_hat."""
    out = {}
    for k in (1, 2):
        pos = np.column_stack([getattr(hedft, f"{ax}{k}_A") for ax in "xyz"])
        vel = np.column_stack([getattr(hedft, f"v{k}_{ax}_Aps") for ax in "xyz"])
        t = hedft.time_ps
        in_win = (t >= window_9a[0]) & (t <= window_9a[1])
        frac = transverse_fraction(pos[in_win], vel[in_win])
        depth = np.linalg.norm(pos, axis=1) - droplet_radius_A
        rho = rho_he_ratio(depth, steepness=GATE_WIDTH_A)
        ex = (t >= EXIT_WINDOW_PS[0]) & (t <= EXIT_WINDOW_PS[1])
        out[k] = {
            "mean_transverse_fraction": float(np.mean(frac)),
            "v_exit_Aps": window_mean(t, getattr(hedft, f"v{k}_magnitude_Aps"), EXIT_WINDOW_PS),
            "rho_hat_max_exit": float(np.max(rho[ex])),
        }
    return out


def score_run(run_dir: Path, hedft, window_9a) -> dict:
    """Exit-bound score of one finished 9 A Tier-0 run directory."""
    from i2_helium_md.simulation.run_directory import RunDirectory
    from scripts.post_processing.tier0_drag_comparison import ensemble_mean_series

    run = RunDirectory(run_dir)
    ion, cfg = run.load_ion(), run.load_cfg()
    if abs(float(ion.time_ps[0]) - float(hedft.time_ps[0])) > 1e-6:
        raise ValueError(f"{run_dir.name}: MD/HeDFT time origins differ")
    if float(ion.time_ps[-1]) < EXIT_WINDOW_PS[1]:
        raise ValueError(f"{run_dir.name}: run ends before the exit window")
    t, _, v_i1, v_i2 = ensemble_mean_series(ion)
    v_pool = 0.5 * (v_i1 + v_i2)   # radial MD: both atoms equivalent
    radius = float(ion.droplet_radii_angstrom[0])
    atoms = hedft_atom_summary(hedft, radius, window_9a)
    for k, a in atoms.items():
        if a["rho_hat_max_exit"] > RHO_INSIDE_MAX:
            raise RuntimeError(
                f"HeDFT atom {k} still inside the surface in the exit window "
                f"(rho_hat {a['rho_hat_max_exit']:.3f} > {RHO_INSIDE_MAX}): STOP, "
                f"return to the user before scoring (G4.12 step 2)")
    less_transverse = min(atoms, key=lambda k: atoms[k]["mean_transverse_fraction"])
    v_exit_md = window_mean(t, v_pool, EXIT_WINDOW_PS)
    v_ref_mean = 0.5 * (hedft.v1_magnitude_Aps + hedft.v2_magnitude_Aps)
    coeffs = cfg.drag_coefficients
    return {
        "run": run_dir.name, "bundle": cfg.drag_bundle, "form": coeffs.form,
        **{f"coef_{k}": float(v) for k, v in coeffs.coefficients.items()},
        "e_bind_eV": float(cfg.binding_energy_I_ion_eV),
        "m_eff_amu": float(cfg.m_eff_amu),
        "droplet_radius_A": radius,
        "less_transverse_atom": less_transverse,
        "transverse_frac_atom1": atoms[1]["mean_transverse_fraction"],
        "transverse_frac_atom2": atoms[2]["mean_transverse_fraction"],
        "rho_hat_max_exit_atom1": atoms[1]["rho_hat_max_exit"],
        "rho_hat_max_exit_atom2": atoms[2]["rho_hat_max_exit"],
        "v_exit_md_Aps": v_exit_md,
        "v_exit_atom1_Aps": atoms[1]["v_exit_Aps"],
        "v_exit_atom2_Aps": atoms[2]["v_exit_Aps"],
        **exit_bound(v_exit_md, atoms[less_transverse]["v_exit_Aps"]),
        "rmse_vs_v1_Aps": rmse_on_grid(t, v_pool, hedft.time_ps, hedft.v1_magnitude_Aps, window_9a),
        "rmse_vs_v2_Aps": rmse_on_grid(t, v_pool, hedft.time_ps, hedft.v2_magnitude_Aps, window_9a),
        "rmse_vs_vmean_Aps": rmse_on_grid(t, v_pool, hedft.time_ps, v_ref_mean, window_9a),
    }


def _window_9a() -> tuple[float, float]:
    from scripts.tier0_common import window_source_dir
    p = load_json(window_source_dir("9A") / "fit_parameters.json")
    return float(p["t_start"]), float(p["t_end"])


def _write_rows(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows)} rows -> {path.relative_to(PROJECT_ROOT).as_posix()}")


def _print_rows(rows: list[dict]) -> None:
    for r in rows:
        print(f"   {r['run']}: v_exit MD {r['v_exit_md_Aps']:.3f} vs atom1 "
              f"{r['v_exit_atom1_Aps']:.3f} / atom2 {r['v_exit_atom2_Aps']:.3f} -> "
              f"margin {r['margin_Aps']:+.3f} {'PASS' if r['pass'] else 'FAIL'}"
              f"{' (at bound)' if r['at_bound'] else ''}; rmse v1 {r['rmse_vs_v1_Aps']:.3f}, "
              f"v2 {r['rmse_vs_v2_Aps']:.3f}, mean {r['rmse_vs_vmean_Aps']:.3f}; "
              f"less-transverse atom {r['less_transverse_atom']}")


def run_anchor() -> None:
    from scripts.tier0_common import run_dir_name
    hedft, window = load_hedft_trajectory(HEDFT_9A_PATH), _window_9a()
    rows = [score_run(RUNS_ROOT / run_dir_name("9A", v, BOUND_N), hedft, window)
            for v in ANCHOR_VARIANTS]
    _print_rows(rows)
    _write_rows(rows, BOUND_DIR / "bound_anchor.csv")


def _run_point(c: float, a: float, e_bind: float, *, mass_amu: float | None = None,
               tag: str | None = None) -> Path:
    from scripts.tier0_common import build_drag_cfg, run_dir_name
    from i2_helium_md.simulation.ion import run_ion_propagation
    from i2_helium_md.simulation.neutral import run_neutral_propagation
    from i2_helium_md.simulation.run_directory import RunDirectory

    cfg = build_drag_cfg("9A", BOUND_VARIANT, num_molecules=BOUND_N,
                         ion_time_ps=BOUND_ION_TIME_PS, dt_ion_ps=BOUND_DT_ION_PS,
                         seed=SEED, coeff_overrides={"a": a, "c": c},
                         e_bind_override=e_bind, m_eff_amu=mass_amu)
    cfg.validate()
    run_dir = RUNS_ROOT / run_dir_name("9A", BOUND_VARIANT, BOUND_N,
                                       tag or f"g412{point_label(c)}")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    return run_dir


def run_bound(c_values) -> None:
    hedft, window = load_hedft_trajectory(HEDFT_9A_PATH), _window_9a()
    recs = {r["c_amu_per_A"]: r for r in load_ridge(c_values)}
    rows = []
    for c in c_values:
        r = recs[float(c)]
        t0 = time.perf_counter()
        run_dir = _run_point(float(c), r["a_amu_per_ps"], r["e_bind_eV"])
        row = {"c_amu_per_A": float(c), **score_run(run_dir, hedft, window)}
        print(f"    c={c:.3f}: {time.perf_counter() - t0:.0f} s", flush=True)
        rows.append(row)
        dump_json({**stamp(), **row}, BOUND_DIR / f"{point_label(float(c))}.json")
    _print_rows(rows)


def run_bound_assemble() -> None:
    rows = [load_json(BOUND_DIR / f"{point_label(c)}.json") for c in C_GRID]
    for r in rows:
        for k in ("date", "branch", "commit"):
            r.pop(k, None)
    _print_rows(rows)
    _write_rows(rows, BOUND_DIR / "bound.csv")


G414_RIDGE_CSV = RIDGE_DIR / "ridge_g414.csv"
G414_COLUMNS = RIDGE_COLUMNS + ("margin_Aps", "bound_pass", "at_bound", "passes_all")


def run_assemble414() -> None:
    """G4.14 Step 0a: the fine points with their bands and 9 A bound."""
    rows = []
    for r in load_ridge(G414_C_FINE):
        b = load_json(BOUND_DIR / f"{point_label(r['c_amu_per_A'])}.json")
        if b["coef_c"] != r["c_amu_per_A"] or b["coef_a"] != r["a_amu_per_ps"]:
            raise AssertionError(f"bound/ridge disagree at c={r['c_amu_per_A']}")
        rows.append({**{k: r[k] for k in RIDGE_COLUMNS},
                     "margin_Aps": b["margin_Aps"], "bound_pass": b["pass"],
                     "at_bound": b["at_bound"],
                     "passes_all": bool(r["bands_pass"] and b["pass"])})
    with open(G414_RIDGE_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=G414_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"   c={r['c_amu_per_A']}: a={r['a_amu_per_ps']:.4f} E={r['e_bind_eV']:.5f} "
              f"rmse={r['rmse_18A_Aps']:.4f} bands={r['bands_pass']} "
              f"margin9A={r['margin_Aps']:+.3f} -> passes_all={r['passes_all']}")
    print(f"  wrote {len(rows)} rows -> {G414_RIDGE_CSV.relative_to(PROJECT_ROOT).as_posix()}")


# =============================================================================
# G11.4 step 2: the ridge and the 9 A bound at m = I+ + 21 He
# =============================================================================
_G11_STAMP = ("THESIS_FINAL_GOALS.md G11.4 step 2 (18 A ridge at m = I+ + 21 He, s_rho = 3.5 A)",
              "THESIS_FINAL_GOALS.md G11.2 / G11.4, agreed 2026-10-09 before any run")


def _committed_point(c: float) -> dict:
    """The committed m = 19 ridge point at ``c`` (G4.12 / G4.14)."""
    return load_json(RIDGE_DIR / f"{point_label(float(c))}.json")


def run_oracle_m21() -> None:
    """The mass-parameter path at m = 19 reproduces the committed c = 1.6 point."""
    ref = _committed_point(1.6)
    setup = _setup(ref["meff_amu"])
    r = evaluate_form_objective(LINEAR_QUADRATIC, {"a": ref["a_amu_per_ps"], "c": 1.6},
                                ref["e_bind_eV"], setup)
    _check("G11 oracle c=1.6 at m=19 via the mass path", r.rmse_Aps, ref["rmse_18A_Aps"])
    print("  ORACLE PASS: the mass-parameter path is the identity at m = 19.")


def run_check_m21(labels=tuple(G11_CHECK_POINTS)) -> None:
    """Committed production points, coefficients unchanged, re-scored at m = 21."""
    setup = _setup(G11_M_EFF_AMU)
    for label in labels:
        ref = _committed_point(G11_CHECK_POINTS[label])
        a, c, e = ref["a_amu_per_ps"], ref["c_amu_per_A"], ref["e_bind_eV"]
        r = evaluate_form_objective(LINEAR_QUADRATIC, {"a": a, "c": c}, e, setup)
        rec = {
            **stamp(), "goal": _G11_STAMP[0], "preregistration": _G11_STAMP[1],
            "label": label, "source": f"g412_ridge/{point_label(c)}.json (m = 19 fit)",
            "c_amu_per_A": c, "a_amu_per_ps": a, "e_bind_eV": e,
            "rmse_18A_Aps": r.rmse_Aps, "escape_fraction": r.escape_fraction,
            "rmse_18A_Aps_m19": ref["rmse_18A_Aps"],
            "bands_pass": bool(r.rmse_Aps <= RMSE_18A_MAX_APS and r.escape_fraction == 1.0),
            "meff_amu": float(setup.cfg_base.m_eff_amu),
        }
        print(f"    {label} (c={c}, a={a:.4f}, E={e:.5f}) at m=21: rmse {r.rmse_Aps:.4f} "
              f"(m=19 {ref['rmse_18A_Aps']:.4f}), escape {r.escape_fraction} -> "
              f"bands {'PASS' if rec['bands_pass'] else 'FAIL'}", flush=True)
        dump_json(rec, G11_RIDGE_DIR / f"check_{label}.json")


def run_ridge_m21(c_values) -> None:
    paths = {c: G11_RIDGE_DIR / f"{point_label(float(c))}.json" for c in c_values}
    existing = [p.name for p in paths.values() if p.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite m = 21 ridge points {existing}")
    setup = _setup(G11_M_EFF_AMU)
    for c, path in paths.items():
        dump_json(_fit_point(setup, float(c), _G11_STAMP), path)


def _g11_point(label: str) -> dict:
    """A G11 point record by label: a float string = refit ridge point, else a check."""
    if label in G11_CHECK_POINTS:
        return load_json(G11_RIDGE_DIR / f"check_{label}.json")
    return load_json(G11_RIDGE_DIR / f"{point_label(float(label))}.json")


def run_bound_m21(labels) -> None:
    hedft, window = load_hedft_trajectory(HEDFT_9A_PATH), _window_9a()
    for label in labels:
        r = _g11_point(label)
        name = label if label in G11_CHECK_POINTS else point_label(float(label))
        t0 = time.perf_counter()
        run_dir = _run_point(r["c_amu_per_A"], r["a_amu_per_ps"], r["e_bind_eV"],
                             mass_amu=G11_M_EFF_AMU, tag=f"g11m21{name}")
        row = {"label": name, "c_amu_per_A": r["c_amu_per_A"],
               **score_run(run_dir, hedft, window)}
        if row["m_eff_amu"] != G11_M_EFF_AMU:
            raise AssertionError(f"{name}: 9 A run at m = {row['m_eff_amu']}")
        print(f"    {name}: {time.perf_counter() - t0:.0f} s", flush=True)
        _print_rows([row])
        dump_json({**stamp(), **row}, G11_BOUND_DIR / f"{name}.json")


G11_COLUMNS = ("label", "c_amu_per_A", "a_amu_per_ps", "e_bind_eV", "rmse_18A_Aps",
               "escape_fraction", "bands_pass", "a_m19", "e_bind_eV_m19", "rmse_18A_Aps_m19",
               "margin_Aps", "margin_Aps_m19", "bound_pass", "at_bound", "passes_all")


def run_assemble_m21() -> None:
    """Refit ridge + unchanged checks, each with its m = 19 counterpart and 9 A bound."""
    rows = []
    labels = [f"{c:g}" for c in C_GRID] + list(G11_CHECK_POINTS)
    for label in labels:
        r = _g11_point(label)
        name = label if label in G11_CHECK_POINTS else point_label(float(label))
        ref = _committed_point(r["c_amu_per_A"])
        b = load_json(G11_BOUND_DIR / f"{name}.json")
        b19 = load_json(BOUND_DIR / f"{point_label(r['c_amu_per_A'])}.json")
        if b["coef_a"] != r["a_amu_per_ps"] or b["e_bind_eV"] != r["e_bind_eV"]:
            raise AssertionError(f"bound/ridge disagree at {name}")
        rows.append({
            "label": name, "c_amu_per_A": r["c_amu_per_A"],
            "a_amu_per_ps": r["a_amu_per_ps"], "e_bind_eV": r["e_bind_eV"],
            "rmse_18A_Aps": r["rmse_18A_Aps"], "escape_fraction": r["escape_fraction"],
            "bands_pass": r["bands_pass"],
            "a_m19": ref["a_amu_per_ps"], "e_bind_eV_m19": ref["e_bind_eV"],
            "rmse_18A_Aps_m19": ref["rmse_18A_Aps"],
            "margin_Aps": b["margin_Aps"], "margin_Aps_m19": b19["margin_Aps"],
            "bound_pass": b["pass"], "at_bound": b["at_bound"],
            "passes_all": bool(r["bands_pass"] and b["pass"]),
        })
    with open(G11_RIDGE_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=G11_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"   {r['label']}: c={r['c_amu_per_A']} a={r['a_amu_per_ps']:.4f} "
              f"(m19 {r['a_m19']:.4f}) E={r['e_bind_eV']:.5f} (m19 {r['e_bind_eV_m19']:.5f}) "
              f"rmse={r['rmse_18A_Aps']:.4f} bands={r['bands_pass']} "
              f"margin9A={r['margin_Aps']:+.3f} (m19 {r['margin_Aps_m19']:+.3f}) "
              f"-> passes_all={r['passes_all']}")
    print(f"  wrote {len(rows)} rows -> {G11_RIDGE_CSV.relative_to(PROJECT_ROOT).as_posix()}")


def main(mode: str = MODE, c_values=C_GRID) -> int:
    print(f"===== G4.12 ridge / exit bound, MODE={mode} =====")
    if mode == "oracle":
        run_oracle()
    elif mode == "ridge":
        run_ridge(c_values)
    elif mode == "assemble":
        run_assemble()
    elif mode == "anchor":
        run_anchor()
    elif mode == "bound":
        run_bound(c_values)
    elif mode == "bound_assemble":
        run_bound_assemble()
    elif mode == "oracle414":
        run_oracle414()
    elif mode == "ridge414":
        run_ridge414(c_values)
    elif mode == "bound414":
        run_bound(c_values)
    elif mode == "assemble414":
        run_assemble414()
    elif mode == "oracle_m21":
        run_oracle_m21()
    elif mode == "check_m21":
        run_check_m21()
    elif mode == "ridge_m21":
        run_ridge_m21(c_values)
    elif mode == "bound_m21":
        run_bound_m21(c_values)
    elif mode == "assemble_m21":
        run_assemble_m21()
    else:
        raise ValueError(f"unknown MODE {mode!r}")
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = args[0] if args else MODE
    if mode == "bound_m21":   # labels: floats = refit ridge points, names = checks
        labels = tuple(args[1:]) or (*G11_CHECK_POINTS, *(f"{c:g}" for c in C_GRID))
        raise SystemExit(main(mode, labels))
    default_c = G414_C_FINE if mode.endswith("414") else C_GRID
    raise SystemExit(main(mode, tuple(float(x) for x in args[1:]) or default_c))
