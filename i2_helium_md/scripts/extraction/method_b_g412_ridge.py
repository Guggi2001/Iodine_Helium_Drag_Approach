"""The 18 A linear + quadratic ridge and its 9 A exit bound at m = I+ + 21 He.

THESIS_FINAL_GOALS.md G4.12 / G11.4. 18 A fixes one number, the drag at
v ~ 3 A/ps, so ``F = a v + c v^2`` along a one-parameter 18 A ridge fits it. This
driver maps that ridge at the s_rho = 3.5 A gate and the 18 A shell mass
(production ``m(n)`` convention, ``G11_M_EFF_AMU``) and bounds each point by the
radial 9 A exit speed.

Units: a [amu/ps], c [amu/A]; each term amu*A/ps^2; gamma = a + c|v| [amu/ps].

Modes (USER SETTINGS ``MODE``; argv override ``python ... ridge_m21 0 0.4``):

* ``"ridge_m21"`` -- fixed-c 18 A fits (``lq_fixed_c``: a >= 0, E_bind free) for
  ``C_GRID`` (or argv); refuses to overwrite an existing point JSON.
* ``"bound_m21"`` -- the 9 A exit bound for argv labels (floats = refit ridge
  points, names = check points; default all).

``ridge_m21.csv`` and the ``check_*.json`` records are frozen: B2 / A5 are the
m = 19 fits re-scored unchanged at m = 21, and those fits are kept only in git
tag ``pre-drag-cleanup``.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_g412_ridge.py
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
import sys
import time

import numpy as np

# =============================================================================
# USER SETTINGS
# =============================================================================
MODE = "ridge_m21"   # "ridge_m21" | "bound_m21"

# Pre-registered ridge grid [amu/A] (G4.12 step 1; FIXED).
C_GRID = (0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.224)

# --- Pre-registered bands / constants (FIXED; do not re-tune) -----------------
RMSE_18A_MAX_APS = 0.19            # G4.9 band; escape fraction must be 1.0
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
BOUND_VARIANT = "LQ"   # lq bundle; {a, c, E} overridden
GATE_OVERRIDES = {"drag_spatial_gate": "erf_independent", "drag_gate_steepness": 3.5}

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    LQ_FIXED_C,
    _git_branch,
    build_case_setup,
    evaluate_form_objective,
    fit_form_trajectory_matching,
)
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU  # noqa: E402
from i2_helium_md.physics.drag import LINEAR_QUADRATIC  # noqa: E402
from i2_helium_md.physics.helium_density import rho_he_ratio  # noqa: E402
from i2_helium_md.postprocess import load_hedft_trajectory  # noqa: E402
from i2_helium_md.presets import CALIBRATION_RIDGE_ROOT, tddft_reference_path  # noqa: E402
from i2_helium_md.simulation.battery import run_neutral_and_ion  # noqa: E402

CASE = "18A"
N = 50                   # ensemble size per forward integration (Tier-0)
SEED = 20260604          # neutral-stage seed (ion drag path is RNG-free)
ION_TIME_PS = 20.0       # covers the [t*, ~14 ps] windows with escape margin
DT_ION_PS = 0.01
HEDFT_9A_PATH = tddft_reference_path("9A")
RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
GATE_WIDTH_A = float(GATE_OVERRIDES["drag_gate_steepness"])

# G11.2: extraction mass = I+ + 21 He under the production m(n) convention
# (constants.py), so the fit mass equals production's n0 = 21 start exactly.
G11_N_HE = 21
G11_M_EFF_AMU = MASS_I_ION_AMU + G11_N_HE * MASS_HE_AMU   # 210.9546 amu
G11_RIDGE_DIR = CALIBRATION_RIDGE_ROOT / "ridge"
G11_BOUND_DIR = CALIBRATION_RIDGE_ROOT / "exit_bound_9A"
# Committed production points re-scored unchanged at m = 21 (label -> c;
# THESIS_FINAL_GOALS G5 Outcome: B2 = c 1.6, A5 = c 2.0).
G11_CHECK_POINTS = {"B2": 1.6, "A5": 2.0}


# =============================================================================
# Pure helpers (unit-tested in tests/test_g412_ridge.py)
# =============================================================================
def load_json(path: Path) -> dict:
    """The JSON document at ``path``."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def dump_json(payload: dict, path: Path) -> None:
    """Write ``payload`` to ``path`` as indented JSON, creating parent folders."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    print(f"  wrote {path}")


def stamp() -> dict:
    """Date and branch stamp for a written record."""
    return {"date": datetime.date.today().isoformat(), "branch": _git_branch(PROJECT_ROOT)}


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


def _point_record(fit, c: float, goal_stamp) -> dict:
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


def _fit_point(setup, c: float, goal_stamp) -> dict:
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
    from scripts.tier0_common import TIER0_WINDOWS_PS
    return TIER0_WINDOWS_PS["9A"]


def _print_rows(rows: list[dict]) -> None:
    for r in rows:
        print(f"   {r['run']}: v_exit MD {r['v_exit_md_Aps']:.3f} vs atom1 "
              f"{r['v_exit_atom1_Aps']:.3f} / atom2 {r['v_exit_atom2_Aps']:.3f} -> "
              f"margin {r['margin_Aps']:+.3f} {'PASS' if r['pass'] else 'FAIL'}"
              f"{' (at bound)' if r['at_bound'] else ''}; rmse v1 {r['rmse_vs_v1_Aps']:.3f}, "
              f"v2 {r['rmse_vs_v2_Aps']:.3f}, mean {r['rmse_vs_vmean_Aps']:.3f}; "
              f"less-transverse atom {r['less_transverse_atom']}")


def point_cfg(c: float, a: float, e_bind: float, *, mass_amu: float) -> SimConfig:
    """Radial 9 A Tier-0 config of one lq law (c [amu/A], a [amu/ps], E_bind [eV])."""
    from scripts.tier0_common import build_drag_cfg

    cfg = build_drag_cfg("9A", BOUND_VARIANT, num_molecules=BOUND_N,
                         ion_time_ps=BOUND_ION_TIME_PS, dt_ion_ps=BOUND_DT_ION_PS,
                         seed=SEED, coeff_overrides={"a": a, "c": c},
                         e_bind_override=e_bind, m_eff_amu=mass_amu)
    cfg.validate()
    return cfg


# =============================================================================
# G11.4 step 2: the ridge and the 9 A bound at m = I+ + 21 He
# =============================================================================
_G11_STAMP = ("THESIS_FINAL_GOALS.md G11.4 step 2 (18 A ridge at m = I+ + 21 He, s_rho = 3.5 A)",
              "THESIS_FINAL_GOALS.md G11.2 / G11.4, agreed 2026-10-09 before any run")


def check_point(label: str) -> dict:
    """The committed check record of a production law (``G11_CHECK_POINTS`` label)."""
    return load_json(G11_RIDGE_DIR / f"check_{label}.json")


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
        return check_point(label)
    return load_json(G11_RIDGE_DIR / f"{point_label(float(label))}.json")


def run_bound_m21(labels) -> None:
    from scripts.tier0_common import run_dir_name

    hedft, window = load_hedft_trajectory(HEDFT_9A_PATH), _window_9a()
    for label in labels:
        r = _g11_point(label)
        name = label if label in G11_CHECK_POINTS else point_label(float(label))
        t0 = time.perf_counter()
        run_dir = RUNS_ROOT / run_dir_name("9A", BOUND_VARIANT, BOUND_N, f"g11m21{name}")
        run_neutral_and_ion(point_cfg(r["c_amu_per_A"], r["a_amu_per_ps"], r["e_bind_eV"],
                                      mass_amu=G11_M_EFF_AMU), run_dir)
        row = {"label": name, "c_amu_per_A": r["c_amu_per_A"],
               **score_run(run_dir, hedft, window)}
        if row["m_eff_amu"] != G11_M_EFF_AMU:
            raise AssertionError(f"{name}: 9 A run at m = {row['m_eff_amu']}")
        print(f"    {name}: {time.perf_counter() - t0:.0f} s", flush=True)
        _print_rows([row])
        dump_json({**stamp(), **row}, G11_BOUND_DIR / f"{name}.json")


def main(mode: str = MODE, c_values=C_GRID) -> int:
    print(f"===== 18 A ridge / 9 A exit bound at m = 21, MODE={mode} =====")
    if mode == "ridge_m21":
        run_ridge_m21(c_values)
    elif mode == "bound_m21":
        run_bound_m21(c_values)
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
    raise SystemExit(main(mode, tuple(float(x) for x in args[1:]) or C_GRID))
