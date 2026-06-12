"""Method-B trajectory-matching drag extraction (joint {a, b, E_bind} fit).

Fits the ``linear_cubic`` drag coefficients **jointly with the effective
droplet binding energy** by minimizing the forward-integrated, from-onset,
in-window trajectory RMSE of the ensemble-mean clean-atom speed ``|v2|``
against the same-smoothed long reference
(``data/reference/drag/<case>/velocity_smoothed/cleaned_data_long.csv``),
plus an escape penalty that steers the optimizer away from trapping binding
depths (``TIER0_FINDINGS.md`` -> "Correct drag traps the ions").

Method spec: ``METHOD_B_trajectory_matching_extraction.md``. Key properties:

* **From-onset objective** -- each evaluation is a full ion-stage forward
  integration from the explosion onset (t=0 on the reference clock), scored
  in-window only via :func:`compare_speed_to_reference`'s ``window=``. The
  pre-t* transient is part of what the coefficients must absorb (production
  usage), per the locked design decision.
* **Deterministic objective** -- the drag ion path is RNG-free (BAOAB at
  ``T_eff=0``, no collisions, no mass attachment), so the neutral checkpoint
  is computed **once** per (case, seed) and reused for every evaluation;
  the objective is then bit-deterministic in ``(a, b, E_bind)``.
* **Exact §6.5.1 pairing per evaluation** -- each trial stamps the candidate
  ``E_bind`` into the coefficient bundle *and* sets
  ``cfg.binding_energy_I_ion_eV`` to the same value, so the config-load guard
  passes silently (no escape hatch needed inside the fit loop).
* **VMI stays held-out** -- it is never part of this objective; it is the
  downstream arbiter (METHOD_B §4).

Units: ``a`` [amu/ps], ``b`` [amu*ps/A^2], ``E_bind`` [eV], speeds [A/ps],
times [ps]; the objective and RMSE are in A/ps.

The module also carries the **§9 shared-form joint refit** machinery
(:func:`evaluate_joint_objective`, :func:`fit_shared_trajectory_matching`
with the ``shared_3param`` / ``shared_pure_cubic`` / ``diagnostic_4param``
variants, :func:`joint_sensitivity_halfwidths`,
:func:`write_shared_fit_parameters`): one size-independent drag law fit
jointly across both cases, equal-weight per-case objective mean, anchors
from the 18 A Method-A bundle only. Driver:
``scripts/extraction/method_b_shared_refit.py``.
"""

from __future__ import annotations

import datetime as _datetime
import json
import subprocess
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import Bounds, minimize

from ..config import SimConfig
from ..physics.drag import DragCoefficients, LINEAR_CUBIC
from ..presets import (
    REFERENCE_DRAG_ROOT,
    single_pulse_N2000_18Angst_drag,
    single_pulse_N2000_drag,
)
from ..postprocess.compare_trajectories import compare_speed_to_reference
from ..postprocess.hedft_loader import load_smoothed_speed_reference
from ..simulation.checkpoint import IonCheckpoint, NeutralCheckpoint
from ..simulation.ion import run_ion_propagation
from ..simulation.neutral import run_neutral_propagation

# Static solvation energy of I+ in the droplet [eV] -- the physical UPPER
# BOUND for the effective binding (TIER0_FINDINGS §"Correct drag traps the
# ions": the dynamical escape barrier is below the static value, never above).
E_BIND_STATIC_EV = 0.308

# Objective sentinel [A/ps] for a non-finite trajectory: large but finite so
# Nelder-Mead can recover instead of crashing.
_NONFINITE_SENTINEL_APS = 1.0e3

_CASE_PRESETS = {
    "9A": single_pulse_N2000_drag,
    "18A": single_pulse_N2000_18Angst_drag,
}


@dataclass(frozen=True)
class ObjectiveResult:
    """One objective evaluation at ``(a, b, e_bind_eV)``.

    ``objective = rmse_Aps + penalty_Aps`` where
    ``penalty_Aps = penalty_weight * (1 - escape_fraction)`` [A/ps]. An
    accepted fit must have ``penalty_Aps == 0`` (full escape).
    """

    a: float                  # amu/ps
    b: float                  # amu*ps/A^2
    e_bind_eV: float          # eV
    objective: float          # A/ps (rmse + penalty; sentinel if non-finite)
    rmse_Aps: float           # A/ps, in-window |v2| RMSE vs smoothed reference
    escape_fraction: float    # in [0, 1], over all 2N ions at the final step
    penalty_Aps: float        # A/ps
    num_overlap_points: int


@dataclass
class CaseSetup:
    """Per-case fit context: base config, cached neutral stage, reference.

    Build via :func:`build_case_setup`. The neutral checkpoint is computed
    once and reused across all objective evaluations (drag coefficients and
    binding only enter the ion stage).
    """

    case: str
    cfg_base: SimConfig
    neutral: NeutralCheckpoint
    t_ref_ps: np.ndarray      # ps, strictly increasing
    ref_speed_Aps: np.ndarray  # A/ps, same length
    window: tuple[float, float]  # ps; [first, last] sample of the reference
    a0: float                 # amu/ps; optimizer normalization anchor
    b0: float                 # amu*ps/A^2; optimizer normalization anchor
    reference_path: Path


@dataclass
class TrajectoryMatchingFit:
    """A completed Method-B fit for one case (best of all starts)."""

    case: str
    best: ObjectiveResult
    starts: list[ObjectiveResult]      # best point of each optimizer start
    n_evaluations: int
    converged: bool
    seed: int
    n_molecules: int
    window: tuple[float, float]
    reference_path: Path
    meff_amu: float
    # Uncertainty (filled by the driver script; see the module docstring of
    # scripts/extraction/method_b_extraction.py for semantics):
    a_err: float | None = None
    b_err: float | None = None
    e_bind_err_eV: float | None = None
    uncertainty_model: str | None = None
    seed_sweep_seeds: list[int] = field(default_factory=list)
    seed_sweep_results: list[dict] = field(default_factory=list)


def build_case_setup(
    case: str,
    *,
    num_molecules: int = 50,
    seed: int = 20260604,
    ion_time_ps: float = 20.0,
    dt_ion_ps: float = 0.01,
    reference_filename: str = "cleaned_data_long.csv",
) -> CaseSetup:
    """Build the per-case fit context (runs the neutral stage ONCE).

    Parameters
    ----------
    case : str
        ``"9A"`` or ``"18A"``.
    num_molecules : int
        Ensemble size for each forward integration. N=50 is the committed
        Tier-0 evidence point (N=50 vs N=2000 in-window RMSE identical to
        0.001 A/ps, ``tests/test_tier0_drag_comparison.py``).
    seed : int
        Neutral-stage RNG seed (the ion drag path is RNG-free). 20260604 is
        the Tier-0 run-generation seed (``scripts/gen_tier0_runs.py``).
    ion_time_ps, dt_ion_ps : float
        Ion-stage duration / timestep [ps]. 20 ps covers the [t*, ~14 ps]
        window with escape margin at the production dt of 0.01 ps.
    reference_filename : str
        File under ``data/reference/drag/<case>/velocity_smoothed/``; the
        extended same-smoothed reference by default.

    Raises
    ------
    ValueError
        On an unknown case.
    FileNotFoundError
        If the smoothed reference is absent (raised by the loader).
    """
    if case not in _CASE_PRESETS:
        raise ValueError(f"case must be one of {sorted(_CASE_PRESETS)}, got {case!r}")

    cfg = _CASE_PRESETS[case](
        num_molecules=num_molecules,
        ion_simulation_time=ion_time_ps,
        dt_ion=dt_ion_ps,
        seed=seed,
    )
    cfg.validate()

    ref_path = (
        REFERENCE_DRAG_ROOT / case / "velocity_smoothed" / reference_filename
    )
    ref = load_smoothed_speed_reference(ref_path)
    window = (float(ref.time_ps[0]), float(ref.time_ps[-1]))
    if window[1] > ion_time_ps:
        raise ValueError(
            f"reference window end {window[1]:.2f} ps exceeds the ion run "
            f"duration {ion_time_ps:.2f} ps; raise ion_time_ps"
        )

    neutral = run_neutral_propagation(cfg)

    coeffs = cfg.drag_coefficients
    return CaseSetup(
        case=case,
        cfg_base=cfg,
        neutral=neutral,
        t_ref_ps=np.asarray(ref.time_ps, dtype=float),
        ref_speed_Aps=np.asarray(ref.speed_Aps, dtype=float),
        window=window,
        a0=float(coeffs.coefficients["a"]),
        b0=float(coeffs.coefficients["b"]),
        reference_path=ref_path,
    )


def _escape_fraction(ion: IonCheckpoint, gate_steepness_A: float) -> float:
    """Fraction of ions (all 2N) escaped at the final stored step.

    Escaped: final depth ``r - R_droplet`` exceeds ``+2 * gate_steepness_A``
    (well past the erf gate, drag off) with non-negative radial velocity.
    A trapped ensemble (R(t) reversing inside the well) scores ~0.
    """
    px = ion.positions_x[:, -1]
    py = ion.positions_y[:, -1]
    pz = ion.positions_z[:, -1]
    vx = ion.velocities_x[:, -1]
    vy = ion.velocities_y[:, -1]
    vz = ion.velocities_z[:, -1]
    r = np.sqrt(px * px + py * py + pz * pz)
    depth = r - np.asarray(ion.droplet_radii_angstrom, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        v_radial = (px * vx + py * vy + pz * vz) / r
    escaped = (depth > 2.0 * gate_steepness_A) & (v_radial >= 0.0)
    return float(np.mean(escaped))


def evaluate_objective(
    theta: Sequence[float],
    setup: CaseSetup,
    *,
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> ObjectiveResult:
    """Evaluate the Method-B objective at ``theta = (a, b, e_bind_eV)``.

    Builds a config that differs from ``setup.cfg_base`` ONLY in the
    ion-stage fields (candidate coefficients + binding, stamped as an exact
    §6.5.1 pair), forward-integrates the ion stage from the cached neutral
    checkpoint, and scores the ensemble-mean ``|v2|`` against the smoothed
    reference over ``setup.window``.

    Parameters
    ----------
    theta : sequence of 3 floats
        ``(a [amu/ps], b [amu*ps/A^2], e_bind [eV])``.
    penalty_weight_Aps : float
        Weight of the ``(1 - escape_fraction)`` penalty [A/ps]. The default
        2.0 sits an order of magnitude above a good in-window RMSE
        (~0.1-0.4 A/ps), making the trapping-basin boundary unambiguous even
        where the windowed RMSE is locally flat.
    run_fn : callable
        ``(cfg, neutral_ckpt) -> IonCheckpoint``-like. Dependency-injected so
        unit tests stub the MD entirely.

    Returns
    -------
    ObjectiveResult
        With the large finite sentinel as ``objective`` if the trajectory or
        score is non-finite (never raises into the optimizer for that case).
    """
    a, b, e_bind = (float(v) for v in theta)

    coeffs = DragCoefficients(
        form=LINEAR_CUBIC,
        coefficients={"a": a, "b": b},
        extraction_mass_model="constant",
        extraction_mass_amu=setup.cfg_base.m_eff_amu,
        extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=e_bind,
    )
    cfg = replace(
        setup.cfg_base,
        drag_coefficients=coeffs,
        binding_energy_I_ion_eV=e_bind,
        # Exact stamped pairing above -> §6.5.1 passes strictly; no hatch.
        allow_unvalidated_binding_pairing=False,
    )
    cfg.validate()

    ion = run_fn(cfg, setup.neutral)

    comparison = compare_speed_to_reference(
        ion,
        atom="I2",
        t_ref_ps=setup.t_ref_ps,
        ref_speed_Aps=setup.ref_speed_Aps,
        window=setup.window,
    )
    escape = _escape_fraction(ion, cfg.drag_gate_steepness)
    rmse = float(comparison.rmse)
    penalty = penalty_weight_Aps * (1.0 - escape)

    if not np.isfinite(rmse) or not np.isfinite(escape):
        objective = _NONFINITE_SENTINEL_APS
        rmse = float("nan") if not np.isfinite(rmse) else rmse
        penalty = float("nan") if not np.isfinite(escape) else penalty
    else:
        objective = rmse + penalty

    return ObjectiveResult(
        a=a,
        b=b,
        e_bind_eV=e_bind,
        objective=objective,
        rmse_Aps=rmse,
        escape_fraction=escape,
        penalty_Aps=penalty,
        num_overlap_points=int(comparison.num_overlap_points),
    )


def fit_trajectory_matching(
    setup: CaseSetup,
    *,
    e_bind_min_eV: float = 0.01,
    e_bind_max_eV: float = E_BIND_STATIC_EV,
    prescan_e_bind_eV: Sequence[float] | None = None,
    extra_starts: bool = True,
    penalty_weight_Aps: float = 2.0,
    maxfev_per_start: int = 400,
    fatol_Aps: float = 1.0e-3,
    xatol_normalized: float = 1.0e-3,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> TrajectoryMatchingFit:
    """Run the joint 3-parameter Nelder-Mead fit for one case.

    Optimizes in normalized coordinates ``x = (a/a0, b/b0, E/0.308)`` (the
    Method-A-era ``{a0, b0}`` from the case's current ``fit_parameters.json``
    anchor the conditioning) with bounds ``a in [0.1, 10]*a0`` (keeps the
    §3.3 dissipativity ``a > 0`` satisfied), ``b in [1e-3, 10]*b0``, and
    ``E in [e_bind_min_eV, e_bind_max_eV]`` (static solvation = upper bound).

    Initialization: a coarse 1-D pre-scan of ``E_bind`` at fixed
    ``(a0, b0)`` picks the best start; with ``extra_starts`` two more starts
    probe the documented drag<->binding degeneracy ridge
    (``(a0, b0, 0.9*E_max)`` and ``(1.3*a0, 0.7*b0, E_prescan)``).

    Returns
    -------
    TrajectoryMatchingFit
        ``best`` is the lowest-objective minimum across starts; ``starts``
        records each start's minimum so the degeneracy spread is visible.
        Uncertainty fields are left ``None`` (the driver script fills them).
    """
    a0, b0 = setup.a0, setup.b0
    scale = np.array([a0, b0, E_BIND_STATIC_EV], dtype=float)
    lower = np.array([0.1, 1.0e-3, e_bind_min_eV / E_BIND_STATIC_EV])
    upper = np.array([10.0, 10.0, e_bind_max_eV / E_BIND_STATIC_EV])

    n_evals = 0

    def objective_x(x: np.ndarray) -> float:
        nonlocal n_evals
        n_evals += 1
        theta = np.asarray(x, dtype=float) * scale
        return evaluate_objective(
            theta, setup, penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn
        ).objective

    # --- E_bind pre-scan at (a0, b0): pick the best start off the bound ---
    if prescan_e_bind_eV is None:
        prescan_e_bind_eV = np.linspace(0.05, 0.305, 8)
    prescan = [
        evaluate_objective(
            (a0, b0, e), setup,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )
        for e in prescan_e_bind_eV
    ]
    n_evals += len(prescan)
    e_start = min(prescan, key=lambda r: r.objective).e_bind_eV
    e_start_n = e_start / E_BIND_STATIC_EV

    x0s = [np.array([1.0, 1.0, e_start_n])]
    if extra_starts:
        x0s.append(np.array([1.0, 1.0, 0.9 * upper[2]]))
        x0s.append(np.array([1.3, 0.7, e_start_n]))
    x0s = [np.clip(x0, lower, upper) for x0 in x0s]

    start_results: list[ObjectiveResult] = []
    converged = True
    for x0 in x0s:
        res = minimize(
            objective_x,
            x0,
            method="Nelder-Mead",
            bounds=Bounds(lower, upper),
            options={
                "fatol": fatol_Aps,
                "xatol": xatol_normalized,
                "maxfev": maxfev_per_start,
            },
        )
        converged = converged and bool(res.success)
        theta_min = np.asarray(res.x, dtype=float) * scale
        start_results.append(
            evaluate_objective(
                theta_min, setup,
                penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
            )
        )
        n_evals += 1

    best = min(start_results, key=lambda r: r.objective)
    return TrajectoryMatchingFit(
        case=setup.case,
        best=best,
        starts=start_results,
        n_evaluations=n_evals,
        converged=converged,
        seed=setup.cfg_base.seed,
        n_molecules=setup.cfg_base.num_molecules,
        window=setup.window,
        reference_path=setup.reference_path,
        meff_amu=float(setup.cfg_base.m_eff_amu),
    )


def _objective_rise_halfwidths(
    objective_at: Callable[[np.ndarray], float],
    theta_best: Sequence[float],
    best_objective: float,
    *,
    rel_rise: float,
    factors: Sequence[float],
) -> tuple[float, ...]:
    """Generic multiplicative-factor half-width scan (sensitivity band core).

    For each parameter, scans multiplicative ``factors`` in both directions
    (other parameters held at the optimum) and records the parameter change
    at which ``objective_at`` first exceeds
    ``(1 + rel_rise) * best_objective``; the half-width is the mean of the
    up/down crossing distances (one side if only one crosses, the largest
    scanned change if neither does -- conservative). A parameter that is
    exactly 0 (a fixed pure-cubic ``a``) cannot be scanned multiplicatively
    and gets half-width 0.0 by convention (it was not a free parameter).
    """
    theta_best = np.asarray(theta_best, dtype=float)
    threshold = (1.0 + rel_rise) * best_objective
    halfwidths = []
    for i in range(theta_best.size):
        if theta_best[i] == 0.0:
            halfwidths.append(0.0)
            continue
        crossings = []
        for direction in (1.0, -1.0):
            crossing = None
            for f in factors:
                factor = f if direction > 0 else 1.0 / f
                theta = theta_best.copy()
                theta[i] = theta_best[i] * factor
                if objective_at(theta) > threshold:
                    crossing = abs(theta[i] - theta_best[i])
                    break
            if crossing is not None:
                crossings.append(crossing)
        if crossings:
            halfwidths.append(float(np.mean(crossings)))
        else:
            # Never crossed within the scan: flat direction -- report the
            # largest scanned change (conservative, visible in review).
            halfwidths.append(float(theta_best[i] * (factors[-1] - 1.0)))
    return tuple(halfwidths)


def sensitivity_halfwidths(
    setup: CaseSetup,
    best: ObjectiveResult,
    *,
    rel_rise: float = 0.10,
    factors: Sequence[float] = (1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0),
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> tuple[float, float, float]:
    """RMSE-sensitivity half-widths for ``(a, b, e_bind)`` around ``best``.

    Multiplicative factor scan (see :func:`_objective_rise_halfwidths` for
    the crossing semantics). Coarse by design: this is a *sensitivity band*
    on a deterministic objective, NOT a statistical confidence interval; it
    is combined with the seed-sweep spread by the driver script and labeled
    ``uncertainty_model = "seed_sweep_std_plus_rmse_sensitivity"``.
    """

    def objective_at(theta: np.ndarray) -> float:
        return evaluate_objective(
            theta, setup, penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn
        ).objective

    return _objective_rise_halfwidths(  # type: ignore[return-value]
        objective_at,
        (best.a, best.b, best.e_bind_eV),
        best.objective,
        rel_rise=rel_rise,
        factors=factors,
    )


def _reference_file_str(reference_path: Path, repo_root: Path) -> str:
    """Repo-relative POSIX path for the provenance stamp (absolute if outside,
    e.g. a test tmp dir)."""
    try:
        return str(reference_path.relative_to(repo_root).as_posix())
    except ValueError:
        return reference_path.as_posix()


def _git_branch(repo_root: Path) -> str:
    """Current git branch, or 'unknown' (provenance stamp must not crash)."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=repo_root, capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def write_fit_parameters(
    fit: TrajectoryMatchingFit,
    out_dir: Path,
    *,
    transverse_contaminated: bool,
) -> Path:
    """Write a Method-B ``fit_parameters.json`` bundle (METHOD_B §6 schema).

    The schema is a superset of the Method-A file, loadable by
    :func:`i2_helium_md.presets.load_drag_coefficients` (which requires the
    Method-B provenance keys when ``extraction_method`` is
    ``"trajectory_matching"``).

    Parameters
    ----------
    fit : TrajectoryMatchingFit
        A completed fit. ``a_err`` / ``b_err`` must be filled (the loader
        requires them); the driver script computes them.
    out_dir : Path
        Target directory (e.g. ``REFERENCE_DRAG_ROOT / case /
        "trajectory_matching"``); created if needed. Method-A artifacts are
        never touched.
    transverse_contaminated : bool
        The standing 9 A flag: True when the reference is the genuinely
        non-radial 9 A case (METHOD_B §3.5/§5); False for clean-radial 18 A.

    Raises
    ------
    ValueError
        If the fit still carries a penalty (trapped ions must never be
        committed) or the uncertainty fields are unfilled.
    """
    if not (fit.best.penalty_Aps == 0.0):
        raise ValueError(
            "refusing to write a trapped fit: escape_fraction="
            f"{fit.best.escape_fraction} (penalty {fit.best.penalty_Aps} "
            "A/ps != 0); the binding<->drag pair must permit full escape"
        )
    if fit.a_err is None or fit.b_err is None or fit.uncertainty_model is None:
        raise ValueError(
            "uncertainty fields (a_err, b_err, uncertainty_model) must be "
            "filled before writing the bundle"
        )

    repo_root = REFERENCE_DRAG_ROOT.parents[2]
    reference_file = _reference_file_str(fit.reference_path, repo_root)
    payload = {
        "extraction_method": "trajectory_matching",
        "form": LINEAR_CUBIC,
        "a": fit.best.a,
        "b": fit.best.b,
        "a_err": fit.a_err,
        "b_err": fit.b_err,
        "uncertainty_model": fit.uncertainty_model,
        "meff_amu": fit.meff_amu,
        "extraction_mass_model": "constant",
        "effective_binding_energy_I_ion_eV": fit.best.e_bind_eV,
        "effective_binding_energy_err_eV": fit.e_bind_err_eV,
        "t_start": fit.window[0],
        "t_end": fit.window[1],
        "reference_file": reference_file,
        "objective": "in_window_rmse_mean_abs_v2_plus_escape_penalty",
        "objective_rmse_Aps": fit.best.rmse_Aps,
        "escape_fraction": fit.best.escape_fraction,
        "n_molecules_fit": fit.n_molecules,
        "seed": fit.seed,
        "seed_sweep_seeds": fit.seed_sweep_seeds,
        "seed_sweep_results": fit.seed_sweep_results,
        "optimizer": {
            "method": "nelder-mead",
            "n_starts": len(fit.starts),
            "evals": fit.n_evaluations,
            "converged": fit.converged,
            "start_minima": [
                {
                    "a": r.a, "b": r.b, "e_bind_eV": r.e_bind_eV,
                    "objective_Aps": r.objective,
                }
                for r in fit.starts
            ],
        },
        "date": _datetime.date.today().isoformat(),
        "branch": _git_branch(repo_root),
        "flags": {
            "radial_projection_convention": True,
            "transverse_contaminated_non_radial_reference": (
                transverse_contaminated
            ),
            "full_window_heldout_window_axis_forfeited": True,
        },
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fit_parameters.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return out_path


# ===========================================================================
# Shared-form joint refit (METHOD_B §9)
# ===========================================================================
# One size-independent drag law {a, b} (+ binding) carried across BOTH cases:
# over-constrained (3 parameters / 2 full trajectories), hence falsifiable --
# unlike the independent per-case fits. Variants per §9.2; Stage 1 (score the
# untouched 9 A prediction from the delivered §8 18 A bundle) needs no new
# machinery -- it is one `evaluate_objective` call by the driver script.

SHARED_3PARAM = "shared_3param"          # fully shared {a, b, E_bind}, a >= 0
SHARED_PURE_CUBIC = "shared_pure_cubic"  # a == 0 fixed, fit {b, E_bind}
DIAGNOSTIC_4PARAM = "diagnostic_4param"  # shared {a, b}, per-case E_bind
JOINT_VARIANTS = (SHARED_3PARAM, SHARED_PURE_CUBIC, DIAGNOSTIC_4PARAM)


@dataclass(frozen=True)
class JointObjectiveResult:
    """One joint objective evaluation at shared ``(a, b)`` + per-case binding.

    ``objective`` is the **equal-weight mean** of the per-case objectives
    (each = in-window RMSE + escape penalty, per :class:`ObjectiveResult`).
    Per-case first, then mean: the windows differ in length, so averaging the
    per-case scalars prevents long-window dominance (METHOD_B §9.2).
    """

    a: float                            # amu/ps (0.0 under pure-cubic)
    b: float                            # amu*ps/A^2
    e_bind_eV: dict[str, float]         # per case; equal under shared variants
    objective: float                    # A/ps, equal-weight mean
    per_case: dict[str, ObjectiveResult]


@dataclass
class SharedTrajectoryMatchingFit:
    """A completed §9 Stage-2 joint fit (best of all starts) for one variant."""

    variant: str
    cases: list[str]
    best: JointObjectiveResult
    starts: list[JointObjectiveResult]  # best point of each optimizer start
    n_evaluations: int                  # joint evaluations (each = 2 ion runs)
    converged: bool
    seed: int
    n_molecules: int
    windows: dict[str, tuple[float, float]]
    reference_paths: dict[str, Path]
    meff_amu: float
    anchors: tuple[float, float]        # (a0, b0), 18A Method-A bundle only
    # Uncertainty (filled by the driver script, same semantics as the
    # per-case fit; the shared refit's band is sensitivity-only -- the §8
    # seed sweep showed the HeDFT-comparison presets are seed-insensitive):
    a_err: float | None = None
    b_err: float | None = None
    e_bind_err_eV: float | None = None
    uncertainty_model: str | None = None


def evaluate_joint_objective(
    a: float,
    b: float,
    e_bind_by_case_eV: dict[str, float],
    setups: dict[str, CaseSetup],
    *,
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> JointObjectiveResult:
    """Evaluate the §9 joint objective: per-case (RMSE + penalty), then mean.

    Parameters
    ----------
    a, b : float
        Shared ``linear_cubic`` coefficients [amu/ps, amu*ps/A^2]. ``a = 0``
        is the legal pure-cubic point (§9.5 guard relaxation).
    e_bind_by_case_eV : dict[str, float]
        Effective binding per case [eV]; equal values realize the shared-
        ``E_bind`` variants, distinct values the 4-parameter diagnostic.
    setups : dict[str, CaseSetup]
        Per-case contexts; keys must match ``e_bind_by_case_eV`` exactly.

    Returns
    -------
    JointObjectiveResult
        With ``objective`` the equal-weight mean over cases (a non-finite
        case contributes its large finite sentinel; never raises into the
        optimizer for that reason).
    """
    if set(e_bind_by_case_eV) != set(setups):
        raise ValueError(
            f"e_bind_by_case_eV cases {sorted(e_bind_by_case_eV)} do not "
            f"match setups {sorted(setups)}"
        )
    per_case = {
        case: evaluate_objective(
            (a, b, e_bind_by_case_eV[case]),
            setup,
            penalty_weight_Aps=penalty_weight_Aps,
            run_fn=run_fn,
        )
        for case, setup in setups.items()
    }
    objective = float(np.mean([r.objective for r in per_case.values()]))
    return JointObjectiveResult(
        a=float(a),
        b=float(b),
        e_bind_eV={c: float(e) for c, e in e_bind_by_case_eV.items()},
        objective=objective,
        per_case=per_case,
    )


def fit_shared_trajectory_matching(
    setups: dict[str, CaseSetup],
    *,
    variant: str,
    anchors: tuple[float, float],
    e_bind_min_eV: float = 0.01,
    e_bind_max_eV: float = E_BIND_STATIC_EV,
    prescan_e_bind_eV: Sequence[float] | None = None,
    extra_starts: bool = True,
    penalty_weight_Aps: float = 2.0,
    maxfev_per_start: int = 400,
    fatol_Aps: float = 1.0e-3,
    xatol_normalized: float = 1.0e-3,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> SharedTrajectoryMatchingFit:
    """Run the §9 Stage-2 joint Nelder-Mead fit across all cases in ``setups``.

    Parameter vector by ``variant`` (normalized coordinates, scale
    ``(a0, b0, 0.308 eV)``):

    * ``shared_3param`` -- ``(a, b, E)``, ``a`` free with **lower bound 0**
      (the §9 decision; the §9.5 guard accepts ``a = 0`` while ``b > 0``);
    * ``shared_pure_cubic`` -- ``(b, E)`` with ``a = 0`` fixed;
    * ``diagnostic_4param`` -- ``(a, b, E_case1, ..., E_caseK)`` in
      ``setups`` insertion order (run only if the primary fails its bands:
      localizes drag-law vs. binding generalization failure).

    Parameters
    ----------
    setups : dict[str, CaseSetup]
        Per-case contexts (insertion order fixes the parameter order). All
        cases must share ``seed``, ``num_molecules``, and ``m_eff_amu``
        (one joint provenance stamp; raises otherwise).
    anchors : tuple[float, float]
        Normalization anchors ``(a0, b0)``. Per §9.3 these must come from
        the **18 A Method-A bundle only** -- the 9 A Method-A artifact is
        provenance-broken and must not condition the fit, which is why the
        per-case ``setup.a0/b0`` are deliberately NOT used here.
    prescan_e_bind_eV : sequence of float, optional
        Coarse shared-``E_bind`` pre-scan grid at the anchor coefficients
        (``(a0, b0)``; ``(0, b0)`` under pure-cubic); default 8 points in
        [0.05, 0.305] eV.
    extra_starts : bool
        Adds the two ridge-probing starts (high-binding, coefficient-tilted)
        mirroring the per-case fit's multi-start scheme.

    Returns
    -------
    SharedTrajectoryMatchingFit
        ``best`` is the lowest-objective minimum across starts; ``starts``
        records each start's minimum so any residual ridge is visible.
        Uncertainty fields are left ``None`` (the driver script fills them).
    """
    if variant not in JOINT_VARIANTS:
        raise ValueError(
            f"variant must be one of {JOINT_VARIANTS}, got {variant!r}"
        )
    if len(setups) < 2:
        raise ValueError(
            f"the shared-form refit is a cross-case fit; need >= 2 cases, "
            f"got {sorted(setups)}"
        )
    cases = list(setups)
    first = setups[cases[0]]
    for case in cases[1:]:
        s = setups[case]
        if (
            s.cfg_base.seed != first.cfg_base.seed
            or s.cfg_base.num_molecules != first.cfg_base.num_molecules
            or s.cfg_base.m_eff_amu != first.cfg_base.m_eff_amu
        ):
            raise ValueError(
                "all cases must share seed, num_molecules, and m_eff_amu for "
                f"one joint provenance stamp; {cases[0]!r} vs {case!r} differ"
            )
    a0, b0 = (float(v) for v in anchors)
    if not (a0 > 0.0 and b0 > 0.0):
        raise ValueError(f"anchors must be positive, got {anchors!r}")

    e_lo_n = e_bind_min_eV / E_BIND_STATIC_EV
    e_hi_n = e_bind_max_eV / E_BIND_STATIC_EV
    if variant == SHARED_3PARAM:
        # a lower bound 0 -- the §9 decision (pure-cubic point reachable).
        lower = np.array([0.0, 1.0e-3, e_lo_n])
        upper = np.array([10.0, 10.0, e_hi_n])
        scale = np.array([a0, b0, E_BIND_STATIC_EV])
    elif variant == SHARED_PURE_CUBIC:
        lower = np.array([1.0e-3, e_lo_n])
        upper = np.array([10.0, e_hi_n])
        scale = np.array([b0, E_BIND_STATIC_EV])
    else:  # DIAGNOSTIC_4PARAM
        k = len(cases)
        lower = np.array([0.0, 1.0e-3] + [e_lo_n] * k)
        upper = np.array([10.0, 10.0] + [e_hi_n] * k)
        scale = np.array([a0, b0] + [E_BIND_STATIC_EV] * k)

    def unpack(x: np.ndarray) -> tuple[float, float, dict[str, float]]:
        theta = np.asarray(x, dtype=float) * scale
        if variant == SHARED_3PARAM:
            return theta[0], theta[1], {c: theta[2] for c in cases}
        if variant == SHARED_PURE_CUBIC:
            return 0.0, theta[0], {c: theta[1] for c in cases}
        return theta[0], theta[1], dict(zip(cases, theta[2:]))

    n_evals = 0

    def evaluate_x(x: np.ndarray) -> JointObjectiveResult:
        nonlocal n_evals
        n_evals += 1
        a, b, e_by_case = unpack(x)
        return evaluate_joint_objective(
            a, b, e_by_case, setups,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )

    # --- shared-E_bind pre-scan at the anchor coefficients ---
    if prescan_e_bind_eV is None:
        prescan_e_bind_eV = np.linspace(0.05, 0.305, 8)
    a_pre = 0.0 if variant == SHARED_PURE_CUBIC else a0
    prescan = [
        evaluate_joint_objective(
            a_pre, b0, {c: float(e) for c in cases}, setups,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )
        for e in prescan_e_bind_eV
    ]
    n_evals += len(prescan)
    e_start = min(prescan, key=lambda r: r.objective).e_bind_eV[cases[0]]
    e_start_n = e_start / E_BIND_STATIC_EV

    if variant == SHARED_3PARAM:
        x0s = [np.array([1.0, 1.0, e_start_n])]
        if extra_starts:
            x0s.append(np.array([1.0, 1.0, 0.9 * e_hi_n]))
            x0s.append(np.array([1.3, 0.7, e_start_n]))
    elif variant == SHARED_PURE_CUBIC:
        x0s = [np.array([1.0, e_start_n])]
        if extra_starts:
            x0s.append(np.array([1.0, 0.9 * e_hi_n]))
            x0s.append(np.array([0.7, e_start_n]))
    else:  # DIAGNOSTIC_4PARAM
        k = len(cases)
        x0s = [np.array([1.0, 1.0] + [e_start_n] * k)]
        if extra_starts:
            x0s.append(np.array([1.0, 1.0] + [0.9 * e_hi_n] * k))
            x0s.append(np.array([1.3, 0.7] + [e_start_n] * k))
    x0s = [np.clip(x0, lower, upper) for x0 in x0s]

    start_results: list[JointObjectiveResult] = []
    converged = True
    for x0 in x0s:
        res = minimize(
            lambda x: evaluate_x(x).objective,
            x0,
            method="Nelder-Mead",
            bounds=Bounds(lower, upper),
            options={
                "fatol": fatol_Aps,
                "xatol": xatol_normalized,
                "maxfev": maxfev_per_start,
            },
        )
        converged = converged and bool(res.success)
        start_results.append(evaluate_x(np.asarray(res.x, dtype=float)))

    best = min(start_results, key=lambda r: r.objective)
    return SharedTrajectoryMatchingFit(
        variant=variant,
        cases=cases,
        best=best,
        starts=start_results,
        n_evaluations=n_evals,
        converged=converged,
        seed=first.cfg_base.seed,
        n_molecules=first.cfg_base.num_molecules,
        windows={c: setups[c].window for c in cases},
        reference_paths={c: setups[c].reference_path for c in cases},
        meff_amu=float(first.cfg_base.m_eff_amu),
        anchors=(a0, b0),
    )


def joint_sensitivity_halfwidths(
    setups: dict[str, CaseSetup],
    best: JointObjectiveResult,
    *,
    rel_rise: float = 0.10,
    factors: Sequence[float] = (1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0),
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> tuple[float, float, float]:
    """RMSE-sensitivity half-widths for ``(a, b, E_shared)`` on the joint
    objective (shared-``E_bind`` variants only; same band semantics as
    :func:`sensitivity_halfwidths`). Under pure-cubic (``best.a == 0``) the
    ``a`` half-width is 0.0 by convention -- ``a`` was fixed, not fitted.

    Raises
    ------
    ValueError
        If ``best`` carries per-case (unequal) bindings -- the 4-parameter
        diagnostic has no single shared ``E`` to scan and never produces a
        production bundle, so its sensitivity is deliberately not defined.
    """
    e_values = sorted(set(best.e_bind_eV.values()))
    if len(e_values) != 1:
        raise ValueError(
            "joint sensitivity is defined only for shared-E_bind variants; "
            f"got per-case bindings {best.e_bind_eV}"
        )
    e_shared = e_values[0]

    def objective_at(theta: np.ndarray) -> float:
        return evaluate_joint_objective(
            theta[0], theta[1], {c: theta[2] for c in setups}, setups,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        ).objective

    return _objective_rise_halfwidths(  # type: ignore[return-value]
        objective_at,
        (best.a, best.b, e_shared),
        best.objective,
        rel_rise=rel_rise,
        factors=factors,
    )


def write_shared_fit_parameters(
    fit: SharedTrajectoryMatchingFit,
    out_dir: Path,
    *,
    stage: str = "stage2_joint",
) -> Path:
    """Write a shared-form ``fit_parameters.json`` bundle (METHOD_B §9.6).

    The schema is a superset of the per-case Method-B bundle and loads
    through :func:`i2_helium_md.presets.load_drag_coefficients` unchanged
    (same loader contract), with the §9.6 provenance additions:
    ``calibration_cases``, ``stage``, ``variant``, the per-case windows /
    references / scores, and the anchor provenance. The loader-required
    ``t_start``/``t_end`` span the union of the per-case windows (the
    per-case windows are in ``windows``); ``reference_file`` lists both
    per-case references.

    Both shared variants are recorded to their own directories by the driver
    (the production-candidate choice is deferred to the preset-rewiring
    decision); the 4-parameter diagnostic is refused -- it has no single
    jointly-validated binding and must never become a loadable bundle.

    Raises
    ------
    ValueError
        On the diagnostic variant, a trapped fit in ANY case, or unfilled
        uncertainty fields.
    """
    if fit.variant == DIAGNOSTIC_4PARAM:
        raise ValueError(
            "refusing to write a loadable bundle for the 4-parameter "
            "diagnostic: per-case bindings are a failure-localization "
            "device, never a production pairing (METHOD_B §9.2)"
        )
    trapped = {
        c: r.escape_fraction
        for c, r in fit.best.per_case.items()
        if r.penalty_Aps != 0.0
    }
    if trapped:
        raise ValueError(
            "refusing to write a trapped fit: escape_fraction by case "
            f"{trapped} (penalty != 0); the binding<->drag pair must permit "
            "full escape in every case"
        )
    if fit.a_err is None or fit.b_err is None or fit.uncertainty_model is None:
        raise ValueError(
            "uncertainty fields (a_err, b_err, uncertainty_model) must be "
            "filled before writing the bundle"
        )
    e_values = sorted(set(fit.best.e_bind_eV.values()))
    assert len(e_values) == 1  # guaranteed: diagnostic refused above
    e_shared = e_values[0]

    repo_root = REFERENCE_DRAG_ROOT.parents[2]
    payload = {
        "extraction_method": "trajectory_matching",
        "form": LINEAR_CUBIC,
        "a": fit.best.a,
        "b": fit.best.b,
        "a_err": fit.a_err,
        "b_err": fit.b_err,
        "uncertainty_model": fit.uncertainty_model,
        "meff_amu": fit.meff_amu,
        "extraction_mass_model": "constant",
        "effective_binding_energy_I_ion_eV": e_shared,
        "effective_binding_energy_err_eV": fit.e_bind_err_eV,
        # Loader-required window keys: union span; per-case in "windows".
        "t_start": min(w[0] for w in fit.windows.values()),
        "t_end": max(w[1] for w in fit.windows.values()),
        "reference_file": [
            _reference_file_str(fit.reference_paths[c], repo_root)
            for c in fit.cases
        ],
        # --- §9.6 shared-form provenance ---
        "calibration_cases": fit.cases,
        "stage": stage,
        "variant": fit.variant,
        "anchors": {
            "a0": fit.anchors[0],
            "b0": fit.anchors[1],
            "provenance": (
                "18A Method-A bundle only (9A Method-A artifact "
                "provenance-broken, drag_migration_log.md 2026-06-11)"
            ),
        },
        "objective": (
            "equal_weight_mean_over_cases_of_in_window_rmse_mean_abs_v2_"
            "plus_escape_penalty"
        ),
        "objective_mean_Aps": fit.best.objective,
        "windows": {c: list(fit.windows[c]) for c in fit.cases},
        "per_case": {
            c: {
                "rmse_Aps": r.rmse_Aps,
                "escape_fraction": r.escape_fraction,
                "objective_Aps": r.objective,
                "num_overlap_points": r.num_overlap_points,
            }
            for c, r in fit.best.per_case.items()
        },
        "n_molecules_fit": fit.n_molecules,
        "seed": fit.seed,
        "optimizer": {
            "method": "nelder-mead",
            "n_starts": len(fit.starts),
            "evals": fit.n_evaluations,
            "converged": fit.converged,
            "start_minima": [
                {
                    "a": r.a, "b": r.b, "e_bind_eV": r.e_bind_eV,
                    "objective_Aps": r.objective,
                }
                for r in fit.starts
            ],
        },
        "date": _datetime.date.today().isoformat(),
        "branch": _git_branch(repo_root),
        "flags": {
            "radial_projection_convention": True,
            # 9A (genuinely non-radial) is INSIDE the Stage-2 joint
            # objective -- the §3.5 accepted risk, carried by the bundle.
            "transverse_contaminated_non_radial_reference": True,
            "full_window_heldout_window_axis_forfeited": True,
            "cross_case_axis_consumed_by_joint_fit": True,
        },
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fit_parameters.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return out_path
