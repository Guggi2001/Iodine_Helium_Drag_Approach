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

The **§10 alternative-form discrimination** machinery
(:func:`evaluate_form_objective`, :func:`fit_shared_form_trajectory_matching`
with the ``lq_shared_3param`` / ``lq_shared_pure_quadratic`` /
``pl_shared_3param`` variants, the single-case Stage-1-analog
:func:`fit_form_trajectory_matching`, :func:`form_sensitivity_halfwidths`,
:func:`write_shared_form_fit_parameters`) runs the ``linear_quadratic`` and
``power_law`` families through the same joint objective against the
pure-cubic incumbent; the ``power_law`` optimizer works in the locked
pivot parameterization ``(gamma_ref, n)`` (§10.4.1) and stamps raw
``{C, n}``. Drivers: ``scripts/extraction/method_b_form_refit_*.py``.
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
from ..physics.drag import (
    DragCoefficients,
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    _REQUIRED_COEFF_KEYS,
)
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


def _run_and_score(
    coeffs: DragCoefficients,
    e_bind: float,
    setup: CaseSetup,
    *,
    penalty_weight_Aps: float,
    run_fn: Callable[..., IonCheckpoint],
) -> tuple[float, float, float, float, int]:
    """Forward-integrate one candidate and score it (the shared objective core).

    Builds a config that differs from ``setup.cfg_base`` ONLY in the
    ion-stage fields (candidate form + coefficients + binding, stamped as an
    exact §6.5.1 pair), forward-integrates the ion stage from the cached
    neutral checkpoint, and scores the ensemble-mean ``|v2|`` against the
    smoothed reference over ``setup.window``.

    Returns ``(objective, rmse_Aps, escape_fraction, penalty_Aps,
    num_overlap_points)`` with the large finite sentinel as ``objective`` if
    the trajectory or score is non-finite (never raises into the optimizer
    for that reason).
    """
    cfg = replace(
        setup.cfg_base,
        drag_form=coeffs.form,
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

    return objective, rmse, escape, penalty, int(comparison.num_overlap_points)


def evaluate_objective(
    theta: Sequence[float],
    setup: CaseSetup,
    *,
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> ObjectiveResult:
    """Evaluate the Method-B objective at ``theta = (a, b, e_bind_eV)``.

    ``linear_cubic`` front-end of :func:`_run_and_score` (the §10 form-phase
    families evaluate through :func:`evaluate_form_objective` instead).

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
    objective, rmse, escape, penalty, n_overlap = _run_and_score(
        coeffs, e_bind, setup,
        penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
    )

    return ObjectiveResult(
        a=a,
        b=b,
        e_bind_eV=e_bind,
        objective=objective,
        rmse_Aps=rmse,
        escape_fraction=escape,
        penalty_Aps=penalty,
        num_overlap_points=n_overlap,
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


# ===========================================================================
# Alternative drag-form discrimination (METHOD_B §10)
# ===========================================================================
# Form-generic analog of the §9 shared machinery: the linear_quadratic and
# power_law families run through the SAME joint objective (per-case RMSE +
# escape penalty, then equal-weight mean), same windows, same N/seed/anchoring
# discipline, compared against the pure-cubic incumbent under the reused §9.4
# bands plus the locked two-threshold T_form scheme (§10.4.1). The §9
# linear_cubic machinery above is the frozen incumbent record and stays
# untouched; this section adds the families BESIDE it.
#
# power_law optimizer parameterization (§10.4.1, locked): the optimizer works
# in the PIVOT space ``(gamma_ref, n)`` with ``gamma_ref = C * v_ref**(n-1)``
# at the locked pivot speed ``V_REF_APS = 3.0`` (passed in via ``anchors``) --
# it axis-aligns the ``log C ~ const - n*log(v)`` matching ridge so the fitted
# n-hat's sensitivity half-width is meaningful. Optimizer-internal only: every
# result and bundle stamps the raw closed-form ``{C, n}``.

LQ_SHARED_3PARAM = "lq_shared_3param"                # {a, c, E}, a >= 0
LQ_SHARED_PURE_QUADRATIC = "lq_shared_pure_quadratic"  # a == 0, fit {c, E}
PL_SHARED_3PARAM = "pl_shared_3param"                # {C, n, E}, n in bounds
FORM_VARIANTS = (LQ_SHARED_3PARAM, LQ_SHARED_PURE_QUADRATIC, PL_SHARED_3PARAM)

_VARIANT_FORM = {
    LQ_SHARED_3PARAM: LINEAR_QUADRATIC,
    LQ_SHARED_PURE_QUADRATIC: LINEAR_QUADRATIC,
    PL_SHARED_3PARAM: POWER_LAW,
}

# Anchor keys required per family (pre-registered constants from the driver's
# USER SETTINGS, METHOD_B §10.4.1 -- never a runtime bundle read; the re-wired
# presets carry a0 = 0, which must not condition any fit).
_VARIANT_ANCHOR_KEYS = {
    LQ_SHARED_3PARAM: ("a0", "c0"),
    LQ_SHARED_PURE_QUADRATIC: ("a0", "c0"),
    PL_SHARED_3PARAM: ("C0", "n0", "v_ref_Aps"),
}


def _pivot_to_C(gamma_ref: float, n: float, v_ref_Aps: float) -> float:
    """Invert the §10.4.1 pivot: ``C = gamma_ref / v_ref**(n-1)``."""
    return float(gamma_ref) / float(v_ref_Aps) ** (float(n) - 1.0)


def _C_to_pivot(C: float, n: float, v_ref_Aps: float) -> float:
    """The §10.4.1 pivot coordinate: ``gamma_ref = C * v_ref**(n-1)``."""
    return float(C) * float(v_ref_Aps) ** (float(n) - 1.0)


@dataclass(frozen=True)
class FormObjectiveResult:
    """One single-case objective evaluation of a form-phase candidate.

    Same objective semantics as :class:`ObjectiveResult`
    (``objective = rmse_Aps + penalty_Aps``); the coefficients are carried as
    the form's raw closed-form dict (``{a, c}`` or ``{C, n}``).
    """

    form: str
    coefficients: dict[str, float]
    e_bind_eV: float
    objective: float          # A/ps (rmse + penalty; sentinel if non-finite)
    rmse_Aps: float
    escape_fraction: float
    penalty_Aps: float
    num_overlap_points: int


@dataclass(frozen=True)
class JointFormObjectiveResult:
    """One joint objective evaluation of a form-phase candidate.

    ``objective`` is the equal-weight mean of the per-case objectives, exactly
    as :class:`JointObjectiveResult` (per-case first prevents long-window
    dominance, §9.2 -- reused unchanged by §10.4).
    """

    form: str
    coefficients: dict[str, float]
    e_bind_eV: dict[str, float]         # per case; equal under shared variants
    objective: float                    # A/ps, equal-weight mean
    per_case: dict[str, FormObjectiveResult]


@dataclass
class SharedFormTrajectoryMatchingFit:
    """A completed §10 form-phase fit (best of all starts) for one variant.

    Produced by :func:`fit_shared_form_trajectory_matching` (Stage 2, >= 2
    cases) or :func:`fit_form_trajectory_matching` (the Stage-1 analog's
    single-case 18 A-only fit -- record-only, never written as a bundle).
    """

    variant: str
    form: str
    cases: list[str]
    best: JointFormObjectiveResult
    starts: list[JointFormObjectiveResult]
    n_evaluations: int                  # joint evaluations (1 ion run/case)
    converged: bool
    seed: int
    n_molecules: int
    windows: dict[str, tuple[float, float]]
    reference_paths: dict[str, Path]
    meff_amu: float
    anchors: dict[str, float]           # pre-registered constants (§10.4.1)
    # Uncertainty (filled by the driver script; sensitivity-only per the §8
    # seed-insensitivity finding, same convention as the §9 shared fit):
    coeff_errs: dict[str, float] | None = None
    e_bind_err_eV: float | None = None
    uncertainty_model: str | None = None


def evaluate_form_objective(
    form: str,
    coefficients: dict[str, float],
    e_bind_eV: float,
    setup: CaseSetup,
    *,
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> FormObjectiveResult:
    """Evaluate the Method-B objective for a form-phase candidate (one case).

    Form-generic front-end of :func:`_run_and_score`: builds the candidate
    bundle of the given ``form`` (raw closed-form ``coefficients``), stamps
    the exact §6.5.1 pairing, forward-integrates, and scores -- identical
    objective semantics to :func:`evaluate_objective`.

    Parameters
    ----------
    form : str
        A realised drag form (``linear_quadratic`` or ``power_law``;
        ``linear_cubic`` is also accepted for cross-checks against
        :func:`evaluate_objective`).
    coefficients : dict[str, float]
        The form's raw coefficients (``{a, c}`` / ``{C, n}`` / ``{a, b}``).
        power_law candidates pass raw ``{C, n}`` -- the pivot is
        optimizer-internal and never reaches this function.
    e_bind_eV : float
        Candidate effective binding [eV].
    """
    coeffs = DragCoefficients(
        form=form,
        coefficients={k: float(v) for k, v in coefficients.items()},
        extraction_mass_model="constant",
        extraction_mass_amu=setup.cfg_base.m_eff_amu,
        extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=float(e_bind_eV),
    )
    objective, rmse, escape, penalty, n_overlap = _run_and_score(
        coeffs, float(e_bind_eV), setup,
        penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
    )
    return FormObjectiveResult(
        form=form,
        coefficients=dict(coeffs.coefficients),
        e_bind_eV=float(e_bind_eV),
        objective=objective,
        rmse_Aps=rmse,
        escape_fraction=escape,
        penalty_Aps=penalty,
        num_overlap_points=n_overlap,
    )


def evaluate_joint_form_objective(
    form: str,
    coefficients: dict[str, float],
    e_bind_by_case_eV: dict[str, float],
    setups: dict[str, CaseSetup],
    *,
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> JointFormObjectiveResult:
    """Joint objective for a form-phase candidate: per-case score, then mean.

    The form-generic analog of :func:`evaluate_joint_objective` (same
    per-case-first equal-weight-mean construction, §9.2/§10.4).
    """
    if set(e_bind_by_case_eV) != set(setups):
        raise ValueError(
            f"e_bind_by_case_eV cases {sorted(e_bind_by_case_eV)} do not "
            f"match setups {sorted(setups)}"
        )
    per_case = {
        case: evaluate_form_objective(
            form, coefficients, e_bind_by_case_eV[case], setup,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )
        for case, setup in setups.items()
    }
    objective = float(np.mean([r.objective for r in per_case.values()]))
    return JointFormObjectiveResult(
        form=form,
        coefficients={k: float(v) for k, v in coefficients.items()},
        e_bind_eV={c: float(e) for c, e in e_bind_by_case_eV.items()},
        objective=objective,
        per_case=per_case,
    )


def _check_form_anchors(variant: str, anchors: dict[str, float]) -> None:
    """Validate the pre-registered anchor dict for a form-phase variant."""
    required = _VARIANT_ANCHOR_KEYS[variant]
    missing = [k for k in required if k not in anchors]
    if missing:
        raise ValueError(
            f"variant {variant!r} requires anchors {required}; missing "
            f"{missing} (pre-registered §10.4.1 constants -- never a runtime "
            f"bundle read)"
        )
    nonpositive = {k: anchors[k] for k in required if not (anchors[k] > 0.0)}
    if nonpositive:
        raise ValueError(
            f"form-phase anchors must be positive, got {nonpositive} (the "
            f"re-wired shared bundle's a0 = 0 must not condition a fit)"
        )


def fit_shared_form_trajectory_matching(
    setups: dict[str, CaseSetup],
    *,
    variant: str,
    anchors: dict[str, float],
    n_bounds: tuple[float, float] = (1.0, 4.0),
    e_bind_min_eV: float = 0.01,
    e_bind_max_eV: float = E_BIND_STATIC_EV,
    prescan_e_bind_eV: Sequence[float] | None = None,
    extra_starts: bool = True,
    penalty_weight_Aps: float = 2.0,
    maxfev_per_start: int = 400,
    fatol_Aps: float = 1.0e-3,
    xatol_normalized: float = 1.0e-3,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> SharedFormTrajectoryMatchingFit:
    """Run the §10 Stage-2 shared joint fit for one form-phase variant.

    Same optimizer discipline as :func:`fit_shared_trajectory_matching`
    (normalized bounded Nelder-Mead, shared-``E_bind`` pre-scan at the anchor
    coefficients, ridge-probing multi-starts), with the family-specific
    parameter vector (normalized coordinates):

    * ``lq_shared_3param`` -- ``(a/a0, c/c0, E/0.308)``, ``a`` lower bound
      **0** (pure-quadratic point reachable; §10.3 guard arm);
    * ``lq_shared_pure_quadratic`` -- ``(c/c0, E/0.308)`` with ``a = 0``
      fixed;
    * ``pl_shared_3param`` -- ``(gamma_ref/gamma_ref0, n, E/0.308)`` in the
      §10.4.1 PIVOT parameterization (``gamma_ref0 = C0*v_ref**(n0-1)``
      computed from the anchors); ``n`` direct-bounded by ``n_bounds``
      (locked ``[1, 4]``: nesting points 2 and 3 interior, ``n >= 1`` keeps
      gamma finite at rest). Results stamp raw ``{C, n}``.

    Parameters
    ----------
    setups : dict[str, CaseSetup]
        Per-case contexts. All cases must share ``seed``, ``num_molecules``,
        and ``m_eff_amu`` (one joint provenance stamp; raises otherwise).
        Requires >= 2 cases; the Stage-1 analog's single-case fit goes
        through :func:`fit_form_trajectory_matching`.
    variant : str
        One of :data:`FORM_VARIANTS`.
    anchors : dict[str, float]
        Pre-registered normalization constants (§10.4.1): ``{a0, c0}`` for
        the lq variants, ``{C0, n0, v_ref_Aps}`` for power_law. Positive;
        validated. NEVER read from a loaded bundle at runtime.
    n_bounds : tuple of 2 floats
        power_law exponent bounds (locked ``(1, 4)``); ignored by the lq
        variants.

    Returns
    -------
    SharedFormTrajectoryMatchingFit
        ``best`` is the lowest-objective minimum across starts; ``starts``
        records each start's minimum so any residual ridge is visible.
        Uncertainty fields are left ``None`` (the driver script fills them).
    """
    if len(setups) < 2:
        raise ValueError(
            f"the shared-form refit is a cross-case fit; need >= 2 cases, "
            f"got {sorted(setups)} (the Stage-1 analog's single-case fit is "
            f"fit_form_trajectory_matching)"
        )
    return _fit_form_engine(
        setups,
        variant=variant,
        anchors=anchors,
        n_bounds=n_bounds,
        e_bind_min_eV=e_bind_min_eV,
        e_bind_max_eV=e_bind_max_eV,
        prescan_e_bind_eV=prescan_e_bind_eV,
        extra_starts=extra_starts,
        penalty_weight_Aps=penalty_weight_Aps,
        maxfev_per_start=maxfev_per_start,
        fatol_Aps=fatol_Aps,
        xatol_normalized=xatol_normalized,
        run_fn=run_fn,
    )


def fit_form_trajectory_matching(
    setup: CaseSetup,
    *,
    variant: str,
    anchors: dict[str, float],
    **kwargs,
) -> SharedFormTrajectoryMatchingFit:
    """Single-case form-phase fit (the §10.4 Stage-1 analog).

    Fits the family's full variant on ONE case (18 A by convention) so the
    untouched other case can be scored as a prediction. Record-only output:
    the §10 decision (2026-06-12 session) is that Stage-1-analog parameters
    are recorded in the stage-1 JSON, never written as a loadable bundle.
    Same engine, same objective, same optimizer discipline as the shared
    fit; the "joint" objective over one case reduces to that case's
    objective.
    """
    return _fit_form_engine(
        {setup.case: setup}, variant=variant, anchors=anchors, **kwargs
    )


def _fit_form_engine(
    setups: dict[str, CaseSetup],
    *,
    variant: str,
    anchors: dict[str, float],
    n_bounds: tuple[float, float] = (1.0, 4.0),
    e_bind_min_eV: float = 0.01,
    e_bind_max_eV: float = E_BIND_STATIC_EV,
    prescan_e_bind_eV: Sequence[float] | None = None,
    extra_starts: bool = True,
    penalty_weight_Aps: float = 2.0,
    maxfev_per_start: int = 400,
    fatol_Aps: float = 1.0e-3,
    xatol_normalized: float = 1.0e-3,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> SharedFormTrajectoryMatchingFit:
    """Shared optimizer engine behind the two public form-phase fitters."""
    if variant not in FORM_VARIANTS:
        raise ValueError(
            f"variant must be one of {FORM_VARIANTS}, got {variant!r}"
        )
    _check_form_anchors(variant, anchors)
    form = _VARIANT_FORM[variant]

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

    e_lo_n = e_bind_min_eV / E_BIND_STATIC_EV
    e_hi_n = e_bind_max_eV / E_BIND_STATIC_EV

    # --- family-specific normalized parameter vector ---
    if variant == LQ_SHARED_3PARAM:
        a0, c0 = anchors["a0"], anchors["c0"]
        scale = np.array([a0, c0, E_BIND_STATIC_EV])
        lower = np.array([0.0, 1.0e-3, e_lo_n])   # a lower bound 0 (§10.3)
        upper = np.array([10.0, 10.0, e_hi_n])

        def unpack(x: np.ndarray) -> tuple[dict[str, float], float]:
            theta = np.asarray(x, dtype=float) * scale
            return {"a": float(theta[0]), "c": float(theta[1])}, float(theta[2])

        anchor_coeffs = {"a": a0, "c": c0}
    elif variant == LQ_SHARED_PURE_QUADRATIC:
        c0 = anchors["c0"]
        scale = np.array([c0, E_BIND_STATIC_EV])
        lower = np.array([1.0e-3, e_lo_n])
        upper = np.array([10.0, e_hi_n])

        def unpack(x: np.ndarray) -> tuple[dict[str, float], float]:
            theta = np.asarray(x, dtype=float) * scale
            return {"a": 0.0, "c": float(theta[0])}, float(theta[1])

        anchor_coeffs = {"a": 0.0, "c": c0}
    else:  # PL_SHARED_3PARAM -- pivot parameterization (§10.4.1)
        C0, n0, v_ref = anchors["C0"], anchors["n0"], anchors["v_ref_Aps"]
        n_lo, n_hi = (float(v) for v in n_bounds)
        if not (n_lo >= 1.0 and n_hi > n_lo):
            raise ValueError(
                f"n_bounds must satisfy 1 <= n_lo < n_hi (gamma finite at "
                f"rest, §10.3), got {n_bounds!r}"
            )
        if not (n_lo <= n0 <= n_hi):
            raise ValueError(
                f"anchor n0={n0} outside n_bounds {n_bounds!r}"
            )
        gamma_ref0 = _C_to_pivot(C0, n0, v_ref)
        # x = (gamma_ref/gamma_ref0, n, E/0.308): n is O(1) in [1, 4], so it
        # is carried raw (scale 1) -- already well-conditioned.
        scale = np.array([gamma_ref0, 1.0, E_BIND_STATIC_EV])
        lower = np.array([1.0e-3, n_lo, e_lo_n])
        upper = np.array([10.0, n_hi, e_hi_n])

        def unpack(x: np.ndarray) -> tuple[dict[str, float], float]:
            theta = np.asarray(x, dtype=float) * scale
            gamma_ref, n = float(theta[0]), float(theta[1])
            return (
                {"C": _pivot_to_C(gamma_ref, n, v_ref), "n": n},
                float(theta[2]),
            )

        anchor_coeffs = {"C": C0, "n": n0}

    n_evals = 0

    def evaluate_x(x: np.ndarray) -> JointFormObjectiveResult:
        nonlocal n_evals
        n_evals += 1
        coefficients, e_bind = unpack(x)
        return evaluate_joint_form_objective(
            form, coefficients, {c: e_bind for c in cases}, setups,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )

    # --- shared-E_bind pre-scan at the anchor coefficients ---
    if prescan_e_bind_eV is None:
        prescan_e_bind_eV = np.linspace(0.05, 0.305, 8)
    prescan = [
        evaluate_joint_form_objective(
            form, anchor_coeffs, {c: float(e) for c in cases}, setups,
            penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
        )
        for e in prescan_e_bind_eV
    ]
    n_evals += len(prescan)
    e_start = min(prescan, key=lambda r: r.objective).e_bind_eV[cases[0]]
    e_start_n = e_start / E_BIND_STATIC_EV

    # --- multi-starts (the §9 discipline adapted per family; conditioning
    # only -- recorded via the optimizer stamp, not pre-registered) ---
    if variant == LQ_SHARED_3PARAM:
        x0s = [np.array([1.0, 1.0, e_start_n])]
        if extra_starts:
            x0s.append(np.array([1.0, 1.0, 0.9 * e_hi_n]))
            x0s.append(np.array([1.3, 0.7, e_start_n]))
    elif variant == LQ_SHARED_PURE_QUADRATIC:
        x0s = [np.array([1.0, e_start_n])]
        if extra_starts:
            x0s.append(np.array([1.0, 0.9 * e_hi_n]))
            x0s.append(np.array([0.7, e_start_n]))
    else:  # PL_SHARED_3PARAM: probe the exponent axis at the nesting points
        n0 = anchors["n0"]
        x0s = [np.array([1.0, n0, e_start_n])]
        if extra_starts:
            x0s.append(np.array([1.0, n0, 0.9 * e_hi_n]))
            x0s.append(np.array([1.0, 2.0, e_start_n]))
            x0s.append(np.array([1.0, 3.0, e_start_n]))
    x0s = [np.clip(x0, lower, upper) for x0 in x0s]

    start_results: list[JointFormObjectiveResult] = []
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
    return SharedFormTrajectoryMatchingFit(
        variant=variant,
        form=form,
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
        anchors={k: float(v) for k, v in anchors.items()},
    )


def form_sensitivity_halfwidths(
    setups: dict[str, CaseSetup],
    best: JointFormObjectiveResult,
    *,
    v_ref_Aps: float | None = None,
    rel_rise: float = 0.10,
    factors: Sequence[float] = (1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0),
    penalty_weight_Aps: float = 2.0,
    run_fn: Callable[..., IonCheckpoint] = run_ion_propagation,
) -> dict[str, float]:
    """RMSE-sensitivity half-widths for a form-phase joint optimum.

    Same crossing semantics as :func:`joint_sensitivity_halfwidths` (via the
    shared :func:`_objective_rise_halfwidths` core), with the family-specific
    parameter vector:

    * ``linear_quadratic`` -- scans ``(a, c, E)``; a fixed ``a = 0``
      (pure-quadratic) gets half-width 0.0 by convention. Returns
      ``{"a": .., "c": .., "e_bind_eV": ..}``.
    * ``power_law`` -- scans in the PIVOT space ``(gamma_ref, n, E)`` (the
      §10.4.1 decision: the pivot axis-aligns the matching ridge so the
      n-hat half-width is the meaningful exponent-identifiability number).
      Requires ``v_ref_Aps``. Returns ``{"gamma_ref": .., "n": ..,
      "e_bind_eV": ..}``; the driver converts ``gamma_ref`` to a ``C`` band
      at fixed ``n`` (``C_err = gamma_ref_err / v_ref**(n-1)``).

    Raises
    ------
    ValueError
        On per-case (unequal) bindings, or a power_law scan without
        ``v_ref_Aps``.
    """
    e_values = sorted(set(best.e_bind_eV.values()))
    if len(e_values) != 1:
        raise ValueError(
            "form sensitivity is defined only for shared-E_bind fits; "
            f"got per-case bindings {best.e_bind_eV}"
        )
    e_shared = e_values[0]

    if best.form == LINEAR_QUADRATIC:
        names = ("a", "c", "e_bind_eV")
        theta_best = (
            best.coefficients["a"], best.coefficients["c"], e_shared,
        )

        def objective_at(theta: np.ndarray) -> float:
            return evaluate_joint_form_objective(
                LINEAR_QUADRATIC,
                {"a": float(theta[0]), "c": float(theta[1])},
                {c: float(theta[2]) for c in setups},
                setups,
                penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
            ).objective

    elif best.form == POWER_LAW:
        if v_ref_Aps is None:
            raise ValueError(
                "power_law sensitivity scans in the pivot space and needs "
                "v_ref_Aps (the locked V_REF_APS pivot speed, §10.4.1)"
            )
        names = ("gamma_ref", "n", "e_bind_eV")
        C, n = best.coefficients["C"], best.coefficients["n"]
        theta_best = (_C_to_pivot(C, n, v_ref_Aps), n, e_shared)

        def objective_at(theta: np.ndarray) -> float:
            gamma_ref, n_trial = float(theta[0]), float(theta[1])
            return evaluate_joint_form_objective(
                POWER_LAW,
                {"C": _pivot_to_C(gamma_ref, n_trial, v_ref_Aps),
                 "n": n_trial},
                {c: float(theta[2]) for c in setups},
                setups,
                penalty_weight_Aps=penalty_weight_Aps, run_fn=run_fn,
            ).objective

    else:
        raise ValueError(
            f"form sensitivity covers the §10 families, got {best.form!r}"
        )

    halfwidths = _objective_rise_halfwidths(
        objective_at,
        theta_best,
        best.objective,
        rel_rise=rel_rise,
        factors=factors,
    )
    return dict(zip(names, halfwidths))


def write_shared_form_fit_parameters(
    fit: SharedFormTrajectoryMatchingFit,
    out_dir: Path,
    *,
    anchor_provenance: str,
    stage: str = "stage2_joint",
) -> Path:
    """Write a §10 form-phase ``fit_parameters.json`` bundle.

    Same loader contract as :func:`write_shared_fit_parameters` (the bundle
    loads through :func:`i2_helium_md.presets.load_drag_coefficients`): the
    form's raw coefficients + per-key ``<k>_err`` bands, the Method-B
    provenance keys, and the §9.6 shared-form provenance superset (variant
    disambiguates the family). power_law bundles additionally record the
    pivot block (``v_ref_Aps``, ``gamma_ref``) for re-derivability.

    Refusals mirror the §9 writer: a single-case (Stage-1 analog) fit is
    refused (record-only by the 2026-06-12 decision), as are trapped fits
    and unfilled uncertainty fields.

    Parameters
    ----------
    anchor_provenance : str
        Human-readable origin of the pre-registered anchors (named-constant
        provenance from the driver's USER SETTINGS, §10.4.1).
    """
    if len(fit.cases) < 2:
        raise ValueError(
            "refusing to write a single-case (Stage-1 analog) bundle: "
            "Stage-1-analog fits are record-only (2026-06-12 decision); "
            "only Stage-2 shared fits become loadable bundles"
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
    coeff_keys = _REQUIRED_COEFF_KEYS[fit.form]
    if (
        fit.coeff_errs is None
        or fit.uncertainty_model is None
        or any(k not in fit.coeff_errs for k in coeff_keys)
    ):
        raise ValueError(
            f"uncertainty fields (coeff_errs for {coeff_keys}, "
            "uncertainty_model) must be filled before writing the bundle"
        )
    e_values = sorted(set(fit.best.e_bind_eV.values()))
    assert len(e_values) == 1  # all §10 variants are shared-E_bind
    e_shared = e_values[0]

    repo_root = REFERENCE_DRAG_ROOT.parents[2]
    payload = {
        "extraction_method": "trajectory_matching",
        "form": fit.form,
        **{k: fit.best.coefficients[k] for k in coeff_keys},
        **{f"{k}_err": fit.coeff_errs[k] for k in coeff_keys},
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
        # --- shared-form provenance (§9.6 superset, §10 variant) ---
        "calibration_cases": fit.cases,
        "stage": stage,
        "variant": fit.variant,
        "anchors": {**fit.anchors, "provenance": anchor_provenance},
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
                    "coefficients": r.coefficients,
                    "e_bind_eV": r.e_bind_eV,
                    "objective_Aps": r.objective,
                }
                for r in fit.starts
            ],
        },
        "date": _datetime.date.today().isoformat(),
        "branch": _git_branch(repo_root),
        "flags": {
            "radial_projection_convention": True,
            # 9A (genuinely non-radial) is INSIDE the joint objective -- the
            # §3.5 accepted risk, carried by the bundle.
            "transverse_contaminated_non_radial_reference": True,
            "full_window_heldout_window_axis_forfeited": True,
            "cross_case_axis_consumed_by_joint_fit": True,
            # §10.4: these comparisons re-use the two seen trajectories --
            # model selection under a pre-registered protocol, NOT fresh
            # held-out validation (the winner's external test is VMI,
            # post-Tier-1).
            "model_selection_on_seen_data": True,
        },
    }
    if fit.form == POWER_LAW:
        v_ref = fit.anchors["v_ref_Aps"]
        C, n = fit.best.coefficients["C"], fit.best.coefficients["n"]
        payload["pivot"] = {
            "v_ref_Aps": v_ref,
            "gamma_ref": _C_to_pivot(C, n, v_ref),
            "note": (
                "optimizer-internal parameterization gamma_ref = "
                "C*v_ref**(n-1) (§10.4.1); the stamped law is the raw "
                "closed-form {C, n}"
            ),
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fit_parameters.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return out_path


def write_per_case_form_fit_parameters(
    fit: SharedFormTrajectoryMatchingFit,
    out_dir: Path,
    *,
    transverse_contaminated: bool,
    anchor_provenance: str,
    v_ref_Aps: float | None = None,
) -> Path:
    """Write a §10.8 per-case (single-case) form ``fit_parameters.json`` bundle.

    The per-case analog of :func:`write_shared_form_fit_parameters`: it writes a
    single-case (9 A-only or 18 A-only) ``linear_quadratic`` / ``power_law`` fit
    as a loadable bundle, filling the missing diagonal that §8 already gave
    ``linear_cubic``. Same loader contract (the bundle loads through
    :func:`i2_helium_md.presets.load_drag_coefficients`): the form's raw
    coefficients + per-key ``<k>_err``, the Method-B provenance keys, and the
    single-case window / reference. ``power_law`` additionally stamps the pivot
    block (``v_ref_Aps``, ``gamma_ref``) for re-derivability.

    These are diagnostic, **calibrated-not-validated** single-case fits: for
    these forms the §9 cross-case held-out axis is already spent, and a
    single-case fit trivially matches its own curve (METHOD_B §5/§8 warning).
    They are NEVER preset-wired -- the incumbent stays ``shared_pure_cubic``
    (§10.8); the honesty flags carry that status in the artifact.

    This is a distinct entry point from :func:`write_shared_form_fit_parameters`
    (which deliberately refuses single-case fits, the 2026-06-12 "Stage-1 analog
    is record-only" decision); that refusal is left intact.

    Parameters
    ----------
    fit : SharedFormTrajectoryMatchingFit
        A completed SINGLE-case form fit (``len(fit.cases) == 1``), as produced
        by :func:`fit_form_trajectory_matching`. ``coeff_errs`` /
        ``uncertainty_model`` must be filled (the driver computes them from the
        sensitivity scan).
    out_dir : Path
        Target directory (created if needed); the §10.8 convention is
        ``REFERENCE_DRAG_ROOT / case / "trajectory_matching" / variant``.
    transverse_contaminated : bool
        The standing 9 A flag: True for the genuinely non-radial 9 A reference,
        False for clean-radial 18 A (same semantics as
        :func:`write_fit_parameters`).
    anchor_provenance : str
        Human-readable origin of the pre-registered §10.4.1 anchors (the anchors
        are conditioning-only; carried for provenance).
    v_ref_Aps : float, optional
        REQUIRED for ``power_law`` -- the locked pivot speed (§10.4.1) stamped in
        the pivot block. Ignored for ``linear_quadratic``.

    Raises
    ------
    ValueError
        If the fit is not single-case, carries an escape penalty (trapped), has
        unfilled uncertainty fields, or is ``power_law`` without ``v_ref_Aps``.
    """
    if len(fit.cases) != 1:
        raise ValueError(
            "write_per_case_form_fit_parameters writes a SINGLE-case bundle "
            f"(the §10.8 per-case diagnostic), got cases {fit.cases}; the "
            "joint Stage-2 fit uses write_shared_form_fit_parameters"
        )
    case = fit.cases[0]
    trapped = {
        c: r.escape_fraction
        for c, r in fit.best.per_case.items()
        if r.penalty_Aps != 0.0
    }
    if trapped:
        raise ValueError(
            "refusing to write a trapped fit: escape_fraction by case "
            f"{trapped} (penalty != 0); the binding<->drag pair must permit "
            "full escape"
        )
    coeff_keys = _REQUIRED_COEFF_KEYS[fit.form]
    if (
        fit.coeff_errs is None
        or fit.uncertainty_model is None
        or any(k not in fit.coeff_errs for k in coeff_keys)
    ):
        raise ValueError(
            f"uncertainty fields (coeff_errs for {coeff_keys}, "
            "uncertainty_model) must be filled before writing the bundle"
        )
    if fit.form == POWER_LAW and v_ref_Aps is None:
        raise ValueError(
            "power_law per-case bundle needs v_ref_Aps (the locked V_REF_APS "
            "pivot speed, §10.4.1) to stamp the pivot block"
        )
    e_values = sorted(set(fit.best.e_bind_eV.values()))
    assert len(e_values) == 1  # single case -> trivially one binding
    e_shared = e_values[0]

    repo_root = REFERENCE_DRAG_ROOT.parents[2]
    r = fit.best.per_case[case]
    payload = {
        "extraction_method": "trajectory_matching",
        "form": fit.form,
        **{k: fit.best.coefficients[k] for k in coeff_keys},
        **{f"{k}_err": fit.coeff_errs[k] for k in coeff_keys},
        "uncertainty_model": fit.uncertainty_model,
        "meff_amu": fit.meff_amu,
        "extraction_mass_model": "constant",
        "effective_binding_energy_I_ion_eV": e_shared,
        "effective_binding_energy_err_eV": fit.e_bind_err_eV,
        "t_start": fit.windows[case][0],
        "t_end": fit.windows[case][1],
        "reference_file": _reference_file_str(
            fit.reference_paths[case], repo_root
        ),
        # --- per-case form provenance (§10.8) ---
        "calibration_case": case,
        "stage": "per_case_diagnostic",
        "variant": fit.variant,
        "anchors": {**fit.anchors, "provenance": anchor_provenance},
        "objective": "in_window_rmse_mean_abs_v2_plus_escape_penalty",
        "objective_rmse_Aps": r.rmse_Aps,
        "escape_fraction": r.escape_fraction,
        "num_overlap_points": r.num_overlap_points,
        "n_molecules_fit": fit.n_molecules,
        "seed": fit.seed,
        "optimizer": {
            "method": "nelder-mead",
            "n_starts": len(fit.starts),
            "evals": fit.n_evaluations,
            "converged": fit.converged,
            "start_minima": [
                {
                    "coefficients": s.coefficients,
                    "e_bind_eV": s.e_bind_eV,
                    "objective_Aps": s.objective,
                }
                for s in fit.starts
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
            # --- §10.8 single-case honesty flags ---
            # Diagnostic only: a single-case fit trivially matches its own
            # curve (cross-case axis spent for these forms); NEVER preset-wired.
            "per_case_calibrated_not_validated": True,
            "cross_case_axis_not_applied": True,
        },
    }
    if fit.form == POWER_LAW:
        C, n = fit.best.coefficients["C"], fit.best.coefficients["n"]
        payload["pivot"] = {
            "v_ref_Aps": v_ref_Aps,
            "gamma_ref": _C_to_pivot(C, n, v_ref_Aps),
            "note": (
                "optimizer-internal parameterization gamma_ref = "
                "C*v_ref**(n-1) (§10.4.1); the stamped law is the raw "
                "closed-form {C, n}. C_err is the fixed-n partial band "
                "C_err = gamma_ref_err / v_ref**(n-1), NOT a marginal "
                "uncertainty (§10.8)"
            ),
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fit_parameters.json"
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
    return out_path
