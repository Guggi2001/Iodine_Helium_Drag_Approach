"""Force-balance drag cross-check on the TDDFT traces (Method A, G4.10).

THESIS_FINAL_GOALS.md G4.10. Method B (``trajectory_matching``) calibrates the
drag on one radial curve (18 A). This module is the independent local check:
at every time point of a TDDFT trace the drag is what remains of the force
balance along the direction of motion of the compared atom (atom 2),

.. math::

    \\gamma(v)\\,\\hat\\rho(r)\\,v
        = F_C\\cos\\theta - \\frac{dU}{dr}\\,(\\hat r\\cdot\\hat v)
          - m\\,\\frac{d|v|}{dt}
    \\qquad [\\text{amu·Å/ps}^2]

with ``theta`` the angle between v2 and r2 - r1, ``r`` the distance of atom 2
from the droplet centre (origin, the I2 centre of mass), ``dU/dr`` the package
solvation-well derivative (positive; the force -dU/dr r_hat points inward) and
``rho_hat`` the erf-complement density gate. Only ``F . v_hat`` changes the
speed, so projecting every force on ``v_hat`` is exact for ``d|v|/dt``; the
legacy scalar form (``cos theta = 1``, no well) books ``F_C (1 - cos theta)``
as drag whenever the motion is not along the pair axis.

Units follow the friction convention (CLAUDE.md): forces in amu·Å/ps², ``a``
[amu/ps], ``b`` [amu·ps/Å²], ``F = a v + b v^3``; mass enters only the
inertial term. No MD is run.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import savgol_filter

from ..physics.constants import EV, EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2, U, coulomb_energy
from ..physics.drag import spatial_gate
from ..physics.potentials import droplet_force
from ..postprocess.hedft_loader import HedftTrajectory

# (amu·Å/ps²) per (eV/Å): F[eV/A] * EV*1e-4 / U  (the leapfrog.py conversion
# with the mass expressed in amu instead of kg). ~= 9648.533.
AMU_A_PER_PS2_PER_EV_PER_A: float = EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 / U

PROJECTIONS: tuple[str, ...] = ("projected", "scalar", "legacy_no_well")
FORMS: tuple[str, ...] = ("linear", "cubic", "linear_cubic")


def _uniform_dt_ps(time_ps: np.ndarray, rtol: float = 1e-3) -> float:
    """Return the sampling step [ps]; raise on non-uniform or too-short axes."""
    t = np.asarray(time_ps, dtype=float)
    if t.ndim != 1 or t.size < 3:
        raise ValueError(f"time axis must be 1-D with >= 3 samples, got shape {t.shape}")
    dt = np.diff(t)
    dt_mean = float(np.mean(dt))
    if dt_mean <= 0.0 or not np.allclose(dt, dt_mean, rtol=rtol, atol=1e-12):
        raise ValueError("time axis must be strictly increasing and uniformly sampled")
    return dt_mean


def _odd_window_points(window_ps: float, dt_ps: float, polyorder: int) -> int:
    """SG window length in samples: nearest odd count to ``window_ps / dt_ps``."""
    if window_ps <= 0.0:
        raise ValueError(f"window_ps must be > 0, got {window_ps}")
    n = int(round(window_ps / dt_ps))
    n = n + 1 if n % 2 == 0 else n
    if n <= polyorder:
        raise ValueError(
            f"window {window_ps} ps = {n} samples is too short for polyorder {polyorder}"
        )
    return n


@dataclass(frozen=True)
class BalanceGeometry:
    """Pair geometry of a TDDFT trace on its native time axis.

    Attributes
    ----------
    time_ps : np.ndarray, shape (T,)
    pair_distance_A : np.ndarray, shape (T,)
        ``|r2 - r1|`` [A] from the raw positions.
    radius_A : np.ndarray, shape (T,)
        ``|r2|`` [A], distance of atom 2 from the droplet centre (origin).
    cos_theta : np.ndarray, shape (T,)
        ``v_hat2 . (r2 - r1)/|r2 - r1|``; NaN where the smoothed speed is 0.
    radial_cos : np.ndarray, shape (T,)
        ``v_hat2 . r_hat2``; NaN where the smoothed speed is 0.
    """

    time_ps: np.ndarray
    pair_distance_A: np.ndarray
    radius_A: np.ndarray
    cos_theta: np.ndarray
    radial_cos: np.ndarray


def balance_geometry(
    traj: HedftTrajectory, *, direction_window_ps: float = 0.5
) -> BalanceGeometry:
    """Build the projection geometry of atom 2 from a HeDFT trajectory.

    The velocity *direction* comes from the atom-2 components, each smoothed
    with an order-1 Savitzky-Golay filter over ``direction_window_ps`` (G4.10:
    0.5 ps); distances come from the raw positions. The droplet centre is the
    origin (the I2 centre of mass at t = 0; G4.10 assumption).

    Parameters
    ----------
    traj : HedftTrajectory
        Loaded TDDFT trace (uniform time axis required).
    direction_window_ps : float
        SG window for the velocity components [ps], > 0.

    Returns
    -------
    BalanceGeometry

    Raises
    ------
    ValueError
        On a non-uniform time axis or a window too short for the filter.
    """
    dt = _uniform_dt_ps(traj.time_ps)
    w = _odd_window_points(direction_window_ps, dt, polyorder=1)
    v = np.stack(
        [
            savgol_filter(np.asarray(c, dtype=float), w, 1, mode="interp")
            for c in (traj.v2_x_Aps, traj.v2_y_Aps, traj.v2_z_Aps)
        ],
        axis=1,
    )
    p1 = np.stack([traj.x1_A, traj.y1_A, traj.z1_A], axis=1).astype(float)
    p2 = np.stack([traj.x2_A, traj.y2_A, traj.z2_A], axis=1).astype(float)
    pair = p2 - p1
    pair_distance = np.linalg.norm(pair, axis=1)
    radius = np.linalg.norm(p2, axis=1)
    speed = np.linalg.norm(v, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        v_hat = np.where(speed[:, None] > 0.0, v / speed[:, None], np.nan)
        cos_theta = np.sum(v_hat * pair / pair_distance[:, None], axis=1)
        radial_cos = np.sum(v_hat * p2 / radius[:, None], axis=1)
    return BalanceGeometry(
        time_ps=np.asarray(traj.time_ps, dtype=float),
        pair_distance_A=pair_distance,
        radius_A=radius,
        cos_theta=cos_theta,
        radial_cos=radial_cos,
    )


def savgol_speed_and_derivative(
    time_ps: np.ndarray,
    speed_Aps: np.ndarray,
    window_ps: float,
    *,
    smooth_speed: bool,
    polyorder: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Speed and its time derivative by Savitzky-Golay.

    Parameters
    ----------
    time_ps, speed_Aps : np.ndarray, shape (T,)
        Uniform time axis [ps] and speed [A/ps].
    window_ps : float
        SG window [ps] (nearest odd sample count).
    smooth_speed : bool
        ``False`` for an already smoothed speed (the committed reference): the
        speed is returned unchanged, only the derivative is filtered.
        ``True`` for a raw speed: the same filter also returns the smoothed
        speed (deriv 0).
    polyorder : int
        SG polynomial order (2: local quadratic, so the derivative is not a
        window-mean slope).

    Returns
    -------
    speed_Aps, dvdt_Aps2 : np.ndarray, shape (T,)
        [A/ps], [A/ps^2].

    Raises
    ------
    ValueError
        On shape mismatch, non-uniform time or a window shorter than the order.
    """
    t = np.asarray(time_ps, dtype=float)
    v = np.asarray(speed_Aps, dtype=float)
    if v.shape != t.shape:
        raise ValueError(f"speed shape {v.shape} != time shape {t.shape}")
    dt = _uniform_dt_ps(t)
    w = _odd_window_points(window_ps, dt, polyorder)
    if w > t.size:
        raise ValueError(f"window {w} samples exceeds the trace length {t.size}")
    dvdt = savgol_filter(v, w, polyorder, deriv=1, delta=dt, mode="interp")
    v_out = savgol_filter(v, w, polyorder, mode="interp") if smooth_speed else v.copy()
    return v_out, dvdt


@dataclass(frozen=True)
class ForceBalanceTrace:
    """Per-time-point force balance along v_hat (all forces amu·Å/ps²).

    ``F_drag_apparent = F_coulomb_along - F_well_along - F_inertia`` is the
    gate-included drag ``gamma*rho_hat*v``; ``F_drag_bulk = F_drag_apparent /
    rho_hat`` on interior points (``rho_hat >= interior_rho_min``), NaN
    elsewhere. ``F_well_along = dU/dr * (r_hat . v_hat)`` (positive when moving
    outward, i.e. the well decelerates).
    """

    time_ps: np.ndarray
    speed_Aps: np.ndarray
    dvdt_Aps2: np.ndarray
    mass_amu: np.ndarray
    depth_A: np.ndarray
    rho_hat: np.ndarray
    cos_theta: np.ndarray
    F_coulomb: np.ndarray
    F_coulomb_along: np.ndarray
    F_well_along: np.ndarray
    F_inertia: np.ndarray
    F_drag_apparent: np.ndarray
    interior: np.ndarray
    F_drag_bulk: np.ndarray


def force_balance_trace(
    geometry: BalanceGeometry,
    time_ps: np.ndarray,
    speed_Aps: np.ndarray,
    dvdt_Aps2: np.ndarray,
    *,
    mass_amu: float | np.ndarray,
    binding_energy_eV: float,
    potential_steepness_A: float,
    droplet_radius_A: float,
    gate_steepness_A: float,
    projection: str = "projected",
    interior_rho_min: float = 0.9,
) -> ForceBalanceTrace:
    """Evaluate the force balance on the speed time axis.

    Parameters
    ----------
    geometry : BalanceGeometry
        From :func:`balance_geometry`; linearly interpolated onto ``time_ps``.
    time_ps, speed_Aps, dvdt_Aps2 : np.ndarray, shape (T,)
        Evaluation axis [ps], speed [A/ps], its derivative [A/ps^2].
    mass_amu : float or np.ndarray, shape (T,)
        Inertial mass [amu] (``m_eff`` or an anchored ``m(t)``).
    binding_energy_eV, potential_steepness_A : float
        Solvation well depth [eV] and width [A] (``droplet_force``).
    droplet_radius_A : float
        Droplet radius [A]; depth = ``r - droplet_radius_A``.
    gate_steepness_A : float
        Density-gate width ``s_rho`` [A].
    projection : {"projected", "scalar", "legacy_no_well"}
        ``projected``: the G4.10 balance. ``scalar``: ``cos theta -> 1``,
        well kept. ``legacy_no_well``: the old Method A, ``F_C - m dv/dt``.
    interior_rho_min : float
        Interior cut on ``rho_hat`` (G4.10: 0.9).

    Returns
    -------
    ForceBalanceTrace

    Raises
    ------
    ValueError
        On an unknown projection, shape mismatches, or an evaluation axis
        outside the geometry's time range.
    """
    if projection not in PROJECTIONS:
        raise ValueError(f"projection must be one of {PROJECTIONS}, got {projection!r}")
    t = np.asarray(time_ps, dtype=float)
    v = np.asarray(speed_Aps, dtype=float)
    a = np.asarray(dvdt_Aps2, dtype=float)
    if not (t.shape == v.shape == a.shape) or t.ndim != 1:
        raise ValueError("time, speed and dvdt must be 1-D arrays of equal shape")
    m = np.broadcast_to(np.asarray(mass_amu, dtype=float), t.shape).copy()
    g = geometry
    if t[0] < g.time_ps[0] or t[-1] > g.time_ps[-1]:
        raise ValueError(
            f"evaluation axis [{t[0]}, {t[-1]}] ps lies outside the geometry "
            f"[{g.time_ps[0]}, {g.time_ps[-1]}] ps"
        )

    pair = np.interp(t, g.time_ps, g.pair_distance_A)
    radius = np.interp(t, g.time_ps, g.radius_A)
    cos_theta = np.interp(t, g.time_ps, g.cos_theta)
    radial_cos = np.interp(t, g.time_ps, g.radial_cos)
    depth = radius - droplet_radius_A

    # F_C = -dU_C/dR = U_C/R for the pure 1/R Coulomb pair [eV/A].
    F_coulomb = coulomb_energy(pair) / EV / pair * AMU_A_PER_PS2_PER_EV_PER_A
    dU_dr = droplet_force(depth, potential_steepness_A, binding_energy_eV)
    F_well_along = dU_dr * radial_cos * AMU_A_PER_PS2_PER_EV_PER_A
    if projection == "projected":
        F_coulomb_along = F_coulomb * cos_theta
    else:
        F_coulomb_along = F_coulomb.copy()
    if projection == "legacy_no_well":
        F_well_along = np.zeros_like(F_well_along)

    F_inertia = m * a
    F_app = F_coulomb_along - F_well_along - F_inertia
    rho_hat = spatial_gate(depth, gate_steepness_A)
    interior = rho_hat >= interior_rho_min
    F_bulk = np.where(interior, F_app / np.where(interior, rho_hat, 1.0), np.nan)
    return ForceBalanceTrace(
        time_ps=t,
        speed_Aps=v,
        dvdt_Aps2=a,
        mass_amu=m,
        depth_A=depth,
        rho_hat=rho_hat,
        cos_theta=cos_theta,
        F_coulomb=F_coulomb,
        F_coulomb_along=F_coulomb_along,
        F_well_along=F_well_along,
        F_inertia=F_inertia,
        F_drag_apparent=F_app,
        interior=interior,
        F_drag_bulk=F_bulk,
    )


@dataclass(frozen=True)
class DragFormFit:
    """Least-squares fit of ``F = a v + b v^3`` (through the origin).

    ``a`` [amu/ps], ``b`` [amu·ps/A^2] (0 where the form pins it);
    ``a_se``/``b_se`` are naive OLS standard errors -- the residuals are
    strongly autocorrelated (dense smoothed samples), so they understate the
    real uncertainty and are recorded only for scale. ``rmse_force`` in
    amu·A/ps^2; ``v_min``/``v_max`` [A/ps] give the speed lever arm.
    """

    form: str
    a: float
    b: float
    a_se: float
    b_se: float
    rmse_force: float
    n_points: int
    v_min: float
    v_max: float


def fit_drag_form(speed_Aps: np.ndarray, force: np.ndarray, form: str) -> DragFormFit:
    """Fit one drag form to (speed, force) samples.

    Parameters
    ----------
    speed_Aps, force : np.ndarray, shape (K,)
        Interior speeds [A/ps] and bulk drag forces [amu·A/ps^2]; non-finite
        pairs are dropped.
    form : {"linear", "cubic", "linear_cubic"}

    Returns
    -------
    DragFormFit

    Raises
    ------
    ValueError
        On an unknown form, shape mismatch, or fewer finite samples than
        parameters + 1.
    """
    if form not in FORMS:
        raise ValueError(f"form must be one of {FORMS}, got {form!r}")
    v = np.asarray(speed_Aps, dtype=float)
    F = np.asarray(force, dtype=float)
    if v.shape != F.shape or v.ndim != 1:
        raise ValueError("speed and force must be 1-D arrays of equal shape")
    ok = np.isfinite(v) & np.isfinite(F)
    v, F = v[ok], F[ok]
    columns = {"linear": [v], "cubic": [v**3], "linear_cubic": [v, v**3]}[form]
    X = np.stack(columns, axis=1)
    k = X.shape[1]
    if v.size < k + 1:
        raise ValueError(f"need >= {k + 1} finite samples for {form}, got {v.size}")
    coef, *_ = np.linalg.lstsq(X, F, rcond=None)
    resid = F - X @ coef
    sigma2 = float(resid @ resid) / (v.size - k)
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    if form == "linear":
        a, b, a_se, b_se = coef[0], 0.0, se[0], 0.0
    elif form == "cubic":
        a, b, a_se, b_se = 0.0, coef[0], 0.0, se[0]
    else:
        a, b, a_se, b_se = coef[0], coef[1], se[0], se[1]
    return DragFormFit(
        form=form,
        a=float(a),
        b=float(b),
        a_se=float(a_se),
        b_se=float(b_se),
        rmse_force=float(np.sqrt(np.mean(resid**2))),
        n_points=int(v.size),
        v_min=float(v.min()),
        v_max=float(v.max()),
    )
