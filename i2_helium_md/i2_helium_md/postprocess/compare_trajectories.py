"""Numerical comparison of MD ion trajectories vs HeDFT reference.

Ports the trajectory comparison block in
``legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
simulation_image_only_trajectories.m`` (lines 97-118):

    dx = data_ion.x_ci(1:Nmol,:) - data_ion.x_ci(1+Nmol:end,:)
    dy = data_ion.y_ci(1:Nmol,:) - data_ion.y_ci(1+Nmol:end,:)
    dz = data_ion.z_ci(1:Nmol,:) - data_ion.z_ci(1+Nmol:end,:)
    dR_mean = mean(sqrt(dx.^2+dy.^2+dz.^2), 1)
    tmax = min(max(tR), max(data_ion.time_i))
    ...
    dR_mean_on_tR = interp1(t_md, d_md, tR_use, 'linear')
    rmse  = sqrt(mean((dR_mean_on_tR - R_use).^2))
    ratio = mean(dR_mean_on_tR ./ R_use)

The same overlap-and-resample contract is reused for the I1/I2 velocity
magnitude comparison, which the MATLAB script plots but does not summarise
numerically.

Out of scope (per CLAUDE.md): plotting, Abel inversion, mass filtering,
Bayesian histograms, pump-probe, CLI wrappers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..simulation.checkpoint import IonCheckpoint
from .hedft_loader import HedftTrajectory


_QUANTITY_DISTANCE = "distance_A"
_QUANTITY_V1 = "v1_magnitude_Aps"
_QUANTITY_V2 = "v2_magnitude_Aps"


@dataclass(frozen=True)
class TrajectoryComparison:
    """Result of comparing one MD-derived series to a HeDFT reference series.

    Attributes
    ----------
    quantity
        Identifier of the compared series (``"distance_A"``,
        ``"v1_magnitude_Aps"``, or ``"v2_magnitude_Aps"``).
    overlap_t_min_ps, overlap_t_max_ps
        Endpoints of the time interval on which both MD and HeDFT have
        samples, in picoseconds.
    num_overlap_points
        Number of HeDFT samples used for the comparison (i.e. the length
        of ``t_overlap_ps``). Non-finite samples after interpolation are
        still included in this count; the RMSE / mean ratio are computed
        on the finite subset.
    rmse
        Root-mean-square error in the units of the compared quantity:
        angstrom for distance, angstrom/ps for velocity magnitude.
    mean_ratio
        Mean of MD/HeDFT over samples where both are finite *and* the
        HeDFT denominator is non-zero. Dimensionless. ``nan`` if every
        finite reference sample is zero.
    t_overlap_ps
        HeDFT time samples that lie inside the overlap.
    md_on_hedft_grid
        MD series linearly interpolated onto ``t_overlap_ps``.
    hedft_on_overlap
        HeDFT reference values on ``t_overlap_ps``.
    """

    quantity: str
    overlap_t_min_ps: float
    overlap_t_max_ps: float
    num_overlap_points: int
    rmse: float
    mean_ratio: float
    t_overlap_ps: np.ndarray
    md_on_hedft_grid: np.ndarray
    hedft_on_overlap: np.ndarray


def compare_distance(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
) -> TrajectoryComparison:
    """Compare mean MD I-I separation against the HeDFT ``R_distance`` curve.

    Mirrors ``simulation_image_only_trajectories.m:97-118``.

    Parameters
    ----------
    ion
        Ion-stage checkpoint produced by the MD pipeline.
    hedft
        HeDFT reference trajectory loaded by
        :func:`load_hedft_trajectory`.

    Returns
    -------
    TrajectoryComparison
        With ``quantity == "distance_A"`` and RMSE in angstrom.
    """
    n = ion.num_molecules
    dx = ion.positions_x[:n] - ion.positions_x[n:]
    dy = ion.positions_y[:n] - ion.positions_y[n:]
    dz = ion.positions_z[:n] - ion.positions_z[n:]
    distance_md = np.mean(np.sqrt(dx * dx + dy * dy + dz * dz), axis=0)

    return _compare_series(
        quantity=_QUANTITY_DISTANCE,
        t_md=np.asarray(ion.time_ps, dtype=float),
        y_md=np.asarray(distance_md, dtype=float),
        t_ref=np.asarray(hedft.time_ps, dtype=float),
        y_ref=np.asarray(hedft.distance_A, dtype=float),
    )


def compare_velocity_magnitude(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    *,
    atom: Literal["I1", "I2"],
) -> TrajectoryComparison:
    """Compare the mean MD speed of I1 (or I2) against the matching HeDFT |v|.

    The MATLAB script overlays per-particle MD velocities with the HeDFT
    velocity magnitudes (``simulation_image_only_trajectories.m:180-202``);
    this function exposes the analogous numerical comparison.

    Parameters
    ----------
    ion
        Ion-stage checkpoint produced by the MD pipeline.
    hedft
        HeDFT reference trajectory loaded by
        :func:`load_hedft_trajectory`.
    atom
        ``"I1"`` selects atoms ``[0, num_molecules)`` and the
        ``v1_magnitude_Aps`` reference column. ``"I2"`` selects atoms
        ``[num_molecules, 2 * num_molecules)`` and ``v2_magnitude_Aps``.

    Returns
    -------
    TrajectoryComparison
        With ``quantity == "v1_magnitude_Aps"`` or
        ``"v2_magnitude_Aps"``; RMSE in angstrom/ps.
    """
    n = ion.num_molecules
    if atom == "I1":
        slc = slice(0, n)
        ref_series = hedft.v1_magnitude_Aps
        quantity = _QUANTITY_V1
    elif atom == "I2":
        slc = slice(n, 2 * n)
        ref_series = hedft.v2_magnitude_Aps
        quantity = _QUANTITY_V2
    else:
        raise ValueError(
            f"atom must be 'I1' or 'I2', got {atom!r}"
        )

    speed_per_atom = np.sqrt(
        ion.velocities_x[slc] ** 2
        + ion.velocities_y[slc] ** 2
        + ion.velocities_z[slc] ** 2
    )
    speed_md = np.mean(speed_per_atom, axis=0)

    return _compare_series(
        quantity=quantity,
        t_md=np.asarray(ion.time_ps, dtype=float),
        y_md=np.asarray(speed_md, dtype=float),
        t_ref=np.asarray(hedft.time_ps, dtype=float),
        y_ref=np.asarray(ref_series, dtype=float),
    )


def _compare_series(
    *,
    quantity: str,
    t_md: np.ndarray,
    y_md: np.ndarray,
    t_ref: np.ndarray,
    y_ref: np.ndarray,
) -> TrajectoryComparison:
    """Resample ``y_md`` onto the HeDFT grid and compute RMSE + mean ratio.

    Parameters mirror the MATLAB block:

        tmax = min(max(tR), max(t_md))
        mask_md = t_md <= tmax;  maskR = tR <= tmax
        dR_mean_on_tR = interp1(t_md, d_md, tR_use, 'linear')
        good = isfinite(dR_mean_on_tR) & isfinite(R_use)
        rmse  = sqrt(mean((dR_mean_on_tR(good) - R_use(good)).^2))
        ratio = mean(dR_mean_on_tR(good) ./ R_use(good))

    Lower bound is the symmetric ``max(min(tR), min(t_md))`` so that the
    interp call never has to extrapolate.
    """
    t_min = float(max(t_md[0], t_ref[0]))
    t_max = float(min(t_md[-1], t_ref[-1]))

    mask = (t_ref >= t_min) & (t_ref <= t_max)
    t_overlap = t_ref[mask]
    if t_overlap.size < 2:
        raise ValueError(
            "ion and hedft time axes do not overlap on at least 2 samples: "
            f"ion=[{t_md[0]}, {t_md[-1]}] ps, "
            f"hedft=[{t_ref[0]}, {t_ref[-1]}] ps, "
            f"overlap_points={int(t_overlap.size)}"
        )

    md_on_grid = np.interp(t_overlap, t_md, y_md)
    ref_on_overlap = y_ref[mask]

    finite = np.isfinite(md_on_grid) & np.isfinite(ref_on_overlap)
    if not np.any(finite):
        raise ValueError(
            f"No finite samples in the overlap window for quantity={quantity!r}"
        )

    diff = md_on_grid[finite] - ref_on_overlap[finite]
    rmse = float(np.sqrt(np.mean(diff * diff)))

    # Mean ratio is defined only where the HeDFT denominator is non-zero;
    # for velocities the reference starts at 0 so a few early samples must
    # be dropped. RMSE keeps the full finite mask above.
    ratio_mask = finite & (ref_on_overlap != 0.0)
    if np.any(ratio_mask):
        mean_ratio = float(
            np.mean(md_on_grid[ratio_mask] / ref_on_overlap[ratio_mask])
        )
    else:
        mean_ratio = float("nan")

    return TrajectoryComparison(
        quantity=quantity,
        overlap_t_min_ps=t_min,
        overlap_t_max_ps=t_max,
        num_overlap_points=int(t_overlap.size),
        rmse=rmse,
        mean_ratio=mean_ratio,
        t_overlap_ps=t_overlap,
        md_on_hedft_grid=md_on_grid,
        hedft_on_overlap=ref_on_overlap,
    )
