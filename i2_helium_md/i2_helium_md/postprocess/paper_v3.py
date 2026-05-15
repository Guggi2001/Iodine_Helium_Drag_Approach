"""Literal helpers for ``post_process_single_pulse_paper_v3.m``.

This module captures the active non-effusive droplet branch of the legacy
MATLAB paper-figure script. It intentionally keeps MATLAB's plotting recipe
visible: rounded mass selection, outside-droplet filtering, detector-projected
velocity, fixed bin edges, moving means, and max-only normalisation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint
from ._smoothing import moving_mean


PAPER_V3_VELOCITY_BIN_WIDTH_APS = 0.05
PAPER_V3_VELOCITY_MAX_APS = 35.0
PAPER_V3_VELOCITY_SMOOTHING_WINDOW = 15
PAPER_V3_PHI_BIN_WIDTH_RAD = 0.05
PAPER_V3_PHI_SMOOTHING_WINDOW = 15


@dataclass(frozen=True)
class PaperV3RadialReference:
    """Experimental radial reference exported from the legacy MATLAB recipe.

    ``velocity_mps`` is in m/s because v3 labels the top panel ``v / m/s``.
    ``signal`` is a 2-D array with one column per signal trace.
    """

    velocity_mps: np.ndarray
    signal: np.ndarray
    signal_labels: tuple[str, ...]
    source_path: Path


@dataclass(frozen=True)
class PaperV3PhiReference:
    """Experimental angular reference exported from the legacy MATLAB recipe."""

    phi_rad: np.ndarray
    signal_arb: np.ndarray
    source_path: Path


@dataclass(frozen=True)
class PaperV3VelocityCurve:
    """MATLAB-v3 simulated projected-velocity curve for one mass selection."""

    mass_amu: float
    bin_centers_Aps: np.ndarray
    bin_centers_mps: np.ndarray
    bin_edges_Aps: np.ndarray
    counts: np.ndarray
    smoothed: np.ndarray
    normalised: np.ndarray
    num_atoms_used: int
    smoothing_window: int


@dataclass(frozen=True)
class PaperV3PhiCurve:
    """MATLAB-v3 simulated phi curve for one mass selection."""

    mass_amu: float
    bin_centers_rad: np.ndarray
    bin_edges_rad: np.ndarray
    counts: np.ndarray
    smoothed: np.ndarray
    normalised: np.ndarray
    num_atoms_used: int
    smoothing_window: int


def matlab_max_normalise(values) -> np.ndarray:
    """Normalise by maximum only, matching ``y = y / max(y)`` in MATLAB."""

    data = np.asarray(values, dtype=float)
    if data.size == 0:
        return data.copy()
    scale = float(np.nanmax(data))
    if not np.isfinite(scale) or scale <= 0.0:
        return np.zeros_like(data, dtype=float)
    return data / scale


def load_paper_v3_radial_reference(path: str | Path) -> PaperV3RadialReference:
    """Load a v3 radial reference CSV.

    Required columns are ``v_mps`` followed by one or more signal columns.
    This supports both a single I+He radial trace and a wide timescan export.
    """

    p = Path(path)
    data = _load_named_csv(p)
    names = data.dtype.names or ()
    if "v_mps" not in names:
        raise ValueError(f"{p.name} must contain a v_mps column")
    signal_labels = tuple(name for name in names if name != "v_mps")
    if not signal_labels:
        raise ValueError(f"{p.name} must contain at least one signal column")
    signal = np.column_stack([np.asarray(data[name], dtype=float) for name in signal_labels])
    return PaperV3RadialReference(
        velocity_mps=np.asarray(data["v_mps"], dtype=float),
        signal=signal,
        signal_labels=signal_labels,
        source_path=p.resolve(),
    )


def load_paper_v3_phi_reference(path: str | Path) -> PaperV3PhiReference:
    """Load a v3 angular reference CSV with ``phi_rad,signal_arb`` columns."""

    p = Path(path)
    data = _load_named_csv(p)
    names = data.dtype.names or ()
    required = ("phi_rad", "signal_arb")
    if not all(name in names for name in required):
        raise ValueError(f"{p.name} must contain columns phi_rad and signal_arb")
    return PaperV3PhiReference(
        phi_rad=np.asarray(data["phi_rad"], dtype=float),
        signal_arb=np.asarray(data["signal_arb"], dtype=float),
        source_path=p.resolve(),
    )


def paper_v3_velocity_curve(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
) -> PaperV3VelocityCurve:
    """Compute the v3 projected-velocity histogram for one mass selection."""

    selected = _paper_v3_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x)[selected]
    vy = np.asarray(ion.velocities_final_y)[selected]
    projected_speed = np.sqrt(vx * vx + vy * vy)
    edges = _velocity_edges()
    counts, edges = np.histogram(projected_speed, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = moving_mean(counts, PAPER_V3_VELOCITY_SMOOTHING_WINDOW)
    return PaperV3VelocityCurve(
        mass_amu=float(mass_amu),
        bin_centers_Aps=centers,
        bin_centers_mps=centers * 100.0,
        bin_edges_Aps=edges,
        counts=counts,
        smoothed=smoothed,
        normalised=matlab_max_normalise(smoothed),
        num_atoms_used=int(np.count_nonzero(selected)),
        smoothing_window=PAPER_V3_VELOCITY_SMOOTHING_WINDOW,
    )


def paper_v3_phi_curve(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
) -> PaperV3PhiCurve:
    """Compute the v3 phi histogram for one mass selection."""

    selected = _paper_v3_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x)[selected]
    vy = np.asarray(ion.velocities_final_y)[selected]
    phi = np.mod(np.arctan2(vy, vx) + np.pi, 2.0 * np.pi)
    edges = _phi_edges()
    counts, edges = np.histogram(phi, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = moving_mean(counts, PAPER_V3_PHI_SMOOTHING_WINDOW)
    return PaperV3PhiCurve(
        mass_amu=float(mass_amu),
        bin_centers_rad=centers,
        bin_edges_rad=edges,
        counts=counts,
        smoothed=smoothed,
        normalised=matlab_max_normalise(smoothed),
        num_atoms_used=int(np.count_nonzero(selected)),
        smoothing_window=PAPER_V3_PHI_SMOOTHING_WINDOW,
    )


def _paper_v3_selection(ion: IonCheckpoint, *, mass_amu: float) -> np.ndarray:
    masses_amu = np.round(np.asarray(ion.mass_final_kg, dtype=float) / U_KG)
    mass_mask = masses_amu == round(float(mass_amu))
    outside = np.concatenate([ion.b_ion_outside, ion.b_ion_outside]).astype(bool)
    selected = mass_mask & outside
    if not np.any(selected):
        raise ValueError(
            f"No atoms passed the paper-v3 mass={mass_amu:.0f} amu "
            "and b_ion_outside filter."
        )
    return selected


def _velocity_edges() -> np.ndarray:
    return np.arange(
        0.0,
        PAPER_V3_VELOCITY_MAX_APS + PAPER_V3_VELOCITY_BIN_WIDTH_APS,
        PAPER_V3_VELOCITY_BIN_WIDTH_APS,
    )


def _phi_edges() -> np.ndarray:
    return np.arange(0.0, 2.0 * np.pi + PAPER_V3_PHI_BIN_WIDTH_RAD, PAPER_V3_PHI_BIN_WIDTH_RAD)


def _load_named_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"paper-v3 reference file not found: {path.resolve()}")
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    if data.dtype.names is None:
        raise ValueError(f"{path.name} must have a header row")
    return np.atleast_1d(data)
