"""Literal helpers for ``post_process_single_pulse_paper_v4.m``.

The active non-effusive v4 branch is a one-panel radial VMI comparison plus
a separate simulated angular pair-covariance plot. This module keeps that
recipe explicit: line-26-to-34 experimental radial references, projected
detector velocity, fixed velocity bins, v4 smoothing, and the no-diagonal-
removal covariance behavior from the MATLAB script.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint
from ._smoothing import moving_mean


PAPER_V4_VELOCITY_BIN_WIDTH_APS = 0.05
PAPER_V4_VELOCITY_MAX_APS = 35.0
PAPER_V4_VELOCITY_SMOOTHING_WINDOW = 40
PAPER_V4_COVARIANCE_BINS = 90


@dataclass(frozen=True)
class PaperV4RadialReference:
    """One v4 experimental radial reference exported from MATLAB."""

    velocity_mps: np.ndarray
    signal_arb: np.ndarray
    label: str
    source_path: Path


@dataclass(frozen=True)
class PaperV4VelocityCurve:
    """MATLAB-v4 simulated projected-velocity curve for one mass selection."""

    mass_amu: float
    bin_centers_Aps: np.ndarray
    bin_centers_mps: np.ndarray
    bin_edges_Aps: np.ndarray
    counts: np.ndarray
    smoothed: np.ndarray
    shifted: np.ndarray
    normalised: np.ndarray
    num_atoms_used: int
    smoothing_window: int


@dataclass(frozen=True)
class PaperV4AngularCovariance:
    """Literal v4 angular pair covariance and scatter-overlay points."""

    counts: np.ndarray
    theta_centers_rad: np.ndarray
    theta_edges_rad: np.ndarray
    theta_pairs_rad: np.ndarray
    num_pairs_used: int


def load_paper_v4_radial_references(directory: str | Path, *, only_include_iplus_he: bool = False) -> list[PaperV4RadialReference]:
    """Load all ``*_radial.csv`` files from a v4 reference directory."""

    root = Path(directory)
    if not root.exists():
        return []
    refs = []
    for path in sorted(root.glob("*_radial.csv")):
        if only_include_iplus_he and not path.stem.startswith("iplus_he"):
            continue
        refs.append(load_paper_v4_radial_reference(path))
    return refs


def load_paper_v4_radial_reference(path: str | Path) -> PaperV4RadialReference:
    """Load a v4 radial reference CSV with ``v_mps,signal_arb`` columns."""

    p = Path(path)
    data = _load_named_csv(p)
    names = data.dtype.names or ()
    if not all(name in names for name in ("v_mps", "signal_arb")):
        raise ValueError(f"{p.name} must contain columns v_mps and signal_arb")
    return PaperV4RadialReference(
        velocity_mps=np.asarray(data["v_mps"], dtype=float),
        signal_arb=np.asarray(data["signal_arb"], dtype=float),
        label=_label_from_reference_name(p),
        source_path=p.resolve(),
    )


def paper_v4_velocity_curve(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
) -> PaperV4VelocityCurve:
    """Compute the v4 projected-velocity histogram for one mass selection."""

    selected = _paper_v4_atom_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x, dtype=float)[selected]
    vy = np.asarray(ion.velocities_final_y, dtype=float)[selected]
    projected_speed = np.sqrt(vx * vx + vy * vy)

    edges = _velocity_edges()
    counts, edges = np.histogram(projected_speed, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = moving_mean(counts, PAPER_V4_VELOCITY_SMOOTHING_WINDOW)
    shifted = smoothed - float(np.nanmin(smoothed))
    return PaperV4VelocityCurve(
        mass_amu=float(mass_amu),
        bin_centers_Aps=centers,
        bin_centers_mps=centers * 100.0,
        bin_edges_Aps=edges,
        counts=counts.astype(int),
        smoothed=smoothed,
        shifted=shifted,
        normalised=max_normalise(shifted),
        num_atoms_used=int(np.count_nonzero(selected)),
        smoothing_window=PAPER_V4_VELOCITY_SMOOTHING_WINDOW,
    )


def paper_v4_angular_pair_covariance(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
    n_theta_bins: int = PAPER_V4_COVARIANCE_BINS,
) -> PaperV4AngularCovariance:
    """Compute the literal simulated angular covariance block from v4."""

    if n_theta_bins < 1:
        raise ValueError(f"n_theta_bins must be >= 1, got {n_theta_bins}")

    n = ion.num_molecules
    atom_ok = _paper_v4_atom_selection(ion, mass_amu=mass_amu, raise_empty=False)
    pair_ok = atom_ok[:n] & atom_ok[n:]

    vx_a = np.asarray(ion.velocities_final_x[:n], dtype=float)[pair_ok]
    vy_a = np.asarray(ion.velocities_final_y[:n], dtype=float)[pair_ok]
    vx_b = np.asarray(ion.velocities_final_x[n:], dtype=float)[pair_ok]
    vy_b = np.asarray(ion.velocities_final_y[n:], dtype=float)[pair_ok]

    theta_a = np.mod(np.arctan2(vx_a, vy_a) + np.pi, 2.0 * np.pi)
    theta_b = np.mod(np.arctan2(vx_b, vy_b) + np.pi, 2.0 * np.pi)
    edges = np.linspace(0.0, 2.0 * np.pi, n_theta_bins + 1)
    counts, _, _ = np.histogram2d(theta_a, theta_b, bins=(edges, edges))
    theta_pairs = np.column_stack([theta_a, theta_b]) if theta_a.size else np.empty((0, 2))
    return PaperV4AngularCovariance(
        counts=counts.astype(float),
        theta_centers_rad=0.5 * (edges[:-1] + edges[1:]),
        theta_edges_rad=edges,
        theta_pairs_rad=theta_pairs,
        num_pairs_used=int(np.count_nonzero(pair_ok)),
    )


def max_normalise(values) -> np.ndarray:
    """Normalise by maximum only, returning zeros for non-positive traces."""

    data = np.asarray(values, dtype=float)
    if data.size == 0:
        return data.copy()
    scale = float(np.nanmax(data))
    if not np.isfinite(scale) or scale <= 0.0:
        return np.zeros_like(data, dtype=float)
    return data / scale


def _paper_v4_atom_selection(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
    raise_empty: bool = True,
) -> np.ndarray:
    masses_amu = np.round(np.asarray(ion.mass_final_kg, dtype=float) / U_KG)
    mass_mask = masses_amu == round(float(mass_amu))
    outside = np.concatenate([ion.b_ion_outside, ion.b_ion_outside]).astype(bool)
    selected = mass_mask & outside
    if raise_empty and not np.any(selected):
        raise ValueError(
            f"No atoms passed the paper-v4 mass={mass_amu:.0f} amu "
            "and b_ion_outside filter."
        )
    return selected


def _velocity_edges() -> np.ndarray:
    return np.arange(
        0.0,
        PAPER_V4_VELOCITY_MAX_APS + PAPER_V4_VELOCITY_BIN_WIDTH_APS,
        PAPER_V4_VELOCITY_BIN_WIDTH_APS,
    )


def _load_named_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"paper-v4 reference file not found: {path.resolve()}")
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    if data.dtype.names is None:
        raise ValueError(f"{path.name} must have a header row")
    return np.atleast_1d(data)


def _label_from_reference_name(path: Path) -> str:
    match = re.fullmatch(
        r"iplus_(gas|drop|he)_([0-9]+)mw_([0-9]+)_radial",
        path.stem,
    )
    if match is None:
        return path.stem
    channel, power, measurement = match.groups()
    channel_label = {
        "gas": "I+ gas",
        "drop": "I+ drop",
        "he": "I+He",
    }[channel]
    return f"{channel_label} {power} mW ({measurement})"
