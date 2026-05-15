"""Literal helpers for ``post_process_single_pulse_paper_IplusHe_comparison.m``.

The Python name ``paper_v2`` is a consistency alias for the legacy MATLAB
I+He comparison script. This module keeps the recipe explicit: MATLAB-exported
radial references, MATLAB-exported processed 2-D VMI images, and a simulated
detector-plane velocity map using the script's nearest-bin counting.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

import numpy as np
from scipy.io import loadmat

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint
from ._smoothing import moving_mean


PAPER_V2_IMAGE_BINS_APS = np.arange(-35.0, 35.0 + 0.2, 0.2)
PAPER_V2_VELOCITY_BIN_WIDTH_APS = 0.05
PAPER_V2_VELOCITY_MAX_APS = 35.0
PAPER_V2_VELOCITY_SMOOTHING_WINDOW = 15
PAPER_V2_PHI_BIN_WIDTH_RAD = 0.05
PAPER_V2_PHI_SMOOTHING_WINDOW = 15
PAPER_V2_RADIAL_REFERENCE_ORDER = (
    "iplus_gas_300mw_43562_radial.csv",
    "iplus_he_high_snr_radial.csv",
    "iplus_he_160mw_43556_radial.csv",
    "iplus_he_600mw_43569_radial.csv",
)


@dataclass(frozen=True)
class PaperV2RadialReference:
    """One paper-v2 experimental radial reference exported from MATLAB."""

    velocity_Aps: np.ndarray
    velocity_mps: np.ndarray
    signal_arb: np.ndarray
    label: str
    source_path: Path


@dataclass(frozen=True)
class PaperV2VMIImageReference:
    """Processed experimental 2-D VMI image exported from MATLAB."""

    vx_Aps: np.ndarray
    vy_Aps: np.ndarray
    vx_mps: np.ndarray
    vy_mps: np.ndarray
    intensity: np.ndarray
    metadata: dict[str, Any]
    source_path: Path


@dataclass(frozen=True)
class PaperV2VMIPolarImageReference:
    """Processed experimental 2-D polar VMI image (phi rows, v_radius cols).

    Mirrors the MATLAB ``res.image_polar`` matrix exported by
    ``data/reference/scripts/export_paper_v2_reference_data.m``.
    """

    phi_rad: np.ndarray
    v_radius_Aps: np.ndarray
    v_radius_mps: np.ndarray
    intensity: np.ndarray
    metadata: dict[str, Any]
    source_path: Path


@dataclass(frozen=True)
class PaperV2PhiReference:
    """One paper-v2 experimental phi reference exported from MATLAB."""

    phi_rad: np.ndarray
    signal_arb: np.ndarray
    source_path: Path


@dataclass(frozen=True)
class PaperV2VelocityMap:
    """MATLAB-v2 simulated detector-plane velocity map.

    ``counts`` intentionally follows the legacy matrix layout:
    ``counts[vx_index, vy_index]``. The paper-v2 plotting script transposes it
    for Matplotlib display so the final figure has physical ``v_x`` on the
    horizontal axis and ``v_y`` on the vertical axis.
    """

    velocity_bins_Aps: np.ndarray
    counts: np.ndarray
    mass_amu: float
    num_atoms_used: int


@dataclass(frozen=True)
class PaperV2VelocityCurve:
    """Paper-v2 simulated projected-velocity curve for one mass selection."""

    mass_amu: float
    bin_centers_Aps: np.ndarray
    bin_centers_mps: np.ndarray
    bin_edges_Aps: np.ndarray
    bin_edges_mps: np.ndarray
    counts: np.ndarray
    smoothed: np.ndarray
    normalised: np.ndarray
    num_atoms_used: int
    smoothing_window: int


@dataclass(frozen=True)
class PaperV2PhiCurve:
    """Paper-v2 simulated phi curve for one mass selection."""

    mass_amu: float
    bin_centers_rad: np.ndarray
    bin_edges_rad: np.ndarray
    counts: np.ndarray
    smoothed: np.ndarray
    normalised: np.ndarray
    num_atoms_used: int
    smoothing_window: int


def load_paper_v2_radial_references(directory: str | Path) -> list[PaperV2RadialReference]:
    """Load the two radial references actually plotted by paper-v2 MATLAB."""

    root = Path(directory)
    if not root.exists():
        return []
    return [
        load_paper_v2_radial_reference(root / filename)
        for filename in PAPER_V2_RADIAL_REFERENCE_ORDER
        if (root / filename).exists()
    ]


def load_paper_v2_radial_reference(path: str | Path) -> PaperV2RadialReference:
    """Load a paper-v2 radial reference CSV.

    Canonical columns are ``v_mps,signal_arb``. Legacy files written with
    ``v_Aps,signal_arb`` are accepted and converted by a factor of 100.
    """

    p = Path(path)
    data = _load_named_csv(p)
    names = data.dtype.names or ()
    if "signal_arb" not in names:
        raise ValueError(f"{p.name} must contain a signal_arb column")
    if "v_mps" in names:
        velocity_Aps = np.asarray(data["v_mps"], dtype=float) / 100.0
    elif "v_Aps" in names:
        velocity_Aps = np.asarray(data["v_Aps"], dtype=float)
    else:
        raise ValueError(f"{p.name} must contain a v_mps (or legacy v_Aps) column")
    return PaperV2RadialReference(
        velocity_Aps=velocity_Aps,
        velocity_mps=velocity_Aps * 100.0,
        signal_arb=np.asarray(data["signal_arb"], dtype=float),
        label=_label_from_reference_name(p),
        source_path=p.resolve(),
    )


def load_paper_v2_vmi_image_reference(path: str | Path) -> PaperV2VMIImageReference:
    """Load a paper-v2 processed VMI image reference from ``.mat`` or ``.npz``.

    Canonical axis fields are ``vx_mps,vy_mps,intensity``. Legacy files
    using ``vx_Aps,vy_Aps,intensity`` (A/ps) are accepted; the loader
    converts the axes to A/ps either way (m/s axes are divided by 100).
    """

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"paper-v2 VMI image reference not found: {p.resolve()}")

    def _pick(container, has_key) -> tuple[np.ndarray, np.ndarray]:
        if has_key("vx_mps") and has_key("vy_mps"):
            vx = np.atleast_1d(np.squeeze(np.asarray(container["vx_mps"], dtype=float))) / 100.0
            vy = np.atleast_1d(np.squeeze(np.asarray(container["vy_mps"], dtype=float))) / 100.0
            return vx, vy
        if has_key("vx_Aps") and has_key("vy_Aps"):
            vx = np.atleast_1d(np.squeeze(np.asarray(container["vx_Aps"], dtype=float)))
            vy = np.atleast_1d(np.squeeze(np.asarray(container["vy_Aps"], dtype=float)))
            return vx, vy
        raise ValueError(
            f"{p.name} must contain (vx_mps, vy_mps) or legacy (vx_Aps, vy_Aps) axes"
        )

    if p.suffix.lower() == ".mat":
        z = loadmat(p)
        if "intensity" not in z:
            raise ValueError(f"{p.name} must contain an intensity field")
        vx, vy = _pick(z, lambda k: k in z)
        intensity = np.asarray(z["intensity"], dtype=float)
    else:
        with np.load(p, allow_pickle=False) as z:
            if "intensity" not in z.files:
                raise ValueError(f"{p.name} must contain an intensity field")
            vx, vy = _pick(z, lambda k: k in z.files)
            intensity = np.asarray(z["intensity"], dtype=float)

    if intensity.ndim != 2:
        raise ValueError(f"{p.name} intensity must be a 2-D array")
    axes_are_vectors = vx.ndim == 1 and vy.ndim == 1
    axes_are_grids = vx.shape == intensity.shape and vy.shape == intensity.shape
    if axes_are_vectors:
        expected_shape = (vy.size, vx.size)
        if intensity.shape != expected_shape:
            raise ValueError(
                f"{p.name} intensity shape {intensity.shape} must match "
                f"(len(vy_Aps), len(vx_Aps)) = {expected_shape}"
            )
    elif not axes_are_grids:
        raise ValueError(
            f"{p.name} vx_Aps and vy_Aps must either be 1-D axis vectors "
            f"or 2-D coordinate grids matching intensity shape {intensity.shape}"
        )
    _validate_axis_range(p, "vx_Aps", vx)
    _validate_axis_range(p, "vy_Aps", vy)

    metadata_path = p.with_suffix(".json")
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    return PaperV2VMIImageReference(
        vx_Aps=vx,
        vy_Aps=vy,
        vx_mps=vx * 100.0,
        vy_mps=vy * 100.0,
        intensity=intensity,
        metadata=metadata,
        source_path=p.resolve(),
    )


def load_paper_v2_vmi_polar_image_reference(
    path: str | Path,
) -> PaperV2VMIPolarImageReference:
    """Load a paper-v2 processed polar VMI image from ``.mat`` or ``.npz``.

    Required fields are ``intensity_polar``, ``phi_rad``, and one of
    ``v_radius_mps`` (canonical, m/s on disk) or legacy ``v_radius_Aps``
    (A/ps). The loader converts to A/ps either way and validates that
    ``intensity_polar.shape == (phi_rad.size, v_radius.size)``.
    """

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"paper-v2 polar VMI image reference not found: {p.resolve()}"
        )

    def _pick_v(container, has_key) -> np.ndarray:
        if has_key("v_radius_mps"):
            return (
                np.atleast_1d(np.squeeze(np.asarray(container["v_radius_mps"], dtype=float)))
                / 100.0
            )
        if has_key("v_radius_Aps"):
            return np.atleast_1d(
                np.squeeze(np.asarray(container["v_radius_Aps"], dtype=float))
            )
        raise ValueError(
            f"{p.name} must contain v_radius_mps (or legacy v_radius_Aps)"
        )

    if p.suffix.lower() == ".mat":
        z = loadmat(p)
        if "intensity_polar" not in z:
            raise ValueError(f"{p.name} must contain an intensity_polar field")
        if "phi_rad" not in z:
            raise ValueError(f"{p.name} must contain a phi_rad field")
        phi = np.atleast_1d(np.squeeze(np.asarray(z["phi_rad"], dtype=float)))
        v_radius_Aps = _pick_v(z, lambda k: k in z)
        intensity = np.asarray(z["intensity_polar"], dtype=float)
    else:
        with np.load(p, allow_pickle=False) as z:
            if "intensity_polar" not in z.files:
                raise ValueError(f"{p.name} must contain an intensity_polar field")
            if "phi_rad" not in z.files:
                raise ValueError(f"{p.name} must contain a phi_rad field")
            phi = np.atleast_1d(np.squeeze(np.asarray(z["phi_rad"], dtype=float)))
            v_radius_Aps = _pick_v(z, lambda k: k in z.files)
            intensity = np.asarray(z["intensity_polar"], dtype=float)

    if intensity.ndim != 2:
        raise ValueError(f"{p.name} intensity_polar must be a 2-D array")
    if phi.ndim != 1 or v_radius_Aps.ndim != 1:
        raise ValueError(
            f"{p.name} phi_rad and v_radius axis must be 1-D vectors"
        )
    expected_shape = (phi.size, v_radius_Aps.size)
    if intensity.shape != expected_shape:
        raise ValueError(
            f"{p.name} intensity_polar shape {intensity.shape} must match "
            f"(len(phi_rad), len(v_radius)) = {expected_shape}"
        )
    _validate_axis_range(p, "phi_rad", phi)
    _validate_axis_range(p, "v_radius_Aps", v_radius_Aps)

    metadata_path = p.with_suffix(".json")
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    return PaperV2VMIPolarImageReference(
        phi_rad=phi,
        v_radius_Aps=v_radius_Aps,
        v_radius_mps=v_radius_Aps * 100.0,
        intensity=intensity,
        metadata=metadata,
        source_path=p.resolve(),
    )


def load_paper_v2_phi_reference(path: str | Path) -> PaperV2PhiReference:
    """Load a paper-v2 phi reference CSV with ``phi_rad,signal_arb`` columns."""

    p = Path(path)
    data = _load_named_csv(p)
    names = data.dtype.names or ()
    if not all(name in names for name in ("phi_rad", "signal_arb")):
        raise ValueError(f"{p.name} must contain columns phi_rad and signal_arb")
    return PaperV2PhiReference(
        phi_rad=np.asarray(data["phi_rad"], dtype=float),
        signal_arb=np.asarray(data["signal_arb"], dtype=float),
        source_path=p.resolve(),
    )


def paper_v2_velocity_map(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
) -> PaperV2VelocityMap:
    """Count selected final ``vx, vy`` samples into MATLAB's nearest bins."""

    selected = _paper_v2_atom_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x, dtype=float)[selected]
    vy = np.asarray(ion.velocities_final_y, dtype=float)[selected]
    bins = PAPER_V2_IMAGE_BINS_APS
    counts = np.zeros((bins.size, bins.size), dtype=float)
    for vx_i, vy_i in zip(vx, vy, strict=True):
        ix = int(np.argmin(np.abs(bins - vx_i)))
        iy = int(np.argmin(np.abs(bins - vy_i)))
        counts[ix, iy] += 1.0
    return PaperV2VelocityMap(
        velocity_bins_Aps=bins.copy(),
        counts=counts,
        mass_amu=float(mass_amu),
        num_atoms_used=int(np.count_nonzero(selected)),
    )


def paper_v2_velocity_curve(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
    vmin_Aps: float = 0.0,
) -> PaperV2VelocityCurve:
    """Compute the paper-v2 projected-velocity histogram."""

    selected = _paper_v2_atom_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x, dtype=float)[selected]
    vy = np.asarray(ion.velocities_final_y, dtype=float)[selected]
    projected = np.sqrt(vx * vx + vy * vy)
    projected = projected[projected > float(vmin_Aps)]
    edges = _velocity_edges()
    counts, edges = np.histogram(projected, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = moving_mean(counts, PAPER_V2_VELOCITY_SMOOTHING_WINDOW)
    return PaperV2VelocityCurve(
        mass_amu=float(mass_amu),
        bin_centers_Aps=centers,
        bin_centers_mps=centers * 100.0,
        bin_edges_Aps=edges,
        bin_edges_mps=edges * 100.0,
        counts=counts.astype(int),
        smoothed=smoothed,
        normalised=max_normalise(smoothed),
        num_atoms_used=int(np.count_nonzero(selected)),
        smoothing_window=PAPER_V2_VELOCITY_SMOOTHING_WINDOW,
    )


def paper_v2_phi_curve(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
) -> PaperV2PhiCurve:
    """Compute the paper-v2 simulated phi curve."""

    selected = _paper_v2_atom_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x, dtype=float)[selected]
    vy = np.asarray(ion.velocities_final_y, dtype=float)[selected]
    phi = np.mod(np.arctan2(vy, vx) + np.pi, 2.0 * np.pi)
    edges = np.arange(0.0, 2.0 * np.pi + PAPER_V2_PHI_BIN_WIDTH_RAD, PAPER_V2_PHI_BIN_WIDTH_RAD)
    counts, edges = np.histogram(phi, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = moving_mean(counts, PAPER_V2_PHI_SMOOTHING_WINDOW)
    return PaperV2PhiCurve(
        mass_amu=float(mass_amu),
        bin_centers_rad=centers,
        bin_edges_rad=edges,
        counts=counts.astype(int),
        smoothed=smoothed,
        normalised=max_normalise(smoothed),
        num_atoms_used=int(np.count_nonzero(selected)),
        smoothing_window=PAPER_V2_PHI_SMOOTHING_WINDOW,
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


def _validate_axis_range(path: Path, axis_name: str, axis: np.ndarray) -> None:
    finite = np.asarray(axis, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size and float(np.nanmax(finite) - np.nanmin(finite)) <= 0.0:
        raise ValueError(
            f"{path.name} {axis_name} axis has zero range. This usually means the "
            "2-D VMI reference was exported with a constant row/column slice; "
            "rerun data/reference/scripts/export_paper_v2_reference_data.m "
            "after the full-grid axis export fix."
        )


def _paper_v2_atom_selection(ion: IonCheckpoint, *, mass_amu: float) -> np.ndarray:
    masses_amu = np.round(np.asarray(ion.mass_final_kg, dtype=float) / U_KG)
    mass_mask = masses_amu == round(float(mass_amu))
    outside = np.concatenate([ion.b_ion_outside, ion.b_ion_outside]).astype(bool)
    selected = mass_mask & outside
    if not np.any(selected):
        raise ValueError(
            f"No atoms passed the paper-v2 mass={mass_amu:.0f} amu "
            "and b_ion_outside filter."
        )
    return selected


def _velocity_edges() -> np.ndarray:
    return np.arange(
        0.0,
        PAPER_V2_VELOCITY_MAX_APS + PAPER_V2_VELOCITY_BIN_WIDTH_APS,
        PAPER_V2_VELOCITY_BIN_WIDTH_APS,
    )


def _load_named_csv(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"paper-v2 reference file not found: {path.resolve()}")
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    if data.dtype.names is None:
        raise ValueError(f"{path.name} must have a header row")
    return np.atleast_1d(data)


def _label_from_reference_name(path: Path) -> str:
    special = {
        "iplus_he_high_snr_radial": "I+He high-SNR",
    }
    if path.stem in special:
        return special[path.stem]
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
