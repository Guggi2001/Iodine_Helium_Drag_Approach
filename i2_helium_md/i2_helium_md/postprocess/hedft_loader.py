"""Loader for HeDFT/TDDFT reference trajectory CSV files.

Reads the normalized 8-column reference format produced from the legacy
MATLAB pipeline:

    Time_ps, V1_mag, V2_mag, V1_z, V2_z, V1_x, V2_x, R_distance

Two reference files exist under ``data/reference/``:

    9A_All_Data.csv    -- 9 angstrom droplet
    18A_All_Data.csv   -- 18 angstrom droplet

This module supersedes the split legacy importers
``9Angstroem/importfile_v2.m``, ``9Angstroem/importfile_R1_R2.m``, and
``18Angstroem/importfile_marti.m``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


_EXPECTED_COLUMNS: tuple[str, ...] = (
    "Time_ps",
    "V1_mag",
    "V2_mag",
    "V1_z",
    "V2_z",
    "V1_x",
    "V2_x",
    "R_distance",
)

# Pattern matches a leading "<int>A_" or "<float>A_" prefix in the filename
# stem, e.g. ``9A_All_Data`` or ``18A_All_Data``.
_RADIUS_FROM_STEM = re.compile(r"^(\d+(?:\.\d+)?)A_", re.IGNORECASE)


@dataclass(frozen=True)
class HedftTrajectory:
    """One HeDFT reference trajectory loaded from an ``*_All_Data.csv`` file.

    All arrays share the same length T = len(time_ps). Units mirror the MD
    pipeline so no conversion is needed at comparison time.

    Attributes
    ----------
    time_ps : np.ndarray, shape (T,)
        Time grid in picoseconds, strictly increasing.
    v1_magnitude_Aps, v2_magnitude_Aps : np.ndarray, shape (T,)
        Velocity magnitudes |v| of atoms I1 and I2 in angstrom/ps.
    v1_z_Aps, v2_z_Aps : np.ndarray, shape (T,)
        z-components of the I1 and I2 velocities, angstrom/ps.
    v1_x_Aps, v2_x_Aps : np.ndarray, shape (T,)
        x-components of the I1 and I2 velocities, angstrom/ps.
    distance_A : np.ndarray, shape (T,)
        I-I separation in angstrom (the legacy ``R_distance`` column).
    droplet_radius_A : float
        9.0 or 18.0 angstrom; inferred from the filename prefix when
        not supplied explicitly.
    source_path : Path
        Absolute path of the file that was loaded.
    """

    time_ps: np.ndarray
    v1_magnitude_Aps: np.ndarray
    v2_magnitude_Aps: np.ndarray
    v1_z_Aps: np.ndarray
    v2_z_Aps: np.ndarray
    v1_x_Aps: np.ndarray
    v2_x_Aps: np.ndarray
    distance_A: np.ndarray
    droplet_radius_A: float
    source_path: Path


def load_hedft_trajectory(
    path: str | Path,
    *,
    droplet_radius_A: float | None = None,
) -> HedftTrajectory:
    """Load an 8-column HeDFT reference CSV into a :class:`HedftTrajectory`.

    Parameters
    ----------
    path
        Path to the CSV. Must exist; ``FileNotFoundError`` otherwise.
    droplet_radius_A
        Optional override. If ``None``, the radius is inferred from the
        filename prefix (e.g. ``9A_All_Data.csv`` -> 9.0). If both
        inference fails *and* this argument is omitted, ``ValueError``
        is raised.

    Returns
    -------
    HedftTrajectory
        Frozen dataclass with the eight column arrays plus metadata.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the header is missing/has wrong columns, the time column is
        not strictly increasing, or the droplet radius cannot be
        determined.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"HeDFT reference file not found: {p.resolve()}"
        )

    # numpy.genfromtxt with names=True normalises the header into a
    # structured array. The expected column names are already valid
    # Python identifiers, so no munging is needed.
    structured = np.genfromtxt(
        p,
        delimiter=",",
        names=True,
        dtype=float,
    )

    actual_columns = tuple(structured.dtype.names or ())
    missing = tuple(c for c in _EXPECTED_COLUMNS if c not in actual_columns)
    extra = tuple(c for c in actual_columns if c not in _EXPECTED_COLUMNS)
    if missing or extra:
        raise ValueError(
            f"HeDFT reference {p.name} has unexpected columns. "
            f"missing={list(missing)}, unexpected={list(extra)}, "
            f"expected={list(_EXPECTED_COLUMNS)}"
        )

    time_ps = np.asarray(structured["Time_ps"], dtype=float)
    if time_ps.ndim != 1 or time_ps.size < 2:
        raise ValueError(
            f"HeDFT reference {p.name} must have at least 2 time samples, "
            f"got shape {time_ps.shape}"
        )
    if not np.all(np.diff(time_ps) > 0.0):
        raise ValueError(
            f"HeDFT reference {p.name} has a non-monotonic Time_ps column"
        )

    radius = _resolve_droplet_radius(p, override=droplet_radius_A)

    return HedftTrajectory(
        time_ps=time_ps,
        v1_magnitude_Aps=np.asarray(structured["V1_mag"], dtype=float),
        v2_magnitude_Aps=np.asarray(structured["V2_mag"], dtype=float),
        v1_z_Aps=np.asarray(structured["V1_z"], dtype=float),
        v2_z_Aps=np.asarray(structured["V2_z"], dtype=float),
        v1_x_Aps=np.asarray(structured["V1_x"], dtype=float),
        v2_x_Aps=np.asarray(structured["V2_x"], dtype=float),
        distance_A=np.asarray(structured["R_distance"], dtype=float),
        droplet_radius_A=float(radius),
        source_path=p.resolve(),
    )


def _resolve_droplet_radius(
    path: Path,
    *,
    override: float | None,
) -> float:
    if override is not None:
        return float(override)

    match = _RADIUS_FROM_STEM.match(path.stem)
    if match is None:
        raise ValueError(
            f"Cannot infer droplet radius from filename {path.name!r}; "
            f"pass droplet_radius_A explicitly. Expected a stem like "
            f"'9A_All_Data' or '18A_All_Data'."
        )
    return float(match.group(1))
