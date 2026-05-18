"""Literal helpers for ``post_process_single_pulse_paper_IplusHe_comparison_cov.m``.

The MATLAB ``_cov`` script extends the IplusHe comparison with two
pair-covariance diagnostics computed from the experimental I+He droplet
VMI measurements (IDs 45668, 45662, 45667). The Python port treats the
experimental covariance matrices as frozen reference data exported once
from MATLAB (see ``data/reference/scripts/export_paper_cov_reference_data.m``)
and adds a simulated counterpart computed from the Python ion checkpoint
so the figure can be a side-by-side experiment / simulation comparison.

The two simulated helpers in this module mirror the
``generate_VMI_covariance_matrices`` recipe of the legacy script for the
angular and radial-speed pair-covariance respectively. The angular
side reuses :func:`paper_v4_angular_pair_covariance` from the v4 module
already; this module only adds the radial-speed counterpart and the
reference dataclass + loader.

Conventions preserved verbatim from ``_cov.m``:

* pair indexing ``[0, N) <-> [N, 2 N)``,
* mass selection by ``round(mass_amu / U)`` plus ``b_ion_outside``,
* radial speed ``v_r = sqrt(vx^2 + vy^2)`` in A/ps,
* diagonal removal ``cov - diag(diag(cov))``,
* 2 x 2 ``movmean`` smoothing applied along each axis.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint
from ._smoothing import moving_mean
from .paper_v2 import max_normalise


PAPER_COV_VELOCITY_BINS = 90
PAPER_COV_VELOCITY_MAX_APS = 30.0
PAPER_COV_SMOOTHING_WINDOW = 2
PAPER_COV_VELOCITY_TRACE_MIN_APS = 4.0
PAPER_COV_VELOCITY_TRACE_MAX_APS = 22.0
PAPER_COV_COLOR_CLIP_FRACTION = 0.7
PAPER_COV_PHI_BINS = 126
PAPER_COV_PHI_SMOOTHING_WINDOW = 15
PAPER_COV_TRACE_SMOOTHING_WINDOW = 3


@dataclass(frozen=True)
class PaperCovExperimentalReference:
    """MATLAB-exported experimental pair covariance for I+He droplet.

    ``cov_angular`` has its diagonal already zeroed; ``cov_radial`` has
    its diagonal zeroed AND has the 2 x 2 ``movmean`` smoothing already
    applied (these are the on-disk states; the loader does no further
    post-processing).
    """

    cov_angular: np.ndarray
    cov_radial: np.ndarray
    theta_centers_rad: np.ndarray
    velocity_centers_Aps: np.ndarray
    velocity_centers_mps: np.ndarray
    metadata: dict[str, Any]
    source_path: Path


@dataclass(frozen=True)
class PhiAngularDistribution:
    """1-D phi histogram normalised so the peak equals one.

    Mirrors the MATLAB ``_cov.m`` overlay at lines 178-205 (experimental
    polar-image mean) and 240-252 (simulated ``atan2(vy, vx) + pi``
    histogram). The two recipes share the same bin grid in
    ``[0, 2 pi)`` so the two curves can be overlaid on a common axis.
    """

    signal_normalised: np.ndarray
    phi_centers_rad: np.ndarray
    num_samples_used: int
    smoothing_window: int


@dataclass(frozen=True)
class RadialPairCovariance:
    """Simulated radial-speed pair covariance ``(v_r_a, v_r_b)``.

    Mirrors the radial branch of ``generate_VMI_covariance_matrices``:
    speed-speed 2-D histogram, optional diagonal zeroing, optional 2 x 2
    moving mean along each axis.
    """

    counts: np.ndarray
    velocity_centers_Aps: np.ndarray
    velocity_centers_mps: np.ndarray
    velocity_edges_Aps: np.ndarray
    num_pairs_used: int
    mass_amu: float
    smoothing_window: int


def load_paper_cov_experimental_reference(
    path: str | Path,
) -> PaperCovExperimentalReference:
    """Load the MATLAB-exported pair covariance reference (``.mat`` or ``.npz``).

    Required keys: ``cov_angular``, ``cov_radial``, ``theta_centers_rad``,
    and one of ``velocity_centers_mps`` (canonical, m/s) or legacy
    ``velocity_centers_Aps`` (A/ps).
    """

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"paper-cov experimental reference not found: {p.resolve()}"
        )

    def _has(container, key, files_attr=None):
        if files_attr is not None:
            return key in files_attr
        return key in container

    def _pick_velocity(container, has) -> np.ndarray:
        if has("velocity_centers_mps"):
            return (
                np.atleast_1d(
                    np.squeeze(np.asarray(container["velocity_centers_mps"], dtype=float))
                )
                / 100.0
            )
        if has("velocity_centers_Aps"):
            return np.atleast_1d(
                np.squeeze(np.asarray(container["velocity_centers_Aps"], dtype=float))
            )
        raise ValueError(
            f"{p.name} must contain velocity_centers_mps (or legacy velocity_centers_Aps)"
        )

    if p.suffix.lower() == ".mat":
        z = loadmat(p)
        has = lambda k: k in z  # noqa: E731
        for required in ("cov_angular", "cov_radial", "theta_centers_rad"):
            if not has(required):
                raise ValueError(f"{p.name} must contain a {required} field")
        cov_angular = np.asarray(z["cov_angular"], dtype=float)
        cov_radial = np.asarray(z["cov_radial"], dtype=float)
        theta_centers = np.atleast_1d(
            np.squeeze(np.asarray(z["theta_centers_rad"], dtype=float))
        )
        velocity_centers_Aps = _pick_velocity(z, has)
    else:
        with np.load(p, allow_pickle=False) as z:
            files = z.files
            has = lambda k: k in files  # noqa: E731
            for required in ("cov_angular", "cov_radial", "theta_centers_rad"):
                if not has(required):
                    raise ValueError(f"{p.name} must contain a {required} field")
            cov_angular = np.asarray(z["cov_angular"], dtype=float)
            cov_radial = np.asarray(z["cov_radial"], dtype=float)
            theta_centers = np.atleast_1d(
                np.squeeze(np.asarray(z["theta_centers_rad"], dtype=float))
            )
            velocity_centers_Aps = _pick_velocity(z, has)

    if cov_angular.ndim != 2 or cov_angular.shape[0] != cov_angular.shape[1]:
        raise ValueError(
            f"{p.name} cov_angular must be a square 2-D array, got shape {cov_angular.shape}"
        )
    if cov_radial.ndim != 2 or cov_radial.shape[0] != cov_radial.shape[1]:
        raise ValueError(
            f"{p.name} cov_radial must be a square 2-D array, got shape {cov_radial.shape}"
        )
    if theta_centers.size != cov_angular.shape[0]:
        raise ValueError(
            f"{p.name} theta_centers_rad length {theta_centers.size} must match "
            f"cov_angular size {cov_angular.shape[0]}"
        )
    if velocity_centers_Aps.size != cov_radial.shape[0]:
        raise ValueError(
            f"{p.name} velocity_centers length {velocity_centers_Aps.size} must match "
            f"cov_radial size {cov_radial.shape[0]}"
        )

    metadata_path = p.with_suffix(".json")
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    return PaperCovExperimentalReference(
        cov_angular=cov_angular,
        cov_radial=cov_radial,
        theta_centers_rad=theta_centers,
        velocity_centers_Aps=velocity_centers_Aps,
        velocity_centers_mps=velocity_centers_Aps * 100.0,
        metadata=metadata,
        source_path=p.resolve(),
    )


def radial_pair_speed_covariance(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
    n_velocity_bins: int = PAPER_COV_VELOCITY_BINS,
    v_max_Aps: float = PAPER_COV_VELOCITY_MAX_APS,
    smoothing_window: int = PAPER_COV_SMOOTHING_WINDOW,
    remove_diagonal: bool = True,
) -> RadialPairCovariance:
    """Simulated 2-D speed-speed pair covariance for mass-selected ion pairs.

    Mirrors the radial branch of ``generate_VMI_covariance_matrices`` in
    the legacy ``_cov.m`` script: for each surviving molecule the two
    fragment speeds ``sqrt(vx^2 + vy^2)`` are binned into a 2-D
    histogram, the diagonal is zeroed (suppressing self-binning), then a
    2 x 2 moving mean is applied along each axis.

    Parameters
    ----------
    ion
        Final-state checkpoint.
    mass_amu
        Mass selection (rounded to integer amu).
    n_velocity_bins
        Number of velocity bins on each axis.
    v_max_Aps
        Upper edge of the velocity range (A/ps).
    smoothing_window
        Movmean window applied along each axis (set 0 or 1 to skip).
    remove_diagonal
        If True, zero the main diagonal before smoothing.
    """

    if n_velocity_bins < 1:
        raise ValueError(f"n_velocity_bins must be >= 1, got {n_velocity_bins}")
    if v_max_Aps <= 0.0:
        raise ValueError(f"v_max_Aps must be > 0, got {v_max_Aps}")

    n = ion.num_molecules
    atom_ok = _paper_cov_atom_selection(ion, mass_amu=mass_amu)
    pair_ok = atom_ok[:n] & atom_ok[n:]

    vx_a = np.asarray(ion.velocities_final_x[:n], dtype=float)[pair_ok]
    vy_a = np.asarray(ion.velocities_final_y[:n], dtype=float)[pair_ok]
    vx_b = np.asarray(ion.velocities_final_x[n:], dtype=float)[pair_ok]
    vy_b = np.asarray(ion.velocities_final_y[n:], dtype=float)[pair_ok]
    v_r_a = np.sqrt(vx_a * vx_a + vy_a * vy_a)
    v_r_b = np.sqrt(vx_b * vx_b + vy_b * vy_b)

    edges = np.linspace(0.0, float(v_max_Aps), n_velocity_bins + 1)
    counts, _, _ = np.histogram2d(v_r_a, v_r_b, bins=(edges, edges))
    counts = counts.astype(float)
    if remove_diagonal:
        np.fill_diagonal(counts, 0.0)
    if smoothing_window and smoothing_window > 1:
        counts = _movmean_2d(counts, smoothing_window)

    centers = 0.5 * (edges[:-1] + edges[1:])
    return RadialPairCovariance(
        counts=counts,
        velocity_centers_Aps=centers,
        velocity_centers_mps=centers * 100.0,
        velocity_edges_Aps=edges,
        num_pairs_used=int(np.count_nonzero(pair_ok)),
        mass_amu=float(mass_amu),
        smoothing_window=int(smoothing_window),
    )


def radial_covariance_trace(
    cov_radial: np.ndarray,
    velocity_centers_Aps: np.ndarray,
    *,
    v_min_Aps: float = PAPER_COV_VELOCITY_TRACE_MIN_APS,
    v_max_Aps: float = PAPER_COV_VELOCITY_TRACE_MAX_APS,
) -> np.ndarray:
    """1-D trace extracted from ``cov_radial`` per the ``_cov.m`` recipe.

    MATLAB ``_cov.m`` (lines 419-429) selects rows whose corresponding
    velocity falls in the window ``[v_min, v_max]`` and sums them
    along axis 0, halving the result for symmetry::

        b_v = (v > v_min) & (v < v_max)
        trace = sum(cov_radial[b_v, :], axis=0) / 2

    Returns a 1-D array of length ``cov_radial.shape[1]``.
    """

    v = np.asarray(velocity_centers_Aps, dtype=float)
    cov = np.asarray(cov_radial, dtype=float)
    if v.size != cov.shape[1]:
        raise ValueError(
            f"velocity_centers length {v.size} must match cov_radial columns {cov.shape[1]}"
        )
    band = (v > float(v_min_Aps)) & (v < float(v_max_Aps))
    if not np.any(band):
        return np.zeros(cov.shape[1], dtype=float)
    return cov[band, :].sum(axis=0) / 2.0


def simulated_phi_distribution(
    ion: IonCheckpoint,
    *,
    mass_amu: float = 131.0,
    n_phi_bins: int = PAPER_COV_PHI_BINS,
    smoothing_window: int = PAPER_COV_PHI_SMOOTHING_WINDOW,
) -> PhiAngularDistribution:
    """Simulated 1-D phi distribution for mass-selected ion final velocities.

    Mirrors lines 240-252 of ``_cov.m``: take the ``b_mass & b_outside``
    atom selection, compute ``phi = atan2(vy, vx) + pi`` wrapped into
    ``[0, 2 pi)``, histogram on edges ``0:2*pi/n_phi_bins:2*pi``
    (defaults to the legacy ``0:0.05:2*pi`` -> 126 bins), apply
    ``movmean(h, 15)``, and normalise by the maximum.
    """

    if n_phi_bins < 1:
        raise ValueError(f"n_phi_bins must be >= 1, got {n_phi_bins}")

    atom_ok = _paper_cov_atom_selection(ion, mass_amu=mass_amu)
    vx = np.asarray(ion.velocities_final_x, dtype=float)[atom_ok]
    vy = np.asarray(ion.velocities_final_y, dtype=float)[atom_ok]
    phi = np.mod(np.arctan2(vy, vx) + np.pi, 2.0 * np.pi)

    edges = np.linspace(0.0, 2.0 * np.pi, n_phi_bins + 1)
    counts, _ = np.histogram(phi, bins=edges)
    counts = counts.astype(float)
    if smoothing_window and smoothing_window > 1:
        counts = moving_mean(counts, smoothing_window)
    centers = 0.5 * (edges[:-1] + edges[1:])

    return PhiAngularDistribution(
        signal_normalised=max_normalise(counts),
        phi_centers_rad=centers,
        num_samples_used=int(np.count_nonzero(atom_ok)),
        smoothing_window=int(smoothing_window),
    )


def covariance_axis_sum_normalised(
    matrix: np.ndarray,
    *,
    axis: int = 0,
    smoothing_window: int = PAPER_COV_TRACE_SMOOTHING_WINDOW,
) -> np.ndarray:
    """1-D axis-sum trace of a pair-covariance matrix, MATLAB-normalised.

    Mirrors the recipe used at lines 498-522 of ``_cov.m`` for both the
    angular and radial pair-covariance traces overlaid on the second
    figure of the legacy script::

        trace = sum(matrix, axis)        # MATLAB ``sum(., 1)`` == numpy axis=0
        trace = movmean(trace, k)
        trace = trace - min(trace)
        trace = trace / max(trace)

    Returns a 1-D ``float`` array. If the post-baseline trace has a
    non-positive max, zeros are returned (matches ``max_normalise``'s
    safety behaviour).
    """

    arr = np.asarray(matrix, dtype=float)
    trace = arr.sum(axis=int(axis))
    if smoothing_window and smoothing_window > 1:
        trace = moving_mean(trace, smoothing_window)
    finite = trace[np.isfinite(trace)]
    if finite.size == 0:
        return np.zeros_like(trace, dtype=float)
    trace = trace - float(np.nanmin(trace))
    peak = float(np.nanmax(trace))
    if not np.isfinite(peak) or peak <= 0.0:
        return np.zeros_like(trace, dtype=float)
    return trace / peak


def _movmean_2d(values: np.ndarray, window: int) -> np.ndarray:
    """Apply :func:`moving_mean` along axis 0 then axis 1.

    Equivalent to the MATLAB pair ``movmean(X, k, 1); movmean(X, k, 2)``.
    """

    arr = np.asarray(values, dtype=float)
    out = np.apply_along_axis(lambda col: moving_mean(col, window), 0, arr)
    out = np.apply_along_axis(lambda row: moving_mean(row, window), 1, out)
    return out


def _paper_cov_atom_selection(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
) -> np.ndarray:
    """Mass + ``b_ion_outside`` filter mirroring the v4 atom selection."""

    masses_amu = np.round(np.asarray(ion.mass_final_kg, dtype=float) / U_KG)
    mass_mask = masses_amu == round(float(mass_amu))
    outside = np.concatenate([ion.b_ion_outside, ion.b_ion_outside]).astype(bool)
    return mass_mask & outside
