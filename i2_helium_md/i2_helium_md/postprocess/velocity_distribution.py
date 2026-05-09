"""Final-velocity histogram + experimental VMI reference loader.

Used by the post-processing plotting script to reproduce the bottom tile
of legacy ``simulation_image.m`` (lines 159-256), which overlays:

    * experimental I+ gas-phase VMI signal,
    * experimental I+He droplet VMI signal,
    * simulated I+He histogram (mass = 131 amu),
    * simulated I+He2 histogram (mass = 135 amu).

Numerical responsibilities live here so the plotting script stays a thin
wrapper around postprocess primitives.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint


_VMI_EXPECTED_COLUMNS: tuple[str, ...] = ("v_Aps", "signal_arb")


@dataclass(frozen=True)
class VmiReference:
    """One experimental VMI reference loaded from a 2-column CSV.

    Attributes
    ----------
    velocity_Aps : np.ndarray, shape (M,)
        Velocity bins in angstrom/ps (the ``v_Aps`` column).
    signal_arb : np.ndarray, shape (M,)
        Signal in arbitrary units (the ``signal_arb`` column). Raw
        unnormalised values; the plotting layer applies its own scaling.
    source_path : Path
        Resolved absolute path of the file that was loaded.
    """

    velocity_Aps: np.ndarray
    signal_arb: np.ndarray
    source_path: Path


def load_vmi_reference(path: str | Path) -> VmiReference:
    """Read a 2-column ``v_Aps,signal_arb`` VMI reference CSV.

    Parameters
    ----------
    path
        Path to the CSV. Both ``data/reference/vmi_iplus_he.csv`` and
        ``data/reference/vmi_iplus_gas.csv`` follow this format.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the header is missing or the columns differ from
        ``("v_Aps", "signal_arb")``.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"VMI reference file not found: {p.resolve()}"
        )

    frame = pd.read_csv(p)
    actual = tuple(frame.columns)
    if actual != _VMI_EXPECTED_COLUMNS:
        raise ValueError(
            f"VMI reference {p.name} has columns {list(actual)}; "
            f"expected {list(_VMI_EXPECTED_COLUMNS)}"
        )

    return VmiReference(
        velocity_Aps=np.asarray(frame["v_Aps"].to_numpy(), dtype=float),
        signal_arb=np.asarray(frame["signal_arb"].to_numpy(), dtype=float),
        source_path=p.resolve(),
    )


@dataclass(frozen=True)
class FinalVelocityHistogram:
    """1-D histogram of |v_final| for atoms passing a mass + outside filter.

    Attributes
    ----------
    bin_centers_Aps : np.ndarray, shape (B,)
        Bin centers in angstrom/ps.
    bin_edges_Aps : np.ndarray, shape (B + 1,)
        Bin edges in angstrom/ps.
    counts : np.ndarray, shape (B,)
        Raw count per bin.
    density : np.ndarray, shape (B,)
        ``counts / bin_width``. Integrates (via trapezoid) to roughly
        ``num_atoms_used``.
    mass_amu : float
        Target atomic mass that was used to filter atoms.
    num_atoms_used : int
        Number of atoms that passed the mass + outside filter.
    """

    bin_centers_Aps: np.ndarray
    bin_edges_Aps: np.ndarray
    counts: np.ndarray
    density: np.ndarray
    mass_amu: float
    num_atoms_used: int


def compute_final_velocity_histogram(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
    num_bins: int = 120,
    v_max_Aps: float = 28.0,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> FinalVelocityHistogram:
    """Bin ``|v_final|`` over atoms whose final mass matches ``mass_amu``.

    Mirrors the simulation block in
    ``legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
    simulation_image.m:213-256``:

        mass_select(j) % e.g. 131 (I+He) or 135 (I+He2)
        b_use = round(data_ion.mass_i(:,end)/u) == mass_select & ...
                data_ion.b_ion_outside == 1
        v_total = sqrt(vx^2+vy^2+vz^2)
        bayes_hist(v_total, x_bins, ...)

    The legacy ``bayes_hist`` adds Bayesian uncertainty bars; we use
    plain ``np.histogram`` because the plotting layer only displays the
    central trace.

    Parameters
    ----------
    ion
        Ion-stage checkpoint. ``velocities_final_x/y/z`` (shape
        ``(2*N,)``), ``mass_final_kg`` (shape ``(2*N,)``), and
        ``b_ion_outside`` (shape ``(N,)``, per molecule) are read.
    mass_amu
        Target atomic mass in amu. Atoms whose
        ``round(mass_final_kg / U_KG)`` is within
        ``mass_tolerance_amu`` of ``mass_amu`` are kept.
    num_bins
        Number of histogram bins between 0 and ``v_max_Aps``.
    v_max_Aps
        Upper edge of the histogram. Atoms with ``|v| > v_max_Aps`` end
        up in the overflow and are excluded from the returned bins
        (matches the MATLAB plot xlim).
    mass_tolerance_amu
        Half-width of the mass-acceptance window in amu.
    require_outside
        If ``True`` (default), only atoms belonging to molecules whose
        ``b_ion_outside`` flag is set contribute. The flag is
        per-molecule; both atoms of an "outside" molecule pass.

    Returns
    -------
    FinalVelocityHistogram

    Raises
    ------
    ValueError
        If no atoms pass the mass + outside filter, or if
        ``v_max_Aps <= 0`` or ``num_bins < 1``.
    """
    if num_bins < 1:
        raise ValueError(f"num_bins must be >= 1, got {num_bins}")
    if v_max_Aps <= 0.0:
        raise ValueError(f"v_max_Aps must be > 0, got {v_max_Aps}")

    n = ion.num_molecules
    mass_amu_per_atom = np.round(np.asarray(ion.mass_final_kg) / U_KG)
    mass_mask = np.abs(mass_amu_per_atom - mass_amu) <= mass_tolerance_amu

    if require_outside:
        outside_per_atom = np.concatenate(
            [ion.b_ion_outside, ion.b_ion_outside]
        ).astype(bool)
        select = mass_mask & outside_per_atom
    else:
        select = mass_mask

    num_used = int(np.count_nonzero(select))
    if num_used == 0:
        raise ValueError(
            f"No atoms passed the mass={mass_amu} amu "
            f"(tolerance {mass_tolerance_amu}) "
            f"{'+ outside' if require_outside else ''} filter; "
            f"check the run or relax the tolerance."
        )

    speed_final = np.sqrt(
        np.asarray(ion.velocities_final_x)[select] ** 2
        + np.asarray(ion.velocities_final_y)[select] ** 2
        + np.asarray(ion.velocities_final_z)[select] ** 2
    )

    edges = np.linspace(0.0, v_max_Aps, num_bins + 1)
    counts, _ = np.histogram(speed_final, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]
    density = counts.astype(float) / bin_width

    return FinalVelocityHistogram(
        bin_centers_Aps=centers,
        bin_edges_Aps=edges,
        counts=counts,
        density=density,
        mass_amu=float(mass_amu),
        num_atoms_used=num_used,
    )


@dataclass(frozen=True)
class BimodalGaussianFit:
    """Sum-of-two-Gaussians decomposition of a 1-D velocity histogram.

    Reproduces the high-velocity tail decomposition of
    ``legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper.m``
    (``nlinfit(v, y, fun, [814, 1000, 1, 2000, 100, 1])`` style six-
    parameter fit reduced here to two centred Gaussians).

    Model: ``f(v) = A1 * exp(-(v - mu1)^2 / (2 sig1^2))
                  + A2 * exp(-(v - mu2)^2 / (2 sig2^2))``.

    All parameters are ``np.nan`` and ``success`` is ``False`` when the
    fit does not converge or the input has fewer non-zero bins than
    parameters.
    """

    amplitude_1: float
    mean_1_Aps: float
    sigma_1_Aps: float
    amplitude_2: float
    mean_2_Aps: float
    sigma_2_Aps: float
    fitted_curve: np.ndarray
    residual: float
    success: bool


def _double_gaussian(
    v: np.ndarray,
    a1: float, mu1: float, s1: float,
    a2: float, mu2: float, s2: float,
) -> np.ndarray:
    return (
        a1 * np.exp(-((v - mu1) ** 2) / (2.0 * s1 * s1))
        + a2 * np.exp(-((v - mu2) ** 2) / (2.0 * s2 * s2))
    )


def bimodal_gaussian_fit(
    histogram: FinalVelocityHistogram,
) -> BimodalGaussianFit:
    """Fit two Gaussians to the histogram density.

    Initial guesses pick the two largest local maxima as ``mu1, mu2`` and
    use one quarter of the histogram support as both ``sigma`` seeds.
    Bounds keep the means inside the bin range and sigmas positive.
    """
    v = np.asarray(histogram.bin_centers_Aps, dtype=float)
    y = np.asarray(histogram.density, dtype=float)
    nan_curve = np.full_like(y, np.nan)

    if v.size < 6 or np.count_nonzero(y > 0.0) < 6:
        return BimodalGaussianFit(
            amplitude_1=float("nan"), mean_1_Aps=float("nan"),
            sigma_1_Aps=float("nan"),
            amplitude_2=float("nan"), mean_2_Aps=float("nan"),
            sigma_2_Aps=float("nan"),
            fitted_curve=nan_curve, residual=float("nan"), success=False,
        )

    v_min = float(v[0])
    v_max = float(v[-1])
    span = v_max - v_min
    sigma_seed = max(span / 8.0, (v[1] - v[0]) if v.size > 1 else 1.0)

    # Two seeded peaks: weighted mean of the lower and upper halves.
    weights = np.clip(y, 0.0, None)
    mid = v_min + 0.5 * span
    lower = v <= mid
    upper = ~lower
    w_low = weights[lower].sum()
    w_up = weights[upper].sum()
    mu1_seed = (
        float(np.average(v[lower], weights=weights[lower]))
        if w_low > 0 else v_min + 0.25 * span
    )
    mu2_seed = (
        float(np.average(v[upper], weights=weights[upper]))
        if w_up > 0 else v_min + 0.75 * span
    )
    a_seed = float(y.max())

    p0 = (a_seed, mu1_seed, sigma_seed, a_seed, mu2_seed, sigma_seed)
    bounds = (
        (0.0, v_min, 1e-3, 0.0, v_min, 1e-3),
        (np.inf, v_max, span, np.inf, v_max, span),
    )

    try:
        popt, _pcov = curve_fit(
            _double_gaussian, v, y,
            p0=p0, bounds=bounds, maxfev=10_000,
        )
    except (RuntimeError, ValueError):
        return BimodalGaussianFit(
            amplitude_1=float("nan"), mean_1_Aps=float("nan"),
            sigma_1_Aps=float("nan"),
            amplitude_2=float("nan"), mean_2_Aps=float("nan"),
            sigma_2_Aps=float("nan"),
            fitted_curve=nan_curve, residual=float("nan"), success=False,
        )

    fit_curve = _double_gaussian(v, *popt)
    residual = float(np.sqrt(np.mean((fit_curve - y) ** 2)))
    a1, mu1, s1, a2, mu2, s2 = (float(p) for p in popt)
    # Order so that amplitude_1 is always the dominant peak.
    if a2 > a1:
        a1, mu1, s1, a2, mu2, s2 = a2, mu2, s2, a1, mu1, s1
    return BimodalGaussianFit(
        amplitude_1=a1, mean_1_Aps=mu1, sigma_1_Aps=s1,
        amplitude_2=a2, mean_2_Aps=mu2, sigma_2_Aps=s2,
        fitted_curve=fit_curve, residual=residual, success=True,
    )
