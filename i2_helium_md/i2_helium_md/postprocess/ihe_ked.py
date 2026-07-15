"""Loaders + sim-side helpers for the I+He_n kinetic-energy reference.

Reads the frozen experimental export under ``data/reference/ihe_ked/``
(exporter: MATLAB ``export_IHe_KED_reference.m``, measurement 17.10.24,
FLIR camera; full provenance in that folder's README/provenance JSON):

* ``IHe_KED_reference.csv`` -- per-fragment first moments <E> of the 3-D
  kinetic-energy distribution for n = 0..17, with the full error model
  (per-point stat + analysis-systematic errors, two correlated fractional
  scale bands, trust flags).
* ``IHe_KED_curves_n{0..4}.csv`` -- trusted curves for n = 0..4: 2-D
  detector projection and 3-D reconstruction, each per unit speed and per
  unit energy, peak-normalized, on a dual ``E_eV``/``v_mps`` axis pair.

Comparison conventions (reference README, binding for consumers):

* compare **mean-to-mean**, never mean-to-peak -- the distributions are
  strongly skewed;
* the MD 3-D ``|v|`` histogram maps onto ``signal_3d_Pv``; the in-plane
  projected speed maps onto ``signal_2d_Pv``;
* mass model ``m(n) = 126.90 + 4.0026 n`` u
  (:func:`i2_helium_md.physics.shell_schedule.complex_mass_amu` -- matches
  the reference's ``massCenter_u_per_e`` column);
* both sides are droplet-rest-frame;
* curves are independently peak-normalized: shapes only, never amplitudes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..physics.constants import EV, U as U_KG
from ..physics.shell_schedule import complex_mass_amu
from ..simulation.checkpoint import IonCheckpoint
from .csv_contract import validate_columns
from .velocity_distribution import select_final_mass_gate


_REFERENCE_COLUMNS: tuple[str, ...] = (
    "n",
    "label",
    "massCenter_u_per_e",
    "N_counts",
    "N_eff",
    "meanKE_eV",
    "modeKE_eV",
    "medianKE_eV",
    "sigmaKE_eV",
    "statErr_meanKE_eV",
    "sysErr_meanKE_eV",
    "calibSyst_frac",
    "conditionSyst_frac",
    "bgOffShift_eV",
    "dominantError",
    "noiseLimited",
)

# Legal dominantError tags (reference COLUMNS.md).
_LEGAL_DOMINANT_ERROR: frozenset[str] = frozenset(
    {"stat", "analysis", "calib", "bg-structural"}
)

# Fragments with trusted full curves (IHe_KED_curves_n{0..4}.csv).
CURVE_N_MAX: int = 4


@dataclass(frozen=True)
class IHeKedReference:
    """Per-fragment <E> reference table (one array entry per fragment n).

    All energies in eV. ``calib_syst_frac`` / ``condition_syst_frac`` are
    fractional **correlated** scale bands (they shift the whole curve, not
    individual points; never fold them into per-point error bars).

    Attributes
    ----------
    n : np.ndarray, shape (Nn,), int
        He count, contiguous ascending from 0 (Nn = 18 in the frozen file).
    label : np.ndarray, shape (Nn,), str
        Fragment name ``I+He_n``.
    mass_center_u : np.ndarray, shape (Nn,), float
        Fragment mass in u: 126.90 + 4.0026 n.
    N_counts, N_eff : np.ndarray, shape (Nn,), float
        Raw gate counts / background-subtracted in-mask counts.
    mean_KE_eV : np.ndarray, shape (Nn,), float
        THE reference values -- first moments of the unsmoothed 3-D P(E).
    mode_KE_eV, median_KE_eV, sigma_KE_eV : np.ndarray, shape (Nn,), float
        Shape descriptors (sigma is a physical width, NOT an error bar).
    stat_err_mean_KE_eV, sys_err_mean_KE_eV : np.ndarray, shape (Nn,), float
        Per-point uncertainties of the mean; combine in quadrature.
    calib_syst_frac, condition_syst_frac : np.ndarray, shape (Nn,), float
        Correlated fractional scale bands (multiply by ``mean_KE_eV``).
    bg_off_shift_eV : np.ndarray, shape (Nn,), float
        Diagnostic-only background-off sensitivity; not an uncertainty.
    dominant_error : np.ndarray, shape (Nn,), str
        Largest error source per row (one of stat/analysis/calib/
        bg-structural).
    noise_limited : np.ndarray, shape (Nn,), bool
        True would flag <E> as an upper bound (0 for all frozen rows).
    source_path : Path
        Resolved path of the loaded CSV.
    """

    n: np.ndarray
    label: np.ndarray
    mass_center_u: np.ndarray
    N_counts: np.ndarray
    N_eff: np.ndarray
    mean_KE_eV: np.ndarray
    mode_KE_eV: np.ndarray
    median_KE_eV: np.ndarray
    sigma_KE_eV: np.ndarray
    stat_err_mean_KE_eV: np.ndarray
    sys_err_mean_KE_eV: np.ndarray
    calib_syst_frac: np.ndarray
    condition_syst_frac: np.ndarray
    bg_off_shift_eV: np.ndarray
    dominant_error: np.ndarray
    noise_limited: np.ndarray
    source_path: Path

    @property
    def point_err_eV(self) -> np.ndarray:
        """Per-point error of the mean: sqrt(stat^2 + sys^2) [eV]."""
        return np.sqrt(
            self.stat_err_mean_KE_eV ** 2 + self.sys_err_mean_KE_eV ** 2
        )

    @property
    def gold_mask(self) -> np.ndarray:
        """Calib-limited ("gold") fragments: a disagreement there is real
        physics beyond the two correlated bands (reference README)."""
        return self.dominant_error == "calib"


def load_ihe_ked_reference(path: str | Path) -> IHeKedReference:
    """Load ``IHe_KED_reference.csv`` into an :class:`IHeKedReference`.

    Parameters
    ----------
    path
        Path to the CSV. Must exist; ``FileNotFoundError`` otherwise.

    Returns
    -------
    IHeKedReference
        Frozen dataclass with the per-fragment KE moments and full error
        model, one array entry per fragment n = 0..17.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        Missing/extra columns; non-finite numeric entries; ``n`` not
        contiguous ascending from 0; a non-positive ``meanKE_eV``; a
        fractional band outside [0, 1]; an unknown ``dominantError`` tag;
        or ``noiseLimited`` not in {0, 1}.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"IHe KED reference file not found: {p.resolve()}"
        )

    frame = pd.read_csv(p)
    validate_columns(
        list(frame.columns),
        _REFERENCE_COLUMNS,
        file_label=f"IHe KED reference {p.name}",
    )

    def _finite(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"IHe KED reference {p.name} column {name!r} has "
                f"non-finite entries: {values}."
            )
        return values

    n_float = _finite("n")
    n_int = np.rint(n_float).astype(int)
    if n_float.size == 0 or np.any(n_float != n_int) or not np.array_equal(
        n_int, np.arange(n_int.size)
    ):
        raise ValueError(
            f"IHe KED reference {p.name} n column must be contiguous "
            f"integers ascending from 0; got {n_float.tolist()}."
        )

    mean_KE_eV = _finite("meanKE_eV")
    if np.any(mean_KE_eV <= 0.0):
        raise ValueError(
            f"IHe KED reference {p.name} has non-positive meanKE_eV: "
            f"{mean_KE_eV[mean_KE_eV <= 0.0]}."
        )

    for frac_name in ("calibSyst_frac", "conditionSyst_frac"):
        frac = _finite(frac_name)
        if np.any((frac < 0.0) | (frac > 1.0)):
            raise ValueError(
                f"IHe KED reference {p.name} {frac_name} outside [0, 1]: "
                f"{frac}."
            )

    dominant = np.asarray(frame["dominantError"].to_numpy(), dtype=str)
    unknown = set(dominant.tolist()) - _LEGAL_DOMINANT_ERROR
    if unknown:
        raise ValueError(
            f"IHe KED reference {p.name} has unknown dominantError tags "
            f"{sorted(unknown)}; legal: {sorted(_LEGAL_DOMINANT_ERROR)}."
        )

    noise = _finite("noiseLimited")
    if not np.all(np.isin(noise, (0.0, 1.0))):
        raise ValueError(
            f"IHe KED reference {p.name} noiseLimited must be 0 or 1; "
            f"got {noise}."
        )

    return IHeKedReference(
        n=n_int,
        label=np.asarray(frame["label"].to_numpy(), dtype=str),
        mass_center_u=_finite("massCenter_u_per_e"),
        N_counts=_finite("N_counts"),
        N_eff=_finite("N_eff"),
        mean_KE_eV=mean_KE_eV,
        mode_KE_eV=_finite("modeKE_eV"),
        median_KE_eV=_finite("medianKE_eV"),
        sigma_KE_eV=_finite("sigmaKE_eV"),
        stat_err_mean_KE_eV=_finite("statErr_meanKE_eV"),
        sys_err_mean_KE_eV=_finite("sysErr_meanKE_eV"),
        calib_syst_frac=_finite("calibSyst_frac"),
        condition_syst_frac=_finite("conditionSyst_frac"),
        bg_off_shift_eV=_finite("bgOffShift_eV"),
        dominant_error=dominant,
        noise_limited=noise.astype(bool),
        source_path=p.resolve(),
    )


_CURVE_COLUMNS: tuple[str, ...] = (
    "E_eV",
    "v_mps",
    "signal_2d_Pv",
    "signal_2d_PE",
    "signal_3d_Pv",
    "signal_3d_PE",
)


@dataclass(frozen=True)
class IHeKedCurve:
    """One trusted per-fragment curve file (n = 0..4).

    Both representations on one shared axis pair. Each signal column is
    independently peak-normalized (smoothed envelope = 1): compare shapes,
    never amplitudes. NaN marks excluded regions (3-D columns only:
    Abel-center spike / low-E cut); matplotlib renders them as gaps.

    Attributes
    ----------
    n : int
        He count of the fragment (0..4).
    E_eV : np.ndarray, shape (M,)
        Kinetic energy axis (mass-independent detector mapping).
    v_mps : np.ndarray, shape (M,)
        Speed axis of the mass-m(n) complex; differs per fragment at
        equal E.
    signal_2d_Pv, signal_2d_PE : np.ndarray, shape (M,)
        2-D detector projection per unit speed / per unit energy. The
        coordinate is the projected in-plane speed (sim counterpart:
        ``compute_final_velocity_histogram(..., projected=True)``).
    signal_3d_Pv, signal_3d_PE : np.ndarray, shape (M,)
        Reconstructed 3-D speed / energy distribution (sim counterpart:
        the 3-D ``|v|`` histogram).
    source_path : Path
        Resolved path of the loaded CSV.
    """

    n: int
    E_eV: np.ndarray
    v_mps: np.ndarray
    signal_2d_Pv: np.ndarray
    signal_2d_PE: np.ndarray
    signal_3d_Pv: np.ndarray
    signal_3d_PE: np.ndarray
    source_path: Path


def load_ihe_ked_curve(directory: str | Path, n: int) -> IHeKedCurve:
    """Load ``IHe_KED_curves_n{n}.csv`` for one fragment n in 0..4.

    Raises
    ------
    ValueError
        ``n`` outside 0..``CURVE_N_MAX``; missing/extra columns; a
        non-finite axis value; a non-ascending axis; or a signal column
        with no finite entries at all.
    FileNotFoundError
        If the curve file does not exist in ``directory``.
    """
    if not (0 <= int(n) <= CURVE_N_MAX):
        raise ValueError(
            f"n must be in 0..{CURVE_N_MAX} (trusted curves); got {n}."
        )
    p = Path(directory) / f"IHe_KED_curves_n{int(n)}.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"IHe KED curve file not found: {p.resolve()}"
        )

    frame = pd.read_csv(p)
    validate_columns(
        list(frame.columns),
        _CURVE_COLUMNS,
        file_label=f"IHe KED curve {p.name}",
    )

    def _axis(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"IHe KED curve {p.name} axis {name!r} has non-finite "
                f"entries."
            )
        if np.any(np.diff(values) <= 0.0):
            raise ValueError(
                f"IHe KED curve {p.name} axis {name!r} must be strictly "
                f"ascending."
            )
        return values

    def _signal(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.isfinite(values).any():
            raise ValueError(
                f"IHe KED curve {p.name} signal {name!r} has no finite "
                f"entries."
            )
        if np.any(np.isinf(values)):
            raise ValueError(
                f"IHe KED curve {p.name} signal {name!r} contains inf."
            )
        return values

    return IHeKedCurve(
        n=int(n),
        E_eV=_axis("E_eV"),
        v_mps=_axis("v_mps"),
        signal_2d_Pv=_signal("signal_2d_Pv"),
        signal_2d_PE=_signal("signal_2d_PE"),
        signal_3d_Pv=_signal("signal_3d_Pv"),
        signal_3d_PE=_signal("signal_3d_PE"),
        source_path=p.resolve(),
    )
