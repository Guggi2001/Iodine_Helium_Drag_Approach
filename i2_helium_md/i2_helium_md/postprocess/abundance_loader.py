"""Loader for the experimental I+He_n abundance reference (Tier-2 Phase E, Slice E3).

Reads ``data/reference/integrated_i_he_abundance.csv`` -- the Tier-2 arbiter:
the measured I+He_n integer-n size distribution the generative mechanism is
falsified against (scored by the E4 Wasserstein comparison).

The contract mirrors :func:`~i2_helium_md.postprocess.hedft_loader.load_hedft_trajectory`:
a frozen-dataclass return, an **order-independent membership** header match
(missing *and* extra columns raise), loud ``FileNotFoundError`` / ``ValueError``,
and ``source_path = p.resolve()``.

Data contract (verified against the real file, 2026-07-03)
----------------------------------------------------------
Columns ``n, label, massCenter_u_per_e, massWindowLower_u_per_e,
massWindowUpper_u_per_e, ionCounts, ionPercent``; 21 rows with ``n`` contiguous
``0..20`` (``n = 0`` = bare I+, labels ``I^+`` .. ``I^+He_20``); ``massCenter``
step 4.0026 u/e (He); ``ionPercent`` sums to 100. The reference's ``massCenter``
uses the integer I mass (127 u at ``n = 0``), distinct from the simulation-side
``complex_mass_amu`` convention (126.90 base) -- E3 loads the reference column
verbatim and never recomputes it, so the two conventions do not collide.

Provenance is **to be completed**: the file ships without a documented producer
(no exporter under ``data/reference/scripts/`` for it; the only in-repo
reference is the consumer ``data/reference/scripts/plotting_histogram.py``). The
verified data contract is recorded in ``data/reference/README.md``; the user
supplies measurement IDs / the producing script when available.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .csv_contract import validate_columns

# Order-independent: membership only. Both missing and extra columns raise.
_EXPECTED_COLUMNS: tuple[str, ...] = (
    "n",
    "label",
    "massCenter_u_per_e",
    "massWindowLower_u_per_e",
    "massWindowUpper_u_per_e",
    "ionCounts",
    "ionPercent",
)

# ionPercent must sum to 100 within this tolerance (the real file sums to
# exactly 100.0; a modest band absorbs any future rounding without masking a
# genuinely broken export).
_PERCENT_SUM: float = 100.0
_PERCENT_SUM_TOL: float = 1e-3


@dataclass(frozen=True)
class HeAbundanceReference:
    """Experimental I+He_n abundance reference on integer-n support.

    Attributes
    ----------
    n : np.ndarray, shape (Nn,), int
        Shell count, contiguous and ascending from 0.
    label : np.ndarray, shape (Nn,), str
        Species label (``I^+``, ``I^+He_1``, ...).
    mass_center_u : np.ndarray, shape (Nn,), float
        Mass-window centre in u/e (the reference's own mass grid).
    ion_counts : np.ndarray, shape (Nn,), float
        Per-species ion intensity (non-negative). Despite the column name,
        the real file's values are fractional (column sum ~ 0.19) -- a
        normalized intensity, **not** raw detector counts; the absolute
        normalization is part of the open provenance item.
    ion_fraction : np.ndarray, shape (Nn,), float
        ``ionPercent`` normalized by its own sum; sums to exactly 1.0. This is
        the distribution the E4 Wasserstein comparison scores against.
    source_path : Path
        Absolute path of the file that was loaded.
    """

    n: np.ndarray
    label: np.ndarray
    mass_center_u: np.ndarray
    ion_counts: np.ndarray
    ion_fraction: np.ndarray
    source_path: Path


def load_he_abundance_reference(path: str | Path) -> HeAbundanceReference:
    """Load the experimental I+He_n abundance CSV into a :class:`HeAbundanceReference`.

    Parameters
    ----------
    path
        Path to the CSV. Must exist; ``FileNotFoundError`` otherwise.

    Returns
    -------
    HeAbundanceReference
        Frozen dataclass with the integer support, labels, mass grid, counts,
        and the sum-normalized ``ion_fraction``.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the header is missing/has extra columns; any numeric column has a
        non-finite entry (blank/corrupted cell); ``n`` is non-integer,
        non-contiguous, or does not ascend from 0; ``ionCounts`` or
        ``ionPercent`` has a negative entry; or ``ionPercent`` does not sum
        to 100 within tolerance.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"I+He_n abundance reference file not found: {p.resolve()}"
        )

    # dtype=None + encoding lets the string ``label`` column coexist with the
    # numeric columns in one structured array; names come from the header.
    structured = np.atleast_1d(
        np.genfromtxt(p, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )

    validate_columns(
        structured.dtype.names,
        _EXPECTED_COLUMNS,
        file_label=f"abundance reference {p.name}",
    )

    def _finite_column(name: str) -> np.ndarray:
        """Load a numeric column, refusing NaN/inf (a blank or corrupted cell
        becomes NaN under genfromtxt, and NaN comparisons are silently False
        in every downstream guard)."""
        values = np.asarray(structured[name], dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"abundance reference {p.name} column {name!r} has non-finite "
                f"entries (NaN/inf -- blank or corrupted cell): {values}."
            )
        return values

    n_float = _finite_column("n")
    if n_float.size == 0:
        raise ValueError(f"abundance reference {p.name} is empty (no rows).")
    if np.any(n_float != np.rint(n_float)):
        raise ValueError(
            f"abundance reference {p.name} has non-integer n values: {n_float}."
        )
    n_int = np.rint(n_float).astype(int)
    if not np.array_equal(n_int, np.arange(n_int.size)):
        raise ValueError(
            f"abundance reference {p.name} n column must be contiguous and "
            f"ascending from 0; got {n_int.tolist()}."
        )

    ion_counts = _finite_column("ionCounts")
    if np.any(ion_counts < 0.0):
        raise ValueError(
            f"abundance reference {p.name} has negative ionCounts: "
            f"{ion_counts[ion_counts < 0.0]}."
        )

    ion_percent = _finite_column("ionPercent")
    if np.any(ion_percent < 0.0):
        raise ValueError(
            f"abundance reference {p.name} has negative ionPercent entries "
            f"(a negative abundance is unphysical even when the sum is 100): "
            f"{ion_percent[ion_percent < 0.0]}."
        )
    percent_sum = float(ion_percent.sum())
    if abs(percent_sum - _PERCENT_SUM) > _PERCENT_SUM_TOL:
        raise ValueError(
            f"abundance reference {p.name} ionPercent must sum to "
            f"{_PERCENT_SUM} (within {_PERCENT_SUM_TOL}); got {percent_sum}."
        )
    ion_fraction = ion_percent / ion_percent.sum()

    return HeAbundanceReference(
        n=n_int,
        label=np.asarray(structured["label"], dtype=str),
        mass_center_u=_finite_column("massCenter_u_per_e"),
        ion_counts=ion_counts,
        ion_fraction=ion_fraction,
        source_path=p.resolve(),
    )
