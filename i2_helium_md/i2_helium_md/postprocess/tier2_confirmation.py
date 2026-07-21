"""Tier-2 confirmation-matrix scorer (§I.11.4 Stage 0 — the T4-absorbed build).

Repo-ization of the T9 leg-scoring conventions (findings §4m–§4r) plus the
experimental scoring layer against the two committed references. Scope:
pure post-processing of finished run dirs — no physics, no config surface.

Conventions (frozen as scored throughout the §I.11 oracle chain):

- ``droplet_retained`` ions are **excluded from every read** (they carry the
  verbatim in-droplet handover state, not a detector arrival — the
  ``DetectionResult.detected_mask`` contract, review fix 2026-07-18);
- ``suppressed`` ions score in **bin 0** (the RQ3 bare-candidate class —
  the never-opened side of the Δ× race, counted not evolved);
- histograms live on the integer support ``0..n_max`` (default ``N_STAR``);
- W1 is the CDF-gap sum on unit rung spacing, computed by
  :func:`~i2_helium_md.postprocess.distribution_compare.wasserstein_integer_support`
  (no second metric implementation);
- the experimental solvated read renormalizes to **n >= 1** (RQ8: the bare
  bin is a channel-branching quantity, not a model target);
- the KE-axis score implements the **committed ihe_ked error model**
  (``data/reference/ihe_ked/COLUMNS.md``): per-point
  ``sqrt(statErr^2 + sysErr^2)`` scatter (optionally widened by the sim
  bin SE), plus the calib and condition fractional bands as two
  *correlated* whole-curve shifts profiled as unit-normal nuisances
  (closed-form 2x2 linear profile).

Twin comparison inputs are the committed pre-registration CSVs written by
``scripts/tier2_h2b_forward_model.py`` (``h2b_leg_*_predictions.csv`` /
``h2b_leg_*_ke.csv``).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ..simulation.detection_stage import DetectionResult, load_detection_result
from .abundance_loader import HeAbundanceReference
from .distribution_compare import wasserstein_integer_support
from .ihe_ked import IHeKedReference
from .size_distribution import N_STAR

__all__ = [
    "ConfirmationDetectionRead",
    "KEPerBin",
    "TwinPredictionRow",
    "TwinKECurve",
    "HistogramScore",
    "KECurveScore",
    "read_confirmation_detection",
    "load_confirmation_run",
    "load_twin_prediction",
    "load_twin_ke_curve",
    "wasserstein_between",
    "solvated_renormalized",
    "score_histogram_vs_reference",
    "score_ke_curve_vs_reference",
]

# The retained class marker (detection-stage vocabulary).
_RETAINED = "droplet_retained"
_SUPPRESSED = "suppressed"

# Twin prediction CSVs are written at 4-decimal rounding over 22 bins, so a
# committed row's histogram may sum to 1 +- ~2e-3. Accept within this band
# and renormalize exactly; anything larger is a wrong/corrupted row.
_TWIN_H_SUM_TOL = 5e-3

_DETECTION_FILENAME = "detection.npz"


# ---------------------------------------------------------------------------
# per-run detection read
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class KEPerBin:
    """Per-shell-bin detected kinetic energy of one scored run.

    Attributes
    ----------
    n : np.ndarray, shape (Nb,), int
        Occupied scored bins, ascending (suppressed ions occupy bin 0).
    mean_eV : np.ndarray, shape (Nb,)
        Mean detected kinetic energy per bin [eV].
    se_eV : np.ndarray, shape (Nb,)
        Standard error of the bin mean (ddof=1; 0.0 for single-member bins —
        the ``fragment_mean_kinetic_energy`` precedent) [eV].
    count : np.ndarray, shape (Nb,), int
        Ions per bin.
    """

    n: np.ndarray
    mean_eV: np.ndarray
    se_eV: np.ndarray
    count: np.ndarray


@dataclass(frozen=True)
class ConfirmationDetectionRead:
    """One run's detection artifact reduced under the §4r scoring conventions.

    Attributes
    ----------
    label : str
        Caller-chosen run label (e.g. ``"dc1"``).
    num_ions : int
        Total ions in the artifact (``2N``).
    num_scored : int
        Ions entering every read (``droplet_retained`` excluded).
    trapped_frac : float
        ``droplet_retained`` fraction of **all** ions.
    suppressed_frac : float
        ``suppressed`` fraction of the **scored** ions.
    n_scored : np.ndarray, shape (num_scored,)
        Scored terminal shell counts (int-valued; suppressed mapped to 0).
    ke_scored_eV : np.ndarray, shape (num_scored,)
        Detected kinetic energies of the scored ions [eV].
    n_values : np.ndarray, shape (n_max+1,), int
        Histogram support ``0..n_max``.
    fraction : np.ndarray, shape (n_max+1,)
        Scored-ensemble fractions (sum exactly 1).
    n_mean : float
        Mean scored ``n`` (suppressed counted at 0).
    n1_frac : float
        ``fraction[1]``.
    """

    label: str
    num_ions: int
    num_scored: int
    trapped_frac: float
    suppressed_frac: float
    n_scored: np.ndarray
    ke_scored_eV: np.ndarray
    n_values: np.ndarray
    fraction: np.ndarray
    n_mean: float
    n1_frac: float

    def ke_by_n(self, *, min_count: int = 1) -> KEPerBin:
        """Per-bin mean detected KE over the scored ions.

        Parameters
        ----------
        min_count
            Drop bins with fewer ions than this (default 1 = keep all
            occupied bins).
        """
        ns: list[int] = []
        means: list[float] = []
        ses: list[float] = []
        counts: list[int] = []
        for n in np.unique(self.n_scored):
            sel = self.ke_scored_eV[self.n_scored == n]
            if sel.size < min_count:
                continue
            ns.append(int(n))
            means.append(float(sel.mean()))
            ses.append(
                float(sel.std(ddof=1) / np.sqrt(sel.size)) if sel.size > 1
                else 0.0
            )
            counts.append(int(sel.size))
        return KEPerBin(
            n=np.asarray(ns, dtype=int),
            mean_eV=np.asarray(means, dtype=float),
            se_eV=np.asarray(ses, dtype=float),
            count=np.asarray(counts, dtype=int),
        )


def read_confirmation_detection(
    detection: DetectionResult,
    *,
    label: str = "",
    n_max: int = N_STAR,
) -> ConfirmationDetectionRead:
    """Reduce one ``DetectionResult`` under the §4r scoring conventions.

    ``droplet_retained`` ions are excluded from every read; ``suppressed``
    ions score in bin 0; the histogram lives on ``0..n_max``.

    Raises
    ------
    ValueError
        If the exclusion empties the ensemble, any scored ``n`` is
        non-integer, or any scored ``n`` exceeds ``n_max``.
    """
    state = np.asarray(detection.state_reason)
    n_det = np.asarray(detection.n_detected, dtype=float)
    ke = np.asarray(detection.E_kin_detected_eV, dtype=float)

    num_ions = int(n_det.size)
    keep = state != _RETAINED
    if not np.any(keep):
        raise ValueError(
            f"run {label!r}: every ion is {_RETAINED} — the exclude policy "
            f"leaves no scored ensemble."
        )

    n_scored = n_det[keep].copy()
    ke_scored = ke[keep].copy()
    state_scored = state[keep]

    if not np.allclose(n_scored, np.round(n_scored)):
        raise ValueError(
            f"run {label!r}: n_detected must be int-valued (v7 convention); "
            f"got non-integer entries."
        )
    n_scored = np.round(n_scored).astype(int)
    if np.any(n_scored < 0) or np.any(n_scored > n_max):
        raise ValueError(
            f"run {label!r}: scored n outside [0, n_max={n_max}]: "
            f"[{n_scored.min()}, {n_scored.max()}]."
        )

    # the §4r convention: the never-opened (suppressed) class is the RQ3
    # bare-candidate — it scores in bin 0, not at its handover n.
    suppressed = state_scored == _SUPPRESSED
    n_scored[suppressed] = 0

    counts = np.bincount(n_scored, minlength=n_max + 1).astype(float)
    fraction = counts / counts.sum()

    return ConfirmationDetectionRead(
        label=label,
        num_ions=num_ions,
        num_scored=int(n_scored.size),
        trapped_frac=float(np.count_nonzero(~keep) / num_ions),
        suppressed_frac=float(np.count_nonzero(suppressed) / n_scored.size),
        n_scored=n_scored,
        ke_scored_eV=ke_scored,
        n_values=np.arange(n_max + 1),
        fraction=fraction,
        n_mean=float(n_scored.mean()),
        n1_frac=float(fraction[1]),
    )


def load_confirmation_run(
    run_dir: str | Path,
    *,
    label: str | None = None,
    n_max: int = N_STAR,
) -> ConfirmationDetectionRead:
    """Load ``detection.npz`` from a finished run dir and reduce it.

    Parameters
    ----------
    run_dir
        Run directory holding a ``detection.npz`` artifact.
    label
        Row label; defaults to the directory name.
    """
    p = Path(run_dir)
    detection_path = p / _DETECTION_FILENAME
    if not detection_path.exists():
        raise FileNotFoundError(
            f"no {_DETECTION_FILENAME} in {p.resolve()} — the confirmation "
            f"scorer reads the detected ensemble only."
        )
    detection = load_detection_result(detection_path)
    return read_confirmation_detection(
        detection, label=p.name if label is None else label, n_max=n_max
    )


# ---------------------------------------------------------------------------
# twin pre-registration CSV loaders
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TwinPredictionRow:
    """One (leg, config) row of a committed twin prediction CSV.

    ``h`` is the twin's detected histogram on ``0..21`` (trapped excluded,
    suppressed in bin 0), renormalized exactly to 1 (the CSV is written at
    4-decimal rounding).
    """

    leg: str
    config: str
    h: np.ndarray
    nbar_det: float
    suppressed_frac: float
    trapped_frac: float
    bare_frac: float
    source_path: Path


@dataclass(frozen=True)
class TwinKECurve:
    """Per-bin mean detected KE of one (leg, config) twin row."""

    leg: str
    config: str
    n: np.ndarray
    weight: np.ndarray
    mean_KE_eV: np.ndarray
    source_path: Path

    def mean_ke_for(self, n: int) -> float | None:
        """The twin's mean KE in bin ``n`` [eV], or None if absent."""
        hit = np.nonzero(self.n == n)[0]
        return float(self.mean_KE_eV[hit[0]]) if hit.size else None


def _select_rows(path: Path, leg: str, config: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [
            row for row in csv.DictReader(fh)
            if row["leg"] == leg and row["config"] == config
        ]
    if not rows:
        raise ValueError(
            f"no (leg={leg!r}, config={config!r}) row in {path.resolve()}."
        )
    return rows


def load_twin_prediction(
    path: str | Path, *, leg: str, config: str, n_max: int = N_STAR
) -> TwinPredictionRow:
    """Load one (leg, config) histogram row from a twin predictions CSV."""
    p = Path(path)
    rows = _select_rows(p, leg, config)
    if len(rows) > 1:
        raise ValueError(
            f"{len(rows)} rows match (leg={leg!r}, config={config!r}) in "
            f"{p.resolve()} — the pre-registration contract is one row."
        )
    row = rows[0]
    h = np.asarray(
        [float(row[f"h{k}"]) for k in range(n_max + 1)], dtype=float
    )
    if np.any(h < 0.0) or not np.all(np.isfinite(h)):
        raise ValueError(
            f"twin histogram has negative/non-finite entries in {p.resolve()}."
        )
    total = float(h.sum())
    if abs(total - 1.0) > _TWIN_H_SUM_TOL:
        raise ValueError(
            f"twin histogram sum {total} is beyond the rounding tolerance "
            f"{_TWIN_H_SUM_TOL} of 1 in {p.resolve()} — wrong row or "
            f"corrupted export."
        )
    return TwinPredictionRow(
        leg=leg,
        config=config,
        h=h / total,
        nbar_det=float(row["nbar_det"]),
        suppressed_frac=float(row["suppressed_frac"]),
        trapped_frac=float(row["trapped_frac"]),
        bare_frac=float(row["bare_frac"]),
        source_path=p,
    )


def load_twin_ke_curve(
    path: str | Path, *, leg: str, config: str
) -> TwinKECurve:
    """Load one (leg, config) per-bin KE curve from a twin KE CSV."""
    p = Path(path)
    rows = _select_rows(p, leg, config)
    n = np.asarray([int(r["n"]) for r in rows], dtype=int)
    order = np.argsort(n)
    if np.unique(n).size != n.size:
        raise ValueError(
            f"duplicate n rows for (leg={leg!r}, config={config!r}) in "
            f"{p.resolve()}."
        )
    return TwinKECurve(
        leg=leg,
        config=config,
        n=n[order],
        weight=np.asarray([float(r["weight"]) for r in rows], dtype=float)[order],
        mean_KE_eV=np.asarray(
            [float(r["mean_KE_eV"]) for r in rows], dtype=float
        )[order],
        source_path=p,
    )


# ---------------------------------------------------------------------------
# distribution comparison helpers
# ---------------------------------------------------------------------------

def wasserstein_between(
    n_a: np.ndarray,
    f_a: np.ndarray,
    n_b: np.ndarray,
    f_b: np.ndarray,
) -> float:
    """W1 between two integer-support fraction vectors.

    Thin adapter onto :func:`wasserstein_integer_support` (the single W1
    implementation — rule 1); both fraction vectors must sum to 1.
    """
    return wasserstein_integer_support(
        SimpleNamespace(n_values=np.asarray(n_a), fraction=np.asarray(f_a)),
        SimpleNamespace(n=np.asarray(n_b), ion_fraction=np.asarray(f_b)),
    )


def solvated_renormalized(
    n_values: np.ndarray,
    fraction: np.ndarray,
    *,
    n_min: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Drop bins below ``n_min`` and renormalize (the RQ8 solvated read).

    Raises
    ------
    ValueError
        If no probability mass survives at ``n >= n_min``.
    """
    n_arr = np.asarray(n_values)
    f_arr = np.asarray(fraction, dtype=float)
    keep = n_arr >= n_min
    mass = float(f_arr[keep].sum())
    if mass <= 0.0:
        raise ValueError(
            f"no solvated mass at n >= {n_min}: the whole distribution sits "
            f"in the bare/suppressed bin."
        )
    return n_arr[keep], f_arr[keep] / mass


@dataclass(frozen=True)
class HistogramScore:
    """Solvated-branch histogram score of one run vs the abundance reference.

    All quantities are on the ``n >= n_min`` renormalized branch (RQ8).
    ``ratio_n1_over_n2`` is ``inf`` when ``sim_n2 == 0 < sim_n1`` and
    ``nan`` when both vanish.
    """

    w1_solvated: float
    sim_n1: float
    sim_n2: float
    ref_n1: float
    ref_n2: float
    ratio_n1_over_n2: float
    n_min: int


def score_histogram_vs_reference(
    read: ConfirmationDetectionRead,
    ref: HeAbundanceReference,
    *,
    n_min: int = 1,
) -> HistogramScore:
    """Score one run's solvated histogram against the abundance reference."""
    sim_n, sim_f = solvated_renormalized(
        read.n_values, read.fraction, n_min=n_min
    )
    ref_n, ref_f = solvated_renormalized(ref.n, ref.ion_fraction, n_min=n_min)

    def _frac(n_arr: np.ndarray, f_arr: np.ndarray, n: int) -> float:
        hit = np.nonzero(n_arr == n)[0]
        return float(f_arr[hit[0]]) if hit.size else 0.0

    sim_n1 = _frac(sim_n, sim_f, 1)
    sim_n2 = _frac(sim_n, sim_f, 2)
    if sim_n2 > 0.0:
        ratio = sim_n1 / sim_n2
    else:
        ratio = float("inf") if sim_n1 > 0.0 else float("nan")
    return HistogramScore(
        w1_solvated=wasserstein_between(sim_n, sim_f, ref_n, ref_f),
        sim_n1=sim_n1,
        sim_n2=sim_n2,
        ref_n1=_frac(ref_n, ref_f, 1),
        ref_n2=_frac(ref_n, ref_f, 2),
        ratio_n1_over_n2=ratio,
        n_min=n_min,
    )


# ---------------------------------------------------------------------------
# KE-curve score — the committed error model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class KECurveScore:
    """Per-n mean-KE score under the committed ihe_ked error model.

    The two correlated bands enter as unit-normal nuisances ``(a, b)``
    scaling the reference curve coherently
    (``r -> r * (1 + a*calib_frac + b*condition_frac)``); the reported
    ``chi2_profiled`` is minimized over ``(a, b)`` with the ``a^2 + b^2``
    prior penalty included. ``n_within_1p25`` is the coarse fallback bar
    (bins with sim/ref ratio within x1.25, both directions).

    ``n1_anchor`` / ``n1_widen_bg`` record the n = 1 comparison convention
    in force (I77 provenance): under ``"median"`` the ``ref_mean_eV`` entry
    at n = 1 holds the substituted ``median_KE_eV`` value.
    """

    n: np.ndarray
    sim_mean_eV: np.ndarray
    sim_se_eV: np.ndarray
    sim_count: np.ndarray
    ref_mean_eV: np.ndarray
    sigma_eV: np.ndarray
    z_profiled: np.ndarray
    chi2_profiled: float
    chi2_unprofiled: float
    a_calib: float
    b_condition: float
    n_points: int
    n_within_1p25: int
    include_sim_se: bool
    n1_anchor: str
    n1_widen_bg: bool


def score_ke_curve_vs_reference(
    read: ConfirmationDetectionRead,
    ref: IHeKedReference,
    *,
    n_min: int = 1,
    n_max: int | None = None,
    min_count: int = 1,
    include_sim_se: bool = True,
    n1_anchor: str = "median",
    n1_widen_bg: bool = False,
) -> KECurveScore:
    """Score one run's per-n mean detected KE against the ihe_ked reference.

    Parameters
    ----------
    read
        The reduced detection read (KE per scored bin).
    ref
        The committed reference (``load_ihe_ked_reference``).
    n_min, n_max
        Scored bin range; defaults to the solvated branch ``n >= 1`` (RQ8)
        up to the reference's own support.
    min_count
        Minimum ions per sim bin (thin bins carry no stable mean).
    include_sim_se
        Widen the per-point sigma by the sim bin SE in quadrature (the sim
        mean at N = 50 carries its own statistical error; reported so the
        choice is visible).
    n1_anchor
        The n = 1 comparison convention (I77: the experimental n = 1 mean
        KE is a two-population mixture mean; the solvated core the drag
        model produces sits at the median). One of:

        - ``"mean"`` — legacy behaviour, byte-identical to the pre-I77
          scorer (the §4s/§4t recorded-χ² regression anchor);
        - ``"median"`` (default) — the n = 1 reference comparison value is
          ``ref.median_KE_eV`` [eV] instead of ``ref.mean_KE_eV``; all
          other bins, the per-point errors, and the two correlated bands
          are unchanged;
        - ``"exclude"`` — drop the n = 1 bin from the fit (treated like
          the bare n = 0 bin; ``n_points`` shrinks by one).
    n1_widen_bg
        Robustness sub-variant (default off; only legal with
        ``n1_anchor="median"``): additionally fold the n = 1
        ``bg_off_shift_eV`` [eV] into that bin's per-point sigma in
        quadrature (a symmetric Gaussian stand-in for the one-sided
        background-off systematic; §4u showed median+bg ≈ exclude).

    Raises
    ------
    ValueError
        If no sim bin overlaps the reference in the scored range, on an
        unknown ``n1_anchor``, or if ``n1_widen_bg`` is combined with an
        anchor other than ``"median"``.
    """
    if n1_anchor not in ("mean", "median", "exclude"):
        raise ValueError(
            f"unknown n1_anchor {n1_anchor!r}; "
            f"expected 'mean', 'median', or 'exclude'."
        )
    if n1_widen_bg and n1_anchor != "median":
        raise ValueError(
            f"n1_widen_bg=True requires n1_anchor='median'; "
            f"got n1_anchor={n1_anchor!r}."
        )

    ke = read.ke_by_n(min_count=min_count)
    hi = int(ref.n.max()) if n_max is None else n_max

    ref_index: dict[int, int] = {int(v): i for i, v in enumerate(ref.n)}
    keep = [
        i for i, n in enumerate(ke.n)
        if n_min <= int(n) <= hi and int(n) in ref_index
        and not (n1_anchor == "exclude" and int(n) == 1)
    ]
    if not keep:
        raise ValueError(
            f"run {read.label!r}: no scorable KE bins in [{n_min}, {hi}] "
            f"overlap the reference support."
        )

    n_sel = ke.n[keep]
    s = ke.mean_eV[keep]
    se = ke.se_eV[keep]
    count = ke.count[keep]
    ridx = np.asarray([ref_index[int(n)] for n in n_sel], dtype=int)
    r = ref.mean_KE_eV[ridx].copy()
    point_err = ref.point_err_eV[ridx]
    c_frac = ref.calib_syst_frac[ridx]
    d_frac = ref.condition_syst_frac[ridx]

    n1_hits = np.nonzero(n_sel == 1)[0]
    if n1_anchor == "median" and n1_hits.size:
        r[n1_hits] = ref.median_KE_eV[ridx[n1_hits]]

    sigma_sq = point_err**2 + (se**2 if include_sim_se else 0.0)
    if n1_widen_bg and n1_hits.size:
        sigma_sq[n1_hits] += ref.bg_off_shift_eV[ridx[n1_hits]] ** 2
    if np.any(sigma_sq <= 0.0):
        raise ValueError(
            "zero per-point sigma — the reference carries no error at a "
            "scored bin; cannot weight the chi^2."
        )
    w = 1.0 / sigma_sq

    # profile the two correlated bands: rho = s - r, basis u1 = c*r, u2 = d*r,
    # minimize sum w*(rho - a*u1 - b*u2)^2 + a^2 + b^2  (closed-form 2x2).
    rho = s - r
    basis = np.stack([c_frac * r, d_frac * r], axis=1)  # (Np, 2)
    lhs = basis.T @ (w[:, None] * basis) + np.eye(2)
    rhs = basis.T @ (w * rho)
    shifts = np.linalg.solve(lhs, rhs)
    resid = rho - basis @ shifts
    chi2_prof = float(np.sum(w * resid**2) + np.sum(shifts**2))
    ratio = s / r

    return KECurveScore(
        n=n_sel,
        sim_mean_eV=s,
        sim_se_eV=se,
        sim_count=count,
        ref_mean_eV=r,
        sigma_eV=np.sqrt(sigma_sq),
        z_profiled=resid / np.sqrt(sigma_sq),
        chi2_profiled=chi2_prof,
        chi2_unprofiled=float(np.sum(w * rho**2)),
        a_calib=float(shifts[0]),
        b_condition=float(shifts[1]),
        n_points=int(n_sel.size),
        n_within_1p25=int(
            np.count_nonzero(np.maximum(ratio, 1.0 / ratio) <= 1.25)
        ),
        include_sim_se=include_sim_se,
        n1_anchor=n1_anchor,
        n1_widen_bg=n1_widen_bg,
    )
