"""Score a pooled battery against the experimental size and energy references.

One row per pool: size distribution, mean KE per fragment size, joint score,
high-n tail and the law that produced it. The physics reads come from
:mod:`i2_helium_md.postprocess.tier2_confirmation`.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from ..presets import REFERENCE_DRAG_ROOT
from ..simulation.detection_stage import load_detection_result
from ..simulation.run_directory import RunDirectory
from .abundance_loader import load_he_abundance_reference
from .ihe_ked import load_ihe_ked_reference
from .tier2_confirmation import (
    deep_bin_ke_ratio,
    midhot_ratio,
    read_confirmation_detection,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
)

# The references are read from the source tree (editable install).
_REFERENCE_ROOT = REFERENCE_DRAG_ROOT.parent
ABUNDANCE_REFERENCE_CSV = (
    _REFERENCE_ROOT / "integrated_i_he_abundance.csv"
)
IHE_KED_REFERENCE_CSV = (
    _REFERENCE_ROOT / "ihe_ked" / "IHe_KED_reference.csv"
)

# KE-curve scoring knobs — the committed §I.11.4.1 / I77 conventions.
KE_N_MIN = 1
MIN_KE_BIN_COUNT = 2
INCLUDE_SIM_SE = True
N1_ANCHOR = "median"

# The design's provisional normalizers (plan §3.5e). Every S is printed under
# these AND under the measured set until the §1.2 deep-bin item is settled.
PROVISIONAL_NORM = {
    "w1_ref": 0.571,
    "midhot_ln": float(np.log(1.15)),
    "deepke_ln": float(np.log(1.15)),
}

# The finals S convention (Block-0-licensed axes, provisional norms).
LICENSED_AXES = ("w1", "midhot")

# Accepted n = 1 share (solvated) and mean fragment size.
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)


def joint_score(
    w1: float,
    midhot: float,
    deepke: float,
    norm: dict[str, float],
    axes: Sequence[str] = ("w1", "midhot", "deepke"),
) -> float:
    """The plan §3.5e joint score. Lower is better; 0 is a perfect cell.

    ``S = W1/w1_ref + |ln midHot|/midhot_ln + |ln deepKE|/deepke_ln``

    Equal weights are a **declared convention, not a truth**: each term is the
    observable's distance from the experiment divided by a scale that makes the
    three commensurable. ``axes`` drops any term Block 0 did not license, so an
    unlicensed observable never silently contributes to a ranking.

    Units: ``w1`` is a Wasserstein distance in He-number bins; ``midhot`` and
    ``deepke`` are dimensionless sim/ref KE ratios (1.0 = experiment). Returns
    NaN if any *scored* axis is NaN.
    """
    total = 0.0
    if "w1" in axes:
        total += w1 / norm["w1_ref"]
    if "midhot" in axes:
        total += abs(math.log(midhot)) / norm["midhot_ln"] if midhot > 0 \
            else float("nan")
    if "deepke" in axes:
        total += abs(math.log(deepke)) / norm["deepke_ln"] if deepke > 0 \
            else float("nan")
    if any(np.isnan(v) for v in (
        w1 if "w1" in axes else 0.0,
        midhot if "midhot" in axes else 1.0,
        deepke if "deepke" in axes else 1.0,
    )):
        return float("nan")
    return float(total)


def licensed_score(row: dict[str, Any]) -> float:
    """Joint score on the licensed axes, provisional norms (plan §3.5e/§3.5f)."""
    return joint_score(float(row["w1_solv"]), float(row["midHot"]),
                       float(row["deepKE"]), PROVISIONAL_NORM,
                       axes=LICENSED_AXES)


def observable_columns(read, abundance_ref, ked_ref) -> dict[str, Any]:
    """Size-distribution and KE observables of one reduced read (pure).

    Parameters
    ----------
    read : ConfirmationDetectionRead
        Reduced detection result.
    abundance_ref, ked_ref
        From ``load_he_abundance_reference`` and ``load_ihe_ked_reference``.

    Returns
    -------
    dict
        ``num_scored`` and ``ke_npts`` (counts); ``trap``, ``trap_bound``,
        ``trap_marg``, ``supp`` (shares); ``nbar_det`` (mean He number);
        ``n1_solv`` (share of n = 1 among n >= 1); ``w1_solv`` (Wasserstein
        distance [He atoms]); ``midHot``, ``deepKE`` (mean KE / reference) with
        their bin counts; ``chi2_med`` (KE curve).
    """
    hist = score_histogram_vs_reference(read, abundance_ref)
    ke = score_ke_curve_vs_reference(
        read, ked_ref, n_min=KE_N_MIN, min_count=MIN_KE_BIN_COUNT,
        include_sim_se=INCLUDE_SIM_SE, n1_anchor=N1_ANCHOR,
    )
    mid = midhot_ratio(read, ked_ref)
    deep = deep_bin_ke_ratio(read, ked_ref)
    return {
        "num_scored": read.num_scored,
        "trap": read.trapped_frac,
        # Decomposed on purpose (plan §3.5b item 2): trap_bound is the exact
        # conservative verdict (physics — the ion cannot escape at any
        # relaxation length), trap_marginal is a MODELLING exclusion under
        # `exclude_all_coupled`, conditional on the drag law and the window.
        # Never merge them in a table; `trap` is their sum.
        "trap_bound": read.trap_bound_frac,
        "trap_marg": read.trap_marginal_frac,
        "supp": read.suppressed_frac,
        "nbar_det": read.n_mean,
        "n1_solv": hist.sim_n1,
        "w1_solv": hist.w1_solvated,
        "midHot": mid.value,
        "midHot_bins": mid.n_bins,
        "deepKE": deep.value,
        "deepKE_bins": deep.n_bins,
        "chi2_med": ke.chi2_profiled,
        # Reported next to chi2 on purpose: the deep-birth cells vacate the
        # low-n bins, so their chi2 is a sum over far fewer terms and a smaller
        # value there is NOT a better KE landing.
        "ke_npts": ke.n_points,
    }


def lowke_columns(read) -> dict[str, Any]:
    """Mean, median and mode of the KE at n = 1 and 2 [eV]; NaN below 5 ions.

    Mode = midpoint of the tallest 40-bin histogram cell on [0, 3.5] eV
    (display diagnostic; the scored quantities are KE1_mean / KE2_mean).
    """
    out: dict[str, Any] = {}
    for n in (1, 2):
        sel = read.ke_scored_eV[read.n_scored == n]
        if sel.size < 5:
            out[f"KE{n}_mean"] = float("nan")
            out[f"KE{n}_med"] = float("nan")
            out[f"KE{n}_mode"] = float("nan")
            continue
        hist, edges = np.histogram(sel, bins=40, range=(0.0, 3.5))
        i = int(np.argmax(hist))
        out[f"KE{n}_mean"] = float(sel.mean())
        out[f"KE{n}_med"] = float(np.median(sel))
        out[f"KE{n}_mode"] = float(0.5 * (edges[i] + edges[i + 1]))
    return out


def score_row(label: str, member: str, pipeline: str, seed, read, ab, ked) -> dict[str, Any]:
    """One score row: ids, observables, low-n KE [eV], gate and joint score.

    Parameters
    ----------
    label, member, pipeline, seed
        Identification columns, passed through.
    read, ab, ked
        As for :func:`observable_columns`.

    Returns
    -------
    dict
        Columns of :func:`observable_columns` and :func:`lowke_columns`, plus
        ``md_gate`` (1 if n = 1 share and mean n lie in their bands) and ``S``
        (:func:`licensed_score`, lower is better).
    """
    row: dict[str, Any] = {"label": label, "member": member,
                           "pipeline": pipeline, "seed": seed}
    row.update(observable_columns(read, ab, ked))
    row.update(lowke_columns(read))
    row["md_gate"] = int(GATE_N1[0] <= float(row["n1_solv"]) <= GATE_N1[1]
                         and GATE_NBAR[0] <= float(row["nbar_det"]) <= GATE_NBAR[1])
    row["S"] = licensed_score(row)
    return row


def tail_columns(read, ked) -> dict[str, float]:
    """High-n tail and deep-bin KE of one reduced read.

    Parameters
    ----------
    read : ConfirmationDetectionRead
        Reads ``n_scored`` [He atoms] and ``ke_scored_eV`` [eV].
    ked
        KE reference with ``n`` and ``mean_KE_eV`` [eV].

    Returns
    -------
    dict
        ``mass_n15plus_solv``: share of n >= 15 among n >= 1 (NaN if none);
        ``n_max``; ``KE10_14_over_ref``: mean KE / reference over n = 10-14,
        bins with more than 3 ions (NaN if none).

    Raises
    ------
    ValueError
        If the read holds no scored ions.
    """
    ns, ke = np.asarray(read.n_scored), np.asarray(read.ke_scored_eV)
    if ns.size == 0:
        raise ValueError("no scored ions: the high-n tail is undefined")
    solv = ns >= 1
    ref = dict(zip(np.asarray(ked.n).astype(int), np.asarray(ked.mean_KE_eV, dtype=float)))
    ratios = [ke[ns == k].mean() / ref[k] for k in range(10, 15) if (ns == k).sum() > 3]
    return {"mass_n15plus_solv": float((ns[solv] >= 15).mean()) if solv.any() else float("nan"),
            "n_max": int(ns.max()),
            "KE10_14_over_ref": float(np.mean(ratios)) if ratios else float("nan")}


def law_columns(cfg) -> dict[str, Any]:
    """Drag law, well, birth-energy and cooling parameters of a drag config.

    Parameters
    ----------
    cfg : SimConfig
        Must carry ``drag_coefficients``.

    Returns
    -------
    dict
        ``drag_form``, ``drag_coefficients`` (law units: a [amu/ps], b [amu ps/A^2],
        c [amu/A], v_c [A/ps]), ``e_bind_eV``, ``E0_eV``, ``tau_ps``, ``sampler``,
        ``gate_A`` [Angstrom].
    """
    co = cfg.drag_coefficients.coefficients
    return {"drag_form": cfg.drag_form, "drag_coefficients": dict(co),
            "e_bind_eV": float(cfg.binding_energy_I_ion_eV),
            "E0_eV": float(cfg.internal_energy_partition_fraction * cfg.coulomb_available_eV),
            "tau_ps": float(cfg.internal_energy_cooling_tau_ps),
            "sampler": cfg.droplet_size_sampler_mode,
            "gate_A": float(cfg.drag_gate_steepness)}


def score_pool(pool_dir: Path, label: str) -> dict[str, float | int | str]:
    """Score one pooled battery directory (``cfg.json`` + ``detection.npz``).

    Parameters
    ----------
    pool_dir : Path
        Pool directory, see :func:`i2_helium_md.simulation.battery.pool_battery`.
    label : str
        Row label, e.g. the law name.

    Returns
    -------
    dict
        Columns of ``g13_decoupled_table.csv``. Energies in eV, times in ps,
        widths in Angstrom; ``KE*``, ``midHot`` and ``deepKE`` are NaN when too
        few ions populate their bins.

    Raises
    ------
    FileNotFoundError
        If a file of the pool or a reference is missing.
    """
    pool_dir = Path(pool_dir)
    cfg = RunDirectory(pool_dir).load_cfg()
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    read = read_confirmation_detection(
        load_detection_result(pool_dir / "detection.npz"), label=label)
    return {
        **score_row(label, "pooled", cfg.droplet_size_sampler_mode, "-", read, ab, ked),
        **tail_columns(read, ked),
        **law_columns(cfg),
        "num_molecules": int(cfg.num_molecules),
        "container": pool_dir.name,
        "ihe_ked_reference": IHE_KED_REFERENCE_CSV.name,
        "self_bound_gate_release": cfg.self_bound_gate_release,
    }


def _fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Fixed-width text table."""
    if not rows:
        return "(no rows)"
    cols = list(rows[0].keys())
    cells = [[_fmt(r[c]) for c in cols] for r in rows]
    widths = [
        max(len(c), *(len(row[i]) for row in cells))
        for i, c in enumerate(cols)
    ]
    lines = [
        "  ".join(c.ljust(w) for c, w in zip(cols, widths)),
        "  ".join("-" * w for w in widths),
    ]
    lines += ["  ".join(v.ljust(w) for v, w in zip(row, widths)) for row in cells]
    return "\n".join(lines)


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write the rows as CSV; returns the resolved path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return p.resolve()
