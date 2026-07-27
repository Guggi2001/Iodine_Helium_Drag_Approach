"""Atlas D2b §4.3 — geometry-grid re-weighting: oracle + corrected-ensemble forecast (zero MD).

The 11-cell Axis A grid is the measured transfer function
{(R, birth law) -> observable vector}. This report re-weights its cells with a
candidate droplet-size density to predict that ensemble's observables without
new MD (plan §4.3). Two products, in order:

1. **Re-weighting oracle (mandatory, §4.3).** Reconstruct the pooled N = 5000
   standing battery from the grid's L3 column (the production birth law)
   weighted by the battery's own *realized* droplet radii. If a known ensemble
   cannot be reconstructed, no reconstruction of an unknown one is admissible.
   The oracle also FIXES the interpolation convention.
2. **Corrected-geometry forecast (the G2 input).** The corrected ensemble
   (`legacy` + `raw` sizes at the preset's own source conditions, Boltzmann
   births) re-weighted onto the L2 column, plus the 2 x 2 decomposition
   counterfactuals (standing/corrected sizes x L3/L2 birth law).

Method conventions (frozen before the first run of this script)
---------------------------------------------------------------
- **Column-matched 1-D re-weighting.** The candidate laws coincide with grid
  columns exactly (production `uniform_volume` m3 = L3; corrected Boltzmann
  313.2 K = L2), so the conditional birth-depth distribution at each R is the
  cell's own by construction and only the R marginal is re-weighted. A 2-D
  (R, depth) interpolation across columns is deliberately NOT attempted: the
  L1/L2 columns are near-degenerate in mean depth but differ in spread, which
  makes depth-based interpolation ill-posed (plan §3.1).
- **Two candidate conventions**, adjudicated by the oracle: ``nearest``
  (midpoint binning in R) and ``linear`` (piecewise-linear hat weights in R).
  Mass outside the support is clamped to the end cells and its fraction is
  reported (no silent caps).
- **Pre-registered admissibility gate** (frozen before results were seen): per
  oracle observable, |reconstruction - recorded| <= max(2 x seed-SD,
  10 % x |recorded|), with the seed-SDs measured from the five N = 1000
  battery members. The convention passing more of the 7 columns is adopted
  (tie -> smaller total |error|/tolerance); the adopted convention must pass
  >= 6 of 7, else the re-weighting method is INADMISSIBLE on this support and
  the forecast below is printed as indicative only. The adopted convention's
  per-observable errors are the *stated interpolation error* that travels
  with every forecast row.
- **Mixture scoring reuses the committed scorer** (rule 1): the mixture is a
  duck-typed read whose sufficient statistics (fraction vector, per-bin KE
  means, fate fractions) are exact weighted combinations of the per-cell
  reads; ``score_histogram_vs_reference`` / ``midhot_ratio`` /
  ``deep_bin_ke_ratio`` consume it unchanged. ``chi2_med`` is deliberately
  ABSENT from mixture rows: its committed error model widens per-point sigma
  by the sim bin SE, an ion-count convention that has no exact weighted
  analogue — omitting it is honest, approximating it silently is not.
- **In-support oracle (declared POST-HOC).** Added after the pre-registered
  oracle returned INADMISSIBLE with 53.8 % of the standing density clamped
  below the R = 26.6 Å support edge: reconstruct the battery's own
  R >= support-min sub-ensemble (a known target fully inside the support)
  from the grid cells. This separates *method* error (column re-weighting on
  covered support) from *support* error (the clamp). It does not overturn the
  pre-registered verdict on the full standing reconstruction; it qualifies
  what that verdict means for the corrected forecast, whose density the
  support covers to ~95 %.
- **Arm B range rows.** Cells with a non-empty ``droplet_retained_marginal``
  class carry the §3.5b NOT-TIGHT bracket, so every mixture touching them is
  reported twice: Arm A (marginals excluded — headline) and Arm B (marginals
  injected at their exact conservative asymptotic KE). The Arm-B construction
  is imported from ``tier2atlas_retained_bracket.py`` (its one sanctioned
  consumer); the marginal fraction is conditional on a cubic law extrapolated
  ~3x beyond its calibration band (§3.5b item 6) — the caveat travels with
  every Arm-B row.

Invocation: edit USER SETTINGS, then

    python scripts/post_processing/tier2atlas_geometry_reweight.py
"""

from __future__ import annotations

import dataclasses
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.physics.constants import (  # noqa: E402
    droplet_radius_bulk_angstrom,
)
from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    ConfirmationDetectionRead,
    KEPerBin,
    deep_bin_ke_ratio,
    load_confirmation_run,
    midhot_ratio,
    read_confirmation_detection,
    score_histogram_vs_reference,
)
from i2_helium_md.sampling.droplet_sizes import (  # noqa: E402
    mean_droplet_size,
    sample_droplet_sizes,
)
from i2_helium_md.simulation.checkpoint import (  # noqa: E402
    load_neutral_checkpoint,
)
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    DetectionResult,
    load_detection_result,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    CELL_LABELS,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RECORDED,
    ORACLE_RUN,
    RUN_SELECTION,
    RUNS_ROOT,
    check_oracle,
    observable_row,
)
from scripts.post_processing.tier2atlas_retained_bracket import (  # noqa: E402
    arm_b_detection,
    marginal_dossier,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# The five standing-battery members: seed-SD yardstick + the realized standing
# droplet-radius density (their neutral checkpoints' sampled radii).
BATTERY_MEMBERS: tuple[str, ...] = tuple(
    f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{i}"
    for i in range(1, 6)
)

# Corrected-ensemble draw: `legacy` + `raw` at the preset's own source
# conditions (read from a grid cell's cfg.json — no number is typed here).
CORRECTED_DENSITY_SOURCE_CELL = "r3l2"
CORRECTED_SAMPLE_N = 200_000
CORRECTED_SAMPLE_SEED = 20260727          # the grid seed; a density draw only

# Pre-registered admissibility gate (see module docstring).
GATE_SD_MULT = 2.0
GATE_REL_TOL = 0.10
GATE_MIN_PASS = 6                          # of the 7 oracle columns
GATE_COLUMNS: tuple[str, ...] = (
    "trap", "supp", "nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE",
)

CONVENTIONS: tuple[str, ...] = ("nearest", "linear")

# ---------------------------------------------------------------------------
# interpolation weights (pure; tested)
# ---------------------------------------------------------------------------


def interpolation_weights(
    radii_A: np.ndarray,
    support_A: np.ndarray,
    convention: str,
) -> tuple[np.ndarray, float, float]:
    """Cell weights for a sampled radius density on a 1-D support.

    Parameters
    ----------
    radii_A : np.ndarray, shape (M,)
        Sampled droplet radii [Å] of the candidate ensemble.
    support_A : np.ndarray, shape (C,)
        Cell radii [Å], strictly ascending.
    convention
        ``"nearest"`` — midpoint binning; ``"linear"`` — piecewise-linear hat
        weights. Both clamp mass outside ``[support[0], support[-1]]`` to the
        end cells.

    Returns
    -------
    (weights, clamped_below, clamped_above)
        ``weights`` sums to 1; the clamp fractions are the sample mass outside
        the support (reported, never silent).
    """
    r = np.asarray(radii_A, dtype=float)
    s = np.asarray(support_A, dtype=float)
    if s.ndim != 1 or s.size < 2 or np.any(np.diff(s) <= 0):
        raise ValueError(f"support must be ascending with >= 2 cells; got {s}")
    if r.size == 0:
        raise ValueError("empty radius sample — no density to re-weight.")

    below = float(np.mean(r < s[0]))
    above = float(np.mean(r > s[-1]))
    w = np.zeros(s.size, dtype=float)
    if convention == "nearest":
        edges = np.concatenate(([-np.inf], (s[:-1] + s[1:]) / 2.0, [np.inf]))
        idx = np.searchsorted(edges, r, side="right") - 1
        w = np.bincount(idx, minlength=s.size).astype(float)
    elif convention == "linear":
        rc = np.clip(r, s[0], s[-1])
        seg = np.clip(np.searchsorted(s, rc, side="right") - 1, 0, s.size - 2)
        t = (rc - s[seg]) / (s[seg + 1] - s[seg])
        np.add.at(w, seg, 1.0 - t)
        np.add.at(w, seg + 1, t)
    else:
        raise ValueError(
            f"unknown convention {convention!r}; expected 'nearest' or "
            f"'linear'."
        )
    return w / w.sum(), below, above


# ---------------------------------------------------------------------------
# weighted mixture of per-cell reads (pure; tested)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MixtureRead:
    """Weighted mixture of per-cell scored reads — scorer duck-type.

    Exposes exactly the surface the committed observable functions consume
    (``n_values`` / ``fraction`` for the histogram score, ``ke_by_n`` for the
    KE band ratios) so the mixture is scored by the SAME code as a single run
    (rule 1). The sufficient statistics are exact weighted combinations; no
    resampling and no replication.

    ``detected_yield`` is the weighted scored share Σ w_c (1 − trap_c) — the
    fraction of source ions the mixture actually scores; at the anchored radii
    it is well below 1 and must be reported next to every mixture row.
    """

    label: str
    trapped_frac: float
    trap_bound_frac: float
    trap_marginal_frac: float
    suppressed_frac: float
    detected_yield: float
    n_values: np.ndarray
    fraction: np.ndarray
    n_mean: float
    n1_frac: float
    _ke_n: np.ndarray
    _ke_mean_eV: np.ndarray
    _ke_count: np.ndarray

    def ke_by_n(self, *, min_count: int = 1) -> KEPerBin:
        """Weighted per-bin mean detected KE of the mixture.

        Only ``min_count=1`` (all occupied bins) is meaningful for a weighted
        mixture — a bin's "count" is a weighted quantity, so any other
        threshold would silently compare apples to oranges. Fails loudly.
        ``se_eV`` is 0 by construction (no weighted-SE convention exists; the
        band ratios never read it, and ``chi2_med`` — which would — is
        deliberately not computed for mixtures).
        """
        if min_count != 1:
            raise ValueError(
                f"MixtureRead supports min_count=1 only (weighted bin counts "
                f"have no ion-count semantics); got {min_count}."
            )
        return KEPerBin(
            n=self._ke_n.copy(),
            mean_eV=self._ke_mean_eV.copy(),
            se_eV=np.zeros_like(self._ke_mean_eV),
            count=self._ke_count.copy(),
        )


def mixture_read(
    reads: list[ConfirmationDetectionRead],
    weights: np.ndarray,
    *,
    label: str,
) -> MixtureRead:
    """Combine per-cell reads into one weighted-ensemble read.

    ``weights`` are per-cell source-ion weights (the candidate density mass on
    each cell); they are renormalized to 1. All fate fractions combine at the
    source-ion level; the scored-ensemble statistics (histogram, per-bin KE)
    combine with each cell's scored share ``1 - trapped_frac`` folded in, so
    detection-conditional observables automatically carry the detected-yield
    weighting.
    """
    if len(reads) != np.asarray(weights).size:
        raise ValueError(
            f"{len(reads)} reads vs {np.asarray(weights).size} weights."
        )
    w = np.asarray(weights, dtype=float)
    if np.any(w < 0) or w.sum() <= 0:
        raise ValueError(f"weights must be non-negative with positive sum; got {w}")
    w = w / w.sum()

    n_values = reads[0].n_values
    for rd in reads[1:]:
        if not np.array_equal(rd.n_values, n_values):
            raise ValueError("all reads must share one histogram support.")

    trap = float(np.sum(w * [rd.trapped_frac for rd in reads]))
    trap_b = float(np.sum(w * [rd.trap_bound_frac for rd in reads]))
    trap_m = float(np.sum(w * [rd.trap_marginal_frac for rd in reads]))

    # scored-ensemble weights: source weight x scored share per cell
    u = w * np.asarray([rd.num_scored / rd.num_ions for rd in reads])
    yield_frac = float(u.sum())
    if yield_frac <= 0.0:
        raise ValueError(f"mixture {label!r}: no scored mass — every "
                         f"contributing cell is fully retained.")
    u_n = u / yield_frac

    fraction = np.zeros_like(reads[0].fraction)
    for uc, rd in zip(u_n, reads):
        fraction = fraction + uc * rd.fraction
    fraction = fraction / fraction.sum()

    supp = float(np.sum(u_n * [rd.suppressed_frac for rd in reads]))
    n_mean = float(np.sum(n_values * fraction))

    # per-bin KE: weighted mean of per-cell bin means, bin weight
    # u_c * count_{c,n} / num_scored_c (exact — means are linear statistics)
    acc: dict[int, list[float]] = {}
    for uc, rd in zip(u_n, reads):
        if uc == 0.0:
            continue   # a zero-weight cell must not open empty KE bins
        kpb = rd.ke_by_n(min_count=1)
        for n_bin, mean_eV, count in zip(kpb.n, kpb.mean_eV, kpb.count):
            beta = uc * float(count) / rd.num_scored
            wsum, wmean, csum = acc.get(int(n_bin), [0.0, 0.0, 0.0])
            acc[int(n_bin)] = [wsum + beta, wmean + beta * float(mean_eV),
                               csum + float(count)]
    ns = np.asarray(sorted(acc), dtype=int)
    ke_mean = np.asarray([acc[n][1] / acc[n][0] for n in ns], dtype=float)
    ke_count = np.asarray([int(round(acc[n][2])) for n in ns], dtype=int)

    return MixtureRead(
        label=label,
        trapped_frac=trap,
        trap_bound_frac=trap_b,
        trap_marginal_frac=trap_m,
        suppressed_frac=supp,
        detected_yield=yield_frac,
        n_values=n_values.copy(),
        fraction=fraction,
        n_mean=n_mean,
        n1_frac=float(fraction[1]),
        _ke_n=ns,
        _ke_mean_eV=ke_mean,
        _ke_count=ke_count,
    )


def mixture_row(
    mix: MixtureRead | ConfirmationDetectionRead, abundance_ref, ked_ref
) -> dict[str, Any]:
    """The committed observable vector of one mixture (chi2_med absent — see
    module docstring). Also accepts a plain per-run read, so a measured
    target ensemble (the in-support oracle) prints in the same columns."""
    hist = score_histogram_vs_reference(mix, abundance_ref)
    mid = midhot_ratio(mix, ked_ref)
    deep = deep_bin_ke_ratio(mix, ked_ref)
    det_yield = (
        mix.detected_yield if isinstance(mix, MixtureRead)
        else mix.num_scored / mix.num_ions
    )
    return {
        "label": mix.label,
        "trap": mix.trapped_frac,
        "trap_bound": mix.trap_bound_frac,
        "trap_marg": mix.trap_marginal_frac,
        "det_yield": det_yield,
        "supp": mix.suppressed_frac,
        "nbar_det": mix.n_mean,
        "n1_solv": hist.sim_n1,
        "w1_solv": hist.w1_solvated,
        "midHot": mid.value,
        "midHot_bins": mid.n_bins,
        "deepKE": deep.value,
        "deepKE_bins": deep.n_bins,
    }


# ---------------------------------------------------------------------------
# densities and columns
# ---------------------------------------------------------------------------


def column_cells(column: str) -> list[tuple[str, float]]:
    """``(label, R [Å])`` of one birth-law column, ascending in R.

    R comes from the authoritative ``cfg.json`` (``single_droplet_size``
    through the shared bulk conversion), never from the run tag.
    """
    out: list[tuple[str, float]] = []
    for label in CELL_LABELS:
        if not label.endswith(column):
            continue
        cfg = json.loads(
            (RUNS_ROOT / RUN_SELECTION[label] / "cfg.json").read_text(
                encoding="utf-8"
            )
        )
        if not cfg["use_single_droplet_size"]:
            raise ValueError(f"{label}: not a fixed-size grid cell.")
        out.append(
            (label, float(droplet_radius_bulk_angstrom(cfg["single_droplet_size"])))
        )
    out.sort(key=lambda item: item[1])
    return out


def standing_radii_A() -> np.ndarray:
    """The standing ensemble's realized droplet radii [Å] — the pooled
    battery's own sampled sizes. The checkpoint stores one value per ion
    (2N; both ions of a molecule share its droplet), so every molecule
    appears twice — an exact duplication that leaves the density unchanged."""
    radii: list[np.ndarray] = []
    for name in BATTERY_MEMBERS:
        ckpt = load_neutral_checkpoint(RUNS_ROOT / name / "neutral.npz")
        radii.append(np.asarray(ckpt.droplet_radii, dtype=float))
        del ckpt
    return np.concatenate(radii)


def _minimal_detection(
    n_detected: np.ndarray,
    state_reason: np.ndarray,
    ke_eV: np.ndarray,
) -> DetectionResult:
    """Scoring-side minimal ``DetectionResult`` for a masked sub-ensemble.

    ``read_confirmation_detection`` consumes exactly three per-ion fields
    (``state_reason``, ``n_detected``, ``E_kin_detected_eV``); the rest are
    zero-filled and the event stream is empty. Never write this to disk — it
    is a scoring container, not a run artifact.
    """
    m = int(np.asarray(n_detected).size)
    zeros = np.zeros(m, dtype=float)
    return DetectionResult(
        num_molecules=max(m // 2, 1),
        t_handover_ps=30.0,
        detection_time_ps=8.53e6,
        n_detected=np.asarray(n_detected, dtype=float),
        E_int_detected_eV=zeros.copy(),
        mass_detected_kg=zeros.copy(),
        vx_detected=zeros.copy(),
        vy_detected=zeros.copy(),
        vz_detected=zeros.copy(),
        E_kin_detected_eV=np.asarray(ke_eV, dtype=float),
        E_pot_detected_eV=zeros.copy(),
        E_dissip_detected_eV=zeros.copy(),
        E_mass_transfer_detected_eV=zeros.copy(),
        state_reason=np.asarray(state_reason),
        event_offsets=np.zeros(m + 1, dtype=np.int64),
        event_time_ps=np.zeros(0, dtype=float),
        event_pre_shed_n=np.zeros(0, dtype=float),
        event_dE_int_eV=np.zeros(0, dtype=float),
        event_dE_bind_fold_eV=np.zeros(0, dtype=float),
        event_dE_mass_transfer_eV=np.zeros(0, dtype=float),
    )


def battery_in_support_target(
    r_min_A: float,
) -> tuple[ConfirmationDetectionRead, np.ndarray, float]:
    """The battery's own R >= r_min sub-ensemble: a known in-support target.

    ``droplet_radii`` is stored per ion (2N, block layout), aligned with the
    detection rows, so the mask is a direct element-wise comparison.
    Returns ``(read, masked_per_ion_radii, kept_ion_fraction)``.
    """
    n_parts: list[np.ndarray] = []
    state_parts: list[np.ndarray] = []
    ke_parts: list[np.ndarray] = []
    radii_parts: list[np.ndarray] = []
    total_ions = 0
    kept_ions = 0
    for name in BATTERY_MEMBERS:
        run_dir = RUNS_ROOT / name
        ckpt = load_neutral_checkpoint(run_dir / "neutral.npz")
        radii = np.asarray(ckpt.droplet_radii, dtype=float)
        del ckpt
        detection = load_detection_result(run_dir / "detection.npz")
        if radii.size != np.asarray(detection.n_detected).size:
            raise ValueError(
                f"{name}: {radii.size} per-ion radii vs "
                f"{np.asarray(detection.n_detected).size} detection rows — "
                f"layout assumption violated."
            )
        mask = radii >= r_min_A
        total_ions += int(mask.size)
        kept_ions += int(mask.sum())
        n_parts.append(np.asarray(detection.n_detected, dtype=float)[mask])
        state_parts.append(np.asarray(detection.state_reason)[mask])
        ke_parts.append(
            np.asarray(detection.E_kin_detected_eV, dtype=float)[mask]
        )
        radii_parts.append(radii[mask])
    target = _minimal_detection(
        np.concatenate(n_parts),
        np.concatenate(state_parts),
        np.concatenate(ke_parts),
    )
    read = read_confirmation_detection(
        target, label=f"battery R>={r_min_A:.1f}"
    )
    return read, np.concatenate(radii_parts), kept_ions / total_ions


def corrected_radii_A() -> tuple[np.ndarray, float, float]:
    """Draw the corrected ensemble's radii [Å]: `legacy` + `raw` ln-normal at
    the preset's own source conditions (G0 decision 1).

    Returns ``(radii, N_mean_drawn, N_mean_exact)`` so the nozzle-correlation
    anchor (⟨N⟩ ≈ 12794 at 40 mbar / 14 K, D0 §15) is checked, not assumed.
    """
    cfg = RunDirectory(
        RUNS_ROOT / RUN_SELECTION[CORRECTED_DENSITY_SOURCE_CELL]
    ).load_cfg()
    draw_cfg = dataclasses.replace(cfg, num_molecules=CORRECTED_SAMPLE_N)
    rng = np.random.default_rng(CORRECTED_SAMPLE_SEED)
    sizes = sample_droplet_sizes(draw_cfg, mode="raw", rng=rng)
    n_mean_exact, _ = mean_droplet_size(cfg.p_source_mbar, cfg.T_source_K)
    return (
        np.asarray(droplet_radius_bulk_angstrom(sizes), dtype=float),
        float(np.mean(sizes)),
        float(n_mean_exact),
    )


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------


def _column_reads(
    column: str, *, arm_b: bool = False
) -> tuple[list[str], np.ndarray, list[ConfirmationDetectionRead], list[str]]:
    """Load one column's cell reads (Arm A, or Arm B where marginals exist).

    Returns ``(labels, support_R, reads, arm_b_cells)`` where ``arm_b_cells``
    lists the cells whose read was actually replaced by the Arm-B injection.
    """
    labels: list[str] = []
    support: list[float] = []
    reads: list[ConfirmationDetectionRead] = []
    replaced: list[str] = []
    for label, R in column_cells(column):
        run_dir = RUNS_ROOT / RUN_SELECTION[label]
        read = load_confirmation_run(run_dir, label=label)
        if arm_b and read.trap_marginal_frac > 0.0:
            dossier = marginal_dossier(run_dir)
            read = read_confirmation_detection(
                arm_b_detection(dossier), label=f"{label}+marg"
            )
            replaced.append(label)
        labels.append(label)
        support.append(R)
        reads.append(read)
    return labels, np.asarray(support), reads, replaced


def _weight_note(labels: list[str], w: np.ndarray, below: float,
                 above: float) -> str:
    parts = ", ".join(f"{lab} {wi:.3f}" for lab, wi in zip(labels, w))
    return (f"    weights: {parts}   "
            f"(clamped below/above support: {below:.1%} / {above:.1%})")


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    print("Atlas D2b §4.3 — geometry-grid re-weighting "
          "(oracle first, then the corrected-ensemble forecast)\n")

    # -- scorer-drift oracle (identical discipline to the grid table) -------
    scorer_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(scorer_row)
    if drift:
        print("*** SCORER ORACLE DRIFT — session invalid ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer oracle OK (pooled battery row reproduced).\n")

    # -- seed-SD yardstick, measured from the five battery members ---------
    member_rows = [
        observable_row(name[-2:], RUNS_ROOT / name, abundance_ref, ked_ref,
                       with_geometry=False)
        for name in BATTERY_MEMBERS
    ]
    seed_sd = {
        col: float(np.std([r[col] for r in member_rows], ddof=1))
        for col in GATE_COLUMNS
    }
    print("seed-SD yardstick (5 x N = 1000 members): "
          + "  ".join(f"{k}={v:.4g}" for k, v in seed_sd.items()) + "\n")

    # -- densities ----------------------------------------------------------
    std_radii = standing_radii_A()
    corr_radii, n_mean_drawn, n_mean_exact = corrected_radii_A()
    for name, r_arr in (("standing (battery realized)", std_radii),
                        ("corrected (legacy+raw draw)", corr_radii)):
        q = np.percentile(r_arr, [5, 50, 95])
        print(f"{name}: M={r_arr.size}, R mean {r_arr.mean():.1f} Å, "
              f"q05/med/q95 {q[0]:.1f}/{q[1]:.1f}/{q[2]:.1f} Å")
    print(f"corrected ⟨N⟩: drawn {n_mean_drawn:.0f} vs nozzle correlation "
          f"{n_mean_exact:.0f} (D0 §15 anchor ≈ 12794)\n")

    # -- re-weighting oracle: standing density on the L3 column -------------
    print("=== Re-weighting oracle: standing density x L3 column vs the "
          "recorded pooled row ===")
    l3_labels, l3_support, l3_reads, _ = _column_reads("l3")
    print(f"L3 support: "
          + ", ".join(f"{l} R={r:.1f}" for l, r in zip(l3_labels, l3_support)))

    # the recorded row is listed first, and format_table takes its columns
    # from the first row — so it must carry the full mixture column set
    # (placeholders where the committed record has no value).
    recorded_row: dict[str, Any] = {
        "label": "pooled(recorded)",
        "trap": ORACLE_RECORDED["trap"],
        "trap_bound": "-",
        "trap_marg": "-",
        "det_yield": "-",
        "supp": ORACLE_RECORDED["supp"],
        "nbar_det": ORACLE_RECORDED["nbar_det"],
        "n1_solv": ORACLE_RECORDED["n1_solv"],
        "w1_solv": ORACLE_RECORDED["w1_solv"],
        "midHot": ORACLE_RECORDED["midHot"],
        "midHot_bins": "-",
        "deepKE": ORACLE_RECORDED["deepKE"],
        "deepKE_bins": "-",
    }
    oracle_rows: list[dict[str, Any]] = [recorded_row]
    verdicts: dict[str, tuple[int, float, dict[str, float]]] = {}
    for convention in CONVENTIONS:
        w, below, above = interpolation_weights(std_radii, l3_support,
                                                convention)
        mix = mixture_read(l3_reads, w, label=f"recon std x L3 ({convention})")
        row = mixture_row(mix, abundance_ref, ked_ref)
        oracle_rows.append(row)
        print(_weight_note(l3_labels, w, below, above))

        n_pass, score = 0, 0.0
        errors: dict[str, float] = {}
        for col in GATE_COLUMNS:
            err = float(row[col]) - float(ORACLE_RECORDED[col])
            tol = max(GATE_SD_MULT * seed_sd[col],
                      GATE_REL_TOL * abs(ORACLE_RECORDED[col]))
            errors[col] = err
            score += abs(err) / tol
            if abs(err) <= tol:
                n_pass += 1
        verdicts[convention] = (n_pass, score, errors)
    print()
    print(format_table(oracle_rows))
    print()

    for convention in CONVENTIONS:
        n_pass, score, errors = verdicts[convention]
        err_txt = "  ".join(
            f"{col} {err:+.3f} ({abs(err) / seed_sd[col]:.1f} SD)"
            for col, err in errors.items()
        )
        print(f"{convention}: {n_pass}/{len(GATE_COLUMNS)} columns pass the "
              f"pre-registered gate (Σ|err|/tol {score:.2f})\n    {err_txt}")

    adopted = max(
        CONVENTIONS,
        key=lambda c: (verdicts[c][0], -verdicts[c][1]),
    )
    adopted_pass = verdicts[adopted][0]
    admissible = adopted_pass >= GATE_MIN_PASS
    print(f"\nadopted convention: {adopted!r} "
          f"({adopted_pass}/{len(GATE_COLUMNS)} pass) -> "
          f"{'ADMISSIBLE' if admissible else 'INADMISSIBLE'} "
          f"(gate: >= {GATE_MIN_PASS}/7; per-column tol = "
          f"max({GATE_SD_MULT} x seed-SD, {GATE_REL_TOL:.0%} of recorded))")
    if not admissible:
        print("Per plan §4.3: if a known ensemble cannot be reconstructed, no "
              "reconstruction of an unknown one is admissible. The forecast "
              "below is printed as INDICATIVE ONLY; the designed remedy is a "
              "below-support confirmation cell (the §4.3 '<= 2 confirmations' "
              "budget), not a re-tune.")
    print()

    # -- in-support oracle (POST-HOC — see module docstring) ----------------
    print("=== In-support oracle (post-hoc): battery R >= support-min "
          "sub-ensemble ===")
    r_min = float(l3_support[0])
    target_read, masked_radii, kept_frac = battery_in_support_target(r_min)
    target_row = mixture_row(target_read, abundance_ref, ked_ref)
    print(f"target: pooled battery restricted to R >= {r_min:.1f} Å "
          f"({kept_frac:.1%} of ions kept — the density mass the support "
          f"actually covers).")
    in_rows: list[dict[str, Any]] = [target_row]
    in_verdicts: dict[str, tuple[int, dict[str, float]]] = {}
    for convention in CONVENTIONS:
        w, below, above = interpolation_weights(masked_radii, l3_support,
                                                convention)
        mix = mixture_read(
            l3_reads, w, label=f"recon in-support ({convention})"
        )
        row = mixture_row(mix, abundance_ref, ked_ref)
        in_rows.append(row)
        print(_weight_note(l3_labels, w, below, above))
        n_pass = 0
        errors: dict[str, float] = {}
        for col in GATE_COLUMNS:
            err = float(row[col]) - float(target_row[col])
            tol = max(GATE_SD_MULT * seed_sd[col],
                      GATE_REL_TOL * abs(float(target_row[col])))
            errors[col] = err
            if abs(err) <= tol:
                n_pass += 1
        in_verdicts[convention] = (n_pass, errors)
    print()
    print(format_table(in_rows))
    print()
    for convention in CONVENTIONS:
        n_pass, errors = in_verdicts[convention]
        err_txt = "  ".join(
            f"{col} {err:+.3f} ({abs(err) / seed_sd[col]:.1f} SD)"
            for col, err in errors.items()
        )
        print(f"{convention}: {n_pass}/{len(GATE_COLUMNS)} in-support columns "
              f"pass\n    {err_txt}")
    print(
        "\nReading: the pre-registered verdict above is on the FULL standing "
        "density and stands. This post-hoc check separates method error from "
        "support error: a pass here says column re-weighting works where the "
        "support covers the density; the full-density failure is then the "
        "53.8 % below-support clamp, not the method. The target and the "
        "grid cells are statistically independent runs, so these errors "
        "include both interpolation bias and one N = 500 seed's scatter.\n"
    )

    # -- the forecast: 2 x 2 decomposition + Arm B ranges --------------------
    print("=== Forecast (adopted convention; chi2_med intentionally absent) "
          "===")
    forecast_rows: list[dict[str, Any]] = [recorded_row]
    for density_name, radii in (("std", std_radii), ("corr", corr_radii)):
        for column in ("l3", "l2"):
            labels, support, reads_a, _ = _column_reads(column)
            w, below, above = interpolation_weights(radii, support, adopted)
            name = f"{density_name} x {column.upper()}"
            mix_a = mixture_read(reads_a, w, label=f"{name} (Arm A)")
            forecast_rows.append(mixture_row(mix_a, abundance_ref, ked_ref))
            print(_weight_note(labels, w, below, above).replace(
                "    weights", f"    {name} weights"))

            has_marginal = any(
                rd.trap_marginal_frac > 0.0 and wi > 1e-6
                for rd, wi in zip(reads_a, w)
            )
            if has_marginal:
                labels_b, support_b, reads_b, replaced = _column_reads(
                    column, arm_b=True
                )
                mix_b = mixture_read(reads_b, w, label=f"{name} (Arm B)")
                forecast_rows.append(mixture_row(mix_b, abundance_ref, ked_ref))
                print(f"    Arm B replaces cells: {replaced}")
    print()
    print(format_table(forecast_rows))
    print(
        "\nReading notes:\n"
        "  - 'corr x L2' is the corrected-geometry forecast (the G2 input); "
        "'corr x L3' and 'std x L2' decompose size vs birth-law "
        "contributions; 'std x L3' is the oracle reconstruction.\n"
        "  - Arm A excludes the droplet_retained_marginal class; Arm B "
        "injects it at exact conservative asymptotic KE. The A/B pair is the "
        "§3.5b NOT-TIGHT bracket propagated to the ensemble level — quote "
        "deepKE and nbar_det as the A–B range, never a single number.\n"
        "  - The marginal fraction is conditional on a cubic law extrapolated "
        "~3x beyond its 9/18 Å calibration band (§3.5b item 6).\n"
        "  - det_yield is the scored share of source ions: the corrected "
        "ensemble's detected subset is small-R/shallow-birth biased, and "
        "every detection-conditional column must be read together with it.\n"
        "  - Interpolation error: the adopted convention's per-observable "
        "oracle errors above travel with every forecast row (coarse 4-point "
        "R support; stated, not implied)."
    )


if __name__ == "__main__":
    main()
