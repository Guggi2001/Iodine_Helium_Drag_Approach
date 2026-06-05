"""Tier-0 TDDFT drag comparison for one finished drag run.

Implements the §6.4 Tier-0 validation step specified in
``TIER0_COMPARISON_spec.md``: score a deterministic, fixed-``m_eff``,
``linear_cubic`` drag run against the TDDFT distance + velocity-magnitude
reference traces **inside the extraction window only**, and report the numbers
that set (or check) the Tier-0 acceptance thresholds.

What this script does
---------------------
1. Loads a finished drag run (``RunDirectory``) and the matching HeDFT
   reference CSV (``9A_All_Data.csv`` / ``18A_All_Data.csv``).
2. Reads the scored window ``[t_start, t_end]`` and ``m_eff`` from the drag
   ``fit_parameters.json`` -- the *same provenance source* as the coefficients
   (spec §3.1), so the scored region can never drift from the extraction's
   actual ``[t*, t_end]``. No hard-coded 2.67 / 8.5.
3. Asserts MD<->HeDFT time-origin alignment (spec §3.2) as a hard safeguard.
4. Scores, **in-window only** (spec §4):
     * distance RMSE (angstrom),
     * mean(I1, I2) velocity-magnitude RMSE (angstrom/ps).
   and reports, as diagnostics (not gates): ``mean_ratio`` and the I1-I2
   velocity split (symmetry sanity, spec §4).
5. Plots the **full** trajectory (plotted-but-not-scored, spec §3 D2) with the
   scored window shaded.
6. Optionally exports the ensemble-mean trajectories
   ``{time_ps, mean_distance_A, mean_speed_I1_Aps, mean_speed_I2_Aps}`` to a
   small CSV -- the committed regression reference for
   ``tests/test_tier0_drag_comparison.py`` (spec §7). These four series fully
   reproduce every scored + diagnostic number, since the comparisons reduce to
   means over molecules.

Interpretation guard (spec §2): a mid-window-good / ends-poorer residual is the
*expected signature of the constant-m_eff approximation*, NOT a ``linear_cubic``
failure. Do not promote ``linear_quadratic`` on end-of-window divergence alone;
a real form problem shows as mid-window shape mismatch.

How to use
----------
Edit USER SETTINGS, then::

    python scripts/post_processing/tier0_drag_comparison.py
"""

from __future__ import annotations

import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
# The finished drag run to score (produced by single_pulse_N2000_drag).
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_drag_tier0_N50"

# HeDFT reference trace for the same case (9A or 18A).
HEDFT_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"

# Drag coefficient directory -- the window [t_start, t_end] and m_eff are read
# from its fit_parameters.json (same provenance as the coefficients).
DRAG_COEFF_DIR = (
    PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "linear_and_cubic"
)

# Optional: write the ensemble-mean series here (the regression reference).
# Set to None to skip the export and only print + plot.
EXPORT_MEAN_SERIES_PATH = (
    PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier0"
    / "md_mean_trajectory.csv"
)

# Show the figure window. Off for headless / batch use.
SHOW_FIGURE = True


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    load_hedft_trajectory,
)
from i2_helium_md.postprocess.compare_trajectories import (  # noqa: E402
    compare_distance,
    compare_velocity_magnitude,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# Tolerance for the MD<->HeDFT time-origin alignment safeguard (spec §3.2).
_ORIGIN_TOL_PS = 1e-6


def read_drag_window(coeff_dir: Path) -> tuple[float, float, float]:
    """Read ``(t_start, t_end, meff_amu)`` from ``coeff_dir/fit_parameters.json``.

    The window and the extraction mass share the coefficients' provenance file
    (spec §3.1), so the scored region travels with the calibration.
    """
    params_path = coeff_dir / "fit_parameters.json"
    data = json.loads(params_path.read_text(encoding="utf-8"))
    for key in ("t_start", "t_end", "meff_amu"):
        if key not in data:
            raise KeyError(
                f"{params_path} is missing required key {key!r}; "
                f"present keys: {sorted(data)}"
            )
    t_start = float(data["t_start"])
    t_end = float(data["t_end"])
    if not t_end > t_start:
        raise ValueError(
            f"degenerate extraction window from {params_path}: "
            f"t_start={t_start}, t_end={t_end}"
        )
    return t_start, t_end, float(data["meff_amu"])


def assert_time_origin_aligned(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    window: tuple[float, float],
    *,
    tol_ps: float = _ORIGIN_TOL_PS,
) -> None:
    """Hard precondition: MD ion-stage t=0 coincides with HeDFT t=0 (spec §3.2).

    A half-ps offset would silently score the wrong segment. The ion stage
    starts at charge-switch (t=0) and the HeDFT trace starts at ionization
    (t=0); both axes are in ps. Also checks the scored window lies inside both
    series so the windowed overlap is not silently clipped.
    """
    t0_md = float(ion.time_ps[0])
    t0_ref = float(hedft.time_ps[0])
    if abs(t0_md - t0_ref) > tol_ps:
        raise ValueError(
            "MD<->HeDFT time origins are not aligned: "
            f"ion.time_ps[0]={t0_md} ps, hedft.time_ps[0]={t0_ref} ps "
            f"(tol={tol_ps} ps). Resolve the offset before scoring (spec §3.2)."
        )

    t_start, t_end = window
    md_lo, md_hi = float(ion.time_ps[0]), float(ion.time_ps[-1])
    ref_lo, ref_hi = float(hedft.time_ps[0]), float(hedft.time_ps[-1])
    lo = max(md_lo, ref_lo)
    hi = min(md_hi, ref_hi)
    if t_start < lo - tol_ps or t_end > hi + tol_ps:
        raise ValueError(
            f"scored window [{t_start}, {t_end}] ps is not contained in the "
            f"MD<->HeDFT overlap [{lo}, {hi}] ps; the windowed RMSE would be "
            "silently clipped (spec §3.2)."
        )


def ensemble_mean_series(
    ion: IonCheckpoint,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(time_ps, mean_distance_A, mean_speed_I1, mean_speed_I2)``.

    These are exactly the molecule-means the comparisons reduce to, so they are
    a lossless reference for every Tier-0 scored + diagnostic quantity.
    """
    n = ion.num_molecules
    dx = ion.positions_x[:n] - ion.positions_x[n:]
    dy = ion.positions_y[:n] - ion.positions_y[n:]
    dz = ion.positions_z[:n] - ion.positions_z[n:]
    mean_distance = np.mean(np.sqrt(dx * dx + dy * dy + dz * dz), axis=0)

    speed = np.sqrt(
        ion.velocities_x**2 + ion.velocities_y**2 + ion.velocities_z**2
    )
    mean_speed_i1 = np.mean(speed[:n], axis=0)
    mean_speed_i2 = np.mean(speed[n:], axis=0)
    return (
        np.asarray(ion.time_ps, dtype=float),
        mean_distance,
        mean_speed_i1,
        mean_speed_i2,
    )


def export_mean_series(
    path: Path,
    ion: IonCheckpoint,
    *,
    provenance: str,
) -> None:
    """Write the four ensemble-mean series to a small, inspectable CSV."""
    t, dist, v1, v2 = ensemble_mean_series(ion)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        f"{provenance}\n"
        "time_ps,mean_distance_A,mean_speed_I1_Aps,mean_speed_I2_Aps"
    )
    np.savetxt(
        path,
        np.column_stack([t, dist, v1, v2]),
        delimiter=",",
        header=header,
        comments="# ",
    )
    print(f"Wrote mean-series reference ({t.size} rows) -> {path}")


def score(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    window: tuple[float, float],
) -> dict:
    """Compute the Tier-0 scored numbers + diagnostics for one run."""
    dist = compare_distance(ion, hedft, window=window)
    v1 = compare_velocity_magnitude(ion, hedft, atom="I1", window=window)
    v2 = compare_velocity_magnitude(ion, hedft, atom="I2", window=window)

    mean_velocity_rmse = 0.5 * (v2.rmse + v2.rmse)
    velocity_split = abs(v1.rmse - v2.rmse)
    return {
        "window": window,
        "n_scored": dist.num_overlap_points,
        "distance_rmse_A": dist.rmse,
        "distance_mean_ratio": dist.mean_ratio,
        "v_I1_rmse_Aps": v1.rmse,
        "v_I2_rmse_Aps": v2.rmse,
        "v_mean_rmse_Aps": mean_velocity_rmse,
        "v_split_Aps": velocity_split,
        "v_I1_mean_ratio": v1.mean_ratio,
        "v_I2_mean_ratio": v2.mean_ratio,
    }


def print_summary(metrics: dict, *, label: str) -> None:
    t_start, t_end = metrics["window"]
    print()
    print(f"===== Tier-0 comparison: {label} =====")
    print(f"  scored window      : [{t_start:.3f}, {t_end:.3f}] ps "
          f"({metrics['n_scored']} HeDFT samples)")
    print(f"  GATE distance RMSE : {metrics['distance_rmse_A']:.4f} A")
    print(f"  GATE mean |v| RMSE : {metrics['v_mean_rmse_Aps']:.4f} A/ps "
          f"(= mean of I1={metrics['v_I1_rmse_Aps']:.4f}, "
          f"I2={metrics['v_I2_rmse_Aps']:.4f})")
    print("  -- diagnostics (reported, not gated) --")
    print(f"  distance mean_ratio: {metrics['distance_mean_ratio']:.4f}")
    print(f"  I1-I2 vel split    : {metrics['v_split_Aps']:.4f} A/ps")
    print(f"  vel mean_ratio I1/I2: {metrics['v_I1_mean_ratio']:.4f} / "
          f"{metrics['v_I2_mean_ratio']:.4f}")
    print("=" * (28 + len(label)))


def build_figure(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    window: tuple[float, float],
):
    """Full-trajectory distance + velocity panels with the scored window shaded.

    The threshold reads only the windowed RMSE; the full trajectory (pre-t*
    transient + post-window divergence) is shown but not scored (spec §3 D2).
    """
    import matplotlib.pyplot as plt

    t_md, dist_md, v1_md, v2_md = ensemble_mean_series(ion)
    t_start, t_end = window

    fig, (ax_d, ax_v) = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)

    ax_d.plot(t_md, dist_md, color="tab:blue", lw=1.5, label="MD mean")
    ax_d.plot(hedft.time_ps, hedft.distance_A, color="black", lw=1.5,
              label="HeDFT / TDDFT")
    ax_d.axvspan(t_start, t_end, color="tab:green", alpha=0.12,
                 label="scored window")
    ax_d.set_ylabel(r"$R_1 - R_2$ / $\mathrm{\AA}$")
    ax_d.legend(frameon=False)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)

    ax_v.plot(t_md, v1_md, color="tab:blue", lw=1.2, label="MD mean |v| I1")
    ax_v.plot(t_md, v2_md, color="tab:cyan", lw=1.2, label="MD mean |v| I2")
    ax_v.plot(hedft.time_ps, hedft.v1_magnitude_Aps, color="black", lw=1.2,
              ls="-", label="HeDFT |v1|")
    ax_v.plot(hedft.time_ps, hedft.v2_magnitude_Aps, color="dimgray", lw=1.2,
              ls=":", label="HeDFT |v2|")
    ax_v.axvspan(t_start, t_end, color="tab:green", alpha=0.12)
    ax_v.set_ylabel(r"$|v|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.set_xlabel("t / ps")
    ax_v.legend(frameon=False, ncol=2)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)

    fig.suptitle(f"Tier-0 drag comparison  (r0={hedft.droplet_radius_A:.0f} A)")
    fig.tight_layout()
    return fig


def main() -> int:
    ion = RunDirectory(RUN_DIR).load_ion()
    hedft = load_hedft_trajectory(HEDFT_PATH)
    t_start, t_end, meff_amu = read_drag_window(DRAG_COEFF_DIR)
    window = (t_start, t_end)

    print(
        f"Loaded drag ion checkpoint: N={ion.num_molecules}, "
        f"time=[{ion.time_ps[0]:.3f}, {ion.time_ps[-1]:.3f}] ps "
        f"({ion.time_ps.size} samples)"
    )
    print(
        f"Loaded HeDFT reference: r0={hedft.droplet_radius_A:.1f} A, "
        f"time=[{hedft.time_ps[0]:.3f}, {hedft.time_ps[-1]:.3f}] ps "
        f"({hedft.time_ps.size} samples)"
    )
    print(f"Drag window from fit_parameters.json: "
          f"[{t_start:.3f}, {t_end:.3f}] ps, m_eff={meff_amu:.6f} amu")

    assert_time_origin_aligned(ion, hedft, window)

    metrics = score(ion, hedft, window)
    print_summary(metrics, label=f"r0={hedft.droplet_radius_A:.0f}A N={ion.num_molecules}")

    if EXPORT_MEAN_SERIES_PATH is not None:
        provenance = (
            f"Tier-0 MD ensemble-mean trajectory. run={RUN_DIR.name}, "
            f"N={ion.num_molecules}, ion_steps={ion.time_ps.size}, "
            f"window=[{t_start},{t_end}] ps, m_eff={meff_amu} amu."
        )
        export_mean_series(EXPORT_MEAN_SERIES_PATH, ion, provenance=provenance)

    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        build_figure(ion, hedft, window)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
