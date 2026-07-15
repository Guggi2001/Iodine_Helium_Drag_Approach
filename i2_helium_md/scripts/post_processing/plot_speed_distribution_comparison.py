"""Interactive comparison: I+He_n per-fragment speed distributions (n = 0-3).

Shows the 2-D and 3-D per-fragment overlay figures on BOTH axis
conventions -- v-axis (the interactive twin of ``plot_run_summary.py``'s
``_section_ihe_ked_curves`` section, called once per representation,
imported from there, not reimplemented) and, additionally, E-axis
(built locally in this script -- ``run_summary`` does not gain these).
Intended for interactive IDE use (e.g. PyCharm): edit the
``USER SETTINGS`` block below, then either run this file directly::

    python scripts/post_processing/plot_speed_distribution_comparison.py

or run it from the IDE. All four figures are shown together with a
single ``plt.show()`` call.

Energy-axis measure change (Task 10)
-------------------------------------
The reference ``signal_{2d,3d}_PE`` columns are per-unit-energy, peak-
normalized densities on the ``E_eV`` axis. The simulation histogram
(``compute_final_velocity_histogram``) is a per-unit-speed density, so
overlaying it on the energy axis requires the measure change
``P(E) dE = P(v) dv`` -> ``P(E) = P(v) * |dv/dE|``. With
``E = 1/2 m v^2``, ``dE/dv = m v``, so ``P(E) ~ P(v) / v`` (the constant
``m`` drops out under peak normalization, so it is never applied).
Smoothing happens FIRST in the uniform-v domain (a boxcar on the
nonuniform E grid produced by ``v -> E`` would not be a proper moving
average), then the 1/v Jacobian, then peak normalization
(``normalise_trace``).
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    compute_final_velocity_histogram,
    energy_eV_of_speed_mps,
    fragment_mean_kinetic_energy,
    load_ihe_ked_curve,
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess._smoothing import (  # noqa: E402
    moving_mean,
    normalise_trace,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

import plot_run_summary as run_summary  # noqa: E402


# Panel count for the 2x2 per-fragment grids (n = 0..N_PANELS-1), shared
# with run_summary._section_ihe_ked_curves's panel convention.
N_PANELS = 4

# Signal-region floor for the simulation energy curves, taken from the
# reference's own convention (provenance JSON / README: curves are excluded
# below E = 0.01 eV and the envelope normalization is defined over the
# signal region). Without it the 1/v Jacobian turns any histogram weight in
# the lowest velocity bins (e.g. gate atoms moving nearly along z, whose
# in-plane speed is ~0) into an unphysical spike at E ~ 0 that would
# dominate the peak normalization and squash the physical peak.
ENERGY_SIGNAL_CUT_EV = 0.01

_ENERGY_FIGURE_TITLES = {
    "3d": ("3-D kinetic-energy distributions P(E) vs I$^+$He$_n$ "
           "reference (peak-normalized, shapes only)"),
    "2d": ("2-D detector projections P(E) vs I$^+$He$_n$ reference "
           "(peak-normalized, shapes only)"),
}


def _build_energy_figure(ion, ked_dir, ked_ref, representation) -> plt.Figure:
    """Per-fragment kinetic-energy-distribution overlays for n = 0..3.

    Energy-axis counterpart of ``run_summary._section_ihe_ked_curves``:
    same 2x2 layout, per-fragment reference raw + movmean-smoothed
    curves, mean-energy markers, and empty-gate annotation, but plotted
    against ``E_eV`` instead of ``v_mps``. The reference columns
    (``signal_{2d,3d}_PE``) are already per-unit-energy; the simulation
    histogram is converted from its native per-unit-speed density via
    the Jacobian ``P(E) ~ P(v) / v`` (module docstring has the full
    derivation).

    Parameters
    ----------
    ion : IonCheckpoint
        Ion-stage checkpoint providing the final-velocity ensemble.
    ked_dir : Path
        Directory holding ``IHe_KED_curves_n{0..4}.csv``.
    ked_ref : IHeKedReference
        Per-fragment <E> reference table (for the exp <E> markers).
    representation : {"2d", "3d"}
        Detector projection vs 3-D reconstruction.

    Returns
    -------
    plt.Figure
        2x2 grid, n = 0..3, constrained layout.

    Raises
    ------
    ValueError
        If ``representation`` is not "2d" or "3d".
    """
    if representation not in ("2d", "3d"):
        raise ValueError(f"representation must be '2d' or '3d', "
                         f"got {representation!r}")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0),
                             constrained_layout=True)
    axes = axes.ravel()
    fig.suptitle(_ENERGY_FIGURE_TITLES[representation])

    for n in range(N_PANELS):
        ax = axes[n]
        curve = load_ihe_ked_curve(ked_dir, n)
        ref_signal = getattr(curve, f"signal_{representation}_PE")
        # Raw export as faint points joined by a fine line; matplotlib
        # breaks the line at the NaN cuts, so gaps stay gaps.
        ax.plot(curve.E_eV, ref_signal, ".-", color="tab:blue",
                markersize=2.5, linewidth=0.5, alpha=0.30,
                label="experiment (raw)")
        ax.plot(curve.E_eV,
                run_summary._nan_aware_moving_mean(
                    ref_signal, run_summary.IHE_KED_EXP_SMOOTHING_WINDOW,
                ),
                color="tab:blue", linewidth=1.4,
                label="experiment (movmean 10)")

        mass_amu = float(complex_mass_amu(n))
        sim_note = None
        try:
            hist = compute_final_velocity_histogram(
                ion, mass_amu=mass_amu,
                num_bins=run_summary.IHE_KED_HIST_NUM_BINS,
                v_max_Aps=run_summary.IHE_KED_HIST_V_MAX_APS,
                projected=(representation == "2d"),
            )
        except ValueError:
            sim_note = "sim: no atoms in gate"
        else:
            # Measure change v -> E: smooth FIRST in the uniform-v domain
            # (a boxcar on the nonuniform E grid would not be a proper
            # moving average), then apply the 1/v Jacobian (P(E) ~ P(v)/v;
            # the constant m drops out under peak normalization), then
            # peak-normalize.
            v_mps = hist.bin_centers_mps
            smoothed_v = moving_mean(
                hist.density, run_summary.HIST_SMOOTHING_WINDOW,
            )
            density_E = smoothed_v / v_mps
            x = energy_eV_of_speed_mps(v_mps, mass_amu)
            # Reference signal-region cut: drop the sub-cut bins BEFORE
            # normalizing so the Jacobian's low-E artifact cannot set the
            # peak (see ENERGY_SIGNAL_CUT_EV).
            keep = x >= ENERGY_SIGNAL_CUT_EV
            x = x[keep]
            y = normalise_trace(density_E[keep])
            ax.plot(x, y, "--", color="tab:red", linewidth=1.4,
                    label=f"simulation (N={hist.num_atoms_used})")

        # <E> markers: reference values are already energies, no
        # conversion needed (unlike the v-axis v(<E>) markers).
        ax.axvline(float(ked_ref.mean_KE_eV[n]), color="tab:blue",
                   linestyle=":", linewidth=1.0,
                   label=r"exp $\langle E\rangle$")
        try:
            sim_mean = fragment_mean_kinetic_energy(ion, n)
        except ValueError:
            pass
        else:
            ax.axvline(sim_mean.mean_KE_eV, color="tab:red",
                       linestyle=":", linewidth=1.0,
                       label=r"sim $\langle E\rangle$")

        if sim_note:
            ax.annotate(sim_note, (0.97, 0.9), xycoords="axes fraction",
                        ha="right", fontsize=8, color="tab:red")
        ax.set(title=f"n = {n}",
               xlim=(0.0, float(np.nanmax(curve.E_eV))),
               ylim=(0.0, 1.25))
        ax.set_xlabel("E / eV")
        ax.set_ylabel("signal / arb. units")
        if n == 0:
            ax.legend(frameon=False, fontsize=7)
    return fig


# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Path to the run directory holding cfg.json + neutral.npz + ion.npz.
RUN_DIR: Path = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"

# Directory holding the frozen I+He_n kinetic-energy reference
# (IHe_KED_reference.csv + IHe_KED_curves_n{0..4}.csv). Required -- this
# script exists specifically to compare against the reference curves, so
# unlike plot_run_summary.py's optional section gating, None is refused
# with a clear error rather than silently skipping the figures.
IHE_KED_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked"
)


def main() -> int:
    """Load the configured run's ion checkpoint and the ihe_ked reference,
    build the 3-D and 2-D per-fragment overlay figures on both the v-axis
    (shared section builder) and the E-axis (local builder), and show all
    four together.

    Returns
    -------
    int
        0 on success (mirrors ``plot_run_summary.main``'s convention).

    Raises
    ------
    ValueError
        If ``IHE_KED_REFERENCE_DIR`` is ``None`` (this script exists to
        compare against the reference) or if the run directory has no ion
        checkpoint.
    """
    if IHE_KED_REFERENCE_DIR is None:
        raise ValueError(
            "IHE_KED_REFERENCE_DIR is None -- this script exists to "
            "compare against the I+He_n kinetic-energy reference curves; "
            "set it to a valid directory holding IHe_KED_reference.csv "
            "and IHe_KED_curves_n{0..4}.csv."
        )
    ked_dir = Path(IHE_KED_REFERENCE_DIR)

    run = RunDirectory(RUN_DIR)
    if not run.has_ion():
        raise ValueError(
            f"{RUN_DIR} has no ion checkpoint -- this is an ion-stage "
            "comparison (I+He_n speed distributions)."
        )
    ion = run.load_ion()

    ked_ref = load_ihe_ked_reference(ked_dir / "IHe_KED_reference.csv")

    run_summary._section_ihe_ked_curves(ion, ked_dir, ked_ref, "3d")
    run_summary._section_ihe_ked_curves(ion, ked_dir, ked_ref, "2d")
    _build_energy_figure(ion, ked_dir, ked_ref, "3d")
    _build_energy_figure(ion, ked_dir, ked_ref, "2d")

    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
