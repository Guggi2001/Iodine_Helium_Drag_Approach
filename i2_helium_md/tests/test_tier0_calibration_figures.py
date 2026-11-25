"""The five Tier-0 calibration figures: run table (no propagation) and thesis layout."""

from __future__ import annotations

import struct
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from scripts.post_processing.thesis_style import TEXTWIDTH_IN  # noqa: E402
from scripts.post_processing.tier0_calibration_figures import (  # noqa: E402
    DPI,
    OUT_DIR,
    RUNS,
)
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    CAL_ALONE_FIGSIZE_IN,
    CAL_FIGSIZE_IN,
    CAL_Y_MARGINS,
    CAL_ZERO_PAD,
    _finish_calibration_axis,
    surface_crossing_time_ps,
)
import numpy as np  # noqa: E402


def _png_size(path):
    with open(path, "rb") as fh:
        return struct.unpack(">II", fh.read(24)[16:24])


def test_the_five_figures():
    names = {r.velocity_png for r in RUNS} | {r.force_png for r in RUNS if r.force_png}
    assert names == {"18A_Stuck_Velocity_R.png", "18A_Stuck_Force_Energy.png",
                     "18A_B2.png", "18A_B2_Force_Energy.png", "9A_B2.png"}


def test_committed_pngs_keep_the_paired_size():
    """Pairs sit side by side in the thesis: half the text width, same size."""
    assert CAL_FIGSIZE_IN[0] == TEXTWIDTH_IN / 2
    paired = tuple(int(side * DPI) for side in CAL_FIGSIZE_IN)
    alone = tuple(int(side * DPI) for side in CAL_ALONE_FIGSIZE_IN)
    for run in RUNS:
        assert _png_size(OUT_DIR / run.velocity_png) == paired
        if run.force_png:
            assert _png_size(OUT_DIR / run.force_png) == (alone if run.force_alone else paired)


def test_only_the_18A_B2_force_figure_stands_alone():
    assert [r.force_png for r in RUNS if r.force_alone] == ["18A_B2_Force_Energy.png"]


def _panel_ylim(**kwargs):
    """y limits of a panel with a line reaching 100 only outside its x limits."""
    fig, ax = plt.subplots()
    ax.plot([0.0, 1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 100.0])
    ax.set_xlim(0.0, 2.0)
    _finish_calibration_axis(ax, "title", "label", "a", **kwargs)
    ylim = ax.get_ylim()
    plt.close(fig)
    return ylim


def test_panel_y_limits_fit_only_the_visible_lines():
    low, high = CAL_Y_MARGINS
    assert _panel_ylim() == pytest.approx((1.0 - 2.0 * low, 3.0 + 2.0 * high))


def test_zero_based_panel_keeps_a_line_at_zero_whole():
    bottom, top = _panel_ylim(from_zero=True)
    assert bottom == pytest.approx(-CAL_ZERO_PAD * top)


def test_stuck_and_b2_runs_get_distinct_directories():
    names = [r.sim.name for r in RUNS]
    assert len(set(names)) == len(names)


def _radial_ion(r_frag1, r_frag2, radius=10.0):
    """Synthetic 2N checkpoint stand-in: atoms on the +x axis at radii r [Å], t = 0, 1, ... ps."""
    r = np.vstack([r_frag1, r_frag2]).astype(float)
    return SimpleNamespace(num_molecules=len(r_frag1), time_ps=np.arange(r.shape[1], dtype=float),
                           positions_x=r, positions_y=np.zeros_like(r),
                           positions_z=np.zeros_like(r),
                           droplet_radii_angstrom=np.full(r.shape[0], radius))


def test_surface_crossing_is_the_median_first_step_at_or_past_the_surface():
    # fragment 1 crosses R_drop = 10 at t = 1, 2, 3 ps -> median 2; fragment 2 at 3 ps
    frag1 = [[5, 10, 12, 13], [5, 8, 11, 13], [5, 6, 9, 10]]
    frag2 = [[5, 6, 7, 10]] * 3
    ion = _radial_ion(frag1, frag2)
    assert surface_crossing_time_ps(ion, atom_index=0) == 2.0
    assert surface_crossing_time_ps(ion, atom_index=4) == 3.0


def test_surface_crossing_is_none_when_most_ions_stay_inside():
    frag1 = [[5, 12, 13], [5, 6, 7], [5, 6, 7]]   # one of three escapes: median = inf
    assert surface_crossing_time_ps(_radial_ion(frag1, frag1), atom_index=0) is None
    with pytest.raises(ValueError):
        surface_crossing_time_ps(_radial_ion(frag1, frag1), atom_index=6)


def test_only_the_18A_B2_force_figure_marks_the_surface_crossing():
    assert [r.force_png for r in RUNS if r.mark_surface] == ["18A_B2_Force_Energy.png"]
