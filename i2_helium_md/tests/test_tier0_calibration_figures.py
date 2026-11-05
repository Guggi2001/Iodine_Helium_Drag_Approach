"""The five Tier-0 calibration figures: run table (no propagation) and thesis layout."""

from __future__ import annotations

import struct

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from scripts.post_processing.thesis_style import TEXTWIDTH_IN  # noqa: E402
from scripts.post_processing.tier0_calibration_figures import (  # noqa: E402
    DPI,
    OUT_DIR,
    RUNS,
    build_cfg,
)
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    CAL_ALONE_FIGSIZE_IN,
    CAL_FIGSIZE_IN,
    CAL_Y_MARGINS,
    CAL_ZERO_PAD,
    _finish_calibration_axis,
)
from scripts.tier0_common import run_dir_name  # noqa: E402


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
    names = [run_dir_name(r.case, "LQ", 50, r.run_tag) for r in RUNS]
    assert len(set(names)) == len(names)


@pytest.mark.parametrize("run", RUNS, ids=lambda r: f"{r.case}_{r.run_tag}")
def test_run_configs(run):
    cfg = build_cfg(run)
    assert cfg.drag_coefficients.coefficients == {"a": 18.007189923460135, "c": 1.6}
    assert cfg.m_eff_amu == 126.90 + 21 * 4.0026
    assert (cfg.num_molecules, cfg.seed) == (50, 20260604)
    assert cfg.ion_simulation_time == run.ion_time_ps
    assert cfg.binding_energy_I_ion_eV == (run.e_bind_eV or 0.11727328616561973)
