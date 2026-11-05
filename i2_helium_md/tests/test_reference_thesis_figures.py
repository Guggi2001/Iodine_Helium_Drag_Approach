"""The ``data/reference/scripts`` exporters of ``data/figures/thesis_rest/``.

``visualize_trajectory_velocity.py`` (TDDFT reference), ``plotting_histogram.py``
(I+He_n abundance) and ``experimental_reference_fragments.py``: the derived
numbers and that the builders draw every series. Agg backend, figures closed, nothing written to disk.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    load_he_abundance_reference,
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from scripts.post_processing.thesis_style import EXP_COLOR, TEXTWIDTH_IN  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_ROOT / "data" / "reference" / "scripts"
N_MAX = 20    # largest n of integrated_i_he_abundance.csv
TIER0_SCORES_CSV = (
    PROJECT_ROOT / "data" / "figures" / "thesis_microscopic_calibration"
    / "tier0_B2_B3_h405_scores.csv"
)


def _load_script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"{name}_under_test", SCRIPTS_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tddft() -> ModuleType:
    return _load_script("visualize_trajectory_velocity")


@pytest.fixture(scope="module")
def histogram() -> ModuleType:
    return _load_script("plotting_histogram")


class TestTddftReference:
    def test_exit_bound_is_the_committed_tier0_value(self, tddft):
        with open(TIER0_SCORES_CSV, newline="", encoding="utf-8") as fh:
            committed = float(next(csv.DictReader(fh))["v_exit_atom1_Aps"])
        bound = tddft.exit_speed_bound_Aps(load_hedft_trajectory(tddft.HEDFT_9A_PATH))
        # same data and window as the Tier-0 scoring; only CSV rounding may differ
        assert bound == pytest.approx(committed, rel=1e-9)

    def test_exit_bound_rejects_an_empty_window(self, tddft, monkeypatch):
        monkeypatch.setattr(tddft, "EXIT_WINDOW_PS", (1.0e6, 2.0e6))
        with pytest.raises(ValueError, match="exit window"):
            tddft.exit_speed_bound_Aps(load_hedft_trajectory(tddft.HEDFT_9A_PATH))

    def test_builder_draws_every_series(self, tddft):
        fig = tddft.build_figure(
            load_hedft_trajectory(tddft.HEDFT_9A_PATH),
            load_hedft_trajectory(tddft.HEDFT_18A_PATH),
            load_smoothed_speed_reference(tddft.SMOOTHED_18A_PATH),
        )
        ax_9A, ax_18A = fig.axes
        assert len(ax_9A.get_lines()) == 3    # two atoms, exit bound
        assert len(ax_18A.get_lines()) == 5   # two atoms, cleaned, two thresholds
        plt.close(fig)


class TestAbundanceHistogram:
    def test_reweighted_bins_exclude_bare_iodine_and_sum_to_100(self, histogram):
        ref = load_he_abundance_reference(histogram.ABUNDANCE_CSV)
        n, percent = histogram.reweighted_percent(ref)
        np.testing.assert_array_equal(n, np.arange(1, N_MAX + 1))
        assert percent.sum() == pytest.approx(100.0)

    def test_builder_draws_one_bar_per_bin(self, histogram):
        ref = load_he_abundance_reference(histogram.ABUNDANCE_CSV)
        fig = histogram.build_figure(ref)
        assert len(fig.axes[0].patches) == N_MAX
        plt.close(fig)

    def test_with_bare_figure_shows_all_bins_on_a_log_axis(self, histogram):
        ref = load_he_abundance_reference(histogram.ABUNDANCE_CSV)
        fig = histogram.build_figure_with_bare(ref)
        ax = fig.axes[0]
        heights = [bar.get_height() for bar in ax.patches]
        assert ax.get_yscale() == "log"
        assert len(heights) == N_MAX + 1
        assert sum(heights) == pytest.approx(100.0)
        assert heights[0] == pytest.approx(ref.ion_fraction[0] * 100.0)
        assert ax.get_ylim()[0] < min(heights)    # every bar rises above the axis floor
        plt.close(fig)

    @pytest.mark.parametrize("builder", ["build_figure", "build_figure_with_bare"])
    def test_figures_are_half_text_width_in_the_experiment_colour(self, histogram, builder):
        ref = load_he_abundance_reference(histogram.ABUNDANCE_CSV)
        fig = getattr(histogram, builder)(ref)
        assert fig.get_figwidth() == pytest.approx(TEXTWIDTH_IN / 2)
        assert fig.axes[0].patches[0].get_facecolor() == pytest.approx(to_rgba(EXP_COLOR))
        plt.close(fig)

    def test_log_axis_rejects_an_empty_bin(self, histogram):
        with pytest.raises(ValueError, match="positive bins"):
            histogram._bar_figure(np.array([0, 1]), np.array([100.0, 0.0]), "y", log_y=True)
        plt.close("all")


class TestExperimentalFragments:
    @pytest.fixture(scope="class")
    def fragments(self) -> ModuleType:
        return _load_script("experimental_reference_fragments")

    @pytest.fixture(scope="class")
    def data(self, fragments) -> tuple[dict, dict]:
        curves = {n: fragments.load_curves(n) for n in fragments.N_FRAGMENTS}
        return curves, fragments.load_reference_means()

    def test_movmean_keeps_the_nan_cuts(self, fragments):
        smoothed = fragments.movmean(np.array([np.nan, 1.0, 3.0, 5.0]), 3)
        assert np.isnan(smoothed[0])
        np.testing.assert_allclose(smoothed[1:], [2.0, 3.0, 4.0])

    def test_fragment_colors_are_the_viridis_samples(self, fragments):
        expected = [plt.get_cmap("viridis")(0.85 * n / 4) for n in fragments.N_FRAGMENTS]
        np.testing.assert_allclose(fragments.FRAGMENT_COLORS, expected)

    @pytest.mark.parametrize("keys, shape, mean_lines", [
        (("2d_v", "2d_E", "3d_v", "3d_E"), (2, 2), (0, 0, 5, 5)),
        (("2d_v", "3d_v"), (1, 2), (0, 5)),
    ])
    def test_every_panel_draws_every_fragment(self, fragments, data, keys, shape, mean_lines):
        fig = fragments.build_figure(keys, shape, (TEXTWIDTH_IN, 3.0), *data)
        n_fragments = len(fragments.N_FRAGMENTS)
        assert [len(ax.get_lines()) for ax in fig.axes] == [
            2 * n_fragments + extra for extra in mean_lines   # raw + smoothed + <E> markers
        ]
        assert len(fig.legends[0].get_texts()) == n_fragments
        plt.close(fig)
