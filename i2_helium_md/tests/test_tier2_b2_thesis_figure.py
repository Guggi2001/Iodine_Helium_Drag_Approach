"""``tier2_b2_thesis_figure.py`` — the B2 thesis figures (B2-1 / B2-2 / B2-3).

Locks recipe parity with ``summary_sections`` (the script re-draws, it must
not re-derive): every plotted array of the new data layer equals the one the
existing detection-summary section draws, up to the Å/ps <-> m/s display
conversion (exact factor 100). Also locks the style merge, the anchor guard
and that the builders render. Synthetic detections only (same fixture as
``test_plot_detection_summary``); committed ihe_ked / abundance references;
Agg backend, figures closed, nothing written to disk.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from tests.test_detected_view import make_rich_detection  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    detected_ensemble_view,
    load_he_abundance_reference,
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "tier2_b2_thesis_figure.py"
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("tier2_b2_thesis_figure_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def sections(mod):
    # imported by the script, so its directory is on sys.path
    import summary_sections
    return summary_sections


@pytest.fixture(scope="module")
def inputs():
    # 6 molecules / 12 ions: solvated n = 1 / 2 / 3 (several each, so the
    # chi^2 min-count passes), one suppressed, one retained; A/ps speeds.
    n = [1.0, 1.0, 2.0, 2.0, 1.0, 3.0, 2.0, 1.0, 6.0, 9.0, 3.0, 1.0]
    reason = ["frozen"] * 8 + ["suppressed", "droplet_retained", "frozen", "frozen"]
    det = make_rich_detection(
        n, reason,
        vx=[9.0, 9.5, 8.0, 8.4, 9.2, 7.5, 8.1, 9.8, 8.0, 99.0, 7.2, 9.1],
        vy=[1.0, 0.5, 1.0, 0.2, 0.3, 0.4, 0.9, 0.1, 6.0, 99.0, 0.6, 0.8],
        vz=[0.1, 0.2, 0.3, 0.1, 0.0, 0.2, 0.1, 0.3, 0.0, 99.0, 0.1, 0.2],
    )
    read = read_confirmation_detection(det, label="synthetic")
    return dict(
        read=read,
        view=detected_ensemble_view(det),
        ked_ref=load_ihe_ked_reference(IHE_KED_DIR / "IHe_KED_reference.csv"),
        abundance=load_he_abundance_reference(ABUNDANCE),
    )


class TestStyleMerge:
    def test_override_wins_and_nested_dicts_merge(self, mod):
        common = {"sim": {"color": "b", "lw": 1.0}, "a": 1}
        st = mod.merged_style({"sim": {"color": "r"}, "a": 2, "new": 3}, common)
        assert st == {"sim": {"color": "r", "lw": 1.0}, "a": 2, "new": 3}
        assert common["sim"]["color"] == "b"   # common untouched


class TestRecipeParity:
    def test_size_matches_section_bars(self, mod, sections, inputs):
        d = mod.size_data(inputs["read"], inputs["abundance"])
        fig = sections._section_detected_size_distribution(inputs["read"], inputs["abundance"])
        try:
            scored, solv = fig.axes[0], fig.axes[1]
            np.testing.assert_allclose([p.get_height() for p in scored.patches],
                                       100.0 * d.frac_scored)
            heights = [p.get_height() for p in solv.patches]
            n_ref = d.n_ref_solv.size
            np.testing.assert_allclose(heights[:n_ref], 100.0 * d.frac_ref_solv)
            np.testing.assert_allclose(heights[n_ref:], 100.0 * d.frac_sim_solv)
        finally:
            plt.close(fig)

    def test_ke_scores_match_section(self, mod, sections, inputs):
        d = mod.ke_data(inputs["read"], inputs["ked_ref"])
        med, _ = sections._detected_ke_scores(inputs["read"], inputs["ked_ref"])
        assert d.chi2_score.chi2_profiled == med.chi2_profiled
        np.testing.assert_array_equal(d.chi2_score.z_profiled, med.z_profiled)
        # scored bins are exactly the section's min-count bins
        ke2 = inputs["read"].ke_by_n(min_count=sections.DETECTED_KE_MIN_BIN_COUNT)
        keep = ke2.n >= sections.DETECTED_KE_N_MIN
        np.testing.assert_array_equal(d.n[d.scored_mask], ke2.n[keep])
        np.testing.assert_allclose(d.mean_eV[d.scored_mask], ke2.mean_eV[keep])

    @pytest.mark.parametrize("representation", ["3d", "2d"])
    def test_pv_matches_section_lines(self, mod, sections, inputs, representation):
        fig = sections._section_ihe_ked_curves(
            inputs["view"], IHE_KED_DIR, inputs["ked_ref"], representation)
        try:
            for n in (1, 2):
                d = mod.pv_data(inputs["view"], IHE_KED_DIR, inputs["ked_ref"], n,
                                representation)
                lines = fig.axes[n].get_lines()
                # section order: exp raw, exp smoothed, sim, exp v(<E>), sim v(<E>)
                raw, smooth, sim = (ln.get_xydata() for ln in lines[:3])
                np.testing.assert_allclose(raw[:, 0], 100.0 * d.v_exp_Aps)
                np.testing.assert_allclose(smooth[:, 1], d.exp_smooth, equal_nan=True)
                np.testing.assert_allclose(sim[:, 0], 100.0 * d.v_sim_Aps)
                np.testing.assert_allclose(sim[:, 1], d.sim_trace)
                np.testing.assert_allclose(lines[3].get_xdata()[0], 100.0 * d.v_mean_E_exp_Aps)
                np.testing.assert_allclose(lines[4].get_xdata()[0], 100.0 * d.v_mean_E_sim_Aps)
        finally:
            plt.close(fig)

    def test_ratio_is_the_band_ratio_ingredient(self, mod, inputs):
        # per-bin ratio = ke_band_ratio's (sim mean / raw ref mean), so the
        # plotted points are the ones midHot / deepKE aggregate
        from i2_helium_md.postprocess.tier2_confirmation import ke_band_ratio
        d = mod.ke_data(inputs["read"], inputs["ked_ref"])
        ref = ke_band_ratio(inputs["read"], inputs["ked_ref"], n_lo=1, n_hi=int(d.n.max()))
        np.testing.assert_allclose(d.ratio[np.isin(d.n, ref.n)], ref.ratio)
        assert d.ratio_err.shape == d.ratio.shape

    def test_tail_mass_is_the_winners_table_column(self, mod, inputs):
        # mass_n15plus_solv = solvated share at n >= 15 (tier2_g5_winners_table);
        # tail_n 3 so the two synthetic n = 3 ions land in it
        d = mod.size_data(inputs["read"], inputs["abundance"], tail_n=3)
        np.testing.assert_allclose(d.sim_tail_mass, d.frac_sim_solv[d.n_sim_solv >= 3].sum())
        assert d.sim_tail_mass > 0.0
        assert 0.0 < d.ref_tail_mass < 1.0

    def test_pv_rejects_unknown_representation(self, mod, inputs):
        with pytest.raises(ValueError, match="representation"):
            mod.pv_data(inputs["view"], IHE_KED_DIR, inputs["ked_ref"], 1, "1d")


class TestAnchorGuard:
    def _csv(self, tmp_path, values, label="B2_pp_c160", container="c"):
        path = tmp_path / "winners.csv"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["label", "container", *values])
            w.writeheader()
            w.writerow({"label": label, "container": container, **values})
        return path

    def test_passes_on_4dp_match_and_raises_on_drift(self, mod, tmp_path):
        path = self._csv(tmp_path, {"nbar_det": 3.90284712, "midHot": 1.0016225})
        mod.check_anchor({"nbar_det": 3.90284, "midHot": 1.00162}, path, "B2_pp_c160", "c")
        with pytest.raises(ValueError, match="nbar_det"):
            mod.check_anchor({"nbar_det": 3.9040, "midHot": 1.00162}, path, "B2_pp_c160")

    def test_raises_on_wrong_container_or_missing_label(self, mod, tmp_path):
        path = self._csv(tmp_path, {"nbar_det": 3.9})
        with pytest.raises(ValueError, match="belongs to"):
            mod.check_anchor({"nbar_det": 3.9}, path, "B2_pp_c160", "other")
        with pytest.raises(ValueError, match="expected one row"):
            mod.check_anchor({"nbar_det": 3.9}, path, "B3")

    def test_anchor_values_cover_quoted_scalars(self, mod, inputs):
        v = mod.anchor_values(mod.size_data(inputs["read"], inputs["abundance"]),
                              mod.ke_data(inputs["read"], inputs["ked_ref"]))
        assert set(v) == {"num_scored", "supp", "nbar_det", "n1_solv", "w1_solv",
                          "mass_n15plus_solv", "midHot", "midHot_bins", "deepKE",
                          "deepKE_bins", "chi2_med", "ke_npts", "KE1_mean"}


class TestBuilders:
    @pytest.mark.parametrize("layout,n_axes", [("1x2", 2), ("2x1", 2), ("solvated_only", 1)])
    def test_b2_1_layouts(self, mod, inputs, layout, n_axes):
        fig = mod.build_b2_1(mod.size_data(inputs["read"], inputs["abundance"]),
                             {**mod.FIG_B2_1, "layout": layout})
        try:
            assert len(fig.axes) == n_axes
        finally:
            plt.close(fig)

    def test_b2_2_split_matches_section(self, mod, sections):
        assert mod.FIG_B2_2["split_n"] == sections.IHE_KED_MEAN_ENERGY_SPLIT_N

    @pytest.mark.parametrize("thin", [False, True])
    @pytest.mark.parametrize("lower", ["ratio", "z"])
    def test_b2_2_three_panels_linear(self, mod, inputs, thin, lower):
        fig = mod.build_b2_2(mod.ke_data(inputs["read"], inputs["ked_ref"]),
                             {**mod.FIG_B2_2, "show_thin_bins": thin, "lower_panel": lower})
        try:
            assert len(fig.axes) == 3
            assert all(ax.get_yscale() == "linear" for ax in fig.axes[:2])
        finally:
            plt.close(fig)

    def test_b2_2_rejects_unknown_lower_panel(self, mod, inputs):
        with pytest.raises(ValueError, match="lower_panel"):
            mod.build_b2_2(mod.ke_data(inputs["read"], inputs["ked_ref"]),
                           {**mod.FIG_B2_2, "lower_panel": "chi2"})

    def test_figures_span_the_text_width(self, mod):
        for fig_style in (mod.FIG_B2_1, mod.FIG_B2_2, mod.FIG_B2_3, mod.FIG_B2_4):
            assert fig_style["figsize"][0] == pytest.approx(mod.TEXTWIDTH_CM / 2.54)

    def test_b2_3_default_is_n1_to_n6_with_figure_legend(self, mod, inputs):
        assert mod.FIG_B2_3["n_values"] == tuple(range(1, 7))
        panels = [mod.pv_data(inputs["view"], IHE_KED_DIR, inputs["ked_ref"], n, "3d")
                  for n in mod.FIG_B2_3["n_values"]]
        fig = mod.build_b2_3(panels)
        try:
            assert len(fig.axes) == 6
            assert len(fig.legends) == 1
            assert all(ax.get_legend() is None for ax in fig.axes)
            # per-n speed range; inner column / upper rows carry no axis labels
            assert fig.axes[0].get_xlim() == mod.FIG_B2_3["v_xlim_by_n"][1]
            assert fig.axes[1].get_ylabel() == "" and fig.axes[0].get_xlabel() == ""
            assert fig.axes[4].get_xlabel() and fig.axes[4].get_ylabel()
        finally:
            plt.close(fig)

    def test_b2_4_is_b2_3_on_the_2d_projection(self, mod):
        diff = {k for k in mod.FIG_B2_3 if mod.FIG_B2_3[k] != mod.FIG_B2_4[k]}
        assert diff == {"stem", "representation"}
        assert mod.FIG_B2_4["representation"] == "2d"

    def test_b2_3_grid_hides_spare_axes_and_figure_legend(self, mod, inputs):
        panels = [mod.pv_data(inputs["view"], IHE_KED_DIR, inputs["ked_ref"], n, "3d")
                  for n in (1, 2, 3)]
        fig = mod.build_b2_3(panels, {**mod.FIG_B2_3, "legend_panel": "figure"})
        try:
            assert len(fig.axes) == 4
            assert not fig.axes[3].axison
            assert len(fig.legends) == 1
        finally:
            plt.close(fig)
