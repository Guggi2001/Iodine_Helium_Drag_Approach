"""``tier2_law_comparison_figure.py`` — the B2 / B3 / h405 discussion figures.

The script re-uses the B2 script's data layer, so recipe parity is locked in
``test_tier2_b2_thesis_figure``; here: the law table, the legend label
plumbing (``sim_label``, default B2), the display offsets and that the three
builders render one series per law. Synthetic detections only (the B2 test
fixture); committed ihe_ked / abundance references; Agg backend, figures
closed, nothing written to disk.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from tests.test_tier2_b2_thesis_figure import inputs  # noqa: E402,F401  (shared fixture)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "tier2_law_comparison_figure.py"
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("tier2_law_comparison_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def laws(mod, inputs):  # noqa: F811
    size = mod.size_data(inputs["read"], inputs["abundance"])
    ke = mod.ke_data(inputs["read"], inputs["ked_ref"])
    return [mod.LawData(key=key, name=name, color=color, marker=marker, dash=dash, slot=slot,
                        size=size, ke=ke, view=inputs["view"])
            for key, (_, name, color, marker, dash, slot) in mod.LAWS.items()]


def _labels(fig) -> list[str]:
    return [t.get_text() for leg in fig.legends for t in leg.get_texts()]


class TestSettings:
    def test_law_table_is_the_three_thesis_pools(self, mod):
        assert list(mod.LAWS) == ["B2", "B3", "h405"]
        assert [v[0] for v in mod.LAWS.values()] == [
            "B2_N6000_pooled", "B3_N6000_pooled", "h405_N6000_pooled"]

    def test_colours_markers_dashes_slots_are_distinct(self, mod):
        for col in (2, 3, 4, 5):
            assert len({v[col] for v in mod.LAWS.values()}) == len(mod.LAWS)

    def test_l1_l2_are_matching_half_width_subfigures(self, mod):
        # side by side in the thesis: half the text width each, same height
        assert mod.FIG_L_1["figsize"] == mod.FIG_L_2["figsize"]
        assert mod.FIG_L_1["figsize"][0] == pytest.approx(mod.TEXTWIDTH_IN / 2)

    @pytest.mark.parametrize("label,short", [
        (r"Simulation (CC$_{5.4}$) $\bar{E}$", "CC$_{5.4}$"),
        ("Simulation, < 2 ions", "Open: < 2 ions"),
        (r"Experiment $\bar{E}$", "Experiment"),
        ("unmapped", "unmapped"),
    ])
    def test_short_legend_labels(self, mod, label, short):
        assert mod._short_label(label, mod.FIG_L_2) == short

    def test_variants_lead_with_b2_at_slot_zero(self, mod):
        assert mod.VARIANTS == {"h405": ("B2", "h405"), "h405_b3": ("B2", "B3", "h405")}
        for keys in mod.VARIANTS.values():
            assert keys[0] == "B2" and mod.LAWS["B2"][5] == 0
            assert set(keys) <= set(mod.LAWS)


class TestBuilders:
    def test_printed_names_come_from_thesis_style(self, mod):
        # B2 / h405 are repo labels; the figures print LQ / CC5.4 (thesis_style.THESIS_NAMES)
        b2 = sys.modules["scripts.post_processing.tier2_b2_thesis_figure"]
        assert b2.COMMON_STYLE["sim_label"] == mod.THESIS_NAMES["B2"] == "LQ"
        assert mod.THESIS_NAMES["B3"] == r"LQ$_\mathrm{eff}$"
        assert {k: v[1] for k, v in mod.LAWS.items()} == mod.THESIS_NAMES

    @pytest.mark.parametrize("keys", [("B2", "h405"), ("B2", "B3", "h405")])
    def test_l_1_draws_one_line_per_law(self, mod, laws, keys):
        sub = [law for law in laws if law.key in keys]
        fig = mod.build_l_1(sub)
        try:
            lin, log = fig.axes
            assert len(lin.lines) == len(sub) and log.get_yscale() == "log"
            labels = _labels(fig)
            assert all(law.name in labels for law in sub) and labels[0] == "Experiment"
            assert (mod.THESIS_NAMES["B3"] in labels) == ("B3" in keys)
        finally:
            plt.close(fig)

    def test_l_2_rejects_an_offset_first_law(self, mod, laws):
        with pytest.raises(ValueError, match="slot 0"):
            mod.build_l_2([laws[2], laws[0]])

    def test_l_2_legend_names_every_law(self, mod, laws):
        fig = mod.build_l_2(laws)
        try:
            top, bottom = fig.axes
            assert top.get_yscale() == "log"
            # panel letters continue L-1's (a) / (b) in the side-by-side figure
            assert [t.get_text() for ax in (top, bottom) for t in ax.texts
                    if t.get_text() in ("(a)", "(b)", "(c)", "(d)")] == ["(c)", "(d)"]
            assert top.get_shared_x_axes().joined(top, bottom)
            labels = _labels(fig)
            for law in laws:
                assert law.name in labels          # short half-width label
        finally:
            plt.close(fig)

    def test_l_3_overlays_every_law(self, mod, laws, inputs):  # noqa: F811
        panels = {law.key: [mod.pv_data(law.view, IHE_KED_DIR, inputs["ked_ref"], n, "3d")
                            for n in mod.FIG_L_3["n_values"]] for law in laws}
        fig = mod.build_l_3(laws, panels)
        try:
            labels = _labels(fig)
            for law in laws:
                assert f"Simulation ({law.name})" in labels
        finally:
            plt.close(fig)
