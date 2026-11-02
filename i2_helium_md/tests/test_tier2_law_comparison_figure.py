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
import numpy as np  # noqa: E402
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
                        size=size, ke=ke, view=inputs["view"], read=inputs["read"])
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
        assert mod.VARIANTS == {"lq": ("B2",), "h405": ("B2", "h405"),
                                "b3": ("B2", "B3")}
        for keys in mod.VARIANTS.values():
            assert keys[0] == "B2" and mod.LAWS["B2"][5] == 0
            assert set(keys) <= set(mod.LAWS)

    def test_figure_variants(self, mod):
        # the n = 0 figures are LQ-only; L-1 to L-3 keep the law pairs
        assert set(mod.FIGURE_VARIANTS) == set(mod.FIGURES)
        for fig_key, suffixes in mod.FIGURE_VARIANTS.items():
            assert set(suffixes) <= set(mod.VARIANTS)
            if fig_key in ("l_4", "l_5"):
                assert suffixes == ("lq",)
            else:
                assert "lq" not in suffixes


class TestBareIonInputs:
    def _cfg(self, tmp_path, scale: float) -> Path:
        path = tmp_path / "cfg.json"
        path.write_text(f'{{"R0_GS_angstrom": 2.666, "E_coulomb_scale": {scale}}}',
                        encoding="utf-8")
        return path

    def test_benchmark_matches_the_documented_value(self, mod, tmp_path):
        # pair budget e^2/(4 pi eps0 2.666 A) = 5.4015 eV, half per fragment
        v = mod.coulomb_benchmarks_eV(self._cfg(tmp_path, 1.0))
        assert v == {"symmetric": pytest.approx(5.4015 / 2, abs=1e-3)}
        mod.check_benchmarks(v)                       # 2.70 to 2 dp

    def test_benchmark_guard_raises_off_the_production_budget(self, mod, tmp_path):
        with pytest.raises(ValueError, match="benchmark mismatch"):
            mod.check_benchmarks(mod.coulomb_benchmarks_eV(self._cfg(tmp_path, 0.8)))

    @pytest.mark.parametrize("n", [0, 1, 2])
    def test_ke_moments_are_the_scored_mean_and_median(self, mod, inputs, n):  # noqa: F811
        read = inputs["read"]
        sel = read.ke_scored_eV[read.n_scored == n]
        mean, median = mod.ke_moments(read, n)
        assert mean == pytest.approx(sel.mean(), rel=1e-12)
        assert median == pytest.approx(np.median(sel), rel=1e-12)

    def test_exp_moments_come_from_the_reference(self, mod, inputs):  # noqa: F811
        # current v2 reference (IHe_KED_reference.csv), not the v1 log numbers
        assert mod.exp_ke_moments(inputs["ked_ref"], 0) == pytest.approx((3.59357, 3.63925))
        assert mod.exp_ke_moments(inputs["ked_ref"], 1) == pytest.approx((1.19751, 1.05270))

    def test_ke_moments_raise_on_an_empty_bin(self, mod, inputs):  # noqa: F811
        with pytest.raises(ValueError, match="n = 19"):
            mod.ke_moments(inputs["read"], 19)


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
            # every law on B2-3's kinetic-energy axis (P(E) via ked_traces)
            ax = fig.axes[0]
            t = mod.ked_traces(panels["h405"][0], "E", mod.FIG_B2_3["e_min_eV"])
            assert ax.get_xlim() == mod.FIG_B2_3["e_xlim_by_n"][1]
            np.testing.assert_array_equal(ax.lines[-1].get_xdata(), t.x_sim)
        finally:
            plt.close(fig)

    def test_l_4_shows_bare_unrenormalised_beside_the_solvated_branch(self, mod, laws,
                                                                      inputs):  # noqa: F811
        fig = mod.build_l_4(laws, inputs["abundance"])
        try:
            bare_ax, solv_ax = fig.axes
            ab = inputs["abundance"]
            heights = sorted(round(p.get_height(), 9) for p in bare_ax.patches)
            assert heights == sorted(round(100.0 * f, 9) for f in ab.ion_fraction)
            hatched = [p for p in bare_ax.patches if p.get_hatch()]
            assert len(hatched) == 1 and hatched[0].get_x() < 0 < hatched[0].get_x() + 1
            assert bare_ax.get_xlim()[0] < 0 < solv_ax.get_xlim()[0]
            assert solv_ax.get_yscale() == "linear"          # L-1 (a)
            # per law: an open n = 0 marker, joined to n = 1 by a line
            zero = [ln for ln in bare_ax.lines if list(ln.get_xdata()) == [0]]
            joint = [ln for ln in bare_ax.lines if list(ln.get_xdata()) == [0, 1]]
            assert len(zero) == len(joint) == len(laws)
            assert any("not scored" in t.get_text() for t in bare_ax.texts)
            assert _labels(fig)[0] == "Experiment"
        finally:
            plt.close(fig)

    def _l5(self, mod, laws, inputs):  # noqa: F811
        panels = {law.key: [mod.pv_data(law.view, IHE_KED_DIR, inputs["ked_ref"], n, "3d")
                            for n in mod.FIG_L_5["n_values"]] for law in laws}
        return mod.build_l_5(laws, panels, inputs["ked_ref"], {"symmetric": 2.7})

    def test_l_5_benchmark_only_in_the_bare_panel(self, mod, laws, inputs):  # noqa: F811
        fig = self._l5(mod, laws, inputs)
        try:
            def vlines(ax):
                return {round(ln.get_xdata()[0], 2) for ln in ax.lines
                        if ln.get_xdata()[0] == ln.get_xdata()[-1] and len(ln.get_xdata()) == 2}
            a, b, c = fig.axes
            assert 2.7 in vlines(a)
            assert 2.7 not in (vlines(b) | vlines(c))
            labels = _labels(fig)
            assert "Mean" in labels and "Median" in labels
            assert sum(lab.startswith("Symmetric") for lab in labels) == 1
            assert [ax.get_xlim() for ax in fig.axes] == [
                mod.FIG_L_5["e_xlim_by_n"][n] for n in mod.FIG_L_5["n_values"]]
        finally:
            plt.close(fig)

    def test_l_5_benchmark_switch(self, mod, laws, inputs, monkeypatch):  # noqa: F811
        monkeypatch.setitem(mod.FIG_L_5, "show_symmetric_line", False)
        fig = self._l5(mod, laws, inputs)
        try:
            assert not any(lab.startswith("Symmetric") for lab in _labels(fig))
        finally:
            plt.close(fig)

    @pytest.mark.parametrize("count,boxed", [(1, True), (3, False)])
    def test_l_5_value_box_only_for_few_series(self, mod, laws, inputs, count,
                                               boxed):  # noqa: F811
        fig = self._l5(mod, laws[:count], inputs)
        try:
            for ax in fig.axes:
                assert any(t.get_text().startswith("Mean / median") for t in ax.texts) == boxed
        finally:
            plt.close(fig)
