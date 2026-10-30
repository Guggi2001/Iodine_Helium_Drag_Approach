"""``tier0_law_figure.py`` — the thesis-style Tier-0 law figure (plot-only).

Runs no MD: it reads the committed Tier-0 series / scores CSVs
(``data/figures/thesis_tier0/``), the TDDFT references and the drag bundles.
Locks the anchor guard (B2 RMSE / 9 Å margin bit-exact vs the committed G11
JSONs), the variant table, the bulk-law coefficients and that the builder
draws one trace per law in every panel. Agg backend, figures closed, nothing
written to disk.
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

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "tier0_law_figure.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("tier0_law_figure_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def data(mod):
    from i2_helium_md.postprocess import load_hedft_trajectory
    from i2_helium_md.postprocess.hedft_loader import load_smoothed_speed_reference
    return dict(
        series=mod.read_series(), scores=mod.read_scores(), coeffs=mod.law_coefficients(),
        ref18=load_smoothed_speed_reference(
            mod.REFERENCE_DRAG_ROOT / "18A" / "velocity_smoothed" / mod.REFERENCE_FILENAME),
        hedft18=load_hedft_trajectory(mod.HEDFT_18A_PATH),
        hedft9=load_hedft_trajectory(mod.HEDFT_9A_PATH),
    )


class TestData:
    def test_variants(self, mod):
        assert mod.VARIANTS == {"tier2_B2_h405": ("B2", "h405"),
                                "tier2_B2_B3_h405": ("B2", "B3", "h405")}

    def test_exit_bound_is_the_committed_atom1_exit_speed(self, mod, data):
        bound = mod.exit_bound_Aps(data["scores"], ("B2", "B3", "h405"))
        assert bound == float(data["scores"]["B2"]["v_exit_atom1_Aps"])
        # margin = law exit mean - bound (the pass criterion the figure annotates)
        for k in ("B2", "B3", "h405"):
            row = data["scores"][k]
            assert float(row["v_exit_9A_Aps"]) - bound == pytest.approx(
                float(row["margin_9A_Aps"]), abs=1e-12)

    def test_anchor_passes_on_committed_csvs(self, mod, data):
        mod.check_anchor(data["scores"], data["series"])

    def test_anchor_rejects_a_drifted_b2_rmse(self, mod, data):
        scores = {k: dict(v) for k, v in data["scores"].items()}
        scores["B2"]["rmse_18A_Aps"] = str(float(scores["B2"]["rmse_18A_Aps"]) + 1e-12)
        with pytest.raises(ValueError, match="RMSE"):
            mod.check_anchor(scores, data["series"])

    def test_coefficients_match_the_scored_laws(self, mod, data):
        # B2 / B3: lq at the LQ_LAWS (a, c); h405: capped cubic at the production cap
        for label, (a, c, _, _) in mod.LQ_LAWS.items():
            assert data["coeffs"][label].coefficients == {"a": a, "c": c}
        h = data["coeffs"]["h405"].coefficients
        assert h["v_c"] == 5.40
        assert f"b {h['b']:.3f}" in data["scores"]["h405"]["drag"]

    def test_bulk_force_is_ungated(self, mod, data):
        v = np.array([1.0, 3.0])
        a, c = mod.LQ_LAWS["B2"][:2]
        f = mod.drag_force(v, np.full_like(v, mod.BULK_DEPTH_A), data["coeffs"]["B2"], mod.GATE_A)
        np.testing.assert_allclose(f, a * v + c * v * v, rtol=1e-12)


class TestBuilder:
    @pytest.mark.parametrize("keys", [("B2", "h405"), ("B2", "B3", "h405")])
    def test_one_trace_per_law_in_every_panel(self, mod, data, keys):
        fig = mod.build_tier0(keys, data["series"], data["scores"], data["coeffs"],
                              data["ref18"], data["hedft18"], data["hedft9"])
        try:
            assert fig.get_size_inches()[1] == pytest.approx(
                mod.FIG_T0["height_by_nlaws"][len(keys)])
            ax18, ax9, axF = fig.axes
            assert len(ax18.lines) == 3 + len(keys)      # raw v1, raw v2, smoothed
            assert len(ax9.lines) == 3 + len(keys)       # raw atoms 1 / 2, exit bound
            assert len(axF.lines) == len(keys)
            # (b) draws the TDDFT atoms exactly like (a)'s raw traces (one legend entry)
            raw18, raw9 = ax18.lines[0], ax9.lines[0]
            assert (raw18.get_color(), raw18.get_linewidth()) == (raw9.get_color(),
                                                                  raw9.get_linewidth())
            labels = [t.get_text() for leg in fig.legends for t in leg.get_texts()]
            assert labels[:3] == [r"TDDFT $|v_1|$, $|v_2|$", r"Smoothed $|v_1|$",
                                  "Exit-speed bound"]
            assert labels[3:] == [mod.LAWS[k][1] for k in keys]
        finally:
            plt.close(fig)

    def test_figure_spans_the_text_width(self, mod):
        for fig_style in (mod.FIG_T0, mod.FIG_RIDGE):
            assert fig_style["figsize"][0] == pytest.approx(mod.TEXTWIDTH_IN)


class TestRidge:
    def test_grid_rows_and_b2_check(self, mod, data):
        ridge = mod.read_ridge(data["coeffs"]["B2"])
        cs = [c for c, _, _ in ridge]
        assert cs == sorted(cs) and 1.6 in cs
        # the 9 Å bound cuts the steepest point only (c 7.224, THESIS_DISCUSSION §0)
        assert [c for c, _, ok in ridge if not ok] == [7.224]

    def test_rejects_a_b2_that_is_not_the_ridge_b2(self, mod, data):
        from dataclasses import replace
        b2 = data["coeffs"]["B2"]
        off = replace(b2, coefficients={**b2.coefficients, "a": b2.coefficients["a"] + 1e-9})
        with pytest.raises(ValueError, match="B2 row"):
            mod.read_ridge(off)

    def test_draws_every_ridge_point_once_per_panel(self, mod, data):
        b2 = data["coeffs"]["B2"]
        ridge = mod.read_ridge(b2)
        fig = mod.build_ridge(ridge, b2, data["ref18"])
        try:
            ax, axz, cbar = fig.axes
            for a in (ax, axz):
                assert len(a.lines) == len(ridge)        # c = 1.6 drawn as B2, not twice
            labels = [t.get_text() for t in ax.get_legend().get_texts()]
            assert labels == ["18 Å trace speeds", "Fails 9 Å bound", "LQ ($c$ = 1.6)"]
        finally:
            plt.close(fig)
