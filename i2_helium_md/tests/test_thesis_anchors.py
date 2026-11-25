"""The thesis runs on disk reproduce the committed numbers the figures quote.

The figure scripts only report a difference; the exact checks are here.
"""

from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md.simulation.thesis_runs import RUNS_ROOT, TABLES_ROOT, member_dir, pool_dir
from scripts.post_processing import (
    anchors,
    tier0_law_figure,
    tier0_thesis_figure,
    tier2_anisotropy_figure,
    tier2_b2_thesis_figure,
    tier2_g13_shed_timing,
    tier2_law_comparison_figure,
    tier2_thesis_observable_table,
)
from scripts.tier0_runs import tier0_runs

LAWS = ("B2", "B3", "h405")
pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")
needs_pools = pytest.mark.skipif(
    not all((pool_dir(law) / "detection.npz").is_file() for law in LAWS),
    reason="thesis pools not on disk")


def copy_with(path, tmp_path, label, column, value):
    """Copy of a committed table with one cell changed."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    hit = [r for r in rows if r["label"] == label]
    assert len(hit) == 1 and column in hit[0]
    hit[0][column] = value
    out = tmp_path / path.name
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    return out


class TestAnchors:
    @pytest.fixture
    def table(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("label,container,nbar_det,midHot\nB2,c,3.90284712,1.0016225\n")
        return path

    def test_equal_to_4_dp_is_no_difference(self, table):
        values = {"nbar_det": 3.90284, "midHot": 1.00162}
        assert anchors.row_differences(values, table, "B2", "c") == []

    def test_a_difference_names_table_law_and_column(self, table):
        bad = anchors.row_differences({"nbar_det": 3.9040, "midHot": 1.00162}, table, "B2")
        assert bad == ["t.csv B2 nbar_det: 3.9040 vs 3.9028"]

    def test_nan_is_a_difference(self, table):
        assert anchors.row_differences({"nbar_det": float("nan")}, table, "B2")

    def test_wrong_container_or_missing_row_raises(self, table):
        with pytest.raises(ValueError, match="belongs to"):
            anchors.row_differences({"nbar_det": 3.9}, table, "B2", "other")
        with pytest.raises(ValueError, match="expected one row"):
            anchors.row_differences({"nbar_det": 3.9}, table, "B3")

    def test_report_prints_each_difference_or_ok(self, capsys):
        anchors.report(["a", "b"], "fine")
        assert capsys.readouterr().out.splitlines() == [
            anchors.PREFIX + "a", anchors.PREFIX + "b"]
        anchors.report([], "fine")
        assert capsys.readouterr().out == "fine\n"


# --- a changed committed value is reported ------------------------------------

def test_g13_reports_law_and_column(tmp_path, monkeypatch):
    row = anchors.committed_row(TABLES_ROOT / "g13_decoupled_table.csv", "B2")
    n_det = np.array([0.0, 1.0])
    read = SimpleNamespace(scored=np.array([True, True]), n_det=n_det,
                           ke_det_eV=np.array([0.0, float(row["KE1_mean"])]))
    bad = tier2_g13_shed_timing.anchor_differences("B2", read)
    assert any("B2 num_scored" in text for text in bad)
    assert any("B2 nbar_det" in text for text in bad)
    assert not any("KE1_mean" in text for text in bad)


def test_anisotropy_reports_the_size():
    series = tier2_anisotropy_figure.BetaSeries(n=np.array([0, 1]), beta=np.array([1.88, 2.1]))
    assert tier2_anisotropy_figure.anchor_differences(series, {0: 1.88, 1: 2.04, 2: 1.97}) == [
        "THESIS_DISCUSSION 13.3 beta n = 1: 2.1 vs 2.04",
        "THESIS_DISCUSSION 13.3 beta n = 2: None vs 1.97"]


def test_tier0_thesis_reports_law_and_quantity():
    def runs(rmse, margin):
        return {"B2": {"rmse": rmse}}, {"B2": {"row": {"margin_Aps": margin}}}

    rmse, margin = tier0_thesis_figure.committed_b2()
    assert tier0_thesis_figure.anchor_differences(*runs(rmse, margin)) == []
    bad = tier0_thesis_figure.anchor_differences(*runs(rmse * (1 + 1e-15), margin + 1.0))
    assert len(bad) == 2 and all("B2" in text for text in bad)
    assert "RMSE" in bad[0] and "margin" in bad[1]


def references(fig):
    return (fig.load_he_abundance_reference(fig.ABUNDANCE_REF_PATH),
            fig.load_ihe_ked_reference(fig.IHE_KED_REFERENCE_DIR / "IHe_KED_reference.csv"))


def b2_differences():
    fig = tier2_b2_thesis_figure
    abundance, ked = references(fig)
    _, read = fig.load_pool()
    return fig.anchor_differences(fig.size_data(read, abundance), fig.ke_data(read, ked))


def law_differences(law):
    fig = tier2_law_comparison_figure
    return fig.anchor_differences(fig.load_law(law, *references(fig)))


# --- the stored runs reproduce the committed numbers ----------------------------

@needs_pools
def test_b2_figure():
    assert b2_differences() == []


@needs_pools
def test_b2_figure_reports_a_changed_value(tmp_path, monkeypatch):
    changed = copy_with(tier2_b2_thesis_figure.WINNERS_CSV, tmp_path, "B2", "midHot", "9.0")
    monkeypatch.setattr(tier2_b2_thesis_figure, "WINNERS_CSV", changed)
    (text,) = b2_differences()
    assert "B2 midHot" in text


@needs_pools
@pytest.mark.parametrize("law", LAWS)
def test_law_comparison(law):
    assert law_differences(law) == []


@needs_pools
def test_law_comparison_reports_a_changed_value(tmp_path, monkeypatch):
    changed = copy_with(tier2_law_comparison_figure.EXTRA_CSV, tmp_path, "B3", "bare_share", "9")
    monkeypatch.setattr(tier2_law_comparison_figure, "EXTRA_CSV", changed)
    (text,) = law_differences("B3")
    assert "B3 bare_share" in text


@needs_pools
def test_anisotropy():
    fig = tier2_anisotropy_figure
    model = fig.model_beta(fig.RUN_DIR, fig.N_MAX_MODEL, fig.MIN_COUNT)
    assert fig.anchor_differences(model, fig.ANCHOR_BETA_MODEL) == []


@needs_pools
def test_observable_table():
    table = tier2_thesis_observable_table
    per_n, extra = table.build_tables()
    assert table.anchor_differences(per_n, extra) == []


@needs_pools
def test_observable_table_reports_a_changed_value(tmp_path, monkeypatch):
    table = tier2_thesis_observable_table
    changed = copy_with(table.SCORE_TABLE_CSV, tmp_path, "h405", "KE2_med", "9.0")
    monkeypatch.setattr(table, "SCORE_TABLE_CSV", changed)
    (text,) = table.anchor_differences(*table.build_tables())
    assert "h405 KE2_med" in text


@pytest.mark.skipif(
    not all((member_dir(law, m) / "stage1_events.npz").is_file()
            for law in LAWS for m in ("s1", "s6")), reason="thesis members not on disk")
@pytest.mark.parametrize("law", LAWS)
def test_g13(law):
    timing = tier2_g13_shed_timing
    assert timing.anchor_differences(law, timing.read_law(law)) == []


@pytest.mark.skipif(not all((RUNS_ROOT / name / "ion.npz").is_file() for name in tier0_runs()),
                    reason="Tier-0 runs not on disk")
def test_tier0_thesis():
    runs18, _ = tier0_thesis_figure.run_18a()
    runs9, _ = tier0_thesis_figure.run_9a()
    assert tier0_thesis_figure.anchor_differences(runs18, runs9) == []


def test_tier0_law_figure():
    fig = tier0_law_figure
    assert fig.anchor_differences(fig.read_scores(), fig.read_series()) == []


def test_missing_run_names_the_driver(tmp_path, monkeypatch):
    """Scripts that read one run directory name the driver that makes it."""
    from scripts.post_processing import legacy_ke_size_figure, legacy_ked_shapes_figure

    for mod, driver in ((tier2_b2_thesis_figure, "gen_thesis_laws.py"),
                        (tier2_anisotropy_figure, "gen_thesis_laws.py"),
                        (legacy_ke_size_figure, "gen_legacy_runs.py"),
                        (legacy_ked_shapes_figure, "gen_legacy_runs.py")):
        monkeypatch.setattr(mod, "RUN_DIR", tmp_path / "gone")
        with pytest.raises(FileNotFoundError, match=driver):
            mod.main()


def test_missing_member_or_pool_names_the_driver(tmp_path, monkeypatch):
    from scripts.post_processing import legacy_continuum_drag_figure, tier2_g6_chord_map

    gone = tmp_path / "runs"
    for mod in (tier2_g13_shed_timing, tier2_g6_chord_map):
        monkeypatch.setattr(mod, "member_dir", lambda law, m: gone / f"{law}_{m}")
        with pytest.raises(FileNotFoundError, match="gen_thesis_laws.py"):
            mod.main()
    monkeypatch.setattr(tier2_thesis_observable_table, "pool_dir", lambda law: gone / law)
    with pytest.raises(FileNotFoundError, match="gen_thesis_laws.py"):
        tier2_thesis_observable_table.main()
    monkeypatch.setattr(tier2_law_comparison_figure, "RUNS_DIR", gone)
    with pytest.raises(FileNotFoundError, match="gen_thesis_laws.py"):
        tier2_law_comparison_figure.main()
    monkeypatch.setattr(legacy_continuum_drag_figure, "RUNS_ROOT", gone)
    with pytest.raises(FileNotFoundError, match="gen_legacy_runs.py"):
        legacy_continuum_drag_figure.main()


def test_unfinished_run_names_the_driver(tmp_path, monkeypatch):
    """A directory without its result file, as an interrupted run leaves it."""
    from scripts.post_processing import (
        legacy_continuum_drag_figure,
        legacy_ke_size_figure,
        legacy_ked_shapes_figure,
        tier2_g6_chord_map,
    )

    empty = tmp_path / "unfinished"
    empty.mkdir()
    laws, legacy = "gen_thesis_laws.py", "gen_legacy_runs.py"
    monkeypatch.setattr(tier2_law_comparison_figure, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(legacy_continuum_drag_figure, "RUNS_ROOT", tmp_path)
    for name in ("B2_N6000_pooled", "legacy_9A_N2000"):
        (tmp_path / name).mkdir()
    for mod in (tier2_g13_shed_timing, tier2_g6_chord_map):
        monkeypatch.setattr(mod, "member_dir", lambda law, m: empty)
    monkeypatch.setattr(tier2_thesis_observable_table, "pool_dir", lambda law: empty)
    for mod in (tier2_b2_thesis_figure, tier2_anisotropy_figure, legacy_ke_size_figure,
                legacy_ked_shapes_figure):
        monkeypatch.setattr(mod, "RUN_DIR", empty)
    for mod, driver in ((tier2_g13_shed_timing, laws), (tier2_g6_chord_map, laws),
                        (tier2_thesis_observable_table, laws), (tier2_b2_thesis_figure, laws),
                        (tier2_anisotropy_figure, laws), (tier2_law_comparison_figure, laws),
                        (legacy_ke_size_figure, legacy), (legacy_ked_shapes_figure, legacy),
                        (legacy_continuum_drag_figure, legacy)):
        with pytest.raises(FileNotFoundError, match=driver):
            mod.main()
