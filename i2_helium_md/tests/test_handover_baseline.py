"""Temporary handover tooling (removed in P6): figure baseline against git HEAD."""

import io
import subprocess
from pathlib import Path

import pytest
from PIL import Image

from scripts.handover.baseline import (
    FigureRow,
    assert_clean,
    parse_figure_table,
    render_all,
    same_content,
)

ROOT = Path(__file__).resolve().parents[1]


def _row(figure, script):
    return FigureRow(figure=Path(figure), script=Path(script), inputs="")


@pytest.fixture
def tmp_git_repo(tmp_path):
    (tmp_path / "data").mkdir()
    (tmp_path / "data/x.csv").write_bytes(b"a\r\n1\r\n")
    (tmp_path / "s_bad.py").write_text("raise SystemExit('boom')\n")
    (tmp_path / "s_noop.py").write_text("")
    (tmp_path / "s_ok.py").write_text(
        "open('data/x.csv', 'wb').write(bytes([97, 10, 49, 10]))")
    (tmp_path / "s_change.py").write_text(
        "open('data/x.csv', 'wb').write(bytes([97, 10, 50, 10]))")
    for args in (["init", "-q"], ["add", "-A"],
                 ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                  "commit", "-qm", "t"]):
        subprocess.run(["git", *args], cwd=tmp_path, check=True, capture_output=True)
    return tmp_path


def test_readme_lists_every_thesis_figure():
    rows = parse_figure_table(ROOT / "data/figures/README.md")
    on_disk = {p.relative_to(ROOT).as_posix()
               for p in (ROOT / "data/figures").glob("thesis_*/*")}
    assert len(rows) == 43
    assert {r.figure.as_posix() for r in rows} == on_disk
    assert all((ROOT / r.script).is_file() for r in rows)


def test_csv_line_endings_do_not_count():
    assert same_content(b"a,b\r\n1,2\r\n", b"a,b\n1,2\n", ".csv")
    assert not same_content(b"a,b\n1,3\n", b"a,b\n1,2\n", ".csv")


def test_png_compares_pixels_not_metadata():
    im = Image.new("RGB", (4, 4), "red")

    def png(dpi):
        buf = io.BytesIO()
        im.save(buf, "PNG", dpi=(dpi, dpi))
        return buf.getvalue()

    assert png(72) != png(200)
    assert same_content(png(72), png(200), ".png")


def test_assert_clean_refuses_dirty_data(tmp_git_repo):
    (tmp_git_repo / "data/x.csv").write_text("changed")
    with pytest.raises(RuntimeError, match="x.csv"):
        assert_clean(tmp_git_repo)


def test_crashing_script_is_recorded_and_others_continue(tmp_git_repo):
    out = render_all([_row("data/y.csv", "s_bad.py"), _row("data/x.csv", "s_ok.py")],
                     tmp_git_repo)
    assert out["data/y.csv"]["status"] == "not_rendered"
    assert "boom" in out["data/y.csv"]["reason"]
    assert out["data/x.csv"]["status"] == "identical"


def test_stale_figure_is_not_called_identical(tmp_git_repo):
    out = render_all([_row("data/x.csv", "s_noop.py")], tmp_git_repo)
    assert out["data/x.csv"]["status"] == "not_rendered"


def test_changed_figure_differs_and_bytes_are_restored(tmp_git_repo):
    out = render_all([_row("data/x.csv", "s_change.py")], tmp_git_repo)
    assert out["data/x.csv"]["status"] == "differs"
    assert (tmp_git_repo / "data/x.csv").read_bytes() == b"a\r\n1\r\n"
