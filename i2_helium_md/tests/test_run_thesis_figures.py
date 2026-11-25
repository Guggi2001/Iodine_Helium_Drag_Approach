"""Figure driver, on a small fake project: copies, restores, reports."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from scripts import run_thesis_figures as driver
from scripts.post_processing.anchors import PREFIX

FIGURES = {"good": "a.csv", "bad": "b.csv", "second": "c.csv", "missing": "d.csv",
           "anchor": "e.csv"}
SOURCES = {
    "good": "open('data/figures/ch/a.csv', 'wb').write(b'x\\n1\\n')",
    "bad": ("open('data/figures/ch/b.csv', 'wb').write(b'half')\n"
            "open('data/figures/ch/input.csv', 'wb').write(b'half')\n"
            "raise SystemExit('boom')"),
    "second": ("open('data/figures/ch/c.csv', 'wb').write("
               "open('data/figures/ch/input.csv', 'rb').read())"),
    "missing": "raise FileNotFoundError('runs/gone is missing; make it with scripts/run_x.py')",
    "anchor": ("open('data/figures/ch/e.csv', 'wb').write(b'old\\n')\n"
               f"print('{PREFIX}B2 nbar_det: 1 vs 2')"),
    "table": "open('data/runs/h2b_forward_model/t.csv', 'wb').write(b'changed')",
    # sees the committed table, not the one another script has just written
    "reader": "print(open('data/runs/h2b_forward_model/t.csv').read())\n"
              "assert open('data/runs/h2b_forward_model/t.csv').read() == 'old\\n'",
}


@pytest.fixture
def project(tmp_path, monkeypatch):
    root = tmp_path / "project"
    figures, tables = root / "data/figures/ch", root / "data/runs/h2b_forward_model"
    for folder in (figures, tables, root / "scripts"):
        folder.mkdir(parents=True)
    for name in (*FIGURES.values(), "input.csv"):
        (figures / name).write_bytes(b"old\r\n")
    (tables / "t.csv").write_bytes(b"old\n")
    for stem, source in SOURCES.items():
        (root / f"scripts/{stem}.py").write_text(source)
    (root / "data/figures/README.md").write_text("".join(
        f"| [`{f}`](ch/{f}) | [`{s}.py`](../../scripts/{s}.py) | x |\n"
        for s, f in FIGURES.items()))
    for args in (["init", "-q"], ["add", "-A"],
                 ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                  "commit", "-qm", "t"]):
        subprocess.run(["git", *args], cwd=root, check=True, capture_output=True)
    for name, value in dict(PROJECT_ROOT=root, README=root / "data/figures/README.md",
                            OUTPUT_ROOT=root / "data/figures_regenerated", SCRIPTS="all",
                            TABLE_SCRIPTS=()).items():
        monkeypatch.setattr(driver, name, value)
    return root


def run(monkeypatch, *stems):
    monkeypatch.setattr(driver, "SCRIPTS", list(stems))
    return driver.main()


def git_status(root):
    return subprocess.run(["git", "status", "--short", "--untracked-files=no"], cwd=root,
                          capture_output=True, text=True).stdout


def test_figure_is_copied_and_the_committed_one_restored(project, monkeypatch, capsys):
    assert run(monkeypatch, "good") == 0
    assert (project / "data/figures_regenerated/ch/a.csv").read_bytes() == b"x\n1\n"
    assert (project / "data/figures/ch/a.csv").read_bytes() == b"old\r\n"
    assert git_status(project) == ""
    assert "differs" in capsys.readouterr().out


def test_failing_script_is_named_and_the_next_one_runs_on_restored_inputs(
        project, monkeypatch, capsys):
    assert run(monkeypatch, "bad", "second") == 1
    out = capsys.readouterr().out
    assert "not rendered" in out and "boom" in out
    assert (project / "data/figures_regenerated/ch/c.csv").read_bytes() == b"old\r\n"
    assert not (project / "data/figures_regenerated/ch/b.csv").exists()
    assert git_status(project) == ""


def test_dirty_figure_is_refused_before_any_script(project, monkeypatch):
    path = project / "data/figures/ch/b.csv"
    path.write_text("mine")
    with pytest.raises(RuntimeError, match="b.csv"):
        run(monkeypatch, "good")
    assert path.read_text() == "mine"
    assert not (project / "data/figures_regenerated").exists()


def test_edited_readme_does_not_block(project, monkeypatch):
    readme = project / "data/figures/README.md"
    readme.write_text(readme.read_text() + "\nmore text\n")
    assert run(monkeypatch, "good") == 0
    assert readme.read_text().endswith("more text\n")


def test_unknown_stem_lists_the_choices(project, monkeypatch):
    with pytest.raises(KeyError, match="second"):
        run(monkeypatch, "nope")


def test_missing_run_shows_its_message(project, monkeypatch, capsys):
    assert run(monkeypatch, "missing") == 1
    assert "make it with scripts/run_x.py" in capsys.readouterr().out


def test_anchor_difference_is_repeated_and_fails(project, monkeypatch, capsys):
    assert run(monkeypatch, "anchor") == 1
    out = capsys.readouterr().out
    assert "identical" in out
    assert "B2 nbar_det: 1 vs 2" in out


def test_all_runs_every_script(project, capsys):
    assert driver.main() == 1
    out = capsys.readouterr().out
    assert all(name in out for name in FIGURES.values())


def test_tables_are_copied_apart_and_restored_after_each_script(project, monkeypatch, capsys):
    monkeypatch.setattr(driver, "TABLE_SCRIPTS",
                        (Path("scripts/table.py"), Path("scripts/reader.py")))
    assert run(monkeypatch, "table", "reader") == 0
    assert (project / "data/figures_regenerated/tables/t.csv").read_bytes() == b"changed"
    assert (project / "data/runs/h2b_forward_model/t.csv").read_bytes() == b"old\n"
    assert "t.csv" in capsys.readouterr().out


def test_failing_table_script_fails_the_driver(project, monkeypatch, capsys):
    (project / "scripts/table.py").write_text("raise SystemExit('no table')")
    monkeypatch.setattr(driver, "TABLE_SCRIPTS", (Path("scripts/table.py"),))
    assert run(monkeypatch, "table") == 1
    assert "no table" in capsys.readouterr().out


def test_stale_copy_is_removed_when_the_figure_is_not_rendered(project, monkeypatch):
    stale = project / "data/figures_regenerated/ch/b.csv"
    stale.parent.mkdir(parents=True)
    stale.write_text("from an earlier run")
    assert run(monkeypatch, "bad") == 1
    assert not stale.exists()


def test_table_scripts_exist():
    root = Path(__file__).resolve().parents[1]
    assert driver.TABLE_SCRIPTS and all((root / s).is_file() for s in driver.TABLE_SCRIPTS)


def test_one_failed_restore_does_not_leave_the_others_modified(project, monkeypatch):
    """E.g. a locked file: the rest is put back, the error names the one that is not."""
    (project / "scripts/good.py").write_text(
        "for f in 'abc':\n    open(f'data/figures/ch/{f}.csv', 'wb').write(b'new')")
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.invalid",
                    "commit", "-qam", "t"], cwd=project, check=True, capture_output=True)
    write = Path.write_bytes

    def locked(self, data):
        if self == project / "data/figures/ch/b.csv":
            raise PermissionError("locked")
        return write(self, data)

    monkeypatch.setattr(Path, "write_bytes", locked)
    with pytest.raises(RuntimeError, match="b.csv"):
        run(monkeypatch, "good")
    assert (project / "data/figures/ch/a.csv").read_bytes() == b"old\r\n"
    assert (project / "data/figures/ch/c.csv").read_bytes() == b"old\r\n"
    assert (project / "data/figures/ch/b.csv").read_bytes() == b"new"


def test_readme_edited_during_the_run_is_kept(project, monkeypatch):
    (project / "scripts/good.py").write_text(
        SOURCES["good"] + "\nopen('data/figures/README.md', 'a').write('mine')")
    assert run(monkeypatch, "good") == 0
    assert (project / "data/figures/README.md").read_text().endswith("mine")


def test_deleted_document_does_not_block(project, monkeypatch):
    (project / "data/notes.md").write_text("x")
    for args in (["add", "-A"], ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                                 "commit", "-qm", "t"]):
        subprocess.run(["git", *args], cwd=project, check=True, capture_output=True)
    (project / "data/notes.md").unlink()
    assert run(monkeypatch, "good") == 0


def test_unreadable_image_is_not_rendered_and_the_next_script_runs(
        project, monkeypatch, capsys):
    png = project / "data/figures/ch/p.png"
    driver.Image.new("RGB", (2, 2)).save(png)
    (project / "scripts/broken.py").write_text(
        "open('data/figures/ch/p.png', 'wb').write(b'not an image')")
    readme = project / "data/figures/README.md"
    readme.write_text("| [`p.png`](ch/p.png) | [`broken.py`](../../scripts/broken.py) | x |\n"
                      + readme.read_text())
    for args in (["add", "-A"], ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                                 "commit", "-qm", "t"]):
        subprocess.run(["git", *args], cwd=project, check=True, capture_output=True)
    assert run(monkeypatch, "broken", "good") == 1
    out = capsys.readouterr().out
    assert "not rendered  data/figures/ch/p.png" in out
    assert (project / "data/figures_regenerated/ch/a.csv").exists()
    assert not (project / "data/figures_regenerated/ch/p.png").exists()
    assert git_status(project) == ""
