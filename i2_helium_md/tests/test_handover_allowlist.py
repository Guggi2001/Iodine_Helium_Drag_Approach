"""Temporary handover tooling (removed in P6): import closure and allowlist checks."""

import dataclasses
from pathlib import Path

import pytest

from i2_helium_md.config import SimConfig
from scripts.handover.allowlist import (
    COLUMNS,
    check_coverage,
    decisions,
    draft_rows,
    import_closure,
    keep_violations,
    option_readers,
    path_loads,
    read_rows,
    write_draft,
)

ROOT = Path(__file__).resolve().parents[1]


# --- imports ---------------------------------------------------------------

#: Package modules to be pruned in P5 that kept scripts still reach through the package.
PRUNED_IN_P5 = {"i2_helium_md/physics/exit_strip.py", "i2_helium_md/physics/state_coupling.py",
                "i2_helium_md/sampling/ce_channels.py"}


def kept_closure():
    """Every file the kept scripts (drivers included) import, themselves included."""
    rows = [r for r in read_rows(ROOT / "docs/handover/ALLOWLIST.csv")
            if r["kind"] == "script" and r["decision"] in ("keep", "rewire")
            and r["item"].endswith(".py")]
    assert "scripts/run_thesis_figures.py" in {r["item"] for r in rows}
    return set().union(*(import_closure(ROOT / r["item"], ROOT) for r in rows))


def test_kept_scripts_reach_no_file_to_be_deleted():
    decide = decisions(ROOT / "docs/handover/ALLOWLIST.csv")
    gone = {p for p in kept_closure() if decide(p) in ("move", "delete")}
    assert gone == PRUNED_IN_P5


def test_kept_scripts_reach_no_study_generator():
    assert {p for p in kept_closure() if Path(p).name.startswith("gen_")} == {
        "scripts/gen_legacy_runs.py", "scripts/gen_thesis_laws.py",
        "scripts/gen_tier0_runs.py"}


def test_closure_reaches_scripts_and_package():
    c = import_closure(ROOT / "scripts/post_processing/tier2_b2_thesis_figure.py", ROOT)
    assert "scripts/post_processing/anchors.py" in c
    assert any(p.startswith("i2_helium_md/") for p in c)


def test_closure_follows_function_level_and_bare_imports(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/a.py").write_text("def f():\n    from scripts import b\n")
    (tmp_path / "scripts/b.py").write_text("import c\n")
    (tmp_path / "scripts/c.py").write_text("")
    assert import_closure(tmp_path / "scripts/a.py", tmp_path) == {
        "scripts/a.py", "scripts/b.py", "scripts/c.py"}


def test_package_reexport_reaches_only_the_defining_module(tmp_path):
    pkg = tmp_path / "i2_helium_md"
    pkg.mkdir()
    (pkg / "__init__.py").write_text("from .x import X\nfrom .y import Y\n")
    (pkg / "x.py").write_text("X = 1\n")
    (pkg / "y.py").write_text("Y = 1\n")
    (tmp_path / "e.py").write_text("from i2_helium_md import X\n")
    assert import_closure(tmp_path / "e.py", tmp_path) == {
        "e.py", "i2_helium_md/__init__.py", "i2_helium_md/x.py"}


# --- draft -----------------------------------------------------------------

def test_draft_lists_figures_scripts_and_options():
    rows = {r["item"]: r for r in draft_rows(ROOT)}
    figures = [r for r in rows.values() if r["kind"] == "figure"]
    assert len(figures) == 43 and all(r["decision"] == "keep" for r in figures)
    for fig in figures:
        script = fig["inputs"].split(";")[0]
        assert Path(fig["item"]).name in rows[script]["needed_by"].split()
    fields = {f.name for f in dataclasses.fields(SimConfig)}
    assert {k for k in rows if k.startswith("SimConfig.")} == {f"SimConfig.{f}" for f in fields}


@pytest.mark.skipif(not (ROOT / "data/runs/B2_N1000_s1/cfg.json").is_file(),
                    reason="thesis runs not on disk")
def test_option_rows_show_the_thesis_values():
    rows = {r["item"]: r for r in draft_rows(ROOT)}
    inputs = rows["SimConfig.drag_form"]["inputs"]
    assert "B2=linear_quadratic" in inputs and "h405=capped_cubic" in inputs


def test_write_draft_never_overwrites(tmp_path):
    path = tmp_path / "ALLOWLIST.csv"
    path.write_text("decided by hand")
    with pytest.raises(FileExistsError):
        write_draft(path, [])
    assert path.read_text() == "decided by hand"


# --- coverage --------------------------------------------------------------

def _write(path, rows):
    write_draft(path, [dict(zip(COLUMNS, (i, k, "", "", d, ""))) for i, k, d in rows])
    return path


def test_reports_uncovered_vanished_duplicate_and_undecided(tmp_path):
    path = _write(tmp_path / "a.csv", [
        ("a.py", "script", "keep"), ("a.py", "script", "delete"), ("gone.py", "script", "keep"),
        ("d/", "doc", ""), ("SimConfig.x", "option", "keep")])
    problems = check_coverage(path, ["a.py", "new.py", "d/r.md"], ["x", "y"])
    assert sorted(problems) == ["duplicate: a.py", "uncovered: SimConfig.y", "uncovered: new.py",
                                "undecided: d/", "vanished: gone.py"]


def test_unknown_decision_and_kind_are_reported(tmp_path):
    path = _write(tmp_path / "a.csv", [("a.py", "script", "maybe"), ("b.py", "thing", "keep")])
    assert check_coverage(path, ["a.py", "b.py"], []) == [
        "bad decision 'maybe': a.py", "bad kind 'thing': b.py"]


def test_keep_violations_name_direct_importers(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/k.py").write_text("import old\n")
    (tmp_path / "scripts/old.py").write_text("")
    path = _write(tmp_path / "a.csv", [("scripts/k.py", "script", "keep"),
                                       ("scripts/", "script", "delete")])
    assert keep_violations(path, tmp_path) == [("scripts/k.py", "scripts/old.py")]


def test_kept_code_needs_only_kept_code():
    """A file that imports deleted code must be marked ``rewire``."""
    path = ROOT / "docs/handover/ALLOWLIST.csv"
    decide = decisions(path)
    assert [p for p in keep_violations(path, ROOT) if decide(p[0]) == "keep"] == []


def test_moved_code_is_checked_for_deleted_imports(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/m.py").write_text("import old\n")
    (tmp_path / "scripts/old.py").write_text("")
    path = _write(tmp_path / "a.csv", [("scripts/m.py", "script", "move"),
                                       ("scripts/old.py", "script", "delete")])
    assert keep_violations(path, tmp_path) == [("scripts/m.py", "scripts/old.py")]


def test_path_loads_count_as_imports(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/k.py").write_text("import importlib\np = 'scripts/old.py'\nq = 'old2'\n")
    (tmp_path / "scripts/plain.py").write_text("p = 'old2'\n")
    (tmp_path / "scripts/old.py").write_text("")
    (tmp_path / "scripts/old2.py").write_text("")
    assert path_loads(tmp_path / "scripts/k.py", tmp_path) == {
        tmp_path / "scripts/old.py", tmp_path / "scripts/old2.py"}
    assert path_loads(tmp_path / "scripts/plain.py", tmp_path) == set()


def test_option_readers_name_kept_files_that_read_a_removed_field(tmp_path):
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/k.py").write_text("getattr(cfg, 'gone')\ncfg.stays\n")
    path = _write(tmp_path / "a.csv", [
        ("scripts/k.py", "script", "keep"), ("SimConfig.gone", "option", "delete"),
        ("SimConfig.stays", "option", "keep")])
    assert option_readers(path, tmp_path) == [("scripts/k.py", "gone")]


def test_kept_code_reads_only_kept_options():
    assert option_readers(ROOT / "docs/handover/ALLOWLIST.csv", ROOT) == []


def test_moved_code_names_what_it_still_needs():
    """A moved script that imports deleted code says so in its note."""
    path = ROOT / "docs/handover/ALLOWLIST.csv"
    note = {r["item"]: r["note"] for r in read_rows(path)}
    decide = decisions(path)
    missing = [(a, b) for a, b in keep_violations(path, ROOT)
               if decide(a) == "move" and Path(b).stem not in note[a]]
    assert missing == []
