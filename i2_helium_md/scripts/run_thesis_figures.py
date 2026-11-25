"""Draw all thesis figures and tables and compare them with the committed ones.

Set the options below, then: ``python scripts/run_thesis_figures.py``.

Needs the runs of ``gen_thesis_laws.py``, ``gen_legacy_runs.py`` and ``gen_tier0_runs.py``. The
figure scripts are the ones in ``data/figures/README.md``; they write in place, so
``data/`` must be clean, and every tracked file there is put back afterwards. What
they wrote goes to ``OUTPUT_ROOT``.

Reports each figure as identical, differs or not rendered, and every quoted number
that differs from the committed one. Exits 1 if a script failed, a figure was not
rendered or a number differs. A figure that only ``differs`` is no failure: the
legacy runs are seeded now, and pixels depend on the matplotlib version.
"""

from __future__ import annotations

import io
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.post_processing.anchors import PREFIX  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================
SCRIPTS: str | list[str] = "all"             # or script names without .py
OUTPUT_ROOT = PROJECT_ROOT / "data" / "figures_regenerated"
# =============================================================================

README = PROJECT_ROOT / "data" / "figures" / "README.md"
#: Scripts that write only tables; they run after the figure scripts.
TABLE_SCRIPTS = (Path("scripts/post_processing/tier2_thesis_observable_table.py"),)
TABLES_DIR = Path("data/runs/h2b_forward_model")
SCRIPT_TIMEOUT_S = 3600

_ROW = re.compile(r"^\|\s*\[`[^`]+`\]\(([^)]+)\)\s*\|\s*\[`[^`]+`\]\(([^)]+)\)\s*\|\s*(.*?)\s*\|$")


@dataclass(frozen=True)
class FigureRow:
    """One thesis figure: paths relative to the project root."""

    figure: Path
    script: Path
    inputs: str


@dataclass
class FigureResult:
    """Outcome of one figure; ``old`` / ``new`` are file contents or None."""

    script: Path
    status: str = "not_rendered"             # "identical" | "differs" | "not_rendered"
    reason: str = ""
    old: bytes | None = None
    new: bytes | None = None


@dataclass
class Report:
    """Outcome of a render: per figure, per table, failed scripts, anchor lines."""

    figures: dict[Path, FigureResult] = field(default_factory=dict)
    tables: dict[Path, str] = field(default_factory=dict)
    failures: dict[Path, str] = field(default_factory=dict)
    anchors: list[str] = field(default_factory=list)


def git(root: Path, *args: str) -> str:
    """Output of one git command in ``root``."""
    return subprocess.run(["git", *args], cwd=root, check=True, capture_output=True,
                          encoding="utf-8").stdout


def parse_figure_table(readme: Path) -> list[FigureRow]:
    """Rows of the figure tables; links are resolved relative to the README."""
    base = readme.parent.resolve()
    root = base.parents[1]
    rows = []
    for line in readme.read_text(encoding="utf-8").splitlines():
        m = _ROW.match(line)
        if m:
            figure, script = ((base / m.group(i)).resolve().relative_to(root) for i in (1, 2))
            rows.append(FigureRow(figure, script, m.group(3)))
    return rows


def assert_clean(root: Path, subdir: str = "data") -> None:
    """Raise if a tracked file under ``subdir`` differs from HEAD (documents aside)."""
    dirty = git(root, "status", "--porcelain", "--untracked-files=no", "--", subdir,
                ":(exclude)*.md")
    if dirty.strip():
        raise RuntimeError(f"tracked files under {subdir}/ are modified; commit them or "
                           f"put them back first:\n{dirty}")


def same_content(new: bytes, old: bytes, suffix: str) -> bool:
    """Equal pixels (images), equal text up to line endings (CSV), else equal bytes."""
    if suffix in (".png", ".jpg"):
        a, b = (Image.open(io.BytesIO(x)).convert("RGB") for x in (new, old))
        return a.size == b.size and a.tobytes() == b.tobytes()
    if suffix == ".csv":
        return new.replace(b"\r\n", b"\n") == old.replace(b"\r\n", b"\n")
    return new == old


def _run_script(script: Path, root: Path) -> tuple[str, str]:
    """Run one script; return (reason of failure or an empty string, stdout)."""
    env = {**os.environ, "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8"}
    try:
        p = subprocess.run([sys.executable, str(root / script)], cwd=root, env=env,
                           capture_output=True, encoding="utf-8", errors="replace",
                           stdin=subprocess.DEVNULL, timeout=SCRIPT_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return f"timeout after {SCRIPT_TIMEOUT_S} s", ""
    reason = "" if p.returncode == 0 else f"exit {p.returncode}: {p.stderr.strip()[-800:]}"
    return reason, p.stdout


def _restore(before: dict[Path, bytes]) -> None:
    """Put every file back; a file that cannot be written does not stop the others."""
    failed = []
    for path, content in before.items():
        try:
            if not path.exists() or path.read_bytes() != content:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
        except OSError as exc:
            failed.append(f"{path}: {exc}")
    if failed:
        raise RuntimeError("not put back; use git checkout on these files:\n"
                           + "\n".join(failed))


def _status(new: bytes, old: bytes | None, suffix: str) -> str:
    return "identical" if old is not None and same_content(new, old, suffix) else "differs"


def _mtimes(paths) -> dict[Path, int]:
    return {p: p.stat().st_mtime_ns for p in paths if p.exists()}


def _save(out_root: Path | None, relative: Path, content: bytes) -> None:
    if out_root is not None:
        (out_root / relative).parent.mkdir(parents=True, exist_ok=True)
        (out_root / relative).write_bytes(content)


def render(rows: list[FigureRow], root: Path, *, out_root: Path | None = None,
           table_scripts: tuple[Path, ...] = ()) -> Report:
    """Run each script once; classify what it wrote and put tracked ``data/`` back.

    A table is put back after every script, so each one compares with the committed
    numbers; figure data stay until the end (``tier0_law_figure`` draws from the CSVs
    of ``tier0_thesis_figure``). After a failed script everything is put back.

    Parameters
    ----------
    rows : list of FigureRow
        Figures to draw.
    root : Path
        Project root, a git work tree with a clean ``data/``.
    out_root : Path, optional
        Copies go to ``<out_root>/<chapter>/`` and ``<out_root>/tables/``.
    table_scripts : tuple of Path
        Scripts without figures, run last.
    """
    assert_clean(root)
    tracked = [root / f for f in git(root, "ls-files", "-z", "--", "data", ":(exclude)*.md"
                                     ).split("\0") if f]
    before = {p: p.read_bytes() for p in tracked}
    tables = [p for p in tracked if p.is_relative_to(root / TABLES_DIR)]
    not_figures = {p: c for p, c in before.items() if not p.is_relative_to(root / "data/figures")}
    if out_root is not None:
        for row in rows:
            (out_root / row.figure.parent.name / row.figure.name).unlink(missing_ok=True)
    report = Report()
    try:
        for script in dict.fromkeys([*(r.script for r in rows), *table_scripts]):
            mine = [r.figure for r in rows if r.script == script]
            mtimes = _mtimes([*(root / f for f in mine), *tables])
            print(f"[figures] {script.as_posix()}", flush=True)
            failure, stdout = _run_script(script, root)
            if failure:
                report.failures[script] = failure
            for line in stdout.splitlines():
                if line.startswith(PREFIX) and line[len(PREFIX):] not in report.anchors:
                    report.anchors.append(line[len(PREFIX):])
            for figure in mine:
                path = root / figure
                result = FigureResult(script, reason=failure, old=before.get(path))
                report.figures[figure] = result
                if not path.exists() or path.stat().st_mtime_ns == mtimes.get(path):
                    result.reason = failure or "script did not write this file"
                    continue
                result.new = path.read_bytes()
                if not failure:
                    try:
                        result.status = _status(result.new, result.old, path.suffix)
                    except OSError as exc:            # Pillow: not an image
                        result.reason = f"unreadable: {exc}"
                        continue
                    _save(out_root, Path(figure.parent.name, figure.name), result.new)
            if not failure:
                for path in tables:
                    if path.exists() and path.stat().st_mtime_ns != mtimes[path]:
                        new = path.read_bytes()
                        report.tables[path.relative_to(root)] = _status(
                            new, before[path], path.suffix)
                        _save(out_root, Path("tables", path.name), new)
            _restore(before if failure else not_figures)
    finally:
        _restore(before)
    return report


def _selection(rows: list[FigureRow]) -> tuple[list[FigureRow], tuple[Path, ...]]:
    """Figure rows and table scripts chosen by ``SCRIPTS``."""
    if SCRIPTS == "all":
        return rows, tuple(TABLE_SCRIPTS)
    known = dict.fromkeys([*(r.script.stem for r in rows), *(s.stem for s in TABLE_SCRIPTS)])
    unknown = [s for s in SCRIPTS if s not in known]
    if unknown:
        raise KeyError(f"unknown script {unknown}; expected some of {list(known)}")
    return ([r for r in rows if r.script.stem in SCRIPTS],
            tuple(s for s in TABLE_SCRIPTS if s.stem in SCRIPTS))


def main() -> int:
    """Draw the selected figures and print the report."""
    rows, table_scripts = _selection(parse_figure_table(README))
    report = render(rows, PROJECT_ROOT, out_root=OUTPUT_ROOT, table_scripts=table_scripts)
    for figure, result in report.figures.items():
        status = result.status.replace("_", " ")
        print(f"{status:13s} {figure.as_posix()}")
        if result.status == "not_rendered":
            print(f"              {result.reason.splitlines()[-1]}")
    for table, status in report.tables.items():
        print(f"{status:13s} {table.as_posix()}")
    for script, reason in report.failures.items():
        print(f"failed        {script.as_posix()}: {reason.splitlines()[-1]}")
    for text in report.anchors:
        print(PREFIX + text)
    not_rendered = [f for f, r in report.figures.items() if r.status == "not_rendered"]
    print(f"{len(report.figures)} figures, {len(not_rendered)} not rendered, "
          f"{len(report.failures)} scripts failed, {len(report.anchors)} numbers differ; "
          f"copies in {OUTPUT_ROOT}")
    return 1 if not_rendered or report.failures or report.anchors else 0


if __name__ == "__main__":
    raise SystemExit(main())
