"""Record which thesis figures re-render identically (temporary; removed in P6).

Runs every script of ``data/figures/README.md`` in place, compares each figure
with its content before the run, and puts every tracked file under ``data/``
back byte for byte. Requires a clean ``data/`` tree.

Usage, from the project root: ``python -m scripts.handover.baseline``
writes ``docs/handover/baseline.json``.
"""

from __future__ import annotations

import hashlib
import json
import sys
from importlib.metadata import version
from pathlib import Path

from scripts.run_thesis_figures import (  # noqa: F401  (tests read them from here)
    FigureRow,
    assert_clean,
    git as _git,
    parse_figure_table,
    render,
    same_content,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
README = PROJECT_ROOT / "data" / "figures" / "README.md"
OUT_PATH = PROJECT_ROOT / "docs" / "handover" / "baseline.json"


def _sha256(content: bytes | None) -> str | None:
    return None if content is None else hashlib.sha256(content).hexdigest()


def render_all(rows: list[FigureRow], root: Path) -> dict[str, dict]:
    """Run each script once and classify its figures; tracked ``data/`` is restored."""
    return {
        figure.as_posix(): {
            "script": r.script.as_posix(), "status": r.status, "reason": r.reason,
            "sha256_head": _sha256(r.old), "sha256_new": _sha256(r.new)}
        for figure, r in render(rows, root).figures.items()}


def _run_dir_state(root: Path) -> dict[str, tuple[int, int]]:
    """Size and mtime of every run-directory file (git-ignored, so not restored)."""
    runs = root / "data" / "runs"
    tracked = set(_git(root, "ls-files", "-z", "--", "data/runs").split("\0"))
    return {rel: (p.stat().st_size, p.stat().st_mtime_ns)
            for p in runs.rglob("*") if p.is_file()
            and (rel := p.relative_to(root).as_posix()) not in tracked}


def main() -> int:
    state = _run_dir_state(PROJECT_ROOT)
    figures = render_all(parse_figure_table(README), PROJECT_ROOT)
    after = _run_dir_state(PROJECT_ROOT)
    record = {
        "head": _git(PROJECT_ROOT, "rev-parse", "HEAD").strip(),
        "versions": {"python": sys.version.split()[0],
                     **{n: version(n) for n in ("numpy", "scipy", "matplotlib", "pillow")}},
        "figures": figures,
        "run_dir_changes": sorted(k for k in state.keys() | after.keys()
                                  if state.get(k) != after.get(k)),
    }
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    for name, rec in figures.items():
        if rec["status"] != "identical":
            print(f"  {rec['status']}: {name}  {rec['reason'][-200:]}")
    print("  run-directory files changed:", record["run_dir_changes"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
