"""Record which thesis figures re-render identically (temporary; removed in P6).

Runs every script of ``data/figures/README.md`` in place, compares each figure
with its content before the run, and puts every tracked file under ``data/``
back byte for byte. Requires a clean ``data/`` tree.

Usage, from the project root: ``python -m scripts.handover.baseline``
writes ``docs/handover/baseline.json``.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parents[2]
README = PROJECT_ROOT / "data" / "figures" / "README.md"
OUT_PATH = PROJECT_ROOT / "docs" / "handover" / "baseline.json"
SCRIPT_TIMEOUT_S = 3600

_ROW = re.compile(r"^\|\s*\[`[^`]+`\]\(([^)]+)\)\s*\|\s*\[`[^`]+`\]\(([^)]+)\)\s*\|\s*(.*?)\s*\|$")


@dataclass(frozen=True)
class FigureRow:
    """One thesis figure: paths relative to the project root."""

    figure: Path
    script: Path
    inputs: str


def _git(root: Path, *args: str) -> str:
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
    """Raise if a tracked file under ``subdir`` differs from HEAD."""
    dirty = _git(root, "status", "--porcelain", "--untracked-files=no", "--", subdir)
    if dirty.strip():
        raise RuntimeError(f"tracked files under {subdir}/ are modified:\n{dirty}")


def same_content(new: bytes, old: bytes, suffix: str) -> bool:
    """Equal pixels (images), equal text up to line endings (CSV), else equal bytes."""
    if suffix in (".png", ".jpg"):
        a, b = (Image.open(io.BytesIO(x)).convert("RGB") for x in (new, old))
        return a.size == b.size and a.tobytes() == b.tobytes()
    if suffix == ".csv":
        return new.replace(b"\r\n", b"\n") == old.replace(b"\r\n", b"\n")
    return new == old


def _run_script(script: Path, root: Path) -> str:
    """Run one script; return an empty string on success, else the reason."""
    env = {**os.environ, "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8"}
    try:
        p = subprocess.run([sys.executable, str(root / script)], cwd=root, env=env,
                           capture_output=True, encoding="utf-8", errors="replace",
                           stdin=subprocess.DEVNULL, timeout=SCRIPT_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return f"timeout after {SCRIPT_TIMEOUT_S} s"
    return "" if p.returncode == 0 else f"exit {p.returncode}: {p.stderr.strip()[-800:]}"


def render_all(rows: list[FigureRow], root: Path) -> dict[str, dict]:
    """Run each script once and classify its figures; tracked ``data/`` is restored."""
    assert_clean(root)
    tracked = [root / f for f in _git(root, "ls-files", "-z", "--", "data").split("\0") if f]
    before = {p: p.read_bytes() for p in tracked}
    out: dict[str, dict] = {}
    try:
        for script in dict.fromkeys(r.script for r in rows):
            mine = [r for r in rows if r.script == script]
            mtimes = {r.figure: (root / r.figure).stat().st_mtime_ns
                      for r in mine if (root / r.figure).exists()}
            print(f"[baseline] {script.as_posix()}", flush=True)
            failure = _run_script(script, root)
            for r in mine:
                path = root / r.figure
                old = before.get(path)
                rec = {"script": script.as_posix(), "status": "not_rendered", "reason": failure,
                       "sha256_head": hashlib.sha256(old).hexdigest() if old else None,
                       "sha256_new": None}
                out[r.figure.as_posix()] = rec
                if not path.exists() or path.stat().st_mtime_ns == mtimes.get(r.figure):
                    rec["reason"] = failure or "script did not write this file"
                    continue
                new = path.read_bytes()
                rec["sha256_new"] = hashlib.sha256(new).hexdigest()
                if not failure:
                    same = old is not None and same_content(new, old, path.suffix)
                    rec["status"] = "identical" if same else "differs"
    finally:
        for path, content in before.items():
            if not path.exists() or path.read_bytes() != content:
                path.write_bytes(content)
    return out


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
