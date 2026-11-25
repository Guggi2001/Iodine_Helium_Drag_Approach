"""Map thesis figures to the code they need and check the allowlist (temporary; removed in P6).

Usage, from the project root::

    python -m scripts.handover.allowlist draft   # once; never overwrites
    python -m scripts.handover.allowlist check
"""

from __future__ import annotations

import ast
import csv
import dataclasses
import subprocess
import sys
import warnings
from collections.abc import Callable
from functools import cache
from pathlib import Path

from i2_helium_md.config import SimConfig
from i2_helium_md.presets import single_pulse_droplet_distribution
from i2_helium_md.simulation.run_directory import RunDirectory
from scripts.handover.baseline import parse_figure_table

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CSV_PATH = PROJECT_ROOT / "docs" / "handover" / "ALLOWLIST.csv"
COLUMNS = ("item", "kind", "needed_by", "inputs", "decision", "note")
KINDS = ("figure", "number", "script", "module", "option", "data", "doc", "test")
DECISIONS = ("keep", "rewire", "move", "delete")
BARE_IMPORT_ROOTS = ("scripts", "data/reference/scripts")
CODE_ROOTS = ("i2_helium_md", "scripts", "data/reference/scripts", "tests")
#: label -> run directory of every simulation the thesis figures read
THESIS_RUNS = {
    "B2": "B2_N1000_s1",
    "B3": "B3_N1000_s1",
    "h405": "h405_N1000_s1",
    "T0_18A": "18A_drag_srho3p5A18A_linear_quadratic_N50_B2m21",
    "T0_18A_stuck": "18A_drag_srho3p5A18A_linear_quadratic_N50_B2m21stuck",
    "T0_9A": "9A_drag_srho3p5A18A_linear_quadratic_N50_B2m21",
    "T0_9A_h405": "9A_drag_shared_pure_cubic_N50_g11fig_h405pp35m21",
    "HS_9A": "9A_hedft_comparison",
    "HS_18A": "18A_hedft_comparison",
}
#: tracked folders that get one row instead of one row per file
FOLDERS = ("../Drag_function/", "../molecular_dynamics_lammps/", "old_data/", ".crossnote/",
           "i2_helium_md.egg-info/", "scripts/cross_reference/", "data/tier2_probe_figures/")


# --- imports ---------------------------------------------------------------

def _parse(source: Path) -> ast.Module:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # invalid escapes in scanned scripts
        return ast.parse(source.read_text(encoding="utf-8-sig"))


@cache
def _bare_modules(root: Path) -> dict[str, Path]:
    return {p.stem: p for d in BARE_IMPORT_ROOTS for p in (root / d).rglob("*.py")}


@cache
def _names(folder: Path) -> frozenset[str]:
    """Exact-case entry names; Windows' ``is_file`` would match ``X.py`` for ``x.py``."""
    return frozenset(p.name for p in folder.iterdir()) if folder.is_dir() else frozenset()


def _module_file(dotted: str, root: Path) -> Path | None:
    """File of a dotted module, or its package ``__init__.py``."""
    base = root.joinpath(*dotted.split("."))
    for p in (base.with_suffix(".py"), base / "__init__.py"):
        if p.name in _names(p.parent):
            return p
    return None


def _absolute(node: ast.ImportFrom, source: Path, root: Path) -> str:
    """Dotted module named by a ``from`` import, with relative levels resolved."""
    if not node.level:
        return node.module or ""
    package = source.relative_to(root).parts[:-node.level]
    return ".".join((*package, *((node.module,) if node.module else ())))


def _targets_of_name(module: str, name: str, root: Path) -> set[Path]:
    """Files needed for ``from module import name``."""
    sub = _module_file(f"{module}.{name}", root)
    if sub:
        return {sub}
    home = _module_file(module, root)
    if home is None:
        bare = _bare_modules(root).get(module)
        return {bare} if bare else set()
    if home.name != "__init__.py":
        return {home}
    found = {home}
    for node in ast.walk(_parse(home)):
        if isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if (alias.asname or alias.name) == name:
                    found |= _targets_of_name(_absolute(node, home, root), alias.name, root)
    return found


def direct_imports(source: Path, root: Path) -> set[Path]:
    """Repo files imported anywhere in ``source`` (function level included).

    A package ``__init__.py`` is listed but not followed: only the modules that
    define the imported names count.
    """
    if source.name == "__init__.py":
        return set()
    found: set[Path] = set()
    for node in ast.walk(_parse(source)):
        if isinstance(node, ast.ImportFrom):
            module = _absolute(node, source, root)
            for alias in node.names:
                found |= _targets_of_name(module, alias.name, root)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                hit = _module_file(alias.name, root) or _bare_modules(root).get(alias.name)
                if hit:
                    found.add(hit)
    return found - {source}


def import_closure(entry: Path, root: Path) -> set[str]:
    """Posix paths, relative to ``root``, of every repo file ``entry`` reaches."""
    seen: set[Path] = set()
    todo = [entry]
    while todo:
        source = todo.pop()
        if source not in seen:
            seen.add(source)
            todo.extend(direct_imports(source, root) - seen)
    return {p.relative_to(root).as_posix() for p in seen}


# --- draft -----------------------------------------------------------------

def repo_files(root: Path) -> list[str]:
    """Tracked and not-ignored files, relative to ``root`` (siblings as ``../``)."""
    out = subprocess.run(["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard",
                          ".."], cwd=root, check=True, capture_output=True, encoding="utf-8")
    return sorted({f for f in out.stdout.split("\0") if f and (root / f).exists()})


def _kind(item: str) -> str:
    top = item.split("/")[0]
    if top == "tests":
        return "test"
    if top == "i2_helium_md":
        return "module"
    if item.endswith((".py", ".m", ".asv")) or top == "scripts" or item in FOLDERS[:2]:
        return "script"
    return "data" if top in ("data", "old_data") else "doc"


def _option_rows(root: Path) -> list[dict]:
    cfgs = {label: RunDirectory(root / "data/runs" / name).load_cfg()
            for label, name in THESIS_RUNS.items()
            if (root / "data/runs" / name / "cfg.json").is_file()}
    # single_pulse_droplet_long: its cfg.json was overwritten on 2026-11-25 and no longer
    # describes ion.npz, so the matching preset stands in.
    cfgs["HS_long_preset"] = single_pulse_droplet_distribution(num_molecules=20000)
    rows = []
    for field in dataclasses.fields(SimConfig):
        shown = "; ".join(f"{label}={str(getattr(cfg, field.name))[:40]}"
                          for label, cfg in cfgs.items())
        rows.append(dict(zip(COLUMNS, (f"SimConfig.{field.name}", "option", "", shown, "", ""))))
    return rows


def draft_rows(root: Path) -> list[dict]:
    """Undecided rows for every file and config field; figures are ``keep``."""
    table = parse_figure_table(root / "data/figures/README.md")
    figure = {r.figure.as_posix(): r for r in table}
    needed: dict[str, set[str]] = {}
    for script in dict.fromkeys(r.script for r in table):
        names = {r.figure.name for r in table if r.script == script}
        for path in import_closure(root / script, root):
            needed.setdefault(path, set()).update(names)
    rows = []
    for f in repo_files(root):
        if f.startswith(FOLDERS):
            continue
        row = figure.get(f)
        rows.append(dict(zip(COLUMNS, (
            f, "figure" if row else _kind(f), " ".join(sorted(needed.get(f, ()))),
            f"{row.script.as_posix()}; {row.inputs}" if row else "", "keep" if row else "", ""))))
    rows += [dict(zip(COLUMNS, (f, _kind(f), "", "", "", ""))) for f in FOLDERS]
    return rows + _option_rows(root)


def write_draft(path: Path, rows: list[dict]) -> None:
    """Write the CSV; raise ``FileExistsError`` if it exists (it is edited by hand)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# --- checks ----------------------------------------------------------------

def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def decisions(path: Path) -> Callable[[str], str]:
    """Lookup file -> decision: exact row first, then the longest folder row."""
    table = {r["item"]: r["decision"] for r in read_rows(path)}
    folders = sorted((i for i in table if i.endswith("/")), key=len, reverse=True)

    def decide(item: str) -> str:
        key = item if item in table else next((f for f in folders if item.startswith(f)), "")
        return table.get(key, "")

    return decide


def check_coverage(csv_path: Path, tracked: list[str], fields: list[str]) -> list[str]:
    """Problems: uncovered, vanished, duplicate, undecided or invalid rows."""
    rows = read_rows(csv_path)
    items = [r["item"] for r in rows]
    problems = [f"duplicate: {i}" for i in sorted({i for i in items if items.count(i) > 1})]
    for r in rows:
        if not r["decision"]:
            problems.append(f"undecided: {r['item']}")
        elif r["decision"] not in DECISIONS:
            problems.append(f"bad decision {r['decision']!r}: {r['item']}")
        if r["kind"] not in KINDS:
            problems.append(f"bad kind {r['kind']!r}: {r['item']}")
    files = {r["item"] for r in rows if r["kind"] not in ("option", "number")}
    folders = tuple(i for i in files if i.endswith("/"))
    known = set(tracked)
    problems += [f"uncovered: {t}" for t in tracked if t not in files and not t.startswith(folders)]
    problems += [f"vanished: {i}" for i in sorted(files) if not i.endswith("/") and i not in known]
    problems += [f"vanished: {i}" for i in folders if not any(t.startswith(i) for t in tracked)]
    options = {r["item"] for r in rows if r["kind"] == "option"}
    wanted = {f"SimConfig.{f}" for f in fields}
    problems += [f"uncovered: {o}" for o in sorted(wanted - options)]
    return problems + [f"vanished: {o}" for o in sorted(options - wanted)]


def path_loads(source: Path, root: Path) -> set[Path]:
    """Repo scripts named in a string of a file that loads code dynamically."""
    tree = _parse(source)
    if not any("importlib" in (alias.name, getattr(node, "module", "") or "")
               or alias.name.startswith("importlib.") or alias.name in ("runpy", "subprocess")
               for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))
               for alias in node.names):
        return set()
    stems = {Path(n.value).stem for n in ast.walk(tree)
             if isinstance(n, ast.Constant) and isinstance(n.value, str) and len(n.value) < 200}
    return {p for s, p in _bare_modules(root).items() if s in stems} - {source}


def option_readers(csv_path: Path, root: Path) -> list[tuple[str, str]]:
    """(file, field) where a ``keep`` file names a config field that is not ``keep``."""
    decide = decisions(csv_path)
    fields = {r["item"].removeprefix("SimConfig.") for r in read_rows(csv_path)
              if r["kind"] == "option" and r["decision"] in ("delete", "rewire")}
    found = []
    for path in sorted(p for d in CODE_ROOTS for p in (root / d).rglob("*.py")):
        name = path.relative_to(root).as_posix()
        if decide(name) == "keep" and "handover" not in name:
            words = set()
            for n in ast.walk(_parse(path)):
                words |= {getattr(n, a) for a in ("attr", "id", "arg")
                          if isinstance(getattr(n, a, None), str)}
                if isinstance(n, ast.Constant) and isinstance(n.value, str):
                    words.add(n.value)
            found += [(name, f) for f in sorted(fields & words)]
    return found


def keep_violations(csv_path: Path, root: Path) -> list[tuple[str, str]]:
    """(importer, imported) pairs where kept, rewired or moved code needs deleted code."""
    decide = decisions(csv_path)
    found = []
    for path in sorted(p for d in CODE_ROOTS for p in (root / d).rglob("*.py")):
        importer = path.relative_to(root).as_posix()
        if decide(importer) in ("keep", "rewire", "move"):
            targets = sorted(p.relative_to(root).as_posix()
                             for p in direct_imports(path, root) | path_loads(path, root))
            found += [(importer, t) for t in targets if decide(t) == "delete"]
    return found


def main(argv: list[str]) -> int:
    if argv == ["draft"]:
        write_draft(CSV_PATH, draft_rows(PROJECT_ROOT))
        return 0
    if argv == ["check"]:
        fields = [f.name for f in dataclasses.fields(SimConfig)]
        problems = check_coverage(CSV_PATH, repo_files(PROJECT_ROOT), fields)
        for p in problems:
            print(p)
        print(f"{len(problems)} problems")
        return int(bool(problems))
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
