"""Compare rerun values with committed ones and report, without stopping.

On another machine or numpy version a rerun may differ in rounding, so a figure
script prints the differences and still draws. ``tests/test_thesis_anchors.py``
holds the exact checks; ``run_thesis_figures.py`` collects the printed lines.
"""

from __future__ import annotations

import csv
from pathlib import Path

#: Start of every reported line.
PREFIX = "differs from committed: "


def committed_row(csv_path: Path, label: str, container: str | None = None) -> dict[str, str]:
    """The one row ``label`` of a committed table.

    Raises
    ------
    ValueError
        Not exactly one such row, or it belongs to another ``container``.
    """
    with open(csv_path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["label"] == label]
    if len(rows) != 1:
        raise ValueError(f"{csv_path.name}: expected one row {label!r}, found {len(rows)}")
    if container is not None and rows[0].get("container") != container:
        raise ValueError(
            f"{label!r} row belongs to {rows[0].get('container')!r}, not {container!r}")
    return rows[0]


def row_differences(values: dict[str, float], csv_path: Path, label: str,
                    container: str | None = None) -> list[str]:
    """Columns whose value differs from the committed row at 4 dp, one text each."""
    row = committed_row(csv_path, label, container)
    return [f"{csv_path.name} {label} {k}: {float(v):.4f} vs {float(row[k]):.4f}"
            for k, v in values.items() if round(float(v), 4) != round(float(row[k]), 4)]


def report(differences: list[str], ok: str) -> None:
    """Print each difference on a ``PREFIX`` line, or ``ok`` if there is none."""
    for text in differences:
        print(PREFIX + text, flush=True)
    if not differences:
        print(ok, flush=True)
