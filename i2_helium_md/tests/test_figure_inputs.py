"""Every run directory a kept script or the figure table names is made by a driver."""

from __future__ import annotations

import csv
import re
from itertools import product
from pathlib import Path

import pytest

from i2_helium_md.presets import LEGACY_THESIS_RUNS, THESIS_MEMBER_SEEDS, THESIS_PRESETS
from i2_helium_md.simulation.thesis_runs import TABLES_ROOT, member_dir, pool_dir
from scripts.tier0_runs import tier0_runs

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "data" / "figures" / "README.md"
# RUN_DIR of gen_legacy_runs is the output of a manual run.
EXEMPT = ("scripts/gen_legacy_runs.py", "scripts/handover/")
RETIRED = re.compile(r"(?:9A|18A|\})_hedft_comparison\b|\bsingle_pulse_droplet(?:_long)?\b")
NAMED = re.compile(r"data/runs/([^\s`\"'|)/,]+(?:,[^\s`\"'|)/,]+)*[^\s`\"'|)/]*)"
                   r"|\"runs\"\s*/\s*\"([^\"]+)\"")


def kept_scripts():
    with open(ROOT / "docs/handover/ALLOWLIST.csv", newline="", encoding="utf-8") as fh:
        items = [r["item"] for r in csv.DictReader(fh)
                 if r["kind"] == "script" and r["decision"] in ("keep", "rewire")]
    return [i for i in items if i.endswith(".py") and (ROOT / i).is_file()
            and not i.startswith(EXEMPT)]


def expand(name):
    """``{a,b}_s1…s3`` -> every name it stands for."""
    parts = re.split(r"(\{[^}]*\})", name)
    choices = [p[1:-1].split(",") if p.startswith("{") else [p] for p in parts]
    out = []
    for combo in product(*choices):
        text = "".join(combo)
        span = re.search(r"s(\d)…s(\d)", text)
        numbers = range(int(span[1]), int(span[2]) + 1) if span else [None]
        out += [text if n is None else text.replace(span[0], f"s{n}") for n in numbers]
    return [t.replace("law", law) for t in out
            for law in (THESIS_PRESETS if "law" in t else ("",))]


def made_by_a_driver():
    names = {pool_dir(law).name for law in THESIS_PRESETS}
    names |= {member_dir(law, m).name for law in THESIS_PRESETS for m in THESIS_MEMBER_SEEDS}
    return names | set(LEGACY_THESIS_RUNS) | set(tier0_runs()) | {TABLES_ROOT.name}


def named_runs(path):
    text = path.read_text(encoding="utf-8")
    return {n for m in NAMED.finditer(text) for n in expand(m[1] or m[2])}


def test_expand():
    assert expand("{B2,h405}_x") == ["B2_x", "h405_x"]
    assert expand("B2_N1000_s1…s3") == ["B2_N1000_s1", "B2_N1000_s2", "B2_N1000_s3"]
    assert "h405_N1000_s6" in expand("{law}_N1000_s1…s6")


def test_readme_names_are_found():
    assert {"B2_N6000_pooled", "B2_N1000_s6", "legacy_droplet_N20000"} <= named_runs(README)


@pytest.mark.parametrize("item", ["data/figures/README.md", *kept_scripts()])
def test_named_runs_are_made_by_a_driver(item):
    path = ROOT / item
    assert named_runs(path) - made_by_a_driver() == set()
    assert not RETIRED.findall(path.read_text(encoding="utf-8"))


def test_continuum_drag_compares_with_the_seeded_presets():
    from scripts.post_processing.legacy_continuum_drag_figure import CASES

    for preset, name, _, _ in CASES.values():
        assert preset == LEGACY_THESIS_RUNS[name]()
