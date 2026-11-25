"""Names of the thesis run directories and the missing-run message."""

from __future__ import annotations

import pytest

from i2_helium_md.presets import THESIS_MEMBER_SEEDS, THESIS_PRESETS
from i2_helium_md.simulation.thesis_runs import (
    RUNS_ROOT,
    TABLES_ROOT,
    member_dir,
    pool_dir,
    require_run,
)


def test_names_of_members_and_pools(tmp_path):
    assert member_dir("B2", "s6").name == "B2_N1000_s6"
    assert pool_dir("h405", tag="_x").name == "h405_x_N6000_pooled"
    assert member_dir("B3", "s1", root=tmp_path, tag="_x") == tmp_path / "B3_x_N1000_s1"
    assert pool_dir("B2").parent == RUNS_ROOT


def test_runs_root_is_the_project_one():
    assert RUNS_ROOT.parts[-3:] == ("i2_helium_md", "data", "runs")
    assert (TABLES_ROOT / "g13_decoupled_table.csv").is_file()


@pytest.mark.parametrize("call", [lambda: member_dir("B9", "s1"), lambda: pool_dir("B9")])
def test_unknown_law_lists_the_choices(call):
    with pytest.raises(KeyError, match="h405"):
        call()


def test_unknown_member_lists_the_choices():
    with pytest.raises(KeyError, match="s6"):
        member_dir("B2", "r6")


def test_missing_run_names_the_path_and_the_driver(tmp_path):
    with pytest.raises(FileNotFoundError, match="gone") as err:
        require_run(tmp_path / "gone", "gen_thesis_laws.py")
    assert "scripts/gen_thesis_laws.py" in str(err.value)
    assert require_run(tmp_path, "gen_thesis_laws.py") == tmp_path


def test_names_equal_the_old_generator_names():
    from scripts.gen_tier2_detfix_postpickup import thesis_pool_dir, thesis_run_dir

    for law in THESIS_PRESETS:
        assert pool_dir(law) == thesis_pool_dir(law)
        for member in THESIS_MEMBER_SEEDS:
            old = "r6" if member == "s6" else member
            assert member_dir(law, member) == thesis_run_dir(law, old)
