"""G4.13 generator + scorer rules (no MD, no run dirs).

- cell labels / the 16 + 1 cell set match the pre-registered grid;
- the pre-registered cfg-diff key sets;
- the frozen Stage-B selection rule and decision rule on synthetic rows.
"""

import pytest

from scripts.gen_tier2_g413_ridge_battery import (
    CELLS, GRID, H405_GATE_LABEL, cell_label, expected_diff_keys,
)
from scripts.post_processing.tier2_g413_table import decide, select_stage_b


def test_cell_set_matches_preregistered_grid():
    assert len(GRID) == 16 and len(CELLS) == 17
    assert cell_label(5.2, 0.37) == "t52e370"
    assert cell_label(5.6, 0.335) == "t56e335"
    assert not CELLS[H405_GATE_LABEL].ridge_law


def test_expected_diff_keys():
    ridge = expected_diff_keys("t48e390")
    assert {"drag_form", "drag_coefficients", "binding_energy_I_ion_eV",
            "internal_energy_cooling_tau_ps",
            "internal_energy_partition_fraction"} <= ridge
    assert expected_diff_keys(H405_GATE_LABEL) == {
        "drag_spatial_gate", "drag_gate_steepness",
        "allow_inconsistent_gate_pairing"}


def _r(cell, gate, midhot, ke1, w1=0.7):
    return {"cell": cell, "md_gate": gate, "midHot": midhot,
            "KE1_mean": ke1, "w1_solv": w1}


def test_select_in_band_max_ke1():
    rows = [_r("a", 1, 1.0, 0.55), _r("b", 1, 1.1, 0.60), _r("c", 0, 1.0, 0.9),
            _r("d", 1, 1.3, 0.70)]
    assert select_stage_b(rows) == "b"


def test_select_nearest_band_then_w1():
    rows = [_r("a", 1, 1.30, 0.6, w1=0.8), _r("b", 1, 1.20, 0.5, w1=0.9),
            _r("c", 1, 0.80, 0.5, w1=0.7)]
    # b and c are both 0.05 from the band; c has the lower W1.
    assert select_stage_b(rows) == "c"


def test_select_none_gated():
    assert select_stage_b([_r("a", 0, 1.0, 0.6)]) is None


@pytest.mark.parametrize("gate,ke1,midhot,out", [
    (1, 0.641, 1.0, "matches_or_beats_h405"),
    (1, 0.639, 1.0, "dominated"),
    (1, 0.70, 1.2, "dominated"),
    (0, 0.70, 1.0, "dominated"),
])
def test_decide(gate, ke1, midhot, out):
    assert decide({"md_gate": gate, "KE1_mean": ke1, "midHot": midhot}) == out
