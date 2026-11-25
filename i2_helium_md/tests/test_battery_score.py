"""Battery score: reproduces the committed decoupled rows of the thesis table."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md.postprocess import battery_score
from i2_helium_md.postprocess.battery_score import (
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    PROVISIONAL_NORM,
    joint_score,
    law_columns,
    lowke_columns,
    score_pool,
    tail_columns,
    write_rows_csv,
)
from i2_helium_md.presets import drag_production_B2

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

RUNS = Path(__file__).resolve().parents[1] / "data" / "runs"
TABLE = RUNS / "h2b_forward_model" / "g13_decoupled_table.csv"


def decoupled_row(label):
    with TABLE.open(newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["label"] == label]
    assert len(rows) == 1 and rows[0]["self_bound_gate_release"] == "decoupled"
    return rows[0]


def is_number(text):
    try:
        float(text)
    except ValueError:
        return False
    return True


@pytest.mark.parametrize("law", ["B2", "B3", "h405"])
def test_score_reproduces_the_committed_decoupled_row(law):
    pool = RUNS / f"{law}_N6000_pooled"
    if not (pool / "detection.npz").is_file():
        pytest.skip("thesis pool not on disk")
    want = decoupled_row(law)
    got = score_pool(pool, law)
    numeric = [k for k, v in want.items() if is_number(v)]
    assert {"n1_solv", "w1_solv", "KE1_mean", "chi2_med", "S"} <= set(numeric)
    for k in numeric:
        assert round(float(got[k]), 4) == round(float(want[k]), 4), k
    assert got["label"] == law
    assert got["container"] == want["container"]
    assert got["self_bound_gate_release"] == "decoupled"


def _read():
    """Synthetic scored read: n = 1 x 10, n = 2 x 5, n = 12 x 4, n = 15 x 1, bare x 2."""
    n = np.repeat([0, 1, 2, 12, 15], [2, 10, 5, 4, 1])
    return SimpleNamespace(n_scored=n, ke_scored_eV=np.where(n == 12, 0.1, 0.05 * n))


def test_tail_columns_of_a_synthetic_read():
    got = tail_columns(_read(), SimpleNamespace(n=np.arange(1, 16), mean_KE_eV=np.full(15, 0.2)))
    assert got == {"mass_n15plus_solv": pytest.approx(1 / 20), "n_max": 15,
                   "KE10_14_over_ref": pytest.approx(0.5)}


def test_tail_columns_without_scored_ions_is_an_error():
    read = SimpleNamespace(n_scored=np.array([]), ke_scored_eV=np.array([]))
    with pytest.raises(ValueError, match="no scored ions"):
        tail_columns(read, SimpleNamespace(n=np.arange(1, 16), mean_KE_eV=np.full(15, 0.2)))


def test_tail_columns_of_bare_ions_only():
    read = SimpleNamespace(n_scored=np.zeros(3), ke_scored_eV=np.ones(3))
    got = tail_columns(read, SimpleNamespace(n=np.arange(1, 16), mean_KE_eV=np.full(15, 0.2)))
    assert math.isnan(got["mass_n15plus_solv"]) and got["n_max"] == 0


def test_lowke_columns_of_a_synthetic_read():
    got = lowke_columns(_read())
    assert got["KE1_mean"] == got["KE1_med"] == pytest.approx(0.05)
    assert got["KE2_mean"] == pytest.approx(0.1)
    assert got["KE1_mode"] == pytest.approx(0.5 * 3.5 / 40)


def test_law_columns_of_b2():
    got = law_columns(drag_production_B2())
    assert (got["E0_eV"], got["tau_ps"], got["gate_A"]) == (pytest.approx(0.3297), 6.4, 3.5)
    assert got["drag_coefficients"]["c"] == 1.6 and got["sampler"] == "post_pickup"


@pytest.mark.parametrize("n1, nbar, gate", [(0.20, 4.0, 1), (0.18, 4.0, 0), (0.20, 4.4, 0)])
def test_score_row_gate_and_joint_score(monkeypatch, n1, nbar, gate):
    monkeypatch.setattr(battery_score, "observable_columns", lambda *a: {
        "n1_solv": n1, "nbar_det": nbar, "w1_solv": 0.571, "midHot": 1.15, "deepKE": 0.1})
    row = battery_score.score_row("x", "pooled", "p", "-", _read(), None, None)
    assert row["md_gate"] == gate
    assert row["S"] == pytest.approx(2.0)      # deepKE is not a licensed axis


def test_references_are_on_disk():
    assert ABUNDANCE_REFERENCE_CSV.is_file() and IHE_KED_REFERENCE_CSV.is_file()


def test_joint_score_uses_only_the_named_axes():
    two = joint_score(0.571, 1.15, 0.5, PROVISIONAL_NORM, axes=("w1", "midhot"))
    assert two == pytest.approx(2.0)
    assert joint_score(0.571, 1.15, 0.5, PROVISIONAL_NORM) > two
    assert math.isnan(joint_score(0.571, 0.0, 1.0, PROVISIONAL_NORM))


def test_write_rows_csv_round_trip(tmp_path):
    path = write_rows_csv(tmp_path / "sub" / "t.csv", [{"label": "a", "x": 1.5}])
    with open(path, newline="", encoding="utf-8") as fh:
        assert list(csv.DictReader(fh)) == [{"label": "a", "x": "1.5"}]
