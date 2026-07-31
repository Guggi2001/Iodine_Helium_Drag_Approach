"""Focused tests for the free-form linear ring scorer's pure helpers.

Band arithmetic and verdict-stamp helpers only — no run directories, no MD
(the scorer's oracles gate the real read at invocation time). Tolerances
exact: these are declarative band checks, not numerical physics.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

from scripts.post_processing.tier2atlas_linring_table import (
    DEEPKE_HOT,
    DEEPKE_NEUTRAL,
    GATE_N1,
    GATE_NBAR,
    KKE_MIN_DELTA_EV,
    SUCCESS_KE1_EV,
    deepke_class,
    hard_gate,
    ke_width_columns,
)


class TestFrozenBands:
    def test_bands_are_the_preregistered_ones(self):
        # Plan §6.1/§6.2 + the g3ring MD-side gate precedent.
        assert GATE_N1 == (0.19, 0.30)
        assert GATE_NBAR == (3.77, 4.37)
        assert KKE_MIN_DELTA_EV == 0.05
        assert SUCCESS_KE1_EV == 0.75
        assert DEEPKE_NEUTRAL == (0.6, 1.6)
        assert DEEPKE_HOT == 2.0


class TestHardGate:
    @pytest.mark.parametrize("n1,nbar,want", [
        (0.19, 3.77, 1), (0.30, 4.37, 1), (0.24, 4.0, 1),   # in incl. edges
        (0.1899, 4.0, 0), (0.3001, 4.0, 0),                 # n1 out
        (0.24, 3.7699, 0), (0.24, 4.3701, 0),               # nbar out
        (0.18, 3.5, 0),                                     # both out
    ])
    def test_boundaries_inclusive(self, n1, nbar, want):
        assert hard_gate(n1, nbar) == want


class TestDeepkeClass:
    @pytest.mark.parametrize("value,want", [
        (0.6, "neutral[0.6,1.6]"), (1.6, "neutral[0.6,1.6]"),
        (1.0, "neutral[0.6,1.6]"),
        (1.7, "strip(1.6,2.0)"), (1.99, "strip(1.6,2.0)"),
        (2.0, "hot(>=2.0)"), (3.4, "hot(>=2.0)"),
        (0.5, "cold(<0.6)"),
        (float("nan"), "nan"),
    ])
    def test_band_edges(self, value, want):
        assert deepke_class(value) == want


class TestKeWidthColumns:
    def _read(self, n, ke):
        return types.SimpleNamespace(
            n_scored=np.asarray(n), ke_scored_eV=np.asarray(ke, dtype=float))

    def test_sd_and_above115_on_the_n1_bin_only(self):
        read = self._read(
            [1, 1, 1, 1, 1, 1, 2, 2],
            [0.5, 0.9, 1.0, 1.1, 1.2, 1.3, 5.0, 5.0],
        )
        cols = ke_width_columns(read)
        sel = np.array([0.5, 0.9, 1.0, 1.1, 1.2, 1.3])
        assert cols["KE1_SD"] == pytest.approx(sel.std(ddof=1), rel=1e-12)
        assert cols["above115_n1"] == pytest.approx(2.0 / 6.0, rel=1e-12)

    def test_under_five_fragments_is_nan(self):
        # Matches the lowke_columns small-bin convention.
        cols = ke_width_columns(self._read([1, 1, 1, 1], [1.0] * 4))
        assert np.isnan(cols["KE1_SD"])
        assert np.isnan(cols["above115_n1"])
