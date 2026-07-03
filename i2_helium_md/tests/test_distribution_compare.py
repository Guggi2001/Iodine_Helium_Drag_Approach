"""Tests for postprocess/distribution_compare.py (Tier-2 Phase E, Slice E4).

Integer-support 1-D Wasserstein W1 = sum_n |F_sim(n) - F_ref(n)| (the CDF-gap
sum) between the E1 simulated size distribution and the E3 experimental
abundance reference, plus the ``validation_histogram_metric`` config-dispatch
(the field's first physics-live reader; chi2/ks refuse). Hand-computed micro
oracles + a scipy cross-check (import-guarded, tests only).
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference
from i2_helium_md.postprocess.distribution_compare import (
    compare_size_distributions,
    wasserstein_integer_support,
)
from i2_helium_md.postprocess.size_distribution import (
    compute_terminal_shell_distribution,
)

_REAL_CSV = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "reference"
    / "integrated_i_he_abundance.csv"
)


def _sim(n_values, fraction):
    """Duck-typed ShellDistribution: only n_values + fraction are read."""
    return SimpleNamespace(
        n_values=np.asarray(n_values, dtype=int),
        fraction=np.asarray(fraction, dtype=float),
    )


def _ref(n, ion_fraction):
    """Duck-typed HeAbundanceReference: only n + ion_fraction are read."""
    return SimpleNamespace(
        n=np.asarray(n, dtype=int),
        ion_fraction=np.asarray(ion_fraction, dtype=float),
    )


# ---------------------------------------------------------------------------
# Hand-computed W1 oracles
# ---------------------------------------------------------------------------
class TestWassersteinOracles:
    def test_identical_is_zero(self):
        sim = _sim([0, 1, 2], [0.2, 0.5, 0.3])
        ref = _ref([0, 1, 2], [0.2, 0.5, 0.3])
        assert wasserstein_integer_support(sim, ref) == pytest.approx(0.0)

    def test_shift_by_one(self):
        # Unit mass moved one rung -> W1 = 1 (moving distance 1).
        sim = _sim([0, 1, 2, 3], [0, 0, 1, 0])
        ref = _ref([0, 1, 2, 3], [0, 0, 0, 1])
        assert wasserstein_integer_support(sim, ref) == pytest.approx(1.0)

    def test_shift_by_two(self):
        sim = _sim([0, 1, 2], [1, 0, 0])
        ref = _ref([0, 1, 2], [0, 0, 1])
        assert wasserstein_integer_support(sim, ref) == pytest.approx(2.0)

    def test_micro_hand_case(self):
        # F_sim=[.5,1,1], F_ref=[0,.5,1]; |gap|=[.5,.5,0] -> W1 = 1.0.
        sim = _sim([0, 1, 2], [0.5, 0.5, 0.0])
        ref = _ref([0, 1, 2], [0.0, 0.5, 0.5])
        assert wasserstein_integer_support(sim, ref) == pytest.approx(1.0)

    def test_n21_meets_zero_filled_reference(self):
        # sim covers 0..21 (legal n=21); ref only 0..20. Union support handles it.
        sim_frac = np.zeros(22)
        sim_frac[21] = 1.0
        ref_frac = np.zeros(21)
        ref_frac[20] = 1.0
        sim = _sim(np.arange(22), sim_frac)
        ref = _ref(np.arange(21), ref_frac)
        assert wasserstein_integer_support(sim, ref) == pytest.approx(1.0)

    def test_symmetric(self):
        sim = _sim([0, 1, 2], [0.7, 0.2, 0.1])
        ref = _ref([0, 1, 2], [0.1, 0.2, 0.7])
        a = wasserstein_integer_support(sim, ref)
        b = wasserstein_integer_support(_sim([0, 1, 2], [0.1, 0.2, 0.7]),
                                        _ref([0, 1, 2], [0.7, 0.2, 0.1]))
        assert a == pytest.approx(b)


# ---------------------------------------------------------------------------
# Support guards
# ---------------------------------------------------------------------------
class TestSupportGuards:
    def test_disjoint_support_raises(self):
        sim = _sim([0, 1, 2], [0.3, 0.3, 0.4])
        ref = _ref([10, 11, 12], [0.3, 0.3, 0.4])
        with pytest.raises(ValueError, match="disjoint"):
            wasserstein_integer_support(sim, ref)

    def test_empty_sim_raises(self):
        with pytest.raises(ValueError, match="empty"):
            wasserstein_integer_support(_sim([], []), _ref([0, 1], [0.5, 0.5]))

    def test_empty_ref_raises(self):
        with pytest.raises(ValueError, match="empty"):
            wasserstein_integer_support(_sim([0, 1], [0.5, 0.5]), _ref([], []))


# ---------------------------------------------------------------------------
# Input validation (garbage in must raise, not score)
# ---------------------------------------------------------------------------
class TestInputValidation:
    def test_counts_instead_of_fraction_rejected(self):
        # ShellDistribution exposes both .counts and .fraction; passing counts
        # by mistake inflates every CDF gap into a plausible-looking but
        # meaningless W1. The sum-to-1 guard catches the misuse.
        sim = _sim([0, 1, 2], [2.0, 1.0, 3.0])
        ref = _ref([0, 1, 2], [0.2, 0.5, 0.3])
        with pytest.raises(ValueError, match="sum"):
            wasserstein_integer_support(sim, ref)

    def test_nan_sim_fraction_rejected(self):
        with pytest.raises(ValueError, match="finite"):
            wasserstein_integer_support(
                _sim([0, 1], [np.nan, 1.0]), _ref([0, 1], [0.5, 0.5])
            )

    def test_nan_ref_fraction_rejected(self):
        with pytest.raises(ValueError, match="finite"):
            wasserstein_integer_support(
                _sim([0, 1], [0.5, 0.5]), _ref([0, 1], [np.nan, 1.0])
            )

    def test_negative_ref_fraction_rejected(self):
        # Sums to 1, so only a sign guard catches it (the E3 negative-percent
        # pathway before that loader gained its own guard).
        with pytest.raises(ValueError, match="negative"):
            wasserstein_integer_support(
                _sim([0, 1], [0.5, 0.5]), _ref([0, 1], [1.5, -0.5])
            )

    def test_duplicate_support_rejected(self):
        # numpy fancy-index scatter is last-write-wins: a duplicate n would
        # silently drop probability mass before the CDF sum.
        with pytest.raises(ValueError, match="duplicate"):
            wasserstein_integer_support(
                _sim([0, 1, 1], [0.5, 0.25, 0.25]), _ref([0, 1], [0.5, 0.5])
            )


# ---------------------------------------------------------------------------
# Real E1 + E3 integration (finite score)
# ---------------------------------------------------------------------------
class TestRealIntegration:
    def test_real_ref_vs_e1_sim_finite(self):
        ref = load_he_abundance_reference(_REAL_CSV)
        # Tiny synthetic E1 distribution: terminal n in {18,19,20,21}.
        ckpt = SimpleNamespace(
            n_shell=np.array([[18.0], [19.0], [20.0], [21.0], [19.0], [20.0]])
        )
        sim = compute_terminal_shell_distribution(ckpt)
        w1 = wasserstein_integer_support(sim, ref)
        assert np.isfinite(w1)
        assert w1 >= 0.0


# ---------------------------------------------------------------------------
# Config-dispatch (validation_histogram_metric activation)
# ---------------------------------------------------------------------------
class TestDispatch:
    def test_wasserstein_matches_direct(self):
        sim = _sim([0, 1, 2], [0.5, 0.5, 0.0])
        ref = _ref([0, 1, 2], [0.0, 0.5, 0.5])
        assert compare_size_distributions(sim, ref, metric="wasserstein") == pytest.approx(
            wasserstein_integer_support(sim, ref)
        )

    def test_chi2_refuses(self):
        sim = _sim([0, 1], [0.5, 0.5])
        ref = _ref([0, 1], [0.5, 0.5])
        with pytest.raises(NotImplementedError, match="chi2"):
            compare_size_distributions(sim, ref, metric="chi2")

    def test_ks_refuses(self):
        sim = _sim([0, 1], [0.5, 0.5])
        ref = _ref([0, 1], [0.5, 0.5])
        with pytest.raises(NotImplementedError, match="ks"):
            compare_size_distributions(sim, ref, metric="ks")

    def test_unknown_metric_raises(self):
        sim = _sim([0, 1], [0.5, 0.5])
        ref = _ref([0, 1], [0.5, 0.5])
        with pytest.raises(ValueError, match="metric"):
            compare_size_distributions(sim, ref, metric="frobenius")

    def test_config_field_flows_through(self):
        # The config field is the dispatcher's intended source: exercise the
        # linkage that retires the rule-2 carry.
        cfg = single_pulse_N2000()
        assert cfg.validation_histogram_metric == "wasserstein"
        sim = _sim([0, 1, 2], [0.5, 0.5, 0.0])
        ref = _ref([0, 1, 2], [0.0, 0.5, 0.5])
        result = compare_size_distributions(
            sim, ref, metric=cfg.validation_histogram_metric
        )
        assert result == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# scipy cross-check (import-guarded, tests only)
# ---------------------------------------------------------------------------
class TestScipyCrossCheck:
    def test_matches_scipy_wasserstein(self):
        stats = pytest.importorskip("scipy.stats")
        cases = [
            ([0, 1, 2, 3], [0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]),
            ([0, 1, 2, 3, 4], [0.2, 0.2, 0.2, 0.2, 0.2], [0.0, 0.1, 0.5, 0.3, 0.1]),
        ]
        for support, fs, fr in cases:
            sim = _sim(support, fs)
            ref = _ref(support, fr)
            mine = wasserstein_integer_support(sim, ref)
            theirs = stats.wasserstein_distance(support, support, fs, fr)
            assert mine == pytest.approx(theirs)
