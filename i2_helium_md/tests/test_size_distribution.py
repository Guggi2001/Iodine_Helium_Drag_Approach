"""Tests for postprocess/size_distribution.py (Tier-2 Phase E, Slice E1).

The E1 extractor turns a finished generative ion run (or an E2 relaxation
result) into the terminal I+He_n integer-n size distribution -- the direct
counterpart of the experimental abundance reference. These tests exercise the
pure extraction against hand-counted oracles and the fail-loud validation arms
(fractional / out-of-range n, empty ensemble, wrong source type). Tiny
synthetic checkpoints only; no figures, no production-sized arrays.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from i2_helium_md.physics.shell_schedule import complex_mass_amu
from i2_helium_md.postprocess.size_distribution import (
    N_STAR,
    ShellDistribution,
    compute_terminal_shell_distribution,
)
from tests.test_checkpoint import _make_ion_checkpoint


# ---------------------------------------------------------------------------
# Duck-typed stubs. E1 reads only the terminal shell-count vector; the sim-end
# mode duck-types on ``n_shell (2N, T)`` (IonCheckpoint), the relaxed mode on
# ``terminal_n (2N,)`` (the future E2 RelaxationResult).
# ---------------------------------------------------------------------------
def _ckpt_stub(terminal_ns, *, num_steps: int = 3) -> SimpleNamespace:
    """IonCheckpoint-like stub: ``n_shell (2N, num_steps)``.

    Earlier columns are pinned to the n* sentinel (21.0) so a bug that reads
    the wrong column instead of the terminal one is caught by the oracles.
    """
    terminal = np.asarray(terminal_ns, dtype=float)
    two_n = terminal.size
    n_shell = np.full((two_n, num_steps), float(N_STAR))
    if two_n:
        n_shell[:, -1] = terminal
    return SimpleNamespace(n_shell=n_shell)


def _relax_stub(terminal_ns) -> SimpleNamespace:
    """RelaxationResult-like stub: exposes ``terminal_n (2N,)``."""
    return SimpleNamespace(terminal_n=np.asarray(terminal_ns, dtype=float))


# ---------------------------------------------------------------------------
# Hand-counted oracle
# ---------------------------------------------------------------------------
class TestHandCount:
    def test_hand_count_matches(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([0, 0, 1, 2, 2, 2]))
        assert dist.counts[0] == 2
        assert dist.counts[1] == 1
        assert dist.counts[2] == 3
        assert np.all(dist.counts[3:] == 0)
        # Terminal column is [0,0,1,2,2,2]; the 21.0 sentinel in earlier
        # columns must NOT leak in.
        assert dist.counts[N_STAR] == 0
        expected = np.array([2, 1, 3] + [0] * (N_STAR - 2), dtype=float) / 6.0
        np.testing.assert_allclose(dist.fraction, expected)
        assert dist.source == "sim_end"

    def test_fraction_sums_to_one(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([0, 1, 2, 3, 4, 5]))
        assert dist.fraction.sum() == pytest.approx(1.0)

    def test_counts_sum_equals_atom_count(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([5, 5, 6, 7]))
        assert int(dist.counts.sum()) == 4


# ---------------------------------------------------------------------------
# Integer support (0..n*), incl. the legal n = 21 endpoint
# ---------------------------------------------------------------------------
class TestIntegerSupport:
    def test_support_is_contiguous_zero_to_nstar(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([10, 11]))
        np.testing.assert_array_equal(dist.n_values, np.arange(N_STAR + 1))
        assert dist.n_values.dtype.kind == "i"
        assert dist.counts.shape == dist.n_values.shape == dist.fraction.shape

    def test_n_equals_nstar_counted_not_clipped(self):
        # n = 21 is a legal simulated terminal outcome (runs start at n* = 21).
        dist = compute_terminal_shell_distribution(_ckpt_stub([21, 21, 20]))
        assert dist.counts[N_STAR] == 2
        assert dist.counts[20] == 1

    def test_custom_n_max(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([2, 3]), n_max=5)
        np.testing.assert_array_equal(dist.n_values, np.arange(6))
        assert dist.counts[2] == 1 and dist.counts[3] == 1

    def test_monotone_falling_envelope_reproduced(self):
        terminal = [0] * 5 + [1] * 3 + [2] * 2 + [3] * 1
        dist = compute_terminal_shell_distribution(_ckpt_stub(terminal))
        env = dist.counts[:4]
        assert env[0] == 5 and env[1] == 3 and env[2] == 2 and env[3] == 1
        assert np.all(np.diff(env) <= 0)  # monotone non-increasing
        assert np.all(dist.counts[4:] == 0)


# ---------------------------------------------------------------------------
# Both input modes
# ---------------------------------------------------------------------------
class TestInputModes:
    def test_relaxed_mode_tag(self):
        dist = compute_terminal_shell_distribution(_relax_stub([0, 1, 1]))
        assert dist.source == "relaxed"
        assert dist.counts[0] == 1 and dist.counts[1] == 2

    def test_both_modes_identical_histogram(self):
        terminal = [1, 2, 2, 3, 3, 3]
        d_sim = compute_terminal_shell_distribution(_ckpt_stub(terminal))
        d_relax = compute_terminal_shell_distribution(_relax_stub(terminal))
        np.testing.assert_array_equal(d_sim.counts, d_relax.counts)
        np.testing.assert_allclose(d_sim.fraction, d_relax.fraction)
        assert d_sim.source == "sim_end" and d_relax.source == "relaxed"

    def test_real_ion_checkpoint_accepted(self):
        # A real v7 IonCheckpoint (n_shell int-valued floats in [14, 22)).
        ckpt = _make_ion_checkpoint(num_molecules=3, num_steps=6)
        dist = compute_terminal_shell_distribution(ckpt)
        assert dist.source == "sim_end"
        assert int(dist.counts.sum()) == 2 * 3
        np.testing.assert_array_equal(dist.n_values, np.arange(N_STAR + 1))
        assert np.all(dist.counts >= 0)
        assert dist.fraction.sum() == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Explicit provenance override (reloaded relaxation.npz case)
# ---------------------------------------------------------------------------
class TestSourceTagOverride:
    """A relaxation checkpoint reloaded from ``relaxation.npz`` is a bona fide
    v7 ``IonCheckpoint``: the duck-typed dispatch alone would infer
    ``"sim_end"`` although its terminal column IS the matched-time relaxed
    ``n``. The explicit ``source_tag`` override is the caller's provenance
    escape hatch (E4/F3 pair the two scores per run on this tag)."""

    def test_override_relaxed_on_checkpoint_input(self):
        dist = compute_terminal_shell_distribution(
            _ckpt_stub([1, 2, 2]), source_tag="relaxed"
        )
        assert dist.source == "relaxed"

    def test_override_sim_end_on_relaxed_input(self):
        dist = compute_terminal_shell_distribution(
            _relax_stub([1, 2, 2]), source_tag="sim_end"
        )
        assert dist.source == "sim_end"

    def test_override_does_not_change_histogram(self):
        terminal = [1, 2, 2, 3]
        d_default = compute_terminal_shell_distribution(_ckpt_stub(terminal))
        d_override = compute_terminal_shell_distribution(
            _ckpt_stub(terminal), source_tag="relaxed"
        )
        np.testing.assert_array_equal(d_default.counts, d_override.counts)

    def test_invalid_override_rejected(self):
        with pytest.raises(ValueError, match="sim_end|relaxed"):
            compute_terminal_shell_distribution(
                _ckpt_stub([1]), source_tag="matched"
            )

    def test_default_inference_unchanged(self):
        assert compute_terminal_shell_distribution(_ckpt_stub([1])).source == "sim_end"


# ---------------------------------------------------------------------------
# Fail-loud validation arms
# ---------------------------------------------------------------------------
class TestValidation:
    def test_fractional_entry_rejected(self):
        with pytest.raises(ValueError, match="integer-valued"):
            compute_terminal_shell_distribution(_relax_stub([1.0, 2.5, 3.0]))

    def test_out_of_range_high_rejected(self):
        with pytest.raises(ValueError, match="range"):
            compute_terminal_shell_distribution(_relax_stub([20, 22]))

    def test_out_of_range_negative_rejected(self):
        with pytest.raises(ValueError, match="range"):
            compute_terminal_shell_distribution(_relax_stub([-1, 5]))

    def test_non_finite_rejected(self):
        with pytest.raises(ValueError, match="finite"):
            compute_terminal_shell_distribution(_relax_stub([np.nan, 3.0]))

    def test_empty_ensemble_sim_end_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            compute_terminal_shell_distribution(SimpleNamespace(n_shell=np.empty((0, 3))))

    def test_zero_steps_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            compute_terminal_shell_distribution(SimpleNamespace(n_shell=np.empty((4, 0))))

    def test_empty_ensemble_relaxed_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            compute_terminal_shell_distribution(SimpleNamespace(terminal_n=np.empty(0)))

    def test_one_d_n_shell_rejected(self):
        with pytest.raises(ValueError, match="2-D"):
            compute_terminal_shell_distribution(SimpleNamespace(n_shell=np.array([1.0, 2.0])))

    def test_bad_source_type_rejected(self):
        with pytest.raises(TypeError, match="n_shell.*terminal_n|terminal_n.*n_shell"):
            compute_terminal_shell_distribution(SimpleNamespace())


# ---------------------------------------------------------------------------
# Rung mass labelling
# ---------------------------------------------------------------------------
class TestMassLabelling:
    def test_mass_amu_matches_complex_mass(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([0, 1, 2]))
        np.testing.assert_allclose(dist.mass_amu, complex_mass_amu(dist.n_values))
        assert dist.mass_amu[0] == pytest.approx(126.90)
        assert dist.mass_amu[1] == pytest.approx(126.90 + 4.0026)

    def test_frozen(self):
        dist = compute_terminal_shell_distribution(_ckpt_stub([0, 1]))
        with pytest.raises(Exception):
            dist.source = "relaxed"  # frozen dataclass


class TestNStarSingleSource:
    def test_n_star_is_the_shared_physics_constant(self):
        # Quality Principle 1: the Langmuir cap has one source of truth
        # (physics/constants.py); size_distribution re-exports it, it does not
        # redeclare it.
        from i2_helium_md.physics import constants

        assert N_STAR == constants.N_STAR
