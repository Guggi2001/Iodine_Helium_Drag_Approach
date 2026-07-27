"""Tests for the D2b §4.3 grid re-weighting machinery.

Covers the pure functions of
``scripts/post_processing/tier2atlas_geometry_reweight.py``:

- ``interpolation_weights`` — the two pre-registered conventions (nearest /
  linear) including the clamped-mass accounting (no silent caps);
- ``mixture_read`` — the weighted-mixture duck-type. The load-bearing
  property: a mixture whose cell weights are proportional to the cells' ion
  counts must equal the read of the concatenated ensemble EXACTLY (float
  tolerance), for the fate fractions, the histogram, and the per-bin KE
  means — that is what makes the committed scorers (histogram score, midHot,
  deepKE) valid on a mixture. ``chi2_med`` is deliberately unsupported
  (weighted SE has no committed convention), enforced by the ``min_count``
  guard test.

All inputs are tiny synthetic ``DetectionResult`` objects (reusing the
``test_tier2_confirmation`` builders); no run dirs, no figures.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from i2_helium_md.postprocess.tier2_confirmation import (
    deep_bin_ke_ratio,
    midhot_ratio,
    read_confirmation_detection,
    score_histogram_vs_reference,
)
from scripts.post_processing.tier2atlas_geometry_reweight import (
    interpolation_weights,
    mixture_read,
)
from tests.test_tier2_confirmation import (
    make_abundance_reference,
    make_detection,
    make_ked_reference,
)


# ---------------------------------------------------------------------------
# interpolation_weights
# ---------------------------------------------------------------------------


class TestInterpolationWeights:
    SUPPORT = np.array([26.6, 34.0, 49.4, 68.3])

    def test_nearest_assigns_by_midpoint(self):
        # midpoints: 30.3, 41.7, 58.85
        radii = np.array([20.0, 30.0, 31.0, 50.0, 70.0])
        w, below, above = interpolation_weights(radii, self.SUPPORT, "nearest")
        np.testing.assert_allclose(w, [0.4, 0.2, 0.2, 0.2])
        assert below == pytest.approx(0.2)   # only the 20.0 sample
        assert above == pytest.approx(0.2)   # only the 70.0 sample

    def test_linear_splits_between_bracketing_cells(self):
        # one sample exactly halfway between the first two support points
        radii = np.array([(26.6 + 34.0) / 2.0])
        w, below, above = interpolation_weights(radii, self.SUPPORT, "linear")
        np.testing.assert_allclose(w, [0.5, 0.5, 0.0, 0.0])
        assert below == 0.0 and above == 0.0

    def test_linear_on_support_point_is_pure(self):
        w, _, _ = interpolation_weights(
            np.array([49.4]), self.SUPPORT, "linear"
        )
        np.testing.assert_allclose(w, [0.0, 0.0, 1.0, 0.0], atol=1e-12)

    def test_clamping_sends_outside_mass_to_end_cells(self):
        radii = np.array([5.0, 10.0, 100.0])
        for convention in ("nearest", "linear"):
            w, below, above = interpolation_weights(
                radii, self.SUPPORT, convention
            )
            np.testing.assert_allclose(w, [2 / 3, 0.0, 0.0, 1 / 3])
            assert below == pytest.approx(2 / 3)
            assert above == pytest.approx(1 / 3)

    def test_weights_sum_to_one(self):
        rng = np.random.default_rng(1)
        radii = rng.uniform(15.0, 80.0, size=500)
        for convention in ("nearest", "linear"):
            w, _, _ = interpolation_weights(radii, self.SUPPORT, convention)
            assert w.sum() == pytest.approx(1.0)

    def test_rejects_unknown_convention_and_bad_support(self):
        with pytest.raises(ValueError, match="unknown convention"):
            interpolation_weights(np.array([30.0]), self.SUPPORT, "cubic")
        with pytest.raises(ValueError, match="ascending"):
            interpolation_weights(
                np.array([30.0]), np.array([34.0, 26.6]), "nearest"
            )
        with pytest.raises(ValueError, match="empty"):
            interpolation_weights(np.array([]), self.SUPPORT, "nearest")


# ---------------------------------------------------------------------------
# mixture_read
# ---------------------------------------------------------------------------


def _cell_a():
    """4 ions: two detected (n=3, n=1), one suppressed, one retained-bound."""
    return make_detection(
        n_detected=[3.0, 7.0, 1.0, 12.0],
        state_reason=["time_exhausted", "suppressed", "time_exhausted",
                      "droplet_retained"],
        ke_eV=[0.5, 0.2, 1.0, 0.0],
    )


def _cell_b():
    """3 ions, all detected (n=3, n=2, n=12 — the deep band is occupied)."""
    return make_detection(
        n_detected=[3.0, 2.0, 12.0],
        state_reason=["time_exhausted", "time_exhausted", "time_exhausted"],
        ke_eV=[0.7, 0.3, 0.05],
    )


def _pooled_ab():
    return make_detection(
        n_detected=[3.0, 7.0, 1.0, 12.0, 3.0, 2.0, 12.0],
        state_reason=["time_exhausted", "suppressed", "time_exhausted",
                      "droplet_retained", "time_exhausted", "time_exhausted",
                      "time_exhausted"],
        ke_eV=[0.5, 0.2, 1.0, 0.0, 0.7, 0.3, 0.05],
    )


class TestMixtureRead:
    def test_ion_count_weights_reproduce_pooled_read_exactly(self):
        read_a = read_confirmation_detection(_cell_a(), label="a")
        read_b = read_confirmation_detection(_cell_b(), label="b")
        pooled = read_confirmation_detection(_pooled_ab(), label="ab")

        # weights proportional to source-ion counts == literal pooling
        mix = mixture_read([read_a, read_b], np.array([4.0, 3.0]),
                           label="mix")

        assert mix.trapped_frac == pytest.approx(pooled.trapped_frac,
                                                 abs=1e-14)
        assert mix.trap_bound_frac == pytest.approx(pooled.trap_bound_frac,
                                                    abs=1e-14)
        assert mix.trap_marginal_frac == pytest.approx(
            pooled.trap_marginal_frac, abs=1e-14
        )
        assert mix.suppressed_frac == pytest.approx(pooled.suppressed_frac,
                                                    abs=1e-14)
        np.testing.assert_allclose(mix.fraction, pooled.fraction, atol=1e-14)
        assert mix.n_mean == pytest.approx(pooled.n_mean, abs=1e-12)
        assert mix.n1_frac == pytest.approx(pooled.n1_frac, abs=1e-14)
        assert mix.detected_yield == pytest.approx(
            pooled.num_scored / pooled.num_ions, abs=1e-14
        )

        kpb_mix = mix.ke_by_n()
        kpb_pool = pooled.ke_by_n()
        np.testing.assert_array_equal(kpb_mix.n, kpb_pool.n)
        np.testing.assert_allclose(kpb_mix.mean_eV, kpb_pool.mean_eV,
                                   atol=1e-14)
        np.testing.assert_array_equal(kpb_mix.count, kpb_pool.count)

    def test_committed_scorers_agree_with_pooled_read(self):
        """The duck-type contract: histogram score and both band ratios give
        the same numbers on the mixture as on the literally pooled read."""
        read_a = read_confirmation_detection(_cell_a(), label="a")
        read_b = read_confirmation_detection(_cell_b(), label="b")
        pooled = read_confirmation_detection(_pooled_ab(), label="ab")
        mix = mixture_read([read_a, read_b], np.array([4.0, 3.0]),
                           label="mix")

        n_ref = np.arange(1, 18)
        abundance = make_abundance_reference(
            n_ref, np.full(n_ref.size, 1.0 / n_ref.size)
        )
        ked = make_ked_reference(n_ref, np.full(n_ref.size, 0.5))

        h_mix = score_histogram_vs_reference(mix, abundance)
        h_pool = score_histogram_vs_reference(pooled, abundance)
        assert h_mix.w1_solvated == pytest.approx(h_pool.w1_solvated,
                                                  abs=1e-12)
        assert h_mix.sim_n1 == pytest.approx(h_pool.sim_n1, abs=1e-14)

        for fn in (midhot_ratio, deep_bin_ke_ratio):
            band_mix = fn(mix, ked)
            band_pool = fn(pooled, ked)
            np.testing.assert_array_equal(band_mix.n, band_pool.n)
            assert band_mix.value == pytest.approx(band_pool.value, abs=1e-12)

    def test_degenerate_weight_reproduces_single_cell(self):
        read_a = read_confirmation_detection(_cell_a(), label="a")
        read_b = read_confirmation_detection(_cell_b(), label="b")
        mix = mixture_read([read_a, read_b], np.array([1.0, 0.0]),
                           label="only-a")
        np.testing.assert_allclose(mix.fraction, read_a.fraction, atol=1e-14)
        assert mix.trapped_frac == pytest.approx(read_a.trapped_frac)
        assert mix.n_mean == pytest.approx(read_a.n_mean)

    def test_marginal_class_fraction_is_weighted(self):
        det = _cell_a()
        reasons = np.asarray(det.state_reason, dtype="<U32").copy()
        reasons[3] = "droplet_retained_marginal"   # the U16 builder truncates
        det_marg = dataclasses.replace(det, state_reason=reasons)
        read_m = read_confirmation_detection(det_marg, label="m")
        read_b = read_confirmation_detection(_cell_b(), label="b")
        assert read_m.trap_marginal_frac == pytest.approx(0.25)

        mix = mixture_read([read_m, read_b], np.array([0.5, 0.5]),
                           label="mix")
        assert mix.trap_marginal_frac == pytest.approx(0.125)
        assert mix.trap_bound_frac == 0.0
        assert mix.trapped_frac == pytest.approx(0.125)

    def test_min_count_guard_fails_loudly(self):
        read_a = read_confirmation_detection(_cell_a(), label="a")
        mix = mixture_read([read_a], np.array([1.0]), label="mix")
        with pytest.raises(ValueError, match="min_count=1 only"):
            mix.ke_by_n(min_count=2)

    def test_input_validation(self):
        read_a = read_confirmation_detection(_cell_a(), label="a")
        read_b = read_confirmation_detection(_cell_b(), label="b")
        with pytest.raises(ValueError, match="reads vs"):
            mixture_read([read_a, read_b], np.array([1.0]), label="mix")
        with pytest.raises(ValueError, match="non-negative"):
            mixture_read([read_a, read_b], np.array([1.0, -0.5]),
                         label="mix")
        with pytest.raises(ValueError, match="positive sum"):
            mixture_read([read_a, read_b], np.array([0.0, 0.0]),
                         label="mix")
