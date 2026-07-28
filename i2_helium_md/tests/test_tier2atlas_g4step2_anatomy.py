"""Tests for the G4 Step 2 Block D W₁-anatomy helpers (plan §3.5f).

Covers `scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py`:

- ``top_share_bins`` — the minimal bin set carrying >= share of W1
  (sorted by (-|gap|, n); ValueError outside 0 < share <= 1);
- ``solvated_gap_profile`` — the identity oracle in unit form: the gap
  profile on the RQ8 solvated branch reproduces the committed
  ``w1_solvated`` to 1e-12.

Synthetic inputs only — no run dirs, no figures (testing rules).
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.postprocess.distribution_compare import CdfGapProfile
from i2_helium_md.postprocess.tier2_confirmation import (
    ConfirmationDetectionRead,
    score_histogram_vs_reference,
)
from scripts.post_processing.tier2atlas_g4step2_w1_anatomy import (
    solvated_gap_profile,
    top_share_bins,
)
from tests.test_tier2_confirmation import make_abundance_reference


def _profile(gaps):
    g = np.asarray(gaps, dtype=float)
    z = np.zeros_like(g)
    return CdfGapProfile(support=np.arange(1, g.size + 1),
                         f_sim=z, f_ref=z, gaps=g)


class TestTopShareBins:
    def test_minimal_set_and_share(self):
        prof = _profile([0.5, -0.3, 0.1, -0.1])   # W1 = 1.0
        bins, reached = top_share_bins(prof, 0.70)
        assert list(bins) == [1, 2]
        assert reached == pytest.approx(0.8)

    def test_tie_broken_by_n(self):
        prof = _profile([0.4, -0.4, 0.2])
        bins, _ = top_share_bins(prof, 0.5)
        assert list(bins) == [1, 2]   # equal |gap|: lower n first

    def test_share_one_returns_all_nonzero(self):
        prof = _profile([0.4, 0.0, -0.6])
        bins, reached = top_share_bins(prof, 1.0)
        assert list(bins) == [1, 3]   # zero-gap bin excluded
        assert reached == pytest.approx(1.0)

    def test_bins_ascending(self):
        prof = _profile([0.1, -0.5, 0.2, -0.4])
        bins, _ = top_share_bins(prof, 0.7)
        assert list(bins) == sorted(bins)

    def test_bad_share_raises(self):
        with pytest.raises(ValueError):
            top_share_bins(_profile([1.0]), 0.0)
        with pytest.raises(ValueError):
            top_share_bins(_profile([1.0]), 1.5)


def _read(n_scored, ke, *, num_ions, n_max=25):
    counts = np.bincount(n_scored, minlength=n_max + 1).astype(float)
    return ConfirmationDetectionRead(
        label="m", num_ions=num_ions, num_scored=len(n_scored),
        trapped_frac=0.0, trap_bound_frac=0.0, trap_marginal_frac=0.0,
        suppressed_frac=0.0,
        n_scored=np.asarray(n_scored), ke_scored_eV=np.asarray(ke),
        n_values=np.arange(n_max + 1), fraction=counts / counts.sum(),
        n_mean=float(np.mean(n_scored)),
        n1_frac=float(counts[1] / counts.sum()),
    )


class TestSolvatedGapIdentity:
    def test_matches_committed_w1_on_synthetic(self):
        # a read with bare (n = 0) mass, so the solvated renormalization is
        # actually exercised; a reference over n = 1..6.
        read = _read([0, 1, 1, 2, 3, 5, 5, 6], [0.1] * 8, num_ions=8)
        ref = make_abundance_reference(
            [1, 2, 3, 4, 5, 6], [0.24, 0.22, 0.20, 0.14, 0.12, 0.08]
        )
        prof = solvated_gap_profile(read, ref)
        committed = score_histogram_vs_reference(read, ref).w1_solvated
        assert prof.w1 == pytest.approx(committed, abs=1e-12)

    def test_signed_direction(self):
        # all sim mass at n = 1 vs all ref mass at n = 2: sim CDF leads (+).
        read = _read([1, 1], [0.1, 0.1], num_ions=2)
        ref = make_abundance_reference([1, 2], [0.0, 1.0])
        prof = solvated_gap_profile(read, ref)
        assert prof.gaps[0] == pytest.approx(1.0)
