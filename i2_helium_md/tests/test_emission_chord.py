"""Straight path from the birth point to the droplet surface along the explosion axis."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.postprocess.emission_chord import CHORD_EDGES_ANGSTROM, emission_chords


def run_with(tmp_path, positions, radius):
    """Run directory with one molecule; first stored column is a decoy."""
    p = np.asarray(positions, dtype=float)
    np.savez(tmp_path / "neutral.npz", droplet_radii=np.full(2, radius),
             **{f"positions_{k}": np.stack([np.full(2, 99.0), p[:, i]], axis=1)
                for i, k in enumerate("xyz")})
    return tmp_path


@pytest.mark.parametrize("positions, radius, chord", [
    ([(10, 0, 0), (-10, 0, 0)], 50.0, 40.0),
    ([(3, 4, 0), (-3, 4, 0)], 5.0, 0.0),
])
def test_chord_of_both_fragments(tmp_path, positions, radius, chord):
    got = emission_chords(run_with(tmp_path, positions, radius))
    np.testing.assert_allclose(got, [chord, chord], atol=1e-12)


def test_off_centre_fragments_get_their_own_chords(tmp_path):
    got = emission_chords(run_with(tmp_path, [(30, 0, 0), (10, 0, 0)], 50.0))
    np.testing.assert_allclose(got, [20.0, 60.0], atol=1e-12)


def test_bin_edges():
    assert CHORD_EDGES_ANGSTROM == (0, 25, 30, 33, 35, 37, 39, 41, 44, 48, 55, 65, 80, 200)
