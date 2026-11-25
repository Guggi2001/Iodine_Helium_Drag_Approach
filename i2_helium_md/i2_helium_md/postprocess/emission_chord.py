"""Path length from the birth point to the droplet surface along the explosion axis."""

from __future__ import annotations

from pathlib import Path

import numpy as np

#: Chord bin edges [Angstrom] of the chord maps.
CHORD_EDGES_ANGSTROM = (0, 25, 30, 33, 35, 37, 39, 41, 44, 48, 55, 65, 80, 200)


def emission_chords(run_dir: Path) -> np.ndarray:
    """Straight-line path [Å] from each fragment's birth point to the droplet
    surface along the Coulomb-explosion axis (partner -> fragment).

    Uses the last stored neutral-stage column (= the first ion-stage column,
    bit-identical); fragment ``i`` pairs with ``i ± num_molecules``.
    Chord = -p·ê + sqrt((p·ê)² - |p|² + R²) for birth position p, unit axis
    ê, droplet radius R. Ignores deflection and the soft erf edge of the
    density gate. Shape (2 * num_molecules,).
    """
    with np.load(run_dir / "neutral.npz", allow_pickle=False) as z:
        R = np.asarray(z["droplet_radii"], dtype=float)
        P = np.stack([z[k][:, -1] for k in ("positions_x", "positions_y",
                                            "positions_z")], axis=1)
    n_mol = R.size // 2
    partner = np.r_[np.arange(n_mol, 2 * n_mol), np.arange(n_mol)]
    axis = P - P[partner]
    axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    pe = np.sum(P * axis, axis=1)
    return -pe + np.sqrt(pe ** 2 - np.sum(P * P, axis=1) + R ** 2)
