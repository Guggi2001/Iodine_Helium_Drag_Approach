"""``DetectedEnsembleView`` — the checkpoint-shaped detected ensemble
(`TIER2_DETECTION_SUMMARY_SPEC.md` §3, as amended 2026-07-22).

The view duck-types exactly the ``IonCheckpoint`` surface the shared
mass-gated diagnostics read (``mass_final_kg`` (2N,),
``velocities_final_{x,y,z}`` (2N,) [A/ps], ``b_ion_outside`` (N,),
``num_molecules``), so every legacy plotting recipe runs verbatim on the
detection ensemble. Forcing rules under test:

* frozen / time_exhausted — ``mass_detected_kg`` pass-through, which is
  exactly ``m(n_detected)`` for the discrete-shed mechanism (the
  n_detected <-> mass-gate equivalence pin);
* suppressed — forced to the bare gate ``m(0)`` (the frozen
  suppressed -> n = 0 convention; they ride at handover n on disk);
* droplet_retained — forced to NaN (matches no gate; per-fragment
  exclusion, incl. pair panels via the pair-AND).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from tests import test_tier2_confirmation as t2c

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.physics.shell_schedule import complex_mass_amu
from i2_helium_md.postprocess.detected_view import (
    DetectedEnsembleView,
    detected_ensemble_view,
)
from i2_helium_md.postprocess.ihe_ked import fragment_mean_kinetic_energy
from i2_helium_md.postprocess.velocity_distribution import (
    compute_final_velocity_histogram,
    select_final_mass_gate,
)


def make_rich_detection(n_detected, state_reason, *, vx, vy, vz):
    """Synthetic ``DetectionResult`` with physical masses and velocities.

    ``mass_detected_kg`` follows the on-disk convention: ``m(n_detected)``
    for every row (suppressed rows ride at their handover ``n``, which IS
    their ``n_detected`` on disk; retained rows hold the verbatim handover
    state). Velocities in A/ps.
    """
    det = t2c.make_detection(
        n_detected=n_detected,
        state_reason=state_reason,
        ke_eV=[0.0] * len(n_detected),
    )
    mass_kg = complex_mass_amu(np.asarray(n_detected, dtype=float)) * U_KG
    return replace(
        det,
        mass_detected_kg=np.asarray(mass_kg, dtype=float),
        vx_detected=np.asarray(vx, dtype=float),
        vy_detected=np.asarray(vy, dtype=float),
        vz_detected=np.asarray(vz, dtype=float),
    )


def _four_ion_detection():
    # 2 molecules / 4 ions: frozen n=2, suppressed handover n=5,
    # droplet_retained handover n=9, time_exhausted n=1.
    return make_rich_detection(
        [2.0, 5.0, 9.0, 1.0],
        ["frozen", "suppressed", "droplet_retained", "time_exhausted"],
        vx=[1.0, 2.0, 3.0, 4.0],
        vy=[0.5, 0.6, 0.7, 0.8],
        vz=[0.0, -0.1, -0.2, -0.3],
    )


class TestViewShapesAndPassthrough:
    def test_shapes_and_velocity_passthrough(self):
        det = _four_ion_detection()
        view = detected_ensemble_view(det)
        assert isinstance(view, DetectedEnsembleView)
        assert view.num_molecules == 2
        assert view.mass_final_kg.shape == (4,)
        assert view.velocities_final_x.tolist() == [1.0, 2.0, 3.0, 4.0]
        assert view.velocities_final_y.tolist() == [0.5, 0.6, 0.7, 0.8]
        assert view.velocities_final_z.tolist() == [0.0, -0.1, -0.2, -0.3]

    def test_outside_flag_is_all_true(self):
        # exclusion is by mass sentinel, never by the molecule-level flag
        view = detected_ensemble_view(_four_ion_detection())
        assert view.b_ion_outside.shape == (2,)
        assert view.b_ion_outside.dtype == bool
        assert view.b_ion_outside.all()

    def test_shape_mismatch_raises(self):
        det = _four_ion_detection()
        bad = replace(det, vx_detected=np.zeros(3))
        with pytest.raises(ValueError, match="shape"):
            detected_ensemble_view(bad)


class TestMassForcingRules:
    def test_frozen_mass_equals_m_of_n_exactly(self):
        # the n_detected <-> mass-gate equivalence pin (exact, no tolerance:
        # the discrete-shed mechanism writes m(n) bit-exactly)
        view = detected_ensemble_view(_four_ion_detection())
        assert view.mass_final_kg[0] == complex_mass_amu(2) * U_KG
        assert view.mass_final_kg[3] == complex_mass_amu(1) * U_KG

    def test_suppressed_forced_to_bare_gate(self):
        det = _four_ion_detection()
        view = detected_ensemble_view(det)
        assert view.mass_final_kg[1] == complex_mass_amu(0) * U_KG
        in_bare = select_final_mass_gate(view, mass_amu=complex_mass_amu(0))
        in_handover = select_final_mass_gate(view, mass_amu=complex_mass_amu(5))
        assert in_bare[1]
        assert not in_handover[1]

    def test_retained_matches_no_gate(self):
        view = detected_ensemble_view(_four_ion_detection())
        assert np.isnan(view.mass_final_kg[2])
        for n in range(0, 22):
            gate = select_final_mass_gate(view, mass_amu=complex_mass_amu(n))
            assert not gate[2], f"retained ion leaked into the n={n} gate"


class TestExistingRecipesRunVerbatim:
    def test_final_velocity_histogram_counts_gated_rows(self):
        # two n=1 detected ions; the retained sibling is invisible
        det = make_rich_detection(
            [1.0, 1.0, 7.0, 1.0],
            ["frozen", "frozen", "droplet_retained", "time_exhausted"],
            vx=[3.0, 4.0, 99.0, 0.0],
            vy=[4.0, 3.0, 99.0, 5.0],
            vz=[0.0, 0.0, 99.0, 0.0],
        )
        view = detected_ensemble_view(det)
        hist = compute_final_velocity_histogram(
            view, mass_amu=complex_mass_amu(1), num_bins=10, v_max_Aps=10.0
        )
        assert hist.num_atoms_used == 3  # |v| = 5 A/ps three times

    def test_fragment_mean_ke_over_detected_ensemble(self):
        det = make_rich_detection(
            [1.0, 1.0, 7.0, 1.0],
            ["frozen", "frozen", "droplet_retained", "time_exhausted"],
            vx=[3.0, 4.0, 99.0, 0.0],
            vy=[4.0, 3.0, 99.0, 5.0],
            vz=[0.0, 0.0, 99.0, 0.0],
        )
        view = detected_ensemble_view(det)
        mean = fragment_mean_kinetic_energy(view, 1)
        assert mean.num_atoms_used == 3
        # all three at |v| = 5 A/ps -> mean E = 1/2 m(1) v^2, in eV; the
        # helper owns the conversion, so pin via the helper's own speed
        assert mean.v_of_mean_E_mps == pytest.approx(5.0e2, rel=1e-12)
