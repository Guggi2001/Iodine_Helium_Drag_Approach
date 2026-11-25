"""Generated rung tables of the dissociation ladder (moved with the functions, P2)."""

from __future__ import annotations

import pytest

from i2_helium_md.physics import dissociation_ladder as ladder
from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum


def test_form_u_rungs_match_d0_of_n_and_pin_sigma21():
    """The base table is the mixture kappa=1 Form-U ladder read through
    d0_of_n(1..N_STAR+11) -- the Form-U cache-height convention (Slice T2).
    Numerical anchor: Sigma(21) = 0.18783720 eV (the probe-program wiring
    oracle, findings §4j)."""
    rungs = ladder.form_u_rungs_eV()
    assert len(rungs) == N_STAR + 11 == 32
    for i, rung in enumerate(rungs, start=1):
        assert rung == pytest.approx(
            d0_of_n(i, picture="statistical_mixture", kappa=1.0)
        )
    assert sum(rungs[:N_STAR]) == pytest.approx(0.18783720, abs=1e-7)
    assert sum(rungs[:N_STAR]) == pytest.approx(
        ladder_cumsum(N_STAR, picture="statistical_mixture", kappa=1.0)
    )


def test_rq4graded_table_tapers_first_three_rungs():
    """rq4graded = the RQ4-graded diagnostic taper 2.2 : 1.5 : 1.3 on rungs
    1-3 of the flat mixture bottom, rest = Form-U verbatim (findings I51)."""
    base = ladder.form_u_rungs_eV()
    graded = ladder.rq4graded_rungs_eV()
    assert len(graded) == len(base)
    assert ladder.RQ4GRADED_TAPER_MULTIPLIERS == (2.2, 1.5, 1.3)
    for i, mult in enumerate(ladder.RQ4GRADED_TAPER_MULTIPLIERS):
        assert graded[i] == pytest.approx(mult * base[i])
    assert graded[3:] == base[3:]
