"""Shared Form-U rung-table builders for the tabulated == form_u oracles.

Review dedup 2026-07-16 (rule 1; the Slice-G test-helper precedent): the
equivalence-oracle feed convention -- the 1-indexed ``D_0(1..length)`` table
at the Form-U cache height ``N_STAR + 11 = 32`` -- lives here once instead of
as a per-file copy, so a change to the convention (length, the ``atleast_1d``
cast) cannot silently leave the per-stage bitwise oracles testing different
feeds.
"""

from __future__ import annotations

import numpy as np

from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.physics.dissociation_ladder import d0_of_n, tabulated_ladder


def form_u_rungs(picture="statistical_mixture", kappa=1.0, length=N_STAR + 11):
    """The Form-U rung table ``D_0(1..length)`` [eV] -- the equivalence-oracle feed."""
    return tuple(
        np.atleast_1d(d0_of_n(np.arange(1, length + 1), picture=picture, kappa=kappa))
    )


def form_u_rungs_for(cfg, length=N_STAR + 11):
    """:func:`form_u_rungs` at ``cfg``'s (picture, kappa)."""
    return form_u_rungs(
        picture=cfg.ladder_electronic_picture,
        kappa=cfg.ladder_steepness,
        length=length,
    )


def form_u_ladder(picture="statistical_mixture", kappa=1.0, length=N_STAR + 11):
    """A :class:`TabulatedLadder` built from the Form-U rungs -- the injectable twin."""
    return tabulated_ladder(form_u_rungs(picture=picture, kappa=kappa, length=length))
