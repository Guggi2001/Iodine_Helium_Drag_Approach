"""Form-U dissociation ladder ``D_0(n)`` and the integrated self-bound gate
threshold ``Sigma(n)`` for the I+He_n complex (Tier-2 Phase-A Slice L).

Pure, stateless energetics: no physics state, no RNG, no integrator, no mass and
no velocity (Phase A is mass-agnostic; the friction/force-coefficient ``gamma`` is
never touched here). This is the Tier-2 analog of the Tier-1a ``shell_schedule.py``
"first, fully-independent build unit" -- everything stateful (the pickup /
evaporation channels, Phase B) and composed (the generative driver, Phase C)
*consumes* this module and lives later.

Encoded form (TIER2_PHASE_A_IMPLEMENTATION_PLAN.md §2.2 (L); MASS §8 R3 / §11)
-----------------------------------------------------------------------------
The single-rung dissociation cost interpolates the picture-keyed first rung
``D_0(1)`` down to the bulk-He floor ``D_floor`` across a sigmoid cliff centred at
``n* + 1/2`` (between the last in-shell atom ``n*`` and the first shell-2 atom)::

    D_0(n) = D_floor + (D_0(1) - D_floor) * (1 - sigma(n)) / (1 - sigma(1))
    sigma(n) = [1 + exp(-kappa * (n - n* - 1/2))]^(-1)
    Sigma(n) = sum_{i=1}^{n} D_0(i)          (the cumulative gate threshold)

The normalisation ``(1 - sigma(1))`` pins ``D_0(1)`` to the sourced first rung
*exactly* for every ``kappa`` (the rung is picture-set; ``kappa`` only shapes the
cliff). ``Sigma(n)`` is the cumulative cost an evaporation must overcome to stay
self-bound; ``Sigma(21)`` is ~kappa-independent (~11% over kappa in [0.3, 5]) and
is therefore pinnable ahead of the kappa fit (the decoupling property).

Electronic picture (the Tier-2 co-fit knob, OQ1 carry-both)
-----------------------------------------------------------
``picture`` selects the first rung ``D_0(1)``:

* ``statistical_mixture`` (default) -- 74.4 cm^-1, the SO-statistical (X2+I1+I0)/3;
* ``x2_only``                       -- 106.9 cm^-1, the X2 exact J=0 ZPE;
* ``cooling_relaxed``               -- live but **provisional**: the arithmetic
  mean of the two pinned rungs. The rule-2 carry is the concrete
  relaxation-weighted blend *value* (pinned at Phase F), not the enum arm --
  the arm is read. Only the ordering mix < cooling_relaxed < X2 is asserted;
  do not read its value as a 4-figure oracle.

Units follow the package convention: energies in eV (``*_eV``); ``n`` and
``kappa`` are dimensionless (``kappa`` is per-unit-``n``). Vectorised scalar-in
-> float-out / array-in -> ndarray-out, mirroring ``shell_schedule.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .constants import (
    D0_1_COOLING_RELAXED_WAVENUMBER,
    D0_1_MIX_WAVENUMBER,
    D0_1_X2_WAVENUMBER,
    D_FLOOR_WAVENUMBER,
    EV_PER_WAVENUMBER,
    N_STAR,
)

# Picture -> first-rung dissociation energy [eV]. Single source of truth for the
# picture keying; ``config.LadderElectronicPicture`` must list the same members.
_FIRST_RUNG_EV: dict[str, float] = {
    "statistical_mixture": D0_1_MIX_WAVENUMBER * EV_PER_WAVENUMBER,
    "x2_only": D0_1_X2_WAVENUMBER * EV_PER_WAVENUMBER,
    "cooling_relaxed": D0_1_COOLING_RELAXED_WAVENUMBER * EV_PER_WAVENUMBER,
}

#: Bulk-He floor [eV] -- the asymptote ``D_0(n)`` relaxes to past the shell edge.
D_FLOOR_EV: float = D_FLOOR_WAVENUMBER * EV_PER_WAVENUMBER

#: Cliff centre on the ``n`` axis (between the last in-shell atom and shell 2).
CLIFF_CENTER: float = N_STAR + 0.5


def first_rung_d0_eV(picture: str = "statistical_mixture") -> float:
    """First-rung dissociation energy ``D_0(1)`` [eV] for ``picture``.

    The picture-keyed sourced anchor; ``D_0(1)`` is ``kappa``-independent by the
    Form-U normalisation. Raises ``ValueError`` on an unknown picture (fail-loud;
    mirrors the ``mass_scenario`` reject arm in ``config.check_drag_config``).
    """
    try:
        return _FIRST_RUNG_EV[picture]
    except KeyError:
        raise ValueError(
            f"unknown ladder electronic picture {picture!r}; expected one of "
            f"{tuple(_FIRST_RUNG_EV)}"
        ) from None


def sigma(n, *, kappa: float):
    """Centred Form-U sigmoid ``sigma(n) = [1 + exp(-kappa*(n - n* - 1/2))]^-1``.

    Dimensionless, monotone non-decreasing in ``n``; ``sigma(n* + 1/2) = 1/2``
    exactly for any ``kappa`` (the cliff centre). Vectorised: scalar-in -> float,
    array-in -> ndarray (same shape).

    Implemented in the algebraically-identical ``0.5*(1 + tanh(x/2))`` form (with
    ``x = kappa*(n - CLIFF_CENTER)``): ``tanh`` saturates to +/-1 instead of
    overflowing, so a sharp cliff (large ``kappa``) stays finite and warning-free
    where the bare ``1/(1+exp(-x))`` would form an ``inf`` intermediate.
    """
    n_arr = np.asarray(n, dtype=float)
    out = 0.5 * (1.0 + np.tanh(0.5 * kappa * (n_arr - CLIFF_CENTER)))
    return float(out) if np.ndim(n) == 0 else out


def d0_of_n(n, *, picture: str = "statistical_mixture", kappa: float, ladder=None):
    """Single-rung dissociation cost ``D_0(n)`` [eV] (Form U or injected table).

    Parameters
    ----------
    n : int or array-like
        Shell occupancy / rung index (>= 1 for a physical rung).
    picture : str, optional
        Electronic picture selecting ``D_0(1)`` (default ``statistical_mixture``).
    kappa : float
        Cliff steepness (per-unit-``n``); the Free Tier-2 knob.
    ladder : TabulatedLadder, optional
        Injected non-parametric ladder (Slice T2, §I.10): when set, the lookup
        delegates to ``ladder.d0_of_n(n)`` and **``picture``/``kappa`` are
        ignored** (the table replaces the Form-U parametrisation entirely).
        ``None`` (default) is the byte-inert Form-U path. Consumers thread this
        kwarg; the config bridge is :func:`resolve_ladder`.

    Returns
    -------
    float or ndarray
        ``D_0(n)`` in eV. On the Form-U path: monotone non-increasing in ``n``
        and bounded to ``[D_floor, D_0(1)]`` for ``n >= 1``. Scalar-in -> float,
        array-in -> ndarray of the same shape. The tabulated path additionally
        requires integer ``n`` within the table (fail-loud).
    """
    if ladder is not None:
        return ladder.d0_of_n(n)
    d1 = first_rung_d0_eV(picture)
    if not (d1 > D_FLOOR_EV):  # defensive: the sourced rungs all satisfy this
        raise ValueError(
            f"D_0(1)={d1!r} eV must exceed the floor D_floor={D_FLOOR_EV!r} eV "
            f"for picture {picture!r}"
        )
    s_n = sigma(n, kappa=kappa)
    s_1 = sigma(1.0, kappa=kappa)
    out = D_FLOOR_EV + (d1 - D_FLOOR_EV) * (1.0 - np.asarray(s_n)) / (1.0 - s_1)
    return float(out) if np.ndim(n) == 0 else out


def _as_integer_occupancy(n, *, context: str) -> np.ndarray:
    """Coerce occupancy ``n`` to an integer ndarray (fail-loud on fractional).

    A cumulative/indexed ladder lookup is defined over integer occupancy only:
    integer-valued floats are accepted and cast; a genuinely fractional ``n``
    raises a loud ``ValueError`` (a bare integer-indexed table lookup would die
    with numpy's cryptic ``IndexError`` instead). Shared by the Form-U
    :func:`ladder_cumsum` and both :class:`TabulatedLadder` lookups so the
    fail-loud contract stays single-source.
    """
    n_arr = np.asarray(n)
    if not np.issubdtype(n_arr.dtype, np.integer):
        if np.any(n_arr != np.floor(n_arr)):
            raise ValueError(f"{context} requires integer occupancy n; got {n!r}")
        n_arr = n_arr.astype(int)
    return n_arr


@lru_cache(maxsize=64)
def _sigma_prefix_table(picture: str, kappa: float, n_top: int) -> np.ndarray:
    """Cached ``Sigma`` prefix table ``cum[k] = Sigma(k)`` for ``k in [0, n_top]``.

    Built exactly as the pre-cache :func:`ladder_cumsum` body did (same rung
    construction, same sequential ``np.cumsum``), so every indexed prefix value
    is **bit-identical** to a fresh computation -- prefix sums do not depend on
    later rungs, so one table serves every request with ``n <= n_top``. This is
    the hot-path fix for the Slice-G driver, which evaluates ``Sigma`` several
    times per step for a ``(picture, kappa)`` pair fixed over the whole run
    (~1e5 steps); the whole output space is a ~22-entry table. Read-only (the
    cached array is shared across callers).

    Cache-coherence hazard: the key is ``(picture, kappa, n_top)`` only, but the
    table is built through this module's global ``d0_of_n`` / ``_FIRST_RUNG_EV``.
    Monkeypatching either *in this module's namespace* after a table is cached
    serves stale values -- patch the consumer module's namespace instead (as the
    suites do), or call ``_sigma_prefix_table.cache_clear()``.
    """
    rungs = np.atleast_1d(
        d0_of_n(np.arange(1, n_top + 1), picture=picture, kappa=kappa)
    )
    cum = np.concatenate(([0.0], np.cumsum(rungs)))
    cum.setflags(write=False)
    return cum


#: Default cached-table height: covers the full first shell (n* = 21) plus the
#: modest overshoot the generative channels can reach, so the production run
#: hits one cache entry. Larger ``n`` transparently builds a taller table.
_SIGMA_TABLE_N_TOP: int = N_STAR + 11


def ladder_cumsum(n, *, picture: str = "statistical_mixture", kappa: float,
                  ladder=None):
    """Cumulative gate threshold ``Sigma(n) = sum_{i=1}^{n} D_0(i)`` [eV].

    ``Sigma(0) = 0``. Consumed by Slice K (occupancy-resolved asymptote) and Slice
    U (self-unbound floor). Vectorised: scalar-in -> float, array-in -> ndarray.
    Negative ``n`` is rejected (fail-loud). Values come from the cached
    :func:`_sigma_prefix_table` (bit-identical to a fresh computation; see there).

    ``ladder`` (Slice T2, §I.10): when set, delegates to
    ``ladder.ladder_cumsum(n)`` and ``picture``/``kappa`` are ignored -- see
    :func:`d0_of_n`.
    """
    if ladder is not None:
        return ladder.ladder_cumsum(n)
    n_arr = np.asarray(n)
    if np.any(n_arr < 0):
        raise ValueError(f"ladder_cumsum requires n >= 0; got {n!r}")
    n_arr = _as_integer_occupancy(n, context="ladder_cumsum")
    n_max = int(n_arr.max(initial=0))
    # cum[k] = Sigma(k): prepend 0 so cum[0] == 0, then index by n. The table is
    # padded to a fixed height so the hot path reuses one cache entry.
    cum = _sigma_prefix_table(picture, float(kappa), max(n_max, _SIGMA_TABLE_N_TOP))
    out = cum[n_arr]
    return float(out) if np.ndim(n) == 0 else out


def gate_threshold(n, *, picture: str = "statistical_mixture", kappa: float,
                   ladder=None):
    """Readable alias for :func:`ladder_cumsum` -- the self-bound gate threshold."""
    return ladder_cumsum(n, picture=picture, kappa=kappa, ladder=ladder)


@dataclass(frozen=True)
class TabulatedLadder:
    """A non-parametric ladder built from an explicit per-rung ``D_0`` table.

    The declared Form-U fallback: a non-monotone / second-cliff histogram that the
    sigmoid cannot represent would be supplied here instead. ``rungs_eV[i-1]`` is
    ``D_0(i)`` (1-indexed rungs); not the default path.
    """

    rungs_eV: tuple[float, ...]

    def d0_of_n(self, n):
        """``D_0(n)`` [eV] by table lookup (1-indexed).

        Raises on out-of-range or fractional ``n`` (integer-valued floats are
        accepted and cast -- the :func:`ladder_cumsum` contract).
        """
        n_arr = _as_integer_occupancy(n, context="TabulatedLadder.d0_of_n")
        if np.any(n_arr < 1) or np.any(n_arr > len(self.rungs_eV)):
            raise ValueError(
                f"n must be in [1, {len(self.rungs_eV)}] for this tabulated "
                f"ladder; got {n!r}"
            )
        table = np.asarray(self.rungs_eV, dtype=float)
        out = table[n_arr - 1]
        return float(out) if np.ndim(n) == 0 else out

    def ladder_cumsum(self, n):
        """``Sigma(n)`` [eV] by cumulative table lookup; ``Sigma(0) = 0``.

        Raises on out-of-range or fractional ``n`` (integer-valued floats are
        accepted and cast -- the :func:`ladder_cumsum` contract).
        """
        n_arr = _as_integer_occupancy(n, context="TabulatedLadder.ladder_cumsum")
        if np.any(n_arr < 0) or np.any(n_arr > len(self.rungs_eV)):
            raise ValueError(
                f"n must be in [0, {len(self.rungs_eV)}] for this tabulated "
                f"ladder; got {n!r}"
            )
        cum = np.concatenate(([0.0], np.cumsum(np.asarray(self.rungs_eV, dtype=float))))
        out = cum[n_arr]
        return float(out) if np.ndim(n) == 0 else out


def tabulated_ladder(rungs_eV) -> TabulatedLadder:
    """Build a :class:`TabulatedLadder` from a 1-indexed per-rung ``D_0`` table [eV].

    The declared Form-U fallback (not the default path). Round-trips the supplied
    rungs and their cumulative sum exactly.
    """
    rungs = tuple(float(x) for x in rungs_eV)
    if len(rungs) == 0:
        raise ValueError("tabulated_ladder requires at least one rung")
    return TabulatedLadder(rungs_eV=rungs)


def resolve_ladder(dissociation_ladder: str, rungs_eV=None):
    """Resolve the config ladder surface to an injectable ladder object (Slice T2).

    The single ``(cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV)`` ->
    ``ladder`` bridge: ``config.check_ladder_config`` calls it at config load and
    the three stages (ion ``biphasic_step`` / relaxation / detection) call it at
    point-of-use, so an unvalidated config view still fails loudly here rather
    than silently running Form-U physics (the pre-T2 relaxation-stage hazard).

    Parameters
    ----------
    dissociation_ladder : str
        The config selector (``"form_u"`` or ``"tabulated"``).
    rungs_eV : sequence of float, optional
        The 1-indexed per-rung ``D_0`` table [eV] (``rungs_eV[i-1] = D_0(i)``);
        the ``cfg.tabulated_ladder_rungs_eV`` payload.

    Returns
    -------
    TabulatedLadder or None
        ``None`` for ``form_u`` (consumers use the Form-U module functions via
        ``picture``/``kappa``); a :class:`TabulatedLadder` for ``tabulated``.

    Raises
    ------
    ValueError
        * ``"tabulated"`` without a table (no silent Form-U fallback);
        * a table under ``"form_u"`` (a silently-ignored table is a stale-intent
          hazard -- the off-diagonal pairing is refused, mirroring the drag
          form/coefficient cross-check);
        * fewer than ``N_STAR`` rungs (``Sigma(n*)`` must be table-covered: the
          solvation split normalises by it). Longer tables cover the pickup
          overshoot band (the Form-U cache height is ``N_STAR + 11``); a runtime
          lookup past the table stays a loud :class:`TabulatedLadder` error --
          the accepted fail-loud convention (discussion 2026-07-16);
        * any non-positive or non-finite rung (a ``D_0 <= 0`` rung is
          unphysical: it opens the RRK bracket at zero cost).
    """
    if dissociation_ladder == "form_u":
        if rungs_eV is not None:
            raise ValueError(
                "tabulated_ladder_rungs_eV is set but dissociation_ladder="
                "'form_u' would silently ignore it; select "
                "dissociation_ladder='tabulated' or drop the table."
            )
        return None
    if dissociation_ladder == "tabulated":
        if rungs_eV is None:
            raise ValueError(
                "dissociation_ladder='tabulated' requires "
                "tabulated_ladder_rungs_eV (the 1-indexed per-rung D_0 table "
                "[eV]); refusing to silently run the Form-U ladder."
            )
        ladder = tabulated_ladder(rungs_eV)
        if len(ladder.rungs_eV) < N_STAR:
            raise ValueError(
                f"tabulated_ladder_rungs_eV needs at least n* = {N_STAR} rungs "
                f"(Sigma(n*) must be table-covered); got "
                f"{len(ladder.rungs_eV)}."
            )
        rungs_arr = np.asarray(ladder.rungs_eV, dtype=float)
        if not np.all(np.isfinite(rungs_arr) & (rungs_arr > 0.0)):
            raise ValueError(
                "tabulated_ladder_rungs_eV entries must be positive finite "
                f"energies [eV]; got {ladder.rungs_eV!r}."
            )
        return ladder
    raise ValueError(
        f"unknown dissociation_ladder {dissociation_ladder!r}; expected one of "
        "('form_u', 'tabulated')"
    )
