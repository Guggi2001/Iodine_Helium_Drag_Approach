"""Depth-graded exit stripping — form functions (Tier-2 (C) design, (A) arm).

The outbound-surface-crossing He knockout of
``TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md`` §3.3. At an outbound crossing with
exit speed ``v_x`` [Å/ps] and mechanical shell ``n_x``, each rung
``j = 1…n_x`` (j = 1 innermost) is independently knocked out with
probability::

    P_knock(j) = P₀(v_x) · G(j)
    P₀(v)      = min(1, (v / v_strip)^a)                 [dimensionless]
    G(j)       = 1 / (1 + exp(−(j − j₀) / w_j))          [dimensionless]

``G`` encodes "outer eager, inner shadowed" (large j = outer rung = more
exposed); ``v_strip`` is Sourced ≈ 9.9 Å/ps (the equal-mass max-transfer
full-strip threshold vs the production rq4graded rungs — the §3.5j CF-3
convention); ``a`` is the Bounded velocity-gate steepness (prior 2, ram
∝ v²). This module is **pure form arithmetic** — no RNG, no state, no mass:
the step operator that draws the Bernoullis and books the ledger lives in
``simulation/ion_propagation_step.exit_strip_step`` (the evaporation /
pickup module split).

Units: velocities Å/ps; ``j``, ``j₀``, ``w_j`` in rungs; probabilities
dimensionless in [0, 1].
"""

from __future__ import annotations

import numpy as np


def strip_velocity_gate(
    v_exit_aps, *, v_strip_aps: float, exponent: float,
) -> np.ndarray | float:
    """Velocity gate ``P₀(v) = min(1, (v / v_strip)^a)`` [dimensionless].

    Parameters
    ----------
    v_exit_aps : array_like
        Outbound exit speed(s) [Å/ps], ≥ 0.
    v_strip_aps : float
        Full-strip reference speed [Å/ps], > 0 (Sourced 9.9).
    exponent : float
        Gate steepness a > 0 (Bounded, prior 2).

    Raises
    ------
    ValueError
        On non-positive ``v_strip_aps`` / ``exponent`` or a negative speed.
    """
    if not (np.isfinite(v_strip_aps) and v_strip_aps > 0):
        raise ValueError(f"v_strip_aps must be finite and > 0, got {v_strip_aps!r}")
    if not (np.isfinite(exponent) and exponent > 0):
        raise ValueError(f"exponent must be finite and > 0, got {exponent!r}")
    v = np.asarray(v_exit_aps, dtype=float)
    if np.any(v < 0):
        raise ValueError("exit speed must be >= 0 everywhere")
    out = np.minimum(1.0, (v / v_strip_aps) ** exponent)
    return float(out) if np.ndim(v_exit_aps) == 0 else out


def strip_depth_grading(
    j, *, protect_j0: float, width_rungs: float,
) -> np.ndarray | float:
    """Depth grading ``G(j) = 1 / (1 + exp(−(j − j₀)/w_j))`` [dimensionless].

    Monotone increasing in ``j``: outer rungs (large j) are eagerly knocked,
    the innermost 1–3 are shadowed (the element §3.5j measured as required).

    Parameters
    ----------
    j : array_like
        Rung index/indices (1 = innermost), ≥ 1.
    protect_j0 : float
        Protection depth j₀ [rungs], ≥ 0 (Bounded, P1 box [1.5, 2]).
    width_rungs : float
        Grading width w_j [rungs], > 0 (Bounded, P1 box [0.5, 1]).

    Raises
    ------
    ValueError
        On a non-positive width, a negative ``j₀``, or a rung index < 1.
    """
    if not (np.isfinite(width_rungs) and width_rungs > 0):
        raise ValueError(f"width_rungs must be finite and > 0, got {width_rungs!r}")
    if not (np.isfinite(protect_j0) and protect_j0 >= 0):
        raise ValueError(f"protect_j0 must be finite and >= 0, got {protect_j0!r}")
    j_arr = np.asarray(j, dtype=float)
    if np.any(j_arr < 1):
        raise ValueError("rung index j must be >= 1 (1 = innermost)")
    out = 1.0 / (1.0 + np.exp(-(j_arr - protect_j0) / width_rungs))
    return float(out) if np.ndim(j) == 0 else out


def strip_knock_probabilities(
    n_x: int,
    v_exit_aps: float,
    *,
    v_strip_aps: float,
    exponent: float,
    protect_j0: float,
    width_rungs: float,
) -> np.ndarray:
    """Per-rung knockout probabilities ``P₀(v_x)·G(j)`` for ``j = 1…n_x``.

    Returns an ``(n_x,)`` array (empty for ``n_x = 0``); element ``j-1``
    is rung ``j``'s probability. Validation propagates from the two form
    functions.
    """
    n_int = int(n_x)
    if n_int < 0:
        raise ValueError(f"n_x must be >= 0, got {n_x!r}")
    if n_int == 0:
        return np.zeros(0, dtype=float)
    j = np.arange(1, n_int + 1, dtype=float)
    p0 = strip_velocity_gate(
        float(v_exit_aps), v_strip_aps=v_strip_aps, exponent=exponent,
    )
    return p0 * np.asarray(
        strip_depth_grading(j, protect_j0=protect_j0, width_rungs=width_rungs)
    )
