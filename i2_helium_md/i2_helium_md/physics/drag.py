"""Gated drag-force physics for I+ in a helium bubble (Slice 1).

TDDFT-calibrated continuous drag model that will replace the hard-sphere
collision model (:mod:`i2_helium_md.physics.collisions`) for the ion stage.
This module is the pure-physics swap-point content: three config-free,
**mass-free** functions plus the coefficient-bundle *type* they consume.
``collisions.py`` is left intact and importable -- this module is **additive
and parallel**, not a replacement.

Friction convention (unified, ``DRAG_PORT_DESIGN_DECISIONS.md`` §1.2)
---------------------------------------------------------------------
``gamma(v)`` is a **force coefficient** with units ``amu/ps``, defined
``gamma(v) = |F_drag(v)| / v``, so the friction force is ``gamma(v) * v`` with
**no leading mass**. The friction *rate* ``gamma/m`` [1/ps] appears only inside
the BAOAB damping exponent ``e^(-gamma*dt/m)`` (Slice 2, not here). Because
``gamma`` is a force coefficient, this module is **fully mass-agnostic**: it
never takes ``m``. Mass enters only at the integrator's O-step, as one explicit
division by ``m(t)``.

Governing equations (primary form ``linear_cubic``, the only form realised)
---------------------------------------------------------------------------
With ``depth = r_atom - r_droplet`` (negative inside the droplet, positive
outside -- the same convention as :func:`i2_helium_md.physics.potentials.droplet_potential`)::

    g(depth)        = 0.5 * (1 - erf(depth / steepness))      # dimensionless in [0, 1]
    F_drag(v,depth) = g(depth) * (a*v + b*v**3)               # amu*A/ps^2
    gamma(v,depth)  = g(depth) * (a + b*v**2)                 # amu/ps

Units (``a`` in amu/ps, ``b`` in amu*ps/A^2) balance to a force; see
``SLICE1_GOALS_gated_drag_module.md`` §3. The drag opposes motion; these
functions return the **positive magnitude** form and the consumer applies
``-F_drag`` along ``v_hat`` at the integrator (Slice 2).

``gamma`` is exposed via its **closed form** ``g*(a + b*v**2)``, *never* via
``|F_drag|/v``: the two are analytically equal, but the division manufactures a
``0/0`` singularity at ``v -> 0`` that the ``linear_cubic`` form does not have
(``gamma -> g*a`` there). This is a physics-definition point, made explicit per
form -- it is exactly the ``v -> 0`` boundary where ``power_law`` (n<0) genuinely
*does* diverge and would need the §3.8 floor.

The other three forms (``linear_quadratic``, ``threshold``, ``power_law``) are
reserved behind the same dispatch and raise :class:`NotImplementedError`; see
``SLICE1_GOALS_gated_drag_module.md`` §5. Adding them later is a branch, not a
signature change.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from scipy.special import erf


# Drag-form tags. Only LINEAR_CUBIC is realised in Slice 1; the rest are
# reserved behind the dispatch (see :func:`_raise_unrealised_form`).
LINEAR_CUBIC = "linear_cubic"
LINEAR_QUADRATIC = "linear_quadratic"
THRESHOLD = "threshold"
POWER_LAW = "power_law"

# Required coefficient keys per form (variable arity by form, §3.8).
_REQUIRED_COEFF_KEYS: dict[str, tuple[str, ...]] = {
    LINEAR_CUBIC: ("a", "b"),            # amu/ps, amu*ps/A^2
    LINEAR_QUADRATIC: ("a", "c"),        # amu/ps, amu/A
    THRESHOLD: ("F_sat", "v0"),          # amu*A/ps^2, A/ps
    POWER_LAW: ("gamma", "n"),           # amu*A^(1-n)*ps^(n-2), dimensionless
}

_VALID_MASS_MODELS = ("constant", "time_resolved")


@dataclass(frozen=True)
class DragCoefficients:
    """Form-tagged, mass-stamped drag-coefficient bundle.

    Slice 1 owns this *type*; the ``SimConfig`` enum surface and the §6.5
    ``mass_scenario`` <-> ``drag_coefficients`` consistency guard are Slice 3.
    The module **consumes** a bundle and never constructs one from config.

    Attributes
    ----------
    form : str
        Drag-form tag, one of ``{LINEAR_CUBIC, LINEAR_QUADRATIC, THRESHOLD,
        POWER_LAW}``. Only ``LINEAR_CUBIC`` is realised in Slice 1.
    coefficients : Mapping[str, float]
        Form-tagged, variable-arity coefficients. For ``LINEAR_CUBIC``:
        ``{"a": <amu/ps>, "b": <amu*ps/A^2>}``.
    extraction_mass_model : str
        How mass was treated during extraction: ``"constant"`` or
        ``"time_resolved"``. The §6.5 guard (Slice 3) reads this; Slice 1
        only carries it.
    extraction_mass_amu : float
        The constant effective mass [amu] the coefficients were extracted
        under (or the reference value for a ``time_resolved`` ``m(t)``).
        **Provenance only** -- this module is mass-agnostic and never uses it
        in any force or gamma evaluation.

    Raises
    ------
    ValueError
        On unknown form, missing required coefficients, invalid mass model,
        or non-positive extraction mass.
    """

    form: str
    coefficients: Mapping[str, float]
    extraction_mass_model: str
    extraction_mass_amu: float

    def __post_init__(self) -> None:
        if self.form not in _REQUIRED_COEFF_KEYS:
            raise ValueError(
                f"unknown drag form {self.form!r}; expected one of "
                f"{sorted(_REQUIRED_COEFF_KEYS)}"
            )
        required = _REQUIRED_COEFF_KEYS[self.form]
        missing = [k for k in required if k not in self.coefficients]
        if missing:
            raise ValueError(
                f"form {self.form!r} requires coefficients {required}; "
                f"missing {missing}"
            )
        if self.extraction_mass_model not in _VALID_MASS_MODELS:
            raise ValueError(
                f"extraction_mass_model must be one of {_VALID_MASS_MODELS}, "
                f"got {self.extraction_mass_model!r}"
            )
        if not (self.extraction_mass_amu > 0):
            raise ValueError(
                f"extraction_mass_amu must be positive, got "
                f"{self.extraction_mass_amu!r}"
            )


def _raise_unrealised_form(form: str) -> None:
    """Raise for any form not realised in Slice 1 (explicit, never silent)."""
    if form in (LINEAR_QUADRATIC, THRESHOLD):
        raise NotImplementedError(
            f"drag form {form!r} is reserved but not yet extracted (no fit "
            f"pass; see DRAG_PORT_DESIGN_DECISIONS.md §3.7). Only "
            f"{LINEAR_CUBIC!r} is realised in Slice 1."
        )
    if form == POWER_LAW:
        raise NotImplementedError(
            f"drag form {POWER_LAW!r} is deferred: it needs the §3.8 "
            f"low-velocity floor and divergent-gamma handling (out of Tier-0 "
            f"scope). Only {LINEAR_CUBIC!r} is realised in Slice 1."
        )
    # DragCoefficients.__post_init__ already rejects unknown forms; defensive.
    raise ValueError(f"unknown drag form {form!r}")


def spatial_gate(depth, steepness: float) -> np.ndarray:
    """Erf-complement spatial gate ``g(depth)`` in ``[0, 1]``.

    ``g(depth) = 0.5 * (1 - erf(depth / steepness))``: 1 deep inside
    (``depth << 0``), 0.5 at the nominal surface (``depth = 0``), 0 outside
    (``depth >> 0``). Smooth and ``C^1`` -- the continuity the discarded sharp
    boolean gate (G1) lacked. Reuses the same erf/steepness machinery as
    :func:`i2_helium_md.physics.potentials.droplet_potential`, *complemented* so
    that drag turns **off** outside the droplet (G4-collapsing-to-G2, §5.5).

    Parameters
    ----------
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom. Negative inside, positive outside.
    steepness : float
        Width of the erf transition in Angstrom (> 0); e.g.
        ``cfg.potential_steepness = 14.2``.

    Returns
    -------
    np.ndarray
        Dimensionless gate factor in ``[0, 1]``, broadcast to the shape of
        ``depth``.

    Raises
    ------
    ValueError
        If ``steepness`` is not positive.
    """
    if not (steepness > 0):
        raise ValueError(f"steepness must be positive, got {steepness!r}")
    depth = np.asarray(depth, dtype=float)
    return 0.5 * (1.0 - erf(depth / steepness))


def drag_force(v, depth, coeffs: DragCoefficients, steepness: float) -> np.ndarray:
    """Gated drag-force magnitude ``F_drag(v, depth)`` [amu*A/ps^2].

    For ``linear_cubic``::

        F_drag(v, depth) = g(depth) * (a*v + b*v**3)

    Returns the **positive magnitude** form (drag opposes motion); the consumer
    applies ``-F_drag`` along ``v_hat`` at the integrator (Slice 2).

    Parameters
    ----------
    v : array_like
        Speed in A/ps. Same shape as (or broadcastable with) ``depth``.
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom (negative inside).
    coeffs : DragCoefficients
        Form-tagged coefficient bundle. Only ``LINEAR_CUBIC`` is realised.
    steepness : float
        Gate width in Angstrom (> 0).

    Returns
    -------
    np.ndarray
        Drag-force magnitude in amu*A/ps^2.

    Raises
    ------
    NotImplementedError
        If ``coeffs.form`` is a reserved/deferred (non-``linear_cubic``) form.
    """
    v = np.asarray(v, dtype=float)
    g = spatial_gate(depth, steepness)
    if coeffs.form == LINEAR_CUBIC:
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        return g * (a * v + b * v**3)
    _raise_unrealised_form(coeffs.form)


def drag_gamma(v, depth, coeffs: DragCoefficients, steepness: float) -> np.ndarray:
    """Gated friction force-coefficient ``gamma(v, depth)`` [amu/ps] (closed form).

    For ``linear_cubic``::

        gamma(v, depth) = g(depth) * (a + b*v**2)

    Exposed via the **closed form**, *not* ``|F_drag|/v``: analytically equal,
    but the division is singular at ``v -> 0`` where ``linear_cubic`` is finite
    (``gamma -> g*a``). Later consumed by *both* the O-step rate ``gamma/m`` and
    the FDT noise amplitude ``sqrt(2*gamma*kB*Teff)``; carries the **same**
    ``g(depth)`` as :func:`drag_force` (hard FDT coupling, §5.2), so the noise
    the future O-step reads is gated consistently with the drag.

    Parameters
    ----------
    v : array_like
        Speed in A/ps. Same shape as (or broadcastable with) ``depth``.
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom (negative inside).
    coeffs : DragCoefficients
        Form-tagged coefficient bundle. Only ``LINEAR_CUBIC`` is realised.
    steepness : float
        Gate width in Angstrom (> 0).

    Returns
    -------
    np.ndarray
        Friction force-coefficient in amu/ps. Finite at ``v = 0``
        (``-> g*a``) for ``linear_cubic``.

    Raises
    ------
    NotImplementedError
        If ``coeffs.form`` is a reserved/deferred (non-``linear_cubic``) form.
    """
    v = np.asarray(v, dtype=float)
    g = spatial_gate(depth, steepness)
    if coeffs.form == LINEAR_CUBIC:
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        return g * (a + b * v**2)
    _raise_unrealised_form(coeffs.form)
