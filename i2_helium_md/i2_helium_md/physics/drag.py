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

Governing equations (realised forms; METHOD_B §10.3 dimensional analysis)
---------------------------------------------------------------------------
With ``depth = r_atom - r_droplet`` (negative inside the droplet, positive
outside -- the same convention as :func:`i2_helium_md.physics.potentials.droplet_potential`)::

    g(depth) = 0.5 * (1 - erf(depth / steepness))    # dimensionless in [0, 1]

    linear_cubic     {a [amu/ps], b [amu*ps/A^2]}:
        F_drag = g * (a*v + b*v**3)                  # amu*A/ps^2
        gamma  = g * (a + b*v**2)                    # amu/ps
    linear_quadratic {a [amu/ps], c [amu/A]}:
        F_drag = g * (a*v + c*v**2)                  # [c*v**2] = amu*A/ps^2 OK
        gamma  = g * (a + c*v)                       # [c*v]    = amu/ps     OK
    power_law        {C [amu*A^(1-n)*ps^(n-2)], n [dimensionless]}:
        F_drag = g * C * v**n                        # amu*A/ps^2 by [C]
        gamma  = g * C * v**(n-1)                    # amu/ps     by [C]

All unit checks balance to a force / force coefficient; see
``SLICE1_GOALS_gated_drag_module.md`` §3 and METHOD_B §10.3. The drag opposes
motion; these functions return the **positive magnitude** form and the consumer
applies ``-F_drag`` along ``v_hat`` at the integrator (Slice 2).

Nesting identities (the §10.3 cross-check obligations, exact by construction):
``power_law(n=2, C=c) == linear_quadratic(a=0, c)`` (pure-quadratic) and
``power_law(n=3, C=b) == linear_cubic(a=0, b)`` (pure-cubic).

``gamma`` is exposed via its **closed form** per form, *never* via
``|F_drag|/v``: the two are analytically equal, but the division manufactures a
``0/0`` singularity at ``v -> 0`` that the closed forms do not have
(``gamma -> g*a`` for ``linear_cubic``/``linear_quadratic``; ``gamma -> 0`` for
``power_law`` with ``n > 1`` and ``-> g*C`` at ``n = 1``). The §3.3 config
guard enforces ``n >= 1``; an ``n < 1`` law would diverge at rest and would
activate the §3.8 ``drag_low_v_floor`` obligation (the floor stays inert).

``threshold`` remains reserved behind the same dispatch and raises
:class:`NotImplementedError`; see ``SLICE1_GOALS_gated_drag_module.md`` §5.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from ._gates import _erf_complement


# Drag-form tags. LINEAR_CUBIC (Slice 1), LINEAR_QUADRATIC and POWER_LAW
# (METHOD_B §10 form phase) are realised; THRESHOLD stays reserved behind the
# dispatch (see :func:`_raise_unrealised_form`).
LINEAR_CUBIC = "linear_cubic"
LINEAR_QUADRATIC = "linear_quadratic"
THRESHOLD = "threshold"
POWER_LAW = "power_law"

# Required coefficient keys per form (variable arity by form, §3.8).
# POWER_LAW's amplitude key is "C" (METHOD_B §10.3/§10.5 raw-{C, n} stamp);
# the historical Method-A export 18A/power/fit_parameters.json keeps its
# legacy "gamma" key -- it is frozen evidence read once for the locked C0
# anchor constant, never loaded through load_drag_coefficients.
_REQUIRED_COEFF_KEYS: dict[str, tuple[str, ...]] = {
    LINEAR_CUBIC: ("a", "b"),            # amu/ps, amu*ps/A^2
    LINEAR_QUADRATIC: ("a", "c"),        # amu/ps, amu/A
    THRESHOLD: ("F_sat", "v0"),          # amu*A/ps^2, A/ps
    POWER_LAW: ("C", "n"),               # amu*A^(1-n)*ps^(n-2), dimensionless
}

# The forms with realised force/gamma branches below -- the single source for
# the loader's form acceptance and the Slice-4 driver's scope guard. THRESHOLD
# is deliberately absent (reserved; out of the METHOD_B §10 form-phase scope).
REALIZED_FORMS: tuple[str, ...] = (LINEAR_CUBIC, LINEAR_QUADRATIC, POWER_LAW)

_VALID_MASS_MODELS = ("constant", "time_resolved")

# How the coefficients were extracted (provenance, METHOD_B doc §6):
# "force_balance" = Method A (direct F_drag-vs-v regression, the original
# pipeline); "trajectory_matching" = Method B (joint {a, b, E_bind} fit by
# forward-integrated trajectory RMSE).
_VALID_EXTRACTION_METHODS = ("force_balance", "trajectory_matching")


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
        POWER_LAW}``. ``THRESHOLD`` is reserved (not realised).
    coefficients : Mapping[str, float]
        Form-tagged, variable-arity coefficients. ``LINEAR_CUBIC``:
        ``{"a": <amu/ps>, "b": <amu*ps/A^2>}``; ``LINEAR_QUADRATIC``:
        ``{"a": <amu/ps>, "c": <amu/A>}``; ``POWER_LAW``:
        ``{"C": <amu*A^(1-n)*ps^(n-2)>, "n": <dimensionless>}``.
    extraction_mass_model : str
        How mass was treated during extraction: ``"constant"`` or
        ``"time_resolved"``. The §6.5 guard (Slice 3) reads this; Slice 1
        only carries it.
    extraction_mass_amu : float
        The constant effective mass [amu] the coefficients were extracted
        under (or the reference value for a ``time_resolved`` ``m(t)``).
        **Provenance only** -- this module is mass-agnostic and never uses it
        in any force or gamma evaluation.
    extraction_method : str
        Provenance: ``"force_balance"`` (Method A, the original direct
        ``F_drag``-vs-``v`` regression -- the default, so legacy bundles load
        unchanged) or ``"trajectory_matching"`` (Method B, the joint
        ``{a, b, E_bind}`` forward-integration fit). Carried for the §6.5.1
        guard and artifact stamping; never used in any force evaluation.
    effective_binding_energy_I_ion_eV : float or None
        The effective droplet binding depth [eV] the coefficients were
        **jointly calibrated with** (Method B, §6.5.1 coupled pair), or
        ``None`` for a legacy bundle whose drag<->binding pairing was never
        jointly validated. **Provenance only** -- consumed by the config-load
        guard, never by this module.

    Raises
    ------
    ValueError
        On unknown form, missing required coefficients, invalid mass model,
        invalid extraction method, non-positive extraction mass, or
        non-positive effective binding energy.
    """

    form: str
    coefficients: Mapping[str, float]
    extraction_mass_model: str
    extraction_mass_amu: float
    extraction_method: str = "force_balance"
    effective_binding_energy_I_ion_eV: float | None = None

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
        if self.extraction_method not in _VALID_EXTRACTION_METHODS:
            raise ValueError(
                f"extraction_method must be one of "
                f"{_VALID_EXTRACTION_METHODS}, got {self.extraction_method!r}"
            )
        if self.effective_binding_energy_I_ion_eV is not None and not (
            self.effective_binding_energy_I_ion_eV > 0
        ):
            raise ValueError(
                f"effective_binding_energy_I_ion_eV must be positive or None, "
                f"got {self.effective_binding_energy_I_ion_eV!r}"
            )


def _raise_unrealised_form(form: str) -> None:
    """Raise for any form not realised (explicit, never silent)."""
    if form == THRESHOLD:
        raise NotImplementedError(
            f"drag form {THRESHOLD!r} is reserved but not realised (explicitly "
            f"out of the METHOD_B §10 form-phase scope; see "
            f"DRAG_PORT_DESIGN_DECISIONS.md §3.7). Realised forms: "
            f"{LINEAR_CUBIC!r}, {LINEAR_QUADRATIC!r}, {POWER_LAW!r}."
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

    Notes
    -----
    Routed through the single-source :func:`i2_helium_md.physics._gates._erf_complement`
    (Slice rho) so this drag gate and the Tier-2 ``rho_He/rho_bulk`` density gate are
    one formula in one place (CLAUDE.md rule 1). No-behaviour-change refactor:
    arithmetic and the ``steepness > 0`` guard are identical to the former in-lined
    form, locked by the parity regression tests.
    """
    return _erf_complement(depth, steepness)


def drag_force(v, depth, coeffs: DragCoefficients, steepness: float) -> np.ndarray:
    """Gated drag-force magnitude ``F_drag(v, depth)`` [amu*A/ps^2].

    Per realised form::

        linear_cubic:     F_drag = g(depth) * (a*v + b*v**3)
        linear_quadratic: F_drag = g(depth) * (a*v + c*v**2)
        power_law:        F_drag = g(depth) * C * v**n

    Returns the **positive magnitude** form (drag opposes motion); the consumer
    applies ``-F_drag`` along ``v_hat`` at the integrator (Slice 2).

    Parameters
    ----------
    v : array_like
        Speed in A/ps. Same shape as (or broadcastable with) ``depth``.
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom (negative inside).
    coeffs : DragCoefficients
        Form-tagged coefficient bundle (``threshold`` is reserved).
    steepness : float
        Gate width in Angstrom (> 0).

    Returns
    -------
    np.ndarray
        Drag-force magnitude in amu*A/ps^2.

    Raises
    ------
    NotImplementedError
        If ``coeffs.form`` is the reserved ``threshold`` form.
    """
    v = np.asarray(v, dtype=float)
    g = spatial_gate(depth, steepness)
    if coeffs.form == LINEAR_CUBIC:
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        return g * (a * v + b * v**3)
    if coeffs.form == LINEAR_QUADRATIC:
        a = float(coeffs.coefficients["a"])
        c = float(coeffs.coefficients["c"])
        return g * (a * v + c * v**2)
    if coeffs.form == POWER_LAW:
        C = float(coeffs.coefficients["C"])
        n = float(coeffs.coefficients["n"])
        return g * C * v**n
    _raise_unrealised_form(coeffs.form)


def drag_gamma(v, depth, coeffs: DragCoefficients, steepness: float) -> np.ndarray:
    """Gated friction force-coefficient ``gamma(v, depth)`` [amu/ps] (closed form).

    Per realised form::

        linear_cubic:     gamma = g(depth) * (a + b*v**2)
        linear_quadratic: gamma = g(depth) * (a + c*v)
        power_law:        gamma = g(depth) * C * v**(n-1)

    Exposed via the **closed form**, *not* ``|F_drag|/v``: analytically equal,
    but the division is singular at ``v -> 0`` where the closed forms are
    regular (``gamma -> g*a`` for ``linear_cubic``/``linear_quadratic``;
    ``gamma -> 0`` for ``power_law`` with ``n > 1``, ``-> g*C`` at ``n = 1`` --
    the §3.3 guard enforces ``n >= 1``, where ``v**(n-1)`` is finite at rest).
    Later consumed by *both* the O-step rate ``gamma/m`` and the FDT noise
    amplitude ``sqrt(2*gamma*kB*Teff)``; carries the **same** ``g(depth)`` as
    :func:`drag_force` (hard FDT coupling, §5.2), so the noise the future
    O-step reads is gated consistently with the drag.

    Parameters
    ----------
    v : array_like
        Speed in A/ps. Same shape as (or broadcastable with) ``depth``.
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom (negative inside).
    coeffs : DragCoefficients
        Form-tagged coefficient bundle (``threshold`` is reserved).
    steepness : float
        Gate width in Angstrom (> 0).

    Returns
    -------
    np.ndarray
        Friction force-coefficient in amu/ps. Finite at ``v = 0`` for every
        realised form within its guard-validated coefficient domain.

    Raises
    ------
    NotImplementedError
        If ``coeffs.form`` is the reserved ``threshold`` form.
    """
    v = np.asarray(v, dtype=float)
    g = spatial_gate(depth, steepness)
    if coeffs.form == LINEAR_CUBIC:
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        return g * (a + b * v**2)
    if coeffs.form == LINEAR_QUADRATIC:
        a = float(coeffs.coefficients["a"])
        c = float(coeffs.coefficients["c"])
        return g * (a + c * v)
    if coeffs.form == POWER_LAW:
        C = float(coeffs.coefficients["C"])
        n = float(coeffs.coefficients["n"])
        # n >= 1 (guard-enforced) keeps v**(n-1) finite at v = 0; numpy's
        # 0.0**0.0 == 1.0 realizes the n = 1 limit gamma -> g*C exactly.
        return g * C * v ** (n - 1.0)
    _raise_unrealised_form(coeffs.form)
