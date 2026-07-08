"""Terminal I+He_n integer-n size-distribution extractor (Tier-2 Phase E, Slice E1).

From a finished generative (``biphasic_energy_gated``) ion run -- or the E2
post-ejection relaxation result -- build the simulated I+He_n size distribution
on integer-n support. This is the direct counterpart of the experimental
abundance reference ``data/reference/integrated_i_he_abundance.csv`` (loaded by
E3), and the object the E4 Wasserstein comparison scores.

The native Poisson spread of the terminal shell count (MASS doc Sec.2.2) is
preserved: the distribution is a histogram over the integer shell count itself,
never a binning of a continuous endpoint.

Three input modes, duck-typed on the terminal shell-count column:

* **sim-end** -- a raw :class:`~i2_helium_md.simulation.checkpoint.IonCheckpoint`
  (any object exposing ``n_shell (2N, T)``). The terminal column ``n_shell[:, -1]``
  is the R5 *upper bound* on the true terminal ``n`` (the cascade continues past
  the 20 ps ion stage; MASS doc Sec.R5). ``source = "sim_end"``.
* **relaxed** -- an E2 ``RelaxationResult`` (any object exposing
  ``terminal_n (2N,)``), the matched-time terminal ``n``. ``source = "relaxed"``.
* **detected** -- a Slice-DS ``DetectionResult`` (``terminal_n`` plus the
  disambiguating ``state_reason``), the detector-arrival terminal ``n`` at the
  Sourced ``t_detect``. ``source = "detected"`` -- the Tier-2 arbitration
  observable; the relaxed read is a convergence diagnostic beside it.

The inference has one blind spot: a relaxation checkpoint **reloaded from
``relaxation.npz``** is itself a v7 ``IonCheckpoint`` and would be inferred
``"sim_end"`` although its terminal column is the matched-time relaxed ``n``.
Callers on that path pass the explicit ``source_tag="relaxed"`` override.

Both iodine fragments of every molecule are counted: the I2 Coulomb explosion
yields two I+ ions, each carrying its own He shell, and each is an independent
I+He_n detection. The extractor therefore histograms all ``2N`` atoms with no
masking.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Langmuir occupancy cap / initial shell count n* (MASS doc Sec.11; the runs
# start at n0 = n* = 21). n = n* is a legal simulated terminal outcome; n > n*
# is forbidden by the cap, so its appearance signals upstream corruption.
# Single source: physics/constants.py (the value the pickup channel enforces
# at run time) -- re-exported here, never redeclared.
from ..physics.constants import N_STAR
from ..physics.shell_schedule import complex_mass_amu

# Legal provenance tags for ShellDistribution.source (the E4/F3 pairing
# discriminator): the R5 sim-end upper bound, the E2 matched-time result, and
# the Slice-DS detector-arrival read (the Tier-2 arbitration observable; the
# relaxed read is demoted to a convergence diagnostic once a detected read
# exists -- TIER2_DETECTION_STAGE_DESIGN.md §3.5).
_SOURCE_TAGS: tuple[str, ...] = ("sim_end", "relaxed", "detected")


@dataclass(frozen=True)
class ShellDistribution:
    """Terminal I+He_n size distribution on integer-n support.

    Attributes
    ----------
    n_values : np.ndarray, shape (Nn,), int
        Contiguous integer support ``0 .. n_max`` (``Nn = n_max + 1``).
    counts : np.ndarray, shape (Nn,), int
        Per-n atom counts; ``counts.sum()`` is the number of ion atoms (2N).
    fraction : np.ndarray, shape (Nn,), float
        ``counts`` normalized by its own sum; sums to 1.
    source : str
        Provenance tag: ``"sim_end"`` (the R5 sim-end upper bound, from a raw
        ``IonCheckpoint``), ``"relaxed"`` (the E2 matched-time terminal ``n``;
        a convergence diagnostic once a detected read exists), or
        ``"detected"`` (the Slice-DS detector-arrival read -- the Tier-2
        arbitration observable). E4/F3 pair scores per run by this tag.
    """

    n_values: np.ndarray
    counts: np.ndarray
    fraction: np.ndarray
    source: str

    @property
    def mass_amu(self) -> np.ndarray:
        """Complex mass ``m(n) = 126.90 + 4.0026 n`` [amu] labelling each rung.

        Convenience for cross-checking a rung against the experimental
        mass-window columns; computed on demand from :attr:`n_values`.
        """
        return complex_mass_amu(self.n_values)

    def moments(self) -> tuple[float, float]:
        """Return ``(mean, spread)`` of the distribution over integer ``n``.

        ``mean = sum n*fraction``; ``spread = sqrt(sum fraction*(n-mean)^2)`` --
        the **population** standard deviation of the size distribution (weights
        already sum to 1). Single source of truth for the terminal-``n`` moments
        reported by both the F3 scoreboard and the staircase-probe report, so the
        two never diverge on the spread convention (CLAUDE.md Principle 1).
        """
        n = self.n_values.astype(float)
        p = self.fraction.astype(float)
        mean = float(np.sum(n * p))
        spread = float(np.sqrt(np.sum(p * (n - mean) ** 2)))
        return mean, spread


def _terminal_n_and_source(source) -> tuple[np.ndarray, str]:
    """Extract the terminal per-atom shell-count vector and the provenance tag.

    Dispatches by attribute presence (not ``isinstance``) so the sim-end mode
    lands with E1 and the relaxed mode is admitted once E2 ships its
    ``RelaxationResult`` -- no import of the not-yet-built type.

    Returns
    -------
    (terminal, source) : (np.ndarray shape (2N,), str)

    Raises
    ------
    ValueError
        If ``n_shell`` is not 2-D, or either input carries no atoms/steps.
    TypeError
        If ``source`` exposes neither ``n_shell`` nor ``terminal_n``.
    """
    n_shell = getattr(source, "n_shell", None)
    if n_shell is not None:
        arr = np.asarray(n_shell, dtype=float)
        if arr.ndim != 2:
            raise ValueError(
                f"n_shell must be 2-D (2N, T); got shape {arr.shape}."
            )
        if arr.shape[0] == 0 or arr.shape[1] == 0:
            raise ValueError(
                f"empty ensemble: n_shell has shape {arr.shape}, nothing to "
                f"histogram."
            )
        return arr[:, -1], "sim_end"

    terminal_n = getattr(source, "terminal_n", None)
    if terminal_n is not None:
        arr = np.asarray(terminal_n, dtype=float).ravel()
        if arr.size == 0:
            raise ValueError(
                "empty ensemble: terminal_n has size 0, nothing to histogram."
            )
        # A DetectionResult also exposes terminal_n; its state_reason array
        # (frozen/suppressed/time_exhausted) is the attribute that separates
        # the detector-arrival read from E2's matched-time relaxed read.
        if getattr(source, "state_reason", None) is not None:
            return arr, "detected"
        return arr, "relaxed"

    raise TypeError(
        "source must expose either 'n_shell' (IonCheckpoint, sim-end mode) or "
        f"'terminal_n' (RelaxationResult relaxed mode / DetectionResult "
        f"detected mode); got {type(source).__name__}."
    )


def compute_terminal_shell_distribution(
    source,
    *,
    n_max: int = N_STAR,
    source_tag: str | None = None,
) -> ShellDistribution:
    """Build the terminal I+He_n integer-n size distribution.

    Parameters
    ----------
    source
        Either an ``IonCheckpoint`` (``n_shell (2N, T)`` -> sim-end upper bound)
        or an E2 ``RelaxationResult`` (``terminal_n (2N,)`` -> matched time).
    n_max
        Upper edge of the integer support (default :data:`N_STAR` = 21). ``n =
        n_max`` is a legal outcome; ``n > n_max`` fails loud.
    source_tag
        Explicit provenance override, one of ``"sim_end"`` / ``"relaxed"`` /
        ``"detected"``. ``None`` (default) infers the tag from the input type
        (a ``DetectionResult``'s ``state_reason`` attribute marks it
        ``"detected"``). **Required for a relaxation checkpoint reloaded from
        ``relaxation.npz``**: that is a bona fide v7 ``IonCheckpoint``, so the
        duck-typed inference would tag its matched-time terminal column
        ``"sim_end"`` and F3 would pair the two scores wrongly -- pass
        ``source_tag="relaxed"`` there.

    Returns
    -------
    ShellDistribution
        Frozen histogram on integer support ``0 .. n_max``.

    Raises
    ------
    ValueError
        Non-finite terminal count; a fractional ``n`` (violates the v7
        int-valued contract -> upstream corruption); an ``n`` outside
        ``[0, n_max]`` (``n > n_max`` violates the Langmuir cap); an empty
        ensemble; or a ``source_tag`` outside the legal tags
        (``"sim_end"`` / ``"relaxed"`` / ``"detected"``).
    TypeError
        If ``source`` is neither input mode.
    """
    terminal, inferred_tag = _terminal_n_and_source(source)
    if source_tag is not None and source_tag not in _SOURCE_TAGS:
        raise ValueError(
            f"source_tag must be one of {_SOURCE_TAGS}; got {source_tag!r}."
        )
    tag = source_tag if source_tag is not None else inferred_tag

    if not np.all(np.isfinite(terminal)):
        raise ValueError(
            "terminal shell counts contain non-finite entries (NaN/inf)."
        )

    # v7 stores n_shell as an int-valued float array. A fractional entry means
    # upstream corruption -- fail loud rather than silently round it away.
    residual = terminal - np.rint(terminal)
    if np.any(residual != 0.0):
        bad = int(np.argmax(np.abs(residual)))
        raise ValueError(
            "n_shell must be integer-valued (v7 contract); found fractional "
            f"entry {terminal[bad]!r} at atom index {bad} -- upstream "
            f"corruption."
        )
    n_int = np.rint(terminal).astype(int)

    if np.any(n_int < 0) or np.any(n_int > n_max):
        lo, hi = int(n_int.min()), int(n_int.max())
        raise ValueError(
            f"terminal shell counts out of range [0, {n_max}]: observed "
            f"[{lo}, {hi}]. n > {n_max} violates the Langmuir cap "
            f"(corruption); n < 0 is unphysical."
        )

    counts = np.bincount(n_int, minlength=n_max + 1).astype(int)
    n_values = np.arange(n_max + 1, dtype=int)
    fraction = counts / counts.sum()

    return ShellDistribution(
        n_values=n_values,
        counts=counts,
        fraction=fraction,
        source=tag,
    )
