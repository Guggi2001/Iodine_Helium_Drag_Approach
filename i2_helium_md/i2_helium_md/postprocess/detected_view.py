"""Checkpoint-shaped view of the detected ensemble.

``DetectedEnsembleView`` duck-types exactly the ``IonCheckpoint`` surface
the shared mass-gated diagnostics read (``select_final_mass_gate`` and
every consumer built on it: final-velocity histograms, ihe_ked fragment
moments, paper-v2/v4 VMI and covariance selections), so the legacy
plotting recipes run **verbatim** on the detection-stage ensemble — rule 1
with zero recipe duplication (`TIER2_DETECTION_SUMMARY_SPEC.md` §3, as
amended 2026-07-22).

Forcing rules (the spec's detected-ensemble conventions):

* ``frozen`` / ``time_exhausted`` — ``mass_detected_kg`` pass-through;
  the discrete-shed mechanism writes ``m(n_detected)`` bit-exactly, so
  the existing exact-mass gate *is* ``n_detected`` selection;
* ``suppressed`` — mass forced to the bare gate ``m(0)``: on disk these
  ions ride at their handover ``n``, but the frozen Tier-2 convention
  scores them at n = 0 (bare I⁺ at the detector, the RQ3 channel); their
  stored detected velocities are kept;
* ``droplet_retained`` — mass forced to NaN: the rows hold the verbatim
  in-droplet handover state and must match **no** gate. Exclusion is
  thereby per-fragment, and pair diagnostics drop mixed pairs through
  their existing pair-AND.

``b_ion_outside`` is all-True: on the detected ensemble the outside
criterion is meaningless (everything that reaches the detector is
outside) and exclusion is carried entirely by the mass sentinel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..physics.constants import U as U_KG
from ..physics.shell_schedule import complex_mass_amu
from ..simulation.detection_stage import DetectionResult

__all__ = ["DetectedEnsembleView", "detected_ensemble_view"]


@dataclass(frozen=True)
class DetectedEnsembleView:
    """Checkpoint-shaped detected ensemble (2N rows, checkpoint layout).

    Attributes
    ----------
    num_molecules : int
        ``N``; per-ion arrays have shape ``(2N,)`` in the checkpoint
        layout (fragment-1 block then fragment-2 block), inherited from
        ``DetectionResult``.
    mass_final_kg : np.ndarray, shape (2N,)
        Gate-defining mass [kg] after the forcing rules (module
        docstring): ``m(n_detected)`` for detected in-band ions,
        ``m(0)`` for suppressed, NaN for droplet-retained.
    velocities_final_x, velocities_final_y, velocities_final_z :
        np.ndarray, shape (2N,)
        Terminal detected velocity components [A/ps], pass-through from
        ``v{x,y,z}_detected``.
    b_ion_outside : np.ndarray, shape (N,), bool
        All True (see module docstring).
    """

    num_molecules: int
    mass_final_kg: np.ndarray
    velocities_final_x: np.ndarray
    velocities_final_y: np.ndarray
    velocities_final_z: np.ndarray
    b_ion_outside: np.ndarray


def detected_ensemble_view(det: DetectionResult) -> DetectedEnsembleView:
    """Build the checkpoint-shaped view of one detection result.

    Parameters
    ----------
    det
        A loaded ``detection.npz`` payload (``load_detection_result``).

    Returns
    -------
    DetectedEnsembleView
        See the class docstring for the forcing rules applied.

    Raises
    ------
    ValueError
        If any per-ion array does not have shape ``(2N,)`` for
        ``N = det.num_molecules``.
    """
    n_ions = 2 * int(det.num_molecules)
    reason = np.asarray(det.state_reason)
    arrays = {
        "state_reason": reason,
        "n_detected": np.asarray(det.n_detected, dtype=float),
        "mass_detected_kg": np.asarray(det.mass_detected_kg, dtype=float),
        "vx_detected": np.asarray(det.vx_detected, dtype=float),
        "vy_detected": np.asarray(det.vy_detected, dtype=float),
        "vz_detected": np.asarray(det.vz_detected, dtype=float),
    }
    for name, arr in arrays.items():
        if arr.shape != (n_ions,):
            raise ValueError(
                f"detection array {name!r} has shape {arr.shape}, expected "
                f"({n_ions},) for num_molecules={det.num_molecules}."
            )

    mass_kg = arrays["mass_detected_kg"].copy()
    mass_kg[reason == "suppressed"] = complex_mass_amu(0) * U_KG
    mass_kg[reason == "droplet_retained"] = np.nan

    return DetectedEnsembleView(
        num_molecules=int(det.num_molecules),
        mass_final_kg=mass_kg,
        velocities_final_x=arrays["vx_detected"],
        velocities_final_y=arrays["vy_detected"],
        velocities_final_z=arrays["vz_detected"],
        b_ion_outside=np.ones(int(det.num_molecules), dtype=bool),
    )
