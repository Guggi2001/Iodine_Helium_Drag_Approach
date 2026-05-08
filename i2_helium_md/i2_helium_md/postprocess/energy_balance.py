"""Recipe functions for the legacy MATLAB post-processing figures.

Reproduces the in-script "live debug" panels of
``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` (energy balance, line 965)
and ``vmi_sim_3d_ion_propa.m`` (energy balance, line 898; temperature
diagnostic, line 883 -- the third column there is the lab-frame
scattering angle, not the COM angle), plus the simulation-side panels
of ``post_process_single_pulse_paper_v3.m`` (azimuthal phi histogram,
line 314; final ion mass spectrum, line 397) that don't require
polar VMI image data.

All functions are pure: they take a ``NeutralCheckpoint`` or
``IonCheckpoint`` and return numpy arrays. No plotting, no I/O.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..physics.constants import U
from ..simulation.checkpoint import IonCheckpoint, NeutralCheckpoint


# ===========================================================================
# Energy totals
# ===========================================================================
@dataclass
class EnergyTotals:
    """Summed-over-atoms energy traces for one stage of the simulation.

    All traces share the same ``time_ps`` axis. ``E_system`` is the
    sum of the other components (matches MATLAB ``E_system =
    E_kin + E_pot + E_dissip [+ E_mass_attach_defect]``); it should be
    ~flat over time within Verlet drift.
    """
    time_ps: np.ndarray
    E_kin_eV: np.ndarray
    E_pot_eV: np.ndarray
    E_dissip_eV: np.ndarray
    E_system_eV: np.ndarray
    # Only populated for the ion stage; ``None`` for neutral.
    E_mass_attach_defect_eV: np.ndarray | None = None


def neutral_energy_totals(ckpt: NeutralCheckpoint) -> EnergyTotals:
    """Sum-over-atoms energy traces for a neutral run.

    Mirrors ``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` line 965 where
    each component is plotted as ``sum(E_*, 1)``. ``E_system`` is
    ``E_kin + E_pot + E_dissip``.
    """
    e_kin = np.sum(ckpt.E_kin_eV, axis=0)
    e_pot = np.sum(ckpt.E_pot_eV, axis=0)
    e_dis = np.sum(ckpt.E_dissip_eV, axis=0)
    return EnergyTotals(
        time_ps=ckpt.time_ps,
        E_kin_eV=e_kin,
        E_pot_eV=e_pot,
        E_dissip_eV=e_dis,
        E_system_eV=e_kin + e_pot + e_dis,
        E_mass_attach_defect_eV=None,
    )


def ion_energy_totals(ckpt: IonCheckpoint) -> EnergyTotals:
    """Per-molecule energy traces for an ion run.

    Mirrors ``vmi_sim_3d_ion_propa.m`` line 898 where each component
    is plotted as ``sum(E_*, 1) / num_molecules``. ``E_system``
    includes ``E_mass_attach_defect`` so the total is conserved up
    to Verlet drift.
    """
    n = float(ckpt.num_molecules)
    e_kin = np.sum(ckpt.E_kin_eV, axis=0) / n
    e_pot = np.sum(ckpt.E_pot_eV, axis=0) / n
    e_dis = np.sum(ckpt.E_dissip_eV, axis=0) / n
    e_def = np.sum(ckpt.E_mass_attach_defect_eV, axis=0) / n
    return EnergyTotals(
        time_ps=ckpt.time_ps,
        E_kin_eV=e_kin,
        E_pot_eV=e_pot,
        E_dissip_eV=e_dis,
        E_system_eV=e_kin + e_pot + e_dis + e_def,
        E_mass_attach_defect_eV=e_def,
    )


# ===========================================================================
# Phi histogram for the paper figure
# ===========================================================================
@dataclass
class PhiHistogram:
    """Azimuthal final-velocity histogram for one mass selection.

    ``counts`` is the raw histogram; ``density`` is normalised to
    integrate to 1. Bin centers are the midpoints of ``bin_edges``.
    """
    bin_centers_rad: np.ndarray
    bin_edges_rad: np.ndarray
    counts: np.ndarray
    density: np.ndarray
    num_atoms_used: int


def phi_histogram(
    ckpt: IonCheckpoint,
    *,
    bin_width_rad: float = 0.05,
    mass_amu: float | None = None,
    mass_tolerance_amu: float = 0.5,
) -> PhiHistogram:
    """Histogram ``arctan2(vy_final, vx_final) + pi`` over [0, 2*pi).

    Mirrors ``post_process_single_pulse_paper_v3.m`` line 314 which
    plots ``phi_sim = atan2(vy_total, vx_total) + pi`` for
    mass-selected ions.

    Parameters
    ----------
    ckpt : IonCheckpoint
    bin_width_rad : float, optional
        Bin width in radians. Default 0.05 rad ~ 2.86 deg.
    mass_amu : float, optional
        If given, restrict to atoms whose ``mass_final_kg`` is within
        ``mass_tolerance_amu`` of this value. ``None`` uses all atoms.
    mass_tolerance_amu : float, optional
        Half-width of the mass selection window in amu.
    """
    mass_final_amu = ckpt.mass_final_kg / U
    if mass_amu is None:
        mask = np.ones(mass_final_amu.shape, dtype=bool)
    else:
        mask = np.abs(mass_final_amu - mass_amu) <= mass_tolerance_amu

    vx = ckpt.velocities_final_x[mask]
    vy = ckpt.velocities_final_y[mask]
    phi = np.mod(np.arctan2(vy, vx) + np.pi, 2.0 * np.pi)

    edges = np.arange(0.0, 2.0 * np.pi + bin_width_rad, bin_width_rad)
    counts, edges = np.histogram(phi, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = counts.sum()
    if total > 0:
        density = counts / (total * bin_width_rad)
    else:
        density = np.zeros_like(counts, dtype=float)

    return PhiHistogram(
        bin_centers_rad=centers,
        bin_edges_rad=edges,
        counts=counts,
        density=density,
        num_atoms_used=int(mask.sum()),
    )


# ===========================================================================
# Mass spectrum for the paper figure
# ===========================================================================
@dataclass
class MassSpectrum:
    """Final ion mass histogram in 1-amu bins."""
    bin_centers_amu: np.ndarray
    bin_edges_amu: np.ndarray
    counts: np.ndarray


def mass_spectrum(
    ckpt: IonCheckpoint,
    *,
    bin_width_amu: float = 1.0,
) -> MassSpectrum:
    """Histogram ``mass_final_kg / U`` with bins centred on integer amu.

    Mirrors ``post_process_single_pulse_paper_v3.m`` line 397 which
    calls ``histogram(data_ion.mass_i(:,end)/u)``. Bin edges are at
    half-integer amu so that 127.0 lands in the bin centred on 127
    rather than on a bin edge.
    """
    masses_amu = ckpt.mass_final_kg / U
    if masses_amu.size == 0:
        return MassSpectrum(
            bin_centers_amu=np.array([], dtype=float),
            bin_edges_amu=np.array([], dtype=float),
            counts=np.array([], dtype=int),
        )

    half = 0.5 * bin_width_amu
    center_lo = float(np.floor(masses_amu.min())) - bin_width_amu
    center_hi = float(np.ceil(masses_amu.max())) + bin_width_amu
    edges = np.arange(
        center_lo - half, center_hi + half + bin_width_amu, bin_width_amu,
    )
    counts, edges = np.histogram(masses_amu, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return MassSpectrum(
        bin_centers_amu=centers,
        bin_edges_amu=edges,
        counts=counts,
    )
