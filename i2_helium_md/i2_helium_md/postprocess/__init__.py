"""Post-processing utilities for MD ion-stage trajectories.

Step 13: HeDFT/TDDFT reference loading and numerical comparison against
``IonCheckpoint`` output. Mirrors the trajectory comparison block in
``legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
simulation_image_only_trajectories.m``.
"""

from .boltzmann_overlay import BoltzmannCurve, boltzmann_population
from .compare_trajectories import (
    NeutralComparison,
    TrajectoryComparison,
    compare_distance,
    compare_neutral_to_hedft,
    compare_velocity_magnitude,
)
from .energy_balance import (
    EnergyTotals,
    MassSpectrum,
    PhiHistogram,
    ion_energy_totals,
    mass_spectrum,
    neutral_energy_totals,
    phi_histogram,
)
from .hedft_loader import HedftTrajectory, load_hedft_trajectory
from .pair_correlation import (
    CovarianceMatrix,
    DistanceHistogram,
    angular_pair_covariance,
    interparticle_distance_histogram,
)
from .paper_v3 import (
    PaperV3PhiCurve,
    PaperV3PhiReference,
    PaperV3RadialReference,
    PaperV3VelocityCurve,
    load_paper_v3_phi_reference,
    load_paper_v3_radial_reference,
    matlab_max_normalise,
    paper_v3_phi_curve,
    paper_v3_velocity_curve,
)
from .paper_v2 import (
    PaperV2PhiCurve,
    PaperV2PhiReference,
    PaperV2RadialReference,
    PaperV2VelocityCurve,
    PaperV2VelocityMap,
    PaperV2VMIImageReference,
    PaperV2VMIPolarImageReference,
    load_paper_v2_radial_reference,
    load_paper_v2_phi_reference,
    load_paper_v2_he2_radial_references,
    load_paper_v2_radial_references,
    load_paper_v2_vmi_image_reference,
    load_paper_v2_vmi_polar_image_reference,
    paper_v2_phi_curve,
    paper_v2_velocity_curve,
    paper_v2_velocity_map,
)
from .paper_v4 import (
    PaperV4AngularCovariance,
    PaperV4RadialReference,
    PaperV4VelocityCurve,
    load_paper_v4_radial_reference,
    load_paper_v4_radial_references,
    paper_v4_angular_pair_covariance,
    paper_v4_velocity_curve,
)
from .polar_velocity import (
    AnisotropyFit,
    BetaCurve,
    PolarHistogram,
    anisotropy_fit,
    beta_of_velocity,
    polar_velocity_histogram,
)
from .time_resolved import RadialEvolution, radial_distribution_evolution
from .velocity_2d import Velocity2DHistogram, velocity_density_2d
from .velocity_distribution import (
    BimodalGaussianFit,
    FinalVelocityHistogram,
    VmiReference,
    bimodal_gaussian_fit,
    compute_final_velocity_histogram,
    load_vmi_reference,
)

__all__ = [
    "HedftTrajectory",
    "load_hedft_trajectory",
    "TrajectoryComparison",
    "compare_distance",
    "compare_velocity_magnitude",
    "NeutralComparison",
    "compare_neutral_to_hedft",
    "VmiReference",
    "load_vmi_reference",
    "FinalVelocityHistogram",
    "compute_final_velocity_histogram",
    "BimodalGaussianFit",
    "bimodal_gaussian_fit",
    "EnergyTotals",
    "neutral_energy_totals",
    "ion_energy_totals",
    "PhiHistogram",
    "phi_histogram",
    "MassSpectrum",
    "mass_spectrum",
    "PolarHistogram",
    "polar_velocity_histogram",
    "AnisotropyFit",
    "anisotropy_fit",
    "BetaCurve",
    "beta_of_velocity",
    "Velocity2DHistogram",
    "velocity_density_2d",
    "DistanceHistogram",
    "interparticle_distance_histogram",
    "CovarianceMatrix",
    "angular_pair_covariance",
    "PaperV3RadialReference",
    "PaperV3PhiReference",
    "PaperV3VelocityCurve",
    "PaperV3PhiCurve",
    "load_paper_v3_radial_reference",
    "load_paper_v3_phi_reference",
    "paper_v3_velocity_curve",
    "paper_v3_phi_curve",
    "matlab_max_normalise",
    "PaperV2RadialReference",
    "PaperV2VMIImageReference",
    "PaperV2VMIPolarImageReference",
    "PaperV2PhiReference",
    "PaperV2VelocityMap",
    "PaperV2VelocityCurve",
    "PaperV2PhiCurve",
    "load_paper_v2_radial_reference",
    "load_paper_v2_phi_reference",
    "load_paper_v2_he2_radial_references",
    "load_paper_v2_radial_references",
    "load_paper_v2_vmi_image_reference",
    "load_paper_v2_vmi_polar_image_reference",
    "paper_v2_velocity_map",
    "paper_v2_velocity_curve",
    "paper_v2_phi_curve",
    "PaperV4RadialReference",
    "PaperV4VelocityCurve",
    "PaperV4AngularCovariance",
    "load_paper_v4_radial_reference",
    "load_paper_v4_radial_references",
    "paper_v4_velocity_curve",
    "paper_v4_angular_pair_covariance",
    "RadialEvolution",
    "radial_distribution_evolution",
    "BoltzmannCurve",
    "boltzmann_population",
]
