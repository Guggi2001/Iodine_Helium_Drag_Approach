"""Post-processing utilities for MD ion-stage trajectories.

Step 13: HeDFT/TDDFT reference loading and numerical comparison against
``IonCheckpoint`` output. Mirrors the trajectory comparison block in
``legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
simulation_image_only_trajectories.m``.
"""

from .compare_trajectories import (
    TrajectoryComparison,
    compare_distance,
    compare_velocity_magnitude,
)
from .hedft_loader import HedftTrajectory, load_hedft_trajectory
from .velocity_distribution import (
    FinalVelocityHistogram,
    VmiReference,
    compute_final_velocity_histogram,
    load_vmi_reference,
)

__all__ = [
    "HedftTrajectory",
    "load_hedft_trajectory",
    "TrajectoryComparison",
    "compare_distance",
    "compare_velocity_magnitude",
    "VmiReference",
    "load_vmi_reference",
    "FinalVelocityHistogram",
    "compute_final_velocity_histogram",
]
