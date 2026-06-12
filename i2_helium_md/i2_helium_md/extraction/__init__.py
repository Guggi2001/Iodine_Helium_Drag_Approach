"""Drag-law extraction (Method B, trajectory matching).

Extraction is its own concern: it *runs* simulations to fit model parameters,
so it belongs neither in ``physics/`` (pure science, no I/O or orchestration)
nor in ``postprocess/`` (analysis of finished runs) nor ``simulation/``
(orchestration of a single run). See
``METHOD_B_trajectory_matching_extraction.md``.
"""

from .trajectory_matching import (
    CaseSetup,
    ObjectiveResult,
    TrajectoryMatchingFit,
    build_case_setup,
    evaluate_objective,
    fit_trajectory_matching,
    sensitivity_halfwidths,
    write_fit_parameters,
)

__all__ = [
    "CaseSetup",
    "ObjectiveResult",
    "TrajectoryMatchingFit",
    "build_case_setup",
    "evaluate_objective",
    "fit_trajectory_matching",
    "sensitivity_halfwidths",
    "write_fit_parameters",
]
