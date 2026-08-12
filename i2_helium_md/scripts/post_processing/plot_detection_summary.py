"""Detection-stage summary for one finished Tier-2 confirmation run.

The detector-facing counterpart of ``plot_run_summary.py``
(`TIER2_DETECTION_SUMMARY_SPEC.md`, delivered 2026-07-22): every panel is
rendered from the **detection stage** (``detection.npz``, t_detect), the
asymptotic state of the Tier-2 physics — under the drag + biphasic mass
mechanism the ion keeps evolving long after the MD window (E2 relaxation
+ the µs RRK cascade), so end-of-window comparisons are handover
diagnostics, not detector observables.

Conventions (spec §3): the recipe sections run **verbatim** on a
checkpoint-shaped ``DetectedEnsembleView`` — retained ions match no gate
(excluded per-fragment), suppressed ions score in the bare n = 0 gate
with their stored detected velocities (the RQ3 bare-peak channel, newly
populated), and the cov panels pair the molecule's two **detected**
fragments (the E1 2N convention carried through the detection stage).
VMI/polar/cov panels carry the Tier-3 caveat: ensemble second moments
stay under-dispersed until the noise tier is built.

This document has **no legacy mode**: a run dir without ``detection.npz``
raises ``FileNotFoundError``.

Edit the ``USER SETTINGS`` block below and run the script (e.g. from
PyCharm)::

    python scripts/post_processing/plot_detection_summary.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    detected_ensemble_view,
    load_he_abundance_reference,
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    knob_columns_from_cfg,
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    STATE_REASONS,
    load_detection_result,
)

from summary_sections import (  # noqa: E402
    _SectionSkipped,
    _section_detected_ke_anatomy,
    _section_detected_ked_mean_energy,
    _section_detected_size_distribution,
    _section_ihe_ked_curves,
    _section_mass_resolved,
    _section_paper_cov_angular,
    _section_paper_cov_phi_distribution,
    _section_paper_cov_radial,
    _section_paper_cov_radial_distribution,
    _section_paper_cov_traces,
    _section_paper_v2_polar,
    _section_paper_v2_vmi,
)

# =============================================================================
# Stage banners (spec §3/§4)
# =============================================================================
_DETECTED_NOTE = (
    "detected ensemble at t_detect; n = 0 = suppressed/bare channel "
    "(RQ3); droplet-retained excluded"
)
_TIER3_NOTE = (
    "detected ensemble at t_detect; ensemble second moments "
    "under-dispersed — Tier-3 noise stubbed; droplet-retained excluded"
)

# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Run directory holding cfg.json + ion.npz + relaxation.npz + detection.npz.Ag
RUN_DIR: Path = (
    PROJECT_ROOT / "data" / "runs"
    / "9A_drag_shared_pure_cubic_N1000_tier2atlas_conf270_g4fh415"
)

# Reference paths (same meanings as in plot_run_summary.py). ``None``
# skips the matching comparison panels.
IHE_KED_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked"
)
ABUNDANCE_REF_PATH: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)
PAPER_V2_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "paper_v2"
)
PAPER_COV_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "paper_cov"
)

# Experimental-panel noise clip fraction (see plot_run_summary.py).
EXPERIMENTAL_NOISE_FLOOR: float = 0.20

# Mass channel for the VMI / polar / cov panels: on the detected ensemble
# the m = 131 amu gate is exactly the detected n = 1 channel.
PAPER_V2_MASS_AMU: float = 131.0

# Output directory for the PDF and per-panel PNGs. ``None`` -> <RUN_DIR>/figures.
OUT_DIR: Path | None = None


# =============================================================================
# Cover / metadata section
# =============================================================================
def _section_detection_metadata(run_dir, cfg, det, read) -> plt.Figure:
    """Stage-identity cover page: times, per-reason census, scored
    totals, and the physics knobs from the authoritative cfg.json
    (the F3 convention — knobs are never parsed from the dir tag)."""
    fig, ax = plt.subplots(figsize=(8.5, 6.0), constrained_layout=True)
    ax.axis("off")
    ax.set_title(f"Detection summary -- {Path(run_dir).name}")

    lines: list[str] = [
        f"run_dir = {run_dir}",
        "stage: DETECTION (terminal detector ensemble)",
        f"  t_handover = {det.t_handover_ps:.1f} ps",
        f"  t_detect   = {det.detection_time_ps:.0f} ps "
        f"({det.detection_time_ps / 1e6:.2f} us)",
        f"ions: 2N = {2 * det.num_molecules}",
    ]
    reason = np.asarray(det.state_reason)
    for name in STATE_REASONS:
        count = int(np.count_nonzero(reason == name))
        lines.append(f"  {name}: {count} "
                     f"({100.0 * count / reason.size:.1f} %)")
    lines.append(
        f"scored (retained excluded): {read.num_scored}/{read.num_ions}"
    )
    lines.append(
        "  suppressed (bare) fraction: "
        f"{100.0 * read.suppressed_frac:.1f} %"
    )
    lines.append(f"  n_mean(det) = {read.n_mean:.3f}")

    if cfg is not None:
        lines.append("knobs (cfg.json):")
        for key, value in knob_columns_from_cfg(cfg).items():
            lines.append(f"  {key} = {value}")
    else:
        lines.append("no cfg.json found")

    ax.text(0.0, 1.0, "\n".join(lines), family="monospace", fontsize=9,
            verticalalignment="top")
    return fig


# =============================================================================
# Wiring
# =============================================================================
def _load_inputs(run_dir: Path):
    """Load (det, read, view, cfg) from a confirmation run dir.

    Raises ``FileNotFoundError`` when ``detection.npz`` is absent — this
    document has no legacy mode (spec §2.2)."""
    run_dir = Path(run_dir)
    det_path = run_dir / "detection.npz"
    if not det_path.exists():
        raise FileNotFoundError(
            f"{det_path} not found — plot_detection_summary.py renders the "
            "detection-stage ensemble and has no legacy mode; for a run "
            "without a detection stage use plot_run_summary.py."
        )
    det = load_detection_result(det_path)
    read = read_confirmation_detection(det, label=run_dir.name)
    view = detected_ensemble_view(det)
    cfg_path = run_dir / "cfg.json"
    cfg = (
        json.loads(cfg_path.read_text(encoding="utf-8"))
        if cfg_path.exists() else None
    )
    return det, read, view, cfg


def _sections(*, run_dir, cfg, det, read, view, ked_dir, ked_ref,
              abundance, paper_v2_ref_dir, paper_cov_ref_dir):
    """The §4 roster, in order: (label, builder) pairs."""
    return [
        ("detection_metadata",
         lambda: _section_detection_metadata(run_dir, cfg, det, read)),
        ("detected_size_distribution",
         lambda: _section_detected_size_distribution(read, abundance)),
        ("detected_ihe_ked_mean_energy",
         lambda: _section_detected_ked_mean_energy(read, ked_ref)),
        ("detected_ke_anatomy",
         lambda: _section_detected_ke_anatomy(read, ked_ref)),
        ("detected_ihe_ked_curves_3d",
         lambda: _section_ihe_ked_curves(
             view, ked_dir, ked_ref, "3d", stage_note=_DETECTED_NOTE)),
        ("detected_ihe_ked_curves_2d",
         lambda: _section_ihe_ked_curves(
             view, ked_dir, ked_ref, "2d", stage_note=_DETECTED_NOTE)),
        ("detected_mass_resolved_velocities",
         lambda: _section_mass_resolved(view, stage_note=_DETECTED_NOTE)),
        ("detected_paper_v2_vmi",
         lambda: _section_paper_v2_vmi(
             view, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR,
             mass_amu=PAPER_V2_MASS_AMU, stage_note=_TIER3_NOTE)),
        ("detected_paper_v2_polar",
         lambda: _section_paper_v2_polar(
             view, paper_v2_ref_dir, EXPERIMENTAL_NOISE_FLOOR,
             mass_amu=PAPER_V2_MASS_AMU, stage_note=_TIER3_NOTE)),
        ("detected_paper_cov_radial_distribution",
         lambda: _section_paper_cov_radial_distribution(
             view, paper_cov_ref_dir, paper_v2_ref_dir,
             mass_amu=PAPER_V2_MASS_AMU, stage_note=_TIER3_NOTE)),
        ("detected_paper_cov_phi_distribution",
         lambda: _section_paper_cov_phi_distribution(
             view, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU,
             stage_note=_TIER3_NOTE)),
        ("detected_paper_cov_angular_pair_cov",
         lambda: _section_paper_cov_angular(
             view, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU,
             stage_note=_TIER3_NOTE)),
        ("detected_paper_cov_radial_pair_cov",
         lambda: _section_paper_cov_radial(
             view, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU,
             stage_note=_TIER3_NOTE)),
        ("detected_paper_cov_pair_cov_traces",
         lambda: _section_paper_cov_traces(
             view, paper_cov_ref_dir, mass_amu=PAPER_V2_MASS_AMU,
             stage_note=_TIER3_NOTE)),
    ]


def main(run_dir: Path | None = None) -> int:
    run_dir = Path(RUN_DIR if run_dir is None else run_dir)
    det, read, view, cfg = _load_inputs(run_dir)

    ked_dir = Path(IHE_KED_REFERENCE_DIR) if IHE_KED_REFERENCE_DIR else None
    ked_ref = (
        load_ihe_ked_reference(ked_dir / "IHe_KED_reference.csv")
        if ked_dir else None
    )
    abundance = (
        load_he_abundance_reference(Path(ABUNDANCE_REF_PATH))
        if ABUNDANCE_REF_PATH else None
    )
    paper_v2_ref_dir = (
        Path(PAPER_V2_REFERENCE_DIR) if PAPER_V2_REFERENCE_DIR else None
    )
    paper_cov_ref_dir = (
        Path(PAPER_COV_REFERENCE_DIR) if PAPER_COV_REFERENCE_DIR else None
    )

    out_dir = Path(OUT_DIR) if OUT_DIR else (run_dir / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "detection_summary.pdf"

    print(f"[detection_summary] writing {pdf_path}")
    with PdfPages(pdf_path) as pdf:
        for label, builder in _sections(
            run_dir=run_dir, cfg=cfg, det=det, read=read, view=view,
            ked_dir=ked_dir, ked_ref=ked_ref, abundance=abundance,
            paper_v2_ref_dir=paper_v2_ref_dir,
            paper_cov_ref_dir=paper_cov_ref_dir,
        ):
            try:
                fig = builder()
            except _SectionSkipped as skip:
                print(f"[detection_summary] skip {label}: {skip}")
                continue
            if fig is None:
                continue
            pdf.savefig(fig)
            fig.savefig(out_dir / f"{label}.png", dpi=150)
            plt.close(fig)
            print(f"[detection_summary] wrote {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
