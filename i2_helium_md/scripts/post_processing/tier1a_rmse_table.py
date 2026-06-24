"""Emit the Tier-1a t* RMSE table for finished run directories.

This is the reporting half of the Tier-1a section-5/9 deliverable. It scores the
fixed null and three ``anchored_discrete`` runs against the 9 A TDDFT reference
over the single smoothed-reference window ``[2.67, ref_end]``. The output is a
diagnostic table only: no ranking, threshold, or scientific verdict is asserted
by code.

Usage::

    python scripts/post_processing/tier1a_rmse_table.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 50
T_STAR_VALUES_PS = [0.5, 5.0, 9.0]
WINDOW_START_PS = 2.67
PLOT_T_STAR_PS = 0.5

SAVE_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier1a_rmse_table.csv"
EXPORT_MEAN_SERIES_DIR = PROJECT_ROOT / "data" / "reference" / "drag" / CASE / "tier1a"

SHOW_FIGURE = True
POSITIONS_FIGURE = False
ENERGY_FIGURE = False
FORCE_FIGURE = False
FORCE_ATOM_INDEX = 1  # atom 2, matching the clean Tier-0 extraction atom


from scripts.tier1a_common import (  # noqa: E402
    tier1a_run_dir_name,
    tier1a_run_tag,
)
from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    SmoothedSpeedReference,
    compare_distance,
    compare_speed_to_reference,
    ion_ledger_closure,
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    _reconstruct_radial_forces,
    ensemble_mean_series,
    export_mean_series,
    plot_energy_analysis,
)


TIER1A_TABLE_COLUMNS = [
    "case",
    "variant",
    "N",
    "run_tag",
    "scenario",
    "t_star_ps",
    "window_start_ps",
    "window_end_ps",
    "n_scored_R",
    "n_scored_v2_smoothed",
    "RMSE_R_raw_A",
    "RMSE_v2_smoothed_Aps",
    "v2_smoothed_mean_ratio",
    "ledger_max_resid_eV",
    "n_sheds",
    "n_shell_start",
    "n_shell_end",
]


@dataclass(frozen=True)
class Tier1aRunRecord:
    """Loaded Tier-1a run plus the metadata needed for tables and figures."""

    label: str
    run_tag: str
    scenario: str
    t_star_ps: float | None
    ion: IonCheckpoint
    cfg: SimConfig | None
    row: dict[str, Any]


def _shell_summary(n_shell: np.ndarray) -> tuple[int, int, int]:
    """Return ``(n_start, n_end, n_sheds)`` from a v6 ``n_shell`` array."""
    if n_shell.ndim != 2 or n_shell.shape[1] == 0:
        raise ValueError(
            f"n_shell must have shape (2N, T) with T > 0, got {n_shell.shape}"
        )
    n_start = int(round(float(n_shell[0, 0])))
    n_end = int(round(float(n_shell[0, -1])))
    return n_start, n_end, max(0, n_start - n_end)


def _run_label(scenario: str, t_star_ps: float | None) -> str:
    if scenario == "fixed":
        return "fixed"
    return f"anchored t*={float(t_star_ps):.1f}"


def _score_ion(
    ion: IonCheckpoint,
    *,
    case: str,
    variant: str,
    n: int,
    run_tag: str,
    scenario: str,
    t_star_ps: float | None,
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window_start_ps: float,
) -> dict[str, Any]:
    """Score one loaded ion checkpoint and return one table row."""
    window_end_ps = float(smoothed.time_ps[-1])
    window = (float(window_start_ps), window_end_ps)

    dist = compare_distance(ion, hedft, window=window)
    v2_smoothed = compare_speed_to_reference(
        ion,
        atom="I2",
        t_ref_ps=smoothed.time_ps,
        ref_speed_Aps=smoothed.speed_Aps,
        window=window,
    )
    closure = ion_ledger_closure(ion)
    n_shell_start, n_shell_end, n_sheds = _shell_summary(ion.n_shell)

    return {
        "case": case,
        "variant": variant,
        "N": int(n),
        "run_tag": run_tag,
        "scenario": scenario,
        "t_star_ps": t_star_ps,
        "window_start_ps": window[0],
        "window_end_ps": window[1],
        "n_scored_R": dist.num_overlap_points,
        "n_scored_v2_smoothed": v2_smoothed.num_overlap_points,
        "RMSE_R_raw_A": dist.rmse,
        "RMSE_v2_smoothed_Aps": v2_smoothed.rmse,
        "v2_smoothed_mean_ratio": v2_smoothed.mean_ratio,
        "ledger_max_resid_eV": closure.max_abs_residual_eV,
        "n_sheds": n_sheds,
        "n_shell_start": n_shell_start,
        "n_shell_end": n_shell_end,
    }


def score_tier1a_run(
    run: RunDirectory | str | Path,
    *,
    case: str,
    variant: str,
    n: int,
    run_tag: str,
    scenario: str,
    t_star_ps: float | None,
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window_start_ps: float = WINDOW_START_PS,
) -> dict[str, Any]:
    """Score one finished Tier-1a run and return one table row.

    The distance metric is raw HeDFT ``R(t)``. The speed metric is the
    same-smoothed I2 speed reference that Method-B minimized against. The ledger
    residual is a wiring-correctness diagnostic, not a physics verdict.
    """
    run_dir = run if isinstance(run, RunDirectory) else RunDirectory(run)
    ion = run_dir.load_ion()
    return _score_ion(
        ion,
        case=case,
        variant=variant,
        n=n,
        run_tag=run_tag,
        scenario=scenario,
        t_star_ps=t_star_ps,
        hedft=hedft,
        smoothed=smoothed,
        window_start_ps=window_start_ps,
    )


def _load_references(
    project_root: Path,
    case: str,
) -> tuple[HedftTrajectory, SmoothedSpeedReference]:
    hedft = load_hedft_trajectory(
        project_root / "data" / "reference" / f"{case}_All_Data.csv"
    )
    smoothed = load_smoothed_speed_reference(
        project_root
        / "data"
        / "reference"
        / "drag"
        / case
        / "velocity_smoothed"
        / "cleaned_data_long.csv"
    )
    return hedft, smoothed


def collect_tier1a_records(
    *,
    project_root: Path = PROJECT_ROOT,
    case: str = CASE,
    variant: str = VARIANT,
    n: int = N,
    t_star_values_ps: Iterable[float] = tuple(T_STAR_VALUES_PS),
    window_start_ps: float = WINDOW_START_PS,
) -> tuple[list[Tier1aRunRecord], HedftTrajectory, SmoothedSpeedReference]:
    """Load the Tier-1a run matrix with table rows and references."""
    hedft, smoothed = _load_references(project_root, case)
    run_root = project_root / "data" / "runs"
    matrix: list[tuple[str, float | None]] = [("fixed", None)]
    matrix.extend(("anchored_discrete", float(t)) for t in t_star_values_ps)

    records: list[Tier1aRunRecord] = []
    for scenario, t_star_ps in matrix:
        run_tag = tier1a_run_tag(scenario, t_star_ps)
        run_dir = run_root / tier1a_run_dir_name(
            case, variant, n, scenario, t_star_ps
        )
        run = RunDirectory(run_dir)
        ion = run.load_ion()
        cfg = run.load_cfg() if run.has_cfg() else None
        row = _score_ion(
            ion,
            case=case,
            variant=variant,
            n=n,
            run_tag=run_tag,
            scenario=scenario,
            t_star_ps=t_star_ps,
            hedft=hedft,
            smoothed=smoothed,
            window_start_ps=window_start_ps,
        )
        records.append(
            Tier1aRunRecord(
                label=_run_label(scenario, t_star_ps),
                run_tag=run_tag,
                scenario=scenario,
                t_star_ps=t_star_ps,
                ion=ion,
                cfg=cfg,
                row=row,
            )
        )
    return records, hedft, smoothed


def collect_tier1a_rows(
    *,
    project_root: Path = PROJECT_ROOT,
    case: str = CASE,
    variant: str = VARIANT,
    n: int = N,
    t_star_values_ps: Iterable[float] = tuple(T_STAR_VALUES_PS),
    window_start_ps: float = WINDOW_START_PS,
) -> list[dict[str, Any]]:
    """Load the Tier-1a run matrix and return all table rows."""
    records, _, _ = collect_tier1a_records(
        project_root=project_root,
        case=case,
        variant=variant,
        n=n,
        t_star_values_ps=t_star_values_ps,
        window_start_ps=window_start_ps,
    )
    return [record.row for record in records]


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Return a plain-text table with the Tier-1a columns."""
    if not rows:
        return "(no Tier-1a rows)"
    widths = {
        col: max(len(col), *(len(_format_value(row[col])) for row in rows))
        for col in TIER1A_TABLE_COLUMNS
    }
    header = "  ".join(col.ljust(widths[col]) for col in TIER1A_TABLE_COLUMNS)
    sep = "  ".join("-" * widths[col] for col in TIER1A_TABLE_COLUMNS)
    body = [
        "  ".join(_format_value(row[col]).ljust(widths[col])
                  for col in TIER1A_TABLE_COLUMNS)
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write table rows to CSV using the stable Tier-1a column order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIER1A_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return p


def build_trajectory_figure(
    records: list[Tier1aRunRecord],
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window: tuple[float, float],
    *,
    selected_t_star_ps: float = PLOT_T_STAR_PS,
    title: str | None = None,
):
    """Build a clear Tier-1a |v2| comparison for one selected anchored t*."""
    import matplotlib.pyplot as plt

    selected = _selected_velocity_records(records, selected_t_star_ps)
    fig, ax_v = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    t_start, t_end = window

    for record in selected:
        t_md, _, _, v2_md = ensemble_mean_series(record.ion)
        ax_v.plot(t_md, v2_md, lw=1.4, ls="--",
                  label=f"{record.label} MD mean |v2|")

    ax_v.plot(
        smoothed.time_ps,
        smoothed.speed_Aps,
        color="black",
        lw=1.6,
        ls=":",
        label="CEEMDAN+SG |v2|",
    )
    ax_v.plot(
        hedft.time_ps,
        hedft.v2_magnitude_Aps,
        color="0.15",
        lw=1.2,
        alpha=0.65,
        label="HeDFT |v2|",
    )
    ax_v.axvspan(t_start, t_end, color="tab:green", alpha=0.12,
                 label="scored window")
    ax_v.set_ylabel(r"$|v|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.set_xlabel("t / ps")
    ax_v.legend(frameon=False, ncol=2)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)

    fig.suptitle(
        title
        or (
            f"Tier-1a |v2| comparison ({hedft.droplet_radius_A:.0f} A, "
            f"anchored t*={selected_t_star_ps:.1f} ps)"
        )
    )
    return fig


def _selected_velocity_records(
    records: list[Tier1aRunRecord],
    selected_t_star_ps: float,
) -> list[Tier1aRunRecord]:
    """Return fixed plus the selected anchored record for the |v2| plot."""
    fixed = [record for record in records if record.scenario == "fixed"]
    anchored = [
        record
        for record in records
        if (
            record.scenario == "anchored_discrete"
            and record.t_star_ps is not None
            and np.isclose(record.t_star_ps, selected_t_star_ps)
        )
    ]
    if len(fixed) != 1:
        raise ValueError(f"expected exactly one fixed Tier-1a record, got {len(fixed)}")
    if len(anchored) != 1:
        available = [
            record.t_star_ps
            for record in records
            if record.scenario == "anchored_discrete"
        ]
        raise ValueError(
            f"expected exactly one anchored record for t*={selected_t_star_ps}; "
            f"available t* values are {available}"
        )
    return [fixed[0], anchored[0]]


def build_positions_figure(records: list[Tier1aRunRecord]):
    """Build mean x/y/z position panels for all Tier-1a runs."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        3, 1, figsize=(9.5, 7.0), sharex=True, constrained_layout=True
    )
    labels = ("x", "y", "z")
    attrs = ("positions_x", "positions_y", "positions_z")
    for ax, label, attr in zip(axes, labels, attrs):
        for record in records:
            arr = getattr(record.ion, attr)
            n = record.ion.num_molecules
            ax.plot(record.ion.time_ps, np.mean(arr[:n], axis=0), lw=1.2,
                    label=f"{record.label} {label}1")
            ax.plot(record.ion.time_ps, np.mean(arr[n:], axis=0), lw=1.2,
                    ls="--", label=f"{record.label} {label}2")
        ax.set_ylabel(fr"{label} / $\mathrm{{\AA}}$")
        ax.legend(frameon=False, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].set_xlabel("t / ps")
    fig.suptitle("Tier-1a mean positions")
    return fig


def build_energy_figures(
    records: list[Tier1aRunRecord],
    window: tuple[float, float],
) -> list:
    """Build the Tier-0-style energy figure for each run with a saved cfg."""
    figures = []
    for record in records:
        if record.cfg is None:
            continue
        fig = plot_energy_analysis(record.ion, record.cfg, window)
        fig.suptitle(f"Tier-1a energy balance: {record.label}")
        figures.append(fig)
    return figures


def build_force_figures(
    records: list[Tier1aRunRecord],
    window: tuple[float, float],
    *,
    atom_index: int = FORCE_ATOM_INDEX,
) -> list:
    """Build radial-projected force figures for each run with a saved cfg."""
    import matplotlib.pyplot as plt

    figures = []
    for record in records:
        if record.cfg is None:
            continue
        drag_r, coul_r, drop_r = _reconstruct_radial_forces(
            record.ion,
            record.cfg,
            atom_index=atom_index,
        )
        net_force = coul_r + drop_r + drag_r
        t_start, t_end = window
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(8.0, 7.0))
        ax_top.axhline(0.0, color="0.6", lw=0.8)
        ax_top.plot(record.ion.time_ps, coul_r, color="tab:red", label="Coulomb")
        ax_top.plot(record.ion.time_ps, drop_r, color="tab:blue",
                    label="droplet (confining)")
        ax_top.plot(record.ion.time_ps, drag_r, color="tab:green", label="drag")
        ax_top.axvspan(t_start, t_end, color="tab:green", alpha=0.12,
                       label="scored window")
        ax_top.set_ylabel(r"$F \cdot \hat{r}$ / $\mathrm{amu}\,\mathrm{\AA}/\mathrm{ps}^2$")
        ax_top.set_title(f"{record.label}: radial-projected forces, atom {atom_index + 1}")
        ax_top.legend(frameon=False, ncol=2)
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        ax_bottom.axhline(0.0, color="0.6", lw=0.8)
        ax_bottom.plot(record.ion.time_ps, net_force, color="tab:purple",
                       label="net force", lw=2)
        ax_bottom.axvspan(t_start, t_end, color="tab:green", alpha=0.12)
        ax_bottom.set_xlabel("t / ps")
        ax_bottom.set_ylabel(r"$F_{\mathrm{net}} \cdot \hat{r}$ / $\mathrm{amu}\,\mathrm{\AA}/\mathrm{ps}^2$")
        ax_bottom.legend(frameon=False)
        ax_bottom.spines["top"].set_visible(False)
        ax_bottom.spines["right"].set_visible(False)
        fig.tight_layout()
        figures.append(fig)
    return figures


def export_tier1a_mean_series(
    directory: str | Path,
    records: list[Tier1aRunRecord],
) -> list[Path]:
    """Write one Tier-0-format mean-series CSV per Tier-1a run."""
    out_dir = Path(directory)
    paths: list[Path] = []
    for record in records:
        path = out_dir / f"{record.run_tag}_mean_trajectory.csv"
        provenance = (
            "Tier-1a MD ensemble-mean trajectory. "
            f"run_tag={record.run_tag}, scenario={record.scenario}, "
            f"t_star_ps={record.t_star_ps}, N={record.ion.num_molecules}, "
            f"ion_steps={record.ion.time_ps.size}."
        )
        export_mean_series(path, record.ion, provenance=provenance)
        paths.append(path)
    return paths


def main() -> int:
    """Print the Tier-1a RMSE table and optionally save it as CSV."""
    records, hedft, smoothed = collect_tier1a_records()
    rows = [record.row for record in records]
    window = (WINDOW_START_PS, float(smoothed.time_ps[-1]))
    print(
        "Tier-1a RMSE table: speed metric is same-smoothed I2; "
        "ledger closure is a wiring diagnostic, not a physics verdict."
    )
    print(format_table(rows))
    if SAVE_CSV_PATH is not None:
        path = write_rows_csv(SAVE_CSV_PATH, rows)
        print(f"Wrote Tier-1a RMSE table -> {path}")
    if EXPORT_MEAN_SERIES_DIR is not None:
        for path in export_tier1a_mean_series(EXPORT_MEAN_SERIES_DIR, records):
            print(f"Wrote Tier-1a mean series -> {path}")
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt

        build_trajectory_figure(
            records,
            hedft,
            smoothed,
            window,
            selected_t_star_ps=PLOT_T_STAR_PS,
        )
        if POSITIONS_FIGURE:
            build_positions_figure(records)
        if ENERGY_FIGURE:
            build_energy_figures(records, window)
        if FORCE_FIGURE:
            build_force_figures(records, window, atom_index=FORCE_ATOM_INDEX)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
