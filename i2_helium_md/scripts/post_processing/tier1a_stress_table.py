"""Emit the Tier-1a onset-strip stress table for finished run directories.

This scorer is intentionally separate from ``tier1a_rmse_table.py``. It keeps
the diagnostic onset-strip stress rows in their own table while still loading
the fixed Tier-1a run for trajectory plotting and mean-series export.

Usage::

    python scripts/post_processing/tier1a_stress_table.py
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
WINDOW_START_PS = 2.67
PLOT_N_FINAL = 0

SAVE_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier1a_stress_table.csv"
EXPORT_MEAN_SERIES_DIR = PROJECT_ROOT / "data" / "reference" / "drag" / CASE / "tier1a_stress"

SHOW_FIGURE = True


from scripts.tier1a_common import (  # noqa: E402
    TIER1A_STRESS_N_FINAL_VALUES,
    TIER1A_STRESS_T_STRIP_PS,
    tier1a_run_dir_name,
    tier1a_run_tag,
    tier1a_stress_run_dir_name,
    tier1a_stress_run_tag,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    SmoothedSpeedReference,
    compare_distance,
    compare_speed_to_reference,
    ion_ledger_closure,
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    ensemble_mean_series,
    export_mean_series,
)


T_STRIP_PS = TIER1A_STRESS_T_STRIP_PS
N_FINAL_VALUES = TIER1A_STRESS_N_FINAL_VALUES


TIER1A_STRESS_TABLE_COLUMNS = [
    "case",
    "variant",
    "N",
    "run_tag",
    "stress_family",
    "n_final_requested",
    "t_strip_ps",
    "window_start_ps",
    "window_end_ps",
    "n_scored_R",
    "n_scored_v2_smoothed",
    "RMSE_R_raw_A",
    "RMSE_v2_smoothed_Aps",
    "v2_smoothed_mean_ratio",
    "ledger_max_resid_eV",
    "n_removed",
    "n_shell_start",
    "n_shell_end",
]


@dataclass(frozen=True)
class Tier1aStressRunRecord:
    """Loaded Tier-1a stress run plus metadata for tables and figures."""

    label: str
    run_tag: str
    n_final: int | None
    t_strip_ps: float | None
    ion: IonCheckpoint
    cfg: SimConfig | None
    row: dict[str, Any]


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


def _shell_summary(n_shell: np.ndarray) -> tuple[int, int, int]:
    """Return ``(n_start, n_end, n_removed)`` from a v6 ``n_shell`` array."""
    if n_shell.ndim != 2 or n_shell.shape[1] == 0:
        raise ValueError(
            f"n_shell must have shape (2N, T) with T > 0, got {n_shell.shape}"
        )
    start_values = np.asarray(n_shell[:, 0], dtype=float)
    end_values = np.asarray(n_shell[:, -1], dtype=float)
    if not np.all(np.isfinite(start_values)) or not np.all(np.isfinite(end_values)):
        raise ValueError("n_shell start/end values must be finite")
    start_rounded = np.rint(start_values)
    end_rounded = np.rint(end_values)
    if not np.allclose(start_values, start_rounded, rtol=0.0, atol=1e-9):
        raise ValueError(
            f"n_shell start values must be integer-valued: {start_values.tolist()}"
        )
    if not np.allclose(end_values, end_rounded, rtol=0.0, atol=1e-9):
        raise ValueError(
            f"n_shell end values must be integer-valued: {end_values.tolist()}"
        )
    starts = start_rounded.astype(int)
    ends = end_rounded.astype(int)
    if not np.all(starts == starts[0]):
        raise ValueError(f"n_shell start values are not uniform: {starts.tolist()}")
    if not np.all(ends == ends[0]):
        raise ValueError(f"n_shell end values are not uniform: {ends.tolist()}")
    n_start = int(starts[0])
    n_end = int(ends[0])
    return n_start, n_end, max(0, n_start - n_end)


def _raise_cfg_mismatch(run_tag: str, field: str, expected: Any, actual: Any) -> None:
    raise ValueError(
        f"{run_tag} cfg metadata mismatch: {field} expected {expected!r}, "
        f"got {actual!r}"
    )


def _raise_ion_mismatch(run_tag: str, field: str, expected: Any, actual: Any) -> None:
    raise ValueError(
        f"{run_tag} ion metadata mismatch: {field} expected {expected!r}, "
        f"got {actual!r}"
    )


def _load_required_cfg(run: RunDirectory, run_tag: str) -> SimConfig:
    if not run.has_cfg():
        raise ValueError(
            f"{run_tag} missing cfg.json in {run.root}; refusing to score from "
            "path metadata alone."
        )
    return run.load_cfg()


def _validate_fixed_cfg(run_tag: str, cfg: SimConfig) -> None:
    if cfg.mass_scenario != "fixed":
        _raise_cfg_mismatch(run_tag, "mass_scenario", "fixed", cfg.mass_scenario)


def _validate_fixed_ion(run_tag: str, ion: IonCheckpoint) -> None:
    if ion.mass_scenario != "fixed":
        _raise_ion_mismatch(run_tag, "ion.mass_scenario", "fixed", ion.mass_scenario)


def _validate_stress_cfg(
    run_tag: str,
    cfg: SimConfig,
    *,
    n_final: int,
    t_strip_ps: float,
) -> None:
    if cfg.mass_scenario != "anchored_discrete":
        _raise_cfg_mismatch(
            run_tag,
            "mass_scenario",
            "anchored_discrete",
            cfg.mass_scenario,
        )
    if cfg.anchor_mode != "onset_strip":
        _raise_cfg_mismatch(run_tag, "anchor_mode", "onset_strip", cfg.anchor_mode)
    if cfg.anchor_n_final != int(n_final):
        _raise_cfg_mismatch(run_tag, "anchor_n_final", int(n_final), cfg.anchor_n_final)
    if not np.isclose(float(cfg.t_star_ps), float(t_strip_ps)):
        _raise_cfg_mismatch(run_tag, "t_star_ps", float(t_strip_ps), cfg.t_star_ps)


def _validate_stress_call_metadata(
    *,
    run_tag: str,
    n_final: int,
    t_strip_ps: float,
) -> tuple[int, float]:
    expected_run_tag = tier1a_stress_run_tag(
        n_final=n_final,
        t_strip_ps=t_strip_ps,
    )
    if run_tag != expected_run_tag:
        raise ValueError(
            f"{run_tag} stress metadata mismatch: run_tag expected "
            f"{expected_run_tag!r} for n_final={n_final!r}, "
            f"t_strip_ps={t_strip_ps!r}"
        )
    return int(n_final), float(t_strip_ps)


def _validate_stress_ion(
    run_tag: str,
    ion: IonCheckpoint,
    *,
    n_final: int,
) -> tuple[int, int, int]:
    if ion.mass_scenario != "anchored_discrete":
        _raise_ion_mismatch(
            run_tag,
            "ion.mass_scenario",
            "anchored_discrete",
            ion.mass_scenario,
        )
    n_shell_start, n_shell_end, n_removed = _shell_summary(ion.n_shell)
    if n_shell_end != int(n_final):
        _raise_ion_mismatch(run_tag, "n_shell_end", int(n_final), n_shell_end)
    return n_shell_start, n_shell_end, n_removed


def _score_tier1a_stress_ion(
    ion: IonCheckpoint,
    *,
    case: str,
    variant: str,
    n: int,
    run_tag: str,
    n_final: int,
    t_strip_ps: float,
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window_start_ps: float,
    shell_summary: tuple[int, int, int] | None = None,
) -> dict[str, Any]:
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
    if shell_summary is None:
        shell_summary = _shell_summary(ion.n_shell)
    n_shell_start, n_shell_end, n_removed = shell_summary

    return {
        "case": case,
        "variant": variant,
        "N": int(n),
        "run_tag": run_tag,
        "stress_family": "onset_strip",
        "n_final_requested": int(n_final),
        "t_strip_ps": float(t_strip_ps),
        "window_start_ps": window[0],
        "window_end_ps": window[1],
        "n_scored_R": dist.num_overlap_points,
        "n_scored_v2_smoothed": v2_smoothed.num_overlap_points,
        "RMSE_R_raw_A": dist.rmse,
        "RMSE_v2_smoothed_Aps": v2_smoothed.rmse,
        "v2_smoothed_mean_ratio": v2_smoothed.mean_ratio,
        "ledger_max_resid_eV": closure.max_abs_residual_eV,
        "n_removed": n_removed,
        "n_shell_start": n_shell_start,
        "n_shell_end": n_shell_end,
    }


def score_tier1a_stress_run(
    run: RunDirectory | str | Path,
    *,
    case: str,
    variant: str,
    n: int,
    run_tag: str,
    n_final: int,
    t_strip_ps: float,
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window_start_ps: float = WINDOW_START_PS,
    cfg: SimConfig | None = None,
) -> dict[str, Any]:
    """Score one finished Tier-1a onset-strip stress run."""
    n_final_value, t_strip_ps_value = _validate_stress_call_metadata(
        run_tag=run_tag,
        n_final=n_final,
        t_strip_ps=t_strip_ps,
    )
    run_dir = run if isinstance(run, RunDirectory) else RunDirectory(run)
    if cfg is None:
        cfg = _load_required_cfg(run_dir, run_tag)
    _validate_stress_cfg(
        run_tag,
        cfg,
        n_final=n_final_value,
        t_strip_ps=t_strip_ps_value,
    )
    ion = run_dir.load_ion(cfg=cfg)
    shell_summary = _validate_stress_ion(
        run_tag,
        ion,
        n_final=n_final_value,
    )
    return _score_tier1a_stress_ion(
        ion,
        case=case,
        variant=variant,
        n=n,
        run_tag=run_tag,
        n_final=n_final_value,
        t_strip_ps=t_strip_ps_value,
        hedft=hedft,
        smoothed=smoothed,
        window_start_ps=window_start_ps,
        shell_summary=shell_summary,
    )


def _run_label(n_final: int | None) -> str:
    if n_final is None:
        return "fixed"
    return f"onset strip n={n_final}"


def collect_tier1a_stress_records(
    *,
    project_root: Path = PROJECT_ROOT,
    case: str = CASE,
    variant: str = VARIANT,
    n: int = N,
    n_final_values: Iterable[int] = tuple(N_FINAL_VALUES),
    t_strip_ps: float = T_STRIP_PS,
    window_start_ps: float = WINDOW_START_PS,
) -> tuple[list[Tier1aStressRunRecord], HedftTrajectory, SmoothedSpeedReference]:
    """Load fixed plus onset-strip stress runs with rows only for stress runs."""
    hedft, smoothed = _load_references(project_root, case)
    run_root = project_root / "data" / "runs"

    fixed_tag = tier1a_run_tag("fixed", None)
    fixed_run = RunDirectory(
        run_root / tier1a_run_dir_name(case, variant, n, "fixed", None)
    )
    fixed_cfg = _load_required_cfg(fixed_run, fixed_tag)
    _validate_fixed_cfg(fixed_tag, fixed_cfg)
    fixed_ion = fixed_run.load_ion(cfg=fixed_cfg)
    _validate_fixed_ion(fixed_tag, fixed_ion)
    records = [
        Tier1aStressRunRecord(
            label=_run_label(None),
            run_tag=fixed_tag,
            n_final=None,
            t_strip_ps=None,
            ion=fixed_ion,
            cfg=fixed_cfg,
            row={},
        )
    ]

    for n_final in n_final_values:
        run_tag = tier1a_stress_run_tag(
            n_final=n_final,
            t_strip_ps=t_strip_ps,
        )
        n_final_value, t_strip_ps_value = _validate_stress_call_metadata(
            run_tag=run_tag,
            n_final=n_final,
            t_strip_ps=t_strip_ps,
        )
        run = RunDirectory(
            run_root
            / tier1a_stress_run_dir_name(
                case,
                variant,
                n,
                n_final=n_final,
                t_strip_ps=t_strip_ps,
            )
        )
        cfg = _load_required_cfg(run, run_tag)
        _validate_stress_cfg(
            run_tag,
            cfg,
            n_final=n_final_value,
            t_strip_ps=t_strip_ps_value,
        )
        ion = run.load_ion(cfg=cfg)
        shell_summary = _validate_stress_ion(
            run_tag,
            ion,
            n_final=n_final_value,
        )
        row = _score_tier1a_stress_ion(
            ion,
            case=case,
            variant=variant,
            n=n,
            run_tag=run_tag,
            n_final=n_final_value,
            t_strip_ps=t_strip_ps_value,
            hedft=hedft,
            smoothed=smoothed,
            window_start_ps=window_start_ps,
            shell_summary=shell_summary,
        )
        records.append(
            Tier1aStressRunRecord(
                label=_run_label(n_final_value),
                run_tag=run_tag,
                n_final=n_final_value,
                t_strip_ps=t_strip_ps_value,
                ion=ion,
                cfg=cfg,
                row=row,
            )
        )

    return records, hedft, smoothed


def _selected_velocity_records(
    records: list[Tier1aStressRunRecord],
    selected_n_final: int,
) -> list[Tier1aStressRunRecord]:
    """Return fixed plus the selected onset-strip record for the |v2| plot."""
    fixed = [record for record in records if record.n_final is None]
    selected = [
        record
        for record in records
        if record.n_final is not None and record.n_final == int(selected_n_final)
    ]
    if len(fixed) != 1:
        raise ValueError(
            f"expected exactly one fixed Tier-1a record, got {len(fixed)}"
        )
    if len(selected) != 1:
        available = [
            record.n_final for record in records if record.n_final is not None
        ]
        raise ValueError(
            f"expected exactly one stress record for n_final={selected_n_final}; "
            f"available n_final values are {available}"
        )
    return [fixed[0], selected[0]]


def build_stress_trajectory_figure(
    records: list[Tier1aStressRunRecord],
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window: tuple[float, float],
    *,
    selected_n_final: int = PLOT_N_FINAL,
    title: str | None = None,
):
    """Build a Tier-1a stress |v2| comparison for one selected endpoint."""
    import matplotlib.pyplot as plt

    selected = _selected_velocity_records(records, selected_n_final)
    fig, ax_v = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    t_start, t_end = window

    for record in selected:
        t_md, _, _, v2_md = ensemble_mean_series(record.ion)
        ax_v.plot(
            t_md,
            v2_md,
            lw=1.4,
            ls="--",
            label=f"{record.label} MD mean |v2|",
        )

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
    ax_v.axvspan(t_start, t_end, color="tab:green", alpha=0.12, label="scored window")
    ax_v.set_ylabel(r"$|v|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.set_xlabel("t / ps")
    ax_v.legend(frameon=False, ncol=2)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)
    fig.suptitle(
        title
        or (
            f"Tier-1a onset-strip stress |v2| comparison "
            f"({hedft.droplet_radius_A:.0f} A, n_final={selected_n_final})"
        )
    )
    return fig


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Return a plain-text table with the Tier-1a stress columns."""
    if not rows:
        return "(no Tier-1a stress rows)"
    widths = {
        col: max(len(col), *(len(_format_value(row[col])) for row in rows))
        for col in TIER1A_STRESS_TABLE_COLUMNS
    }
    header = "  ".join(
        col.ljust(widths[col]) for col in TIER1A_STRESS_TABLE_COLUMNS
    )
    sep = "  ".join("-" * widths[col] for col in TIER1A_STRESS_TABLE_COLUMNS)
    body = [
        "  ".join(
            _format_value(row[col]).ljust(widths[col])
            for col in TIER1A_STRESS_TABLE_COLUMNS
        )
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write stress table rows to CSV using the stable column order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIER1A_STRESS_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return p


def export_tier1a_stress_mean_series(
    directory: str | Path,
    records: list[Tier1aStressRunRecord],
) -> list[Path]:
    """Write one Tier-0-format mean-series CSV per loaded stress comparison run."""
    out_dir = Path(directory)
    paths: list[Path] = []
    for record in records:
        path = out_dir / f"{record.run_tag}_mean_trajectory.csv"
        provenance = (
            "Tier-1a onset-strip stress MD ensemble-mean trajectory. "
            f"run_tag={record.run_tag}, n_final={record.n_final}, "
            f"t_strip_ps={record.t_strip_ps}, N={record.ion.num_molecules}, "
            f"ion_steps={record.ion.time_ps.size}."
        )
        export_mean_series(path, record.ion, provenance=provenance)
        paths.append(path)
    return paths


def main() -> int:
    """Print the Tier-1a onset-strip stress table and optional artifacts."""
    records, hedft, smoothed = collect_tier1a_stress_records()
    rows = [record.row for record in records if record.row]
    window = (WINDOW_START_PS, float(smoothed.time_ps[-1]))
    print(
        "Tier-1a onset-strip stress table: speed metric is same-smoothed I2; "
        "fixed run is included only for plotting/export comparison."
    )
    print(format_table(rows))
    if SAVE_CSV_PATH is not None:
        path = write_rows_csv(SAVE_CSV_PATH, rows)
        print(f"Wrote Tier-1a stress table -> {path}")
    if EXPORT_MEAN_SERIES_DIR is not None:
        for path in export_tier1a_stress_mean_series(
            EXPORT_MEAN_SERIES_DIR,
            records,
        ):
            print(f"Wrote Tier-1a stress mean series -> {path}")
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt

        build_stress_trajectory_figure(
            records,
            hedft,
            smoothed,
            window,
            selected_n_final=PLOT_N_FINAL,
        )
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
