"""Run the single-pulse neutral + ion simulation pipeline.

This is the public Step-12 entry point. It deliberately stays narrow:
configure the canonical single-pulse preset, run neutral propagation, run ion
propagation from the neutral checkpoint, and save artifacts through
``RunDirectory``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from i2_helium_md import single_pulse_N2000  # noqa: E402
from i2_helium_md.simulation.ion import (  # noqa: E402
    DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    run_ion_propagation,
)
from i2_helium_md.simulation.neutral import (  # noqa: E402
    DEFAULT_MAX_CHECKPOINT_BYTES,
    run_neutral_propagation,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Run the single-pulse neutral + ion MD pipeline.",
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="Output run directory for cfg.json, neutral.npz, and ion.npz.",
    )
    parser.add_argument(
        "--num-molecules",
        type=int,
        default=None,
        help="Override the single_pulse_N2000 molecule count.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Override the random seed. Omit for fresh randomness.",
    )
    parser.add_argument(
        "--neutral-max-bytes",
        type=int,
        default=DEFAULT_MAX_CHECKPOINT_BYTES,
        help="Neutral checkpoint memory budget before storage striding.",
    )
    parser.add_argument(
        "--ion-max-bytes",
        type=int,
        default=DEFAULT_MAX_CHECKPOINT_BYTES_ION,
        help="Ion checkpoint memory budget before storage striding.",
    )
    parser.add_argument(
        "--ion-simulation-time",
        type=float,
        default=None,
        help="Override ion simulation time in ps, mainly for smoke runs.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting an existing cfg.json, neutral.npz, or ion.npz.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print driver-level progress from the neutral and ion stages.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the single-pulse pipeline and return a process exit code."""
    args = build_parser().parse_args(argv)
    run = RunDirectory(args.run_dir)

    _check_outputs(run, force=args.force)

    overrides = {}
    if args.num_molecules is not None:
        overrides["num_molecules"] = args.num_molecules
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.ion_simulation_time is not None:
        overrides["ion_simulation_time"] = args.ion_simulation_time

    cfg = single_pulse_N2000(**overrides)
    cfg.validate()

    if args.force:
        run.save_cfg(cfg)

    print(f"run directory: {run.root}")
    print("running neutral propagation")
    neutral = run_neutral_propagation(
        cfg,
        run_dir=run,
        max_bytes=args.neutral_max_bytes,
        verbose=args.verbose,
    )

    print("running ion propagation")
    ion = run_ion_propagation(
        cfg,
        neutral,
        run_dir=run,
        max_bytes=args.ion_max_bytes,
        verbose=args.verbose,
    )

    print("single-pulse run complete")
    print(f"cfg:     {run.cfg_path}")
    print(f"neutral: {run.neutral_path} ({neutral.time_ps.size} stored steps)")
    print(f"ion:     {run.ion_path} ({ion.time_ps.size} stored steps)")
    return 0


def _check_outputs(run: RunDirectory, *, force: bool) -> None:
    """Refuse to overwrite existing run artifacts unless requested."""
    existing = [
        path
        for path in (run.cfg_path, run.neutral_path, run.ion_path)
        if path.exists()
    ]
    if existing and not force:
        existing_text = ", ".join(str(path) for path in existing)
        raise SystemExit(
            "run directory already contains output artifacts; use --force "
            f"to overwrite them: {existing_text}"
        )


if __name__ == "__main__":
    raise SystemExit(main())
