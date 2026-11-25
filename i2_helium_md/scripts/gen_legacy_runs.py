"""Generate legacy hard-sphere runs: the thesis runs, or one run of your own.

Open this file, edit the values in USER SETTINGS, then run it from your editor or
with:

    python scripts/gen_legacy_runs.py

Every run directory holds ``cfg.json``, ``neutral.npz`` and ``ion.npz``.

* ``INPUT_PRESET = "thesis_legacy"`` makes the three seeded runs behind the legacy
  figures (``presets.LEGACY_THESIS_RUNS``) in ``data/runs/<name>``; the name of one of
  them makes just that one. Size, seed and directory follow from the name. A finished
  run is skipped; a partial one or one with another model is refused. The runs of the
  submitted figures were unseeded, so these agree with them statistically only.
* Any other preset makes one run in ``RUN_DIR`` with your settings. Start with
  ``RUN_SIZE = "smoke"``.
"""

from __future__ import annotations

from pathlib import Path
import sys


# Absolute path to this project folder. Use this for all default paths so the
# script writes to the same place whether it is launched from PyCharm, a
# terminal in the repo root, or a terminal inside scripts/.
PROJECT_ROOT = Path(__file__).resolve().parents[1]


# =============================================================================
# USER SETTINGS
# =============================================================================
# This is the section you should edit for normal use. You do not need command
# line arguments.

# "thesis_legacy"                       all three thesis runs below
# "legacy_9A_N2000"                     thesis run, HeDFT comparison 9 A
# "legacy_18A_N2000"                    thesis run, HeDFT comparison 18 A
# "legacy_droplet_N20000"               thesis run, droplet size distribution
#
# Migrated MATLAB input files, for your own runs:
# "single_pulse_N2000"                  HeDFT comparison 9 A
# "single_pulse_N2000_18Angst"          HeDFT comparison 18 A
# "single_pulse_droplet_distribution"   droplet size distribution, as the experiment
INPUT_PRESET = "thesis_legacy"

# Choose "smoke", "custom", or "production".
#
# "smoke"      quick check that the full pipeline works
# "custom"     uses NUM_MOLECULES, SEED, and ION_TIME_PS below
# "production" uses the selected INPUT_PRESET exactly, unless you override
#              values in the PRODUCTION OVERRIDES section
#
# Thesis runs: "production", or "smoke" for NUM_MOLECULES and ION_TIME_PS in
# data/runs/smoke.
RUN_SIZE = "production"

# Where your own run is written (thesis runs have fixed names). Change the final folder name for each run
# you want to keep. This intentionally points to the project-level data/runs
# folder, not to scripts/data/runs and not to a top-level results folder.
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "legacy_custom"

# If False, the script stops when RUN_DIR already contains outputs. This
# prevents accidental overwrites. Set True only when you intentionally want to
# rerun into the same folder.
OVERWRITE_EXISTING_RUN = False

# Custom/smoke settings. These are ignored by production mode unless you copy
# them into the production overrides below.
NUM_MOLECULES = 10
SEED = 123
ION_TIME_PS = 0.1

# Print internal driver details such as storage stride and checkpoint size.
VERBOSE = True


# =============================================================================
# PRODUCTION OVERRIDES
# =============================================================================
# Leave these as None to use the selected INPUT_PRESET exactly (required for thesis runs).

# Set one of these values only if you want a production-like run with a specific
# override.
PRODUCTION_NUM_MOLECULES = None
PRODUCTION_SEED = None
PRODUCTION_ION_TIME_PS = None

PRODUCTION_CROSSSECTION = None
PRODUCTION_BINDING_ENERGY = None
PRODUCTION_MASS_ATTACH_PROBABILITY = None


# =============================================================================
# ADVANCED SETTINGS
# =============================================================================
# The neutral and ion drivers automatically store a strided trajectory if a full
# checkpoint would be too large. Most users should leave these budgets alone.
NEUTRAL_MAX_BYTES = 300_000_000
ION_MAX_BYTES = 1_000_000_000


# =============================================================================
# IMPORT SETUP
# =============================================================================
# Running this file directly puts scripts/ on sys.path. Add the repository root
# so the local i2_helium_md package imports reliably even if it is not installed.
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md import (  # noqa: E402
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
    single_pulse_droplet_distribution,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.presets import LEGACY_THESIS_RUNS  # noqa: E402
from i2_helium_md.simulation.battery import run_neutral_and_ion  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT  # noqa: E402

OUTPUT_ROOT = RUNS_ROOT
PRODUCTION_OVERRIDES = (
    "PRODUCTION_NUM_MOLECULES", "PRODUCTION_SEED", "PRODUCTION_ION_TIME_PS",
    "PRODUCTION_CROSSSECTION", "PRODUCTION_BINDING_ENERGY",
    "PRODUCTION_MASS_ATTACH_PROBABILITY",
)


PRESET_BUILDERS = {
    "single_pulse_N2000": single_pulse_N2000,
    "single_pulse_N2000_18Angst": single_pulse_N2000_18Angst,
    "single_pulse_droplet_distribution": single_pulse_droplet_distribution,
}


def main() -> int:
    """Make the run or runs selected at the top of this file."""
    if INPUT_PRESET == "thesis_legacy" or INPUT_PRESET in LEGACY_THESIS_RUNS:
        for run_dir, thesis_cfg in thesis_runs().items():
            run_neutral_and_ion(thesis_cfg, run_dir)
        return 0
    cfg = build_config()
    run = RunDirectory(RUN_DIR)

    refuse_accidental_overwrite(run)
    print_run_header(cfg, run)

    # Validate before writing files. This catches impossible settings early.
    cfg.validate()

    # RunDirectory normally keeps the first cfg.json it wrote. For an explicit
    # overwrite run, replace cfg.json before the drivers save checkpoints so
    # all three output files describe the same rerun.
    if OVERWRITE_EXISTING_RUN:
        run.save_cfg(cfg)

    print("[1/2] Neutral propagation")
    neutral = run_neutral_propagation(
        cfg,
        run_dir=run,
        max_bytes=NEUTRAL_MAX_BYTES,
        verbose=VERBOSE,
    )

    print("[2/2] Ion propagation")
    ion = run_ion_propagation(
        cfg,
        neutral,
        run_dir=run,
        max_bytes=ION_MAX_BYTES,
        verbose=VERBOSE,
    )

    print_run_footer(run, neutral_steps=neutral.time_ps.size, ion_steps=ion.time_ps.size)
    return 0


def thesis_runs() -> dict[Path, SimConfig]:
    """The selected thesis runs: directory -> config.

    Raises ``SystemExit`` naming a setting that contradicts the fixed size, seed and
    directory of a thesis run.
    """
    set_overrides = [name for name in PRODUCTION_OVERRIDES if globals()[name] is not None]
    if set_overrides:
        raise SystemExit(f"{', '.join(set_overrides)} must be None for a thesis run")
    if OVERWRITE_EXISTING_RUN:
        raise SystemExit("OVERWRITE_EXISTING_RUN must be False for a thesis run")
    if RUN_SIZE not in ("production", "smoke"):
        raise SystemExit(
            f"RUN_SIZE must be 'production' or 'smoke' for a thesis run, got {RUN_SIZE!r}")
    names = list(LEGACY_THESIS_RUNS) if INPUT_PRESET == "thesis_legacy" else [INPUT_PRESET]
    if RUN_SIZE == "smoke":
        size = {"num_molecules": NUM_MOLECULES, "ion_simulation_time": ION_TIME_PS}
        return {OUTPUT_ROOT / "smoke" / n: LEGACY_THESIS_RUNS[n](**size) for n in names}
    return {OUTPUT_ROOT / n: LEGACY_THESIS_RUNS[n]() for n in names}


def build_config():
    """Create the SimConfig for the selected input preset and RUN_SIZE."""
    preset_builder = get_preset_builder(INPUT_PRESET)

    if RUN_SIZE == "smoke":
        return preset_builder(
            num_molecules=NUM_MOLECULES,
            seed=SEED,
            ion_simulation_time=ION_TIME_PS,
        )

    if RUN_SIZE == "custom":
        return preset_builder(
            num_molecules=NUM_MOLECULES,
            seed=SEED,
            ion_simulation_time=ION_TIME_PS,
        )

    if RUN_SIZE == "production":
        overrides = {}
        if PRODUCTION_NUM_MOLECULES is not None:
            overrides["num_molecules"] = PRODUCTION_NUM_MOLECULES
        if PRODUCTION_SEED is not None:
            overrides["seed"] = PRODUCTION_SEED
        if PRODUCTION_ION_TIME_PS is not None:
            overrides["ion_simulation_time"] = PRODUCTION_ION_TIME_PS
        if PRODUCTION_CROSSSECTION is not None:
            overrides["geometric_scattering_crosssection_Iplus"] = PRODUCTION_CROSSSECTION
        if PRODUCTION_BINDING_ENERGY is not None:
            overrides["binding_energy_I_ion_eV"] = PRODUCTION_BINDING_ENERGY
        if PRODUCTION_MASS_ATTACH_PROBABILITY is not None:
            overrides["mass_attach_probability"] = PRODUCTION_MASS_ATTACH_PROBABILITY
        return preset_builder(**overrides)

    raise ValueError(
        f"unknown RUN_SIZE {RUN_SIZE!r}; use 'smoke', 'custom', or 'production'"
    )


def get_preset_builder(name: str):
    """Return the preset builder for a migrated MATLAB input file name."""
    try:
        return PRESET_BUILDERS[name]
    except KeyError as exc:
        choices = ", ".join(["thesis_legacy", *LEGACY_THESIS_RUNS, *sorted(PRESET_BUILDERS)])
        raise ValueError(
            f"unknown INPUT_PRESET {name!r}; use one of: {choices}"
        ) from exc


def refuse_accidental_overwrite(run: RunDirectory) -> None:
    """Stop before overwriting existing run artifacts unless allowed."""
    existing = [
        path
        for path in (run.cfg_path, run.neutral_path, run.ion_path)
        if path.exists()
    ]
    if not existing or OVERWRITE_EXISTING_RUN:
        return

    existing_text = "\n  ".join(str(path) for path in existing)
    raise SystemExit(
        "RUN_DIR already contains simulation output:\n"
        f"  {existing_text}\n\n"
        "Choose a new RUN_DIR, or set OVERWRITE_EXISTING_RUN = True if you "
        "really want to replace these files."
    )


def print_run_header(cfg, run: RunDirectory) -> None:
    """Print the resolved settings before the expensive work starts."""
    print("")
    print("Single-pulse iodine/helium simulation")
    print("=====================================")
    print(f"input preset:  {INPUT_PRESET}")
    print(f"run size:      {RUN_SIZE}")
    print(f"run directory: {run.root}")
    print(f"molecules:     {cfg.num_molecules}")
    print(f"seed:          {cfg.seed}")
    print(f"ion time:      {cfg.ion_simulation_time} ps")
    print("")


def print_run_footer(run: RunDirectory, *, neutral_steps: int, ion_steps: int) -> None:
    """Print output paths and the minimal code needed to load the run."""
    print("")
    print("Run complete")
    print("============")
    print(f"Config:  {run.cfg_path}")
    print(f"Neutral: {run.neutral_path} ({neutral_steps} stored time points)")
    print(f"Ion:     {run.ion_path} ({ion_steps} stored time points)")
    print("")
    print("Load this run later with:")
    print("")
    print("    from i2_helium_md.simulation.run_directory import RunDirectory")
    print(f"    run = RunDirectory({str(run.root)!r})")
    print("    cfg = run.load_cfg()")
    print("    neutral = run.load_neutral()")
    print("    ion = run.load_ion()")


if __name__ == "__main__":
    raise SystemExit(main())
