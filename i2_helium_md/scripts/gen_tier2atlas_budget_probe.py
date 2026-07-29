"""Atlas §3.5k — the budget-slope probe at the h405 pins (3 × N = 1000).

Runs the registered probe of `TIER2_SENSITIVITY_ATLAS_PLAN.md` §3.5k: each
cell is the committed `g4fh405` config **verbatim** — corrected geometry,
(v_c 5.5, τ 4.4), seed 20260729 (the finals seed, CRN-paired) — with only
the per-fragment Coulomb budget moved to a measured gas-phase CE channel
via `E_coulomb_scale` (fixed birth geometry R₀ = 2.666 Å), the
`coulomb_available_eV` stamp moved with it, and E₀ handled per cell:

* **bud226k** (2.26 eV/frag, E₀ pinned 0.405)  — kinematic line, down,
* **bud411k** (4.11 eV/frag, E₀ pinned 0.405)  — kinematic line, up,
* **bud411**  (4.11 eV/frag, f_int 0.15 kept)  — the (B) wiring contrast.

The `k` cells pin E₀ by rescaling `internal_energy_partition_fraction`
(f_int = 0.405/budget), isolating the kinematic KER → KE₁ transfer from
the E₀-gate action; bud411 keeps the standing f_int so E₀ scales with the
budget (the (B) proportional wiring, E₀ = 0.6165).

**Oracles (§1.4, enforced before any MD and by --dry-run):**

1. **cfg-diff** — each cell's cfg diffs against the committed `g4fh405`
   `cfg.json` in exactly the registered field set and nothing else
   (all three fields predate the reference, so the diff runs directly).
2. **kinematic unit oracle** — `scale·14.39964548/(2·R₀)` equals the
   registered budget and `f_int·stamp` equals the registered E₀, both
   to 1e-9; R₀ asserted at 2.666 (the budget moves through the scale,
   never through the geometry).

Atlas stance: instrument runs — nothing here moves `finc1v725`; the
BP-KILL consequence (honest-residual branch at S_k < 0.2) is
pre-registered in plan §3.5k.

Usage::

    python scripts/gen_tier2atlas_budget_probe.py --dry-run   # oracles only
    python scripts/gen_tier2atlas_budget_probe.py             # all 3 cells
    python scripts/gen_tier2atlas_budget_probe.py bud411k     # selected cells
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path
import sys
from typing import NamedTuple, Optional

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    _SPEC_BY_LABEL as _FINALS_SPECS,
    N,
    SEED,
    build_cell,
    finals_run_dir_name,
    verify_corrected_geometry,
)
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"

# The committed h405 run this probe is CRN-paired against (its cfg.json is
# the cfg-diff oracle reference; its seed is inherited via build_cell).
H405_RUN = finals_run_dir_name("h405")

# The Coulomb constant as hard-coded in physics/interactions.py (the single
# package source; restated here as an oracle literal so the registered budget
# arithmetic is checkable without running the potential).
COULOMB_EV_ANGSTROM = 14.39964548
R0_PRODUCTION_ANGSTROM = 2.666
E_FRAG_REF_EV = COULOMB_EV_ANGSTROM / R0_PRODUCTION_ANGSTROM / 2.0  # 2.7006

E0_PIN_EV = 0.405        # the h405 onset (f_int 0.15 x stamp 2.7), pinned in k-cells
F_INT_STANDING = 0.15    # the standing fraction, kept in the wiring cell

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz",
                       "relaxation.npz", "detection.npz")

_KINEMATIC_FIELDS = ("E_coulomb_scale", "coulomb_available_eV")
_PIN_FIELDS = _KINEMATIC_FIELDS + ("internal_energy_partition_fraction",)


class BudgetCell(NamedTuple):
    """One §3.5k probe cell."""

    label: str
    budget_eV: float          # registered per-fragment kinematic budget
    pin_E0: bool              # True: f_int rescaled to hold E0 = 0.405
    role: str

    @property
    def scale(self) -> float:
        return self.budget_eV / E_FRAG_REF_EV

    @property
    def f_int(self) -> float:
        return (E0_PIN_EV / self.budget_eV) if self.pin_E0 else F_INT_STANDING

    @property
    def E0_eV(self) -> float:
        return self.f_int * self.budget_eV

    @property
    def expected_diff(self) -> tuple[str, ...]:
        return _PIN_FIELDS if self.pin_E0 else _KINEMATIC_FIELDS


RING: tuple[BudgetCell, ...] = (
    BudgetCell("bud226k", 2.26, True, "kinematic line, down (measured CE channel)"),
    BudgetCell("bud411k", 4.11, True, "kinematic line, up (measured fast channel)"),
    BudgetCell("bud411", 4.11, False, "(B) proportional wiring contrast"),
)
_SPEC_BY_LABEL = {c.label: c for c in RING}


def budget_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N, run_tag=f"tier2atlas_conf270_{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_budget_cell(cell: BudgetCell) -> SimConfig:
    """The committed h405 config with only the registered budget fields moved."""
    cfg = build_cell(_FINALS_SPECS["h405"])
    replacements: dict[str, float] = {
        "E_coulomb_scale": cell.scale,
        "coulomb_available_eV": cell.budget_eV,
    }
    if cell.pin_E0:
        replacements["internal_energy_partition_fraction"] = cell.f_int
    cfg = dataclasses.replace(cfg, **replacements)
    cfg.validate()
    return cfg


def verify_budget_cfg_oracle(cfg: SimConfig, cell: BudgetCell) -> None:
    """§3.5k cfg-diff oracle: exactly the registered field set, nothing else."""
    ref_path = PROJECT_ROOT / "data" / "runs" / H405_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref_path, context=f"budget {cell.label}")
    if sorted(diff) != sorted(cell.expected_diff):
        raise AssertionError(
            f"[{cell.label}] cfg diff vs committed h405 is {sorted(diff)}; "
            f"pre-registered exactly {sorted(cell.expected_diff)}."
        )


def verify_budget_unit_oracle(cfg: SimConfig, cell: BudgetCell) -> None:
    """§3.5k kinematic unit oracle: budget and E0 arithmetic to 1e-9."""
    if float(cfg.R0_GS_angstrom) != R0_PRODUCTION_ANGSTROM:
        raise AssertionError(
            f"[{cell.label}] R0_GS {cfg.R0_GS_angstrom!r} != "
            f"{R0_PRODUCTION_ANGSTROM} — the budget must move through "
            "E_coulomb_scale, never the birth geometry."
        )
    per_frag = (float(cfg.E_coulomb_scale) * COULOMB_EV_ANGSTROM
                / float(cfg.R0_GS_angstrom) / 2.0)
    if abs(per_frag - cell.budget_eV) > 1e-9:
        raise AssertionError(
            f"[{cell.label}] kinematic budget {per_frag!r} != registered "
            f"{cell.budget_eV!r} (scale {cfg.E_coulomb_scale!r})."
        )
    if abs(float(cfg.coulomb_available_eV) - cell.budget_eV) > 1e-9:
        raise AssertionError(
            f"[{cell.label}] stamp {cfg.coulomb_available_eV!r} != registered "
            f"budget {cell.budget_eV!r}."
        )
    onset = (float(cfg.internal_energy_partition_fraction)
             * float(cfg.coulomb_available_eV))
    if abs(onset - cell.E0_eV) > 1e-9:
        raise AssertionError(
            f"[{cell.label}] E0 = f_int*stamp = {onset!r} != registered "
            f"{cell.E0_eV!r}."
        )


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _verify_cell(cell: BudgetCell) -> SimConfig:
    cfg = build_budget_cell(cell)
    verify_corrected_geometry(cfg, _FINALS_SPECS["h405"])
    verify_budget_cfg_oracle(cfg, cell)
    verify_budget_unit_oracle(cfg, cell)
    return cfg


def _run_one(cell: BudgetCell) -> str:
    cfg = _verify_cell(cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / budget_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] budget={cell.budget_eV} eV/frag "
          f"scale={cell.scale:.6f} f_int={cell.f_int:.6f} "
          f"E0={cell.E0_eV:.4f} ({cell.role}) seed={cfg.seed}", flush=True)
    print(f"[{cell.label}] neutral ...", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{cell.label}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(f"[{cell.label}] relaxation (E2) ...", flush=True)
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")
    print(f"[{cell.label}] detection ...", flush=True)
    detect = run_detection_stage(
        relax.checkpoint, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    nd = float(detect.n_detected[dm].mean()) if dm.any() else float("nan")
    print(f"[{cell.label}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size})",
          flush=True)
    return cell.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[BudgetCell]:
    if not labels:
        return list(RING)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in RING]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas 3.5k budget probe: {len(cells)} cell(s), N={N}, seed={SEED} "
          f"(= the finals seed, CRN-paired), h405 pins, "
          f"E_frag(scale=1)={E_FRAG_REF_EV:.6f} eV, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    if args.dry_run:
        for cell in cells:
            _verify_cell(cell)
            print(f"[{cell.label}] cfg + kinematic oracles OK "
                  f"budget={cell.budget_eV} scale={cell.scale:.6f} "
                  f"f_int={cell.f_int:.6f} E0={cell.E0_eV:.4f} -> "
                  f"{budget_run_dir_name(cell.label)}", flush=True)
        print("dry run complete: all cfg + kinematic oracles pass.", flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(cells) == 1:
        for cell in cells:
            _run_one(cell)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=concurrency) as pool:
        futures = {pool.submit(_run_one_by_label, cell.label): cell.label
                   for cell in cells}
        for future in as_completed(futures):
            future.result()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
