"""Atlas §3.5i — the s(n) drag-state-coupling probe at the h405 pins (3 × N = 1000).

Runs the registered probe of `TIER2_DRAG_STATE_COUPLING_DESIGN.md` §8: each
cell is the committed `g4fh405` config **verbatim** — corrected geometry,
(v_c 5.5, τ 4.4, E₀ 0.405), seed 20260729 (the finals seed, CRN-paired) —
with only the drag-state-coupling surface switched on
(`drag_state_coupling = "shell_area"`, R_core = 3.2 Å) and ρ_shell swept
across its **Bounded physical range** (the parameter's uncertainty interval,
not a fit axis):

* **sa22** (ρ = 0.0218, bulk)    — strong-coupling end,
* **sa30** (ρ = 0.030, prior)    — the geometric prior,
* **sa44** (ρ = 0.0436, 2×bulk)  — weak-coupling end.

(Direction: lower shell density → fluffier dressed object → larger
dressed/bare contrast → stronger coupling.)

**Oracles (§1.4, enforced before any MD and by --dry-run):**

1. **cfg-diff** — each cell's cfg diffs against the committed `g4fh405`
   `cfg.json` in exactly the three new coupling fields and nothing else.
   The reference predates the fields, so (the finals-generator
   post-reference precedent) the diff runs on a copy pinned back to the
   dataclass defaults and the real values are asserted field-by-field;
   the corrected-geometry stamps are re-checked (the finals guard, reused).
2. **unit oracle** — `s(n_ref) = 1` exactly, and the design §8 registered
   s-table (post-erratum 2026-07-29: exact closure arithmetic) reproduces
   to 3 decimals per cell, the §3 prior table additionally at sa30.

Atlas stance: instrument runs — nothing here moves `finc1v725`; nothing
adopts at probe level (design §8 success shape → pooled-battery candidate).

Usage::

    python scripts/gen_tier2atlas_sn_probe.py --dry-run   # oracles only
    python scripts/gen_tier2atlas_sn_probe.py             # all 3 cells
    python scripts/gen_tier2atlas_sn_probe.py sa30        # selected cells
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
from i2_helium_md.physics.state_coupling import (  # noqa: E402
    derive_n_ref_amu,
    effective_radius_angstrom,
    shell_area_state_factor,
)
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

R_CORE_ANGSTROM = 3.2   # design §8: fixed at the Bounded prior for all cells

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz",
                       "relaxation.npz", "detection.npz")

# The three new config fields (the pre-registered cfg-diff surface).
_COUPLING_FIELDS = ("drag_state_coupling", "state_coupling_R_core_angstrom",
                    "state_coupling_rho_shell_per_A3")


class SnCell(NamedTuple):
    """One §3.5i probe cell."""

    label: str
    rho_shell_per_A3: float
    role: str
    # Registered unit-oracle values (design §8 table, post-erratum):
    # (R_eff(19) [A] to 2 dp, s(1) to 3 dp, s(0) to 3 dp).
    reff19_2dp: float
    s1_3dp: float
    s0_3dp: float


RING: tuple[SnCell, ...] = (
    SnCell("sa22", 0.0218, "strong-coupling end (bulk density)",
           6.22, 0.321, 0.265),
    SnCell("sa30", 0.030, "the geometric prior",
           5.69, 0.366, 0.317),
    SnCell("sa44", 0.0436, "weak-coupling end (2x bulk)",
           5.15, 0.428, 0.386),
)
_SPEC_BY_LABEL = {c.label: c for c in RING}

# Design §3 prior table (post-erratum; asserted at sa30, the prior cell).
_PRIOR_S_TABLE_3DP = {21: 1.057, 19: 1.000, 14: 0.850, 8: 0.650,
                      5: 0.538, 2: 0.412, 1: 0.366, 0: 0.317}


def sn_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N, run_tag=f"tier2atlas_conf270_{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_sn_cell(cell: SnCell) -> SimConfig:
    """The committed h405 config with only the state coupling switched on."""
    cfg = build_cell(_FINALS_SPECS["h405"])
    cfg = dataclasses.replace(
        cfg,
        drag_state_coupling="shell_area",
        state_coupling_R_core_angstrom=R_CORE_ANGSTROM,
        state_coupling_rho_shell_per_A3=cell.rho_shell_per_A3,
    )
    cfg.validate()
    return cfg


def verify_sn_cfg_oracle(cfg: SimConfig, cell: SnCell) -> None:
    """§3.5i cfg-diff oracle: exactly the three coupling fields, nothing else.

    The committed h405 ``cfg.json`` predates the coupling fields, and the
    shared diff helper (correctly) refuses to certify a post-reference field
    off its default. The diff therefore runs on a copy pinned back to the
    dataclass defaults (the finals-generator ``droplet_size_sampler_mode``
    precedent) — it must then be EMPTY — and the real coupling values are
    asserted explicitly here.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / H405_RUN / "cfg.json"
    fields = SimConfig.__dataclass_fields__
    cfg_for_diff = dataclasses.replace(
        cfg, **{name: fields[name].default for name in _COUPLING_FIELDS},
    )
    diff = cfg_diff_vs_reference(cfg_for_diff, ref_path,
                                 context=f"sn {cell.label}")
    if diff:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs committed h405 (coupling fields "
            f"pinned) is {sorted(diff)}; pre-registered exactly [] — a "
            f"non-coupling pin drifted."
        )
    if cfg.drag_state_coupling != "shell_area":
        raise AssertionError(f"[{cell.label}] coupling enum not shell_area.")
    if float(cfg.state_coupling_R_core_angstrom) != R_CORE_ANGSTROM:
        raise AssertionError(
            f"[{cell.label}] R_core {cfg.state_coupling_R_core_angstrom!r} "
            f"!= registered {R_CORE_ANGSTROM!r}."
        )
    if float(cfg.state_coupling_rho_shell_per_A3) != cell.rho_shell_per_A3:
        raise AssertionError(
            f"[{cell.label}] rho_shell "
            f"{cfg.state_coupling_rho_shell_per_A3!r} != cell value "
            f"{cell.rho_shell_per_A3!r}."
        )


def verify_sn_unit_oracle(cfg: SimConfig, cell: SnCell) -> None:
    """§8 unit oracle: s(n_ref) = 1 exactly + the registered s-table to 3 dp."""
    n_ref = derive_n_ref_amu(cfg.drag_coefficients.extraction_mass_amu)
    if n_ref != 19:
        raise AssertionError(
            f"[{cell.label}] n_ref derived {n_ref} != 19 (the bundle stamp "
            "changed?)"
        )
    kwargs = dict(
        R_core_angstrom=cfg.state_coupling_R_core_angstrom,
        rho_shell_per_A3=cfg.state_coupling_rho_shell_per_A3,
    )
    s_ref = float(shell_area_state_factor(float(n_ref), n_ref=n_ref, **kwargs))
    if s_ref != 1.0:
        raise AssertionError(f"[{cell.label}] s(n_ref) = {s_ref!r} != 1.0 exact.")
    reff19 = float(effective_radius_angstrom(float(n_ref), **kwargs))
    checks = [("R_eff(19)", round(reff19, 2), cell.reff19_2dp)]
    for n_val, want in ((1, cell.s1_3dp), (0, cell.s0_3dp)):
        got = float(shell_area_state_factor(float(n_val), n_ref=n_ref, **kwargs))
        checks.append((f"s({n_val})", round(got, 3), want))
    if cell.label == "sa30":
        for n_val, want in _PRIOR_S_TABLE_3DP.items():
            got = float(
                shell_area_state_factor(float(n_val), n_ref=n_ref, **kwargs)
            )
            checks.append((f"prior s({n_val})", round(got, 3), want))
    bad = [f"{name}: computed {got} vs registered {want}"
           for name, got, want in checks if got != want]
    if bad:
        raise AssertionError(
            f"[{cell.label}] unit oracle vs the registered s-table failed: "
            + "; ".join(bad)
        )


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _verify_cell(cell: SnCell) -> SimConfig:
    cfg = build_sn_cell(cell)
    verify_corrected_geometry(cfg, _FINALS_SPECS["h405"])
    verify_sn_cfg_oracle(cfg, cell)
    verify_sn_unit_oracle(cfg, cell)
    return cfg


def _run_one(cell: SnCell) -> str:
    cfg = _verify_cell(cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / sn_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] rho_shell={cell.rho_shell_per_A3} "
          f"R_core={R_CORE_ANGSTROM} ({cell.role}) seed={cfg.seed}",
          flush=True)
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


def _select(labels: list[str]) -> list[SnCell]:
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
    print(f"Atlas 3.5i s(n) probe: {len(cells)} cell(s), N={N}, seed={SEED} "
          f"(= the finals seed, CRN-paired), h405 pins, "
          f"R_core={R_CORE_ANGSTROM}, concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    if args.dry_run:
        for cell in cells:
            _verify_cell(cell)
            print(f"[{cell.label}] cfg + unit oracles OK "
                  f"rho_shell={cell.rho_shell_per_A3} -> "
                  f"{sn_run_dir_name(cell.label)}", flush=True)
        print("dry run complete: all cfg + unit oracles pass.", flush=True)
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
