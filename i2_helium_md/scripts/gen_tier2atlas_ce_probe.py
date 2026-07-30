"""Atlas §3.5l — the (C) CE-channel + exit-strip probe at the h405 pins.

Runs the REGISTERED probe of `TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` §8/§9
(registration FROZEN 2026-07-29, user-approved bands): each cell is the
committed `g4fh405` config **verbatim** — corrected geometry, (v_c 5.5,
τ 4.4, E₀ 0.405), seed 20260729 (the finals seed: A-only is CRN-paired
against the committed h405 on every pre-existing stream; C-full/B-only
break per-molecule pairing at the source per design PC-5) — with only the
(C) surfaces switched on per cell:

* **cfull**  — mixture (frozen weights (0.30, 0.50, 0.20)) + strip (P1 box
  center), f_int,Q3 = 0.0985 (the pinned-E₀ bracket end the P3 forecast
  rides).
* **aonly**  — single channel at the retiring 2.7006 scalar budget + strip
  on — attributes (A) in isolation against the §3.5j bracket.
* **bonly**  — mixture on, strip off — attributes (B); expected to
  reproduce the §3.5k complementarity (fast branch parks suppressed).
* **cq3hi**  — the coupling arm: cfull with f_int,Q3 = 0.15 (the full-
  proportional bracket end) — the midHot/supp trade measured in-mixture.

Frozen inputs (design §8): weights (0.30, 0.50, 0.20); f 0.80;
σ (single, Q2, Q3) = (0.42, 0.31, 0.55) eV; E_single 0.53 eV; strip box
center a 2 / j₀ 1.75 / w_j 0.75; v_strip 9.9 Sourced; ε_carry 0.025
(mid of [0, 0.05]). `drag_state_coupling = "off"` everywhere (the STOPPED
§3.5i axis; the cfg-diff oracle catches stray activation).

**Oracles (§1.4, enforced before any MD and by --dry-run):**

1. **cfg-diff** — each cell's cfg diffs against the committed `g4fh405`
   `cfg.json` in exactly zero fields once the post-reference (C) fields are
   pinned back to the dataclass defaults (the sn-probe precedent); the real
   (C) values are asserted field-by-field per cell, and the corrected-
   geometry stamps are re-checked (the finals guard, reused).
2. **unit oracle** — E_ref reproduces 2.7006 (4 dp), the channel means
   reproduce 2.16/4.32 (2 dp) at f = 0.80, the strip forms reproduce
   P₀(10.16) = 1 (the committed §3.5j n = 1-ender exit speed is above
   v_strip) and G(j₀) = 1/2 exactly.

Scoring: `scripts/post_processing/tier2atlas_ce_probe_table.py` (CP-1..8).

Atlas stance: instrument runs — nothing here moves `finc1v725`; nothing
adopts at probe level (PC-2 sequencing governs a pass).

Usage::

    python scripts/gen_tier2atlas_ce_probe.py --dry-run   # oracles only
    python scripts/gen_tier2atlas_ce_probe.py             # all 4 cells
    python scripts/gen_tier2atlas_ce_probe.py cfull       # selected cells
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
from i2_helium_md.physics.exit_strip import (  # noqa: E402
    strip_depth_grading,
    strip_velocity_gate,
)
from i2_helium_md.sampling.ce_channels import (  # noqa: E402
    E_REF_PER_ION_EV,
    ce_channel_means_eV,
)
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS — the frozen registration inputs (design §8; do not tune)
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"

H405_RUN = finals_run_dir_name("h405")

FROZEN_WEIGHTS = (0.30, 0.50, 0.20)      # (w_single, w_Q2, w_Q3)
FROZEN_F = 0.80
FROZEN_SIGMAS = (0.42, 0.31, 0.55)       # (single, Q2, Q3) [eV]
FROZEN_E_SINGLE = 0.53                   # [eV]
F_INT_Q2 = 0.15                          # standing (h405 pin on the dominant channel)
F_INT_SINGLE = 0.15                      # default (weakly identified)
F_INT_Q3_PINNED = 0.0985                 # the pinned-E0 bracket end (P3 basis)
F_INT_Q3_FULL = 0.15                     # the full-proportional bracket end

STRIP_BOX_CENTER = dict(
    exit_strip_v_ref=9.9,                # Sourced
    exit_strip_exponent=2.0,             # a (box [1.5, 2.5])
    exit_strip_protect_j0=1.75,          # j0 (box [1.5, 2])
    exit_strip_width_rungs=0.75,         # w_j (box [0.5, 1])
    exit_strip_carry_eV=0.025,           # eps (box [0, 0.05] mid)
)

# The committed §3.5j n = 1-ender exit speed (unit-oracle anchor).
V_N1_ENDER_APS = 10.16

DEFAULT_CONCURRENCY = 4
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True

_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz",
                       "relaxation.npz", "detection.npz")

# The (C) config fields (the pre-registered cfg-diff surface: pinned back to
# the dataclass defaults for the diff, asserted explicitly per cell).
_CE_FIELDS = (
    "ce_channel_mode", "ce_channel_weights", "ce_fraction_f",
    "ce_channel_sigma_eV", "ce_single_ker_eV", "ce_q3_partner_mask",
    "ce_internal_energy_partition_fractions",
    "exit_strip_mode", "exit_strip_v_ref", "exit_strip_exponent",
    "exit_strip_protect_j0", "exit_strip_width_rungs", "exit_strip_carry_eV",
)


class CeCell(NamedTuple):
    """One §9 probe cell."""

    label: str
    channels: bool
    strip: bool
    f_int_q3: float
    role: str


RING: tuple[CeCell, ...] = (
    CeCell("cfull", True, True, F_INT_Q3_PINNED,
           "C-full: mixture + strip, f_int,Q3 pinned end"),
    CeCell("aonly", False, True, F_INT_Q3_PINNED,
           "A-only: retiring 2.7006 scalar budget + strip"),
    CeCell("bonly", True, False, F_INT_Q3_PINNED,
           "B-only: mixture, no strip (§3.5k complementarity check)"),
    CeCell("cq3hi", True, True, F_INT_Q3_FULL,
           "coupling arm: C-full at the full-proportional end"),
)
_SPEC_BY_LABEL = {c.label: c for c in RING}


def ce_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N, run_tag=f"tier2atlas_conf270_ce{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_ce_cell(cell: CeCell) -> SimConfig:
    """The committed h405 config with only this cell's (C) surfaces on."""
    cfg = build_cell(_FINALS_SPECS["h405"])
    extra: dict[str, object] = {}
    if cell.channels:
        extra.update(
            ce_channel_mode="sampled",
            ce_channel_weights=FROZEN_WEIGHTS,
            ce_fraction_f=FROZEN_F,
            ce_channel_sigma_eV=FROZEN_SIGMAS,
            ce_single_ker_eV=FROZEN_E_SINGLE,
            ce_q3_partner_mask=True,
            ce_internal_energy_partition_fractions=(
                F_INT_SINGLE, F_INT_Q2, cell.f_int_q3,
            ),
        )
    if cell.strip:
        extra.update(exit_strip_mode="depth_graded", **STRIP_BOX_CENTER)
    cfg = dataclasses.replace(cfg, **extra)
    cfg.validate()
    return cfg


def verify_ce_cfg_oracle(cfg: SimConfig, cell: CeCell) -> None:
    """cfg-diff oracle: exactly the (C) fields differ, nothing else.

    The committed h405 ``cfg.json`` predates the (C) fields, so the diff
    runs on a copy pinned back to the dataclass defaults (sn-probe
    precedent) — it must then be EMPTY — and the real values are asserted
    explicitly here (incl. the STOPPED s(n) axis staying off).
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / H405_RUN / "cfg.json"
    fields = SimConfig.__dataclass_fields__
    cfg_for_diff = dataclasses.replace(
        cfg, **{name: fields[name].default for name in _CE_FIELDS},
    )
    diff = cfg_diff_vs_reference(cfg_for_diff, ref_path,
                                 context=f"ce {cell.label}")
    if diff:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs committed h405 ((C) fields pinned) "
            f"is {sorted(diff)}; pre-registered exactly [] — a non-(C) pin "
            "drifted."
        )
    if cfg.drag_state_coupling != "off":
        raise AssertionError(
            f"[{cell.label}] STOPPED s(n) axis activated "
            f"({cfg.drag_state_coupling!r}) — standing exclusion violated."
        )
    # The real (C) values, field by field.
    if cell.channels:
        checks = [
            ("ce_channel_mode", cfg.ce_channel_mode, "sampled"),
            ("ce_channel_weights", tuple(cfg.ce_channel_weights),
             FROZEN_WEIGHTS),
            ("ce_fraction_f", float(cfg.ce_fraction_f), FROZEN_F),
            ("ce_channel_sigma_eV", tuple(cfg.ce_channel_sigma_eV),
             FROZEN_SIGMAS),
            ("ce_single_ker_eV", float(cfg.ce_single_ker_eV),
             FROZEN_E_SINGLE),
            ("ce_q3_partner_mask", bool(cfg.ce_q3_partner_mask), True),
            ("ce_internal_energy_partition_fractions",
             tuple(cfg.ce_internal_energy_partition_fractions),
             (F_INT_SINGLE, F_INT_Q2, cell.f_int_q3)),
        ]
    else:
        checks = [("ce_channel_mode", cfg.ce_channel_mode, "off")]
    if cell.strip:
        checks.append(("exit_strip_mode", cfg.exit_strip_mode, "depth_graded"))
        checks.extend(
            (name, float(getattr(cfg, name)), value)
            for name, value in STRIP_BOX_CENTER.items()
        )
    else:
        checks.append(("exit_strip_mode", cfg.exit_strip_mode, "off"))
    bad = [f"{name}: {got!r} != registered {want!r}"
           for name, got, want in checks if got != want]
    if bad:
        raise AssertionError(
            f"[{cell.label}] (C) field assertion failed: " + "; ".join(bad)
        )


def verify_ce_unit_oracle() -> None:
    """Unit oracle: E_ref / channel means / strip-form landmark values."""
    checks = [
        ("E_ref [eV]", round(E_REF_PER_ION_EV, 4), 2.7006),
    ]
    e_s, e_q2, e_q3 = ce_channel_means_eV(
        fraction_f=FROZEN_F, single_ker_eV=FROZEN_E_SINGLE,
    )
    checks += [("E_single", e_s, FROZEN_E_SINGLE),
               ("E_Q2 (2 dp)", round(e_q2, 2), 2.16),
               ("E_Q3 (2 dp)", round(e_q3, 2), 4.32)]
    p0 = strip_velocity_gate(
        V_N1_ENDER_APS,
        v_strip_aps=STRIP_BOX_CENTER["exit_strip_v_ref"],
        exponent=STRIP_BOX_CENTER["exit_strip_exponent"],
    )
    checks.append(("P0(10.16)", float(p0), 1.0))
    g_mid = strip_depth_grading(
        STRIP_BOX_CENTER["exit_strip_protect_j0"],
        protect_j0=STRIP_BOX_CENTER["exit_strip_protect_j0"],
        width_rungs=STRIP_BOX_CENTER["exit_strip_width_rungs"],
    )
    checks.append(("G(j0)", float(g_mid), 0.5))
    bad = [f"{name}: computed {got} vs registered {want}"
           for name, got, want in checks if got != want]
    if bad:
        raise AssertionError("unit oracle failed: " + "; ".join(bad))


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _verify_cell(cell: CeCell) -> SimConfig:
    cfg = build_ce_cell(cell)
    verify_corrected_geometry(cfg, _FINALS_SPECS["h405"])
    verify_ce_cfg_oracle(cfg, cell)
    return cfg


def _run_one(cell: CeCell) -> str:
    cfg = _verify_cell(cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / ce_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{cell.label}] channels={'on' if cell.channels else 'off'} "
          f"strip={'on' if cell.strip else 'off'} f_int,Q3={cell.f_int_q3} "
          f"({cell.role}) seed={cfg.seed}", flush=True)
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
    stripped = int(np.count_nonzero(relax.checkpoint.ce_strip_count > 0))
    print(f"[{cell.label}] done -> {run_dir} (n_detect_mean={nd:.2f}, "
          f"droplet_retained={int(np.count_nonzero(~dm))}/{dm.size}, "
          f"stripped_ions={stripped})", flush=True)
    return cell.label


def _run_one_by_label(label: str) -> str:
    return _run_one(_SPEC_BY_LABEL[label])


def _select(labels: list[str]) -> list[CeCell]:
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
    print(f"Atlas (C) probe: {len(cells)} cell(s), N={N}, seed={SEED} "
          f"(= the finals seed; A-only CRN-paired vs committed h405, PC-5), "
          f"h405 pins, drag_state_coupling=off, "
          f"concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_ce_unit_oracle()
    print("unit oracle OK (E_ref/means/strip landmarks reproduce).",
          flush=True)

    if args.dry_run:
        for cell in cells:
            _verify_cell(cell)
            print(f"[{cell.label}] cfg oracle OK -> "
                  f"{ce_run_dir_name(cell.label)}", flush=True)
        print("dry run complete: unit + cfg oracles pass.", flush=True)
        return 0

    concurrency = max(1, min(args.concurrency, len(cells)))
    if concurrency == 1 or len(cells) == 1:
        for cell in cells:
            _run_one(cell)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=concurrency) as pool:
        futures = {pool.submit(_run_one_by_label, c.label): c.label
                   for c in cells}
        for fut in as_completed(futures):
            label = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{label}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(cells)} cell(s) failed.", flush=True)
        return 1
    print(f"all {len(cells)} cells complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
