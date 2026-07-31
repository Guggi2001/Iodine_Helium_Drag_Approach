"""Atlas free-form linear sweep — the CRN MD confirmation ring vs h405.

Runs the `TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN` §6.1 ring: **8 cells ×
N = 500, one fresh shared seed 20260731 (CRN across every cell including
the capped incumbent partner)** at the G2-adopted corrected geometry —
7 `pure_linear` cells from the arm-1 twin basin (sub-plateau core
a 27.5–35) plus the h405 capped partner (v_c 5.5 / τ 4.4 / E₀ 0.405,
the g4finals pins verbatim).

Cells (plan §6.1 table; per-cell knobs = (form, a | v_c, E_bind, τ, E₀),
everything else the standing finc1v725 pins at the corrected geometry):

* lr1  a 27.5 / eb0482 / τ 4.8 / E₀ 0.35 — KE₁ optimum, Φ 0.64
* lr2  a 27.5 / eb1168 / τ 4.8 / E₀ 0.35 — well axis at the optimum
* lr3  a 30.0 / eb0482 / τ 4.8 / E₀ 0.36 — sub-plateau W₁-best
* lr4  a 30.0 / eb1168 / τ 4.8 / E₀ 0.36 — well axis, second point
* lr5  a 32.5 / eb0482 / τ 4.8 / E₀ 0.36 — mid-core
* lr6  a 35.0 / eb0482 / τ 4.8 / E₀ 0.37 — W₁ optimum (a*), Φ 0.81
* lr7  a 27.5 / eb0482 / τ 4.8 / E₀ 0.37 — E₀ axis at fixed chord
* h405p  capped_cubic v_c 5.5 / eb1168 / τ 4.4 / E₀ 0.405 — CRN partner

**LR-P1 oracle (enforced before any MD and by --dry-run):** the 7 frozen
lin twin rows must reproduce **string-exact** from the committed
`atlas_linsweep.csv`, the h405 partner anchors on the committed
`atlas_ke1_authority.csv` pooled row (gate-on-committed-artifacts rule),
and every cell's cfg must diff against the standing battery reference in
exactly its pre-registered key set.

Provenance posture (the first free-form MD coefficients of the program):
the `pure_linear` bundle is stamped `extraction_method="free_form"` with
**no** `effective_binding_energy_I_ion_eV` (a free-form parameter was
never jointly calibrated with any binding), so §6.5.1 warn-fires on
every lin cell under the documented `allow_unvalidated_binding_pairing`
hatch; §6.5 mass pairing warn-fires as in every biphasic run.

Retained policy: `exclude_all_coupled` headline (Arm A); the scorer
(`scripts/post_processing/tier2atlas_linring_table.py`) re-scores every
cell under the §3.5b Arm-B marginal-injection end per plan §6.2 item 3.

Atlas stance: instrument runs — nothing here moves `finc1v725`/h405.

Usage::

    python scripts/gen_tier2atlas_linring.py --dry-run   # oracles + cfg guards
    python scripts/gen_tier2atlas_linring.py             # all 8 cells
    python scripts/gen_tier2atlas_linring.py lr1 h405p   # selected cells
"""

from __future__ import annotations

import argparse
import csv
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


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 500                         # fragments per cell (plan §6.1, N adjudicated)
SEED = 20260731                 # ONE fresh seed shared by every cell (CRN)
ION_TIME_PS = 30.0
DT_ION_PS = 0.01

BUDGET_EV = 2.70
S_EFF = 8.0
F_RET = 0.1
LAMBDA0_PER_PS = 0.9
COOLING_GATE = "density_scaled"
P_TAIL = -1.0                   # h405 partner only (capped_cubic tail)

# The G2-adopted corrected geometry (G0 spec — no new numbers).
SIZE_PRIOR = "legacy"
SIZE_SAMPLER_MODE = "raw"
BIRTH_LAW = "boltzmann"
BIRTH_MARGIN_ANGSTROM = 0.0
PARENT_WELL_K = 313.2

# Arm-A headline policy; the scorer adds the Arm-B injection end (§6.2 item 3).
RETAINED_POLICY = "exclude_all_coupled"
SHED_CONVENTION = "co_moving"
SHELL_MODEL = "density_tied"
PARTITION_LAW = "sigma_proportional"
R0_GS_ANGSTROM = 2.666
E_COULOMB_SCALE = 1.0
RELAXATION_TIME_PS = 8000.0
RELAXATION_FORCES = "coulomb"
RELAXATION_DISSIPATION = "landau_gated_drag"
V_LIMIT_M_PER_S = 58.0
DETECTION_TIME_PS = 8.53e6

# Exact E_bind values behind the well tags (§3.5c convention; eb1168 is the
# shared_pure_cubic bundle stamp value — set explicitly on lin cells because
# the free_form bundle carries no stamp to inherit it from).
EBIND_EV = {"eb0482": 0.0482, "eb1168": 0.11675778353879479}

# Reference run for the cfg-diff guard: a standing-point battery member.
REFERENCE_RUN = "9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s1"

# Keys every cell must differ in vs the reference (N + seed + the corrected
# geometry + the interim retained policy + this ring's τ/E₀/chord moves).
# `droplet_size_sampler_mode` handling as in gen_tier2atlas_g3ring.py.
_BASE_DIFF_KEYS = frozenset({
    "num_molecules", "seed",
    "droplet_size_prior",
    "birth_position_law", "initial_position_margin_angstrom",
    "binding_energy_molecule_K",
    "detection_droplet_retained_policy",
    "internal_energy_cooling_tau_ps",
    "internal_energy_partition_fraction",
    "drag_coefficients",
})
_LIN_DIFF_KEYS = frozenset({
    "drag_form", "allow_unvalidated_binding_pairing",
})

# Standing finc1v725 knob values (context for the diff-key comments above:
# every ring τ/E₀ differs from these, so both knobs sit in _BASE_DIFF_KEYS).
STANDING_TAU_PS = 3.2
STANDING_E0_EV = 0.27

# The committed twin-scan record the LR-P1 oracle re-reads (7 lin cells).
LINSWEEP_CSV = Path("data/runs/h2b_forward_model/atlas_linsweep.csv")
# The committed Step-0 authority record anchoring the h405 partner.
KE1AUTH_CSV = Path("data/runs/h2b_forward_model/atlas_ke1_authority.csv")

DEFAULT_CONCURRENCY = 3
SKIP_COMPLETED_RUNS = True
OVERWRITE_EXISTING_RUN = True


class LinRingCell(NamedTuple):
    """One §6.1 ring cell (lin: ``a`` set, ``v_c`` None; partner: inverse)."""

    label: str
    a: Optional[float]          # amu/ps (pure_linear cells)
    v_c: Optional[float]        # A/ps   (capped partner)
    eb_tag: str
    tau_ps: float
    e0_eV: float
    role: str


RING_MATRIX: tuple[LinRingCell, ...] = (
    LinRingCell("lr1", 27.5, None, "eb0482", 4.8, 0.35, "KE1 optimum, Phi 0.64"),
    LinRingCell("lr2", 27.5, None, "eb1168", 4.8, 0.35, "well axis at the optimum"),
    LinRingCell("lr3", 30.0, None, "eb0482", 4.8, 0.36, "sub-plateau W1-best"),
    LinRingCell("lr4", 30.0, None, "eb1168", 4.8, 0.36, "well axis, second point"),
    LinRingCell("lr5", 32.5, None, "eb0482", 4.8, 0.36, "mid-core"),
    LinRingCell("lr6", 35.0, None, "eb0482", 4.8, 0.37, "W1 optimum (a*), Phi 0.81"),
    LinRingCell("lr7", 27.5, None, "eb0482", 4.8, 0.37, "E0 axis at fixed chord"),
    LinRingCell("h405p", None, 5.5, "eb1168", 4.4, 0.405, "CRN incumbent partner"),
)
_SPEC_BY_LABEL = {c.label: c for c in RING_MATRIX}

# The frozen lin twin rows (plan §3.1 basin; committed atlas_linsweep.csv
# strings, compared EXACTLY): label -> (trapped_frac, suppressed_frac,
# nbar_det, n1_solv, w1_solv, n1_ke_eV, ke2_eV, deepke, gate).
TWIN_ROWS: dict[str, tuple[str, ...]] = {
    "lr1": ("0.0013", "0.0767", "4.488", "0.1988", "0.6729",
            "0.8922", "0.827", "3.375", "1"),
    "lr2": ("0.0031", "0.074", "4.507", "0.1965", "0.6897",
            "0.8546", "0.7885", "2.8414", "1"),
    "lr3": ("0.003", "0.1067", "4.423", "0.2083", "0.5871",
            "0.817", "0.7491", "2.7296", "1"),
    "lr4": ("0.0076", "0.1036", "4.413", "0.2067", "0.624",
            "0.7801", "0.7114", "2.2312", "1"),
    "lr5": ("0.007", "0.0954", "4.686", "0.1958", "0.6093",
            "0.7725", "0.7019", "2.3022", "1"),
    "lr6": ("0.0128", "0.1242", "4.625", "0.2008", "0.5857",
            "0.7027", "0.6322", "1.8634", "1"),
    "lr7": ("0.0013", "0.1685", "3.838", "0.2397", "0.5836",
            "0.8366", "0.7733", "3.1395", "0"),
}
_TWIN_ROW_COLS = ("trapped_frac", "suppressed_frac", "nbar_det", "n1_solv",
                  "w1_solv", "n1_ke_eV", "ke2_eV", "deepke", "gate")

# The frozen h405 partner anchor (committed atlas_ke1_authority.csv pooled
# row, string-exact): md_KE1, md_KE2, twin_KE1, twin_KE2, twin_n1_solv,
# twin_nbar, twin_w1.
H405_ANCHOR: tuple[str, ...] = (
    "0.6406", "0.5493", "0.6236", "0.5407", "0.2083", "4.484", "0.7025",
)
_H405_ANCHOR_COLS = ("md_KE1", "md_KE2", "twin_KE1", "twin_KE2",
                     "twin_n1_solv", "twin_nbar", "twin_w1")


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_md_confirmation import rq4graded_rungs_eV  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import build_biphasic_cfg, cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

_REQUIRED_ARTIFACTS: tuple[str, ...] = (
    "cfg.json", "neutral.npz", "ion.npz", "relaxation.npz", "detection.npz",
)


def verify_twin_preregistration() -> None:
    """LR-P1 first half: frozen rows reproduce string-exact from the committed
    artifacts. Fails loud — a mismatch means the frozen table or a committed
    CSV drifted, and the ring must not launch."""
    lin_path = PROJECT_ROOT / LINSWEEP_CSV
    with open(lin_path, newline="", encoding="utf-8") as fh:
        idx = {
            (r["arm"], r["a"], r["E_bind_tag"], r["tau_ps"], r["E0_eV"]): r
            for r in csv.DictReader(fh)
        }
    for cell in RING_MATRIX:
        if cell.a is None:
            continue
        key = ("lin", f"{cell.a:.1f}", cell.eb_tag, str(cell.tau_ps),
               str(cell.e0_eV))
        row = idx.get(key)
        if row is None:
            raise AssertionError(
                f"[{cell.label}] LR-P1 FAILED: no committed linsweep row at "
                f"{key}."
            )
        got = tuple(row[col] for col in _TWIN_ROW_COLS)
        if got != TWIN_ROWS[cell.label]:
            raise AssertionError(
                f"[{cell.label}] LR-P1 FAILED: committed linsweep row {got} "
                f"!= frozen pre-registration {TWIN_ROWS[cell.label]}."
            )

    auth_path = PROJECT_ROOT / KE1AUTH_CSV
    anchor = None
    with open(auth_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["group"] == "h405bat" and r["label"] == "pooled":
                anchor = tuple(r[col] for col in _H405_ANCHOR_COLS)
                break
    if anchor is None:
        raise AssertionError(
            "[h405p] LR-P1 FAILED: no h405bat pooled row in the committed "
            "authority CSV."
        )
    if anchor != H405_ANCHOR:
        raise AssertionError(
            f"[h405p] LR-P1 FAILED: committed authority row {anchor} != "
            f"frozen pre-registration {H405_ANCHOR}."
        )
    print(f"LR-P1 twin oracle PASSED: all {len(TWIN_ROWS)} frozen lin rows + "
          "the h405 anchor reproduce string-exact from the committed "
          "artifacts.", flush=True)


def ring_run_dir_name(label: str) -> str:
    """Run-dir basename in the atlas namespace (no ``_tier2_`` substring)."""
    name = run_dir_name(
        CASE, VARIANT, N,
        run_tag=f"tier2atlas_conf{int(round(BUDGET_EV * 100))}_linr{label}",
    )
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name


def build_cell(cell: LinRingCell) -> SimConfig:
    """One ring cell's validated config: corrected geometry + this cell's
    (form, a | v_c, E_bind, tau, E0), the standing pins otherwise."""
    if cell.a is not None:
        drag_form = "pure_linear"
        drag_overrides: dict[str, float] = {"a": cell.a}
    else:
        drag_form = "capped_cubic"
        drag_overrides = {"v_c": cell.v_c, "p_tail": P_TAIL}
    cfg = build_biphasic_cfg(
        CASE, VARIANT,
        num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=SEED, lambda0_per_ps=LAMBDA0_PER_PS,
        f_int=cell.e0_eV / BUDGET_EV, f_ret=F_RET, tau_ps=cell.tau_ps,
        evap_rrk_dof=S_EFF, coulomb_available_eV=BUDGET_EV,
        relaxation_time_ps=RELAXATION_TIME_PS,
        relaxation_forces=RELAXATION_FORCES,
        cooling_spatial_gate=COOLING_GATE, drag_form=drag_form,
        drag_coefficient_overrides=drag_overrides,
        dissociation_ladder="tabulated",
        tabulated_ladder_rungs_eV=rq4graded_rungs_eV(),
        R0_GS_angstrom=R0_GS_ANGSTROM, E_coulomb_scale=E_COULOMB_SCALE,
        single_initial_position=False,
        detection_time_ps=DETECTION_TIME_PS,
        birth_position_law=BIRTH_LAW,
        initial_position_margin_angstrom=BIRTH_MARGIN_ANGSTROM,
        detection_droplet_retained_policy=RETAINED_POLICY,
        evaporation_shed_convention=SHED_CONVENTION,
        initial_shell_model=SHELL_MODEL,
        internal_energy_partition_law=PARTITION_LAW,
        droplet_size_prior=SIZE_PRIOR,
        use_single_droplet_size=False,
        droplet_size_sampler_mode=SIZE_SAMPLER_MODE,
        binding_energy_molecule_K=PARENT_WELL_K,
    )
    extra: dict[str, object] = dict(
        relaxation_dissipation=RELAXATION_DISSIPATION,
        v_limit_m_per_s=V_LIMIT_M_PER_S,
    )
    if cell.a is not None:
        # Free-form provenance: a was twin-selected, never extracted, and
        # never jointly calibrated with any binding — the stamp goes to
        # None and the documented §6.5.1 hatch is set (plan §6.1; the
        # coefficients-only replace keeps the pure_linear form + a).
        extra["drag_coefficients"] = dataclasses.replace(
            cfg.drag_coefficients,
            extraction_method="free_form",
            effective_binding_energy_I_ion_eV=None,
        )
        extra["binding_energy_I_ion_eV"] = EBIND_EV[cell.eb_tag]
        extra["allow_unvalidated_binding_pairing"] = True
    cfg = dataclasses.replace(cfg, **extra)
    cfg.validate()
    return cfg


def verify_corrected_geometry(cfg: SimConfig, cell: LinRingCell) -> None:
    """The G0/G2 corrected-geometry stamps + this ring's form/provenance
    posture, checked field by field."""
    if (cfg.droplet_size_prior != "legacy" or cfg.use_single_droplet_size
            or cfg.droplet_size_sampler_mode != "raw"):
        raise AssertionError(
            f"[{cell.label}] corrected size law broken: prior="
            f"{cfg.droplet_size_prior!r}, single={cfg.use_single_droplet_size}, "
            f"mode={cfg.droplet_size_sampler_mode!r}."
        )
    if (cfg.birth_position_law != "boltzmann"
            or cfg.initial_position_margin_angstrom != 0.0
            or cfg.binding_energy_molecule_K != PARENT_WELL_K):
        raise AssertionError(
            f"[{cell.label}] corrected birth law broken: law="
            f"{cfg.birth_position_law!r}, margin="
            f"{cfg.initial_position_margin_angstrom}, well="
            f"{cfg.binding_energy_molecule_K}."
        )
    if cfg.detection_droplet_retained_policy != RETAINED_POLICY:
        raise AssertionError(f"[{cell.label}] retained policy drifted.")
    if cell.a is not None:
        coeffs = cfg.drag_coefficients
        if (cfg.drag_form != "pure_linear" or coeffs.form != "pure_linear"
                or float(coeffs.coefficients["a"]) != cell.a):
            raise AssertionError(
                f"[{cell.label}] lin chord broken: form={cfg.drag_form!r}, "
                f"coefficients={dict(coeffs.coefficients)!r}."
            )
        if (coeffs.extraction_method != "free_form"
                or coeffs.effective_binding_energy_I_ion_eV is not None):
            raise AssertionError(
                f"[{cell.label}] free-form provenance posture broken: "
                f"method={coeffs.extraction_method!r}, stamp="
                f"{coeffs.effective_binding_energy_I_ion_eV!r}."
            )
        if cfg.binding_energy_I_ion_eV != EBIND_EV[cell.eb_tag]:
            raise AssertionError(f"[{cell.label}] well not stamped.")
        if not cfg.allow_unvalidated_binding_pairing:
            raise AssertionError(f"[{cell.label}] §6.5.1 hatch not set.")
    else:
        coeffs = cfg.drag_coefficients
        if (cfg.drag_form != "capped_cubic"
                or float(coeffs.coefficients["v_c"]) != cell.v_c
                or float(coeffs.coefficients["p_tail"]) != P_TAIL):
            raise AssertionError(
                f"[{cell.label}] h405 partner chord broken: form="
                f"{cfg.drag_form!r}, coefficients="
                f"{dict(coeffs.coefficients)!r}."
            )
        stamped = coeffs.effective_binding_energy_I_ion_eV
        if cfg.binding_energy_I_ion_eV != stamped:
            raise AssertionError(
                f"[{cell.label}] E_bind {cfg.binding_energy_I_ion_eV!r} != "
                f"bundle stamp {stamped!r} on the bundle-paired partner."
            )


def expected_diff_keys(cell: LinRingCell) -> frozenset[str]:
    """The pre-registered cfg-diff key set for this cell vs the reference."""
    keys = set(_BASE_DIFF_KEYS)
    if cell.a is not None:
        keys |= _LIN_DIFF_KEYS
        if EBIND_EV[cell.eb_tag] != EBIND_EV["eb1168"]:
            keys.add("binding_energy_I_ion_eV")
    return frozenset(keys)


def verify_against_reference(cfg: SimConfig, cell: LinRingCell) -> list[str]:
    """Diff vs the standing battery reference; only the pre-registered keys
    may differ (LR-P1 second half — an unexpected key = a drifted pin).

    ``droplet_size_sampler_mode`` handling as in gen_tier2atlas_g3ring.py:
    the reference cfg.json predates the field, so the diff runs on a copy
    pinned to the default; the real ``"raw"`` value is asserted by
    :func:`verify_corrected_geometry`.
    """
    ref_path = PROJECT_ROOT / "data" / "runs" / REFERENCE_RUN / "cfg.json"
    cfg_for_diff = dataclasses.replace(
        cfg, droplet_size_sampler_mode="post_pickup"
    )
    diff = cfg_diff_vs_reference(cfg_for_diff, ref_path, context=cell.label)
    expected = expected_diff_keys(cell)
    if set(diff) != expected:
        raise AssertionError(
            f"[{cell.label}] cfg diff vs {REFERENCE_RUN} is {sorted(diff)}; "
            f"pre-registered exactly {sorted(expected)}."
        )
    return diff


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _run_one(cell: LinRingCell) -> str:
    cfg = build_cell(cell)
    verify_corrected_geometry(cfg, cell)
    verify_against_reference(cfg, cell)
    run_dir = PROJECT_ROOT / "data" / "runs" / ring_run_dir_name(cell.label)
    if run_dir.exists():
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{cell.label}] skip (complete) -> {run_dir}", flush=True)
            return cell.label
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(f"{run_dir} exists and is not complete.")
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    chord = (f"a={cell.a}" if cell.a is not None
             else f"v_c={cell.v_c} p_tail={P_TAIL}")
    print(f"[{cell.label}] {cfg.drag_form} {chord} {cell.eb_tag} "
          f"tau={cell.tau_ps} E0={cell.e0_eV} ({cell.role}) seed={cfg.seed}",
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


def _select(labels: list[str]) -> list[LinRingCell]:
    if not labels:
        return list(RING_MATRIX)
    unknown = sorted(set(labels) - set(_SPEC_BY_LABEL))
    if unknown:
        raise SystemExit(f"unknown label(s) {unknown}; expected among "
                         f"{[c.label for c in RING_MATRIX]}")
    return [_SPEC_BY_LABEL[label] for label in labels]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cells = _select(args.labels)
    print(f"Atlas free-form linear MD ring (plan §6.1): {len(cells)} cell(s), "
          f"N={N}, seed={SEED} (fresh, shared CRN incl. the h405 partner), "
          f"corrected geometry ({SIZE_PRIOR}+{SIZE_SAMPLER_MODE}, {BIRTH_LAW} "
          f"{PARENT_WELL_K} K), concurrency={args.concurrency}"
          f"{' [DRY RUN]' if args.dry_run else ''} "
          "(instrument runs -- standing point does not move).", flush=True)

    verify_twin_preregistration()

    if args.dry_run:
        for cell in cells:
            cfg = build_cell(cell)
            verify_corrected_geometry(cfg, cell)
            diff = verify_against_reference(cfg, cell)
            chord = (f"a={cell.a}" if cell.a is not None
                     else f"v_c={cell.v_c}")
            print(f"[{cell.label}] cfg OK {cfg.drag_form} {chord} "
                  f"{cell.eb_tag} tau={cell.tau_ps} E0={cell.e0_eV} -> "
                  f"{ring_run_dir_name(cell.label)}", flush=True)
            print(f"[{cell.label}]   diff vs standing: {sorted(diff)}",
                  flush=True)
        print("dry run complete: twin oracle + all cfg guards pass.",
              flush=True)
        return 0

    concurrency = max(1, args.concurrency)
    if concurrency == 1 or len(cells) == 1:
        for cell in cells:
            _run_one(cell)
        return 0

    from concurrent.futures import ProcessPoolExecutor, as_completed
    workers = min(concurrency, len(cells))
    print(f"launching {len(cells)} cells through a {workers}-slot pool ...",
          flush=True)
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(_run_one_by_label, c.label): c.label
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
