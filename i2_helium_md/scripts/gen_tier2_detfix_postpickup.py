"""Droplet-size sampler A/B at h405 — the production battery re-run with ``post_pickup``.

The production pipeline (``gen_tier2_detfix_battery.py``) samples droplet sizes
with ``droplet_size_sampler_mode="raw"``: plain log-normal draws from the source
distribution (40 bar / 14 K / 5 µm, ⟨N⟩ = 12794, δ = 0.625). The G0-1 decision
(D0 §15.5–15.6) chose ``raw`` from a quantile match to the parent document's
quoted droplet range; the ``post_pickup`` arm (single-I₂ pickup conditioning +
pickup-induced evaporation — the legacy MATLAB production sampler) was never
measured at the standing point. This battery measures it.

Members: the six production members (``s1``–``s5``, ``r6``) with the same seed
and pipeline, **one field flipped**: ``droplet_size_sampler_mode="post_pickup"``.

**Not CRN.** ``post_pickup`` oversamples the log-normal (≥ 10 000 draws) and
runs the pickup Monte Carlo on the same ``default_rng(cfg.seed)`` stream before
the birth positions and velocities are drawn, so every downstream draw differs
from the production twin. The comparison is battery-vs-battery (pooled
N = 6000 each), not paired per ion.

Oracles before any MD:

1. ``verify_member_vs_production`` — the cfg diffs against the committed
   production member's ``cfg.json`` in exactly :data:`AB_DIFF_KEYS`.
2. ``verify_handover_alignment`` (re-used) — the last stored ion column is the
   final state at t_h = 500 ps.
3. ``verify_postpickup_sizes`` — the sampler realizes the pickup-weighted
   ensemble (sample mean within :data:`POSTPICKUP_MEAN_BAND` He of the
   expected ≈ 16.4k), checked on the actual cfg and seed before launch.

Usage::

    python scripts/gen_tier2_detfix_postpickup.py --dry-run   # oracles only
    python scripts/gen_tier2_detfix_postpickup.py             # all 6, 2 at a time
    python scripts/gen_tier2_detfix_postpickup.py s1 r6       # selected

Pool afterwards (figures containers, N = 5000 each, s1–s5; see
:data:`POOL_MEMBERS` for why r6 is out)::

    python scripts/build_pooled_detection_container.py --target detfix_h405_postpickup
    python scripts/build_pooled_detection_container.py --target detfix_h405_s1s5

**lq A5 cell (THESIS_FINAL_GOALS G5, 2026-10-08).** ``--cell lq`` runs the
same flip on the thesis production law — the committed G4.14 A5 member
(``gen_tier2_g414_tauline.confirm_cells``: lq c 2.0, a 16.717, E_bind 0.1175,
3.5 Å gate, τ 6.4, E₀ 0.3155) with ``droplet_size_sampler_mode="post_pickup"``
and, for B1/B1'/B2, E₀ or τ re-placed. Oracle 1 then diffs against the committed A5
member's ``cfg.json`` (exactly the sampler, plus E₀ / τ when re-placed); oracles 2
and 3 are unchanged. ion.npz is dropped once detection.npz reloads with scored
ions (the scorer reads cfg.json + neutral.npz + detection.npz)::

    python scripts/gen_tier2_detfix_postpickup.py --cell lq --dry-run
    python scripts/gen_tier2_detfix_postpickup.py --cell lq          # A: s1-s5
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b1     # B1: E0 column, s1
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b1tau  # B1': tau column, s1
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b1c    # B1c: ridge c column, s1
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b1w    # B1w: 0.0482 well, a column, s1
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b3     # B3: s1-s5 + r6 at a = B3_A
    python scripts/gen_tier2_detfix_postpickup.py --cell lq --b2     # B2: s1-s5 + r6 at the B2 cell

Pools: ``--target g5_lq_s1s5`` (raw A5 arm, N = 5000), ``g5_lqpp`` (A) and,
once the B2 cell is set, ``g5_lqppb2`` (N = 6000). Scorer:
``scripts/post_processing/tier2_g5_postpickup_table.py``.

**Thesis production (THESIS_FINAL_GOALS G13, 2026-10-26).** ``--b2`` / ``--b3``
now run under ``self_bound_gate_release="decoupled"`` (:data:`THESIS_GATE_RELEASE`)
into ``data/runs/B2_N1000_s1`` ... ``_s6`` (``B3_…``; r6 -> s6), ion.npz kept;
the A5 oracle is unchanged and the release is the one field added on top.
The committed ``never`` B2 / B3 batteries keep their long names (history).
Pools: ``gen_tier2_g13_decoupled.py --pool`` -> ``B2_N6000_pooled`` / ``B3_N6000_pooled``.

Max two concurrent MD runs (``--concurrency`` is capped at 2).
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
from pathlib import Path
import sys
from typing import Optional

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MAX_BYTES_ION, MEMBER_SEEDS, T_HANDOVER_PS, build_member,
    member_run_dir_name as production_run_dir_name, verify_handover_alignment,
)
from scripts.gen_tier2_g413_ridge_battery import DRAG_FIELDS, LQ_BUNDLE, _E0_KEY, _TAU_KEY  # noqa: E402
from scripts.gen_tier2_g414_tauline import (  # noqa: E402
    ION_DROPPED_MARKER, LINE_CSV, Cell, build as a5_build,
    confirm_cells as a5_confirm_cells, drop_ion_npz, run_dir_for as a5_run_dir_for,
)
from scripts.gen_tier2atlas_g4finals import BUDGET_EV, CASE, N, VARIANT  # noqa: E402
from scripts.tier0_common import build_drag_cfg, run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes  # noqa: E402
from i2_helium_md.simulation.detection_stage import run_detection_stage  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

# =============================================================================
# USER SETTINGS
# =============================================================================

SAMPLER_MODE = "post_pickup"
#: The only field allowed to differ from the committed production member.
AB_DIFF_KEYS = frozenset({"droplet_size_sampler_mode"})
#: Expected post-pickup mean ⟨N⟩ at 40 bar / 14 K (D0 §15.2: 16.3k–16.7k
#: realized; raw is 12.6k–12.9k). The band separates the two arms by > 3k He.
POSTPICKUP_MEAN_BAND = (15_000.0, 18_000.0)

#: The A/B pool (both arms, N = 5000). r6 is excluded: its post_pickup run
#: tripped the detection-stage marginal-partner safeguard (D4, first firing:
#: one slow escaper, ion 610, 56 Å outside the surface at 0.12 Å/ps, still
#: helium-coupled at t_h = 500 ps while its partner had escaped). Re-running
#: r6 at a longer t_h would make the arms differ in two fields; user decision
#: 2026-10-02: compare s1–s5 against production's s1–s5 (`pooled_new`).
POOL_MEMBERS = ("s1", "s2", "s3", "s4", "s5")

# lq A5 cell (G5). B1: E0 column at s1 around the projected re-placement
# (h405 A/B: post_pickup raised nbar by +0.72; dnbar/dE0 ~ -25 /eV, so
# +0.03 eV ~ -0.75); the A run at s1 is the column's E0 = 0.3155 point.
B1_E0_EV = (0.3355, 0.3455, 0.3555)
B1_MEMBER = "s1"
# B1' (user, 2026-10-08): E0 re-placement trades KE (E0 is taken from the
# Coulomb budget), so nbar is restored with the cooling time instead: tau
# column at s1, A5 E0. Estimate from the committed g414_tradeoff line:
# dnbar/dtau ~ -1.6 /ps at fixed E0 -> tau ~ 6.85 ps for -0.7 (beyond the
# swept 4.0-6.6 ps). The A run at s1 is the tau = 6.4 point.
B1_TAU_PS = (6.8, 7.2, 7.6)
# B1c (user, 2026-10-08): tau and E0 both trade KE (D0 §15.9: they change the
# chord -> (n, KE) map); the drag does too, in the other direction. Lower c
# on the TDDFT ridge raises KE at matched nbar (atlas_g414_filter.csv), so
# post_pickup + E0 re-placement + lower c may hold both. Per c: the ridge
# point (committed g412_ridge JSON), tau 6.4, E0 = the committed twin line's
# nbar 3.89 value (interpolated in nbar from 3.85 / 4.05 at tau 6.4) plus the
# post_pickup E0 shift measured in B1 (s1: 0.3155 -> 0.3368 eV for nbar 3.87).
B1C_C = (1.6, 1.8)
B1C_TAU_PS = 6.4
B1C_NBAR = 3.89
PP_E0_SHIFT_EV = 0.0213
#: B2 E0 [eV] / tau [ps]: set (user) from B1 / B1' before ``--b2``; None =
#: the A5 cell's own value.
# B2 (user, 2026-10-08): the B1c c = 1.6 cell (s1: nbar 3.94, midHot 0.994,
# W1 0.93, S 1.67) at all six production members, pooled N = 6000.
B2_C: Optional[float] = 1.6
B2_E0_EV: Optional[float] = 0.3297     # = b1c_e0(1.6); checked in main
B2_TAU_PS: Optional[float] = None
# B1w (user, 2026-10-09): the G10 high-tail recipe on the B2 cell -- shallow
# escape well E_bind 0.0482 eV with a raised to restore midHot (G10 on the
# raw A5 cell: a 16.72 -> 19.5 undid the well's midHot +0.135). Effective
# model, not Tier-0-consistent (G10.5). E0 per a holds nbar at the B2 value:
# B2 E0 + B1W_WELL_E0_EV (the well's nbar rise, raw pools 3.886 -> 4.05, at
# dnbar/dE0 ~ -31 /eV) + B1W_DE0_DA * (a - ridge a) (the committed twin's
# nbar-matched dE0/da at 0.0482, atlas_g415_summary.csv a 18.5 -> 20.5).
B1W_C = 1.6
B1W_E_BIND_EV = 0.0482
B1W_A = (19.8, 20.8, 21.8)
B1W_WELL_E0_EV = 0.0053
B1W_DE0_DA = 0.002185
#: B3 (user, 2026-10-09): the B1w a 19.8 cell (s1: midHot 1.004, n >= 15
#: 0.0085, W1 0.767) at all six production members, pooled N = 6000.
B3_A = 19.8
#: Thesis production (THESIS_FINAL_GOALS G13, user 2026-10-26): B2 / B3
#: (``--b2`` / ``--b3``) and the re-landed h405 contrast law
#: (``gen_tier2_h405_pp35.py --scan-pool``) run with the self-bound gate lifted
#: outside the droplet, into plain run dirs ``<law>_N1000_s1`` ... ``_s6``
#: (member r6 -> s6; every other detail is in the run's cfg.json), ion.npz
#: kept. The exploratory columns (A5, B1, B1', B1c, B1w) keep the committed
#: ``never`` release and their long run-dir names.
THESIS_GATE_RELEASE = "decoupled"
LQ_DROP_REASON = ("G5 lq post_pickup: disk (T: full); scorer reads cfg.json + "
                  "neutral.npz + detection.npz")

MAX_CONCURRENCY = 2
SKIP_COMPLETED_RUNS = True
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", "detection.npz")


def _select(labels: list[str], default: tuple[str, ...] = tuple(MEMBER_SEEDS)) -> list[str]:
    if not labels:
        return list(default)
    unknown = sorted(set(labels) - set(MEMBER_SEEDS))
    if unknown:
        raise SystemExit(f"unknown member(s) {unknown}; expected among "
                         f"{list(MEMBER_SEEDS)}")
    return labels


def member_run_dir_name(member: str) -> str:
    """Run-dir basename; ``pp`` marks the post_pickup arm of the A/B."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}")
    return run_dir_name(CASE, VARIANT, N,
                        run_tag=f"detfix_conf270_h405pp{member}_th500")


def production_run_dir(member: str) -> Path:
    """The committed production member this A/B member mirrors."""
    return PROJECT_ROOT / "data" / "runs" / production_run_dir_name(member)


def build_postpickup_member(member: str) -> SimConfig:
    """The production member's cfg with only the sampler mode flipped."""
    cfg = dataclasses.replace(build_member(member),
                              droplet_size_sampler_mode=SAMPLER_MODE)
    cfg.validate()
    return cfg


def verify_member_vs_production(member: str) -> list[str]:
    """Oracle 1: the cfg differs from the committed production cfg.json in
    exactly :data:`AB_DIFF_KEYS`. Returns the diff keys."""
    cfg = build_postpickup_member(member)
    diff = cfg_diff_vs_reference(cfg, production_run_dir(member) / "cfg.json",
                                 context=f"post_pickup {member}")
    if set(diff) != AB_DIFF_KEYS:
        raise AssertionError(
            f"[{member}] cfg vs production differs in {sorted(diff)}; "
            f"pre-registered exactly {sorted(AB_DIFF_KEYS)}.")
    return diff


def verify_postpickup_sizes(cfg: SimConfig) -> float:
    """Oracle 3: the sampler realizes the pickup-weighted ensemble on this
    cfg and seed (same call the neutral stage makes). Returns the sample
    mean ⟨N⟩ [He atoms]."""
    sizes = sample_droplet_sizes(cfg, mode=cfg.droplet_size_sampler_mode,
                                 rng=np.random.default_rng(cfg.seed))
    mean_N = float(np.mean(sizes))
    lo, hi = POSTPICKUP_MEAN_BAND
    if sizes.size != cfg.num_molecules or not lo <= mean_N <= hi:
        raise AssertionError(
            f"post_pickup sizes: {sizes.size} draws, mean {mean_N:.0f} He "
            f"(expected {cfg.num_molecules} draws, mean in [{lo:.0f}, {hi:.0f}]).")
    return mean_N


def a5_member_cell(member: str):
    """The committed G4.14 A5 confirmation cell of a production member."""
    return next(c for c in a5_confirm_cells() if c.label.endswith(f"conf{member}"))


#: Drag-law keys a ridge-point change moves (c, its ridge a and E_bind).
LAW_DIFF_KEYS = frozenset({"drag_coefficients", "binding_energy_I_ion_eV"})


def b1c_e0(c: float) -> float:
    """E0 [eV] for the B1c cell at ridge c: the committed twin line
    (atlas_g414_twin_line.csv, tau B1C_TAU_PS) linearly interpolated in nbar
    to B1C_NBAR, plus PP_E0_SHIFT_EV; rounded to 4 dp. Raises if the line
    rows are missing."""
    with open(LINE_CSV, newline="", encoding="utf-8") as fh:
        pts = sorted((float(r["nbar_target_md"]), float(r["E0_eV"]))
                     for r in csv.DictReader(fh)
                     if float(r["c"]) == c and float(r["tau_ps"]) == B1C_TAU_PS)
    if len(pts) != 2:
        raise AssertionError(f"twin line: {len(pts)} rows at c {c}, tau {B1C_TAU_PS}")
    (n0, e0), (n1, e1) = pts
    return round(e0 + (B1C_NBAR - n0) / (n1 - n0) * (e1 - e0) + PP_E0_SHIFT_EV, 4)


def b1w_e0(a: float) -> float:
    """E0 [eV] for the B1w cell at drag a [amu/ps] (see B1W_* constants)."""
    from scripts.gen_tier2_g414_tauline import ridge_point
    return round(b1c_e0(B1W_C) + B1W_WELL_E0_EV
                 + B1W_DE0_DA * (a - ridge_point(B1W_C)["a"]), 4)


def lq_label(member: str, e0_eV: Optional[float] = None,
             tau_ps: Optional[float] = None, c: Optional[float] = None,
             a: Optional[float] = None, e_bind_eV: Optional[float] = None) -> str:
    """e.g. ``g5lqpps1`` (A5 cell), ``g5lqppe3455s1`` (E0 0.3455 eV),
    ``g5lqppt68s1`` (tau 6.8 ps), ``g5lqppc160e3297s1`` (ridge c 1.6),
    ``g5lqppc160a208w0482e3411s1`` (c 1.6, a 20.8, E_bind 0.0482)."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}")
    tag = ("" if c is None else f"c{round(c * 100):03d}") + (
        "" if a is None else f"a{round(a * 10):03d}") + (
        "" if e_bind_eV is None else f"w{round(e_bind_eV * 1e4):04d}") + (
        "" if e0_eV is None else f"e{round(e0_eV * 1e4):04d}") + (
        "" if tau_ps is None else f"t{round(tau_ps * 10):d}")
    return f"g5lqpp{tag}{member}"


def lq_run_dir(member: str, e0_eV: Optional[float] = None,
               tau_ps: Optional[float] = None, c: Optional[float] = None,
               a: Optional[float] = None, e_bind_eV: Optional[float] = None) -> Path:
    label = lq_label(member, e0_eV, tau_ps, c, a, e_bind_eV)
    name = run_dir_name(CASE, LQ_BUNDLE, N, run_tag=f"detfix_conf270_{label}_th500")
    return PROJECT_ROOT / "data" / "runs" / name


def build_lq_member(member: str, e0_eV: Optional[float] = None,
                    tau_ps: Optional[float] = None, c: Optional[float] = None,
                    a: Optional[float] = None, e_bind_eV: Optional[float] = None,
                    ) -> SimConfig:
    """The A5 member cfg with the sampler flipped and, optionally, E0 [eV],
    the cooling time tau [ps], the ridge point c [amu/Å] (with its committed
    a and E_bind) and / or an off-ridge a [amu/ps] / E_bind [eV] re-placed
    (the B1w effective law; not Tier-0-consistent)."""
    cell = a5_member_cell(member)
    if c is not None:
        cell = Cell(cell.label, c, cell.tau_ps, cell.nbar_target, cell.e0_eV, cell.seed)
    base = a5_build(cell)
    repl: dict = {"droplet_size_sampler_mode": SAMPLER_MODE}
    if a is not None or e_bind_eV is not None:
        coeffs = base.drag_coefficients.coefficients
        drag = build_drag_cfg(
            CASE, LQ_BUNDLE, num_molecules=base.num_molecules,
            ion_time_ps=base.ion_simulation_time, dt_ion_ps=base.dt_ion, seed=base.seed,
            coeff_overrides={"a": coeffs["a"] if a is None else float(a), "c": coeffs["c"]},
            e_bind_override=(base.binding_energy_I_ion_eV if e_bind_eV is None
                             else float(e_bind_eV)))
        repl.update({k: getattr(drag, k) for k in DRAG_FIELDS})
        if repl["drag_gate_steepness"] != base.drag_gate_steepness:
            raise AssertionError(f"bundle gate {repl['drag_gate_steepness']} != "
                                 f"{base.drag_gate_steepness}")
    if e0_eV is not None:
        repl[_E0_KEY] = float(e0_eV) / BUDGET_EV
    if tau_ps is not None:
        repl[_TAU_KEY] = float(tau_ps)
    cfg = dataclasses.replace(base, **repl)
    cfg.validate()
    return cfg


def verify_lq_member(member: str, e0_eV: Optional[float] = None,
                     tau_ps: Optional[float] = None, c: Optional[float] = None,
                     a: Optional[float] = None, e_bind_eV: Optional[float] = None,
                     ) -> list[str]:
    """Oracle 1 (lq): the cfg differs from the committed A5 member's cfg.json
    in exactly the sampler (+ E0 / tau / drag law when re-placed). Returns
    the diff keys."""
    job = (member, e0_eV, tau_ps, c, a, e_bind_eV)
    cfg = build_lq_member(*job)
    ref = a5_run_dir_for(a5_member_cell(member)) / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref, context=f"lq post_pickup {member}")
    law = c is not None or a is not None or e_bind_eV is not None
    want = (set(AB_DIFF_KEYS) | ({_E0_KEY} if e0_eV is not None else set())
            | ({_TAU_KEY} if tau_ps is not None else set())
            | (set(LAW_DIFF_KEYS) if law else set()))
    if set(diff) != want:
        raise AssertionError(f"[{lq_label(*job)}] cfg vs A5 differs in "
                             f"{sorted(diff)}; pre-registered exactly {sorted(want)}.")
    return diff


def run_lq(job: tuple) -> str:
    """lq A5 post_pickup member: guard, neutral -> ion -> detection, drop ion.npz.
    ``job`` = (member, e0_eV, tau_ps[, c[, a[, e_bind_eV]]])."""
    job = (*job, None, None, None)[:6]
    label = lq_label(*job)
    verify_lq_member(*job)
    cfg = build_lq_member(*job)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = lq_run_dir(*job)
    done = (all((run_dir / n).exists() for n in ("cfg.json", "neutral.npz", "detection.npz"))
            and ((run_dir / "ion.npz").exists() or (run_dir / ION_DROPPED_MARKER).exists()))
    if SKIP_COMPLETED_RUNS and done:
        drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
        print(f"[{label}] skip (complete) -> {run_dir.name}", flush=True)
        return label
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{label}] seed={cfg.seed} E0 {getattr(cfg, _E0_KEY) * BUDGET_EV:.4f} eV "
          f"tau {getattr(cfg, _TAU_KEY)} ps sampler={SAMPLER_MODE} ion {steps} steps "
          f"(stride {stride}) -> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(f"[{label}] ion.npz last stored column at {t_last} ps, "
                             f"not t_h {T_HANDOVER_PS} ps.")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    del ion, neutral
    drop_ion_npz(run_dir, reason=LQ_DROP_REASON)
    dm = det.detected_mask
    print(f"[{label}] done: scored {int(dm.sum())}/{dm.size}, retained "
          f"{int(np.count_nonzero(~dm))}", flush=True)
    return label


#: Thesis lq laws: name -> job tail after the member (e0, tau, c[, a, e_bind]).
THESIS_LQ_JOBS: dict[str, tuple] = {
    "B2": (B2_E0_EV, B2_TAU_PS, B2_C),
    "B3": (b1w_e0(B3_A), None, B1W_C, B3_A, B1W_E_BIND_EV),
}


def thesis_member_label(member: str) -> str:
    """Thesis run-dir member label: ``s1``-``s5`` unchanged, ``r6`` -> ``s6``."""
    if member not in MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}")
    return "s6" if member == "r6" else member


def thesis_run_dir(law: str, member: str) -> Path:
    """Thesis production member dir, e.g. ``data/runs/B2_N1000_s1``."""
    return PROJECT_ROOT / "data" / "runs" / f"{law}_N{N}_{thesis_member_label(member)}"


def thesis_pool_dir(law: str) -> Path:
    """Thesis pooled figures container, e.g. ``data/runs/B2_N6000_pooled``."""
    return PROJECT_ROOT / "data" / "runs" / f"{law}_N{N * len(MEMBER_SEEDS)}_pooled"


def build_thesis_lq(law: str, member: str) -> SimConfig:
    """Thesis B2 / B3 member: the lq post_pickup cell under
    :data:`THESIS_GATE_RELEASE`."""
    cfg = dataclasses.replace(build_lq_member(member, *THESIS_LQ_JOBS[law]),
                              self_bound_gate_release=THESIS_GATE_RELEASE)
    cfg.validate()
    return cfg


def run_thesis_member(cfg: SimConfig, run_dir: Path, label: str) -> str:
    """neutral -> ion (to t_h) -> detection into ``run_dir``, ion.npz kept,
    plus the RNG-free Stage-I event log ``stage1_events.npz`` (§16).

    A complete dir is skipped if its cfg.json equals ``cfg`` and refused
    otherwise (never silently accepts a different model)."""
    steps, stride = verify_handover_alignment(cfg)
    if all((run_dir / f).exists() for f in _REQUIRED_ARTIFACTS):
        if RunDirectory(run_dir).load_cfg() != cfg:
            raise RuntimeError(f"[{label}] existing {run_dir.name} has a different cfg")
        print(f"[{label}] skip (complete) -> {run_dir.name}", flush=True)
        return label
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{label}] seed {cfg.seed}, {steps} steps (stride {stride}) -> {run_dir.name}",
          flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION, record_events=True)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(f"[{label}] ion.npz last column at {t_last} ps, not t_h")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    print(f"[{label}] done: scored {int(det.detected_mask.sum())}/{det.detected_mask.size}",
          flush=True)
    return label


def run_thesis_lq(job: tuple[str, str]) -> str:
    """Thesis B2 / B3 member: A5 oracle, then :func:`run_thesis_member`.
    ``job`` = (law, member)."""
    law, member = job
    verify_lq_member(member, *THESIS_LQ_JOBS[law])
    return run_thesis_member(build_thesis_lq(law, member), thesis_run_dir(law, member),
                             f"{law} {member}")


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def run_one(member: str) -> str:
    """Guard, then neutral -> ion (to t_h) -> detection (skip path + closure)."""
    verify_member_vs_production(member)
    cfg = build_postpickup_member(member)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = PROJECT_ROOT / "data" / "runs" / member_run_dir_name(member)
    if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
        print(f"[{member}] skip (complete) -> {run_dir}", flush=True)
        return member
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{member}] seed={cfg.seed} sampler={SAMPLER_MODE} ion {steps} steps "
          f"(stride {stride}) to t_h={T_HANDOVER_PS} ps -> {run_dir}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{member}] ion ...", flush=True)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise AssertionError(
            f"[{member}] ion.npz last stored column at {t_last} ps, not t_h "
            f"{T_HANDOVER_PS} ps.")
    print(f"[{member}] detection (seeded from ion.npz at {t_last} ps) ...",
          flush=True)
    detect = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    dm = detect.detected_mask
    print(f"[{member}] done: scored {int(dm.sum())}/{dm.size}, retained "
          f"{int(np.count_nonzero(~dm))}", flush=True)
    return member


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", nargs="*")
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cell", choices=("h405", "lq"), default="h405",
                        help="h405 (the 2026-10-02 A/B) or the lq A5 thesis cell (G5)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--b1", action="store_true", help="lq: the B1 E0 column at s1")
    group.add_argument("--b1tau", action="store_true", help="lq: the B1' tau column at s1")
    group.add_argument("--b1c", action="store_true",
                       help="lq: the B1c ridge-c column at s1 (E0 re-placed per c)")
    group.add_argument("--b1w", action="store_true",
                       help="lq: the B1w shallow-well a column at s1 (c 1.6, E_bind 0.0482)")
    group.add_argument("--b3", action="store_true",
                       help="lq: s1-s5 + r6 at the B3 cell (B1w, a = B3_A)")
    group.add_argument("--b2", action="store_true",
                       help="lq: s1-s5 + r6 at the B2 cell (B2_C / B2_E0_EV / B2_TAU_PS)")
    args = parser.parse_args(argv)
    if args.cell == "h405" and (args.b1 or args.b1tau or args.b1c or args.b1w or args.b2 or args.b3):
        raise SystemExit("--b1 / --b1tau / --b1c / --b1w / --b2 / --b3 apply to --cell lq only")

    if args.cell == "h405":
        members = _select(args.labels)
        for member in members:
            diff = verify_member_vs_production(member)
            cfg = build_postpickup_member(member)
            steps, stride = verify_handover_alignment(cfg)
            mean_N = verify_postpickup_sizes(cfg)
            print(f"[{member}] guards PASSED: seed {MEMBER_SEEDS[member]}, diff vs "
                  f"production {diff}, {steps} steps, stride {stride}, "
                  f"post_pickup <N> {mean_N:.0f} -> {member_run_dir_name(member)}",
                  flush=True)
        fn, jobs = run_one, members
    elif args.b2 or args.b3:
        # Thesis production: decoupled release, plain run dirs (THESIS_GATE_RELEASE).
        if args.b2 and B2_E0_EV is None and B2_TAU_PS is None and B2_C is None:
            raise SystemExit("B2 cell not set: choose from B1 / B1' / B1c first")
        if args.b2 and B2_C is not None and B2_E0_EV != b1c_e0(B2_C):
            raise SystemExit(f"B2_E0_EV {B2_E0_EV} != b1c_e0({B2_C}) = {b1c_e0(B2_C)}")
        law = "B2" if args.b2 else "B3"
        jobs = [(law, m) for m in _select(args.labels, tuple(MEMBER_SEEDS))]
        for _, m in jobs:
            diff = verify_lq_member(m, *THESIS_LQ_JOBS[law])
            cfg = build_thesis_lq(law, m)
            steps, stride = verify_handover_alignment(cfg)
            mean_N = verify_postpickup_sizes(cfg)
            print(f"[{law} {m}] guards PASSED: seed {cfg.seed}, diff vs A5 {diff} + "
                  f"self_bound_gate_release={THESIS_GATE_RELEASE}, {steps} steps, "
                  f"stride {stride}, post_pickup <N> {mean_N:.0f} -> "
                  f"{thesis_run_dir(law, m).name}", flush=True)
        fn = run_thesis_lq
    else:
        if args.b1:
            jobs = [(m, e0, None) for m in _select(args.labels, (B1_MEMBER,))
                    for e0 in B1_E0_EV]
        elif args.b1w:
            jobs = [(m, b1w_e0(a), None, B1W_C, a, B1W_E_BIND_EV)
                    for m in _select(args.labels, (B1_MEMBER,)) for a in B1W_A]
        elif args.b1c:
            jobs = [(m, b1c_e0(c), None, c) for m in _select(args.labels, (B1_MEMBER,))
                    for c in B1C_C]
        elif args.b1tau:
            jobs = [(m, None, t) for m in _select(args.labels, (B1_MEMBER,))
                    for t in B1_TAU_PS]
        else:
            jobs = [(m, None, None) for m in _select(args.labels, POOL_MEMBERS)]
        for job in jobs:
            diff = verify_lq_member(*job)
            cfg = build_lq_member(*job)
            steps, stride = verify_handover_alignment(cfg)
            mean_N = verify_postpickup_sizes(cfg)
            print(f"[{lq_label(*job)}] guards PASSED: seed {cfg.seed}, E0 "
                  f"{getattr(cfg, _E0_KEY) * BUDGET_EV:.4f} eV, tau {getattr(cfg, _TAU_KEY)} ps, "
                  f"diff vs A5 {diff}, {steps} steps, stride {stride}, post_pickup <N> "
                  f"{mean_N:.0f} -> {lq_run_dir(*job).name}", flush=True)
        fn = run_lq
    if args.dry_run:
        return 0

    workers = max(1, min(args.concurrency, MAX_CONCURRENCY, len(jobs)))
    if workers == 1:
        for job in jobs:
            fn(job)
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(fn, j): j for j in jobs}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    if failures:
        print(f"{failures}/{len(jobs)} run(s) failed.", flush=True)
        return 1
    print(f"all {len(jobs)} runs complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
