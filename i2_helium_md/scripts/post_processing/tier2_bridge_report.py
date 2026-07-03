"""Phase-D bridge report (Slice Z): generative vs anchored/TDDFT at 9 Å / 0.80 eV.

Reads the finished ``tier2_bridge_biphasic`` run directory and emits the
bridge diagnostic (plan §3 comparison + report step):

* emergent mean ``n(t)`` overlaid on the anchored staircase *family*
  ``build_shell_schedule(t★)`` for t★ ∈ {0.5, 5.0, 9.0} ps (t★ = 5.0 primary);
* ``R(t)`` / ``|v(t)|`` vs the frozen 9 Å TDDFT reference
  (``compare_distance`` / ``compare_velocity_magnitude``);
* reconstructed per-ion ``t×`` vs the §2.2 closed form
  ``τ·ln(f_int·E_avail/Σ(21))`` (sharp), the GAH25 5–6.5 ps prior (±factor 2)
  and the [1, 15] ps sanity band;
* the ``Π(t)`` regime trace (mean ± envelope; corrected freeze-side
  expectation — Π = 0 at gate-open, Π < 1 throughout, Π → 0 at exit);
* the 5-term ledger residual (``ion_ledger_closure``).

Overlay figures go into ``<run_dir>/figures/``; the numeric summary is printed
and written to ``<run_dir>/bridge_summary.txt``. The verdict is **editorial**
(reported, not auto-adjudicated) and lives in
``docs/drag_port/Tier2/TIER2_PHASE_D_BRIDGE_FINDINGS.md``.

Usage::

    python scripts/post_processing/tier2_bridge_report.py
"""

from __future__ import annotations

from pathlib import Path
import sys

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

RUN_DIR = None  # None -> data/runs/<tier2_bridge_run_dir_name(CASE, VARIANT, N)>
HEDFT_CSV_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"

T_STAR_FAMILY_PS = (0.5, 5.0, 9.0)
PRIMARY_T_STAR_PS = 5.0

# GAH25 shell-1 t0 prior (Na+ number; I+ is Rb+-like -> +/- factor 2 expected,
# factor 10 is the genuine flag) and the MASS section-6.11 sanity band.
GAH25_PRIOR_PS = (5.0, 6.5)
T_X_SANITY_BAND_PS = (1.0, 15.0)

SHOW_FIGURES = False
SAVE_FIGURES = True


import matplotlib  # noqa: E402

if not SHOW_FIGURES:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scripts.tier2_common import tier2_bridge_run_dir_name  # noqa: E402
from i2_helium_md.physics.dissociation_ladder import ladder_cumsum  # noqa: E402
from i2_helium_md.physics.shell_schedule import build_shell_schedule  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    compare_distance,
    compare_velocity_magnitude,
    ion_ledger_closure,
    load_hedft_trajectory,
)
from i2_helium_md.postprocess.bridge_diagnostics import (  # noqa: E402
    crossing_time_ps,
    mean_shell_count,
    regime_parameter,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _fmt_band(band: tuple[float, float]) -> str:
    return f"[{band[0]:g}, {band[1]:g}] ps"


def build_summary_lines(cfg, ion, hedft) -> tuple[list[str], dict]:
    """Compute every bridge number and return (report lines, plot payload)."""
    time_ps = ion.time_ps
    lines: list[str] = []
    payload: dict = {}

    lines.append("Tier-2 Phase-D bridge report (Slice Z)")
    lines.append(
        f"run: {CASE} {VARIANT} N={N} scenario={cfg.mass_scenario} "
        f"E_avail={cfg.coulomb_available_eV:.2f} eV"
    )
    lines.append(
        f"pinned point: lambda0={cfg.pickup_rate_coefficient:g}/ps "
        f"f_int={cfg.internal_energy_partition_fraction:g} "
        f"f_ret={cfg.internal_energy_retained_fraction:g} "
        f"tau={cfg.internal_energy_cooling_tau_ps:g} ps "
        f"kappa={cfg.ladder_steepness:g} picture={cfg.ladder_electronic_picture}"
    )
    lines.append("")

    # --- 1. emergent mean n(t) vs the anchored staircase family -------------
    mean_n = mean_shell_count(ion.n_shell)
    payload["mean_n"] = mean_n
    payload["schedules"] = {
        t_star: build_shell_schedule(t_star).n_of_t(time_ps)
        for t_star in T_STAR_FAMILY_PS
    }
    lines.append("[mean n(t) vs anchored 21->19->14 staircase family]")
    probe_times = [0.0, 5.0, 10.0, 20.0, float(time_ps[-1])]
    for t_probe in probe_times:
        idx = int(np.argmin(np.abs(time_ps - t_probe)))
        stair = ", ".join(
            f"t*={t_star:.1f}: {payload['schedules'][t_star][idx]:d}"
            for t_star in T_STAR_FAMILY_PS
        )
        lines.append(
            f"  t={time_ps[idx]:5.1f} ps  mean n = {mean_n[idx]:6.2f}   "
            f"anchored ({stair})"
        )
    lines.append(
        f"  terminal: mean n = {mean_n[-1]:.2f}, per-ion range "
        f"[{ion.n_shell[:, -1].min():.0f}, {ion.n_shell[:, -1].max():.0f}] "
        f"(anchored terminal: 14)"
    )
    lines.append("")

    # --- 2. R(t) / |v(t)| vs the 9 A TDDFT trace ----------------------------
    lines.append("[R(t) / |v(t)| vs 9 A TDDFT (compare_* on the overlap window)]")
    comp_R = compare_distance(ion, hedft)
    payload["comp_R"] = comp_R
    lines.append(
        f"  R(t):  RMSE = {comp_R.rmse:.3f} A, mean MD/TDDFT = "
        f"{comp_R.mean_ratio:.3f}, window "
        f"[{comp_R.overlap_t_min_ps:.2f}, {comp_R.overlap_t_max_ps:.2f}] ps "
        f"({comp_R.num_overlap_points} pts)"
    )
    payload["comp_v"] = {}
    for atom in ("I1", "I2"):
        comp_v = compare_velocity_magnitude(ion, hedft, atom=atom)
        payload["comp_v"][atom] = comp_v
        lines.append(
            f"  |v_{atom}|: RMSE = {comp_v.rmse:.3f} A/ps, mean MD/TDDFT = "
            f"{comp_v.mean_ratio:.3f}"
        )
    lines.append("")

    # --- 3. t_x: per-ion reconstruction vs the sharp closed form ------------
    t_x = crossing_time_ps(
        ion.E_int_eV, ion.n_shell, time_ps,
        picture=cfg.ladder_electronic_picture, kappa=cfg.ladder_steepness,
    )
    payload["t_x"] = t_x
    sigma21 = float(
        ladder_cumsum(
            21, picture=cfg.ladder_electronic_picture, kappa=cfg.ladder_steepness
        )
    )
    e0 = cfg.internal_energy_partition_fraction * cfg.coulomb_available_eV
    t_x_closed = cfg.internal_energy_cooling_tau_ps * np.log(e0 / sigma21)
    payload["t_x_closed"] = t_x_closed
    n_nan = int(np.isnan(t_x).sum())
    lines.append("[t_x reconstruction (sharp oracle, plan §2.2)]")
    lines.append(
        f"  Sigma(21) = {sigma21:.4f} eV, E_int(0) = f_int*E_avail = {e0:.3f} eV"
    )
    lines.append(
        f"  closed form tau*ln(E0/Sigma(21)) = {t_x_closed:.3f} ps"
    )
    if n_nan:
        lines.append(f"  WARNING: {n_nan} ions never cross in-window (NaN)")
    finite = t_x[~np.isnan(t_x)]
    if finite.size:
        agree = bool(np.all(finite == finite[0]))
        lines.append(
            f"  reconstructed t_x = {finite[0]:.3f} ps; all ions agree: "
            f"{agree} (spread {finite.max() - finite.min():.4f} ps)"
        )
        lines.append(
            f"  offset vs closed form = {finite[0] - t_x_closed:+.4f} ps "
            f"(dt discreteness bound: one stored step)"
        )
        in_sanity = T_X_SANITY_BAND_PS[0] <= finite[0] <= T_X_SANITY_BAND_PS[1]
        in_prior = GAH25_PRIOR_PS[0] <= finite[0] <= GAH25_PRIOR_PS[1]
        in_prior_x2 = (
            GAH25_PRIOR_PS[0] / 2.0 <= finite[0] <= GAH25_PRIOR_PS[1] * 2.0
        )
        lines.append(
            f"  sanity band {_fmt_band(T_X_SANITY_BAND_PS)}: "
            f"{'inside' if in_sanity else 'OUTSIDE'}; GAH25 prior "
            f"{_fmt_band(GAH25_PRIOR_PS)}: {'inside' if in_prior else 'outside'}"
            f"; prior +/- factor 2: {'inside' if in_prior_x2 else 'OUTSIDE'}"
        )
    lines.append("")

    # --- 4. Pi(t) regime trace (corrected freeze-side expectation) ----------
    pi = regime_parameter(ion, cfg)
    payload["pi"] = pi
    pi_mean = pi.mean(axis=0)
    payload["pi_mean"] = pi_mean
    lines.append("[Pi(t) = lambda(n)*f_ret*tau regime trace]")
    lines.append(
        f"  Pi at gate-open (t=0): {pi_mean[0]:.4f} (expected 0: full-shell "
        f"Langmuir cap)"
    )
    lines.append(
        f"  max mean Pi = {pi_mean.max():.4f} at t = "
        f"{time_ps[int(np.argmax(pi_mean))]:.2f} ps; max per-ion Pi = "
        f"{pi.max():.4f}"
    )
    lines.append(
        f"  terminal mean Pi = {pi_mean[-1]:.4f} (expected -> 0 at exit)"
    )
    lines.append(
        f"  freeze side (Pi < 1 throughout): {bool(np.all(pi < 1.0))} "
        f"(a computed Pi > 1 at this condition is itself a flag)"
    )
    lines.append("")

    # --- 5. 5-term invariant -------------------------------------------------
    closure = ion_ledger_closure(ion)
    payload["closure"] = closure
    rel = closure.max_abs_residual_eV / abs(closure.E_system_eV[0])
    lines.append("[5-term invariant (ion_ledger_closure)]")
    lines.append(
        f"  max |residual| = {closure.max_abs_residual_eV:.3e} eV "
        f"({rel * 100:.4f}% of E_system(0) = {closure.E_system_eV[0]:.4f} eV)"
    )
    return lines, payload


def emit_figures(payload, ion, out_dir: Path) -> list[Path]:
    """Write the three bridge overlay figures; return the saved paths."""
    time_ps = ion.time_ps
    saved: list[Path] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    # (a) mean n(t) + anchored staircase family.
    fig, ax = plt.subplots(figsize=(8, 5))
    for t_star, stair in payload["schedules"].items():
        primary = t_star == PRIMARY_T_STAR_PS
        ax.step(
            time_ps, stair, where="post",
            lw=2.0 if primary else 1.0,
            alpha=1.0 if primary else 0.45,
            label=f"anchored t*={t_star:.1f} ps" + (" (primary)" if primary else ""),
        )
    ax.plot(time_ps, payload["mean_n"], "k-", lw=2.0, label="generative mean n(t)")
    ax.set_xlabel("time [ps]")
    ax.set_ylabel("He-shell count n")
    ax.set_title("Bridge: emergent mean n(t) vs anchored 21→19→14 family")
    ax.legend()
    fig.tight_layout()
    path = out_dir / "bridge_mean_n_overlay.png"
    fig.savefig(path, dpi=150)
    saved.append(path)

    # (b) R(t) and |v(t)| vs TDDFT on the overlap grid.
    fig, (ax_r, ax_v) = plt.subplots(1, 2, figsize=(12, 5))
    comp_R = payload["comp_R"]
    ax_r.plot(comp_R.t_overlap_ps, comp_R.hedft_on_overlap, label="TDDFT R(t)")
    ax_r.plot(comp_R.t_overlap_ps, comp_R.md_on_hedft_grid, label="MD mean R(t)")
    ax_r.set_xlabel("time [ps]")
    ax_r.set_ylabel("I–I distance [Å]")
    ax_r.set_title(f"R(t): RMSE {comp_R.rmse:.2f} Å")
    ax_r.legend()
    for atom, comp_v in payload["comp_v"].items():
        ax_v.plot(comp_v.t_overlap_ps, comp_v.hedft_on_overlap, ls="--",
                  label=f"TDDFT |v_{atom}|")
        ax_v.plot(comp_v.t_overlap_ps, comp_v.md_on_hedft_grid,
                  label=f"MD |v_{atom}|")
    ax_v.set_xlabel("time [ps]")
    ax_v.set_ylabel("|v| [Å/ps]")
    ax_v.set_title("|v(t)| vs 9 Å TDDFT")
    ax_v.legend()
    fig.tight_layout()
    path = out_dir / "bridge_kinematics.png"
    fig.savefig(path, dpi=150)
    saved.append(path)

    # (c) Pi(t) mean +/- envelope, with the t_x marker and the Pi = 1 line.
    fig, ax = plt.subplots(figsize=(8, 5))
    pi = payload["pi"]
    ax.fill_between(
        time_ps, pi.min(axis=0), pi.max(axis=0), alpha=0.3,
        label="per-ion envelope",
    )
    ax.plot(time_ps, payload["pi_mean"], lw=2.0, label="mean Π(t)")
    ax.axhline(1.0, color="r", ls=":", label="Π = 1 (freeze/shed boundary)")
    t_x = payload["t_x"]
    finite = t_x[~np.isnan(t_x)]
    if finite.size:
        ax.axvline(finite[0], color="k", ls="--", lw=1.0,
                   label=f"t× = {finite[0]:.2f} ps")
    ax.set_xlabel("time [ps]")
    ax.set_ylabel("Π = λ(n)·f_ret·τ")
    ax.set_title("Bridge: regime order parameter (freeze side expected)")
    ax.legend()
    fig.tight_layout()
    path = out_dir / "bridge_regime_pi.png"
    fig.savefig(path, dpi=150)
    saved.append(path)

    if SHOW_FIGURES:
        plt.show()
    plt.close("all")
    return saved


def main() -> int:
    run_dir = (
        Path(RUN_DIR)
        if RUN_DIR is not None
        else PROJECT_ROOT / "data" / "runs"
        / tier2_bridge_run_dir_name(CASE, VARIANT, N)
    )
    run = RunDirectory(run_dir)
    cfg = run.load_cfg()
    ion = run.load_ion(cfg)
    hedft = load_hedft_trajectory(HEDFT_CSV_PATH)

    lines, payload = build_summary_lines(cfg, ion, hedft)

    if SAVE_FIGURES:
        saved = emit_figures(payload, ion, run_dir / "figures")
        lines.append("")
        lines.append("[figures]")
        lines.extend(f"  {p}" for p in saved)

    text = "\n".join(lines)
    print(text)
    (run_dir / "bridge_summary.txt").write_text(text, encoding="utf-8")
    print(f"\nsummary written -> {run_dir / 'bridge_summary.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
