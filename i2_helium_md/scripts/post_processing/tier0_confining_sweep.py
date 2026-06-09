"""Confining-potential parameter sweep for the t*-seeded Tier-0 harness.

Studies how the droplet **confining potential** shapes the seeded I+ trajectory.
The ion-stage confining force is
``droplet_force(depth, steepness=cfg.potential_steepness,
binding_energy=cfg.binding_energy_I_ion_eV)``
(``physics/leapfrog.py`` ``_droplet_acceleration`` -> ``physics/potentials.py``).
Both knobs are config fields, already overridable via the preset builder
(``presets.py`` ends each builder with ``replace(cfg, **overrides)``), so this
script just loops over a small set of overrides and overlays the resulting
trajectories against the TDDFT reference.

Two facts make the overlay clean:

* The seed (``R``, ``|v|`` at t*) is parameter-independent, so every sweep curve
  starts at the **same** t* state -- divergence afterward is the pure influence
  of the confining parameter.
* ``potential_steepness`` is independent of the drag gate
  (``drag_gate_steepness`` is a separate field), so sweeping the confining width
  does not move the drag region. The confining influence is isolated.

The reduced MD series (``R(t)``, ``|v2|(t)``) reuse
``tier0_drag_comparison.ensemble_mean_series`` -- the same molecule-means the
single-run figure plots -- and the forward integration reuses
``tier0_tstar_seeded_comparison.run_case`` with its new ``overrides=`` argument.

How to use
----------
Edit USER SETTINGS, then::

    python scripts/post_processing/tier0_confining_sweep.py
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
CASE = "9A"  # "9A" or "18A"
# A from-onset run for the same case, used only to read the droplet radius the
# real pipeline assigns to this droplet size (so the spatial gate matches).
ONSET_RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_drag_tier0_N50"
# Each entry: label -> SimConfig field overrides applied on top of the preset.
# "baseline" (empty override) reproduces the committed single-run t*-seeded run.
SWEEP = {
    "baseline":   {},
    "deep well":  {"binding_energy_I_ion_eV": 0.6},
    "shallow":    {"binding_energy_I_ion_eV": 0.15},
    "stiff wall": {"potential_steepness": 8.0},
}
SHOW_FIGURE = True


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# This script's own directory is on sys.path when run directly, so the sibling
# harness modules import cleanly.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import tier0_drag_comparison as T  # noqa: E402
import tier0_tstar_seeded_comparison as TS  # noqa: E402


def run_sweep(case, onset_run_dir, sweep):
    """Run one t*-seeded trajectory per override set; collect series + metrics.

    Returns a list of ``dict`` records with the label, the override applied, the
    MD series ``(t, distance, |v2|)``, the reference ``hedft``, the scored
    ``window``, and the ``metrics`` from ``T.score``.
    """
    records = []
    hedft = None
    window = None
    for label, overrides in sweep.items():
        print(f"\n----- sweep: {label}  overrides={overrides or '{}'} -----")
        ion, hedft, window, _drop, _cfg, metrics = TS.run_case(
            case, onset_run_dir, overrides=overrides,
        )
        t_md, dist_md, _v1_md, v2_md = T.ensemble_mean_series(ion)
        records.append({
            "label": label,
            "overrides": overrides,
            "t": t_md,
            "distance": dist_md,
            "v2": v2_md,
            "metrics": metrics,
        })
    return records, hedft, window


def print_sweep_table(records):
    """Print in-window scored numbers per sweep entry (quantitative influence)."""
    print("\n===== confining sweep: in-window RMSE vs reference =====")
    print(f"  {'label':<12} {'dist RMSE / A':>14} {'|v2| RMSE / (A/ps)':>20}")
    for r in records:
        m = r["metrics"]
        print(f"  {r['label']:<12} {m['distance_rmse_A']:>14.4f} "
              f"{m['v_I2_rmse_Aps']:>20.4f}")
    print("=" * 52)


def build_sweep_figure(records, hedft, window):
    """Overlay R(t) (top) and |v2|(t) (bottom) for every sweep entry vs reference."""
    import matplotlib.pyplot as plt

    t_start, t_end = window
    fig, (ax_d, ax_v) = plt.subplots(
        2, 1, figsize=(8.0, 6.0), sharex=True, constrained_layout=True,
    )

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, r in enumerate(records):
        c = colors[i % len(colors)]
        ax_d.plot(r["t"], r["distance"], color=c, lw=1.5, label=r["label"])
        ax_v.plot(r["t"], r["v2"], color=c, lw=1.5, label=r["label"])

    # Reference (black) + shaded scored window, matching build_figure conventions.
    ax_d.plot(hedft.time_ps, hedft.distance_A, color="black", lw=1.5,
              label="HeDFT / TDDFT")
    ax_d.axvspan(t_start, t_end, color="tab:green", alpha=0.12,
                 label="scored window")
    ax_d.set_ylabel(r"$R_1 - R_2$ / $\mathrm{\AA}$")
    ax_d.legend(frameon=False, ncol=2)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)

    ax_v.plot(hedft.time_ps, hedft.v2_magnitude_Aps, color="black", lw=1.5,
              ls=":", label="HeDFT |v2|")
    ax_v.axvspan(t_start, t_end, color="tab:green", alpha=0.12)
    ax_v.set_ylabel(r"$|v_2|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax_v.set_xlabel("t / ps")
    ax_v.legend(frameon=False, ncol=2)
    ax_v.spines["top"].set_visible(False)
    ax_v.spines["right"].set_visible(False)

    fig.suptitle(f"{CASE} t*-seeded -- confining-potential sweep")
    return fig


def main() -> int:
    # Keep the sibling harness's CASE-dependent state (titles, etc.) in sync.
    TS.CASE = CASE
    records, hedft, window = run_sweep(CASE, ONSET_RUN_DIR, SWEEP)
    print_sweep_table(records)
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        build_sweep_figure(records, hedft, window)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
