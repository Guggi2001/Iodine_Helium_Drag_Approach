"""Tier-0 same-smoothed consistency comparison.

The one remaining Tier-0 task (``tier0_comparison_tasks_left.md`` §4.1,
``TIER0_FINDINGS.md``). Tier 0 is an *internal-consistency check*: does the
drag implementation + BAOAB driver forward-integrate the extracted gamma back
into the (denoised) trajectory gamma was fit to?

The 9 A residual (~0.86 A/ps in the clean band) is hypothesised to be the
1.2 ps **bubble-mode oscillation** present in the *raw* ``9A_All_Data.csv`` but
removed by the extraction's CEEMDAN+SG before gamma was fit. This script scores
the MD clean-atom speed |v2| against BOTH:

  * the **raw** reference |v2| (``compare_velocity_magnitude``), and
  * the **same-smoothed** reference |v2| (the CEEMDAN+SG-denoised
    ``data/reference/drag/<case>/velocity_smoothed/cleaned_data.csv``),

over the **cleaned reference's own window** (9 A -> [2.67, 6.0],
18 A -> [4.54, 8.0]), for BOTH comparison modes:

  * **from-onset** -- the production-representative MD run, integrated from the
    ion-stage explosion onset (loaded from a finished drag run directory);
  * **t\\*-seeded** -- the clean-form diagnostic that seeds the ion at t* with
    the reference state (reuses ``tier0_tstar_seeded_comparison.run_case``).

It reports same-smoothed *alongside* raw, never instead of it (the raw number is
the harsher honest cross-check). If the same-smoothed RMSE collapses to
~18 A-class while the raw stays elevated, the residual was bubble-mode variance
the law never claimed to reproduce -> Tier 0 is a clean consistency pass.

IMPORTANT (provenance): the 9 A drag law was re-extracted on the clean window
([2.67, 6.0], a=24.88/b=2.085 in fit_parameters.json). Any from-onset run scored
here must be **regenerated with the current gamma** (``scripts/gen_tier0_runs.py
9A 50``); the committed ``md_mean_trajectory_N50.csv`` predates the
re-extraction and must not be reused. The t*-seeded mode runs the preset inline,
so it is automatically on the current gamma.

How to use
----------
Regenerate the current-gamma runs first::

    python scripts/gen_tier0_runs.py 9A 50
    python scripts/gen_tier0_runs.py 18A 50

then edit USER SETTINGS and::

    python scripts/post_processing/tier0_same_smoothed_comparison.py
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
CASE = "9A"  # "9A" or "18A"

# A current-gamma from-onset drag run for this case (regenerated via
# scripts/gen_tier0_runs.py). Also used by the t*-seeded mode to read the
# droplet radius the real pipeline assigns.
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "9A_drag_tier0_N50"

SHOW_FIGURE = True


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
# This script's own directory is on sys.path when run directly, so the sibling
# Tier-0 harness modules import cleanly.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import tier0_drag_comparison as T  # noqa: E402
import tier0_tstar_seeded_comparison as TS  # noqa: E402
from i2_helium_md.postprocess.compare_trajectories import (  # noqa: E402
    compare_speed_to_reference,
    compare_velocity_magnitude,
)
from i2_helium_md.postprocess.hedft_loader import (  # noqa: E402
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _paths(case: str):
    raw = PROJECT_ROOT / "data" / "reference" / f"{case}_All_Data.csv"
    cleaned = (
        PROJECT_ROOT / "data" / "reference" / "drag" / case
        / "velocity_smoothed" / "cleaned_data.csv"
    )
    coeff = PROJECT_ROOT / "data" / "reference" / "drag" / case / "linear_and_cubic"
    return raw, cleaned, coeff


def score_mode(ion, hedft, cleaned, window) -> dict:
    """Score one MD trajectory against raw |v2| and same-smoothed |v2|.

    Returns the clean-atom (I2) raw and smoothed RMSE plus the mean(I1,I2) raw
    number (continuity with prior findings) -- all on the same window.
    """
    raw_i2 = compare_velocity_magnitude(ion, hedft, atom="I2", window=window)
    raw_i1 = compare_velocity_magnitude(ion, hedft, atom="I1", window=window)
    smoothed_i2 = compare_speed_to_reference(
        ion,
        atom="I2",
        t_ref_ps=cleaned.time_ps,
        ref_speed_Aps=cleaned.speed_Aps,
        window=window,
    )
    return {
        "n_raw": raw_i2.num_overlap_points,
        "n_smoothed": smoothed_i2.num_overlap_points,
        "raw_i2_rmse": raw_i2.rmse,
        "smoothed_i2_rmse": smoothed_i2.rmse,
        "raw_mean_rmse": 0.5 * (raw_i1.rmse + raw_i2.rmse),
    }


def print_block(label: str, m: dict) -> None:
    delta = m["raw_i2_rmse"] - m["smoothed_i2_rmse"]
    print(f"  [{label}]")
    print(f"    raw |v2| RMSE       : {m['raw_i2_rmse']:.4f} A/ps "
          f"({m['n_raw']} raw samples)")
    print(f"    smoothed |v2| RMSE  : {m['smoothed_i2_rmse']:.4f} A/ps "
          f"({m['n_smoothed']} cleaned samples)")
    print(f"    raw - smoothed      : {delta:+.4f} A/ps")
    print(f"    raw mean(I1,I2) RMSE: {m['raw_mean_rmse']:.4f} A/ps "
          "(continuity diagnostic)")


def run_case(case: str, run_dir: Path):
    raw_path, cleaned_path, coeff = _paths(case)
    hedft = load_hedft_trajectory(raw_path)
    cleaned = load_smoothed_speed_reference(cleaned_path)
    window = (float(cleaned.time_ps[0]), float(cleaned.time_ps[-1]))

    t_start, t_end, meff = T.read_drag_window(coeff)
    print()
    print(f"===== Tier-0 same-smoothed comparison: {case} =====")
    print(f"  fit_parameters window : [{t_start:.3f}, {t_end:.3f}] ps, "
          f"m_eff={meff:.4f} amu")
    print(f"  cleaned-ref window    : [{window[0]:.3f}, {window[1]:.3f}] ps "
          f"(scored window)")
    if abs(window[0] - t_start) > 1e-6 or window[1] > t_end + 1e-6:
        print("  NOTE: cleaned-ref window differs from fit_parameters window; "
              "scoring on the cleaned extent.")

    # from-onset mode -- the production-representative MD run.
    ion_onset = RunDirectory(run_dir).load_ion()
    print(f"  from-onset run        : {run_dir.name}, N={ion_onset.num_molecules}, "
          f"t=[{ion_onset.time_ps[0]:.2f}, {ion_onset.time_ps[-1]:.2f}] ps")
    m_onset = score_mode(ion_onset, hedft, cleaned, window)

    # t*-seeded mode -- the clean-form diagnostic (current gamma inline).
    ion_seeded, _, _, _, _ = TS.run_case(case, run_dir)
    m_seeded = score_mode(ion_seeded, hedft, cleaned, window)

    print()
    print_block("from-onset", m_onset)
    print_block("t*-seeded", m_seeded)
    print("=" * (38 + len(case)))
    return cleaned, hedft, window, ion_onset, ion_seeded


def build_figure(case, cleaned, hedft, window, ion_onset, ion_seeded):
    import matplotlib.pyplot as plt

    def md_v2(ion):
        n = ion.num_molecules
        speed = np.sqrt(
            ion.velocities_x[n:] ** 2
            + ion.velocities_y[n:] ** 2
            + ion.velocities_z[n:] ** 2
        )
        return np.asarray(ion.time_ps, float), np.mean(speed, axis=0)

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(hedft.time_ps, hedft.v2_magnitude_Aps, color="0.7", lw=1.0,
            label="raw |v2| (TDDFT)")
    ax.plot(cleaned.time_ps, cleaned.speed_Aps, color="black", lw=1.8,
            label="same-smoothed |v2| (CEEMDAN+SG)")
    t_o, v_o = md_v2(ion_onset)
    t_s, v_s = md_v2(ion_seeded)
    ax.plot(t_o, v_o, color="tab:blue", lw=1.2, label="MD |v2| from-onset")
    ax.plot(t_s, v_s, color="tab:red", lw=1.2, ls="--", label="MD |v2| t*-seeded")
    ax.axvspan(window[0], window[1], color="tab:green", alpha=0.10,
               label="scored window")
    ax.set_xlim(window[0] - 1.0, window[1] + 1.0)
    ax.set_xlabel("t / ps")
    ax.set_ylabel(r"$|v_2|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax.set_title(f"Tier-0 same-smoothed comparison ({case})")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def main() -> int:
    cleaned, hedft, window, ion_onset, ion_seeded = run_case(CASE, RUN_DIR)
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt
        build_figure(CASE, cleaned, hedft, window, ion_onset, ion_seeded)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
