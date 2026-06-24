"""Plot the Tier-1a anchored He-shell schedule n_bar(t) with its shed events.

Schedule diagnostic (not a run post-processor): draws the piecewise-linear shell
count ``n_bar(t)`` from :func:`i2_helium_md.physics.shell_schedule.build_shell_schedule`
for each swept onset ``t*`` over ``t in [0, 20]`` ps, marking the 7 cold-shed
events. It visualizes directly that the two segment-1 sheds move with ``t*`` while
the five segment-2 sheds ``{10.4..13.6} ps`` are ``t*``-independent
(``TIER1A_IMPLEMENTATION_PLAN.md`` Sec.4 / Sec.10). It reads no run directory and
no ``SimConfig``; it lives here only because ``post_processing/`` owns the plot
entry points.

Run with::

    python scripts/post_processing/plot_shell_schedule.py

Output is interactive (``plt.show()``); a PNG copy is also written next to this
script's configured output path.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
T_STAR_VALUES_PS = (0.5, 5.0, 9.0)       # Tier-1a onset sweep
TIME_RANGE_PS = (0.0, 20.0)              # plot window
OUTPUT_PNG = PROJECT_ROOT / "data" / "figures" / "shell_schedule.png"


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from i2_helium_md.physics.shell_schedule import build_shell_schedule  # noqa: E402


def main() -> int:
    schedules = {t_star: build_shell_schedule(t_star) for t_star in T_STAR_VALUES_PS}
    for t_star, sched in schedules.items():
        times = ", ".join(f"{e.time_ps:.2f}" for e in sched.events)
        print(f"t*={t_star:>4} ps : 7 sheds at [{times}] ps")

    fig = _build_figure(schedules)
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=150)
    print(f"wrote {OUTPUT_PNG}")

    plt.show()
    return 0


def _build_figure(schedules) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
    t = np.linspace(TIME_RANGE_PS[0], TIME_RANGE_PS[1], 2001)

    for t_star, sched in schedules.items():
        # Physical integer shell count n(t) -- the staircase (foregrounded).
        (line,) = ax.plot(
            t, sched.n_of_t(t), linewidth=1.8,
            drawstyle="steps-post", label=rf"$n(t)$, $t^*={t_star}$ ps",
        )
        color = line.get_color()
        # The continuous anchored loss curve n_bar(t) -- thin guide whose
        # half-integer crossings define the steps (NOT the physical count).
        ax.plot(t, sched.n_bar(t), linewidth=0.8, linestyle="--", alpha=0.35, color=color)
        # Mark each shed at its fire time (n drops n_before -> n_after).
        ev_t = [e.time_ps for e in sched.events]
        ev_n = [e.n_after for e in sched.events]
        ax.scatter(ev_t, ev_n, s=24, color=color, zorder=3)

    ax.set_title(r"Tier-1a anchored He-shell count $n(t)$ (staircase) and $\bar{n}(t)$ guide")
    ax.set_xlabel("t / ps")
    ax.set_ylabel("shell count")
    ax.set_xlim(*TIME_RANGE_PS)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig


if __name__ == "__main__":
    raise SystemExit(main())
