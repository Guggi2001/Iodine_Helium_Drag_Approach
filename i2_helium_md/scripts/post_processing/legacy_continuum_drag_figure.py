"""Thesis figure: the legacy hard-sphere collision model is a continuum drag law.

Thesis section ``sec:legacy_continuum`` (figure ``fig:hard_sphere_is_continuum``,
``Pictures/fig_continuum_drag.png``). Two panels, one curve per legacy preset
(9 A: sigma_0 = 2500 A^2, 18 A: sigma_0 = 1600 A^2, exponent e = -2):

(a) Mean free path of the Mode-3 sampling kernel,
    ``lambda(v) = 1 / (n sigma(v)) = v^2 / (n sigma_0)``, against the He
    nearest-neighbour spacing ``n^(-1/3)``. Over the probed speeds lambda stays
    well below one He spacing: the "collisions" are not resolved encounters.
(b) The resulting stopping force. Measured per speed bin from the committed
    legacy runs as ``F_eff = sum(dE_dissip) / sum(path)`` over in-droplet steps
    (error bar: cluster standard error over molecules), against the closed form
    of the kernel,
    ``F = n sigma_0 * (m v^2 / 2) * 2 rho / (1 + rho)^2 * v^-2``
    (rho = m_ion / m_He), which is independent of v for e = -2.

Both panels shade the path-weighted (in-droplet) speed density of each run, i.e.
the speed range the legacy ion stage actually probes.

Inputs: ``data/runs/{9A,18A}_hedft_comparison`` (legacy hard-sphere ion runs);
their ``cfg.json`` must match the presets ``single_pulse_N2000`` /
``single_pulse_N2000_18Angst`` (checked, fails loudly).

Writes ``data/figures/thesis_legacy/fig_continuum_drag.{png,pdf}`` and the
plotted points ``fig_continuum_drag_points.csv``.

Usage::

    python scripts/post_processing/legacy_continuum_drag_figure.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
OUT_DIR = PROJECT_ROOT / "data" / "figures" / "thesis_legacy"
SHOW_FIGURE = False
SPEED_EDGES_APS = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0])
MIN_BIN_PATH_A = 100.0      # drop F_eff bins with less in-droplet path than this [A]
DENSITY_BINS = 120          # bins of the shaded path-weighted speed density
DENSITY_HEIGHT = 0.30       # shaded density peak, as a fraction of the axes height

scaling_factor = 0.9

import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

from i2_helium_md.physics.constants import DENSITY_DROPLET, EV, U  # noqa: E402
from i2_helium_md.presets import single_pulse_N2000, single_pulse_N2000_18Angst  # noqa: E402

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
STEM = "fig_continuum_drag"
#: eV/A -> amu A/ps^2
EV_PER_A_TO_AMU_A_PS2 = (EV / 1e-10) / (U * 1e-10 / 1e-24)
_PLASMA = matplotlib.colormaps["plasma"]
#: label -> (preset, run directory, colour, marker); colours sampled from plasma
CASES = {
    "9 Å": (single_pulse_N2000(), "9A_hedft_comparison", _PLASMA(0.15), "o"),
    "18 Å": (single_pulse_N2000_18Angst(), "18A_hedft_comparison", _PLASMA(0.68), "s"),
}
BG, AXIS, TEXT = "white", "#b5b2aa", "#4a4a4a"


def _check_run_matches_preset(run_dir: Path, preset) -> None:
    """The committed run must carry the preset's collision settings."""
    cfg = json.loads((run_dir / "cfg.json").read_text(encoding="utf-8"))
    for key in ("geometric_scattering_crosssection_Iplus", "sigma_ion_exponent",
                "scatter_mass_ion_amu", "hard_sphere_collision_mode", "R0_GS_angstrom"):
        if cfg[key] != getattr(preset, key):
            raise AssertionError(f"{run_dir.name}: {key} {cfg[key]!r} != preset {getattr(preset, key)!r}")


def measure_run(run_dir: Path) -> dict:
    """Per-step in-droplet path, dissipated energy and pre-step speed of one legacy run."""
    d = np.load(run_dir / "ion.npz")
    x, y, z = d["positions_x"], d["positions_y"], d["positions_z"]
    ds = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2 + np.diff(z) ** 2)
    inside = (np.sqrt(x[:, 1:] ** 2 + y[:, 1:] ** 2 + z[:, 1:] ** 2)
              - d["droplet_radii_angstrom"][:, None]) < 0.0
    del x, y, z
    speed = np.sqrt(d["velocities_x"] ** 2 + d["velocities_y"] ** 2
                    + d["velocities_z"] ** 2)[:, :-1]
    dE = np.diff(d["E_dissip_eV"])
    n_mol = int(d["num_molecules"])
    # Atoms i and i + N belong to the same molecule (one cluster for the SE).
    molecule = np.broadcast_to((np.arange(2 * n_mol) % n_mol)[:, None], ds.shape)
    m_ion = float(d["mass_kg"][0] / U)
    return dict(ds=ds[inside], dE=dE[inside], speed=speed[inside],
                molecule=molecule[inside], n_mol=n_mol, m_ion_amu=m_ion)


def binned_force(m: dict) -> list[dict]:
    """F_eff = sum(dE) / sum(path) per speed bin, with a cluster SE over molecules."""
    idx = np.digitize(m["speed"], SPEED_EDGES_APS) - 1
    rows = []
    for b in range(len(SPEED_EDGES_APS) - 1):
        k = idx == b
        L = m["ds"][k].sum()
        if L < MIN_BIN_PATH_A:
            continue
        E_mol = np.bincount(m["molecule"][k], weights=m["dE"][k], minlength=m["n_mol"])
        L_mol = np.bincount(m["molecule"][k], weights=m["ds"][k], minlength=m["n_mol"])
        F = E_mol.sum() / L
        # Ratio-estimator (delta-method) standard error over molecule clusters.
        resid = E_mol - F * L_mol
        se = np.sqrt(m["n_mol"] / (m["n_mol"] - 1) * np.sum(resid ** 2)) / L
        rows.append(dict(v_lo=SPEED_EDGES_APS[b], v_hi=SPEED_EDGES_APS[b + 1],
                         v_mid=0.5 * (SPEED_EDGES_APS[b] + SPEED_EDGES_APS[b + 1]),
                         path_A=L, F=F * EV_PER_A_TO_AMU_A_PS2,
                         F_se=se * EV_PER_A_TO_AMU_A_PS2))
    return rows


def closed_form_force(preset, m_ion_amu: float) -> float:
    """Kernel stopping force n sigma_0 <dE>/E (m/2) for e = -2 [amu A/ps^2]."""
    if preset.sigma_ion_exponent != -2.0:
        raise ValueError("closed form is velocity independent only for e = -2")
    rho = m_ion_amu / preset.scatter_mass_ion_amu
    return (DENSITY_DROPLET * preset.geometric_scattering_crosssection_Iplus
            * 0.5 * m_ion_amu * 2.0 * rho / (1.0 + rho) ** 2)


def _shade_density(ax, m: dict, color: str) -> None:
    """Path-weighted speed density, peak at DENSITY_HEIGHT of the axes height."""
    h, e = np.histogram(m["speed"], bins=DENSITY_BINS,
                        range=(SPEED_EDGES_APS[0], SPEED_EDGES_APS[-1]), weights=m["ds"])
    c = 0.5 * (e[1:] + e[:-1])
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    y = DENSITY_HEIGHT * h / h.max()
    ax.fill_between(c, 0.0, y, color=color, alpha=0.15, lw=0, transform=trans, zorder=0)
    #ax.plot(c, y, color=color, alpha=0.55, lw=1.0, transform=trans, zorder=0)


def _style(ax) -> None:
    ax.set_facecolor(BG)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=TEXT, labelsize=11)
    ax.grid(axis="y", color="#dddbd5", lw=0.8)
    ax.set_axisbelow(True)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("===== legacy hard-sphere model: continuum drag figure =====")
    data = {}
    for label, (preset, run_name, color, marker) in CASES.items():
        run_dir = RUNS_ROOT / run_name
        _check_run_matches_preset(run_dir, preset)
        m = measure_run(run_dir)
        F_cf = closed_form_force(preset, m["m_ion_amu"])
        rows = binned_force(m)
        data[label] = dict(preset=preset, m=m, rows=rows, F_cf=F_cf, color=color, marker=marker)
        print(f"  {label}: sigma_0 {preset.geometric_scattering_crosssection_Iplus:g} A^2, "
              f"closed form {F_cf:.1f}, measured "
              + ", ".join(f"{r['v_mid']:.2f}:{r['F']:.1f}" for r in rows))

    fig, (axL, axF) = plt.subplots(1, 2, figsize=(10.0*scaling_factor, 4.95*scaling_factor))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.875, bottom=0.28, wspace=0.28)
    spacing = DENSITY_DROPLET ** (-1.0 / 3.0)
    v = np.linspace(SPEED_EDGES_APS[0], SPEED_EDGES_APS[-1], 400)

    # (a) mean free path of the sampling kernel
    _style(axL)
    axL.set_yscale("log")
    for label, d in data.items():
        _shade_density(axL, d["m"], d["color"])
        lam = v ** 2 / (DENSITY_DROPLET * d["preset"].geometric_scattering_crosssection_Iplus)
        axL.plot(v, lam, color=d["color"], lw=2.6)
    axL.axhline(spacing, color="#808080", ls="--", lw=1.8, dashes=(4, 3))
    axL.text(0.65, spacing * 1.12, rf"He spacing  $n^{{-1/3}}$ = {spacing:.2f} Å",
             color=TEXT, fontsize=11, va="bottom")
    s18 = data["18 Å"]["preset"].geometric_scattering_crosssection_Iplus
    s9 = data["9 Å"]["preset"].geometric_scattering_crosssection_Iplus
    axL.text(7.0, 7.0 ** 2 / (DENSITY_DROPLET * s18) * 1.12, "18 Å", fontsize=12, va="bottom")
    axL.text(8.9, 8.9 ** 2 / (DENSITY_DROPLET * s9) / 1.2, "9 Å", fontsize=12,
             ha="right", va="top")
    axL.set(xlim=(SPEED_EDGES_APS[0], SPEED_EDGES_APS[-1]), ylim=(8e-3, 9.0))
    axL.set_xlabel(r"Ion speed  $v$ / $\mathrm{\AA\,ps^{-1}}$", fontsize=12.5)
    axL.set_ylabel(r"Mean free path  $\lambda$  / Å", fontsize=12.5)
    axL.set_title("(a) Mean free path",
                  fontsize=13, loc="left", x=0.0)

    # (b) resulting stopping force
    _style(axF)
    for label, d in data.items():
        _shade_density(axF, d["m"], d["color"])
        axF.axhline(d["F_cf"], color=d["color"], ls="--", lw=2.2, dashes=(5, 3), alpha=0.9)
        axF.text(SPEED_EDGES_APS[-1], d["F_cf"] + 6, f"{d['F_cf']:.0f}", color=TEXT,
                 fontsize=11, ha="right", va="bottom")
        r = d["rows"]
        axF.errorbar([x["v_mid"] for x in r], [x["F"] for x in r],
                     yerr=[x["F_se"] for x in r], color=d["color"], marker=d["marker"],
                     ms=6, ls="none", capsize=4, elinewidth=1.6, capthick=1.6, zorder=3)
    axF.set(xlim=(SPEED_EDGES_APS[0], SPEED_EDGES_APS[-1]),
            ylim=(0, 1.15 * max(d["F_cf"] for d in data.values())))
    axF.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    axF.set_xlabel(r"Ion speed  $v$ / $\mathrm{\AA\,ps^{-1}}$", fontsize=12.5)
    axF.set_ylabel(r"$F_{\rm eff}$  / amu Å ps$^{-2}$", fontsize=12.5)
    axF.set_title("(b) Stopping force",
                  fontsize=13, loc="left", x=0.0)

    handles = [Line2D([], [], color=d["color"], lw=2.6, ls="--", marker=d["marker"], ms=6,
                      label=f"{label} preset  "
                      rf"($\sigma_0$ = {d['preset'].geometric_scattering_crosssection_Iplus:.0f} Å$^2$)")
               for label, d in data.items()]
    handles += [Line2D([], [], color="#808080", ls="--", lw=2.2, dashes=(4, 3),
                       label=r"Analytic mean force  $F = n\,\sigma_0\,\mu^2/m_{\rm He}$"),
                Patch(facecolor="#d9d9d9", edgecolor="#a6a6a6", lw=1.0, label="Path-weighted speed density (shaded)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=11, labelcolor=TEXT, columnspacing=3.0, bbox_to_anchor=(0.5, 0.01))

    png, pdf = OUT_DIR / f"{STEM}.png", OUT_DIR / f"{STEM}.pdf"
    fig.savefig(png, dpi=200, facecolor=BG)
    fig.savefig(pdf, facecolor=BG)
    with open(OUT_DIR / f"{STEM}_points.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(("preset", "sigma0_A2", "v_lo_Aps", "v_hi_Aps", "path_A",
                    "F_eff_amu_A_ps2", "F_eff_se", "F_closed_form"))
        for label, d in data.items():
            for r in d["rows"]:
                w.writerow((label.replace(" Å", "A"),
                            d["preset"].geometric_scattering_crosssection_Iplus,
                            r["v_lo"], r["v_hi"], r["path_A"], r["F"], r["F_se"], d["F_cf"]))
    for p in (png, pdf):
        print(f"  wrote {p.relative_to(PROJECT_ROOT).as_posix()}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
