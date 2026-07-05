"""Emit the Tier-2 identifiability + regime-determination report (Slice F4).

This is the campaign's **central deliverable**: it turns the F3 size-distribution
scoreboard into a **per-quantity identifiability statement**, not a point fit. The
single observable (the terminal I+He_n size distribution) is asked to carry the
whole CALIBRATION_MAP "Tier-2 load-bearing" set -- 8+ quantities -- and F4 reports,
quantity by quantity, whether that observable can actually *separate* them:

* **ladder shape (kappa)** -- the W1-vs-kappa co-fit landscape per electronic
  picture; identifiable only if the landscape has a well-separated *interior*
  minimum. An argmin at a grid edge is reported ``not_bracketed`` (the optimum is
  not enclosed -> extend the grid / the RRK-dof mechanism-level OQ), per the
  Phase-D bridge finding that the under-shed is kinetic and pushes toward the
  sharp-cliff end.
* **electronic picture** -- separation of the best-W1-per-picture across the
  co-fit; identifiable only if the pictures separate beyond a flatness floor.
* **tau** -- the post-crossing tail: terminal-n *insensitivity* is confirmed
  (``insensitive``) or flagged (``sensitive``, a live calibration dimension) over
  the Stage-2 tau sweep. The Phase-D bridge expects the sensitive arm at 0.80 eV.
* **terminal regime** -- the shell-retaining <-> total-strip axis (E5 t*/Pi
  spine): which regime the campaign realizes and which regime's runs score the
  lowest W1.

Reported deterministically (no computation -- fixed by scope / mechanism):

* **f_ret** -- ``not_identifiable``: the 9 Aa-only campaign (2026-07-03 scope
  decision) has no 18 Aa density-contrast leg, the only route to f_ret. Held at
  its prior.
* **lambda_attach** -- ``not_identifiable``: at 9 Aa / in-window, pickup is
  structurally dead (Pi <= 0.005, Phase-D bridge finding #4); held fixed.
* **f_int** -- ``not_identifiable`` by construction: a timing-only knob (moves t*,
  not the cascade budget), so the terminal-n observable is insensitive to it.
* **nu, s** -- ``held_fixed`` (RRK prefactor / dof not swept this campaign); s
  carries the s<->kappa coupling caveat (a kappa fit is conditional on the fixed
  RRK dof).

**Reported, not auto-adjudicated** (the Tier-1a reporting-gate stance): no code
asserts a fidelity pass/fail. Reporting a factual W1 *minimum* (argmin kappa /
best picture / best-W1 regime) is a summary, not a verdict; the flatness /
bracketing / sensitivity thresholds below are reader aids, not acceptance gates.

The report is **sectioned per budget** (0.80 eV validation vs 2.70 eV
production carry different identifiability meaning; at F4's first run only 0.80
exists, 2.70 arrives with F5). Input is the F3 typed records
(``collect_tier2_records``), never the lossy CSV.

Usage::

    python scripts/post_processing/tier2_identifiability_report.py
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any
import warnings


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# USER SETTINGS
# =============================================================================

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
REFERENCE_PATH = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

# Report only runs at this scenario budget [eV] (None -> every complete tier2 run,
# sectioned per budget).
BUDGET_EV: float | None = None

SAVE_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier2_identifiability.csv"

# --- Reader-aid thresholds (NOT acceptance gates; reported, not adjudicated) ---

# tau is a "live" (sensitive) dimension when terminal n moves by more than this
# many He atoms across the swept tau band; below it, terminal-n insensitivity is
# confirmed. Half a He atom of terminal-mean movement is the smallest shift worth
# calling physical for an integer-n size distribution.
TAU_SENSITIVITY_THRESHOLD_HE = 0.5

# A W1-vs-kappa landscape is "flat" (degenerate -> kappa not identifiable) when
# its relative W1 range (max-min)/min falls below this. 10% keeps a real but
# shallow minimum identifiable while rejecting numerical noise.
KAPPA_FLAT_FRACTION = 0.10

# Electronic pictures are "separated" when the relative spread of their
# best-over-kappa W1 exceeds this fraction of the best. Same 10% floor.
PICTURE_SEPARATION_FRACTION = 0.10


from i2_helium_md.physics.constants import N_STAR  # noqa: E402

from scripts.post_processing.tier2_size_distribution_table import (  # noqa: E402
    collect_tier2_records,
)


# --- status vocabulary -------------------------------------------------------
STATUS_IDENTIFIABLE = "identifiable"
STATUS_NOT_IDENTIFIABLE = "not_identifiable"
STATUS_NOT_BRACKETED = "not_bracketed"
STATUS_INSUFFICIENT = "insufficient_data"
STATUS_HELD_FIXED = "held_fixed"
STATUS_SENSITIVE = "sensitive"
STATUS_INSENSITIVE = "insensitive"

# Statuses that mean "the observable pins this quantity" for the headline count.
_SEPARATED_STATUSES = frozenset({STATUS_IDENTIFIABLE, STATUS_SENSITIVE})

# The CALIBRATION_MAP "Tier-2 load-bearing" set, computed quantities first.
LOAD_BEARING_QUANTITIES = [
    "ladder_shape_kappa",
    "electronic_picture",
    "tau",
    "terminal_regime",
    "f_int",
    "f_ret",
    "lambda_attach",
    "nu",
    "s",
]

IDENTIFIABILITY_COLUMNS = [
    "budget_eV",
    "quantity",
    "status",
    "value_or_range",
    "evidence",
    "notes",
]


# ---------------------------------------------------------------------------
# grouping dataclasses
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class KappaLandscape:
    """W1-vs-kappa at one fixed (budget, picture, f_int, f_ret, tau) cell."""

    budget_eV: float
    picture: str
    tau_ps: float | None
    points: tuple[tuple[float, float], ...]  # sorted (kappa, W1_matched)

    @property
    def w1_min(self) -> float:
        return min(w for _, w in self.points)


@dataclass(frozen=True)
class TauSweep:
    """terminal-n vs tau at one fixed (budget, picture, kappa, f_int, f_ret) cell."""

    budget_eV: float
    picture: str
    kappa: float
    points: tuple[tuple[float, float], ...]  # sorted (tau_ps, n_terminal_mean)


# ---------------------------------------------------------------------------
# cell keys
# ---------------------------------------------------------------------------
def _rk(value: float | None) -> float | None:
    """Round a float knob for stable dict-key grouping (None passes through)."""
    return None if value is None else round(float(value), 6)


def _kappa_cell_key(row: dict[str, Any]) -> tuple:
    return (
        _rk(row["budget_eV"]),
        row["picture"],
        _rk(row["f_int"]),
        _rk(row["f_ret"]),
        _rk(row["tau_ps"]),
    )


def _tau_cell_key(row: dict[str, Any]) -> tuple:
    return (
        _rk(row["budget_eV"]),
        row["picture"],
        _rk(row["kappa"]),
        _rk(row["f_int"]),
        _rk(row["f_ret"]),
    )


# ---------------------------------------------------------------------------
# assessment primitives
# ---------------------------------------------------------------------------
def assess_kappa_landscape(
    points: list[tuple[float, float]] | tuple[tuple[float, float], ...],
) -> tuple[str, float, float, str]:
    """Classify a W1-vs-kappa landscape.

    ``points`` is ``(kappa, W1)`` pairs. Returns
    ``(status, argmin_kappa, w1_min, notes)``. Order of checks: too few points
    (insufficient) -> flat (degenerate, not identifiable) -> argmin at a grid edge
    (not bracketed) -> interior minimum (identifiable). Flatness is checked before
    bracketing so a degenerate landscape is never reported as a real edge optimum.
    """
    pts = sorted(points, key=lambda kw: kw[0])
    if len(pts) < 2:
        # Guard before the argmin so an empty landscape degrades to the documented
        # insufficient status instead of raising on min() over an empty range.
        if not pts:
            return STATUS_INSUFFICIENT, float("nan"), float("nan"), "no kappa points"
        (only_kappa, only_w1) = pts[0]
        return STATUS_INSUFFICIENT, only_kappa, only_w1, "single kappa point"

    kappas = [k for k, _ in pts]
    w1s = [w for _, w in pts]
    argmin_i = min(range(len(pts)), key=lambda i: w1s[i])
    argmin_kappa = kappas[argmin_i]
    w1_min = w1s[argmin_i]

    w1_max = max(w1s)
    if w1_min > 0.0:
        rel_range = (w1_max - w1_min) / w1_min
        flat = rel_range < KAPPA_FLAT_FRACTION
    else:
        flat = w1_max == 0.0  # all-zero W1: degenerate
    if flat:
        return (
            STATUS_NOT_IDENTIFIABLE,
            argmin_kappa,
            w1_min,
            f"flat landscape (relative W1 range < {KAPPA_FLAT_FRACTION:.0%}): "
            "the observable does not resolve kappa",
        )

    if argmin_i in (0, len(pts) - 1):
        return (
            STATUS_NOT_BRACKETED,
            argmin_kappa,
            w1_min,
            f"minimum at grid edge kappa={argmin_kappa:g}: optimum not bracketed "
            "-> extend the grid / RRK-dof mechanism-level OQ",
        )

    return (
        STATUS_IDENTIFIABLE,
        argmin_kappa,
        w1_min,
        f"well-separated interior minimum at kappa={argmin_kappa:g}",
    )


def assess_picture_separation(
    best_w1_by_picture: dict[str, float],
) -> tuple[str, str | None, str]:
    """Classify how well the electronic pictures separate in best-over-kappa W1.

    Returns ``(status, best_picture, notes)``. Needs >= 2 pictures; a relative
    spread below :data:`PICTURE_SEPARATION_FRACTION` is reported not identifiable.
    """
    if len(best_w1_by_picture) < 2:
        only = next(iter(best_w1_by_picture), None)
        return STATUS_INSUFFICIENT, only, "fewer than two pictures scored"

    best_picture = min(best_w1_by_picture, key=best_w1_by_picture.get)
    w_best = best_w1_by_picture[best_picture]
    w_worst = max(best_w1_by_picture.values())
    if w_best > 0.0:
        rel = (w_worst - w_best) / w_best
    else:
        rel = float("inf") if w_worst > 0.0 else 0.0
    if rel < PICTURE_SEPARATION_FRACTION:
        return (
            STATUS_NOT_IDENTIFIABLE,
            best_picture,
            f"pictures separate by < {PICTURE_SEPARATION_FRACTION:.0%} in best W1",
        )
    return (
        STATUS_IDENTIFIABLE,
        best_picture,
        f"best picture {best_picture!r} separates by {rel:.0%} in W1",
    )


def assess_tau_sensitivity(
    points: list[tuple[float, float]] | tuple[tuple[float, float], ...],
) -> tuple[str, float, str]:
    """Classify terminal-n sensitivity to tau over a sweep.

    ``points`` is ``(tau_ps, n_terminal_mean)`` pairs. Returns
    ``(status, spread_he, notes)`` where ``spread_he`` is the terminal-n range.
    """
    if len(points) < 2:
        return STATUS_INSUFFICIENT, 0.0, "single tau point: sweep not yet run"
    n_means = [n for _, n in points]
    spread = max(n_means) - min(n_means)
    if spread > TAU_SENSITIVITY_THRESHOLD_HE:
        return (
            STATUS_SENSITIVE,
            spread,
            f"terminal n moves {spread:.2f} He across tau "
            f"(> {TAU_SENSITIVITY_THRESHOLD_HE} He): live calibration dimension",
        )
    return (
        STATUS_INSENSITIVE,
        spread,
        f"terminal n flat to {spread:.2f} He across tau "
        f"(<= {TAU_SENSITIVITY_THRESHOLD_HE} He): insensitivity confirmed",
    )


def assess_regime(
    rows: list[dict[str, Any]],
) -> tuple[dict[str, int], str | None, str]:
    """Aggregate the E5 regime labels and report the lowest-W1 regime.

    Returns ``(label_counts, best_regime, notes)``. ``best_regime`` is the regime
    label of the run with the smallest ``W1_matched`` (a factual report of which
    regime the observable favours, not a fidelity verdict).
    """
    counts = dict(Counter(r["regime_label"] for r in rows))
    scored = [r for r in rows if _is_number(r.get("W1_matched"))]
    best_regime = None
    if scored:
        best_regime = min(scored, key=lambda r: r["W1_matched"])["regime_label"]
    reachable = sum(1 for r in rows if r.get("total_strip_reachable"))
    notes = (
        f"regime counts {counts}; total_strip_reachable in {reachable}/{len(rows)} "
        "runs (E5 t*/Pi spine)"
    )
    return counts, best_regime, notes


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value == value


# ---------------------------------------------------------------------------
# grouping builders
# ---------------------------------------------------------------------------
def _accumulate_cell(
    cells: dict[tuple, dict[float, float]],
    order: list[tuple],
    key: tuple,
    sub_key: float,
    value: float,
    *,
    axis: str,
) -> None:
    """Insert ``value`` at ``cells[key][sub_key]``, warning loudly on a real drop.

    Two complete run dirs sharing the same knob cell *and* the same swept-axis
    value (a stale re-run left in place, or an undeclared replicate) would collapse
    last-write-wins. Rather than drop the earlier W1/terminal-n silently, warn when
    the incoming value differs from what is already there (identical values are a
    harmless deterministic replay). Discovery order is sorted, so which value is
    "last" is deterministic -- the concern is the silent loss, not nondeterminism.
    """
    bucket = cells.get(key)
    if bucket is None:
        bucket = cells[key] = {}
        order.append(key)
    if sub_key in bucket and bucket[sub_key] != value:
        warnings.warn(
            f"duplicate {axis} cell {key} at {axis}={sub_key:g}: keeping {value:g}, "
            f"dropping prior {bucket[sub_key]:g} (stale or replicate run dir?)",
            RuntimeWarning,
            stacklevel=3,
        )
    bucket[sub_key] = value


def build_kappa_landscapes(rows: list[dict[str, Any]]) -> list[KappaLandscape]:
    """Group scoreboard rows into W1-vs-kappa landscapes (cells with >= 2 kappa).

    A landscape holds one fixed (budget, picture, f_int, f_ret, tau) cell with the
    ladder steepness swept -- i.e. the Stage-1 co-fit cell. A cell with a single
    kappa is not a landscape and is dropped. Non-finite ``W1_matched`` is skipped
    (defensive: W1 is fail-loud upstream in E4, so this never fires in practice --
    the guard keeps this builder consistent with ``_best_w1_by_picture`` /
    ``assess_regime``).
    """
    cells: dict[tuple, dict[float, float]] = {}
    order: list[tuple] = []
    for row in rows:
        w1 = row["W1_matched"]
        if not _is_number(w1):
            continue
        _accumulate_cell(
            cells, order, _kappa_cell_key(row), float(row["kappa"]), float(w1),
            axis="kappa",
        )

    landscapes: list[KappaLandscape] = []
    for key in order:
        kappa_w1 = cells[key]
        if len(kappa_w1) < 2:
            continue
        budget, picture, _f_int, _f_ret, tau = key
        pts = tuple(sorted(kappa_w1.items()))
        landscapes.append(
            KappaLandscape(
                budget_eV=budget, picture=picture, tau_ps=tau, points=pts
            )
        )
    return landscapes


def build_tau_sweeps(rows: list[dict[str, Any]]) -> list[TauSweep]:
    """Group scoreboard rows into terminal-n-vs-tau sweeps (cells with >= 2 tau).

    Non-finite ``n_terminal_mean`` is skipped (defensive, as in
    :func:`build_kappa_landscapes`): a NaN would otherwise make the tau spread NaN
    and, since ``NaN > threshold`` is False, silently report tau ``insensitive``
    -- a false negative on a live calibration dimension.
    """
    cells: dict[tuple, dict[float, float]] = {}
    order: list[tuple] = []
    for row in rows:
        if row.get("tau_ps") is None:
            continue
        n_mean = row["n_terminal_mean"]
        if not _is_number(n_mean):
            continue
        _accumulate_cell(
            cells, order, _tau_cell_key(row), float(row["tau_ps"]), float(n_mean),
            axis="tau",
        )

    sweeps: list[TauSweep] = []
    for key in order:
        tau_n = cells[key]
        if len(tau_n) < 2:
            continue
        budget, picture, kappa, _f_int, _f_ret = key
        pts = tuple(sorted(tau_n.items()))
        sweeps.append(
            TauSweep(budget_eV=budget, picture=picture, kappa=kappa, points=pts)
        )
    return sweeps


def _best_w1_by_picture(landscapes: list[KappaLandscape]) -> dict[str, float]:
    """Lowest co-fit W1 per electronic picture, taken over the kappa landscapes.

    Restricting the picture comparison to landscape cells (>= 2 kappa) makes it a
    genuine kappa x picture co-fit: single-kappa probe runs and the engineered
    ``total_strip`` / Stage-2 tau-sweep cells (distinct f_int/tau -> single-kappa
    cells -> not landscapes) cannot contaminate a picture's best W1. This keeps the
    electronic-picture verdict consistent with the kappa verdict (same source) and
    guarantees the globally best picture always owns a landscape for the kappa
    assessment to run on (no silent fallback to a different picture).
    """
    best: dict[str, float] = {}
    for ls in landscapes:
        w = ls.w1_min
        if ls.picture not in best or w < best[ls.picture]:
            best[ls.picture] = w
    return best


# ---------------------------------------------------------------------------
# the checklist (one budget)
# ---------------------------------------------------------------------------
def _kappa_checklist_entry(
    budget: float,
    landscapes: list[KappaLandscape],
    best_picture: str | None,
) -> tuple[dict[str, Any], float | None]:
    """kappa verdict from the co-fit landscape under the globally best picture.

    Returns ``(entry, best_kappa)``. ``best_picture`` is derived from the same
    ``landscapes`` (:func:`_best_w1_by_picture`), so whenever it is set it owns a
    landscape here -- the assessment is never silently attributed to a different
    picture. ``best_kappa`` (the argmin, passed on to the tau entry so tau is read
    at the co-fit optimum cell) is ``None`` unless a real minimum was located.
    """
    if not landscapes:
        return (
            _entry(
                budget,
                "ladder_shape_kappa",
                STATUS_INSUFFICIENT,
                "-",
                "no kappa sweep on disk",
                "need >= 2 kappa at a fixed (picture, tau) cell",
            ),
            None,
        )
    chosen = [ls for ls in landscapes if ls.picture == best_picture]
    # `or landscapes` covers only best_picture is None (no picture chosen); a set
    # best_picture always has a landscape, so no cross-picture fallback occurs.
    landscape = min(chosen or landscapes, key=lambda ls: ls.w1_min)
    status, argmin_kappa, w1_min, notes = assess_kappa_landscape(landscape.points)
    grid = ", ".join(f"{k:g}" for k, _ in landscape.points)
    best_kappa = argmin_kappa if status != STATUS_INSUFFICIENT else None
    entry = _entry(
        budget,
        "ladder_shape_kappa",
        status,
        f"kappa*={argmin_kappa:g}",
        f"picture={landscape.picture}, W1_min={w1_min:.4g}, kappa grid [{grid}]",
        notes + " (fit conditional on fixed RRK dof s; see 's')",
    )
    return entry, best_kappa


def _picture_checklist_entry(
    budget: float, best_w1_by_picture: dict[str, float]
) -> tuple[dict[str, Any], str | None]:
    status, best_picture, notes = assess_picture_separation(best_w1_by_picture)
    value = best_picture if best_picture is not None else "-"
    evidence = "; ".join(
        f"{p}:W1={w:.4g}" for p, w in sorted(best_w1_by_picture.items())
    )
    return (
        _entry(budget, "electronic_picture", status, value, evidence, notes),
        best_picture,
    )


def _tau_checklist_entry(
    budget: float,
    rows: list[dict[str, Any]],
    best_picture: str | None,
    best_kappa: float | None,
) -> dict[str, Any]:
    sweeps = build_tau_sweeps(rows)
    if not sweeps:
        return _entry(
            budget,
            "tau",
            STATUS_INSUFFICIENT,
            "-",
            "no tau sweep on disk (single tau per cell)",
            "Stage-2 tau sweep not yet run",
        )
    # tau identifiability is a statement about the cell the campaign realizes, so
    # prefer the sweep at the co-fit optimum (best picture AND best kappa), then any
    # sweep under the best picture, then the first available -- never just insertion
    # order across arbitrary kappa cells.
    best_kr = _rk(best_kappa) if best_kappa is not None else None

    def _priority(s: TauSweep) -> int:
        pic_ok = best_picture is None or s.picture == best_picture
        kap_ok = best_kr is not None and s.kappa == best_kr
        if pic_ok and kap_ok:
            return 0
        if pic_ok:
            return 1
        return 2

    sweep = min(sweeps, key=lambda s: (_priority(s), s.kappa))
    status, spread, notes = assess_tau_sensitivity(sweep.points)
    band = ", ".join(f"{t:g}" for t, _ in sweep.points)
    off_optimum = "" if _priority(sweep) == 0 else (
        " (no sweep at the co-fit optimum cell; nearest reported)"
    )
    return _entry(
        budget,
        "tau",
        status,
        f"dn={spread:.2f} He",
        f"picture={sweep.picture}, kappa={sweep.kappa:g}, tau band [{band}] ps",
        notes + off_optimum,
    )


def _regime_checklist_entry(
    budget: float, rows: list[dict[str, Any]]
) -> dict[str, Any]:
    counts, best_regime, notes = assess_regime(rows)
    # Identifiable as a determination when both regimes actually appear (the axis
    # is exercised); otherwise the campaign realizes a single regime (reported).
    status = STATUS_IDENTIFIABLE if len(counts) >= 2 else STATUS_INSUFFICIENT
    value = best_regime if best_regime is not None else "-"
    return _entry(budget, "terminal_regime", status, value, notes,
                  "lowest-W1 regime reported; shell-retaining <-> total-strip axis")


def _static_entries(budget: float) -> list[dict[str, Any]]:
    """The deterministic (scope-/mechanism-fixed) checklist entries."""
    return [
        _entry(
            budget,
            "f_int",
            STATUS_NOT_IDENTIFIABLE,
            "pinned (floor)",
            "timing-only knob",
            "timing-only: moves t* not the cascade budget; terminal n insensitive "
            "to f_int by construction -> not separable by this observable",
        ),
        _entry(
            budget,
            "f_ret",
            STATUS_NOT_IDENTIFIABLE,
            "prior 0.1",
            "9 Aa-only campaign",
            "the only route to f_ret is the 9/18 Aa density contrast; no 18 Aa "
            "shell reference exists (2026-07-03 scope), so f_ret is held, not fit",
        ),
        _entry(
            budget,
            "lambda_attach",
            STATUS_NOT_IDENTIFIABLE,
            "held fixed",
            "9 Aa-only campaign",
            "pickup is structurally dead in-window at 9 Aa (Pi <= 0.005, Phase-D "
            "bridge finding #4), so the size distribution carries ~no lambda "
            "sensitivity here",
        ),
        _entry(
            budget,
            "nu",
            STATUS_HELD_FIXED,
            "held fixed",
            "RRK prefactor not swept",
            "not a campaign knob; would be retuned only under an explicit "
            "mechanism-level OQ, never silently",
        ),
        _entry(
            budget,
            "s",
            STATUS_HELD_FIXED,
            "held fixed (s=59)",
            "RRK dof not swept",
            "s and ladder shape kappa are entangled (s<->kappa coupling): the "
            "kappa fit is conditional on this fixed RRK dof convention",
        ),
    ]


def _entry(
    budget: float,
    quantity: str,
    status: str,
    value_or_range: Any,
    evidence: str,
    notes: str,
) -> dict[str, Any]:
    return {
        "budget_eV": budget,
        "quantity": quantity,
        "status": status,
        "value_or_range": value_or_range,
        "evidence": evidence,
        "notes": notes,
    }


def build_identifiability_rows(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Assemble the per-quantity identifiability checklist, sectioned by budget.

    ``rows`` are F3 scoreboard rows (any mix of budgets). Returns one checklist row
    per (budget, load-bearing quantity), in :data:`LOAD_BEARING_QUANTITIES` order
    within each budget. Computed quantities (kappa/picture/tau/regime) are assessed
    from the rows; the rest are deterministic scope/mechanism statements.
    """
    by_budget: dict[float, list[dict[str, Any]]] = {}
    order: list[float] = []
    for row in rows:
        b = round(float(row["budget_eV"]), 6)
        if b not in by_budget:
            by_budget[b] = []
            order.append(b)
        by_budget[b].append(row)

    out: list[dict[str, Any]] = []
    for budget in sorted(order):
        budget_rows = by_budget[budget]
        # Build the co-fit landscapes once; kappa and picture share this source so
        # their verdicts stay consistent, and the best picture always owns a
        # landscape for the kappa assessment.
        landscapes = build_kappa_landscapes(budget_rows)
        picture_entry, best_picture = _picture_checklist_entry(
            budget, _best_w1_by_picture(landscapes)
        )
        kappa_entry, best_kappa = _kappa_checklist_entry(
            budget, landscapes, best_picture
        )
        entries = {
            "ladder_shape_kappa": kappa_entry,
            "electronic_picture": picture_entry,
            "tau": _tau_checklist_entry(
                budget, budget_rows, best_picture, best_kappa
            ),
            "terminal_regime": _regime_checklist_entry(budget, budget_rows),
        }
        for static in _static_entries(budget):
            entries[static["quantity"]] = static
        out.extend(entries[q] for q in LOAD_BEARING_QUANTITIES)
    return out


# ---------------------------------------------------------------------------
# text report
# ---------------------------------------------------------------------------
def _format_cell(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _headline(checklist: list[dict[str, Any]]) -> str:
    total = len(checklist)
    separated = sum(1 for r in checklist if r["status"] in _SEPARATED_STATUSES)
    return (
        f"Identifiability headline: {separated} of {total} load-bearing "
        "quantities are separated by the size distribution alone "
        "(identifiable/sensitive); the rest are held, timing-only, out of the "
        "9 Aa scope, or not yet swept."
    )


def format_report(
    rows: list[dict[str, Any]],
    *,
    checklist: list[dict[str, Any]] | None = None,
) -> str:
    """Render the full identifiability report as headless text, per budget.

    Pass a pre-built ``checklist`` (from :func:`build_identifiability_rows`) to
    avoid re-running the whole assessment; otherwise it is built from ``rows``.
    """
    if checklist is None:
        checklist = build_identifiability_rows(rows)
    if not checklist:
        return "(no Tier-2 rows: nothing to assess)"

    budgets = sorted({r["budget_eV"] for r in checklist})
    lines: list[str] = [
        "Tier-2 identifiability + regime-determination report.",
        "Reported, not auto-adjudicated: minima/regime are factual summaries; "
        "flatness/bracketing/sensitivity thresholds are reader aids, not gates.",
        "",
    ]
    for budget in budgets:
        section = [r for r in checklist if r["budget_eV"] == budget]
        lines.append(f"=== budget {budget:g} eV ===")
        lines.append(_headline(section))
        lines.append("")
        lines.append(_format_checklist_table(section))
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _format_checklist_table(section: list[dict[str, Any]]) -> str:
    cols = ["quantity", "status", "value_or_range", "evidence", "notes"]
    widths = {
        c: max(len(c), *(len(_format_cell(r[c])) for r in section)) for c in cols
    }
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    sep = "  ".join("-" * widths[c] for c in cols)
    body = [
        "  ".join(_format_cell(r[c]).ljust(widths[c]) for c in cols)
        for r in section
    ]
    return "\n".join([header, sep, *body])


def write_report_csv(path: str | Path, checklist_rows: list[dict[str, Any]]) -> Path:
    """Write identifiability checklist rows to CSV in the stable column order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=IDENTIFIABILITY_COLUMNS)
        writer.writeheader()
        writer.writerows(checklist_rows)
    return p


# ---------------------------------------------------------------------------
# collector wrapper + operator entry
# ---------------------------------------------------------------------------
def collect_identifiability_rows(
    *,
    runs_root: str | Path = RUNS_ROOT,
    reference_path: str | Path = REFERENCE_PATH,
    budget_eV: float | None = BUDGET_EV,
    n_max: int = N_STAR,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Score the campaign (via F3) and assemble the identifiability checklist.

    Returns ``(scoreboard_rows, checklist_rows)``: the F3 rows the assessment ran
    on and the per-(budget, quantity) checklist. Imports the F3 *typed* records,
    never the lossy CSV.
    """
    records, _ = collect_tier2_records(
        runs_root=runs_root,
        reference_path=reference_path,
        budget_eV=budget_eV,
        n_max=n_max,
    )
    scoreboard_rows = [record.row for record in records]
    checklist_rows = build_identifiability_rows(scoreboard_rows)
    return scoreboard_rows, checklist_rows


def main() -> int:
    """Print the Tier-2 identifiability report and optionally save the CSV."""
    scoreboard_rows, checklist_rows = collect_identifiability_rows(
        runs_root=RUNS_ROOT,
        reference_path=REFERENCE_PATH,
        budget_eV=BUDGET_EV,
    )
    print(f"Assessed {len(scoreboard_rows)} campaign run(s) under {RUNS_ROOT}.")
    print(format_report(scoreboard_rows, checklist=checklist_rows))
    if SAVE_CSV_PATH is not None and checklist_rows:
        path = write_report_csv(SAVE_CSV_PATH, checklist_rows)
        print(f"Wrote Tier-2 identifiability report -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
