"""Names of the thesis run directories under ``data/runs``."""

from __future__ import annotations

from pathlib import Path

from ..presets import THESIS_MEMBER_SEEDS, THESIS_PRESETS

#: Run directories (git-ignored) of the source tree (editable install).
RUNS_ROOT = Path(__file__).resolve().parents[2] / "data" / "runs"
#: Committed tables of the thesis numbers.
TABLES_ROOT = RUNS_ROOT / "h2b_forward_model"
#: Molecules per thesis member.
PRODUCTION_N = 1000


def _check_law(law: str) -> None:
    if law not in THESIS_PRESETS:
        raise KeyError(f"unknown law {law!r}; expected one of {list(THESIS_PRESETS)}")


def member_dir(law: str, member: str, *, root: Path = RUNS_ROOT, tag: str = "") -> Path:
    """Run directory of one battery member.

    Parameters
    ----------
    law : str
        A key of ``THESIS_PRESETS``.
    member : str
        A key of ``THESIS_MEMBER_SEEDS``.
    root : Path
        Parent directory.
    tag : str
        Suffix of the law name, for variants of a thesis law.

    Raises
    ------
    KeyError
        Unknown law or member; the message lists the choices.
    """
    _check_law(law)
    if member not in THESIS_MEMBER_SEEDS:
        raise KeyError(f"unknown member {member!r}; expected one of {list(THESIS_MEMBER_SEEDS)}")
    return root / f"{law}{tag}_N{PRODUCTION_N}_{member}"


def pool_dir(law: str, *, root: Path = RUNS_ROOT, tag: str = "") -> Path:
    """Directory of the pool of all members; arguments as for :func:`member_dir`."""
    _check_law(law)
    return root / f"{law}{tag}_N{PRODUCTION_N * len(THESIS_MEMBER_SEEDS)}_pooled"


def require_run(path: Path, driver: str) -> Path:
    """Return ``path``, or name the script that makes it.

    Parameters
    ----------
    path : Path
        Run directory, or a file in one.
    driver : str
        File name of the script under ``scripts/`` that writes it.

    Raises
    ------
    FileNotFoundError
        ``path`` does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; make it with scripts/{driver}")
    return path
