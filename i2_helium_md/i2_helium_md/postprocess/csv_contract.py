"""Shared CSV column-contract validation for the postprocess reference loaders.

Both reference loaders (:mod:`~i2_helium_md.postprocess.hedft_loader`,
:mod:`~i2_helium_md.postprocess.abundance_loader`) enforce the same
**order-independent membership** header match: missing *and* extra columns
raise. Factored here so the contract lives once (Quality Principle 1); a
future Phase-F loader reuses this instead of copying the block a third time.
"""

from __future__ import annotations

from collections.abc import Sequence


def validate_columns(
    actual_columns: Sequence[str] | None,
    expected_columns: Sequence[str],
    *,
    file_label: str,
) -> None:
    """Raise unless ``actual_columns`` matches ``expected_columns`` by membership.

    Order-independent: columns may appear in any order, but every expected
    column must be present and no unexpected column may appear.

    Parameters
    ----------
    actual_columns
        Column names found in the file (e.g. ``structured.dtype.names``;
        ``None`` is treated as no columns).
    expected_columns
        The loader's expected column names.
    file_label
        Human-readable file identifier for the error message, e.g.
        ``f"HeDFT reference {p.name}"``.

    Raises
    ------
    ValueError
        Listing the missing and unexpected columns.
    """
    actual = tuple(actual_columns or ())
    expected = tuple(expected_columns)
    missing = tuple(c for c in expected if c not in actual)
    extra = tuple(c for c in actual if c not in expected)
    if missing or extra:
        raise ValueError(
            f"{file_label} has unexpected columns. "
            f"missing={list(missing)}, unexpected={list(extra)}, "
            f"expected={list(expected)}"
        )
