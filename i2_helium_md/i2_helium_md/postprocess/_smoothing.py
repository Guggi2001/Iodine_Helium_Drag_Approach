"""Internal smoothing helpers shared by post-processing scripts.

Both ``plot_experimental_comparison.py`` and ``plot_paper_v3.py``
need MATLAB-style ``movmean`` (centred moving mean with shortened
endpoint windows) plus a baseline-subtract / unit-max trace
normaliser. Keeping them here avoids the script-to-script imports
that the legacy MATLAB code resorts to and gives unit tests a single
place to pin behaviour.
"""

from __future__ import annotations

import numpy as np


def moving_mean(values, window: int) -> np.ndarray:
    """Centred moving mean with shortened endpoint windows (MATLAB ``movmean``).

    Parameters
    ----------
    values : array-like
        1-D input.
    window : int
        Window size in samples. ``window <= 1`` returns a copy.

    Returns
    -------
    np.ndarray
        Same length as input.
    """
    data = np.asarray(values, dtype=float)
    if window <= 1 or data.size == 0:
        return data.copy()

    half_left = (window - 1) // 2
    half_right = window // 2
    result = np.empty_like(data, dtype=float)
    for idx in range(data.size):
        start = max(0, idx - half_left)
        stop = min(data.size, idx + half_right + 1)
        result[idx] = data[start:stop].mean()
    return result


def normalise_trace(values) -> np.ndarray:
    """Baseline-subtract and scale a trace to unit maximum for plotting.

    Returns zeros if the input is constant.
    """
    data = np.asarray(values, dtype=float)
    shifted = data - data.min()
    scale = shifted.max()
    if scale <= 0.0:
        return np.zeros_like(shifted)
    return shifted / scale
