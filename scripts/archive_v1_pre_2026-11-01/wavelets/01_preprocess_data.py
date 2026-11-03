#!/usr/bin/env python3
"""Preprocess one cell's trace for continuous wavelet analysis."""
from __future__ import annotations

import os
import re
from typing import Any

import numpy as np
import pandas as pd

MIN_RESAMPLED_POINTS = 100
TIME_COLUMNS = ("Time", "time", "Timestamp", "timestamp", "t", "Time_s", "Time (s)")


def _time_seconds(values: pd.Series, column: str) -> np.ndarray:
    """Convert supported numeric, Excel, duration, and datetime values to seconds.

    Numeric values below one are treated as Excel fractions of a day, following
    the MATLAB preprocessing helper. Other numeric values are treated as
    seconds. Text clock values support both ``mm:ss[.fraction]`` and
    ``hh:mm:ss[.fraction]`` forms, with optional AM/PM on the latter.
    """
    if pd.api.types.is_datetime64_any_dtype(values):
        seconds = (values - values.dropna().iloc[0]).dt.total_seconds().to_numpy(float)
        return seconds
    if pd.api.types.is_timedelta64_dtype(values):
        return (values - values.dropna().iloc[0]).dt.total_seconds().to_numpy(float)

    def relative_seconds(seconds: np.ndarray) -> np.ndarray:
        seconds = np.asarray(seconds, dtype=float)
        finite = np.flatnonzero(np.isfinite(seconds))
        if finite.size:
            seconds = seconds - seconds[finite[0]]
        return seconds

    text = values.astype("string").str.strip().replace("", pd.NA)
    # CSV readers usually infer numeric columns, but also handle numeric text.
    numeric = pd.to_numeric(text, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
    nonempty = text.notna().to_numpy()
    if np.isfinite(numeric[nonempty]).all() and nonempty.any():
        largest = float(np.max(numeric[nonempty]))
        if largest < 1:
            numeric *= 86400.0  # Excel fraction of one day
        else:
            unit = column.lower()
            if "ms" in unit:
                numeric /= 1000.0
            elif "min" in unit:
                numeric *= 60.0
            else:
                declared_unit = os.environ.get("CALCIUM_TIME_UNIT", "s").lower()
                if declared_unit in {"ms", "millisecond", "milliseconds"}:
                    numeric /= 1000.0
                elif declared_unit in {"min", "minute", "minutes"}:
                    numeric *= 60.0
                elif declared_unit not in {"s", "sec", "seconds"}:
                    raise ValueError(f"Unsupported CALCIUM_TIME_UNIT: {declared_unit}")
        return relative_seconds(numeric)

    result = np.full(len(text), np.nan, dtype=float)
    clock_pattern = re.compile(
        r"^(?P<a>\d+):(?P<b>\d{1,2}(?:\.\d+)?)(?::(?P<c>\d{1,2}(?:\.\d+)?))?"
        r"(?:\s*(?P<ampm>AM|PM))?$", re.IGNORECASE
    )
    clock_ok = True
    for index, item in enumerate(text):
        if pd.isna(item) or not item:
            continue
        match = clock_pattern.fullmatch(str(item))
        if match is None:
            clock_ok = False
            break
        a, b, c, ampm = match.group("a", "b", "c", "ampm")
        if c is None:  # mm:ss[.fraction], minutes may exceed 59
            second = float(b)
            if second >= 60:
                clock_ok = False
                break
            result[index] = int(a) * 60.0 + second
        else:  # hh:mm:ss[.fraction]
            if "." in b:
                clock_ok = False
                break
            hour, minute, second = int(a), int(b), float(c)
            if minute >= 60 or second >= 60:
                clock_ok = False
                break
            if ampm:
                if hour < 1 or hour > 12:
                    clock_ok = False
                    break
                hour = hour % 12 + (12 if ampm.upper() == "PM" else 0)
            result[index] = hour * 3600.0 + minute * 60.0 + second
    if clock_ok and np.isfinite(result).any():
        return relative_seconds(result)

    # Duration parser handles day-qualified durations and other pandas forms.
    duration = pd.to_timedelta(text, errors="coerce")
    if duration.notna().any() and duration.notna().sum() == text.notna().sum():
        return relative_seconds(duration.dt.total_seconds().to_numpy(float))

    # Finally accept calendar date/time strings, including fractional seconds
    # and AM/PM forms. Subtracting the first valid timestamp makes elapsed time.
    try:
        parsed = pd.to_datetime(text, errors="coerce", format="mixed")
    except (TypeError, ValueError):  # pandas versions before format="mixed"
        parsed = pd.to_datetime(text, errors="coerce")
    if parsed.notna().sum() < text.notna().sum():
        parsed = pd.to_datetime(text, errors="coerce")
    if parsed.notna().any() and parsed.notna().sum() == text.notna().sum():
        seconds = parsed.astype("int64").to_numpy(float) / 1e9
        seconds[parsed.isna().to_numpy()] = np.nan
        return relative_seconds(seconds)

    examples = ", ".join(str(v) for v in text.dropna().head(5))
    raise ValueError(f"Cannot parse time column {column!r}; example values: {examples}")


def preprocess_data(rows: pd.DataFrame) -> dict[str, Any]:
    """Return uniformly resampled ΔF/F, sampling rate, coordinates, and gap data.

    A gap is any interval greater than twice the median interval. Segments are
    linearly resampled independently, retaining endpoints. Segments shorter
    than 100 resampled observations are discarded.
    """
    if "ID" in rows.columns and rows["ID"].nunique(dropna=True) > 1:
        raise ValueError("preprocess_data received rows for multiple cell IDs")
    time_col = next((name for name in TIME_COLUMNS if name in rows.columns), None)
    if time_col is None:
        raise ValueError(f"No time column found; expected one of {TIME_COLUMNS}")
    intensity_col = next((c for c in rows.columns if c.lower() == "meanintensity"), None)
    if intensity_col is None:
        intensity_col = next((c for c in rows.columns
                              if "mean" in c.lower() and "intens" in c.lower()), None)
    if intensity_col is None:
        raise ValueError("No MeanIntensity-like column found")
    missing = [c for c in ("PositionX", "PositionY") if c not in rows.columns]
    if missing:
        raise ValueError(f"Missing spatial columns: {missing}")

    time = _time_seconds(rows[time_col], time_col)
    trace = pd.to_numeric(rows[intensity_col], errors="coerce").to_numpy(float)
    x = pd.to_numeric(rows["PositionX"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(rows["PositionY"], errors="coerce").to_numpy(float)
    valid = np.isfinite(time) & np.isfinite(trace) & np.isfinite(x) & np.isfinite(y)
    time, trace, x, y = time[valid], trace[valid], x[valid], y[valid]
    if time.size < 2:
        raise ValueError("Fewer than two complete observations")
    if np.any(np.diff(time) < 0):
        raise ValueError("Time values are not monotonic; trace was not reordered")

    # Match MATLAB unique(..., 'stable'): retain the first observation at a time.
    keep = np.r_[True, np.diff(time) != 0]
    time, trace, x, y = time[keep], trace[keep], x[keep], y[keep]
    dt = np.diff(time)
    dt_med = float(np.median(dt))
    if not np.isfinite(dt_med) or dt_med <= 0:
        raise ValueError("Invalid median sampling interval")
    gap_indices = np.flatnonzero(dt > 2 * dt_med)
    starts = np.r_[0, gap_indices + 1]
    ends = np.r_[gap_indices + 1, time.size]

    time_parts, trace_parts, x_parts, y_parts = [], [], [], []
    for start, end in zip(starts, ends):
        if end - start < 2:
            continue
        t_seg, f_seg, x_seg, y_seg = time[start:end], trace[start:end], x[start:end], y[start:end]
        expected_n = max(2, int(round((t_seg[-1] - t_seg[0]) / dt_med)) + 1)
        if expected_n < MIN_RESAMPLED_POINTS:
            continue
        t_uniform = np.linspace(t_seg[0], t_seg[-1], expected_n)
        time_parts.append(t_uniform)
        trace_parts.append(np.interp(t_uniform, t_seg, f_seg))
        x_parts.append(np.interp(t_uniform, t_seg, x_seg))
        y_parts.append(np.interp(t_uniform, t_seg, y_seg))
    if not time_parts:
        raise ValueError(f"No gap-free segment had at least {MIN_RESAMPLED_POINTS} resampled points")

    t_uniform = np.concatenate(time_parts)
    fluorescence = np.concatenate(trace_parts)
    baseline = max(float(np.percentile(fluorescence, 10)), np.finfo(float).eps)
    return {
        "time": t_uniform,
        "dFF": (fluorescence - baseline) / baseline,
        "Fs": 1.0 / dt_med,
        "X": np.concatenate(x_parts),
        "Y": np.concatenate(y_parts),
        "gaps": gap_indices,
        "dt": np.diff(t_uniform),
        "dt_med": dt_med,
        "skip": False,
    }
