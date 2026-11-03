#!/usr/bin/env python3
"""Preprocess one cell's trace for continuous wavelet analysis (v2).

Matches the PDT paper's preprocessing wherever the paper defines it (Methods
2.1, "Image segmentation and cell tracking" and "Identification of active
cells"). The v1 scripts are archived in scripts/archive_v1_pre_2026-11-01.

* Gaps (paper rules):
  - Tracking dropouts of up to 3 missing frames (TrackBinaries maximum
    3-frame gap; frame-index step <= CALCIUM_MAX_MISSING_FRAMES+1) are bridged
    by linear interpolation.
  - Refocus pauses split the trace into segments that are processed
    separately. The data point just before and just after each pause is
    ignored (CALCIUM_TRIM_AT_SPLIT, default 1). A pause is a frame-to-frame
    time jump above 1.5x the embryo's frame interval while the frame counter
    (ND.T) advances by one.
  - Larger frame-index jumps (> 3 missing frames) also split.
  - Without ND.T, the fallback is a time rule: split at intervals
    > CALCIUM_GAP_SPLIT_SECONDS (default 10 s).
* Inclusion: a cell needs >= 100 observed frames in usable segments (the
  paper's 100-frame tracking criterion; CALCIUM_MIN_CELL_FRAMES). A segment
  is usable only if it is longer than twice the wavelet cone of influence at
  the lowest frequency (~122 s at 0.02 Hz), because shorter segments have no
  valid low-frequency samples.
* Baseline: F0 = average of the lower ten percent of the segment's
  intensities (paper baseline definition (a)), computed per segment because
  the paper processes the segments around refocus events separately.
  dF/F = (F - F0) / F0.
* Detrending: each segment's dF/F has its least-squares line subtracted
  (CALCIUM_DETREND=linear, the default; "none" disables it) to remove
  photobleaching drift before the wavelet transform.
* Time base: exported timestamps are rounded to whole seconds. Each frame
  gets a de-jittered time from an embryo-wide linear fit of time on ND.T
  within pause-free runs. Columns mixing hh:mm:ss, mm:ss.fff and plain
  seconds are parsed value by value.
"""
from __future__ import annotations

import os
import re
from typing import Any

import numpy as np
import pandas as pd

TIME_COLUMNS = ("Time", "time", "Timestamp", "timestamp", "t", "Time_s", "Time (s)")
GAP_SPLIT_SECONDS = float(os.environ.get("CALCIUM_GAP_SPLIT_SECONDS", "10"))  # fallback only
MAX_MISSING_FRAMES = int(os.environ.get("CALCIUM_MAX_MISSING_FRAMES", "3"))
BASELINE_FRACTION = float(os.environ.get("CALCIUM_BASELINE_FRACTION", "0.10"))
TRIM_AT_SPLIT = int(os.environ.get("CALCIUM_TRIM_AT_SPLIT", "1"))
MIN_CELL_FRAMES = int(os.environ.get("CALCIUM_MIN_CELL_FRAMES", "100"))
DETREND = os.environ.get("CALCIUM_DETREND", "linear").lower()  # linear | none
if DETREND not in {"linear", "none"}:
    raise ValueError(f"CALCIUM_DETREND must be 'linear' or 'none', not {DETREND!r}")

_CLOCK = re.compile(
    r"^(?P<a>\d+):(?P<b>\d{1,2}(?:\.\d+)?)(?::(?P<c>\d{1,2}(?:\.\d+)?))?"
    r"(?:\s*(?P<ampm>AM|PM))?$", re.IGNORECASE)


def _clock_seconds(item: str) -> float | None:
    """Parse mm:ss[.f] or hh:mm:ss[.f] [AM|PM]; return None if not a clock value.

    Also accepts the export artefact ``mm:ss.-ff`` (e.g. ``00:00.-01``), which
    some first frames use for a small negative offset. In plain-seconds exports
    the same frame appears as ``-0.001``.
    """
    negative = re.fullmatch(r"(\d+):(\d{1,2})\.-(\d+)", item)
    if negative:
        return int(negative[1]) * 60.0 + int(negative[2]) - float("0." + negative[3])
    match = _CLOCK.fullmatch(item)
    if match is None:
        return None
    a, b, c, ampm = match.group("a", "b", "c", "ampm")
    if c is None:  # mm:ss[.fraction]; minutes may exceed 59
        second = float(b)
        return int(a) * 60.0 + second if second < 60 else None
    if "." in b:
        return None
    hour, minute, second = int(a), int(b), float(c)
    if minute >= 60 or second >= 60:
        return None
    if ampm:
        if not 1 <= hour <= 12:
            return None
        hour = hour % 12 + (12 if ampm.upper() == "PM" else 0)
    return hour * 3600.0 + minute * 60.0 + second


def _relative(seconds: np.ndarray) -> np.ndarray:
    seconds = np.asarray(seconds, dtype=float)
    finite = np.flatnonzero(np.isfinite(seconds))
    return seconds - seconds[finite[0]] if finite.size else seconds


def _time_seconds(values: pd.Series, column: str) -> np.ndarray:
    """Convert numeric, clock-text, mixed, duration or datetime values to elapsed seconds."""
    if pd.api.types.is_datetime64_any_dtype(values) or pd.api.types.is_timedelta64_dtype(values):
        return (values - values.dropna().iloc[0]).dt.total_seconds().to_numpy(float)

    text = values.astype("string").str.strip().replace("", pd.NA)
    nonempty = text.notna().to_numpy()
    numeric = pd.to_numeric(text, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
    if nonempty.any() and np.isfinite(numeric[nonempty]).all():
        if float(np.max(numeric[nonempty])) < 1:
            numeric *= 86400.0  # Excel fraction of one day
        else:
            unit = column.lower()
            declared = os.environ.get("CALCIUM_TIME_UNIT", "s").lower()
            if "ms" in unit or declared in {"ms", "millisecond", "milliseconds"}:
                numeric /= 1000.0
            elif "min" in unit or declared in {"min", "minute", "minutes"}:
                numeric *= 60.0
            elif declared not in {"s", "sec", "seconds"}:
                raise ValueError(f"Unsupported CALCIUM_TIME_UNIT: {declared}")
        return _relative(numeric)

    # Clock text, or a per-value mixture of clock text and plain seconds. One
    # previous-student export (190724SP_E2) mixes hh:mm:ss, mm:ss.fff and seconds
    # within single tracks; all three share the same elapsed-time origin.
    result = np.full(len(text), np.nan, dtype=float)
    parsed_all = True
    for index, item in enumerate(text):
        if pd.isna(item):
            continue
        value = _clock_seconds(str(item))
        if value is None and np.isfinite(numeric[index]):
            value = float(numeric[index])  # plain seconds
        if value is None:
            parsed_all = False
            break
        result[index] = value
    if parsed_all and np.isfinite(result).any():
        return _relative(result)

    duration = pd.to_timedelta(text, errors="coerce")
    if duration.notna().any() and duration.notna().sum() == text.notna().sum():
        return _relative(duration.dt.total_seconds().to_numpy(float))
    try:
        parsed = pd.to_datetime(text, errors="coerce", format="mixed", utc=True)
    except (TypeError, ValueError):
        parsed = pd.to_datetime(text, errors="coerce", utc=True)
    if parsed.notna().any() and parsed.notna().sum() == text.notna().sum():
        seconds = parsed.astype("int64").to_numpy(float) / 1e9
        seconds[parsed.isna().to_numpy()] = np.nan
        return _relative(seconds)
    examples = ", ".join(str(v) for v in text.dropna().head(5))
    raise ValueError(f"Cannot parse time column {column!r}; example values: {examples}")


def embryo_frame_times(table: pd.DataFrame) -> dict[str, np.ndarray] | None:
    """Embryo-wide de-jittered row times and refocus-pause locations.

    All cells in one frame (``ND.T``) share an acquisition time. Times are
    regressed linearly on the frame index within each pause-free run. A pause
    is a step whose time per frame exceeds 1.5x the typical frame interval.
    Returns ``{"time": seconds per row, "pause_frames": first frame after each
    pause, "frame_dt": typical interval}``, or None if ``ND.T`` is absent.
    """
    time_col = next((name for name in TIME_COLUMNS if name in table.columns), None)
    if time_col is None or "ND.T" not in table.columns:
        return None
    seconds = _time_seconds(table[time_col], time_col)
    frame = pd.to_numeric(table["ND.T"], errors="coerce").to_numpy(float)
    ok = np.isfinite(seconds) & np.isfinite(frame)
    per_frame = (pd.DataFrame({"f": frame[ok], "t": seconds[ok]})
                 .groupby("f")["t"].median().sort_index())
    f, t = per_frame.index.to_numpy(float), per_frame.to_numpy(float)
    if f.size < 3:
        return None
    step = np.diff(t) / np.diff(f)
    positive = step[step > 0]
    if positive.size == 0:
        return None
    typical = float(np.median(positive))
    breaks = np.flatnonzero((step > 1.5 * typical) | (step <= 0)) + 1
    fitted, slopes = t.copy(), []
    for start, end in zip(np.r_[0, breaks], np.r_[breaks, f.size]):
        if end - start >= 3:
            slope, intercept = np.polyfit(f[start:end], t[start:end], 1)
            fitted[start:end] = intercept + slope * f[start:end]
            slopes.append((end - start, slope))
    frame_dt = (sum(n * b for n, b in slopes) / sum(n for n, _ in slopes)) if slopes else typical
    lookup = dict(zip(f, fitted))
    out = np.full(len(table), np.nan)
    out[ok] = [lookup[v] for v in frame[ok]]
    return {"time": out, "pause_frames": f[breaks], "frame_dt": float(frame_dt)}


def min_segment_seconds(lowest_frequency: float, coi_factor: float) -> float:
    """Shortest segment with at least one cone-of-influence-valid sample at f_min."""
    return 2.0 * coi_factor / lowest_frequency


def lower_fraction_mean(values: np.ndarray, fraction: float = BASELINE_FRACTION) -> float:
    """Average of the lowest ``fraction`` of values (paper baseline definition (a))."""
    ordered = np.sort(np.asarray(values, dtype=float))
    n = max(1, int(np.ceil(fraction * ordered.size)))
    return float(np.mean(ordered[:n]))


def preprocess_data(rows: pd.DataFrame, lowest_frequency: float, coi_factor: float,
                    frame_info: dict | None = None) -> dict[str, Any]:
    """Return per-segment uniformly resampled dF/F traces and QC metadata.

    ``coi_factor`` is the cone-of-influence e-folding time x frequency (s*Hz)
    from the spectral stage. ``frame_info`` comes from
    :func:`embryo_frame_times`, with ``time`` already subset to ``rows``.
    Without it, the time column is parsed and the time-based fallback gap rule
    is used.
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

    use_frames = frame_info is not None and "ND.T" in rows.columns
    if use_frames:
        time = np.asarray(frame_info["time"], dtype=float)
        frame = pd.to_numeric(rows["ND.T"], errors="coerce").to_numpy(float)
    else:
        time = _time_seconds(rows[time_col], time_col)
        frame = np.full(time.size, np.nan)
    trace = pd.to_numeric(rows[intensity_col], errors="coerce").to_numpy(float)
    x = pd.to_numeric(rows["PositionX"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(rows["PositionY"], errors="coerce").to_numpy(float)
    valid = np.isfinite(time) & np.isfinite(trace) & np.isfinite(x) & np.isfinite(y)
    if use_frames:
        valid &= np.isfinite(frame)
    time, frame, trace, x, y = time[valid], frame[valid], trace[valid], x[valid], y[valid]
    if time.size < 2:
        raise ValueError("Fewer than two complete observations")
    # Merged multi-layer tracks can hold gap-filling rows out of time order.
    # Each row carries its own frame time, so sort (stable) rather than reject.
    order = np.argsort(time, kind="stable")
    time, frame, trace, x, y = time[order], frame[order], trace[order], x[order], y[order]
    keep = np.r_[True, np.diff(time) != 0]  # first observation at duplicated times
    time, frame, trace, x, y = time[keep], frame[keep], trace[keep], x[keep], y[keep]
    span_seconds = float(time[-1] - time[0])
    dt = np.diff(time)

    if use_frames:
        frame_dt = float(frame_info["frame_dt"])
        fstep = np.diff(frame)
        pauses = np.asarray(frame_info["pause_frames"], dtype=float)
        # A pause lies between consecutive observations if its first post-pause
        # frame is in (frame[i], frame[i+1]].
        pauses = np.sort(pauses)
        crosses_pause = (np.searchsorted(pauses, frame[1:], side="right")
                         > np.searchsorted(pauses, frame[:-1], side="right"))
        split_mask = crosses_pause | (fstep > MAX_MISSING_FRAMES + 1)
        bridged_mask = ~split_mask & (fstep > 1)
    else:
        bridged = dt[dt <= GAP_SPLIT_SECONDS]
        if bridged.size == 0:
            raise ValueError("No inter-frame interval within the gap tolerance")
        single = bridged[bridged < 1.5 * np.median(bridged)]
        frame_dt = float(np.mean(single if single.size else bridged))
        split_mask = dt > GAP_SPLIT_SECONDS
        bridged_mask = ~split_mask & (dt > 1.5 * frame_dt)
    splits = np.flatnonzero(split_mask)
    starts = np.r_[0, splits + 1]
    ends = np.r_[splits + 1, time.size]

    min_seconds = min_segment_seconds(lowest_frequency, coi_factor)
    segments, frames_used = [], 0
    for index, (start, end) in enumerate(zip(starts, ends)):
        if index > 0:
            start += TRIM_AT_SPLIT  # point just after a refocus pause
        if index < len(starts) - 1:
            end -= TRIM_AT_SPLIT  # point just before a refocus pause
        if end - start < 2:
            continue
        t_seg, f_seg = time[start:end], trace[start:end]
        duration = float(t_seg[-1] - t_seg[0])
        if duration <= min_seconds:
            continue
        n = max(2, int(round(duration / frame_dt)) + 1)
        t_uniform = np.linspace(t_seg[0], t_seg[-1], n)
        f_uniform = np.interp(t_uniform, t_seg, f_seg)
        baseline = max(lower_fraction_mean(f_seg), np.finfo(float).eps)
        dff = (f_uniform - baseline) / baseline
        if DETREND == "linear":
            # Remove slow photobleaching drift; inside the cone of influence this
            # leaves the 0.005-0.08 Hz band essentially unchanged.
            dff = dff - np.polyval(np.polyfit(t_uniform - t_uniform[0], dff, 1),
                                   t_uniform - t_uniform[0])
        segments.append({
            "time": t_uniform,
            "dFF": dff,
            "dt": duration / (n - 1),
            "X": np.interp(t_uniform, t_seg, x[start:end]),
            "Y": np.interp(t_uniform, t_seg, y[start:end]),
            "baseline": baseline,
        })
        frames_used += end - start
    if not segments:
        raise ValueError(f"No segment longer than {min_seconds:.0f} s after splitting")
    if frames_used < MIN_CELL_FRAMES:
        raise ValueError(f"Only {frames_used} observed frames in usable segments "
                         f"(< {MIN_CELL_FRAMES})")
    return {
        "segments": segments,
        "frame_dt": frame_dt,
        "span_seconds": span_seconds,
        "frames_observed": int(time.size),
        "frames_used": int(frames_used),
        "n_bridged_gaps": int(bridged_mask.sum()),
        "n_split_gaps": int(splits.size),
    }
