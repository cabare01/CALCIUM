#!/usr/bin/env python3
"""Assemble a cell-level result row from preprocessing and spectral stages (v2)."""
from __future__ import annotations

import numpy as np

RESULT_COLUMNS = [
    "CellID", "RunID", "SourceFile", "DominantFreqHz", "SpectralCentroid",
    "SpectralSpread", "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
    "Duration_min", "MeanX", "MeanY",
    # v2 QC columns
    "UsedDuration_min", "FrameInterval_s", "FramesObserved", "FramesUsed",
    "NSegments", "NBridgedGaps", "NSplitGaps",
]


def build_result_row(data: dict, spectral_feat: dict, cell_id: object,
                     source_file: str, run_id: str) -> dict:
    """Equivalent of MATLAB buildResultRow plus segment/gap QC fields.

    MeanX/MeanY are averaged over the resampled samples of the used segments.
    Duration_min is the elapsed span of the whole track, as in v1.
    """
    segments = data["segments"]
    x = np.concatenate([s["X"] for s in segments])
    y = np.concatenate([s["Y"] for s in segments])
    used = sum(float(s["time"][-1] - s["time"][0]) for s in segments)
    return {
        "CellID": str(cell_id),
        "RunID": run_id,
        "SourceFile": str(source_file),
        **spectral_feat,
        "Duration_min": data["span_seconds"] / 60.0,
        "MeanX": float(np.mean(x)),
        "MeanY": float(np.mean(y)),
        "UsedDuration_min": used / 60.0,
        "FrameInterval_s": data["frame_dt"],
        "FramesObserved": data["frames_observed"],
        "FramesUsed": data["frames_used"],
        "NSegments": len(segments),
        "NBridgedGaps": data["n_bridged_gaps"],
        "NSplitGaps": data["n_split_gaps"],
    }
