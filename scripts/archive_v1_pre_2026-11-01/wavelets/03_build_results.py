#!/usr/bin/env python3
"""Assemble a cell-level result row from preprocessing and spectral stages."""
from __future__ import annotations

from datetime import datetime

import numpy as np


def build_result_row(data: dict, spectral_feat: dict, cell_id: object,
                     source_file: str) -> dict:
    """Equivalent of MATLAB buildResultRow."""
    return {
        "CellID": str(cell_id),
        "RunID": datetime.now().strftime("%Y-%m-%d_%H%M%S"),
        "SourceFile": str(source_file),
        **spectral_feat,
        "Duration_min": float((data["time"][-1] - data["time"][0]) / 60.0),
        "MeanX": float(np.mean(data["X"])),
        "MeanY": float(np.mean(data["Y"])),
    }
