#!/usr/bin/env python3
"""Extract common-axis CWT power and spectral descriptors using PyWavelets."""
from __future__ import annotations

import numpy as np

WAVELET = "cmor1.5-1.0"


def spectral_features(data: dict, frequency_axis: np.ndarray
                      ) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """Equivalent of MATLAB extractSpectralFeatures for a PyWavelets CWT.

    Mean squared coefficient magnitude is calculated at each requested
    frequency, which directly places every cell on the shared frequency axis.
    """
    import pywt

    axis = np.asarray(frequency_axis, dtype=float).reshape(-1)
    if axis.size == 0 or np.any(~np.isfinite(axis)) or np.any(axis <= 0):
        raise ValueError("frequency_axis must contain finite positive frequencies")
    result = {name: float("nan") for name in (
        "DominantFreqHz", "SpectralCentroid", "SpectralSpread",
        "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
    )}
    signal = np.asarray(data["dFF"], dtype=float)
    if signal.size < 100:
        return result, {"Frequency": axis, "MeanPower": np.zeros_like(axis)}
    fs = float(data["Fs"])
    dt = 1.0 / fs
    wavelet = pywt.ContinuousWavelet(WAVELET)
    scales = pywt.central_frequency(wavelet) / (axis * dt)
    coefficients, _ = pywt.cwt(signal, scales, wavelet, sampling_period=dt, method="fft")
    mean_power = np.maximum(np.mean(np.abs(coefficients) ** 2, axis=1), 0.0)
    total_power = float(np.sum(mean_power))
    eps = np.finfo(float).eps
    if total_power > 0 and np.isfinite(total_power):
        probabilities = mean_power / (total_power + eps)
        peak_index = int(np.argmax(mean_power))
        centroid = float(np.sum(axis * mean_power) / (total_power + eps))
        result.update({
            "DominantFreqHz": float(axis[peak_index]),
            "SpectralCentroid": centroid,
            "SpectralSpread": float(np.sqrt(np.sum((axis - centroid) ** 2 * mean_power) /
                                             (total_power + eps))),
            "SpectralEntropy": float(-np.sum(probabilities * np.log2(probabilities + eps))),
            "TotalBandPower": total_power,
            "PeakPowerFraction": float(mean_power[peak_index] / (total_power + eps)),
        })
    return result, {"Frequency": axis, "MeanPower": mean_power}
