#!/usr/bin/env python3
"""Scale-rectified, cone-of-influence-masked CWT spectral descriptors (v2).

Changes from v1:

* Power is rectified by scale (Liu, Liang & Weisberg 2007, J. Atmos. Oceanic
  Technol. 24:2093). PyWavelets coefficient power for a pure sinusoid grows in
  proportion to scale, which pushed most dominant frequencies to the lowest
  bin in v1. ``|W|^2 / scale`` gives equal-amplitude sinusoids equal power at
  every frequency.
* Each continuous segment is transformed separately. The time-averaged power
  at each frequency uses only samples inside the cone of influence (COI): at
  least ``sqrt(B) * C / f`` seconds from a segment edge for wavelet
  ``cmorB-C``. Valid samples are pooled across segments.
* Frequency axis (2026-11-02 revision): 0.005-0.08 Hz, 80 linear bins.
  - Top: cmor1.5-1.0 has a frequency spread of sigma = f/(pi*sqrt(2B)) ~ 0.18 f,
    so a wavelet centred at f reaches ~1.37 f at 2 sigma. 0.08 Hz keeps that
    below the lowest Nyquist frequency in the data set (0.111 Hz at the
    ~4.5 s frame interval of three embryos). Measured calcium events last
    ~8 frames (~30 s) and carry essentially no power above 0.06 Hz.
  - Bottom: event rates in these embryos are ~0.005-0.02 Hz, so the floor
    was lowered from 0.02 Hz to 0.005 Hz so slow activity is measured rather
    than piling up at the band edge. Cost: ~93% of template-mappable cells are
    kept (vs ~97% at f_min >= 0.01 Hz); the extra exclusions are cells with
    < ~11.5 min of usable tracking, too short to resolve a 200 s period.
    Each segment is linearly detrended in 01_preprocess_data.py to remove
    photobleaching drift.
  Override with CALCIUM_FMIN, CALCIUM_FMAX and CALCIUM_NFREQ.
* Minimum data rule: a cell is kept only if, at every frequency, its pooled
  cone-of-influence-valid time is at least CALCIUM_MIN_VALID_CYCLES (default 1)
  full cycles. The lowest frequency is the binding one (1 cycle at 0.005 Hz =
  200 s). Cells that fail are listed in skipped_cells.csv with the reason,
  instead of appearing in the results table with all-NaN features.
"""
from __future__ import annotations

import math
import os
import re

import numpy as np

WAVELET = os.environ.get("CALCIUM_PYWT_WAVELET", "cmor1.5-1.0")
FMIN = float(os.environ.get("CALCIUM_FMIN", "0.005"))
FMAX = float(os.environ.get("CALCIUM_FMAX", "0.08"))
NFREQ = int(os.environ.get("CALCIUM_NFREQ", "80"))
MIN_VALID_CYCLES = float(os.environ.get("CALCIUM_MIN_VALID_CYCLES", "1"))
FEATURES = ("DominantFreqHz", "SpectralCentroid", "SpectralSpread",
            "SpectralEntropy", "TotalBandPower", "PeakPowerFraction")


def frequency_axis() -> np.ndarray:
    return np.linspace(FMIN, FMAX, NFREQ)


def coi_factor(wavelet_name: str = WAVELET) -> float:
    """COI e-folding time x frequency (s*Hz) for a PyWavelets cmorB-C wavelet.

    psi(t) ~ exp(-t^2/B) exp(i 2 pi C t). The amplitude envelope has
    sigma = sqrt(B/2), so the Torrence & Compo COI is sqrt(2)*sigma*scale =
    sqrt(B)*scale wavelet-time units. In seconds this is sqrt(B)*C/f.
    """
    match = re.fullmatch(r"cmor(?P<b>[\d.]+)-(?P<c>[\d.]+)", wavelet_name)
    if match is None:
        raise ValueError(f"COI is defined here only for cmorB-C wavelets, not {wavelet_name}")
    return math.sqrt(float(match["b"])) * float(match["c"])


def mean_power_spectrum(segments: list[dict], axis: np.ndarray
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return COI-masked, scale-rectified mean power, valid-sample counts and
    valid seconds per frequency."""
    import pywt

    wavelet = pywt.ContinuousWavelet(WAVELET)
    center = pywt.central_frequency(wavelet)
    coi_seconds = coi_factor() / axis
    power_sum = np.zeros(axis.size)
    counts = np.zeros(axis.size, dtype=int)
    valid_seconds = np.zeros(axis.size)
    for segment in segments:
        signal = np.asarray(segment["dFF"], dtype=float)
        dt = float(segment["dt"])
        scales = center / (axis * dt)
        coefficients, _ = pywt.cwt(signal, scales, wavelet, sampling_period=dt, method="fft")
        power = np.abs(coefficients) ** 2 / scales[:, None]
        elapsed = np.arange(signal.size) * dt
        edge_distance = np.minimum(elapsed, elapsed[-1] - elapsed)
        valid = edge_distance[None, :] >= coi_seconds[:, None]
        power_sum += np.where(valid, power, 0.0).sum(axis=1)
        counts += valid.sum(axis=1)
        valid_seconds += valid.sum(axis=1) * dt
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_power = np.where(counts > 0, power_sum / np.maximum(counts, 1), np.nan)
    return mean_power, counts, valid_seconds


def spectral_features(data: dict, axis: np.ndarray
                      ) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """Six spectral descriptors with the same definitions as v1 and the MATLAB method."""
    axis = np.asarray(axis, dtype=float).reshape(-1)
    if axis.size == 0 or np.any(~np.isfinite(axis)) or np.any(axis <= 0):
        raise ValueError("frequency axis must contain finite positive frequencies")
    result = {name: float("nan") for name in FEATURES}
    mean_power, counts, valid_seconds = mean_power_spectrum(data["segments"], axis)
    spectrum = {"Frequency": axis, "MeanPower": mean_power, "ValidSamples": counts}
    # Every frequency must rest on enough data; otherwise the cell is skipped
    # (pipeline.py records the reason) rather than kept with NaN features.
    needed = MIN_VALID_CYCLES / axis
    short = valid_seconds < needed
    if np.any(short):
        i = int(np.flatnonzero(short)[0])
        raise ValueError(
            f"Only {valid_seconds[i]:.0f} s of cone-of-influence-valid data at "
            f"{axis[i]:.4f} Hz (< {needed[i]:.0f} s = {MIN_VALID_CYCLES:g} cycle(s))")
    mean_power = np.maximum(mean_power, 0.0)
    total_power = float(np.sum(mean_power))
    eps = np.finfo(float).eps
    if total_power > 0 and np.isfinite(total_power):
        probabilities = mean_power / (total_power + eps)
        peak_index = int(np.argmax(mean_power))
        centroid = float(np.sum(axis * mean_power) / (total_power + eps))
        result.update({
            "DominantFreqHz": float(axis[peak_index]),
            "SpectralCentroid": centroid,
            "SpectralSpread": float(np.sqrt(np.sum((axis - centroid) ** 2 * mean_power)
                                            / (total_power + eps))),
            "SpectralEntropy": float(-np.sum(probabilities * np.log2(probabilities + eps))),
            "TotalBandPower": total_power,
            "PeakPowerFraction": float(mean_power[peak_index] / (total_power + eps)),
        })
    return result, spectrum
