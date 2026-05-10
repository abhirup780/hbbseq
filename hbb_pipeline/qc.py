from __future__ import annotations

import statistics
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hbb_pipeline.models import TraceData

_MIN_PEAKS       = 200    # traces shorter than this are skipped
_EDGE_WINDOW     = 100    # bases at each end used for dye-blob and drop-off checks
_DYE_BLOB_RATIO  = 5.0    # early-peak / mid-median ratio threshold
_NOISE_RATIO     = 0.15   # avg-min / avg-max ratio threshold (SNR proxy)
_DROPOFF_RATIO   = 0.10   # tail-avg / head-avg ratio threshold


def evaluate_trace_artifacts(trace: "TraceData") -> list[str]:
    """Evaluate raw trace data for common artifacts.
    
    Checks for:
    1. Dye blobs (massive spikes in the first 100 bases)
    2. High baseline noise (SNR proxy)
    3. Signal drop-off (late trace degradation)
    
    Returns a list of human-readable warning strings.
    """
    warnings = []
    
    if not trace.channels or not trace.peak_positions:
        return warnings
        
    n_peaks = len(trace.peak_positions)
    if n_peaks < _MIN_PEAKS:
        return warnings
        
    # Pre-calculate max and min intensities at each peak position
    max_peaks = []
    min_peaks = []
    
    for scan_pos in trace.peak_positions:
        intensities = []
        for base in "ACGT":
            ch = trace.channels.get(base, [])
            if scan_pos < len(ch):
                intensities.append(ch[scan_pos])
            else:
                intensities.append(0)
        max_peaks.append(max(intensities) if intensities else 0)
        min_peaks.append(min(intensities) if intensities else 0)
        
    avg_peak = statistics.mean(max_peaks)
    if avg_peak == 0:
        return warnings
        
    # 1. Dye blob: spike in first _EDGE_WINDOW bases vs median of the rest
    first_max  = max(max_peaks[:_EDGE_WINDOW])
    rest_median = statistics.median(max_peaks[_EDGE_WINDOW:])
    if rest_median > 0 and (first_max / rest_median) > _DYE_BLOB_RATIO:
        warnings.append("Massive early signal peak detected (possible dye blob, primer dimer, or excess template). Early sequence may be obscured.")

    # 2. High baseline noise: avg-min / avg-max ratio
    avg_min = statistics.mean(min_peaks)
    if (avg_min / avg_peak) > _NOISE_RATIO:
        warnings.append("High baseline noise detected. Signal-to-noise ratio is poor, which may cause false positive heterozygous calls.")

    # 3. Signal drop-off: tail avg < _DROPOFF_RATIO of head avg
    first_avg = statistics.mean(max_peaks[:_EDGE_WINDOW])
    last_avg  = statistics.mean(max_peaks[-_EDGE_WINDOW:])
    if first_avg > 0 and (last_avg / first_avg) < _DROPOFF_RATIO:
        warnings.append("Severe signal drop-off detected toward the end of the trace.")
        
    return warnings
