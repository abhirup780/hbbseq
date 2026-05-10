from hbb_pipeline.models import TraceData
from hbb_pipeline.qc import evaluate_trace_artifacts

def _make_trace(peak_positions: list[int], max_intensities: list[int], min_intensities: list[int]) -> TraceData:
    channels = {
        "A": [],
        "C": [],
        "G": [],
        "T": []
    }
    # Populate the "A" channel with the max intensity and "C", "G", "T" with min intensities
    for pos, max_val, min_val in zip(peak_positions, max_intensities, min_intensities):
        while len(channels["A"]) <= pos:
            channels["A"].append(0)
            channels["C"].append(0)
            channels["G"].append(0)
            channels["T"].append(0)
        channels["A"][pos] = max_val
        channels["C"][pos] = min_val
        channels["G"][pos] = min_val
        channels["T"][pos] = min_val
        
    return TraceData(
        sequence="A" * len(peak_positions),
        phred_scores=[20] * len(peak_positions),
        channels=channels,
        peak_positions=peak_positions,
        sample_name="test_qc"
    )

def test_dye_blob_detection() -> None:
    peaks = list(range(10, 2500, 10))  # 249 peaks
    max_vals = [1000] * 249
    min_vals = [50] * 249
    
    # Spike at index 50 (in first 100 bases)
    max_vals[50] = 6000
    
    trace = _make_trace(peaks, max_vals, min_vals)
    warnings = evaluate_trace_artifacts(trace)
    assert any("dye blob" in w for w in warnings)

def test_high_baseline_noise() -> None:
    peaks = list(range(10, 2500, 10))  # 249 peaks
    max_vals = [1000] * 249
    # Very high minimums (200 / 1000 = 0.20 > 0.15 threshold)
    min_vals = [200] * 249
    
    trace = _make_trace(peaks, max_vals, min_vals)
    warnings = evaluate_trace_artifacts(trace)
    assert any("baseline noise" in w for w in warnings)

def test_signal_drop_off() -> None:
    peaks = list(range(10, 2500, 10))  # 249 peaks
    max_vals = [1000] * 249
    min_vals = [50] * 249
    
    # Drop off last 100 bases
    for i in range(149, 249):
        max_vals[i] = 50
        
    trace = _make_trace(peaks, max_vals, min_vals)
    warnings = evaluate_trace_artifacts(trace)
    assert any("drop-off" in w for w in warnings)

def test_clean_trace() -> None:
    peaks = list(range(10, 2500, 10))  # 249 peaks
    max_vals = [1000] * 249
    min_vals = [50] * 249
    
    trace = _make_trace(peaks, max_vals, min_vals)
    warnings = evaluate_trace_artifacts(trace)
    assert not warnings
