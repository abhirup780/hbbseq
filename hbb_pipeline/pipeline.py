"""Shared trace-processing core used by both the Streamlit app and the CLI.

Accepts an already-parsed TraceData and returns the trimmed, (optionally
reverse-complemented) trace together with its alignment.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from hbb_pipeline.alignment import (
    align_to_reference,
    apply_iupac_symbols,
    reverse_complement_trace,
)
from hbb_pipeline.parsing import trim_by_quality

if TYPE_CHECKING:
    from hbb_pipeline.models import AlignedRead, TraceData
    from hbb_pipeline.reference import HBBReference

logger = logging.getLogger(__name__)


def process_trace(
    raw: "TraceData",
    ref: "HBBReference",
    is_reverse: bool,
    min_phred: int = 20,
    het_ratio: float = 0.25,
) -> tuple["TraceData", "AlignedRead | None", int, int, str | None]:
    """Apply IUPAC coding, quality trim, optional RC, and align to reference.

    Args:
        raw: Parsed (untrimmed) TraceData.
        ref: Loaded HBBReference.
        is_reverse: If True, reverse-complement after trimming.
        min_phred: Minimum Phred quality used for Mott trimming threshold.
        het_ratio: Secondary-peak ratio forwarded to the aligner.

    Returns:
        (processed_trace, aligned_read, trim_start, trim_end, error_msg)
        error_msg is None on success; aligned_read is None on alignment failure.
    """
    try:
        iupac = apply_iupac_symbols(raw, cutoff=het_ratio)
        threshold = 10 ** (-min_phred / 10.0)
        trimmed, (trim_s, trim_e) = trim_by_quality(iupac, threshold=threshold)
        if is_reverse:
            trimmed = reverse_complement_trace(trimmed)
        aligned = align_to_reference(trimmed, ref)
        return trimmed, aligned, trim_s, trim_e, None
    except Exception as exc:
        logger.warning("Trace processing failed for %s: %s", raw.sample_name, exc)
        return raw, None, 0, len(raw.sequence), str(exc)
