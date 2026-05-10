"""Streamlit UI for HBBseq — β-Globin Sanger variant analysis pipeline."""

from __future__ import annotations

import hashlib
import logging
import sys
import tempfile
from pathlib import Path

from PIL import Image, ImageDraw

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent))

from plots import plot_chromatogram, plot_coverage_map
from hbb_pipeline.alignment import (
    align_to_reference,
    apply_iupac_symbols,
    build_consensus,
    reverse_complement_trace,
)
from hbb_pipeline.heterozygosity import detect_het_indel_breakpoint
from hbb_pipeline.coordinates import CoordinateTranslator
from hbb_pipeline.models import TraceData
from hbb_pipeline.parsing import parse_abi
from hbb_pipeline.pipeline import process_trace as _process_trace_core
from hbb_pipeline.qc import evaluate_trace_artifacts
from hbb_pipeline.reference import HBBReference
from hbb_pipeline.reporting import generate_report, render_markdown_report
from hbb_pipeline.known_variants import lookup_variant
from hbb_pipeline.variants import call_variants_from_alignment

logger = logging.getLogger(__name__)
_DEFAULT_REF = Path(__file__).parent / "reference" / "HBB_reference.fasta"

# Bump this whenever pipeline code changes that affect cached results.
_PIPELINE_VERSION = "35"


# ---------------------------------------------------------------------------
# Pipeline helpers
# ---------------------------------------------------------------------------

def _slice_trace(trace: TraceData, s: int, e: int) -> TraceData:
    return TraceData(
        sequence=trace.sequence[s:e],
        phred_scores=trace.phred_scores[s:e],
        channels=trace.channels,
        peak_positions=trace.peak_positions[s:e],
        sample_name=trace.sample_name,
        is_reverse_complemented=trace.is_reverse_complemented,
    )


def _cds_coverage(consensus: str, ref: HBBReference) -> str:
    covered = sum(
        1 for gs, ge in ref.exons
        for pos in range(gs, ge)
        if pos < len(consensus) and consensus[pos] not in ("?", "N")
    )
    total = sum(ge - gs for gs, ge in ref.exons)
    return f"{100 * covered / total:.1f}%" if total else "N/A"


def _resolve_idx(genomic_pos: int, aligned) -> int:
    from hbb_pipeline.variants import _resolve_trace_index
    return _resolve_trace_index(genomic_pos, aligned)  # -1 = not covered


def _process_trace(
    raw_bytes: bytes,
    name: str,
    tmpdir: str,
    ref: HBBReference,
    is_reverse: bool,
    min_phred: int = 20,
    het_ratio: float = 0.25,
) -> tuple[TraceData | None, TraceData | None, object | None, int, int, str | None]:
    """Parse bytes into a TraceData, then delegate to the shared processing core.

    Returns (raw, trimmed, aligned, trim_s, trim_e, error_msg).
    On failure returns (None, None, None, 0, 0, error_message).
    """
    try:
        path = Path(tmpdir) / name
        path.write_bytes(raw_bytes)
        raw = parse_abi(path)
    except Exception as exc:
        logger.warning("Trace %s failed to parse: %s", name, exc)
        return None, None, None, 0, 0, f"{type(exc).__name__}: {exc}"

    trimmed, aligned, s, e, err = _process_trace_core(
        raw, ref, is_reverse=is_reverse, min_phred=min_phred, het_ratio=het_ratio,
    )
    if err:
        return None, None, None, 0, 0, err
    return raw, trimmed, aligned, s, e, None


def _run_variant_pipeline(
    ref: HBBReference,
    translator: CoordinateTranslator,
    fwd_trimmed: TraceData | None,
    rev_rc: TraceData | None,
    fwd_aligned,
    rev_aligned,
    het_ratio: float = 0.25,
) -> tuple[list, str]:
    """Run consensus → variant calling → merge.

    HET SNV detection relies on IUPAC codes already embedded in the traces
    by apply_iupac_symbols() (called in _process_trace before alignment).

    Returns (merged_variants, consensus_str).
    """
    consensus, cons_quals = build_consensus(fwd_aligned, rev_aligned, ref)

    alignment_variants = call_variants_from_alignment(
        consensus, cons_quals, ref, translator,
        fwd_trimmed, rev_rc, fwd_aligned, rev_aligned,
        het_ratio=het_ratio,
    )
    return alignment_variants, consensus


# ---------------------------------------------------------------------------
# Cached pipeline entry points
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Running analysis…")
def run_pipeline(
    fwd_bytes: bytes,
    rev_bytes: bytes,
    fwd_name: str,
    rev_name: str,
    ref_path_str: str,
    min_phred: int,
    het_ratio: float,
) -> dict:
    """Paired-trace pipeline with graceful single-strand degradation.

    If one trace fails to parse or align, the analysis continues on the
    remaining strand and a warning is surfaced in the result.
    """
    try:
        ref = HBBReference(Path(ref_path_str))
        translator = CoordinateTranslator(ref)
        warnings: list[str] = []

        with tempfile.TemporaryDirectory() as tmpdir:
            fwd_raw, fwd_trimmed, fwd_aligned, fwd_s, fwd_e, fwd_err = \
                _process_trace(fwd_bytes, fwd_name, tmpdir, ref, is_reverse=False, min_phred=min_phred, het_ratio=het_ratio)
            rev_raw, rev_rc, rev_aligned, rev_s, rev_e, rev_err = \
                _process_trace(rev_bytes, rev_name, tmpdir, ref, is_reverse=True, min_phred=min_phred, het_ratio=het_ratio)

            if fwd_err:
                warnings.append(f"Forward trace failed ({fwd_err}) — analysis uses reverse strand only")
            elif fwd_raw:
                warnings.extend([f"[Forward] {w}" for w in evaluate_trace_artifacts(fwd_raw)])
                
            if rev_err:
                warnings.append(f"Reverse trace failed ({rev_err}) — analysis uses forward strand only")
            elif rev_raw:
                warnings.extend([f"[Reverse] {w}" for w in evaluate_trace_artifacts(rev_raw)])

            if fwd_aligned is None and rev_aligned is None:
                return {"error": "Both traces failed to process.\n" + "\n".join(warnings)}

            merged, consensus = _run_variant_pipeline(
                ref, translator,
                fwd_trimmed, rev_rc,
                fwd_aligned, rev_aligned,
                het_ratio=het_ratio,
            )

            # Het INDEL suspicion: scan each raw (untrimmed) IUPAC trace for a
            # sustained clean→noisy peak transition — the hallmark of a het INDEL.
            # Must use the untrimmed trace because Mott trimming removes the very
            # low-quality noisy region that is the primary evidence for the breakpoint.
            for strand_label, strand_raw, is_rev in (
                ("Forward", fwd_raw, False),
                ("Reverse", rev_raw, True),
            ):
                if strand_raw is None:
                    continue
                try:
                    iupac_full = apply_iupac_symbols(strand_raw, cutoff=het_ratio)
                    if is_rev:
                        iupac_full = reverse_complement_trace(iupac_full)
                    iupac_aligned = align_to_reference(iupac_full, ref)
                    if iupac_aligned is None:
                        continue
                    bp = detect_het_indel_breakpoint(iupac_full, iupac_aligned, translator)
                except Exception:
                    bp = None
                if bp is not None:
                    warnings.append(
                        f"⚠ Possible heterozygous INDEL ({strand_label}): mixed peaks "
                        f"detected from approximately c.{bp} onwards in the {strand_label.lower()} "
                        f"trace. This may indicate a heterozygous insertion or deletion "
                        f"in one allele. Please review the trace manually."
                    )

            qc: dict = {}
            if fwd_raw is not None:
                qc["mean_phred_fwd"] = round(
                    sum(fwd_raw.phred_scores) / max(len(fwd_raw.phred_scores), 1), 1)
                qc["usable_length_fwd"] = fwd_e - fwd_s
            else:
                qc["mean_phred_fwd"] = "N/A (trace failed)"
                qc["usable_length_fwd"] = "N/A"
            if rev_raw is not None:
                qc["mean_phred_rev"] = round(
                    sum(rev_raw.phred_scores) / max(len(rev_raw.phred_scores), 1), 1)
                qc["usable_length_rev"] = rev_e - rev_s
            else:
                qc["mean_phred_rev"] = "N/A (trace failed)"
                qc["usable_length_rev"] = "N/A"
            qc["cds_coverage_pct"] = _cds_coverage(consensus, ref)
            for w in warnings:
                logger.warning("[paired] %s", w)
            if warnings:
                qc["analysis_warnings"] = "; ".join(warnings)

            sample_id = fwd_name.replace(".ab1", "").replace(".abi", "")
            report = generate_report(merged, qc, sample_id, consensus)
            md = render_markdown_report(report)

        return {
            "report": report,
            "markdown": md,
            "fwd_trace": fwd_trimmed,
            "rev_rc": rev_rc,
            "fwd_aligned": fwd_aligned,
            "rev_aligned": rev_aligned,
            "fwd_phred_raw": fwd_raw.phred_scores if fwd_raw else None,
            "rev_phred_raw": rev_raw.phred_scores if rev_raw else None,
            "warnings": warnings,
            "error": None,
        }
    except Exception as exc:
        import traceback as _tb
        logger.error("Pipeline error: %s", _tb.format_exc())
        return {"error": f"{type(exc).__name__}: {exc}"}


@st.cache_data(show_spinner="Running analysis…")
def run_pipeline_single(
    trace_bytes: bytes,
    trace_name: str,
    is_reverse: bool,
    ref_path_str: str,
    min_phred: int,
    het_ratio: float,
) -> dict:
    """Single-trace pipeline — one forward or one reverse file."""
    try:
        ref = HBBReference(Path(ref_path_str))
        translator = CoordinateTranslator(ref)

        with tempfile.TemporaryDirectory() as tmpdir:
            raw, trace_proc, aligned, s, e, err = \
                _process_trace(trace_bytes, trace_name, tmpdir, ref, is_reverse=is_reverse, min_phred=min_phred, het_ratio=het_ratio)

            if err:
                return {"error": f"Trace failed to process: {err}"}
                
            if raw:
                strand_str = "Reverse" if is_reverse else "Forward"
                # Note: 'warnings' was not initialized here; fixed by adding it
                single_warnings: list[str] = [f"[{strand_str}] {w}" for w in evaluate_trace_artifacts(raw)]
            else:
                single_warnings = []

            if is_reverse:
                fwd_trimmed, rev_rc = None, trace_proc
                fwd_aligned, rev_aligned = None, aligned
            else:
                fwd_trimmed, rev_rc = trace_proc, None
                fwd_aligned, rev_aligned = aligned, None

            merged, consensus = _run_variant_pipeline(
                ref, translator,
                fwd_trimmed, rev_rc,
                fwd_aligned, rev_aligned,
                het_ratio=het_ratio,
            )

            # Het INDEL suspicion check — use untrimmed IUPAC trace so that the
            # low-quality noisy region after the breakpoint is not Mott-trimmed away.
            try:
                iupac_full = apply_iupac_symbols(raw, cutoff=het_ratio)
                if is_reverse:
                    iupac_full = reverse_complement_trace(iupac_full)
                iupac_aligned = align_to_reference(iupac_full, ref)
                bp = detect_het_indel_breakpoint(iupac_full, iupac_aligned, translator) if iupac_aligned else None
            except Exception:
                bp = None
            if bp is not None:
                single_warnings.append(
                    f"⚠ Possible heterozygous INDEL: mixed peaks detected from "
                    f"approximately c.{bp} onwards. This may indicate a heterozygous "
                    f"insertion or deletion in one allele. Please review the trace manually."
                )

            strand_label = "reverse" if is_reverse else "forward"
            qc: dict = {
                f"mean_phred_{strand_label}": round(
                    sum(raw.phred_scores) / max(len(raw.phred_scores), 1), 1),
                f"usable_length_{strand_label}": e - s,
                "cds_coverage_pct": _cds_coverage(consensus, ref),
                "analysis_mode": f"Single strand ({strand_label})",
            }
            for w in single_warnings:
                logger.warning("[single] %s", w)
            if single_warnings:
                qc["analysis_warnings"] = "; ".join(single_warnings)

            sample_id = trace_name.replace(".ab1", "").replace(".abi", "")
            report = generate_report(merged, qc, sample_id, consensus)
            md = render_markdown_report(report)

        return {
            "report": report,
            "markdown": md,
            "fwd_trace": fwd_trimmed,
            "rev_rc": rev_rc,
            "fwd_aligned": fwd_aligned,
            "rev_aligned": rev_aligned,
            "fwd_phred_raw": raw.phred_scores if not is_reverse else None,
            "rev_phred_raw": raw.phred_scores if is_reverse else None,
            "warnings": [
                f"Single-strand analysis ({strand_label} only). All variant calls require confirmation.",
                *single_warnings,
            ],
            "error": None,
        }
    except Exception as exc:
        import traceback as _tb
        logger.error("Pipeline error (single): %s", _tb.format_exc())
        return {"error": f"{type(exc).__name__}: {exc}"}


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------

def _inject_css() -> None:
    st.markdown("""<style>
/* ── Typography & base ───────────────────────────────────────────────── */
.app-title {
    font-size: 3rem; font-weight: 800; letter-spacing: -1.5px;
    color: #F1F5F9; margin: 0 0 6px; line-height: 1;
}
.app-title .title-hbb { color: #F1F5F9; }
.app-title .title-seq {
    color: #3B82F6; font-weight: 600; letter-spacing: -0.5px;
}
.app-subtitle {
    font-size: 0.75rem; color: #7B8FA8; margin: 0 0 28px;
    letter-spacing: 0.5px; font-weight: 600; text-transform: uppercase;
}
.sidebar-label {
    font-size: 0.66rem; font-weight: 700; letter-spacing: 1.5px;
    text-transform: uppercase; color: #7B8FA8; margin: 18px 0 6px; display: block;
}
.section-rule {
    border: none; border-top: 1px solid rgba(255,255,255,0.07); margin: 20px 0 16px;
}

/* ── Alert banners ───────────────────────────────────────────────────── */
.alert {
    border-radius: 5px; padding: 11px 16px; margin: 3px 0 8px;
    border-left: 3px solid; font-size: 0.855rem; line-height: 1.6;
}
.alert-warn { border-left-color: #92400E; background: rgba(120,80,20,0.10); color: #A8957A; }
.alert-pass { border-left-color: #10B981; background: rgba(16,185,129,0.09); color: #6EE7B7; }
.alert-info { border-left-color: #3B82F6; background: rgba(59,130,246,0.09); color: #93C5FD; }

/* ── Variant cards ───────────────────────────────────────────────────── */
.vcard {
    border-radius: 10px;
    padding: 22px 26px 20px 22px;
    margin: 0 0 16px;
    border-left: 4px solid #3F4F60;
    background: rgba(255,255,255,0.040);
    transition: background 0.18s;
}
.vcard:hover { background: rgba(255,255,255,0.058); }
.vcard.pathogenic {
    border-left-color: #EF4444;
    background: rgba(239,68,68,0.045);
    box-shadow: 0 0 0 1px rgba(239,68,68,0.12), 0 4px 20px rgba(239,68,68,0.06);
}
.vcard.pathogenic:hover { background: rgba(239,68,68,0.07); }
.vcard.review {
    border-left-color: #F59E0B;
    background: rgba(245,158,11,0.04);
}
.vcard.benign      { border-left-color: #3B82F6; }
.vcard.unclassified { border-left-color: #475569; }

.vcard-top { display: flex; align-items: baseline; flex-wrap: wrap; gap: 14px; margin-bottom: 6px; }
.vcard-hgvs {
    font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', ui-monospace, monospace;
    font-size: 1.38rem; font-weight: 700; color: #F1F5F9; letter-spacing: -0.3px;
}
.vcard-protein {
    font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', ui-monospace, monospace;
    font-size: 1rem; color: #64748B; font-weight: 400;
}
.vcard-name {
    font-size: 0.97rem; color: #94A3B8; margin: 5px 0 16px; font-weight: 500;
    letter-spacing: 0.1px;
}
.vcard-badges { display: flex; flex-wrap: wrap; gap: 7px; align-items: center; }

/* ── Badges ──────────────────────────────────────────────────────────── */
.badge {
    display: inline-flex; align-items: center;
    font-size: 0.69rem; font-weight: 700; letter-spacing: 0.7px;
    padding: 3px 10px; border-radius: 4px; text-transform: uppercase; white-space: nowrap;
}
.badge-path   { background: rgba(239,68,68,0.18);  color: #FCA5A5; border: 1px solid rgba(239,68,68,0.35); }
.badge-benign { background: rgba(59,130,246,0.18); color: #93C5FD; border: 1px solid rgba(59,130,246,0.35); }
.badge-mod    { background: rgba(245,158,11,0.18); color: #FCD34D; border: 1px solid rgba(245,158,11,0.35); }
.badge-review { background: rgba(245,158,11,0.15); color: #FCD34D; border: 1px solid rgba(245,158,11,0.32); }
.badge-neutral{ background: rgba(148,163,184,0.12); color: #64748B; border: 1px solid rgba(148,163,184,0.2); }
.badge-hom    { background: rgba(234,88,12,0.15);  color: #FB923C; border: 1px solid rgba(234,88,12,0.32); }
.badge-het    { background: rgba(139,92,246,0.15); color: #A78BFA; border: 1px solid rgba(139,92,246,0.32); }

/* ── Result summary line ─────────────────────────────────────────────── */
.result-header {
    font-size: 0.72rem; font-weight: 700; letter-spacing: 1.5px; text-transform: uppercase;
    margin: 4px 0 22px; padding-bottom: 12px;
    border-bottom: 1px solid rgba(255,255,255,0.07);
}
.result-header-count { color: #94A3B8; }
.result-header-sample { color: #475569; font-weight: 500; }

/* ── Metric boxes ────────────────────────────────────────────────────── */
[data-testid="metric-container"] {
    background: rgba(255,255,255,0.035) !important;
    border: 1px solid rgba(255,255,255,0.09) !important;
    border-radius: 7px !important;
    padding: 14px 18px !important;
}
[data-testid="stMetricLabel"] {
    font-size: 0.72rem !important; font-weight: 600 !important;
    letter-spacing: 0.8px !important; text-transform: uppercase !important;
    color: #64748B !important;
}
[data-testid="stMetricValue"] { color: #E2E8F0 !important; font-size: 1.35rem !important; }
[data-testid="stMetricDelta"] { font-size: 0.78rem !important; }

/* ── Sidebar ─────────────────────────────────────────────────────────── */
section[data-testid="stSidebar"] {
    border-right: 1px solid rgba(255,255,255,0.07);
}
/* Widget labels and help text — lift contrast from Streamlit's default */
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] .stWidgetLabel,
section[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {
    color: #C4CFDC !important;
}
section[data-testid="stSidebar"] .stFileUploaderDropzoneInstructions,
section[data-testid="stSidebar"] small,
section[data-testid="stSidebar"] .stCaption {
    color: #64748B !important;
}

/* ── Tab bar ─────────────────────────────────────────────────────────── */
[data-baseweb="tab-list"] {
    gap: 0;
    border-bottom: 1px solid rgba(255,255,255,0.08) !important;
}
[data-baseweb="tab"] {
    font-size: 0.83rem !important; font-weight: 500 !important;
    padding: 9px 18px !important; color: #64748B !important;
}
[data-baseweb="tab"]:hover { color: #CBD5E1 !important; }
[aria-selected="true"][data-baseweb="tab"] { color: #E2E8F0 !important; }
[data-baseweb="tab-highlight"] {
    background-color: #3B82F6 !important;
    height: 2px !important;
}

/* ── Rejected-calls table ────────────────────────────────────────────── */
.rejected-row {
    display: flex; align-items: center; padding: 7px 4px;
    border-bottom: 1px solid rgba(255,255,255,0.05); font-size: 0.83rem; color: #64748B;
}
.rejected-hgvs { font-family: ui-monospace, monospace; text-decoration: line-through; margin-right: 8px; }

/* ── Variant pills (st.pills selector) ──────────────────────────────── */
[data-testid="stPills"] {
    gap: 8px !important;
    flex-wrap: wrap !important;
}
[data-testid="stPills"] button {
    font-family: 'JetBrains Mono', 'Fira Code', ui-monospace, monospace !important;
    font-size: 0.82rem !important;
    font-weight: 500 !important;
    padding: 6px 16px !important;
    border-radius: 20px !important;
    border: 1px solid rgba(255,255,255,0.14) !important;
    background: rgba(255,255,255,0.05) !important;
    color: #94A3B8 !important;
    transition: all 0.15s ease !important;
    letter-spacing: 0.2px !important;
}
[data-testid="stPills"] button:hover {
    border-color: rgba(59,130,246,0.5) !important;
    background: rgba(59,130,246,0.1) !important;
    color: #CBD5E1 !important;
}
[data-testid="stPills"] button[aria-pressed="true"],
[data-testid="stPills"] button[aria-selected="true"] {
    background: rgba(59,130,246,0.18) !important;
    border-color: #3B82F6 !important;
    color: #93C5FD !important;
    font-weight: 600 !important;
}

/* ── Reject / dismiss: float ✕ button into card top-right corner ─────── */
/* Hook: target horizontal blocks whose first column contains our .vcard HTML */
div[data-testid="stHorizontalBlock"]:has(.vcard),
div[data-testid="stHorizontalBlock"]:has(.alert-warn) {
    align-items: flex-start !important;
    gap: 0 !important;
}
div[data-testid="stHorizontalBlock"]:has(.vcard) > div:first-child,
div[data-testid="stHorizontalBlock"]:has(.alert-warn) > div:first-child {
    flex: 1 1 auto !important;
    max-width: 100% !important;
}
div[data-testid="stHorizontalBlock"]:has(.vcard) > div:last-child,
div[data-testid="stHorizontalBlock"]:has(.alert-warn) > div:last-child {
    flex: 0 0 auto !important;
    margin-left: -52px;
    padding-top: 9px;
    z-index: 20;
    position: relative;
}
div[data-testid="stHorizontalBlock"]:has(.vcard) > div:last-child button,
div[data-testid="stHorizontalBlock"]:has(.alert-warn) > div:last-child button {
    padding: 3px 9px !important;
    min-height: 28px !important;
    width: 32px !important;
    background: rgba(22,22,32,0.72) !important;
    border: 1px solid rgba(255,255,255,0.13) !important;
    color: rgba(255,255,255,0.38) !important;
    font-size: 0.75rem !important;
    border-radius: 6px !important;
    line-height: 1 !important;
    backdrop-filter: blur(6px) !important;
}
div[data-testid="stHorizontalBlock"]:has(.vcard) > div:last-child button:hover,
div[data-testid="stHorizontalBlock"]:has(.alert-warn) > div:last-child button:hover {
    border-color: rgba(255,255,255,0.32) !important;
    color: rgba(255,255,255,0.72) !important;
    background: rgba(50,50,68,0.88) !important;
}

/* ── Variant review blocks (st.container border=True) ───────────────── */
[data-testid="stVerticalBlockBorderWrapper"] {
    border-color: rgba(255,255,255,0.09) !important;
    border-radius: 12px !important;
    margin-bottom: 22px !important;
    overflow: hidden !important;
    padding: 0 !important;
}
[data-testid="stVerticalBlockBorderWrapper"] > div {
    padding: 0 !important;
}
/* Card inside block: flush top, no radius (container clips it), thicker indicator */
[data-testid="stVerticalBlockBorderWrapper"] .vcard {
    margin: 0 !important;
    border-radius: 0 !important;
    border-left-width: 5px !important;
    border-bottom: 1px solid rgba(255,255,255,0.07) !important;
    padding: 20px 24px 18px 20px !important;
}
[data-testid="stVerticalBlockBorderWrapper"] .vcard.pathogenic {
    box-shadow: none !important;
}
/* Chromatogram sub-section padding */
[data-testid="stVerticalBlockBorderWrapper"] > div > div[data-testid="stVerticalBlock"] > *:not([data-testid="stMarkdownContainer"]:first-child) {
    padding-left: 6px;
    padding-right: 6px;
}
.chrom-section-pad { padding: 14px 20px 16px; }

/* ── Sidebar brand ───────────────────────────────────────────────────── */
.sidebar-brand { padding: 6px 0 20px; }
.sb-logo {
    font-size: 1.7rem; font-weight: 800; letter-spacing: -1px;
    line-height: 1; margin-bottom: 4px;
}
.sb-hbb { color: #F1F5F9; }
.sb-seq { color: #3B82F6; font-weight: 600; }
.sb-sub {
    font-size: 0.60rem; color: #475569;
    letter-spacing: 1.5px; text-transform: uppercase; font-weight: 600;
}

/* ── Landing page ────────────────────────────────────────────────────── */
[data-testid="stMarkdownContainer"]:has(.landing-wrap) {
    width: 100%;
    text-align: center;
}
.landing-wrap {
    width: 100%; margin: 52px 0 0;
    display: flex; flex-direction: column; align-items: center;
}
.landing-title {
    font-size: 5.2rem; font-weight: 800; letter-spacing: -4px;
    line-height: 1; margin-bottom: 18px; color: #F1F5F9;
    text-align: center; width: 100%;
}
.landing-title .lt-seq { color: #3B82F6; font-weight: 600; letter-spacing: -1px; }
.landing-tagline {
    font-size: 1.05rem; color: #64748B; line-height: 1.75;
    margin-bottom: 52px; font-weight: 400;
    text-align: center; width: 100%;
}
.landing-step {
    background: rgba(255,255,255,0.025);
    border: 1px solid rgba(255,255,255,0.07);
    border-radius: 14px; padding: 28px 22px 26px; text-align: center;
    transition: border-color 0.2s;
}
.landing-step:hover { border-color: rgba(59,130,246,0.25); }
.landing-step-num {
    font-size: 2rem; font-weight: 800; letter-spacing: -1px;
    color: rgba(59,130,246,0.25); line-height: 1; margin-bottom: 16px;
}
.landing-step-title {
    font-size: 0.95rem; font-weight: 600; color: #E2E8F0; margin-bottom: 8px;
}
.landing-step-desc { font-size: 0.82rem; color: #475569; line-height: 1.65; }
/* ── Full chromatogram tab section headers ───────────────────────────── */
.trace-section-header {
    display: flex; align-items: center; gap: 9px;
    font-size: 0.78rem; font-weight: 600; letter-spacing: 0.5px;
    text-transform: uppercase; color: #64748B;
    padding: 6px 2px 4px; margin-bottom: -6px;
    border-bottom: 1px solid rgba(255,255,255,0.06);
}
.trace-section-dot {
    width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0;
}
.dot-fwd { background: #3B82F6; }
.dot-rev { background: #A78BFA; }

</style>""", unsafe_allow_html=True)


def _variant_card_html(v, sig: str | None) -> str:
    """Build an HTML variant card (no Q score, badge-based metadata)."""
    hgvs    = f"HBB:{v.hgvs_c}"
    protein = v.hgvs_p_hgvs or ""
    zyg     = v.zygosity.value.capitalize()
    name    = v.known_variant_name or ""
    region  = v.region

    reason     = v.review_reason or ""
    is_fwd_only = "Forward-strand only" in reason
    is_rev_only = "Reverse-strand only" in reason

    if is_fwd_only or is_rev_only or v.requires_manual_review:
        card_cls = "review"
    elif sig == "Pathogenic":
        card_cls = "pathogenic"
    elif sig == "Benign":
        card_cls = "benign"
    else:
        card_cls = "unclassified"

    if sig == "Pathogenic":
        sig_html = '<span class="badge badge-path">Pathogenic</span>'
    elif sig == "Benign":
        sig_html = '<span class="badge badge-benign">Benign</span>'
    elif sig == "Modifier":
        sig_html = '<span class="badge badge-mod">Modifier</span>'
    else:
        sig_html = '<span class="badge badge-neutral">Unclassified</span>'

    if is_fwd_only:
        review_html = '<span class="badge badge-review">Forward strand only</span>'
    elif is_rev_only:
        review_html = '<span class="badge badge-review">Reverse strand only</span>'
    elif v.requires_manual_review:
        review_html = '<span class="badge badge-review">Manual review</span>'
    else:
        review_html = ""

    zyg_cls = (
        "badge-hom" if zyg == "Homozygous"
        else "badge-het" if zyg == "Heterozygous"
        else "badge-neutral"
    )

    protein_html = (
        f'<div style="font-family:ui-monospace,monospace;font-size:0.9rem;'
        f'color:#475569;margin:2px 0 0;letter-spacing:0.1px;">{protein}</div>'
    ) if protein else ""
    name_html = (
        f'<div class="vcard-name">{name}</div>'
    ) if name else '<div style="margin-bottom:14px;"></div>'

    return (
        f'<div class="vcard {card_cls}">'
        f'<div class="vcard-hgvs">{hgvs}</div>'
        f'{protein_html}'
        f'{name_html}'
        f'<div class="vcard-badges">'
        f'<span class="badge {zyg_cls}">{zyg}</span>'
        f'<span class="badge badge-neutral">{region}</span>'
        f'{sig_html}{review_html}'
        f'</div></div>'
    )


def _sig_for_variant(v) -> str | None:
    """Return clinical significance string for a variant, or None."""
    key = v.hgvs_c[2:] if v.hgvs_c.startswith("c.") else v.hgvs_c
    kv = lookup_variant(key)
    return kv.clinical_significance if kv else None


def _render_variant_block(
    v,
    sig: str | None,
    fwd_trace,
    rev_rc,
    fwd_aligned,
    rev_aligned,
) -> None:
    """Render one variant's card + chromatograms in an integrated bordered block."""
    gpos = v.ref_pos_genomic
    with st.container(border=True):
        st.markdown(_variant_card_html(v, sig), unsafe_allow_html=True)
        st.markdown('<div class="chrom-section-pad">', unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        with c1:
            st.markdown('<span class="sidebar-label">Forward strand</span>',
                        unsafe_allow_html=True)
            if fwd_trace is None:
                st.markdown(
                    '<div class="alert alert-info" style="margin:6px 0 0">Forward trace not available.</div>',
                    unsafe_allow_html=True)
            else:
                tidx = _resolve_idx(gpos, fwd_aligned)
                if tidx < 0:
                    st.markdown(
                        '<div class="alert alert-info" style="margin:6px 0 0">Not covered by forward read.</div>',
                        unsafe_allow_html=True)
                else:
                    fig = plot_chromatogram(
                        fwd_trace,
                        max(0, tidx - 20),
                        min(len(fwd_trace.sequence), tidx + 21),
                        highlight_pos=tidx,
                    )
                    st.plotly_chart(fig, key=f"chrom_fwd_{gpos}", use_container_width=True)
        with c2:
            st.markdown('<span class="sidebar-label">Reverse complement</span>',
                        unsafe_allow_html=True)
            if rev_rc is None:
                st.markdown(
                    '<div class="alert alert-info" style="margin:6px 0 0">Reverse complement trace not available.</div>',
                    unsafe_allow_html=True)
            else:
                tidx_r = _resolve_idx(gpos, rev_aligned)
                if tidx_r < 0:
                    st.markdown(
                        '<div class="alert alert-info" style="margin:6px 0 0">Not covered by reverse complement read.</div>',
                        unsafe_allow_html=True)
                else:
                    fig = plot_chromatogram(
                        rev_rc,
                        max(0, tidx_r - 20),
                        min(len(rev_rc.sequence), tidx_r + 21),
                        highlight_pos=tidx_r,
                    )
                    st.plotly_chart(fig, key=f"chrom_rev_{gpos}", use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Favicon
# ---------------------------------------------------------------------------

def _make_favicon() -> Image.Image:
    """Blue rounded-square with a white bold H — renders cleanly at all favicon sizes."""
    s = 64
    blue = (59, 130, 246, 255)
    white = (255, 255, 255, 255)

    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    # Background: blue rounded square
    try:
        draw.rounded_rectangle([0, 0, s - 1, s - 1], radius=13, fill=blue)
    except AttributeError:  # Pillow < 8.2
        draw.rectangle([0, 0, s - 1, s - 1], fill=blue)

    # H shape — two vertical bars + crossbar
    pad, bar, mid_h = 14, 10, 4
    mid = s // 2
    draw.rectangle([pad,          pad, pad + bar,      s - pad], fill=white)  # left bar
    draw.rectangle([s - pad - bar, pad, s - pad,        s - pad], fill=white)  # right bar
    draw.rectangle([pad,          mid - mid_h, s - pad, mid + mid_h], fill=white)  # crossbar

    return img


# ---------------------------------------------------------------------------
# Main UI
# ---------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(
        page_title="HBBseq",
        page_icon=_make_favicon(),
        layout="wide",
        initial_sidebar_state="expanded",
    )
    _inject_css()
    # Override Streamlit Cloud's "· Streamlit" title suffix
    st.components.v1.html(
        "<script>window.parent.document.title = 'HBBseq';</script>",
        height=0,
    )

    # ── Sidebar ───────────────────────────────────────────────────────────────
    with st.sidebar:
        st.markdown(
            '<div class="sidebar-brand">'
            '<div class="sb-logo"><span class="sb-hbb">HBB</span><span class="sb-seq">seq</span></div>'
            '<div class="sb-sub">Beta-Globin &nbsp;·&nbsp; Sanger Analysis</div>'
            '</div>',
            unsafe_allow_html=True,
        )
        st.markdown('<hr class="section-rule" style="margin-top:0">', unsafe_allow_html=True)
        st.markdown('<span class="sidebar-label">Input Mode</span>', unsafe_allow_html=True)
        mode = st.radio(
            "mode", ["Paired (Fwd + Rev)", "Single Strand"],
            horizontal=True, label_visibility="collapsed",
        )

        st.markdown('<span class="sidebar-label">Trace Files</span>', unsafe_allow_html=True)
        if mode == "Paired (Fwd + Rev)":
            fwd_file = st.file_uploader("Forward trace", type=["ab1", "abi"],
                                        label_visibility="visible")
            rev_file = st.file_uploader("Reverse trace", type=["ab1", "abi"],
                                        label_visibility="visible")
            single_file = None
            is_reverse  = False
            can_run     = fwd_file is not None and rev_file is not None
        else:
            strand_dir = st.radio("Strand", ["Forward", "Reverse"], horizontal=True)
            is_reverse  = strand_dir == "Reverse"
            single_file = st.file_uploader(
                f"{'Reverse' if is_reverse else 'Forward'} trace",
                type=["ab1", "abi"], label_visibility="visible",
            )
            fwd_file = rev_file = None
            can_run = single_file is not None

        st.markdown('<hr class="section-rule">', unsafe_allow_html=True)
        run_btn = st.button(
            "Run Analysis", type="primary",
            use_container_width=True, disabled=not can_run,
        )
        if not can_run:
            st.caption(
                "Upload both trace files to run." if mode == "Paired (Fwd + Rev)"
                else "Upload a trace file to run."
            )

        st.markdown('<hr class="section-rule">', unsafe_allow_html=True)
        st.markdown('<span class="sidebar-label">Analysis Parameters</span>', unsafe_allow_html=True)
        min_phred = st.slider("Trim quality (Phred)", 10, 30, 20)
        het_ratio = st.slider("HET peak ratio", 0.15, 0.50, 0.25, step=0.01)

    ref_path_str = str(_DEFAULT_REF)

    if not can_run:
        st.markdown(
            '<div class="landing-wrap">'
            '<div class="landing-title"><span>HBB</span><span class="lt-seq">seq</span></div>'
            '<div class="landing-tagline">'
            'Variant detection and HGVS annotation for <em>HBB</em> Sanger sequencing data.'
            '</div>'
            '</div>',
            unsafe_allow_html=True,
        )
        c1, c2, c3 = st.columns(3, gap="large")
        c1.markdown(
            '<div class="landing-step">'
            '<div class="landing-step-num">01</div>'
            '<div class="landing-step-title">Upload traces</div>'
            '<div class="landing-step-desc">Add your forward and reverse .ab1 files using the sidebar panel</div>'
            '</div>',
            unsafe_allow_html=True,
        )
        c2.markdown(
            '<div class="landing-step">'
            '<div class="landing-step-num">02</div>'
            '<div class="landing-step-title">Run analysis</div>'
            '<div class="landing-step-desc">Dual-strand alignment, IUPAC variant calling, and HGVS annotation in seconds</div>'
            '</div>',
            unsafe_allow_html=True,
        )
        c3.markdown(
            '<div class="landing-step">'
            '<div class="landing-step-num">03</div>'
            '<div class="landing-step-title">Review results</div>'
            '<div class="landing-step-desc">Variant cards, interactive chromatograms, QC metrics, and a downloadable report</div>'
            '</div>',
            unsafe_allow_html=True,
        )
        return

    # ── Cache invalidation ────────────────────────────────────────────────────
    if mode == "Paired (Fwd + Rev)":
        files_key = hashlib.md5(fwd_file.getvalue() + rev_file.getvalue()).hexdigest()
    else:
        files_key = hashlib.md5(single_file.getvalue() + bytes([int(is_reverse)])).hexdigest()

    if st.session_state.get("_files_key") != files_key:
        st.session_state["_files_key"] = files_key
        st.session_state.pop("result", None)
        st.session_state.pop("dismissed_warnings", None)
        st.session_state.pop("rejected_variants", None)
    if st.session_state.get("_pipeline_v") != _PIPELINE_VERSION:
        st.session_state["_pipeline_v"] = _PIPELINE_VERSION
        st.session_state.pop("result", None)
        st.session_state.pop("dismissed_warnings", None)
        st.session_state.pop("rejected_variants", None)
        run_pipeline.clear()
        run_pipeline_single.clear()

    if run_btn:
        if mode == "Paired (Fwd + Rev)":
            result = run_pipeline(
                fwd_file.getvalue(), rev_file.getvalue(),
                fwd_file.name, rev_file.name,
                ref_path_str, min_phred, het_ratio,
            )
        else:
            result = run_pipeline_single(
                single_file.getvalue(), single_file.name,
                is_reverse, ref_path_str, min_phred, het_ratio,
            )
        st.session_state["result"] = result

    if "result" not in st.session_state:
        file_rows = ""
        if fwd_file:
            file_rows += (
                f'<div style="display:flex;align-items:baseline;gap:14px;margin-bottom:10px;">'
                f'<span style="font-size:0.65rem;font-weight:700;letter-spacing:1.5px;'
                f'text-transform:uppercase;color:#3B82F6;width:64px;flex-shrink:0;">Forward</span>'
                f'<span style="font-family:ui-monospace,monospace;font-size:1rem;color:#CBD5E1;">'
                f'{fwd_file.name}</span></div>'
            )
        if rev_file:
            file_rows += (
                f'<div style="display:flex;align-items:baseline;gap:14px;margin-bottom:10px;">'
                f'<span style="font-size:0.65rem;font-weight:700;letter-spacing:1.5px;'
                f'text-transform:uppercase;color:#A78BFA;width:64px;flex-shrink:0;">Reverse</span>'
                f'<span style="font-family:ui-monospace,monospace;font-size:1rem;color:#CBD5E1;">'
                f'{rev_file.name}</span></div>'
            )
        if single_file:
            label = "Reverse" if is_reverse else "Forward"
            colour = "#A78BFA" if is_reverse else "#3B82F6"
            file_rows += (
                f'<div style="display:flex;align-items:baseline;gap:14px;margin-bottom:10px;">'
                f'<span style="font-size:0.65rem;font-weight:700;letter-spacing:1.5px;'
                f'text-transform:uppercase;color:{colour};width:64px;flex-shrink:0;">{label}</span>'
                f'<span style="font-family:ui-monospace,monospace;font-size:1rem;color:#CBD5E1;">'
                f'{single_file.name}</span></div>'
            )
        st.markdown(
            f'<div style="max-width:560px;margin:100px auto;text-align:left;">'
            f'<div style="font-size:0.65rem;font-weight:700;letter-spacing:2.5px;'
            f'text-transform:uppercase;color:#3B82F6;margin-bottom:24px;">Files loaded</div>'
            f'{file_rows}'
            f'<div style="margin-top:28px;padding-top:24px;'
            f'border-top:1px solid rgba(255,255,255,0.07);'
            f'font-size:0.9rem;color:#475569;">'
            f'Click <strong style="color:#94A3B8;font-weight:600;">Run Analysis</strong> '
            f'in the sidebar to begin.</div>'
            f'</div>',
            unsafe_allow_html=True,
        )
        return

    result = st.session_state["result"]

    if result.get("error"):
        st.markdown(
            f'<div class="alert alert-warn"><strong>Analysis failed</strong> — '
            f'{result["error"]}<br><small style="opacity:0.7">Check the application '
            f'logs for details.</small></div>',
            unsafe_allow_html=True,
        )
        return

    report    = result["report"]
    fwd_trace = result["fwd_trace"]
    rev_rc    = result["rev_rc"]
    warnings  = result.get("warnings", [])

    # Shared filtered variant list — UTR background calls without a known name
    # are excluded everywhere (Summary, Variant Review, QC map) for consistency.
    all_sorted = sorted(report.variants, key=lambda x: x.ref_pos_genomic)
    filtered_variants = [
        v for v in all_sorted
        if v.known_variant_name or v.region not in ("5UTR", "3UTR")
    ]
    utr_noise_count = len(all_sorted) - len(filtered_variants)

    tab_sum, tab_chrom, tab_qc, tab_trace, tab_report = st.tabs(
        ["Summary", "Variant Review", "Quality Control", "Chromatogram", "Report"]
    )

    # ── Summary ───────────────────────────────────────────────────────────────
    with tab_sum:
        rejected_v: set = st.session_state.setdefault("rejected_variants", set())

        # Dismissible warnings — inside the tab to avoid tab-reset on rerun
        dismissed_w: set = st.session_state.setdefault("dismissed_warnings", set())
        for i, w in enumerate(warnings):
            if i in dismissed_w:
                continue
            col_msg, col_x = st.columns([0.97, 0.03])
            col_msg.markdown(f'<div class="alert alert-warn">{w}</div>',
                             unsafe_allow_html=True)
            col_x.button("✕", key=f"dismiss_w_{i}", help="Dismiss warning",
                         on_click=dismissed_w.add, args=(i,))

        display_variants = [v for v in filtered_variants if v.hgvs_c not in rejected_v]
        n = len(display_variants)

        if n == 0:
            msg = "No variants detected — sequence is wild type at all covered positions."
            if utr_noise_count:
                msg = (
                    f"No pathogenic variants detected. "
                    f"{utr_noise_count} background UTR haplotype variant(s) excluded from display."
                )
            st.markdown(f'<div class="alert alert-pass">{msg}</div>',
                        unsafe_allow_html=True)
        else:
            count_txt = f"{n} variant{'s' if n > 1 else ''} detected"
            st.markdown(
                f'<div class="result-header">'
                f'<span class="result-header-count">{count_txt}</span>'
                f'<span class="result-header-sample"> &nbsp;—&nbsp; {report.sample_id}</span>'
                f'</div>',
                unsafe_allow_html=True,
            )

            # Pre-compute significance for all display variants
            disp_sigs = {v.hgvs_c: _sig_for_variant(v) for v in display_variants}

            for v in display_variants:
                sig = disp_sigs[v.hgvs_c]
                col_card, col_x = st.columns([0.97, 0.03])
                col_card.markdown(_variant_card_html(v, sig), unsafe_allow_html=True)
                col_x.button(
                    "✕", key=f"reject_{v.ref_pos_genomic}",
                    help="Reject this variant call",
                    on_click=rejected_v.add, args=(v.hgvs_c,),
                )

        # Rejected calls — collapsible recovery
        if rejected_v:
            rejected_display = [v for v in all_sorted if v.hgvs_c in rejected_v]
            with st.expander(f"Rejected calls ({len(rejected_display)})  —  click to restore"):
                for v in rejected_display:
                    col_r, col_undo = st.columns([0.88, 0.12])
                    desc = v.hgvs_c
                    if v.known_variant_name:
                        desc += f"  ·  {v.known_variant_name}"
                    desc += f"  ·  {v.zygosity.value}  ·  {v.region}"
                    col_r.markdown(
                        f'<div class="rejected-row"><span class="rejected-hgvs">{desc}</span></div>',
                        unsafe_allow_html=True,
                    )
                    col_undo.button("Restore", key=f"restore_{v.ref_pos_genomic}",
                                    on_click=rejected_v.discard, args=(v.hgvs_c,))


    # ── Quality Control ───────────────────────────────────────────────────────
    with tab_qc:
        qc     = report.qc_metrics
        fwd_al = result["fwd_aligned"]
        rev_al = result["rev_aligned"]

        show_var_ticks = st.checkbox(
            "Show variant positions on gene map", value=True, key="cov_show_variants"
        )
        fig_cov = plot_coverage_map(
            fwd_al.reference_start if fwd_al else None,
            fwd_al.reference_end   if fwd_al else None,
            rev_al.reference_start if rev_al else None,
            rev_al.reference_end   if rev_al else None,
            variants=filtered_variants if show_var_ticks else None,
        )
        st.plotly_chart(fig_cov, use_container_width=True)

        st.markdown('<hr class="section-rule">', unsafe_allow_html=True)
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Quality (Fwd)", qc.get("mean_phred_fwd", "—"), help="Mean Phred score, forward read")
        c2.metric("Quality (Rev)", qc.get("mean_phred_rev", "—"), help="Mean Phred score, reverse read")
        c3.metric("CDS Coverage",  qc.get("cds_coverage_pct", "—"), help="Fraction of coding exons covered")
        c4.metric("Usable bp (Fwd)", f"{qc.get('usable_length_fwd', '—')} bp", help="Bases after quality trimming")
        c5.metric("Usable bp (Rev)", f"{qc.get('usable_length_rev', '—')} bp", help="Bases after quality trimming")

        import plotly.graph_objects as go
        fwd_phred = result.get("fwd_phred_raw")
        rev_phred = result.get("rev_phred_raw")
        if fwd_phred or rev_phred:
            st.markdown('<hr class="section-rule">', unsafe_allow_html=True)
            fig_q = go.Figure()
            if fwd_phred:
                fig_q.add_trace(go.Histogram(x=fwd_phred, name="Forward",
                                             opacity=0.75, nbinsx=40,
                                             marker_color="#3B82F6"))
            if rev_phred:
                fig_q.add_trace(go.Histogram(x=rev_phred, name="Reverse",
                                             opacity=0.75, nbinsx=40,
                                             marker_color="#EF4444"))
            fig_q.update_layout(
                barmode="overlay", height=260,
                title=dict(text="Phred score distribution", font=dict(size=13, color="#94A3B8")),
                xaxis=dict(title="Phred score", color="#64748B",
                           gridcolor="rgba(100,116,139,0.2)", zeroline=False),
                yaxis=dict(title="Base count", color="#64748B",
                           gridcolor="rgba(100,116,139,0.2)", zeroline=False),
                plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
                font=dict(color="#94A3B8"),
                legend=dict(orientation="h", x=1, xanchor="right", y=1,
                            bgcolor="rgba(0,0,0,0)"),
                margin=dict(l=10, r=10, t=40, b=10),
            )
            st.plotly_chart(fig_q, use_container_width=True)

    # ── Variant Review (chromatogram) ─────────────────────────────────────────
    with tab_chrom:
        if not filtered_variants:
            st.markdown('<div class="alert alert-info">No variants detected — nothing to review.</div>',
                        unsafe_allow_html=True)
        else:
            sorted_vars = filtered_variants
            pill_labels = [v.hgvs_c for v in sorted_vars]

            # Pre-select pathogenic variants; fall back to all if none are pathogenic
            default_sel = [v.hgvs_c for v in sorted_vars if _sig_for_variant(v) == "Pathogenic"]
            if not default_sel:
                default_sel = pill_labels[:]

            st.markdown('<span class="sidebar-label">Select variants to review</span>',
                        unsafe_allow_html=True)
            sel_labels = st.pills(
                "variants", pill_labels,
                selection_mode="multi",
                default=default_sel,
                key="vr_pill_sel",
                label_visibility="collapsed",
            )
            st.markdown('<hr class="section-rule">', unsafe_allow_html=True)

            if not sel_labels:
                st.markdown(
                    '<div class="alert alert-info">Select one or more variants above to view chromatograms.</div>',
                    unsafe_allow_html=True)
            else:
                selected_vars = [v for v in sorted_vars if v.hgvs_c in sel_labels]
                fwd_al = result["fwd_aligned"]
                rev_al = result["rev_aligned"]

                _SHOW_LIMIT = 5
                visible = selected_vars[:_SHOW_LIMIT]
                hidden  = selected_vars[_SHOW_LIMIT:]

                for v in visible:
                    _render_variant_block(v, _sig_for_variant(v),
                                          fwd_trace, rev_rc, fwd_al, rev_al)

                if hidden:
                    with st.expander(
                        f"Show {len(hidden)} more variant{'s' if len(hidden) > 1 else ''} "
                        f"— likely artefacts or UTR haplotype calls",
                        expanded=False,
                    ):
                        for v in hidden:
                            _render_variant_block(v, _sig_for_variant(v),
                                                  fwd_trace, rev_rc, fwd_al, rev_al)

    # ── Chromatogram viewer ───────────────────────────────────────────────────
    with tab_trace:
        import plotly.graph_objects as go

        _TV_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#999999", "T": "#d62728"}

        def _render_full_trace(trace: TraceData) -> None:
            channels = trace.channels
            n_scans  = max((len(ch) for ch in channels.values() if ch), default=1)
            stride   = max(1, n_scans // 3000)

            all_sampled: list[int] = []
            for ch in channels.values():
                if ch:
                    all_sampled.extend(ch[i] for i in range(0, len(ch), stride) if i < len(ch))
            if all_sampled:
                all_sampled.sort()
                y_cap = all_sampled[min(int(len(all_sampled) * 0.99), len(all_sampled) - 1)]
            else:
                y_cap = 1000
            y_range = [0, max(y_cap * 1.08, 100)]

            fig = go.Figure()
            for base, color in _TV_COLORS.items():
                ch = channels.get(base, [])
                if not ch:
                    continue
                xs = list(range(0, n_scans, stride))
                ys = [ch[i] for i in xs if i < len(ch)]
                fig.add_trace(go.Scattergl(
                    x=xs[:len(ys)], y=ys, mode="lines", name=base,
                    line=dict(color=color, width=1),
                ))

            peaks      = trace.peak_positions
            seq        = trace.sequence
            label_step = max(1, len(seq) // 50)
            for i in range(0, len(seq), label_step):
                if i >= len(peaks):
                    break
                b = seq[i]
                fig.add_annotation(
                    x=peaks[i], y=-0.06, xref="x", yref="paper",
                    text=f"<b>{b}</b>", showarrow=False,
                    font=dict(color=_TV_COLORS.get(b.upper(), "#888888"), size=9),
                )

            _AXIS = dict(color="#64748B", gridcolor="rgba(100,116,139,0.18)",
                         zerolinecolor="rgba(100,116,139,0.25)")
            fig.update_layout(
                height=360,
                margin=dict(l=10, r=10, t=10, b=50),
                xaxis=dict(
                    title="Scan position",
                    rangeslider=dict(visible=True, thickness=0.07,
                                     bgcolor="rgba(30,30,40,0.95)",
                                     bordercolor="rgba(100,116,139,0.35)",
                                     borderwidth=1),
                    **_AXIS,
                ),
                yaxis=dict(title="Intensity", range=y_range, fixedrange=False, **_AXIS),
                font=dict(color="#94A3B8"),
                legend=dict(orientation="h", x=1, xanchor="right", y=1, yanchor="bottom",
                            bgcolor="rgba(0,0,0,0)"),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
            )
            st.plotly_chart(fig, use_container_width=True)

        if fwd_trace is None and rev_rc is None:
            st.markdown('<div class="alert alert-info">No trace data available.</div>',
                        unsafe_allow_html=True)
        else:
            if fwd_trace is not None:
                st.markdown(
                    '<div class="trace-section-header">'
                    '<span class="trace-section-dot dot-fwd"></span>'
                    'Forward strand'
                    '</div>',
                    unsafe_allow_html=True,
                )
                _render_full_trace(fwd_trace)
            else:
                st.markdown('<div class="alert alert-info">Forward trace not available.</div>',
                            unsafe_allow_html=True)

            if rev_rc is not None:
                st.markdown(
                    '<div class="trace-section-header" style="margin-top:18px">'
                    '<span class="trace-section-dot dot-rev"></span>'
                    'Reverse complement'
                    '</div>',
                    unsafe_allow_html=True,
                )
                _render_full_trace(rev_rc)
            else:
                st.markdown('<div class="alert alert-info">Reverse complement trace not available.</div>',
                            unsafe_allow_html=True)

    # ── Report ────────────────────────────────────────────────────────────────
    with tab_report:
        st.markdown(result["markdown"])
        st.download_button(
            "Download report (.md)",
            data=result["markdown"],
            file_name=f"{report.sample_id}_HBB_report.md",
            mime="text/markdown",
        )


if __name__ == "__main__":
    main()
