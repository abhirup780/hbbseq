"""Plotly chromatogram visualizations (UI layer — no hbb_pipeline imports here except models).

All functions return plotly.graph_objects.Figure so they can be embedded in
Streamlit with st.plotly_chart() or exported to HTML.

Color scheme (colorblind-safe):
    A = #2ca02c (green)
    C = #1f77b4 (blue)
    G = #999999 (mid-grey, visible on both light and dark backgrounds)
    T = #d62728 (red)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import plotly.graph_objects as go

if TYPE_CHECKING:
    from hbb_pipeline.models import TraceData, Variant

_CHANNEL_COLORS = {
    "A": "#2ca02c",
    "C": "#1f77b4",
    "G": "#999999",
    "T": "#d62728",
}


def plot_chromatogram(
    trace: "TraceData",
    start: int,
    end: int,
    highlight_pos: int | None = None,
    trim_start: int | None = None,
    trim_end: int | None = None,
    max_points: int = 1500,
) -> go.Figure:
    """Plot a region of a Sanger chromatogram as four overlaid line traces.

    The x-axis is raw scan index (intensity sample number); base calls are
    annotated at their peak_positions.  A range slider is added for zooming.

    Args:
        trace: Parsed TraceData object.
        start: First *base index* (0-based) to display.
        end: Last base index (0-based, exclusive).
        highlight_pos: If given, draw a vertical dashed line at this base index.
        trim_start: If given, shade bases before this index as trimmed (green boundary).
        trim_end: If given, shade bases from this index onward as trimmed (red boundary).

    Returns:
        A plotly Figure with four channel traces and base annotations.

    Example:
        fig = plot_chromatogram(trace, 0, len(trace.sequence))
        fig.show()
    """
    # Determine scan range from peak positions of the requested base window
    peak_positions = trace.peak_positions
    clipped_start = max(0, min(start, len(peak_positions) - 1))
    clipped_end = max(clipped_start + 1, min(end, len(peak_positions)))

    scan_start = peak_positions[clipped_start]
    scan_end = peak_positions[clipped_end - 1]

    # Add a small margin (10% of window) for visual breathing room
    margin = max(10, (scan_end - scan_start) // 10)
    x_min = max(0, scan_start - margin)
    x_max = scan_end + margin

    fig = go.Figure()

    # Compute stride so the browser never receives more than max_points per channel
    raw_len = x_max - x_min + 1
    stride = max(1, raw_len // max_points)

    # Add one line per channel
    for base in ("A", "C", "G", "T"):
        channel = trace.channels.get(base, [])
        if not channel:
            continue
        x_slice = list(range(x_min, min(x_max + 1, len(channel)), stride))
        y_slice = channel[x_min: x_max + 1: stride]
        fig.add_trace(
            go.Scatter(
                x=x_slice,
                y=y_slice,
                mode="lines",
                name=base,
                line=dict(color=_CHANNEL_COLORS[base], width=1.5),
            )
        )

    # Base call annotations at peak positions
    for i in range(clipped_start, clipped_end):
        if i >= len(peak_positions) or i >= len(trace.sequence):
            break
        px = peak_positions[i]
        base = trace.sequence[i]
        color = _CHANNEL_COLORS.get(base.upper(), "#888888")
        fig.add_annotation(
            x=px,
            y=-0.08,
            xref="x",
            yref="paper",
            text=f"<b>{base}</b>",
            showarrow=False,
            font=dict(color=color, size=10),
        )

    # Trim region shading — grey out what gets discarded
    if trim_start is not None and 0 < trim_start < len(peak_positions):
        ts_scan = peak_positions[max(0, trim_start - 1)]
        fig.add_vrect(
            x0=x_min, x1=ts_scan,
            fillcolor="rgba(120,120,120,0.18)", line_width=0,
        )
        fig.add_vline(
            x=ts_scan, line_color="#2ca02c", line_width=2,
            annotation_text="trim start", annotation_position="top right",
            annotation_font_color="#2ca02c",
        )

    if trim_end is not None and 0 < trim_end <= len(peak_positions):
        te_scan = peak_positions[min(trim_end, len(peak_positions)) - 1]
        fig.add_vrect(
            x0=te_scan, x1=x_max,
            fillcolor="rgba(120,120,120,0.18)", line_width=0,
        )
        fig.add_vline(
            x=te_scan, line_color="#d62728", line_width=2,
            annotation_text="trim end", annotation_position="top left",
            annotation_font_color="#d62728",
        )

    # Highlight position (vertical dashed line)
    if highlight_pos is not None and 0 <= highlight_pos < len(peak_positions):
        hl_x = peak_positions[highlight_pos]
        fig.add_vline(
            x=hl_x,
            line_dash="dash",
            line_color="orange",
            line_width=2,
            annotation_text=f"pos {highlight_pos}",
            annotation_position="top right",
        )

    fig.update_layout(
        title=dict(text=f"Chromatogram — {trace.sample_name}", font=dict(size=14)),
        xaxis=dict(
            title="Scan",
            rangeslider=dict(visible=True),
            range=[x_min, x_max],
            color="#aaaaaa",
            gridcolor="rgba(120,120,120,0.2)",
        ),
        yaxis=dict(title="Intensity", color="#aaaaaa", gridcolor="rgba(120,120,120,0.2)"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
        margin=dict(b=60),
        height=400,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#aaaaaa"),
    )

    return fig


def _genomic_to_c(pos: int) -> str:
    """Convert 0-based genomic position to HGVS c. label (HBB-specific)."""
    if pos < 978:
        return f"c.{pos - 978}"           # e.g. c.-165
    if pos < 1070:                         # Exon 1
        return f"c.{pos - 977}"
    if pos < 1200:                         # IVS-I
        d_from = pos - 1070 + 1
        d_to   = 1200 - pos
        return f"c.92+{d_from}" if d_from <= d_to else f"c.93-{d_to}"
    if pos < 1423:                         # Exon 2
        return f"c.{pos - 1107}"
    if pos < 2023:                         # IVS-II
        d_from = pos - 1423 + 1
        d_to   = 2023 - pos
        return f"c.315+{d_from}" if d_from <= d_to else f"c.316-{d_to}"
    if pos < 2152:                         # Exon 3
        return f"c.{pos - 1707}"
    return f"c.*{pos - 2151}"             # 3' flank


def _short_variant_label(v: "Variant") -> str:  # type: ignore[name-defined]
    """Compact HGVS label for a variant tick."""
    label = v.hgvs_c
    return label if len(label) <= 12 else label[:11] + "…"


def plot_coverage_map(
    fwd_ref_start: int | None,
    fwd_ref_end: int | None,
    rev_ref_start: int | None,
    rev_ref_end: int | None,
    variants: "list[Variant] | None" = None,  # type: ignore[name-defined]
) -> go.Figure:
    """Gene-structure diagram with forward and/or reverse read coverage overlaid."""

    TOTAL     = 2750
    DISP_START = 700   # left edge of displayed region (captures all promoter variants)
    DISP_END   = 2480  # right edge (past Exon 3 end + primer buffer)
    Y_GENE    = 2.10   # centre of gene cartoon
    H_EXON    = 0.32   # half-height for exon blocks
    H_UTR     = 0.16   # half-height for UTR blocks
    Y_FWD     = 1.10   # centre of forward coverage bar
    Y_REV     = 0.38   # centre of reverse coverage bar
    H_COVER   = 0.12
    GENE_TOP  = Y_GENE + H_EXON

    # ── Variant tier layout ───────────────────────────────────────────────────
    # Greedy lane assignment: each variant goes into the lowest tier whose last
    # occupant is ≥ MIN_SEP genomic bases away.  Vertical labels (textangle=-90)
    # take no horizontal space so the only collision axis is Y.
    MIN_SEP   = 140
    TIER_STEP = 0.52   # Y gap between consecutive tiers
    STEM_BASE = GENE_TOP + 0.18  # bottom of first tier's stem top

    show_variants = variants is not None and len(variants) > 0
    sorted_vars: list = []
    var_tiers:   list[int] = []
    n_tiers = 0

    if show_variants:
        sorted_vars = sorted(variants, key=lambda v: v.ref_pos_genomic)
        tier_last: list[float] = []
        for v in sorted_vars:
            gpos = v.ref_pos_genomic
            placed = False
            for i, last in enumerate(tier_last):
                if gpos - last >= MIN_SEP:
                    tier_last[i] = gpos
                    var_tiers.append(i)
                    placed = True
                    break
            if not placed:
                var_tiers.append(len(tier_last))
                tier_last.append(float(gpos))
        n_tiers = max(var_tiers) + 1

    # Label space above the highest tier stem (~1.1 Y units for vertical text)
    Y_TOP = STEM_BASE + n_tiers * TIER_STEP + 1.15 if show_variants else GENE_TOP + 0.5

    shapes: list[dict] = []
    annotations: list[dict] = []

    # ── Gene cartoon ─────────────────────────────────────────────────────────
    _SEGMENTS = [
        ("5′ flank",  0,    978,  "#1E3A5F", H_UTR,  False),
        ("Exon 1",   978,  1070, "#2563EB", H_EXON, True),
        ("IVS-I",   1070, 1200, None,      0,      False),
        ("Exon 2",  1200,  1423, "#2563EB", H_EXON, True),
        ("IVS-II",  1423, 2023, None,      0,      False),
        ("Exon 3",  2023,  2152, "#2563EB", H_EXON, True),
        ("3′ flank", 2152, 2750, "#1E3A5F", H_UTR,  False),
    ]

    for label, s, e, color, h, is_exon in _SEGMENTS:
        if color is None:
            # Intron — draw as a thin horizontal line (classic gene diagram style)
            shapes.append(dict(
                type="line", x0=s, x1=e, y0=Y_GENE, y1=Y_GENE,
                line=dict(color="#4B5563", width=2.5),
            ))
            annotations.append(dict(
                x=(s + e) / 2, y=Y_GENE + 0.04,
                text=label, showarrow=False,
                font=dict(size=10, color="#4B5563"),
                xanchor="center", yanchor="bottom",
                bgcolor="rgba(10,15,26,0.0)",
            ))
        else:
            shapes.append(dict(
                type="rect", x0=s, x1=e,
                y0=Y_GENE - h, y1=Y_GENE + h,
                fillcolor=color, line_width=0,
            ))
            vis_s = max(s, DISP_START)
            vis_e = min(e, DISP_END)
            if vis_e - vis_s > 50:
                annotations.append(dict(
                    x=(vis_s + vis_e) / 2, y=Y_GENE,
                    text=label, showarrow=False,
                    font=dict(size=10, color="white" if is_exon else "#93C5FD"),
                    xanchor="center", yanchor="middle",
                ))

    # ── Exon/intron boundary guide lines ─────────────────────────────────────
    for pos in (978, 1070, 1200, 1423, 2023, 2152):
        shapes.append(dict(
            type="line", x0=pos, x1=pos, y0=0, y1=GENE_TOP,
            line=dict(color="rgba(255,255,255,0.06)", width=1, dash="dot"),
        ))

    # ── Coverage bars ─────────────────────────────────────────────────────────
    def _coverage_bar(y: float, start: int, end: int, color: str, label: str) -> None:
        shapes.append(dict(
            type="rect", x0=start, x1=end,
            y0=y - H_COVER, y1=y + H_COVER,
            fillcolor=color, opacity=0.82, line_width=0,
        ))
        for xpos in (start, end):
            annotations.append(dict(
                x=xpos, y=y - H_COVER - 0.06,
                text=_genomic_to_c(xpos), showarrow=False,
                font=dict(size=10, color=color),
                xanchor="center", yanchor="top",
            ))
        annotations.append(dict(
            xref="paper", x=-0.02, y=y,
            text=f"<b>{label}</b>", showarrow=False,
            font=dict(size=12, color=color),
            xanchor="right", yanchor="middle",
        ))

    def _missing_bar(y: float, label: str) -> None:
        annotations.append(dict(
            xref="paper", x=-0.02, y=y,
            text=f"<b>{label}</b>", showarrow=False,
            font=dict(size=12, color="#4B5563"),
            xanchor="right", yanchor="middle",
        ))
        annotations.append(dict(
            x=(DISP_START + DISP_END) / 2, y=y,
            text="not available", showarrow=False,
            font=dict(size=12, color="#4B5563"),
            xanchor="center", yanchor="middle",
        ))

    if fwd_ref_start is not None and fwd_ref_end is not None:
        _coverage_bar(Y_FWD, fwd_ref_start, fwd_ref_end, "#22C55E", "Fwd")
    else:
        _missing_bar(Y_FWD, "Fwd")

    if rev_ref_start is not None and rev_ref_end is not None:
        _coverage_bar(Y_REV, rev_ref_start, rev_ref_end, "#F87171", "Rev")
    else:
        _missing_bar(Y_REV, "Rev")

    # ── Variant lollipops ─────────────────────────────────────────────────────
    dot_xs: list[float] = []
    dot_ys: list[float] = []
    dot_colors: list[str] = []

    for v, tier in zip(sorted_vars, var_tiers):
        gpos     = v.ref_pos_genomic
        is_path  = getattr(v, "known_variant_name", None) is not None
        color    = "#EF4444" if is_path else "#64748B"
        stem_top = STEM_BASE + tier * TIER_STEP

        shapes.append(dict(
            type="line", x0=gpos, x1=gpos,
            y0=GENE_TOP, y1=stem_top,
            line=dict(color=color, width=1.5),
        ))
        dot_xs.append(gpos)
        dot_ys.append(stem_top)
        dot_colors.append(color)

        short = _short_variant_label(v)
        annotations.append(dict(
            x=gpos, y=stem_top + 0.06,
            text=f"<b>{short}</b>" if is_path else short,
            showarrow=False,
            font=dict(size=11, color=color),
            xanchor="center", yanchor="bottom",
            textangle=-90,
        ))

    # ── Figure ────────────────────────────────────────────────────────────────
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[DISP_START, DISP_END], y=[0, Y_TOP],
        mode="markers", marker=dict(opacity=0),
        hoverinfo="skip", showlegend=False,
    ))
    if dot_xs:
        fig.add_trace(go.Scatter(
            x=dot_xs, y=dot_ys,
            mode="markers",
            marker=dict(size=6, color=dot_colors, line_width=0),
            hoverinfo="skip", showlegend=False,
        ))

    # Gene + transcript label — bottom-right, anchors HGVS coordinates to RefSeq
    annotations.append(dict(
        xref="paper", yref="paper", x=1.0, y=0.0,
        text="<i>HBB</i> · NM_000518.5", showarrow=False,
        font=dict(size=9, color="#334155"),
        xanchor="right", yanchor="top",
    ))

    height = 260 + n_tiers * 45
    fig.update_layout(
        shapes=shapes,
        annotations=annotations,
        height=height,
        margin=dict(l=55, r=20, t=8, b=30),
        xaxis=dict(
            range=[DISP_START, DISP_END],
            showticklabels=False,
            showgrid=False, zeroline=False,
        ),
        yaxis=dict(range=[0, Y_TOP], showticklabels=False,
                   showgrid=False, zeroline=False),
        font=dict(color="#94A3B8"),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )
    return fig
