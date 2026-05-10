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
    return f"c.*{pos - 2151}"             # 3'UTR


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
    """Gene-structure diagram with forward and/or reverse read coverage overlaid.

    Pass None for either strand's coordinates to omit that bar (single-strand mode
    or when a trace failed to align).  Pass a list of Variant objects to overlay
    tick marks above the gene cartoon; pass None to hide them.

    Draws the HBB locus (5'UTR → Exon1 → IVS-I → Exon2 → IVS-II → Exon3 → 3'UTR)
    as a classic gene cartoon, then shows the trimmed+aligned span of each read
    as a coloured bar underneath.

    Coordinates are 0-based genomic positions matching reference.py _EXPECTED.
    """
    # ── Gene structure (hardcoded to match reference.py _EXPECTED) ─────────────
    # (label, start, end, fill_color, bar_height_fraction)
    _SEGMENTS = [
        ("5'UTR",  0,    978,  "#aec7e8", 0.30),
        ("Exon 1", 978,  1070, "#1f77b4", 0.60),
        ("IVS-I",  1070, 1200, "#c7c7c7", 0.12),
        ("Exon 2", 1200, 1423, "#1f77b4", 0.60),
        ("IVS-II", 1423, 2023, "#c7c7c7", 0.12),
        ("Exon 3", 2023, 2152, "#1f77b4", 0.60),
        ("3'UTR",  2152, 2750, "#aec7e8", 0.30),
    ]

    TOTAL   = 2750
    show_variants = variants is not None and len(variants) > 0
    Y_TOP   = 4.60 if show_variants else 3.20   # top of plot area
    Y_GENE  = 2.55   # centre of gene row
    Y_FWD   = 1.55   # centre of forward coverage row
    Y_REV   = 0.55   # centre of reverse coverage row
    H_GENE  = 0.40   # half-height for full exon blocks
    H_COVER = 0.22   # half-height for coverage bars
    GENE_TOP = Y_GENE + H_GENE  # top edge of gene cartoon

    F_SEGMENT = 11   # gene segment label size
    F_EDGE    = 11   # coverage bar edge c. label size
    F_ROW     = 12   # Fwd/Rev row label size
    F_VAR     = 12   # variant tick label size

    shapes: list[dict] = []
    annotations: list[dict] = []

    # Gene cartoon
    for label, s, e, color, frac in _SEGMENTS:
        h = H_GENE * frac / 0.60  # scale to fraction
        shapes.append(dict(
            type="rect", x0=s, x1=e,
            y0=Y_GENE - h, y1=Y_GENE + h,
            fillcolor=color, line_width=0,
        ))
        # Only annotate segments wide enough to fit text
        if (e - s) > 60:
            text_color = "white" if "Exon" in label else "#444444"
            annotations.append(dict(
                x=(s + e) / 2, y=Y_GENE,
                text=label, showarrow=False,
                font=dict(size=F_SEGMENT, color=text_color),
                xanchor="center", yanchor="middle",
            ))

    def _coverage_bar(y: float, start: int, end: int, color: str, label: str) -> None:
        shapes.append(dict(
            type="rect", x0=start, x1=end,
            y0=y - H_COVER, y1=y + H_COVER,
            fillcolor=color, opacity=0.75, line_width=0,
        ))
        # c. coordinate labels at bar edges
        for xpos, anchor in ((start, "center"), (end, "center")):
            annotations.append(dict(
                x=xpos, y=y - H_COVER - 0.08,
                text=_genomic_to_c(xpos), showarrow=False,
                font=dict(size=F_EDGE, color=color),
                xanchor=anchor, yanchor="top",
            ))
        # Row label on the left
        annotations.append(dict(
            x=-30, y=y,
            text=f"<b>{label}</b>", showarrow=False,
            font=dict(size=F_ROW, color=color),
            xanchor="right", yanchor="middle",
        ))

    if fwd_ref_start is not None and fwd_ref_end is not None:
        _coverage_bar(Y_FWD, fwd_ref_start, fwd_ref_end, "#2ca02c", "Fwd")
    else:
        annotations.append(dict(
            x=-30, y=Y_FWD, text="<b>Fwd</b>", showarrow=False,
            font=dict(size=F_ROW, color="#888888"), xanchor="right", yanchor="middle",
        ))
        annotations.append(dict(
            x=TOTAL / 2, y=Y_FWD, text="not available", showarrow=False,
            font=dict(size=F_ROW, color="#888888"), xanchor="center", yanchor="middle",
        ))

    if rev_ref_start is not None and rev_ref_end is not None:
        _coverage_bar(Y_REV, rev_ref_start, rev_ref_end, "#d62728", "Rev")
    else:
        annotations.append(dict(
            x=-30, y=Y_REV, text="<b>Rev</b>", showarrow=False,
            font=dict(size=F_ROW, color="#888888"), xanchor="right", yanchor="middle",
        ))
        annotations.append(dict(
            x=TOTAL / 2, y=Y_REV, text="not available", showarrow=False,
            font=dict(size=F_ROW, color="#888888"), xanchor="center", yanchor="middle",
        ))

    # Vertical guide lines at exon/intron boundaries
    for pos in (978, 1070, 1200, 1423, 2023, 2152):
        shapes.append(dict(
            type="line", x0=pos, x1=pos, y0=0, y1=GENE_TOP,
            line=dict(color="#dddddd", width=1, dash="dot"),
        ))

    # ── Variant tick marks ────────────────────────────────────────────────────
    dot_xs: list[float] = []
    dot_ys: list[float] = []
    dot_colors: list[str] = []

    if show_variants:
        # Two label rows to avoid overlap: alternate when consecutive variants
        # are within MIN_SEP genomic positions.
        MIN_SEP     = 180
        Y_TICK_LOW  = GENE_TOP + 0.20   # tick tip — row A
        Y_TICK_HIGH = GENE_TOP + 0.70   # tick tip — row B (staggered)
        Y_LBL_LOW   = Y_TICK_LOW  + 0.12
        Y_LBL_HIGH  = Y_TICK_HIGH + 0.12

        sorted_vars = sorted(variants, key=lambda v: v.ref_pos_genomic)
        last_pos = -9999
        use_high = False

        for v in sorted_vars:
            gpos    = v.ref_pos_genomic
            is_path = getattr(v, "known_variant_name", None) is not None
            # Determine row
            if gpos - last_pos < MIN_SEP:
                use_high = not use_high
            else:
                use_high = False
            last_pos = gpos

            tick_y = Y_TICK_HIGH if use_high else Y_TICK_LOW
            lbl_y  = Y_LBL_HIGH  if use_high else Y_LBL_LOW
            color  = "#ff6b6b" if is_path else "#bbbbbb"

            # Tick stem from gene top to tick_y
            shapes.append(dict(
                type="line", x0=gpos, x1=gpos,
                y0=GENE_TOP, y1=tick_y,
                line=dict(color=color, width=2),
            ))
            # Collect dot position — rendered as scatter markers (pixel-space, stays round)
            dot_xs.append(gpos)
            dot_ys.append(tick_y)
            dot_colors.append(color)

            # Label
            short = _short_variant_label(v)
            annotations.append(dict(
                x=gpos, y=lbl_y,
                text=f"<b>{short}</b>" if is_path else short,
                showarrow=False,
                font=dict(size=F_VAR, color=color),
                xanchor="center", yanchor="bottom",
                textangle=-40,
            ))

    fig = go.Figure()
    # Invisible scatter to set axis range (shapes don't drive axes)
    fig.add_trace(go.Scatter(x=[0, TOTAL], y=[0, Y_TOP],
                             mode="markers", marker=dict(opacity=0)))
    # Lollipop heads — pixel-space markers stay perfectly round regardless of axis scaling
    if dot_xs:
        fig.add_trace(go.Scatter(
            x=dot_xs, y=dot_ys,
            mode="markers",
            marker=dict(size=9, color=dot_colors, line_width=0),
            hoverinfo="skip",
            showlegend=False,
        ))

    fig.update_layout(
        shapes=shapes,
        annotations=annotations,
        height=320 if show_variants else 260,
        margin=dict(l=55, r=20, t=35, b=20),
        xaxis=dict(
            range=[-80, TOTAL + 50],
            tickvals=[0, 978, 1070, 1200, 1423, 2023, 2152, 2750],
            ticktext=["0", "978\n(E1)", "1070\n(I1)", "1200\n(E2)",
                      "1423\n(I2)", "2023\n(E3)", "2152\n(3')", "2750"],
            tickfont=dict(size=10),
            showgrid=False,
            zeroline=False,
            color="rgba(150,150,150,1)",
        ),
        yaxis=dict(range=[0, Y_TOP], showticklabels=False,
                   showgrid=False, zeroline=False),
        font=dict(color="#aaaaaa"),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
        title=dict(text="HBB gene coverage", font=dict(size=13)),
    )
    return fig
