"""CLI entry point for the HBB Sanger analysis pipeline.

Usage:
    hbb-analyze run fwd.ab1 rev.ab1 --reference reference/HBB_reference.fasta --out report.md
    hbb-analyze run fwd.ab1 rev.ab1 --json out.json
    hbb-analyze validate-reference reference/HBB_reference.fasta
"""

from __future__ import annotations

import json
import logging
import sys
import tempfile
from pathlib import Path

import typer

# Ensure hbb_pipeline importable when run from hbb_sanger/
sys.path.insert(0, str(Path(__file__).parent))

app = typer.Typer(
    name="hbb-analyze",
    help="HBB Sanger sequencing analysis — beta-thalassemia & sickle cell variant detection.",
    add_completion=False,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

_DEFAULT_REF = Path(__file__).parent / "reference" / "HBB_reference.fasta"


@app.command()
def run(
    fwd: Path = typer.Argument(..., help="Forward .ab1 trace file"),
    rev: Path = typer.Argument(..., help="Reverse .ab1 trace file"),
    reference: Path = typer.Option(_DEFAULT_REF, "--reference", "-r", help="Reference FASTA"),
    out: Path | None = typer.Option(None, "--out", "-o", help="Output Markdown report path"),
    json_out: Path | None = typer.Option(None, "--json", help="Output JSON report path"),
    min_phred: int = typer.Option(20, "--min-phred", help="Minimum Phred for quality trimming"),
) -> None:
    """Run the full HBB analysis pipeline on a forward + reverse trace pair."""
    from hbb_pipeline.alignment import build_consensus
    from hbb_pipeline.coordinates import CoordinateTranslator
    from hbb_pipeline.parsing import parse_abi
    from hbb_pipeline.pipeline import process_trace
    from hbb_pipeline.qc import evaluate_trace_artifacts
    from hbb_pipeline.reference import HBBReference, ReferenceValidationError
    from hbb_pipeline.reporting import generate_report, render_markdown_report
    from hbb_pipeline.variants import call_variants_from_alignment

    # Validate reference
    typer.echo(f"Loading reference: {reference}")
    try:
        ref = HBBReference(reference)
    except ReferenceValidationError as exc:
        typer.echo(f"ERROR: Reference validation failed: {exc}", err=True)
        raise typer.Exit(1)
    typer.echo(f"  ✓ Reference OK — {ref.length} bp, CDS {len(ref.cds)} bp")

    translator = CoordinateTranslator(ref)

    # Parse and process traces
    typer.echo(f"Parsing: {fwd.name} (fwd), {rev.name} (rev)")
    from hbb_pipeline.parsing import InvalidTraceFileError

    try:
        fwd_raw = parse_abi(fwd)
        rev_raw = parse_abi(rev)
    except InvalidTraceFileError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(1)

    fwd_warnings = evaluate_trace_artifacts(fwd_raw)
    for w in fwd_warnings:
        typer.echo(f"  [Fwd Warning] {w}", err=True)
    rev_warnings = evaluate_trace_artifacts(rev_raw)
    for w in rev_warnings:
        typer.echo(f"  [Rev Warning] {w}", err=True)
    all_warnings = [f"[Forward] {w}" for w in fwd_warnings] + [f"[Reverse] {w}" for w in rev_warnings]

    fwd_iupac, fwd_aligned, fs, fe, fwd_err = process_trace(fwd_raw, ref, is_reverse=False, min_phred=min_phred)
    rev_iupac, rev_aligned, rs, re, rev_err = process_trace(rev_raw, ref, is_reverse=True,  min_phred=min_phred)

    if fwd_err or rev_err:
        typer.echo(f"ERROR: Trace processing failed — fwd={fwd_err} rev={rev_err}", err=True)
        raise typer.Exit(1)

    typer.echo(f"  Fwd trimmed: {fe - fs} bp  |  Rev trimmed: {re - rs} bp")

    typer.echo("Aligning and calling variants…")
    consensus, cons_quals = build_consensus(fwd_aligned, rev_aligned, ref)
    merged = call_variants_from_alignment(
        consensus, cons_quals, ref, translator,
        fwd_iupac, rev_iupac, fwd_aligned, rev_aligned,
    )
    typer.echo(f"  {len(merged)} variant(s) called")

    # Build QC
    qc = {
        "mean_phred_fwd": round(sum(fwd_raw.phred_scores) / max(len(fwd_raw.phred_scores), 1), 1),
        "mean_phred_rev": round(sum(rev_raw.phred_scores) / max(len(rev_raw.phred_scores), 1), 1),
        "usable_length_fwd": fe - fs,
        "usable_length_rev": re - rs,
        "analysis_warnings": all_warnings,
    }

    sample_id = fwd.stem
    report = generate_report(merged, qc, sample_id, consensus)

    # Output
    md = render_markdown_report(report)

    if out:
        out.write_text(md, encoding="utf-8")
        typer.echo(f"  Markdown report: {out}")
    else:
        typer.echo("\n" + md)

    if json_out:
        json_out.write_text(
            report.model_dump_json(indent=2),
            encoding="utf-8",
        )
        typer.echo(f"  JSON report: {json_out}")

    # Print summary of known variants
    known = [v for v in merged if v.known_variant_name]
    if known:
        typer.echo("\n=== Clinically Significant Findings ===")
        for v in known:
            typer.echo(
                f"  {v.hgvs_c}  {v.known_variant_name}  [{v.zygosity.value}]  "
                f"p.: {v.hgvs_p_hgvs or '—'} / {v.hgvs_p_legacy or '—'}"
            )

    flagged = [v for v in merged if v.requires_manual_review]
    if flagged:
        typer.echo(f"\n⚠  {len(flagged)} variant(s) require manual review.")


@app.command(name="validate-reference")
def validate_reference(
    reference: Path = typer.Argument(..., help="Path to the case-annotated reference FASTA"),
) -> None:
    """Validate the structure of a case-annotated HBB reference FASTA."""
    from hbb_pipeline.reference import HBBReference, ReferenceValidationError

    typer.echo(f"Validating: {reference}")
    try:
        ref = HBBReference(reference)
    except ReferenceValidationError as exc:
        typer.echo(f"FAIL: {exc}", err=True)
        raise typer.Exit(1)

    typer.echo(f"✓ PASS: {ref.length} bp")
    typer.echo(f"  5' flank:  {ref.utr5[1] - ref.utr5[0]} bp")
    typer.echo(f"  Exon 1:    {ref.exons[0][1] - ref.exons[0][0]} bp  (starts ATG: {ref.upper_seq[978:981]})")
    typer.echo(f"  Intron 1:  {ref.introns[0][1] - ref.introns[0][0]} bp")
    typer.echo(f"  Exon 2:    {ref.exons[1][1] - ref.exons[1][0]} bp")
    typer.echo(f"  Intron 2:  {ref.introns[1][1] - ref.introns[1][0]} bp")
    typer.echo(f"  Exon 3:    {ref.exons[2][1] - ref.exons[2][0]} bp  (ends TAA: {ref.upper_seq[2149:2152]})")
    typer.echo(f"  3' flank:  {ref.utr3[1] - ref.utr3[0]} bp")
    typer.echo(f"  CDS:       {len(ref.cds)} bp → 147 aa β-globin")


if __name__ == "__main__":
    app()
