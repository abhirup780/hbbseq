"""Basic tests for reporting.py — generate_report and render_markdown_report."""
from datetime import datetime

from hbb_pipeline.models import ClinicalReport, Variant, Zygosity
from hbb_pipeline.reporting import generate_report, render_markdown_report


def _make_variant(**kwargs) -> Variant:
    defaults = dict(
        ref_pos_genomic=1006,
        ref_allele="A",
        alt_allele="T",
        zygosity=Zygosity.HET,
        variant_type="SNV",
        hgvs_c="c.20A>T",
        hgvs_p_hgvs="p.Glu7Val",
        hgvs_p_legacy="p.Glu6Val",
        region="exon1",
        phred_support=35,
        secondary_peak_ratio_fwd=0.45,
        secondary_peak_ratio_rev=0.42,
        called_by=["alignment"],
        known_variant_name="HbS",
        requires_manual_review=False,
    )
    defaults.update(kwargs)
    return Variant(**defaults)


def test_generate_report_fields() -> None:
    v = _make_variant()
    report = generate_report([v], {"mean_phred_fwd": 32.1}, "SAMPLE01", "ACGT")
    assert report.sample_id == "SAMPLE01"
    assert len(report.variants) == 1
    assert isinstance(report.analysis_timestamp, datetime)
    assert report.qc_metrics["mean_phred_fwd"] == 32.1


def test_generate_report_empty_variants() -> None:
    report = generate_report([], {}, "EMPTY", "")
    assert report.variants == []


def test_render_contains_sample_id() -> None:
    report = generate_report([], {}, "MYSAMPLE", "")
    md = render_markdown_report(report)
    assert "MYSAMPLE" in md


def test_render_contains_disclaimer() -> None:
    report = generate_report([], {}, "X", "")
    md = render_markdown_report(report)
    assert "research use only" in md.lower()


def test_render_variant_row() -> None:
    v = _make_variant()
    report = generate_report([v], {}, "S1", "ACGT")
    md = render_markdown_report(report)
    assert "c.20A>T" in md
    assert "HbS" in md
    assert "p.Glu7Val" in md


def test_render_no_variants_message() -> None:
    report = generate_report([], {}, "S1", "")
    md = render_markdown_report(report)
    assert "No variants detected" in md


def test_render_manual_review_flag() -> None:
    v = _make_variant(requires_manual_review=True, review_reason="Forward-strand only")
    report = generate_report([v], {}, "S1", "ACGT")
    md = render_markdown_report(report)
    assert "Forward-strand only" in md


def test_render_known_variant_section() -> None:
    v = _make_variant(known_variant_name="HbS")
    report = generate_report([v], {}, "S1", "ACGT")
    md = render_markdown_report(report)
    assert "Clinically Significant" in md
    assert "HbS" in md
