# HBBseq

**Sanger sequencing analysis for beta-thalassaemia and sickle cell disease.**  
Automated variant calling, HGVS annotation, zygosity classification, and interactive trace visualisation — in a browser, with no bioinformatics expertise required.

[![CI](https://github.com/abhirup780/hbbseq/actions/workflows/ci.yml/badge.svg)](https://github.com/abhirup780/hbbseq/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://hbb-sanger.streamlit.app/)

---

## Overview

HBBseq is a research-grade analysis pipeline for **HBB** (beta-globin) Sanger sequencing data. It accepts a forward and reverse `.ab1` trace pair, performs dual-strand alignment against the HBB locus, and produces a structured variant report with clinical annotation.

The tool is designed for research laboratories that use Sanger sequencing to screen for haemoglobin disorders. It removes the manual steps of visually inspecting chromatograms for known pathogenic sites and cross-referencing variant databases.

**Live demo:** [https://hbb-sanger.streamlit.app/](https://hbb-sanger.streamlit.app/)

---

## Features

| | Capability | Detail |
|---|---|---|
| **Detection** | SNV calling | Pre-alignment IUPAC ambiguity coding captures heterozygous calls that the basecaller suppresses |
| | Indel calling | Alignment-gap-based calling on dual-strand consensus |
| | Zygosity | Secondary peak ratio analysis across both strands (threshold: 0.25) |
| **Annotation** | HGVS notation | Programmatic c. and p. coordinates — NM_000518.5 reference |
| | Known variant lookup | 1,087-entry registry sourced from [Ithanet](https://www.ithanet.eu/), keyed by HGVS c. string |
| | Clinical significance | Pathogenic / Benign / Modifier classification with population frequency data |
| **Quality** | Mott trimming | Q20-based read-end trimming before alignment |
| | Artifact detection | Dye blob, high-noise, and signal dropoff warnings |
| | Manual review flags | Single-strand calls and sub-threshold signals flagged automatically |
| **Output** | Interactive UI | Dual chromatogram viewer, gene coverage map, per-variant cards |
| | Report export | Markdown report and structured JSON |

---

## Quickstart

### Web application

```bash
pip install -r requirements.txt
streamlit run app.py
```

Upload a forward and reverse `.ab1` file → click **Run Analysis**.

### Command-line interface

```bash
# Full analysis — write Markdown report
python cli.py run forward.ab1 reverse.ab1 --out report.md

# Full analysis — write JSON output
python cli.py run forward.ab1 reverse.ab1 --json result.json

# Validate reference FASTA structure
python cli.py validate-reference reference/HBB_reference.fasta
```

---

## How it works

```
Forward .ab1                   Reverse .ab1
     │                               │
     └──────────┬────────────────────┘
                │
     apply_iupac_symbols()
     Heterozygous positions → IUPAC ambiguity codes
     Phred score preserved through quality trim
                │
     trim_by_quality()
     Mott's algorithm, Q20 threshold
                │
     align_to_reference()
     Biopython PairwiseAligner, local alignment
     2,750 bp HBB locus reference (NM_000518.5)
                │
     evaluate_trace_artifacts()
     Dye blobs · noise ratio · signal dropoff
                │
     build_consensus()
     Dual-strand merge with quality reconciliation
                │
     call_variants_from_alignment()
     IUPAC codes  → heterozygous SNV calls
     Alignment gaps → indel calls
     Known variant lookup → clinical annotation
                │
     generate_report()
     HGVS c./p. · zygosity · manual review flags
     Markdown + JSON output
```

---

## Variant registry

The bundled registry (`hbb_pipeline/known_variants.py`) contains **1,087 HBB variants** sourced from [Ithanet](https://www.ithanet.eu/), keyed by HGVS c. notation. Commonly screened variants include:

| Variant | HGVS | Disease | Populations |
|---|---|---|---|
| HbS | c.20A>T | Sickle cell disease | Sub-Saharan African, Mediterranean |
| HbC | c.19G>A | Sickle cell / HbC disease | West African |
| IVS1-110 | c.93-21G>A | Beta-thalassaemia | Mediterranean, South Asian |
| IVS1-1 | c.92+1G>A | Beta-thalassaemia | Mediterranean, Middle Eastern |
| Cd39 | c.118C>T | Beta-thalassaemia | Mediterranean (Sardinian, Italian) |
| Cd6 (-A) | c.20delA | Beta-thalassaemia | South Asian, Southeast Asian |
| IVSI-5 | c.92+5G>C | Beta-thalassaemia | South Asian, Mediterranean |

The registry is regenerated from the Ithanet CSV export using `generate_variants.py`.

---

## Reference sequence

The bundled `reference/HBB_reference.fasta` is a **case-annotated 2,750 bp** sequence spanning the complete HBB locus. Case annotation (upper/lower) encodes functional regions and drives coordinate translation.

| Region | Coordinates | Length |
|---|---|---|
| 5′ flank | 0–977 | 978 bp |
| Exon 1 | 978–1071 | 94 bp |
| Intron 1 | 1072–1541 | 470 bp |
| Exon 2 | 1542–1694 | 153 bp |
| Intron 2 | 1695–2148 | 454 bp |
| Exon 3 | 2149–2279 | 131 bp |
| 3′ flank | 2280–2749 | 470 bp |

HGVS coordinates follow **NM_000518.5** (CDS start = genomic position 978, Met = codon 1).

---

## Project structure

```
hbbseq/
├── app.py                   # Streamlit web application
├── cli.py                   # Command-line interface (Typer)
├── plots.py                 # Plotly chromatogram and coverage visualisations
├── requirements.txt         # Runtime dependencies
├── requirements-dev.txt     # Development dependencies (pytest)
│
├── hbb_pipeline/            # Core analysis engine
│   ├── alignment.py         # Pairwise alignment and IUPAC pre-coding
│   ├── coordinates.py       # Genomic ↔ HGVS coordinate translation
│   ├── heterozygosity.py    # Secondary peak detection and zygosity
│   ├── known_variants.py    # Ithanet-sourced pathogenic variant registry
│   ├── models.py            # Pydantic v2 data models
│   ├── parsing.py           # ABI trace parser and Mott quality trimming
│   ├── pipeline.py          # Shared trace-processing core
│   ├── qc.py                # Trace artifact detection
│   ├── reference.py         # HBB reference FASTA loader and validator
│   ├── reporting.py         # Markdown and JSON report generation
│   └── variants.py          # Variant calling and HGVS annotation
│
├── reference/
│   └── HBB_reference.fasta  # Case-annotated 2,750 bp HBB locus reference
│
├── generate_variants.py     # Utility: regenerate known_variants.py from Ithanet CSV
│
└── tests/                   # pytest test suite
```

---

## Development

```bash
# Install with dev dependencies
pip install -r requirements-dev.txt

# Run tests
pytest --tb=short -q
```

---

## Disclaimer

This tool is intended **for research use only**. It has not been validated as an in vitro diagnostic device. Variant calls — particularly those flagged for manual review — must be confirmed by an accredited clinical laboratory before informing any clinical or reproductive decision.

---

## License

MIT © 2026 Abhirup Sarkar
