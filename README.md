# HBBseq

> Sanger sequencing analysis for **beta-thalassaemia** and **sickle cell disease** — variant calling, HGVS annotation, zygosity, and trace visualisation in a single web app.

[![CI](https://github.com/abhirup780/hbb-sanger/actions/workflows/ci.yml/badge.svg)](https://github.com/abhirup780/hbb-sanger/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Streamlit](https://img.shields.io/badge/Streamlit-Community%20Cloud-FF4B4B?logo=streamlit)](https://streamlit.io/cloud)

---

## Live Demo

**Try it online — no installation required:**  
[https://hbb-sanger.streamlit.app/](https://hbb-sanger.streamlit.app/)

---

## Features

| Capability | Detail |
|---|---|
| **SNV detection** | Pre-alignment IUPAC coding — catches HET calls the basecaller silences |
| **Indel detection** | Alignment-gap-based calling on dual-strand consensus |
| **HGVS annotation** | Programmatic c. and p. notation — no hardcoded variant lists |
| **Zygosity** | Dual-strand secondary peak ratios (threshold 0.25) |
| **Known variants** | 1087-entry Ithanet registry with automatic lookup and clinical significance |
| **Trace viewer** | Interactive Plotly chromatograms — forward and reverse simultaneously |
| **Coverage map** | Gene locus cartoon with per-variant lollipop markers (HGVS-labelled) |
| **QC report** | Phred scores, usable read length, artifact warnings |
| **Export** | Markdown report + JSON for downstream processing |

---

## Quickstart

### Web app

```bash
pip install -r requirements.txt
streamlit run app.py
```

Upload a forward + reverse `.ab1` pair → click **Run Analysis**.

### Command-line

```bash
pip install -r requirements.txt

# Run full analysis
python cli.py run forward.ab1 reverse.ab1 --out report.md
python cli.py run forward.ab1 reverse.ab1 --json result.json

# Validate reference FASTA
python cli.py validate-reference reference/HBB_reference.fasta
```

---

## Project structure

```
hbb-sanger/
├── app.py                  # Streamlit web app entry point
├── cli.py                  # Command-line interface (Typer)
├── plots.py                # Plotly chromatogram and coverage visualisations
├── requirements.txt
│
├── hbb_pipeline/           # Core analysis engine
│   ├── alignment.py        # Pairwise alignment + IUPAC pre-coding
│   ├── coordinates.py      # Genomic ↔ HGVS coordinate translation
│   ├── heterozygosity.py   # Secondary peak detection and zygosity
│   ├── known_variants.py   # Ithanet-sourced pathogenic variant registry (1087 entries)
│   ├── models.py           # Pydantic data models
│   ├── parsing.py          # ABI trace parser + Mott quality trimming
│   ├── pipeline.py         # Shared trace-processing core (app + CLI)
│   ├── qc.py               # Trace artifact detection (dye blobs, noise, dropoff)
│   ├── reference.py        # HBB reference FASTA loader + validator
│   ├── reporting.py        # Markdown / JSON report generation
│   └── variants.py         # Variant calling + HGVS annotation
│
├── reference/
│   └── HBB_reference.fasta # Case-annotated 2750 bp HBB locus reference
│
├── generate_variants.py    # Dev tool: regenerate known_variants.py from Ithanet CSV
│
└── tests/                  # pytest test suite
```

---

## How it works

```
ABI file (forward + reverse)
   │
   ├─ apply_iupac_symbols()       ← Secondary peak scan, IUPAC coding (pre-trim)
   │    Basecaller-silenced HET positions → IUPAC ambiguity codes
   │    Phred boosted to Q20 so Mott trimmer preserves them
   │
   ├─ trim_by_quality()           ← Mott's modified algorithm (Q20 default)
   │
   ├─ align_to_reference()        ← Biopython PairwiseAligner, local mode
   │
   ├─ build_consensus()           ← Dual-strand merge with quality reconciliation
   │
   ├─ evaluate_trace_artifacts()  ← Dye blob, noise ratio, signal dropoff checks
   │
   └─ call_variants_from_alignment()  ← SNV + indel calling from consensus
        IUPAC codes → HET SNV calls
        Alignment gaps → indel calls
        Known variant lookup → clinical annotation
        generate_report() → Markdown + JSON output
```

---

## Reference

The bundled `reference/HBB_reference.fasta` is a **case-annotated 2750 bp** sequence spanning the complete HBB locus:

| Region | Coordinates |
|---|---|
| 5′ flank | 0–977 |
| Exon 1 | 978–1071 |
| Intron 1 | 1072–1541 |
| Exon 2 | 1542–1694 |
| Intron 2 | 1695–2148 |
| Exon 3 | 2149–2279 |
| 3′ flank | 2280–2749 |

HGVS coordinates are relative to **NM_000518.5** (CDS start = genomic position 978).

---

## Known variant registry

`hbb_pipeline/known_variants.py` contains 1087 HBB variants sourced from [Ithanet](https://www.ithanet.eu/), keyed by HGVS c. notation. Each entry carries:

- Common name (e.g. HbS, IVS1-110, Cd39)
- Associated disease(s)
- Clinical significance (Pathogenic / Benign / Modifier)
- Affected populations (where recorded)

The registry is regenerated from the Ithanet CSV export using `generate_variants.py`.

---

## Disclaimer

> **For research use only.** This tool is not validated for clinical diagnostic use. All variant calls, particularly those flagged for manual review, must be confirmed by a certified clinical laboratory before any clinical decision is made.

---

## License

MIT © 2026 Abhirup Sarkar
