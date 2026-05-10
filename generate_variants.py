"""One-shot script: generate hbb_pipeline/known_variants.py from Ithanet CSV."""
import csv
import re
import sys

CSV_PATH = "HB_Variants_HGVS(ithanet).csv"
OUT_PATH = "hbb_pipeline/known_variants.py"

POPULATIONS = {
    "20A>T":          ["Sub-Saharan African", "Mediterranean", "Middle Eastern", "South Asian"],
    "19G>A":          ["West African"],
    "79G>A":          ["South Asian", "Southeast Asian"],
    "118C>T":         ["Mediterranean (Sardinian, Italian)"],
    "47G>A":          ["South Asian", "Southeast Asian"],
    "27dupG":         ["South Asian (Indian subcontinent)"],
    "126_129delCTTT": ["South Asian", "East Asian (Chinese)"],
    "92+1G>A":        ["Mediterranean", "Middle Eastern"],
    "92+5G>C":        ["South Asian (Indian subcontinent)", "Mediterranean"],
    "93-21G>A":       ["Mediterranean", "South Asian"],
    "315+1G>A":       ["Mediterranean", "Middle Eastern"],
    "316-197C>T":     ["East Asian (Chinese)", "Southeast Asian"],
    "316-2A>G":       ["Mediterranean"],
    "-78A>G":         ["South Asian (Indian subcontinent)", "African American"],
    "-79A>G":         ["African American", "Chinese"],
    "-87C>G":         ["Mediterranean", "South Asian"],
    "-88C>T":         ["African American"],
    "-140C>T":        ["Mediterranean", "Middle Eastern", "South Asian"],
}


def clean_phenotype(s: str) -> str:
    s = s.replace("?-thalassaemia", "beta-thalassaemia")
    s = s.replace("?-chain variant", "beta-chain variant")
    s = s.replace("Hb F levels", "elevated Hb F")
    return "" if s == "N/A" else s


def func_to_sig(f: str) -> str:
    return {"Causative": "Pathogenic", "Neutral": "Benign", "Modifier": "Modifier"}.get(f, f)


def build_name(common: str, hb: str) -> str:
    hb = hb.strip()
    common = common.strip()
    if hb and hb not in ("N/A", "?", ""):
        return f"{hb} ({common})"
    return common


def pystr(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


with open(CSV_PATH, encoding="latin-1") as f:
    rows = list(csv.DictReader(f))

entries: list[tuple] = []
seen_keys: set[str] = set()

for r in rows:
    hgvs = r["HGVS Name"].strip()
    if not hgvs.startswith("HBB:c."):
        continue
    if "[" in hgvs or ";" in hgvs:
        continue
    key = hgvs[len("HBB:c."):].split(" | ")[0].strip()
    if not key or "NG_" in key or key in seen_keys:
        continue
    seen_keys.add(key)
    entries.append((
        key,
        build_name(r["Common Name"], r["Hb Name"]),
        clean_phenotype(r["Phenotype"]),
        func_to_sig(r["Functionality"]),
        POPULATIONS.get(key, []),
    ))


def sort_key(e: tuple) -> tuple:
    m = re.match(r"^(-?\d+)", e[0])
    return (int(m.group(1)) if m else 99999, e[0])


entries.sort(key=sort_key)
print(f"Generating {len(entries)} entries...", file=sys.stderr)

HEADER = '''\
"""Lookup table of HBB variants for clinical annotation.

Sourced from Ithanet (ithanet.eu) HBB variant database ({count} entries).
Keys are HGVS c. strings (without the "c." prefix).
Population data retained for key clinically important variants.

Do not edit by hand — regenerate with generate_variants.py.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class KnownVariant:
    name: str
    disease: str
    clinical_significance: str  # "Pathogenic" | "Benign" | "Modifier"
    common_populations: list[str]


# Key: hgvs_c string (no "c." prefix)
KNOWN_VARIANTS: dict[str, KnownVariant] = {{
'''

FOOTER = '''\
}


def lookup_variant(hgvs_c: str) -> KnownVariant | None:
    """Look up a variant by c. coordinate string (without "c." prefix).

    Example:
        kv = lookup_variant("20A>T")
        assert kv is not None
        print(kv.name)  # "HbS (CD 6 GAG>GTG [Glu>Val])"
    """
    return KNOWN_VARIANTS.get(hgvs_c)
'''

with open(OUT_PATH, "w", encoding="utf-8") as out:
    out.write(HEADER.format(count=len(entries)))
    for key, name, disease, sig, pops in entries:
        out.write(
            f'    "{pystr(key)}": KnownVariant('
            f'name="{pystr(name)}", '
            f'disease="{pystr(disease)}", '
            f'clinical_significance="{sig}", '
            f'common_populations={repr(pops)}),\n'
        )
    out.write(FOOTER)

print(f"Written to {OUT_PATH}", file=sys.stderr)
