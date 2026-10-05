"""Annotation / cross-reference validation for Human-GEM.

Two kinds of check over the annotation tables (model/metabolites.tsv,
model/reactions.tsv, model/genes.tsv), aimed at catching curation mistakes:

  * Format validation. Each row must have as many fields as the header, and each
    external-database identifier must match the format of its namespace (KEGG,
    ChEBI, HMDB, PubChem, MetaNetX, Rhea, LipidMaps, EHMN, HepatoNET1, Reactome,
    TCDB, Ensembl, UniProt, NCBI Gene). Ensembl ids are unversioned. Freeform
    namespaces such as BiGG and Recon3D are not format-checked.
  * Cross-compartment consistency. The same metabolite (same metsNoComp) in
    different compartments is the same chemical, so it should carry the same
    cross-references. When two compartments give different non-empty values for a
    namespace, that is almost always a mistake.

Findings are written to data/testResults/qc_annotation_issues.csv and a coverage
summary is printed, so a pull request that introduces a broken or inconsistent
cross-reference is visible in the committed diff. This is a report and does not
fail the build.

Usage:
    python code/test/annotationTest.py
"""

import csv
import re
import sys

METABOLITES_TSV = "model/metabolites.tsv"
REACTIONS_TSV = "model/reactions.tsv"
GENES_TSV = "model/genes.tsv"
ISSUES_CSV = "data/testResults/qc_annotation_issues.csv"

# Namespace -> regex a single identifier must match (values may be ';'-separated).
# KEGG ids may be compounds (C), drugs (D) or glycans (G).
MET_PATTERNS = {
    "metKEGGID": r"[CDG]\d{5}",
    "metChEBIID": r"CHEBI:\d+",
    "metHMDBID": r"HMDB\d+",
    "metPubChemID": r"\d+",
    "metMetaNetXID": r"MNXM\d+|WATER",  # WATER is MetaNetX's id for H2O
    "metSeedID": r"cpd\d{5}",
    "metLipidMapsID": r"LM[A-Z]{2}\w+",
    "metEHMNID": r"C[A-Z]\d+",
    "metHepatoNET1ID": r"HC\d+",
}
RXN_PATTERNS = {
    "rxnKEGGID": r"R\d{5}",
    "rxnMetaNetXID": r"MNXR\d+",
    "rxnSeedID": r"rxn\d{5}",
    "rxnRheaID": r"RHEA:\d+",
    "rxnRheaMasterID": r"RHEA:\d+",
    "rxnREACTOMEID": r"R-HSA-\d+(?:\.\d+)?",
    "rxnHepatoNET1ID": r"r\d+",
    "rxnTCDBID": r"\d+\.[A-Z]\.\d+(?:\.\d+)*",
}
# Ensembl ids are stored without a version suffix (ENST00000375991, not ENST00000375991.9).
GENE_PATTERNS = {
    "geneENSTID": r"ENST\d{11}",
    "geneENSPID": r"ENSP\d{11}",
    "geneUniProtID": r"[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z0-9]{2}\d){1,2}",
    "geneEntrezID": r"\d+",
}
# Namespaces that should be identical across a metabolite's compartments.
CROSS_COMPARTMENT_COLS = [
    "metKEGGID", "metChEBIID", "metHMDBID", "metPubChemID", "metMetaNetXID", "metLipidMapsID", "metSeedID",
]


def _parts(value: str) -> set[str]:
    return {p.strip() for p in (value or "").split(";") if p.strip()}


def _read_table(path: str, id_col: str, issues: list) -> list:
    """Rows of an annotation table; a row with a different number of fields than the
    header is reported, because tools that read the table by position misplace its values."""
    with open(path, newline="", encoding="utf-8") as fh:
        lines = list(csv.reader(fh, delimiter="\t"))
    header, rows = lines[0], []
    for fields in lines[1:]:
        if not fields:
            continue
        if len(fields) != len(header):
            issues.append((fields[0], id_col, f"malformed: row has {len(fields)} fields, header has {len(header)}"))
        rows.append(dict(zip(header, fields)))
    return rows


def _check_format(rows: list, id_col: str, patterns: dict, issues: list) -> dict:
    total = len(rows)
    coverage = {}
    for column, pattern in patterns.items():
        rx = re.compile(f"^(?:{pattern})$")
        n_nonempty = 0
        for row in rows:
            value = (row.get(column) or "").strip()
            if not value:
                continue
            n_nonempty += 1
            for part in _parts(value):
                if not rx.match(part):
                    issues.append((row[id_col], column, f"malformed: {part}"))
        coverage[column] = (n_nonempty, total)
    return coverage


def _check_cross_compartment(rows: list, issues: list) -> None:
    groups: dict[str, list] = {}
    for row in rows:
        groups.setdefault(row.get("metsNoComp", ""), []).append(row)
    for base, members in groups.items():
        if not base or len(members) < 2:
            continue
        for column in CROSS_COMPARTMENT_COLS:
            values = [_parts(m.get(column, "")) for m in members]
            non_empty = [v for v in values if v]
            if len(non_empty) > 1 and any(v != non_empty[0] for v in non_empty):
                distinct = sorted({";".join(sorted(v)) for v in non_empty})
                issues.append((base, column, "inconsistent across compartments: " + " | ".join(distinct)))


def main() -> int:
    issues: list = []
    mets = _read_table(METABOLITES_TSV, "mets", issues)
    rxns = _read_table(REACTIONS_TSV, "rxns", issues)
    genes = _read_table(GENES_TSV, "genes", issues)

    coverage = {}
    coverage.update(_check_format(mets, "mets", MET_PATTERNS, issues))
    coverage.update(_check_format(rxns, "rxns", RXN_PATTERNS, issues))
    coverage.update(_check_format(genes, "genes", GENE_PATTERNS, issues))
    _check_cross_compartment(mets, issues)
    issues.sort()

    with open(ISSUES_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["id", "column", "issue"])
        writer.writerows(issues)

    n_malformed = sum(1 for _, _, issue in issues if issue.startswith("malformed"))
    n_inconsistent = sum(1 for _, _, issue in issues if issue.startswith("inconsistent"))
    print("Cross-reference coverage (non-empty / total):")
    for column, (n_nonempty, total) in coverage.items():
        pct = (100 * n_nonempty / total) if total else 0
        print(f"  {column}: {n_nonempty}/{total} ({pct:.0f}%)")
    print(f"Malformed cross-reference identifiers: {n_malformed}")
    print(f"Cross-references inconsistent across compartments: {n_inconsistent}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
