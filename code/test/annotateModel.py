"""Attach database cross-references from the annotation tables to a cobra model.

The committed YAML model carries only ids and names; the cross-references live in
model/metabolites.tsv, model/reactions.tsv and model/genes.tsv. MEMOTE reads
annotations from the SBML export, so without them it scores every annotation test
0% even though Human-GEM has the cross-references. This enriches an in-memory model
- used only to build the temporary SBML that MEMOTE runs on, never committed - so
the MEMOTE annotation scores reflect the identifiers Human-GEM actually carries.

Only registry (identifiers.org) namespaces MEMOTE can validate are attached, and
values are normalised to each namespace's expected form (e.g. Rhea drops its
"RHEA:" prefix; KEGG metabolite ids are split into compound/glycan/drug by prefix).
Legacy-only columns (EHMN, HepatoNET1, Recon3D, HMR2, Ratcon) are skipped.

Usage:
    import annotateModel
    annotateModel.enrich(model)   # mutates the model in place
"""

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
METABOLITES_TSV = REPO_ROOT / "model" / "metabolites.tsv"
REACTIONS_TSV = REPO_ROOT / "model" / "reactions.tsv"
GENES_TSV = REPO_ROOT / "model" / "genes.tsv"

# TSV column -> (identifiers.org namespace, prefix to strip so the value conforms).
MET_MAP = {
    "metBiGGID": ("bigg.metabolite", None),
    "metHMDBID": ("hmdb", None),
    "metChEBIID": ("chebi", None),            # values already "CHEBI:1234"
    "metPubChemID": ("pubchem.compound", None),
    "metLipidMapsID": ("lipidmaps", None),
    "metMetaNetXID": ("metanetx.chemical", None),
    "metSeedID": ("seed.compound", None),
}
RXN_MAP = {
    "rxnKEGGID": ("kegg.reaction", None),
    "rxnBiGGID": ("bigg.reaction", None),
    "rxnMetaNetXID": ("metanetx.reaction", None),
    "rxnRheaID": ("rhea", "RHEA:"),           # identifiers.org rhea wants the bare number
    "rxnREACTOMEID": ("reactome", None),
    "rxnTCDBID": ("tcdb", None),
}
GENE_MAP = {
    "geneUniProtID": ("uniprot", None),
    "geneEntrezID": ("ncbigene", None),
}
# KEGG metabolite ids span three namespaces, told apart by their first letter.
KEGG_MET_NS = {"C": "kegg.compound", "G": "kegg.glycan", "D": "kegg.drug"}


def _read(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _values(raw: str, strip_prefix) -> list[str]:
    parts = [p.strip() for p in (raw or "").split(";") if p.strip()]
    if strip_prefix:
        parts = [p[len(strip_prefix):] if p.startswith(strip_prefix) else p for p in parts]
    return parts


def _add(ann: dict, namespace: str, value: str) -> None:
    """Add one identifier under a namespace, keeping a list when there are several."""
    cur = ann.get(namespace)
    if cur is None:
        ann[namespace] = value
    elif isinstance(cur, list):
        if value not in cur:
            cur.append(value)
    elif cur != value:
        ann[namespace] = [cur, value]


def _apply(entity, row: dict, mapping: dict) -> bool:
    ann = dict(entity.annotation or {})
    before = len(ann)
    for column, (namespace, strip) in mapping.items():
        for value in _values(row.get(column, ""), strip):
            _add(ann, namespace, value)
    entity.annotation = ann
    return len(ann) > before


def enrich(model) -> dict:
    """Attach cross-references to every metabolite/reaction/gene that has a table
    row. Returns {'metabolites': n, 'reactions': n, 'genes': n} annotated."""
    mets = {r["mets"]: r for r in _read(METABOLITES_TSV)}
    rxns = {r["rxns"]: r for r in _read(REACTIONS_TSV)}
    genes = {r["genes"]: r for r in _read(GENES_TSV)}
    counts = {"metabolites": 0, "reactions": 0, "genes": 0}

    for met in model.metabolites:
        row = mets.get(met.id)
        if not row:
            continue
        touched = _apply(met, row, MET_MAP)
        ann = dict(met.annotation or {})
        for value in _values(row.get("metKEGGID", ""), None):
            namespace = KEGG_MET_NS.get(value[:1])
            if namespace:
                _add(ann, namespace, value)
                touched = True
        met.annotation = ann
        counts["metabolites"] += touched

    for rxn in model.reactions:
        row = rxns.get(rxn.id)
        if row and _apply(rxn, row, RXN_MAP):
            counts["reactions"] += 1

    for gene in model.genes:
        ann = dict(gene.annotation or {})
        # The gene id is itself the Ensembl gene identifier.
        if gene.id.startswith("ENSG"):
            _add(ann, "ensembl", gene.id)
        gene.annotation = ann
        touched = bool(ann)
        row = genes.get(gene.id)
        if row:
            touched = _apply(gene, row, GENE_MAP) or touched
        counts["genes"] += touched

    return counts
