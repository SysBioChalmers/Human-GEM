"""Merge the annotation tables (cross-references + SBO) onto a Human-GEM model.

The YAML model stores only inline fields (`eccodes`, `metFrom`, `smiles`); the
external identifiers live in `reactions.tsv` / `metabolites.tsv` / `genes.tsv`.
This reads those tables and writes the ids onto each cobra entity's
`annotation` dict (namespace -> list of ids), then adds SBO terms via
raven-toolbox's canonical `add_sbo_terms`.

This is the packaged home of the logic that also lives in `code/annotateGEM.py`
(the CI helper); that module is intended to become a thin shim importing from
here so there is a single source of truth.
"""
from __future__ import annotations

from pathlib import Path

import cobra
import pandas as pd
from raven_toolbox.annotation import add_sbo_terms

# TSV column -> identifiers.org namespace (from annotateGEM.m id2miriam).
RXN_ID2MIRIAM = {
    "rxnKEGGID": "kegg.reaction",
    "rxnBiGGID": "bigg.reaction",
    "rxnREACTOMEID": "reactome",
    "rxnRecon3DID": "vmhreaction",
    "rxnMetaNetXID": "metanetx.reaction",
    "rxnTCDBID": "tcdb",
    "rxnRheaID": "rhea",
    "rxnRheaMasterID": "rhea",
}
MET_ID2MIRIAM = {
    "metBiGGID": "bigg.metabolite",
    "metKEGGID": "kegg.compound",
    "metHMDBID": "hmdb",
    "metChEBIID": "chebi",
    "metPubChemID": "pubchem.compound",
    "metLipidMapsID": "lipidmaps",
    "metRecon3DID": "vmhmetabolite",
    "metMetaNetXID": "metanetx.chemical",
    "metSeedID": "seed.compound",
}
GENE_ID2MIRIAM = {
    "genes": "ensembl",
    "geneENSTID": "ensembl",
    "geneENSPID": "ensembl",
    "geneUniProtID": "uniprot",
    "geneSymbols": "hgnc.symbol",
    "geneEntrezID": "ncbigene",
}

_SBO_GENE = "SBO:0000243"  # gene; add_sbo_terms covers reactions/metabolites, not genes
_BIOMASS_RXN_NAME = "Generic human cell biomass reaction"


def _read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV annotation table as text (empty cells become ``""``)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _split_ids(cell: str) -> list[str]:
    """Split a ``";"``-separated annotation cell into clean, non-empty ids."""
    return [part.strip() for part in str(cell).split(";") if part.strip()]


def _chebi(ids: list[str]) -> list[str]:
    """Ensure every ChEBI id carries the ``CHEBI:`` prefix."""
    return [i if i.upper().startswith("CHEBI:") else f"CHEBI:{i}" for i in ids]


def _rhea(ids: list[str]) -> list[str]:
    """Strip the ``RHEA:`` prefix (it must not appear in the identifiers.org URL)."""
    return [i[5:] if i.upper().startswith("RHEA:") else i for i in ids]


def _apply_row(annotation: dict, row: pd.Series, id2miriam: dict) -> None:
    """Add the mapped id columns of ``row`` to ``annotation`` (namespace -> list)."""
    for column, namespace in id2miriam.items():
        if column not in row:
            continue
        ids = _split_ids(row[column])
        if namespace == "chebi":
            ids = _chebi(ids)
        elif namespace == "rhea":
            ids = _rhea(ids)
        if not ids:
            continue
        merged = list(annotation.get(namespace, []))
        merged.extend(ids)
        # Dedupe, preserving first-seen order (columns can share a namespace).
        annotation[namespace] = list(dict.fromkeys(merged))


def annotate_gem(
    model: cobra.Model,
    model_dir: str | Path,
    *,
    types: tuple[str, ...] = ("rxn", "met", "gene"),
) -> cobra.Model:
    """Merge the TSV cross-references and SBO terms into ``model`` in place.

    Returns the same ``model`` object (pass a copy to keep an un-annotated one).
    """
    model_dir = Path(model_dir)

    if "met" in types:
        mets = _read_tsv(model_dir / "metabolites.tsv").set_index("mets")
        for met in model.metabolites:
            if met.id in mets.index:
                _apply_row(met.annotation, mets.loc[met.id], MET_ID2MIRIAM)

    if "rxn" in types:
        rxns = _read_tsv(model_dir / "reactions.tsv").set_index("rxns")
        for rxn in model.reactions:
            if rxn.id in rxns.index:
                _apply_row(rxn.annotation, rxns.loc[rxn.id], RXN_ID2MIRIAM)

    if "gene" in types:
        genes = _read_tsv(model_dir / "genes.tsv").set_index("genes")
        for gene in model.genes:
            if gene.id in genes.index:
                row = genes.loc[gene.id].copy()
                row["genes"] = gene.id  # the gene id itself is an ensembl id
                _apply_row(gene.annotation, row, GENE_ID2MIRIAM)
            gene.annotation["sbo"] = _SBO_GENE

    if "met" in types or "rxn" in types:
        add_sbo_terms(model, biomass_rxn_name=_BIOMASS_RXN_NAME)

    return model
