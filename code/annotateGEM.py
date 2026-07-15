"""Attach cross-reference (MIRIAM) and SBO annotation to a Human-GEM model.

Python/raven-toolbox port of ``code/annotateGEM.m``.

The YAML model stores only the inline fields (``eccodes``, ``metFrom``,
``smiles``); the full set of external identifiers lives in the annotation
tables ``model/reactions.tsv``, ``model/metabolites.tsv`` and
``model/genes.tsv``. This module reads those tables and writes the identifiers
onto each cobra entity's ``annotation`` dict (namespace -> list of ids). SBO terms
for metabolites and reactions come from raven_toolbox's canonical ``add_sbo_terms``
(classifying exchange/demand/sink, transport, biomass, simple chemical, ...); genes,
which that helper does not cover, get SBO:0000243 here. The exported SBML / Excel /
txt then carry the annotation, while the YAML and ``.mat`` files stay
annotation-light (their cross-references remain the TSV tables), exactly as the
MATLAB release flow produces them.

Used by ``code/io/increaseHumanGEMVersion.py``; can also be run standalone to
inspect the merge:

    python code/annotateGEM.py            # counts only, model not modified on disk
"""
from __future__ import annotations

from pathlib import Path

import cobra
import pandas as pd
from raven_toolbox.annotation import add_sbo_terms

# Map TSV column names to identifiers.org namespaces (from annotateGEM.m id2miriam).
_RXN_ID2MIRIAM = {
    "rxnKEGGID": "kegg.reaction",
    "rxnBiGGID": "bigg.reaction",
    "rxnREACTOMEID": "reactome",
    "rxnRecon3DID": "vmhreaction",
    "rxnMetaNetXID": "metanetx.reaction",
    "rxnTCDBID": "tcdb",
    "rxnRheaID": "rhea",
    "rxnRheaMasterID": "rhea",
}
_MET_ID2MIRIAM = {
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
_GENE_ID2MIRIAM = {
    "genes": "ensembl",
    "geneENSTID": "ensembl",
    "geneENSPID": "ensembl",
    "geneUniProtID": "uniprot",
    "geneSymbols": "hgnc.symbol",
    "geneEntrezID": "ncbigene",
}

# Reaction and metabolite SBO terms come from raven_toolbox.annotation.add_sbo_terms
# (the canonical assignment); it does not cover genes, so gene SBO is set here.
_SBO_GENE = "SBO:0000243"  # gene
# Human-GEM's biomass (objective) reaction, so add_sbo_terms tags it SBO:0000629.
_BIOMASS_RXN_NAME = "Generic human cell biomass reaction"


def _read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV annotation table as text (empty cells become ``""``)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _split_ids(cell: str) -> list[str]:
    """Split a ``";"``-separated annotation cell into clean, non-empty ids."""
    return [part.strip() for part in str(cell).split(";") if part.strip()]


def _chebi(ids: list[str]) -> list[str]:
    """Ensure every ChEBI id has the ``CHEBI:`` prefix (annotateGEM.m rule)."""
    out = []
    for i in ids:
        out.append(i if i.upper().startswith("CHEBI:") else f"CHEBI:{i}")
    return out


def _rhea(ids: list[str]) -> list[str]:
    """Strip the ``RHEA:`` prefix; it must not appear in the identifiers.org URL."""
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
        # Dedupe while preserving first-seen order (columns can share a namespace).
        annotation[namespace] = list(dict.fromkeys(merged))


def annotate_gem(
    model: cobra.Model,
    model_dir: str | Path,
    *,
    types: tuple[str, ...] = ("rxn", "met", "gene"),
) -> cobra.Model:
    """Merge the TSV cross-references and SBO terms into ``model`` in place.

    Parameters
    ----------
    model
        Model whose reactions/metabolites/genes carry Human-GEM ids.
    model_dir
        Directory holding ``reactions.tsv`` / ``metabolites.tsv`` / ``genes.tsv``.
    types
        Which annotation classes to add (``"rxn"``, ``"met"``, ``"gene"``).

    Returns
    -------
    cobra.Model
        The same ``model`` object, now annotated. Pass a copy if the caller
        needs to keep an un-annotated version (the release keeps the YAML/.mat
        exports annotation-light this way).
    """
    model_dir = Path(model_dir)

    if "met" in types:
        mets = _read_tsv(model_dir / "metabolites.tsv").set_index("mets")
        for met in model.metabolites:
            if met.id in mets.index:
                _apply_row(met.annotation, mets.loc[met.id], _MET_ID2MIRIAM)

    if "rxn" in types:
        rxns = _read_tsv(model_dir / "reactions.tsv").set_index("rxns")
        for rxn in model.reactions:
            if rxn.id in rxns.index:
                _apply_row(rxn.annotation, rxns.loc[rxn.id], _RXN_ID2MIRIAM)

    if "gene" in types:
        genes = _read_tsv(model_dir / "genes.tsv").set_index("genes")
        for gene in model.genes:
            if gene.id in genes.index:
                row = genes.loc[gene.id].copy()
                row["genes"] = gene.id  # the gene id itself is an ensembl id
                _apply_row(gene.annotation, row, _GENE_ID2MIRIAM)
            # add_sbo_terms below covers metabolites and reactions, not genes.
            gene.annotation["sbo"] = _SBO_GENE

    # SBO terms for metabolites and reactions: the canonical raven-toolbox
    # assignment (exchange/demand/sink, transport, biomass, simple chemical, ...).
    if "met" in types or "rxn" in types:
        add_sbo_terms(model, biomass_rxn_name=_BIOMASS_RXN_NAME)

    return model


def _main() -> int:
    from raven_toolbox.io import read_yaml_model

    repo_root = Path(__file__).resolve().parents[1]
    model_dir = repo_root / "model"
    model = read_yaml_model(model_dir / "Human-GEM.yml")
    annotate_gem(model, model_dir)
    n_rxn = sum(1 for r in model.reactions if any(k != "sbo" for k in r.annotation))
    n_met = sum(1 for m in model.metabolites if m.annotation)
    n_gene = sum(1 for g in model.genes if g.annotation)
    print(f"annotated reactions (cross-refs): {n_rxn}/{len(model.reactions)}")
    print(f"annotated metabolites:            {n_met}/{len(model.metabolites)}")
    print(f"annotated genes:                  {n_gene}/{len(model.genes)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
