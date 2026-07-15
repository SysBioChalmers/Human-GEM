"""Attach cross-reference (MIRIAM) and SBO annotation to a Human-GEM model.

Python/raven-toolbox port of ``code/annotateGEM.m``.

The YAML model stores only the inline fields (``eccodes``, ``metFrom``,
``smiles``); the full set of external identifiers lives in the annotation
tables ``model/reactions.tsv``, ``model/metabolites.tsv`` and
``model/genes.tsv``. This module reads those tables and writes the identifiers
onto each cobra entity's ``annotation`` dict (namespace -> list of ids), and
assigns an SBO term to every reaction (classified into
biochemical/transport/exchange/demand/sink/biomass), metabolite (simple chemical)
and gene. The exported SBML / Excel / txt then carry the annotation, while the YAML
and ``.mat`` files stay annotation-light (their cross-references remain the TSV
tables), exactly as the MATLAB release flow produces them.

Used by ``code/io/increaseHumanGEMVersion.py``; can also be run standalone to
inspect the merge:

    python code/annotateGEM.py            # counts only, model not modified on disk
"""
from __future__ import annotations

from pathlib import Path

import cobra
import pandas as pd

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

# Reaction SBO terms (annotateGEM.m). Precedence low -> high: default, biomass,
# transport, then boundary (a later match overrides an earlier one). Boundary
# reactions are split into exchange/demand/sink, matching how cobra (and hence
# MEMOTE) classifies them, so each gets the SBO term its check expects.
_SBO_DEFAULT = "SBO:0000176"    # biochemical reaction
_SBO_BIOMASS = "SBO:0000629"    # biomass production
_SBO_TRANSPORT = "SBO:0000185"  # translocation reaction
_SBO_EXCHANGE = "SBO:0000627"   # exchange reaction
_SBO_DEMAND = "SBO:0000628"     # demand reaction
_SBO_SINK = "SBO:0000632"       # sink reaction
# Metabolite and gene SBO terms (one each; no classification needed).
_SBO_METABOLITE = "SBO:0000247"  # simple chemical
_SBO_GENE = "SBO:0000243"        # gene


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


def _is_transport(rxn: cobra.Reaction) -> bool:
    """True if a metabolite name appears in more than one compartment (RAVEN
    getTransportRxns): the reaction moves a species across compartments."""
    comps_by_name: dict[str, set[str]] = {}
    for met in rxn.metabolites:
        comps_by_name.setdefault(met.name, set()).add(met.compartment)
    return any(len(comps) > 1 for comps in comps_by_name.values())


def _boundary_sbo(model: cobra.Model) -> dict:
    """Map each boundary reaction to its SBO term (exchange / demand / sink).

    Uses cobra's own ``exchanges`` / ``demands`` / ``sinks`` classification, which
    MEMOTE also uses, so the assigned term matches the check MEMOTE will apply. If
    cobra cannot classify (e.g. it fails to find an external compartment, or the
    model type lacks these properties), fall back to treating every boundary
    reaction as an exchange, as the original port did."""
    try:
        exchanges, demands, sinks = model.exchanges, model.demands, model.sinks
    except Exception:  # noqa: BLE001 - any classification failure -> safe fallback
        exchanges, demands, sinks = model.boundary, [], []
    sbo = {}
    for rxn in exchanges:
        sbo[rxn.id] = _SBO_EXCHANGE
    for rxn in demands:
        sbo[rxn.id] = _SBO_DEMAND
    for rxn in sinks:
        sbo[rxn.id] = _SBO_SINK
    # Any boundary reaction cobra did not place lands as an exchange.
    for rxn in model.boundary:
        sbo.setdefault(rxn.id, _SBO_EXCHANGE)
    return sbo


def _assign_sbo(model: cobra.Model) -> None:
    """Set ``annotation['sbo']`` on every reaction (annotateGEM.m SBO logic)."""
    boundary_sbo = _boundary_sbo(model)
    for rxn in model.reactions:
        sbo = _SBO_DEFAULT
        is_biomass = "biomass" in rxn.id.lower() or any(
            met.name.lower() == "biomass" and coeff > 0
            for met, coeff in rxn.metabolites.items()
        )
        if is_biomass:
            sbo = _SBO_BIOMASS
        if _is_transport(rxn):
            sbo = _SBO_TRANSPORT
        if rxn.id in boundary_sbo:
            sbo = boundary_sbo[rxn.id]
        rxn.annotation["sbo"] = sbo


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
            met.annotation["sbo"] = _SBO_METABOLITE

    if "rxn" in types:
        rxns = _read_tsv(model_dir / "reactions.tsv").set_index("rxns")
        for rxn in model.reactions:
            if rxn.id in rxns.index:
                _apply_row(rxn.annotation, rxns.loc[rxn.id], _RXN_ID2MIRIAM)
        _assign_sbo(model)

    if "gene" in types:
        genes = _read_tsv(model_dir / "genes.tsv").set_index("genes")
        for gene in model.genes:
            if gene.id in genes.index:
                row = genes.loc[gene.id].copy()
                row["genes"] = gene.id  # the gene id itself is an ensembl id
                _apply_row(gene.annotation, row, _GENE_ID2MIRIAM)
            gene.annotation["sbo"] = _SBO_GENE

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
