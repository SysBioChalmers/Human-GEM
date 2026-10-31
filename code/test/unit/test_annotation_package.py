"""code/annotateGEM.py and human_gem.load_model() share one annotation implementation."""

import sys
from pathlib import Path

import cobra

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # code/, for annotateGEM.py
import annotateGEM  # noqa: E402  (puts the repository root, with human_gem, on the path)
import human_gem.annotation  # noqa: E402


def test_code_and_package_use_the_same_functions():
    assert annotateGEM.annotate_gem is human_gem.annotation.annotate_gem
    assert annotateGEM.add_inchis is human_gem.annotation.add_inchis
    assert human_gem.annotation.RXN_ID2MIRIAM["rxnSeedID"] == "seed.reaction"


def _tables(tmp_path):
    (tmp_path / "reactions.tsv").write_text(
        "rxns\trxnKEGGID\trxnSeedID\trxnRheaID\nR1\tR00001\trxn00001\tRHEA:10000\n", encoding="utf-8")
    (tmp_path / "metabolites.tsv").write_text(
        "mets\tmetChEBIID\tmetInChI\nM1c\t15422\tInChI=1S/X\n", encoding="utf-8")
    (tmp_path / "genes.tsv").write_text(
        "genes\tgeneSymbols\tgeneUniProtID\nENSG00000000001\tABC1\tP00001\n", encoding="utf-8")
    return tmp_path


def _model():
    model = cobra.Model("m")
    reaction = cobra.Reaction("R1")
    reaction.add_metabolites({cobra.Metabolite("M1c", compartment="c"): -1})
    reaction.gene_reaction_rule = "ENSG00000000001"
    model.add_reactions([reaction])
    return model


def test_annotate_gem_merges_every_table(tmp_path):
    model = human_gem.annotation.annotate_gem(_model(), _tables(tmp_path))
    rxn, met, gene = model.reactions.R1, model.metabolites.M1c, model.genes.ENSG00000000001
    assert rxn.annotation["seed.reaction"] == ["rxn00001"]
    assert rxn.annotation["rhea"] == ["10000"]
    assert met.annotation["chebi"] == ["CHEBI:15422"]
    assert met.annotation["inchi"] == ["InChI=1S/X"]
    assert met.notes["inchis"] == "InChI=1S/X"
    assert gene.annotation["hgnc.symbol"] == ["ABC1"]
    assert gene.annotation["sbo"] == "SBO:0000243"


def test_add_inchis_adds_only_the_inchi(tmp_path):
    model = human_gem.annotation.add_inchis(_model(), _tables(tmp_path))
    assert model.metabolites.M1c.notes["inchis"] == "InChI=1S/X"
    assert not model.metabolites.M1c.annotation
