"""gpr_redundancy flags GPR clauses that cannot change the value of the rule."""

import cobra
import pytest

import qcModelChecks
from qcModelChecks import gpr_redundancy


def test_plain_rules_have_no_issue():
    assert gpr_redundancy("") == []
    assert gpr_redundancy("G1") == []
    assert gpr_redundancy("G1 or G2 or G3") == []
    assert gpr_redundancy("(G1 and G2) or (G1 and G3)") == []
    assert gpr_redundancy("G1 and G2 and G3") == []


def test_gene_listed_twice_in_an_or_rule():
    assert gpr_redundancy("G1 or G1 or G2") == ["repeated clause G1"]


def test_gene_listed_twice_in_a_complex():
    assert gpr_redundancy("G1 and G2 and G1") == ["repeated gene G1"]


def test_clause_absorbed_by_a_smaller_clause():
    assert gpr_redundancy("(G1 and G2 and G3) or (G1 and G2)") == ["absorbed clause G1 and G2 and G3"]
    assert gpr_redundancy("G1 or (G1 and G2)") == ["absorbed clause G1 and G2"]


def test_parentheses_around_the_whole_rule_are_looked_through():
    assert gpr_redundancy("((G1 and G2) or (G1 and G2))") == ["repeated clause G1 and G2"]
    assert gpr_redundancy("(G1 or G2)") == []
    assert gpr_redundancy("(G1 and G2) or (G3)") == []


def test_nested_or_inside_and_is_not_expanded():
    assert gpr_redundancy("G1 and (G2 or G3)") == []
    assert gpr_redundancy("(G1 and (G2 or G3)) or G1") == []


@pytest.fixture
def sanity_csv(tmp_path, monkeypatch):
    genes = tmp_path / "genes.tsv"
    genes.write_text("genes\nG1\nG2\n", encoding="utf-8")
    monkeypatch.setattr(qcModelChecks, "GENES_TSV", str(genes))
    path = tmp_path / "qc_reaction_sanity.csv"
    monkeypatch.setattr(qcModelChecks, "REACTION_SANITY_CSV", str(path))
    return path


def test_reaction_sanity_reports_a_redundant_gpr(sanity_csv):
    model = cobra.Model("test")
    reaction = cobra.Reaction("MAR00001", name="a reaction", lower_bound=0, upper_bound=1000)
    reaction.add_metabolites({cobra.Metabolite("MAM00001c", compartment="c"): -1,
                              cobra.Metabolite("MAM00002c", compartment="c"): 1})
    model.add_reactions([reaction])
    reaction.gene_reaction_rule = "G1 or G1 or G2"
    assert qcModelChecks.check_reaction_sanity(model) == 1
    assert "redundant_gpr:repeated clause G1" in sanity_csv.read_text(encoding="utf-8")
