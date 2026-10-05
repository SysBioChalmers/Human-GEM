"""generateAnimalGEM: ortholog filtering, GPR rewriting, species network, annotation merge.

Everything here runs on small hand-built tables and models, without Human-GEM or a solver.
"""

import cobra
import pandas as pd
import pytest

import generateAnimalGEM as gen


def _write_alliance(path, rows):
    header = "fromGeneId\tfromSymbol\ttoGeneId\ttoSymbol\tbest\tbestReverse\tmethodCount\ttotalMethodCount"
    path.write_text("\n".join([header] + ["\t".join(map(str, r)) + "\t10" for r in rows]) + "\n")


def test_orthologs_single_hit_is_kept(tmp_path):
    f = tmp_path / "o.tsv"
    _write_alliance(f, [("H1", "A", "M1", "a", "Yes", "No", 5)])
    assert gen.read_alliance_orthologs(f).values.tolist() == [["A", "a"]]


def test_orthologs_neither_best_is_dropped(tmp_path):
    f = tmp_path / "o.tsv"
    _write_alliance(f, [("H1", "A", "M1", "a", "No", "No", 9), ("H2", "B", "M2", "b", "Yes", "Yes", 3)])
    assert gen.read_alliance_orthologs(f).values.tolist() == [["B", "b"]]


def test_orthologs_multiple_hits_prefer_best_both_ways(tmp_path):
    f = tmp_path / "o.tsv"
    _write_alliance(f, [("H1", "A", "M1", "a1", "Yes", "No", 9), ("H1", "A", "M2", "a2", "Yes", "Yes", 2),
                        ("H1", "A", "M3", "a3", "Yes", "Yes", 1)])
    assert gen.read_alliance_orthologs(f).values.tolist() == [["A", "a2"], ["A", "a3"]]


def test_orthologs_multiple_hits_fall_back_to_method_count(tmp_path):
    f = tmp_path / "o.tsv"
    _write_alliance(f, [("H1", "A", "M1", "a1", "Yes", "No", 2), ("H1", "A", "M2", "a2", "No", "Yes", 7)])
    assert gen.read_alliance_orthologs(f).values.tolist() == [["A", "a2"]]


@pytest.mark.parametrize("rule, mapping, expected", [
    ("G1", {"G1": ["a"]}, "a"),
    ("G1 and G2", {"G1": ["a"], "G2": ["b"]}, "a and b"),
    ("G1 or G2", {"G1": ["a"]}, "a"),                            # unmapped gene leaves the group
    ("G1 and G2", {"G1": ["a"]}, "a"),                           # complex keeps mapped subunits
    ("G1", {}, ""),
    ("G1 or G2", {"G1": ["a", "b"], "G2": ["b"]}, "a or b"),     # duplicates collapse
    ("(G1 or G2) and G3", {"G1": ["a"], "G2": ["b"], "G3": ["c"]}, "(a or b) and c"),
    ("G1 and G2", {"G1": ["a", "b"], "G2": ["c"]}, "(a or b) and c"),
    ("", {"G1": ["a"]}, ""),
])
def test_map_gene_reaction_rule(rule, mapping, expected):
    assert gen.map_gene_reaction_rule(rule, mapping) == expected


def _template():
    m = cobra.Model("t")
    mets = {k: cobra.Metabolite(k, compartment="c") for k in "xyz"}
    m.add_metabolites(list(mets.values()))
    for rid, rule, (a, b) in [("R1", "G1", "xy"), ("R2", "G2 or G3", "yz"), ("R3", "", "zx")]:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.add_metabolites({mets[a]: -1, mets[b]: 1})
        r.gene_reaction_rule = rule
    return m


def test_draft_rewrites_rules_and_drops_reactions_that_lose_their_genes():
    draft = gen.build_ortholog_draft(_template(), {"G2": ["b"], "G3": ["c"]})
    assert {r.id for r in draft.reactions} == {"R2", "R3"}      # R1 lost G1, R3 never had a rule
    assert draft.reactions.R2.gene_reaction_rule == "b or c"
    assert {g.id for g in draft.genes} == {"b", "c"}


def test_species_network_adds_mets_and_reactions():
    m = _template()
    mets = pd.DataFrame({"mets": ["MAM1c"], "metNames": ["foo"], "metFormulas": ["C2"],
                         "metCharges": ["0"], "compartments": ["c"]})
    rxns = pd.DataFrame({"rxns": ["MAR1"], "equations": ["foo[c] => foo[c]"],
                         "subSystems": ["S"], "grRules": ["g1"], "lb": ["0"], "ub": ["1000"],
                         "eccodes": ["1.1.1.1;2.2.2.2"], "rxnReferences": ["PMID:1"],
                         "rxnConfidenceScores": ["3"], "rxnNames": ["a rxn"]})
    m.metabolites.x.name = "bar"
    rxns.loc[0, "equations"] = "bar[c] => foo[c]"
    assert gen.add_species_network(m, rxns, mets) == ["MAR1"]
    r = m.reactions.MAR1
    assert (r.lower_bound, r.upper_bound) == (0, 1000)
    assert r.gene_reaction_rule == "g1"
    assert r.annotation["ec-code"] == ["1.1.1.1", "2.2.2.2"]
    assert r.notes["confidence_score"] == 3
    assert r.metabolites[m.metabolites.MAM1c] == 1


def test_species_network_refuses_existing_ids():
    m = _template()
    mets = pd.DataFrame({"mets": ["x"], "metNames": ["n"], "metFormulas": [""], "metCharges": ["0"],
                         "compartments": ["c"]})
    rxns = pd.DataFrame({"rxns": ["MAR1"], "equations": ["n[c] => n[c]"], "subSystems": [""], "grRules": [""]})
    with pytest.raises(ValueError, match="already in the model"):
        gen.add_species_network(m, rxns, mets)


def test_merge_annotation_appends_species_rows_and_follows_model_order():
    human = pd.DataFrame({"rxns": ["MAR1", "MAR2", "MAR3"], "rxnKEGGID": ["K1", "K2", "K3"]})
    species = pd.DataFrame({"rxns": ["MAR9"], "rxnKEGGID": ["K9"], "lb": ["0"]})
    out = gen.merge_annotation(human, species, "rxns", ["MAR9", "MAR1"])
    assert out.values.tolist() == [["MAR9", "K9"], ["MAR1", "K1"]]
    assert list(out.columns) == ["rxns", "rxnKEGGID"]


def test_merge_annotation_rejects_unknown_components():
    human = pd.DataFrame({"rxns": ["MAR1"], "rxnKEGGID": ["K1"]})
    with pytest.raises(ValueError, match="no annotation row"):
        gen.merge_annotation(human, human.iloc[0:0], "rxns", ["MAR1", "MAR7"])


def test_stamp_metadata(tmp_path):
    m = cobra.Model()
    repo = gen.AnimalRepo("Mouse", tmp_path)
    gen.stamp_metadata(m, repo, "1.9.0", "2026-10-05")
    assert m.id == "Mouse-GEM"
    assert m.notes["metaData"]["version"] == "1.9.0" == m.notes["version"]
    assert m.notes["metaData"]["date"] == "2026-10-05"
    assert m.notes["metaData"]["taxonomy"] == "10090"


def test_repo_paths(tmp_path):
    repo = gen.AnimalRepo("Fruitfly", tmp_path)
    assert repo.orthologs.name == "human2FruitflyOrthologs.tsv"
    assert repo.specific_rxns.name == "fruitflySpecificRxns.tsv"
    assert repo.specific_mets.name == "fruitflySpecificMets.tsv"


def _biomass_model():
    m = cobra.Model("h")
    lipoyl = cobra.Metabolite("L", name="[protein]-N6-(lipoyl)lysine", compartment="m")
    acid = cobra.Metabolite("A", name="lipoic acid", compartment="c")
    pool = cobra.Metabolite("P", name="cofactors and vitamins", compartment="c")
    r22, r65 = cobra.Reaction("MAR00022"), cobra.Reaction("MAR10065")
    m.add_reactions([r22, r65])
    r22.add_metabolites({lipoyl: -1, pool: 1})
    r65.add_metabolites({acid: -0.5, pool: 1})
    return m


def test_fix_lipoyl_biomass_swaps_lipoyl_lysine_for_lipoic_acid():
    m = _biomass_model()
    assert gen.fix_lipoyl_biomass(m) is True
    assert {x.id: c for x, c in m.reactions.MAR00022.metabolites.items()} == {"A": -1, "P": 1}
    assert gen.fix_lipoyl_biomass(m) is False       # nothing left to swap


def test_fix_lipoyl_biomass_leaves_other_models_alone():
    assert gen.fix_lipoyl_biomass(cobra.Model("x")) is False
