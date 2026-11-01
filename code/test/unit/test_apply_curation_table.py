"""applyCurationTable on small model files."""

import pytest

import applyCurationTable as act

YML = """!!omap
- metaData: !!omap
  - id: HumanGEM
- metabolites:
  - !!omap
    - id: MAM00001c
    - name: glycine
    - compartment: c
    - charge: 0.0
    - formula: C2H5NO2
    - metFrom: HMRdatabase
- reactions:
  - !!omap
    - id: MAR00001
    - name: transport of glycine
    - metabolites: !!omap
      - MAM00001c: -1.0
    - lower_bound: 0.0
    - upper_bound: 1000.0
    - subsystem:
      - Transport reactions
    - annotation: !!omap
      - ec-code:
        - 3.1.1.4
    - references: PMID:1
    - confidence_score: 0.0
  - !!omap
    - id: MAR00002
    - name: glycine kinase
    - metabolites: !!omap
      - MAM00001c: -1.0
    - lower_bound: 0.0
    - upper_bound: 1000.0
    - gene_reaction_rule: ENSG00000000001
    - subsystem:
      - Glycine metabolism
    - rxnFrom: HMRdatabase
    - confidence_score: 0.0
    - notes: "quoted: kept as written"
- genes:
  - !!omap
    - id: ENSG00000000001
    - name: GLYK
- compartments: !!omap
  - c: Cytosol
"""
RXN_TSV = "rxns\trxnKEGGID\trxnMetaNetXID\nMAR00001\t\tMNXR1\nMAR00002\tR00001\tMNXR2\n"
MET_TSV = "mets\tmetsNoComp\tmetChEBIID\tmetSmiles\nMAM00001c\tMAM00001\tCHEBI:1\tNCC(=O)O\n"
GENE_TSV = "genes\tgeneSymbols\nENSG00000000001\tGLYK\n"


@pytest.fixture
def root(tmp_path):
    (tmp_path / "model").mkdir()
    (tmp_path / "model" / "Human-GEM.yml").write_text(YML, encoding="utf-8")
    (tmp_path / "model" / "reactions.tsv").write_text(RXN_TSV, encoding="utf-8")
    (tmp_path / "model" / "metabolites.tsv").write_text(MET_TSV, encoding="utf-8")
    (tmp_path / "model" / "genes.tsv").write_text(GENE_TSV, encoding="utf-8")
    return tmp_path


def run(root, rows):
    rows = [dict(zip(("entity", "field", "before", "after"), r)) for r in rows]
    changed, problems, files = act.apply(rows, root)
    for f in files.values():
        f.save()
    return changed, problems


def yml(root):
    return (root / "model" / "Human-GEM.yml").read_text(encoding="utf-8")


def test_name_needing_quotes_is_double_quoted(root):
    changed, problems = run(root, [("MAR00002", "name", "glycine kinase", "kinase: glycine")])
    assert problems == [] and changed["name"] == 1
    assert '    - name: "kinase: glycine"' in yml(root)


def test_ec_removed_with_its_annotation_block(root):
    run(root, [("MAR00001", "ec-code", "3.1.1.4", "")])
    text = yml(root)
    assert "3.1.1.4" not in text and "annotation" not in text.split("MAR00002")[0]
    assert text.count("\n") == YML.count("\n") - 3


def test_ec_added_before_the_references_or_rxnfrom(root):
    run(root, [("MAR00002", "ec-code", "", "2.7.1.31;2.7.1.-")])
    lines = yml(root).split("\n")
    i = lines.index("    - annotation: !!omap", lines.index("    - id: MAR00002"))
    assert lines[i - 2:i + 5] == ["    - subsystem:", "      - Glycine metabolism", "    - annotation: !!omap",
                                  "      - ec-code:", "        - 2.7.1.31", "        - 2.7.1.-",
                                  "    - rxnFrom: HMRdatabase"]


def test_rxnfrom_added_before_confidence_score(root):
    run(root, [("MAR00001", "rxnFrom", "", "HumanGEM")])
    text = yml(root)
    assert "    - references: PMID:1\n    - rxnFrom: HumanGEM\n    - confidence_score: 0.0" in text


def test_metabolite_charge_formula_and_tables(root):
    changed, problems = run(root, [
        ("MAM00001c", "charge", "0", "-1"), ("MAM00001c", "formula", "C2H5NO2", "C2H4NO2"),
        ("MAM00001c", "metSmiles", "NCC(=O)O", "NCC(=O)[O-]"), ("MAR00001", "rxnKEGGID", "", "R1"),
        ("ENSG00000000001", "geneSymbols", "GLYK", "GLYCTK")])
    assert problems == []
    assert "    - charge: -1.0\n    - formula: C2H4NO2" in yml(root)
    assert "MAR00001\tR1\tMNXR1" in (root / "model" / "reactions.tsv").read_text(encoding="utf-8")
    assert "NCC(=O)[O-]" in (root / "model" / "metabolites.tsv").read_text(encoding="utf-8")


def test_untouched_lines_stay_as_written(root):
    run(root, [("MAR00002", "name", "glycine kinase", "glycerate kinase")])
    assert '    - notes: "quoted: kept as written"' in yml(root)
    assert yml(root).replace("glycerate kinase", "glycine kinase") == YML


def test_a_value_that_differs_from_before_is_reported(root):
    changed, problems = run(root, [("MAR00002", "name", "something else", "x")])
    assert problems == ["MAR00002 name: 'glycine kinase', expected 'something else'"]
    assert yml(root) == YML


def test_a_value_already_set_is_skipped(root):
    changed, problems = run(root, [("MAR00002", "rxnMetaNetXID", "MNXR9", "MNXR2")])
    assert problems == [] and not changed


def test_gene_reaction_rule_written_with_symbols(root):
    changed, problems = run(root, [("MAR00002", "gene_reaction_rule", "GLYK", "GLYK or GLYK")])
    assert problems == [] and changed["gene_reaction_rule"] == 1
    assert "    - gene_reaction_rule: ENSG00000000001 or ENSG00000000001\n" in yml(root)


def test_gene_reaction_rule_removed(root):
    run(root, [("MAR00002", "gene_reaction_rule", "GLYK", "")])
    assert "gene_reaction_rule" not in yml(root)


def test_gene_reaction_rule_with_unknown_symbol_is_an_error(root):
    with pytest.raises(ValueError):
        run(root, [("MAR00002", "gene_reaction_rule", "GLYK", "NOPE")])


def test_metabolites_replaced_in_the_given_order(root):
    changed, problems = run(root, [("MAR00001", "metabolites", "MAM00001c:-1.0", "MAM00001c:-1.0;MAM00001c:2.0")])
    assert problems == [] and changed["metabolites"] == 1
    assert "    - metabolites: !!omap\n      - MAM00001c: -1.0\n      - MAM00001c: 2.0\n    - lower_bound" in yml(root)


def test_references_added_before_rxnfrom(root):
    run(root, [("MAR00002", "references", "", "doi:10.1/x")])
    assert "    - references: doi:10.1/x\n    - rxnFrom: HMRdatabase" in yml(root)


def test_unused_gene_is_deleted_from_both_files(root):
    changed, problems = run(root, [("MAR00002", "gene_reaction_rule", "GLYK", ""),
                                   ("ENSG00000000001", "gene", "present", "")])
    assert problems == [] and changed["gene"] == 1
    assert "ENSG00000000001" not in yml(root)
    assert "ENSG00000000001" not in (root / "model" / "genes.tsv").read_text(encoding="utf-8")


def test_gene_still_in_a_rule_is_kept(root):
    changed, problems = run(root, [("ENSG00000000001", "gene", "present", "")])
    assert problems == ["ENSG00000000001 gene: still in the rule of MAR00002"]
    assert "ENSG00000000001" in yml(root)


def test_edits_that_change_the_line_count_leave_later_entries_findable(root):
    changed, problems = run(root, [
        ("MAR00001", "references", "PMID:1", ""),
        ("MAR00001", "gene_reaction_rule", "", "GLYK"),
        ("MAR00002", "references", "", "doi:10.1/x"),
        ("MAR00002", "name", "glycine kinase", "glycerate kinase")])
    assert problems == [] and sum(changed.values()) == 4
    text = yml(root)
    assert "    - references: PMID:1" not in text
    assert "    - name: glycerate kinase" in text and "    - references: doi:10.1/x\n    - rxnFrom" in text


def test_rule_row_is_skipped_when_a_gene_of_the_old_rule_was_deleted(root):
    run(root, [("MAR00002", "gene_reaction_rule", "GLYK", ""), ("ENSG00000000001", "gene", "present", "")])
    changed, problems = run(root, [("MAR00002", "gene_reaction_rule", "GLYK or GONE", "")])
    assert problems == [] and not changed


def test_bounds_changed_as_floats(root):
    changed, problems = run(root, [("MAR00001", "lower_bound", "0", "-1000"), ("MAR00001", "upper_bound", "1000.0", "500")])
    assert problems == [] and changed["lower_bound"] == 1 and changed["upper_bound"] == 1
    assert "    - lower_bound: -1000.0\n    - upper_bound: 500.0\n" in yml(root)


def test_bound_that_differs_from_before_is_reported(root):
    changed, problems = run(root, [("MAR00001", "lower_bound", "-5", "-1")])
    assert problems == ["MAR00001 lower_bound: '0.0', expected '-5'"]
