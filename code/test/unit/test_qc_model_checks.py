"""qcModelChecks: name consistency, identifier format and the growth medium."""

import cobra
import pytest

import qcModelChecks


@pytest.fixture
def csv_path(tmp_path, monkeypatch):
    """Send the check's CSV to a temporary file instead of data/testResults."""
    path = tmp_path / "qc_name_consistency.csv"
    monkeypatch.setattr(qcModelChecks, "NAME_CONSISTENCY_CSV", str(path))
    return path


def _model(metabolites, reaction_name="a reaction"):
    """A one-reaction model holding the given (id, name) metabolites."""
    model = cobra.Model("test")
    reaction = cobra.Reaction("MAR00001")
    reaction.name = reaction_name
    reaction.add_metabolites(
        {cobra.Metabolite(mid, name=name, compartment=mid[-1]): 1.0 for mid, name in metabolites}
    )
    model.add_reactions([reaction])
    return model


def test_consistent_names_produce_no_rows(csv_path):
    model = _model([("MAM00001c", "pyruvate"), ("MAM00001m", "pyruvate")])
    assert qcModelChecks.check_name_consistency(model) == []


def test_missing_metabolite_name_is_flagged(csv_path):
    model = _model([("MAM00001c", "")])
    assert ("metabolite", "MAM00001c", "missing name") in qcModelChecks.check_name_consistency(model)


def test_missing_reaction_name_is_flagged(csv_path):
    model = _model([("MAM00001c", "pyruvate")], reaction_name="")
    assert ("reaction", "MAR00001", "missing name") in qcModelChecks.check_name_consistency(model)


def test_name_differing_across_compartments_is_flagged(csv_path):
    model = _model([("MAM02766c", "pristanic acid"), ("MAM02766m", "(2S)-pristanic acid")])
    issues = [row for row in qcModelChecks.check_name_consistency(model) if row[1] == "MAM02766"]
    assert len(issues) == 1
    kind, _base, issue = issues[0]
    assert kind == "metabolite"
    assert issue.startswith("name differs across compartments")
    # Both spellings are named, so the report says what to reconcile.
    assert "pristanic acid" in issue and "(2S)-pristanic acid" in issue


def test_the_same_name_on_two_different_compounds_is_not_flagged(csv_path):
    """Distinct base identifiers may legitimately share a name (isomers, pseudo-metabolites)."""
    model = _model([("MAM00001c", "retinal"), ("MAM20002c", "retinal")])
    assert qcModelChecks.check_name_consistency(model) == []


def test_rows_are_written_to_the_csv(csv_path):
    model = _model([("MAM00001c", "")])
    qcModelChecks.check_name_consistency(model)
    contents = csv_path.read_text(encoding="utf-8").splitlines()
    assert contents[0] == "kind,id,issue"
    assert any("MAM00001c" in line for line in contents[1:])


@pytest.fixture
def tables(tmp_path, monkeypatch):
    """Model tables and deprecated lists in a temporary directory."""
    def write(name, header, ids):
        path = tmp_path / name
        empty = "\t" * header.count("\t")  # the other columns stay empty
        path.write_text(header + "\n" + "".join(f"{i}{empty}\n" for i in ids), encoding="utf-8")
        return str(path)

    def setup(rxns, mets, genes):
        monkeypatch.setattr(qcModelChecks, "REACTIONS_TSV", write("reactions.tsv", "rxns\tspontaneous", rxns))
        monkeypatch.setattr(qcModelChecks, "METABOLITES_TSV", write("metabolites.tsv", "mets", mets))
        monkeypatch.setattr(qcModelChecks, "GENES_TSV", write("genes.tsv", "genes", genes))
        monkeypatch.setattr(qcModelChecks, "DEPRECATED_RXN_TSV", write("depr.tsv", "rxns", []))
        monkeypatch.setattr(qcModelChecks, "DEPRECATED_MET_TSV", write("depm.tsv", "mets", []))
        monkeypatch.setattr(qcModelChecks, "ANNOTATION_CONSISTENCY_CSV", str(tmp_path / "out.csv"))
    return setup


def _id_model(rxn_id, met_id, gene_id):
    model = cobra.Model("ids")
    reaction = cobra.Reaction(rxn_id)
    reaction.add_metabolites({cobra.Metabolite(met_id, compartment=met_id[-1]): -1})
    reaction.gene_reaction_rule = gene_id
    model.add_reactions([reaction])
    return model


def _format_issues(issues):
    return sorted((kind, i) for kind, i, issue in issues if issue.startswith(qcModelChecks.ID_FORMAT_ISSUE))


def test_human_gem_identifiers_pass_the_format_check(tables):
    tables(["MAR00001"], ["MAM00001c"], ["ENSG00000000419"])
    issues = qcModelChecks.check_annotation_consistency(_id_model("MAR00001", "MAM00001c", "ENSG00000000419"))
    assert _format_issues(issues) == []


def test_six_digit_and_malformed_identifiers_are_flagged(tables):
    tables(["MAR020190"], ["MAM0001c"], ["ENSG0000000041"])
    issues = qcModelChecks.check_annotation_consistency(_id_model("MAR020190", "MAM0001c", "ENSG0000000041"))
    assert _format_issues(issues) == [
        ("gene", "ENSG0000000041"), ("metabolite", "MAM0001c"), ("reaction", "MAR020190")]


def test_growth_medium_opens_the_task_inputs_capped(tmp_path, monkeypatch):
    tasks = tmp_path / "tasks.txt"
    tasks.write_text(
        "\tID\tDESCRIPTION\tIN\tOUT\n"
        "\tER\tother task\tlactate[e]\t\n"
        "\tGR\tgrowth\tglucose[e]\tbiomass[e]\n"
        "\t\t\tO2[e]\t\n"
        "\t\t\tfructose[e]\t\n", encoding="utf-8")
    inputs = qcModelChecks._medium_inputs(str(tasks))
    assert inputs == ["glucose[e]", "O2[e]", "fructose[e]"]  # the GR task's inputs only
    monkeypatch.setattr(qcModelChecks, "_medium_inputs", lambda: inputs)
    model = cobra.Model("medium")
    for mid, name in (("glc_e", "glucose"), ("o2_e", "O2"), ("lac_e", "lactate")):
        exchange = cobra.Reaction(f"EX_{mid}", lower_bound=-1000, upper_bound=1000)
        exchange.add_metabolites({cobra.Metabolite(mid, name=name, compartment="e"): -1})
        model.add_reactions([exchange])
    with model:
        unresolved = qcModelChecks.apply_growth_medium(model)
        bounds = {r.id: r.lower_bound for r in model.reactions}
    assert unresolved == ["fructose[e]"]
    assert bounds["EX_lac_e"] == 0  # not a GR input: uptake closed
    expected = {"glucose": "EX_glc_e", "O2": "EX_o2_e"}
    for name, rid in expected.items():
        cap = 1000.0 if name in qcModelChecks.MEDIUM_UNCAPPED else qcModelChecks.MEDIUM_UPTAKE_CAP
        assert bounds[rid] == -cap
