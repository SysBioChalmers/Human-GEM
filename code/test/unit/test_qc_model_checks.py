"""check_name_consistency flags missing names and names that disagree by compartment."""

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
