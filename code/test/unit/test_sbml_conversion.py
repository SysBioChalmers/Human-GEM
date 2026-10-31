"""The SBML gate compares names, subsystems (SBML groups) and reaction notes."""

import cobra
from cobra.core import Group

import testSbmlConversion as ts


def _model(subsystem="Glycolysis", name="hexokinase", notes=None):
    model = cobra.Model("m")
    reaction = cobra.Reaction("R1", name=name)
    reaction.add_metabolites({cobra.Metabolite("A", compartment="c"): -1})
    reaction.subsystem = subsystem
    reaction.notes = dict(notes or {})
    model.add_reactions([reaction])
    return model


def _reloaded(group_names, name="hexokinase", notes=None):
    model = _model(subsystem="", name=name, notes=notes)
    model.add_groups([Group(f"g{i}", name=n, members=[model.reactions.R1]) for i, n in enumerate(group_names)])
    return model


def test_subsystems_read_back_from_groups():
    assert ts.compare_subsystems(_model(), _reloaded(["Glycolysis"])) == []
    assert ts.compare_subsystems(_model(), _reloaded([])) == ["reaction R1 subsystems ['Glycolysis'] -> []"]
    two = _model(subsystem=["Glycolysis", "Transport"])
    assert ts.compare_subsystems(two, _reloaded(["Transport", "Glycolysis"])) == []


def test_names_are_compared():
    assert ts.compare_names(_model(), _reloaded(["Glycolysis"], name="HK")) == [
        "reaction R1 name 'hexokinase' -> 'HK'"]


def test_notes_are_compared_under_their_sbml_labels():
    notes = {"references": "PMID:1", "confidence_score": 2.0}
    written = {"AUTHORS": "PMID:1", "Confidence Level": "2.0"}
    assert ts.compare_notes(_model(notes=notes), _reloaded([], notes=written)) == []
    assert ts.compare_notes(_model(notes=notes), _reloaded([], notes={"AUTHORS": "PMID:1"})) == [
        "reaction R1 note 'confidence_score' missing"]


def test_escaped_characters_in_read_back_notes_are_resolved():
    notes = {"note": "13'-carboxychromanols"}
    written = {"NOTES": "13&apos;-carboxychromanols"}
    assert ts.compare_notes(_model(notes=notes), _reloaded([], notes=written)) == []
