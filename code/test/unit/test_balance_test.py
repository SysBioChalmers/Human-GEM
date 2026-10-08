"""balanceTest skips the lumped reactions, and keeps R and X as elements."""

import cobra
import pytest

import balanceTest


def _rxn(rid, mets, subsystem=None):
    rxn = cobra.Reaction(rid, lower_bound=0, upper_bound=1000)
    rxn.add_metabolites(mets)
    if subsystem is not None:
        rxn.subsystem = subsystem
    return rxn


def _met(mid, formula, charge=0):
    m = cobra.Metabolite(mid, compartment="c")
    m.formula, m.charge = formula, charge
    return m


@pytest.mark.parametrize("subsystem", ["Pool reactions", "Artificial reactions",
                                       ["Pool reactions"], ["Artificial reactions", "x"]])
def test_lumped_subsystems_are_skipped(subsystem):
    rxn = _rxn("MAR00001", {}, subsystem)
    assert balanceTest._subsystems(rxn) & balanceTest.LUMPED_SUBSYSTEMS


@pytest.mark.parametrize("subsystem", ["Glycolysis", "", None, [], ["Glycolysis"]])
def test_other_subsystems_are_not_skipped(subsystem):
    rxn = _rxn("MAR00001", {}, subsystem)
    assert not balanceTest._subsystems(rxn) & balanceTest.LUMPED_SUBSYSTEMS


def test_subsystems_accepts_a_string_or_a_list():
    assert balanceTest._subsystems(_rxn("R", {}, "Pool reactions")) == {"Pool reactions"}
    assert balanceTest._subsystems(_rxn("R", {}, ["a", "b"])) == {"a", "b"}
    assert balanceTest._subsystems(_rxn("R", {}, None)) == set()


def test_R_and_X_are_balanced_as_elements():
    """A reaction that balances in C and H but not in R is still an imbalance."""
    a, b = _met("a_c", "C2H4R"), _met("b_c", "C2H4")
    rxn = _rxn("MAR00002", {a: -1, b: 1}, "Glycolysis")
    cobra.Model("m").add_reactions([rxn])
    assert rxn.check_mass_balance().get("R") == -1

    c, d = _met("c_c", "X"), _met("d_c", "X")
    ok = _rxn("MAR00003", {c: -1, d: 1}, "Glycolysis")
    cobra.Model("m2").add_reactions([ok])
    assert not ok.check_mass_balance()
