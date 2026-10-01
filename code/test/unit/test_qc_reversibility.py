"""reversibilityTest.check_reversibility flags thermodynamically very unlikely directions."""

import importlib.util
from pathlib import Path

import cobra
import pytest

import reversibilityTest

HEADER = "reaction\tcompartments\ttype\tdGm_kJ_per_mol\tsd_kJ_per_mol\tstoichiometry_hash\n"


@pytest.fixture
def files(tmp_path, monkeypatch):
    monkeypatch.setattr(reversibilityTest, "REVERSIBILITY_CSV", str(tmp_path / "qc_reversibility.csv"))
    return tmp_path


def _model(lb, ub, two_by_two=False):
    """A -> B, or A + C -> B + D with two_by_two."""
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00001", lower_bound=lb, upper_bound=ub)
    mets = {cobra.Metabolite("MAM00001c", compartment="c"): -1,
            cobra.Metabolite("MAM00002c", compartment="c"): 1}
    if two_by_two:
        mets.update({cobra.Metabolite("MAM00003c", compartment="c"): -1,
                     cobra.Metabolite("MAM00004c", compartment="c"): 1})
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    return model


def _check(files, model, dg, sd, exceptions="", stoich_hash=None, kind="single"):
    rxn = model.reactions.MAR00001
    h = stoich_hash or reversibilityTest.stoichiometry_hash({m.id: c for m, c in rxn.metabolites.items()})
    (files / "dg.tsv").write_text(HEADER + f"MAR00001\tc\t{kind}\t{dg}\t{sd}\t{h}\n")
    (files / "exc.tsv").write_text("reaction\treason\n" + exceptions)
    return reversibilityTest.check_reversibility(model, str(files / "dg.tsv"), str(files / "exc.tsv"))


def _verdicts(rows):
    return [r[7] for r in rows]


def test_reversible_reaction_with_large_dg_is_impossible(files):
    assert _verdicts(_check(files, _model(-1000, 1000), -60, 5)) == ["impossible"]


def test_irreversible_in_the_downhill_direction_passes(files):
    assert _check(files, _model(0, 1000), -60, 5) == []
    assert _check(files, _model(-1000, 0), 60, 5) == []


def test_irreversible_in_the_uphill_direction_is_impossible(files):
    assert _verdicts(_check(files, _model(0, 1000), 60, 5)) == ["impossible"]


def test_uncertainty_is_subtracted(files):
    # 60 - 2 * 15 = 30: questionable, not impossible
    assert _verdicts(_check(files, _model(-1000, 1000), -60, 15)) == ["questionable"]
    assert _check(files, _model(-1000, 1000), -60, 25) == []


def test_warning_is_normalised_for_stoichiometry(files):
    # a 25 kJ/mol margin needs a ~24000-fold change for A -> B, but ~150-fold for A + C -> B + D
    assert _verdicts(_check(files, _model(-1000, 1000), -27, 1)) == ["questionable"]
    assert _check(files, _model(-1000, 1000, two_by_two=True), -27, 1) == []


def test_many_reactants_need_more_than_40_kj_to_alarm(files):
    # A + 4 C -> B + 4 D (N = 10): 50 kJ/mol is a 2.2-fold change per reactant, no flag
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00001", lower_bound=-1000, upper_bound=1000)
    rxn.add_metabolites({cobra.Metabolite("MAM00001c", compartment="c"): -1,
                         cobra.Metabolite("MAM00003c", compartment="c"): -4,
                         cobra.Metabolite("MAM00002c", compartment="c"): 1,
                         cobra.Metabolite("MAM00004c", compartment="c"): 4})
    model.add_reactions([rxn])
    assert _check(files, model, -52, 1) == []
    assert _verdicts(_check(files, model, -100, 1)) == ["impossible"]


def test_water_and_protons_do_not_count_for_the_index():
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00005")
    rxn.add_metabolites({cobra.Metabolite("MAM00001c", compartment="c"): -1,
                         cobra.Metabolite("MAM02040c", compartment="c"): -1,
                         cobra.Metabolite("MAM02039c", compartment="c"): 2,
                         cobra.Metabolite("MAM00002c", compartment="c"): 1})
    model.add_reactions([rxn])
    assert reversibilityTest._abs_sum_coefficients(rxn) == 2


def test_proxy_estimate_is_at_most_questionable(files):
    assert _verdicts(_check(files, _model(-1000, 1000), -90, 5, kind="single+proxy")) == ["questionable"]


def test_listed_exception_does_not_fail(files):
    rows = _check(files, _model(-1000, 1000), -60, 5, exceptions="MAR00001\treverse electron transfer\n")
    assert _verdicts(rows) == ["exception"]
    assert "reverse electron transfer" in rows[0][8]


def test_changed_stoichiometry_is_outdated(files):
    assert _verdicts(_check(files, _model(-1000, 1000), -60, 5, stoich_hash="0000000000")) == ["outdated"]


def test_estimator_hashes_the_same_way():
    path = Path(reversibilityTest.__file__).resolve().parents[1] / "qc" / "estimateReactionDeltaG.py"
    spec = importlib.util.spec_from_file_location("estimateReactionDeltaG", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stoich = {"MAM00002c": 1.0, "MAM00001c": -1.0, "MAM00003c": -2.5}
    assert module.stoichiometry_hash(stoich) == reversibilityTest.stoichiometry_hash(stoich)


def _o2_model(lb, ub, extra=None):
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00002", lower_bound=lb, upper_bound=ub)
    mets = {cobra.Metabolite("MAM00001c", compartment="c"): -1,
            cobra.Metabolite("MAM02630c", compartment="c"): -1,  # O2
            cobra.Metabolite("MAM00002c", compartment="c"): 1}
    mets.update(extra or {})
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    return model


def test_oxygenase_running_backward_is_impossible(files):
    rows = reversibilityTest.check_reversibility(_o2_model(-1000, 1000), str(files / "none.tsv"), str(files / "none.tsv"))
    assert _verdicts(rows) == ["impossible"]
    assert "release O2" in rows[0][8]


def test_irreversible_oxygenase_passes(files):
    assert reversibilityTest.check_reversibility(_o2_model(0, 1000), str(files / "none.tsv"), str(files / "none.tsv")) == []


def test_catalase_releasing_o2_passes(files):
    # 2 H2O2 -> O2 + 2 H2O, written backward and allowed to run that way
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00003", lower_bound=-1000, upper_bound=0)
    rxn.add_metabolites({cobra.Metabolite("MAM02630c", compartment="c"): -1,
                         cobra.Metabolite("MAM02040c", compartment="c"): -2,
                         cobra.Metabolite("MAM02041c", compartment="c"): 2})
    model.add_reactions([rxn])
    assert reversibilityTest.check_reversibility(model, str(files / "none.tsv"), str(files / "none.tsv")) == []


def test_o2_transport_is_ignored(files):
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00004", lower_bound=-1000, upper_bound=1000)
    rxn.add_metabolites({cobra.Metabolite("MAM02630c", compartment="c"): -1,
                         cobra.Metabolite("MAM02630m", compartment="m"): 1})
    model.add_reactions([rxn])
    assert reversibilityTest.check_reversibility(model, str(files / "none.tsv"), str(files / "none.tsv")) == []
