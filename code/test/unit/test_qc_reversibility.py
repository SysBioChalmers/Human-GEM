"""check_reversibility judges reaction bounds against stored ΔG'm estimates."""

import importlib.util
from pathlib import Path

import cobra
import pytest

import qcModelChecks

HEADER = "reaction\tcompartments\ttype\tdGm_kJ_per_mol\tsd_kJ_per_mol\tstoichiometry_hash\n"


@pytest.fixture
def files(tmp_path, monkeypatch):
    monkeypatch.setattr(qcModelChecks, "REVERSIBILITY_CSV", str(tmp_path / "qc_reversibility.csv"))
    return tmp_path


def _model(lb, ub):
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00001", lower_bound=lb, upper_bound=ub)
    rxn.add_metabolites({cobra.Metabolite("MAM00001c", compartment="c"): -1,
                         cobra.Metabolite("MAM00002c", compartment="c"): 1})
    model.add_reactions([rxn])
    return model


def _check(files, model, dg, sd, exceptions="", stoich_hash=None, kind="single"):
    rxn = model.reactions.MAR00001
    h = stoich_hash or qcModelChecks.stoichiometry_hash({m.id: c for m, c in rxn.metabolites.items()})
    (files / "dg.tsv").write_text(HEADER + f"MAR00001\tc\t{kind}\t{dg}\t{sd}\t{h}\n")
    (files / "exc.tsv").write_text("reaction\treason\n" + exceptions)
    return qcModelChecks.check_reversibility(model, str(files / "dg.tsv"), str(files / "exc.tsv"))


def _verdicts(rows):
    return [r[6] for r in rows]


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


def test_proxy_estimate_is_at_most_questionable(files):
    assert _verdicts(_check(files, _model(-1000, 1000), -90, 5, kind="single+proxy")) == ["questionable"]


def test_listed_exception_does_not_fail(files):
    rows = _check(files, _model(-1000, 1000), -60, 5, exceptions="MAR00001\treverse electron transfer\n")
    assert _verdicts(rows) == ["exception"]
    assert "reverse electron transfer" in rows[0][7]


def test_changed_stoichiometry_is_outdated(files):
    assert _verdicts(_check(files, _model(-1000, 1000), -60, 5, stoich_hash="0000000000")) == ["outdated"]


def test_estimator_hashes_the_same_way():
    path = Path(qcModelChecks.__file__).resolve().parents[1] / "qc" / "estimateReactionDeltaG.py"
    spec = importlib.util.spec_from_file_location("estimateReactionDeltaG", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stoich = {"MAM00002c": 1.0, "MAM00001c": -1.0, "MAM00003c": -2.5}
    assert module.stoichiometry_hash(stoich) == qcModelChecks.stoichiometry_hash(stoich)


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
    rows = qcModelChecks.check_reversibility(_o2_model(-1000, 1000), str(files / "none.tsv"), str(files / "none.tsv"))
    assert _verdicts(rows) == ["impossible"]
    assert "release O2" in rows[0][7]


def test_irreversible_oxygenase_passes(files):
    assert qcModelChecks.check_reversibility(_o2_model(0, 1000), str(files / "none.tsv"), str(files / "none.tsv")) == []


def test_catalase_releasing_o2_passes(files):
    # 2 H2O2 -> O2 + 2 H2O, written backward and allowed to run that way
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00003", lower_bound=-1000, upper_bound=0)
    rxn.add_metabolites({cobra.Metabolite("MAM02630c", compartment="c"): -1,
                         cobra.Metabolite("MAM02040c", compartment="c"): -2,
                         cobra.Metabolite("MAM02041c", compartment="c"): 2})
    model.add_reactions([rxn])
    assert qcModelChecks.check_reversibility(model, str(files / "none.tsv"), str(files / "none.tsv")) == []


def test_o2_transport_is_ignored(files):
    model = cobra.Model("test")
    rxn = cobra.Reaction("MAR00004", lower_bound=-1000, upper_bound=1000)
    rxn.add_metabolites({cobra.Metabolite("MAM02630c", compartment="c"): -1,
                         cobra.Metabolite("MAM02630m", compartment="m"): 1})
    model.add_reactions([rxn])
    assert qcModelChecks.check_reversibility(model, str(files / "none.tsv"), str(files / "none.tsv")) == []
