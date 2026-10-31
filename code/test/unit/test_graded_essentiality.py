"""The guards, the reaction-order runs and the growth medium of gradedEssentiality.py."""

import sys
import types

import cobra
import pytest
from cobra.io import model_to_dict
from raven_toolbox.tasks.tasklist import Task

import gradedEssentiality as ge


def _toy_model():
    """glucose[e] <-> glucose[c] -> biomass; gene g1 on the uptake, g2 on the biomass."""
    model = cobra.Model("toy")
    glc_e = cobra.Metabolite("glc_e", name="glucose", compartment="e")
    glc_c = cobra.Metabolite("glc_c", name="glucose", compartment="c")
    ex = cobra.Reaction("EX_glc", lower_bound=-1000, upper_bound=1000)
    ex.add_metabolites({glc_e: -1})
    up = cobra.Reaction("UP", lower_bound=0, upper_bound=1000)
    up.add_metabolites({glc_e: -1, glc_c: 1})
    up.gene_reaction_rule = "g1"
    bio = cobra.Reaction("MAR13082", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({glc_c: -1})
    bio.gene_reaction_rule = "g2 or g3"
    model.add_reactions([ex, up, bio])
    model.objective = "MAR13082"
    return model


def test_missing_tasks_are_identified_by_id_and_description():
    tasks = [Task(id="BS", description="a"), Task(id="BS", description="b"), Task(id="GR", description="growth")]
    feasible = [Task(id="BS", description="a"), Task(id="BS", description="b")]
    assert [t.description for t in ge.missing_tasks(tasks, feasible)] == ["growth"]
    assert ge.missing_tasks(tasks, tasks) == []


def test_growth_ratios_on_the_task_medium():
    task = Task(id="GR", description="growth", inputs=[("glucose[e]", 0, 1000)])
    ratios, wild_type, n_media, unresolved = ge.growth_ratios(_toy_model(), task)
    assert wild_type == pytest.approx(1000)
    assert n_media == 1 and unresolved == []
    assert ratios["g1"] == 0.0          # the only uptake
    assert ratios["g2"] == pytest.approx(1.0)  # an isozyme takes over


def test_no_growth_gives_no_ratios():
    task = Task(id="GR", description="growth", inputs=[("oxygen[e]", 0, 1000)])
    ratios, wild_type, _n, unresolved = ge.growth_ratios(_toy_model(), task)
    assert ratios == {} and wild_type == 0.0
    assert unresolved == ["oxygen[e]"]


def test_reorder_model_keeps_the_model():
    model = _toy_model()
    assert ge.reorder_model(model, 0) is model
    shuffled = ge.reorder_model(model, 3)
    assert sorted(r.id for r in shuffled.reactions) == sorted(r.id for r in model.reactions)
    for rxn in model.reactions:
        other = shuffled.reactions.get_by_id(rxn.id)
        assert {m.id: c for m, c in other.metabolites.items()} == {m.id: c for m, c in rxn.metabolites.items()}
        assert other.gene_reaction_rule == rxn.gene_reaction_rule
        assert other.bounds == rxn.bounds
    assert shuffled.slim_optimize() == pytest.approx(model.slim_optimize())


def _long_model(n=12):
    model = cobra.Model("long")
    for i in range(n):
        rxn = cobra.Reaction(f"R{i:02d}", lower_bound=0, upper_bound=1000)
        rxn.add_metabolites({cobra.Metabolite(f"m{i:02d}_c", compartment="c"): -1})
        rxn.gene_reaction_rule = f"g{i:02d}"
        model.add_reactions([rxn])
    return model


def test_reorder_model_changes_the_order_and_nothing_else():
    model = _long_model()
    shuffled = ge.reorder_model(model, 1)
    assert [r.id for r in shuffled.reactions] != [r.id for r in model.reactions]
    assert [m.id for m in shuffled.metabolites] != [m.id for m in model.metabolites]

    def by_id(data):
        return {key: sorted(data[key], key=lambda entry: entry["id"]) if isinstance(data[key], list) else data[key]
                for key in data}

    assert by_id(model_to_dict(shuffled, sort=False)) == by_id(model_to_dict(model, sort=False))


def test_reorder_model_is_reproducible():
    a = [r.id for r in ge.reorder_model(_toy_model(), 7).reactions]
    b = [r.id for r in ge.reorder_model(_toy_model(), 7).reactions]
    assert a == b


def test_unstable_cells_compare_any_task_and_viability_calls():
    run0 = {"DLD1": {"g1": {"tasks": ["GR"]}, "g2": {"tasks": ["BS"]}, "g3": {"tasks": []}}}
    run1 = {"DLD1": {"g1": {"tasks": ["ER"]}, "g2": {"tasks": ["GR"]}}}
    cells = ge.unstable_cells([run0, run1], ["DLD1"], ["g1", "g2", "g3"])
    # g1: viability in both (GR vs ER) -> stable; g2: capability-only vs viability ->
    # unstable; g3: in the model in run 0, absent in run 1, non-essential both -> stable.
    assert cells == [("g2", "DLD1", ["BS", "GR"])]
    text = ge.unstable_csv(cells, 2, {"g2": "G2"})
    assert text.splitlines() == ["genes,geneSymbol,cellLine,order0_tasks,order1_tasks", "g2,G2,DLD1,BS,GR"]


def test_stability_table_gives_the_mcc_range():
    rows0 = [{"cellLine": "DLD1", "allTaskMCC": 0.30, "viabilityMCC": 0.35}, {"cellLine": "all"}]
    rows1 = [{"cellLine": "DLD1", "allTaskMCC": 0.32, "viabilityMCC": 0.35}, {"cellLine": "all"}]
    text = ge.stability_markdown([rows0, rows1], [("g2", "DLD1", ["BS", "GR"])], ["DLD1"])
    assert "| DLD1 | 0.3000-0.3200 | 0.3500 | 1 |" in text


def test_environment_section_names_versions_and_warned_lines():
    metas = {
        ("d0", "DLD1"): {"environment": {"raven-toolbox": "3.0.0b1"}, "ftinit_warnings": []},
        ("d0", "RPE1"): {"environment": {"raven-toolbox": "3.0.0b1"}, "reaction_order": 0,
                         "ftinit_warnings": ["ftINIT step 0 hit the 900.0s time limit"]},
    }
    text = ge.environment_markdown(metas, ["DLD1", "RPE1"])
    assert "- raven-toolbox 3.0.0b1" in text
    assert "- RPE1 (reaction order 0): 1 warning(s)" in text
    assert "DLD1 (" not in text


def test_run_environment_lists_cobra():
    assert ge.run_environment()["cobra"] == cobra.__version__


GROWTH = Task(id="GR", description="growth", inputs=[("glucose[e]", 0, 1000)])
OTHER = Task(id="BS", description="a capability")


def _stub_run(monkeypatch, tmp_path, *, feasible, context_task=GROWTH):
    """Patch everything main() reads or solves with Gurobi, for a one-line scan."""
    monkeypatch.setitem(sys.modules, "gurobipy", types.SimpleNamespace(setParam=lambda *a: None))
    monkeypatch.setattr(ge, "parse_task_list", lambda path: [GROWTH, OTHER])
    monkeypatch.setattr(ge, "_read_rnaseq", lambda path: (["L1"], {"L1": {}}))
    monkeypatch.setattr(ge, "read_yaml_model", lambda path: _toy_model())
    monkeypatch.setattr(ge, "_use_gurobi_solver", lambda model: None)
    monkeypatch.setattr(ge, "_dispose_gurobi_session", lambda: None)
    monkeypatch.setattr(ge, "_load_or_build_prep", lambda model, tasks, cache: types.SimpleNamespace(tasks=feasible))
    monkeypatch.setattr(ge, "_build_logged", lambda prep, model, expression: _toy_model())
    monkeypatch.setattr(ge, "find_task_essential_categories", lambda *a, **k: {})
    if context_task is not GROWTH:
        monkeypatch.setattr(ge, "_growth_task", lambda feasible: context_task)
    return ["--checkpoint-dir", str(tmp_path / "ckpt"), "--processes", "1"]


def test_main_fails_when_step1_drops_an_essential_task(monkeypatch, tmp_path):
    argv = _stub_run(monkeypatch, tmp_path, feasible=[OTHER])  # GR dropped
    assert ge.main(argv) == 1
    assert not list((tmp_path / "ckpt").glob("graded_*.json"))


def test_main_fails_without_writing_a_checkpoint_when_the_context_model_cannot_grow(monkeypatch, tmp_path):
    no_medium = Task(id="GR", description="growth", inputs=[("oxygen[e]", 0, 1000)])
    argv = _stub_run(monkeypatch, tmp_path, feasible=[GROWTH, OTHER], context_task=no_medium)
    assert ge.main(argv) == 1
    assert not (tmp_path / "ckpt" / "graded_L1.json").exists()
