"""A gene outside a context model, or outside the model, counts as non-essential there."""

import cobra
import pytest

import diffGeneEssentiality

LINES = ["DLD1", "GBM", "HCT116", "HELA", "RPE1"]


def _calls(tasks, growth=None):
    return {t: {"tasks": frozenset(tasks), "growth": growth} for t in LINES}


def _diff(base, head, hart=None):
    return diffGeneEssentiality.diff_gene("ENSG0", LINES, base, head, hart or {}, growth_tolerance=0.01)


def test_gene_new_to_the_model_that_blocks_growth_is_a_consensus_gain():
    result = _diff({}, _calls({"GR"}))
    assert result["viability_verdict"] == "consensus gained"
    assert result["viability_count"] == "5/5"


def test_gene_that_enters_the_context_models_is_scored():
    base = {t: None for t in LINES}  # in the model, but in no context model
    result = _diff(base, _calls({"GR"}))
    assert result["viability_verdict"] == "consensus gained"
    assert result["membership_changed_lines"] == LINES


def test_removed_essential_gene_is_a_consensus_loss():
    result = _diff(_calls({"GR"}), {})
    assert result["viability_verdict"] == "consensus lost"


def test_lines_outside_the_context_model_on_both_sides_are_not_counted():
    base = {t: None for t in LINES}
    head = {**base, "DLD1": {"tasks": frozenset({"GR"}), "growth": 0.0}}
    result = _diff(base, head)
    assert result["viability_verdict"] == "isolated (likely noise)"
    assert result["viability_count"] == "1/1"


def test_removed_gene_is_scored_against_hart():
    assert _diff(_calls({"GR"}), {}, {t: False for t in LINES})["hart_verdict"] == "improvement"
    assert _diff(_calls({"GR"}), {}, {t: True for t in LINES})["hart_verdict"] == "regression"


def test_growth_on_one_side_only_is_not_a_growth_shift():
    result = _diff(_calls({"GR"}, growth=0.0), {})
    assert result["growth_shift_lines"] == []


def test_removed_gene_with_a_capability_role_only_is_a_capability_loss():
    result = _diff(_calls({"BS"}, growth=1.0), {})
    assert result["any_verdict"] == "consensus lost"
    assert result["viability_verdict"] == "stable"


def test_lines_leaving_the_context_model_add_to_real_flips():
    head = {
        "DLD1": None,
        "GBM": None,
        "HCT116": {"tasks": frozenset(), "growth": 1.0},
        "HELA": {"tasks": frozenset({"GR"}), "growth": 0.0},
        "RPE1": {"tasks": frozenset({"GR"}), "growth": 0.0},
    }
    result = _diff(_calls({"GR"}, growth=0.0), head)
    assert result["viability_verdict"] == "consensus lost"
    assert result["viability_count"] == "3/5"
    assert result["membership_changed_lines"] == ["DLD1", "GBM"]


def _model(gene_ids):
    """One reaction A -> B per gene, so every gene neighbours every changed reaction."""
    model = cobra.Model("m")
    a, b = cobra.Metabolite("A"), cobra.Metabolite("B")
    for gene_id in gene_ids:
        rxn = cobra.Reaction(f"R_{gene_id}")
        rxn.add_metabolites({a: -1, b: 1})
        rxn.gene_reaction_rule = gene_id
        model.add_reactions([rxn])
    return model


def _matrix(path, rows):
    header = "genes,geneSymbol," + ",".join(f"{t}_class,{t}_tasks,{t}_growth" for t in LINES)
    lines = [header] + [f"{g},{g}," + ",".join(f"TN,{tasks},1.000000" for _ in LINES) for g, tasks in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def report(tmp_path, monkeypatch):
    """build_report on two in-memory models, without Hart 2015 calls."""
    models = {}
    monkeypatch.setattr(diffGeneEssentiality, "read_yaml_model", lambda path: models[str(path)])
    monkeypatch.setattr(diffGeneEssentiality, "bayes_factors", lambda *args, **kwargs: {})

    def run(base_genes, head_genes, base_rows, head_rows):
        models["base"], models["head"] = _model(base_genes), _model(head_genes)
        base_matrix = _matrix(tmp_path / "base.csv", base_rows)
        head_matrix = _matrix(tmp_path / "head.csv", head_rows)
        return diffGeneEssentiality.build_report("base", "head", base_matrix, head_matrix)

    return run


def test_genes_removed_from_or_new_to_the_model_are_labelled(report):
    _summary, detail, detail_csv = report(
        ["G1", "G2"], ["G1", "G3"], [("G1", ""), ("G2", "GR")], [("G1", ""), ("G3", "BS")],
    )
    assert "G2: 5/5 lines, removed from the model, knockout blocked growth" in detail
    assert "G3: 5/5 lines, biosynthesis -- new to the model, required" in detail
    assert ",removed," in detail_csv
    assert ",added," in detail_csv


def test_gene_missing_from_a_matrix_older_than_its_model_is_not_compared(report):
    _summary, detail, _csv = report(
        ["G1", "G2"], ["G1", "G2", "G3"], [("G1", "")], [("G1", ""), ("G2", "GR"), ("G3", "")],
    )
    assert "**Not compared** (missing from a gene-essentiality matrix older than its model): G2." in detail
    assert "removed from the model" not in detail


def test_matrices_without_a_shared_cell_line_are_refused(report, tmp_path):
    other = tmp_path / "other.csv"
    other.write_text("genes,geneSymbol,X_class,X_tasks,X_growth\nG1,G1,TN,,1.000000\n", encoding="utf-8")
    report(["G1"], ["G1"], [("G1", "")], [("G1", "")])
    with pytest.raises(ValueError, match="share no cell line"):
        diffGeneEssentiality.build_report("base", "head", tmp_path / "base.csv", other)
