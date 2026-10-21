"""Metrics, symbol matching and the gene universe of the Hart 2015 comparison.

Every expected value below is worked out by hand in the comment next to it.
"""

import math

import pytest

import evaluateHart2015Essentiality as ev

LINES = ["DLD1", "GBM"]


# --- metrics -------------------------------------------------------------------------

def test_auroc_counts_pairs_and_averages_ties():
    # positives score 0.9 and 0.5, negatives 0.5 and 0.1. Pairs (pos > neg): 0.9 beats
    # both, 0.5 ties 0.5 (1/2) and beats 0.1: (2 + 1.5) / 4.
    assert ev.auroc([(0.9, 1), (0.5, 1), (0.5, 0), (0.1, 0)]) == pytest.approx(0.875)


def test_auroc_needs_both_classes():
    assert math.isnan(ev.auroc([(0.1, 1), (0.2, 1)]))


def test_auprc_takes_tied_scores_together():
    # Scores 1.0 (pos, neg) then 0.5 (pos). At 1.0: recall 1/2, precision 1/2.
    # At 0.5: recall 1, precision 2/3. AP = 0.5 * 0.5 + 0.5 * 2/3.
    assert ev.auprc([(1.0, 1), (1.0, 0), (0.5, 1)]) == pytest.approx(0.25 + 1 / 3)


def test_auprc_does_not_depend_on_the_order_of_ties():
    a = ev.auprc([(1.0, 1), (1.0, 0), (0.5, 1)])
    b = ev.auprc([(1.0, 0), (1.0, 1), (0.5, 1)])
    assert a == b


def test_mcc():
    # TP 2, TN 3, FP 1, FN 1: (6 - 1) / sqrt(3 * 3 * 4 * 4) = 5 / 12.
    assert ev._mcc(dict(TP=2, TN=3, FP=1, FN=1)) == pytest.approx(5 / 12)
    assert math.isnan(ev._mcc(dict(TP=0, TN=3, FP=0, FN=0)))


# --- symbol matching -----------------------------------------------------------------

HGNC = {"RPF1": "ENSG_RPF1", "NEDD4": "ENSG_NEDD4", "GMPPA": "ENSG_SELF"}


def test_own_symbol_wins_over_aliases():
    names = ev.hart_names("ENSG_NEDD4", {"NEDD4"}, {"RPF1"}, {"NEDD4", "RPF1"}, HGNC)
    assert names == {"NEDD4"}


def test_alias_that_is_another_genes_symbol_is_not_used():
    # NEDD4's own symbol unscored: RPF1 is the current symbol of another gene.
    assert ev.hart_names("ENSG_NEDD4", {"NEDD4"}, {"RPF1"}, {"RPF1"}, HGNC) == set()


def test_alias_that_is_no_current_symbol_is_used():
    assert ev.hart_names("ENSG_X", {"NEWNAME"}, {"OLDNAME"}, {"OLDNAME"}, HGNC) == {"OLDNAME"}


def test_alias_that_is_the_genes_own_current_symbol_is_used():
    assert ev.hart_names("ENSG_SELF", {"OUTDATED"}, {"GMPPA"}, {"GMPPA"}, HGNC) == {"GMPPA"}


# --- the gene universe ---------------------------------------------------------------

@pytest.fixture
def hart(monkeypatch):
    """Four genes scored by Hart in both lines: g1, g2 fitness genes, g3, g4 not.
    Thresholds: DLD1 3.57, GBM 3.20."""
    bf = {"g1": 10.0, "g2": 10.0, "g3": -5.0, "g4": -5.0}
    monkeypatch.setattr(ev, "bayes_factors", lambda genes, tissues, **kw: {t: dict(bf) for t in tissues})
    return bf


def _row(rows, line):
    return next(r for r in rows if r["cellLine"] == line)


def test_gene_outside_the_context_model_is_scored_as_non_essential(hart):
    # g2 is not in either context model: it counts as a false negative, with growth 1.0.
    record = {
        "g1": {"tasks": ["GR"], "growth": 0.0},
        "g3": {"tasks": [], "growth": 1.0},
        "g4": {"tasks": ["BS"], "growth": 1.0},
    }
    per_tissue = {t: record for t in LINES}
    rows, matrix = ev.evaluate_graded(per_tissue, LINES, ["g1", "g2", "g3", "g4"],
                                      symbol_of={"g2": "G2"})
    dld1 = _row(rows, "DLD1")
    assert dld1["genes"] == 4
    # all tasks: TP g1, FN g2, TN g3, FP g4 -> MCC 0.
    assert dld1["allTaskFP"] == 1 and dld1["allTaskMCC"] == 0.0
    # viability: TP g1, FN g2, TN g3, TN g4 -> (1*2 - 0*1) / sqrt(1 * 2 * 2 * 3) = 2 / sqrt(12).
    assert dld1["viabilityMCC"] == round(2 / math.sqrt(12), 4)
    # growth scores 1 - ratio: g1 1.0 (pos), g2 0.0 (pos), g3 0.0, g4 0.0 (neg).
    # AUROC: g1 beats both negatives (2), g2 ties both (1): 3 / 4.
    assert dld1["growthAUROC"] == 0.75
    assert dld1["baseRate"] == 0.5
    assert dld1["capabilityOnly"] == 1
    # The matrix marks the absent gene with "." tasks and growth, and the class it is scored with.
    g2_row = next(line for line in matrix.splitlines() if line.startswith("g2,"))
    assert g2_row == "g2,G2,FN,.,.,FN,.,."


def test_all_row_needs_the_call_in_every_line(hart):
    # g1 essential in both lines (TP); g3 essential only in DLD1 (not essential in all: TN).
    per_tissue = {
        "DLD1": {"g1": {"tasks": ["GR"], "growth": 0.0}, "g3": {"tasks": ["GR"], "growth": 0.0}},
        "GBM": {"g1": {"tasks": ["GR"], "growth": 0.0}},
    }
    rows, _ = ev.evaluate_graded(per_tissue, LINES, ["g1", "g2", "g3", "g4"], symbol_of={})
    every = _row(rows, "all")
    assert every["genes"] == 4
    assert every["allTaskFP"] == 0
    # TP g1, FN g2, TN g3 and g4: (1*2 - 0) / sqrt(1 * 2 * 2 * 3).
    assert every["allTaskMCC"] == round(2 / math.sqrt(12), 4)


def test_summary_states_the_universe(hart):
    rows, _ = ev.evaluate_graded({t: {} for t in LINES}, LINES, ["g1"], symbol_of={})
    text = ev.summary_to_markdown(rows)
    assert "outside the cell line's context model counts as non-essential" in text
    assert "| cellLine | genes |" in text
