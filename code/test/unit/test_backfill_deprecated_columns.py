"""backfillDeprecatedColumns on small deprecated lists."""

import backfillDeprecatedColumns as bdc


def test_backfill_adds_columns_and_fills_empty_cells(tmp_path):
    (tmp_path / "model").mkdir()
    (tmp_path / "data" / "deprecatedIdentifiers").mkdir(parents=True)
    (tmp_path / "model" / "reactions.tsv").write_text("rxns\trxnMetaNetXID\trxnSeedID\n", encoding="utf-8")
    (tmp_path / "model" / "metabolites.tsv").write_text("mets\tmetSmiles\n", encoding="utf-8")
    dep = tmp_path / "data" / "deprecatedIdentifiers"
    (dep / "deprecatedReactions.tsv").write_text("rxns\trxnMetaNetXID\nMAR1\t\nMAR2\tMNXR2\n", encoding="utf-8")
    (dep / "deprecatedMetabolites.tsv").write_text("mets\nMAM1c\n", encoding="utf-8")
    rows = [{"list": "reactions", "id": "MAR1", "column": "rxnSeedID", "value": "rxn1"},
            {"list": "reactions", "id": "MAR1", "column": "rxnMetaNetXID", "value": "MNXR1"},
            {"list": "metabolites", "id": "MAM9c", "column": "metSmiles", "value": "C"}]
    added, filled, problems, lists = bdc.backfill(rows, tmp_path)
    for table in lists.values():
        table.save()
    assert problems == [] and added == {"reactions": ["rxnSeedID"], "metabolites": ["metSmiles"]}
    assert (dep / "deprecatedReactions.tsv").read_text(encoding="utf-8") == (
        "rxns\trxnMetaNetXID\trxnSeedID\nMAR1\tMNXR1\trxn1\nMAR2\tMNXR2\t\n")
    assert (dep / "deprecatedMetabolites.tsv").read_text(encoding="utf-8") == "mets\tmetSmiles\nMAM1c\t\n"
