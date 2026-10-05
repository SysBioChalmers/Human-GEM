"""annotationTest flags rows with the wrong number of fields and malformed gene identifiers."""

import annotationTest


def _write(path, lines):
    path.write_text("\n".join("\t".join(fields) for fields in lines) + "\n", encoding="utf-8")
    return str(path)


def test_short_row_is_flagged(tmp_path):
    path = _write(tmp_path / "genes.tsv", [
        ["genes", "geneENSTID", "geneUniProtID"],
        ["ENSG00000000001", "ENST00000000001", "P10109"],
        ["ENSG00000000002", "ENST00000000002"],
    ])
    issues = []
    rows = annotationTest._read_table(path, "genes", issues)
    assert len(rows) == 2
    assert issues == [("ENSG00000000002", "genes", "malformed: row has 2 fields, header has 3")]


def test_versioned_ensembl_ids_are_malformed():
    rows = [{"genes": "ENSG00000135070", "geneENSTID": "ENST00000375991.9;ENST00000326094",
             "geneENSPID": "ENSP00000365159"}]
    issues = []
    annotationTest._check_format(rows, "genes", annotationTest.GENE_PATTERNS, issues)
    assert issues == [("ENSG00000135070", "geneENSTID", "malformed: ENST00000375991.9")]


def test_uniprot_and_entrez_formats():
    rows = [
        {"genes": "ENSG00000000001", "geneUniProtID": "P10109", "geneEntrezID": "2230"},
        {"genes": "ENSG00000000002", "geneUniProtID": "A0A024RBG1", "geneEntrezID": "81689"},
        {"genes": "ENSG00000000003", "geneUniProtID": "ISCA1", "geneEntrezID": "iron"},
    ]
    issues = []
    annotationTest._check_format(rows, "genes", annotationTest.GENE_PATTERNS, issues)
    assert sorted(issues) == [
        ("ENSG00000000003", "geneEntrezID", "malformed: iron"),
        ("ENSG00000000003", "geneUniProtID", "malformed: ISCA1"),
    ]
