"""buildReport's counting and delta helpers, which turn the result files into rows."""

import buildReport


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def test_count_csv_counts_rows_not_lines(tmp_path):
    path = _write(tmp_path / "c.csv", "kind,id\nmetabolite,MAM1\ngene,ENSG1\n")
    assert buildReport._count_csv(path) == 2


def test_count_csv_applies_the_predicate(tmp_path):
    path = _write(tmp_path / "c.csv", "kind,id\nmetabolite,MAM1\ngene,ENSG1\n")
    assert buildReport._count_csv(path, lambda row: row["kind"] == "gene") == 1


def test_count_csv_of_a_header_only_file_is_zero(tmp_path):
    assert buildReport._count_csv(_write(tmp_path / "c.csv", "kind,id\n")) == 0


def test_count_csv_of_a_missing_file_is_none(tmp_path):
    """None means "not measured", which the report shows as pending rather than 0."""
    assert buildReport._count_csv(tmp_path / "absent.csv") is None


def test_distinct_csv_counts_unique_values(tmp_path):
    path = _write(tmp_path / "d.csv", "group,reaction\n1,MAR1\n1,MAR2\n2,MAR3\n")
    assert buildReport._distinct_csv(path, "group") == 2


def test_status_map_skips_the_header_and_blank_lines(tmp_path):
    _write(tmp_path / "qc_status.tsv", "check\tresult\ngrowth\t12.5\n\nyamllint\tpass\n")
    assert buildReport._status_map(tmp_path) == {"growth": "12.5", "yamllint": "pass"}


def test_growth_reads_the_value_from_the_status_file(tmp_path):
    _write(tmp_path / "qc_status.tsv", "check\tresult\ngrowth\t12.5\n")
    assert buildReport._growth(tmp_path) == 12.5


def test_growth_is_none_when_unset_or_unparsable(tmp_path):
    assert buildReport._growth(tmp_path) is None
    _write(tmp_path / "qc_status.tsv", "check\tresult\ngrowth\tnot-a-number\n")
    assert buildReport._growth(tmp_path) is None


def test_a_rising_count_is_a_regression():
    delta, icon, regression, fatal = buildReport._icon(3, 1, "count")
    assert delta == "+2"
    assert icon == ":x:"
    assert regression is True
    assert fatal is False


def test_a_falling_count_is_not_a_regression():
    delta, icon, regression, _fatal = buildReport._icon(1, 3, "count")
    assert delta == "-2"
    assert regression is False
    assert icon == ":warning:"          # non-zero count still warns


def test_an_unchanged_zero_count_is_clean():
    delta, icon, regression, _fatal = buildReport._icon(0, 0, "count")
    assert delta == "0"
    assert icon == ":white_check_mark:"
    assert regression is False


def test_no_growth_is_fatal():
    _delta, icon, _regression, fatal = buildReport._icon(0.0, 1.0, "growth")
    assert icon == ":x:"
    assert fatal is True


def test_growth_is_not_fatal():
    _delta, icon, _regression, fatal = buildReport._icon(12.5, 12.5, "growth")
    assert icon == ":white_check_mark:"
    assert fatal is False


def test_a_missing_base_reads_as_new():
    delta, _icon, _regression, _fatal = buildReport._icon(2, None, "count")
    assert delta == "new"
