"""qcStatus upserts one key at a time into a shared file."""

import qcStatus


def _use_temp_status(tmp_path, monkeypatch):
    path = tmp_path / "qc_status.tsv"
    monkeypatch.setattr(qcStatus, "STATUS_FILE", path)
    return path


def test_missing_file_reads_as_empty(tmp_path, monkeypatch):
    path = _use_temp_status(tmp_path, monkeypatch)
    assert qcStatus.read_status(path) == {}
    assert qcStatus.get_status("growth") == ""


def test_set_then_get(tmp_path, monkeypatch):
    _use_temp_status(tmp_path, monkeypatch)
    qcStatus.set_status("growth", "123.4")
    assert qcStatus.get_status("growth") == "123.4"


def test_second_key_keeps_the_first(tmp_path, monkeypatch):
    _use_temp_status(tmp_path, monkeypatch)
    qcStatus.set_status("roundtrip_cobra", "pass")
    qcStatus.set_status("roundtrip_raven", "pass")
    assert qcStatus.read_status() == {"roundtrip_cobra": "pass", "roundtrip_raven": "pass"}


def test_overwriting_a_key_does_not_duplicate_it(tmp_path, monkeypatch):
    path = _use_temp_status(tmp_path, monkeypatch)
    qcStatus.set_status("yamllint", "fail")
    qcStatus.set_status("yamllint", "pass")
    assert qcStatus.get_status("yamllint") == "pass"
    assert sum(1 for line in path.read_text().splitlines() if line.startswith("yamllint")) == 1


def test_keys_are_sorted_and_the_header_is_written_once(tmp_path, monkeypatch):
    path = _use_temp_status(tmp_path, monkeypatch)
    for key in ("yamllint", "growth", "roundtrip_sbml"):
        qcStatus.set_status(key, "pass")
    lines = path.read_text().splitlines()
    assert lines[0] == "check\tresult"
    keys = [line.split("\t")[0] for line in lines[1:]]
    assert keys == sorted(keys)


def test_cli_get_and_set(tmp_path, monkeypatch, capsys):
    _use_temp_status(tmp_path, monkeypatch)
    assert qcStatus.main(["tasks_essential", "0/57"]) == 0
    assert qcStatus.main(["--get", "tasks_essential"]) == 0
    assert capsys.readouterr().out.strip() == "0/57"


def test_cli_rejects_wrong_argument_count(tmp_path, monkeypatch):
    _use_temp_status(tmp_path, monkeypatch)
    assert qcStatus.main(["only-one-argument"]) == 2
