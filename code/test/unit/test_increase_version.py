"""The release stamps its version into the package and citation files."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "io"))
import increaseHumanGEMVersion as release  # noqa: E402


def _tree(root):
    (root / "human_gem").mkdir()
    (root / "pyproject.toml").write_text('[project]\nname = "human-gem"\nversion = "2.1.0"\n', encoding="utf-8")
    (root / "human_gem" / "__init__.py").write_text('__version__ = "2.1.0"\n', encoding="utf-8")
    (root / "citation.cff").write_text("title: Human-GEM\nversion: 2.1.0\ndate-released: 2026-09-30\n",
                                       encoding="utf-8")
    return root


def test_stamp_version_files(tmp_path):
    release.stamp_version_files("2.2.0", "2026-11-01", _tree(tmp_path))
    assert 'version = "2.2.0"' in (tmp_path / "pyproject.toml").read_text()
    assert '__version__ = "2.2.0"' in (tmp_path / "human_gem" / "__init__.py").read_text()
    cff = (tmp_path / "citation.cff").read_text()
    assert "version: 2.2.0" in cff and "date-released: 2026-11-01" in cff


def test_stamp_refuses_a_file_without_its_line_and_writes_nothing(tmp_path):
    _tree(tmp_path)
    (tmp_path / "citation.cff").write_text("title: Human-GEM\n", encoding="utf-8")
    with pytest.raises(ValueError, match="citation.cff"):
        release.stamp_version_files("2.2.0", "2026-11-01", tmp_path)
    assert 'version = "2.1.0"' in (tmp_path / "pyproject.toml").read_text()


def test_the_repository_files_carry_one_version():
    root = Path(__file__).resolve().parents[3]
    if not (root / "pyproject.toml").is_file():
        pytest.skip("not a full checkout")
    import re
    found = {
        re.search(r'(?m)^version = "([^"]*)"', (root / "pyproject.toml").read_text()).group(1),
        re.search(r'(?m)^__version__ = "([^"]*)"', (root / "human_gem" / "__init__.py").read_text()).group(1),
        re.search(r"(?m)^version: (.*)$", (root / "citation.cff").read_text()).group(1).strip(),
    }
    assert len(found) == 1, found
