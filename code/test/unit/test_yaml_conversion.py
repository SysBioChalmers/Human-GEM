"""The YAML round-trip compares the first write with the committed file."""

import testYamlConversion as ty

YAML = """!!omap
- metaData: !!omap
  - id: HumanGEM
  - version: ""
- metabolites:
  - !!omap
    - id: MAM00001c
    - name: "NO"
    - charge: 0.0
    - metFrom: HMRdatabase
- reactions:
  - !!omap
    - id: MAR00001
    - metabolites: !!omap
      - MAM00001c: -1.0
    - rxnFrom: HMRdatabase
- genes: []
- compartments: !!omap
  - c: cytosol
"""


def test_load_sections_reads_omaps(tmp_path):
    path = tmp_path / "m.yml"
    path.write_text(YAML, encoding="utf-8")
    sections = ty.load_sections(path)
    assert sections["metabolites"]["MAM00001c"]["name"] == "NO"
    assert sections["reactions"]["MAR00001"]["metabolites"] == {"MAM00001c": -1.0}
    assert sections["metaData"]["metaData"]["id"] == "HumanGEM"
    assert sections["compartments"]["compartments"] == {"c": "cytosol"}


def _sections():
    return {
        "metaData": {"metaData": {"id": "HumanGEM"}},
        "metabolites": {"M1": {"id": "M1", "charge": 0.0, "metFrom": "HMR"}},
        "reactions": {"R1": {"id": "R1", "lower_bound": -1000.0, "references": "PMID:1"}},
    }


def test_known_cobra_losses_are_counted_not_failed():
    written = {"metabolites": {"M1": {"id": "M1", "charge": 0}},
               "reactions": {"R1": {"id": "R1", "lower_bound": -1000}}}
    unexpected, expected = ty.compare_first_read(_sections(), written, ty.COBRA_DROPS)
    assert unexpected == []
    assert expected == {"metaData (whole section)": 1, "metabolites.metFrom": 1, "reactions.references": 1}


def test_any_other_loss_or_change_fails():
    written = {"metaData": {"metaData": {"id": "HumanGEM"}},
               "metabolites": {"M1": {"id": "M1", "charge": 1.0, "metFrom": "HMR"}},
               "reactions": {}}
    unexpected, _ = ty.compare_first_read(_sections(), written, {})
    assert "metabolites M1.charge: 0.0 -> 1.0" in unexpected
    assert "reactions R1: missing after the first read" in unexpected
