"""The annotation merge writes each KEGG id of metKEGGID to the namespace of its KEGG database."""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
for path in (ROOT, ROOT / "code"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
try:
    from human_gem import annotation as module
except ImportError:
    import annotateGEM as module

MET_ID2MIRIAM = getattr(module, "MET_ID2MIRIAM", None) or module._MET_ID2MIRIAM


def test_kegg_ids_go_to_compound_drug_and_glycan_namespaces():
    annotation = {}
    row = pd.Series({"metKEGGID": "C00031;D00009;G00054", "metChEBIID": "4167"})
    module._apply_row(annotation, row, MET_ID2MIRIAM)
    assert annotation["kegg.compound"] == ["C00031"]
    assert annotation["kegg.drug"] == ["D00009"]
    assert annotation["kegg.glycan"] == ["G00054"]
    assert annotation["chebi"] == ["CHEBI:4167"]


def test_compound_only_cell_is_unchanged():
    annotation = {}
    module._apply_row(annotation, pd.Series({"metKEGGID": "C00002"}), MET_ID2MIRIAM)
    assert annotation == {"kegg.compound": ["C00002"]}
