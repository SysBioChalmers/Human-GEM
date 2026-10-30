"""Attach cross-reference (MIRIAM) and SBO annotation to a Human-GEM model.

The implementation is ``human_gem/annotation.py`` (the module ``human_gem.load_model()``
uses); this file makes it importable from ``code/`` for the release
(``code/io/increaseHumanGEMVersion.py``) and the CI checks, which run the scripts as
files. Run standalone, it prints how much of the model the tables annotate:

    python code/annotateGEM.py            # counts only, model not modified on disk
"""
from __future__ import annotations

import sys
from pathlib import Path

# The human_gem package sits at the repository root, one level above code/.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from human_gem.annotation import (  # noqa: E402,F401
    GENE_ID2MIRIAM,
    MET_ID2MIRIAM,
    RXN_ID2MIRIAM,
    add_inchis,
    annotate_gem,
)


def _main() -> int:
    from raven_toolbox.io import read_yaml_model

    repo_root = Path(__file__).resolve().parents[1]
    model_dir = repo_root / "model"
    model = read_yaml_model(model_dir / "Human-GEM.yml")
    annotate_gem(model, model_dir)
    n_rxn = sum(1 for r in model.reactions if any(k != "sbo" for k in r.annotation))
    n_met = sum(1 for m in model.metabolites if m.annotation)
    n_gene = sum(1 for g in model.genes if g.annotation)
    n_inchi = sum(1 for m in model.metabolites if m.notes.get("inchis"))
    print(f"annotated reactions (cross-refs): {n_rxn}/{len(model.reactions)}")
    print(f"annotated metabolites:            {n_met}/{len(model.metabolites)}")
    print(f"annotated genes:                  {n_gene}/{len(model.genes)}")
    print(f"metabolites with an InChI:        {n_inchi}/{len(model.metabolites)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
