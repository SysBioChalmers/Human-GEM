"""Model SBML round-trip check for Human-GEM.

``model/Human-GEM.xml`` is regenerated from ``model/Human-GEM.yml`` at release time
by ``code/io/increaseHumanGEMVersion.py``, with the TSV cross-references and SBO terms
merged in (see annotateGEM.py). Nothing reads it back afterwards, so a loss in the
SBML writer or reader reaches a release unnoticed. This writes the annotated model to
SBML, reads it back, and compares the two models, which catches annotation losses and
identifier rewrites for whatever the model currently holds.

The comparison runs on the annotated model rather than the plain one because the
annotated model is what is published as SBML.

Usage:
    python code/test/testSbmlConversion.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from cobra.io import read_sbml_model, write_sbml_model

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = REPO_ROOT / "model"

# annotateGEM.py sits in code/, next to this file's parent.
sys.path.insert(0, str(REPO_ROOT / "code"))
from annotateGEM import annotate_gem  # noqa: E402

from raven_toolbox.comparison import diff_models  # noqa: E402
from raven_toolbox.io import read_yaml_model  # noqa: E402


def normalise_annotations(model) -> None:
    """Reduce every annotation to the set of identifiers it holds.

    Two representations differ without the annotation differing. annotate_gem stores a
    single cross-reference as a one-element tuple, while the SBML reader returns a
    plain string for it; and an identifier listed twice on one entity is written once
    by the SBML writer. Both sides are therefore reduced to sorted, de-duplicated
    strings, which keeps the comparison about annotation content: a dropped, added or
    rewritten identifier still differs.
    """
    for entity in (*model.reactions, *model.metabolites, *model.genes):
        annotation = entity.annotation
        for key, value in list(annotation.items()):
            if isinstance(value, (list, tuple, set)):
                items = sorted({str(item) for item in value})
                annotation[key] = items[0] if len(items) == 1 else items
            else:
                annotation[key] = str(value)


def main() -> int:
    model = annotate_gem(read_yaml_model(MODEL_DIR / "Human-GEM.yml"), MODEL_DIR)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "Human-GEM.xml"
        write_sbml_model(model, str(path))
        reloaded = read_sbml_model(str(path))
    normalise_annotations(model)
    normalise_annotations(reloaded)
    report = diff_models(model, reloaded)
    print(report)
    return 0 if report.equal else 1


if __name__ == "__main__":
    raise SystemExit(main())
