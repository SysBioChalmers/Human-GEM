"""Model SBML round-trip check for Human-GEM.

``model/Human-GEM.xml`` is regenerated from ``model/Human-GEM.yml`` at release time
by ``code/io/increaseHumanGEMVersion.py``: the annotated model (TSV cross-references
and SBO terms merged in, see annotateGEM.py) goes through raven-toolbox's
``export_for_git``, which also writes the subsystems as SBML groups and renames the
reaction notes to the labels RAVEN reads. Nothing reads the file back afterwards, so
a loss in that writer or in the SBML reader reaches a release unnoticed. This check
writes the SBML exactly as the release does, validates it with libSBML, reads it back
and compares it with the annotated model:

  * identifiers, stoichiometry, bounds, objective, GPRs, formulas, charges,
    compartments and annotations (raven-toolbox ``diff_models``);
  * reaction, metabolite and gene names;
  * every reaction's subsystems, read back from the SBML groups;
  * the reaction notes (note, references, confidence score).

The comparison runs on the annotated model rather than the plain one because the
annotated model is what is published as SBML.

Usage:
    python code/test/testSbmlConversion.py
"""

from __future__ import annotations

import html
import sys
import tempfile
from pathlib import Path

from cobra.io import validate_sbml_model

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = REPO_ROOT / "model"

# annotateGEM.py sits in code/, next to this file's parent.
sys.path.insert(0, str(REPO_ROOT / "code"))
from annotateGEM import annotate_gem  # noqa: E402

from raven_toolbox.comparison import diff_models  # noqa: E402
from raven_toolbox.io import export_for_git, read_yaml_model  # noqa: E402

# Reaction notes as the YAML holds them -> the label export_for_git writes them under.
NOTE_LABELS = {
    "note": "NOTES",
    "references": "AUTHORS",
    "confidence_score": "Confidence Level",
}

# Differences listed per check before the rest are only counted.
MAX_LISTED = 20


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


def _subsystems(reaction) -> set[str]:
    value = reaction.subsystem
    names = value if isinstance(value, (list, tuple, set)) else [value]
    return {str(name).strip() for name in names if name and str(name).strip()}


def _same_value(a, b) -> bool:
    """Equal as numbers when both are numbers, else as stripped strings.

    cobra's SBML reader takes a note's text from the XML without resolving character
    references, so an apostrophe written as ``&apos;`` comes back as ``&apos;``; the
    read-back value is unescaped before the comparison.
    """
    b = html.unescape(b) if isinstance(b, str) else b
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a).strip() == str(b).strip()


def compare_names(model, reloaded) -> list[str]:
    diffs = []
    for kind in ("reactions", "metabolites", "genes"):
        other = getattr(reloaded, kind)
        for entity in getattr(model, kind):
            if other.has_id(entity.id) and (entity.name or "") != (other.get_by_id(entity.id).name or ""):
                diffs.append(f"{kind[:-1]} {entity.id} name {entity.name!r} -> "
                             f"{other.get_by_id(entity.id).name!r}")
    return diffs


def compare_subsystems(model, reloaded) -> list[str]:
    """Each reaction's subsystems against the SBML groups it is a member of."""
    groups_of: dict[str, set[str]] = {}
    for group in reloaded.groups:
        for member in group.members:
            groups_of.setdefault(member.id, set()).add(str(group.name).strip())
    diffs = []
    for reaction in model.reactions:
        expected = _subsystems(reaction)
        found = groups_of.get(reaction.id, set())
        if expected != found:
            diffs.append(f"reaction {reaction.id} subsystems {sorted(expected)} -> {sorted(found)}")
    return diffs


def compare_notes(model, reloaded) -> list[str]:
    """The reaction notes, under the labels the SBML export gives them."""
    diffs = []
    for reaction in model.reactions:
        if not reloaded.reactions.has_id(reaction.id):
            continue
        written = reloaded.reactions.get_by_id(reaction.id).notes or {}
        for key, value in (reaction.notes or {}).items():
            if key not in NOTE_LABELS or value in (None, ""):
                continue
            label = NOTE_LABELS[key]
            if label not in written:
                diffs.append(f"reaction {reaction.id} note {key!r} missing")
            elif not _same_value(value, written[label]):
                diffs.append(f"reaction {reaction.id} note {key!r} {value!r} -> {written[label]!r}")
    return diffs


def _report(label: str, diffs: list[str]) -> bool:
    if not diffs:
        print(f"{label}: identical.")
        return True
    print(f"::error::{label}: {len(diffs)} difference(s).")
    for line in diffs[:MAX_LISTED]:
        print(f"  {line}")
    if len(diffs) > MAX_LISTED:
        print(f"  ... and {len(diffs) - MAX_LISTED} more")
    return False


def main() -> int:
    model = annotate_gem(read_yaml_model(MODEL_DIR / "Human-GEM.yml"), MODEL_DIR)
    with tempfile.TemporaryDirectory() as tmp:
        export_for_git(model, tmp, prefix="Human-GEM", formats=("xml",), sub_dirs=False)
        reloaded, errors = validate_sbml_model(str(Path(tmp) / "Human-GEM.xml"))
    # Warnings (for example units) are listed; errors and fatal errors fail the check.
    validation = [f"{kind}: {message}" for kind, messages in errors.items()
                  if not kind.endswith("WARNING") for message in messages]
    n_warnings = sum(len(messages) for kind, messages in errors.items() if kind.endswith("WARNING"))
    if n_warnings:
        print(f"libSBML validation: {n_warnings} warning(s).")
    ok = _report("libSBML validation", validation)
    if reloaded is None:
        print("::error::The SBML file could not be read back.")
        return 1

    ok &= _report("Subsystems (SBML groups)", compare_subsystems(model, reloaded))
    ok &= _report("Names", compare_names(model, reloaded))
    ok &= _report("Reaction notes", compare_notes(model, reloaded))
    normalise_annotations(model)
    normalise_annotations(reloaded)
    report = diff_models(model, reloaded)
    print(report)
    ok &= report.equal
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
