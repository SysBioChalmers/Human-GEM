"""Cut a new Human-GEM release: bump the version and regenerate the model exports.

Python/raven-toolbox port of ``code/io/increaseHumanGEMVersion.m``.

Reads ``model/Human-GEM.yml``, checks it against the annotation tables, then
regenerates every export in ``model/``:

* ``Human-GEM.yml`` and ``Human-GEM.mat`` from the plain model (cross-references
  stay in the TSV tables);
* ``Human-GEM.xml`` (SBML), ``Human-GEM.xlsx`` and ``Human-GEM.txt`` from a copy
  that has the TSV cross-references and SBO terms merged in (see annotateGEM.py).

``build`` additionally writes ``version.txt``, stamps the version into the model
metadata, and fills the ``{{nRXN}}`` / ``{{nMET}}`` / ``{{nGENE}}`` / ``{{DATE}}``
placeholders in ``README.md``. Those placeholders are the state ``develop`` keeps;
only the release branch and ``main`` carry the substituted values.

The release notes for a version live in ``docs/releaseNotes/<version>.md`` and are
written before the release branch is cut. ``validate`` refuses a version whose notes
file is missing, so the release workflow fails before it creates a branch rather than
partway through the build.

Usage:
    python code/io/increaseHumanGEMVersion.py validate --version 2.1.0
    python code/io/increaseHumanGEMVersion.py build --version 2.1.0
    python code/io/increaseHumanGEMVersion.py export   # regenerate exports only
"""
from __future__ import annotations

import argparse
import datetime
import subprocess
import sys
from pathlib import Path

import cobra
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = REPO_ROOT / "model"
VERSION_TXT = REPO_ROOT / "version.txt"
RELEASE_NOTES_DIR = REPO_ROOT / "docs" / "releaseNotes"

# Make the sibling code/annotateGEM.py importable regardless of the caller's cwd.
sys.path.insert(0, str(REPO_ROOT / "code"))
from annotateGEM import annotate_gem  # noqa: E402

from raven_toolbox.io import export_for_git, read_yaml_model  # noqa: E402

# model attribute <-> TSV file <-> id column, for the consistency check.
_ID_TABLES = (
    ("reactions", "reactions.tsv", "rxns"),
    ("metabolites", "metabolites.tsv", "mets"),
    ("genes", "genes.tsv", "genes"),
)


def _parse_version(text: str) -> tuple[int, int, int]:
    parts = text.strip().split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise SystemExit(f"not a x.y.z version: {text.strip()!r}")
    major, minor, patch = (int(p) for p in parts)
    return major, minor, patch


def _legal_bumps(old: tuple[int, int, int]) -> dict[str, tuple[int, int, int]]:
    major, minor, patch = old
    return {
        "major": (major + 1, 0, 0),
        "minor": (major, minor + 1, 0),
        "patch": (major, minor, patch + 1),
    }


def _current_version_text() -> str:
    """The currently released version.

    ``version.txt`` is tracked on ``main`` only (see branch-hygiene.yml), while
    ``validate`` runs from a develop checkout when cutting a release. Read the local
    file when the checkout has one (``main``, or an already-cut release branch);
    otherwise read it from ``origin/main``.
    """
    if VERSION_TXT.is_file():
        return VERSION_TXT.read_text(encoding="utf-8")
    subprocess.run(
        ["git", "fetch", "--depth=1", "origin", "main"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    )
    result = subprocess.run(
        ["git", "show", "origin/main:version.txt"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    )
    return result.stdout


def release_notes_path(version: str) -> Path:
    return RELEASE_NOTES_DIR / f"{version}.md"


def _check_version(version: str) -> None:
    """Error unless ``version`` is a legal increment with release notes present."""
    old = _parse_version(_current_version_text())
    new = _parse_version(version)
    allowed = _legal_bumps(old)
    if new not in allowed.values():
        options = ", ".join(".".join(map(str, v)) for v in allowed.values())
        raise SystemExit(
            f"{version} is not a major, minor or patch increment of "
            f"{'.'.join(map(str, old))}. Expected one of: {options}."
        )

    notes = release_notes_path(version)
    if not notes.is_file() or not notes.read_text(encoding="utf-8").strip():
        raise SystemExit(
            f"{notes.relative_to(REPO_ROOT).as_posix()} is missing or empty. "
            "Write the release notes before cutting the release branch."
        )


def _check_tsv_consistency(model: cobra.Model) -> None:
    """Error if any model id is missing from its TSV table, or vice versa."""
    problems = []
    for attr, fname, col in _ID_TABLES:
        table = pd.read_csv(MODEL_DIR / fname, sep="\t", dtype=str, keep_default_na=False)
        tsv_ids = set(table[col])
        model_ids = {entity.id for entity in getattr(model, attr)}
        only_model = sorted(model_ids - tsv_ids)
        only_tsv = sorted(tsv_ids - model_ids)
        if only_model:
            problems.append(f"in model.{attr} but not {fname}: {only_model}")
        if only_tsv:
            problems.append(f"in {fname} but not model.{attr}: {only_tsv}")
    if problems:
        raise ValueError("Model / TSV mismatch:\n  " + "\n  ".join(problems))


def _set_version(model: cobra.Model, new_version: str) -> None:
    """Write the version into the metaData block that write_yaml_model emits."""
    notes = model.notes or {}
    meta = dict(notes.get("metaData") or {})
    meta["version"] = new_version          # metaData wins in write_yaml_model
    notes["metaData"] = meta
    notes["version"] = new_version
    model.notes = notes


def _update_readme(model: cobra.Model) -> None:
    readme = REPO_ROOT / "README.md"
    content = readme.read_text(encoding="utf-8")
    today = datetime.date.today().isoformat()
    for token, value in (
        ("{{DATE}}", today),
        ("{{nRXN}}", str(len(model.reactions))),
        ("{{nMET}}", str(len(model.metabolites))),
        ("{{nGENE}}", str(len(model.genes))),
    ):
        content = content.replace(token, value)
    readme.write_text(content, encoding="utf-8")


def _export(model: cobra.Model) -> None:
    """Write every derived model file from ``model``.

    The plain formats (yml, mat) keep their cross-references in the TSV tables; the
    annotated formats (xml, xlsx, txt) carry the merged TSV cross-references and SBO
    terms (see annotateGEM.py). export_for_git also (re)writes model/dependencies.txt.
    varname pins the .mat struct name to "humanGEM".
    """
    export_for_git(model, MODEL_DIR, prefix="Human-GEM",
                   formats=("yml", "mat"), sub_dirs=False, varname="humanGEM")
    export_for_git(annotate_gem(model.copy(), MODEL_DIR), MODEL_DIR,
                   prefix="Human-GEM", formats=("xml", "xlsx", "txt"), sub_dirs=False)


def cmd_validate(args: argparse.Namespace) -> int:
    """Check the version and its release notes. No side effects."""
    _check_version(args.version)
    print(f"{args.version} is a legal increment of {_current_version_text().strip()}, "
          f"and {release_notes_path(args.version).relative_to(REPO_ROOT).as_posix()} is ready.")
    return 0


def cmd_build(args: argparse.Namespace) -> int:
    """Stamp the version, regenerate the exports and fill the README placeholders."""
    _check_version(args.version)
    model = read_yaml_model(MODEL_DIR / "Human-GEM.yml")
    _set_version(model, args.version)
    _check_tsv_consistency(model)
    _export(model)
    VERSION_TXT.write_text(args.version, encoding="utf-8")
    _update_readme(model)
    print(f"Human-GEM bumped to {args.version}")
    return 0


def cmd_export(_args: argparse.Namespace) -> int:
    """Regenerate the exports from the committed model, leaving the version alone."""
    model = read_yaml_model(MODEL_DIR / "Human-GEM.yml")
    _check_tsv_consistency(model)
    _export(model)
    print("Exports regenerated, version unchanged.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Bump the Human-GEM version and regenerate the model exports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate", help="check the version increment and release notes")
    p.add_argument("--version", required=True, help="new version, e.g. 2.1.0")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("build", help="stamp the version and refresh the release files")
    p.add_argument("--version", required=True, help="new version, e.g. 2.1.0")
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("export", help="regenerate the model exports only")
    p.set_defaults(func=cmd_export)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
