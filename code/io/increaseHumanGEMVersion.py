"""Cut a new Human-GEM release: bump the version and regenerate the model exports.

Python/raven-toolbox port of ``code/io/increaseHumanGEMVersion.m``.

Reads ``model/Human-GEM.yml``, checks it against the annotation tables, then
regenerates every export in ``model/``:

* ``Human-GEM.yml`` and ``Human-GEM.mat`` from the plain model (cross-references
  stay in the TSV tables);
* ``Human-GEM.xml`` (SBML), ``Human-GEM.xlsx`` and ``Human-GEM.txt`` from a copy
  that has the TSV cross-references and SBO terms merged in (see annotateGEM.py).

Outside test mode it also refuses to run off ``main``, bumps ``version.txt`` and
fills the ``{{nRXN}}`` / ``{{nMET}}`` / ``{{nGENE}}`` / ``{{DATE}}`` placeholders
in ``README.md``.

Usage:
    python code/io/increaseHumanGEMVersion.py {major|minor|patch}
    python code/io/increaseHumanGEMVersion.py patch --test   # regenerate only
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

# Make the sibling code/annotateGEM.py importable regardless of the caller's cwd.
sys.path.insert(0, str(REPO_ROOT / "code"))
from annotateGEM import annotate_gem  # noqa: E402

from raven_toolbox.io import export_for_git, read_yaml_model, write_yaml_model  # noqa: E402

# model attribute <-> TSV file <-> id column, for the consistency check.
_ID_TABLES = (
    ("reactions", "reactions.tsv", "rxns"),
    ("metabolites", "metabolites.tsv", "mets"),
    ("genes", "genes.tsv", "genes"),
)


def _bump(old: str, bump_type: str) -> str:
    """Return the ``major``/``minor``/``patch`` increment of a ``x.y.z`` string."""
    parts = [int(p) for p in old.strip().split(".")]
    if len(parts) != 3:
        raise ValueError(f"version.txt is not x.y.z: {old!r}")
    major, minor, patch = parts
    if bump_type == "major":
        major, minor, patch = major + 1, 0, 0
    elif bump_type == "minor":
        minor, patch = minor + 1, 0
    elif bump_type == "patch":
        patch += 1
    else:
        raise ValueError('bump_type must be "major", "minor" or "patch"')
    return f"{major}.{minor}.{patch}"


def _current_branch() -> str:
    out = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    )
    return out.stdout.strip()


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


def increase_human_gem_version(bump_type: str, test: bool = False) -> str | None:
    """Regenerate the model exports and (unless ``test``) bump the version.

    Returns the new version string, or ``None`` in test mode.
    """
    version_file = REPO_ROOT / "version.txt"
    new_version = None

    if not test:
        branch = _current_branch()
        if branch != "main":
            raise RuntimeError(f"not on main (current branch: {branch})")
        new_version = _bump(version_file.read_text(encoding="utf-8"), bump_type)

    model = read_yaml_model(MODEL_DIR / "Human-GEM.yml")

    if not test:
        _set_version(model, new_version)

    _check_tsv_consistency(model)

    # Plain exports keep their cross-references in the TSV tables. YAML via
    # raven-toolbox; MATLAB via cobra with the "humanGEM" variable name (which
    # export_for_git would not set - it uses cobra's default, the model id).
    write_yaml_model(model, MODEL_DIR / "Human-GEM.yml")
    cobra.io.save_matlab_model(model, str(MODEL_DIR / "Human-GEM.mat"), varname="humanGEM")

    # Annotated exports (xml/xlsx/txt) carry the merged TSV cross-references and SBO
    # terms (see annotateGEM.py); export_for_git also (re)writes dependencies.txt.
    export_for_git(annotate_gem(model.copy(), MODEL_DIR), MODEL_DIR,
                   prefix="Human-GEM", formats=("xml", "xlsx", "txt"), sub_dirs=False)

    if not test:
        version_file.write_text(new_version, encoding="utf-8")
        _update_readme(model)
        print(f"Human-GEM bumped to {new_version}")
    else:
        print("Test run: exports regenerated, version unchanged.")

    return new_version


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bump_type", choices=("major", "minor", "patch"),
                        help="which part of the semantic version to increment")
    parser.add_argument("--test", action="store_true",
                        help="regenerate the exports without checking the branch or "
                             "bumping the version (may run on a development branch)")
    args = parser.parse_args(argv)
    increase_human_gem_version(args.bump_type, test=args.test)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
