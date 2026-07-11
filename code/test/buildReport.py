"""Build one consolidated model-quality report for the pull-request comment.

Combines the three result sets so the model's quality shows up in a single
comment instead of several:

  * Model QC checks (the Model QC workflow): growth, metabolite completeness,
    reaction bound/GPR sanity, annotation cross-references and the MEMOTE score,
    as a compact status table with the change versus the target branch.
  * MACAW and mass/charge balance (the QC-tests workflow): its committed summary.
  * Gene essentiality, Hart 2015 (the gene-essentiality workflow): the per-cell-
    line metrics table.

Each of those workflows calls this script and posts the result to the same
comment identifier, rebuilding the whole comment from the committed result files,
so the last one to finish shows the complete picture. A result set a workflow has
not committed yet for this pull request shows as pending.

Icons: white_check_mark = fine / no regression, warning = changed vs the target
branch (review it), x = failed, new = no baseline to compare against,
hourglass = not committed yet for this pull request.

Usage:
    BASE_RESULTS_DIR=<dir> BASE_REF=<branch> COMMIT_SHA=<sha> \
        python code/test/buildReport.py
"""

import csv
import os
import re
import sys
from pathlib import Path

RESULTS = Path("data/testResults")
SUMMARY_MD = RESULTS / "model_qc_summary.md"
BASE_DIR = os.environ.get("BASE_RESULTS_DIR", "")
BASE_REF = os.environ.get("BASE_REF", "the target branch")
COMMIT_SHA = os.environ.get("COMMIT_SHA", "")

# Each row: (label, metric key, kind). kind sets which direction is "good":
#   count  -> lower is better (more findings is a regression)
#   score  -> higher is better
#   growth -> pass/fail (must be positive)
ROWS = [
    ("Growth (biomass producible)", "growth", "growth"),
    ("Metabolites missing formula", "missing_formula", "count"),
    ("Metabolites missing charge", "missing_charge", "count"),
    ("Reaction bound / GPR issues", "reaction_issues", "count"),
    ("Malformed cross-references", "malformed", "count"),
    ("Cross-refs inconsistent across compartments", "inconsistent", "count"),
    ("MEMOTE score (%)", "memote", "score"),
]


def _count_csv(path: Path, predicate=None) -> int | None:
    if not path.exists():
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        return sum(1 for row in csv.DictReader(fh) if predicate is None or predicate(row))


def _growth(directory: Path) -> float | None:
    try:
        return float((directory / "qc_growth.txt").read_text(encoding="utf-8").strip())
    except (FileNotFoundError, ValueError):
        return None


def _memote_score(directory: Path) -> float | None:
    path = directory / "memote_score.md"
    if not path.exists():
        return None
    match = re.search(r"Total score:\s*([\d.]+)\s*%", path.read_text(encoding="utf-8"))
    return float(match.group(1)) if match else None


def _metrics(directory: Path) -> dict:
    completeness = directory / "qc_metabolite_completeness.csv"
    annotation = directory / "qc_annotation_issues.csv"
    return {
        "growth": _growth(directory),
        "missing_formula": _count_csv(completeness, lambda r: r.get("missing_formula") == "yes"),
        "missing_charge": _count_csv(completeness, lambda r: r.get("missing_charge") == "yes"),
        "reaction_issues": _count_csv(directory / "qc_reaction_sanity.csv"),
        "malformed": _count_csv(annotation, lambda r: r.get("issue", "").startswith("malformed")),
        "inconsistent": _count_csv(annotation, lambda r: r.get("issue", "").startswith("inconsistent")),
        "memote": _memote_score(directory),
    }


def _format_value(value, kind: str) -> str:
    if kind == "growth":
        return f"{value:.3g}"
    if kind == "score":
        return f"{value:.1f}"
    return str(int(value))


def _delta_and_icon(current, base, kind: str) -> tuple[str, str]:
    if kind == "growth":
        # A pass/fail gate: show the verdict even without a baseline to compare against.
        icon = ":white_check_mark:" if current > 1e-6 else ":x:"
        if base is None:
            return "new", icon
        change = current - base
        return (f"{change:+.3g}" if abs(change) > 1e-9 else "0"), icon
    if base is None:
        return "new", ":new:"
    if kind == "score":
        change = current - base
        icon = ":white_check_mark:" if change >= -1e-9 else ":warning:"
        return (f"{change:+.1f}" if abs(change) > 1e-9 else "0"), icon
    change = int(current) - int(base)  # count: lower is better
    icon = ":white_check_mark:" if change <= 0 else ":warning:"
    return (f"{change:+d}" if change != 0 else "0"), icon


def _qc_section() -> tuple[list[str], str]:
    """Return (table lines, overall header) for the Model QC checks."""
    current = _metrics(RESULTS)
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    base = _metrics(Path(BASE_DIR)) if have_base else dict.fromkeys(current)

    table, changed, failed = [], 0, False
    for label, key, kind in ROWS:
        value = current[key]
        if value is None:  # not committed yet for this pull request
            table.append(f"| {label} | pending | | :hourglass_flowing_sand: |")
            continue
        delta, icon = _delta_and_icon(value, base.get(key), kind)
        changed += icon == ":warning:"
        failed = failed or icon == ":x:"
        table.append(f"| {label} | {_format_value(value, kind)} | {delta} | {icon} |")

    if failed:
        header = ":x: **A check failed.** See the detailed reports below."
    elif not have_base:
        header = ":information_source: First run for this comparison; no target-branch baseline yet."
    elif changed:
        header = f":warning: **{changed} Model QC check(s) changed** compared to `{BASE_REF}`."
    else:
        header = f":white_check_mark: **Model QC checks fine**, no changes compared to `{BASE_REF}`."
    return table, header


def _macaw_balance_section() -> str:
    """The MACAW / mass-and-charge-balance summary written by the QC-tests workflow."""
    path = RESULTS / "qc_summary.md"
    if path.exists() and path.read_text(encoding="utf-8").strip():
        return path.read_text(encoding="utf-8").strip()
    return "_Not yet run for this pull request._"


def _gene_essentiality_section() -> str:
    """Render the Hart 2015 per-cell-line metrics table."""
    path = RESULTS / "gene-essential.csv"
    if not path.exists():
        return "_Not yet run for this pull request._"
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if r]
    if len(rows) < 2:
        return "_No gene-essentiality results._"
    header, *body = rows
    lines = ["| " + " | ".join(header) + " |",
             "| " + " | ".join("---" for _ in header) + " |"]
    lines += ["| " + " | ".join(cell for cell in row) + " |" for row in body]
    return "\n".join(lines)


def main() -> int:
    table, header = _qc_section()

    lines = [
        "## Model quality report",
        "",
        header,
        "",
        "### Model QC checks",
        "",
        f"| Check | Result | &Delta; vs `{BASE_REF}` | |",
        "| --- | ---: | ---: | :---: |",
        *table,
        "",
        "### MACAW and mass/charge balance",
        "",
        _macaw_balance_section(),
        "",
        "### Gene essentiality (Hart 2015)",
        "",
        _gene_essentiality_section(),
        "",
        "Per-finding detail is committed to `data/testResults/` "
        "(`qc_metabolite_completeness.csv`, `qc_reaction_sanity.csv`, "
        "`qc_annotation_issues.csv`, `macaw_results.csv`, `balance_results.csv`, "
        "`gene-essential.csv`); the full MEMOTE result is uploaded as a build artifact.",
    ]
    if COMMIT_SHA:
        lines += ["", f"Results for commit {COMMIT_SHA[:7]}."]

    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
