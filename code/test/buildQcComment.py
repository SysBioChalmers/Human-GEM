"""Build a compact, visual Model QC summary for the pull-request comment.

Reads the current QC result files in data/testResults/ and, when available, the
target-branch copies in the directory named by the BASE_RESULTS_DIR environment
variable, then writes a Markdown summary to data/testResults/model_qc_summary.md.
The summary is a per-check status table: current value, change versus the target
branch, and an icon for a quick visual check. The detailed per-finding output
stays in the committed CSVs, which is where you look to find out why something
changed or failed.

Icons: white_check_mark = fine / no regression, warning = changed vs the target
branch (review it), x = failed, new = no baseline to compare against.

Usage:
    BASE_RESULTS_DIR=<dir> BASE_REF=<branch> COMMIT_SHA=<sha> \
        python code/test/buildQcComment.py
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
    path = directory / "qc_growth.txt"
    try:
        return float(path.read_text(encoding="utf-8").strip())
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
    if value is None:
        return "n/a"
    if kind == "growth":
        return f"{value:.3g}"
    if kind == "score":
        return f"{value:.1f}"
    return str(int(value))


def _delta_and_icon(current, base, kind: str) -> tuple[str, str]:
    if current is None:
        return "n/a", ":question:"
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


def main() -> int:
    current = _metrics(RESULTS)
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    base = _metrics(Path(BASE_DIR)) if have_base else dict.fromkeys(current)

    memote_pending = bool(os.environ.get("MEMOTE_PENDING"))
    table, changed, failed = [], 0, False
    for label, key, kind in ROWS:
        if key == "memote" and memote_pending:
            # MEMOTE runs in its own (slow) job; show it as running for now.
            table.append(f"| {label} | running | | :hourglass_flowing_sand: |")
            continue
        delta, icon = _delta_and_icon(current[key], base.get(key), kind)
        changed += icon == ":warning:"
        failed = failed or icon == ":x:"
        table.append(f"| {label} | {_format_value(current[key], kind)} | {delta} | {icon} |")

    if failed:
        header = ":x: **A check failed.** See the detailed reports below."
    elif not have_base:
        header = ":information_source: First run for this comparison; no target-branch baseline yet."
    elif changed:
        header = f":warning: **{changed} check(s) changed** compared to `{BASE_REF}`."
    else:
        header = f":white_check_mark: **All checks fine**, no changes compared to `{BASE_REF}`."

    lines = [
        "## Model QC checks",
        "",
        header,
        "",
        f"| Check | Result | &Delta; vs `{BASE_REF}` | |",
        "| --- | ---: | ---: | :---: |",
        *table,
        "",
        "Detailed findings are committed to `data/testResults/`: "
        "`qc_metabolite_completeness.csv`, `qc_reaction_sanity.csv`, `qc_annotation_issues.csv`. "
        "The full MEMOTE result is uploaded as a build artifact. Look there to see which "
        "reactions, metabolites or identifiers changed.",
    ]
    if memote_pending:
        lines += ["", "_MEMOTE runs in a separate, slower job and will fill in its "
                  "row of this comment when it finishes (a fast core subset on every "
                  "pull request, the full suite on pull requests to `main`)._"]
    if COMMIT_SHA:
        lines += ["", f"Results for commit {COMMIT_SHA[:7]}."]

    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
