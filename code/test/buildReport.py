"""Build one consolidated model-quality report for the pull-request comment.

Combines the result sets so the model's quality shows up in a single comment
instead of several:

  * Model QC checks (the Model QC workflow): growth, metabolite completeness,
    reaction bound/GPR sanity, annotation cross-references and the MEMOTE score,
    as a compact status table with the change versus the target branch.
  * MACAW and mass/charge balance (the QC-tests workflow): dead-end / duplicate
    reaction findings and mass/charge imbalances, in the same status-table style.
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
QC_ROWS = [
    ("Growth (biomass producible)", "growth", "growth"),
    ("Metabolites missing formula", "missing_formula", "count"),
    ("Metabolites missing charge", "missing_charge", "count"),
    ("Reaction bound / GPR issues", "reaction_issues", "count"),
    ("Malformed cross-references", "malformed", "count"),
    ("Cross-refs inconsistent across compartments", "inconsistent", "count"),
    ("MEMOTE score (%)", "memote", "score"),
]

# MACAW and mass/charge balance, counted from the committed result CSVs so the
# summary carries the same value / delta / icon as the Model QC checks.
MB_ROWS = [
    ("Reactions flagged by MACAW dead-end test", "dead_end", "count"),
    ("Reactions flagged as MACAW duplicates", "duplicates", "count"),
    ("Mass-imbalanced reactions", "mass_imbalance", "count"),
    ("Charge-imbalanced reactions", "charge_imbalance", "count"),
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


def _qc_metrics(directory: Path) -> dict:
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


_DUP_COLS = ("duplicate_test_exact", "duplicate_test_directions", "duplicate_test_coefficients")


def _mb_metrics(directory: Path) -> dict:
    macaw = directory / "macaw_results.csv"
    balance = directory / "balance_results.csv"
    return {
        # dead_end_test is "ok" when fine, otherwise a reason (a metabolite list or
        # "only when going forwards/backwards"): a non-ok value is a flagged reaction.
        "dead_end": _count_csv(macaw, lambda r: r.get("dead_end_test", "") not in ("ok", "")),
        "duplicates": _count_csv(
            macaw, lambda r: any(r.get(c, "") not in ("ok", "N/A", "") for c in _DUP_COLS)
        ),
        "mass_imbalance": _count_csv(balance, lambda r: r.get("mass_imbalance", "").strip() != ""),
        "charge_imbalance": _count_csv(balance, lambda r: r.get("charge_imbalance", "").strip() != ""),
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


def _table(rows, current: dict, base: dict) -> tuple[list[str], int, bool]:
    """Build the status-table body for a set of rows. Returns (lines, changed, failed)."""
    lines, changed, failed = [], 0, False
    for label, key, kind in rows:
        value = current.get(key)
        if value is None:  # not committed yet for this pull request
            lines.append(f"| {label} | pending | | :hourglass_flowing_sand: |")
            continue
        delta, icon = _delta_and_icon(value, base.get(key), kind)
        changed += icon == ":warning:"
        failed = failed or icon == ":x:"
        lines.append(f"| {label} | {_format_value(value, kind)} | {delta} | {icon} |")
    return lines, changed, failed


def _header(changed: int, failed: bool, have_base: bool, what: str) -> str:
    if failed:
        return f":x: **A {what} check failed.** See the detailed reports below."
    if not have_base:
        return ":information_source: First run for this comparison; no target-branch baseline yet."
    if changed:
        return f":warning: **{changed} {what} check(s) changed** compared to `{BASE_REF}`."
    return f":white_check_mark: **{what} checks fine**, no changes compared to `{BASE_REF}`."


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
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    base_dir = Path(BASE_DIR) if have_base else None

    qc_cur = _qc_metrics(RESULTS)
    qc_base = _qc_metrics(base_dir) if have_base else dict.fromkeys(qc_cur)
    qc_table, qc_changed, qc_failed = _table(QC_ROWS, qc_cur, qc_base)

    mb_cur = _mb_metrics(RESULTS)
    mb_base = _mb_metrics(base_dir) if have_base else dict.fromkeys(mb_cur)
    mb_table, mb_changed, mb_failed = _table(MB_ROWS, mb_cur, mb_base)

    delta_head = f"| Check | Result | &Delta; vs `{BASE_REF}` | |"
    delta_sep = "| --- | ---: | ---: | :---: |"

    lines = [
        "## Model quality report",
        "",
        _header(qc_changed, qc_failed, have_base, "Model QC"),
        "",
        "### Model QC checks",
        "",
        delta_head,
        delta_sep,
        *qc_table,
        "",
        "### MACAW and mass/charge balance",
        "",
        _header(mb_changed, mb_failed, have_base, "MACAW / balance"),
        "",
        delta_head,
        delta_sep,
        *mb_table,
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
