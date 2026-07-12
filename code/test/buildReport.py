"""Build one model-quality report for the pull-request comment.

Turns the committed result files under data/testResults/ into a single comment
that leads with a one-line verdict (blocked / regressions / running / clean) and
then three status tables:

  * Build gates: the checks that fail the build (duplicate keys, empty reactions,
    model vs annotation-table consistency, growth). A red X here means the merge
    is blocked; the offending entries are in the linked CSV.
  * Model QC reports: quality metrics tracked with a delta versus the target
    branch (completeness, reaction sanity, exact duplicates, unused entities,
    cross-references, MEMOTE score).
  * MACAW and mass/charge balance: dead-end / duplicate / imbalance counts.

Each result set is stamped with the commit it was computed for (data/testResults/
qc_<group>.sha). If that does not match this pull request's head, the set has not
re-run for the current commit and its rows show as pending instead of showing
stale numbers as if they were current.

Icons: white_check_mark = fine / gate passing / no regression, warning = a report
metric worsened vs the target branch, x = a gate is failing, new = no baseline to
compare against, hourglass = not (re-)run for this commit yet.

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

# A result set is produced by one workflow job and stamped with the commit it ran
# on. buildReport marks a set pending when its stamp does not match this commit.
GROUP_STAMPS = {"checks": "qc_checks.sha", "memote": "qc_memote.sha", "macaw": "qc_macaw.sha"}

# Each row: (label, key, kind, group, detail_file). kind sets the good direction:
#   gate   -> 0 is passing, anything > 0 fails the build (icon x)
#   count  -> lower is better (a rise vs the target branch is a regression)
#   score  -> higher is better
#   growth -> pass/fail (must be positive)
GATE_ROWS = [
    ("Duplicate `!!omap` keys", "dup_keys", "gate", "checks", "qc_duplicate_keys.csv"),
    ("Reactions with no metabolites", "empty_rxn", "gate", "checks", "qc_empty_reactions.csv"),
    ("Model / annotation-table inconsistencies", "annot_consistency", "gate", "checks",
     "qc_annotation_consistency.csv"),
    ("Growth (biomass producible)", "growth", "growth", "checks", "qc_growth_blockers.csv"),
]
REPORT_ROWS = [
    ("Metabolites missing formula", "missing_formula", "count", "checks", "qc_metabolite_completeness.csv"),
    ("Metabolites missing charge", "missing_charge", "count", "checks", "qc_metabolite_completeness.csv"),
    ("Reaction bound / GPR issues", "reaction_issues", "count", "checks", "qc_reaction_sanity.csv"),
    ("Exact-duplicate reaction groups", "dup_reactions", "count", "checks", "qc_duplicate_reactions.csv"),
    ("Unused metabolites", "unused_met", "count", "checks", "qc_unused_entities.csv"),
    ("Unused genes", "unused_gene", "count", "checks", "qc_unused_entities.csv"),
    ("Malformed cross-references", "malformed", "count", "checks", "qc_annotation_issues.csv"),
    ("Cross-refs inconsistent across compartments", "inconsistent", "count", "checks", "qc_annotation_issues.csv"),
    ("MEMOTE score (%)", "memote", "score", "memote", "memote_score.md"),
]
MB_ROWS = [
    ("Reactions flagged by MACAW dead-end test", "dead_end", "count", "macaw", "macaw_results.csv"),
    ("Reactions flagged as MACAW duplicates", "duplicates", "count", "macaw", "macaw_results.csv"),
    ("Mass-imbalanced reactions", "mass_imbalance", "count", "macaw", "balance_results.csv"),
    ("Charge-imbalanced reactions", "charge_imbalance", "count", "macaw", "balance_results.csv"),
]
_DUP_COLS = ("duplicate_test_exact", "duplicate_test_directions", "duplicate_test_coefficients")


def _count_csv(path: Path, predicate=None) -> int | None:
    if not path.exists():
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        return sum(1 for row in csv.DictReader(fh) if predicate is None or predicate(row))


def _distinct_csv(path: Path, column: str) -> int | None:
    if not path.exists():
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        return len({row[column] for row in csv.DictReader(fh) if row.get(column)})


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
    unused = directory / "qc_unused_entities.csv"
    macaw = directory / "macaw_results.csv"
    balance = directory / "balance_results.csv"
    return {
        # gates
        "dup_keys": _count_csv(directory / "qc_duplicate_keys.csv"),
        "empty_rxn": _count_csv(directory / "qc_empty_reactions.csv"),
        "annot_consistency": _count_csv(directory / "qc_annotation_consistency.csv"),
        "growth": _growth(directory),
        # model QC reports
        "missing_formula": _count_csv(completeness, lambda r: r.get("missing_formula") == "yes"),
        "missing_charge": _count_csv(completeness, lambda r: r.get("missing_charge") == "yes"),
        "reaction_issues": _count_csv(directory / "qc_reaction_sanity.csv"),
        "dup_reactions": _distinct_csv(directory / "qc_duplicate_reactions.csv", "group"),
        "unused_met": _count_csv(unused, lambda r: r.get("kind") == "metabolite"),
        "unused_gene": _count_csv(unused, lambda r: r.get("kind") == "gene"),
        "malformed": _count_csv(annotation, lambda r: r.get("issue", "").startswith("malformed")),
        "inconsistent": _count_csv(annotation, lambda r: r.get("issue", "").startswith("inconsistent")),
        "memote": _memote_score(directory),
        # MACAW / balance
        "dead_end": _count_csv(macaw, lambda r: r.get("dead_end_test", "") not in ("ok", "")),
        "duplicates": _count_csv(macaw, lambda r: any(r.get(c, "") not in ("ok", "N/A", "") for c in _DUP_COLS)),
        "mass_imbalance": _count_csv(balance, lambda r: r.get("mass_imbalance", "").strip() != ""),
        "charge_imbalance": _count_csv(balance, lambda r: r.get("charge_imbalance", "").strip() != ""),
    }


def _stale_groups() -> set[str]:
    """Groups whose committed stamp does not match this pull request's head commit."""
    if not COMMIT_SHA:
        return set()
    stale = set()
    for group, stamp in GROUP_STAMPS.items():
        path = RESULTS / stamp
        # Only a present-but-different stamp counts as stale; a missing stamp is
        # unknown, so leave it to the file-absence (pending) handling below.
        if path.exists() and path.read_text(encoding="utf-8").strip() != COMMIT_SHA:
            stale.add(group)
    return stale


def _format_value(value, kind: str) -> str:
    if kind == "growth":
        return f"{value:.3g}"
    if kind == "score":
        return f"{value:.1f}"
    return str(int(value))


def _delta_and_icon(current, base, kind: str):
    """Return (delta_text, icon, worse). worse flags a report regression."""
    if kind == "gate":
        icon = ":x:" if current > 0 else ":white_check_mark:"
        if base is None:
            return "new", icon, False
        change = int(current) - int(base)
        return (f"{change:+d}" if change else "0"), icon, False
    if kind == "growth":
        icon = ":white_check_mark:" if current > 1e-6 else ":x:"
        if base is None:
            return "new", icon, False
        change = current - base
        return (f"{change:+.3g}" if abs(change) > 1e-6 else "0"), icon, False
    if base is None:
        return "new", ":new:", False
    if kind == "score":
        change = current - base
        worse = change < -1e-9
        return (f"{change:+.1f}" if abs(change) > 1e-9 else "0"), (":warning:" if worse else ":white_check_mark:"), worse
    change = int(current) - int(base)  # count: lower is better
    worse = change > 0
    return (f"{change:+d}" if change != 0 else "0"), (":warning:" if worse else ":white_check_mark:"), worse


def _table(rows, current: dict, base: dict, stale: set[str]):
    """Return (lines, n_failed_gates, n_worse, n_pending)."""
    lines, failed, worse, pending = [], 0, 0, 0
    for label, key, kind, group, detail in rows:
        value = current.get(key)
        if value is None or group in stale:
            lines.append(f"| {label} | pending | | :hourglass_flowing_sand: |")
            pending += 1
            continue
        delta, icon, is_worse = _delta_and_icon(value, base.get(key), kind)
        if icon == ":x:":
            failed += 1
        worse += is_worse
        cell = f"{_format_value(value, kind)}"
        if detail and value and kind != "growth":
            cell = f"[{cell}]({detail})"  # link the count to its detail file
        lines.append(f"| {label} | {cell} | {delta} | {icon} |")
    return lines, failed, worse, pending


def _gene_essentiality_section() -> str:
    path = RESULTS / "gene-essential.csv"
    if not path.exists():
        return "_Not yet run for this pull request._"
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if r]
    if len(rows) < 2:
        return "_No gene-essentiality results._"
    header, *body = rows
    lines = ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    lines += ["| " + " | ".join(row) + " |" for row in body]
    return "\n".join(lines)


def main() -> int:
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    base_dir = Path(BASE_DIR) if have_base else None
    stale = _stale_groups()

    current = _metrics(RESULTS)
    base = _metrics(base_dir) if have_base else {}

    gate_tbl, gate_fail, _, gate_pending = _table(GATE_ROWS, current, base, stale)
    rep_tbl, _, rep_worse, rep_pending = _table(REPORT_ROWS, current, base, stale)
    mb_tbl, _, mb_worse, mb_pending = _table(MB_ROWS, current, base, stale)

    worse = rep_worse + mb_worse
    pending = gate_pending + rep_pending + mb_pending

    if gate_fail:
        verdict = (f":x: **Merge blocked: {gate_fail} build gate(s) failing.** "
                   "Open the linked CSV(s) for the exact entries to fix.")
    elif worse:
        verdict = f":warning: **{worse} metric(s) worse than `{BASE_REF}`.** Review the changes below."
    elif pending and not have_base:
        verdict = ":information_source: First run for this comparison; results are still coming in."
    elif pending:
        verdict = f":hourglass_flowing_sand: **Some results are still running** ({pending} pending); the rest are unchanged vs `{BASE_REF}`."
    elif not have_base:
        verdict = ":information_source: First run for this comparison; no target-branch baseline yet."
    else:
        verdict = f":white_check_mark: **All checks pass, no regressions vs `{BASE_REF}`.**"

    head = f"| Check | Result | &Delta; vs `{BASE_REF}` | |"
    sep = "| --- | ---: | ---: | :---: |"
    lines = [
        "## Model quality report",
        "",
        verdict,
        "",
        "### Build gates",
        "_A red :x: blocks the merge; the count links to the CSV listing what to fix._",
        "",
        head, sep, *gate_tbl,
        "",
        "### Model QC reports",
        "",
        head, sep, *rep_tbl,
        "",
        "### MACAW and mass/charge balance",
        "",
        head, sep, *mb_tbl,
        "",
        "### Gene essentiality (Hart 2015)",
        "",
        _gene_essentiality_section(),
        "",
        "Per-finding detail is in the linked CSVs under `data/testResults/`; "
        "the full MEMOTE result is uploaded as a build artifact.",
    ]
    if COMMIT_SHA:
        lines += ["", f"Results for commit {COMMIT_SHA[:7]}."]

    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
