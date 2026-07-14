"""Build one model-quality report for the pull-request comment.

Turns the result files under data/testResults/ into a single comment that leads
with a one-line verdict and then three status tables (structural checks, model QC
reports, MACAW/balance). Groups still being computed on this run are passed in the
RUNNING_GROUPS environment variable and their rows show as *running* (hourglass);
the workflow calls this once with "all" before anything has run, once with
"memote" while the slow MEMOTE snapshot is still going, and once with nothing when
everything is in. No stamp files are involved.

Icon rule (per row):
  * growth: white_check_mark if the model grows, x if it cannot (blocks the merge).
  * MEMOTE score: warning if the score dropped versus the target branch, else
    white_check_mark (a non-zero score is good).
  * every other (count) metric: x if the count rose versus the target branch (a
    regression this pull request introduced), warning if the count is non-zero
    (a pre-existing finding, non-blocking), white_check_mark if it is zero.
  * hourglass: the group is still running on this pull request.

Only two conditions fail the build: the model cannot load (duplicate `!!omap`
keys) or cannot grow. Everything else is reported but does not block; a red x
just flags a regression for review.

Usage:
    RUNNING_GROUPS=<all|memote|...> BASE_RESULTS_DIR=<dir> BASE_REF=<branch> \
        RESULTS_URL_BASE=<url> python code/test/buildReport.py
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
# e.g. https://github.com/OWNER/REPO/blob/<branch>/data/testResults - used to link
# each finding count to its CSV. Empty when run locally (then counts are plain text).
URL_BASE = os.environ.get("RESULTS_URL_BASE", "").rstrip("/")

# Groups whose results are still being computed on this run; their rows show as
# "running". The workflow passes this on each call - "all" before anything has run,
# "memote" while the (slow) MEMOTE snapshot is still going, empty once everything is
# in - so no commit or stamp file is needed to track freshness.
ALL_GROUPS = {"checks", "memote", "macaw"}
_running = os.environ.get("RUNNING_GROUPS", "")
RUNNING = set(ALL_GROUPS) if _running.strip() == "all" else {g.strip() for g in _running.split(",") if g.strip()}

# (label, key, kind, group, detail_file)
STRUCTURAL_ROWS = [
    ("Duplicate `!!omap` keys", "dup_keys", "count", "checks", "qc_duplicate_keys.csv"),
    ("Reactions with no metabolites", "empty_rxn", "count", "checks", "qc_empty_reactions.csv"),
    ("Model / annotation-table inconsistencies", "annot_consistency", "count", "checks",
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
        "dup_keys": _count_csv(directory / "qc_duplicate_keys.csv"),
        "empty_rxn": _count_csv(directory / "qc_empty_reactions.csv"),
        "annot_consistency": _count_csv(directory / "qc_annotation_consistency.csv"),
        "growth": _growth(directory),
        "missing_formula": _count_csv(completeness, lambda r: r.get("missing_formula") == "yes"),
        "missing_charge": _count_csv(completeness, lambda r: r.get("missing_charge") == "yes"),
        "reaction_issues": _count_csv(directory / "qc_reaction_sanity.csv"),
        "dup_reactions": _distinct_csv(directory / "qc_duplicate_reactions.csv", "group"),
        "unused_met": _count_csv(unused, lambda r: r.get("kind") == "metabolite"),
        "unused_gene": _count_csv(unused, lambda r: r.get("kind") == "gene"),
        "malformed": _count_csv(annotation, lambda r: r.get("issue", "").startswith("malformed")),
        "inconsistent": _count_csv(annotation, lambda r: r.get("issue", "").startswith("inconsistent")),
        "memote": _memote_score(directory),
        "dead_end": _count_csv(macaw, lambda r: r.get("dead_end_test", "") not in ("ok", "")),
        "duplicates": _count_csv(macaw, lambda r: any(r.get(c, "") not in ("ok", "N/A", "") for c in _DUP_COLS)),
        "mass_imbalance": _count_csv(balance, lambda r: r.get("mass_imbalance", "").strip() != ""),
        "charge_imbalance": _count_csv(balance, lambda r: r.get("charge_imbalance", "").strip() != ""),
    }


def _icon(value, base, kind):
    """Return (delta_text, icon, regression, fatal)."""
    if kind == "growth":
        grows = value > 1e-6
        icon = ":white_check_mark:" if grows else ":x:"
        if base is None:
            return "new", icon, False, (not grows)
        change = value - base
        return (f"{change:+.3g}" if abs(change) > 1e-6 else "0"), icon, False, (not grows)
    if kind == "score":  # higher is better; a drop is a (non-blocking) warning
        if base is None:
            return "new", (":white_check_mark:" if value > 0 else ":warning:"), False, False
        change = value - base
        dropped = change < -1e-9
        return (f"{change:+.1f}" if abs(change) > 1e-9 else "0"), (":warning:" if dropped else ":white_check_mark:"), dropped, False
    # count: rose vs base -> regression (x); non-zero -> pre-existing (warning); zero -> ok
    if base is None:
        return "new", (":warning:" if value > 0 else ":white_check_mark:"), False, False
    change = int(value) - int(base)
    if change > 0:
        return f"+{change}", ":x:", True, False
    if value > 0:
        return ("0" if change == 0 else str(change)), ":warning:", False, False
    return ("0" if change == 0 else str(change)), ":white_check_mark:", False, False


def _cell(value, kind, detail) -> str:
    text = f"{value:.3g}" if kind == "growth" else (f"{value:.1f}" if kind == "score" else str(int(value)))
    # link a positive count (or a growth failure) to its CSV, if we know the repo URL
    linkable = (kind == "count" and value) or (kind == "growth" and value <= 1e-6)
    if URL_BASE and detail and linkable:
        return f"[{text}]({URL_BASE}/{detail})"
    return text


def _table(rows, current: dict, base: dict):
    lines, regressions, warnings, pending = [], 0, 0, 0
    fatal = False
    for label, key, kind, group, detail in rows:
        value = current.get(key)
        is_pending = value is None or group in RUNNING
        if is_pending:
            lines.append(f"| {label} | _running_ | | :hourglass_flowing_sand: |")
            pending += 1
            continue
        delta, icon, regression, row_fatal = _icon(value, base.get(key), kind)
        fatal = fatal or row_fatal or (key == "dup_keys" and value > 0)
        regressions += regression
        warnings += icon == ":warning:"
        lines.append(f"| {label} | {_cell(value, kind, detail)} | {delta} | {icon} |")
    return lines, regressions, warnings, pending, fatal


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
    current = _metrics(RESULTS)
    base = _metrics(Path(BASE_DIR)) if have_base else {}

    st_tbl, st_reg, st_warn, st_pend, fatal = _table(STRUCTURAL_ROWS, current, base)
    rp_tbl, rp_reg, rp_warn, rp_pend, _ = _table(REPORT_ROWS, current, base)
    mb_tbl, mb_reg, mb_warn, mb_pend, _ = _table(MB_ROWS, current, base)

    regressions = st_reg + rp_reg + mb_reg
    warnings = st_warn + rp_warn + mb_warn
    pending = st_pend + rp_pend + mb_pend

    if fatal:
        verdict = ":x: **Merge blocked: the model cannot be loaded or cannot grow.** See the Structural checks table."
    elif regressions:
        extra = f" ({pending} check(s) still running)" if pending else ""
        verdict = f":x: **{regressions} regression(s) vs `{BASE_REF}`** (this pull request increased a finding count){extra}. Review the :x: rows."
    elif pending:
        verdict = f":hourglass_flowing_sand: **{pending} check(s) still running.** The rest are unchanged vs `{BASE_REF}`."
    elif not have_base:
        verdict = ":information_source: First run for this comparison; no target-branch baseline yet."
    elif warnings:
        verdict = f":warning: **{warnings} pre-existing finding(s), no regressions vs `{BASE_REF}`.** Non-blocking."
    else:
        verdict = f":white_check_mark: **All checks clean, no regressions vs `{BASE_REF}`.**"

    head = f"| Check | Result | &Delta; vs `{BASE_REF}` | |"
    sep = "| --- | ---: | ---: | :---: |"
    lines = [
        "## Model quality report",
        "",
        verdict,
        "",
        "### Structural checks",
        "_Duplicate keys (model unloadable) and no growth block the merge; the other rows are non-blocking._",
        "",
        head, sep, *st_tbl,
        "",
        "### Model QC reports",
        "",
        head, sep, *rp_tbl,
        "",
        "### MACAW and mass/charge balance",
        "",
        head, sep, *mb_tbl,
        "",
        "### Gene essentiality (Hart 2015)",
        "",
        _gene_essentiality_section(),
        "",
        ":x: = a count rose vs the target branch (regression) &middot; "
        ":warning: = a pre-existing non-zero finding (non-blocking) &middot; "
        ":hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.",
    ]
    if COMMIT_SHA:
        lines += ["", f"Results for commit {COMMIT_SHA[:7]}."]

    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
