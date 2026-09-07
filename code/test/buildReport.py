"""Build one model-quality report for the pull-request comment.

Turns the result files under data/testResults/ into a single comment that leads
with a one-line verdict and then the status tables (model checks, MACAW/balance,
model-file/metabolic tasks, MEMOTE, gene essentiality). Each check name links to
its explanation in this folder's README. Groups still being computed on this run
are passed in the RUNNING_GROUPS environment variable and their rows show as
*running* (hourglass);
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


def _slug(label: str) -> str:
    """GitHub heading-anchor slug for a table label. Mirrors GitHub's algorithm
    (lowercase, drop punctuation, spaces to hyphens) so a label links to the
    same-named section in this folder's README."""
    s = label.lower().replace("`", "")
    s = re.sub(r"[^\w\s-]", "", s)
    return s.strip().replace(" ", "-")


def _labelled(label: str) -> str:
    """The test name, linked to its explanation in the testResults README when the
    repo URL is known (in CI); plain text when run locally."""
    return f"[{label}]({URL_BASE}/README.md#{_slug(label)})" if URL_BASE else label

# (label, key, kind, group, detail_file)
# Structural gates and the model-QC reports share one table: the split between them
# was arbitrary (growth next to unused genes). The two gates (duplicate keys, growth)
# lead the table; every other row is a non-blocking report. Each label links to the
# matching section in this folder's README (see _labelled).
MODEL_ROWS = [
    ("Duplicate `!!omap` keys", "dup_keys", "count", "checks", "qc_duplicate_keys.csv"),
    ("Growth (biomass producible)", "growth", "growth", "checks", "qc_growth_blockers.csv"),
    ("Reactions with no metabolites", "empty_rxn", "count", "checks", "qc_empty_reactions.csv"),
    ("Model / annotation-table inconsistencies", "annot_consistency", "count", "checks",
     "qc_annotation_consistency.csv"),
    ("Removed reactions or metabolites not deprecated", "removed_not_deprecated", "count", "checks",
     "qc_deprecation_completeness.csv"),
    ("Metabolites missing formula", "missing_formula", "count", "checks", "qc_metabolite_completeness.csv"),
    ("Metabolites missing charge", "missing_charge", "count", "checks", "qc_metabolite_completeness.csv"),
    ("Reaction bound / GPR issues", "reaction_issues", "count", "checks", "qc_reaction_sanity.csv"),
    ("Exact-duplicate reaction groups", "dup_reactions", "count", "checks", "qc_duplicate_reactions.csv"),
    ("Unused metabolites", "unused_met", "count", "checks", "qc_unused_entities.csv"),
    ("Unused genes", "unused_gene", "count", "checks", "qc_unused_entities.csv"),
    ("Malformed cross-references", "malformed", "count", "checks", "qc_annotation_issues.csv"),
    ("Cross-refs inconsistent across compartments", "inconsistent", "count", "checks", "qc_annotation_issues.csv"),
]
MB_ROWS = [
    ("Reactions flagged by MACAW dead-end test", "dead_end", "count", "macaw", "macaw_results.csv"),
    ("Reactions flagged as MACAW duplicates", "duplicates", "count", "macaw", "macaw_results.csv"),
    ("Mass-imbalanced reactions", "mass_imbalance", "count", "macaw", "balance_results.csv"),
    ("Charge-imbalanced reactions", "charge_imbalance", "count", "macaw", "balance_results.csv"),
    ("Structure vs formula/charge inconsistencies", "structure_inconsistent", "count", "macaw",
     "qc_structure_consistency.csv"),
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


def _status_map(directory: Path) -> dict:
    """The combined qc_status.tsv as {check: result} ({} if absent). Holds the
    one-line checks (round-trip, yamllint, metabolic tasks) and the growth value."""
    path = directory / "qc_status.tsv"
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] != "check":
            out[parts[0]] = parts[1]
    return out


def _growth(directory: Path) -> float | None:
    try:
        return float(_status_map(directory)["growth"])
    except (KeyError, ValueError):
        return None


# memote_score.md is split into these two sections (see memoteSnapshot.py); each is
# parsed and compared only against the same section on the base branch, so a subset
# score is never diffed against a full-suite score.
MEMOTE_CORE = "Core subset"
MEMOTE_FULL = "Full suite"


def _memote_meta(directory: Path, title: str):
    """Parse one section of memote_score.md ->
    (total, mode, {section: score}, [(section, test, score)]), or None if that section
    is absent or not yet computed (a placeholder with no total)."""
    path = directory / "memote_score.md"
    if not path.exists():
        return None
    full = path.read_text(encoding="utf-8")
    m = re.search(rf"^## {re.escape(title)}\s*$(.*?)(?=^## |\Z)", full, re.M | re.S)
    if m:
        text = m.group(1)
    elif title == MEMOTE_CORE:
        text = full  # back-compat: an older single-section file is the core subset
    else:
        return None
    total = re.search(r"Total score:\s*([\d.]+)\s*%", text)
    if not total:
        return None
    mode = re.search(r"Mode:\s*(.+?)\.", text)
    sections = {m.group(1): float(m.group(2))
                for m in re.finditer(r"^\| (\w+) \| ([\d.]+)% \|$", text, re.M)}
    detailed = [(s, t, sc) for s, t, sc in re.findall(r"^\| (.+?) \| (.+?) \| ([\d.]+)% \|$", text, re.M)]
    return (float(total.group(1)), mode.group(1) if mode else "", sections, detailed)


def _score_delta(cur, base) -> str:
    if cur is None or base is None:
        return ""
    d = cur - base
    if abs(d) < 0.05:
        return "0"
    return f"{d:+.1f} {':warning:' if d < 0 else ':white_check_mark:'}"


def _memote_section(current: Path, base: Path | None) -> str:
    if "memote" in RUNNING:
        return "_running_ &middot; :hourglass_flowing_sand:"
    core = _memote_meta(current, MEMOTE_CORE)
    if core is None:
        return "_running_ &middot; :hourglass_flowing_sand:"
    total, mode, sections, detailed = core
    base_ok = base and base.exists()
    b = _memote_meta(base, MEMOTE_CORE) if base_ok else None
    b_total, b_sections = (b[0], b[2]) if b else (None, {})
    lines = [f"**Total score: {total:.1f}%** ({mode}) &nbsp; {_score_delta(total, b_total)}".rstrip(), ""]
    if sections:
        lines += ["| Section | Score | &Delta; vs base |", "| --- | ---: | ---: |"]
        lines += [f"| {sec} | {sc:.1f}% | {_score_delta(sc, b_sections.get(sec))} |"
                  for sec, sc in sections.items()]
    if detailed:
        lines += ["", "<details><summary>Per-test scores</summary>", "",
                  "| Section | Test | Score |", "| --- | --- | ---: |"]
        lines += [f"| {s} | {t} | {sc}% |" for s, t, sc in detailed]
        lines += ["", "</details>"]

    # Full suite: shown only if a /run memote result is committed. Compared to the
    # full-suite section on the base branch, never to the subset score above.
    full = _memote_meta(current, MEMOTE_FULL)
    if full is not None:
        bf = _memote_meta(base, MEMOTE_FULL) if base_ok else None
        bf_total = bf[0] if bf else None
        lines += ["", f"**Full suite: {full[0]:.1f}%** &nbsp; {_score_delta(full[0], bf_total)} "
                  "&middot; _from the last_ `/run memote`.".rstrip()]
    else:
        lines += ["", "_Full suite not run for this commit; comment_ `/run memote` _to add it._"]
    return "\n".join(lines)


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
        "removed_not_deprecated": _count_csv(directory / "qc_deprecation_completeness.csv"),
        "growth": _growth(directory),
        "missing_formula": _count_csv(completeness, lambda r: r.get("missing_formula") == "yes"),
        "missing_charge": _count_csv(completeness, lambda r: r.get("missing_charge") == "yes"),
        "reaction_issues": _count_csv(directory / "qc_reaction_sanity.csv"),
        "dup_reactions": _distinct_csv(directory / "qc_duplicate_reactions.csv", "group"),
        "unused_met": _count_csv(unused, lambda r: r.get("kind") == "metabolite"),
        "unused_gene": _count_csv(unused, lambda r: r.get("kind") == "gene"),
        "malformed": _count_csv(annotation, lambda r: r.get("issue", "").startswith("malformed")),
        "inconsistent": _count_csv(annotation, lambda r: r.get("issue", "").startswith("inconsistent")),
        "dead_end": _count_csv(macaw, lambda r: r.get("dead_end_test", "") not in ("ok", "")),
        "duplicates": _count_csv(macaw, lambda r: any(r.get(c, "") not in ("ok", "N/A", "") for c in _DUP_COLS)),
        "mass_imbalance": _count_csv(balance, lambda r: r.get("mass_imbalance", "").strip() != ""),
        "charge_imbalance": _count_csv(balance, lambda r: r.get("charge_imbalance", "").strip() != ""),
        # the CSV lists only the inconsistent metabolites, so its row count is the metric
        "structure_inconsistent": _count_csv(directory / "qc_structure_consistency.csv"),
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
            lines.append(f"| {_labelled(label)} | _running_ | | :hourglass_flowing_sand: |")
            pending += 1
            continue
        delta, icon, regression, row_fatal = _icon(value, base.get(key), kind)
        fatal = fatal or row_fatal or (key == "dup_keys" and value > 0)
        regressions += regression
        warnings += icon == ":warning:"
        lines.append(f"| {_labelled(label)} | {_cell(value, kind, detail)} | {delta} | {icon} |")
    return lines, regressions, warnings, pending, fatal


def _model_integrity_section() -> str:
    """Round-trip, YAML lint and metabolic-task pass/fail from the shared qc_status.tsv
    the workflow writes. That file is committed, so it is present at checkout with
    stale values from a previous run; it is only refreshed once the checks in the
    early "checks" phase have run. While that phase is still going ("checks" in
    RUNNING) show every row as running rather than the stale committed value; a
    missing key likewise means the check has not finished yet."""
    checks = [
        ("YAML round-trip (cobrapy)", "roundtrip_cobra"),
        ("YAML round-trip (RAVEN)", "roundtrip_raven"),
        ("YAML lint", "yamllint"),
        ("Essential metabolic tasks", "tasks_essential"),
        ("Verification metabolic tasks", "tasks_verification"),
    ]
    out = ["| Check | Result | |", "| --- | ---: | :---: |"]
    pending = "checks" in RUNNING
    status = {} if pending else _status_map(RESULTS)
    for label, name in checks:
        val = status.get(name, "")
        if not val:
            out.append(f"| {_labelled(label)} | _running_ | :hourglass_flowing_sand: |")
        elif "/" in val:                       # tasks: "failed/total"
            failed, total = val.split("/")[:2]
            ok = int(failed) == 0
            out.append(f"| {_labelled(label)} | {total + ' passed' if ok else failed + ' failed'} | "
                       f"{':white_check_mark:' if ok else ':x:'} |")
        else:                                  # round-trip / lint: pass|fail
            ok = val.lower() == "pass"
            out.append(f"| {_labelled(label)} | {val} | {':white_check_mark:' if ok else ':x:'} |")
    return "\n".join(out)


def _gene_essentiality_section() -> str:
    # Gene essentiality takes hours and is not run on every pull request. Its result
    # file (gene-essential.csv) is committed and persists across pull requests, so it
    # would be stale here - the result is shown in its own comment when run instead.
    return ("_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` "
            "_to run it on this pull request; the result posts as its own comment._")


def main() -> int:
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    current = _metrics(RESULTS)
    base = _metrics(Path(BASE_DIR)) if have_base else {}

    md_tbl, md_reg, md_warn, md_pend, fatal = _table(MODEL_ROWS, current, base)
    mb_tbl, mb_reg, mb_warn, mb_pend, _ = _table(MB_ROWS, current, base)

    regressions = md_reg + mb_reg
    warnings = md_warn + mb_warn
    pending = md_pend + mb_pend

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
        "_Each check name links to its explanation in the "
        f"[testResults README]({URL_BASE}/README.md)._" if URL_BASE else "",
        "",
        "### Model checks",
        "_Duplicate keys (model unloadable) and no growth block the merge; every other row "
        "is a non-blocking report._",
        "",
        head, sep, *md_tbl,
        "",
        "### MACAW and mass/charge balance",
        "",
        head, sep, *mb_tbl,
        "",
        "### Model file and metabolic tasks",
        "",
        _model_integrity_section(),
        "",
        f"### {_labelled('MEMOTE')}",
        "",
        _memote_section(RESULTS, Path(BASE_DIR) if BASE_DIR else None),
        "",
        "_The score above is the fast core subset. Comment_ `/run memote` "
        "_to run the full suite on this pull request; the score updates here when it finishes._",
        "",
        f"### {_labelled('Gene essentiality (Hart 2015)')}",
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
