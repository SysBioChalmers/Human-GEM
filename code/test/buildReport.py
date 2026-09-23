"""Build the model-quality report: a full detail file and a condensed pull-request comment.

data/testResults/model_qc_summary.md is the full, per-check breakdown (every row,
MEMOTE's per-section and per-test scores, everything) -- committed, so it is always
browsable on GitHub regardless of size. data/testResults/model_qc_comment.md is what
actually gets posted as the pull-request comment: one row per section (not per check),
naming only the checks that need attention (a new regression, or a pre-existing,
non-blocking finding) and linking to the full file for everything else. A check that is
clean contributes only to its section's count, never its own line in the comment.

Groups still being computed on this run are passed in the RUNNING_GROUPS environment
variable and show as *running*; the workflow calls this once with "all" before anything
has run, once with "memote" while the slow MEMOTE snapshot is still going, and once with
nothing when everything is in. No stamp files are involved.

Icon rule (per check, both files):
  * growth: white_check_mark if the model grows, x if it cannot (blocks the merge).
  * MEMOTE score: warning if the score dropped versus the target branch, else
    white_check_mark (a non-zero score is good).
  * every other (count) metric: x if the count rose versus the target branch (a
    regression this pull request introduced), warning if the count is non-zero
    (a pre-existing finding, non-blocking), white_check_mark if it is zero.
  * hourglass: the group is still running on this pull request.

Only two conditions fail the build: the model cannot load (duplicate `!!omap` keys) or
cannot grow. Everything else is reported but does not block; a red x just flags a
regression for review.

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
FULL_MD = RESULTS / "model_qc_summary.md"
COMMENT_MD = RESULTS / "model_qc_comment.md"
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
    repo URL is known (in CI); plain text when run locally. Full-detail file only --
    the condensed comment links whole sections to the full file instead (see module
    docstring)."""
    return f"[{label}]({URL_BASE}/README.md#{_slug(label)})" if URL_BASE else label


# (label, key, kind, group, detail_file)
# Structural gates and the model-QC reports share one section: the split between them
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
TASK_CHECKS = [
    ("YAML round-trip (cobrapy)", "roundtrip_cobra"),
    ("YAML round-trip (RAVEN)", "roundtrip_raven"),
    ("YAML lint", "yamllint"),
    ("Essential metabolic tasks", "tasks_essential"),
    ("Verification metabolic tasks", "tasks_verification"),
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


def _compute_rows(rows_spec, current: dict, base: dict) -> list[dict]:
    """Evaluate every (label, key, kind, group, detail) row once. Both the full-detail
    table and the condensed comment's section row/callouts are built from this list,
    so the two can never disagree with each other."""
    computed = []
    for label, key, kind, group, detail in rows_spec:
        value = current.get(key)
        if value is None or group in RUNNING:
            computed.append({"label": label, "kind": kind, "detail": detail, "key": key, "pending": True,
                              "icon": ":hourglass_flowing_sand:"})
            continue
        delta, icon, regression, fatal = _icon(value, base.get(key), kind)
        fatal = fatal or (key == "dup_keys" and value > 0)
        computed.append({
            "label": label, "kind": kind, "detail": detail, "key": key, "pending": False,
            "value": value, "delta": delta, "icon": icon, "regression": regression, "fatal": fatal,
        })
    return computed


def _full_table(computed: list[dict]) -> list[str]:
    """Every row, unabbreviated. Full-detail file only."""
    lines = []
    for r in computed:
        if r["pending"]:
            lines.append(f"| {_labelled(r['label'])} | _running_ | | :hourglass_flowing_sand: |")
        else:
            lines.append(f"| {_labelled(r['label'])} | {_cell(r['value'], r['kind'], r['detail'])} | "
                          f"{r['delta']} | {r['icon']} |")
    return lines


def _section_summary(computed: list[dict]) -> dict:
    """One section's condensed-comment row, plus the specific non-clean rows to name
    beneath it. A clean row contributes only to the section's count."""
    regressions = [r for r in computed if not r["pending"] and r.get("regression")]
    warnings = [r for r in computed if not r["pending"] and r["icon"] == ":warning:"]
    pending = [r for r in computed if r["pending"]]
    fatal = any(r.get("fatal") for r in computed)
    if fatal:
        status = ":x: **blocked**"
    elif regressions:
        status = f":x: **{len(regressions)}** regression(s)"
    elif pending:
        status = f":hourglass_flowing_sand: **{len(pending)}** running"
    elif warnings:
        status = f":warning: **{len(warnings)}** pre-existing"
    else:
        status = ":white_check_mark: all clean"
    return {"status": status, "fatal": fatal, "regressions": regressions, "warnings": warnings, "pending": pending}


def _callout(rows: list[dict], *, verb: str) -> str:
    """One compact line naming specific rows (e.g. 'Pre-existing, non-blocking: Cross-refs
    inconsistent across compartments 3, ...'). Only called for rows that need naming --
    a clean row is already fully accounted for by the section row's count."""
    if not rows:
        return ""
    items = [f"{r['label']} {_cell(r['value'], r['kind'], r['detail'])}" for r in rows]
    return f"{verb}: " + ", ".join(items) + "."


def _task_rows() -> list[dict]:
    """Round-trip, YAML lint and metabolic-task pass/fail from the shared qc_status.tsv
    the workflow writes. That file is committed, so it is present at checkout with
    stale values from a previous run; it is only refreshed once the checks in the
    early "checks" phase have run. While that phase is still going ("checks" in
    RUNNING) every row shows as running rather than the stale committed value; a
    missing key likewise means the check has not finished yet.

    Every one of these is a merge gate (see the README), unlike most of the model
    checks, so a failure here is always named, never folded into a bare count.
    """
    pending_group = "checks" in RUNNING
    status = {} if pending_group else _status_map(RESULTS)
    rows = []
    for label, name in TASK_CHECKS:
        val = status.get(name, "")
        if not val:
            rows.append({"label": label, "pending": True, "icon": ":hourglass_flowing_sand:", "result": "_running_"})
        elif "/" in val:  # tasks: "failed/total"
            failed, total = val.split("/")[:2]
            ok = int(failed) == 0
            rows.append({"label": label, "pending": False, "ok": ok,
                         "result": f"{total} passed" if ok else f"{failed} failed",
                         "icon": ":white_check_mark:" if ok else ":x:"})
        else:  # round-trip / lint: pass|fail
            ok = val.lower() == "pass"
            rows.append({"label": label, "pending": False, "ok": ok, "result": val,
                         "icon": ":white_check_mark:" if ok else ":x:"})
    return rows


def _task_summary(rows: list[dict]) -> dict:
    failed = [r for r in rows if not r["pending"] and not r["ok"]]
    pending = [r for r in rows if r["pending"]]
    if failed:
        status = f":x: **{len(failed)}** failed"
    elif pending:
        status = f":hourglass_flowing_sand: **{len(pending)}** running"
    else:
        status = ":white_check_mark: all pass"
    return {"status": status, "failed": failed, "pending": pending}


def _memote_data(current: Path, base: Path | None) -> dict:
    """Parsed MEMOTE state, shared by the full-file section and the condensed row."""
    if "memote" in RUNNING:
        return {"pending": True}
    core = _memote_meta(current, MEMOTE_CORE)
    if core is None:
        return {"pending": True}
    total, mode, sections, detailed = core
    base_ok = bool(base and base.exists())
    b = _memote_meta(base, MEMOTE_CORE) if base_ok else None
    b_total, b_sections = (b[0], b[2]) if b else (None, {})
    full = _memote_meta(current, MEMOTE_FULL)
    bf = _memote_meta(base, MEMOTE_FULL) if base_ok else None
    return {
        "pending": False, "total": total, "mode": mode, "sections": sections, "detailed": detailed,
        "b_total": b_total, "b_sections": b_sections, "full": full, "bf_total": bf[0] if bf else None,
    }


def _memote_full_section(data: dict) -> str:
    if data["pending"]:
        return "_running_ &middot; :hourglass_flowing_sand:"
    lines = [f"**Total score: {data['total']:.1f}%** ({data['mode']}) &nbsp; "
             f"{_score_delta(data['total'], data['b_total'])}".rstrip(), ""]
    if data["sections"]:
        lines += ["| Section | Score | &Delta; vs base |", "| --- | ---: | ---: |"]
        lines += [f"| {sec} | {sc:.1f}% | {_score_delta(sc, data['b_sections'].get(sec))} |"
                  for sec, sc in data["sections"].items()]
    if data["detailed"]:
        lines += ["", "<details><summary>Per-test scores</summary>", "",
                  "| Section | Test | Score |", "| --- | --- | ---: |"]
        lines += [f"| {s} | {t} | {sc}% |" for s, t, sc in data["detailed"]]
        lines += ["", "</details>"]
    if data["full"] is not None:
        lines += ["", f"**Full suite: {data['full'][0]:.1f}%** &nbsp; "
                  f"{_score_delta(data['full'][0], data['bf_total'])} "
                  "&middot; _from the last_ `/run memote`.".rstrip()]
    else:
        lines += ["", "_Full suite not run for this commit; comment_ `/run memote` _to add it._"]
    return "\n".join(lines)


def _memote_summary(data: dict) -> dict:
    if data["pending"]:
        return {"status": ":hourglass_flowing_sand: running", "dropped": False}
    dropped = data["b_total"] is not None and (data["total"] - data["b_total"]) < -0.05
    icon = ":warning:" if dropped else ":white_check_mark:"
    return {"status": f"{icon} **{data['total']:.1f}%**", "dropped": dropped}


def _gene_essentiality_full_section() -> str:
    # Gene essentiality takes hours and is not run on every pull request. Its result
    # file (gene-essential.csv) is committed and persists across pull requests, so it
    # would be stale here - the result is shown in its own comment when run instead.
    return ("_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` "
            "_to run it on this pull request; the result posts as its own comment._")


def _gene_essentiality_summary() -> dict:
    has_run = (RESULTS / "gene-essential.csv").exists()
    return {"status": "see the gene-essentiality comment" if has_run else "_not run_"}


def _gates_line(current: dict, base: dict) -> str:
    """One-line status of the two merge gates (duplicate keys, growth), appended to the
    verdict so they are visible without following any link."""
    if "checks" in RUNNING:
        return ""
    parts = []
    for label, key, kind in (("duplicate keys", "dup_keys", "count"), ("growth", "growth", "growth")):
        value = current.get(key)
        if value is None:
            return ""
        _, icon, _, _ = _icon(value, base.get(key), kind)
        text = f"{value:.3g}" if kind == "growth" else str(int(value))
        parts.append(f"{label} **{text}** {icon}")
    return " Gates: " + ", ".join(parts) + "."


def main() -> int:
    have_base = bool(BASE_DIR) and Path(BASE_DIR).exists()
    current = _metrics(RESULTS)
    base = _metrics(Path(BASE_DIR)) if have_base else {}

    network_rows = _compute_rows(MODEL_ROWS + MB_ROWS, current, base)
    network_summary = _section_summary(network_rows)
    task_rows = _task_rows()
    task_summary = _task_summary(task_rows)
    memote_data = _memote_data(RESULTS, Path(BASE_DIR) if BASE_DIR else None)
    memote_summary = _memote_summary(memote_data)
    gene_summary = _gene_essentiality_summary()

    fatal = network_summary["fatal"] or bool(task_summary["failed"])
    regressions = len(network_summary["regressions"])
    pending = len(network_summary["pending"]) + len(task_summary["pending"]) + (1 if memote_data["pending"] else 0)
    warnings = len(network_summary["warnings"])

    if fatal:
        verdict = ":x: **Merge blocked: the model cannot be loaded or cannot grow, or a build gate failed.**"
    elif regressions:
        extra = f" ({pending} check(s) still running)" if pending else ""
        verdict = f":x: **{regressions}** regression(s) vs `{BASE_REF}`{extra}. Review the row(s) below."
    elif pending:
        verdict = f":hourglass_flowing_sand: **{pending}** check(s) still running. The rest are unchanged vs `{BASE_REF}`."
    elif not have_base:
        verdict = ":information_source: First run for this comparison; no target-branch baseline yet."
    elif warnings:
        verdict = f":warning: **{warnings}** pre-existing finding(s), no regressions vs `{BASE_REF}`. Non-blocking."
    else:
        verdict = f":white_check_mark: All checks clean, no regressions vs `{BASE_REF}`."
    verdict += _gates_line(current, base)

    # --- condensed comment: one row per section, only non-clean checks named ---
    full_url = f"{URL_BASE}/model_qc_summary.md" if URL_BASE else ""
    section_table = [
        "| Section | Status |",
        "| --- | --- |",
        f"| Model &amp; network checks ({len(network_rows)}) | {network_summary['status']} |",
        f"| Model file &amp; metabolic tasks ({len(task_rows)} gates) | {task_summary['status']} |",
        f"| MEMOTE | {memote_summary['status']} |",
        f"| Gene essentiality | {gene_summary['status']} |",
    ]
    callouts = [
        _callout(network_summary["regressions"], verb="Regression(s)"),
        (f"Gate failure(s): " + ", ".join(f"{r['label']} ({r['result']})" for r in task_summary["failed"]) + "."
         if task_summary["failed"] else ""),
        _callout(network_summary["warnings"], verb="Pre-existing, non-blocking") if not regressions else "",
    ]
    callouts = [c for c in callouts if c]

    comment_lines = [
        "## Model quality report",
        "",
        verdict,
        "",
        *section_table,
        "",
        *callouts,
    ]
    if callouts:
        comment_lines.append("")
    if full_url:
        comment_lines.append(
            f"[Full report]({full_url}) -- every check, MEMOTE's per-test breakdown, "
            f"and how to run gene essentiality."
        )
    if COMMIT_SHA:
        comment_lines += ["", f"Results for commit {COMMIT_SHA[:7]}."]
    COMMENT_MD.write_text("\n".join(comment_lines) + "\n", encoding="utf-8")

    # --- full detail file: every row, unabbreviated ---
    head = f"| Check | Result | &Delta; vs `{BASE_REF}` | |"
    sep = "| --- | ---: | ---: | :---: |"
    task_head = "| Check | Result | |"
    task_sep = "| --- | ---: | :---: |"
    full_lines = [
        "## Model quality report -- full detail",
        "",
        "_This is the full, per-check breakdown behind the pull-request comment's summary table._ "
        + (f"_Row names link to their explanation in the [testResults README]({URL_BASE}/README.md)._"
           if URL_BASE else ""),
        "",
        "### Model & network checks",
        "_Duplicate keys (model unloadable) and no growth block the merge; every other row "
        "is a non-blocking report._",
        "",
        head, sep, *_full_table(network_rows),
        "",
        "### Model file and metabolic tasks",
        "",
        task_head, task_sep,
        *[f"| {_labelled(r['label'])} | {r['result']} | {r['icon']} |" for r in task_rows],
        "",
        f"### {_labelled('MEMOTE')}",
        "",
        _memote_full_section(memote_data),
        "",
        "_The score above is the fast core subset. Comment_ `/run memote` "
        "_to run the full suite on this pull request; the score updates here when it finishes._",
        "",
        f"### {_labelled('Gene essentiality (Hart 2015)')}",
        "",
        _gene_essentiality_full_section(),
        "",
        ":x: = a count rose vs the target branch (regression) &middot; "
        ":warning: = a pre-existing non-zero finding (non-blocking) &middot; "
        ":hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.",
    ]
    FULL_MD.write_text("\n".join(full_lines) + "\n", encoding="utf-8")

    print("\n".join(comment_lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
