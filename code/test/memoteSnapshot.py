"""Run the MEMOTE test suite and record the total score for tracking.

MEMOTE only computes a total score when it builds a report; the raw `memote run`
result JSON contains per-test results but no score. So this runs MEMOTE in-process
through its Python API and then builds the snapshot report to obtain the score,
rather than shelling out to `memote run` (which would leave us without a score).

Two modes, chosen by the MEMOTE_SUBSET environment variable:

  * Core subset (MEMOTE_SUBSET set): skips the three flux-variability / loopless
    tests that dominate runtime on a genome-scale model, so it finishes within a
    per-pull-request budget. This is what runs on every pull request.
  * Full suite (MEMOTE_SUBSET unset): the complete suite. Much slower (the skipped
    tests do FVA or a loopless MILP over every reaction), so it is reserved for
    pull requests that target the main branch.

Writes the total score to data/testResults/memote_score.md (diff-friendly) and the
scored result JSON to memote_result.json in the repository root, which the workflow
uploads as a build artifact (it is not committed, to avoid bloating the repository).
memote_score.md keeps a "Core subset" and a "Full suite" section; each run rewrites
only its own section, so a routine subset run never overwrites a committed full-suite
score (see _write_section).

Before exporting the SBML, the model is enriched with the database cross-references
from the annotation tables (see annotateModel.py) so MEMOTE's annotation tests score
against the identifiers Human-GEM actually carries. The enriched model exists only in
memory for the temporary SBML; it is never committed.

Set GRB_LICENSE_FILE (a full Gurobi licence) to run with Gurobi; the genome-scale
MILPs are impractical with GLPK. Without it the script falls back to the default
solver.

Usage:
    python code/test/memoteSnapshot.py
"""

import json
import os
import sys
import tempfile

import cobra
import memote.suite.api as api
from memote.suite.reporting import ReportConfiguration, SnapshotReport

import annotateModel

MODEL_FILE = "model/Human-GEM.yml"
RESULT_JSON = "memote_result.json"  # repo root -> uploaded as artifact, not committed
SCORE_MD = "data/testResults/memote_score.md"

# memote_score.md holds two independent sections. The fast core subset runs on every
# pull request and the full suite runs on demand (/run memote); they share the file
# but must not overwrite each other, so each run rewrites only its own section and
# leaves the other intact. A routine subset run therefore never erases a previously
# committed full-suite score, and buildReport compares each section only against the
# same section on the base branch (never subset vs full).
CORE_TITLE = "Core subset"
FULL_TITLE = "Full suite"

# The tests that dominate MEMOTE runtime on a genome-scale model. Two groups:
#  * consistency: MILP / flux-variability / per-metabolite optimisation over the
#    whole model (stoichiometric consistency, energy cycles, blocked reactions,
#    open-bound producibility, ...). This module alone took ~32 min in CI.
#  * matrix: rank / null-space of the stoichiometric matrix, which is O(n^3) and
#    intractable on a genome-scale matrix.
# The core subset skips these so it finishes within a per-pull-request budget; the
# full suite (pull requests to main) runs them.
SLOW_TESTS = [
    # consistency (test_consistency.py)
    "test_stoichiometric_consistency",
    "test_unconserved_metabolites",
    "test_inconsistent_min_stoichiometry",
    "test_detect_energy_generating_cycles",
    "test_find_stoichiometrically_balanced_cycles",
    "test_blocked_reactions",
    "test_find_reactions_unbounded_flux_default_condition",
    "test_find_metabolites_not_produced_with_open_bounds",
    "test_find_metabolites_not_consumed_with_open_bounds",
    # matrix (test_matrix.py)
    "test_number_independent_conservation_relations",
    "test_matrix_rank",
    "test_degrees_of_freedom",
]


def _total_score(scored: dict) -> float | None:
    """Pull the total score (0-1) out of a scored MEMOTE result, tolerating layout."""
    candidates = [scored.get("score")]
    cards = scored.get("cards")
    if isinstance(cards, dict):
        candidates.append(cards.get("score"))
    for card in candidates:
        if isinstance(card, dict) and isinstance(card.get("total_score"), (int, float)):
            return float(card["total_score"])
    if isinstance(scored.get("total_score"), (int, float)):
        return float(scored["total_score"])
    return None


def _section_rows(scored: dict) -> list[tuple[str, float]]:
    """Best-effort per-section scores for the summary table."""
    score = scored.get("score")
    sections = score.get("sections") if isinstance(score, dict) else None
    rows = []
    for section in sections or []:
        name = section.get("section") or section.get("title")
        value = section.get("score")
        if name is not None and isinstance(value, (int, float)):
            rows.append((str(name), float(value)))
    return rows


def _test_metric(scored: dict, test_id: str) -> float | None:
    """The 0-1 metric of a single MEMOTE test (parametrised tests are averaged)."""
    test = (scored.get("tests") or {}).get(test_id)
    if not isinstance(test, dict):
        return None
    metric = test.get("metric")
    if isinstance(metric, (int, float)):
        return float(metric)
    if isinstance(metric, dict):
        values = [v for v in metric.values() if isinstance(v, (int, float))]
        return sum(values) / len(values) if values else None
    return None


def _test_title(scored: dict, test_id: str) -> str:
    test = (scored.get("tests") or {}).get(test_id) or {}
    return str(test.get("title") or test_id)


def _detailed_rows(scored: dict, config) -> list[tuple[str, str, float]]:
    """Best-effort per-test scores grouped by section: (section, test, metric).

    Uses the scoring configuration's section -> cases mapping to place each test,
    and the per-test metric from the scored result. Returns [] if the layout is
    not as expected, so the caller can fall back to the section summary.
    """
    try:
        sections = config["cards"]["scored"]["sections"]
    except (KeyError, TypeError, AttributeError):
        return []
    if not isinstance(sections, dict):
        return []
    rows: list[tuple[str, str, float]] = []
    for section_id, section in sections.items():
        if not isinstance(section, dict):
            continue
        title = str(section.get("title") or section_id)
        for case in section.get("cases") or []:
            metric = _test_metric(scored, case)
            if metric is not None:
                rows.append((title, _test_title(scored, case), metric))
    if rows:
        return rows
    # Fallback: the section -> cases mapping was not where we expected it. List every
    # scored test flat, so the detail is still available even if less tidy.
    tests = scored.get("tests")
    if isinstance(tests, dict):
        for test_id in tests:
            metric = _test_metric(scored, test_id)
            if metric is not None:
                rows.append(("", _test_title(scored, test_id), metric))
    return rows


def _load_sections(path: str) -> dict:
    """Existing memote_score.md as {section_title: body_text}. Empty if absent."""
    sections: dict[str, str] = {}
    if not os.path.exists(path):
        return sections
    current, buf = None, []
    for line in open(path, encoding="utf-8").read().splitlines():
        if line.startswith("## "):
            if current is not None:
                sections[current] = "\n".join(buf).strip("\n")
            current, buf = line[3:].strip(), []
        elif current is not None:
            buf.append(line)
    if current is not None:
        sections[current] = "\n".join(buf).strip("\n")
    return sections


def _placeholder(title: str) -> str:
    if title == FULL_TITLE:
        return "_Not run for this commit. Comment_ `/run memote` _to populate this section._"
    return "_Not yet computed for this commit._"


def _write_section(this_title: str, body: str) -> None:
    """Rewrite only this run's section, preserving the other one (or a placeholder)."""
    sections = _load_sections(SCORE_MD)
    sections[this_title] = body
    out = ["# MEMOTE snapshot", ""]
    for title in (CORE_TITLE, FULL_TITLE):
        out += [f"## {title}", "", sections.get(title) or _placeholder(title), ""]
    with open(SCORE_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out).rstrip() + "\n")


def main() -> int:
    subset = bool(os.environ.get("MEMOTE_SUBSET"))
    skip = SLOW_TESTS if subset else None
    kind = "core subset" if subset else "full suite"

    # Use Gurobi when a full licence is configured; the genome-scale consistency
    # and FVA MILPs are impractical with GLPK.
    if os.environ.get("GRB_LICENSE_FILE"):
        cobra.Configuration().solver = "gurobi"

    # memote reads an SBML model, so convert the canonical YAML model to a
    # temporary SBML file first (memote fails on a .yml directly).
    model = cobra.io.load_yaml_model(MODEL_FILE)

    # The YAML model has only ids and names; attach the database cross-references
    # from the annotation tables so MEMOTE's annotation tests see them. This mutates
    # the in-memory model only - the SBML written below is temporary and the enriched
    # model is never committed.
    counts = annotateModel.enrich(model)
    print(f"Annotated for MEMOTE (temporary): {counts['metabolites']} metabolites, "
          f"{counts['reactions']} reactions, {counts['genes']} genes.", flush=True)

    sbml_path = os.path.join(tempfile.gettempdir(), "human-gem.xml")
    cobra.io.write_sbml_model(model, sbml_path)

    model_obj, sbml_ver, _ = api.validate_model(sbml_path)
    print(f"Running MEMOTE ({kind})"
          + (f", skipping {len(SLOW_TESTS)} slow tests" if skip else ""), flush=True)
    _, result = api.test_model(
        model_obj, sbml_version=sbml_ver, results=True,
        skip=skip, solver_timeout=120,
    )

    # The score is only computed when a report is built: compute_score() runs in
    # the SnapshotReport constructor and writes result["score"]. Build it directly
    # so we do not depend on the report renderer's output format.
    config = ReportConfiguration.load()
    try:
        scored = SnapshotReport(result=result, configuration=config).result
    except Exception as exc:  # noqa: BLE001
        print(f"::warning::MEMOTE scoring failed ({exc}); reporting the raw result.")
        scored = result

    with open(RESULT_JSON, "w", encoding="utf-8") as fh:
        json.dump(scored, fh, default=str)

    print("Scored MEMOTE result top-level keys:", sorted(scored.keys()), flush=True)

    total = _total_score(scored)
    lines = [f"Mode: {kind}."]
    if subset:
        lines.append(f"Skipped (slow) tests: {', '.join(SLOW_TESTS)}.")
    lines.append("")
    if total is None:
        print("::warning::Could not find the MEMOTE total score in the result.")
        lines.append("Total score: unavailable (see workflow log).")
    else:
        pct = total * 100 if total <= 1 else total
        print(f"MEMOTE total score ({kind}): {pct:.1f}%", flush=True)
        lines.append(f"**Total score: {pct:.1f}%**")
        rows = _section_rows(scored)
        if rows:
            lines += ["", "### Section scores", "", "| Section | Score |", "| --- | ---: |"]
            lines += [f"| {name} | {value * 100:.1f}% |" for name, value in rows]

        try:
            detailed = _detailed_rows(scored, config)
        except Exception as exc:  # noqa: BLE001 - detail is optional, never fail on it
            print(f"::warning::Could not build the detailed MEMOTE scores ({exc}).")
            detailed = []
        if detailed:
            lines += ["", "### Detailed scores", "", "| Section | Test | Score |",
                      "| --- | --- | ---: |"]
            lines += [f"| {section} | {test} | {metric * 100:.1f}% |"
                      for section, test, metric in detailed]

    # Rewrite only this run's section (core subset or full suite), keeping the other.
    _write_section(CORE_TITLE if subset else FULL_TITLE, "\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
