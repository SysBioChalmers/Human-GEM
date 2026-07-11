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

MODEL_FILE = "model/Human-GEM.yml"
RESULT_JSON = "memote_result.json"  # repo root -> uploaded as artifact, not committed
SCORE_MD = "data/testResults/memote_score.md"

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
    lines = ["# MEMOTE snapshot", "", f"Mode: {kind}."]
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
            lines += ["", "| Section | Score |", "| --- | --- |"]
            lines += [f"| {name} | {value * 100:.1f}% |" for name, value in rows]

    with open(SCORE_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
