"""Write the pull-request comment that reports a finished full MEMOTE run.

Called by memote-full.yml when the full suite (started with ``/run memote``) ends. It
reads the "Full suite" section of data/testResults/memote_score.md, compares its total
with the same section on the base branch, and prints a short Markdown comment: whether
the run finished, the score and how it changed, a link to the full report, and the
per-test scores in a collapsed table.

Usage:
    python code/test/memoteComment.py --outcome success --base-dir DIR \\
        --report-url URL --run-url URL [--base-ref develop]

``--outcome`` is the outcome of the MEMOTE step (success, failure, ...). A run that
ends without a full-suite score for the current model version, e.g. because it hit the
time limit, is reported as not finished.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import buildReport
import qcStatus

RESULTS = Path(__file__).resolve().parents[2] / "data" / "testResults"


def _change(total: float, base: float | None, base_ref: str) -> str:
    """How the score compares with the base branch, in a few words."""
    if base is None:
        return f"no full-suite score on `{base_ref}` to compare with"
    delta = total - base
    if abs(delta) < 0.05:
        return f":white_check_mark: unchanged from {base:.1f}% on `{base_ref}`"
    if delta > 0:
        return f"{buildReport.IMPROVED} improved from {base:.1f}% on `{base_ref}` (+{delta:.1f})"
    return f":warning: dropped from {base:.1f}% on `{base_ref}` ({delta:.1f})"


def comment(outcome: str, base_dir: Path | None, report_url: str, run_url: str,
            base_ref: str) -> str:
    full = buildReport._memote_meta(RESULTS, buildReport.MEMOTE_FULL)
    current = bool(full and full[4] and full[4] == qcStatus.model_version())
    run = f"[run]({run_url})" if run_url else "run"

    if not current:
        reason = "failed" if outcome == "failure" else "did not finish in time"
        return (f":x: The **full MEMOTE suite** {reason}, so there is no new score "
                f"(see the {run}).")

    total, _, _, detailed, _ = full
    base = None
    if base_dir and base_dir.exists():
        base_full = buildReport._memote_meta(base_dir, buildReport.MEMOTE_FULL)
        base = base_full[0] if base_full else None
    lines = [
        f":microscope: The **full MEMOTE suite** has finished ({run}).",
        "",
        f"**Score: {total:.1f}%**, {_change(total, base, base_ref)}. "
        f"[Full report]({report_url})",
    ]
    if detailed:
        lines += ["", "<details><summary>Per-test scores</summary>", "",
                  "| Section | Test | Score |", "| --- | --- | ---: |"]
        lines += [f"| {section} | {test} | {score}% |" for section, test, score in detailed]
        lines += ["", "</details>"]
    lines += ["", "_The Model QC comment is updated with this score too._"]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outcome", default="success")
    parser.add_argument("--base-dir", type=Path)
    parser.add_argument("--report-url", default="")
    parser.add_argument("--run-url", default="")
    parser.add_argument("--base-ref", default="the base branch")
    args = parser.parse_args(argv)
    print(comment(args.outcome, args.base_dir, args.report_url, args.run_url, args.base_ref))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
