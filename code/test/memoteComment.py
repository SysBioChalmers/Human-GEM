"""Write the pull-request comment that reports a finished full MEMOTE run.

Called by memote-full.yml when the full suite (started with ``/run memote``) ends. It
reads the "Full suite" section of data/testResults/memote_score.md, compares its total
with the same section on the base branch, and prints a short Markdown comment: whether
the run finished, the score and how it changed, a link to the full report, and the
per-test scores in a collapsed table.

Usage:
    python code/test/memoteComment.py --result finished --version HASH --base-dir DIR \\
        --report-url URL --run-url URL [--base-ref develop]

``--result`` is how the MEMOTE step ended (finished, timeout, failed, or skipped when it
did not run), and ``--version`` the model version it scored (qcStatus.model_version()
at its start). Both come from the step itself, so the comment describes this run even
if the branch moved on while it ran.
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


NOT_FINISHED = {
    "timeout": "did not finish in time",
    "failed": "failed",
    "failure": "failed",
    "skipped": "did not start",
}


def comment(result: str, version: str, base_dir: Path | None, report_url: str,
            run_url: str, base_ref: str) -> str:
    run = f"[run]({run_url})" if run_url else "run"
    full = buildReport._memote_meta(RESULTS, buildReport.MEMOTE_FULL)
    scored = result == "finished" and full is not None and full[4] == version
    if not scored:
        reason = NOT_FINISHED.get(result, "left no score")
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
    if version != qcStatus.model_version():
        lines += ["", "_The branch changed during the run; this score is for the model "
                  "it started on._"]
    lines += ["", "_The Model QC comment is updated with this score too._"]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--result", default="finished")
    parser.add_argument("--version", default="")
    parser.add_argument("--base-dir", type=Path)
    parser.add_argument("--report-url", default="")
    parser.add_argument("--run-url", default="")
    parser.add_argument("--base-ref", default="the base branch")
    args = parser.parse_args(argv)
    version = args.version or qcStatus.model_version()
    print(comment(args.result, version, args.base_dir, args.report_url, args.run_url,
                  args.base_ref))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
