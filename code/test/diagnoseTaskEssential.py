"""Diagnose why find_task_essential_reactions (gene-essentiality Step 1a) hangs on
some models: time the call per task and report the slow ones.

find_task_essential_reactions processes tasks sequentially (pfba -> candidate
reactions -> FVA over candidates) and sets no solver time limit, so a single task
whose LP is hard for the model can stall the whole step. This script caps each
solve with a solver time limit, logs the wall-clock time of every task (printing a
"starting" line first so a genuine hang still names the culprit), and stops at an
overall budget so it always finishes and reports.

Usage:
    SOLVER_TIMEOUT=30 DIAG_BUDGET=2000 python code/test/diagnoseTaskEssential.py
"""

import os
import sys
import time

import cobra

from raven_toolbox.tasks.check import find_task_essential_reactions
from raven_toolbox.tasks.tasklist import parse_task_list

MODEL_FILE = "model/Human-GEM.yml"
ESSENTIAL_TASKS = "data/metabolicTasks/metabolicTasks_Essential.txt"
SOLVER_TIMEOUT = int(os.environ.get("SOLVER_TIMEOUT", "30"))   # per-solve cap, seconds
DIAG_BUDGET = int(os.environ.get("DIAG_BUDGET", "2000"))       # overall cap, seconds


def main() -> int:
    if os.environ.get("GRB_LICENSE_FILE"):
        cobra.Configuration().solver = "gurobi"

    model = cobra.io.load_yaml_model(MODEL_FILE)
    try:
        model.solver.configuration.timeout = SOLVER_TIMEOUT
        print(f"Solver: {model.solver.interface.__name__}, per-solve timeout {SOLVER_TIMEOUT}s",
              flush=True)
    except Exception as exc:  # noqa: BLE001
        print(f"::warning::could not set solver timeout ({exc}); solves are uncapped", flush=True)

    tasks = parse_task_list(ESSENTIAL_TASKS)
    print(f"Timing find_task_essential_reactions for {len(tasks)} tasks "
          f"(overall budget {DIAG_BUDGET}s)", flush=True)

    rows = []
    start = time.perf_counter()
    for i, task in enumerate(tasks, 1):
        if time.perf_counter() - start > DIAG_BUDGET:
            print(f"Budget reached; stopping before task {i}/{len(tasks)} ({task.id})", flush=True)
            break
        # Print before the call so a true hang still identifies the culprit task.
        print(f"[{i}/{len(tasks)}] starting {task.id} ...", flush=True)
        t0 = time.perf_counter()
        try:
            result = find_task_essential_reactions(model, [task])
            dt = time.perf_counter() - t0
            n = len(result.reactions)
            print(f"[{i}/{len(tasks)}] {task.id}: {dt:8.1f}s  ({n} essential rxns)", flush=True)
            rows.append((dt, task.id, n))
        except Exception as exc:  # noqa: BLE001
            dt = time.perf_counter() - t0
            print(f"[{i}/{len(tasks)}] {task.id}: FAILED after {dt:.1f}s: "
                  f"{type(exc).__name__}: {exc}", flush=True)
            rows.append((dt, task.id, -1))

    rows.sort(reverse=True)
    print("\n=== slowest tasks ===", flush=True)
    for dt, tid, n in rows[:12]:
        print(f"  {dt:8.1f}s  {tid}  ({'FAILED' if n < 0 else str(n) + ' rxns'})", flush=True)
    total = time.perf_counter() - start
    print(f"\nTimed {len(rows)}/{len(tasks)} tasks in {total:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
