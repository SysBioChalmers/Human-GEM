"""Test metabolic tasks against Human-GEM using the Python RAVEN toolbox.

Python port of code/test/testMetabolicTasks.m. Reads the YAML model, runs the
requested metabolic task list and fails (non-zero exit) if any task does not
pass. Boundary handling that the MATLAB version obtained via addBoundaryMets is
provided here by check_tasks(..., close_boundaries=True), which closes existing
exchange/sink/demand reactions so that inputs and outputs are defined purely by
the tasks (as RAVEN assumes).

Usage:
    python code/test/testMetabolicTasks.py {essential|verification}
"""

import sys
from pathlib import Path

from raven_toolbox.io import read_yaml_model
from raven_toolbox.tasks import check_tasks

# Repository root: this file is code/test/testMetabolicTasks.py
REPO_ROOT = Path(__file__).resolve().parents[2]

TASK_FILES = {
    "essential": REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt",
    "verification": REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_VerifyModel.txt",
}


def main(task_type: str) -> int:
    if task_type not in TASK_FILES:
        print(f"::error::Unknown task type '{task_type}'. Use 'essential' or 'verification'.")
        return 1

    model = read_yaml_model(REPO_ROOT / "model" / "Human-GEM.yml")
    # Use GLPK (bundled with cobra) for the task-feasibility LPs. This keeps the
    # test free of any commercial-solver licence: if gurobipy happens to be
    # installed, optlang would otherwise auto-select the size-limited Gurobi
    # licence, which rejects a genome-scale model.
    model.solver = "glpk"
    results = check_tasks(model, TASK_FILES[task_type], close_boundaries=True)

    failed = [r for r in results if not r.passed]
    if failed:
        for r in failed:
            reason = r.error if r.error else "task not satisfied"
            print(f"::error::Failed task [{r.id}] {r.description}: {reason}")
        print(f"::error::Failed in {task_type} tasks ({len(failed)}/{len(results)} failed).")
        return 1

    print(f"Succeeded with {task_type} tasks ({len(results)} passed).")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("::error::Task type is missing!")
        sys.exit(1)
    sys.exit(main(sys.argv[1]))
