"""Test metabolic tasks against Human-GEM using the Python RAVEN toolbox.

Python port of code/test/testMetabolicTasks.m. Reads the YAML model, runs the
requested metabolic task list and fails (non-zero exit) if any task does not
pass. Boundary handling that the MATLAB version obtained via addBoundaryMets is
provided here by check_tasks(..., close_boundaries=True), which closes existing
exchange/sink/demand reactions so that inputs and outputs are defined purely by
the tasks (as RAVEN assumes).

Before the tasks run, every task file in data/metabolicTasks/ is checked against the
rules RAVEN's checkTasks applies before it solves anything, which would otherwise abort
a RAVEN run on the whole file:

  * every IN/OUT metabolite (``name[compartment]``) exists in the model; ``ALLMETS``
    and ``ALLMETSIN[compartment]`` are accepted;
  * no metabolite is both an input and an output of one task with a nonzero lower
    bound on both sides (checkTasks: "The IN LB and OUT LB ... cannot be nonzero for
    the same metabolite").

The problems go to the log, and pass or fail to qc_status.tsv (``task_files``).

The ``essential`` run also runs the ``GR`` task with the generic cell components
``MAR00021`` as the biomass reaction (``MAR13082`` blocked), as the animal-GEM
generator does when it gap-fills (code/animalGEM/generateAnimalGEM.py); its result
goes to qc_status.tsv as ``tasks_gr_mar00021``.

Usage:
    python code/test/testMetabolicTasks.py {essential|verification|all}
"""

import sys
from pathlib import Path

from raven_toolbox.io import read_yaml_model
from raven_toolbox.tasks import check_tasks
from raven_toolbox.tasks.tasklist import parse_task_list

import qcStatus

# Repository root: this file is code/test/testMetabolicTasks.py
REPO_ROOT = Path(__file__).resolve().parents[2]

TASK_DIR = REPO_ROOT / "data" / "metabolicTasks"

# The human biomass reaction and the generic cell components the animal-GEM generator
# makes the objective instead.
HUMAN_BIOMASS_RXN = "MAR13082"
BIOMASS_COMPONENTS_RXN = "MAR00021"

TASK_FILES = {
    "essential": REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt",
    "verification": REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_VerifyModel.txt",
}


def _split_token(token: str) -> tuple[str, str | None]:
    """``name[compartment]`` -> (name, compartment); no brackets -> (token, None)."""
    if token.endswith("]") and "[" in token:
        name, compartment = token[:-1].rsplit("[", 1)
        return name.strip(), compartment.strip()
    return token.strip(), None


def task_file_problems(model, path: Path) -> list[str]:
    """Problems in one task file that make RAVEN's checkTasks stop (see the module docstring)."""
    # RAVEN upper-cases names and compartments before matching, so both are compared without case.
    known = {(m.name.upper(), m.compartment.lower()) for m in model.metabolites if m.name}
    compartments = {c.lower() for c in model.compartments}
    problems = []
    for task in parse_task_list(path):
        label = f"{path.name} [{task.id}] {task.description}"
        for side, items in (("IN", task.inputs), ("OUT", task.outputs)):
            for token, _lb, _ub in items:
                name, compartment = _split_token(token)
                if name.upper() == "ALLMETS" and compartment is None:
                    continue
                if name.upper() == "ALLMETSIN" and compartment and compartment.lower() in compartments:
                    continue
                if compartment is None or (name.upper(), compartment.lower()) not in known:
                    problems.append(f"{label}: {side} metabolite {token!r} is not in the model")
        nonzero_in = {token.upper() for token, lb, _ub in task.inputs if lb != 0}
        nonzero_out = {token.upper() for token, lb, _ub in task.outputs if lb != 0}
        for token in sorted(nonzero_in & nonzero_out):
            problems.append(f"{label}: {token} is an input and an output, both with a nonzero lower bound")
    return problems


def check_task_files(model, paths) -> int:
    problems = [p for path in paths for p in task_file_problems(model, path)]
    for problem in problems:
        print(f"::error::{problem}")
    print(f"Task files: {len(problems)} problem(s) in {len(paths)} file(s).")
    qcStatus.set_status("task_files", "fail" if problems else "pass")
    return 1 if problems else 0


def check_growth_on_components(model) -> int:
    """The ``GR`` task with ``MAR00021`` instead of ``MAR13082`` as the biomass reaction."""
    for rid in (HUMAN_BIOMASS_RXN, BIOMASS_COMPONENTS_RXN):
        if rid not in model.reactions:
            print(f"::error::{rid} is not in the model.")
            qcStatus.set_status("tasks_gr_mar00021", "fail")
            return 1
    growth = [t for t in parse_task_list(TASK_FILES["essential"]) if t.id == "GR"]
    test = model.copy()
    test.reactions.get_by_id(HUMAN_BIOMASS_RXN).bounds = (0, 0)
    test.reactions.get_by_id(BIOMASS_COMPONENTS_RXN).upper_bound = 1000
    results = check_tasks(test, growth, close_boundaries=True)
    passed = bool(results) and all(r.passed for r in results)
    if passed:
        print(f"GR task with {BIOMASS_COMPONENTS_RXN} as the biomass reaction: passed.")
    else:
        print(f"::error::The GR task fails with {BIOMASS_COMPONENTS_RXN} as the biomass reaction "
              f"({HUMAN_BIOMASS_RXN} blocked), which the animal-GEM generator relies on.")
    qcStatus.set_status("tasks_gr_mar00021", "pass" if passed else "fail")
    return 0 if passed else 1


def _check_one(model, task_type: str) -> int:
    results = check_tasks(model, TASK_FILES[task_type], close_boundaries=True)
    failed = [r for r in results if not r.passed]
    if failed:
        for r in failed:
            reason = r.error if r.error else "task not satisfied"
            print(f"::error::Failed task [{r.id}] {r.description}: {reason}")
        print(f"::error::Failed in {task_type} tasks ({len(failed)}/{len(results)} failed).")
    else:
        print(f"Succeeded with {task_type} tasks ({len(results)} passed).")
    # one-line status for the QC comment, in the shared qc_status.tsv
    qcStatus.set_status(f"tasks_{task_type}", f"{len(failed)}/{len(results)}")
    return 1 if failed else 0


def main(task_type: str) -> int:
    # "all" (the default in CI) runs both task lists in one job; a single type is
    # still accepted for ad-hoc runs.
    types = list(TASK_FILES) if task_type == "all" else [task_type]
    if any(t not in TASK_FILES for t in types):
        print(f"::error::Unknown task type '{task_type}'. Use 'essential', 'verification' or 'all'.")
        return 1

    model = read_yaml_model(REPO_ROOT / "model" / "Human-GEM.yml")
    # Use GLPK (bundled with cobra) for the task-feasibility LPs. This keeps the
    # test free of any commercial-solver licence: if gurobipy happens to be
    # installed, optlang would otherwise auto-select the size-limited Gurobi
    # licence, which rejects a genome-scale model.
    model.solver = "glpk"
    status = check_task_files(model, sorted(TASK_DIR.glob("metabolicTasks_*.txt")))
    status = max(status, *(_check_one(model, t) for t in types))
    if "essential" in types:
        status = max(status, check_growth_on_components(model))
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "all"))
