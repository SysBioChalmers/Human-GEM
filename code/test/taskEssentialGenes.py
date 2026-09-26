"""Identify genes essential for metabolic tasks in a model.

Python port of the RAVEN checkTasksGenes / getTaskEssentialGenes logic used by
the gene-essentiality workflow. A gene is essential for a task when knocking it
out turns an otherwise-feasible task infeasible.

Performance
-----------
The model is copied exactly once. Every task and gene-knockout test is then run
on that single base model inside cobra's ``with model:`` context manager, which
reverts reaction-bound, temporary-reaction and objective changes on exit. The
only changes the context manager does not track are the direct edits to the
metabolite mass-balance constraint bounds that RAVEN's apply_task_constraints
makes, so those are restored explicitly. This matters because copying a
Gurobi-backed cobra model serialises it to a temporary LP file and reads it back
(the "Read LP format model from file ..." lines); doing that per task or per
gene would dominate the runtime.

A gene knockout can only affect a task if at least one reaction it disables
actually carries flux in a feasible solution of that task. For each task we
compute one parsimonious (pFBA) flux distribution up front and only re-test a
gene against the (few) tasks whose flux-carrying reactions it disables.
Reactions essential for a task appear in every feasible solution, so this filter
never drops a real essential gene.

Two genes that knock out the identical set of reactions (subunits of the same
AND-linked complex, e.g. a mitochondrial ETC complex) get the identical
feasibility verdict for every task, so the scan tests each distinct knockout
once and copies the result to every gene that shares it, rather than solving it
once per gene. Distinct knockouts are also independent of each other, so
``processes`` splits them across a process pool the same way
:func:`raven_toolbox.tasks.check.find_task_essential_reactions` splits tasks:
each worker receives one copy of the base model at pool startup, not per
knockout.
"""

from __future__ import annotations

import concurrent.futures as cf
import multiprocessing
from collections import defaultdict
from collections.abc import Callable, Iterable

import cobra
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import pfba

from raven_toolbox.tasks.check import apply_task_constraints, task_name_maps
from raven_toolbox.tasks.tasklist import Task, parse_task_list

_TOL = 1e-8


def _as_tasks(tasks: str | Iterable[Task]) -> list[Task]:
    if isinstance(tasks, (str, bytes)) or hasattr(tasks, "__fspath__"):
        return parse_task_list(tasks)
    return list(tasks)


def _prepare_base(model: cobra.Model) -> cobra.Model:
    """Copy the model once and close its boundary reactions (as check_tasks does)."""
    base = model.copy()
    for rxn in base.boundary:
        rxn.bounds = (0.0, 0.0)
    return base


def _pin_single_threaded(model: cobra.Model) -> None:
    """Force one thread per LP solve (Gurobi-specific, a no-op on other backends).

    The scan re-solves thousands of small LPs that differ from each other by only a
    few reaction bounds. Gurobi's default automatic/concurrent method spins up several
    threads for a problem this size regardless, which adds synchronisation overhead to
    every one of those solves without shortening any single one of them. Pinning to one
    thread removes that overhead and leaves the rest of the machine's cores free for
    ``processes`` to use across knockouts instead of within one.
    """
    try:
        model.solver.problem.Params.Threads = 1
    except Exception:  # noqa: BLE001 - gurobipy optional; no-op without it
        pass


def _set_constraint_bounds(constraint, lb: float, ub: float) -> None:
    """Set an optlang constraint's bounds without a transient lb > ub."""
    if lb > constraint.ub:
        constraint.ub = ub
        constraint.lb = lb
    else:
        constraint.lb = lb
        constraint.ub = ub


def _gene_disabled_reactions(model: cobra.Model) -> dict[str, set[str]]:
    """Map gene id -> ids of reactions disabled when only that gene is knocked out.

    A reaction is disabled by a single-gene knockout when its GPR evaluates to
    False with that gene removed (accounting for isozymes / complexes).
    """
    mapping: dict[str, set[str]] = defaultdict(set)
    for rxn in model.reactions:
        gene_ids = [g.id for g in rxn.genes]
        if not gene_ids:
            continue
        gpr = rxn.gpr
        for gid in gene_ids:
            if not gpr.eval([gid]):
                mapping[gid].add(rxn.id)
    return mapping


def _restore_constraints(base: cobra.Model, met_ids: Iterable[str], saved: dict) -> None:
    """Restore the mass-balance bounds of ``met_ids`` from the ``saved`` snapshot."""
    for mid in met_ids:
        if mid in base.constraints:
            _set_constraint_bounds(base.constraints[mid], *saved[mid])


def _task_flux_set(base, task, name_to_id, comp_to_ids, original_ids, saved) -> set[str] | None:
    """Reactions of ``base`` carrying flux in a pFBA solution of ``task``.

    Runs on ``base`` in place (no copy) and restores it. Returns ``None`` if the
    task is malformed or infeasible on the intact model (it does not pass, so it
    defines no essential genes).
    """
    with base:
        task_mets, error = apply_task_constraints(base, task, name_to_id, comp_to_ids)
        try:
            if error is not None:
                return None
            try:
                fluxes = pfba(base).fluxes
            except OptimizationError:
                return None
            return {rid for rid in original_ids if abs(fluxes.get(rid, 0.0)) > _TOL}
        finally:
            _restore_constraints(base, task_mets or saved, saved)


def _task_feasible_without(
    base, task, name_to_id, comp_to_ids, knockout_rxns, saved
) -> bool:
    """Is ``task`` still feasible on ``base`` with ``knockout_rxns`` forced to zero?

    Runs on ``base`` in place (no copy) and restores it. The knockout is applied
    after the task constraints so it always wins over any reaction bound the task
    itself changes.
    """
    with base:
        task_mets, error = apply_task_constraints(base, task, name_to_id, comp_to_ids)
        try:
            if error is not None:
                return False
            for rid in knockout_rxns:
                if rid in base.reactions:
                    base.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            base.slim_optimize()
            return base.solver.status == "optimal"
        finally:
            _restore_constraints(base, task_mets or saved, saved)


def _prepare_scan(model: cobra.Model, tasks: str | Iterable[Task], emit):
    """Shared setup for the gene scans.

    Copies the model once, builds the task name maps, snapshots the mass-balance
    bounds and computes one parsimonious flux distribution per feasible task.
    Returns ``(base, passing, gene_disabled, name_to_id, comp_to_ids, saved)``.
    """
    tasks = _as_tasks(tasks)
    base = _prepare_base(model)
    _pin_single_threaded(base)
    name_to_id, comp_to_ids = task_name_maps(base)
    original_ids = {r.id for r in base.reactions}
    # Snapshot every mass-balance bound once so a task application can be reverted
    # even when apply_task_constraints errors after a partial modification.
    saved = {m.id: (base.constraints[m.id].lb, base.constraints[m.id].ub) for m in base.metabolites}

    # One parsimonious flux distribution per feasible, non-should-fail task. Each task
    # is paired with its own flux set (by position, not id): the task list reuses a few
    # ids across many distinct tasks (57 tasks under 5 ids: ER/BS/SU/IC/GR), so keying
    # flux sets by task.id let same-id tasks overwrite each other, applying the wrong
    # task's flux filter and dropping genes essential for the overwritten tasks.
    testable = [t for t in tasks if not t.should_fail]
    emit(f"computing flux distributions for {len(testable)} tasks")
    passing: list[tuple[Task, set[str]]] = []
    for task in testable:
        flux_set = _task_flux_set(base, task, name_to_id, comp_to_ids, original_ids, saved)
        if flux_set is not None:
            passing.append((task, flux_set))

    gene_disabled = _gene_disabled_reactions(base)
    emit(f"{len(passing)}/{len(testable)} tasks feasible; "
         f"scanning {len(gene_disabled)} candidate genes")
    return base, passing, gene_disabled, name_to_id, comp_to_ids, saved


def _test_knockout(
    base, passing, name_to_id, comp_to_ids, saved, disabled, all_categories, *, stop_early
) -> tuple[set[str], int]:
    """Task ids broken by forcing ``disabled`` to zero, and how many solves that took.

    With ``stop_early`` this stops at the first broken task, which is all that is
    needed to call the knockout essential and is much cheaper. Without it every
    matching task is tested, so the full set of broken task ids (= categories) is
    known; it also stops once every category has been broken once, since testing
    further tasks cannot add a new category.
    """
    broken_categories: set[str] = set()
    solves = 0
    for task, flux_set in passing:
        # The knockout can only matter if it hits a reaction carrying flux in this
        # task's solution; otherwise that solution survives the knockout.
        if not (disabled & flux_set):
            continue
        solves += 1
        if not _task_feasible_without(base, task, name_to_id, comp_to_ids, disabled, saved):
            broken_categories.add(task.id)
            if stop_early or broken_categories >= all_categories:
                break
    return broken_categories, solves


# Set once per worker process by _init_worker; module-level so ProcessPoolExecutor's
# workers (which import this module fresh rather than inheriting parent state on
# 'spawn') can reach it without re-pickling the shared arguments per knockout.
_WORKER: dict = {}


def _init_worker(base, passing, name_to_id, comp_to_ids, saved, all_categories, stop_early) -> None:
    _pin_single_threaded(base)
    _WORKER.update(
        base=base, passing=passing, name_to_id=name_to_id, comp_to_ids=comp_to_ids,
        saved=saved, all_categories=all_categories, stop_early=stop_early,
    )


def _knockout_worker(disabled: frozenset[str]) -> tuple[frozenset[str], set[str], int]:
    categories, solves = _test_knockout(
        _WORKER["base"], _WORKER["passing"], _WORKER["name_to_id"], _WORKER["comp_to_ids"],
        _WORKER["saved"], disabled, _WORKER["all_categories"], stop_early=_WORKER["stop_early"],
    )
    return disabled, categories, solves


def _scan(base, passing, gene_disabled, name_to_id, comp_to_ids, saved, emit, *, stop_early, processes=1):
    """Knock out each gene and record the id of every task the knockout breaks.

    Genes that disable the identical set of reactions (subunits of the same AND-linked
    complex) get the identical verdict for every task, so they are grouped and each
    distinct reaction set is tested once, then applied to every gene in its group. With
    ``processes`` > 1 the distinct knockouts are independent of each other and are split
    across a :class:`~concurrent.futures.ProcessPoolExecutor`, each worker holding its
    own copy of ``base`` (sent once at pool startup, not per knockout).
    """
    total = len(gene_disabled)
    all_categories = {task.id for task, _flux_set in passing}

    groups: dict[frozenset[str], list[str]] = defaultdict(list)
    for gene_id, disabled in gene_disabled.items():
        if disabled:
            groups[frozenset(disabled)].append(gene_id)
    n_shared = total - len(groups)
    if n_shared:
        emit(f"{len(groups)} distinct knockouts among {total} candidate genes "
             f"({n_shared} share a reaction set with another gene, e.g. multi-subunit complexes)")

    broken: dict[str, set[str]] = {}
    solves = 0
    scanned = 0

    def record(disabled_key: frozenset[str], categories: set[str]) -> None:
        nonlocal scanned
        for gene_id in groups[disabled_key]:
            scanned += 1
            if categories:
                broken[gene_id] = set(categories)
            if scanned % 250 == 0 or scanned == total:
                emit(f"scanned {scanned}/{total} genes, {solves} solves, {len(broken)} essential")

    if processes <= 1 or len(groups) <= 1:
        for disabled_key in groups:
            categories, n_solves = _test_knockout(
                base, passing, name_to_id, comp_to_ids, saved, disabled_key, all_categories,
                stop_early=stop_early,
            )
            solves += n_solves
            record(disabled_key, categories)
    else:
        # 'spawn', not the platform default ('fork' on Linux): forking duplicates the
        # parent's already-running Gurobi environment (and its internal license/logging
        # threads) into every worker, which is exactly the non-fork-safe state that made
        # cobra's default parallel FVA deadlock on Linux CI runners (the reason
        # estimateEssentialGenes.estimate_essential_genes and gradedEssentiality.main
        # both pin FVA to processes=1). 'spawn' starts each worker as a fresh interpreter
        # that only creates its own Gurobi environment inside _init_worker, after the
        # fork boundary, so there is nothing shared to deadlock on.
        ctx = multiprocessing.get_context("spawn")
        with cf.ProcessPoolExecutor(
            max_workers=min(processes, len(groups)),
            mp_context=ctx,
            initializer=_init_worker,
            initargs=(base, passing, name_to_id, comp_to_ids, saved, all_categories, stop_early),
        ) as pool:
            for disabled_key, categories, n_solves in pool.map(_knockout_worker, groups):
                solves += n_solves
                record(disabled_key, categories)
    return broken


def find_task_essential_genes(
    model: cobra.Model,
    tasks: str | Iterable[Task],
    *,
    log: Callable[[str], None] | None = None,
    processes: int = 1,
) -> set[str]:
    """Return the set of gene ids essential for at least one task in ``model``.

    ``model`` is a context-specific model; ``tasks`` is a parsed task list or a
    path to a task-list file. Boundary reactions are closed internally so that
    task inputs/outputs define the exchange, exactly as in check_tasks.

    ``log`` is an optional callable used to report progress; when omitted the
    function is silent. ``processes`` (default 1, sequential) splits the distinct
    gene knockouts across a process pool; see :func:`_scan`.
    """
    emit = log or (lambda _msg: None)
    scan_args = _prepare_scan(model, tasks, emit)
    return set(_scan(*scan_args, emit, stop_early=True, processes=processes))


def find_task_essential_categories(
    model: cobra.Model,
    tasks: str | Iterable[Task],
    *,
    log: Callable[[str], None] | None = None,
    processes: int = 1,
) -> dict[str, set[str]]:
    """Map gene id -> ids of the tasks it is essential for, for every essential gene.

    The essential-task list labels each task with its category, so the returned ids
    are the categories ``ER`` (energy and redox), ``IC`` (internal conversions),
    ``SU`` (substrate utilization), ``BS`` (biosynthesis of products) and ``GR``
    (growth). This distinguishes a gene needed for *viability* (``GR``/``ER``) from
    one needed only for a *capability* the network should have (``SU``/``BS``/``IC``),
    which :func:`find_task_essential_genes` cannot express.

    Unlike :func:`find_task_essential_genes` this cannot stop at the first broken
    task, so it is several times slower; it is meant for analysis, not for the
    per-pull-request run. ``processes`` (default 1, sequential) splits the distinct
    gene knockouts across a process pool; see :func:`_scan`.
    """
    emit = log or (lambda _msg: None)
    scan_args = _prepare_scan(model, tasks, emit)
    return _scan(*scan_args, emit, stop_early=False, processes=processes)
