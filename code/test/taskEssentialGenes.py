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
"""

from __future__ import annotations

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


def find_task_essential_genes(
    model: cobra.Model,
    tasks: str | Iterable[Task],
    *,
    log: Callable[[str], None] | None = None,
) -> set[str]:
    """Return the set of gene ids essential for at least one task in ``model``.

    ``model`` is a context-specific model; ``tasks`` is a parsed task list or a
    path to a task-list file. Boundary reactions are closed internally so that
    task inputs/outputs define the exchange, exactly as in check_tasks.

    ``log`` is an optional callable used to report progress; when omitted the
    function is silent.
    """
    emit = log or (lambda _msg: None)
    tasks = _as_tasks(tasks)
    base = _prepare_base(model)
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
    total = len(gene_disabled)
    emit(f"{len(passing)}/{len(testable)} tasks feasible; scanning {total} candidate genes")

    essential_genes: set[str] = set()
    solves = 0
    for i, (gene_id, disabled) in enumerate(gene_disabled.items(), start=1):
        if not disabled:
            continue
        for task, flux_set in passing:
            # The knockout can only matter if it hits a reaction carrying flux in
            # this task's solution; otherwise that solution survives the knockout.
            if not (disabled & flux_set):
                continue
            solves += 1
            if not _task_feasible_without(base, task, name_to_id, comp_to_ids, disabled, saved):
                essential_genes.add(gene_id)
                break
        if i % 250 == 0 or i == total:
            emit(f"scanned {i}/{total} genes, {solves} solves, {len(essential_genes)} essential")

    return essential_genes
