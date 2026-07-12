"""Generate cell-line-specific models and estimate essential genes.

Python port of code/test/estimateEssentialGenes.m. For each cell line in the
Hart 2015 RNA-seq dataset a context-specific model is built with classic tINIT
(raven_toolbox.init.get_init_model), matching the tINIT2 method the original
MATLAB workflow used, and the genes essential for the essential metabolic tasks
are identified in each model.

Why classic tINIT with the prep_init_model feasibility fixes: forcing the
task-essential reactions to carry flux makes the genome-scale INIT MILP
infeasible unless the template is first prepared the way prep_init_model does for
ftINIT. Two steps are needed and are applied here directly to the classic path:
orient each task-essential reaction irreversibly in its required direction, and
rescale the stoichiometry (rescale_for_init) so a single fixed big-M (100) is
valid across all reactions. After the MILP, redundant low-expression isozyme
genes are pruned with remove_low_score_genes, matching getINITModel2's
removeGenes=true step, so unexpressed alternatives in an OR rule do not mask the
genes that are actually essential for the expressed isozyme.

Gene identifiers
----------------
The MATLAB version relabelled the model genes with their symbols up front. Here
the Ensembl gene ids are kept throughout (they match the RNA-seq identifiers, so
expression scoring is exact) and are mapped to symbols only at the end, for the
comparison against the Hart 2015 table, which uses gene symbols.

Note: building genome-scale tINIT models solves a MILP per cell line and is
computationally heavy; Gurobi is required.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path

import cobra
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions

from raven_toolbox.init import (
    gene_scores_from_expression,
    get_init_model,
    remove_low_score_genes,
)
from raven_toolbox.init.prep import rescale_for_init
from raven_toolbox.tasks.check import find_task_essential_reactions
from raven_toolbox.tasks.tasklist import parse_task_list

from taskEssentialGenes import find_task_essential_genes

# Repository root: this file is code/test/estimateEssentialGenes.py
REPO_ROOT = Path(__file__).resolve().parents[2]
RNASEQ_FILE = REPO_ROOT / "data" / "datasets" / "Hart2015_RNAseq.txt"
ESSENTIAL_TASKS = REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"

# tINIT scoring threshold: genes expressed above 1 (TPM) score positive, below
# score negative (RAVEN's threshold = 1 in the MATLAB workflow).
EXPRESSION_THRESHOLD = 1.0

# tINIT MILP parameters. big_m = 100 is RAVEN's fixed value, valid once the
# stoichiometry is rescaled; mip_gap and time_limit bound each solve.
BIG_M = 100.0
MIP_GAP = 0.001
TIME_LIMIT = 1800.0


def _log(message: str) -> None:
    """Emit a progress line to stderr so it streams into the CI log."""
    print(message, file=sys.stderr, flush=True)


def _set_solver_verbosity(enabled: bool) -> None:
    """Toggle Gurobi's console output (the MILP branch-and-bound log).

    Enabled only around the tINIT solve, so its progress (incumbent, bound, gap
    over time) streams into the CI log, while the copy-heavy phases stay quiet and
    do not flood the log with "Read LP format model from file" banners.
    """
    try:
        import gurobipy
        gurobipy.setParam("OutputFlag", 1 if enabled else 0)
    except Exception:  # noqa: BLE001 - gurobipy optional; no-op without it
        pass


def _read_rnaseq(rnaseq_file: str | Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Read the RNA-seq table into per-cell-line expression maps.

    Returns ``(tissues, expression)`` where ``tissues`` is the column order and
    ``expression[tissue]`` maps gene id -> expression level.
    """
    table = pd.read_table(rnaseq_file)
    gene_col = table.columns[0]
    tissues = [c for c in table.columns[1:]]
    genes = table[gene_col].astype(str)
    expression = {
        tissue: dict(zip(genes, pd.to_numeric(table[tissue], errors="coerce").fillna(0.0)))
        for tissue in tissues
    }
    return tissues, expression


def _gene_symbol_map(model: cobra.Model) -> dict[str, str]:
    """Map gene id -> upper-case symbol (gene.name), falling back to the id."""
    return {g.id: (g.name or g.id).upper() for g in model.genes}


def _orient_forward(rxn: cobra.Reaction, direction: int) -> None:
    """Make ``rxn`` carry flux only in its forced direction (irreversible forward).

    Same primitive prep_init_model applies before the ftINIT MILP: a reaction
    forced in the reverse direction is flipped so its forced direction is forward,
    then it is made irreversible (lower bound >= 0).
    """
    if direction < 0:
        rxn.add_metabolites({m: -2 * c for m, c in rxn.metabolites.items()})
        rxn.bounds = (-rxn.upper_bound, -rxn.lower_bound)
    rxn.lower_bound = max(rxn.lower_bound, 0.0)


def estimate_essential_genes(
    model: cobra.Model,
    rnaseq_file: str | Path = RNASEQ_FILE,
    task_file: str | Path = ESSENTIAL_TASKS,
    *,
    big_m: float = BIG_M,
    mip_gap: float = MIP_GAP,
    time_limit: float | None = TIME_LIMIT,
) -> dict:
    """Build tINIT models per cell line and find their task-essential genes.

    Returns a dict with:
        tissues                  cell-line names, in RNA-seq column order
        gene_ids                 sorted list of all reference-model gene ids
        symbol_of                {gene id: upper-case gene symbol}
        essential_ids_by_tissue  {tissue: set of gene ids essential for any task}
        context_models           {tissue: cobra.Model} the tINIT models
    """
    # Run flux-variability analysis single-process. cobra's default parallel FVA
    # forks worker processes, but Gurobi's environment is not fork-safe, so the
    # workers deadlock non-deterministically on Linux CI runners (the step hangs
    # until the job timeout, leaving orphaned python processes). Serial FVA is a
    # little slower but reliable.
    cobra.Configuration().processes = 1

    tasks = parse_task_list(task_file)
    tissues, expression = _read_rnaseq(rnaseq_file)
    symbol_of = _gene_symbol_map(model)

    # Once-per-template preparation, expression-independent, so done a single time:
    #   1. discover the task-essential reactions and their forced directions;
    #   2. orient them irreversibly forward so the MILP forces them with a lower bound;
    #   3. rescale the stoichiometry so the fixed big-M is valid;
    #   4. drop blocked reactions once (keeping essentials) so get_init_model does not
    #      repeat that sweep per cell line.
    _log(f"Step 1a: finding task-essential reactions ({len(tasks)} tasks, one FVA each) ...")
    essential = find_task_essential_reactions(model, tasks).reactions  # {rxn id: +1/-1}
    essential_ids = set(essential)
    _log(f"Step 1a done: {len(essential_ids)} task-essential reactions")

    _log("Step 1b: orienting essential reactions and rescaling stoichiometry ...")
    prepared = model.copy()
    for rxn_id, direction in essential.items():
        _orient_forward(prepared.reactions.get_by_id(rxn_id), direction)
    rescale_for_init(prepared)

    _log(f"Step 1c: removing blocked reactions from {len(prepared.reactions)} reactions ...")
    blocked = set(find_blocked_reactions(prepared, open_exchanges=True))
    prepared.remove_reactions(sorted(blocked - essential_ids), remove_orphans=True)
    _log(f"Step 1 done: template ready with {len(prepared.reactions)} reactions "
         f"({len(essential_ids)} task-essential)")

    _log(f"Step 2: build {len(tissues)} cell-line-specific models with classic tINIT and "
         f"estimate essential genes")
    essential_ids_by_tissue: dict[str, set[str]] = {}
    context_models: dict[str, cobra.Model] = {}
    for i, tissue in enumerate(tissues, start=1):
        _log(f"  Cell line {i}/{len(tissues)}: {tissue}")
        gene_scores = gene_scores_from_expression(expression[tissue], EXPRESSION_THRESHOLD)

        _log(f"    {tissue}: running tINIT ...")
        context = _build_context_model(prepared, model, gene_scores, essential_ids,
                                       big_m, mip_gap, time_limit)
        context.id = tissue
        context_models[tissue] = context
        _log(f"    {tissue}: model has {len(context.reactions)} reactions, "
             f"{len(context.genes)} genes; finding essential genes ...")

        essential_gene_ids = find_task_essential_genes(
            context, tasks, log=lambda m, t=tissue: _log(f"    {t}: {m}")
        )
        essential_ids_by_tissue[tissue] = set(essential_gene_ids)
        _log(f"    {tissue}: {len(essential_gene_ids)} essential genes")

    return {
        "tissues": tissues,
        "gene_ids": sorted(symbol_of),
        "symbol_of": symbol_of,
        "essential_ids_by_tissue": essential_ids_by_tissue,
        "context_models": context_models,
    }


def _build_context_model(
    prepared: cobra.Model,
    reference: cobra.Model,
    gene_scores: Mapping[str, float],
    essential_ids: set[str],
    big_m: float,
    mip_gap: float,
    time_limit: float | None,
) -> cobra.Model:
    """Run classic tINIT on the prepared template and return the context model.

    The MILP runs on the oriented/rescaled ``prepared`` model; blocked reactions
    were already removed, so ``remove_dead_ends`` is off. The reactions it keeps
    are mapped back onto the unscaled ``reference`` model by id, so the returned
    context model has the original stoichiometry.
    """
    _set_solver_verbosity(True)  # stream the tINIT MILP branch-and-bound log as progress
    try:
        result = get_init_model(
            prepared, gene_scores=gene_scores, essential_rxns=essential_ids,
            remove_dead_ends=False, big_m=big_m, mip_gap=mip_gap, time_limit=time_limit,
        )
    finally:
        _set_solver_verbosity(False)  # quiet again for the copy-heavy gene-essentiality scan
    kept = {r.id for r in result.model.reactions}
    context = reference.copy()
    context.remove_reactions([r for r in context.reactions if r.id not in kept], remove_orphans=True)
    # Prune redundant low-expression isozyme genes, as MATLAB getINITModel2 does with
    # removeGenes=true: keep the highest-scoring isozyme in each OR rule so unexpressed
    # alternatives do not mask genes that are essential for the expressed isozyme.
    context, removed = remove_low_score_genes(context, gene_scores)
    _log(f"    pruned {len(removed)} redundant low-score genes")
    return context
