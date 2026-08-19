"""Graded, task-scoped gene-essentiality analysis (see issue #1076).

This replaces the earlier all-task binary scoring, which called a gene essential when
its knockout made *any* of the 57 tasks in ``metabolicTasks_Essential.txt`` infeasible
and compared that single boolean with the Hart 2015 fitness genes. The task list mixes
two kinds of task:

  * **viability**: ``GR`` growth and ``ER`` energy/redox, what a proliferating cell
    must do;
  * **capability**: ``SU`` substrate utilization, ``BS`` biosynthesis of products and
    ``IC`` internal conversions, what the network should be *able* to do.

Hart 2015 measures proliferation fitness in rich medium, so a gene essential only for
a capability task is counted as a false positive. The clearest case is the ETF complex
(ETFA/ETFB/ETFDH), essential only for the ``SU`` beta-oxidation task and the four
``BS`` phospholipid tasks while its biomass growth ratio is 1.0.

For every cell line this module reports:

  * the task categories each gene is essential for, so capability-only essentiality is
    labelled instead of silently counted as a proliferation prediction;
  * the single-gene-deletion biomass growth ratio under the ``GR`` task's Ham's medium,
    a continuous score that needs no task threshold at all;
  * both binary scorings (all tasks, and viability tasks only) against Hart, plus the
    threshold-free AUROC/AUPRC of the growth ratio against Hart's Bayes Factors.

This is the gene-essentiality entry point of the ``/run gene-essentiality`` workflow,
replacing the earlier all-task binary scoring. It needs Gurobi and takes hours, so it is
not part of the per-pull-request checks. It overwrites the two usual artifacts in
``data/testResults/``: ``gene-essential.csv`` (per gene and cell line) and
``gene-essential_summary.md`` (the metric table posted to the pull request). Usage:

    python code/test/gradedEssentiality.py [--checkpoint-dir DIR]

``--checkpoint-dir`` stores one JSON file per finished cell line and skips the cell
lines already present, so an interrupted run resumes without repeating hours of work.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import cobra
from cobra.flux_analysis import single_gene_deletion

from raven_toolbox.io import read_yaml_model
from raven_toolbox.tasks.tasklist import Task, parse_task_list

from estimateEssentialGenes import (
    BIG_M,
    ESSENTIAL_TASKS,
    MIP_GAP_ABS,
    RNASEQ_FILE,
    TIME_LIMIT,
    # Module-private helpers of the model-building step, so the cell-line models are
    # built exactly as before; only the scoring of the knockouts changed.
    _build_context_model,
    _gene_symbol_map,
    _prep_human_model_for_ftinit,
    _read_rnaseq,
)
from evaluateHart2015Essentiality import evaluate_graded, summary_to_markdown
from taskEssentialGenes import find_task_essential_categories

# Repository root: this file is code/test/gradedEssentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"
RESULTS_DIR = REPO_ROOT / "data" / "testResults"
MATRIX_CSV = RESULTS_DIR / "gene-essential.csv"
SUMMARY_MD = RESULTS_DIR / "gene-essential_summary.md"

# Biomass reaction optimised for the growth ratio (the [GR] task's output).
BIOMASS_REACTION = "MAR13082"

def _log(message: str) -> None:
    print(message, flush=True)


def media_exchange_ids(model: cobra.Model, task: Task) -> tuple[set[str], list[str]]:
    """Boundary-reaction ids for the metabolites a task allows as input.

    Task inputs are written as ``name[compartment]`` (for example ``glucose[e]``), so
    they are resolved by metabolite name and compartment. Returns the matched exchange
    ids and the input tokens that could not be resolved.
    """
    exchange_ids: set[str] = set()
    unresolved: list[str] = []
    for token, _lb, _ub in task.inputs:
        name = token.rsplit("[", 1)[0].strip()
        compartment = token.rsplit("[", 1)[1].rstrip("]") if "[" in token else None
        matched = {
            rxn.id
            for met in model.metabolites
            if met.name and met.name.lower() == name.lower()
            and (compartment is None or met.compartment == compartment)
            for rxn in met.reactions if rxn.boundary
        }
        if matched:
            exchange_ids |= matched
        else:
            unresolved.append(token)
    return exchange_ids, unresolved


def growth_ratios(
    model: cobra.Model,
    growth_task: Task,
    *,
    biomass_id: str = BIOMASS_REACTION,
) -> tuple[dict[str, float], float, int, list[str]]:
    """Single-gene-deletion biomass growth ratio under the growth task's medium.

    The medium is opened exactly on the growth task's inputs and biomass is maximised
    (rather than forced to a fixed flux as the task itself does), so each knockout gets
    a continuous score: 1.0 means fully dispensable for growth, 0.0 growth-lethal.

    Returns ``(ratios, wild_type_flux, n_media_exchanges, unresolved_inputs)``. The
    ratios are empty when the wild type cannot grow on that medium.
    """
    m = model.copy()
    exchange_ids, unresolved = media_exchange_ids(m, growth_task)
    exchange_ids &= {rxn.id for rxn in m.exchanges}
    m.medium = {rxn_id: 1000.0 for rxn_id in exchange_ids}
    m.objective = biomass_id
    wild_type = m.slim_optimize()
    if wild_type is None or math.isnan(wild_type) or wild_type <= 1e-6:
        flux = 0.0 if wild_type is None or math.isnan(wild_type) else wild_type
        return {}, flux, len(exchange_ids), unresolved

    deletions = single_gene_deletion(m, processes=1)
    ratios: dict[str, float] = {}
    for _, row in deletions.iterrows():
        ids = row["ids"]
        gene_id = next(iter(ids)) if isinstance(ids, (frozenset, set)) else ids
        growth = row.get("growth")
        if growth is None or (isinstance(growth, float) and math.isnan(growth)):
            growth = 0.0
        ratios[gene_id] = growth / wild_type
    return ratios, wild_type, len(exchange_ids), unresolved


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Graded, task-scoped gene-essentiality analysis")
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=None,
        help="directory for per-cell-line JSON checkpoints; finished cell lines are skipped",
    )
    args = parser.parse_args(argv)
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    # Silence Gurobi's per-copy "Read LP format model from file ..." banner.
    import gurobipy
    gurobipy.setParam("OutputFlag", 0)
    # Serial FVA: Gurobi's environment is not fork-safe (see estimateEssentialGenes).
    cobra.Configuration().processes = 1

    started = time.time()
    model = read_yaml_model(MODEL_FILE)
    model.solver = "gurobi"
    tasks = parse_task_list(ESSENTIAL_TASKS)
    growth_task = next(task for task in tasks if task.id == "GR")
    tissues, expression = _read_rnaseq(RNASEQ_FILE)
    symbol_of = _gene_symbol_map(model)

    def checkpoint_path(tissue: str) -> Path | None:
        return args.checkpoint_dir / f"graded_{tissue}.json" if args.checkpoint_dir else None

    per_tissue: dict[str, dict[str, dict]] = {}
    pending = [t for t in tissues if not (checkpoint_path(t) and checkpoint_path(t).exists())]

    prep = None
    if pending:
        _log("Step 1: prepHumanModelForftINIT (clean + prepINITModel) ...")
        prep = _prep_human_model_for_ftinit(model, tasks)
        _log(f"Step 1 done: reference {len(prep.ref_model.reactions)} reactions, "
             f"{len(prep.essential_rxns)} task-essential; {len(prep.tasks)} feasible tasks")

    for index, tissue in enumerate(tissues, start=1):
        cached = checkpoint_path(tissue)
        if cached and cached.exists():
            _log(f"Cell line {index}/{len(tissues)}: {tissue} (from checkpoint)")
            per_tissue[tissue] = json.loads(cached.read_text())
            continue

        _log(f"Cell line {index}/{len(tissues)}: {tissue}")
        context = _build_context_model(
            prep, model, expression[tissue], BIG_M, MIP_GAP_ABS, TIME_LIMIT
        )
        context.id = tissue
        _log(f"  {tissue}: model has {len(context.reactions)} reactions, "
             f"{len(context.genes)} genes; scanning task categories ...")
        categories = find_task_essential_categories(
            context, prep.tasks, log=lambda message, t=tissue: _log(f"    {t}: {message}")
        )
        _log(f"  {tissue}: {len(categories)} genes essential for at least one task")

        ratios, wild_type, n_media, unresolved = growth_ratios(context, growth_task)
        if not ratios:
            _log(f"  {tissue}: WARNING no growth on the task medium (flux {wild_type:.3g}), "
                 f"growth ratios unavailable")
        if unresolved:
            _log(f"  {tissue}: WARNING unresolved medium inputs: {unresolved}")
        lethal = sum(1 for value in ratios.values() if value < 1e-6)
        _log(f"  {tissue}: wild-type biomass {wild_type:.3g} on {n_media} medium exchanges; "
             f"{lethal} growth-lethal genes")

        per_tissue[tissue] = {
            gene.id: {"tasks": sorted(categories.get(gene.id, ())), "growth": ratios.get(gene.id)}
            for gene in context.genes
        }
        if cached:
            cached.write_text(json.dumps(per_tissue[tissue]))
            _log(f"  {tissue}: checkpoint written to {cached}")

    rows, matrix_csv = evaluate_graded(
        per_tissue, tissues, sorted(symbol_of), symbol_of=symbol_of
    )
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    MATRIX_CSV.write_text(matrix_csv, encoding="utf-8")
    SUMMARY_MD.write_text(summary_to_markdown(rows), encoding="utf-8")
    _log(summary_to_markdown(rows))
    _log(f"Wrote {MATRIX_CSV} and {SUMMARY_MD} in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
