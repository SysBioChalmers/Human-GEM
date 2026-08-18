"""Graded, task-scoped gene-essentiality analysis (see issue #1076).

The standard run (``geneEssentiality.py``) calls a gene essential when its knockout
makes *any* of the 57 tasks in ``metabolicTasks_Essential.txt`` infeasible, and
compares that binary call with the Hart 2015 fitness genes. That task list mixes two
kinds of task:

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

This is an analysis entry point, not part of the per-pull-request checks: it needs
Gurobi and takes several hours. Two artifacts are written to ``data/testResults/``:
``gene-essential_graded.csv`` (per gene and cell line) and
``gene-essential_graded_summary.md`` (the metric table). Usage:

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
from collections.abc import Iterable, Mapping
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
    # Module-private helpers of the standard pipeline, reused here so the graded
    # analysis builds exactly the same context models as geneEssentiality.py.
    _build_context_model,
    _gene_symbol_map,
    _prep_human_model_for_ftinit,
    _read_rnaseq,
)
from evaluateHart2015Essentiality import BF_THRESHOLDS, bayes_factors
from taskEssentialGenes import find_task_essential_categories

# Repository root: this file is code/test/gradedEssentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"
RESULTS_DIR = REPO_ROOT / "data" / "testResults"
GRADED_CSV = RESULTS_DIR / "gene-essential_graded.csv"
GRADED_SUMMARY = RESULTS_DIR / "gene-essential_graded_summary.md"

# Task categories representing viability of a proliferating cell, and those
# representing a capability the reconstruction should have. Together they are the five
# category ids used by metabolicTasks_Essential.txt.
VIABILITY_CATEGORIES = frozenset({"GR", "ER"})
CAPABILITY_CATEGORIES = frozenset({"SU", "BS", "IC"})

# Biomass reaction optimised for the growth ratio (the [GR] task's output).
BIOMASS_REACTION = "MAR13082"

SUMMARY_COLUMNS = [
    "cellLine", "allTaskMCC", "allTaskFP", "viabilityMCC", "viabilityFP",
    "growthAUROC", "growthAUPRC", "baseRate", "capabilityOnly",
]


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


def auroc(scored: Iterable[tuple[float, int]]) -> float:
    """Area under the ROC curve for ``(score, label)`` pairs, ties rank-averaged."""
    pairs = list(scored)
    n_pos = sum(label for _score, label in pairs)
    n_neg = len(pairs) - n_pos
    if not n_pos or not n_neg:
        return math.nan
    order = sorted(pairs, key=lambda pair: pair[0])
    ranks = [0.0] * len(order)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and order[j + 1][0] == order[i][0]:
            j += 1
        average = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = average
        i = j + 1
    rank_sum = sum(rank for rank, (_score, label) in zip(ranks, order) if label)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def auprc(scored: Iterable[tuple[float, int]]) -> float:
    """Average precision (area under the precision-recall curve).

    Preferred over AUROC as the headline number here because essential genes are the
    minority class, so precision-recall reflects the useful operating range. Compare it
    with the base rate (the share of positives), which is the random-guess level.
    """
    pairs = sorted(scored, key=lambda pair: -pair[0])
    n_pos = sum(label for _score, label in pairs)
    if not n_pos:
        return math.nan
    true_pos = false_pos = 0
    average_precision = 0.0
    previous_recall = 0.0
    for _score, label in pairs:
        if label:
            true_pos += 1
        else:
            false_pos += 1
        recall = true_pos / n_pos
        if recall > previous_recall:
            average_precision += (recall - previous_recall) * (true_pos / (true_pos + false_pos))
            previous_recall = recall
    return average_precision


def _mcc(counts: Mapping[str, int]) -> float:
    tp, tn, fp, fn = counts["TP"], counts["TN"], counts["FP"], counts["FN"]
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn) - (fp * fn)) / denominator if denominator else math.nan


def _confusion_class(predicted: bool, fitness: bool) -> str:
    return ("TP" if fitness else "FP") if predicted else ("FN" if fitness else "TN")


def evaluate_graded(
    per_tissue: Mapping[str, Mapping[str, dict]],
    tissues: list[str],
    gene_ids: Iterable[str],
    *,
    symbol_of: Mapping[str, str],
) -> tuple[list[dict], str]:
    """Score the graded predictions against Hart 2015.

    ``per_tissue`` maps cell line -> gene id -> ``{"tasks": [category ids], "growth":
    ratio or None}``. Returns the per-cell-line summary rows and the per-gene CSV text.
    """
    gene_ids = list(gene_ids)
    hart_bf = bayes_factors(gene_ids, tissues, symbol_of=symbol_of)

    rows: list[dict] = []
    pooled: list[tuple[float, int]] = []
    for tissue in tissues:
        predictions = per_tissue.get(tissue, {})
        by_gene = hart_bf.get(tissue, {})
        threshold = BF_THRESHOLDS[tissue]
        all_counts = dict(TP=0, TN=0, FP=0, FN=0)
        viability_counts = dict(TP=0, TN=0, FP=0, FN=0)
        graded: list[tuple[float, int]] = []
        capability_only = 0
        for gene_id, record in predictions.items():
            categories = set(record["tasks"])
            essential_any = bool(categories)
            essential_viability = bool(categories & VIABILITY_CATEGORIES)
            if essential_any and not essential_viability:
                capability_only += 1
            if gene_id not in by_gene:
                continue
            fitness = by_gene[gene_id] > threshold
            all_counts[_confusion_class(essential_any, fitness)] += 1
            viability_counts[_confusion_class(essential_viability, fitness)] += 1
            growth = record.get("growth")
            if growth is not None:
                graded.append((1.0 - growth, int(fitness)))
        pooled += graded
        base_rate = (sum(label for _s, label in graded) / len(graded)) if graded else math.nan
        rows.append({
            "cellLine": tissue,
            "allTaskMCC": round(_mcc(all_counts), 4),
            "allTaskFP": all_counts["FP"],
            "viabilityMCC": round(_mcc(viability_counts), 4),
            "viabilityFP": viability_counts["FP"],
            "growthAUROC": round(auroc(graded), 4),
            "growthAUPRC": round(auprc(graded), 4),
            "baseRate": round(base_rate, 4),
            "capabilityOnly": capability_only,
        })
    if pooled:
        rows.append({
            "cellLine": "all",
            "allTaskMCC": "", "allTaskFP": "", "viabilityMCC": "", "viabilityFP": "",
            "growthAUROC": round(auroc(pooled), 4),
            "growthAUPRC": round(auprc(pooled), 4),
            "baseRate": round(sum(label for _s, label in pooled) / len(pooled), 4),
            "capabilityOnly": "",
        })

    header = ["genes", "geneSymbol"]
    for tissue in tissues:
        header += [f"{tissue}_class", f"{tissue}_tasks", f"{tissue}_growth"]
    lines = [",".join(header)]
    for gene_id in gene_ids:
        row = [gene_id, symbol_of.get(gene_id, "")]
        for tissue in tissues:
            record = per_tissue.get(tissue, {}).get(gene_id)
            if record is None:
                row += [".", ".", "."]
                continue
            categories = set(record["tasks"])
            essential_any = bool(categories)
            essential_viability = bool(categories & VIABILITY_CATEGORIES)
            by_gene = hart_bf.get(tissue, {})
            if gene_id in by_gene:
                fitness = by_gene[gene_id] > BF_THRESHOLDS[tissue]
                old = _confusion_class(essential_any, fitness)
                new = _confusion_class(essential_viability, fitness)
            else:
                old = "P" if essential_any else "."
                new = "P" if essential_viability else "."
            # Mark the cells whose call the viability scoping actually changes, so a
            # meaningful reclassification (for example ETF moving FP -> TN) stands out.
            growth = record.get("growth")
            row += [
                new if old == new else f"{old}->{new}",
                "|".join(sorted(categories)),
                "" if growth is None else f"{growth:.4f}",
            ]
        lines.append(",".join(row))
    return rows, "\n".join(lines) + "\n"


def summary_to_markdown(rows: list[dict]) -> str:
    """Render the per-cell-line summary as a Markdown table."""
    out = [
        "### Graded gene essentiality vs Hart 2015 (task-scope analysis)",
        "",
        "| " + " | ".join(SUMMARY_COLUMNS) + " |",
        "| " + " | ".join("---" for _ in SUMMARY_COLUMNS) + " |",
    ]
    for row in rows:
        out.append("| " + " | ".join(str(row.get(column, "")) for column in SUMMARY_COLUMNS) + " |")
    return "\n".join(out) + "\n"


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
    GRADED_CSV.write_text(matrix_csv, encoding="utf-8")
    GRADED_SUMMARY.write_text(summary_to_markdown(rows), encoding="utf-8")
    _log(summary_to_markdown(rows))
    _log(f"Wrote {GRADED_CSV} and {GRADED_SUMMARY} in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
