"""Generate cell-line-specific models with ftINIT and estimate essential genes.

Python port of the DevelopWBM ``Gen_ftINIT_models`` Human2 pipeline
(github.com/LiLabTsinghua/DevelopWBM). For each cell line in the Hart 2015 RNA-seq
dataset a context-specific model is built with **ftINIT** (not classic tINIT /
getINITModel2), matching how DevelopWBM produces the published gene-essentiality
numbers, and the genes essential for the essential metabolic tasks are identified
in each model.

The MATLAB flow reproduced here is::

    prepData = prepHumanModelForftINIT(model, false, tasks, reactions.tsv);
    models{i} = ftINIT(prepData, tissue, [], [], data_struct, {}, ...
                       getHumanGEMINITSteps('1+1'), true, true);   % removeGenes, useScoresForTasks
    m = addBoundaryMets(m);
    [~, essentialGenes] = checkTasksGenes(m, [], false, false, true, taskStruct);

Mapping to raven-toolbox:

  * ``prepHumanModelForftINIT`` -> remove drug / amino-acid-triplet reactions and
    MAR13081, load spontaneous reactions from ``model/reactions.tsv``, then
    :func:`raven_toolbox.init.prep_init_model` (== RAVEN ``prepINITModel``) once.
  * ``ftINIT(..., '1+1', removeGenes=true, useScoresForTasks=true)`` ->
    :func:`raven_toolbox.init.ftinit` with ``series='1+1'``, ``gene_scores`` supplied
    (prunes negative-scoring genes, == ``removeGenes``) and ``fill_gaps=True`` with
    score-weighted task gap-filling (== ``useScoresForTasks``).
  * ``checkTasksGenes(..., getEssential=true)`` -> ``find_task_essential_genes``.

Gene identifiers
----------------
The MATLAB version relabelled the model genes with their symbols up front. Here
the Ensembl gene ids are kept throughout (they match the RNA-seq identifiers, so
expression scoring is exact) and are mapped to symbols only at the end, for the
comparison against the Hart 2015 table, which uses gene symbols.

Note: building genome-scale ftINIT models solves a MILP per cell line and is
computationally heavy; Gurobi is required.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path

import cobra
import pandas as pd

from raven_toolbox.init import (
    ftinit,
    gene_scores_from_expression,
    prep_init_model,
    score_reactions_from_genes,
)
from raven_toolbox.tasks.tasklist import parse_task_list

from taskEssentialGenes import find_task_essential_genes

# ftINIT gap-filling copies genome-scale cobra models via recursive deepcopy, which can
# exceed Python's default 1000-frame recursion limit on a model this size.
sys.setrecursionlimit(10000)

# Repository root: this file is code/test/estimateEssentialGenes.py
REPO_ROOT = Path(__file__).resolve().parents[2]
RNASEQ_FILE = REPO_ROOT / "data" / "datasets" / "Hart2015_RNAseq.txt"
ESSENTIAL_TASKS = REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"
REACTIONS_TSV = REPO_ROOT / "model" / "reactions.tsv"

# Extracellular compartment abbreviation in Human-GEM v2.0.0 (RAVEN 's' in Human1).
EXT_COMP = "e"

# tINIT scoring threshold: genes expressed above 1 (TPM) score positive, below
# score negative (RAVEN's threshold = 1 in the MATLAB workflow).
EXPRESSION_THRESHOLD = 1.0

# ftINIT MILP parameters. big_m = 100 is RAVEN's fixed value, valid once the merged
# stoichiometry is rescaled (prep_init_model does this); each per-step solve is bounded
# by an absolute MIP gap and a time limit for tractability on the genome-scale model.
BIG_M = 100.0
MIP_GAP_ABS = 10.0
TIME_LIMIT = 1800.0

# prepHumanModelForftINIT: reactions that can "always be on" and are ignored while
# scoring (protein creation/degradation + metabolite-pooling reactions). The commented
# radical reactions in the MATLAB source are intentionally left out (not spontaneous).
CUSTOM_RXNS_TO_IGNORE = (
    # protein reactions creation/degradation
    "MAR05155", "MAR05156", "MAR05161", "MAR05167", "MAR05168", "MAR05169", "MAR05170",
    "MAR05171", "MAR05172", "MAR05174", "MAR05260", "MAR05262", "MAR05264", "MAR05266",
    "MAR05267", "MAR05268", "MAR05269", "MAR05270", "MAR05271", "MAR05273", "MAR05275",
    "MAR05277", "MAR05279", "MAR05281", "MAR05283", "MAR05291", "MAR09817", "MAR09818",
    # reactions that just pool metabolites
    "MAR00011", "MAR00012", "MAR00477", "MAR05233", "MAR05234", "MAR05238", "MAR05239",
    "MAR05243", "MAR05244", "MAR05247", "MAR09022", "MAR00015", "MAR00016", "MAR00017",
    "MAR10033", "MAR10035", "MAR10036", "MAR10037", "MAR10038", "MAR10062", "MAR10063",
    "MAR10064", "MAR10065", "MAR13082",
)

# The duplicate complex-I reaction with ROS: removed outright (prepHumanModelForftINIT).
DUPLICATE_COMPLEX_I_RXN = "MAR13081"

# removeDrugReactions: metabolites whose name contains one of these drugs (case
# insensitive); every reaction touching such a metabolite is removed.
_DRUG_NAMES = (
    "pravastatin", "gliclazide", "atorvastatin", "fluvastatin", "fluvastain",
    "fluvstatin", "simvastatin", "cyclosporine", "acetaminophen", "cerivastatin",
    "tacrolimus", "ibuprofen", "lovastatin", "losartan", "nifedipine", "pitavastatin",
    "rosuvastatin", "torasemide", "midazolam",
)

# getAATripletReactions: a metabolite is an amino-acid triplet/doublet if its name
# splits on '-' into amino-acyl residues followed by a terminal amino acid (or is the
# single peptide 'Argtyrval'). Reactions touching such metabolites are removed. Matching
# is exact and case-sensitive, as in the MATLAB (strcmp on the '-'-split name parts).
_AMINO_ACIDS_CAP = frozenset({
    "Alanine", "Arginine", "Asparagine", "Aspartate", "Cysteine", "Glutamine",
    "Glutamate", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine",
    "Methionine", "Metheonine", "Phenylalanine", "Proline", "Serine", "Threonine",
    "Tryptophan", "Tyrosine", "Valine",
})
_AMINO_ACYLS = frozenset({
    "Alanyl", "Alaninyl", "Alanine", "Arginyl", "Asparaginyl", "Aspartyl", "Cystyl",
    "Cystinyl", "Cysteinyl", "Glutaminyl", "Glutamyl", "Glutamatsyl", "Glycyl",
    "Histidyl", "Histidinyl", "Isoleucyl", "Isolecyl", "Leucyl", "Lysyl", "Lysine",
    "Methionyl", "Methioninyl", "Phenylalanyl", "Phenylalanine", "Phenylalaninyl",
    "Prolyl", "Seryl", "Threonyl", "Tryptophanyl", "Tyrosyl", "Tyrosinyl", "Valyl",
})
_SINGLE_TRIPLET_NAMES = frozenset({"Argtyrval"})


def _log(message: str) -> None:
    """Emit a progress line to stderr so it streams into the CI log."""
    print(message, file=sys.stderr, flush=True)


def _set_solver_verbosity(enabled: bool) -> None:
    """Toggle Gurobi's console output (the MILP branch-and-bound log).

    Enabled only around the ftINIT solves, so their progress (incumbent, bound, gap
    over time) streams into the CI log, while the copy-heavy phases stay quiet and do
    not flood the log with "Read LP format model from file" banners.
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


def _drug_reactions(model: cobra.Model) -> list[str]:
    """Reactions touching a drug metabolite (RAVEN ``removeDrugReactions``)."""
    drug_mets = [
        m for m in model.metabolites
        if any(d in (m.name or "").lower() for d in _DRUG_NAMES)
    ]
    return sorted({r.id for m in drug_mets for r in m.reactions})


def _is_aa_triplet_metabolite(name: str) -> bool:
    """Whether a metabolite name is an amino-acid triplet/doublet peptide."""
    parts = name.split("-")
    if len(parts) == 1:
        return parts[0] in _SINGLE_TRIPLET_NAMES
    return all(p in _AMINO_ACYLS for p in parts[:-1]) and parts[-1] in _AMINO_ACIDS_CAP


def _aa_triplet_reactions(model: cobra.Model) -> list[str]:
    """Reactions touching an amino-acid-triplet metabolite (``getAATripletReactions``).

    Ports the ``onlyRxnsWithoutGPRs = false`` call used by prepHumanModelForftINIT:
    every reaction touching such a metabolite is returned, regardless of its GPR.
    """
    aa_mets = [m for m in model.metabolites if _is_aa_triplet_metabolite(m.name or "")]
    return sorted({r.id for m in aa_mets for r in m.reactions})


def _load_spontaneous_reactions(tsv_path: str | Path) -> list[str]:
    """Reaction ids flagged spontaneous in ``reactions.tsv`` (spontaneous == 1).

    ``importTsvFile`` returns the column as text when the file has no quoted fields
    (Human-GEM v2.0.0+), so the value is coerced to numeric (see Human-GEM #1020).
    """
    table = pd.read_table(tsv_path, dtype=str)
    spont = pd.to_numeric(table["spontaneous"], errors="coerce").fillna(0)
    return table.loc[spont == 1, "rxns"].astype(str).tolist()


def _prep_human_model_for_ftinit(
    model: cobra.Model, tasks, *, essential_cache_path: str | Path | None = None
):
    """Port of ``prepHumanModelForftINIT``: clean the model and build ftINIT prepData.

    Removes drug reactions, amino-acid-triplet reactions and the duplicate complex-I
    reaction, loads the spontaneous reactions, then runs the once-per-template
    :func:`prep_init_model` (== RAVEN ``prepINITModel``). ``prep_init_model`` discovers
    the task-essential reactions, orients and merges, and rescales the stoichiometry.
    """
    m = model.copy()
    to_remove = set(_drug_reactions(m))
    _log(f"  removeDrugReactions: {len(to_remove)} reactions")
    aa = set(_aa_triplet_reactions(m))
    _log(f"  getAATripletReactions: {len(aa)} reactions")
    to_remove |= aa
    to_remove |= {DUPLICATE_COMPLEX_I_RXN} & {r.id for r in m.reactions}
    m.remove_reactions(sorted(to_remove), remove_orphans=True)
    _log(f"  cleaned template: {len(m.reactions)} reactions, {len(m.genes)} genes")

    spontaneous = [r for r in _load_spontaneous_reactions(REACTIONS_TSV)
                   if r in {rx.id for rx in m.reactions}]
    custom = [r for r in CUSTOM_RXNS_TO_IGNORE if r in {rx.id for rx in m.reactions}]
    _log(f"  spontaneous={len(spontaneous)} custom-ignore={len(custom)}; running prepINITModel ...")

    return prep_init_model(
        m, tasks, ext_comp=EXT_COMP, spontaneous=spontaneous, custom=custom,
        essential_cache_path=essential_cache_path,
    )


def estimate_essential_genes(
    model: cobra.Model,
    rnaseq_file: str | Path = RNASEQ_FILE,
    task_file: str | Path = ESSENTIAL_TASKS,
    *,
    big_m: float = BIG_M,
    mip_gap_abs: float = MIP_GAP_ABS,
    time_limit: float | None = TIME_LIMIT,
    essential_cache_path: str | Path | None = None,
) -> dict:
    """Build ftINIT models per cell line and find their task-essential genes.

    Returns a dict with:
        tissues                  cell-line names, in RNA-seq column order
        gene_ids                 sorted list of all reference-model gene ids
        symbol_of                {gene id: upper-case gene symbol}
        essential_ids_by_tissue  {tissue: set of gene ids essential for any task}
        context_models           {tissue: cobra.Model} the ftINIT models
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

    # Once-per-template preparation (expression-independent): clean the model and run
    # prepINITModel, which discovers the task-essential reactions, orients/merges, and
    # rescales the stoichiometry. This is the ~30-minute step and is shared by every
    # cell line (DevelopWBM builds prepData once, then calls ftINIT per tissue).
    _log("Step 1: prepHumanModelForftINIT (clean + prepINITModel) ...")
    prep = _prep_human_model_for_ftinit(model, tasks, essential_cache_path=essential_cache_path)
    _log(f"Step 1 done: reference {len(prep.ref_model.reactions)} reactions, "
         f"merged {len(prep.min_model.reactions)}, {len(prep.essential_rxns)} task-essential; "
         f"{len(prep.tasks)} feasible tasks")

    _log(f"Step 2: build {len(tissues)} cell-line-specific models with ftINIT '1+1' "
         f"and estimate essential genes")
    essential_ids_by_tissue: dict[str, set[str]] = {}
    context_models: dict[str, cobra.Model] = {}
    for i, tissue in enumerate(tissues, start=1):
        _log(f"  Cell line {i}/{len(tissues)}: {tissue}")
        context = _build_context_model(
            prep, model, expression[tissue], big_m, mip_gap_abs, time_limit
        )
        context.id = tissue
        context_models[tissue] = context
        _log(f"    {tissue}: model has {len(context.reactions)} reactions, "
             f"{len(context.genes)} genes; finding essential genes ...")

        essential_gene_ids = find_task_essential_genes(
            context, prep.tasks, log=lambda m, t=tissue: _log(f"    {t}: {m}")
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
    prep,
    reference: cobra.Model,
    expression: Mapping[str, float],
    big_m: float,
    mip_gap_abs: float,
    time_limit: float | None,
) -> cobra.Model:
    """Run ftINIT ('1+1') on the shared prepData for one cell line.

    Scores the reference reactions from this cell line's expression, then runs the
    multi-step ftINIT MILP. ``gene_scores`` is passed so negative-scoring genes are
    pruned from the GPRs (RAVEN ``removeGenes=true``); ``fill_gaps=True`` restores the
    minimum score-weighted set of reactions needed for every task (``useScoresForTasks``).
    """
    gene_scores = gene_scores_from_expression(expression, EXPRESSION_THRESHOLD)
    rxn_scores = score_reactions_from_genes(prep.ref_model, gene_scores)

    _set_solver_verbosity(True)  # stream the ftINIT MILP branch-and-bound log as progress
    try:
        context = ftinit(
            prep, rxn_scores, gene_scores=gene_scores, series="1+1", fill_gaps=True,
            big_m=big_m, mip_gap_abs=mip_gap_abs, time_limit=time_limit,
        )
    finally:
        _set_solver_verbosity(False)  # quiet again for the copy-heavy gene-essentiality scan
    return context
