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

Building every cell line is independent given the same Step-1 output (the once-per-run
:func:`_prep_human_model_for_ftinit` result: a cleaned, task-essential-annotated,
merged-and-scaled ``PrepData``, identical for every cell line since it is built before
any cell line's expression data is applied), so a run can be split across several
parallel processes or CI jobs by cell line:

    python code/test/gradedEssentiality.py --checkpoint-dir DIR \
        --shard-index I --shard-count N --prep-cache PREP_PATH

with ``I`` from ``0`` to ``N - 1``. Step 1 is expensive (tens of minutes) but produces
the same result regardless of which cell line asks for it, so it belongs to
``--prep-cache`` rather than to any one shard: build it once with ``--prep-only``,

    python code/test/gradedEssentiality.py --prep-only --prep-cache PREP_PATH

then point every shard at that same ``PREP_PATH`` (a pickled ``PrepData``, e.g. shared
between CI jobs as an uploaded/downloaded artifact) so each one loads it instead of
recomputing it. Without ``--prep-cache`` a shard falls back to computing its own Step 1,
as before. Each shard builds only the cell lines where ``index % N == I``, writes their
checkpoints, and exits without writing the final report. Once every shard has finished:

    python code/test/gradedEssentiality.py --checkpoint-dir DIR --aggregate-only

reads every cell line's checkpoint (Gurobi is not needed for this step) and writes the
final ``gene-essential.csv`` / ``gene-essential_summary.md``, exactly as a single
unsharded run would have.

Building a cell line's context-specific model (:func:`_build_context_model`, a genome-scale
MILP) and scanning it for task-essential genes each hold their own Gurobi session for the
whole step. A shared license's concurrent-session limit is best respected by never letting
more cell lines run one of these steps at once than the license allows, so ``--model-dir``
splits the two into a barrier: build every shard's context model and pickle it there first,

    python code/test/gradedEssentiality.py --build-models-only --model-dir MODEL_DIR \
        --shard-index I --shard-count N --prep-cache PREP_PATH

then, once every shard's model-building run has finished (so none of those Gurobi sessions
are still open), point the normal run at the same ``--model-dir`` to load each cell line's
model instead of rebuilding it:

    python code/test/gradedEssentiality.py --checkpoint-dir DIR --model-dir MODEL_DIR \
        --shard-index I --shard-count N --prep-cache PREP_PATH

Without ``--model-dir`` a shard builds its own context models inline, as before.

Within one cell line, the distinct gene knockouts are independent of each other (see
:func:`taskEssentialGenes.find_task_essential_categories`) and are split across
``--processes`` worker processes, one machine's cores rather than one CI shard; it
defaults to the machine's core count.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pickle
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


def _shard(tissues: list[str], shard_index: int, shard_count: int) -> list[str]:
    """This shard's cell lines: every ``shard_count``-th one, starting at ``shard_index``."""
    return tissues if shard_count <= 1 else tissues[shard_index::shard_count]


def _load_or_build_prep(model: cobra.Model, tasks, prep_cache: Path | None):
    """Step 1's ``PrepData``: loaded from ``prep_cache`` if present, else built and cached.

    Step 1 (:func:`_prep_human_model_for_ftinit`) is expensive but expression-independent,
    so its result is identical for every cell line; ``prep_cache`` lets several shards
    share one build of it instead of each running Step 1 for itself.
    """
    if prep_cache and prep_cache.exists():
        _log(f"Step 1: loading shared prepData from {prep_cache}")
        with open(prep_cache, "rb") as fh:
            return pickle.load(fh)

    _log("Step 1: prepHumanModelForftINIT (clean + prepINITModel) ...")
    prep = _prep_human_model_for_ftinit(model, tasks)
    _log(f"Step 1 done: reference {len(prep.ref_model.reactions)} reactions, "
         f"{len(prep.essential_rxns)} task-essential; {len(prep.tasks)} feasible tasks")

    if prep_cache:
        prep_cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = prep_cache.with_name(prep_cache.name + ".part")
        with open(tmp, "wb") as fh:
            pickle.dump(prep, fh)
        tmp.replace(prep_cache)
        _log(f"Step 1: wrote shared prepData to {prep_cache}")
    return prep


def _use_gurobi_solver(model: cobra.Model, *, retries: int = 5, initial_delay: float = 30.0) -> None:
    """Switch ``model`` onto the Gurobi solver, retrying past a transient session cap.

    Acquiring a Gurobi WLS session is a round trip against the license server's own
    session table, and a sibling shard's session does not always clear from that table
    the instant its process exits -- a finished build-test shard's session blocked the
    very next job (aggregate) from acquiring its own with "Too many sessions, 5 active
    sessions for a baseline of 2" moments after the 5-way build-test matrix finished. A
    concurrent shard's session clearing a minute or two into a wait is far more likely
    than the license actually being exhausted, so back off and retry a bounded number
    of times before giving up.

    A failed attempt can leave gurobipy's process-wide default environment in a
    half-initialised state (the crash is inside its own construction), so a retry that
    reuses it blind could behave differently from the first attempt for reasons that
    have nothing to do with the license server; best-effort dispose it before each
    retry so every attempt starts from the same clean slate.

    Logs how many attempts and how long acquisition took, and the CPU count / Gurobi
    thread count once acquired -- separating "time spent getting a session" from "time
    spent solving" and surfacing possible CPU oversubscription across concurrent
    shards, both invisible from the outside while a job is still running.
    """
    import gurobipy

    started = time.time()
    delay = initial_delay
    for attempt in range(1, retries + 1):
        try:
            model.solver = "gurobi"
            break
        except gurobipy.GurobiError as exc:
            if attempt == retries:
                raise
            _log(f"Gurobi session unavailable ({exc}); retrying in {delay:.0f}s "
                 f"({attempt}/{retries}) ...")
            try:
                gurobipy.disposeDefaultEnv()
            except Exception:
                pass  # nothing to clean up, or the env never got that far -- fine either way
            time.sleep(delay)
            delay *= 2

    elapsed = time.time() - started
    threads = model.solver.problem.Params.Threads
    _log(f"Gurobi session acquired in {elapsed:.1f}s ({attempt} attempt(s)); "
         f"{os.cpu_count()} CPU(s) visible, solver Threads={threads}")


def _dispose_gurobi_session(*, timeout: float = 20.0) -> None:
    """Best-effort: close this process's default Gurobi environment right away.

    A cobra Model's Gurobi-backed solver otherwise only closes its WLS session when
    the interpreter exits, and even then the license server does not always notice the
    dropped connection immediately -- it can keep counting the session as active for a
    while, which is exactly what let a finished shard's session block the next job's
    session request (see :func:`_use_gurobi_solver`). Call this once every Gurobi-backed
    object this process created has gone out of scope (the caller should first drop its
    own references, e.g. reassign the built context model and the template model to
    ``None``) so the session is released now instead of on the license server's own
    timeout.

    ``disposeDefaultEnv()`` releasing a WLS session is a network round trip to Gurobi's
    license server, unlike the everyday case of a local/named-user license closing
    instantly -- an already-finished, already-written-to-disk shard must never be stuck
    "in progress" waiting on that call, so it runs on a background thread with a bounded
    join instead of inline. A timeout (or any other error) is logged and swallowed: the
    process is about to exit either way, which releases the session itself, just later.
    """
    import gc
    import threading

    gc.collect()
    done = threading.Event()
    error: list[BaseException] = []

    def _dispose() -> None:
        try:
            import gurobipy
            gurobipy.disposeDefaultEnv()
        except BaseException as exc:  # noqa: BLE001 -- reported on the main thread below
            error.append(exc)
        finally:
            done.set()

    thread = threading.Thread(target=_dispose, daemon=True)
    thread.start()
    if not done.wait(timeout):
        _log(f"WARNING: Gurobi session release did not finish within {timeout:.0f}s; "
             f"moving on (the process exiting will release it anyway).")
    elif error:
        _log(f"WARNING: could not release the Gurobi session cleanly: {error[0]}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Graded, task-scoped gene-essentiality analysis")
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=None,
        help="directory for per-cell-line JSON checkpoints; finished cell lines are skipped",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        default=0,
        help="0-based index of this run's slice of cell lines; requires --checkpoint-dir",
    )
    parser.add_argument(
        "--shard-count",
        type=int,
        default=1,
        help="split cell lines into this many slices; this run builds index I where "
             "I %% shard-count == shard-index, then exits without the final report",
    )
    parser.add_argument(
        "--aggregate-only",
        action="store_true",
        help="skip building; combine every cell line's --checkpoint-dir checkpoint "
             "(all must already be present) into the final report",
    )
    parser.add_argument(
        "--prep-cache",
        type=Path,
        default=None,
        help="pickle path for the shared, expression-independent Step-1 PrepData: "
             "loaded from here if present, else computed and written here",
    )
    parser.add_argument(
        "--prep-only",
        action="store_true",
        help="build --prep-cache (required) and exit without building any cell line",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="directory holding one pickled context-specific model per cell line: written "
             "there by --build-models-only, and loaded from there instead of rebuilding if "
             "present otherwise",
    )
    parser.add_argument(
        "--build-models-only",
        action="store_true",
        help="build and pickle this shard's cell-line context models into --model-dir "
             "(required), then exit without gene-essentiality scanning",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=os.cpu_count() or 1,
        help="worker processes for one cell line's gene-knockout scan "
             "(default: this machine's core count)",
    )
    args = parser.parse_args(argv)
    if args.aggregate_only and not args.checkpoint_dir:
        parser.error("--aggregate-only requires --checkpoint-dir")
    if args.shard_count > 1 and not args.build_models_only and not args.checkpoint_dir:
        parser.error("--shard-count > 1 requires --checkpoint-dir (unless --build-models-only)")
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must be in [0, --shard-count) and --shard-count >= 1")
    if args.prep_only and not args.prep_cache:
        parser.error("--prep-only requires --prep-cache")
    if args.prep_only and args.aggregate_only:
        parser.error("--prep-only and --aggregate-only are mutually exclusive")
    if args.build_models_only and not args.model_dir:
        parser.error("--build-models-only requires --model-dir")
    if args.build_models_only and (args.prep_only or args.aggregate_only):
        parser.error("--build-models-only and --prep-only/--aggregate-only are mutually exclusive")
    if args.checkpoint_dir:
        args.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    tasks = parse_task_list(ESSENTIAL_TASKS)
    growth_task = next(task for task in tasks if task.id == "GR")
    tissues, expression = _read_rnaseq(RNASEQ_FILE)
    # cobra defaults to Gurobi automatically whenever gurobipy is importable (its
    # preferred-solver order tries "gurobi" before "glpk"), so building the model below
    # would open a Gurobi WLS session even for --aggregate-only, which never solves
    # anything and does not need one -- exactly what let the aggregate job's model load
    # collide with a still-active build-test session. Force glpk up front; the branches
    # that actually need Gurobi's MILP solver switch to it explicitly below.
    cobra.Configuration().solver = "glpk"
    model = read_yaml_model(MODEL_FILE)
    symbol_of = _gene_symbol_map(model)

    def checkpoint_path(tissue: str) -> Path | None:
        return args.checkpoint_dir / f"graded_{tissue}.json" if args.checkpoint_dir else None

    def model_path(tissue: str) -> Path | None:
        return args.model_dir / f"context_{tissue}.pkl" if args.model_dir else None

    if args.prep_only:
        _use_gurobi_solver(model)
        import gurobipy
        gurobipy.setParam("OutputFlag", 0)
        cobra.Configuration().processes = 1
        _load_or_build_prep(model, tasks, args.prep_cache)
        _log(f"Wrote {args.prep_cache} in {time.time() - started:.0f}s")
        model = None
        _dispose_gurobi_session()
        return 0

    if args.build_models_only:
        _use_gurobi_solver(model)
        import gurobipy
        gurobipy.setParam("OutputFlag", 0)
        cobra.Configuration().processes = 1

        shard = _shard(tissues, args.shard_index, args.shard_count)
        prep = _load_or_build_prep(model, tasks, args.prep_cache)
        args.model_dir.mkdir(parents=True, exist_ok=True)

        context = None
        for index, tissue in enumerate(shard, start=1):
            target = model_path(tissue)
            if target.exists():
                _log(f"Cell line {index}/{len(shard)}: {tissue} model (from cache)")
                continue

            _log(f"Cell line {index}/{len(shard)}: {tissue}")
            context = _build_context_model(
                prep, model, expression[tissue], BIG_M, MIP_GAP_ABS, TIME_LIMIT
            )
            context.id = tissue
            _log(f"  {tissue}: model has {len(context.reactions)} reactions, "
                 f"{len(context.genes)} genes")
            tmp = target.with_name(target.name + ".part")
            with open(tmp, "wb") as fh:
                pickle.dump(context, fh)
            tmp.replace(target)
            _log(f"  {tissue}: wrote context model to {target}")

        _log(f"Shard {args.shard_index}/{args.shard_count} model-build done "
             f"in {time.time() - started:.0f}s")
        context = prep = model = None
        _dispose_gurobi_session()
        return 0

    if args.aggregate_only:
        missing = [t for t in tissues if not (checkpoint_path(t) and checkpoint_path(t).exists())]
        if missing:
            _log(f"ERROR: no checkpoint for {len(missing)}/{len(tissues)} cell line(s): {missing}")
            return 1
        per_tissue = {t: json.loads(checkpoint_path(t).read_text()) for t in tissues}
    else:
        _use_gurobi_solver(model)
        # Silence Gurobi's per-copy "Read LP format model from file ..." banner.
        import gurobipy
        gurobipy.setParam("OutputFlag", 0)
        # Serial FVA: Gurobi's environment is not fork-safe (see estimateEssentialGenes).
        cobra.Configuration().processes = 1

        shard = _shard(tissues, args.shard_index, args.shard_count)
        pending = [t for t in shard if not (checkpoint_path(t) and checkpoint_path(t).exists())]

        per_tissue = {}
        prep = None
        if pending:
            prep = _load_or_build_prep(model, tasks, args.prep_cache)

        for index, tissue in enumerate(shard, start=1):
            cached = checkpoint_path(tissue)
            if cached and cached.exists():
                _log(f"Cell line {index}/{len(shard)}: {tissue} (from checkpoint)")
                per_tissue[tissue] = json.loads(cached.read_text())
                continue

            _log(f"Cell line {index}/{len(shard)}: {tissue}")
            cached_model = model_path(tissue)
            if cached_model and cached_model.exists():
                _log(f"  {tissue}: loading context model from {cached_model}")
                with open(cached_model, "rb") as fh:
                    context = pickle.load(fh)
            else:
                context = _build_context_model(
                    prep, model, expression[tissue], BIG_M, MIP_GAP_ABS, TIME_LIMIT
                )
                context.id = tissue
            _log(f"  {tissue}: model has {len(context.reactions)} reactions, "
                 f"{len(context.genes)} genes; scanning task categories ...")
            categories = find_task_essential_categories(
                context, prep.tasks, processes=args.processes,
                log=lambda message, t=tissue: _log(f"    {t}: {message}"),
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

        context = prep = model = None
        _dispose_gurobi_session()
        if args.shard_count > 1:
            _log(f"Shard {args.shard_index}/{args.shard_count} done in {time.time() - started:.0f}s; "
                 f"run --aggregate-only once every shard has finished.")
            return 0

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
