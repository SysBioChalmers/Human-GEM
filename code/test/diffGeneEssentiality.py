"""Scope gene-essentiality comparison to the reactions a model diff actually touches.

The global metrics in gene-essential_summary.md (allTaskMCC, growthAUROC, ...) are
computed over every gene Hart 2015 scored, so a pull request that edits a handful of
reactions is swamped by ftINIT's run-to-run non-determinism in the other ~3000 genes.
This script instead:

  1. Loads the model at two git refs (base and head) with raven-toolbox and finds the
     reactions whose id, gene_reaction_rule or stoichiometry differ (added, removed or
     modified). The union of their GPR genes on either side is the "direct" gene set --
     this also catches a gene silently dropped from a shared reaction's rule, such as
     the flavoenzyme genes that used to ride along the ETF:ubiquinone-oxidoreductase
     OR-rule before #1028 (see issue #1076 / the etf-ubiquinone coupling fix).
  2. Expands one hop to "neighbor" genes: genes of other reactions (unchanged
     themselves) that share a metabolite with a changed reaction, e.g. everything else
     drawing on the ubiquinone/ubiquinol pool.
  3. Diffs the two already-committed gene-essential.csv matrices (no ftINIT/Gurobi
     re-run) for just this gene set. A gene is reported only when its predicted
     essentiality moved in the same direction in a majority of the five cell lines --
     a flip reproduced across most cell lines is far more likely a real network effect
     than a single cell line's ftINIT noise.
  4. Splits those consistent changes into two kinds, because they are not evaluated the
     same way:
       - a change to whether the knockout blocks *growth* (the GR/ER "viability" task
         categories) is compared against the Hart 2015 experimental fitness screen, so
         each one is labelled an improvement (moved the prediction closer to the
         experimental call) or a regression (moved it further away);
       - a change to whether the knockout blocks some other *capability* (SU/BS/IC --
         the network being able to do something, not the cell being able to grow) is
         reported but not scored against Hart, which only measures proliferation. This
         is expected to fire for genes intentionally coupled into a pathway (e.g. the
         ETF complex becoming essential for beta-oxidation) and is not by itself a
         concern.

Needs raven-toolbox to read the two model versions and the Hart 2015 table already in
data/datasets/; does not need Gurobi or a fresh ftINIT run.

Usage:
    python code/test/diffGeneEssentiality.py \\
        --base-model /path/to/base/Human-GEM.yml \\
        --base-matrix /path/to/base/gene-essential.csv \\
        [--head-model model/Human-GEM.yml] [--head-matrix data/testResults/gene-essential.csv] \\
        [--growth-tolerance 0.01] [--out-csv PATH] [--out-summary PATH]

The base model/matrix are typically pulled from the target branch with ``git show``,
for example:

    git show develop:model/Human-GEM.yml > base/Human-GEM.yml
    git show develop:data/testResults/gene-essential.csv > base/gene-essential.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from raven_toolbox.io import read_yaml_model

from evaluateHart2015Essentiality import BF_THRESHOLDS, bayes_factors

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HEAD_MODEL = REPO_ROOT / "model" / "Human-GEM.yml"
DEFAULT_HEAD_MATRIX = REPO_ROOT / "data" / "testResults" / "gene-essential.csv"

# GR/ER are what a proliferating cell must do and are the only categories comparable
# to Hart 2015 (a proliferation screen); SU/BS/IC are network capabilities the
# reconstruction should have, which Hart does not measure. See gradedEssentiality.py.
VIABILITY_CATEGORIES = frozenset({"GR", "ER"})
CATEGORY_LABELS = {
    "GR": "growth (biomass production)",
    "ER": "energy/redox balance",
    "SU": "substrate utilization",
    "BS": "biosynthesis",
    "IC": "internal conversion",
}

# A gene needs to flip in the same direction in at least this many compared cell
# lines to be reported as a consistent change rather than an isolated one.
CONSENSUS_MIN_LINES = 3


def _describe_categories(categories) -> str:
    return ", ".join(CATEGORY_LABELS.get(c, c) for c in sorted(categories)) if categories else "none"


def changed_reactions(base_model, head_model) -> dict[str, dict]:
    """Reactions added, removed, or modified (GPR or stoichiometry) between the two models.

    Returns ``{reaction_id: {"kind": "added"|"removed"|"modified", "genes": {gene_id, ...}}}``.
    ``genes`` is the union of the reaction's GPR genes on whichever side(s) it exists,
    so a gene dropped from (or added to) a surviving reaction's rule is included.
    """
    base_ids = {r.id for r in base_model.reactions}
    head_ids = {r.id for r in head_model.reactions}
    changed: dict[str, dict] = {}

    for rid in sorted(head_ids - base_ids):
        rxn = head_model.reactions.get_by_id(rid)
        changed[rid] = {"kind": "added", "genes": {g.id for g in rxn.genes}}

    for rid in sorted(base_ids - head_ids):
        rxn = base_model.reactions.get_by_id(rid)
        changed[rid] = {"kind": "removed", "genes": {g.id for g in rxn.genes}}

    for rid in sorted(base_ids & head_ids):
        b_rxn = base_model.reactions.get_by_id(rid)
        h_rxn = head_model.reactions.get_by_id(rid)
        b_mets = {m.id: coeff for m, coeff in b_rxn.metabolites.items()}
        h_mets = {m.id: coeff for m, coeff in h_rxn.metabolites.items()}
        if b_rxn.gene_reaction_rule != h_rxn.gene_reaction_rule or b_mets != h_mets:
            changed[rid] = {
                "kind": "modified",
                "genes": {g.id for g in b_rxn.genes} | {g.id for g in h_rxn.genes},
            }

    return changed


def neighbor_genes(head_model, changed: dict[str, dict], direct_genes: set[str]) -> set[str]:
    """One-hop genes: other reactions sharing a metabolite with a changed reaction.

    Only reactions still present in ``head_model`` seed the expansion (a removed
    reaction's old metabolite neighborhood is not followed), and genes already in
    ``direct_genes`` are excluded.
    """
    seed_metabolites: set[str] = set()
    for rid, info in changed.items():
        if info["kind"] == "removed" or not head_model.reactions.has_id(rid):
            continue
        rxn = head_model.reactions.get_by_id(rid)
        seed_metabolites |= {m.id for m in rxn.metabolites}

    changed_rxn_ids = set(changed)
    neighbors: set[str] = set()
    for met_id in seed_metabolites:
        met = head_model.metabolites.get_by_id(met_id)
        for rxn in met.reactions:
            if rxn.id in changed_rxn_ids:
                continue
            neighbors |= {g.id for g in rxn.genes}
    return neighbors - direct_genes


def read_matrix(path: Path) -> tuple[list[str], dict[str, str], dict[str, dict[str, dict | None]]]:
    """Parse a gene-essential.csv into ``(tissues, symbol_of, per_gene)``.

    ``per_gene[gene_id][tissue]`` is ``{"tasks": frozenset[str], "growth": float | None}``,
    or ``None`` when the gene was not part of that cell line's context-specific model.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        tissues = [c[: -len("_class")] for c in columns if c.endswith("_class")]
        symbol_of: dict[str, str] = {}
        per_gene: dict[str, dict[str, dict | None]] = {}
        for row in reader:
            gene_id = row["genes"]
            symbol_of[gene_id] = row.get("geneSymbol", "") or ""
            record: dict[str, dict | None] = {}
            for tissue in tissues:
                tasks_raw = (row.get(f"{tissue}_tasks") or "").strip()
                if tasks_raw == ".":
                    record[tissue] = None
                    continue
                growth_raw = (row.get(f"{tissue}_growth") or "").strip()
                record[tissue] = {
                    "tasks": frozenset(tasks_raw.split("|")) if tasks_raw else frozenset(),
                    "growth": None if growth_raw in ("", ".") else float(growth_raw),
                }
            per_gene[gene_id] = record
    return tissues, symbol_of, per_gene


def _direction(before: bool, after: bool) -> str | None:
    if before == after:
        return None
    return "gained" if after else "lost"


def _hart_call(direction: str, fitness: bool | None) -> str | None:
    """Whether a viability flip moved the prediction toward or away from Hart 2015.

    ``fitness`` is Hart's experimental call (True = fitness/essential gene in this cell
    line), or ``None`` when Hart did not score this gene in this cell line. A flip to
    "essential" that Hart backs up (or a flip to "non-essential" that Hart backs up
    too) is an improvement; the opposite direction is a regression.
    """
    if fitness is None:
        return None
    predicted_essential_after = direction == "gained"
    return "improvement" if predicted_essential_after == fitness else "regression"


def diff_gene(
    gene_id: str,
    tissues: list[str],
    base_record: dict[str, dict | None],
    head_record: dict[str, dict | None],
    hart_fitness: dict[str, bool],
    *,
    growth_tolerance: float,
) -> dict:
    """Per-cell-line comparison for one gene, plus consensus verdicts across lines."""
    per_line: dict[str, dict] = {}
    any_directions: list[str] = []
    viability_directions: list[str] = []
    hart_calls: list[str] = []
    gained_categories: set[str] = set()
    lost_categories: set[str] = set()
    growth_shift_lines: list[str] = []
    n_compared = 0

    for tissue in tissues:
        before = base_record.get(tissue)
        after = head_record.get(tissue)
        line: dict = {
            "before_tasks": "|".join(sorted(before["tasks"])) if before else ".",
            "after_tasks": "|".join(sorted(after["tasks"])) if after else ".",
            "before_growth": before["growth"] if before else None,
            "after_growth": after["growth"] if after else None,
        }
        if before is None or after is None:
            line["membership_changed"] = (before is None) != (after is None)
            per_line[tissue] = line
            continue

        line["membership_changed"] = False
        n_compared += 1
        any_dir = _direction(bool(before["tasks"]), bool(after["tasks"]))
        via_dir = _direction(
            bool(before["tasks"] & VIABILITY_CATEGORIES),
            bool(after["tasks"] & VIABILITY_CATEGORIES),
        )
        if any_dir:
            any_directions.append(any_dir)
            gained_categories |= after["tasks"] - before["tasks"]
            lost_categories |= before["tasks"] - after["tasks"]
        if via_dir:
            viability_directions.append(via_dir)
            call = _hart_call(via_dir, hart_fitness.get(tissue))
            if call:
                hart_calls.append(call)
            line["hart_fitness"] = hart_fitness.get(tissue)
            line["hart_call"] = call

        b_growth, a_growth = before["growth"], after["growth"]
        if b_growth is not None and a_growth is not None and abs(a_growth - b_growth) > growth_tolerance:
            growth_shift_lines.append(tissue)

        line["any_flip"] = any_dir
        line["viability_flip"] = via_dir
        per_line[tissue] = line

    def verdict(directions: list[str]) -> tuple[str, str]:
        if not directions:
            return "stable", ""
        gained = directions.count("gained")
        lost = directions.count("lost")
        # Denominator is every line comparable on both sides, not just the lines that
        # flipped -- "3/3" would otherwise look identical whether 3/5 or 3/3 lines moved.
        if max(gained, lost) >= CONSENSUS_MIN_LINES:
            return ("consensus gained" if gained >= lost else "consensus lost"), f"{max(gained, lost)}/{n_compared}"
        return "isolated (likely noise)", f"{len(directions)}/{n_compared}"

    any_verdict, any_count = verdict(any_directions)
    viability_verdict, viability_count = verdict(viability_directions)

    if not hart_calls:
        hart_verdict = "not scored by Hart 2015" if viability_directions else "n/a"
    elif len(set(hart_calls)) == 1:
        hart_verdict = hart_calls[0]
    else:
        hart_verdict = "mixed"

    return {
        "gene": gene_id,
        "per_line": per_line,
        "any_verdict": any_verdict,
        "any_count": any_count,
        "viability_verdict": viability_verdict,
        "viability_count": viability_count,
        "hart_verdict": hart_verdict,
        "gained_categories": gained_categories,
        "lost_categories": lost_categories,
        "growth_shift_lines": growth_shift_lines,
        "membership_changed_lines": [t for t, line in per_line.items() if line.get("membership_changed")],
    }


def build_report(
    base_model_path: Path,
    head_model_path: Path,
    base_matrix_path: Path,
    head_matrix_path: Path,
    *,
    growth_tolerance: float = 0.01,
    base_label: str = "the target branch",
    csv_url: str = "",
) -> tuple[str, str]:
    """Returns ``(summary_text, detail_csv_text)``.

    ``csv_url`` is the committed blob URL of the per-gene detail CSV this call is about
    to write (e.g. from ``--out-csv`` under ``data/testResults/``), linked from the
    summary instead of duplicating any of that detail into the summary text itself.
    """
    print(f"Loading base model from {base_model_path} ...", file=sys.stderr)
    base_model = read_yaml_model(base_model_path)
    print(f"Loading head model from {head_model_path} ...", file=sys.stderr)
    head_model = read_yaml_model(head_model_path)

    changed = changed_reactions(base_model, head_model)
    direct = set().union(*(info["genes"] for info in changed.values())) if changed else set()
    neighbors = neighbor_genes(head_model, changed, direct)

    tissues_b, symbol_b, per_gene_b = read_matrix(base_matrix_path)
    tissues_h, symbol_h, per_gene_h = read_matrix(head_matrix_path)
    tissues = [t for t in tissues_h if t in tissues_b] or tissues_h
    symbol_of = {**symbol_b, **symbol_h}

    scope = [(g, "direct") for g in sorted(direct)] + [(g, "neighbor") for g in sorted(neighbors)]
    scope_gene_ids = [g for g, _ in scope]

    print(f"Loading Hart 2015 fitness calls for {len(scope_gene_ids)} gene(s) in scope ...", file=sys.stderr)
    hart_bf = bayes_factors(scope_gene_ids, tissues, symbol_of=symbol_of)
    hart_fitness_by_gene = {
        gene_id: {
            tissue: hart_bf.get(tissue, {}).get(gene_id, None) is not None
            and hart_bf[tissue][gene_id] > BF_THRESHOLDS[tissue]
            for tissue in tissues
            if gene_id in hart_bf.get(tissue, {})
        }
        for gene_id in scope_gene_ids
    }

    results = []
    for gene_id, hop in scope:
        base_record = per_gene_b.get(gene_id, {})
        head_record = per_gene_h.get(gene_id, {})
        if not base_record and not head_record:
            continue  # gene not in either matrix (e.g. not built into any cell-line model)
        result = diff_gene(
            gene_id, tissues, base_record, head_record, hart_fitness_by_gene.get(gene_id, {}),
            growth_tolerance=growth_tolerance,
        )
        result["hop"] = hop
        result["symbol"] = symbol_of.get(gene_id, "")
        results.append(result)

    growth_relevant = [r for r in results if r["viability_verdict"].startswith("consensus")]
    capability_only = [
        r for r in results
        if r not in growth_relevant and r["any_verdict"].startswith("consensus")
    ]
    consensus_ids = {id(r) for r in growth_relevant + capability_only}
    isolated = [
        r for r in results if id(r) not in consensus_ids and (
            r["any_verdict"].startswith("isolated") or r["viability_verdict"].startswith("isolated")
        )
    ]
    stable = [r for r in results if id(r) not in consensus_ids and r not in isolated]

    n_added = sum(1 for i in changed.values() if i["kind"] == "added")
    n_removed = sum(1 for i in changed.values() if i["kind"] == "removed")
    n_modified = sum(1 for i in changed.values() if i["kind"] == "modified")
    summary_lines = [
        "### Gene essentiality: effect of this change",
        "",
        f"**{len(changed)}** reaction(s) changed vs `{base_label}` "
        f"(**{n_added}** added, **{n_removed}** removed, **{n_modified}** modified). "
        f"Checked **{len(direct)}** gene(s) in those reactions plus **{len(neighbors)}** more "
        f"that share a metabolite with one of them.",
        "",
        f"_A gene counts as changed below only if the flip agrees in direction across at least "
        f"{CONSENSUS_MIN_LINES} of the 5 cell-line models; fewer than that is grouped as likely "
        f"noise from ftINIT's run-to-run variability instead._",
        "",
    ]

    def _grouped_rows(rows: list[dict], key_fn) -> list[tuple[tuple, list[str]]]:
        """Group genes sharing the same outcome (same key) into one row's gene list.

        Several genes commonly move together for the same reason (e.g. the three ETF
        complex subunits), and listing each on its own line just repeats the same
        three cells three times; one row with a comma-joined gene list reads the same
        information without the repetition.
        """
        groups: dict[tuple, list[str]] = {}
        for r in rows:
            groups.setdefault(key_fn(r), []).append(r["symbol"] or r["gene"])
        return sorted(groups.items(), key=lambda kv: sorted(kv[1]))

    def _growth_key(r: dict) -> tuple[str, str, str]:
        direction = "gained" if r["viability_verdict"].endswith("gained") else "lost"
        change_text = "knockout now blocks growth" if direction == "gained" else "knockout no longer blocks growth"
        icon = {
            "improvement": ":white_check_mark: correct",
            "regression": ":x: wrong",
            "mixed": ":warning: mixed across lines",
        }.get(r["hart_verdict"], f":question: {r['hart_verdict']}")
        return r["viability_count"], change_text, icon

    def _capability_key(r: dict) -> tuple[str, str, str]:
        direction = "gained" if r["any_verdict"].endswith("gained") else "lost"
        change_text = "now required" if direction == "gained" else "no longer required"
        categories = r["gained_categories"] | r["lost_categories"]
        return r["any_count"], _describe_categories(categories), change_text

    def _named_detail(rows: list[dict], key_fn) -> str:
        """Genes named with what changed for them, grouped by identical outcome --
        only called for the two small, actionable categories (growth-relevant and
        capability-only); noise and unchanged genes are counted, never named here."""
        groups = _grouped_rows(rows, key_fn)
        parts = [f"{', '.join(sorted(genes))} ({', '.join(str(k) for k in key)})" for key, genes in groups]
        return "; ".join(parts)

    csv_ref = f"[gene-essential-diff.csv]({csv_url})" if csv_url else "the full detail CSV"

    table_rows = [
        ("Growth effect changed (vs Hart 2015)", len(growth_relevant), _named_detail(growth_relevant, _growth_key)),
        ("Other role changed (not Hart-comparable)", len(capability_only), _named_detail(capability_only, _capability_key)),
        (f"Likely noise (<{CONSENSUS_MIN_LINES}/5 lines)", len(isolated), f"see {csv_ref}" if isolated else ""),
        ("No change", len(stable), ""),
    ]
    summary_lines += ["| Category | Genes | Detail |", "| --- | --- | --- |"]
    summary_lines += [f"| {cat} | **{n}** | {detail or '--'} |" for cat, n, detail in table_rows]
    summary_lines.append("")
    summary_lines.append(f"Full per-gene, per-line detail (every gene checked, not just the ones named above): {csv_ref}.")
    summary = "\n".join(summary_lines) + "\n"

    detail_header = [
        "gene", "symbol", "hop", "any_verdict", "any_count", "viability_verdict", "viability_count", "hart_verdict",
    ]
    for tissue in tissues:
        detail_header += [f"{tissue}_before_tasks", f"{tissue}_after_tasks", f"{tissue}_before_growth", f"{tissue}_after_growth"]
    detail_rows = [",".join(detail_header)]
    for r in sorted(results, key=lambda r: (r["hop"], r["gene"])):
        row = [
            r["gene"], r["symbol"], r["hop"], r["any_verdict"], r["any_count"],
            r["viability_verdict"], r["viability_count"], r["hart_verdict"],
        ]
        for tissue in tissues:
            line = r["per_line"].get(tissue, {})
            row += [
                str(line.get("before_tasks", ".")),
                str(line.get("after_tasks", ".")),
                "" if line.get("before_growth") is None else f"{line['before_growth']:.6f}",
                "" if line.get("after_growth") is None else f"{line['after_growth']:.6f}",
            ]
        detail_rows.append(",".join(row))
    detail_csv = "\n".join(detail_rows) + "\n"

    return summary, detail_csv


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base-model", type=Path, required=True, help="Human-GEM.yml checked out at the target branch")
    parser.add_argument("--base-matrix", type=Path, required=True, help="gene-essential.csv committed on the target branch")
    parser.add_argument("--head-model", type=Path, default=DEFAULT_HEAD_MODEL, help="Human-GEM.yml for this branch (default: model/Human-GEM.yml)")
    parser.add_argument("--head-matrix", type=Path, default=DEFAULT_HEAD_MATRIX, help="gene-essential.csv for this branch (default: data/testResults/gene-essential.csv)")
    parser.add_argument("--growth-tolerance", type=float, default=0.01, help="minimum |growth ratio delta| to flag independent of task-based flips (default: 0.01)")
    parser.add_argument("--base-label", type=str, default="the target branch", help="name shown for the base branch, e.g. 'develop' (default: 'the target branch')")
    parser.add_argument("--csv-url", type=str, default="", help="committed blob URL of --out-csv, linked from the summary instead of duplicating its content there")
    parser.add_argument("--out-summary", type=Path, default=None, help="write the summary text here instead of only stdout")
    parser.add_argument("--out-csv", type=Path, default=None, help="write the per-gene per-line detail CSV here")
    args = parser.parse_args(argv)

    summary, detail_csv = build_report(
        args.base_model, args.head_model, args.base_matrix, args.head_matrix,
        growth_tolerance=args.growth_tolerance, base_label=args.base_label, csv_url=args.csv_url,
    )
    print(summary)
    if args.out_summary:
        args.out_summary.write_text(summary, encoding="utf-8")
    if args.out_csv:
        args.out_csv.write_text(detail_csv, encoding="utf-8")
    else:
        print(detail_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
