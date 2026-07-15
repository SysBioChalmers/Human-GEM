"""Gene-essentiality CI entry point.

Python port of the inline MATLAB in .github/workflows/gene-essentiality.yml.
Builds cell-line-specific tINIT models from the Hart 2015 RNA-seq data, predicts
gene essentiality for the essential metabolic tasks, and writes two artifacts to
data/testResults/:

  * gene-essential.csv          per-gene prediction-vs-Hart2015 table: one row per
                                model gene, one column per cell line holding the
                                confusion class (TP/TN/FP/FN, or P/. when Hart did not
                                score the gene). Every gene is kept, so run-to-run diffs
                                show only changed cells (gaining FN/FP is a regression).
  * gene-essential_summary.md   the summary statistics of the comparison against
                                the Hart 2015 CRISPR fitness screen, as a Markdown
                                table (used in the PR comment).

Requires Gurobi (the tINIT MILP is genome-scale). Usage:
    python code/test/geneEssentiality.py
"""

import sys
from pathlib import Path

from raven_toolbox.io import read_yaml_model

from estimateEssentialGenes import estimate_essential_genes
from evaluateHart2015Essentiality import (
    essentiality_matrix_csv,
    evaluate_hart2015_essentiality,
    experimental_status,
    results_to_markdown,
)

# Repository root: this file is code/test/geneEssentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"
RESULTS_DIR = REPO_ROOT / "data" / "testResults"
MATRIX_CSV = RESULTS_DIR / "gene-essential.csv"
SUMMARY_MD = RESULTS_DIR / "gene-essential_summary.md"


def main() -> int:
    # Silence Gurobi's console output. Copying a genome-scale model (which
    # raven-toolbox does per task while preparing the template, and again per
    # cell line) round-trips it through a temporary LP file, and Gurobi prints a
    # "Read LP format model from file ..." banner every time, flooding the log.
    # Setting OutputFlag on the default environment before any model is created
    # suppresses that noise while keeping our own progress lines.
    import gurobipy
    gurobipy.setParam("OutputFlag", 0)

    model = read_yaml_model(MODEL_FILE)
    # tINIT solves a genome-scale MILP, which needs Gurobi.
    model.solver = "gurobi"

    egenes = estimate_essential_genes(model)
    tissues = egenes["tissues"]
    symbol_of = egenes["symbol_of"]
    essential_ids_by_tissue = egenes["essential_ids_by_tissue"]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Hart 2015 experimental fitness genes mapped into model-gene (Ensembl) space, per
    # cell line (see issue #970): Hart's gene symbols are matched to model genes through
    # their symbols and aliases, so essential gene ids are compared directly rather than
    # collapsed to a single symbol.
    experimental = experimental_status(egenes["gene_ids"], tissues, symbol_of=symbol_of)

    # Detailed per-gene prediction-vs-Hart2015 table (the diffable artifact).
    MATRIX_CSV.write_text(
        essentiality_matrix_csv(
            egenes["gene_ids"], symbol_of, tissues, essential_ids_by_tissue, experimental
        ),
        encoding="utf-8",
    )

    # Summary statistics against the Hart 2015 data.
    results = evaluate_hart2015_essentiality(
        egenes["gene_ids"], tissues, essential_ids_by_tissue, symbol_of=symbol_of
    )
    SUMMARY_MD.write_text(results_to_markdown(results), encoding="utf-8")

    # Echo the summary to the log.
    print(results.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
