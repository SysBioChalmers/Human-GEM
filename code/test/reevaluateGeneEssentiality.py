"""Regenerate the gene-essentiality summary from the committed per-gene matrix.

Reads data/testResults/gene-essential.csv (the per-gene matrix written by
gradedEssentiality.py) and recomputes data/testResults/gene-essential_summary.md by
comparing to the Hart 2015 data. This does NOT re-run the expensive ftINIT /
essential-gene solve, so it finishes in seconds. Use it after changing only the
evaluation logic (for example the gene-ID matching from issue #970, or the split
between viability and capability tasks from issue #1076), instead of re-running the
full gene-essentiality workflow.

The matrix carries, per cell line, the task categories each knockout breaks and the
biomass growth ratio, which is everything the scoring needs; only the Hart comparison
is redone here.

Needs neither Gurobi nor raven-toolbox.

Usage:
    python code/test/reevaluateGeneEssentiality.py
"""

import csv
import sys
from pathlib import Path

from evaluateHart2015Essentiality import evaluate_graded, summary_to_markdown

# Repository root: this file is code/test/reevaluateGeneEssentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MATRIX_CSV = REPO_ROOT / "data" / "testResults" / "gene-essential.csv"
SUMMARY_MD = REPO_ROOT / "data" / "testResults" / "gene-essential_summary.md"


def read_matrix(path: Path) -> tuple[list[str], list[str], dict[str, str], dict[str, dict[str, dict]]]:
    """Parse the per-gene matrix back into the structure ``evaluate_graded`` expects.

    Returns ``(tissues, gene_ids, symbol_of, per_tissue)``. A gene missing from a
    cell-line model is written as ``.`` and is left out of that cell line.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        tissues = [c[: -len("_class")] for c in columns if c.endswith("_class")]
        gene_ids: list[str] = []
        symbol_of: dict[str, str] = {}
        per_tissue: dict[str, dict[str, dict]] = {t: {} for t in tissues}
        for row in reader:
            gene_id = row["genes"]
            gene_ids.append(gene_id)
            symbol_of[gene_id] = row.get("geneSymbol", "") or ""
            for tissue in tissues:
                tasks = (row.get(f"{tissue}_tasks") or "").strip()
                if tasks == ".":
                    continue  # gene not in this cell-line model
                growth = (row.get(f"{tissue}_growth") or "").strip()
                per_tissue[tissue][gene_id] = {
                    "tasks": tasks.split("|") if tasks else [],
                    "growth": None if growth in ("", ".") else float(growth),
                }
    return tissues, gene_ids, symbol_of, per_tissue


def main() -> int:
    if not MATRIX_CSV.exists():
        print(f"::error::{MATRIX_CSV} not found; run the gene-essentiality workflow first.")
        return 1

    tissues, gene_ids, symbol_of, per_tissue = read_matrix(MATRIX_CSV)
    rows, _matrix_csv = evaluate_graded(per_tissue, tissues, gene_ids, symbol_of=symbol_of)
    summary = summary_to_markdown(rows)
    SUMMARY_MD.write_text(summary, encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
