"""Regenerate the gene-essentiality summary from the committed per-gene matrix.

Reads data/testResults/gene-essential.csv (the per-gene essentiality matrix) and
recomputes data/testResults/gene-essential_summary.md by comparing to the Hart
2015 data. This does NOT re-run the expensive ftINIT / essential-gene solve, so
it finishes in seconds. Use it after changing only the evaluation logic (for
example the gene-ID matching from issue #970), instead of re-running the full
gene-essentiality workflow.

Needs neither Gurobi nor raven-toolbox.

Usage:
    python code/test/reevaluateGeneEssentiality.py
"""

import csv
import sys
from pathlib import Path

from evaluateHart2015Essentiality import evaluate_hart2015_essentiality, results_to_markdown

# Repository root: this file is code/test/reevaluateGeneEssentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MATRIX_CSV = REPO_ROOT / "data" / "testResults" / "gene-essential.csv"
SUMMARY_MD = REPO_ROOT / "data" / "testResults" / "gene-essential_summary.md"


def main() -> int:
    if not MATRIX_CSV.exists():
        print(f"::error::{MATRIX_CSV} not found; run the gene-essentiality workflow first.")
        return 1

    with open(MATRIX_CSV, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fixed = {"genes", "geneSymbol"}
        tissues = [c for c in (reader.fieldnames or []) if c not in fixed]
        gene_ids: list[str] = []
        symbol_of: dict[str, str] = {}
        essential_ids_by_tissue: dict[str, set[str]] = {t: set() for t in tissues}
        for row in reader:
            gid = row["genes"]
            gene_ids.append(gid)
            symbol_of[gid] = row.get("geneSymbol", "") or ""
            for t in tissues:
                if (row.get(t) or "").strip().lower() == "yes":
                    essential_ids_by_tissue[t].add(gid)

    results = evaluate_hart2015_essentiality(
        gene_ids, tissues, essential_ids_by_tissue, symbol_of=symbol_of
    )
    SUMMARY_MD.write_text(results_to_markdown(results), encoding="utf-8")
    print(results.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
