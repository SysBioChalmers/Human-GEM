"""Compare predicted gene essentiality with Hart 2015 experimental fitness data.

Python port of code/test/evaluateHart2015Essentiality.m. Given the genes predicted
essential in each cell-line-specific model (see estimateEssentialGenes.py), this
computes TP/TN/FP/FN and the derived metrics (accuracy, sensitivity, specificity,
F1, MCC) against the Hart et al. 2015 CRISPR fitness screen.

Gene identifier matching (see issue #970)
-----------------------------------------
The experimental table (Hart2015_TableS2.xlsx) uses gene symbols, while the model
uses Ensembl gene ids. Mapping each model gene to its single symbol loses matches,
because the symbol the model carries is often not the one Hart uses (for example
the model's ACP3 vs Hart's ACPP, both aliases of ENSG00000014257). Following the
recommendation in issue #970, the comparison is done in Ensembl space: each Hart
symbol is matched to a model gene through that gene's symbol *and* aliases (from
model/genes.tsv), so alias differences no longer drop true matches. This agrees
with the issue author's manually curated Hart-to-Ensembl mapping at Jaccard
0.96-0.98 per cell line.

Experimental data: Table S2 from Hart 2015. By definition in that study essential
genes are a subset of "fitness" genes; a gene is a fitness gene in a cell line
when its Bayes Factor exceeds the 5% FDR threshold reported in the paper.
"""

from __future__ import annotations

import csv
import math
import re
import zipfile
from collections.abc import Iterable, Mapping
from pathlib import Path

import pandas as pd

# Repository root: this file is code/test/evaluateHart2015Essentiality.py
REPO_ROOT = Path(__file__).resolve().parents[2]
HART_TABLE_S2 = REPO_ROOT / "data" / "datasets" / "Hart2015_TableS2.xlsx"
GENES_TSV = REPO_ROOT / "model" / "genes.tsv"

# Bayes-Factor thresholds (5% FDR) per cell line, from the Hart 2015 supporting
# information. Keys match the upper-cased BF_* column names of Table S2.
BF_THRESHOLDS = {
    "HCT116": 1.57,
    "HELA": 15.47,
    "GBM": 3.20,
    "RPE1": 6.84,
    "DLD1": 3.57,
}

# Metric columns, in the order written to the summary.
RESULT_COLUMNS = [
    "cellLine", "TP", "TN", "FP", "FN",
    "accuracy", "sensitivity", "specificity", "F1", "MCC",
]


def _read_xlsx_rows(path: str | Path) -> tuple[list[str], list[dict[str, object]]]:
    """Minimal single-sheet .xlsx reader (shared strings + numeric cells).

    Returns ``(headers, rows)`` where each row is a dict keyed by the first-row
    headers, numeric cells as floats, text cells as strings, empty cells as None.
    This is all the Hart table needs and avoids an openpyxl dependency, which also
    lets the summary be regenerated (reevaluateGeneEssentiality.py) without it.
    """
    with zipfile.ZipFile(path) as z:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = z.read("xl/sharedStrings.xml").decode("utf-8", "replace")
            shared = ["".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
                      for si in re.findall(r"<si>(.*?)</si>", ss, re.S)]
        sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8", "replace")

    def col_index(ref: str) -> int:
        letters = re.match(r"[A-Z]+", ref).group(0)
        n = 0
        for ch in letters:
            n = n * 26 + (ord(ch) - 64)
        return n - 1

    raw_rows: list[dict[int, object]] = []
    for row_xml in re.findall(r"<row[^>]*>.*?</row>", sheet, re.S):
        cells: dict[int, object] = {}
        # Parse the cell's r (reference) and t (type) attributes independently of
        # their order, since a style attribute (s="..") can sit between them.
        for m in re.finditer(r"<c\b([^>]*)>(?:<v>(.*?)</v>)?</c>", row_xml):
            attrs, v = m.group(1), m.group(2)
            ref = re.search(r'r="([A-Z]+)\d+"', attrs)
            if ref is None:
                continue
            ci = col_index(ref.group(1))
            type_attr = re.search(r't="(\w+)"', attrs)
            if v is None:
                cells[ci] = None
            elif type_attr is not None and type_attr.group(1) == "s":
                cells[ci] = shared[int(v)]
            else:
                try:
                    cells[ci] = float(v)
                except ValueError:
                    cells[ci] = v
        raw_rows.append(cells)

    if not raw_rows:
        return [], []
    width = max((max(c) + 1) if c else 0 for c in raw_rows)
    headers = [str(raw_rows[0].get(i, "")) for i in range(width)]
    rows = [{headers[i]: c.get(i) for i in range(width)} for c in raw_rows[1:]]
    return headers, rows


def _load_hart2015(table_path: str | Path = HART_TABLE_S2) -> dict[str, dict[str, set[str]]]:
    """Load experimental fitness genes per cell line from Hart 2015 Table S2.

    Returns ``{cell_line: {"scored": set(symbols), "essential": set(symbols)}}``
    for the five cell lines. Symbols are upper case. The ``all`` category is
    derived later, in Ensembl space.
    """
    headers, rows = _read_xlsx_rows(table_path)

    # Remove the duplicated MARCH1 / MARCH2 rows (drop the last occurrence of each),
    # matching the MATLAB implementation.
    for dup_gene in ("MARCH1", "MARCH2"):
        idxs = [i for i, r in enumerate(rows) if str(r.get("Gene") or "").upper() == dup_gene]
        if idxs:
            del rows[idxs[-1]]

    # The five cell-line columns are the BF_* headers whose cell line has a threshold
    # (this excludes BF_a375_GeCKo and BF_hct116_shRNA).
    bf_headers = [
        h for h in headers
        if h.upper().startswith("BF_") and h.upper().replace("BF_", "") in BF_THRESHOLDS
    ]
    result: dict[str, dict[str, set[str]]] = {}
    for header in bf_headers:
        cell_line = header.upper().replace("BF_", "")
        threshold = BF_THRESHOLDS[cell_line]
        scored: set[str] = set()
        essential: set[str] = set()
        for r in rows:
            bf = r.get(header)
            if isinstance(bf, float) and not math.isnan(bf):
                gene = str(r.get("Gene") or "").upper()
                scored.add(gene)
                if bf > threshold:
                    essential.add(gene)
        result[cell_line] = {"scored": scored, "essential": essential}
    return result


def _model_gene_names(genes_tsv: str | Path = GENES_TSV) -> dict[str, set[str]]:
    """Map model gene id -> set of upper-case symbols and aliases from genes.tsv.

    Both the primary ``geneSymbols`` and the semicolon-separated ``geneAliases``
    are included so a Hart symbol matches whichever alias the model gene carries.
    """
    names: dict[str, set[str]] = {}
    with open(genes_tsv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene_id = row["genes"]
            entry = set()
            symbol = (row.get("geneSymbols") or "").strip().upper()
            if symbol:
                entry.add(symbol)
            aliases = row.get("geneAliases") or ""
            entry |= {a.strip().upper() for a in aliases.split(";") if a.strip()}
            names[gene_id] = entry
    return names


def _safe_div(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else math.nan


def _metrics(tp: int, tn: int, fp: int, fn: int) -> dict[str, float]:
    sensitivity = _safe_div(tp, tp + fn)
    specificity = _safe_div(tn, tn + fp)
    accuracy = _safe_div(tp + tn, tp + tn + fp + fn)
    f1 = _safe_div(2 * tp, 2 * tp + fp + fn)
    mcc_denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = _safe_div((tp * tn) - (fp * fn), mcc_denom)
    return {
        "accuracy": accuracy,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "F1": f1,
        "MCC": mcc,
    }


def experimental_status(
    ref_gene_ids: Iterable[str],
    tissues: list[str],
    *,
    symbol_of: Mapping[str, str] | None = None,
    genes_tsv: str | Path = GENES_TSV,
    table_path: str | Path = HART_TABLE_S2,
) -> dict[str, dict[str, set[str]]]:
    """Map Hart 2015 fitness data into model-gene (Ensembl) space, per cell line.

    Returns ``{cell_line: {"scored": set(gene_id), "essential": set(gene_id)}}`` (plus an
    ``"all"`` key = intersection across cell lines). A model gene is *scored* / *essential*
    in a cell line when any of its symbols or aliases is scored / a fitness gene there.
    """
    ref = set(ref_gene_ids)
    hart = _load_hart2015(table_path)
    gene_names = _model_gene_names(genes_tsv)

    def names_for(gene_id: str) -> set[str]:
        entry = set(gene_names.get(gene_id, ()))
        if symbol_of and symbol_of.get(gene_id):
            entry.add(symbol_of[gene_id].upper())
        return entry

    names = {gene_id: names_for(gene_id) for gene_id in ref}
    status: dict[str, dict[str, set[str]]] = {}
    for cell_line in tissues:
        scored_syms = hart[cell_line]["scored"]
        essential_syms = hart[cell_line]["essential"]
        status[cell_line] = {
            "scored": {g for g in ref if names[g] & scored_syms},
            "essential": {g for g in ref if names[g] & essential_syms},
        }
    if tissues:
        status["all"] = {
            "scored": set.intersection(*(status[c]["scored"] for c in tissues)),
            "essential": set.intersection(*(status[c]["essential"] for c in tissues)),
        }
    return status


def evaluate_hart2015_essentiality(
    ref_gene_ids: Iterable[str],
    tissues: list[str],
    essential_ids_by_tissue: Mapping[str, Iterable[str]],
    *,
    symbol_of: Mapping[str, str] | None = None,
    genes_tsv: str | Path = GENES_TSV,
    table_path: str | Path = HART_TABLE_S2,
) -> pd.DataFrame:
    """Evaluate predicted essentiality against Hart 2015 fitness data in Ensembl space.

    Parameters
    ----------
    ref_gene_ids
        All model gene ids (Ensembl); the universe of genes.
    tissues
        Cell-line names, in output order (e.g. ["DLD1", "GBM", "HCT116", "HELA",
        "RPE1"]). Must correspond to the Hart 2015 cell lines.
    essential_ids_by_tissue
        Mapping cell-line -> set of model gene ids predicted essential for any task.
    symbol_of
        Optional gene id -> symbol map, added to the genes.tsv names as a fallback.

    Returns a DataFrame with one row per cell line plus an "all" row, and columns
    RESULT_COLUMNS.
    """
    ref = set(ref_gene_ids)
    status = experimental_status(
        ref_gene_ids, tissues, symbol_of=symbol_of, genes_tsv=genes_tsv, table_path=table_path
    )
    exp_scored = {c: status[c]["scored"] for c in status}
    exp_essential = {c: status[c]["essential"] for c in status}

    pred = {cell_line: set(essential_ids_by_tissue.get(cell_line, ())) & ref for cell_line in tissues}
    pred["all"] = set.intersection(*pred.values()) if pred else set()

    rows = []
    for cell_line in [*tissues, "all"]:
        if cell_line not in exp_essential:
            continue
        model_essential = pred[cell_line]
        model_non_essential = ref - model_essential

        experimentally_essential = exp_essential[cell_line]
        experimentally_non_essential = exp_scored[cell_line] - experimentally_essential

        tp = len(model_essential & experimentally_essential)
        tn = len(model_non_essential & experimentally_non_essential)
        fp = len(model_essential & experimentally_non_essential)
        fn = len(model_non_essential & experimentally_essential)

        row = {"cellLine": cell_line, "TP": tp, "TN": tn, "FP": fp, "FN": fn}
        row.update(_metrics(tp, tn, fp, fn))
        rows.append(row)

    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def results_to_markdown(results: pd.DataFrame) -> str:
    """Render the summary metrics as a Markdown table (counts as integers, metric
    columns to 4 significant figures)."""
    metric_cols = ("accuracy", "sensitivity", "specificity", "F1", "MCC")
    lines = [
        "### Gene essentiality vs Hart 2015 fitness genes",
        "",
        "| " + " | ".join(RESULT_COLUMNS) + " |",
        "| " + " | ".join("---" for _ in RESULT_COLUMNS) + " |",
    ]
    for _, r in results.iterrows():
        cells = [
            str(r["cellLine"]),
            str(int(r["TP"])), str(int(r["TN"])), str(int(r["FP"])), str(int(r["FN"])),
        ]
        cells.extend(f"{r[col]:.4g}" for col in metric_cols)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def essentiality_matrix_csv(
    gene_ids: list[str],
    symbol_of: Mapping[str, str],
    tissues: list[str],
    essential_ids_by_tissue: Mapping[str, Iterable[str]],
    experimental: Mapping[str, Mapping[str, Iterable[str]]],
) -> str:
    """Build the per-gene prediction-vs-Hart2015 table as CSV text.

    One row per gene, one column per cell line, holding the confusion class of the
    prediction against the Hart 2015 fitness call:

        TP  predicted essential  & Hart fitness gene       (correct)
        TN  predicted non-essential & Hart non-fitness     (correct)
        FP  predicted essential  & Hart non-fitness        (wrong: over-predicted)
        FN  predicted non-essential & Hart fitness gene    (wrong: missed)
        P   predicted essential, gene not scored by Hart   (unvalidated positive)
        .   predicted non-essential, gene not scored       (uninformative)

    Because Hart's truth is fixed, a cell changing between two runs is exactly a
    prediction change, with its direction built in: gaining ``FN``/``FP`` is a
    regression, gaining ``TP``/``TN`` an improvement. Every model gene is kept (one
    stable row per gene, in the given order), so diffs show only changed cells.
    """
    pred = {t: set(essential_ids_by_tissue.get(t, ())) for t in tissues}
    ess = {t: set(experimental.get(t, {}).get("essential", ())) for t in tissues}
    scored = {t: set(experimental.get(t, {}).get("scored", ())) for t in tissues}

    def klass(gid: str, t: str) -> str:
        predicted = gid in pred[t]
        if gid not in scored[t]:
            return "P" if predicted else "."
        fitness = gid in ess[t]
        return ("TP" if fitness else "FP") if predicted else ("FN" if fitness else "TN")

    header = ["genes", "geneSymbol", *tissues]
    lines = [",".join(header)]
    for gid in gene_ids:
        row = [gid, symbol_of.get(gid, ""), *(klass(gid, t) for t in tissues)]
        lines.append(",".join(row))
    return "\n".join(lines) + "\n"
