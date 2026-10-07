"""Pair SLC7 light chains with their SLC3 heavy chain in the transporter GPRs.

Two stages:
  stage1  removes SLC7A6 from the b0,+ rule, drops a heavy chain that appears as a
          standalone OR alternative, and drops SLC7A9 from the bile-acid/GSH exchanges.
  stage2  stage1 plus AND-ing each SLC7A5-A11 light chain with its heavy chain.

Edits only the gene_reaction_rule lines of model/Human-GEM.yml, in place, so the rest
of the file (and its CRLF line endings) is untouched.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

YML = Path("model/Human-GEM.yml")
GENES = Path("model/genes.tsv")

LIGHT_HEAVY = {
    "SLC7A5": "SLC3A2", "SLC7A6": "SLC3A2", "SLC7A7": "SLC3A2", "SLC7A8": "SLC3A2",
    "SLC7A9": "SLC3A1", "SLC7A10": "SLC3A2", "SLC7A11": "SLC3A2",
}
HEAVY = {"SLC3A1", "SLC3A2"}
BILE = {"MAR01879", "MAR01892", "MAR01893"}


def load_symbols():
    sym = {}
    with GENES.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym[row["genes"]] = row["geneSymbols"]
    return sym, {v: k for k, v in sym.items()}


def split_top(expr, op):
    """Split on a top-level ' <op> ', leaving bracketed sub-expressions intact."""
    pat, out, buf, depth, i = f" {op} ", [], "", 0, 0
    while i < len(expr):
        ch = expr[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if depth == 0 and expr.startswith(pat, i):
            out.append(buf)
            buf = ""
            i += len(pat)
            continue
        buf += ch
        i += 1
    out.append(buf)
    return out


def unwrap(term):
    term = term.strip()
    while term.startswith("(") and term.endswith(")") and split_top(term[1:-1], "or") and len(
            split_top(term[1:-1], "or")) >= 1 and term.count("(") == term.count(")"):
        inner = term[1:-1]
        if inner.count("(") != inner.count(")"):
            break
        term = inner.strip()
        if not (term.startswith("(") and term.endswith(")")):
            break
    return term


def parse_flat(rule):
    """Return [[ENSG, ...], ...] for a flat OR-of-ANDs rule, or None if not flat."""
    terms = []
    for raw in split_top(rule, "or"):
        term = unwrap(raw)
        if len(split_top(term, "or")) != 1:
            return None
        genes = []
        for piece in split_top(term, "and"):
            piece = piece.strip()
            if not re.fullmatch(r"ENSG\d+", piece):
                return None
            genes.append(piece)
        terms.append(genes)
    return terms


def render(terms):
    parts = []
    multi = len(terms) > 1
    for genes in terms:
        joined = " and ".join(sorted(genes))
        parts.append(f"({joined})" if multi and len(genes) > 1 else joined)
    return " or ".join(parts)


def curate(rxn_id, terms, sym, rev, stage):
    out, dropped = [], []
    for genes in terms:
        names = [sym.get(g, g) for g in genes]
        if len(genes) == 1 and names[0] in HEAVY:
            dropped.append(f"standalone {names[0]}")
            continue
        keep = list(genes)
        if "SLC7A9" in names and "SLC7A6" in names:
            keep = [g for g in keep if sym.get(g) != "SLC7A6"]
            dropped.append("SLC7A6 from b0,+ complex")
        if rxn_id in BILE:
            keep = [g for g in keep if sym.get(g) != "SLC7A9"]
            dropped.append("SLC7A9 from bile-acid/GSH exchange")
        if not keep:
            continue
        if stage == 2:
            kn = [sym.get(g, g) for g in keep]
            for name in list(kn):
                heavy = LIGHT_HEAVY.get(name)
                if heavy and heavy not in kn:
                    keep.append(rev[heavy])
                    kn.append(heavy)
        out.append(keep)
    seen, uniq = set(), []
    for genes in out:
        key = tuple(sorted(genes))
        if key not in seen:
            seen.add(key)
            uniq.append(genes)
    return uniq, dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, choices=(1, 2), required=True)
    ap.add_argument("--report", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    sym, rev = load_symbols()
    targets = set(LIGHT_HEAVY) | HEAVY
    target_ids = {rev[n] for n in targets if n in rev}

    raw = YML.read_bytes().decode("utf-8")
    lines = raw.split("\n")  # keeps the \r at the end of each line

    rxn_id, changes, refused = None, [], []
    for idx, line in enumerate(lines):
        body = line.rstrip("\r")
        m = re.match(r"^    - id: (MA[RM]\d+)$", body)
        if m:
            rxn_id = m.group(1)
            continue
        m = re.match(r"^(    - gene_reaction_rule: )(.*)$", body)
        if not m:
            continue
        prefix, rule = m.groups()
        if not rule or not (set(re.findall(r"ENSG\d+", rule)) & target_ids):
            continue
        terms = parse_flat(rule)
        if terms is None:
            refused.append((rxn_id, rule))
            continue
        new_terms, dropped = curate(rxn_id, terms, sym, rev, args.stage)
        if not new_terms:
            refused.append((rxn_id, f"EMPTIED: {rule}"))
            continue
        new_rule = render(new_terms)
        if new_rule != rule:
            changes.append((rxn_id, rule, new_rule, "; ".join(sorted(set(dropped)))))
            lines[idx] = prefix + new_rule + ("\r" if line.endswith("\r") else "")

    if refused:
        print(f"REFUSED {len(refused)} rule(s) that are not a flat OR of ANDs:", file=sys.stderr)
        for rid, rule in refused:
            print(f"  {rid}: {rule[:140]}", file=sys.stderr)
        sys.exit(1)

    def pretty(rule):
        return re.sub(r"ENSG\d+", lambda mm: sym.get(mm.group(0), mm.group(0)), rule)

    print(f"stage {args.stage}: {len(changes)} GPRs changed")
    if args.report:
        with args.report.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["reaction", "before", "after", "before_symbols", "after_symbols", "note"])
            for rid, old, new, note in changes:
                w.writerow([rid, old, new, pretty(old), pretty(new), note])
        print(f"report written to {args.report}")

    if not args.dry_run:
        YML.write_bytes("\n".join(lines).encode("utf-8"))
        print(f"{YML} updated")


if __name__ == "__main__":
    main()
