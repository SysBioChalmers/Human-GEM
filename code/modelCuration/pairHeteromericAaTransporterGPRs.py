"""Pair SLC7 light chains with their SLC3 heavy chain in the transporter GPRs.

Three stages, each including the previous ones:
  stage1  removes SLC7A6 from the b0,+ rule, drops a heavy chain that appears as a
          standalone OR alternative, and drops SLC7A9 from the bile-acid/GSH exchanges.
  stage2  AND-s each SLC7A5-A11 light chain with its heavy chain.
  stage3  matches the light chains to the substrates and the Na+ coupling of each
          amino acid exchange:
          - a Na+-coupled exchange of two non-cationic amino acids loses SLC7A5 and
            SLC7A6, because neither couples a neutral/neutral exchange to Na+;
          - an exchange without Na+ that moves lysine, arginine or ornithine loses
            SLC7A5, because LAT1 has no cationic substrate;
          - a Na+-coupled exchange of a cationic amino acid that carries SLC7A6 gains
            SLC7A7 (y+LAT1), the second y+L light chain;
          - MAR09402 (kynurenine/phenylalanine exchange) gains SLC7A8 (LAT2).
          A rule that loses all its genes is removed from the reaction.

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
SODIUM = "MAM02519"
AMINO_ACIDS = {
    "glycine", "alanine", "serine", "threonine", "cysteine", "methionine", "valine",
    "leucine", "isoleucine", "proline", "phenylalanine", "tyrosine", "tryptophan",
    "asparagine", "glutamine", "histidine", "citrulline", "homoserine",
    "lysine", "arginine", "ornithine",
}
CATIONIC = {"lysine", "arginine", "ornithine"}
EXTRA_TERMS = {"MAR09402": ["SLC7A8", "SLC3A2"]}
# references added with a gene: y+LAT1 function (Torrents 1998) and LAT2-dependent
# kynurenine uptake (Liu 2018, Zhou 2023)
SLC7A7_REFS = ["doi:10.1074/jbc.273.49.32437"]
EXTRA_REFS = {"MAR09402": ["doi:10.1016/j.ccell.2018.02.005", "doi:10.1038/s41590-023-01662-3"]}


def classify_exchange(mets, met_names):
    """Return (sodium_coupled, cationic) for an amino acid exchange, else None.

    An exchange moves exactly two amino acids across the membrane in opposite
    directions, optionally with Na+ entering the cell.
    """
    moved = {}
    for met, coef in mets.items():
        moved.setdefault(met[:-1], {})[met[-1]] = coef
    sodium = moved.pop(SODIUM, None)
    crossing = {b: c for b, c in moved.items() if len(c) == 2}
    if len(crossing) != 2 or len(crossing) != len(moved):
        return None
    names = {met_names.get(b) for b in crossing}
    if not names <= AMINO_ACIDS:
        return None
    if {c.get("e", 0) for c in crossing.values()} != {-1.0, 1.0}:
        return None
    sodium_coupled = sodium is not None and sodium.get("e") == -1.0 and sodium.get("c") == 1.0
    return sodium_coupled, bool(names & CATIONIC)


def match_substrates(rxn_id, terms, sym, rev, exchange):
    """Apply the stage-3 substrate rules to the terms of one rule."""
    out, notes = [], []

    def names_of(genes):
        return {sym.get(g, g) for g in genes}

    if exchange is not None:
        sodium_coupled, cationic = exchange
        for genes in terms:
            names = names_of(genes)
            if sodium_coupled and not cationic and names & {"SLC7A5", "SLC7A6"}:
                notes.append("SLC7A5/SLC7A6 from Na+-coupled neutral exchange")
                continue
            if not sodium_coupled and cationic and "SLC7A5" in names:
                notes.append("SLC7A5 from cationic exchange")
                continue
            out.append(genes)
        if sodium_coupled and cationic and any("SLC7A6" in names_of(g) for g in out) \
                and not any("SLC7A7" in names_of(g) for g in out):
            out.append([rev["SLC7A7"], rev["SLC3A2"]])
            notes.append("SLC7A7 added to y+L exchange")
    else:
        out = [list(g) for g in terms]
    extra = EXTRA_TERMS.get(rxn_id)
    if extra and not any(names_of(g) == set(extra) for g in out):
        out.append([rev[n] for n in extra])
        notes.append(f"{extra[0]} added")
    return out, notes


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
        if stage >= 2:
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
    ap.add_argument("--stage", type=int, choices=(1, 2, 3), required=True)
    ap.add_argument("--report", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    sym, rev = load_symbols()
    targets = set(LIGHT_HEAVY) | HEAVY
    target_ids = {rev[n] for n in targets if n in rev}

    raw = YML.read_bytes().decode("utf-8")
    lines = raw.split("\n")  # keeps the \r at the end of each line

    rxn_id, changes, refused, removed = None, [], [], set()
    met_names, mets, pending_refs = {}, {}, {}
    for idx, line in enumerate(lines):
        body = line.rstrip("\r")
        m = re.match(r"^    - id: (MAR\d+|MAM\d+[a-z])$", body)
        if m:
            rxn_id, mets = m.group(1), {}
            continue
        m = re.match(r"^(    - references: )(.*)$", body)
        if m and rxn_id in pending_refs:
            refs = m.group(2).split(";")
            refs += [r for r in pending_refs.pop(rxn_id) if r not in refs]
            lines[idx] = m.group(1) + ";".join(refs) + ("\r" if line.endswith("\r") else "")
            continue
        m = re.match(r"^    - name: (.*)$", body)
        if m and rxn_id and rxn_id.startswith("MAM"):
            met_names.setdefault(rxn_id[:-1], m.group(1).strip().strip('"'))
            continue
        m = re.match(r"^      - (MAM\d+[a-z]): (\S+)$", body)
        if m:
            mets[m.group(1)] = float(m.group(2))
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
        if args.stage >= 3:
            new_terms, notes = match_substrates(
                rxn_id, new_terms, sym, rev, classify_exchange(mets, met_names))
            dropped += notes
        if not new_terms and args.stage < 3:
            refused.append((rxn_id, f"EMPTIED: {rule}"))
            continue
        new_rule = render(new_terms)
        if new_rule != rule:
            changes.append((rxn_id, rule, new_rule, "; ".join(sorted(set(dropped)))))
            if new_rule:
                lines[idx] = prefix + new_rule + ("\r" if line.endswith("\r") else "")
            else:
                removed.add(idx)
            refs = EXTRA_REFS.get(rxn_id, []) + (
                SLC7A7_REFS if "SLC7A7 added to y+L exchange" in dropped else [])
            if refs:
                pending_refs[rxn_id] = refs

    lines = [line for idx, line in enumerate(lines) if idx not in removed]

    if pending_refs:
        refused += [(rid, "no references field for the added gene's reference")
                    for rid in pending_refs]
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
            w = csv.writer(fh, delimiter="\t" if args.report.suffix == ".tsv" else ",",
                           lineterminator="\n")
            w.writerow(["reaction", "before", "after", "before_symbols", "after_symbols", "note"])
            for rid, old, new, note in changes:
                w.writerow([rid, old, new, pretty(old), pretty(new), note])
        print(f"report written to {args.report}")

    if not args.dry_run:
        YML.write_bytes("\n".join(lines).encode("utf-8"))
        print(f"{YML} updated")


if __name__ == "__main__":
    main()
