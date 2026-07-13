"""Phase 2: regenerate protonation-mismatched structures to the model's pH-7.3
microspecies.

For every metabolite whose stored SMILES has the right skeleton but the wrong
protonation state (the "protonation" class from identity.py), enumerate
protonation states with dimorphite-dl and keep the one whose element formula AND
net charge exactly match the curated model formula/charge. The model's
formula/charge is authoritative (the curated pH-7.3 convention); we only pick the
protonation of the correct skeleton that reproduces it.

Writes data/annotation/protonation_resolved.tsv with, per distinct structure:
mets, name, model_formula, model_charge, old_smiles, new_smiles, new_inchi,
status (resolved / ambiguous / unresolved). Applying it to metabolites.tsv is a
separate, reviewed step.
"""

import csv
import re
import sys
from collections import Counter

from dimorphite_dl import protonate_smiles
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

IDENTITY = "data/annotation/metabolite_identity.tsv"
OUT = "data/annotation/protonation_resolved.tsv"
_EL = re.compile(r"([A-Z][a-z]?)(\d*)")


def elems(formula: str) -> dict:
    d = {}
    for sym, num in _EL.findall(formula.strip().rstrip("+-")):
        if sym:
            d[sym] = d.get(sym, 0) + (int(num) if num else 1)
    return d


def chem_score(mol) -> int:
    """Reward the standard pH-7.3 microspecies: negative charge on carboxylate /
    phosphate / sulfate oxygens and ammonium on aliphatic amines; penalize
    alkoxide / phenolate (deprotonated alcohols) and other implausible states."""
    s = 0
    for a in mol.GetAtoms():
        sym, q = a.GetSymbol(), a.GetFormalCharge()
        if sym == "O" and q == -1:
            nb = a.GetNeighbors()
            if not nb:
                continue
            n = nb[0]
            if n.GetSymbol() in ("P", "S"):
                s += 1  # phosphate / sulfate oxygen
            elif n.GetSymbol() == "C":
                carbonyl = any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() == "O"
                               for b in n.GetBonds())
                s += 1 if carbonyl else -2  # carboxylate good, alkoxide/phenolate bad
        elif sym == "N" and q == 1:
            amide = any(any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(nb).GetSymbol() == "O"
                            for b in nb.GetBonds()) for nb in a.GetNeighbors() if nb.GetSymbol() == "C")
            s += -2 if (amide or a.GetIsAromatic()) else 1  # ammonium good, acyl/aromatic N+ bad
    return s


def resolve(smiles: str, tgt_elems: dict, tgt_charge: int):
    """Return (new_smiles, status). status in resolved/ambiguous/unresolved.

    Enumerate protonation states of the skeleton, keep those whose element
    formula and net charge match the model, then pick the most chemically
    standard one by chem_score. resolved = a unique top score; ambiguous = a tie
    (deterministic pick, flagged); unresolved = no state reproduces the model
    formula/charge (usually a model-side error)."""
    try:
        variants = set(protonate_smiles(smiles, ph_min=1.0, ph_max=13.0, precision=1.0, max_variants=64))
    except Exception:
        return "", "unresolved"
    scored = []
    for v in variants:
        m = Chem.MolFromSmiles(v)
        if m is None:
            continue
        if Chem.GetFormalCharge(m) == tgt_charge and elems(rdMolDescriptors.CalcMolFormula(m)) == tgt_elems:
            scored.append((chem_score(m), Chem.MolToSmiles(m)))
    if not scored:
        return "", "unresolved"
    best = max(s for s, _ in scored)
    top = sorted(smi for s, smi in scored if s == best)
    if len(top) == 1:
        return top[0], "resolved"
    return top[0], "ambiguous"


def main() -> int:
    rows = [r for r in csv.DictReader(open(IDENTITY, encoding="utf-8"), delimiter="\t")
            if r["consistency"] == "protonation"]
    seen = {}
    for r in rows:
        seen.setdefault((r["stored_smiles"], r["model_formula"], r["model_charge"]), r)
    uniq = list(seen.values())
    print(f"distinct protonation structures: {len(uniq)}", flush=True)

    out_rows, tally = [], Counter()
    for i, r in enumerate(uniq):
        new_smiles, status = resolve(r["stored_smiles"], elems(r["model_formula"]), int(r["model_charge"]))
        new_inchi = ""
        if new_smiles:
            mol = Chem.MolFromSmiles(new_smiles)
            new_inchi = Chem.MolToInchi(mol) if mol else ""
        tally[status] += 1
        out_rows.append({
            "mets": r["mets"], "name": r["name"],
            "model_formula": r["model_formula"], "model_charge": r["model_charge"],
            "old_smiles": r["stored_smiles"], "new_smiles": new_smiles,
            "new_inchi": new_inchi, "status": status,
        })
        if (i + 1) % 200 == 0:
            print(f"  ...{i+1}/{len(uniq)}", file=sys.stderr, flush=True)

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"resolution: {dict(tally)}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
