"""Build the canonical metabolite identity table and audit structure vs
formula/charge consistency (Phase 0-1 of the annotation overhaul).

For every metabolite it gathers:
  - model name, formula, charge (from model/Human-GEM.yml)
  - stored SMILES / InChI / InChIKey (from data/modelCuration/metabolites_SMILES_Inchi.tsv)
  - RDKit-derived canonical SMILES, Hill formula, net charge, InChI, InChIKey

then classifies each metabolite by how the stored structure relates to the
curated formula/charge:

  ok            structure formula and charge match the model
  protonation   same heavy-atom skeleton, differ only in H count and charge
                (structure is a different protonation state than the model)
  generic       R-group / polymer: model formula contains R, or SMILES has a
                dummy atom (*). Cannot be verified against a concrete structure.
  formula_error heavy-atom composition differs and it is not generic (a real
                discrepancy: the stored structure is the wrong molecule)
  no_structure  no SMILES stored
  unparseable   RDKit cannot parse the SMILES

Run from a Human-GEM checkout:
    python code/annotation/identity.py
Writes data/annotation/metabolite_identity.tsv and prints a summary.
"""

import csv
import re
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

YAML = Path("model/Human-GEM.yml")
STRUCT = Path("data/modelCuration/metabolites_SMILES_Inchi.tsv")
OUT = Path("data/annotation/metabolite_identity.tsv")

_ELEM = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> Counter:
    """Element -> count for a Hill formula string (charge sign ignored)."""
    c = Counter()
    if not formula:
        return c
    for sym, num in _ELEM.findall(formula.strip().rstrip("+-")):
        if sym:
            c[sym] += int(num) if num else 1
    return c


def load_model_metabolites() -> dict:
    """id -> {name, formula, charge} from the YAML metabolites section."""
    mets, cur, insec = {}, {}, False
    for line in YAML.open(encoding="utf-8"):
        if line.startswith("- ") and not line.startswith("- !!omap"):
            insec = line.strip() == "- metabolites:"
            continue
        if not insec:
            continue
        s = line.strip()
        if s == "- !!omap":
            if cur.get("id"):
                mets[cur["id"]] = cur
            cur = {}
            continue
        m = re.match(r'- (\w+):\s*"?(.*?)"?\s*$', s)
        if m:
            cur[m.group(1)] = m.group(2)
    if cur.get("id"):
        mets[cur["id"]] = cur
    return mets


def load_structures() -> dict:
    """id -> {SMILES, inchi, inchikey}."""
    out = {}
    with STRUCT.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["mets"]] = r
    return out


def classify(model_formula, model_charge, smiles):
    """Return (category, rdkit_formula, rdkit_charge, canon_smiles, inchikey)."""
    if not smiles:
        return "no_structure", "", "", "", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "unparseable", "", "", "", ""
    # generic representation: an R group (dummy atom) in the SMILES or an "R" in
    # the model formula means there is no concrete structure to check against.
    generic = "*" in smiles or "R" in (model_formula or "")
    rd_formula = rdMolDescriptors.CalcMolFormula(mol)  # e.g. C2H3O2-
    rd_charge = Chem.GetFormalCharge(mol)
    canon = Chem.MolToSmiles(mol)
    try:
        ik = Chem.MolToInchiKey(mol)
    except Exception:
        ik = ""
    m_el = parse_formula(model_formula)
    r_el = parse_formula(rd_formula)
    try:
        m_charge = int(model_charge)
    except (TypeError, ValueError):
        m_charge = None
    # compare heavy-atom (non-H) composition
    m_heavy = {k: v for k, v in m_el.items() if k != "H"}
    r_heavy = {k: v for k, v in r_el.items() if k != "H"}
    if generic:
        cat = "generic"
    elif m_heavy != r_heavy:
        cat = "formula_error"
    elif m_el == r_el and m_charge == rd_charge:
        cat = "ok"
    else:
        # same heavy atoms; differ in H and/or charge -> protonation variant
        cat = "protonation"
    return cat, rd_formula, str(rd_charge), canon, ik


def main() -> int:
    model = load_model_metabolites()
    struct = load_structures()
    rows, cats = [], Counter()
    for mid, m in model.items():
        smi = (struct.get(mid, {}) or {}).get("SMILES", "").strip()
        cat, rdf, rdc, canon, ik = classify(m.get("formula", ""), m.get("charge", ""), smi)
        cats[cat] += 1
        rows.append({
            "mets": mid,
            "name": m.get("name", ""),
            "model_formula": m.get("formula", ""),
            "model_charge": m.get("charge", ""),
            "stored_smiles": smi,
            "rdkit_formula": rdf,
            "rdkit_charge": rdc,
            "rdkit_inchikey": ik,
            "canonical_smiles": canon,
            "consistency": cat,
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    total = len(rows)
    print(f"metabolites: {total}")
    for cat, n in cats.most_common():
        print(f"  {cat:14} {n:5}  ({100*n/total:.1f}%)")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
