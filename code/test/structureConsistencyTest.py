"""QC test: metabolite formula/charge (model/Human-GEM.yml) must be consistent
with the SMILES/InChI stored in model/metabolites.tsv.

Formula and charge live in the YAML (they define the model's mass and charge
balance); the structures (SMILES, InChI) live in metabolites.tsv alongside the
other metabolite annotations. Because they are stored separately, they can drift
apart. This test parses each metabolite's SMILES with RDKit, derives its formula
and net charge, and compares them to the curated formula/charge. It also checks
that the stored SMILES and InChI describe the same structure.

The pH-7.3 microspecies convention means a metabolite's stored structure should
carry the same net charge as its `charge` field (acetate is C2H3O2 / -1, so its
SMILES must be the anion, not neutral acetic acid).

Categories per metabolite:
  ok             SMILES formula and charge match the model
  protonation    same heavy-atom skeleton, structure is a different protonation
                 state than the model (inconsistent: fix the SMILES/InChI or the charge)
  formula_error  heavy-atom composition disagrees (the wrong molecule is stored)
  smiles_inchi   SMILES and InChI disagree with each other
  generic        R-group / polymer, no concrete structure to check
  no_structure   no SMILES stored

Writes data/testResults/qc_structure_consistency.csv (the inconsistent ones) and
prints a summary. Report only: it does not fail the build. The workflow tracks
the counts with a delta versus the target branch, so a pull request that
introduces a new inconsistency is visible.

    python code/test/structureConsistencyTest.py
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, "code/annotation")
from identity import classify, load_model_metabolites  # noqa: E402

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MET_TSV = Path("model/metabolites.tsv")
OUT = Path("data/testResults/qc_structure_consistency.csv")
INCONSISTENT = {"protonation", "formula_error", "smiles_inchi"}


def load_structures() -> dict:
    """mets -> {metSmiles, metInChI} from metabolites.tsv."""
    with MET_TSV.open(encoding="utf-8", newline="") as fh:
        return {r["mets"]: r for r in csv.DictReader(fh, delimiter="\t")}


def inchikey_from_inchi(inchi: str) -> str:
    if not inchi:
        return ""
    try:
        return Chem.InchiToInchiKey(inchi)
    except Exception:
        return ""


def main() -> int:
    model = load_model_metabolites()
    struct = load_structures()
    rows, counts = [], {}
    for mid, m in model.items():
        s = struct.get(mid, {}) or {}
        smiles = s.get("metSmiles", "").strip()
        inchi = s.get("metInChI", "").strip()
        cat, rdf, rdc, _canon, ik_smiles = classify(m.get("formula", ""), m.get("charge", ""), smiles)
        # cross-check SMILES against InChI when both are present and concrete
        if cat == "ok" and inchi:
            ik_inchi = inchikey_from_inchi(inchi)
            if ik_inchi and ik_smiles and ik_inchi != ik_smiles:
                cat = "smiles_inchi"
        counts[cat] = counts.get(cat, 0) + 1
        if cat in INCONSISTENT:
            rows.append({
                "mets": mid,
                "name": m.get("name", ""),
                "issue": cat,
                "model_formula": m.get("formula", ""),
                "model_charge": m.get("charge", ""),
                "smiles_formula": rdf,
                "smiles_charge": rdc,
                "metSmiles": smiles,
            })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["mets", "name", "issue", "model_formula",
                                           "model_charge", "smiles_formula", "smiles_charge", "metSmiles"])
        w.writeheader()
        w.writerows(rows)
    total = sum(counts.values())
    n_bad = sum(counts.get(c, 0) for c in INCONSISTENT)
    print(f"metabolites: {total}")
    for cat in ("ok", "protonation", "formula_error", "smiles_inchi", "generic", "no_structure"):
        if cat in counts:
            print(f"  {cat:14} {counts[cat]:5}")
    print(f"\nformula/charge inconsistent with structure: {n_bad}; see {OUT}")
    if n_bad:
        print(f"::warning::{n_bad} metabolite(s) whose SMILES/InChI disagree with their formula/charge.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
