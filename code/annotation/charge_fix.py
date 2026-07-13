"""Find metabolites whose curated charge is inconsistent with their structure and
propose the corrected charge.

A metabolite whose stored structure has the right skeleton but whose formula/charge
cannot be reproduced by any protonation state of that skeleton has a model-side
error. The usual case is an anion formula carried at charge 0 (a carboxylate stored
as C..H..O.. / 0), which is a physically impossible (half-integer DBE) neutral.

For each such metabolite this finds the protonation state whose element formula
matches the model formula exactly and reads off its charge: that is the charge the
formula actually implies. It proposes new_charge and the matching structure.

Writes data/annotation/charge_proposals.tsv (mets, name, model_formula,
old_charge, new_charge, new_smiles). Report only.
"""

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, "code/annotation")
from identity import classify, load_model_metabolites, parse_formula  # noqa: E402

from dimorphite_dl import protonate_smiles  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

MET_TSV = Path("model/metabolites.tsv")
OUT = Path("data/annotation/charge_proposals.tsv")


def elems(formula):
    return parse_formula(formula)


# atomic numbers for the elements that occur in the model
_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "P": 15, "S": 16,
      "Cl": 17, "K": 19, "Ca": 20, "Mn": 25, "Fe": 26, "Co": 27, "Cu": 29, "Zn": 30,
      "Se": 34, "Br": 35, "Mo": 42, "I": 53, "As": 33, "R": 0, "X": 0}


def closed_shell(el, charge):
    """A species is closed-shell iff its total electron count (sum Z - charge) is
    even; an odd count is a radical, i.e. an impossible formula/charge combination."""
    try:
        total = sum(_Z.get(k, None) * v for k, v in el.items())
    except TypeError:
        return None  # unknown element, cannot judge
    return (total - charge) % 2 == 0


def physiological_state(smiles):
    """The structure's dominant pH-7.3 microspecies as (smiles, elems, charge)."""
    try:
        out = list(protonate_smiles(smiles, ph_min=7.3, ph_max=7.3, max_variants=1))
    except Exception:
        return None
    if not out:
        return None
    m = Chem.MolFromSmiles(out[0])
    if m is None:
        return None
    return Chem.MolToSmiles(m), elems(rdMolDescriptors.CalcMolFormula(m)), Chem.GetFormalCharge(m)


# names that denote genuine odd-electron (radical) species; their odd electron
# count is intentional, not an error, so the closed-shell test must not touch them
_RADICAL = re.compile(r"radical|semiquinone|monodehydro|semidehydro|\bnitroxide\b", re.I)


def is_radical(name):
    return bool(_RADICAL.search(name or ""))


def resolve_case(smiles, model_elems, model_charge):
    """Decide whether a formula/charge inconsistency is a charge error or a formula
    error by comparing to the physiological state. Returns (status, new_charge,
    new_formula, new_smiles)."""
    phys = physiological_state(smiles)
    if phys is None:
        return "no_phys", None, None, ""
    p_smiles, p_elems, p_charge = phys
    model_ok = closed_shell(model_elems, model_charge)
    # Only override curated formula/charge when it is provably invalid (a radical).
    # If the model is a valid closed-shell species, keep it and fix the structure.
    if model_ok is False:
        if p_elems == model_elems:
            return "charge_fix", p_charge, None, p_smiles      # formula right, charge wrong
        if p_charge == model_charge:
            return "formula_fix", None, _hill(p_elems), p_smiles  # charge right, H count wrong
        return "model_invalid_both", p_charge, _hill(p_elems), p_smiles
    if (p_elems, p_charge) == (model_elems, model_charge):
        return "consistent_phys", None, None, p_smiles
    # model is a valid species that differs from the physiological guess: trust the
    # curated formula/charge, the stored structure is what needs correcting
    return "structure_fix", None, None, ""


def _hill(el):
    order = ["C", "H"] + sorted(k for k in el if k not in ("C", "H"))
    return "".join(f"{k}{el[k] if el[k] > 1 else ''}" for k in order if el.get(k))


def main() -> int:
    model = load_model_metabolites()
    struct = {r["mets"]: r for r in csv.DictReader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t")}
    seen, props = {}, []
    for mid, m in model.items():
        smi = (struct.get(mid, {}) or {}).get("metSmiles", "").strip()
        cat = classify(m.get("formula", ""), m.get("charge", ""), smi)[0]
        if cat != "protonation":
            continue
        key = (smi, m.get("formula", ""), m.get("charge", ""))
        if key in seen:
            continue
        seen[key] = 1
        try:
            old_charge = int(m.get("charge", ""))
        except ValueError:
            old_charge = None
        if is_radical(m.get("name", "")):
            props.append((mid, m.get("name", ""), m.get("formula", ""), str(old_charge),
                          "", "", "", "radical_skip"))
            continue
        status, new_charge, new_formula, new_smiles = resolve_case(
            smi, elems(m.get("formula", "")), old_charge)
        props.append((mid, m.get("name", ""), m.get("formula", ""), str(old_charge),
                      "" if new_charge is None else str(new_charge),
                      new_formula or "", new_smiles, status))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["mets", "name", "model_formula", "old_charge", "new_charge",
                    "new_formula", "new_smiles", "status"])
        w.writerows(props)
    from collections import Counter
    tally = Counter(p[-1] for p in props)
    print(f"distinct inconsistent structures: {len(props)}  {dict(tally)}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
