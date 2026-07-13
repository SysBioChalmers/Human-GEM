"""Apply the model-invalid formula/charge corrections from charge_fix.py.

Only the cases whose current formula/charge is a physically impossible radical
(closed_shell == False) are applied: charge_fix (formula right, charge wrong),
formula_fix (charge right, H count wrong) and model_invalid_both. For each, the
metabolite's formula and/or charge in model/Human-GEM.yml is corrected and its
metSmiles/metInChI in metabolites.tsv is set to the matching structure. Because
these change mass/charge balance, run balanceTest.py before and after.

  python code/annotation/apply_charge.py            # dry run
  python code/annotation/apply_charge.py --write
"""

import csv
import re
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

YAML = Path("model/Human-GEM.yml")
MET_TSV = Path("model/metabolites.tsv")
PROP = Path("data/annotation/charge_proposals.tsv")
APPLY = {"charge_fix", "formula_fix", "model_invalid_both"}


def main(write: bool) -> int:
    props = [r for r in csv.DictReader(PROP.open(encoding="utf-8"), delimiter="\t") if r["status"] in APPLY]
    # key by metsNoComp + current formula/charge so every compartment is updated
    fix = {}
    for r in props:
        noc = r["mets"][:-1] if r["mets"][-1].isalpha() else r["mets"]
        new_smi = r["new_smiles"]
        inchi = Chem.MolToInchi(Chem.MolFromSmiles(new_smi)) if new_smi else ""
        fix[(noc, r["model_formula"], r["old_charge"])] = {
            "formula": r["new_formula"] or r["model_formula"],
            "charge": r["new_charge"] if r["new_charge"] else r["old_charge"],
            "smiles": new_smi, "inchi": inchi,
        }

    # --- metabolites.tsv: current formula/charge live in the YAML, so read them there
    tsv = list(csv.reader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t"))
    th = tsv[0]
    ti = {c: i for i, c in enumerate(th)}

    # parse YAML metabolite blocks
    lines = YAML.read_text(encoding="utf-8").splitlines()
    cur = {}
    blocks = []  # (start_idx, dict of field->line_idx, id, formula, charge)
    insec = False
    for i, ln in enumerate(lines):
        if ln.startswith("- ") and not ln.startswith("- !!omap"):
            insec = ln.strip() == "- metabolites:"
        if not insec:
            continue
        s = ln.strip()
        if s == "- !!omap":
            cur = {"_fields": {}}
            blocks.append(cur)
            continue
        m = re.match(r'- (\w+):\s*"?(.*?)"?\s*$', s)
        if m and cur:
            cur["_fields"][m.group(1)] = i
            cur[m.group(1)] = m.group(2)

    changed_mets = {}
    n_yaml = 0
    for b in blocks:
        mid = b.get("id", "")
        noc = mid[:-1] if mid and mid[-1].isalpha() else mid
        key = (noc, b.get("formula", ""), b.get("charge", ""))
        if key not in fix:
            continue
        f = fix[key]
        n_yaml += 1
        changed_mets[mid] = f
        if write:
            fi = b["_fields"]
            fl = lines[fi["formula"]]
            indent = fl[:len(fl) - len(fl.lstrip())]
            lines[fi["formula"]] = f'{indent}- formula: "{f["formula"]}"'
            lines[fi["charge"]] = f'{indent}- charge: {f["charge"]}'

    for r in tsv[1:]:
        f = changed_mets.get(r[ti["mets"]])
        if f:
            r[ti["metSmiles"]] = f["smiles"]
            r[ti["metInChI"]] = f["inchi"]

    print(f"distinct structures: {len(fix)}   metabolite entries updated: {n_yaml}")
    if write:
        YAML.write_text("\n".join(lines) + "\n", encoding="utf-8")
        with MET_TSV.open("w", encoding="utf-8", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(tsv)
        print("wrote model/Human-GEM.yml and model/metabolites.tsv")
    else:
        print("(dry run; pass --write to apply)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main("--write" in sys.argv))
