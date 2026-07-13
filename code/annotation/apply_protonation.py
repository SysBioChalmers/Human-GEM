"""Apply the resolved pH-7.3 structures to model/metabolites.tsv.

Reads data/annotation/protonation_resolved.tsv and, for every metabolite whose
stored metSmiles matches a resolved structure (same old SMILES, same model
formula/charge), replaces metSmiles and metInChI with the charged microspecies.
Only status == "resolved" entries are applied; ambiguous and unresolved are left
untouched for review. Report only unless run with --write.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, "code/annotation")
from identity import load_model_metabolites  # noqa: E402

MET_TSV = Path("model/metabolites.tsv")
RESOLVED = Path("data/annotation/protonation_resolved.tsv")


def main(write: bool) -> int:
    model = load_model_metabolites()
    resolved = {}
    for r in csv.DictReader(RESOLVED.open(encoding="utf-8"), delimiter="\t"):
        if r["status"] == "resolved":
            resolved[(r["old_smiles"], r["model_formula"], r["model_charge"])] = (r["new_smiles"], r["new_inchi"])

    with MET_TSV.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    i_mets, i_smi, i_inchi = header.index("mets"), header.index("metSmiles"), header.index("metInChI")

    n = 0
    for r in rows[1:]:
        mid = r[i_mets]
        m = model.get(mid)
        if not m:
            continue
        key = (r[i_smi], m.get("formula", ""), m.get("charge", ""))
        if key in resolved:
            new_smi, new_inchi = resolved[key]
            if new_smi and r[i_smi] != new_smi:
                r[i_smi] = new_smi
                r[i_inchi] = new_inchi
                n += 1

    print(f"resolved structures: {len(resolved)}  metabolite rows updated: {n}")
    if write:
        with MET_TSV.open("w", encoding="utf-8", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
        print(f"wrote {MET_TSV}")
    else:
        print("(dry run; pass --write to apply)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main("--write" in sys.argv))
