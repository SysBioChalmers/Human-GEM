"""Harmonize metabolite cross-references across compartments and drop malformed
MetaNetX ids.

A metabolite has one structure, so every compartment must carry the same
cross-references. Applying cross-references per compartment let different
compartments (whose stored structures differ slightly) receive different ids,
which annotationTest flags. This collapses each metabolite's cross-references to a
single value per column: when the target branch already had one consistent value
it is kept (preserving curation); otherwise the most common current value is used.
metMetaNetXID values that are not of the form MNXM<digits> (e.g. the named id
WATER) are dropped first.
"""

import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

MET_TSV = Path("model/metabolites.tsv")
COLS = ["metBiGGID", "metKEGGID", "metHMDBID", "metChEBIID", "metPubChemID",
        "metLipidMapsID", "metMetaNetXID", "metSeedID"]
MNXM_OK = re.compile(r"^MNXM\d+$")


def noc(mid):
    return mid[:-1] if mid and mid[-1].isalpha() else mid


def clean_metanetx(val):
    parts = [p for p in val.split(";") if MNXM_OK.match(p)]
    return ";".join(dict.fromkeys(parts))


def develop_values():
    """metsNoComp -> {col -> set(distinct non-empty values on develop)}."""
    try:
        text = subprocess.run(["git", "show", "origin/develop:model/metabolites.tsv"],
                              capture_output=True, text=True, check=True).stdout
    except Exception:
        return {}
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    h = rows[0]
    idx = {c: i for i, c in enumerate(h)}
    dev = defaultdict(lambda: defaultdict(set))
    for r in rows[1:]:
        for col in COLS:
            if col in idx and idx[col] < len(r) and r[idx[col]]:
                dev[noc(r[0])][col].add(r[idx[col]])
    return dev


def main() -> int:
    dev = develop_values()
    rows = list(csv.reader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t"))
    h = rows[0]
    idx = {c: i for i, c in enumerate(h)}
    groups = defaultdict(list)
    for r in rows[1:]:
        groups[noc(r[0])].append(r)

    n_malformed = n_harmonized = 0
    for nc, grp in groups.items():
        for col in COLS:
            if col not in idx:
                continue
            i = idx[col]
            # drop malformed MetaNetX ids first
            if col == "metMetaNetXID":
                for r in grp:
                    cleaned = clean_metanetx(r[i])
                    if cleaned != r[i]:
                        n_malformed += 1
                        r[i] = cleaned
            vals = [r[i] for r in grp if r[i]]
            if len(set(vals)) <= 1:
                continue  # already consistent (or all empty)
            dev_vals = dev.get(nc, {}).get(col, set())
            consensus = sorted(Counter(vals).items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            # metMetaNetXID was rewritten wholesale by the drift update, so always
            # unify it. For the other columns only fix inconsistencies this PR
            # introduced (develop was consistent or empty); leave develop's own
            # cross-compartment inconsistencies alone.
            if col == "metMetaNetXID":
                canon = consensus
            elif len(dev_vals) == 1:
                canon = next(iter(dev_vals))          # keep the curated value
            elif not dev_vals:
                canon = consensus
            else:
                continue                              # pre-existing on develop, not ours
            for r in grp:
                r[i] = canon
            n_harmonized += 1

    with MET_TSV.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"malformed MetaNetX ids cleaned: {n_malformed}")
    print(f"metabolite/column pairs harmonized across compartments: {n_harmonized}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
