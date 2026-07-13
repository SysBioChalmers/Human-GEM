"""Phase 4 apply (#967): write structure-verified reaction cross-references.

Reuses reaction_verify to match each reaction to MetaNetX by structural signature,
then selects one best identifier per namespace and emits three categories, applied
as separate commits like the metabolite cross-references:

  add      empty cell gains the best verified id
  drift    rxnMetaNetXID: a deprecated MNXR replaced by the current one
  correct  an existing id the structure contradicts (review, not auto-applied)

  python code/annotation/apply_reaction.py <scratch>
  python code/annotation/apply_reaction.py <scratch> --write add
  python code/annotation/apply_reaction.py <scratch> --write drift
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, "code/annotation")
import reaction_verify as rv  # noqa: E402

RXN_TSV = Path("model/reactions.tsv")
PROPOSAL = Path("data/annotation/reaction_proposals.tsv")

COL2DB = {
    "rxnKEGGID": ("kegg.reaction",),
    "rxnBiGGID": ("bigg.reaction",),
    "rxnRheaID": ("rhea",),
}


def best(candidates, col):
    c = Counter(candidates)
    if not c:
        return ""
    if col == "rxnRheaID":                       # Rhea quartet -> master (lowest)
        return "RHEA:" + min(candidates, key=lambda x: int(x) if x.isdigit() else 1 << 30)
    if col == "rxnBiGGID":
        def rank(i):
            return (0 if not i.startswith("R_") else 1, -c[i], len(i), i)
        return sorted(c, key=rank)[0]
    return sorted(c, key=lambda i: (-c[i], i))[0]


def build_proposals(scratch):
    rv.SCRATCH = Path(scratch)
    mnx2key = rv.load_mnx_keys()
    sig2mnxr, mnxr2xref = rv.build_mnxr_index(mnx2key)
    metkey = rv.model_met_keys()
    rxns = rv.parse_model_reactions()
    rrows = {r["rxns"]: r for r in csv.DictReader(RXN_TSV.open(encoding="utf-8", newline=""), delimiter="\t")}

    props = []
    for rid, r in rxns.items():
        sub, prod, ok = defaultdict(float), defaultdict(float), True
        for met, coeff in r["mets"].items():
            noc = met[:-1] if met and met[-1].isalpha() else met
            if noc in ("MAM02039", "MAM02040"):
                continue
            k = metkey.get(met)
            if not k:
                ok = False
                break
            (sub if coeff < 0 else prod)[k] += abs(coeff)
        if not ok:
            continue
        sig = rv.signature(dict(sub), dict(prod))
        mnxrset = sig2mnxr.get(sig) if sig else None
        if not mnxrset:
            continue
        cand = defaultdict(list)
        for mnxr in mnxrset:
            for db, ids in mnxr2xref.get(mnxr, {}).items():
                cand[db].extend(ids)
        rr = rrows.get(rid, {})
        # MetaNetX drift / add
        have_mnx = rv.norm(rr.get("rxnMetaNetXID", ""))
        best_mnxr = sorted(mnxrset, key=lambda m: (-len(mnxr2xref.get(m, {})), m))[0]
        if have_mnx and not (have_mnx & mnxrset):
            props.append((rid, "rxnMetaNetXID", "drift", ";".join(sorted(have_mnx)), best_mnxr))
        elif not have_mnx:
            props.append((rid, "rxnMetaNetXID", "add", "", best_mnxr))
        # other namespaces
        for col, dbs in COL2DB.items():
            pool = []
            for db in dbs:
                pool.extend(cand.get(db, []))
            if not pool:
                continue
            pick = best(pool, col)
            have = rv.norm(rr.get(col, ""))
            vset = {("RHEA:" + x if col == "rxnRheaID" and not x.startswith("RHEA:") else x) for x in pool}
            if col == "rxnRheaID":
                have = {("RHEA:" + x if not x.startswith("RHEA:") else x) for x in have}
            if not have:
                props.append((rid, col, "add", "", pick))
            elif not (have & vset):
                props.append((rid, col, "correct", ";".join(sorted(have)), pick))
    return props


def apply_mode(props, mode):
    rows = list(csv.reader(RXN_TSV.open(encoding="utf-8", newline=""), delimiter="\t"))
    header = rows[0]
    idx = {c: i for i, c in enumerate(header)}
    by = defaultdict(list)
    for rid, col, action, old, new in props:
        if action == mode and new:
            by[rid].append((col, new))
    ri = {r[0]: r for r in rows[1:]}
    n = 0
    for rid, changes in by.items():
        r = ri.get(rid)
        if not r:
            continue
        while len(r) < len(header):
            r.append("")
        for col, new in changes:
            r[idx[col]] = new
            n += 1
    with RXN_TSV.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"applied mode={mode}: {n} cells in {len(by)} reactions")


def main():
    scratch = sys.argv[1]
    props = build_proposals(scratch)
    PROPOSAL.parent.mkdir(parents=True, exist_ok=True)
    with PROPOSAL.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rxns", "column", "action", "old", "new"])
        w.writerows(props)
    for (a, c), n in sorted(Counter((p[2], p[1]) for p in props).items()):
        print(f"  {a:8} {c:16} {n}")
    if "--write" in sys.argv:
        apply_mode(props, sys.argv[sys.argv.index("--write") + 1])


if __name__ == "__main__":
    main()
