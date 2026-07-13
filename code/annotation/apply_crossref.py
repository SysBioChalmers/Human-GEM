"""Phase 3 apply: turn the structure-based cross-reference audit into concrete,
reviewable changes to model/metabolites.tsv.

Selects ONE best identifier per metabolite and namespace from the structure-
verified candidates (MetaNetX hub), and emits three change categories so they can
be committed separately:

  add      an empty cell gains the best verified id
  drift    metMetaNetXID: a deprecated MNXM is replaced by the current one
  correct  an existing id the structure contradicts (for review, not auto-applied)

Best-id selection: prefer identifiers carried by more of the matched MNXM
(consensus ~ primary), then a namespace format preference (modern HMDB, KEGG
compound over drug/glycan, MetaNetX id whose full InChIKey matches the model
exactly), then lexicographic.

Usage (dry run writes only the proposal file):
  python code/annotation/apply_crossref.py <scratch_with_mnx_tsvs>
  python code/annotation/apply_crossref.py <scratch> --write add
  python code/annotation/apply_crossref.py <scratch> --write drift
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MET_TSV = Path("model/metabolites.tsv")
PROPOSAL = Path("data/annotation/crossref_proposals.tsv")

COL2DB = {
    "metKEGGID": ("kegg.compound", "keggC", "kegg.drug", "keggD", "kegg.glycan", "keggG"),
    "metChEBIID": ("chebi", "CHEBI"),
    "metBiGGID": ("bigg.metabolite",),
    "metHMDBID": ("hmdb",),
    "metLipidMapsID": ("lipidmaps",),
    "metSeedID": ("seed.compound",),
}


def pragmatic(ik):
    return ik.rsplit("-", 1)[0] if ik else ""


def model_inchikey(smiles):
    if not smiles:
        return ""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return Chem.MolToInchiKey(m)
    except Exception:
        return ""


def fmt_rank(idv, col):
    if col == "metHMDBID":
        return 0 if (idv.startswith("HMDB") and len(idv) == 11) else 1  # HMDB0000000
    if col == "metKEGGID":
        return {"C": 0, "D": 1, "G": 2}.get(idv[:1], 3)
    return 0


def normalize(idv, col):
    if col == "metChEBIID":
        return idv if idv.startswith("CHEBI:") else "CHEBI:" + idv
    if col == "metHMDBID" and idv.startswith("HMDB") and len(idv) < 11:
        return "HMDB" + idv[4:].zfill(7)  # HMDB03450 -> HMDB0003450
    return idv


def best(candidates, col):
    """candidates: iterable of ids (with repeats across MNXM). Return best one."""
    c = Counter(normalize(x, col) for x in candidates)
    if not c:
        return ""
    return sorted(c, key=lambda i: (-c[i], fmt_rank(i, col), i))[0]


def load_mnx(scratch):
    key2mnx = defaultdict(set)
    mnx2key, mnx2full = {}, {}
    for line in (scratch / "mnx_struct.tsv").open(encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1]:
            key2mnx[pragmatic(p[1])].add(p[0])
            mnx2key[p[0]] = pragmatic(p[1])
            mnx2full[p[0]] = p[1]
    mnx2xref = defaultdict(lambda: defaultdict(list))
    for line in (scratch / "mnx_xref.tsv").open(encoding="utf-8"):
        src, _, mnxm = line.rstrip("\n").partition("\t")
        db, _, extid = src.partition(":")
        if mnxm.strip() and db and extid:
            mnx2xref[mnxm.strip()][db].append(extid)
    return key2mnx, mnx2key, mnx2full, mnx2xref


def build_proposals(scratch):
    key2mnx, mnx2key, mnx2full, mnx2xref = load_mnx(scratch)
    rows = list(csv.DictReader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t"))
    props = []  # (mets, column, action, old, new)
    for r in rows:
        ik = model_inchikey(r.get("metSmiles", "").strip())
        if not ik:
            continue
        pk = pragmatic(ik)
        mnxset = key2mnx.get(pk, set())
        if not mnxset:
            continue
        # candidates per db (with multiplicity across matched MNXM)
        cand = defaultdict(list)
        for mnxm in mnxset:
            for db, ids in mnx2xref.get(mnxm, {}).items():
                cand[db].extend(ids)
        # metMetaNetXID drift: best current MNXM (prefer exact full-InChIKey match)
        exact = sorted(m for m in mnxset if mnx2full.get(m) == ik)
        best_mnx = exact[0] if exact else sorted(mnxset, key=lambda m: (-len(mnx2xref.get(m, {})), m))[0]
        have_mnx = {x.strip() for x in (r.get("metMetaNetXID") or "").split(";") if x.strip()}
        if have_mnx and not any(mnx2key.get(h) == pk for h in have_mnx):
            props.append((r["mets"], "metMetaNetXID", "drift", ";".join(sorted(have_mnx)), best_mnx))
        elif not have_mnx:
            props.append((r["mets"], "metMetaNetXID", "add", "", best_mnx))
        # other namespaces: best id into empty cells (add); flag contradicted (correct)
        for col, dbs in COL2DB.items():
            pool = []
            for db in dbs:
                pool.extend(cand.get(db, []))
            if not pool:
                continue
            pick = best(pool, col)
            have = {x.strip() for x in (r.get(col) or "").split(";") if x.strip()}
            vset = {normalize(x, col) for x in pool}
            if not have:
                props.append((r["mets"], col, "add", "", pick))
            elif not (have & vset):
                props.append((r["mets"], col, "correct", ";".join(sorted(have)), pick))
    return props


def main():
    scratch = Path(sys.argv[1])
    mode = sys.argv[sys.argv.index("--write") + 1] if "--write" in sys.argv else None
    cols = None
    if "--cols" in sys.argv:
        cols = set(sys.argv[sys.argv.index("--cols") + 1].split(","))
    props = build_proposals(scratch)
    PROPOSAL.parent.mkdir(parents=True, exist_ok=True)
    with PROPOSAL.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["mets", "column", "action", "old", "new"])
        w.writerows(props)
    by = Counter((p[2], p[1]) for p in props)
    print(f"proposals: {len(props)}  wrote {PROPOSAL}")
    for (action, col), n in sorted(by.items()):
        print(f"  {action:8} {col:16} {n}")
    if mode:
        apply_mode(props, mode, cols)


def apply_mode(props, mode, cols=None):
    rows = list(csv.reader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t"))
    header = rows[0]
    by_met = defaultdict(list)
    touched = set()
    for m, col, action, old, new in props:
        if action == mode and new and (cols is None or col in cols):
            by_met[m].append((col, new))
            touched.add(col)
    # create a new column (e.g. metSeedID for ModelSEED) only if it is being written
    for col in sorted(touched):
        if col not in header:
            header.append(col)
            for r in rows[1:]:
                r.append("")
    idx = {c: i for i, c in enumerate(header)}
    ri = {r[0]: r for r in rows[1:]}
    n = 0
    for m, changes in by_met.items():
        r = ri.get(m)
        if not r:
            continue
        while len(r) < len(header):
            r.append("")
        for col, new in changes:
            r[idx[col]] = new
            n += 1
    with MET_TSV.open("w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"applied mode={mode}: {n} cells in {len(by_met)} metabolites")


if __name__ == "__main__":
    main()
