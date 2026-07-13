"""Phase 3: verify and extend metabolite cross-references against structure.

Uses MetaNetX as the hub. chem_prop maps MNXM -> InChIKey; chem_xref maps MNXM ->
every other database identifier (KEGG, ChEBI, BiGG, HMDB, LipidMaps, ModelSEED).
Each model metabolite's InChIKey is computed from its metSmiles and matched to
MNXM on the pragmatic key (the InChIKey without its final protonation character,
so a charged model species matches a neutral database entry). From the matched
MNXM we read the structure-verified identifiers for every namespace and compare
them to what metabolites.tsv currently holds:

  confirmed  existing id is among the structure-verified ids
  wrong      existing id is NOT structure-verified (candidate to correct/deprecate)
  new        a structure-verified id is missing from the model

Inputs (produced by streaming MetaNetX, see scratchpad extract):
  mnx_struct.tsv : MNXM  InChIKey  formula  charge
  mnx_xref.tsv   : source(db:extid)  MNXM

Report only; writes data/annotation/crossref_audit.tsv.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MET_TSV = Path("model/metabolites.tsv")
SCRATCH = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
MNX_STRUCT = SCRATCH / "mnx_struct.tsv"
MNX_XREF = SCRATCH / "mnx_xref.tsv"
OUT = Path("data/annotation/crossref_audit.tsv")

# model column  ->  MetaNetX chem_xref source prefix(es). Both the identifiers.org
# style (kegg.compound) and the legacy style (keggC) are included where they carry
# the same extid format, for maximum coverage.
COL2DB = {
    "metKEGGID": ("kegg.compound", "keggC", "kegg.drug", "keggD", "kegg.glycan", "keggG"),
    "metChEBIID": ("chebi", "CHEBI"),
    "metBiGGID": ("bigg.metabolite",),
    "metHMDBID": ("hmdb",),
    "metLipidMapsID": ("lipidmaps",),
    "metMetaNetXID": ("metanetx.chemical",),  # handled specially: the MNXM itself
    "metSeedID": ("seed.compound",),  # ModelSEED (#968); column may not exist yet
}


def pragmatic(inchikey: str) -> str:
    """InChIKey without the final protonation character (protonation-invariant)."""
    return inchikey.rsplit("-", 1)[0] if inchikey else ""


def model_inchikey(smiles: str) -> str:
    if not smiles:
        return ""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return Chem.MolToInchiKey(m)
    except Exception:
        return ""


def load_mnx():
    """pragmatic-key -> set(MNXM); MNXM -> pragmatic key; MNXM -> {db: set(extid)}."""
    key2mnx = defaultdict(set)
    mnx2key = {}
    with MNX_STRUCT.open(encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                k = pragmatic(p[1])
                key2mnx[k].add(p[0])
                mnx2key[p[0]] = k
    mnx2xref = defaultdict(lambda: defaultdict(set))
    with MNX_XREF.open(encoding="utf-8") as fh:
        for line in fh:
            src, _, mnxm = line.rstrip("\n").partition("\t")
            mnxm = mnxm.strip()
            db, _, extid = src.partition(":")
            if mnxm and db and extid:
                mnx2xref[mnxm][db].add(extid)
    return key2mnx, mnx2key, mnx2xref


def norm(val: str) -> set:
    return {x.strip() for x in (val or "").split(";") if x.strip()}


def main() -> int:
    key2mnx, mnx2key, mnx2xref = load_mnx()
    print(f"MetaNetX: {len(key2mnx)} distinct structures, {len(mnx2xref)} cross-referenced MNXM")

    with MET_TSV.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    from collections import Counter
    tally = Counter()
    mets_gaining = defaultdict(set)   # col -> set of mets that would gain an id
    audit = []
    matched = 0
    for r in rows:
        ik = model_inchikey(r.get("metSmiles", "").strip())
        if not ik:
            continue
        pk = pragmatic(ik)
        mnxset = key2mnx.get(pk, set())
        if not mnxset:
            tally["no_mnx_match"] += 1
            continue
        matched += 1
        verified = defaultdict(set)
        for mnxm in mnxset:
            for db, ids in mnx2xref.get(mnxm, {}).items():
                verified[db].update(ids)
        for col, dbs in COL2DB.items():
            have = norm(r.get(col, ""))
            if col == "metMetaNetXID":
                # verify each existing MNXM by its own structure (handles ID drift:
                # an older MNXM that still resolves to this structure is confirmed).
                vset = {m for m in mnxset}
                for hid in have:
                    ok = mnx2key.get(hid) == pk
                    tally[f"{col}:{'confirmed' if ok else 'wrong'}"] += 1
                    if not ok:
                        audit.append((r["mets"], col, "wrong", hid, ";".join(sorted(vset)[:5])))
                if not have and vset:
                    mets_gaining[col].add(r["mets"])
                    audit.append((r["mets"], col, "new", "", ";".join(sorted(vset)[:5])))
                continue
            vset = set()
            for db in dbs:
                vset |= verified.get(db, set())
            if col == "metChEBIID":
                vset = {("CHEBI:" + x) if not x.startswith("CHEBI:") else x for x in vset}
            for hid in have:
                if hid in vset:
                    status = "confirmed"
                elif vset:                       # MetaNetX links this structure to
                    status = "wrong"             # other ids of this db, but not this one
                else:
                    status = "unverified"        # no id of this db for this structure
                tally[f"{col}:{status}"] += 1
                if status in ("wrong", "unverified"):
                    audit.append((r["mets"], col, status, hid, ";".join(sorted(vset)[:8])))
            if not have and vset:
                mets_gaining[col].add(r["mets"])
                audit.append((r["mets"], col, "new", "", ";".join(sorted(vset)[:8])))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["mets", "column", "status", "existing_id", "structure_verified"])
        w.writerows(audit)

    print(f"metabolites structure-matched to MetaNetX: {matched}  (no match: {tally['no_mnx_match']})")
    print(f"{'namespace':16} {'confirmed':>9} {'wrong':>6} {'unverif':>8} {'mets_gaining':>12}")
    for col in COL2DB:
        print(f"  {col:14} {tally[col+':confirmed']:>9} {tally[col+':wrong']:>6} "
              f"{tally[col+':unverified']:>8} {len(mets_gaining[col]):>12}")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
