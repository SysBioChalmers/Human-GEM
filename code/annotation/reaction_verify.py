"""Phase 4 (#967): verify and extend reaction cross-references against structure.

Each reaction is reduced to a structural signature: the multiset of participant
InChIKeys (pragmatic key, protonation-invariant) with their coefficients on each
side, with protons and water removed and the two sides held direction-invariantly.
MetaNetX reac_prop gives the same signature for every MNXR (its equation is in
MNXM, which chem_prop maps to InChIKeys), and reac_xref maps MNXR to KEGG / BiGG /
Rhea / MetaNetX reaction ids. Matching the model's signature to MNXR therefore
verifies the reaction's cross-references and proposes structure-verified new ones.

Inputs (scratch): mnx_struct.tsv, mnxr_eq.tsv (MNXR eq), mnxr_xref.tsv (src MNXR).
Report only; writes data/annotation/reaction_audit.tsv.
"""

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

YAML = Path("model/Human-GEM.yml")
MET_TSV = Path("model/metabolites.tsv")
OUT = Path("data/annotation/reaction_audit.tsv")
SCRATCH = Path(sys.argv[1] if len(sys.argv) > 1 else ".")

SKIP_MNX = {"MNXM01", "MNXM1", "WATER", "MNXM2"}      # H+, H2O
COL2DB = {
    "rxnKEGGID": ("kegg.reaction",),
    "rxnBiGGID": ("bigg.reaction",),
    "rxnMetaNetXID": ("metanetx.reaction",),
    "rxnRheaID": ("rhea",),
}


def pragmatic(ik):
    return ik.rsplit("-", 1)[0] if ik else ""


def signature(sub: dict, prod: dict):
    """direction-invariant structural signature, or None if degenerate/empty."""
    s = frozenset(sub.items())
    p = frozenset(prod.items())
    if not s or not p or s == p:
        return None
    return frozenset((s, p))


def load_mnx_keys():
    mnx2key = {}
    for line in (SCRATCH / "mnx_struct.tsv").open(encoding="utf-8"):
        a = line.rstrip("\n").split("\t")
        if len(a) >= 2 and a[1]:
            mnx2key[a[0]] = pragmatic(a[1])
    return mnx2key


_TOKEN = re.compile(r"([\d.]+)\s+(\S+?)@\S+")


def parse_eq(eq, mnx2key):
    """MetaNetX equation -> (sub dict, prod dict) of key->coeff, or None if any
    non-skip participant lacks a structure key."""
    if "=" not in eq:
        return None
    left, _, right = eq.partition(" = ")
    out = []
    for side in (left, right):
        d = defaultdict(float)
        for coeff, mnxm in _TOKEN.findall(side):
            if mnxm in SKIP_MNX:
                continue
            k = mnx2key.get(mnxm)
            if not k:
                return None            # unmapped participant -> cannot match
            d[k] += float(coeff)
        out.append(dict(d))
    return out[0], out[1]


def build_mnxr_index(mnx2key):
    sig2mnxr = defaultdict(set)
    for line in (SCRATCH / "mnxr_eq.tsv").open(encoding="utf-8"):
        mnxr, _, eq = line.rstrip("\n").partition("\t")
        parsed = parse_eq(eq, mnx2key)
        if not parsed:
            continue
        sig = signature(*parsed)
        if sig:
            sig2mnxr[sig].add(mnxr)
    mnxr2xref = defaultdict(lambda: defaultdict(list))
    for line in (SCRATCH / "mnxr_xref.tsv").open(encoding="utf-8"):
        src, _, mnxr = line.rstrip("\n").partition("\t")
        db, _, extid = src.partition(":")
        if mnxr.strip() and db and extid:
            mnxr2xref[mnxr.strip()][db].append(extid)
    return sig2mnxr, mnxr2xref


def model_met_keys():
    keys = {}
    hplus_water = set()
    for r in csv.DictReader(MET_TSV.open(encoding="utf-8", newline=""), delimiter="\t"):
        smi = r.get("metSmiles", "").strip()
        m = Chem.MolFromSmiles(smi) if smi else None
        if m is not None:
            try:
                keys[r["mets"]] = pragmatic(Chem.MolToInchiKey(m))
            except Exception:
                pass
    return keys


def parse_model_reactions():
    """id -> {met: coeff} from the YAML reactions section."""
    rxns, cur, insec, inmet = {}, None, False, False
    for ln in YAML.read_text(encoding="utf-8").splitlines():
        if ln.startswith("- ") and not ln.startswith("- !!omap"):
            insec = ln.strip() == "- reactions:"
        if not insec:
            continue
        s = ln.strip()
        if s == "- !!omap":
            cur = {"mets": {}}
            inmet = False
            continue
        if cur is None:
            continue
        mid = re.match(r'- id:\s*"?(.*?)"?\s*$', s)
        if mid:
            cur["id"] = mid.group(1)
            rxns[cur["id"]] = cur
            continue
        if s.startswith("- metabolites:"):
            inmet = True
            continue
        if inmet:
            mm = re.match(r'- (MAM\w+):\s*(-?[\d.]+)\s*$', s)
            if mm:
                cur["mets"][mm.group(1)] = float(mm.group(2))
            elif s.startswith("- "):
                inmet = False
    return rxns


def norm(v):
    return {x.strip().replace("RHEA:", "") for x in (v or "").split(";") if x.strip()}


def main() -> int:
    mnx2key = load_mnx_keys()
    sig2mnxr, mnxr2xref = build_mnxr_index(mnx2key)
    print(f"MetaNetX reactions with a structural signature: {len(sig2mnxr)}")
    metkey = model_met_keys()
    rxns = parse_model_reactions()
    print(f"model reactions: {len(rxns)}  metabolites with a key: {len(metkey)}")

    tally = Counter()
    rxn_gain = defaultdict(set)
    audit = []
    matched = 0
    for rid, r in rxns.items():
        sub, prod = defaultdict(float), defaultdict(float)
        ok = True
        for met, coeff in r["mets"].items():
            noc = met[:-1] if met and met[-1].isalpha() else met
            if noc in ("MAM02039", "MAM02040"):   # H+, H2O
                continue
            k = metkey.get(met)
            if not k:
                ok = False
                break
            (sub if coeff < 0 else prod)[k] += abs(coeff)
        if not ok:
            tally["no_structure"] += 1
            continue
        sig = signature(dict(sub), dict(prod))
        if not sig:
            tally["degenerate"] += 1
            continue
        mnxrset = sig2mnxr.get(sig)
        if not mnxrset:
            tally["no_match"] += 1
            continue
        matched += 1
        verified = defaultdict(set)
        for mnxr in mnxrset:
            verified["metanetx.reaction"].add(mnxr)
            for db, ids in mnxr2xref.get(mnxr, {}).items():
                verified[db].update(ids)
        rrow = rrows.get(rid, {})
        for col, dbs in COL2DB.items():
            have = norm(rrow.get(col, ""))
            vset = set()
            for db in dbs:
                vset |= verified.get(db, set())
            for hid in have:
                if hid in vset:
                    tally[f"{col}:confirmed"] += 1
                elif vset:
                    tally[f"{col}:wrong"] += 1
                    audit.append((rid, col, "wrong", hid, ";".join(sorted(vset)[:6])))
                else:
                    tally[f"{col}:unverified"] += 1
            if not have and vset:
                rxn_gain[col].add(rid)
                audit.append((rid, col, "new", "", ";".join(sorted(vset)[:6])))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rxns", "column", "status", "existing_id", "structure_verified"])
        w.writerows(audit)
    print(f"reactions structure-matched: {matched}  "
          f"(no_structure {tally['no_structure']}, no_match {tally['no_match']}, degenerate {tally['degenerate']})")
    print(f"{'namespace':16}{'confirmed':>10}{'wrong':>7}{'unverif':>9}{'rxn_gaining':>12}")
    for col in COL2DB:
        print(f"  {col:14}{tally[col+':confirmed']:>10}{tally[col+':wrong']:>7}"
              f"{tally[col+':unverified']:>9}{len(rxn_gain[col]):>12}")
    print(f"\nwrote {OUT}")
    return 0


rrows = {r["rxns"]: r for r in csv.DictReader(Path("model/reactions.tsv").open(encoding="utf-8", newline=""), delimiter="\t")}

if __name__ == "__main__":
    raise SystemExit(main())
