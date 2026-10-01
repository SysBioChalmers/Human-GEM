"""Estimate the Gibbs energy of Human-GEM reactions with eQuilibrator.

Writes data/thermodynamics/reactionDeltaG.tsv, the table that the reversibility check
in code/test/qcModelChecks.py reads. The estimate is ΔG'm, the transformed Gibbs
energy with every reactant at 1 mM, at the pH of the reaction's compartment, an ionic
strength of 0.25 M and pMg 3, from component contribution (eQuilibrator).

Each metabolite is matched to an eQuilibrator compound through its cross-references
(MetaNetX, KEGG, ChEBI, HMDB, LipidMaps, BiGG, SEED), then its InChI, then its InChI
key, and last by name (a close match, NAME_SCORE or better). A match is only accepted
when the compound has the same formula as the model metabolite, once both are taken to
their neutral form (hydrogen minus charge), so that a cross-reference to a related,
generic, or differently oxidised compound is not used. A name match can still pick a
stereo- or positional isomer; isomers differ little in formation energy, which is well
within the margins the reversibility check uses.

A metabolite that is not in eQuilibrator at all can still be estimated through a proxy:
a compound with the same neutral formula whose structure is close to the model's
(Tanimoto similarity of Morgan fingerprints of at least PROXY_SIMILARITY), typically a
double-bond or chain-position isomer. Reactions estimated with a proxy are marked
"proxy"; the check reports them but never fails on them.

Reactions are skipped when:
  * a metabolite has no accepted match (often generic or lumped metabolites, and
    proteins such as thioredoxin or cytochrome c);
  * the same compound, or H+, occurs in two compartments (a transport, or a
    proton-coupled reaction such as ATP synthase: its ΔG depends on the membrane
    potential and the pH gradient, not only on chemistry);
  * FAD/FADH2, FMN/FMNH2 or ubiquinone/ubiquinol is the redox couple: these are bound
    to their enzymes or sit in the membrane, and their redox potential differs from
    the free compound that eQuilibrator describes (acyl-CoA dehydrogenases would all
    look uphill, and succinate dehydrogenase could not run backwards);
  * the reaction is not balanced in mass (hydrogen included) or charge in the model,
    so its ΔG would describe a different reaction.
Reactions spanning compartments without transporting a compound (such as GPD2) are
estimated at the pH of the first compartment and marked "multi".

Needs equilibrator-api (and rdkit for the InChI key lookups); it is not run in CI. Run
it after curation that adds or changes reactions, and commit the table.

Usage:
    python code/qc/estimateReactionDeltaG.py [--model model/Human-GEM.yml]
        [--out data/thermodynamics/reactionDeltaG.tsv] [--log FILE]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]

# Compartment pH, as used by the published thermodynamic analyses of Human-GEM.
PH = {"c": 7.2, "n": 7.2, "r": 7.2, "m": 8.0, "x": 7.0, "g": 6.6, "l": 5.0, "e": 7.4}
IONIC_STRENGTH = "0.25M"
NAME_SCORE = 0.9
PROXY_SIMILARITY = 0.9
P_MG = 3.0

NAMESPACES = (
    ("metMetaNetXID", "metanetx.chemical:"),
    ("metKEGGID", "kegg:"),
    ("metChEBIID", "chebi:"),
    ("metHMDBID", "hmdb:"),
    ("metLipidMapsID", "lipidmaps:"),
    ("metBiGGID", "bigg.metabolite:"),
    ("metSeedID", "seed.compound:"),
)


def stoichiometry_hash(stoich: dict[str, float]) -> str:
    """Short hash of a reaction's stoichiometry, to tell when an estimate is stale."""
    text = ";".join(f"{m}:{stoich[m]:g}" for m in sorted(stoich))
    return hashlib.sha1(text.encode()).hexdigest()[:10]


ATOMIC_NUMBER = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "P": 15, "S": 16, "Cl": 17,
                 "K": 19, "Ca": 20, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Se": 34,
                 "Br": 35, "Mo": 42, "I": 53}
BOUND_COUPLES = ({"FAD", "FADH2"}, {"FMN", "FMNH2"}, {"ubiquinone", "ubiquinol"})


def _atoms(formula: str) -> Counter | None:
    """Atom counts of a formula, or None if it is generic (R, X, ...) or unparsable."""
    if not formula or re.search(r"[^A-Za-z0-9]", formula):
        return None
    bag: Counter = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if el not in ATOMIC_NUMBER:
            return None
        bag[el] += int(n or 1)
    return bag


def _neutral(bag: Counter, charge: float) -> Counter:
    """Atom counts of the neutral form: hydrogen minus the charge."""
    out = Counter(bag)
    out["H"] -= charge
    return +out if all(v >= 0 for v in out.values()) else out


def load_model(path: Path):
    top = dict(yaml.load(path.read_text(), Loader=yaml.CSafeLoader))
    mets = {}
    for m in top["metabolites"]:
        m = dict(m)
        charge = m.get("charge")
        mets[m["id"]] = {"name": str(m.get("name", "")), "formula": str(m.get("formula") or ""),
                         "charge": float(charge) if charge is not None and charge != "" else None,
                         "compartment": m["compartment"]}
    rxns = []
    for r in top["reactions"]:
        r = dict(r)
        rxns.append({"id": r["id"], "stoich": {k: float(v) for k, v in dict(r["metabolites"]).items()},
                     "lb": float(r["lower_bound"]), "ub": float(r["upper_bound"])})
    return mets, rxns


class Matcher:
    def __init__(self, cc, annotations: dict[str, dict]):
        self.cc = cc
        self.ann = annotations
        self.cache: dict[str, tuple] = {}
        try:
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            self.chem = Chem
        except ImportError:
            self.chem = None

    def _candidates(self, row: dict):
        for key, prefix in NAMESPACES:
            for value in (row.get(key) or "").split(";"):
                value = value.strip()
                if not value:
                    continue
                if key == "metChEBIID":
                    value = value.replace("CHEBI:", "")
                yield key, lambda p=prefix, v=value: self.cc.get_compound(p + v)
        inchi = (row.get("metInChI") or "").strip()
        if inchi:
            yield "InChI", lambda: self.cc.get_compound_by_inchi(inchi)
        if self.chem is not None:
            key = None
            try:
                mol = self.chem.MolFromInchi(inchi) if inchi else (
                    self.chem.MolFromSmiles(row["metSmiles"]) if row.get("metSmiles") else None)
                key = self.chem.MolToInchiKey(mol) if mol is not None else None
            except Exception:  # noqa: BLE001 - malformed structures are skipped
                key = None
            if key:
                ccache = self.cc.ccache
                yield "InChIKey", lambda: (ccache.search_compound_by_inchi_key(key) or [None])[0]

    def _fingerprint(self, inchi: str = "", smiles: str = ""):
        if self.chem is None:
            return None
        from rdkit.Chem import rdFingerprintGenerator
        try:
            mol = self.chem.MolFromInchi(inchi) if inchi else self.chem.MolFromSmiles(smiles)
        except Exception:  # noqa: BLE001
            return None
        if mol is None:
            return None
        if not hasattr(self, "_morgan"):
            self._morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return self._morgan.GetFingerprint(mol)

    def _formula_index(self):
        """Neutral formula -> [(compound id, InChI)] for every eQuilibrator compound with a
        group vector. Built once, on the first proxy lookup."""
        if hasattr(self, "_index"):
            return self._index
        import pickle

        from sqlalchemy import text
        index: dict = {}
        rows = self.cc.ccache.session.execute(text(
            "SELECT id, atom_bag, inchi FROM compounds WHERE atom_bag IS NOT NULL AND group_vector IS NOT NULL"))
        for cid, blob, inchi in rows:
            bag = Counter(pickle.loads(blob))
            electrons = bag.pop("e-", None)
            if electrons is None or not inchi or any(el not in ATOMIC_NUMBER for el in bag):
                continue
            charge = sum(ATOMIC_NUMBER[el] * n for el, n in bag.items()) - electrons
            key = tuple(sorted(_neutral(bag, charge).items()))
            index.setdefault(key, []).append((cid, inchi))
        self._index = index
        return index

    def proxy(self, met: str, formula: str, charge: float | None):
        """(compound, similarity) of the closest same-formula compound, or (None, 0)."""
        from rdkit import DataStructs
        atoms = _atoms(formula)
        row = self.ann.get(met, {})
        if atoms is None or charge is None:
            return None, 0.0
        target = self._fingerprint(row.get("metInChI") or "", row.get("metSmiles") or "")
        if target is None:
            return None, 0.0
        best, best_sim = None, 0.0
        for cid, inchi in self._formula_index().get(tuple(sorted(_neutral(atoms, charge).items())), []):
            fp = self._fingerprint(inchi)
            if fp is None:
                continue
            sim = DataStructs.TanimotoSimilarity(target, fp)
            if sim > best_sim:
                best, best_sim = cid, sim
        if best is None or best_sim < PROXY_SIMILARITY:
            return None, best_sim
        return self.cc.ccache.get_compound_by_internal_id(best), best_sim

    def _by_name(self, name: str):
        try:
            hits = self.cc.ccache.search(name)
        except Exception:  # noqa: BLE001 - no hit raises ValueError
            return []
        return [compound for compound, score in hits if score >= NAME_SCORE]

    def match(self, met: str, formula: str, charge: float | None, name: str = ""):
        """(compound, how) for a metabolite, or (None, reason).

        The match is shared by all compartments of a metabolite (same base id).
        """
        base = met[:-1]
        if base in self.cache:
            return self.cache[base]
        atoms = _atoms(formula)
        want = None if atoms is None or charge is None else _neutral(atoms, charge)
        result = (None, "generic formula or no charge" if want is None else "no match")
        row = self.ann.get(met, {})
        if want is not None:
            candidates = list(self._candidates(row))
            if name:
                candidates += [("name", lambda c=c: c) for c in self._by_name(name)]
            for how, get in candidates:
                try:
                    compound = get()
                except Exception:  # noqa: BLE001 - lookups of unknown ids raise
                    compound = None
                if compound is None:
                    continue
                bag = Counter(compound.atom_bag or {})
                electrons = bag.pop("e-", None)
                if electrons is None or any(el not in ATOMIC_NUMBER for el in bag):
                    continue
                cpd_charge = sum(ATOMIC_NUMBER[el] * n for el, n in bag.items()) - electrons
                if _neutral(bag, cpd_charge) == want:
                    result = (compound, how)
                    break
                result = (None, "formula differs")
        if result[0] is None and result[1] in ("no match", "formula differs"):
            compound, similarity = self.proxy(met, formula, charge)
            if compound is not None:
                result = (compound, f"proxy ({similarity:.2f})")
        self.cache[base] = result
        return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model", type=Path, default=ROOT / "model" / "Human-GEM.yml")
    parser.add_argument("--annotation", type=Path, default=ROOT / "model" / "metabolites.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "thermodynamics" / "reactionDeltaG.tsv")
    parser.add_argument("--log", type=Path, help="also write why each reaction was skipped")
    args = parser.parse_args(argv)

    from equilibrator_api import Q_, ComponentContribution, Reaction

    mets, rxns = load_model(args.model)
    with args.annotation.open() as fh:
        ann = {r["mets"]: r for r in csv.DictReader(fh, delimiter="\t")}
    cc = ComponentContribution()
    cc.ionic_strength = Q_(IONIC_STRENGTH)
    cc.p_mg = Q_(P_MG)
    matcher = Matcher(cc, ann)

    rows, log = [], []
    for n, rxn in enumerate(rxns):
        stoich = rxn["stoich"]
        chem = {m: v for m, v in stoich.items() if mets[m]["name"] != "H+"}
        bases = [m[:-1] for m in chem]
        if len(chem) < 2:
            log.append((rxn["id"], "exchange or single metabolite"))
            continue
        if len(set(bases)) < len(bases):
            log.append((rxn["id"], "transport"))
            continue
        if len({mets[m]["compartment"] for m in stoich if mets[m]["name"] == "H+"}) > 1:
            log.append((rxn["id"], "proton-coupled across a membrane"))
            continue
        names = {mets[m]["name"] for m in stoich}
        if any(couple <= names for couple in BOUND_COUPLES):
            log.append((rxn["id"], "bound redox couple"))
            continue
        net: Counter = Counter()
        charge = 0.0
        balanced = True
        for m, v in stoich.items():
            atoms = _atoms(mets[m]["formula"])
            if atoms is None or mets[m]["charge"] is None:
                balanced = False
                break
            for el, k in atoms.items():
                net[el] += k * v
            charge += mets[m]["charge"] * v
        if balanced and (any(abs(x) > 1e-9 for x in net.values()) or abs(charge) > 1e-9):
            log.append((rxn["id"], "not balanced in the model"))
            continue
        sparse, missing, proxied = {}, [], False
        for m, v in chem.items():
            compound, how = matcher.match(m, mets[m]["formula"], mets[m]["charge"], mets[m]["name"])
            if compound is None:
                missing.append(f"{m} ({how})")
                continue
            proxied = proxied or how.startswith("proxy")
            sparse[compound] = sparse.get(compound, 0) + v
        if missing:
            log.append((rxn["id"], "unmatched: " + ", ".join(missing)))
            continue
        compartments = sorted({mets[m]["compartment"] for m in stoich})
        first = mets[next(iter(stoich))]["compartment"]
        cc.p_h = Q_(PH.get(first, 7.2))
        try:
            dg = cc.physiological_dg_prime(Reaction(sparse))
            mean, sd = dg.value.m_as("kJ/mol"), dg.error.m_as("kJ/mol")
        except Exception as exc:  # noqa: BLE001
            log.append((rxn["id"], f"eQuilibrator error: {str(exc)[:60]}"))
            continue
        kind = "multi" if len(compartments) > 1 else "single"
        if proxied:
            kind += "+proxy"
        rows.append((rxn["id"], "".join(compartments), kind,
                     f"{mean:.1f}", f"{sd:.1f}", stoichiometry_hash(stoich)))
        if n % 1000 == 0:
            print(n, len(rows), file=sys.stderr, flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["reaction", "compartments", "type", "dGm_kJ_per_mol", "sd_kJ_per_mol", "stoichiometry_hash"])
        w.writerows(sorted(rows))
    if args.log:
        with args.log.with_suffix(".matches.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            for base, (compound, how) in sorted(matcher.cache.items()):
                w.writerow([base, how, compound.id if compound is not None else ""])
        with args.log.open("w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(log)
    print(f"estimated {len(rows)} of {len(rxns)} reactions; skipped {len(log)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
