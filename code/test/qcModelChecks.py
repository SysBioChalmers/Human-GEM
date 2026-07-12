"""Consolidated model quality-control checks for Human-GEM.

A single entry point for the structural checks. Each check writes a detailed,
diff-friendly CSV under data/testResults/ so that whatever is wrong is spelled
out in a committed file (not only in the workflow log), and buildReport.py turns
the same files into the pull-request comment.

Checks split into two kinds:

  * Gates (fail the build): a model that has these problems is broken.
      - duplicate keys inside an !!omap entry (cobra cannot load the model);
      - reactions with no metabolites;
      - the model and its annotation tables (reactions.tsv / metabolites.tsv /
        genes.tsv) disagree, or a deprecated identifier is used;
      - the model cannot produce biomass under its default constraints (the
        blocking biomass precursors are written out so it can be fixed).
  * Reports (do not fail): quality metrics tracked with a delta versus the
      target branch.
      - metabolites missing a formula or a charge;
      - reaction bound / GPR sanity;
      - exact-duplicate reactions (same stoichiometry);
      - metabolites and genes not used by any reaction.

Usage:
    python code/test/qcModelChecks.py
"""

import csv
import sys
from collections import Counter, defaultdict

import cobra
import yaml

MODEL_FILE = "model/Human-GEM.yml"
GENES_TSV = "model/genes.tsv"
REACTIONS_TSV = "model/reactions.tsv"
METABOLITES_TSV = "model/metabolites.tsv"
DEPRECATED_RXN_TSV = "data/deprecatedIdentifiers/deprecatedReactions.tsv"
DEPRECATED_MET_TSV = "data/deprecatedIdentifiers/deprecatedMetabolites.tsv"

RESULTS = "data/testResults"
DUP_KEYS_CSV = f"{RESULTS}/qc_duplicate_keys.csv"
DUP_RXN_CSV = f"{RESULTS}/qc_duplicate_reactions.csv"
EMPTY_RXN_CSV = f"{RESULTS}/qc_empty_reactions.csv"
ANNOTATION_CONSISTENCY_CSV = f"{RESULTS}/qc_annotation_consistency.csv"
UNUSED_CSV = f"{RESULTS}/qc_unused_entities.csv"
COMPLETENESS_CSV = f"{RESULTS}/qc_metabolite_completeness.csv"
REACTION_SANITY_CSV = f"{RESULTS}/qc_reaction_sanity.csv"
GROWTH_TXT = f"{RESULTS}/qc_growth.txt"
GROWTH_BLOCKERS_CSV = f"{RESULTS}/qc_growth_blockers.csv"

GROWTH_TOLERANCE = 1e-6

# Pseudo-metabolites (generic class sinks and biomass pools) intrinsically
# have no molecular formula, so they are excluded from the completeness report.
PSEUDO_METABOLITES = {
    "MAM10001", "MAM10002", "MAM10003",              # steroids, xenobiotics, arachidonate derivatives
    "MAM10012", "MAM10013", "MAM10014", "MAM10015",  # cofactor/protein/lipid/metabolite pool_biomass
}


def _tsv_column(path: str, column: str) -> list[str]:
    with open(path, newline="", encoding="utf-8") as fh:
        return [row[column] for row in csv.DictReader(fh, delimiter="\t")]


def _write_csv(path: str, header: list[str], rows: list) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)


# --------------------------------------------------------------------------- #
# Gate: duplicate keys inside an !!omap entry (cobra cannot load such a model)
# --------------------------------------------------------------------------- #
def check_no_duplicate_keys(model_file: str) -> list[tuple]:
    """Return [(entry_id, scope, key, first_line, dup_line)] for duplicate !!omap keys.

    Every metabolite, reaction and gene is an !!omap; a repeated key inside one
    (two 'name' fields, or the same metabolite twice in a stoichiometry) is
    written and re-read happily by RAVEN but makes cobra.io.load_yaml_model raise
    a bare AssertionError with no location. This scan names the entry, key and
    line numbers.
    """
    stack: list[dict] = []
    expect: list[str] = []
    firsts: list[list] = []
    dups: list[tuple] = []

    def owner() -> str:
        for ctx in reversed(stack):
            if ctx["type"] == "omap" and ctx["id"]:
                return ctx["id"]
        return "(top level)"

    with open(model_file, encoding="utf-8") as fh:
        for event in yaml.parse(fh):
            kind = type(event).__name__
            line = event.start_mark.line + 1
            if kind == "SequenceStartEvent":
                is_omap = event.tag == "tag:yaml.org,2002:omap"
                stack.append({"type": "omap" if is_omap else "seq", "keys": {}, "id": None})
            elif kind == "SequenceEndEvent":
                stack.pop()
            elif kind == "MappingStartEvent":
                stack.append({"type": "map"})
                expect.append("key")
                firsts.append([None, None])
            elif kind == "MappingEndEvent":
                stack.pop()
                first_key, first_val = firsts.pop()
                expect.pop()
                if stack and stack[-1]["type"] == "omap" and first_key is not None:
                    parent = stack[-1]
                    if first_key in parent["keys"]:
                        scope = "field" if parent["id"] else "member"
                        dups.append((owner(), scope, first_key, parent["keys"][first_key], line))
                    else:
                        parent["keys"][first_key] = line
                    if first_key == "id" and parent["id"] is None:
                        parent["id"] = first_val
            elif kind == "ScalarEvent" and stack and stack[-1]["type"] == "map":
                i = len(expect) - 1
                if expect[i] == "key":
                    if firsts[i][0] is None:
                        firsts[i][0] = event.value
                    expect[i] = "val"
                else:
                    if firsts[i][1] is None:
                        firsts[i][1] = event.value
                    expect[i] = "key"

    _write_csv(DUP_KEYS_CSV, ["entry", "scope", "key", "first_line", "duplicate_line"],
               sorted(dups))
    return dups


# --------------------------------------------------------------------------- #
# Gate: reactions with no metabolites
# --------------------------------------------------------------------------- #
def check_empty_reactions(model: cobra.Model) -> list[str]:
    empty = sorted(r.id for r in model.reactions if len(r.metabolites) == 0)
    _write_csv(EMPTY_RXN_CSV, ["reaction"], [(rid,) for rid in empty])
    return empty


# --------------------------------------------------------------------------- #
# Gate: the model and its annotation tables must agree; no deprecated ids
# --------------------------------------------------------------------------- #
def check_annotation_consistency(model: cobra.Model) -> list[tuple]:
    """Compare model ids against reactions/metabolites/genes.tsv and the deprecated
    lists. Returns [(kind, id, issue)]."""
    issues: list[tuple] = []

    def compare(kind: str, in_model: set[str], in_tsv: set[str], deprecated: set[str]):
        for used in sorted(in_model & deprecated):
            issues.append((kind, used, "deprecated identifier used"))
        for missing in sorted(in_model - in_tsv):
            issues.append((kind, missing, f"in model but not in the {kind} table"))
        for orphan in sorted(in_tsv - in_model):
            issues.append((kind, orphan, f"in the {kind} table but not in the model"))

    compare("reaction",
            {r.id for r in model.reactions}, set(_tsv_column(REACTIONS_TSV, "rxns")),
            set(_tsv_column(DEPRECATED_RXN_TSV, "rxns")))
    compare("metabolite",
            {m.id for m in model.metabolites}, set(_tsv_column(METABOLITES_TSV, "mets")),
            set(_tsv_column(DEPRECATED_MET_TSV, "mets")))
    # genes.tsv has no deprecated list; only check presence in both directions.
    compare("gene",
            {g.id for g in model.genes}, set(_tsv_column(GENES_TSV, "genes")), set())

    # The 'spontaneous' column in reactions.tsv must be numeric (RAVEN reads it).
    def _numeric(value: str) -> bool:
        try:
            float(value)
            return True
        except ValueError:
            return False

    non_numeric = [v for v in _tsv_column(REACTIONS_TSV, "spontaneous") if v.strip() and not _numeric(v)]
    if non_numeric:
        issues.append(("reactions.tsv", "spontaneous",
                       f"{len(non_numeric)} non-numeric value(s) in the spontaneous column"))

    _write_csv(ANNOTATION_CONSISTENCY_CSV, ["kind", "id", "issue"], issues)
    return issues


# --------------------------------------------------------------------------- #
# Report: exact-duplicate reactions (identical stoichiometry)
# --------------------------------------------------------------------------- #
def check_duplicate_reactions(model: cobra.Model) -> int:
    """Group reactions that share an identical metabolite -> coefficient mapping.

    This is the strict "truly identical" duplicate. Near-duplicates (reverse
    direction, different coefficients, different electron carriers) are the
    remit of the MACAW duplicate_test report, which keeps its own detail.
    """
    by_signature: dict[frozenset, list[str]] = defaultdict(list)
    for rxn in model.reactions:
        signature = frozenset((met.id, coeff) for met, coeff in rxn.metabolites.items())
        if signature:
            by_signature[signature].append(rxn.id)

    rows: list[tuple] = []
    group = 0
    for members in by_signature.values():
        if len(members) > 1:
            group += 1
            example = model.reactions.get_by_id(members[0]).build_reaction_string()
            for rid in sorted(members):
                rows.append((group, rid, example))
    rows.sort()
    _write_csv(DUP_RXN_CSV, ["group", "reaction", "equation"], rows)
    return group


# --------------------------------------------------------------------------- #
# Report: metabolites / genes not used by any reaction
# --------------------------------------------------------------------------- #
def check_unused_entities(model: cobra.Model) -> tuple[int, int]:
    used_mets = {m.id for r in model.reactions for m in r.metabolites}
    used_genes = {g.id for r in model.reactions for g in r.genes}
    rows = [("metabolite", m.id) for m in model.metabolites if m.id not in used_mets]
    rows += [("gene", g.id) for g in model.genes if g.id not in used_genes]
    rows.sort()
    _write_csv(UNUSED_CSV, ["kind", "id"], rows)
    n_met = sum(1 for k, _ in rows if k == "metabolite")
    n_gene = sum(1 for k, _ in rows if k == "gene")
    return n_met, n_gene


# --------------------------------------------------------------------------- #
# Report: metabolites missing a formula or a charge
# --------------------------------------------------------------------------- #
def check_metabolite_completeness(model: cobra.Model) -> tuple[int, int]:
    rows = []
    for met in model.metabolites:
        if met.id[:-1] in PSEUDO_METABOLITES:
            continue
        missing_formula = not (met.formula or "").strip()
        missing_charge = met.charge is None
        if missing_formula or missing_charge:
            rows.append((met.id, met.name or "",
                         "yes" if missing_formula else "", "yes" if missing_charge else ""))
    rows.sort()
    _write_csv(COMPLETENESS_CSV, ["metabolite", "name", "missing_formula", "missing_charge"], rows)
    return sum(1 for r in rows if r[2]), sum(1 for r in rows if r[3])


# --------------------------------------------------------------------------- #
# Report: reaction bound / GPR sanity
# --------------------------------------------------------------------------- #
def check_reaction_sanity(model: cobra.Model) -> int:
    documented = set(_tsv_column(GENES_TSV, "genes"))
    rows = []
    for rxn in model.reactions:
        issues = []
        lb, ub = rxn.lower_bound, rxn.upper_bound
        if lb > ub:
            issues.append("lb>ub")
        if lb < -1000 or ub > 1000:
            issues.append("|bound|>1000")
        undocumented = sorted(g.id for g in rxn.genes if g.id not in documented)
        if undocumented:
            issues.append("undocumented_genes:" + ";".join(undocumented))
        if rxn.boundary and rxn.gene_reaction_rule.strip():
            issues.append("boundary_with_gpr")
        if issues:
            rows.append((rxn.id, rxn.name or "", ";".join(issues)))
    rows.sort()
    _write_csv(REACTION_SANITY_CSV, ["reaction", "name", "issues"], rows)
    return len(rows)


# --------------------------------------------------------------------------- #
# Gate: growth, with the blocking biomass precursors when it fails
# --------------------------------------------------------------------------- #
def check_growth(model: cobra.Model) -> float:
    value = model.slim_optimize()
    return float(value) if value is not None else float("nan")


def write_growth_blockers(model: cobra.Model) -> list[str]:
    """When the model cannot grow, list the biomass precursors it cannot make.

    For each reactant of the objective (biomass) reaction, add a temporary demand
    and see whether any positive flux can be driven through it. Precursors that
    cannot be produced are what to fix. Runs only on a growth failure.
    """
    blockers: list[str] = []
    try:
        objective = [r for r in model.reactions if r.objective_coefficient != 0]
        precursors = sorted({m.id for r in objective for m, c in r.metabolites.items() if c < 0})
        for met_id in precursors:
            with model:
                met = model.metabolites.get_by_id(met_id)
                demand = model.add_boundary(met, type="demand")
                model.objective = demand
                flux = model.slim_optimize()
            if flux is None or abs(flux) <= GROWTH_TOLERANCE:
                blockers.append(met_id)
    except Exception as exc:  # noqa: BLE001 - diagnostics must never mask the growth failure
        print(f"::warning::Could not compute biomass-precursor blockers ({exc}).")
    _write_csv(GROWTH_BLOCKERS_CSV, ["blocked_biomass_precursor"], [(b,) for b in blockers])
    return blockers


def main() -> int:
    gate_failed = False

    # 1. Duplicate keys must run before the cobra load: if they exist, the load
    #    below is the bare AssertionError we are trying to explain.
    duplicates = check_no_duplicate_keys(MODEL_FILE)
    if duplicates:
        for entry_id, scope, key, first_line, dup_line in duplicates:
            print(f"::error::Duplicate {scope} '{key}' in entry {entry_id} "
                  f"({MODEL_FILE} lines {first_line} and {dup_line}).")
        print(f"Found {len(duplicates)} duplicate key(s); see {DUP_KEYS_CSV}. "
              f"cobra cannot load a model with duplicate keys in an !!omap block.")
        # Cannot load the model, so cannot run the rest; fail now.
        return 1

    model = cobra.io.load_yaml_model(MODEL_FILE)
    # GLPK for the growth LP so this does not depend on a Gurobi licence.
    model.solver = "glpk"

    # Gates
    empty = check_empty_reactions(model)
    if empty:
        print(f"::error::{len(empty)} reaction(s) have no metabolites; see {EMPTY_RXN_CSV}.")
        gate_failed = True

    annotation = check_annotation_consistency(model)
    if annotation:
        print(f"::error::{len(annotation)} model/annotation-table inconsistency(ies); "
              f"see {ANNOTATION_CONSISTENCY_CSV}.")
        gate_failed = True

    growth = check_growth(model)
    grows = growth == growth and growth > GROWTH_TOLERANCE  # not NaN and positive
    with open(GROWTH_TXT, "w", encoding="utf-8") as fh:
        fh.write(f"{growth:.6g}\n")
    if not grows:
        blockers = write_growth_blockers(model)
        print(f"::error::Model cannot produce biomass under its default constraints "
              f"({len(blockers)} blocked precursor(s); see {GROWTH_BLOCKERS_CSV}).")
        gate_failed = True
    else:
        _write_csv(GROWTH_BLOCKERS_CSV, ["blocked_biomass_precursor"], [])

    # Reports (never fail the build)
    n_formula, n_charge = check_metabolite_completeness(model)
    n_reaction_issues = check_reaction_sanity(model)
    n_dup_rxn = check_duplicate_reactions(model)
    n_unused_met, n_unused_gene = check_unused_entities(model)

    print(f"Metabolites missing a formula: {n_formula}")
    print(f"Metabolites missing a charge: {n_charge}")
    print(f"Reactions with bound/GPR issues: {n_reaction_issues}")
    print(f"Exact-duplicate reaction groups: {n_dup_rxn}")
    print(f"Unused metabolites / genes: {n_unused_met} / {n_unused_gene}")
    print(f"Growth (max biomass, default constraints): {growth:.4g} "
          f"({'ok' if grows else 'NO GROWTH'})")

    return 1 if gate_failed else 0


if __name__ == "__main__":
    sys.exit(main())
