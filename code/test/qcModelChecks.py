"""Model quality-control checks for Human-GEM.

Three checks that complement the MACAW, balance and MEMOTE tests:

  * Metabolite formula and charge completeness. A metabolite with no chemical
    formula or no charge is silently skipped by the mass/charge balance test, so
    tracking these keeps that test honest.
  * Reaction GPR and bounds sanity. Flux bounds must satisfy lb <= ub and stay
    within the standard +/-1000 range; gene-protein-reaction rules must reference
    valid gene identifiers.
  * Growth sanity. The model must be able to produce biomass (objective reaction)
    under its default constraints.

The completeness and sanity findings are written as diff-friendly CSVs and are
reports (they do not fail the build), following the balance-test convention. The
growth check is a functional gate: a model that cannot grow is broken, so that
one fails the build.

Usage:
    python code/test/qcModelChecks.py
"""

import csv
import sys

import cobra

MODEL_FILE = "model/Human-GEM.yml"
GENES_TSV = "model/genes.tsv"
COMPLETENESS_CSV = "data/testResults/qc_metabolite_completeness.csv"
REACTION_SANITY_CSV = "data/testResults/qc_reaction_sanity.csv"
GROWTH_TXT = "data/testResults/qc_growth.txt"
GROWTH_TOLERANCE = 1e-6


def _documented_genes() -> set[str]:
    """Gene ids listed in model/genes.tsv."""
    with open(GENES_TSV, newline="", encoding="utf-8") as fh:
        return {row["genes"] for row in csv.DictReader(fh, delimiter="\t")}


def check_metabolite_completeness(model: cobra.Model) -> tuple[int, int]:
    """Write metabolites missing a formula or charge; return (n_formula, n_charge)."""
    rows = []
    for met in model.metabolites:
        missing_formula = not (met.formula or "").strip()
        missing_charge = met.charge is None
        if missing_formula or missing_charge:
            rows.append((met.id, met.name or "",
                         "yes" if missing_formula else "",
                         "yes" if missing_charge else ""))
    rows.sort()
    with open(COMPLETENESS_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["metabolite", "name", "missing_formula", "missing_charge"])
        writer.writerows(rows)
    n_formula = sum(1 for r in rows if r[2])
    n_charge = sum(1 for r in rows if r[3])
    return n_formula, n_charge


def check_reaction_sanity(model: cobra.Model) -> int:
    """Write reactions with bound or GPR problems; return the number found."""
    documented = _documented_genes()
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
    with open(REACTION_SANITY_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["reaction", "name", "issues"])
        writer.writerows(rows)
    return len(rows)


def check_growth(model: cobra.Model) -> float:
    """Return the maximum biomass flux under the model's default constraints."""
    value = model.slim_optimize()
    return float(value) if value is not None else float("nan")


def main() -> int:
    model = cobra.io.load_yaml_model(MODEL_FILE)
    # Use GLPK for the single growth LP so this check does not depend on a Gurobi
    # licence (gurobipy is installed for MEMOTE, and its bundled licence would
    # otherwise reject this genome-scale model).
    model.solver = "glpk"

    n_formula, n_charge = check_metabolite_completeness(model)
    n_reaction_issues = check_reaction_sanity(model)
    growth = check_growth(model)
    grows = growth == growth and growth > GROWTH_TOLERANCE  # not NaN and positive

    # Persist the growth value so the comment builder can compare it to the target branch.
    with open(GROWTH_TXT, "w", encoding="utf-8") as fh:
        fh.write(f"{growth:.6g}\n")

    print(f"Metabolites missing a formula: {n_formula}")
    print(f"Metabolites missing a charge: {n_charge}")
    print(f"Reactions with bound/GPR issues: {n_reaction_issues}")
    print(f"Growth (max biomass, default constraints): {growth:.4g}"
          f" ({'ok' if grows else 'NO GROWTH'})")

    if not grows:
        print("::error::Model cannot produce biomass under its default constraints.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
