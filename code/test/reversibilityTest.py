"""Flag reaction directions that thermodynamics makes very unlikely.

The check catches the most obviously wrong directionality only. A reaction that is
not flagged is not thereby shown to have the right reversibility: most reactions
have no estimate, and a modest ΔG can be overcome by concentrations or coupling.

Two criteria, each applied to every direction a reaction's bounds leave open:

  * ΔG'm, the Gibbs energy with all reactants at 1 mM, estimated with eQuilibrator
    (data/thermodynamics/reactionDeltaG.tsv, made by
    code/qc/estimateReactionDeltaG.py). The margin is ΔG'm in the open uphill
    direction minus twice its uncertainty.
      - warning ("questionable", reported): the reversibility index of Noor et al.
        (2012) is above ln(1000), i.e. the reactant concentrations would have to
        change more than 1000-fold to run this direction. ln Γ = (2/N)·margin/RT,
        with N the sum of the absolute coefficients excluding H2O and H+, as in
        eQuilibrator and the MEMOTE thermodynamics test.
      - alarm ("impossible", fails the check): a warning whose margin is also above
        IMPOSSIBLE_DG (40 kJ/mol, a 10^7-fold mass-action ratio). Between 1 uM and
        10 mM a reactant shifts ΔG by at most about 6 kJ/mol per order of magnitude,
        so no physiological concentrations reach it. Requiring the warning as well
        keeps reactions with many reactants, where small changes in each add up, from
        raising the alarm. Estimates that rely on an isomer proxy can only warn.
  * Oxygen: a reaction that can release O2 without consuming H2O2 or superoxide is
    an alarm. Reducing O2 is so exergonic that only the disproportionation of
    reactive oxygen species (catalase, superoxide dismutase) releases O2 in human
    cells; any other O2-releasing direction runs an oxygenase or oxidase backwards.
    This needs no ΔG estimate, so it also covers reactions eQuilibrator cannot
    estimate.

A reaction listed with a reason in data/thermodynamics/reversibilityExceptions.tsv
is reported as "exception" instead of raising the alarm. An estimate whose stored
stoichiometry hash differs from the reaction is reported as "outdated" and not
judged; CI refreshes the estimates of changed reactions before running this check.

Writes data/testResults/qc_reversibility.csv and exits 1 when an alarm is raised.

Usage:
    python code/test/reversibilityTest.py
"""

import csv
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path

import cobra

MODEL_FILE = "model/Human-GEM.yml"
DELTA_G_TSV = "data/thermodynamics/reactionDeltaG.tsv"
EXCEPTIONS_TSV = "data/thermodynamics/reversibilityExceptions.tsv"
REVERSIBILITY_CSV = "data/testResults/qc_reversibility.csv"

IMPOSSIBLE_DG = 40.0  # kJ/mol, margin above which a direction raises the alarm
LN_GAMMA_WARNING = math.log(1000)  # reversibility index above which a direction warns
RT = 8.314462618e-3 * 298.15  # kJ/mol, the temperature of the estimates

WATER_BASE = "MAM02040"
PROTON_BASE = "MAM02039"
O2_BASE = "MAM02630"
ROS_BASES = ("MAM02041", "MAM02631")  # H2O2, superoxide

HEADER = ["reaction", "name", "lower_bound", "upper_bound", "dGm_kJ_per_mol", "sd_kJ_per_mol",
          "ln_gamma", "verdict", "note"]


def stoichiometry_hash(stoichiometry: dict[str, float]) -> str:
    """Short hash of a reaction's stoichiometry (metabolite id -> coefficient).

    reactionDeltaG.tsv stores it next to each estimate, so an estimate is only used
    while the reaction is unchanged. code/qc/estimateReactionDeltaG.py computes it the
    same way.
    """
    text = ";".join(f"{m}:{stoichiometry[m]:g}" for m in sorted(stoichiometry))
    return hashlib.sha1(text.encode()).hexdigest()[:10]


def _read_tsv(path: str) -> list[dict]:
    if not Path(path).exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _abs_sum_coefficients(rxn: cobra.Reaction) -> float:
    return sum(abs(c) for m, c in rxn.metabolites.items() if m.id[:-1] not in (WATER_BASE, PROTON_BASE))


def _oxygen_release(rxn: cobra.Reaction) -> str:
    """The direction in which the reaction can release O2 without consuming H2O2 or
    superoxide ("forward"/"backward"), or "" if it cannot."""
    by_base: dict[str, float] = defaultdict(float)
    for met, coef in rxn.metabolites.items():
        by_base[met.id[:-1]] += coef
    if sum(1 for met in rxn.metabolites if met.id[:-1] == O2_BASE) != 1 or len(rxn.metabolites) < 2:
        return ""  # no O2, or an O2 transport
    for sign, is_open, direction in ((1, rxn.upper_bound > 0, "forward"), (-1, rxn.lower_bound < 0, "backward")):
        if is_open and sign * by_base[O2_BASE] > 0 and not any(sign * by_base[r] < 0 for r in ROS_BASES):
            return direction
    return ""


def check_reversibility(model: cobra.Model, delta_g_tsv: str = DELTA_G_TSV,
                        exceptions_tsv: str = EXCEPTIONS_TSV) -> list[tuple]:
    """Reactions whose open directions thermodynamics makes very unlikely.

    Returns [(reaction, name, lower_bound, upper_bound, dGm, sd, ln_gamma, verdict, note)]
    with verdict "impossible" (alarm), "questionable" (warning), "exception" (an alarm
    listed in exceptions_tsv) or "outdated" (estimate made for another stoichiometry).
    """
    exceptions = {r["reaction"]: r.get("reason", "") for r in _read_tsv(exceptions_tsv)}
    findings: dict[str, dict] = {}

    for est in _read_tsv(delta_g_tsv):
        rid = est["reaction"]
        if rid not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rid)
        dg, sd = float(est["dGm_kJ_per_mol"]), float(est["sd_kJ_per_mol"])
        row = {"rxn": rxn, "dg": dg, "sd": sd, "ln_gamma": "", "notes": []}
        if stoichiometry_hash({m.id: c for m, c in rxn.metabolites.items()}) != est["stoichiometry_hash"]:
            row.update(verdict="outdated", notes=["the reaction changed since its ΔG was estimated"])
            findings[rid] = row
            continue
        # dg is for the reaction as written: positive means forward is uphill
        uphill_open = rxn.upper_bound > 0 if dg > 0 else rxn.lower_bound < 0
        margin = abs(dg) - 2 * sd
        n = _abs_sum_coefficients(rxn)
        ln_gamma = 2 * margin / (n * RT) if n else 0.0
        if not uphill_open or ln_gamma <= LN_GAMMA_WARNING:
            continue
        proxy = "proxy" in est.get("type", "")
        direction = "forward" if dg > 0 else "backward"
        note = f"can run {direction}, where ΔG'm is {abs(dg):+.0f} ± {sd:.0f} kJ/mol"
        if proxy:
            note += " (estimated with an isomer proxy)"
        row.update(ln_gamma=round(ln_gamma, 1), notes=[note],
                   verdict="impossible" if margin > IMPOSSIBLE_DG and not proxy else "questionable")
        findings[rid] = row

    for rxn in model.reactions:
        direction = _oxygen_release(rxn)
        if not direction:
            continue
        note = f"can release O2 running {direction} without consuming H2O2 or superoxide"
        row = findings.setdefault(rxn.id, {"rxn": rxn, "dg": "", "sd": "", "ln_gamma": "", "notes": []})
        row["notes"].append(note)
        row["verdict"] = "impossible"

    rows = []
    for rid, row in findings.items():
        verdict, notes = row["verdict"], row["notes"]
        if verdict == "impossible" and rid in exceptions:
            verdict, notes = "exception", notes + [exceptions[rid]]
        rxn = row["rxn"]
        rows.append((rid, rxn.name or "", rxn.lower_bound, rxn.upper_bound, row["dg"], row["sd"],
                     row["ln_gamma"], verdict, "; ".join(notes)))
    rows.sort()
    Path(REVERSIBILITY_CSV).parent.mkdir(parents=True, exist_ok=True)
    with open(REVERSIBILITY_CSV, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)
    return rows


def main() -> int:
    model = cobra.io.load_yaml_model(MODEL_FILE)
    rows = check_reversibility(model)
    impossible = [r for r in rows if r[7] == "impossible"]
    for r in impossible:
        print(f"::error::Reaction {r[0]} {r[8]}; make it irreversible in the other direction, "
              f"or list it with the evidence in {EXCEPTIONS_TSV}.")
    flagged = {r[0] for r in rows if r[7] == "exception"}
    for rid in sorted({r["reaction"] for r in _read_tsv(EXCEPTIONS_TSV)} - flagged):
        print(f"::notice::{rid} is listed in {EXCEPTIONS_TSV} but raises no alarm any more; "
              f"remove it from the list.")
    for verdict in ("impossible", "questionable", "exception", "outdated"):
        print(f"Reversibility {verdict}: {sum(1 for r in rows if r[7] == verdict)}")
    print("These flag only the most obviously wrong directions; an unflagged reaction is "
          "not thereby shown to have the right reversibility.")
    return 1 if impossible else 0


if __name__ == "__main__":
    sys.exit(main())
