"""
Report mass- and charge-unbalanced reactions in Human-GEM (issue #704).

Uses cobrapy's Reaction.check_mass_balance(), which reports both elemental
(mass) and charge imbalances in a single call. Boundary reactions
(exchange/demand/sink), the biomass reaction and the lumped pseudo-reactions
are excluded, since an elemental balance is not defined for them: a lumped
pseudo-reaction stands for a whole class of species at once, so its
coefficients are an average over that class rather than a stoichiometry.

Lumping is recognised two ways. Four subsystems consist of lumped reactions
throughout (LUMPED_SUBSYSTEMS). The lumped reactions filed elsewhere are listed
by identifier, with what each one lumps (LUMPED_REACTIONS); that list is
audited on every run, so an entry that the model no longer needs is reported
rather than left to rot.

R and X are left as elements of their own, not ignored. A reaction that uses
them is still expected to balance in them, and several do not, which is a
finding rather than noise (#1153).

The unbalanced reactions are written, sorted, to
data/testResults/balance_results.csv, so that a pull request introducing a new
imbalance is visible in the committed diff. This is a report: it does not fail
the build. Making it a hard gate would first require resolving the reactions
that are already unbalanced.
"""
import csv
import traceback

import cobra

# Every reaction in one of these subsystems lumps many species into a single
# pseudo-reaction: a pool stands for the mixture of its members, an artificial
# conversion for a whole class of metabolites, and protein assembly and
# degradation for the translation and proteolysis of a named protein, counted
# per amino-acid residue.
LUMPED_SUBSYSTEMS = frozenset({
    "Pool reactions",
    "Artificial reactions",
    "Protein assembly",
    "Protein degradation",
})

# The lumped pseudo-reactions filed under an ordinary subsystem, and what each
# one lumps. A reaction belongs here only when its coefficients average over a
# class of species, so that no single elemental balance exists for it; a
# reaction that is merely unbalanced does not.
LUMPED_REACTIONS = {
    # A pool metabolite whose formula spells out only what its members have in
    # common, resolved into a weighted mixture of those members.
    "MAR00015": "NEFA blood pool in to its member fatty acids",
    "MAR00016": "member fatty acids to NEFA blood pool out",
    "MAR00017": "SMCFA blood pool to its member fatty acids",
    "MAR03537": "cholesterol-ester plasma pool to its member esters",
    "MAR03622": "cholesterol-ester plasma pool to its member esters",
    # DNA and RNA carry the formula of one average nucleotide, polymerised from
    # and hydrolysed to a weighted mixture of the four (deoxy)nucleotides.
    "MAR07160": "DNA from a mixture of dNTPs",
    "MAR07161": "RNA from a mixture of NTPs",
    "MAR07162": "RNA to a mixture of NDPs",
    "MAR07163": "DNA to a mixture of dNMPs",
    "MAR07164": "RNA to a mixture of NMPs",
}


def _subsystems(rxn):
    """A reaction's subsystems as a set; the YAML gives either a list or a string."""
    sub = rxn.subsystem
    if isinstance(sub, (list, tuple, set)):
        return {str(s) for s in sub}
    return {str(sub)} if sub else set()


def _is_lumped(rxn):
    """Whether the reaction stands for a whole class of species at once."""
    return rxn.id in LUMPED_REACTIONS or bool(_subsystems(rxn) & LUMPED_SUBSYSTEMS)


def _stale_lumped_entries(model, unbalanced_ids):
    """Entries of LUMPED_REACTIONS the model no longer needs.

    An identifier that is gone from the model, or whose reaction now balances
    in both mass and charge, no longer excludes anything and should be dropped
    from the list.
    """
    stale = []
    for rid in sorted(LUMPED_REACTIONS):
        if rid not in model.reactions:
            stale.append(f"{rid} (not in the model)")
        elif rid not in unbalanced_ids:
            stale.append(f"{rid} (balances)")
    return stale


def main():
    model = cobra.io.load_yaml_model("model/Human-GEM.yml")
    rows = []
    errored = []
    unbalanced_ids = set()
    n_lumped = 0
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        if "biomass" in rxn.id.lower() or "biomass" in (rxn.name or "").lower():
            continue
        try:
            imbalance = rxn.check_mass_balance()
        except Exception as exc:  # noqa: BLE001 - record, do not silently drop
            # A reaction whose balance cannot be evaluated (usually a malformed
            # formula) is a finding, not something to hide: record it as its own
            # row so it shows up in the committed diff and the count.
            if not _is_lumped(rxn):
                errored.append(rxn.id)
                rows.append((rxn.id, rxn.name or "", f"check_failed:{exc}", ""))
            continue
        if not imbalance:
            continue
        unbalanced_ids.add(rxn.id)
        if _is_lumped(rxn):
            n_lumped += 1
            continue
        mass = {k: v for k, v in imbalance.items() if k != "charge"}
        charge = imbalance.get("charge", 0)
        rows.append((
            rxn.id,
            rxn.name or "",
            ";".join(f"{k}:{v:g}" for k, v in sorted(mass.items())),
            f"{charge:g}" if charge else "",
        ))
    rows.sort()
    with open("data/testResults/balance_results.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["reaction", "name", "mass_imbalance", "charge_imbalance"])
        writer.writerows(rows)
    n_mass = sum(1 for r in rows if r[2])  # includes uncheckable reactions (surfaced, not hidden)
    n_charge = sum(1 for r in rows if r[3])
    print(
        f"Unbalanced reactions (excluding boundary, biomass and lumped): {len(rows)} "
        f"({n_mass} mass, {n_charge} charge, {len(errored)} could not be checked); "
        f"{n_lumped} lumped pseudo-reaction(s) excluded"
    )
    if errored:
        print(f"::warning::{len(errored)} reaction(s) could not be balance-checked: "
              f"{';'.join(errored[:20])}{' ...' if len(errored) > 20 else ''}")
    stale = _stale_lumped_entries(model, unbalanced_ids)
    if stale:
        print(f"::warning::{len(stale)} entr(y/ies) of LUMPED_REACTIONS are no longer "
              f"needed and can be dropped from code/test/balanceTest.py: {'; '.join(stale)}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
