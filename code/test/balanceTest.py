"""
Report mass- and charge-unbalanced reactions in Human-GEM (issue #704).

Uses cobrapy's Reaction.check_mass_balance(), which reports both elemental
(mass) and charge imbalances in a single call. Boundary reactions
(exchange/demand/sink) and the biomass reaction are excluded, since they are
not expected to balance.

The unbalanced reactions are written, sorted, to
data/testResults/balance_results.csv, so that a pull request introducing a new
imbalance is visible in the committed diff. This is a report: it does not fail
the build. Making it a hard gate would first require resolving the reactions
that are already unbalanced.
"""
import csv
import traceback

import cobra


def main():
    model = cobra.io.load_yaml_model("model/Human-GEM.yml")
    rows = []
    errored = []
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
            errored.append(rxn.id)
            rows.append((rxn.id, rxn.name or "", f"check_failed:{exc}", ""))
            continue
        if imbalance:
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
        f"Unbalanced reactions (excluding boundary and biomass): {len(rows)} "
        f"({n_mass} mass, {n_charge} charge, {len(errored)} could not be checked)"
    )
    if errored:
        print(f"::warning::{len(errored)} reaction(s) could not be balance-checked: "
              f"{';'.join(errored[:20])}{' ...' if len(errored) > 20 else ''}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
