## Model QC checks

:white_check_mark: **All checks fine**, no changes compared to `develop`.

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Growth (biomass producible) | 125 | new | :white_check_mark: |
| Metabolites missing formula | 7 | new | :new: |
| Metabolites missing charge | 0 | new | :new: |
| Reaction bound / GPR issues | 0 | new | :new: |
| Malformed cross-references | 59 | new | :new: |
| Cross-refs inconsistent across compartments | 59 | new | :new: |
| MEMOTE score (%) | n/a | n/a | :question: |

Detailed findings are committed to `data/testResults/`: `qc_metabolite_completeness.csv`, `qc_reaction_sanity.csv`, `qc_annotation_issues.csv`. The full MEMOTE result is uploaded as a build artifact. Look there to see which reactions, metabolites or identifiers changed.

Results for commit d95c81d.
