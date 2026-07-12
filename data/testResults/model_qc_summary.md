## Model quality report

:information_source: First run for this comparison; no target-branch baseline yet.

### Build gates
_A red :x: blocks the merge; the count links to the CSV listing what to fix._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Duplicate `!!omap` keys | 0 | new | :white_check_mark: |
| Reactions with no metabolites | 0 | new | :white_check_mark: |
| Model / annotation-table inconsistencies | 0 | new | :white_check_mark: |
| Growth (biomass producible) | 125 | new | :white_check_mark: |

### Model QC reports

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Metabolites missing formula | [7](qc_metabolite_completeness.csv) | new | :new: |
| Metabolites missing charge | 0 | new | :new: |
| Reaction bound / GPR issues | 0 | new | :new: |
| Exact-duplicate reaction groups | 0 | new | :new: |
| Unused metabolites | 0 | new | :new: |
| Unused genes | 0 | new | :new: |
| Malformed cross-references | [59](qc_annotation_issues.csv) | new | :new: |
| Cross-refs inconsistent across compartments | [59](qc_annotation_issues.csv) | new | :new: |
| MEMOTE score (%) | [20.2](memote_score.md) | new | :new: |

### MACAW and mass/charge balance

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Reactions flagged by MACAW dead-end test | [2510](macaw_results.csv) | new | :new: |
| Reactions flagged as MACAW duplicates | [377](macaw_results.csv) | new | :new: |
| Mass-imbalanced reactions | [87](balance_results.csv) | new | :new: |
| Charge-imbalanced reactions | [240](balance_results.csv) | new | :new: |

### Gene essentiality (Hart 2015)

| cellLine | TP | TN | FP | FN | accuracy | sensitivity | specificity | F1 | MCC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 38 | 2158 | 60 | 280 | 0.8659 | 0.1195 | 0.9729 | 0.1827 | 0.1588 |
| GBM | 34 | 2137 | 64 | 300 | 0.8564 | 0.1018 | 0.9709 | 0.1574 | 0.1276 |
| HCT116 | 47 | 2181 | 55 | 308 | 0.8599 | 0.1324 | 0.9754 | 0.2057 | 0.1906 |
| HELA | 32 | 2241 | 70 | 250 | 0.8766 | 0.1135 | 0.9697 | 0.1667 | 0.1332 |
| RPE1 | 15 | 2179 | 83 | 258 | 0.8655 | 0.05495 | 0.9633 | 0.08086 | 0.02935 |
| all | 7 | 2379 | 95 | 112 | 0.9202 | 0.05882 | 0.9616 | 0.06335 | 0.02199 |

Per-finding detail is in the linked CSVs under `data/testResults/`; the full MEMOTE result is uploaded as a build artifact.
