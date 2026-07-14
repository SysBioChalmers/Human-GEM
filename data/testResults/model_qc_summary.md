## Model quality report

:warning: **6 pre-existing finding(s), no regressions vs `develop`.** Non-blocking.

### Structural checks
_Duplicate keys (model unloadable) and no growth block the merge; the other rows are non-blocking._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Duplicate `!!omap` keys | 0 | 0 | :white_check_mark: |
| Reactions with no metabolites | 0 | 0 | :white_check_mark: |
| Model / annotation-table inconsistencies | 0 | 0 | :white_check_mark: |
| Growth (biomass producible) | 125 | 0 | :white_check_mark: |

### Model QC reports

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Metabolites missing formula | 0 | 0 | :white_check_mark: |
| Metabolites missing charge | 0 | 0 | :white_check_mark: |
| Reaction bound / GPR issues | 0 | 0 | :white_check_mark: |
| Exact-duplicate reaction groups | 0 | 0 | :white_check_mark: |
| Unused metabolites | 0 | 0 | :white_check_mark: |
| Unused genes | 0 | 0 | :white_check_mark: |
| Malformed cross-references | 0 | 0 | :white_check_mark: |
| Cross-refs inconsistent across compartments | [3](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/qc_annotation_issues.csv) | -6 | :warning: |
| MEMOTE score (%) | 20.2 | 0 | :white_check_mark: |

### MACAW and mass/charge balance

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Reactions flagged by MACAW dead-end test | [2510](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/macaw_results.csv) | 0 | :warning: |
| Reactions flagged as MACAW duplicates | [377](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/macaw_results.csv) | 0 | :warning: |
| Mass-imbalanced reactions | [87](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/balance_results.csv) | 0 | :warning: |
| Charge-imbalanced reactions | [234](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/balance_results.csv) | -6 | :warning: |
| Structure vs formula/charge inconsistencies | [397](https://github.com/SysBioChalmers/Human-GEM/blob/annotation-overhaul/data/testResults/qc_structure_consistency.csv) | new | :warning: |

### Gene essentiality (Hart 2015)

| cellLine | TP | TN | FP | FN | accuracy | sensitivity | specificity | F1 | MCC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 38 | 2158 | 60 | 280 | 0.8659 | 0.1195 | 0.9729 | 0.1827 | 0.1588 |
| GBM | 34 | 2137 | 64 | 300 | 0.8564 | 0.1018 | 0.9709 | 0.1574 | 0.1276 |
| HCT116 | 47 | 2181 | 55 | 308 | 0.8599 | 0.1324 | 0.9754 | 0.2057 | 0.1906 |
| HELA | 32 | 2241 | 70 | 250 | 0.8766 | 0.1135 | 0.9697 | 0.1667 | 0.1332 |
| RPE1 | 15 | 2179 | 83 | 258 | 0.8655 | 0.05495 | 0.9633 | 0.08086 | 0.02935 |
| all | 7 | 2379 | 95 | 112 | 0.9202 | 0.05882 | 0.9616 | 0.06335 | 0.02199 |

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
