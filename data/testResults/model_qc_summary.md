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
| Cross-refs inconsistent across compartments | [3](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/qc_annotation_issues.csv) | 0 | :warning: |
| MEMOTE score (%) | 20.2 | 0 | :white_check_mark: |

### MACAW and mass/charge balance

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Reactions flagged by MACAW dead-end test | [2510](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/macaw_results.csv) | 0 | :warning: |
| Reactions flagged as MACAW duplicates | [377](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/macaw_results.csv) | 0 | :warning: |
| Mass-imbalanced reactions | [87](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/balance_results.csv) | 0 | :warning: |
| Charge-imbalanced reactions | [234](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/balance_results.csv) | 0 | :warning: |
| Structure vs formula/charge inconsistencies | [397](https://github.com/SysBioChalmers/Human-GEM/blob/worktree-matlab-to-python-workflows/data/testResults/qc_structure_consistency.csv) | 0 | :warning: |

### Model file and metabolic tasks

| Check | Result | |
| --- | ---: | :---: |
| YAML round-trip (cobrapy) | pass | :white_check_mark: |
| YAML round-trip (RAVEN) | pass | :white_check_mark: |
| YAML lint | pass | :white_check_mark: |
| Essential metabolic tasks | 69 passed | :white_check_mark: |
| Verification metabolic tasks | 21 passed | :white_check_mark: |

### Gene essentiality (Hart 2015)

_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` _to run it on this pull request; the result posts as its own comment._

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
