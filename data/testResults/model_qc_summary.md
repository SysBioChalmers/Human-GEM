## Model quality report

:x: **2 regression(s) vs `develop`** (this pull request increased a finding count). Review the :x: rows.

### Structural checks
_Duplicate keys (model unloadable) and no growth block the merge; the other rows are non-blocking._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Duplicate `!!omap` keys | 0 | 0 | :white_check_mark: |
| Reactions with no metabolites | 0 | 0 | :white_check_mark: |
| Model / annotation-table inconsistencies | [2](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_annotation_consistency.csv) | +2 | :x: |
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
| Cross-refs inconsistent across compartments | [3](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_annotation_issues.csv) | 0 | :warning: |

### MACAW and mass/charge balance

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Reactions flagged by MACAW dead-end test | [2511](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv) | +1 | :x: |
| Reactions flagged as MACAW duplicates | [377](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv) | 0 | :warning: |
| Mass-imbalanced reactions | [87](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv) | 0 | :warning: |
| Charge-imbalanced reactions | [234](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv) | 0 | :warning: |
| Structure vs formula/charge inconsistencies | [397](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_structure_consistency.csv) | 0 | :warning: |

### Model file and metabolic tasks

| Check | Result | |
| --- | ---: | :---: |
| YAML round-trip (cobrapy) | pass | :white_check_mark: |
| YAML round-trip (RAVEN) | pass | :white_check_mark: |
| YAML lint | pass | :white_check_mark: |
| Essential metabolic tasks | 57 passed | :white_check_mark: |
| Verification metabolic tasks | 21 passed | :white_check_mark: |

### MEMOTE

**Total score: 20.2%** (core subset) &nbsp; 0

| Section | Score | &Delta; vs base |
| --- | ---: | ---: |
| consistency | 42.4% | 0 |
| annotation_met | 25.0% | 0 |
| annotation_rxn | 25.0% | 0 |
| annotation_gene | 0.0% | 0 |
| annotation_sbo | 0.0% | 0 |

<details><summary>Per-test scores</summary>

| Section | Test | Score |
| --- | --- | ---: |
| Consistency | Stoichiometric Consistency | 100.0% |
| Consistency | Mass Balance | 0.8% |
| Consistency | Charge Balance | 2.1% |
| Consistency | Metabolite Connectivity | 0.0% |
| Consistency | Unbounded Flux In Default Medium | 100.0% |
| Annotation - Metabolites | Presence of Metabolite Annotation | 100.0% |
| Annotation - Metabolites | Metabolite Annotations Per Database | 100.0% |
| Annotation - Metabolites | Metabolite Annotation Conformity Per Database | 100.0% |
| Annotation - Metabolites | Uniform Metabolite Identifier Namespace | 0.0% |
| Annotation - Reactions | Presence of Reaction Annotation | 100.0% |
| Annotation - Reactions | Reaction Annotations Per Database | 100.0% |
| Annotation - Reactions | Reaction Annotation Conformity Per Database | 100.0% |
| Annotation - Reactions | Uniform Reaction Identifier Namespace | 0.0% |
| Annotation - Genes | Presence of Gene Annotation | 100.0% |
| Annotation - Genes | Gene Annotations Per Database | 100.0% |
| Annotation - Genes | Gene Annotation Conformity Per Database | 100.0% |
| Annotation - SBO Terms | Metabolite General SBO Presence | 100.0% |
| Annotation - SBO Terms | Metabolite SBO:0000247 Presence | 100.0% |
| Annotation - SBO Terms | Reaction General SBO Presence | 100.0% |
| Annotation - SBO Terms | Metabolic Reaction SBO:0000176 Presence | 100.0% |
| Annotation - SBO Terms | Transport Reaction SBO:0000185 Presence | 100.0% |
| Annotation - SBO Terms | Exchange Reaction SBO:0000627 Presence | 100.0% |
| Annotation - SBO Terms | Demand Reaction SBO:0000628 Presence | 100.0% |
| Annotation - SBO Terms | Sink Reactions SBO:0000632 Presence | 100.0% |
| Annotation - SBO Terms | Gene General SBO Presence | 100.0% |
| Annotation - SBO Terms | Gene SBO:0000243 Presence | 100.0% |
| Annotation - SBO Terms | Biomass Reactions SBO:0000629 Presence | 100.0% |

</details>

_The score above is the fast core subset. Comment_ `/run memote` _to run the full suite on this pull request; the score updates here when it finishes._

### Gene essentiality (Hart 2015)

_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` _to run it on this pull request; the result posts as its own comment._

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
