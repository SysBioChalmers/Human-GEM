## Model quality report

:warning: **6 pre-existing finding(s), no regressions vs `develop`.** Non-blocking. Gates: duplicate keys 0 :white_check_mark:, growth 125 :white_check_mark:.

_Row names match the headings in the [testResults README](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/README.md)._

### Model & network checks
_Duplicate keys (model unloadable) and no growth block the merge; every other row is a non-blocking report._

<details><summary>18 more check(s) unchanged vs `develop` (12 clean, 6 pre-existing finding(s)) -- show</summary>

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| Duplicate `!!omap` keys | 0 | 0 | :white_check_mark: |
| Growth (biomass producible) | 125 | 0 | :white_check_mark: |
| Reactions with no metabolites | 0 | 0 | :white_check_mark: |
| Model / annotation-table inconsistencies | 0 | 0 | :white_check_mark: |
| Removed reactions or metabolites not deprecated | 0 | 0 | :white_check_mark: |
| Metabolites missing formula | 0 | 0 | :white_check_mark: |
| Metabolites missing charge | 0 | 0 | :white_check_mark: |
| Reaction bound / GPR issues | 0 | 0 | :white_check_mark: |
| Exact-duplicate reaction groups | 0 | 0 | :white_check_mark: |
| Unused metabolites | 0 | 0 | :white_check_mark: |
| Unused genes | 0 | 0 | :white_check_mark: |
| Malformed cross-references | 0 | 0 | :white_check_mark: |
| Cross-refs inconsistent across compartments | [3](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_annotation_issues.csv) | 0 | :warning: |
| Reactions flagged by MACAW dead-end test | [2510](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv) | 0 | :warning: |
| Reactions flagged as MACAW duplicates | [377](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv) | 0 | :warning: |
| Mass-imbalanced reactions | [87](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv) | 0 | :warning: |
| Charge-imbalanced reactions | [234](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv) | 0 | :warning: |
| Structure vs formula/charge inconsistencies | [397](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_structure_consistency.csv) | 0 | :warning: |

</details>

### Model file and metabolic tasks

All 5 pass: YAML round-trip (cobrapy, RAVEN), YAML lint, 57 essential + 21 verification tasks. :white_check_mark:

### [MEMOTE](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/README.md#memote)

**Total score: 63.2%** (core subset) &nbsp; 0

| Section | Score | &Delta; vs base |
| --- | ---: | ---: |
| consistency | 42.4% | 0 |
| annotation_met | 73.0% | 0 |
| annotation_rxn | 72.7% | 0 |
| annotation_gene | 46.7% | 0 |
| annotation_sbo | 81.7% | 0 |

<details><summary>Per-test scores</summary>

| Section | Test | Score |
| --- | --- | ---: |
| Consistency | Stoichiometric Consistency | 100.0% |
| Consistency | Mass Balance | 0.8% |
| Consistency | Charge Balance | 2.1% |
| Consistency | Metabolite Connectivity | 0.0% |
| Consistency | Unbounded Flux In Default Medium | 100.0% |
| Annotation - Metabolites | Presence of Metabolite Annotation | 0.0% |
| Annotation - Metabolites | Metabolite Annotations Per Database | 62.3% |
| Annotation - Metabolites | Metabolite Annotation Conformity Per Database | 45.8% |
| Annotation - Metabolites | Uniform Metabolite Identifier Namespace | 0.0% |
| Annotation - Reactions | Presence of Reaction Annotation | 0.0% |
| Annotation - Reactions | Reaction Annotations Per Database | 75.9% |
| Annotation - Reactions | Reaction Annotation Conformity Per Database | 33.3% |
| Annotation - Reactions | Uniform Reaction Identifier Namespace | 0.0% |
| Annotation - Genes | Presence of Gene Annotation | 0.0% |
| Annotation - Genes | Gene Annotations Per Database | 80.0% |
| Annotation - Genes | Gene Annotation Conformity Per Database | 80.0% |
| Annotation - SBO Terms | Metabolite General SBO Presence | 0.0% |
| Annotation - SBO Terms | Metabolite SBO:0000247 Presence | 0.1% |
| Annotation - SBO Terms | Reaction General SBO Presence | 0.0% |
| Annotation - SBO Terms | Metabolic Reaction SBO:0000176 Presence | 0.0% |
| Annotation - SBO Terms | Transport Reaction SBO:0000185 Presence | 0.7% |
| Annotation - SBO Terms | Exchange Reaction SBO:0000627 Presence | 0.0% |
| Annotation - SBO Terms | Demand Reaction SBO:0000628 Presence | 100.0% |
| Annotation - SBO Terms | Sink Reactions SBO:0000632 Presence | 100.0% |
| Annotation - SBO Terms | Gene General SBO Presence | 0.0% |
| Annotation - SBO Terms | Gene SBO:0000243 Presence | 0.0% |
| Annotation - SBO Terms | Biomass Reactions SBO:0000629 Presence | 0.0% |

</details>

**Full suite: 64.2%** &nbsp; 0 &middot; _from the last_ `/run memote`.

_The score above is the fast core subset. Comment_ `/run memote` _to run the full suite on this pull request; the score updates here when it finishes._

### [Gene essentiality (Hart 2015)](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/README.md#gene-essentiality-hart-2015)

_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` _to run it on this pull request; the result posts as its own comment._

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
