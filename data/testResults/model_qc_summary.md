## Model quality report

:warning: **6 pre-existing finding(s), no regressions vs `develop`.** Non-blocking.

_Each check name links to its explanation in the [testResults README](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md)._

### Model checks
_Duplicate keys (model unloadable) and no growth block the merge; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Duplicate `!!omap` keys](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#duplicate-omap-keys) | 0 | 0 | :white_check_mark: |
| [Growth (biomass producible)](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#growth-biomass-producible) | 125 | 0 | :white_check_mark: |
| [Reactions with no metabolites](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#reactions-with-no-metabolites) | 0 | 0 | :white_check_mark: |
| [Model / annotation-table inconsistencies](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#model--annotation-table-inconsistencies) | 0 | 0 | :white_check_mark: |
| [Removed reactions or metabolites not deprecated](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#removed-reactions-or-metabolites-not-deprecated) | 0 | 0 | :white_check_mark: |
| [Metabolites missing formula](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#metabolites-missing-formula) | 0 | 0 | :white_check_mark: |
| [Metabolites missing charge](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#metabolites-missing-charge) | 0 | 0 | :white_check_mark: |
| [Reaction bound / GPR issues](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#reaction-bound--gpr-issues) | 0 | 0 | :white_check_mark: |
| [Exact-duplicate reaction groups](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#exact-duplicate-reaction-groups) | 0 | 0 | :white_check_mark: |
| [Unused metabolites](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#unused-metabolites) | 0 | 0 | :white_check_mark: |
| [Unused genes](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#unused-genes) | 0 | 0 | :white_check_mark: |
| [Malformed cross-references](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#malformed-cross-references) | 0 | 0 | :white_check_mark: |
| [Cross-refs inconsistent across compartments](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#cross-refs-inconsistent-across-compartments) | [3](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/qc_annotation_issues.csv) | 0 | :warning: |

### MACAW and mass/charge balance

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Reactions flagged by MACAW dead-end test](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#reactions-flagged-by-macaw-dead-end-test) | [2510](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/macaw_results.csv) | 0 | :warning: |
| [Reactions flagged as MACAW duplicates](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#reactions-flagged-as-macaw-duplicates) | [377](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/macaw_results.csv) | 0 | :warning: |
| [Mass-imbalanced reactions](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#mass-imbalanced-reactions) | [87](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/balance_results.csv) | 0 | :warning: |
| [Charge-imbalanced reactions](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#charge-imbalanced-reactions) | [234](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/balance_results.csv) | 0 | :warning: |
| [Structure vs formula/charge inconsistencies](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#structure-vs-formulacharge-inconsistencies) | [397](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/qc_structure_consistency.csv) | 0 | :warning: |

### Model file and metabolic tasks

| Check | Result | |
| --- | ---: | :---: |
| [YAML round-trip (cobrapy)](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#yaml-round-trip-cobrapy) | pass | :white_check_mark: |
| [YAML round-trip (RAVEN)](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#yaml-round-trip-raven) | pass | :white_check_mark: |
| [YAML lint](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#yaml-lint) | pass | :white_check_mark: |
| [Essential metabolic tasks](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#essential-metabolic-tasks) | 57 passed | :white_check_mark: |
| [Verification metabolic tasks](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#verification-metabolic-tasks) | 21 passed | :white_check_mark: |

### [MEMOTE](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#memote)

**Total score: 63.2%** (core subset) &nbsp; 0

| Section | Score | &Delta; vs base |
| --- | ---: | ---: |
| consistency | 42.4% | 0 |
| annotation_met | 73.0% | 0 |
| annotation_rxn | 72.7% | 0 |
| annotation_gene | 46.6% | -0.1 :warning: |
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
| Annotation - Genes | Gene Annotations Per Database | 80.1% |
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

### [Gene essentiality (Hart 2015)](https://github.com/SysBioChalmers/Human-GEM/blob/fix/ca-influx-channels/data/testResults/README.md#gene-essentiality-hart-2015)

_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` _to run it on this pull request; the result posts as its own comment._

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
