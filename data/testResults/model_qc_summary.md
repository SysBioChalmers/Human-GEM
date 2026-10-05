## Model quality report -- full detail

_This is the full, per-check breakdown behind the pull-request comment's summary table._ _Row names link to their explanation in the [testResults README](README.md)._

### Model & network checks
_Duplicate keys (model unloadable) and no growth block the merge; every other row is a non-blocking report._

| Check | Result | &Delta; vs `develop` | |
| --- | ---: | ---: | :---: |
| [Duplicate `!!omap` keys](README.md#duplicate-omap-keys) | 0 | 0 | :white_check_mark: |
| [Growth (biomass producible)](README.md#growth-biomass-producible) | 1.9 | 0 | :white_check_mark: |
| [Reactions with no metabolites](README.md#reactions-with-no-metabolites) | 0 | 0 | :white_check_mark: |
| [Model / annotation-table inconsistencies](README.md#model--annotation-table-inconsistencies) | 0 | 0 | :white_check_mark: |
| [Removed reactions or metabolites not deprecated](README.md#removed-reactions-or-metabolites-not-deprecated) | 0 | 0 | :white_check_mark: |
| [Metabolites missing formula](README.md#metabolites-missing-formula) | 0 | 0 | :white_check_mark: |
| [Metabolites missing charge](README.md#metabolites-missing-charge) | 0 | 0 | :white_check_mark: |
| [Reaction bound / GPR issues](README.md#reaction-bound--gpr-issues) | 0 | 0 | :white_check_mark: |
| [Naming issues (missing or inconsistent)](README.md#naming-issues-missing-or-inconsistent) | 0 | 0 | :white_check_mark: |
| [Exact-duplicate reaction groups](README.md#exact-duplicate-reaction-groups) | 0 | 0 | :white_check_mark: |
| [Reactions split across compartments](README.md#reactions-split-across-compartments) | [26](qc_split_compartments.csv) | 0 | :warning: |
| [Unused metabolites](README.md#unused-metabolites) | 0 | 0 | :white_check_mark: |
| [Unused genes](README.md#unused-genes) | 0 | 0 | :white_check_mark: |
| [Malformed cross-references](README.md#malformed-cross-references) | [1](qc_annotation_issues.csv) | +1 | :x: |
| [Cross-refs inconsistent across compartments](README.md#cross-refs-inconsistent-across-compartments) | 0 | 0 | :white_check_mark: |
| [Reactions flagged by MACAW dead-end test](README.md#reactions-flagged-by-macaw-dead-end-test) | [1111](macaw_results.tsv) | 0 | :warning: |
| [Reactions flagged as MACAW duplicates](README.md#reactions-flagged-as-macaw-duplicates) | [340](macaw_results.tsv) | 0 | :warning: |
| [Mass-imbalanced reactions](README.md#mass-imbalanced-reactions) | [69](balance_results.csv) | 0 | :warning: |
| [Charge-imbalanced reactions](README.md#charge-imbalanced-reactions) | [172](balance_results.csv) | 0 | :warning: |
| [Structure vs formula/charge inconsistencies](README.md#structure-vs-formulacharge-inconsistencies) | [57](qc_structure_consistency.csv) | 0 | :warning: |

### Model file and metabolic tasks

| Check | Result | |
| --- | ---: | :---: |
| [YAML round-trip (cobrapy)](README.md#yaml-round-trip-cobrapy) | pass | :white_check_mark: |
| [YAML round-trip (RAVEN)](README.md#yaml-round-trip-raven) | pass | :white_check_mark: |
| [SBML round-trip](README.md#sbml-round-trip) | pass | :white_check_mark: |
| [YAML lint](README.md#yaml-lint) | pass | :white_check_mark: |
| [Essential metabolic tasks](README.md#essential-metabolic-tasks) | 57 passed | :white_check_mark: |
| [Verification metabolic tasks](README.md#verification-metabolic-tasks) | 21 passed | :white_check_mark: |

### [MEMOTE](README.md#memote)

**Total score: 63.8%** (core subset) &nbsp; 0

| Section | Score | &Delta; vs base |
| --- | ---: | ---: |
| consistency | 42.5% | 0 |
| annotation_met | 77.7% | 0 |
| annotation_rxn | 76.4% | 0 |
| annotation_gene | 46.7% | 0 |
| annotation_sbo | 81.8% | 0 |

<details><summary>Per-test scores</summary>

| Section | Test | Score |
| --- | --- | ---: |
| Consistency | Mass Balance | 99.3% |
| Consistency | Charge Balance | 98.2% |
| Consistency | Metabolite Connectivity | 100.0% |
| Annotation - Metabolites | Presence of Metabolite Annotation | 100.0% |
| Annotation - Metabolites | Metabolite Annotations Per Database | 47.4% |
| Annotation - Metabolites | Metabolite Annotation Conformity Per Database | 63.3% |
| Annotation - Metabolites | Uniform Metabolite Identifier Namespace | 100.0% |
| Annotation - Reactions | Presence of Reaction Annotation | 100.0% |
| Annotation - Reactions | Reaction Annotations Per Database | 27.8% |
| Annotation - Reactions | Reaction Annotation Conformity Per Database | 77.8% |
| Annotation - Reactions | Uniform Reaction Identifier Namespace | 100.0% |
| Annotation - Genes | Presence of Gene Annotation | 100.0% |
| Annotation - Genes | Gene Annotations Per Database | 20.0% |
| Annotation - Genes | Gene Annotation Conformity Per Database | 20.0% |
| Annotation - SBO Terms | Metabolite General SBO Presence | 100.0% |
| Annotation - SBO Terms | Metabolite SBO:0000247 Presence | 99.9% |
| Annotation - SBO Terms | Reaction General SBO Presence | 100.0% |
| Annotation - SBO Terms | Metabolic Reaction SBO:0000176 Presence | 100.0% |
| Annotation - SBO Terms | Transport Reaction SBO:0000185 Presence | 99.5% |
| Annotation - SBO Terms | Exchange Reaction SBO:0000627 Presence | 100.0% |
| Annotation - SBO Terms | Gene General SBO Presence | 100.0% |
| Annotation - SBO Terms | Gene SBO:0000243 Presence | 100.0% |
| Annotation - SBO Terms | Biomass Reactions SBO:0000629 Presence | 100.0% |

</details>

_Full suite: 64.9%, from an earlier model version; comment_ `/run memote` _to update it._

_The total above is the fast core subset, run on every push. Comment_ `/run memote` _to run the full suite on this pull request; the summary shows it while the model is unchanged._

### [Gene essentiality (Hart 2015)](README.md#gene-essentiality-hart-2015)

_Not run automatically (it takes hours). Comment_ `/run gene-essentiality` _to run it on this pull request; the result posts as its own comment._

:x: = a count rose vs the target branch (regression) &middot; :warning: = a pre-existing non-zero finding (non-blocking) &middot; :hourglass_flowing_sand: = still running. Counts link to the CSV listing the exact entries.
