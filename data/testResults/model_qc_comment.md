## Model quality report

:warning: **5** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **5** pre-existing |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/model_qc_summary.md) |

**Improved:**
- Reactions flagged as MACAW duplicates: [**342**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/macaw_results.tsv)

**Pre-existing:**
- Cross-refs inconsistent across compartments: [**2**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/qc_annotation_issues.csv)
- Reactions flagged by MACAW dead-end test: [**1126**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**70**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**198**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**315**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
