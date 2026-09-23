## Model quality report

:warning: **6** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **125** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **6** pre-existing |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.2%** |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/model_qc_summary.md) |

**Pre-existing:**
- Cross-refs inconsistent across compartments: [**3**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/qc_annotation_issues.csv)
- Reactions flagged by MACAW dead-end test: [**2510**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/macaw_results.csv)
- Reactions flagged as MACAW duplicates: [**377**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/macaw_results.csv)
- Mass-imbalanced reactions: [**87**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**234**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**397**](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
