## Model quality report

:warning: **4** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **125** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **4** pre-existing |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/model_qc_summary.md) |

**Improved:**
- Reactions flagged by MACAW dead-end test: [**1135**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/macaw_results.tsv)
- Structure vs formula/charge inconsistencies: [**314**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/qc_structure_consistency.csv)

**Pre-existing:**
- Cross-refs inconsistent across compartments: [**2**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/qc_annotation_issues.csv)
- Reactions flagged as MACAW duplicates: [**377**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**87**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**234**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/balance_results.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
