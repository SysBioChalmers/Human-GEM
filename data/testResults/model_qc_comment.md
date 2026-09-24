## Model quality report

:warning: **5** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **125** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **5** pre-existing |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.5%** |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/model_qc_summary.md) |

**Improved:**
- Structure vs formula/charge inconsistencies: [**316**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/qc_structure_consistency.csv)

**Pre-existing:**
- Cross-refs inconsistent across compartments: [**3**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/qc_annotation_issues.csv)
- Reactions flagged by MACAW dead-end test: [**1141**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/macaw_results.tsv)
- Reactions flagged as MACAW duplicates: [**377**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**87**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**234**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/reaction-xrefs/data/testResults/balance_results.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
