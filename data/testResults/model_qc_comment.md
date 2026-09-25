## Model quality report

:x: **1** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **125** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (19) | :x: **1** regression(s) |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/880-compartment-check/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Model / annotation-table inconsistencies: [**1**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/880-compartment-check/data/testResults/qc_annotation_consistency.csv)

**Improved:**
- Reactions flagged by MACAW dead-end test: [**1135**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/880-compartment-check/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**86**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/880-compartment-check/data/testResults/balance_results.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
