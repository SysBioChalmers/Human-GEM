## Model quality report

:x: **1** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (19) | :x: **1** regression(s) |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-metabolites/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Reactions flagged by MACAW dead-end test: [**1111**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-metabolites/data/testResults/macaw_results.tsv)

**Improved:**
- Charge-imbalanced reactions: [**172**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-metabolites/data/testResults/balance_results.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
