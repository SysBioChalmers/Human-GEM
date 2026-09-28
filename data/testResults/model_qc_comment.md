## Model quality report

:x: **2** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (19) | :x: **2** regression(s) |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Reactions split across compartments: [**27**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/qc_split_compartments.csv)
- Reactions flagged by MACAW dead-end test: [**1107**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/macaw_results.tsv)

**Improved:**
- Reactions flagged as MACAW duplicates: [**338**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/macaw_results.tsv)
- Charge-imbalanced reactions: [**175**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**57**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
