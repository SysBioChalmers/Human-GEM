## Model quality report

:warning: **6** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:, direction alarms **0** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (23) | :warning: **6** pre-existing |
| Model file &amp; metabolic tasks (6) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/model_qc_summary.md) |

**Improved:**
- Structure vs formula/charge inconsistencies: [**52**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/qc_structure_consistency.csv)

**Pre-existing:**
- Reaction directions: warning: [**76**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/qc_reversibility.csv)
- Mass-imbalanced reactions: [**69**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**172**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/balance_results.csv)
- Reactions flagged as MACAW duplicates: [**340**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/macaw_results.tsv)
- Reactions flagged by MACAW dead-end test: [**1134**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/macaw_results.tsv)
- Reactions split across compartments: [**26**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1153-balance-first-tranche/data/testResults/qc_split_compartments.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression

_Reaction directions: the alarm and the warning catch only directions that thermodynamics makes very unlikely. A reaction without an alarm or warning is not thereby shown to have the right reversibility._
