## Model quality report

:warning: **4** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:, direction alarms **0** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (23) | :warning: **4** pre-existing |
| Model file &amp; metabolic tasks (6) | :white_check_mark: all pass |
| MEMOTE | :sparkles: **63.9%** (core subset) |
| Full report | [model_qc_summary.md](model_qc_summary.md) |

**Improved:**
- Mass-imbalanced reactions: [**66**](balance_results.csv)
- Charge-imbalanced reactions: [**149**](balance_results.csv)
- Reactions flagged by MACAW dead-end test: [**1392**](macaw_results.tsv)

**Pre-existing:**
- Reaction directions: warning: [**74**](qc_reversibility.csv)
- Structure vs formula/charge inconsistencies: [**52**](qc_structure_consistency.csv)
- Reactions flagged as MACAW duplicates: [**340**](macaw_results.tsv)
- Reactions split across compartments: [**26**](qc_split_compartments.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression

_Reaction directions: the alarm and the warning catch only directions that thermodynamics makes very unlikely. A reaction without an alarm or warning is not thereby shown to have the right reversibility._
