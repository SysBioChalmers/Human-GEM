## Model quality report

:x: **1** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:, direction alarms **0** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (24) | :x: **1** regression(s) |
| Model file &amp; metabolic tasks (8) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](model_qc_summary.md) |

**Regression(s):**
- Reaction bound / GPR issues: [**1**](qc_reaction_sanity.csv)

**Improved:**
- Mass-imbalanced reactions: [**65**](balance_results.csv)
- Charge-imbalanced reactions: [**145**](balance_results.csv)
- Structure vs formula/charge inconsistencies: [**51**](qc_structure_consistency.csv)
- Reactions flagged as MACAW duplicates: [**338**](macaw_results.tsv)
- Reactions flagged by MACAW dead-end test: [**1376**](macaw_results.tsv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression

_Reaction directions: the alarm and the warning catch only directions that thermodynamics makes very unlikely. A reaction without an alarm or warning is not thereby shown to have the right reversibility._
