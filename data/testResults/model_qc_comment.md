## Model quality report

:warning: **9** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:, direction alarms **0** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (24) | :warning: **9** pre-existing |
| Model file &amp; metabolic tasks (8) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](model_qc_summary.md) |

**Pre-existing:**
- Reaction directions: warning: [**70**](qc_reversibility.csv)
- Mass-imbalanced reactions: [**65**](balance_results.csv)
- Charge-imbalanced reactions: [**145**](balance_results.csv)
- Structure vs formula/charge inconsistencies: [**51**](qc_structure_consistency.csv)
- Reactions flagged as MACAW duplicates: [**338**](macaw_results.tsv)
- Twin reactions with a missing or conflicting EC number: [**52**](qc_twin_reactions.csv)
- Reactions flagged by MACAW dead-end test: [**1380**](macaw_results.tsv)
- Reactions split across compartments: [**26**](qc_split_compartments.csv)
- Unused genes: [**2**](qc_unused_entities.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression

_Reaction directions: the alarm and the warning catch only directions that thermodynamics makes very unlikely. A reaction without an alarm or warning is not thereby shown to have the right reversibility._
