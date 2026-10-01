## Model quality report

:x: **1** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:, direction alarms **0** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (23) | :x: **1** regression(s) |
| Model file &amp; metabolic tasks (6) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/feat/1083-reversibility/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Reactions flagged by MACAW dead-end test: [**1137**](https://github.com/SysBioChalmers/Human-GEM/blob/feat/1083-reversibility/data/testResults/macaw_results.tsv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression

_Reaction directions: the alarm and the warning catch only directions that thermodynamics makes very unlikely. A reaction without an alarm or warning is not thereby shown to have the right reversibility._
