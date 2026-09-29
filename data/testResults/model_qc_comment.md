## Model quality report

:x: **Merge blocked: the model cannot be loaded or cannot grow, or a build gate failed.** Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (20) | :warning: **6** pre-existing |
| Model file &amp; metabolic tasks (6) | :x: **1** failed |
| MEMOTE | :sparkles: **63.9%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/model_qc_summary.md) |

**Gate failure(s):**
- SBML round-trip: fail

**Pre-existing:**
- Reactions split across compartments: [**26**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/qc_split_compartments.csv)
- Reactions flagged by MACAW dead-end test: [**1111**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/macaw_results.tsv)
- Reactions flagged as MACAW duplicates: [**340**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**69**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**172**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**57**](https://github.com/SysBioChalmers/Human-GEM/blob/chore/raven3-yaml-format/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
