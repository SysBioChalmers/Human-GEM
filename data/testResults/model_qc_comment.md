## Model quality report

:warning: **3** pre-existing finding(s), no regressions vs `main`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (20) | :warning: **3** pre-existing |
| Model file &amp; metabolic tasks (6) | :white_check_mark: all pass |
| MEMOTE | :sparkles: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/model_qc_summary.md) |

**Improved:**
- Cross-refs inconsistent across compartments: **0**
- Mass-imbalanced reactions: [**69**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**172**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**57**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/qc_structure_consistency.csv)

**Pre-existing:**
- Reactions split across compartments: [**26**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/qc_split_compartments.csv)
- Reactions flagged by MACAW dead-end test: [**1111**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/macaw_results.tsv)
- Reactions flagged as MACAW duplicates: [**340**](https://github.com/SysBioChalmers/Human-GEM/blob/release/2.1.0/data/testResults/macaw_results.tsv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `main` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
