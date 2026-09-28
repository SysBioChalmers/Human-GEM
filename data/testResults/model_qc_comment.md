## Model quality report

:x: **2** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (19) | :x: **2** regression(s) |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1000-ca-atpases/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Model / annotation-table inconsistencies: [**9**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1000-ca-atpases/data/testResults/qc_annotation_consistency.csv)
- Reactions flagged as MACAW duplicates: [**340**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1000-ca-atpases/data/testResults/macaw_results.tsv)

**Improved:**
- Reactions split across compartments: [**26**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1000-ca-atpases/data/testResults/qc_split_compartments.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
