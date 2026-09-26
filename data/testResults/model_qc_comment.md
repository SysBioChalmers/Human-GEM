## Model quality report

:x: **Merge blocked: the model cannot be loaded or cannot grow, or a build gate failed.** Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **4** pre-existing |
| Model file &amp; metabolic tasks (5) | :x: **1** failed |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/model_qc_summary.md) |

**Gate failure(s):**
- Verification metabolic tasks: 8 failed

**Improved:**
- Reactions flagged by MACAW dead-end test: [**1130**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/macaw_results.tsv)
- Structure vs formula/charge inconsistencies: [**315**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/qc_structure_consistency.csv)

**Pre-existing:**
- Cross-refs inconsistent across compartments: [**2**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/qc_annotation_issues.csv)
- Reactions flagged as MACAW duplicates: [**377**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/macaw_results.tsv)
- Mass-imbalanced reactions: [**71**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/balance_results.csv)
- Charge-imbalanced reactions: [**233**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/balance_results.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
