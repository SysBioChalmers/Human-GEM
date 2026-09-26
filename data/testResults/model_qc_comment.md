## Model quality report

:x: **1** regression(s) vs `develop`. Review the row(s) below. Gates: duplicate keys **0** :white_check_mark:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :x: **1** regression(s) |
| Model file &amp; metabolic tasks (5) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1085-metabolite-curation/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Removed reactions or metabolites not deprecated: [**1**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1085-metabolite-curation/data/testResults/qc_deprecation_completeness.csv)

**Improved:**
- Cross-refs inconsistent across compartments: **0**
- Charge-imbalanced reactions: [**186**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1085-metabolite-curation/data/testResults/balance_results.csv)
- Structure vs formula/charge inconsistencies: [**70**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1085-metabolite-curation/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
