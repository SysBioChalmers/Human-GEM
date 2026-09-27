## Model quality report

:x: **Merge blocked: the model cannot be loaded or cannot grow, or a build gate failed.** Gates: duplicate keys **6** :x:, growth **1.9** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :x: **blocked** |
| Model file &amp; metabolic tasks (5) | :x: **2** failed |
| MEMOTE | :white_check_mark: **63.8%** (core subset) |
| Full report | [model_qc_summary.md](https://github.com/SysBioChalmers/Human-GEM/blob/fix/retinoate-4hydroxy-13cis/data/testResults/model_qc_summary.md) |

**Regression(s):**
- Duplicate `!!omap` keys: [**6**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/retinoate-4hydroxy-13cis/data/testResults/qc_duplicate_keys.csv)

**Gate failure(s):**
- YAML round-trip (cobrapy): fail
- YAML round-trip (RAVEN): fail

**Improved:**
- Structure vs formula/charge inconsistencies: [**310**](https://github.com/SysBioChalmers/Human-GEM/blob/fix/retinoate-4hydroxy-13cis/data/testResults/qc_structure_consistency.csv)

:white_check_mark: unchanged &middot; :sparkles: improved vs `develop` &middot; :warning: pre-existing, non-blocking &middot; :x: regression
