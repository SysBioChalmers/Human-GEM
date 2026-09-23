## Model quality report

:warning: **6** pre-existing finding(s), no regressions vs `develop`. Non-blocking. Gates: duplicate keys **0** :white_check_mark:, growth **125** :white_check_mark:.

| Section | Status |
| --- | --- |
| Model &amp; network checks (18) | :warning: **6** pre-existing |
| Model file &amp; metabolic tasks (5 gates) | :white_check_mark: all pass |
| MEMOTE | :white_check_mark: **63.2%** |
| Gene essentiality | see the gene-essentiality comment |

Pre-existing, non-blocking: Cross-refs inconsistent across compartments [3](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_annotation_issues.csv), Reactions flagged by MACAW dead-end test [2510](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv), Reactions flagged as MACAW duplicates [377](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/macaw_results.csv), Mass-imbalanced reactions [87](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv), Charge-imbalanced reactions [234](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/balance_results.csv), Structure vs formula/charge inconsistencies [397](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/qc_structure_consistency.csv).

[Full report](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/model_qc_summary.md) -- every check, MEMOTE's per-test breakdown, and how to run gene essentiality.
