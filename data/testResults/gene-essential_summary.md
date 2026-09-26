### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3566 | 130 | 0.3936 | 88 | 0.668 | 0.3349 | 0.168 | 51 |
| GBM | 0.3246 | 139 | 0.3626 | 98 | 0.6627 | 0.3121 | 0.162 | 48 |
| HCT116 | 0.3706 | 131 | 0.3748 | 103 | 0.6704 | 0.3364 | 0.1803 | 41 |
| HELA | 0.312 | 169 | 0.311 | 140 | 0.6715 | 0.2749 | 0.1461 | 40 |
| RPE1 | 0.2619 | 170 | 0.3114 | 123 | 0.6591 | 0.2606 | 0.1358 | 50 |
| all |  |  |  |  | 0.6661 | 0.3033 | 0.1587 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **0** gene(s) directly affected plus **0** more that share a metabolite with one of them.

**No change:** 0 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remove-mar01364/data/testResults/gene-essential-diff.csv).
