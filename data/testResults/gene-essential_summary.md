### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3494 | 139 | 0.3852 | 96 | 0.6671 | 0.3267 | 0.1671 | 52 |
| GBM | 0.3234 | 143 | 0.3602 | 102 | 0.6624 | 0.3093 | 0.1614 | 48 |
| HCT116 | 0.3736 | 129 | 0.3782 | 101 | 0.6731 | 0.3403 | 0.1797 | 41 |
| HELA | 0.3158 | 166 | 0.319 | 134 | 0.6688 | 0.2768 | 0.1456 | 43 |
| RPE1 | 0.2642 | 167 | 0.3145 | 120 | 0.6597 | 0.2628 | 0.1361 | 50 |
| all |  |  |  |  | 0.6662 | 0.3033 | 0.1582 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **0** gene(s) directly affected plus **2** more that share a metabolite with one of them.

**No change:** 2 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/feat/dihydroorotate-transport/data/testResults/gene-essential-diff.csv).
