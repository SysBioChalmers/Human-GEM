### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3444 | 142 | 0.384 | 96 | 0.6691 | 0.3277 | 0.1676 | 55 |
| GBM | 0.3271 | 140 | 0.3586 | 103 | 0.661 | 0.3071 | 0.1615 | 44 |
| HCT116 | 0.3723 | 130 | 0.3751 | 103 | 0.6707 | 0.3358 | 0.1797 | 40 |
| HELA | 0.3124 | 169 | 0.3151 | 137 | 0.6703 | 0.2752 | 0.1457 | 43 |
| RPE1 | 0.2544 | 176 | 0.3035 | 128 | 0.6593 | 0.2573 | 0.1361 | 51 |
| all |  |  |  |  | 0.6659 | 0.3005 | 0.1584 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **3** gene(s) directly affected plus **0** more that share a metabolite with one of them.

**No change:** 3 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/collagen-secretion/data/testResults/gene-essential-diff.csv).
