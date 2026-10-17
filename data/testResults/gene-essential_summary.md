### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3382 | 150 | 0.3629 | 110 | 0.6649 | 0.3111 | 0.1668 | 50 |
| GBM | 0.3153 | 149 | 0.3415 | 114 | 0.6579 | 0.2955 | 0.1612 | 42 |
| HCT116 | 0.3715 | 131 | 0.3652 | 110 | 0.6691 | 0.3278 | 0.179 | 34 |
| HELA | 0.3075 | 174 | 0.3009 | 149 | 0.6645 | 0.2634 | 0.145 | 36 |
| RPE1 | 0.2594 | 171 | 0.2977 | 132 | 0.6539 | 0.2497 | 0.1358 | 42 |
| all |  |  |  |  | 0.6621 | 0.2896 | 0.1578 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **31** gene(s) directly affected plus **529** more that share a metabolite with one of them.

**Likely noise** (5 gene(s), <3/5 lines): CA13, FPGS, GGH, PXMP2, SLC36A1.

**No change:** 555 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
