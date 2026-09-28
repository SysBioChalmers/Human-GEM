### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.348 | 142 | 0.3749 | 102 | 0.6674 | 0.3201 | 0.1671 | 50 |
| GBM | 0.3156 | 143 | 0.3406 | 112 | 0.6554 | 0.2944 | 0.1618 | 37 |
| HCT116 | 0.3671 | 134 | 0.3648 | 110 | 0.6689 | 0.3281 | 0.1796 | 37 |
| HELA | 0.3117 | 173 | 0.303 | 147 | 0.6648 | 0.2651 | 0.1453 | 38 |
| RPE1 | 0.2578 | 172 | 0.287 | 140 | 0.6456 | 0.241 | 0.1363 | 35 |
| all |  |  |  |  | 0.6606 | 0.2894 | 0.1583 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **29** gene(s) directly affected plus **1783** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- AKR1A1, GNPAT: 3/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (10 gene(s), <3/5 lines): FPGS, GGH, H6PD, LIPA, NAGK, PISD, SLC17A5, SLC27A5, SLC35D1, SLC36A1.

**No change:** 1800 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
