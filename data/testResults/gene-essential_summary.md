### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3562 | 133 | 0.3939 | 90 | 0.6693 | 0.334 | 0.1675 | 52 |
| GBM | 0.3246 | 142 | 0.3555 | 105 | 0.6606 | 0.3059 | 0.1614 | 44 |
| HCT116 | 0.3734 | 130 | 0.3777 | 102 | 0.6718 | 0.3372 | 0.1791 | 41 |
| HELA | 0.3103 | 171 | 0.3114 | 140 | 0.6659 | 0.2712 | 0.1455 | 42 |
| RPE1 | 0.2524 | 178 | 0.2996 | 131 | 0.6584 | 0.2546 | 0.136 | 50 |
| all |  |  |  |  | 0.6651 | 0.3 | 0.1582 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **5** gene(s) directly affected plus **1626** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- CYP3A4: 5/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (12 gene(s), <3/5 lines): AKR1A1, CA13, GNPAT, H6PD, IDH1, LIPA, MPC1, MPC2, NAGK, SLC17A5, SLC17A7, SLC27A5.

**No change:** 1618 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/retinoate-4hydroxy-13cis/data/testResults/gene-essential-diff.csv).
