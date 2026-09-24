### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3645 | 134 | 0.3796 | 101 | 0.6731 | 0.329 | 0.1685 | 45 |
| GBM | 0.288 | 169 | 0.3286 | 123 | 0.6514 | 0.2879 | 0.1625 | 53 |
| HCT116 | 0.3775 | 131 | 0.3774 | 104 | 0.6783 | 0.3449 | 0.1807 | 41 |
| HELA | 0.3283 | 166 | 0.3251 | 134 | 0.6706 | 0.2823 | 0.147 | 45 |
| RPE1 | 0.2566 | 183 | 0.3095 | 122 | 0.6623 | 0.2653 | 0.1358 | 70 |
| all |  |  |  |  | 0.6674 | 0.3022 | 0.1592 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **3** gene(s) directly affected plus **1463** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- MAT2A: 4/5 lines, knockout now blocks growth -- :warning: mixed across lines

**No change:** 1465 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/sam-synthesis-gpr/data/testResults/gene-essential-diff.csv).
