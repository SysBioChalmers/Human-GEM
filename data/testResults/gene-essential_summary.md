### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3531 | 134 | 0.3774 | 99 | 0.669 | 0.3229 | 0.1666 | 44 |
| GBM | 0.3161 | 148 | 0.3397 | 115 | 0.6575 | 0.295 | 0.1618 | 40 |
| HCT116 | 0.3684 | 133 | 0.3663 | 109 | 0.6693 | 0.3292 | 0.1796 | 37 |
| HELA | 0.3193 | 163 | 0.314 | 138 | 0.67 | 0.2734 | 0.1455 | 36 |
| RPE1 | 0.2598 | 170 | 0.2958 | 133 | 0.654 | 0.2502 | 0.1365 | 40 |
| all |  |  |  |  | 0.664 | 0.2941 | 0.1583 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **110** gene(s) directly affected plus **2340** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- ACAA2: 5/5 lines, knockout now blocks growth -- :x: wrong
- AKR1A1, CYP3A4, GNPAT: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct

**Likely noise** (20 gene(s), <3/5 lines): CA13, CA5B, CRAT, CYP2J2, CYP3A5, ELOVL1, HSD17B7, IDH1, LIPA, MPC1, MPC2, NAGK, PISD, PXMP2, RDH5, SGPL1, SLC17A5, SLC22A5, SLC25A17, SLC27A5.

**No change:** 2426 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/review-v2.0.1-curation/data/testResults/gene-essential-diff.csv).
