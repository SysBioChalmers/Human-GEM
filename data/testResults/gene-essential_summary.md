### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3444 | 142 | 0.3777 | 100 | 0.6651 | 0.3217 | 0.1676 | 51 |
| GBM | 0.3109 | 153 | 0.3328 | 123 | 0.656 | 0.2907 | 0.1615 | 36 |
| HCT116 | 0.3657 | 135 | 0.3662 | 109 | 0.6681 | 0.3296 | 0.1797 | 39 |
| HELA | 0.3068 | 174 | 0.3088 | 142 | 0.6686 | 0.2708 | 0.1457 | 43 |
| RPE1 | 0.2597 | 171 | 0.2982 | 132 | 0.6579 | 0.2537 | 0.1362 | 42 |
| all |  |  |  |  | 0.663 | 0.2932 | 0.1584 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **5** gene(s) directly affected plus **1551** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- UGP2: 5/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (15 gene(s), <3/5 lines): AKR1A1, ALDH4A1, ASRGL1, FPGS, GNPAT, HSD17B7, MPC1, MPC2, NAGK, PXMP2, RDH5, SLC17A5, SLC25A2, SLC27A5, SLC36A1.

**No change:** 1540 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/galactose-1-phosphate-duplicates/data/testResults/gene-essential-diff.csv).
