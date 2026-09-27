### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3417 | 144 | 0.3792 | 99 | 0.6703 | 0.3266 | 0.1677 | 54 |
| GBM | 0.3203 | 140 | 0.3566 | 102 | 0.6588 | 0.3057 | 0.1616 | 44 |
| HCT116 | 0.3695 | 132 | 0.372 | 105 | 0.6698 | 0.3342 | 0.1799 | 40 |
| HELA | 0.3123 | 169 | 0.3112 | 140 | 0.669 | 0.2728 | 0.1458 | 40 |
| RPE1 | 0.253 | 177 | 0.3019 | 129 | 0.6556 | 0.2542 | 0.1365 | 51 |
| all |  |  |  |  | 0.6647 | 0.2987 | 0.1586 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **9** gene(s) directly affected plus **1986** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- HSD17B7: 3/5 lines, knockout now blocks growth -- :information_source: not scored by Hart 2015

**Likely noise** (11 gene(s), <3/5 lines): AKR1A1, GNPAT, LIPA, MPC1, MPC2, NAGK, PISD, PXMP2, RDH5, SLC17A5, SLC27A5.

**No change:** 1983 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/880-compartment-check/data/testResults/gene-essential-diff.csv).
