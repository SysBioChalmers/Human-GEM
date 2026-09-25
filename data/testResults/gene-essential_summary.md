### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3822 | 125 | 0.3982 | 89 | 0.6779 | 0.343 | 0.1686 | 51 |
| GBM | 0.3046 | 155 | 0.3379 | 119 | 0.6566 | 0.2939 | 0.1622 | 41 |
| HCT116 | 0.3772 | 130 | 0.3774 | 103 | 0.6777 | 0.3435 | 0.1809 | 41 |
| HELA | 0.3174 | 174 | 0.3152 | 140 | 0.6686 | 0.2744 | 0.1462 | 47 |
| RPE1 | 0.2384 | 196 | 0.2825 | 145 | 0.6496 | 0.2415 | 0.1363 | 56 |
| all |  |  |  |  | 0.6663 | 0.298 | 0.1591 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **18** gene(s) directly affected plus **404** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- CA5B: 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- DCK: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct
- MPC1, MPC2: 3/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (6 gene(s), <3/5 lines): AKR1A1, ALDH6A1, CLYBL, DUT, PC, SLC25A15.

**No change:** 412 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/feat/spontaneous-diffusion/data/testResults/gene-essential-diff.csv).
