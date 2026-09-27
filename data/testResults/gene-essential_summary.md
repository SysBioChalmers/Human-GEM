### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3443 | 142 | 0.3792 | 99 | 0.6641 | 0.3217 | 0.1677 | 52 |
| GBM | 0.3194 | 146 | 0.3554 | 105 | 0.6611 | 0.3059 | 0.1616 | 48 |
| HCT116 | 0.3683 | 133 | 0.3706 | 106 | 0.6697 | 0.3326 | 0.1798 | 40 |
| HELA | 0.3179 | 164 | 0.3176 | 135 | 0.6717 | 0.2774 | 0.1458 | 40 |
| RPE1 | 0.2512 | 179 | 0.2982 | 132 | 0.6579 | 0.2537 | 0.1362 | 50 |
| all |  |  |  |  | 0.6647 | 0.2982 | 0.1585 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **2** gene(s) directly affected plus **405** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- GPT2: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct

**Likely noise** (16 gene(s), <3/5 lines): AKR1A1, ALDH4A1, ASRGL1, BCKDHA, BCKDHB, CA5B, DECR2, FPGS, GOT2, HOGA1, MPC1, MPC2, NAGS, PC, SHMT2, SLC25A2.

**No change:** 390 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/418-modifymodel-curations/data/testResults/gene-essential-diff.csv).
