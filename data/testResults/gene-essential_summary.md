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

Checked **1** gene(s) directly affected plus **334** more that share a metabolite with one of them.

**No change:** 335 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/992-cyp2e1-nadph/data/testResults/gene-essential-diff.csv).
