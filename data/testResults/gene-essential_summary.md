### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3393 | 149 | 0.3644 | 109 | 0.6626 | 0.311 | 0.1669 | 50 |
| GBM | 0.3152 | 149 | 0.3386 | 116 | 0.6549 | 0.2925 | 0.1613 | 40 |
| HCT116 | 0.3701 | 132 | 0.3666 | 109 | 0.6671 | 0.327 | 0.1791 | 36 |
| HELA | 0.3075 | 174 | 0.3046 | 146 | 0.6651 | 0.2657 | 0.1449 | 39 |
| RPE1 | 0.2592 | 171 | 0.2963 | 133 | 0.656 | 0.25 | 0.136 | 41 |
| all |  |  |  |  | 0.661 | 0.2893 | 0.1579 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **16** gene(s) directly affected plus **1887** more that share a metabolite with one of them.

**Other role changed** (not Hart-comparable):
- PXMP2: 3/5 lines, substrate utilization -- now required

**Likely noise** (11 gene(s), <3/5 lines): AKR1A1, GNPAT, LIPA, MPC1, MPC2, NAGK, PISD, RDH5, SLC17A5, SLC25A17, SLC27A5.

**No change:** 1891 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
