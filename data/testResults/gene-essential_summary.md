### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3456 | 144 | 0.3659 | 108 | 0.6633 | 0.3121 | 0.1669 | 46 |
| GBM | 0.3172 | 142 | 0.3423 | 111 | 0.6561 | 0.2952 | 0.1614 | 37 |
| HCT116 | 0.3701 | 132 | 0.3637 | 111 | 0.6689 | 0.3269 | 0.1791 | 34 |
| HELA | 0.3077 | 177 | 0.3021 | 148 | 0.6669 | 0.2652 | 0.1449 | 41 |
| RPE1 | 0.2549 | 175 | 0.2937 | 135 | 0.6506 | 0.246 | 0.136 | 43 |
| all |  |  |  |  | 0.6612 | 0.2891 | 0.1579 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **22** gene(s) directly affected plus **1834** more that share a metabolite with one of them.

**Likely noise** (14 gene(s), <3/5 lines): AKR1A1, CA13, CA5B, CYP2R1, GNPAT, IDH1, MPC1, MPC2, PC, RDH5, SLC17A5, SLC25A17, SLC25A20, SLC27A5.

**No change:** 1842 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
