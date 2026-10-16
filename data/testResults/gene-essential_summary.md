### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3337 | 151 | 0.3615 | 111 | 0.6645 | 0.3101 | 0.1668 | 49 |
| GBM | 0.3156 | 146 | 0.3365 | 115 | 0.654 | 0.2912 | 0.1613 | 38 |
| HCT116 | 0.3728 | 130 | 0.3667 | 109 | 0.6694 | 0.3289 | 0.179 | 34 |
| HELA | 0.3123 | 174 | 0.3058 | 143 | 0.6678 | 0.2673 | 0.1443 | 44 |
| RPE1 | 0.2648 | 166 | 0.2937 | 135 | 0.6502 | 0.246 | 0.1359 | 34 |
| all |  |  |  |  | 0.6613 | 0.2888 | 0.1577 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **66** gene(s) directly affected plus **1903** more that share a metabolite with one of them.

**Likely noise** (15 gene(s), <3/5 lines): CA5B, FPGS, GGH, IDH1, MPC1, MPC2, MTHFD1, PC, PISD, RDH5, SGPL1, SLC17A5, SLC27A5, SLC35D1, SLC36A1.

**No change:** 1954 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
