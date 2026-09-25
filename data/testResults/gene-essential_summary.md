### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3643 | 129 | 0.3952 | 89 | 0.6704 | 0.3369 | 0.1681 | 50 |
| GBM | 0.3046 | 158 | 0.3346 | 119 | 0.6551 | 0.2917 | 0.1619 | 47 |
| HCT116 | 0.3747 | 130 | 0.3762 | 102 | 0.6716 | 0.3389 | 0.1805 | 42 |
| HELA | 0.3109 | 168 | 0.3121 | 140 | 0.6688 | 0.2731 | 0.1456 | 38 |
| RPE1 | 0.2616 | 174 | 0.3101 | 124 | 0.6588 | 0.2596 | 0.1358 | 55 |
| all |  |  |  |  | 0.6648 | 0.2995 | 0.1586 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **1** gene(s) directly affected plus **871** more that share a metabolite with one of them.

**Other role changed** (not Hart-comparable):
- PCYT1A: 3/5 lines, biosynthesis -- no longer required

**Likely noise** (17 gene(s), <3/5 lines): AGXT, ALDH7A1, ALDOC, FPGS, GAPDH, GNPAT, GOT2, GPT2, GRHPR, LDHB, PCYT2, PSPH, SLC25A15, SLC25A2, SLC25A38, SLC26A6, SLC36A1.

**No change:** 854 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/feat/pxmp2-peroxisomal-pore/data/testResults/gene-essential-diff.csv).
