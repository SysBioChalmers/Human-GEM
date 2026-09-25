### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3819 | 126 | 0.3959 | 91 | 0.676 | 0.3404 | 0.1682 | 50 |
| GBM | 0.2938 | 170 | 0.3309 | 124 | 0.6564 | 0.2912 | 0.1622 | 55 |
| HCT116 | 0.3794 | 137 | 0.3758 | 106 | 0.6698 | 0.3367 | 0.181 | 49 |
| HELA | 0.3331 | 165 | 0.3292 | 131 | 0.6732 | 0.285 | 0.1467 | 48 |
| RPE1 | 0.2393 | 199 | 0.2887 | 140 | 0.6548 | 0.2482 | 0.1363 | 67 |
| all |  |  |  |  | 0.666 | 0.2998 | 0.1592 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **2** gene(s) directly affected plus **1311** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- DCK: 3/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (23 gene(s), <3/5 lines): ACOT2, AKR1A1, ALDH4A1, ASRGL1, CYP27A1, DGKZ, DUT, FPGS, GNPAT, GRHPR, HSD17B7, IPMK, LDHB, MLYCD, MPC1, MPC2, MTHFD1, PRDX6, SLC25A1, SLC25A15, SLC25A2, SLC27A2, SLC36A1.

**No change:** 1289 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/feat/succinate-export/data/testResults/gene-essential-diff.csv).
