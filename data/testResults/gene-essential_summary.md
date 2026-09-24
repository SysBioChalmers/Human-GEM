### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3613 | 134 | 0.3764 | 101 | 0.6712 | 0.3265 | 0.1684 | 45 |
| GBM | 0.2893 | 168 | 0.33 | 122 | 0.6518 | 0.2888 | 0.1624 | 53 |
| HCT116 | 0.3776 | 131 | 0.3775 | 104 | 0.6783 | 0.3448 | 0.1806 | 41 |
| HELA | 0.325 | 166 | 0.3216 | 134 | 0.6685 | 0.2798 | 0.1469 | 45 |
| RPE1 | 0.2528 | 183 | 0.3054 | 122 | 0.6599 | 0.2625 | 0.1357 | 70 |
| all |  |  |  |  | 0.6663 | 0.3009 | 0.1591 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **4** gene(s) directly affected plus **441** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- PRDX6: 3/5 lines, knockout now blocks growth -- :x: wrong

**Likely noise** (26 gene(s), <3/5 lines): AKR1A1, ALDH4A1, ASRGL1, BCKDHA, BCKDHB, CA5B, CLYBL, CYP27A1, DCK, DUT, ETFA, ETFB, ETFDH, FPGS, GOT2, HSD3B1, LDHB, MLYCD, MPC1, MPC2, NAGS, PC, SHMT2, SLC16A1, SLC25A1, SLC25A2.

**No change:** 418 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/sqor-chdh-prodh2-ubiquinone/data/testResults/gene-essential-diff.csv).
