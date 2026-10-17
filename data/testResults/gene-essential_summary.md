### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3409 | 146 | 0.3681 | 105 | 0.6654 | 0.3147 | 0.1666 | 51 |
| GBM | 0.3181 | 142 | 0.3383 | 112 | 0.6536 | 0.2923 | 0.1611 | 37 |
| HCT116 | 0.3642 | 135 | 0.3587 | 113 | 0.6646 | 0.3216 | 0.1788 | 35 |
| HELA | 0.3326 | 147 | 0.3267 | 121 | 0.6663 | 0.281 | 0.1447 | 38 |
| RPE1 | 0.2357 | 191 | 0.2688 | 153 | 0.6405 | 0.2301 | 0.1355 | 41 |
| all |  |  |  |  | 0.6584 | 0.2874 | 0.1576 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **6** gene(s) directly affected plus **1614** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- NAGK: 3/5 lines, knockout now blocks growth -- :x: wrong
- SLC17A5: 4/5 lines, knockout now blocks growth -- :x: wrong
- SLC27A5: 4/5 lines, knockout now blocks growth -- :warning: mixed across lines

**Likely noise** (16 gene(s), <3/5 lines): AKR1A1, ASRGL1, FPGS, GGH, GNPAT, GRHPR, HAO1, LIPA, MPC1, MPC2, RDH5, SGPL1, SLC25A15, SLC25A2, SLC26A6, SLC36A1.

**No change:** 1601 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
