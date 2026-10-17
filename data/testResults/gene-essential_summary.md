### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3508 | 140 | 0.3751 | 102 | 0.6675 | 0.3199 | 0.1668 | 48 |
| GBM | 0.322 | 141 | 0.3439 | 110 | 0.6588 | 0.2976 | 0.1612 | 38 |
| HCT116 | 0.3729 | 130 | 0.3653 | 110 | 0.6702 | 0.3288 | 0.1789 | 33 |
| HELA | 0.3326 | 152 | 0.324 | 125 | 0.67 | 0.2812 | 0.145 | 40 |
| RPE1 | 0.2649 | 166 | 0.2951 | 134 | 0.6517 | 0.2468 | 0.1358 | 35 |
| all |  |  |  |  | 0.6639 | 0.2952 | 0.1578 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **1** gene(s) directly affected plus **1623** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- SLC17A5: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct
- SLC27A5: 3/5 lines, knockout no longer blocks growth -- :warning: mixed across lines

**Other role changed** (not Hart-comparable):
- IDH1: 3/5 lines, biosynthesis, substrate utilization -- no longer required
- SGPL1: 3/5 lines, biosynthesis -- no longer required

**Likely noise** (14 gene(s), <3/5 lines): ACY1, AKR1A1, ASRGL1, CA13, FPGS, GGH, GNPAT, HIBCH, LIPA, MPC1, MPC2, NAGK, SLC25A2, SLC36A1.

**No change:** 1606 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
