### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1162, 2026-10-18.

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3492 | 137 | 0.3806 | 97 | 0.6669 | 0.3225 | 0.1665 | 49 |
| GBM | 0.3194 | 141 | 0.3458 | 107 | 0.6589 | 0.2989 | 0.1611 | 41 |
| HCT116 | 0.3788 | 124 | 0.3749 | 102 | 0.6708 | 0.3348 | 0.1788 | 35 |
| HELA | 0.3093 | 167 | 0.3031 | 142 | 0.6614 | 0.2631 | 0.1447 | 36 |
| RPE1 | 0.2608 | 167 | 0.2946 | 132 | 0.65 | 0.2464 | 0.1356 | 38 |
| all |  |  |  |  | 0.6618 | 0.2928 | 0.1576 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **5** gene(s) directly affected plus **1590** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- SEC11A (removed from the model), SEC11C (removed from the model), SPCS1 (removed from the model): 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- SPCS2 (removed from the model): 5/5 lines, knockout no longer blocks growth -- :information_source: not scored by Hart 2015
- SPCS3 (removed from the model): 5/5 lines, knockout no longer blocks growth -- :x: wrong

**Other role changed** (not Hart-comparable):
- SGPL1: 3/5 lines, biosynthesis -- now required

**Likely noise** (16 gene(s), <3/5 lines): ACY1, AKR1A1, ASRGL1, FPGS, GGH, GNPAT, HIBCH, LIPA, MPC1, MPC2, MTHFD1, NAGK, SLC17A5, SLC25A2, SLC27A5, SLC36A1.

**No change:** 1573 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
