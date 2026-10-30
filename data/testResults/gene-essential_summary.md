### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1168, 2026-10-20.

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3478 | 134 | 0.3756 | 99 | 0.6656 | 0.3219 | 0.1694 | 43 |
| GBM | 0.3214 | 142 | 0.3425 | 114 | 0.6556 | 0.2959 | 0.1631 | 34 |
| HCT116 | 0.3707 | 127 | 0.3619 | 110 | 0.6689 | 0.3293 | 0.1815 | 29 |
| HELA | 0.3214 | 157 | 0.3198 | 135 | 0.6688 | 0.2765 | 0.1458 | 31 |
| RPE1 | 0.2647 | 160 | 0.2936 | 133 | 0.6506 | 0.2474 | 0.1371 | 29 |
| all |  |  |  |  | 0.662 | 0.2942 | 0.1596 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **53** gene(s) directly affected plus **2230** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- AKR1A1, GNPAT, SLC22A5: 3/5 lines, knockout no longer blocks growth -- :warning: mixed across lines
- COQ4 (added to the model): 5/5 lines, knockout now blocks growth -- :warning: mixed across lines
- CRAT: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct
- ELOVL1: 4/5 lines, knockout no longer blocks growth -- :sparkles: correct
- FADS2, H6PD, MMAB, PRUNE1: 5/5 lines, knockout now blocks growth -- :x: wrong
- HPDL (added to the model): 5/5 lines, knockout now blocks growth -- :x: wrong
- MMAA (removed from the model): 5/5 lines, knockout no longer blocks growth -- :sparkles: correct

**Other role changed** (not Hart-comparable):
- ACAA1, CPT1A, CPT2, EHHADH: 5/5 lines, biosynthesis, internal conversion, substrate utilization -- no longer required
- DECR2, IDH1: 3/5 lines, biosynthesis, substrate utilization -- now required
- ECI2: 4/5 lines, biosynthesis, internal conversion, substrate utilization -- no longer required
- HADH, HSD17B10: 5/5 lines, biosynthesis, substrate utilization -- no longer required
- MYO5B (removed from the model): 5/5 lines, biosynthesis, substrate utilization -- no longer required

**Likely noise** (20 gene(s), <3/5 lines): ACAT1, CA5B, FPGS, GGH, GUK1, LIPA, MPC1, MPC2, NAGK, PC, PISD, PXMP2, RDH5, SLC17A5, SLC25A17, SLC25A20, SLC27A5, SLC36A1, SLC37A4, SOAT1.

**No change:** 2241 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
