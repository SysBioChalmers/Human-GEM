### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1177, 2026-10-21.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2654 | 0.3428 | 135 | 0.3604 | 101 | 0.6453 | 0.2662 | 0.1266 | 43 |
| GBM | 2654 | 0.2926 | 141 | 0.3128 | 110 | 0.6255 | 0.2436 | 0.1319 | 37 |
| HCT116 | 2714 | 0.3527 | 130 | 0.3453 | 112 | 0.645 | 0.2728 | 0.14 | 30 |
| HELA | 2716 | 0.3128 | 162 | 0.3062 | 138 | 0.647 | 0.2264 | 0.1105 | 35 |
| RPE1 | 2654 | 0.2284 | 173 | 0.2563 | 139 | 0.6157 | 0.1942 | 0.1093 | 36 |
| all | 2651 | 0.3036 | 172 | 0.3241 | 145 | 0.6359 | 0.2401 | 0.1237 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **305** gene(s) directly affected plus **2298** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- ACAA2: 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- COQ4 (added to the model): 5/5 lines, knockout now blocks growth -- :warning: mixed across lines
- COQ5, HMGCS1: 5/5 lines, knockout now blocks growth -- :warning: mixed across lines
- CRAT, ELOVL1, SLC22A5: 4/5 lines, knockout no longer blocks growth -- :sparkles: correct
- FADS2, H6PD, MMAB, MTRR, PRUNE1: 5/5 lines, knockout now blocks growth -- :x: wrong
- HPDL (added to the model): 5/5 lines, knockout now blocks growth -- :x: wrong
- MMAA (removed from the model): 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- SLC37A4 (added to the model): 5/5 lines, knockout now blocks growth -- :information_source: not scored by Hart 2015
- SLC37A4 (removed from the model): 3/5 lines, knockout no longer blocks growth -- :information_source: not scored by Hart 2015

**Other role changed** (not Hart-comparable):
- ACAA1, CPT1A, CPT2, EHHADH: 5/5 lines, biosynthesis, internal conversion, substrate utilization -- no longer required
- DECR2, IDH1: 3/5 lines, biosynthesis, substrate utilization -- now required
- ECI2: 4/5 lines, biosynthesis, internal conversion, substrate utilization -- no longer required
- HADH, HSD17B10: 5/5 lines, biosynthesis, substrate utilization -- no longer required
- MYO5B (removed from the model): 5/5 lines, biosynthesis, substrate utilization -- no longer required

**Likely noise** (27 gene(s), <3/5 lines): AKR1A1, ALDH1L2, CA5B, CP, DMGDH, FPGS, GGH, GNPAT, GUK1, LIPA, MPC1, MPC2, MTHFD1, NAGK, PC, PDPR, PISD, PXMP2, RPIA, SGPL1, SLC17A5, SLC17A7, SLC25A17, SLC25A20, SLC27A5, SLC36A1, SOAT1.

**No change:** 2550 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
