### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1177, 2026-10-30.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2654 | 0.339 | 136 | 0.3604 | 101 | 0.6439 | 0.2648 | 0.1266 | 43 |
| GBM | 2654 | 0.2899 | 141 | 0.3111 | 109 | 0.6238 | 0.2425 | 0.1319 | 38 |
| HCT116 | 2714 | 0.3546 | 126 | 0.3514 | 107 | 0.645 | 0.2761 | 0.14 | 30 |
| HELA | 2716 | 0.311 | 164 | 0.3031 | 141 | 0.6465 | 0.225 | 0.1105 | 34 |
| RPE1 | 2654 | 0.2309 | 170 | 0.2603 | 135 | 0.62 | 0.1984 | 0.1093 | 37 |
| all | 2651 | 0.3054 | 170 | 0.3252 | 144 | 0.6359 | 0.241 | 0.1237 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **304** gene(s) directly affected plus **2284** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- ACAA2: 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- COQ5, HMGCS1: 5/5 lines, knockout now blocks growth -- :warning: mixed across lines
- MTRR: 5/5 lines, knockout now blocks growth -- :x: wrong
- SLC25A20: 3/5 lines, knockout now blocks growth -- :x: wrong
- SLC37A4 (added to the model): 5/5 lines, knockout now blocks growth -- :information_source: not scored by Hart 2015
- SLC37A4 (removed from the model): 5/5 lines, knockout no longer blocks growth -- :information_source: not scored by Hart 2015

**Likely noise** (24 gene(s), <3/5 lines): ACAT1, AKR1A1, ALDH1L2, CA5B, CRAT, DMGDH, FPGS, GGH, GNPAT, KL, MPC1, MPC2, PC, PDPR, RDH5, RPIA, SGPL1, SLC17A5, SLC17A7, SLC22A5, SLC25A17, SLC27A5, SLC36A1, SOAT1.

**No change:** 2557 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
