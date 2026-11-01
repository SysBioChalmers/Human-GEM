### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1183, 2026-11-01.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2651 | 0.3405 | 136 | 0.3609 | 104 | 0.6456 | 0.2648 | 0.126 | 39 |
| GBM | 2651 | 0.2894 | 142 | 0.3094 | 111 | 0.6238 | 0.2412 | 0.1316 | 37 |
| HCT116 | 2711 | 0.3541 | 130 | 0.3479 | 111 | 0.6464 | 0.2747 | 0.1394 | 31 |
| HELA | 2713 | 0.3117 | 164 | 0.3038 | 141 | 0.6473 | 0.2252 | 0.1102 | 34 |
| RPE1 | 2651 | 0.2348 | 168 | 0.2635 | 134 | 0.6215 | 0.1993 | 0.1083 | 36 |
| all | 2648 | 0.3052 | 172 | 0.3258 | 145 | 0.6371 | 0.2408 | 0.1231 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **12** gene(s) directly affected plus **2176** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- CRAT: 4/5 lines, knockout now blocks growth -- :x: wrong
- FADS1: 5/5 lines, knockout now blocks growth -- :x: wrong

**Other role changed** (not Hart-comparable):
- CPT2: 5/5 lines, biosynthesis, substrate utilization -- now required

**Likely noise** (12 gene(s), <3/5 lines): CA13, CA5B, FPGS, GGH, MPC1, MPC2, PC, PISD, RDH5, SLC25A17, SLC25A20, SLC36A1.

**No change:** 2173 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
