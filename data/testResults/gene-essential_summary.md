### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1184, 2026-11-01.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2652 | 0.3427 | 134 | 0.3632 | 100 | 0.644 | 0.2653 | 0.1259 | 42 |
| GBM | 2652 | 0.2905 | 141 | 0.3118 | 109 | 0.6235 | 0.2419 | 0.1316 | 38 |
| HCT116 | 2712 | 0.3585 | 126 | 0.3527 | 107 | 0.6451 | 0.2758 | 0.1394 | 31 |
| HELA | 2714 | 0.3118 | 164 | 0.3038 | 138 | 0.6433 | 0.2234 | 0.1102 | 38 |
| RPE1 | 2652 | 0.2315 | 171 | 0.2607 | 136 | 0.6183 | 0.1963 | 0.1086 | 37 |
| all | 2649 | 0.307 | 170 | 0.3269 | 144 | 0.635 | 0.2401 | 0.1232 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- HELA (reaction order 0): 1 warning(s)
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **50** gene(s) directly affected plus **1839** more that share a metabolite with one of them.

**Likely noise** (5 gene(s), <3/5 lines): CA5B, KL, MTHFD1, PC, RDH5.

**No change:** 1884 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
