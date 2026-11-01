### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1184, 2026-11-01.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2652 | 0.3405 | 136 | 0.3619 | 101 | 0.6439 | 0.2644 | 0.1259 | 43 |
| GBM | 2652 | 0.2943 | 140 | 0.3147 | 109 | 0.6247 | 0.2437 | 0.1316 | 37 |
| HCT116 | 2712 | 0.3585 | 126 | 0.3553 | 107 | 0.6465 | 0.2776 | 0.1394 | 30 |
| HELA | 2714 | 0.3118 | 164 | 0.3049 | 140 | 0.6463 | 0.2251 | 0.1102 | 35 |
| RPE1 | 2652 | 0.2341 | 168 | 0.2638 | 133 | 0.6185 | 0.1982 | 0.1086 | 37 |
| all | 2649 | 0.3061 | 171 | 0.3258 | 145 | 0.6362 | 0.2415 | 0.1232 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **51** gene(s) directly affected plus **1828** more that share a metabolite with one of them.

**Other role changed** (not Hart-comparable):
- CA5B: 3/5 lines, biosynthesis, internal conversion -- now required

**Likely noise** (5 gene(s), <3/5 lines): CA13, FPGS, GGH, SLC25A17, SLC36A1.

**No change:** 1873 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
