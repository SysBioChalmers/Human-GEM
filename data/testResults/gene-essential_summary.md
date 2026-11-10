### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1189, 2026-11-10.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2651 | 0.3422 | 137 | 0.3619 | 101 | 0.6441 | 0.2644 | 0.126 | 45 |
| GBM | 2651 | 0.2922 | 142 | 0.3122 | 111 | 0.6247 | 0.2422 | 0.1316 | 37 |
| HCT116 | 2711 | 0.3563 | 128 | 0.3503 | 109 | 0.6447 | 0.2741 | 0.1394 | 31 |
| HELA | 2713 | 0.3089 | 167 | 0.3018 | 143 | 0.6455 | 0.2224 | 0.1102 | 35 |
| RPE1 | 2651 | 0.2313 | 172 | 0.2615 | 136 | 0.6187 | 0.1963 | 0.1083 | 38 |
| all | 2648 | 0.3034 | 174 | 0.3248 | 146 | 0.6357 | 0.2394 | 0.1231 |  |

### Run environment

- cobra 0.32.1, gurobipy 14.0.0, optlang 1.9.1, python 3.11.18, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **5** gene(s) directly affected plus **1383** more that share a metabolite with one of them.

**Likely noise** (8 gene(s), <3/5 lines): CA13, FPGS, MPC1, MPC2, MTHFD1, PISD, RDH5, SLC36A1.

**No change:** 1380 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
