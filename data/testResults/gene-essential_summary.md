### Graded gene essentiality vs Hart 2015 (task-scope analysis)

Run: pull request #1181, 2026-11-01.

Scored over every model gene that Hart 2015 scored in the cell line (`genes`); a gene outside the cell line's context model counts as non-essential, with growth ratio 1.0. The `all` row: essential in every cell line against a fitness gene in every cell line, over the genes Hart scored in every cell line (Human1 used every Hart gene).

| cellLine | genes | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 2654 | 0.339 | 136 | 0.3591 | 102 | 0.6452 | 0.2659 | 0.1266 | 42 |
| GBM | 2654 | 0.2947 | 139 | 0.314 | 109 | 0.6256 | 0.2444 | 0.1319 | 36 |
| HCT116 | 2714 | 0.3571 | 126 | 0.3514 | 107 | 0.645 | 0.2761 | 0.14 | 31 |
| HELA | 2716 | 0.311 | 164 | 0.3021 | 139 | 0.6453 | 0.2244 | 0.1105 | 37 |
| RPE1 | 2654 | 0.2301 | 171 | 0.2593 | 136 | 0.6197 | 0.1978 | 0.1093 | 37 |
| all | 2651 | 0.3027 | 173 | 0.3231 | 146 | 0.6363 | 0.2414 | 0.1237 |  |

### Run environment

- cobra 0.32.1, gurobipy 13.0.3, optlang 1.9.1, python 3.11.17, raven-toolbox 3.0.0b1, raven-toolbox commit ee0f8868bcc4a6f907913151a1e54dcbb4f2a160

ftINIT steps that kept an unproven incumbent (time limit or unconverged tie-break):
- RPE1 (reaction order 0): 1 warning(s)

### Gene essentiality: effect of this change (vs `develop`)

Checked **7** gene(s) directly affected plus **51** more that share a metabolite with one of them.

**Likely noise** (1 gene(s), <3/5 lines): RDH5.

**No change:** 57 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
