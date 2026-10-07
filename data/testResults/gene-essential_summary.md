### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3381 | 150 | 0.3644 | 109 | 0.665 | 0.3121 | 0.1669 | 51 |
| GBM | 0.3231 | 140 | 0.3453 | 109 | 0.6589 | 0.2988 | 0.1614 | 38 |
| HCT116 | 0.3672 | 135 | 0.3604 | 114 | 0.6668 | 0.3231 | 0.1786 | 34 |
| HELA | 0.3141 | 168 | 0.3069 | 144 | 0.6655 | 0.2667 | 0.145 | 35 |
| RPE1 | 0.2486 | 181 | 0.2873 | 140 | 0.6496 | 0.2428 | 0.136 | 44 |
| all |  |  |  |  | 0.6613 | 0.2887 | 0.1578 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **31** gene(s) directly affected plus **529** more that share a metabolite with one of them.

**Likely noise** (1 gene(s), <3/5 lines): LIPA.

**No change:** 559 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
