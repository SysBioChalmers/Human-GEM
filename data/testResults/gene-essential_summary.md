### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3458 | 138 | 0.3694 | 103 | 0.6654 | 0.3187 | 0.1693 | 44 |
| GBM | 0.3177 | 145 | 0.345 | 110 | 0.6599 | 0.3002 | 0.163 | 42 |
| HCT116 | 0.3723 | 128 | 0.3679 | 106 | 0.6717 | 0.3347 | 0.1814 | 35 |
| HELA | 0.3078 | 172 | 0.3013 | 147 | 0.6641 | 0.2647 | 0.1457 | 36 |
| RPE1 | 0.2487 | 181 | 0.2865 | 141 | 0.6512 | 0.244 | 0.1373 | 43 |
| all |  |  |  |  | 0.6625 | 0.292 | 0.1596 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **23** gene(s) directly affected plus **1612** more that share a metabolite with one of them.

**Likely noise** (10 gene(s), <3/5 lines): AKR1A1, ELOVL1, GNPAT, LIPA, NAGK, PISD, RDH5, SGPL1, SLC17A5, SLC27A5.

**No change:** 1625 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](gene-essential-diff.csv).
