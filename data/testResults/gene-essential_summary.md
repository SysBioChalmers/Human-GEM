### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3479 | 139 | 0.3806 | 98 | 0.667 | 0.3255 | 0.168 | 50 |
| GBM | 0.3254 | 141 | 0.3567 | 104 | 0.6608 | 0.3074 | 0.162 | 44 |
| HCT116 | 0.3732 | 129 | 0.3763 | 102 | 0.6711 | 0.3374 | 0.1803 | 40 |
| HELA | 0.3142 | 164 | 0.3148 | 137 | 0.6679 | 0.2744 | 0.146 | 37 |
| RPE1 | 0.2556 | 176 | 0.3034 | 129 | 0.6569 | 0.2543 | 0.1357 | 50 |
| all |  |  |  |  | 0.6647 | 0.2997 | 0.1587 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **6** gene(s) directly affected plus **157** more that share a metabolite with one of them.

**Likely noise** (1 gene(s), <3/5 lines): SLC22A5.

**No change:** 162 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/transport-stoichiometry/data/testResults/gene-essential-diff.csv).
