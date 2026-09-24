### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.38 | 125 | 0.3958 | 89 | 0.6796 | 0.3449 | 0.1684 | 51 |
| GBM | 0.3073 | 161 | 0.3579 | 103 | 0.66 | 0.3112 | 0.1624 | 70 |
| HCT116 | 0.3687 | 137 | 0.3706 | 108 | 0.6743 | 0.3386 | 0.181 | 43 |
| HELA | 0.3238 | 164 | 0.3203 | 135 | 0.6665 | 0.2778 | 0.1469 | 41 |
| RPE1 | 0.2496 | 186 | 0.304 | 126 | 0.6585 | 0.2584 | 0.1358 | 68 |
| all |  |  |  |  | 0.668 | 0.3059 | 0.1592 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **0** gene(s) directly affected plus **0** more that share a metabolite with one of them.

**No change:** 0 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/gene-essential-diff.csv).
