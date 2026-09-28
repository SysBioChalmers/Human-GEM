### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3368 | 148 | 0.3654 | 108 | 0.6649 | 0.3137 | 0.1676 | 50 |
| GBM | 0.3166 | 143 | 0.3415 | 112 | 0.6564 | 0.295 | 0.1615 | 37 |
| HCT116 | 0.3683 | 133 | 0.3676 | 108 | 0.6708 | 0.3324 | 0.1798 | 38 |
| HELA | 0.3145 | 167 | 0.31 | 141 | 0.6696 | 0.272 | 0.1458 | 37 |
| RPE1 | 0.2575 | 173 | 0.2968 | 133 | 0.6553 | 0.2506 | 0.1363 | 43 |
| all |  |  |  |  | 0.6634 | 0.293 | 0.1585 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **35** gene(s) directly affected plus **729** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- ACADSB, AMACR, CYP4F2, HADHA, HADHB: 5/5 lines, knockout now blocks growth -- :x: wrong
- CYP4F11: 4/5 lines, knockout no longer blocks growth -- :sparkles: correct
- ECHS1: 4/5 lines, knockout now blocks growth -- :x: wrong
- ETFA, ETFB, ETFDH: 3/5 lines, knockout now blocks growth -- :x: wrong
- HSD17B7: 3/5 lines, knockout no longer blocks growth -- :information_source: not scored by Hart 2015

**Likely noise** (7 gene(s), <3/5 lines): ACADL, AKR1A1, LIPA, MPC1, MPC2, SGPL1, SLC17A7.

**No change:** 746 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/738-888-vitamin-e-beta-oxidation/data/testResults/gene-essential-diff.csv).
