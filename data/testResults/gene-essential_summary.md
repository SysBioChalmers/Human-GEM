### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3456 | 141 | 0.3745 | 102 | 0.667 | 0.3203 | 0.1676 | 48 |
| GBM | 0.3192 | 141 | 0.3445 | 110 | 0.6562 | 0.2968 | 0.1614 | 37 |
| HCT116 | 0.363 | 137 | 0.3588 | 114 | 0.6684 | 0.3253 | 0.1799 | 36 |
| HELA | 0.3101 | 171 | 0.3051 | 145 | 0.6701 | 0.2704 | 0.1458 | 37 |
| RPE1 | 0.251 | 179 | 0.289 | 139 | 0.6553 | 0.2476 | 0.1364 | 43 |
| all |  |  |  |  | 0.6632 | 0.2919 | 0.1585 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **33** gene(s) directly affected plus **2158** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- HSD17B7: 3/5 lines, knockout now blocks growth -- :information_source: not scored by Hart 2015

**Other role changed** (not Hart-comparable):
- SLC25A17: 3/5 lines, biosynthesis, internal conversion, substrate utilization -- no longer required

**Likely noise** (11 gene(s), <3/5 lines): AKR1A1, GNPAT, LIPA, MPC1, MPC2, NAGK, RDH5, SGPL1, SLC17A5, SLC17A7, SLC27A5.

**No change:** 2178 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/1085-metabolite-curation/data/testResults/gene-essential-diff.csv).
