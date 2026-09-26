### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3427 | 143 | 0.3821 | 97 | 0.6686 | 0.3266 | 0.1681 | 55 |
| GBM | 0.3231 | 135 | 0.3624 | 96 | 0.6606 | 0.3118 | 0.1619 | 45 |
| HCT116 | 0.3719 | 130 | 0.3748 | 103 | 0.6707 | 0.3362 | 0.1803 | 40 |
| HELA | 0.3109 | 170 | 0.3135 | 138 | 0.6684 | 0.2745 | 0.146 | 43 |
| RPE1 | 0.263 | 169 | 0.3128 | 122 | 0.6644 | 0.2641 | 0.1358 | 50 |
| all |  |  |  |  | 0.6663 | 0.3025 | 0.1587 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **70** gene(s) directly affected plus **2319** more that share a metabolite with one of them.

**Other role changed** (not Hart-comparable):
- ACAD9: 4/5 lines, biosynthesis, substrate utilization -- no longer required
- DECR2, IDH1: 3/5 lines, biosynthesis, substrate utilization -- now required
- ECH1: 5/5 lines, biosynthesis, substrate utilization -- no longer required

**Likely noise** (27 gene(s), <3/5 lines): AKR1A1, CA13, CA5B, CLYBL, CRAT, CROT, CYP2J2, ELOVL1, GNPAT, H6PD, HSD17B7, LIPA, MPC1, MPC2, NAGK, PC, PGK1, PISD, PXMP2, SGPL1, SLC17A5, SLC22A5, SLC25A17, SLC25A20, SLC27A5, SLC37A4, SOAT1.

**No change:** 2358 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/duplicate-metabolite-ids/data/testResults/gene-essential-diff.csv).
