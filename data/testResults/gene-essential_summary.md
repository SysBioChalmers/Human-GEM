### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3531 | 138 | 0.3764 | 101 | 0.6651 | 0.319 | 0.1671 | 47 |
| GBM | 0.3183 | 141 | 0.3422 | 111 | 0.6572 | 0.2956 | 0.1616 | 36 |
| HCT116 | 0.3659 | 135 | 0.3634 | 111 | 0.6685 | 0.327 | 0.1795 | 37 |
| HELA | 0.3116 | 173 | 0.3029 | 147 | 0.6634 | 0.2643 | 0.1454 | 38 |
| RPE1 | 0.2677 | 163 | 0.2986 | 131 | 0.6564 | 0.252 | 0.1363 | 35 |
| all |  |  |  |  | 0.6621 | 0.2915 | 0.1582 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **83** gene(s) directly affected plus **2385** more that share a metabolite with one of them.

**Other role changed** (not Hart-comparable):
- PCYT1A: 3/5 lines, biosynthesis -- now required

**Likely noise** (25 gene(s), <3/5 lines): AKR1A1, ALDH4A1, ATP8A1, CA5B, CRAT, ELOVL1, FPGS, GGH, GNPAT, HSD17B7, IDH1, LIPA, MPC1, MPC2, NAGK, PC, PCYT2, PLD1, PXMP2, SGPL1, SLC17A5, SLC22A5, SLC27A5, SLC36A1, SLC37A4.

**No change:** 2442 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/remaining-reactions/data/testResults/gene-essential-diff.csv).
