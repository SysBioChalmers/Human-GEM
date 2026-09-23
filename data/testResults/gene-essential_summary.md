### Graded gene essentiality vs Hart 2015 (task-scope analysis)

| cellLine | allTaskMCC | allTaskFP | viabilityMCC | viabilityFP | growthAUROC | growthAUPRC | baseRate | capabilityOnly |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DLD1 | 0.3464 | 155 | 0.3885 | 103 | 0.6907 | 0.3346 | 0.1636 | 63 |
| GBM | 0.3261 | 154 | 0.3625 | 112 | 0.6815 | 0.3201 | 0.1612 | 50 |
| HCT116 | 0.3772 | 131 | 0.3925 | 96 | 0.6811 | 0.3546 | 0.1811 | 48 |
| HELA | 0.3304 | 164 | 0.3456 | 120 | 0.6833 | 0.2936 | 0.1434 | 56 |
| RPE1 | 0.2715 | 184 | 0.331 | 126 | 0.6766 | 0.272 | 0.1353 | 63 |
| all |  |  |  |  | 0.6821 | 0.3156 | 0.1577 |  |

### Gene essentiality: effect of this change (vs `develop`)

Checked **9** gene(s) directly affected plus **425** more that share a metabolite with one of them.

**Growth effect changed** (vs Hart 2015):
- CRAT: 3/5 lines, knockout no longer blocks growth -- :sparkles: correct

**Other role changed** (not Hart-comparable):
- ECH1: 3/5 lines, biosynthesis, substrate utilization -- no longer required
- ETFA, ETFB, ETFDH: 5/5 lines, biosynthesis, substrate utilization -- now required

**Likely noise** (32 gene(s), <3/5 lines): ACAA1, ACAD9, ACADVL, ACOT2, ACOX1, ADH5, AKR1A1, ALDH3A2, ALDH4A1, ALDH6A1, ALDH7A1, BCKDHA, BCKDHB, CA5B, CHDH, COQ5, CPT2, CYP1B1, DECR1, GLYCTK, GNPAT, GOT2, GRHPR, HSD17B4, MLYCD, NAT8L, PEMT, PGS1, PISD, SARDH, SLC25A1, SLC27A2.

**No change:** 397 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/fix/etf-ubiquinone-coupling/data/testResults/gene-essential-diff.csv).
