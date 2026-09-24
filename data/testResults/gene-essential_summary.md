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

Checked **2848** gene(s): every gene in the matrix, because this pull request changes the algorithm, not the model.

**Growth effect changed** (vs Hart 2015):
- ABCD1, IPMK: 4/5 lines, knockout no longer blocks growth -- :sparkles: correct
- ACAA1, COQ5: 4/5 lines, knockout no longer blocks growth -- :warning: mixed across lines
- CRAT, ELOVL1: 5/5 lines, knockout now blocks growth -- :x: wrong
- LCAT, PTS: 5/5 lines, knockout no longer blocks growth -- :sparkles: correct
- PGS1: 3/5 lines, knockout now blocks growth -- :warning: mixed across lines
- PTPMT1: 3/5 lines, knockout now blocks growth -- :sparkles: correct
- SLC22A5: 5/5 lines, knockout now blocks growth -- :warning: mixed across lines

**Other role changed** (not Hart-comparable):
- ACADVL: 4/5 lines, substrate utilization -- now required
- ALDH4A1: 4/5 lines, biosynthesis, growth (biomass production), substrate utilization -- now required
- ECH1: 3/5 lines, biosynthesis, substrate utilization -- now required
- GAPDH, GPI, PGK1: 3/5 lines, substrate utilization -- now required
- HOGA1: 3/5 lines, biosynthesis, growth (biomass production), substrate utilization -- now required

**Likely noise** (102 gene(s), <3/5 lines): ABAT, ACAD9, ACAT2, ACOT7, ACOX1, ACY1, ACY3, ADH5, ADI1, AFMID, AGXT, AKR1A1, AKR1C3, ALDH3A2, ALDH6A1, ALDH7A1, ALDOC, ALG5, AMD1, AMT, APIP, ASRGL1, BCAT2, BHMT2, CA5B, CHDH, CLYBL, CPT1A, CPT2, CYP4F2, CYP4F3, DCK, DECR1, DGKZ, DUT, ECHS1, ENOPH1, ESD, ETFA, ETFB, ETFDH, FPGS, GCSH, GGH, GLDC, GNPAT, GOT1, GOT2, GPT2, GUK1, H6PD, HAAO, HIBCH, HIF1A, HMGCS1, HPRT1, HSD17B4, HSD3B1, IDH1, KYNU, LDHB, MLYCD, MPC1, MPC2, MRI1, MTAP, MTHFD1, MTHFR, MTR, NADSYN1, NAGS, NAT8L, NMNAT3, PC, PCYT1A, PDPR, PEMT, PISD, PNP, PRDX6, PSPH, QPRT, RPIA, SARDH, SDSL, SGPL1, SHMT2, SLC16A1, SLC25A15, SLC25A17, SLC25A2, SLC25A20, SLC25A36, SLC25A38, SLC29A1, SLC35D1, SLC36A1, SLC37A4, SMOX, TALDO1, TDO2, XDH.

**No change:** 2728 gene(s).

Full per-gene, per-line detail (every gene checked, not just the ones named above): [gene-essential-diff.csv](https://github.com/SysBioChalmers/Human-GEM/blob/test/gene-essentiality-gurobi-ties/data/testResults/gene-essential-diff.csv).
