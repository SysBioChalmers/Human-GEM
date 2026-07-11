# Test results

The file here contains results from the [MACAW](https://github.com/Devlin-Moyer/macaw) `dead_end_test` and `duplicate_test` tests, from a mass and charge balance report, and from cell-line specific gene essentiality prediction based on the [Hart _et al._ (2015)](https://doi.org/10.1016/j.cell.2015.11.015) dataset.

The test results shown here were obtained by the GitHub Actions run in:

- **PR #1046** (QC)
- **PR #973** (gene essentiality)

The results will be updated by any subsequent pull request. Summary results are shown as a comment in the corresponding pull request.

### MACAW: `dead_end_test`
Looks for metabolites in Human-GEM that can only be produced by all reactions they participate in or only consumed, then identifies all reactions that are prevented from sustaining steady-state fluxes because of each of these dead-end metabolites. The simplest case of a dead-end metabolite is one that only participates in a single reaction. Also flags all reversible reactions that can only carry fluxes in a single direction because one of their metabolites can either only be consumed or only be produced by all other reactions it participates in.

### MACAW: `duplicate_test`
Identifies sets of reactions that may be duplicates of each other because they:

- Involve exactly the same metabolites with exactly the same stoichiometric coefficients (but potentially different associated genes).
- Involve exactly the same metabolites, but go in different directions and/or some are reversible and some are not.
- Involve exactly the same metabolites, but with different stoichiometric coefficients.
- Represent the oxidation and/or reduction of the same metabolite, but use different electron acceptors/donors from the given list of pairs of oxidized and reduced forms of various electron carriers (e.g. NAD(H), NADP(H), FAD(H2), ubiquinone/ubiquinol, cytochromes).

It is possible for a single reaction to fit in multiple of the above categories. There are sometimes cases where sets of reactions that fall into one of the above categories are completely legitimate representations of real biochemistry (e.g. separate irreversible reactions for importing vs exporting the same metabolite because two different transporters encoded by different genes are each responsible for transporting that metabolite in only one direction, enzymes that can use NAD(H) or NADP(H) interchangeably to catalyze the same redox reaction), but reactions that meet these criteria are generally worth close examination to ensure that they should actually all exist as separate reactions.

### Mass and charge balance
Reports the reactions whose elemental (mass) or charge sums do not balance, using cobrapy's `check_mass_balance()`. Boundary reactions (exchange/demand/sink) and the biomass reaction are excluded, as they are not expected to balance. The unbalanced reactions are written to `balance_results.csv`, so a pull request that introduces a new imbalance is visible in the committed diff. This is a report and does not fail the build.

### Cell-line specific gene essentiality
Evaluate gene essentiality predictions in 5 cell-line specific GEMs with experimental fitness data gathered from the [Hart _et al._ (2015)](https://doi.org/10.1016/j.cell.2015.11.015).

Cell-line specific GEMs are constructed with tINIT2 for DLD1, GBM, HCT116, HeLa and RPE1 cell lines. Then, the `metabolicTasks_Essential.txt` list of tasks is used to identify essential genes in each of these models. The predicted gene essentiality is compared to results from a high-throughput CRISPR-Cas9 screen for identifying genes that affect fitness. Only the summary statistics of this comparison are kept.

### Model QC checks
A set of lightweight model quality-control checks, run by the `Model QC checks` workflow:

- `qc_metabolite_completeness.csv`: metabolites without a chemical formula or without a charge. Such metabolites are silently skipped by the mass and charge balance test, so tracking them keeps that test meaningful.
- `qc_reaction_sanity.csv`: reactions with invalid flux bounds (`lb > ub` or outside the standard +/-1000 range) or GPR issues (genes not annotated in `genes.tsv`, or a boundary reaction with a gene rule).
- `qc_annotation_issues.csv`: cross-reference problems in the annotation tables. Identifiers whose format does not match their namespace (KEGG, ChEBI, HMDB, PubChem, MetaNetX, Rhea, LipidMaps, EHMN, HepatoNET1, Reactome, TCDB), and metabolites whose cross-references are inconsistent across compartments (the same metabolite should carry the same identifiers in every compartment).
- `memote_score.md`: the total score from the [MEMOTE](https://memote.readthedocs.io) test suite, tracked so a pull request that changes it is visible in the diff. MEMOTE is split by cost: every pull request runs a fast core subset (it skips the flux-variability, stoichiometric-consistency-MILP and matrix-rank tests that dominate runtime on a genome-scale model), and pull requests to `main` run the complete suite. The scored MEMOTE result is uploaded as a build artifact.

The workflow also verifies the model can produce biomass under its default constraints (a growth sanity check), which is the one check that fails the build if it does not hold. The fast checks and MEMOTE run as two separate jobs, so the quick checks report without waiting for the (much slower) MEMOTE snapshot.

All of the results here are combined into a single pull-request comment (`model_qc_summary.md`): a compact status table for the model QC checks (current value, change compared to the target branch, and an icon for a quick visual check), followed by the MACAW and mass/charge balance summary and the gene-essentiality metrics. Each workflow rebuilds and posts that comment from the committed result files, so a result set that has not finished yet shows as pending. The per-finding detail stays in the CSVs above, which is where you look to find out what changed or why something failed.
