# MACAW results

The file here contains results from the [MACAW](https://github.com/Devlin-Moyer/macaw) `dead_end_test` and `duplicate_test` tests.

The test results shown here were obtained by the GitHub Actions run in **PR #883**, and will be updated by any subsequent PR. Summary results are shown as a comment in the corresponding PR.

### `dead_end_test`
Looks for metabolites in Human-GEM that can only be produced by all reactions they participate in or only consumed, then identifies all reactions that are prevented from sustaining steady-state fluxes because of each of these dead-end metabolites. The simplest case of a dead-end metabolite is one that only participates in a single reaction. Also flags all reversible reactions that can only carry fluxes in a single direction because one of their metabolites can either only be consumed or only be produced by all other reactions it participates in.

### `duplicate_test`
Identifies sets of reactions that may be duplicates of each other because they:

- Involve exactly the same metabolites with exactly the same stoichiometric coefficients (but potentially different associated genes).
- Involve exactly the same metabolites, but go in different directions and/or some are reversible and some are not.
- Involve exactly the same metabolites, but with different stoichiometric coefficients.
- Represent the oxidation and/or reduction of the same metabolite, but use different electron acceptors/donors from the given list of pairs of oxidized and reduced forms of various electron carriers (e.g. NAD(H), NADP(H), FAD(H2), ubiquinone/ubiquinol, cytochromes).

It is possible for a single reaction to fit in multiple of the above categories. There are sometimes cases where sets of reactions that fall into one of the above categories are completely legitimate representations of real biochemistry (e.g. separate irreversible reactions for importing vs exporting the same metabolite because two different transporters encoded by different genes are each responsible for transporting that metabolite in only one direction, enzymes that can use NAD(H) or NADP(H) interchangeably to catalyze the same redox reaction), but reactions that meet these criteria are generally worth close examination to ensure that they should actually all exist as separate reactions.