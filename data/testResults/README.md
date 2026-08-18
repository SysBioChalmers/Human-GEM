# Test results

This folder holds the committed quality-control (QC) results for Human-GEM. Every
pull request re-runs the checks, commits the updated files here, and posts a summary
as a single comment on the pull request (`model_qc_summary.md`). Each check name in
that comment links to the matching explanation in section 2 below.

The page has three parts:

1. [Where the current results come from](#1-where-the-current-results-come-from) - which pull request produced each file.
2. [What each check means](#2-what-each-check-means) - the tests shown in the pull-request comment.
3. [Files in this folder](#3-files-in-this-folder) - what every file here contains.

## 1. Where the current results come from

Result files are regenerated and committed by GitHub Actions. Most are produced
together by the **Model QC** workflow, which runs on every pull request. The full
MEMOTE suite and the gene-essentiality prediction take hours, so they run on demand
only (by commenting `/run memote` or `/run gene-essentiality`) and update just their
own files. The pull request in each row is the one whose run last wrote those files.

| Result file(s) | Produced by | Last updated by |
| --- | --- | --- |
| `qc_duplicate_keys.csv`, `qc_empty_reactions.csv`, `qc_annotation_consistency.csv`, `qc_deprecation_completeness.csv`, `qc_metabolite_completeness.csv`, `qc_reaction_sanity.csv`, `qc_duplicate_reactions.csv`, `qc_unused_entities.csv`, `qc_growth_blockers.csv` | `qcModelChecks.py` | **PR #1077** (model QC checks) |
| `qc_annotation_issues.csv` | `annotationTest.py` | **PR #1077** (model QC checks) |
| `qc_status.tsv` (round-trip, YAML lint, metabolic tasks, growth) | `testYamlConversion.py`, `testMetabolicTasks.py`, `action-yamllint`, `qcModelChecks.py` (via `qcStatus.py`) | **PR #1077** (model QC checks) |
| `macaw_results.csv`, `balance_results.csv`, `qc_structure_consistency.csv` | `macawTests.py`, `balanceTest.py`, `structureConsistencyTest.py` | **PR #1077** (MACAW and balance) |
| `memote_score.md` | `memoteSnapshot.py` (fast subset every PR; full suite via `/run memote`) | **PR #1077** (MEMOTE) |
| `gene-essential.csv`, `gene-essential_summary.md` | `geneEssentiality.py` via `/run gene-essentiality` | **PR #1027** (gene essentiality) |

## 2. What each check means

The headings below match the check names in the pull-request comment exactly, so a
name in the comment links straight to its explanation.

In the comment, a count links to the CSV that lists the exact entries, and the icon
reads: :white_check_mark: the count is zero (or the model grows / the score is
non-zero); :warning: a non-zero but pre-existing finding that this pull request did
not make worse (non-blocking); :x: a count that rose versus the target branch (a
regression this pull request introduced), or a failed gate.

### Model checks

Structural integrity and per-entity quality, all from `qcModelChecks.py` unless
noted. Two rows are build gates (a finding blocks the merge); the rest are reports.

#### Duplicate `!!omap` keys
**Gate.** Duplicate keys inside one metabolite/reaction/gene `!!omap` entry (two
`name` fields, or the same metabolite listed twice in a stoichiometry). RAVEN reads
and rewrites these, but `cobra.io.load_yaml_model` then raises a bare
`AssertionError` with no location, so the model cannot be loaded. The CSV names the
entry, key and line numbers.

#### Growth (biomass producible)
**Gate.** Whether the model can produce biomass under its default constraints
(`slim_optimize`). When it cannot, `qc_growth_blockers.csv` lists the biomass
precursors that cannot be made, which are what to fix.

#### Reactions with no metabolites
Reactions whose stoichiometry is empty. Such a reaction does nothing and usually
signals a broken edit.

#### Model / annotation-table inconsistencies
The model and its annotation tables (`reactions.tsv` / `metabolites.tsv` /
`genes.tsv`) must list the same identifiers. Flags identifiers in the model but not
the table (or the reverse), any deprecated identifier still used in the model, and a
non-numeric value in the `spontaneous` column of `reactions.tsv`.

#### Removed reactions or metabolites not deprecated
Human-GEM retires identifiers rather than deleting them, so a removed identifier
stays resolvable. This flags reactions or metabolites that are present on the target
branch but gone from this pull request's model and were **not** added to
`deprecatedReactions.tsv` / `deprecatedMetabolites.tsv`. A non-zero count means an
identifier was dropped without being moved to a deprecated list. (Comparison needs
the target-branch model tables, so it is reported only in CI.)

#### Metabolites missing formula
Metabolites with no chemical formula. They are silently skipped by the mass-balance
test, so tracking them keeps that test meaningful. Generic pool/class
pseudo-metabolites, which have no formula by design, are excluded.

#### Metabolites missing charge
Metabolites with no charge, for the same reason as the formula check.

#### Reaction bound / GPR issues
Reactions with invalid flux bounds (`lb > ub`, or a bound outside +/-1000) or
gene-rule problems (a gene not annotated in `genes.tsv`, or a boundary reaction that
carries a gene rule).

#### Exact-duplicate reaction groups
Groups of two or more reactions with **identical** stoichiometry (same metabolites
and same coefficients). This is the strict "truly identical" case; near-duplicates
(reverse direction, different coefficients, different electron carriers) are the
remit of the MACAW duplicate test below.

#### Unused metabolites
Metabolites not used by any reaction in the model.

#### Unused genes
Genes not referenced by any reaction's gene rule.

#### Malformed cross-references
From `annotationTest.py`. Cross-references in the annotation tables whose format does
not match their namespace (KEGG, ChEBI, HMDB, PubChem, MetaNetX, Rhea, LipidMaps,
EHMN, HepatoNET1, Reactome, TCDB).

#### Cross-refs inconsistent across compartments
From `annotationTest.py`. The same metabolite in different compartments carries
different cross-references, which should agree.

### MACAW and mass/charge balance

Network-level checks from [MACAW](https://github.com/Devlin-Moyer/macaw), the mass
and charge balance report, and the structure-vs-formula check.

#### Reactions flagged by MACAW dead-end test
Reactions prevented from carrying steady-state flux because one of their metabolites
can only ever be produced, or only consumed, by every reaction it takes part in (the
simplest case being a metabolite in a single reaction). Also flags reversible
reactions that can therefore run in only one direction.

#### Reactions flagged as MACAW duplicates
Sets of reactions that may be duplicates because they involve the same metabolites
(with the same or different coefficients or directions), or represent the same
oxidation/reduction using different electron carriers. Some are legitimate; the flag
means "worth checking", not "certainly wrong".

#### Mass-imbalanced reactions
Reactions whose elemental sums do not balance, from cobrapy's `check_mass_balance()`.
Boundary reactions (exchange/demand/sink) and biomass are excluded, since they are
not expected to balance.

#### Charge-imbalanced reactions
Reactions whose charge sums do not balance, with the same exclusions as above.

#### Structure vs formula/charge inconsistencies
From `structureConsistencyTest.py`. Metabolites whose structure (SMILES/InChI in
`metabolites.tsv`) implies a formula or charge that disagrees with the formula/charge
carried in the model YAML.

### Model file and metabolic tasks

Whether the model file survives conversion and satisfies the curated task lists.

#### YAML round-trip (cobrapy)
The model is loaded and re-written with cobrapy and must come back unchanged; a
failure means the YAML does not survive a cobrapy round-trip. **Gate.**

#### YAML round-trip (RAVEN)
The same round-trip through the RAVEN toolbox. **Gate.**

#### YAML lint
`yamllint` over `model/` (line-length rule disabled). **Gate.**

#### Essential metabolic tasks
The number of `metabolicTasks_Essential.txt` tasks the model passes, from
`testMetabolicTasks.py`. Any failure blocks the merge. **Gate.**

#### Verification metabolic tasks
The number of verification tasks the model passes, from `testMetabolicTasks.py`. Any
failure blocks the merge. **Gate.**

### MEMOTE
The total score, plus per-section and per-test scores, from the
[MEMOTE](https://memote.readthedocs.io) suite (`memoteSnapshot.py`). Every pull
request runs a fast core subset (skipping the flux-variability,
stoichiometric-consistency MILP and matrix-rank tests that dominate runtime).
Comment `/run memote` to run the full suite; the score then updates in place. Higher
is better, so the comment warns only when the score drops versus the target branch.

Before running, the model is enriched with the cross-references and SBO terms from
the annotation tables (the canonical `code/annotateGEM.py` helper), so the annotation
tests score against the identifiers Human-GEM actually carries rather than the bare
ids in the YAML. The enriched model is used only to build the temporary SBML MEMOTE
reads; it is not committed. The score is stored in two sections, `Core subset` and
`Full suite`.

### Gene essentiality (Hart 2015)
Gene-essentiality predictions in five cell-line-specific GEMs (DLD1, GBM, HCT116,
HeLa, RPE1), built with tINIT2 and evaluated against the CRISPR-Cas9 fitness screen
of [Hart _et al._ (2015)](https://doi.org/10.1016/j.cell.2015.11.015). This takes
hours and is not run on every pull request; comment `/run gene-essentiality` to run
it, and the result posts as its own comment. Only the summary statistics of the
comparison are kept here.

A gene counts as essential here when its knockout breaks *any* of the 57 essential
tasks, which mixes viability tasks (`GR` growth, `ER` energy and redox) with capability
tasks (`SU` substrate utilization, `BS` biosynthesis, `IC` internal conversions). Genes
essential only for a capability task, such as the ETF complex for beta-oxidation, are
therefore counted as false positives against a proliferation screen. `gradedEssentiality.py`
analyses that scoping: it reports the task categories behind every call together with a
continuous biomass growth ratio, and scores both against Hart. It is run manually and its
output is not committed here. See issue #1076.

## 3. Files in this folder

| File | Contents |
| --- | --- |
| `model_qc_summary.md` | The rendered pull-request comment (built by `buildReport.py` from the files below). Not a test itself. |
| `qc_status.tsv` | One key/value line each for the round-trip, YAML-lint and metabolic-task results and the growth value (see `qcStatus.py`). |
| `qc_duplicate_keys.csv` | Duplicate `!!omap` keys: entry, scope, key, first and duplicate line numbers. |
| `qc_growth_blockers.csv` | Biomass precursors that cannot be produced; empty when the model grows. |
| `qc_empty_reactions.csv` | Reactions with no metabolites. |
| `qc_annotation_consistency.csv` | Model-vs-annotation-table mismatches, deprecated-identifier use, and `spontaneous`-column problems: `kind, id, issue`. |
| `qc_deprecation_completeness.csv` | Reactions/metabolites removed since the target branch but not added to a deprecated list: `kind, id, issue`. |
| `qc_metabolite_completeness.csv` | Metabolites missing a formula and/or a charge: `metabolite, name, missing_formula, missing_charge`. |
| `qc_reaction_sanity.csv` | Reactions with bound or GPR issues: `reaction, name, issues`. |
| `qc_duplicate_reactions.csv` | Exact-duplicate reaction groups: `group, reaction, equation`. |
| `qc_unused_entities.csv` | Metabolites and genes used by no reaction: `kind, id`. |
| `qc_annotation_issues.csv` | Malformed and cross-compartment-inconsistent cross-references. |
| `qc_structure_consistency.csv` | Metabolites whose structure disagrees with the model formula/charge. |
| `macaw_results.csv` | Full MACAW output (dead-end and duplicate tests) per reaction. |
| `balance_results.csv` | Mass- and charge-imbalanced reactions. |
| `memote_score.md` | MEMOTE scores in two sections, core subset and full suite (see the MEMOTE explanation above). |
| `gene-essential.csv` | Per-gene essentiality matrix across the five cell-line models. |
| `gene-essential_summary.md` | Summary statistics of the gene-essentiality comparison against Hart 2015. |
| `README.md` | This file. |
</content>
</invoke>
