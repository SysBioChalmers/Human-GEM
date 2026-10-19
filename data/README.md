# Human-GEM data

This directory contains datasets, metabolic task files, and log files that support some functionality of Human-GEM and its associated scripts, and/or were used for the generation and curation of Human-GEM. The data subdirectories and their contents are briefly summarized below.

### datasets
Data the gene-essentiality analysis (`code/test/gradedEssentiality.py`) compares the model with.
- `Hart2015_RNAseq.txt`, `Hart2015_TableS2.xlsx`: expression of the five cell lines and the CRISPR fitness screen of [Hart _et al._ (2015)](https://doi.org/10.1016/j.cell.2015.11.015).
- `HGNC_symbols_Hart2015.tsv`: the gene symbols of Hart 2015 Table S2 that are current [HGNC](https://www.genenames.org) symbols, with their Ensembl gene id (HGNC complete set of 2026-10-09); outdated Hart symbols are not listed. A Hart symbol is matched to a model gene through an alias only when the alias is not another gene's current symbol.

### metabolicTasks
Contains metabolic task files used for evaluating the functionality of Human-GEM.
- `metabolicTasks_Essential.txt`: Metabolic tasks that are necessary for cell viability.
- `metabolicTasks_Full.txt`: A longer list of metabolic tasks which contains some functions that may not be relevant for all cell or tissue types. However, all of these tasks should pass for the generic Human-GEM.
- `metabolicTasks_VerifyModel.txt`: Metabolic tasks designed to verify that the model is properly mass and energy balanced, and exhibits a biologically meaningful solution space.

### modelCuration
Contains curation-related data files used for making changes to the Human-GEM model. These model curation data files help to improve transparency of changes made to the model.

Note that many of the data files associated with model curation were outdated and not maintained, and have been moved to the `.deprecated` folder in the root directory of this repository.

### deprecatedIdentifiers
Contains the `tsv` files with reaction and metabolite annotation, identical to that in `/model/` for those identifiers that have been used in Human-GEM and subsequently removed. The identifiers are kept here so as to avoid accidentally reusing them in the future.
### thermodynamics
Data for the reversibility check in the Model QC workflow (see `data/testResults/README.md`).
- `reactionDeltaG.tsv`: ΔG'm (all reactants at 1 mM, at the compartment's pH) of the reactions that eQuilibrator can estimate, with its uncertainty and a hash of the stoichiometry it was estimated for. Made with `code/qc/estimateReactionDeltaG.py`; the Model QC workflow estimates new and changed reactions on every pull request and commits the updated table.
- `reversibilityExceptions.tsv`: reactions that may keep a direction the check would call impossible, each with the reason.
