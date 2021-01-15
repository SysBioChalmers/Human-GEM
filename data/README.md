# Human-GEM data

This directory contains datasets, metabolic task files, and log files that support some functionality of Human-GEM and its associated scripts, and/or were used for the generation and curation of Human-GEM. The data subdirectories and their contents are briefly summarized below.

### Ensembl
Contains the Ensembl identifier mapping file, `ensembl_ID_mapping.tsv`, used to translate between different types of gene, protein, and transcript IDs (ENSG, ENSP, ENST, gene symbols, UniProt, NCBI).

### metabolicTasks
Contains metabolic task files used for evaluating the functionality of Human-GEM.
- `metabolicTasks_Essential.xlsx`: Metabolic tasks that are necessary for cell viability.
- `metabolicTasks_Full.xlsx`: A longer list of metabolic tasks which contains some functions that may not be relevant for all cell or tissue types. However, all of these tasks should pass for the generic Human-GEM.
- `metabolicTasks_VerifyModel.xlsx`: Metabolic tasks designed to verify that the model is properly mass and energy balanced, and exhibits a biologically meaningful solution space.

