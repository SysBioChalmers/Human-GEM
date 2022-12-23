# Model curation scripts

This directory contains curation-related scripts and functions used to make changes to the Human-GEM repository. These curation scripts help to improve transparency of changes made to the model when the number of changes is too large to view practically.

- `getCompFromUniprotCellAtlas.py`: Code for collecting subcellular localization information for existing metabolic enzymes from **Swissprot** and **Cell Atlas** ([HPA](https://www.proteinatlas.org/search/has_protein_data_in%3ACell)), and then incorporating the compartment info from both sources into `genes.tsv`.
