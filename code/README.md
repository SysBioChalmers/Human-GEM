# Human-GEM code

This directory contains functions and scripts that facilitate work with the Human-GEM model and other content in the repository. This directory is organized as follows:

```
code
├── GPRs
├── io
├── misc
├── modelCuration
│   ├── GPRs
│   │   └── EnzymeComplexes
│   ├── MetAssociation
│   ├── RxnAssociation
│   └── modelIntegration
├── qc
└── tINIT
```

### GPRs
Functions related to gene-transcript-protein-reaction (GTPR) associations in the model. The functions primarily involve the Human-GEM `grRules` field, facilitating processes such as cleaning or simpliyfing the rules, or converting between gene, transcript, and protein IDs.

### io
Functions associated with input/output of files into and out of MATLAB, such as the yaml-formatted Human-GEM, or documentation of model changes.

### misc
Code used for technical purposes, typically to augment missing or inconvenient functionalities of MATLAB.

### modelCuration
*DEPRECATED - WILL BE REMOVED*
Contains the curation scripts/functions that were used to help keep track of changes made to the Human-GEM `.mat` file because it is binary, and/or were implemented for a one-time curation action. As such, their only remaining purpose is for transparency and re-tracing the steps of the curation process.

Since these scripts will forever remain in past versions of this repository but are no longer in use, all code in this directory is considered deprecated and will soon be **removed from future repository versions**.

### qc
Functions to help with quality control (QC) of Human-GEM, such as checking for duplicate reactions or mass imbalances.

### tINIT
Functions associated with the **t**ask‐driven **I**ntegrative **N**etwork **I**nference for **T**issues (tINIT) algorithm. These functions include updates to the [original algorithm](https://www.ncbi.nlm.nih.gov/pubmed/24646661), as described in the [Human-GEM publication](https://stke.sciencemag.org/lookup/doi/10.1126/scisignal.aaz1482). Note that the updated tINIT implementation included here still requires many functions from the [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN).




