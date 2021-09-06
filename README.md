# Human-GEM: The generic genome-scale metabolic model of _Homo sapiens_
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-1-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

[![Join the chat at https://gitter.im/SysBioChalmers/Human-GEM](https://badges.gitter.im/SysBioChalmers/Human-GEM.svg)](https://gitter.im/SysBioChalmers/Human-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM.svg)](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM)
[![DOI](https://zenodo.org/badge/105752644.svg)](https://zenodo.org/badge/latestdoi/105752644)

### Brief Model Description

This repository contains the latest version of Human-GEM, a human genome-scale metabolic model.

### Cite us:

If you use Human1 in your research, please cite:  

 > J. L. Robinson, P. Kocabasÿ, H. Wang, P.-E. Cholley, et al. An atlas of human metabolism. _Sci. Signal._ 13, eaaz1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)
 
Starting from Human-GEM v1.5.0, all the releases are also archived in [Zenodo](https://doi.org/10.5281/zenodo.4099692) from which specific version can be cited if used.

If you use Mouse1, Rat1, Zebrafish1, Fruitfly1 or Worm1 in your research, please cite:   

  > H. Wang, J. L. Robinson, P. Kocabasÿ, J. Gustafsson, M. Anton, P.-E. Cholley, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. _PNAS_ 118, e2102344118 (2021). [doi.org/10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118)



### Model Keywords

**Utilisation:** predictive simulation, multi-omics integrative analysis, model template  
**Field:** metabolic-network reconstruction  
**Type of Model:** reconstruction, curated  
**Model Source:** HPA, HMR2, iHsa, iHepatocytes2322, Recon3D  
**Omic Source:** genomics, proteomics  
**Taxonomy:** _Homo sapiens_  
**Metabolic System:** general metabolism  
**Condition:** generic metabolism  


### Model Overview

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
|_Homo sapiens_ |   HMR2, Recon3D, iHsa|    13081|  8371|  3625|


### Administration

This repository is administered by Jonathan Robinson ([@JonathanRob](https://github.com/jonathanrob)) and Hao Wang ([@Hao-Chalmers](https://github.com/hao-chalmers)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


## User Guide

Detailed instructions on the installation and use of the Human-GEM model and repository can be found in the [Human-GEM user guide](https://sysbiochalmers.github.io/Human-GEM-guide/).


## Installation

### Required Software
* A functional MATLAB installation (MATLAB 7.3 and higher).
* The [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN).
* The [COBRA toolbox](https://github.com/opencobra/cobratoolbox) (not necessary for most functionality).


### Dependencies - Recommended Software
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.


### Installation Instructions
* Clone the [main branch](https://github.com/SysBioChalmers/Human-GEM/tree/main) of this repository, or [download the latest release](https://github.com/SysBioChalmers/Human-GEM/releases/latest).
* Add the directory to your MATLAB path (instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)).



## Model Files

The model is available as `.xml`, `.xlsx`, `.txt`, `.yml`, and `.mat` in the `model/` directory. Note that only the `.yml` version is available on branches other than `main` (e.g., `develop`), to facilitate tracking of model changes.



## Usage

#### Loading/saving the model

`Human-GEM.mat` (Recommended if on `main` branch)
* Load and save using the built-in MATLAB `load()` and `save()` functions.

`Human-GEM.yml` (Recommended if on `develop` or other branches)
* Load using the `importYaml.m` function (in `code/io/`)
* Save using the `exportYaml.m` function (in `code/io/`)

`Human-GEM.xml` (SBML format)
* Load using the `importModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))
* Save using the `exportModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN))


## Reaction, Metabolite, and Gene Annotations

Additional annotation information and external identifiers for Human-GEM reactions, metabolites, and genes are provided as `tsv` files in the `model/` directory. The structure of the `tsv` files are tabulated below.

* `reactions.tsv` content:

fieldname      |  annotation             |
---------------|------------------------ |
rxns           |identical to `model.rxns`|
rxnKEGGID      |KEGG reaction ID        |
rxnBiGGID      |BiGG reaction ID        |
rxnEHMNID      |EHMN reaction ID        |
rxnHepatoNET1ID|HepatoNET1 reaction ID  |
rxnREACTOMEID  |REACTOME ID             |
rxnRecon3DID   |Recon3D reaction ID     |
rxnMetaNetXID  |MetaNetX reaction ID    |
rxnHMR2ID      |HMR2 reaction ID        |
rxnRatconID    |Ratcon reaction ID      |
rxnTCDBID      |TCDB ID                 |
spontaneous    |Spontaneous status      |
rxnMAID        |MA reaction ID          |
rxnRheaID      |Rhea ID                 |
rxnRheaMasterID|Master Rhea ID          |


* `metabolites.tsv` content:

fieldname      |  annotation             |
---------------|-------------------------|
mets           |identical to `model.mets`|
metsNoComp     |`model.mets` without compartment suffix|
metBiGGID      |BiGG metabolite ID       |
metKEGGID      |KEGG metabolite ID       |
metHMDBID      |HMDB ID                  |
metChEBIID     |ChEBI ID                 |
metPubChemID   |PubChem ID               |
metLipidMapsID |LipidMaps ID             |
metEHMNID      |EHMN metabolite ID       |
metHepatoNET1ID|HepatoNET1 metabolite ID |
metRecon3DID   |Recon3D metabolite ID    |
metMetaNetXID  |MetaNetX metabolite ID   |
metHMR2ID      |HMR2 metabolite ID       |
metMAID        |MA metabolite ID         |


* `genes.tsv` content:

fieldname     |  annotation          |
--------------|----------------------|
genes         |Ensembl gene ID       |
geneENSTID    |Ensembl transcript ID |
geneENSPID    |Ensembl protein ID |
geneUniProtID |UniProt ID            |
geneSymbols   |Gene Symbol           |
geneEntrezID  |NCBI Entrez ID        |
geneNames     |Gene Name             |
geneAliases   |Alias Names           |


To import/export this annotation data to/from MATLAB, use the `importTsvFile` and `exportTsvFile` functions, respectively.



## Websites

- [Metabolic Atlas](https://metabolicatlas.org/) enables visualization and exploration of Human-GEM content.
- The [Human-GEM user guide](https://sysbiochalmers.github.io/Human-GEM-guide/) provides detailed instructions and examples for using the Human-GEM model and repository.



## Metabolic Maps

A collection of manually curated 2D metabolic maps associated with Human-GEM are stored in the [Human-maps repository](https://github.com/SysBioChalmers/Human-maps). These maps can be downloaded from the repository or explored interactively using [Metabolic Atlas](https://metabolicatlas.org/explore/map-viewer/human1).



## Contributing

Contributions are always welcome! Please read the [contribution guidelines](https://github.com/SysBioChalmers/Human-GEM/blob/main/.github/CONTRIBUTING.md) to get started.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://www.helmholtz-muenchen.de/icb/institute/staff/staff/ma/5122/index.html"><img src="https://avatars.githubusercontent.com/u/18048784?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Daniel Weindl</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Adweindl" title="Bug reports">🐛</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!