# Human-GEM: The generic genome-scale metabolic model of _Homo sapiens_

[![Join the chat at https://gitter.im/SysBioChalmers/Human-GEM](https://badges.gitter.im/SysBioChalmers/Human-GEM.svg)](https://gitter.im/SysBioChalmers/Human-GEM?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM.svg)](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM)
[![DOI](https://zenodo.org/badge/105752644.svg)](https://zenodo.org/badge/latestdoi/105752644)<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-24-success.svg)](#contributors)
<!-- ALL-CONTRIBUTORS-BADGE:END --> 

### Brief Model Description

This repository contains the latest version of Human-GEM, a human genome-scale metabolic model. We encourage [contributions](#contributing).

### Cite us:

If you use Human1 in your research, please cite:  

 > J. L. Robinson, P. Kocabas�, H. Wang, P.-E. Cholley, et al. An atlas of human metabolism. _Sci. Signal._ 13, eaaz1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)
 
Starting from Human-GEM v1.5.0, all the releases are also archived in [Zenodo](https://doi.org/10.5281/zenodo.4099692) from which specific version can be cited if used.

If you use Mouse1, Rat1, Zebrafish1, Fruitfly1, or Worm1 in your research, please cite:   

  > H. Wang, J. L. Robinson, P. Kocabas�, J. Gustafsson, M. Anton, P.-E. Cholley, et al. Genome-scale metabolic network reconstruction of model animals as a platform for translational research. _PNAS_ 118, e2102344118 (2021). [doi.org/10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118)



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
|_Homo sapiens_ |   HMR2, Recon3D, iHsa|    13070|  8369|  3067|


## Contributing

Contributions are always welcome! Read more about the project's philosophy in our [wiki](https://github.com/SysBioChalmers/Human-GEM/wiki) or have a look at the [Contributing guidelines](https://github.com/SysBioChalmers/Human-GEM/blob/main/.github/CONTRIBUTING.md) before starting.


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

Additional annotation information and external identifiers for Human-GEM reactions, metabolites, and genes are provided as `tsv` files in the `model/` directory (`reactions.tsv`, `metabolites.tsv`, and `genes.tsv`, respectively).  

To import/export this annotation data to/from MATLAB, use the `importTsvFile` and `exportTsvFile` functions, respectively.


## Websites

- [Metabolic Atlas](https://metabolicatlas.org/) enables visualization and exploration of Human-GEM content.
- The [Human-GEM user guide](https://sysbiochalmers.github.io/Human-GEM-guide/) provides detailed instructions and examples for using the Human-GEM model and repository.


## Metabolic Maps

A collection of manually curated 2D metabolic maps associated with Human-GEM are stored in the [Human-maps repository](https://github.com/SysBioChalmers/Human-maps). These maps can be downloaded from the repository or explored interactively using [Metabolic Atlas](https://metabolicatlas.org/explore/map-viewer/human1).


## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center"><a href="https://github.com/ANiknejad"><img src="https://avatars.githubusercontent.com/u/2682520?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Anne Niknejad</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AANiknejad" title="Bug reports">🐛</a> <a href="#content-ANiknejad" title="Content">🖋</a></td>
      <td align="center"><a href="https://github.com/avlant"><img src="https://avatars.githubusercontent.com/u/5329888?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Avlant</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Aavlant" title="Bug reports">🐛</a> <a href="#content-avlant" title="Content">🖋</a></td>
      <td align="center"><a href="https://github.com/BenjaSanchez"><img src="https://avatars.githubusercontent.com/u/9384349?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Benjamín Sánchez</b></sub></a><br /><a href="#question-BenjaSanchez" title="Answering Questions">💬</a></td>
      <td align="center"><a href="https://github.com/dweindl"><img src="https://avatars.githubusercontent.com/u/18048784?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Daniel Weindl</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Adweindl" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://github.com/edkerk"><img src="https://avatars.githubusercontent.com/u/7326655?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Eduard Kerkhoven</b></sub></a><br /><a href="#question-edkerk" title="Answering Questions">💬</a></td>
      <td align="center"><a href="https://orcid.org/0000-0001-7475-0136"><img src="https://avatars.githubusercontent.com/u/21077367?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Hao Wang</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AHao-Chalmers" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=Hao-Chalmers" title="Code">💻</a> <a href="#data-Hao-Chalmers" title="Data">🔣</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=Hao-Chalmers" title="Documentation">📖</a> <a href="#ideas-Hao-Chalmers" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-Hao-Chalmers" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#platform-Hao-Chalmers" title="Packaging/porting to new platform">📦</a> <a href="#projectManagement-Hao-Chalmers" title="Project Management">📆</a> <a href="#question-Hao-Chalmers" title="Answering Questions">💬</a> <a href="#research-Hao-Chalmers" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3AHao-Chalmers" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=Hao-Chalmers" title="Tests">⚠️</a> <a href="#talk-Hao-Chalmers" title="Talks">📢</a></td>
      <td align="center"><a href="https://jonathanrob.github.io"><img src="https://avatars.githubusercontent.com/u/22366558?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Jonathan Robinson</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AJonathanRob" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=JonathanRob" title="Code">💻</a> <a href="#data-JonathanRob" title="Data">🔣</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=JonathanRob" title="Documentation">📖</a> <a href="#ideas-JonathanRob" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-JonathanRob" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#platform-JonathanRob" title="Packaging/porting to new platform">📦</a> <a href="#projectManagement-JonathanRob" title="Project Management">📆</a> <a href="#question-JonathanRob" title="Answering Questions">💬</a> <a href="#research-JonathanRob" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3AJonathanRob" title="Reviewed Pull Requests">👀</a> <a href="#tutorial-JonathanRob" title="Tutorials">✅</a> <a href="#talk-JonathanRob" title="Talks">📢</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/10344158?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Jorge Ferreira</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ajorgemlferreira" title="Bug reports">🐛</a></td>
    </tr>
    <tr>
      <td align="center"><a href="https://github.com/CadavidJoseL"><img src="https://avatars.githubusercontent.com/u/62765618?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Jose L. Cadavid</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3ACadavidJoseL" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://github.com/mpagni12"><img src="https://avatars.githubusercontent.com/u/45748199?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Marco Pagni</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ampagni12" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://orcid.org/0000-0002-7753-9042"><img src="https://avatars.githubusercontent.com/u/23480589?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Mihail Anton</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Amihai-sysbio" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=mihai-sysbio" title="Code">💻</a> <a href="#ideas-mihai-sysbio" title="Ideas, Planning, & Feedback">🤔</a> <a href="#infra-mihai-sysbio" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Amihai-sysbio" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=mihai-sysbio" title="Tests">⚠️</a> <a href="#talk-mihai-sysbio" title="Talks">📢</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/26245751?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Pierre-Etienne Cholley</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apecholleyc" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=pecholleyc" title="Code">💻</a> <a href="#content-pecholleyc" title="Content">🖋</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Apecholleyc" title="Reviewed Pull Requests">👀</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/2399043?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Pierre-Etienne Cholley</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apecholley" title="Bug reports">🐛</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/8766764?v=4?s=80" width="80px;" alt=""/><br /><sub><b>PkiwiBird</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3APkiwiBird" title="Bug reports">🐛</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/38076281?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Pranas Grigaitis</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apranasag" title="Bug reports">🐛</a> <a href="#content-pranasag" title="Content">🖋</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=pranasag" title="Code">💻</a></td>
      <td align="center"><img src="https://avatars.githubusercontent.com/u/32029599?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Pınar Kocabaş</b></sub><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apinarkocabas" title="Bug reports">🐛</a></td>
    </tr>
    <tr>
      <td align="center"><a href="https://github.com/Rasools"><img src="https://avatars.githubusercontent.com/u/22166601?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Rasool Saghaleyni</b></sub></a><br /><a href="#ideas-Rasools" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center"><a href="https://github.com/cherkaos"><img src="https://avatars.githubusercontent.com/u/4625396?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Sarah Cherkaoui</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Acherkaos" title="Bug reports">🐛</a> <a href="#content-cherkaos" title="Content">🖋</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=cherkaos" title="Code">💻</a></td>
      <td align="center"><a href="https://github.com/simas232"><img src="https://avatars.githubusercontent.com/u/11994076?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Simonas Marcišauskas</b></sub></a><br /><a href="#question-simas232" title="Answering Questions">💬</a></td>
      <td align="center"><a href="https://github.com/TunahanCakir"><img src="https://avatars.githubusercontent.com/u/71440332?v=4?s=80" width="80px;" alt=""/><br /><sub><b>TunahanCakir</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3ATunahanCakir" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://github.com/XuhangLi"><img src="https://avatars.githubusercontent.com/u/41695293?v=4?s=80" width="80px;" alt=""/><br /><sub><b>Xuhang Li</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AXuhangLi" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://github.com/johan-gson"><img src="https://avatars.githubusercontent.com/u/32481323?v=4?s=80" width="80px;" alt=""/><br /><sub><b>johan-gson</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ajohan-gson" title="Bug reports">🐛</a> <a href="#content-johan-gson" title="Content">🖋</a></td>
      <td align="center"><a href="https://orcid.org/0000-0003-3947-488X"><img src="https://avatars.githubusercontent.com/u/3072880?v=4?s=80" width="80px;" alt=""/><br /><sub><b>smoretti</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Asmoretti" title="Bug reports">🐛</a></td>
      <td align="center"><a href="https://github.com/stairs"><img src="https://avatars.githubusercontent.com/u/6586371?v=4?s=80" width="80px;" alt=""/><br /><sub><b>stairs</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Astairs" title="Bug reports">🐛</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->