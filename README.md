# Human-GEM: The generic genome-scale metabolic model of _Homo sapiens_
[![GitHub version](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM.svg)](https://badge.fury.io/gh/sysbiochalmers%2FHuman-GEM)
[![DOI](https://zenodo.org/badge/105752644.svg)](https://zenodo.org/badge/latestdoi/105752644)<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-46-success.svg)](#contributors)
<!-- ALL-CONTRIBUTORS-BADGE:END --> 

### Brief model description
This repository contains the latest version of Human-GEM, a human genome-scale metabolic model. We encourage [contributions](#contributing).

### Cite us:
If you use Human2 in your research, please cite:  
 > Luo J, Wang H, Moyer D, Guo Z, Robinson JL, Gustafsson J, Anton M, Chen Y, Kerkhoven EJ, Nielsen J, Li F. Reconstruction of human metabolic models with large language models. _PNAS_ 123.15:e2516511123 (2026). [doi:10.1073/pnas.2516511123](https://doi.org/10.1073/pnas.2516511123)

If you use Human1 in your research, please cite:  
 > Robinson JL, et al. An atlas of human metabolism. _Sci. Signal._ 13, eaaz1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)
 
Starting from Human-GEM v1.5.0, all the releases are also archived in [Zenodo](https://doi.org/10.5281/zenodo.4099692) from which specific version can be cited if used.

### Model keywords
**Utilisation:** predictive simulation, multi-omics integrative analysis, model template  
**Field:** metabolic-network reconstruction  
**Type of Model:** reconstruction, curated  
**Model Source:** HPA, HMR2, iHsa, iHepatocytes2322, Recon3D  
**Omic Source:** genomics, proteomics  
**Taxonomy:** _Homo sapiens_  
**Metabolic System:** general metabolism  
**Condition:** generic metabolism  

### Model overview
|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|:-----:|
|_Homo sapiens_ |   HMR2, Recon3D, iHsa|    {{nRXN}}|  {{nMET}}|  {{nGENE}}|

## Contributing
Contributions are always welcome! Read more about the project's philosophy in our [wiki](https://github.com/SysBioChalmers/Human-GEM/wiki) or have a look at the [Contributing guidelines](https://github.com/SysBioChalmers/Human-GEM/blob/main/.github/CONTRIBUTING.md) before starting.

## User guide
Detailed instructions on the installation and use of the Human-GEM model and repository can be found in the [Human-GEM user guide](https://sysbiochalmers.github.io/Human-GEM-guide/).

# Installation
## Required software
### Basic user
If you want to use the model for your own model simulations, you can use any software that accepts **SBML L3V1 FBCv3** formatted model files. This includes any of the following:

#### MATLAB-based
* [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) v2.10.3 or later (recommended, see [Installation instructions](https://github.com/SysBioChalmers/RAVEN/wiki/Installation#installation-instructions))
* [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)

#### Python-based
* [cobrapy](https://github.com/opencobra/cobrapy)  

Please see the installation instructions for each software package.

### Developer
#### MATLAB-based  
If you want to contribute to the development of Human-GEM, or otherwise want to run any of the [provided](https://github.com/SysBioChalmers/Human-GEM/tree/main/code) MATLAB functions, then the following software is required:
* [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) v2.10.3 or later (recommended, see [Installation instructions](https://github.com/SysBioChalmers/RAVEN/wiki/Installation#installation-instructions))

#### Python-based  
You can also contribution to the development of Human-GEM via python (e.g. cobrapy), even if you would not be able to run any of the model-specific MATLAB functions. To curate the model, you can still edit `Human-GEM.yml`, `genes.tsv`, `metabolites.tsv` and `reactions.tsv`, all located in the `model/` folder.

### Recommended solver
* When performing simulations with Human-GEM, you are encouraged to use [Gurobi Optimizer](https://www.gurobi.com/academia/academic-program-and-licenses/).

## Installation instructions
* Clone the [main branch](https://github.com/SysBioChalmers/Human-GEM/tree/main) of this repository, or [download the latest release](https://github.com/SysBioChalmers/Human-GEM/releases/latest).
* Add the directory to your MATLAB path either by using the lines below or manually (instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)).
```matlab
% Replace "/my/path/" with the actual path to the Human-GEM folder
cd /my/path/Human-GEM/code
% This will add the relevant paths to the path variable in MATLAB
HumanGEMInstaller.install
% It is also possible to remove Human-GEM from the MATLAB path using
HumanGEMInstaller.uninstall
```

## Model files
The model is available as `.xml`, `.xlsx`, `.txt`, `.yml`, and `.mat` in the `model/` directory.

The model development is based on the `yml` file (to facilitate tracking of model changes), so that this is the only model format available on non-`main` branches (such as `develop`). See also the [location of annotation information](#reaction-metabolite-and-gene-annotations).

## Usage
#### Loading/saving the model in MATLAB
`Human-GEM.mat` (Recommended if on `main` branch):
* Load and save using the built-in MATLAB `load()` and `save()` functions.

`Human-GEM.yml` (Recommended if on `develop` or other branches):
* Load using the `readYAMLmodel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)).
* Save using the `writeYAMLmodel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)).

`Human-GEM.xml` (SBML format):
* Load using the `importModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)).
* Save using the `exportModel.m` function (from [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)).

## Reaction, metabolite, and gene annotations
By default, additional annotation information and external identifiers for Human-GEM reactions, metabolites, and genes are **not** kept in the `yml` file, but rather provided as `tsv` files in the `model/` directory (`reactions.tsv`, `metabolites.tsv`, and `genes.tsv`, respectively). The `annotateGEM` function can add those annotations in MATLAB, while direct import/export of this annotation data is done through the `importTsvFile` and `exportTsvFile` functions, respectively.

## Websites
- [Metabolic Atlas](https://metabolicatlas.org/) enables visualization and exploration of Human-GEM content.
- The [Human-GEM user guide](https://sysbiochalmers.github.io/Human-GEM-guide/) provides detailed instructions and examples for using the Human-GEM model and repository.

## Metabolic maps
A collection of manually curated 2D metabolic maps associated with Human-GEM are stored in the [Human-maps repository](https://github.com/SysBioChalmers/Human-maps). These maps can be downloaded from the repository or explored interactively using [Metabolic Atlas](https://metabolicatlas.org/explore/map-viewer/human1).

## Contributors
Human-GEM is built by everyone who reports an error, brings the evidence for a correction, writes code or reviews a change. If you open an issue that leads to a fix, you are added here; see [Crediting contributors](.github/CONTRIBUTING.md#crediting-contributors).

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

### Core team

<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/edkerk"><img src="https://avatars.githubusercontent.com/u/7326655?v=4?s=64" width="64px;" alt="Eduard Kerkhoven"/><br /><sub><b>Eduard Kerkhoven</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Aedkerk" title="Bug reports">🐛</a> <a href="#research-edkerk" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=edkerk" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Aedkerk" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=edkerk" title="Documentation">📖</a> <a href="#ideas-edkerk" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/feiranl"><img src="https://avatars.githubusercontent.com/u/32157802?v=4?s=64" width="64px;" alt="feiranl"/><br /><sub><b>feiranl</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Afeiranl" title="Bug reports">🐛</a> <a href="#research-feiranl" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=feiranl" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Afeiranl" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=feiranl" title="Documentation">📖</a> <a href="#ideas-feiranl" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://orcid.org/0000-0002-7111-1360"><img src="https://avatars.githubusercontent.com/u/67491919?v=4?s=64" width="64px;" alt="Jiahao Luo"/><br /><sub><b>Jiahao Luo</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AJHL-452b" title="Bug reports">🐛</a> <a href="#research-JHL-452b" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=JHL-452b" title="Code">💻</a> <a href="#ideas-JHL-452b" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://orcid.org/0000-0002-6997-9531"><img src="https://avatars.githubusercontent.com/u/33460176?v=4?s=64" width="64px;" alt="Devlin Moyer"/><br /><sub><b>Devlin Moyer</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3ADevlin-Moyer" title="Bug reports">🐛</a> <a href="#research-Devlin-Moyer" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=Devlin-Moyer" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3ADevlin-Moyer" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://orcid.org/0000-0001-7475-0136"><img src="https://avatars.githubusercontent.com/u/21077367?v=4?s=64" width="64px;" alt="Hao Wang"/><br /><sub><b>Hao Wang</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ahaowang-bioinfo" title="Bug reports">🐛</a> <a href="#research-haowang-bioinfo" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=haowang-bioinfo" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Ahaowang-bioinfo" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=haowang-bioinfo" title="Documentation">📖</a> <a href="#ideas-haowang-bioinfo" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://jonathanrob.github.io"><img src="https://avatars.githubusercontent.com/u/22366558?v=4?s=64" width="64px;" alt="Jonathan Robinson"/><br /><sub><b>Jonathan Robinson</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AJonathanRob" title="Bug reports">🐛</a> <a href="#research-JonathanRob" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=JonathanRob" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3AJonathanRob" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=JonathanRob" title="Documentation">📖</a> <a href="#ideas-JonathanRob" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/pecholleyc"><img src="https://avatars.githubusercontent.com/u/26245751?v=4?s=64" width="64px;" alt="Pierre-Etienne Cholley"/><br /><sub><b>Pierre-Etienne Cholley</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apecholleyc" title="Bug reports">🐛</a> <a href="#research-pecholleyc" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=pecholleyc" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Apecholleyc" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://orcid.org/0000-0002-7753-9042"><img src="https://avatars.githubusercontent.com/u/23480589?v=4?s=64" width="64px;" alt="Mihail Anton"/><br /><sub><b>Mihail Anton</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Amihai-sysbio" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=mihai-sysbio" title="Code">💻</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Amihai-sysbio" title="Reviewed Pull Requests">👀</a> <a href="#ideas-mihai-sysbio" title="Ideas, Planning, & Feedback">🤔</a></td>
    </tr>
  </tbody>
</table>

### Curators and developers

<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/johan-gson"><img src="https://avatars.githubusercontent.com/u/32481323?v=4?s=48" width="48px;" alt="johan-gson"/><br /><sub><b>johan-gson</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ajohan-gson" title="Bug reports">🐛</a> <a href="#research-johan-gson" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=johan-gson" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/PkiwiBird"><img src="https://avatars.githubusercontent.com/u/8766764?v=4?s=48" width="48px;" alt="PkiwiBird"/><br /><sub><b>PkiwiBird</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3APkiwiBird" title="Bug reports">🐛</a> <a href="#research-PkiwiBird" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=PkiwiBird" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/pranasag"><img src="https://avatars.githubusercontent.com/u/38076281?v=4?s=48" width="48px;" alt="Pranas Grigaitis"/><br /><sub><b>Pranas Grigaitis</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apranasag" title="Bug reports">🐛</a> <a href="#research-pranasag" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=pranasag" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/pinarkocabas"><img src="https://avatars.githubusercontent.com/u/32029599?v=4?s=48" width="48px;" alt="Pınar Kocabaş"/><br /><sub><b>Pınar Kocabaş</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Apinarkocabas" title="Bug reports">🐛</a> <a href="#research-pinarkocabas" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/pulls?q=is%3Apr+reviewed-by%3Apinarkocabas" title="Reviewed Pull Requests">👀</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/cherkaos"><img src="https://avatars.githubusercontent.com/u/4625396?v=4?s=48" width="48px;" alt="Sarah Cherkaoui"/><br /><sub><b>Sarah Cherkaoui</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Acherkaos" title="Bug reports">🐛</a> <a href="#research-cherkaos" title="Research">🔬</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=cherkaos" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/ANiknejad"><img src="https://avatars.githubusercontent.com/u/2682520?v=4?s=48" width="48px;" alt="Anne Niknejad"/><br /><sub><b>Anne Niknejad</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3AANiknejad" title="Bug reports">🐛</a> <a href="#research-ANiknejad" title="Research">🔬</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/avlant"><img src="https://avatars.githubusercontent.com/u/5329888?v=4?s=48" width="48px;" alt="Avlant"/><br /><sub><b>Avlant</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Aavlant" title="Bug reports">🐛</a> <a href="#research-avlant" title="Research">🔬</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://www.sriramlab.org/"><img src="https://avatars.githubusercontent.com/u/43740320?v=4?s=48" width="48px;" alt="Carolina Chung, PhD"/><br /><sub><b>Carolina Chung, PhD</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Achechung" title="Bug reports">🐛</a> <a href="#research-chechung" title="Research">🔬</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/h-escoffier"><img src="https://avatars.githubusercontent.com/u/85628846?v=4?s=48" width="48px;" alt="Hugues Esc_"/><br /><sub><b>Hugues Esc_</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ah-escoffier" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=h-escoffier" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/CadavidJoseL"><img src="https://avatars.githubusercontent.com/u/62765618?v=4?s=48" width="48px;" alt="Jose L. Cadavid"/><br /><sub><b>Jose L. Cadavid</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3ACadavidJoseL" title="Bug reports">🐛</a> <a href="#research-CadavidJoseL" title="Research">🔬</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/liamkelley93"><img src="https://avatars.githubusercontent.com/u/9171240?v=4?s=48" width="48px;" alt="Liam"/><br /><sub><b>Liam</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Aliamkelley93" title="Bug reports">🐛</a> <a href="https://github.com/SysBioChalmers/Human-GEM/commits?author=liamkelley93" title="Code">💻</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/mpagni12"><img src="https://avatars.githubusercontent.com/u/45748199?v=4?s=48" width="48px;" alt="Marco Pagni"/><br /><sub><b>Marco Pagni</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Ampagni12" title="Bug reports">🐛</a> <a href="#research-mpagni12" title="Research">🔬</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/Rasools"><img src="https://avatars.githubusercontent.com/u/22166601?v=4?s=48" width="48px;" alt="Rasool Saghaleyni"/><br /><sub><b>Rasool Saghaleyni</b></sub></a><br /><a href="#research-Rasools" title="Research">🔬</a> <a href="#ideas-Rasools" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://orcid.org/0000-0003-3947-488X"><img src="https://avatars.githubusercontent.com/u/3072880?v=4?s=48" width="48px;" alt="smoretti"/><br /><sub><b>smoretti</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Asmoretti" title="Bug reports">🐛</a> <a href="#research-smoretti" title="Research">🔬</a></td>
      <td align="center" valign="top" width="12.5%"><a href="https://github.com/stairs"><img src="https://avatars.githubusercontent.com/u/6586371?v=4?s=48" width="48px;" alt="stairs"/><br /><sub><b>stairs</b></sub></a><br /><a href="https://github.com/SysBioChalmers/Human-GEM/issues?q=author%3Astairs" title="Bug reports">🐛</a> <a href="#research-stairs" title="Research">🔬</a></td>
    </tr>
  </tbody>
</table>

### Error reports and ideas

<p>
<a href="https://github.com/BenjaSanchez" title="Benjamín Sánchez: ideas, planning, & feedback"><img src="https://avatars.githubusercontent.com/u/9384349?v=4?s=32" width="32px;" alt="Benjamín Sánchez"/></a>
<a href="https://github.com/Christoff1993" title="Christoff1993: bug reports"><img src="https://avatars.githubusercontent.com/u/95428150?v=4?s=32" width="32px;" alt="Christoff1993"/></a>
<a href="https://github.com/dweindl" title="Daniel Weindl: bug reports"><img src="https://avatars.githubusercontent.com/u/18048784?v=4?s=32" width="32px;" alt="Daniel Weindl"/></a>
<a href="https://www.oxinabox.net/" title="Frames White: bug reports"><img src="https://avatars.githubusercontent.com/u/5127634?v=4?s=32" width="32px;" alt="Frames White"/></a>
<a href="https://github.com/hchapman1" title="hchapman1: ideas, planning, & feedback"><img src="https://avatars.githubusercontent.com/u/127259422?v=4?s=32" width="32px;" alt="hchapman1"/></a>
<a href="https://github.com/ina999111" title="Ina Maltais-Payette: bug reports"><img src="https://avatars.githubusercontent.com/u/182758843?v=4?s=32" width="32px;" alt="Ina Maltais-Payette"/></a>
<a href="https://github.com/IVANDOMENZAIN" title="Iván Domenzain: ideas, planning, & feedback"><img src="https://avatars.githubusercontent.com/u/26483972?v=4?s=32" width="32px;" alt="Iván Domenzain"/></a>
<a href="https://github.com/dagl1" title="Jelle Bonthuis: bug reports"><img src="https://avatars.githubusercontent.com/u/24440380?v=4?s=32" width="32px;" alt="Jelle Bonthuis"/></a>
<a href="https://github.com/jorgemlferreira" title="Jorge Ferreira: bug reports"><img src="https://avatars.githubusercontent.com/u/10344158?v=4?s=32" width="32px;" alt="Jorge Ferreira"/></a>
<a href="https://github.com/juliette-cooke" title="Juliette: bug reports"><img src="https://avatars.githubusercontent.com/u/90753730?v=4?s=32" width="32px;" alt="Juliette"/></a>
<a href="https://github.com/jreimertz" title="Justin Reimertz: bug reports"><img src="https://avatars.githubusercontent.com/u/119626411?v=4?s=32" width="32px;" alt="Justin Reimertz"/></a>
<a href="https://github.com/manas-kohli" title="manas-kohli: bug reports"><img src="https://avatars.githubusercontent.com/u/5028979?v=4?s=32" width="32px;" alt="manas-kohli"/></a>
<a href="https://mmtf.dev" title="Mehmet Efe Akça: bug reports"><img src="https://avatars.githubusercontent.com/u/13402668?v=4?s=32" width="32px;" alt="Mehmet Efe Akça"/></a>
<a href="https://github.com/mfcesur" title="Müberra Fatma Cesur: bug reports"><img src="https://avatars.githubusercontent.com/u/84909981?v=4?s=32" width="32px;" alt="Müberra Fatma Cesur"/></a>
<a href="https://uni-tuebingen.de/en/206112" title="Nantia Leonidou: bug reports"><img src="https://avatars.githubusercontent.com/u/28543617?v=4?s=32" width="32px;" alt="Nantia Leonidou"/></a>
<a href="https://github.com/zrf492174" title="Runfeng Zhang: bug reports"><img src="https://avatars.githubusercontent.com/u/103118149?v=4?s=32" width="32px;" alt="Runfeng Zhang"/></a>
<a href="https://github.com/simas232" title="Simonas Marcišauskas: ideas, planning, & feedback"><img src="https://avatars.githubusercontent.com/u/11994076?v=4?s=32" width="32px;" alt="Simonas Marcišauskas"/></a>
<a href="https://github.com/szkabel" title="szkabel: bug reports"><img src="https://avatars.githubusercontent.com/u/14231791?v=4?s=32" width="32px;" alt="szkabel"/></a>
<a href="https://github.com/tom-hobbs" title="Thomas Hobbs: bug reports"><img src="https://avatars.githubusercontent.com/u/186342818?v=4?s=32" width="32px;" alt="Thomas Hobbs"/></a>
<a href="https://github.com/TunahanCakir" title="TunahanCakir: bug reports"><img src="https://avatars.githubusercontent.com/u/71440332?v=4?s=32" width="32px;" alt="TunahanCakir"/></a>
<a href="https://github.com/HanWeishang" title="Weishang Han: bug reports"><img src="https://avatars.githubusercontent.com/u/104333317?v=4?s=32" width="32px;" alt="Weishang Han"/></a>
<a href="https://github.com/wshao1" title="wshao1: bug reports"><img src="https://avatars.githubusercontent.com/u/50877702?v=4?s=32" width="32px;" alt="wshao1"/></a>
<a href="https://github.com/XuhangLi" title="Xuhang Li: bug reports"><img src="https://avatars.githubusercontent.com/u/41695293?v=4?s=32" width="32px;" alt="Xuhang Li"/></a>
</p>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
