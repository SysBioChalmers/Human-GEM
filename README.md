# Human-GEM: The generic genome-scale metabolic model of _Homo sapiens_

- Brief Model Description

This repository contains the latest version of Human-GEM, a human genome-scale metabolic model.

- Abstract:

Genome-scale metabolic models (GEMs) are valuable tools to study metabolism, and provide a scaffold for integrative analysis of -omics data. Researchers have developed increasingly comprehensive human GEMs, but the disconnect among different model sources and versions impedes further progress. We therefore integrated and extensively curated the most recent human metabolic models to construct a consensus GEM, Human-GEM. Human-GEM was created using a version-controlled, open-source model development framework to enable community-driven curation and refinement. This framework allows Human-GEM to be an evolving shared resource for future studies of human health and disease.

- Model KeyWords:

**GEM Category:** Species; **Utilisation:** Predictive simulation; **Field:** Metabolic-network reconstruction; **Type of Model:** Reconstruction; **Model Source:** HPA, HMR2, iHsa, iHepatocytes2322, Recon3D; **Omic Source:** Proteomics; **Taxonomy:** _Homo sapiens_; **Metabolic System:** General Metabolism; **Condition:** Generic metabolism;

- Reference:

Article under consideration.

- Pubmed ID: n/a

- Last update: 2019-06-12


- The model contains:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|
|_Homo sapiens_ |	HMR2, Recon3D|	13520|	10103|	3628|




This repository is administered by Jonathan L. Robinson ([@JonathanRob](https://github.com/jonathanrob)) and Hao Wang ([@Hao-Chalmers](https://github.com/hao-chalmers)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.



## Installation

### Required Software:
* A functional Matlab installation (MATLAB 7.3 and higher).
* The [RAVEN toolbox](https://github.com/SysBioChalmers/RAVEN).
* The [COBRA toolbox](https://github.com/opencobra/cobratoolbox).


### Dependencies - Recommended Software:
* The libSBML MATLAB API (version [5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/) is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.


### Installation Instructions
* Clone the model from [master](https://github.com/SysBioChalmers/) branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers)
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)


## Contributors
- [Jonathan L. Robinson](https://www.chalmers.se/en/Staff/Pages/jonrob.aspx), Chalmers University of Technology, Gothenburg Sweden
- [Pınar Kocabaş](https://www.chalmers.se/en/staff/Pages/kocabas.aspx), Chalmers University of Technology, Gothenburg Sweden
- [Pierre-Etienne Cholley](https://www.chalmers.se/en/staff/Pages/cholley.aspx), Chalmers University of Technology, Gothenburg Sweden
- [Avlant Nilsson](https://www.chalmers.se/en/staff/Pages/avlant-nilsson.aspx), Chalmers University of Technology, Gothenburg Sweden
- [Hao Wang](https://www.chalmers.se/en/staff/Pages/hao-wang.aspx), Chalmers University of Technology, Gothenburg Sweden
