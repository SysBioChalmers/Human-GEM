# human-GEM: The generic genome-scale metabolic model of _Homo sapiens_

- Brief Model Description

This repository contains the latest version of metabolic model human-GEM, which is a genome-scale model of the generic human cell.

- Abstract:

Human genome-scale metabolic models (GEMs) are important tools for the study of human health and diseases, by providing a scaffold upon which many different types of data can be analyzed in an integrative manner. Despite advancements in the size and complexity of human GEMs, there currently exists the need for a highly accurate, standardized, and manually curated model of human metabolism. Here, we update the previous version of the Human Metabolic Reaction (HMR2) model to human-GEM, in a process that includes the addition of new reactions, metabolites, and genes, followed by extensive curation. The objective of human-GEM is to serve as a community maintained “gold-standard” of generic human GEMs, further enabling integrative and mechanistic studies of human metabolism.

- Model KeyWords:

**GEM Category:** Species; **Utilisation:** Predictive simulation; **Field:** Metabolic-network reconstruction; **Type of Model:** Reconstruction; **Model Source:** HPA, HMR2, iHepatocytes2322, Recon3D; **Omic Source:** Proteomics; **Taxonomy:** _Homo sapiens_; **Metabolic System:** General Metabolism; **Condition:** Generic metabolism;

- Reference: n/a

- Pubmed ID: n/a

- Last update: 2019-04-03


- The model contains:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|
|_Homo sapiens_ |	HMR2|	13520|	10103|	3627|



This repository is administered by Jonathan L. Robinson ([@JonathanRob](https://github.com/jonathanrob)) and Hao Wang ([@Hao-Chalmers](https://github.com/hao-chalmers)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology



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
- [Hao Wang](https://www.chalmers.se/en/staff/Pages/hao-wang.aspx), Chalmers University of Technology, Gothenburg Sweden
