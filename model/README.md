# Human-GEM model and annotation files

This directory contains the Human-GEM model and annotation files.


## Model

The model is available as `.xml`, `.xlsx`, `.txt`, `.yml`, and `.mat`. Note that only the `.yml` version is available on branches other than `main` (e.g., `develop`), to facilitate tracking of model changes.


## Reaction, Metabolite, and Gene Annotations

Additional annotation information and external identifiers for Human-GEM reactions, metabolites, and genes are provided as `tsv` files. The structure of the `tsv` files are tabulated below.

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
geneENSPID    |Ensembl protein ID    |
geneUniProtID |UniProt ID            |
geneSymbols   |Gene Symbol           |
geneEntrezID  |NCBI Entrez ID        |
geneNames     |Gene Name             |
geneAliases   |Alias Names           |


To import/export this annotation data to/from MATLAB, use the `importTsvFile` and `exportTsvFile` functions, respectively.
