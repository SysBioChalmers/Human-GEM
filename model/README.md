# Human-GEM model and annotation files

This directory contains the Human-GEM model and annotation files.


## Model

The model is available as `.xml`, `.xlsx`, `.txt`, `.yml`, and `.mat`. Note that only the `.yml` version is available on branches other than `main` (e.g., `develop`), to facilitate tracking of model changes.


## Charge and protonation

A metabolite's formula and charge are those of its major microspecies at pH 7.3, the form that MetaNetX and the ChEBI "major microspecies at pH 7.3" entries give: carboxylic acids and phosphates are ionised (pyruvate C3H3O3 −1, ATP −4, acyl-CoAs −4), amines are protonated (choline +1, lysine +1) and amino acids are zwitterions. Reactions are balanced for mass and charge with H+ in these forms. When a new metabolite is added, or a formula or charge changed, take the pH 7.3 form.

Some metabolites do not follow this yet, for example the glycine- and taurine-conjugated bile acids and the steroid glucuronides, which are written as neutral acids. The cross-references and structures in `metabolites.tsv` can name the neutral parent compound instead of the pH 7.3 form (for example a ChEBI id of the acid rather than of its anion), so a charge read from a cross-reference is not proof of the model charge. Ammonia appears both as NH3 and as NH4+ (`MAM02578`, `MAM02579`).


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
rxnSeedID      |ModelSEED reaction ID   |
rxnHMR2ID      |HMR2 reaction ID        |
rxnRatconID    |Ratcon reaction ID      |
rxnTCDBID      |TCDB ID                 |
spontaneous    |1 if the reaction needs no enzyme: non-enzymatic chemistry, or diffusion (transport between cytosol and nucleus through nuclear pores, and gases and water across membranes). ftINIT always keeps these reactions in context-specific models.|
rxnRheaID      |Rhea ID of the reaction's direction: left-to-right or right-to-left for an irreversible reaction, bidirectional for a reversible one|
rxnRheaMasterID|Rhea master ID, the undirected reaction of rxnRheaID|
rxnRetired     |Earlier identifiers of the reaction (for example `HMR_` IDs, or a `MAR` ID merged into this reaction), separated by `;`|


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
metSmiles      |SMILES structure         |
metInChI       |InChI structure          |
metSeedID      |ModelSEED compound ID    |
metHMR2ID      |HMR2 metabolite ID       |
metRetired     |Earlier identifiers of the metabolite (for example `m` IDs from HMR, or a `MAM` ID merged into this metabolite), separated by `;`|


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
compartments  |Subcellular locations of the gene product, separated by `;`|
compDataSource|Sources of the subcellular locations (for example SwissProt, DeepLoc2, CellAtlas), separated by `;`|


To import/export this annotation data to/from MATLAB, use the `importTsvFile` and `exportTsvFile` functions, respectively.

All three tables are tab-separated, and a field is enclosed in double quotes (`"`) only when its value contains a tab, a double quote or a newline. No field currently does, so the files carry no quotes. Every column holds text, including the all-digit identifier columns such as `geneEntrezID`.
