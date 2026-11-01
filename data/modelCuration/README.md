# Model curation data

This directory contains curation-related data files that support changes to the Human-GEM model and improve the transparency of those changes.

- `metabolite_name_synonyms.tsv` stores mapped synonymous metabolite names used in Human-GEM and other databases.
- `heteromericAaTransporterSubstrateGPRs.tsv` lists the GPRs that `code/modelCuration/pairHeteromericAaTransporterGPRs.py --stage 3` changes to match the SLC7 light chains to the substrates and Na+ coupling of each amino acid exchange (see #1142).
- `deprecatedIdentifierColumns.tsv`: per-entity table of the change "deprecated lists keep the SEED, SMILES and InChI columns" (`code/modelCuration/backfillDeprecatedColumns.py`).
- `transportNameDirection.tsv`: per-entity table of the change "transport names follow the direction of the reaction" (`code/modelCuration/applyCurationTable.py`).
- `geneSymbolsCompartments.tsv`: per-entity table of the change "current HGNC symbols and compartments in genes.tsv" (`code/modelCuration/applyCurationTable.py`).
- `acylCoAStructures.tsv`: per-entity table of the change "acyl-CoA structures carry the configuration of coenzyme A" (`code/modelCuration/applyCurationTable.py`).
- `structureIonisation.tsv`: per-entity table of the change "structures without impossible ionisation (cAMP, cGMP, AICAR, APS)" (`code/modelCuration/applyCurationTable.py`).
- `aminoAcidConfiguration.tsv`: per-entity table of the change "L amino-acid residues and (R)-carnitine in structures drawn as D" (`code/modelCuration/applyCurationTable.py`).
- `structuresOfAnotherCompound.tsv`: per-entity table of the change "structures of metabolites that showed another compound" (`code/modelCuration/applyCurationTable.py`).
- `isomerStructures.tsv`: per-entity table of the change "stereo- and geometric isomers drawn as the isomer the name gives" (`code/modelCuration/applyCurationTable.py`).
- `transBetaOxidationStructures.tsv`: per-entity table of the change "trans beta-oxidation intermediates drawn E, and their names" (`code/modelCuration/applyCurationTable.py`).
- `octadecadienoicAcidStructures.tsv`: per-entity table of the change "(6Z,9Z)-octadecadienoic acid and its carnitine are not linoelaidic" (`code/modelCuration/applyCurationTable.py`).
- `fattyAcidLigaseIds.tsv`: per-entity table of the change "fatty-acid ligases and FATP uptakes carry their own ids" (`code/modelCuration/applyCurationTable.py`).
- `atpSynthaseComplexIIds.tsv`: per-entity table of the change "ATP synthase and complex I carry vectorial ids" (`code/modelCuration/applyCurationTable.py`).
- `acylcarnitineKinaseAcetoneIds.tsv`: per-entity table of the change "C22:4 acylcarnitine, PI4K/PI4P and acetone transport ids" (`code/modelCuration/applyCurationTable.py`).
- `metaboliteIdsOfAnotherCompound.tsv`: per-entity table of the change "metabolite ids that belong to another compound" (`code/modelCuration/applyCurationTable.py`).
- `amideAndPterinStructures.tsv`: per-entity table of the change "amides and 2-amino-4-oxo-pterins drawn as such" (`code/modelCuration/applyCurationTable.py`).
- `crossReferencesAcrossCompartments.tsv`: per-entity table of the change "cross-references shared by all compartments of a metabolite" (`code/modelCuration/applyCurationTable.py`).
- `inchiFromRelease2.tsv`: per-entity table of the change "restore the released v2.0.0 InChIs lost or made less specific" (`code/modelCuration/applyCurationTable.py`).
- `metaNetXSuccessors.tsv`: per-entity table of the change "current MetaNetX ids for deprecated ones" (`code/modelCuration/applyCurationTable.py`).
- `rheaDirectionalIds.tsv`: per-entity table of the change "rxnRheaID holds the directional Rhea id" (`code/modelCuration/applyCurationTable.py`).
- `rxnFromHumanGEM.tsv`: per-entity table of the change "rxnFrom HumanGEM on reactions added in Human-GEM" (`code/modelCuration/applyCurationTable.py`).
- `ecFromGprGene.tsv`: per-entity table of the change "EC numbers taken from a GPR gene that describe another reaction" (`code/modelCuration/applyCurationTable.py`).
- `ecSplitSuccessors.tsv`: per-entity table of the change "split EC numbers keep the successor that fits the reaction" (`code/modelCuration/applyCurationTable.py`).
- `ecWrongPredecessor.tsv`: per-entity table of the change "EC numbers carried from older wrong ones to valid successors" (`code/modelCuration/applyCurationTable.py`).
- `ecNotFitting.tsv`: per-entity table of the change "older EC numbers that do not fit the reaction" (`code/modelCuration/applyCurationTable.py`).
- `ecSubstrateRange.tsv`: per-entity table of the change "EC numbers outside the substrate range of the enzyme" (`code/modelCuration/applyCurationTable.py`).
- `ecPeroxisomalIsomerases.tsv`: per-entity table of the change "duplicate EC on `MAR04633`; peroxisomal enoyl-CoA isomerase ECs" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesFromWrongSource.tsv`: per-entity table of the change "reaction names from a wrong EC, the GPR gene or another enzyme" (`code/modelCuration/applyCurationTable.py`).
- `namesDehydrogenaseBileAcid.tsv`: per-entity table of the change "names and ECs of acyl-CoA dehydrogenase and bile-acid reactions" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesFromChemistry.tsv`: per-entity table of the change "names derived from chemistry that contradict the reaction" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesAfterRenames.tsv`: per-entity table of the change "names after metabolite renames; restore abbreviations" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesRestored.tsv`: per-entity table of the change "restore overwritten reaction names and expand truncated ones" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesForIdentifiers.tsv`: per-entity table of the change "name the reactions whose name is their identifier" (`code/modelCuration/applyCurationTable.py`).
- `exchangeNamesCargoDirection.tsv`: per-entity table of the change "exchange names give the direction of each cargo" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesIdentifying.tsv`: per-entity table of the change "reaction names that do not identify the reaction" (`code/modelCuration/applyCurationTable.py`).
- `reactionNamesFollowRenames.tsv`: per-entity table of the change "reaction names follow the metabolite renames" (`code/modelCuration/applyCurationTable.py`).
- `ecOnTransports.tsv`: per-entity table of the change "remove EC numbers that do not describe a transport" (`code/modelCuration/applyCurationTable.py`).
- `crossReferencesOfAnotherReaction.tsv`: per-entity table of the change "cross-references that describe another reaction" (`code/modelCuration/applyCurationTable.py`).
- `atpDrivenTransportNames.tsv`: per-entity table of the change "name ATP-driven transports apart from passive ones" (`code/modelCuration/applyCurationTable.py`).
- `ecShortChainAcylcarnitines.tsv`: per-entity table of the change "EC of the C4 and C5 acylcarnitine formations" (`code/modelCuration/applyCurationTable.py`).
- `rxnFromHuman1Curation.tsv`: per-entity table of the change "rxnFrom on the reactions added in the Human1 curation" (`code/modelCuration/applyCurationTable.py`).

Static dumps of external sources that were kept here during the Human1 series (SwissProt / Cell Atlas / DeepLoc2 subcellular locations, Rhea reaction associations, and metabolite SMILES/InChI) have been removed for the Human2 release. They had gone stale and are better queried fresh from the source when needed. They remain in the git history and in the last Human1 release, [v1.19.0](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.19.0).
- `acox2BileAcidOxidase.tsv`: per-entity table of the change "ACOX2 oxidises 25(S)THCA-CoA with FAD, and the FADH2 reoxidation has a descriptive name" (`code/modelCuration/applyCurationTable.py`).
- `carnitineAcyltransferaseGprs.tsv`: per-entity table of the change "carnitine acyltransferase GPRs follow the chain lengths each enzyme is active on" (`code/modelCuration/applyCurationTable.py`). The group column names the decision: 2a (mitochondrial C16 to C20), 2d (C2 to C4), Q48 (C21 and C22), Q50 (C24 and C26) and Q47 (ER, C16 and longer).
- `fattyAcidDesaturaseGprs.tsv`: per-entity table of the change "FADS1 and FADS2 rules follow the regiospecificity, FADS3 is removed" (`code/modelCuration/applyCurationTable.py`).
- `diphosphataseGprs.tsv`: per-entity table of the change "PRUNE1 is removed from the diphosphate phosphohydrolases" (`code/modelCuration/applyCurationTable.py`).
- `acox2ReferenceStereoisomer.tsv`: per-entity table of the change "MAR01638 cites the 25S specificity of the trihydroxycoprostanoyl-CoA oxidase" (`code/modelCuration/applyCurationTable.py`).
