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

Static dumps of external sources that were kept here during the Human1 series (SwissProt / Cell Atlas / DeepLoc2 subcellular locations, Rhea reaction associations, and metabolite SMILES/InChI) have been removed for the Human2 release. They had gone stale and are better queried fresh from the source when needed. They remain in the git history and in the last Human1 release, [v1.19.0](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.19.0).
