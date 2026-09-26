# Model curation data

This directory contains curation-related data files that support changes to the Human-GEM model and improve the transparency of those changes.

- `metabolite_name_synonyms.tsv` stores mapped synonymous metabolite names used in Human-GEM and other databases.
- `metaboliteNameChEBIdiff.tsv` compares each metabolite name with its ChEBI id's preferred name and synonyms (see #1037), categorised to guide name-curation batches.

Static dumps of external sources that were kept here during the Human1 series (SwissProt / Cell Atlas / DeepLoc2 subcellular locations, Rhea reaction associations, and metabolite SMILES/InChI) have been removed for the Human2 release. They had gone stale and are better queried fresh from the source when needed. They remain in the git history and in the last Human1 release, [v1.19.0](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.19.0).
