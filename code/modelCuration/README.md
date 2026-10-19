# Model curation scripts

This directory holds scripts used to make curation changes to Human-GEM when the number of changes is too large to view or reproduce practically. Keeping the scripts here improves the transparency of how the model was changed.

The one-off curation scripts from the Human1 series (2018-2024) have been removed for the Human2 release to keep the repository lean. They remain in the git history and in the last Human1 release, [v1.19.0](https://github.com/SysBioChalmers/Human-GEM/releases/tag/v1.19.0). New curation scripts should be added here.

## Scripts

Each script edits the model files in place and takes `--dry-run`. The per-entity table of a scripted change is in `data/modelCuration/`, and the commit that made the change gives the command.
- `backfillDeprecatedColumns.py --table data/modelCuration/<name>.tsv` gives `data/deprecatedIdentifiers/deprecatedReactions.tsv` and `deprecatedMetabolites.tsv` the columns of `model/reactions.tsv` and `model/metabolites.tsv` that they lack, and fills cells from a table (list, id, column, value). `--dry-run` reports what would change.
