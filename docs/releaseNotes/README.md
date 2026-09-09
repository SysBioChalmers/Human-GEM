# Release notes

One file per release, named after the version: `docs/releaseNotes/2.1.0.md`.

The file is the body of the GitHub release, published verbatim by the `finalize`
job of `.github/workflows/release.yml`. Write it on `develop` before cutting the
release branch: `increaseHumanGEMVersion.py validate` refuses a version whose notes
file is missing or empty, so the release workflow stops before it creates a branch.

Because the notes are committed, they travel with the release branch and are part
of the release pull request, where they can be reviewed alongside the model changes.

## Structure

The releases published so far open with a short summary of the themes, then list
the changes per pull request. A file that follows that shape:

```markdown
### Short summary:

This release contains three sets of changes:
- Structural model changes: #929, #934
- Non-structural model changes, mostly annotation and naming: #1019, #1029
- Repository maintenance: #1027, #1049

### Longer details

- Fix (PR #929)
  - Recurate the gene-reaction rules where the wrong isoenzymes had been assigned
    - MAR00249, MAR02269, MAR02455
- Feat (PR #995)
  - Introduce MAR20191 as the mitochondrial counterpart of MAR03890
    - Add reaction MAR20191 and gene `ENSG00000135821`
```

List reaction, metabolite and gene identifiers for structural changes. For large
batches (a few hundred entities or more), state the scope and the count instead of
listing every identifier, and point at the pull request for the full diff.
