# Human-GEM annotation overhaul — plan

Consolidates issues #964 (metabolite cross-references), #967 (reaction
cross-references), #968 (ModelSEED / MEMOTE), and #626 / #618 (duplicate
metabolites). Genes are handled separately. This is a complete annotation
overhaul, not a merge of the submitted spreadsheets.

## Decisions (agreed)

1. **Protonation matching is pragmatic.** Match a model metabolite to an
   external database entry on skeleton + stereochemistry; accept that the
   database entry may be a different protonation state (most databases store
   the neutral form). Record the database identifier even when protonation
   differs.
2. **Correct wrong annotations**, do not only add. Existing cross-references
   that fail structural verification are corrected or moved to the deprecated
   list, not kept.
3. **No structural changes in this release.** The coming release changes
   annotations only (cross-references, and the SMILES / InChI / formula /
   charge fields of existing metabolites). Adding or removing metabolites and
   reactions, including duplicate merges (#626 / #618), is deferred to the
   following release, which already has several open PRs. Duplicates are
   flagged in this release, not resolved.
4. **Cheminformatics runs locally** (RDKit / OpenBabel), not only in CI.
5. **Correctness first.** Get existing annotations right; grow ModelSEED
   coverage gradually as one of the databases checked, rather than bulk-adding
   it to lift the MEMOTE score.

## Key facts established

- **Protonation reference: dominant ionic microspecies at pH ~7.3** (HMR2
  inheritance). Measured: ATP -4, ADP -3, Pi (HPO4) -2, acetate -1, citrate -3,
  palmitate -1. Only 39% of metabolites are neutral. Formula and charge encode
  this convention.
- **The stored SMILES / InChI are not consistent with formula / charge.** The
  structure file (from #728) holds mostly neutral or arbitrary-protonation
  structures, while formula / charge encode the pH-7.3 microspecies. Example:
  acetate SMILES `CC(=O)O` (neutral, C2H4O2) versus model C2H3O2 / -1. Making
  SMILES, InChI, formula and charge mutually consistent at the pH-7.3
  convention is a core deliverable.
- **Structure coverage:** ~2,973 of 4,163 distinct metabolites (71%) have a
  computable InChIKey and are structure-verifiable. The rest are generic
  R-group lipids and glycans (no InChIKey) and need name / formula / manual
  handling.
- **Submitted matches are skeleton-level.** Independent PubChem verification of
  the #964 additions (n=40): 50% exact, 12% protonation-only, 30% wrong
  stereochemistry, 8% wrong compound. About 38% are not directly usable.
- **Cross-references live only in the .tsv files** (`metabolites.tsv`,
  `reactions.tsv`); the model YAML carries only id / name / formula / charge.
  The .tsv files are the source of truth, exported to SBML/MIRIAM on release,
  which is what MEMOTE scores. MEMOTE improvement is a side effect.
- **Structure-confirmed duplicate metabolites:** 69 groups / 141 metabolites
  by identical full InChIKey (cleaner than the reference-ID lists in
  #626 / #618, which include R-group false positives). Flag only, this release.

## Coverage now

Metabolites (4,163 distinct): MetaNetX 75%, BiGG 63%, KEGG 45%, ChEBI 43%,
PubChem 43%, HMDB 26%, LipidMaps 18%.

Reactions (12,854): Recon3D 96%, HMR2 62%, BiGG 58%, MetaNetX 53%, KEGG 19%,
Rhea 6%, Reactome 2%.

## Principles

1. Structure-first identity. Anchor each metabolite on its InChIKey (with an
   explicit protonation layer) and each reaction on the multiset of participant
   InChIKeys plus stoichiometry. Cross-references hang off that identity.
2. Verify, do not trust. Every cross-reference, existing or new, is classified
   against the target database's own structure: exact / protonation-variant
   (ok) / stereo-different (review) / mismatch (reject) / dead identifier.
3. Add-only for empty cells; correct-or-deprecate for wrong ones. Never blank a
   valid identifier.
4. Move identifiers, do not delete (deprecated lists).
5. Reproducible scripts plus CI gates, not one-off spreadsheets.

## Phases

- **Phase 0 — Foundations.** Document the pH-7.3 convention. Build one canonical
  identity table (id -> InChI, InChIKey layers, formula, charge, standardized
  structure). Download reference structure tables once (MetaNetX chem_prop,
  ChEBI, ModelSEED compounds, LipidMaps, HMDB, PubChem batch) for local joins.
- **Phase 1 — Audit existing annotations (report only).** Run every current
  metabolite and reaction cross-reference through the verifier. Per-namespace
  audit CSVs: confirmed / protonation-variant / wrong-stereo / wrong-compound /
  dead. Quantify how wrong the current annotation is before adding anything.
- **Phase 2 — Structure / formula / charge consistency (this release).** Make
  SMILES, InChI, formula and charge mutually consistent at pH 7.3 for every
  metabolite that has a parseable structure. Correct formula errors; regenerate
  charged InChI / InChIKey. No metabolites added or removed.
- **Phase 3 — Metabolite cross-references (this release).** Apply only verified,
  protonation-aware matches; add-only for empties, correct/deprecate for
  audited-wrong. Feed the filtered good subset of #964 and #968-metabolites.
  Re-derive MetaNetX from structure rather than applying #964's blank column.
  Priority: HMDB / PubChem / ChEBI / LipidMaps, then MetaNetX / ModelSEED, then
  KEGG (name/formula, lower confidence).
- **Phase 4 — Reaction cross-references (this release).** With metabolites clean,
  match reactions by canonical signature against KEGG / BiGG / MetaNetX / Rhea /
  ModelSEED. Apply filtered #967 and #968-reactions. Reconcile EC and pathway
  fields (ties to the EC curation in #1051).
- **Phase 5 — Balance reconciliation.** Re-run mass/charge balance. Fix
  metabolite-formula errors that cause imbalance; separate those from genuinely
  wrong reactions. Enforce proton balance at pH 7.3.
- **Phase 6 — Duplicate report (defer resolution).** Emit the structure-confirmed
  duplicate list for the next (structural) release. Do not merge here.

## Incremental CI (design goal)

Full-model checks are slow. Add fast, diff-scoped checks that validate only what
a pull request changed: for each changed metabolite, structure / formula / charge
consistency and cross-reference structural agreement; for each changed reaction,
mass/charge balance and cross-reference signature. This reuses the delta-vs-base
machinery already in the QC workflows so a change is checked without re-running
the whole model.

## Coordination

- Open PR #1057 corrects malformed and cross-compartment metabolite
  cross-references; overlaps Phase 1 and should be reconciled, not collided.
- Local branch `feat/modelseed-annotation` (#968) reserializes the whole tsv
  (line-ending churn); treat as input, not a base to build on.
- Structural PRs (#653, #1040-#1047, #1011, #1030) belong to the next release.
