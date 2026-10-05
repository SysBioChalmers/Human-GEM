"""Regenerate an animal GEM (Mouse, Rat, Worm, Fruitfly, Zebrafish) from Human-GEM.

Python/raven-toolbox replacement for the ``masterScript<Species>GEM.m`` scripts that
each Animal-GEM repository carried (see Human-GEM issue #421). The five scripts differed
only in the species name, so the species is an argument here and everything
species-specific is read from the Animal-GEM repository:

    <repo>/data/human2<Species>Orthologs.tsv    Alliance of Genome Resources orthologs
    <repo>/data/<species>SpecificRxns.tsv        reactions that Human-GEM does not have
    <repo>/data/<species>SpecificMets.tsv        metabolites that those reactions need
    <repo>/model/reactions.tsv, metabolites.tsv  annotation tables (rewritten)
    <repo>/model/<Species>-GEM.{yml,mat,xml}    the model (rewritten)
    <repo>/version.txt                           the model version

Steps, in order:

1. Take Human-GEM as the template and rewrite every GPR through the ortholog table
   (Ensembl id -> human symbol -> species gene). Reactions whose GPR becomes empty are
   removed.
2. Add the species-specific metabolites and reactions.
3. Gap-fill the essential metabolic tasks from Human-GEM, so the model can grow on Ham's
   medium. Gap-filled reactions carry no GPR and a note saying so.
4. Stamp the version and date into the model metadata, merge the annotation tables, and
   write the model files.

The script reads and writes nothing outside the two repositories and takes no input from
the terminal, so a developer and a CI job run the same command:

    python code/animalGEM/generateAnimalGEM.py Mouse
    python code/animalGEM/generateAnimalGEM.py Mouse --repo ../Mouse-GEM --version 1.9.0

The Animal-GEM repository defaults to ``Mouse-GEM`` next to this repository; ``--repo`` or
``ANIMAL_GEM_REPO`` points elsewhere (a CI checkout, say). The gap-filling MILP wants a
real MILP solver on a genome-scale model: ``--solver gurobi`` (or any solver that cobra
knows) selects one, otherwise cobra's default is used.

Anything that is not step 4's version/date is deterministic for a given template, input
tables, solver and seed.
"""
from __future__ import annotations

import argparse
import ast
import datetime
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import cobra
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

# Make the sibling code/annotateGEM.py importable regardless of the caller's cwd.
sys.path.insert(0, str(REPO_ROOT / "code"))
from annotateGEM import annotate_gem  # noqa: E402

from raven_toolbox.init import fill_tasks  # noqa: E402
from raven_toolbox.io import export_for_git, read_yaml_model  # noqa: E402
from raven_toolbox.manipulation.add import add_reactions_from_equations  # noqa: E402
from raven_toolbox.tasks import check_tasks, parse_task_list  # noqa: E402

ESSENTIAL_TASKS = REPO_ROOT / "data" / "metabolicTasks" / "metabolicTasks_Essential.txt"

# NCBI taxonomy ids, for the model metadata.
TAXONOMY = {
    "Mouse": "10090",
    "Rat": "10116",
    "Worm": "6239",
    "Fruitfly": "7227",
    "Zebrafish": "7955",
}

# Human-GEM reactions the gap-filling reset relies on: the human biomass reaction is
# blocked and the objective moved to the generic cell components, which every
# eukaryote is assumed to make.
HUMAN_BIOMASS_RXN = "MAR13082"
BIOMASS_COMPONENTS_RXN = "MAR00021"

GAPFILL_NOTE = "reaction added by gap filling"

# Human-GEM releases in which MAR00021 cannot carry flux on Ham's medium (issue #1140):
# its cofactor pool MAR00022 consumes [protein]-N6-(lipoyl)lysine, which nothing makes
# since the lipoylation curation. Remove this, and fix_lipoyl_biomass, once MAR00021 is
# curated in Human-GEM.
LIPOYL_BIOMASS_RELEASES = {"2.0.0", "2.0.1", "2.1.0"}
LIPOYL_LYSINE = "[protein]-N6-(lipoyl)lysine"


@dataclass(frozen=True)
class AnimalRepo:
    """Where the species-specific inputs and outputs of one Animal-GEM repository live."""

    species: str
    root: Path

    @property
    def model_id(self) -> str:
        return f"{self.species}-GEM"

    @property
    def data_dir(self) -> Path:
        return self.root / "data"

    @property
    def model_dir(self) -> Path:
        return self.root / "model"

    @property
    def orthologs(self) -> Path:
        return self.data_dir / f"human2{self.species}Orthologs.tsv"

    @property
    def specific_rxns(self) -> Path:
        return self.data_dir / f"{self.species.lower()}SpecificRxns.tsv"

    @property
    def specific_mets(self) -> Path:
        return self.data_dir / f"{self.species.lower()}SpecificMets.tsv"

    @property
    def version_file(self) -> Path:
        return self.root / "version.txt"


def read_tsv(path: Path) -> pd.DataFrame:
    """Read a tab-separated table as text; empty cells become ``""``."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


# --- orthologs ---------------------------------------------------------------------

def read_alliance_orthologs(path: Path, count_best: bool = True) -> pd.DataFrame:
    """Reduce an Alliance of Genome Resources ortholog table to human/species symbol pairs.

    1. With ``count_best``, drop pairs that are neither best forward nor best reverse.
    2. Keep every human gene with a single remaining hit.
    3. For the others, keep the hits that are both best forward and best reverse.
    4. If none is, keep the hit supported by the most methods (the first on a tie).

    Returns a frame with columns ``from`` (human symbol) and ``to`` (species symbol).
    """
    table = read_tsv(path)
    if count_best:
        table = table[~((table["best"] == "No") & (table["bestReverse"] == "No"))]
    table = table.assign(methodCount=pd.to_numeric(table["methodCount"], errors="coerce"))

    keep = []
    for _, hits in table.groupby("fromGeneId", sort=False):
        if len(hits) == 1:
            keep.append(hits)
            continue
        both = hits[(hits["best"] == "Yes") & (hits["bestReverse"] == "Yes")]
        keep.append(both if len(both) else hits.sort_values("methodCount", ascending=False,
                                                            kind="stable").head(1))
    kept = pd.concat(keep) if keep else table.iloc[0:0]
    pairs = kept[["fromSymbol", "toSymbol"]].rename(columns={"fromSymbol": "from", "toSymbol": "to"})
    return pairs.reset_index(drop=True)


def ensembl_to_species_genes(genes_tsv: Path, orthologs: pd.DataFrame) -> dict[str, list[str]]:
    """Map each Human-GEM gene (Ensembl id) to its species orthologs, via the gene symbol."""
    by_symbol: dict[str, list[str]] = {}
    for human, species in zip(orthologs["from"], orthologs["to"]):
        if species not in by_symbol.setdefault(human, []):
            by_symbol[human].append(species)

    mapping: dict[str, list[str]] = {}
    for gene, symbols in zip(*(read_tsv(genes_tsv)[c] for c in ("genes", "geneSymbols"))):
        found: list[str] = []
        for symbol in (s.strip() for s in symbols.split(";") if s.strip()):
            found.extend(g for g in by_symbol.get(symbol, []) if g not in found)
        if found:
            mapping[gene] = found
    return mapping


# --- GPR rewriting -----------------------------------------------------------------

def _map_node(node: ast.AST, mapping: dict[str, list[str]]):
    """Rewrite a GPR node. Returns ``None`` (nothing left), a gene id, or ``(op, children)``."""
    if isinstance(node, ast.Name):
        genes = mapping.get(node.id, [])
        if not genes:
            return None
        return genes[0] if len(genes) == 1 else ("or", list(genes))
    if isinstance(node, ast.BoolOp):
        op = "and" if isinstance(node.op, ast.And) else "or"
        children = [c for c in (_map_node(v, mapping) for v in node.values) if c is not None]
        return _collapse(op, children)
    raise ValueError(f"unsupported GPR element: {ast.dump(node)}")


def _collapse(op: str, children: list):
    """Flatten nested same-operator groups, drop duplicates, unwrap single children."""
    flat: list = []
    for child in children:
        parts = child[1] if isinstance(child, tuple) and child[0] == op else [child]
        flat.extend(p for p in parts if p not in flat)
    if not flat:
        return None
    return flat[0] if len(flat) == 1 else (op, flat)


def _render(expr, parent: str | None = None) -> str:
    if expr is None:
        return ""
    if isinstance(expr, str):
        return expr
    op, children = expr
    text = f" {op} ".join(_render(c, op) for c in children)
    return f"({text})" if parent is not None and parent != op else text


def map_gene_reaction_rule(rule: str, mapping: dict[str, list[str]]) -> str:
    """Rewrite one GPR through ``mapping`` (gene -> list of replacement genes).

    A gene with several replacements becomes an ``or`` of them. A gene with none is
    dropped from its group, whether ``and`` or ``or``, so a complex keeps the subunits
    that have an ortholog. Duplicate genes and redundant parentheses are removed. A rule
    with no gene left comes back as ``""``.
    """
    if not rule.strip():
        return ""
    return _render(_map_node(cobra.core.gene.GPR.from_string(rule).body, mapping))


def build_ortholog_draft(template: cobra.Model, gene_map: dict[str, list[str]]) -> cobra.Model:
    """Copy of ``template`` with every GPR rewritten through ``gene_map``.

    Reactions that had a GPR and lose all of it are removed, together with metabolites
    only they used. Reactions that never had a GPR (spontaneous, exchange) stay.
    """
    draft = template.copy()
    draft.id = ""
    draft.name = ""
    draft.notes = {}
    for entity in (*draft.reactions, *draft.metabolites):
        entity.notes.pop("rxnFrom", None)
        entity.notes.pop("metFrom", None)

    lost = []
    for rxn in draft.reactions:
        old = rxn.gene_reaction_rule
        if not old:
            continue
        new = map_gene_reaction_rule(old, gene_map)
        if new:
            rxn.gene_reaction_rule = new
        else:
            lost.append(rxn)
    draft.remove_reactions(lost, remove_orphans=True)
    prune_unused_genes(draft)
    return draft


def prune_unused_genes(model: cobra.Model) -> None:
    """Remove genes that no reaction uses."""
    unused = [g for g in model.genes if not g.reactions]
    if unused:
        cobra.manipulation.remove_genes(model, unused, remove_reactions=False)


# --- species-specific network ------------------------------------------------------

def _number(value: str, default: float) -> float:
    return default if value.strip() in ("", "NaN", "nan") else float(value)


def add_species_network(model: cobra.Model, rxns: pd.DataFrame, mets: pd.DataFrame) -> list[str]:
    """Add the species-specific metabolites and reactions; returns the new reaction ids.

    ``mets`` needs ``mets, metNames, metFormulas, metCharges, compartments``; ``rxns``
    needs ``rxns, equations, subSystems, grRules`` and may have ``lb, ub, rxnNames,
    eccodes, rxnReferences, rxnConfidenceScores``. Equations name their metabolites as
    ``name[compartment]``. Ids that the model already has are an error rather than an
    overwrite.
    """
    for table, needed in ((mets, ("mets", "metNames", "metFormulas", "metCharges", "compartments")),
                          (rxns, ("rxns", "equations", "subSystems", "grRules"))):
        missing = [c for c in needed if c not in table.columns]
        if missing:
            raise ValueError(f"species-specific table is missing column(s) {missing}")
    clash = [m for m in mets["mets"] if m in model.metabolites]
    clash += [r for r in rxns["rxns"] if r in model.reactions]
    if clash:
        raise ValueError(f"already in the model, cannot be added: {clash[:10]}")
    if not rxns["equations"].map(lambda e: "[" in e and "]" in e).all():
        raise ValueError('equation metabolites must be written "name[compartment]"')

    model.add_metabolites([
        cobra.Metabolite(row.mets, name=row.metNames, formula=row.metFormulas,
                         charge=int(_number(row.metCharges, 0)), compartment=row.compartments)
        for row in mets.itertuples()
    ])

    specs = []
    for row in rxns.to_dict("records"):
        spec = {"id": row["rxns"], "equation": row["equations"],
                "gene_reaction_rule": row["grRules"]}
        if row.get("rxnNames"):
            spec["name"] = row["rxnNames"]
        if row.get("lb", "") != "" or row.get("ub", "") != "":
            spec["bounds"] = (_number(row.get("lb", ""), -1000.0), _number(row.get("ub", ""), 1000.0))
        specs.append(spec)
    added = add_reactions_from_equations(model, specs, mets_by="name", allow_new_mets=False)

    for rxn, row in zip(added, rxns.to_dict("records")):
        codes = [c.strip() for c in row.get("eccodes", "").split(";") if c.strip()]
        if codes:
            rxn.annotation["ec-code"] = codes
        # Subsystems are a list on every Human-GEM reaction; keep the new ones alike.
        rxn.subsystem = [row["subSystems"]] if row["subSystems"] else []
        if row.get("rxnReferences"):
            rxn.notes["references"] = row["rxnReferences"]
        if row.get("rxnConfidenceScores", "") != "":
            rxn.notes["confidence_score"] = _number(row["rxnConfidenceScores"], float("nan"))
    return [r.id for r in added]


# --- gap-filling -------------------------------------------------------------------

def fix_lipoyl_biomass(model: cobra.Model) -> bool:
    """Let the cofactor pool of ``MAR00021`` use lipoic acid instead of lipoyl-lysine.

    In place; returns whether anything changed. ``MAR00022`` makes the "cofactors and
    vitamins" pool that ``MAR00021`` consumes. It takes ``[protein]-N6-(lipoyl)lysine``,
    which Human-GEM 2.0.0-2.1.0 cannot make, so the generic biomass is blocked. The
    current human cofactor pool (``MAR10065``) takes free lipoic acid, which Ham's medium
    supplies; ``MAR00022`` is changed to match. Only the template of a model in
    ``LIPOYL_BIOMASS_RELEASES`` is passed here; the curated fix belongs in Human-GEM.
    """
    if "MAR00022" not in model.reactions or "MAR10065" not in model.reactions:
        return False
    lipoyl = [m for m in model.reactions.MAR00022.metabolites if m.name == LIPOYL_LYSINE]
    acid = [m for m in model.reactions.MAR10065.metabolites if m.name == "lipoic acid"]
    if len(lipoyl) != 1 or len(acid) != 1:
        return False
    coefficient = model.reactions.MAR00022.metabolites[lipoyl[0]]
    model.reactions.MAR00022.add_metabolites({lipoyl[0]: -coefficient, acid[0]: coefficient})
    return True


def reset_biomass(model: cobra.Model) -> None:
    """Block the human biomass reaction and make the generic cell components the objective."""
    for rid in (HUMAN_BIOMASS_RXN, BIOMASS_COMPONENTS_RXN):
        if rid not in model.reactions:
            raise ValueError(f"{rid} is not in the model; is the template Human-GEM?")
    model.reactions.get_by_id(HUMAN_BIOMASS_RXN).bounds = (0, 0)
    components = model.reactions.get_by_id(BIOMASS_COMPONENTS_RXN)
    components.upper_bound = 1000
    model.objective = {components: 1}


def fill_essential_tasks(model: cobra.Model, reference: cobra.Model, tasks, *,
                         reset: bool = True, time_limit: float | None = None) -> list[str]:
    """Gap-fill ``model`` in place from ``reference`` until the essential tasks pass.

    Returns the ids of the added reactions. They keep no GPR (the reference genes are not
    the species' genes) and are marked in their notes. Raises if a task still fails.
    """
    if len(set(model.reactions.list_attr("id")) & set(reference.reactions.list_attr("id"))) \
            < 0.5 * len(model.reactions):
        raise ValueError("the model shares under half of its reactions with the reference; "
                         "is the reference the template it was built from?")
    reference = reference.copy()
    if reset:
        reset_biomass(model)
        reset_biomass(reference)

    options = {} if time_limit is None else {"time_limit": time_limit}
    filled = fill_tasks(model, reference, tasks, **options)
    if filled.failed_tasks:
        hint = (f" The reference cannot run them with {BIOMASS_COMPONENTS_RXN} as the biomass "
                "reaction (Human-GEM issue #1140); --keep-biomass keeps the human biomass instead."
                if reset else "")
        raise RuntimeError(f"gap-filling could not satisfy tasks: {sorted(set(filled.failed_tasks))}.{hint}")

    new = [rid for rid in filled.added_reactions if rid not in model.reactions]
    _add_reactions_from(model, filled.model, new)
    for rid in new:
        rxn = model.reactions.get_by_id(rid)
        rxn.gene_reaction_rule = ""
        note = rxn.notes.get("note", "").rstrip(";")
        rxn.notes["note"] = f"{note};{GAPFILL_NOTE}" if note else GAPFILL_NOTE
    prune_unused_genes(model)

    failed = [r.id for r in check_tasks(model, tasks, close_boundaries=True) if not r.passed]
    if failed:
        raise RuntimeError(f"tasks still fail after gap-filling: {failed}")
    return new


def _add_reactions_from(target: cobra.Model, source: cobra.Model, rxn_ids: list[str]) -> None:
    """Add ``rxn_ids`` from ``source`` to ``target`` along with the metabolites they need."""
    mets = {m.id: m for rid in rxn_ids for m in source.reactions.get_by_id(rid).metabolites}
    target.add_metabolites([m.copy() for mid, m in mets.items() if mid not in target.metabolites])
    for rid in rxn_ids:
        src = source.reactions.get_by_id(rid)
        rxn = cobra.Reaction(src.id, name=src.name, subsystem=src.subsystem,
                             lower_bound=src.lower_bound, upper_bound=src.upper_bound)
        target.add_reactions([rxn])
        rxn.add_metabolites({target.metabolites.get_by_id(m.id): c for m, c in src.metabolites.items()})
        rxn.annotation = dict(src.annotation)
        rxn.notes = dict(src.notes)


# --- annotation tables -------------------------------------------------------------

def merge_annotation(human: pd.DataFrame, species: pd.DataFrame, id_col: str,
                     ids: list[str]) -> pd.DataFrame:
    """Annotation table of the animal model: Human-GEM's columns, species rows appended.

    Species-specific columns Human-GEM has no counterpart for are not carried over, and
    Human-GEM columns the species table lacks stay empty for its rows. Rows follow
    ``ids``; an id that neither table knows is an error.
    """
    extra = species.reindex(columns=human.columns, fill_value="")
    merged = pd.concat([human, extra], ignore_index=True).drop_duplicates(id_col, keep="last")
    unknown = sorted(set(ids) - set(merged[id_col]))
    if unknown:
        raise ValueError(f"model components with no annotation row: {unknown[:10]}")
    return merged.set_index(id_col).loc[ids].reset_index()[list(human.columns)]


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, quoting=0, lineterminator="\n")


# --- metadata and output -----------------------------------------------------------

def stamp_metadata(model: cobra.Model, repo: AnimalRepo, version: str, date: str) -> None:
    """Set the model id, name, version, date, taxonomy and source URL."""
    model.id = model.name = repo.model_id
    model.notes = {
        "version": version,
        "metaData": {
            "version": version,
            "date": date,
            "taxonomy": TAXONOMY.get(repo.species, ""),
            "sourceUrl": f"https://github.com/SysBioChalmers/{repo.model_id}",
        },
    }


def write_outputs(model: cobra.Model, repo: AnimalRepo, rxn_table: pd.DataFrame,
                  met_table: pd.DataFrame, formats: tuple[str, ...]) -> None:
    """Write the annotation tables and the model files into ``repo.model_dir``."""
    out = repo.model_dir
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(rxn_table, out / "reactions.tsv")
    write_tsv(met_table, out / "metabolites.tsv")

    plain = [f for f in formats if f in ("yml", "mat")]
    annotated = [f for f in formats if f not in ("yml", "mat")]
    if plain:
        export_for_git(model, out, prefix=repo.model_id, formats=plain, sub_dirs=False,
                       varname=f"{repo.species.lower()}GEM")
    if annotated:
        # SBML ids cannot contain a dash. Cross-references and SBO terms are merged into
        # the SBML/Excel/txt copies only, as in Human-GEM.
        copy = annotate_gem(model.copy(), out, types=("rxn", "met"))
        copy.id = repo.model_id.replace("-", "")
        export_for_git(copy, out, prefix=repo.model_id, formats=annotated, sub_dirs=False)


# --- driver ------------------------------------------------------------------------

def generate_animal_gem(species: str, repo_dir: Path, *, version: str | None = None,
                        date: str | None = None, human_model: Path | None = None,
                        formats: tuple[str, ...] = ("yml", "mat", "xml"),
                        time_limit: float | None = None,
                        reset_objective: bool = True) -> cobra.Model:
    """Generate ``species``-GEM from Human-GEM and write it into ``repo_dir``.

    ``version`` is written to ``version.txt`` and the model metadata; without it the
    repository's current ``version.txt`` is kept and stamped. ``date`` defaults to today.
    ``reset_objective`` swaps the human biomass reaction for the generic cell components
    before gap-filling (see :func:`reset_biomass`); without it Human-GEM's own biomass
    reaction stays the objective.
    """
    repo = AnimalRepo(species, Path(repo_dir).resolve())
    for path in (repo.orthologs, repo.specific_rxns, repo.specific_mets):
        if not path.is_file():
            raise FileNotFoundError(path)
    if version is None:
        if not repo.version_file.is_file():
            raise FileNotFoundError(f"{repo.version_file}; pass --version")
        version = repo.version_file.read_text(encoding="utf-8").strip()
    date = date or datetime.date.today().isoformat()
    human_model = Path(human_model or REPO_ROOT / "model" / "Human-GEM.yml")

    print(f"Reading template {human_model}", flush=True)
    template = read_yaml_model(human_model)
    template_version = (template.notes.get("metaData") or {}).get("version")
    if reset_objective and template_version in LIPOYL_BIOMASS_RELEASES \
            and fix_lipoyl_biomass(template):
        print(f"Human-GEM {template_version}: MAR00022 uses lipoic acid instead of "
              f"{LIPOYL_LYSINE} (see issue #1140)", flush=True)

    orthologs = read_alliance_orthologs(repo.orthologs)
    gene_map = ensembl_to_species_genes(REPO_ROOT / "model" / "genes.tsv", orthologs)
    print(f"{len(orthologs)} ortholog pairs, {len(gene_map)} Human-GEM genes with an ortholog", flush=True)

    model = build_ortholog_draft(template, gene_map)
    print(f"Ortholog draft: {len(model.reactions)} of {len(template.reactions)} reactions", flush=True)

    rxns, mets = read_tsv(repo.specific_rxns), read_tsv(repo.specific_mets)
    added = add_species_network(model, rxns, mets)
    print(f"Added {len(added)} species-specific reactions, {len(mets)} metabolites", flush=True)

    filled = fill_essential_tasks(model, template, parse_task_list(ESSENTIAL_TASKS),
                                  reset=reset_objective, time_limit=time_limit)
    print(f"Gap-filled {len(filled)} reactions", flush=True)

    for gene in model.genes:
        gene.name = gene.id
    stamp_metadata(model, repo, version, date)

    rxn_table = merge_annotation(read_tsv(REPO_ROOT / "model" / "reactions.tsv"), rxns, "rxns",
                                 [r.id for r in model.reactions])
    met_table = merge_annotation(read_tsv(REPO_ROOT / "model" / "metabolites.tsv"), mets, "mets",
                                 [m.id for m in model.metabolites])
    write_outputs(model, repo, rxn_table, met_table, formats)
    if not repo.version_file.is_file() or repo.version_file.read_text(encoding="utf-8").strip() != version:
        repo.version_file.write_text(version, encoding="utf-8")
    print(f"{repo.model_id} {version} written to {repo.model_dir}", flush=True)
    return model


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Regenerate an animal GEM from Human-GEM and its species-specific data.")
    parser.add_argument("species", choices=sorted(TAXONOMY), help="species, e.g. Mouse")
    parser.add_argument("--repo", type=Path,
                        default=os.environ.get("ANIMAL_GEM_REPO"),
                        help="Animal-GEM repository (default: ../<Species>-GEM, or $ANIMAL_GEM_REPO)")
    parser.add_argument("--version", help="new model version, also written to version.txt "
                                          "(default: keep the repository's version.txt)")
    parser.add_argument("--date", help="model date, YYYY-MM-DD (default: today)")
    parser.add_argument("--human-model", type=Path, help="template (default: model/Human-GEM.yml)")
    parser.add_argument("--formats", default="yml,mat,xml",
                        help="comma-separated output formats: yml, mat, xml, xlsx, txt")
    parser.add_argument("--solver", help="cobra solver for the gap-filling MILP (default: cobra's)")
    parser.add_argument("--time-limit", type=float, help="seconds per gap-filling MILP")
    parser.add_argument("--keep-biomass", action="store_true",
                        help="keep Human-GEM's biomass reaction instead of resetting the objective "
                             f"to the generic cell components ({BIOMASS_COMPONENTS_RXN})")
    args = parser.parse_args(argv)

    if args.solver:
        cobra.Configuration().solver = args.solver
    repo = args.repo or REPO_ROOT.parent / f"{args.species}-GEM"
    generate_animal_gem(args.species, repo, version=args.version, date=args.date,
                        human_model=args.human_model, time_limit=args.time_limit,
                        reset_objective=not args.keep_biomass,
                        formats=tuple(f.strip() for f in args.formats.split(",") if f.strip()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
