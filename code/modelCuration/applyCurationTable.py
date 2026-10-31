"""Apply a per-entity table of annotation changes to the model files.

    python code/modelCuration/applyCurationTable.py --table data/modelCuration/<name>.tsv [--dry-run]

The table is tab-separated with the columns entity, field, before and after, one row per
changed value; further columns (a note, a source) are ignored. entity is a reaction
(MAR), metabolite (MAM) or gene (ENSG) id, and field is one of:

  - name (reactions, metabolites and genes), formula and charge (metabolites), ec-code and
    rxnFrom (reactions), in model/Human-GEM.yml; several EC numbers are joined by ";";
  - a column of model/reactions.tsv, model/metabolites.tsv or model/genes.tsv.

An empty value stands for an absent field or an empty cell. A row is applied when the
current value equals before. A value that already equals after is left as it is; any other
value is reported and nothing is written. The files are edited in place, line by line, so
every line that the table does not change stays byte-identical. --dry-run reports what
would change and writes nothing.
"""
import argparse
import csv
import io
import re
import sys
from collections import Counter
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

YML = Path("model/Human-GEM.yml")
TSV = {"MAR": Path("model/reactions.tsv"), "MAM": Path("model/metabolites.tsv"),
       "ENSG": Path("model/genes.tsv")}
SECTION = {"MAR": "reactions", "MAM": "metabolites", "ENSG": "genes"}
YAML_FIELDS = {"MAR": {"name", "ec-code", "rxnFrom"}, "MAM": {"name", "formula", "charge"},
               "ENSG": {"name"}}
# Field order of a reaction entry, as raven_toolbox.io.write_yaml_model writes it.
REACTION_ORDER = ["id", "name", "metabolites", "lower_bound", "upper_bound", "gene_reaction_rule",
                  "objective_coefficient", "subsystem", "annotation", "references", "rxnFrom",
                  "confidence_score", "notes"]

# The scalars that write_yaml_model quotes (raven_toolbox.io.yaml._needs_quote): those
# that would read back as another type, and those that are not valid plain scalars.
_RESOLVES_TO_OTHER_TYPE = re.compile(
    r"^(true|True|TRUE|false|False|FALSE)$"
    r"|^[-+]?[0-9][0-9_]*\.[0-9_]*([eE][-+]?[0-9]+)?$"
    r"|^[-+]?[0-9][0-9_]*[eE][-+]?[0-9]+$"
    r"|^[-+]?\.[0-9_]+([eE][-+][0-9]+)?$"
    r"|^[-+]?\.(inf|Inf|INF)$"
    r"|^\.(nan|NaN|NAN)$"
    r"|^[-+]?0b[01_]+$"
    r"|^[-+]?0o?[0-7_]+$"
    r"|^[-+]?[0-9_]+$"
    r"|^[-+]?0x[0-9a-fA-F_]+$"
    r"|^(~|[Nn]ull|NULL)$"
    r"|^[0-9]{4}-[0-9]{2}-[0-9]{2}$"
)
_WS = (" ", "\t", "\n", "\r")


def needs_quote(text: str) -> bool:
    if not text or _RESOLVES_TO_OTHER_TYPE.match(text):
        return True
    first = text[0]
    if first in "#,[]{}&*!|>'\"%@`":
        return True
    followed_by_ws = len(text) == 1 or text[1] in _WS
    if first in "?:-" and followed_by_ws:
        return True
    if first in _WS or text[-1] in _WS or "\n" in text or "\r" in text:
        return True
    for i in range(1, len(text)):
        if text[i] == ":" and (i == len(text) - 1 or text[i + 1] in _WS):
            return True
        if text[i] == "#" and text[i - 1] in _WS:
            return True
    return False


_DUMPER = YAML(typ="rt")
_DUMPER.width = 1_000_000
_LOADER = YAML(typ="safe")


def render(value: str) -> str:
    """A scalar as write_yaml_model writes it: plain, or double-quoted when needed."""
    if not needs_quote(value):
        return value
    out = io.StringIO()
    _DUMPER.dump({"k": DoubleQuotedScalarString(value)}, out)
    return out.getvalue().rstrip("\n")[len("k: "):]


def parse(raw: str) -> str:
    """The value of a scalar as written in the YAML file."""
    value = _LOADER.load("k: " + raw)["k"]
    return "" if value is None else str(value)


def charge_text(value: str) -> str:
    return repr(float(value))


def same_charge(a: str, b: str) -> bool:
    return (a == "" and b == "") or (a != "" and b != "" and float(a) == float(b))


# ------------------------------------------------------------------------------- YAML
class Model:
    """model/Human-GEM.yml as lines, with the line range of every entry."""

    def __init__(self, path: Path):
        self.path = path
        self.lines = path.read_text(encoding="utf-8").split("\n")
        self.index()

    def index(self) -> None:
        self.entries = {}
        section, start, entry_id = None, None, None
        for i, line in enumerate(self.lines + ["- end:"]):
            if line.startswith("- "):
                if entry_id:
                    self.entries[(section, entry_id)] = (start, i)
                section, start, entry_id = line[2:].split(":")[0], None, None
            elif line == "  - !!omap":
                if entry_id:
                    self.entries[(section, entry_id)] = (start, i)
                start, entry_id = i, None
            elif start is not None and entry_id is None and line.startswith("    - id: "):
                entry_id = parse(line[len("    - id: "):])

    def fields(self, section: str, entry: str) -> dict:
        """field -> (first line, end line) of each field block of an entry."""
        start, end = self.entries[(section, entry)]
        out, name, first = {}, None, None
        for i in range(start + 1, end):
            line = self.lines[i]
            if line.startswith("    - ") and not line.startswith("    - - "):
                if name:
                    out[name] = (first, i)
                name, first = line[len("    - "):].split(":")[0], i
        if name:
            out[name] = (first, end)
        return out

    def get(self, section: str, entry: str, field: str) -> str:
        blocks = self.fields(section, entry)
        if field == "ec-code":
            if "annotation" not in blocks:
                return ""
            first, end = blocks["annotation"]
            codes, inside = [], False
            for line in self.lines[first + 1:end]:
                if line.startswith("      - "):
                    inside = line == "      - ec-code:"
                elif inside and line.startswith("        - "):
                    codes.append(parse(line[len("        - "):]))
            return ";".join(codes)
        if field not in blocks:
            return ""
        line = self.lines[blocks[field][0]]
        return parse(line.split(": ", 1)[1]) if ": " in line else ""

    def set(self, section: str, entry: str, field: str, value: str) -> None:
        n_lines = len(self.lines)
        blocks = self.fields(section, entry)
        if field == "ec-code":
            self._set_ec(blocks, entry, value)
        else:
            text = charge_text(value) if field == "charge" else render(value)
            if field in blocks:
                first, end = blocks[field]
                assert end == first + 1, (entry, field)
                if value == "":
                    del self.lines[first]
                else:
                    self.lines[first] = f"    - {field}: {text}"
            elif value != "":
                self.lines.insert(self._position(blocks, entry, field), f"    - {field}: {text}")
        if len(self.lines) != n_lines:
            self.index()

    def _set_ec(self, blocks: dict, entry: str, value: str) -> None:
        codes = [c for c in value.split(";") if c]
        new = ["      - ec-code:"] + [f"        - {render(c)}" for c in codes]
        if "annotation" in blocks:
            first, end = blocks["annotation"]
            inner = self.lines[first + 1:end]
            keep, inside = [], False
            for line in inner:
                if line.startswith("      - "):
                    inside = line == "      - ec-code:"
                    if inside:
                        continue
                if not inside:
                    keep.append(line)
            body = (new if codes else []) + keep
            self.lines[first:end] = (["    - annotation: !!omap"] + body) if body else []
        elif codes:
            pos = self._position(blocks, entry, "annotation")
            self.lines[pos:pos] = ["    - annotation: !!omap"] + new

    def _position(self, blocks: dict, entry: str, field: str) -> int:
        later = REACTION_ORDER[REACTION_ORDER.index(field) + 1:]
        starts = [blocks[f][0] for f in later if f in blocks]
        if starts:
            return min(starts)
        return max(end for _, end in blocks.values())

    def save(self) -> None:
        self.path.write_text("\n".join(self.lines), encoding="utf-8")


# -------------------------------------------------------------------------------- TSV
class Table:
    """A tab-separated annotation table; a changed row is rewritten with every column."""

    def __init__(self, path: Path):
        self.path = path
        self.lines = path.read_text(encoding="utf-8").split("\n")
        self.header = self.lines[0].split("\t")
        self.rows = {line.split("\t", 1)[0]: i for i, line in enumerate(self.lines[1:], 1) if line}

    def get(self, key: str, column: str) -> str:
        values = self.lines[self.rows[key]].split("\t")
        index = self.header.index(column)
        return values[index] if index < len(values) else ""

    def set(self, key: str, column: str, value: str) -> None:
        values = self.lines[self.rows[key]].split("\t")
        values += [""] * (len(self.header) - len(values))
        values[self.header.index(column)] = value
        self.lines[self.rows[key]] = "\t".join(values)

    def save(self) -> None:
        self.path.write_text("\n".join(self.lines), encoding="utf-8")


# ------------------------------------------------------------------------------- main
def kind_of(entity: str) -> str:
    for prefix in ("MAR", "MAM", "ENSG"):
        if entity.startswith(prefix):
            return prefix
    raise ValueError(f"{entity}: not a reaction, metabolite or gene id")


def read_table(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
    missing = {"entity", "field", "before", "after"} - set(rows[0] if rows else {})
    if missing:
        raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
    return rows


def apply(rows: list[dict], root: Path = Path(".")) -> tuple[Counter, list[str], dict]:
    """Apply the rows to the files under root; returns the counts of changed values per
    field, the rows that do not fit the files, and the edited files (not yet saved)."""
    model = Model(root / YML)
    tables = {k: Table(root / p) for k, p in TSV.items()}
    changed, problems = Counter(), []
    for row in rows:
        entity, field = row["entity"], row["field"]
        before, after = row["before"], row["after"]
        kind = kind_of(entity)
        if field in YAML_FIELDS[kind]:
            section = SECTION[kind]
            if (section, entity) not in model.entries:
                problems.append(f"{entity}: not in {YML}")
                continue
            current = model.get(section, entity, field)
            same = same_charge if field == "charge" else str.__eq__
            if same(current, after):
                continue
            if not same(current, before):
                problems.append(f"{entity} {field}: {current!r}, expected {before!r}")
                continue
            model.set(section, entity, field, after)
        else:
            table = tables[kind]
            if field not in table.header:
                problems.append(f"{entity}: {field} is not a field of {YML} or a column of {table.path}")
                continue
            if entity not in table.rows:
                problems.append(f"{entity}: no row in {table.path}")
                continue
            current = table.get(entity, field)
            if current == after:
                continue
            if current != before:
                problems.append(f"{entity} {field}: {current!r}, expected {before!r}")
                continue
            table.set(entity, field, after)
        changed[field] += 1
    return changed, problems, {"model": model, **tables}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    rows = read_table(args.table)
    changed, problems, files = apply(rows)
    if problems:
        print(f"{len(problems)} row(s) do not fit the model files; nothing written:", file=sys.stderr)
        for line in problems:
            print(f"  {line}", file=sys.stderr)
        return 1
    for field, n in sorted(changed.items()):
        print(f"{field}: {n} value(s) {'to change' if args.dry_run else 'changed'}")
    print(f"{sum(changed.values())} of {len(rows)} row(s) {'would apply' if args.dry_run else 'applied'}")
    if not args.dry_run:
        for f in files.values():
            f.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
