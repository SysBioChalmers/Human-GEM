"""Give the deprecated-identifier lists the columns of the model annotation tables, and fill
values from a table.

    python code/modelCuration/backfillDeprecatedColumns.py --table data/modelCuration/<name>.tsv [--dry-run]

data/deprecatedIdentifiers/deprecatedReactions.tsv and deprecatedMetabolites.tsv gain every
column of model/reactions.tsv and model/metabolites.tsv that they lack, at the same position,
empty in every row. The table (tab-separated: list, id, column, value; list is reactions or
metabolites) then fills cells of these lists: an empty cell takes the value, a cell that
already holds it is left as it is, and a cell with another value is reported and nothing is
written. Rows of the table for ids that are not in the list are skipped. --dry-run reports
what would change and writes nothing.
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

LISTS = {
    "reactions": (Path("data/deprecatedIdentifiers/deprecatedReactions.tsv"), Path("model/reactions.tsv")),
    "metabolites": (Path("data/deprecatedIdentifiers/deprecatedMetabolites.tsv"), Path("model/metabolites.tsv")),
}


class Table:
    def __init__(self, path: Path):
        self.path = path
        lines = path.read_text(encoding="utf-8").split("\n")
        assert lines[-1] == "", f"{path} does not end with a newline"
        self.header = lines[0].split("\t")
        self.rows = [line.split("\t") for line in lines[1:-1]]
        for values in self.rows:
            values += [""] * (len(self.header) - len(values))
        self.index = {values[0]: values for values in self.rows}

    def add_columns(self, header: list[str]) -> list[str]:
        """Insert the columns of header that this table lacks, at their position in header."""
        missing = [c for c in header if c not in self.header]
        unknown = [c for c in self.header if c not in header]
        if unknown:
            raise ValueError(f"{self.path}: columns {unknown} are not in the model table")
        for column in missing:
            pos = header.index(column)
            self.header.insert(pos, column)
            for values in self.rows:
                values.insert(pos, "")
        if self.header != header:
            raise ValueError(f"{self.path}: column order differs from the model table")
        return missing

    def save(self) -> None:
        lines = ["\t".join(self.header)] + ["\t".join(values) for values in self.rows]
        self.path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def backfill(rows: list[dict], root: Path = Path(".")) -> tuple[dict, Counter, list[str], dict]:
    lists, added = {}, {}
    for name, (path, model_path) in LISTS.items():
        lists[name] = Table(root / path)
        model_header = (root / model_path).read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
        added[name] = lists[name].add_columns(model_header)
    filled, problems = Counter(), []
    for row in rows:
        table = lists[row["list"]]
        values = table.index.get(row["id"])
        if values is None:
            continue
        if row["column"] not in table.header:
            problems.append(f"{row['id']}: no column {row['column']} in {table.path}")
            continue
        pos = table.header.index(row["column"])
        if values[pos] == row["value"]:
            continue
        if values[pos]:
            problems.append(f"{row['id']} {row['column']}: {values[pos]!r}, table gives {row['value']!r}")
            continue
        values[pos] = row["value"]
        filled[(row["list"], row["column"])] += 1
    return added, filled, problems, lists


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    with args.table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
    added, filled, problems, lists = backfill(rows)
    if problems:
        print(f"{len(problems)} row(s) do not fit the deprecated lists; nothing written:", file=sys.stderr)
        for line in problems:
            print(f"  {line}", file=sys.stderr)
        return 1
    for name, columns in added.items():
        print(f"{LISTS[name][0]}: {'would gain' if args.dry_run else 'gained'} {columns or 'no columns'}")
    for (name, column), n in sorted(filled.items()):
        print(f"{name} {column}: {n} value(s) {'to fill' if args.dry_run else 'filled'}")
    if not args.dry_run:
        for table in lists.values():
            table.save()
    return 0


if __name__ == "__main__":
    sys.exit(main())
