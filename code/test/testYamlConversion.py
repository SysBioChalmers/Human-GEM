"""Model YAML round-trip check for Human-GEM, with and without the RAVEN toolbox.

Verifies that reading and writing model/Human-GEM.yml keeps the model, for the
requested toolchain, so users of either cobrapy or the RAVEN toolbox can rely on it.
The model is read and written (tmp1), then tmp1 is read and written again (tmp2). Two
comparisons follow:

  * first read: tmp1 against the committed file, parsed, entry by entry and field by
    field (number formatting and key order aside). This sees a field the reader drops
    or changes when it reads the committed file;
  * re-read: tmp1 and tmp2 must be byte-for-byte identical, so that a written file
    reads back to itself.

raven-toolbox must keep every field. cobrapy does not read the Human-GEM extensions:
it drops the ``metaData`` block, ``metFrom``, ``rxnFrom``, ``references`` and
``confidence_score`` (:data:`COBRA_DROPS`), and writes the model id, name and its own
format version as top-level entries (:data:`COBRA_ADDS`). Those are listed as a note;
any other loss or change fails, as does a field cobrapy starts to keep. A field that
is absent on one side and empty on the other (an empty GPR, say) counts as the same.

Usage:
    python code/test/testYamlConversion.py [--tool cobra|raven-toolbox]

With no --tool it runs the cobra round-trip and, if raven-toolbox is importable,
the raven-toolbox round-trip as well.
"""

import argparse
import sys
import tempfile
from collections import Counter
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"

# Fields cobrapy does not read from a Human-GEM YAML: section -> fields ("*" = the
# whole section).
COBRA_DROPS = {
    "metaData": {"*"},
    "metabolites": {"metFrom"},
    "reactions": {"rxnFrom", "references", "confidence_score"},
}

# Top-level entries cobrapy writes that the committed file keeps in metaData.
COBRA_ADDS = {"id", "name", "version"}

# Differences listed per comparison before the rest are only counted.
MAX_LISTED = 20


def _plain(value):
    """Turn the PyYAML form of ``!!omap`` (a list of 1-item pairs) into dicts, recursively."""
    if isinstance(value, list):
        if value and all(isinstance(item, tuple) and len(item) == 2 for item in value):
            return {key: _plain(item) for key, item in value}
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    return value


def load_sections(path: Path) -> dict:
    """``{section: {entry id: {field: value}}}``; ``metaData`` and ``compartments``
    become a single entry each, keyed by the section name."""
    with open(path, encoding="utf-8") as fh:
        top = _plain(yaml.load(fh, Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader)))
    sections: dict = {}
    for name, content in top.items():
        if isinstance(content, list):
            sections[name] = {entry.get("id"): entry for entry in content}
        elif isinstance(content, dict):
            sections[name] = {name: content}
        else:  # a top-level scalar, such as cobrapy's model id
            sections[name] = {name: {"value": content}}
    return sections


def _same(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b


def _empty(value) -> bool:
    return value is None or value == "" or value == [] or value == {}


def compare_first_read(
    committed: dict, written: dict, allowed: dict, added: set = frozenset(),
) -> tuple[list[str], Counter]:
    """Differences between the committed and the first-written model.

    Returns ``(unexpected, expected)``: each unexpected difference as a line, and a count
    per ``section.field`` of the losses ``allowed`` lists and the top-level entries
    ``added`` lists.
    """
    unexpected: list[str] = []
    expected: Counter = Counter()
    known_added = (set(written) - set(committed)) & set(added)
    for section in sorted(known_added):
        expected[f"{section} (top-level entry written)"] += 1
    for section in sorted((set(committed) | set(written)) - known_added):
        dropped_fields = allowed.get(section, set())
        old, new = committed.get(section, {}), written.get(section, {})
        if "*" in dropped_fields:
            if section in written:
                unexpected.append(f"{section}: expected to be dropped, but written")
            elif old:
                expected[f"{section} (whole section)"] += 1
            continue
        for entry_id in sorted(set(old) - set(new), key=str):
            unexpected.append(f"{section} {entry_id}: missing after the first read")
        for entry_id in sorted(set(new) - set(old), key=str):
            unexpected.append(f"{section} {entry_id}: added by the first read")
        for entry_id in sorted(set(old) & set(new), key=str):
            a, b = old[entry_id], new[entry_id]
            for field in sorted(set(a) | set(b)):
                if field in dropped_fields:
                    if field in b:
                        unexpected.append(f"{section} {entry_id}.{field}: expected to be dropped, but written")
                    elif a.get(field) not in (None, ""):
                        expected[f"{section}.{field}"] += 1
                elif field not in b:
                    if not _empty(a[field]):
                        unexpected.append(f"{section} {entry_id}.{field}: dropped ({a[field]!r})")
                elif field not in a:
                    if not _empty(b[field]):
                        unexpected.append(f"{section} {entry_id}.{field}: added ({b[field]!r})")
                elif not _same(a[field], b[field]):
                    unexpected.append(f"{section} {entry_id}.{field}: {a[field]!r} -> {b[field]!r}")
    return unexpected, expected


def _roundtrip(read, write, allowed: dict, added: set = frozenset()) -> bool:
    """Read -> write (tmp1) -> read -> write (tmp2); report both comparisons."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp1 = Path(tmpdir) / "roundtrip_1.yml"
        tmp2 = Path(tmpdir) / "roundtrip_2.yml"
        write(read(MODEL_FILE), tmp1)
        write(read(tmp1), tmp2)
        stable = tmp1.read_text(encoding="utf-8") == tmp2.read_text(encoding="utf-8")
        unexpected, expected = compare_first_read(load_sections(MODEL_FILE), load_sections(tmp1), allowed, added)

    for key, count in sorted(expected.items()):
        print(f"  not read (known): {key} in {count} entries")
    if unexpected:
        print(f"::error::first read differs from the committed file in {len(unexpected)} place(s):")
        for line in unexpected[:MAX_LISTED]:
            print(f"  {line}")
        if len(unexpected) > MAX_LISTED:
            print(f"  ... and {len(unexpected) - MAX_LISTED} more")
    if not stable:
        print("::error::the re-exported model differs from the first export.")
    return stable and not unexpected


def _cobra_roundtrip() -> bool:
    import cobra  # cobrapy only; no RAVEN toolbox needed
    return _roundtrip(
        lambda p: cobra.io.load_yaml_model(str(p)),
        lambda m, p: cobra.io.save_yaml_model(m, str(p)),
        COBRA_DROPS,
        COBRA_ADDS,
    )


def _raven_roundtrip() -> bool:
    from raven_toolbox.io import read_yaml_model, write_yaml_model
    return _roundtrip(read_yaml_model, write_yaml_model, {})


TOOLS = {"cobra": _cobra_roundtrip, "raven-toolbox": _raven_roundtrip}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", choices=list(TOOLS),
                        help="round-trip with this toolchain only (default: cobra, "
                             "plus raven-toolbox if it is installed)")
    args = parser.parse_args()

    if args.tool:
        tools = [args.tool]
    else:
        tools = ["cobra"]
        try:
            import raven_toolbox  # noqa: F401
            tools.append("raven-toolbox")
        except ImportError:
            print("raven-toolbox not installed; running the cobra round-trip only.")

    failed = False
    for tool in tools:
        print(f"{tool} round-trip:")
        try:
            lossless = TOOLS[tool]()
        except Exception as exc:  # noqa: BLE001 - report any conversion problem
            print(f"::error::{tool} round-trip raised an error: {exc}")
            failed = True
            continue
        if lossless:
            print(f"{tool} round-trip: the conversion kept the model.")
        else:
            print(f"::error::{tool} round-trip: the model changed; see above.")
            failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
