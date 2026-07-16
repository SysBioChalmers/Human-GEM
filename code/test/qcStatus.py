"""Read and write the combined QC status file (data/testResults/qc_status.tsv).

Several fast checks each produce a single scalar - a pass/fail, a failed/total
count, or the growth value. Rather than commit a separate one-line file per check,
they all upsert into one key/value TSV:

    check	result
    growth	123.4
    roundtrip_cobra	pass
    roundtrip_raven	pass
    tasks_essential	0/57
    tasks_verification	0/21
    yamllint	pass

Upsert (read, set the one key, rewrite) keeps it order-independent and rerun-safe,
and because the QC steps run sequentially in one job there is no contention. Keys
are a fixed set, so nothing stale accumulates.

CLI (used by the workflow's shell steps):
    python code/test/qcStatus.py <key> <value>     # set one key
    python code/test/qcStatus.py --get <key>       # print one value (empty if unset)
"""

import sys
from pathlib import Path

STATUS_FILE = Path(__file__).resolve().parents[2] / "data" / "testResults" / "qc_status.tsv"
_HEADER = ("check", "result")


def read_status(path: Path = STATUS_FILE) -> dict:
    """Return the status file as a {check: result} dict ({} if it does not exist)."""
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] == _HEADER[0]:  # header row
            continue
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def set_status(key: str, value: str) -> None:
    """Upsert one key and rewrite the file (header + keys sorted for a stable diff)."""
    data = read_status()
    data[key] = str(value)
    STATUS_FILE.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(_HEADER)] + [f"{k}\t{data[k]}" for k in sorted(data)]
    STATUS_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")


def get_status(key: str) -> str:
    return read_status().get(key, "")


def main(argv: list[str]) -> int:
    if len(argv) == 2 and argv[0] == "--get":
        print(get_status(argv[1]))
        return 0
    if len(argv) == 2:
        set_status(argv[0], argv[1])
        return 0
    print("usage: qcStatus.py <key> <value> | --get <key>", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
