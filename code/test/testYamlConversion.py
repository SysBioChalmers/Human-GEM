"""Model YAML round-trip check for Human-GEM, with and without the RAVEN toolbox.

Verifies that reading and writing model/Human-GEM.yml is lossless (the committed
file is canonical) for the requested toolchain, so users of either cobrapy or the
RAVEN toolbox can rely on it. The model is read, written (tmp1), read back and
written again (tmp2); if read/write is lossless, tmp1 and tmp2 are byte-for-byte
identical (write preserves the stored metaData verbatim, so nothing needs to be
stripped before comparison).

Usage:
    python code/test/testYamlConversion.py [--tool cobra|raven-toolbox]

With no --tool it runs the cobra round-trip and, if raven-toolbox is importable,
the raven-toolbox round-trip as well.
"""

import argparse
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"


def _roundtrip(read, write) -> bool:
    """Read -> write -> read -> write; return True if the two written files match."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp1 = Path(tmpdir) / "roundtrip_1.yml"
        tmp2 = Path(tmpdir) / "roundtrip_2.yml"
        write(read(MODEL_FILE), tmp1)
        write(read(tmp1), tmp2)
        return tmp1.read_text(encoding="utf-8") == tmp2.read_text(encoding="utf-8")


def _cobra_roundtrip() -> bool:
    import cobra  # cobrapy only; no RAVEN toolbox needed
    return _roundtrip(
        lambda p: cobra.io.load_yaml_model(str(p)),
        lambda m, p: cobra.io.save_yaml_model(m, str(p)),
    )


def _raven_roundtrip() -> bool:
    from raven_toolbox.io import read_yaml_model, write_yaml_model
    return _roundtrip(read_yaml_model, write_yaml_model)


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
        try:
            lossless = TOOLS[tool]()
        except Exception as exc:  # noqa: BLE001 - report any conversion problem
            print(f"::error::{tool} round-trip raised an error: {exc}")
            failed = True
            continue
        if lossless:
            print(f"{tool} round-trip: the conversion was lossless.")
        else:
            print(f"::error::{tool} round-trip: the re-exported model differs from the export.")
            failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
