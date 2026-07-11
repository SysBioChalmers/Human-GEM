"""Test YAML import/export round-trip for Human-GEM using the Python RAVEN toolbox.

Python port of code/test/testYamlConversion.m. Verifies that the read/write YAML
conversion does not change the model content.

The MATLAB version compared the in-memory model read from the original file with
a model obtained after one export/import cycle. Here we test the equivalent
property as an idempotency check on the serialisation: the model is read, written
(tmp1), read back, and written again (tmp2). If read/write is lossless, tmp1 and
tmp2 are byte-for-byte identical. Because write_yaml_model preserves the stored
metaData verbatim (it does not regenerate the date), no field needs to be
stripped before comparison.

Usage:
    python code/test/testYamlConversion.py
"""

import sys
import tempfile
from pathlib import Path

from raven_toolbox.io import read_yaml_model, write_yaml_model

# Repository root: this file is code/test/testYamlConversion.py
REPO_ROOT = Path(__file__).resolve().parents[2]
MODEL_FILE = REPO_ROOT / "model" / "Human-GEM.yml"


def main() -> int:
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp1 = Path(tmpdir) / "roundtrip_1.yml"
            tmp2 = Path(tmpdir) / "roundtrip_2.yml"

            # Read original, write it out, read the result and write it again.
            model = read_yaml_model(MODEL_FILE)
            write_yaml_model(model, tmp1)

            reimported = read_yaml_model(tmp1)
            write_yaml_model(reimported, tmp2)

            content1 = tmp1.read_text(encoding="utf-8")
            content2 = tmp2.read_text(encoding="utf-8")

        if content1 != content2:
            print("::error::Re-imported model is different from export")
            return 1
    except Exception as exc:  # noqa: BLE001 - match MATLAB catch-all behaviour
        print(f"::error::There are problems during the conversion import and export: {exc}")
        return 1

    print("The conversion was successful.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
