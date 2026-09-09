"""Make the QC scripts importable from the tests.

The scripts in code/test are run as files by the workflow, not installed as a
package, so the directory holding them is put on sys.path here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
