"""Load the Human-GEM model in Python as an annotated cobra model."""
from __future__ import annotations

from pathlib import Path

import cobra

from .annotation import annotate_gem

# This package sits at the repository root, next to the model/ directory.
_DEFAULT_MODEL_DIR = Path(__file__).resolve().parents[1] / "model"


def load_model(
    model_dir: str | Path | None = None,
    *,
    annotate: bool = True,
) -> cobra.Model:
    """Load Human-GEM as an annotated cobrapy model.

    Parameters
    ----------
    model_dir
        Directory holding ``Human-GEM.yml`` and the annotation tables
        (``reactions.tsv`` / ``metabolites.tsv`` / ``genes.tsv``). Defaults to
        the ``model/`` directory of the Human-GEM repository this package lives
        in.
    annotate
        If ``True`` (default), merge the cross-references and SBO terms from the
        annotation tables onto the model. If ``False``, return the bare model as
        read from the YAML (formula / charge / eccodes only).

    Returns
    -------
    cobra.Model
        The Human-GEM model, annotated unless ``annotate=False``.
    """
    from raven_toolbox.io import read_yaml_model

    model_dir = Path(model_dir) if model_dir is not None else _DEFAULT_MODEL_DIR
    yml = model_dir / "Human-GEM.yml"
    if not yml.is_file():
        raise FileNotFoundError(
            f"Human-GEM.yml not found in {model_dir}. Pass model_dir=... pointing at "
            "a Human-GEM model/ directory (the packaged distribution bundles it)."
        )

    model = read_yaml_model(str(yml))
    if annotate:
        annotate_gem(model, model_dir)
    return model
