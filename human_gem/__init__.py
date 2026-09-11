"""Python interface to the Human-GEM genome-scale metabolic model.

Loads Human-GEM as a cobrapy model with its full cross-reference and SBO
annotations, on top of `raven-toolbox`. The heavy machinery (YAML I/O, ftINIT,
metabolic tasks, gap-filling) lives in raven-toolbox; this package adds the
Human-GEM-specific glue - most importantly loading the model with the
annotation tables (`reactions.tsv` / `metabolites.tsv` / `genes.tsv`) merged
onto it, which a bare `cobra.io.load_yaml_model` does not do.

    import human_gem
    model = human_gem.load_model()          # annotated cobra.Model
"""
from __future__ import annotations

from .io import load_model

__version__ = "0.1.0"
__all__ = ["load_model"]
