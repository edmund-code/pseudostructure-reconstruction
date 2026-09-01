"""pseudospace: reusable helpers extracted from the kidney pseudospace workflows.

Repo-local package (no installation required) imported directly because the
notebooks run from the repository root. Import submodules explicitly, e.g.::

    from pseudospace.markers import annotate_clusters_by_de
    from pseudospace.trajectory import recompute_subset_dpt

See ``analysis/mouse_only_pseudospace.py`` and ``analysis/notebooks/`` for the
current workflows that drive these functions.
"""

__version__ = "0.1.0"
