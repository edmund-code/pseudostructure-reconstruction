"""Synthetic checks for the human-to-mouse ortholog boundary."""

import numpy as np
import pandas as pd
import anndata as ad

from pseudospace.cross_species import build_one_to_one_ortholog_map, map_human_to_mouse_space


def test_one_to_one_map_excludes_ambiguous_pairs_and_accepts_override():
    table = pd.DataFrame(
        {
            "human_symbol": ["H1", "H2", "H2", "H3"],
            "mouse_symbol": ["M1", "M2", "M3", "M4"],
        }
    )
    mapping = build_one_to_one_ortholog_map(table, overrides={"H2": "M2"})
    assert dict(zip(mapping.human_symbol, mapping.mouse_symbol)) == {"H1": "M1", "H2": "M2", "H3": "M4"}
    assert mapping.loc[mapping.human_symbol.eq("H2"), "mapping_status"].item() == "curated_override"


def test_human_matrix_is_aggregated_in_shared_mouse_space():
    mapping = pd.DataFrame(
        {
            "human_symbol": ["H1", "H2"],
            "mouse_symbol": ["M1", "M1"],
            "mapping_status": ["curated_override", "curated_override"],
        }
    )
    human = ad.AnnData(
        np.array([[1.0, 2.0, 9.0], [3.0, 4.0, 8.0]]),
        var=pd.DataFrame(index=["H1", "H2", "UNMAPPED"]),
    )
    converted, report = map_human_to_mouse_space(human, mapping)
    assert converted.var_names.tolist() == ["M1"]
    values = converted.X.toarray() if hasattr(converted.X, "toarray") else np.asarray(converted.X)
    np.testing.assert_allclose(values, [[3.0], [7.0]])
    assert report.mapping_status.tolist() == ["mapped", "mapped", "unmapped"]
