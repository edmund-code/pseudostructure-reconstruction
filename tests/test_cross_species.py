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


def test_one_to_one_map_resolves_consensus_support():
    table = pd.DataFrame(
        {
            "human_symbol": ["SLC12A1", "SLC12A1", "SLC12A1", "GENE_X", "GENE_NOISE"],
            "mouse_symbol": ["Slc12a1", "Slc12a2", "Slc12a3", "Gene_x", "Gene_noise"],
            "support": [
                "Ensembl,HomoloGene,NCBI,HGNC,OMA",
                "OrthoMCL",
                "OrthoMCL",
                "Ensembl,NCBI,HGNC",
                "OrthoMCL",  # only 1 database, below min_support=3
            ],
        }
    )
    mapping = build_one_to_one_ortholog_map(table, min_support=3)
    assert dict(zip(mapping.human_symbol, mapping.mouse_symbol)) == {
        "GENE_X": "Gene_x",
        "SLC12A1": "Slc12a1",
    }


def test_auto_discover_cluster_identities():
    from pseudospace.markers import auto_discover_cluster_identities

    adata = ad.AnnData(
        np.zeros((6, 5)),
        obs=pd.DataFrame(
            {
                'leiden_coarse': ['0', '0', '1', '1', '2', '2'],
                'comparison_species': ['human', 'mouse', 'human', 'mouse', 'human', 'mouse'],
            }
        ),
        var=pd.DataFrame(index=['Slc5a2', 'Slc12a1', 'Vim', 'Cryab', 'Aqp2']),
    )

    panel_scores = pd.DataFrame(
        [
            {'PT-S1': 1.2, 'PT-S2': 0.9, 'cTAL': -0.2, 'Stroma': 0.0},  # Cluster 0 -> PT
            {'PT-S1': -0.5, 'PT-S2': -0.3, 'cTAL': 2.1, 'Stroma': -0.1},  # Cluster 1 -> cTAL
            {'PT-S1': 0.1, 'PT-S2': -0.1, 'cTAL': -0.4, 'Stroma': 1.5},  # Cluster 2 -> Stroma
        ],
        index=['0', '1', '2'],
    )

    top_de_df = pd.DataFrame(
        {
            'group': ['0', '0', '1', '1', '2', '2'],
            'names': ['Slc5a2', 'Gatm', 'Slc12a1', 'Umod', 'Vim', 'Sparc'],
            'scores': [10.0, 9.0, 15.0, 14.0, 8.0, 7.0],
        }
    )

    labels, summary = auto_discover_cluster_identities(adata, panel_scores, top_de_df)

    assert labels == {'0': 'PT', '1': 'cTAL', '2': 'Stroma'}
    assert summary.loc[summary['Cluster'] == '0', 'Fate'].item() == 'Retained'
    assert summary.loc[summary['Cluster'] == '1', 'Fate'].item() == 'Retained'
    assert summary.loc[summary['Cluster'] == '2', 'Fate'].item() == 'Filtered'
    assert summary.loc[summary['Cluster'] == '0', 'Human Cells'].item() == 1
    assert summary.loc[summary['Cluster'] == '0', 'Mouse Cells'].item() == 1


def test_plot_coarse_cluster_visualizations(tmp_path):
    from pathlib import Path
    from pseudospace.markers import auto_discover_cluster_identities, plot_coarse_cluster_visualizations

    adata = ad.AnnData(
        np.zeros((6, 5)),
        obs=pd.DataFrame(
            {
                'leiden_coarse': ['0', '0', '1', '1', '9', '9'],
                'comparison_species': ['human', 'mouse', 'human', 'mouse', 'human', 'mouse'],
                'celltype_primary_segment': ['PT-S1', 'PT-S2', 'cTAL', 'cTAL', 'PT-S1', 'CCD'],
            }
        ),
        var=pd.DataFrame(index=['Slc5a2', 'Slc12a1', 'Vim', 'Cryab', 'Aqp2']),
    )

    panel_scores = pd.DataFrame(
        [
            {'PT-S1': 1.2, 'PT-S2': 0.9, 'cTAL': -0.2, 'CCD': -0.5, 'Stroma': 0.0},
            {'PT-S1': -0.5, 'PT-S2': -0.3, 'cTAL': 2.1, 'CCD': -0.2, 'Stroma': -0.1},
            {'PT-S1': 0.3, 'PT-S2': -0.1, 'cTAL': -0.1, 'CCD': 0.3, 'Stroma': 0.0},
        ],
        index=['0', '1', '9'],
    )

    top_de_df = pd.DataFrame(
        {
            'group': ['0', '0', '1', '1', '9', '9'],
            'names': ['Slc5a2', 'Gatm', 'Slc12a1', 'Umod', 'Cryab', 'Sgk1'],
            'scores': [10.0, 9.0, 15.0, 14.0, 8.0, 7.0],
            'logfoldchanges': [2.0, 1.5, 3.0, 2.5, 3.2, 1.8],
            'pct_nz_group': [0.9, 0.8, 0.95, 0.85, 0.94, 0.77],
        }
    )

    labels, summary = auto_discover_cluster_identities(adata, panel_scores, top_de_df)
    generated = plot_coarse_cluster_visualizations(adata, panel_scores, top_de_df, summary, tmp_path)

    assert len(generated) >= 4
    for file_path in generated:
        assert Path(file_path).is_file()


