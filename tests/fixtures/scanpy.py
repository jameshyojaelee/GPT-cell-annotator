from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from gpt_cell_annotator.scanpy import MARKER_DB_COLUMNS


def marker_dataframe() -> pd.DataFrame:
    records = [
        {
            "source": "Demo",
            "cell_type": "B cell",
            "ontology_id": "CL:0000236",
            "gene_symbol": "MS4A1",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "Demo",
            "reference": "",
            "evidence_score": "high",
        },
        {
            "source": "Demo",
            "cell_type": "CD4 T cell",
            "ontology_id": "CL:0000624",
            "gene_symbol": "CD3E",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "Demo",
            "reference": "",
            "evidence_score": "high",
        },
    ]
    return pd.DataFrame.from_records(records, columns=MARKER_DB_COLUMNS)


def synthetic_adata() -> ad.AnnData:
    matrix = np.zeros((6, 4), dtype=float)
    matrix[:3, 0] = 5.0  # cluster 0 expresses gene0 (MS4A1)
    matrix[3:, 1] = 5.0  # cluster 1 expresses gene1 (CD3E)
    obs = pd.DataFrame(
        {"cluster": ["0", "0", "0", "1", "1", "1"]},
        index=[f"cell_{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["MS4A1", "CD3E", "GNLY", "LYZ"])
    adata = ad.AnnData(matrix, obs=obs, var=var)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata
