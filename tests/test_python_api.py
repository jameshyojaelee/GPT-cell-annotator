from __future__ import annotations

import pandas as pd

from backend.llm.annotator import Annotator
from gpt_cell_annotator import annotate_markers


def _marker_db() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source": "Demo",
                "cell_type": "B cell",
                "ontology_id": "CL:0000236",
                "gene_symbol": "MS4A1",
                "species": "Homo sapiens",
                "tissue": "Blood",
                "evidence": "",
                "reference": "",
                "evidence_score": "high",
            }
        ]
    )


def test_annotate_markers_returns_report() -> None:
    annotator = Annotator(force_mock=True)
    result = annotate_markers(
        {"0": ["MS4A1", "CD79A"]},
        species="Homo sapiens",
        marker_db=_marker_db(),
        annotator=annotator,
        min_marker_overlap=1,
    )

    assert result.report.summary.total_clusters == 1
    assert result.report.clusters[0].annotation["primary_label"]
