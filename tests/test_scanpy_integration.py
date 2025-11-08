from __future__ import annotations

from pathlib import Path

import anndata as ad
import pytest
import sys

TESTS_DIR = Path(__file__).parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from backend.llm.annotator import Annotator
from config.settings import get_settings
from gpt_cell_annotator.scanpy import (
    BatchOptions,
    DiskAnnotationCache,
    GuardrailConfig,
    annotate_anndata,
    annotate_from_markers,
    annotate_rank_genes,
    report_to_dataframe,
    validate_anndata,
)
from gpt_cell_annotator.scanpy import main as scanpy_main
from fixtures.scanpy import marker_dataframe, synthetic_adata


@pytest.fixture()
def marker_db():
    return marker_dataframe()


@pytest.fixture()
def adata():
    return synthetic_adata()


def test_annotate_anndata_adds_columns(adata, marker_db):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    annotator = Annotator(settings=settings, force_mock=True)

    result = annotate_anndata(
        adata.copy(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotator=annotator,
    )

    assert "gptca_label" in result.adata.obs.columns
    assert "gptca_proposed_label" in result.adata.obs.columns
    labels = {label for label in result.adata.obs["gptca_label"] if label}
    assert labels  # at least one label populated
    assert result.report.dataset.summary.total_clusters == 2
    assert len(result.report.dataset.clusters) == 2
    assert result.stats["total_clusters"] == 2
    assert result.report.offline_mode is True


def test_cli_roundtrip(tmp_path: Path, monkeypatch, adata, marker_db):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1

    input_path = tmp_path / "demo.h5ad"
    adata.copy().write(input_path)

    class DummyCache:
        def load(self, *, frame=None, source_path=None, resolver=None, columns):
            return marker_db.copy()

    monkeypatch.setattr("gpt_cell_annotator.scanpy._marker_cache", lambda: DummyCache())

    summary_csv = tmp_path / "summary.csv"
    json_report = tmp_path / "report.json"
    output_path = tmp_path / "annotated.h5ad"

    exit_code = scanpy_main(
        [
            "annotate",
            str(input_path),
            "--cluster-key",
            "cluster",
            "--species",
            "Homo sapiens",
            "--output",
            str(output_path),
            "--summary-csv",
            str(summary_csv),
            "--json-report",
            str(json_report),
            "--offline",
        ]
    )

    assert exit_code == 0
    annotated = ad.read_h5ad(output_path)
    assert "gptca_label" in annotated.obs.columns
    assert "gptca_proposed_label" in annotated.obs.columns
    assert summary_csv.exists()
    assert json_report.exists()


def test_cli_validate_only(tmp_path: Path, monkeypatch, adata, marker_db):
    settings = get_settings()
    settings.validation_min_marker_overlap = 1
    annotated = adata.copy()
    annotated.obs["label"] = ["B cell"] * annotated.n_obs
    input_path = tmp_path / "demo.h5ad"
    annotated.write(input_path)

    class DummyCache:
        def load(self, *, frame=None, source_path=None, resolver=None, columns):
            return marker_db.copy()

    monkeypatch.setattr("gpt_cell_annotator.scanpy._marker_cache", lambda: DummyCache())

    json_report = tmp_path / "validate.json"

    exit_code = scanpy_main(
        [
            "annotate",
            str(input_path),
            "--cluster-key",
            "cluster",
            "--species",
            "Homo sapiens",
            "--label-column",
            "label",
            "--validate-only",
            "--json-report",
            str(json_report),
        ]
    )

    assert exit_code == 0
    assert json_report.exists()


def test_annotate_from_markers(marker_db):
    markers = {"0": ["MS4A1", "CD79A"], "1": ["CD3E", "CD3D"]}
    result = annotate_from_markers(
        markers,
        species="Homo sapiens",
        marker_db=marker_db,
        batch_options=BatchOptions(chunk_size=1),
    )

    assert result.report.dataset.summary.total_clusters == 2
    labels = {cluster.annotation["primary_label"] for cluster in result.report.dataset.clusters}
    assert any(label != "Unknown or Novel" for label in labels)


def test_annotate_rank_genes_wrapper(marker_db):
    rankings = {
        "0": ["MS4A1", "CD79A", "CD74"],
        "1": ["CD3E", "CD3D", "CD2"],
    }
    result = annotate_rank_genes(
        rankings,
        species="Homo sapiens",
        marker_db=marker_db,
    )
    assert result.report.dataset.summary.total_clusters == 2


def test_validate_anndata_guardrail_override(adata, marker_db):
    sample = adata.copy()
    sample.obs["gpt_label"] = ["B cell"] * sample.n_obs

    report_default = validate_anndata(
        sample.copy(),
        "cluster",
        species="Homo sapiens",
        label_column="gpt_label",
        marker_db=marker_db,
    )
    assert report_default.dataset.summary.flagged_clusters <= 2

    report_strict = validate_anndata(
        sample.copy(),
        "cluster",
        species="Homo sapiens",
        label_column="gpt_label",
        marker_db=marker_db,
        guardrails=GuardrailConfig(min_marker_overlap=5),
    )
    assert report_strict.dataset.summary.flagged_clusters >= report_default.dataset.summary.flagged_clusters


def test_disk_cache_reduces_batches(monkeypatch, adata, marker_db, tmp_path: Path):
    call_counter = {"count": 0}
    cache = DiskAnnotationCache(tmp_path / "cache")

    def fake_annotate(self, payload, context):
        call_counter["count"] += 1
        return {str(item["cluster_id"]): {"primary_label": "B cell"} for item in payload}

    monkeypatch.setattr("backend.llm.annotator.Annotator.annotate_batch", fake_annotate)

    result_first = annotate_anndata(
        adata.copy(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotation_cache=cache,
        batch_options=BatchOptions(chunk_size=1),
    )
    assert result_first.report.llm_batches == 2

    call_counter["count"] = 0
    result_second = annotate_anndata(
        adata.copy(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotation_cache=cache,
        batch_options=BatchOptions(chunk_size=1),
    )
    assert call_counter["count"] == 0
    assert result_second.report.cache_hits == 2


def test_report_to_dataframe_handles_scanpy_dataset(marker_db, adata):
    annotator = Annotator(force_mock=True)
    result = annotate_anndata(
        adata.copy(),
        "cluster",
        species="Homo sapiens",
        marker_db=marker_db,
        annotator=annotator,
    )
    frame = report_to_dataframe(result.report)
    assert set(frame.columns) >= {"cluster_id", "primary_label", "warnings"}
