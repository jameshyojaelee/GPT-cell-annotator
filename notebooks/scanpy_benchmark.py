"""Performance smoke benchmark for the Scanpy integration.

Generates a synthetic 50k-cell AnnData object, runs `annotate_anndata`
with and without the disk cache, and prints runtime plus peak RSS.
"""

from __future__ import annotations

import json
import resource
import tempfile
import time
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.llm.annotator import Annotator
from gpt_cell_annotator.scanpy import (
    BatchOptions,
    DiskAnnotationCache,
    GuardrailConfig,
    ScanpyAnnotationResult,
    annotate_anndata,
)


def _build_marker_db() -> pd.DataFrame:
    records = [
        {
            "source": "Benchmark",
            "cell_type": "B cell",
            "ontology_id": "CL:0000236",
            "gene_symbol": "MS4A1",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "benchmark",
            "reference": "",
            "evidence_score": "high",
        },
        {
            "source": "Benchmark",
            "cell_type": "CD4 T cell",
            "ontology_id": "CL:0000624",
            "gene_symbol": "CD3E",
            "species": "Homo sapiens",
            "tissue": "Blood",
            "evidence": "benchmark",
            "reference": "",
            "evidence_score": "high",
        },
    ]
    return pd.DataFrame.from_records(records)


def _synthetic_adata(*, cells: int = 50_000, genes: int = 200, clusters: int = 40) -> ad.AnnData:
    rng = np.random.default_rng(seed=42)
    gene_names = [f"GENE{i:03d}" for i in range(genes)]
    cluster_ids = np.arange(clusters)
    obs_clusters = rng.choice(cluster_ids, size=cells)
    X = rng.poisson(lam=1.0, size=(cells, genes)).astype(np.float32)

    obs = pd.DataFrame({"cluster": obs_clusters.astype(str)}, index=[f"cell_{i}" for i in range(cells)])
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X, obs=obs, var=var)

    # Populate rank_genes_groups with deterministic top genes per cluster to skip expensive computations.
    names = {}
    for cluster in cluster_ids:
        offset = (cluster * 5) % genes
        top_genes = [gene_names[(offset + i) % genes] for i in range(10)]
        names[str(cluster)] = top_genes
    adata.uns["rank_genes_groups"] = {"names": names}
    return adata


def _rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    # ru_maxrss is kilobytes on Linux.
    return usage.ru_maxrss / 1024.0


def benchmark_run(
    *,
    adata: ad.AnnData,
    marker_db: pd.DataFrame,
    cache: DiskAnnotationCache | None,
) -> tuple[ScanpyAnnotationResult, float, float]:
    annotator = Annotator(force_mock=True)
    start = time.perf_counter()
    result = annotate_anndata(
        adata,
        "cluster",
        species="Homo sapiens",
        top_n_markers=10,
        marker_db=marker_db,
        annotator=annotator,
        compute_rankings=False,
        batch_options=BatchOptions(chunk_size=64),
        guardrails=GuardrailConfig(min_marker_overlap=1),
        annotation_cache=cache,
    )
    duration = time.perf_counter() - start
    rss = _rss_mb()
    return result, duration, rss


def main(output: Path | None = None) -> None:
    marker_db = _build_marker_db()
    adata = _synthetic_adata()

    rows: list[dict[str, float | int | str]] = []

    result_cold, duration_cold, rss_cold = benchmark_run(
        adata=adata.copy(),
        marker_db=marker_db,
        cache=None,
    )
    rows.append(
        {
            "mode": "cold",
            "duration_sec": round(duration_cold, 3),
            "peak_rss_mb": round(rss_cold, 1),
            "cache_hits": result_cold.report.cache_hits,
            "llm_batches": result_cold.report.llm_batches,
        }
    )

    with tempfile.TemporaryDirectory(prefix="gca-scanpy-benchmark-") as tmpdir:
        cache = DiskAnnotationCache(Path(tmpdir))
        # Warm run to populate cache
        benchmark_run(
            adata=adata.copy(),
            marker_db=marker_db,
            cache=cache,
        )
        result_cached, duration_cached, rss_cached = benchmark_run(
            adata=adata.copy(),
            marker_db=marker_db,
            cache=cache,
        )
        rows.append(
            {
                "mode": "cached",
                "duration_sec": round(duration_cached, 3),
                "peak_rss_mb": round(rss_cached, 1),
                "cache_hits": result_cached.report.cache_hits,
                "llm_batches": result_cached.report.llm_batches,
            }
        )

    print("Scanpy Benchmark (50k cells)")
    for row in rows:
        print(
            f" - {row['mode']:>6}: {row['duration_sec']:>6.3f}s | "
            f"RSS {row['peak_rss_mb']:>6.1f} MB | "
            f"cache_hits={row['cache_hits']} | llm_batches={row['llm_batches']}"
        )

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        print(f"\nWrote benchmark telemetry to {output}")


if __name__ == "__main__":
    main()
