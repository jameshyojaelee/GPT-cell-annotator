# Scanpy Integration Guide

Use GPT Cell Annotator directly in Scanpy pipelines via the Python API or the `gca scanpy` CLI.

## Setup

```bash
pip install "gpt-cell-annotator[scanpy]"
gca build-db --offline  # materialise bundled marker assets
```

## Annotate in a notebook

```python
import scanpy as sc
from gpt_cell_annotator import BatchOptions, GuardrailConfig, annotate_anndata

adata = sc.read_h5ad("data/demo/pbmc_demo.h5ad")

result = annotate_anndata(
    adata,
    cluster_key="leiden",
    species="Homo sapiens",
    tissue="Peripheral blood",
    top_n_markers=5,
    batch_options=BatchOptions(chunk_size=16),
    guardrails=GuardrailConfig(min_marker_overlap=2),
)

result.adata.obs[["gptca_label", "gptca_status"]].head()
result.report.summary
```

The helper will compute `rank_genes_groups` when missing, normalise markers (including ortholog mapping), call the annotator (mock or live), apply validation guardrails, and return a `ScanpyAnnotationResult` with:

- `adata`: updated AnnData object with `gptca_*` columns
- `annotations`: list of per-cluster annotations
- `report`: `ScanpyDatasetReport` containing guardrail stats and warnings

## CLI workflow

```bash
gca scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --species "Homo sapiens" \
  --chunk-size 24 \
  --cache-dir ~/.cache/gca/annotations \
  --json-report reports/pbmc_report.json \
  --offline

# Guardrail-only validation of existing labels
gca scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --label-column curated_label \
  --species "Homo sapiens" \
  --validate-only
```

Key flags:
- `--marker-db` to point at a custom Parquet DB (defaults to cached assets or `GCA_MARKER_DB_PATH`)
- `--preset` for quick species/tissue presets (`human_pbmc`, `mouse_pbmc`, `mouse_brain`)
- `--cache-dir` / `GCA_CACHE_DIR` to persist annotation responses between runs
- `--guardrail-min-overlap`, `--guardrail-force-unknown` to override thresholds without editing env vars
- `--offline` to force the heuristic annotator (no network calls)

## Guardrails and retrieval

Environment variables (or `config.settings.Settings`) control validation and retrieval:

- `VALIDATION_MIN_MARKER_OVERLAP` (default 2) downgrades labels with low support.
- `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` adjust confidence tiers.
- `RAG_ENABLED`, `RAG_TOP_K`, `RAG_MIN_OVERLAP` gate retrieval augmentation.
- `SYNONYM_ENABLE_ORTHOLOGS` toggles ortholog mapping; `ORTHOLOG_MAPPING_PATH` overrides the mapping file.

## Troubleshooting

- Missing `marker_db.parquet` → rerun `gca build-db` or set `--marker-db`.
- `scanpy` import errors → ensure the scanpy extra is installed in the active environment.
- Empty markers in offline mode → confirm the CSV `markers` column parses correctly (comma/JSON lists).
- Unexpected mock mode → set `OPENAI_API_KEY` before launching notebooks/CLI, or pass `force_mock=True` explicitly for offline runs.
