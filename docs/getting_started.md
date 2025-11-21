# Getting Started

This guide gets you from install to your first annotations with the CLI, Scanpy helpers, and the R wrapper.

## 1) Environment setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install "gpt-cell-annotator[scanpy]"  # includes CLI + Scanpy extras
```

If you prefer Poetry, run `poetry install --extras "scanpy" --extras "dev"`.

## 2) Offline demo (CLI)

```bash
# Build the bundled marker knowledge base
gca build-db --offline

# Annotate demo PBMC markers and save the JSON report
gca annotate data/demo/pbmc_markers.csv \
  --species "Homo sapiens" \
  --offline \
  --out-json demo_annotations.json
```

The JSON report matches the schema used by the Scanpy helpers and the R wrapper. Add `--out-csv` to emit a flat summary.

## 3) Scanpy workflow

```bash
gca scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --species "Homo sapiens" \
  --cache-dir ~/.cache/gca/annotations \
  --json-report reports/pbmc_report.json
```

Use `--validate-only` to run guardrails on existing labels without calling the LLM, and `--offline` to force the heuristic annotator.

The repo bundles `data/demo/pbmc_demo.h5ad` with Leiden clustering precomputed. To rebuild it yourself, install `python-igraph` and `leidenalg` and run the snippet in the README under “Bundled demos.”

## 4) Programmatic usage (Python)

```python
from gpt_cell_annotator.scanpy import annotate_anndata

result = annotate_anndata(
    "data/demo/pbmc_demo.h5ad",
    cluster_key="leiden",
    species="Homo sapiens",
    cache_dir="~/.cache/gca/annotations",
)

result.adata.obs[["gptca_label", "gptca_status"]].head()
print(result.report_dict())
```

The helper reuses cached marker databases and returns a structured report alongside the updated AnnData object.

## 5) Seurat / R usage

```r
library(gptcellannotator)
cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
gptca_config_set(cfg)

markers <- read.csv(system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator"))
annotations <- gptca_annotate_markers(markers, species = "Homo sapiens")
annotations$clusters
```

Set `offline = FALSE` with a `base_url` and `api_key` only if you need to target an external REST service.

## 6) Validation guardrails

Key knobs (env vars or `config.settings.Settings` fields):

- `VALIDATION_MIN_MARKER_OVERLAP` (default 2) — minimum shared markers required to treat a label as supported.
- `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` — bump confidence tiers based on overlap.
- `RAG_ENABLED` / `RAG_TOP_K` — control retrieval grounding of model responses.

Tune these before running the CLI or Scanpy helpers to adjust conservatism. Offline mode still applies guardrails against the bundled marker DB.

## 7) Troubleshooting

- Missing `marker_db.parquet` → rerun `gca build-db` or set `GCA_MARKER_DB_PATH`.
- `scanpy` import errors → install the scanpy extra.
- Empty annotations in offline mode → ensure markers column parses correctly (`markers` column in CSV with comma/JSON lists).

See `docs/faq.md` for more scenarios.
