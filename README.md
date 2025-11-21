# GPT-Cell-Annotator

GPT-Cell-Annotator is a CLI + Scanpy toolkit with an R wrapper for labeling single-cell RNA-seq clusters. It pairs a curated marker knowledge base with prompt-engineered LLM calls, optional retrieval, and validation guardrails. Everything ships in the CLI with offline demo assets; Scanpy helpers and the R wrapper reuse the same engine. No backend service is required—offline mode is the default.

## What it does
- Annotate marker CSVs from the CLI (`gca annotate`) with offline mock mode by default (live OpenAI calls when `OPENAI_API_KEY` is set).
- Build and ship the marker knowledge base locally (`gca build-db`) with checksum verification.
- Run end-to-end Scanpy workflows (`gca scanpy annotate`) with caching, presets, and guardrails.
- Normalise markers across species/orthologs; validate LLM suggestions against curated databases.
- Bridge to Seurat via the R package (`gptcellannotator`) which shells out to the CLI by default.

## Quick start

```bash
# Install Python package + Scanpy helpers
pip install "gpt-cell-annotator[scanpy]"

# Offline demo (uses bundled assets + heuristic mock annotator)
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json

# Rebuild marker DB artifacts locally
gca build-db --offline --output-dir ~/.cache/gca/db

# Annotate AnnData directly
gca scanpy annotate data/demo/pbmc_demo.h5ad \
  --cluster-key leiden \
  --species "Homo sapiens" \
  --chunk-size 16 \
  --cache-dir ~/.cache/gca/annotations \
  --json-report reports/pbmc_report.json
```

Offline mode disables network calls and leans on the heuristic annotator plus bundled assets. Live LLM calls require `OPENAI_API_KEY`.

### Bundled demos
- Marker CSV: `data/demo/pbmc_markers.csv` (used in `gca annotate` examples).
- AnnData: `data/demo/pbmc_demo.h5ad` with Leiden clustering precomputed for `gca scanpy annotate`.
  - To regenerate: install `python-igraph` and `leidenalg`, then run:
    ```bash
    python - <<'PY'
    import scanpy as sc
    from pathlib import Path
    adata = sc.datasets.pbmc3k()
    sc.pp.recipe_zheng17(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, key_added="leiden")
    out = Path("data/demo/pbmc_demo.h5ad")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out, compression="gzip")
    PY
    ```

## Scanpy + guardrails in Python

```python
from gpt_cell_annotator.scanpy import annotate_anndata
from config.settings import get_settings

settings = get_settings()
settings.validation_min_marker_overlap = 2

result = annotate_anndata(
    "data/demo/pbmc_demo.h5ad",
    cluster_key="leiden",
    species="Homo sapiens",
    cache_dir="~/.cache/gca/annotations",
    settings=settings,
)

print(result.stats)          # cache hits, batches, guardrail thresholds
print(result.report_dict())  # structured JSON-ready payload
```

## Seurat (R) via the CLI

```r
library(gptcellannotator)

# CLI-first: runs `gca annotate` under the hood
cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
gptca_config_set(cfg)

markers <- read.csv(system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator"))
annotations <- gptca_annotate_markers(markers, species = "Homo sapiens")
annotations$clusters
```

Set `offline = FALSE` and provide `base_url`/`api_key` only if you need the optional REST service.

## Architecture (brief)
- **LLM layer:** prompt-engineered OpenAI calls with schema validation; mock annotator for offline mode.
- **Retrieval:** optional top-k marker lookup against local databases to ground model responses.
- **Validation & guardrails:** overlap thresholds, downgrade to Unknown, and structured warnings per cluster.
- **Assets:** bundled demo marker tables, ortholog maps, and prompts materialised under `~/.cache/gpt-cell-annotator`.

## Documentation
- Installation and offline setup: `docs/install.md`
- Guided CLI walkthroughs: `docs/getting_started.md`
- Scanpy integration details: `docs/scanpy_integration.md`
- Seurat integration (R): `docs/seurat_integration.md`
- FAQs and configuration: `docs/faq.md`
