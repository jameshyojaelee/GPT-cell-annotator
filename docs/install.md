# Installation Guide

## Prerequisites
- Python 3.11 with `pip >= 23`
- Optional: Scanpy stack for AnnData workflows (`pip install "gpt-cell-annotator[scanpy]"`)
- Optional: R 4.3+ with Pandoc (for the Seurat vignette)

## Python installation

```bash
# Core CLI + engine
pip install gpt-cell-annotator

# Add Scanpy helpers (recommended for AnnData workflows)
pip install "gpt-cell-annotator[scanpy]"

# Developer tooling
pip install "gpt-cell-annotator[dev]"
```

The first CLI call will materialise bundled assets under `~/.cache/gpt-cell-annotator`. Override with:

- `GPT_CELL_ANNOTATOR_HOME=/path/to/cache` to relocate all assets
- `GPT_CELL_ANNOTATOR_DATA_DIR=/path/to/db` to point at an existing marker database
- `GCA_CACHE_DIR=/path/to/cache` for Scanpy annotation caches
- `GCA_MARKER_DB_PATH=/path/to/marker_db.parquet` to use a custom database

## CLI quick checks

```bash
# Offline demo (no network required)
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json

# Build marker DB artifacts locally
gca build-db --offline --output-dir ~/.cache/gca/db

# Scanpy helper (requires the scanpy extra)
gca scanpy annotate data/demo/pbmc_demo.h5ad --cluster-key leiden --species "Homo sapiens"
```

## R / Seurat installation

```r
if (!require("pak")) install.packages("pak")
pak::pkg_install("github::jameshyojaelee/CellAnnot-GPT@main?subdir=clients/r/gptcellannotator")

# Ensure the Python CLI is available for offline mode
pip install gpt-cell-annotator
```

The R package shells out to `gca annotate` by default. Set `offline = FALSE` and supply `base_url`/`api_key` only when using an external REST service.

## Working from source

```bash
git clone https://github.com/jameshyojaelee/CellAnnot-GPT.git
cd CellAnnot-GPT
python -m venv .venv && source .venv/bin/activate
pip install -e ".[scanpy,dev]"
```

Run `gca annotate data/demo/pbmc_markers.csv --offline` to confirm the environment.
