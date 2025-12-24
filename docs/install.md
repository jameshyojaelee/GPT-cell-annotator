# Installation Guide

## Prerequisites
- Python 3.9+ with `pip >= 23`
- Optional: R 4.3+ with Pandoc (for the Seurat vignette)

## Python installation (no PyPI required)

```bash
# Install directly from GitHub
pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"

# Or install from a local build artifact
#   poetry build
#   pip install dist/gpt_cell_annotator-*.whl

# Editable/development install
pip install -e ".[dev]"
```

The first CLI call will materialise bundled assets under `~/.cache/gpt-cell-annotator`. Override with:

- `GPT_CELL_ANNOTATOR_HOME=/path/to/cache` to relocate assets
- `GPT_CELL_ANNOTATOR_DATA_DIR=/path/to/db` to point at an existing marker database
- `GCA_MARKER_DB_PATH=/path/to/marker_db.parquet` to use a custom database

## CLI quick checks

```bash
# Offline demo (no network required)
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json
```

## R / Seurat installation (CLI-first)

```r
if (!require("pak")) install.packages("pak")
pak::pkg_install("github::jameshyojaelee/CellAnnot-GPT@main?subdir=clients/r/gptcellannotator")

# Ensure the Python CLI is available for offline mode
pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"
```

The R package shells out to `gca annotate` by default. Set `offline = FALSE` and supply `OPENAI_API_KEY` only when using live LLM calls.

## Working from source

```bash
git clone https://github.com/jameshyojaelee/CellAnnot-GPT.git
cd CellAnnot-GPT
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

Run `gca annotate data/demo/pbmc_markers.csv --offline` to confirm the environment.
