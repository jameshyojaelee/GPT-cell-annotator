# GPT-Cell-Annotator

GPT-Cell-Annotator is a minimal CLI toolkit (with Python + R wrappers) for labeling single-cell RNA-seq clusters from marker gene lists. It pairs a curated marker knowledge base with prompt-engineered LLM calls and lightweight validation guardrails. Everything ships with offline demo assets and works without a backend; live LLM calls are opt-in.

## What it does
- Annotate marker CSVs from the CLI (`gca annotate`); offline mock mode is default, live OpenAI calls run when `OPENAI_API_KEY` is set.
- Validate labels against the bundled marker DB and downgrade uncertain clusters to `Unknown or Novel`.
- Provide a small Python API (`annotate_markers`) for in-notebook workflows.
- Bridge to Seurat via the R package (`gptcellannotator`) which shells out to the CLI.

## Quick start

1) Install (GitHub; no PyPI needed):

```bash
pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"
```

Or install from a local build artifact:

```bash
# 1) poetry build
# 2) pip install dist/gpt_cell_annotator-*.whl
```

For editable dev installs:

```bash
pip install -e ".[dev]"
```

2) Try the offline demo (uses bundled assets + heuristic mock annotator):

```bash
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json
```

3) Annotate from Python:

```python
from gpt_cell_annotator import annotate_markers

clusters = {
    "0": ["MS4A1", "CD79A"],
    "1": ["CD3D", "CD3E"],
}

result = annotate_markers(clusters, species="Homo sapiens")
print(result.report_dict())
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

Set `offline = FALSE` and provide `OPENAI_API_KEY` only if you need live LLM calls.

## Bundled demos
- Marker CSV: `data/demo/pbmc_markers.csv` (used in `gca annotate` examples).

Bundled assets live under `src/gpt_cell_annotator/_assets` and are materialised to
`~/.cache/gpt-cell-annotator` automatically when you run the CLI. The top-level `data/`
directory is reserved for locally generated artifacts.

## Architecture (brief)
- **LLM layer:** prompt-engineered OpenAI calls with schema validation; mock annotator for offline mode.
- **Validation & guardrails:** overlap thresholds, downgrade to Unknown, and structured warnings per cluster.
- **Assets:** bundled demo marker tables and a small marker DB materialised under `~/.cache/gpt-cell-annotator`.

## Documentation
- Installation and offline setup: `docs/install.md`
- Guided CLI walkthroughs: `docs/getting_started.md`
- Seurat integration (R): `docs/seurat_integration.md`
- FAQs and configuration: `docs/faq.md`
- Development notes: `docs/development.md`
