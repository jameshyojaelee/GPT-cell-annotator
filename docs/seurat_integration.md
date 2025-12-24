# Seurat Integration Guide

Use the `gptcellannotator` R package to run GPT Cell Annotator from Seurat pipelines. The package shells out to the Python CLI so offline runs match the Python behaviour; live LLM calls are opt-in via `OPENAI_API_KEY`.

## Prerequisites
- R >= 4.3 with dev headers for `xml2`
- Pandoc >= 2.11 (bundled with RStudio)
- Seurat >= 5.0
- Python CLI installed and on `PATH` (`pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"`)

## Install the R package

```r
if (!require("pak")) install.packages("pak")
pak::pkg_install("github::jameshyojaelee/CellAnnot-GPT@main?subdir=clients/r/gptcellannotator")
```

Helper scripts:
- `Rscript clients/r/scripts/install_github.R [ref]`
- `Rscript clients/r/scripts/install_local_tarball.R`

## Quickstart workflow

```r
library(Seurat)
library(gptcellannotator)

# CLI-first: uses `gca annotate` under the hood
cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
gptca_config_set(cfg)

pbmc <- readRDS("pbmc_seurat.rds")
markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)

annotations <- gptca_annotate_seurat(
  pbmc,
  markers = markers,
  species = "Homo sapiens",
  return_validated = TRUE
)

pbmc <- gptca_add_metadata(pbmc, annotations)
gptca_plot_umap(pbmc)
```

Set `offline = FALSE` and supply `OPENAI_API_KEY` if you want live LLM calls.

## Outputs
- `annotations$clusters`: tidy tibble with labels, ontology IDs, status (`supported`, `flagged`, `unknown`), and warnings.
- `pbmc@meta.data`: added columns such as `gptca_label`, `gptca_confidence`, `gptca_status`, and validation notes.
- `gptca_plot_umap()`: ggplot overlay of annotations; flagged clusters are outlined for review.

## Troubleshooting

| Symptom | Resolution |
| --- | --- |
| `gca` binary not found | Install the Python package and ensure `Sys.which("gca")` returns a path; otherwise set `GPTCA_CLI_PATH`. |
| Vignette build fails (Pandoc) | Install Pandoc >= 2.11 or point `RSTUDIO_PANDOC` at your installation. |
| Empty annotations | Increase marker counts (`top_n`) and confirm positive markers from `FindAllMarkers()`. |
| Want live LLM calls | Set `OPENAI_API_KEY` in your shell before running R. |
