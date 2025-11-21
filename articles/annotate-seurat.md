# Annotate Seurat Clusters with GPT Cell Annotator

## Overview

`gptcellannotator` provides an end-to-end bridge from Seurat marker
discovery to the GPT Cell Annotator engine. It is CLI-first (offline by
default) and can call an optional REST endpoint when provided with a
base URL and API key. This vignette walks through preparing marker
tables, invoking the annotator, applying validation guardrails, and
visualising results on a UMAP embedding.

## Setup

``` r
library(gptcellannotator)
library(Seurat)
library(tibble)
```

Configure the CLI path once per session. Environment variables
(`GPTCA_CLI_PATH`) take precedence, but explicit configuration is
clearer in reproducible scripts.

``` r
cfg <- gptca_config(
  offline = TRUE,
  cli_path = Sys.which("gca")
)
gptca_config_set(cfg)
```

## Load markers

Start with a Seurat object that has undergone clustering. Here we load
the demo PBMC object distributed with Seurat, then compute positive
markers per cluster via `FindAllMarkers`. The function returns a long
data frame with columns `gene`, `cluster`, `avg_log2FC`, and more.

``` r
pbmc <- readRDS("path/to/pbmc.rds")
markers <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25)
```

The first rows look like:

``` r
markers %>% as_tibble() %>% dplyr::group_by(cluster) %>% dplyr::slice_head(n = 3)
```

## Submit markers for annotation

[`gptca_annotate_seurat()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_annotate.md)
handles conversion to the backend payload, including limiting the top
markers per cluster and sorting by `avg_log2FC`.

``` r
annotations <- gptca_annotate_seurat(
  pbmc,
  markers = markers,
  species = "Homo sapiens",
  tissue = "Peripheral blood",
  return_validated = TRUE,
  top_n = 25
)
```

The result is a `gptca_annotation` object with convenient accessors:

``` r
annotations$clusters
annotations$summary
```

- `clusters` – tidy tibble with per-cluster labels, ontology IDs,
  validation status, and collected warnings.
- `summary` – dataset-level support metrics from the validation report.
- `raw` – full JSON payload if you need access to every nested field.

## Merge results into Seurat metadata

Annotator outputs should become part of your object metadata so
downstream plots and differential analyses can use the curated labels.
[`gptca_add_metadata()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_annotate.md)
adds a set of columns (`gptca_label`, `gptca_confidence`,
`gptca_status`, `gptca_warnings`) keyed by the existing clustering
column.

``` r
pbmc <- gptca_add_metadata(pbmc, annotations)
head(pbmc@meta.data[, c("seurat_clusters", "gptca_label", "gptca_status")])
```

Clusters flagged during validation will have `gptca_status` set to
`"flagged"` or `"unknown"`, and the warnings column retains the backend
rationale.

## Visualise on UMAP

[`gptca_plot_umap()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_annotate.md)
renders a ggplot built from the stored embeddings. Flagged clusters are
outlined in black to emphasise groups requiring manual review.

``` r
gptca_plot_umap(pbmc)
```

Because the plot is built with `ggplot2` you can further customise
scales, themes, or overlays using standard syntax.

## Optional REST workflow

To target an external REST deployment instead of the CLI, supply a base
URL and API key:

``` r
rest_cfg <- gptca_config(
  offline = FALSE,
  base_url = "https://api.gpt-cell-annotator.org",
  api_key = Sys.getenv("GPTCA_API_KEY")
)
annotations <- gptca_annotate_markers(
  markers,
  species = "Homo sapiens",
  config = rest_cfg
)
```

The CLI path is automatically discovered from `GPTCA_CLI_PATH` or the
system `PATH`. Both modes apply the same validation guardrails and
return the same schema.

## Troubleshooting

- **`gca` binary not found** – install the Python package and ensure
  `Sys.which("gca")` resolves the CLI; override with `GPTCA_CLI_PATH`.
- **Validation disagreements** – inspect
  `annotations$clusters$validation` for missing markers or contradictory
  evidence. You can also access the raw JSON via `annotations$raw`.
- **REST-specific errors (401/429)** – only when `offline = FALSE`.
  Check API credentials and retry settings in
  [`gptca_config()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_config.md).

For more recipes and automation tips, visit the pkgdown site once
published or open an issue on GitHub.
