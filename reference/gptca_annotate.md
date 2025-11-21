# Annotate clusters with GPT Cell Annotator

High-level helpers to send marker tables to GPT Cell Annotator, merge
validated labels into Seurat metadata, and visualise annotations on UMAP
embeddings.

## Usage

``` r
gptca_annotate_markers(
  markers,
  species = NULL,
  tissue = NULL,
  return_validated = TRUE,
  marker_limit = 50,
  offline = NULL,
  fallback = TRUE,
  config = gptca_config_get(),
  cluster_column = NULL,
  gene_column = NULL,
  markers_column = NULL,
  score_column = NULL,
  top_n = NULL
)

gptca_annotate_seurat(
  object,
  markers = NULL,
  species = NULL,
  tissue = NULL,
  return_validated = TRUE,
  marker_limit = 50,
  cluster_column = "cluster",
  gene_column = "gene",
  score_column = NULL,
  top_n = 25,
  offline = NULL,
  fallback = TRUE,
  config = gptca_config_get()
)

gptca_add_metadata(
  object,
  annotation,
  cluster_column = "seurat_clusters",
  label_col = "gptca_label",
  confidence_col = "gptca_confidence",
  status_col = "gptca_status",
  warnings_col = "gptca_warnings"
)

gptca_plot_umap(
  object,
  label_col = "gptca_label",
  status_col = "gptca_status",
  reduction = "umap",
  flagged_status = c("flagged", "unknown"),
  point_size = 0.6,
  alpha = 0.9,
  palette = NULL
)
```

## Arguments

- markers:

  Marker definitions: named list of gene symbols, a data frame with
  \`cluster\` and \`gene\` columns (e.g. `FindAllMarkers` output), or a
  data frame with \`cluster_id\` and list-column \`markers\`.

- species:

  Optional species hint included in the dataset context.

- tissue:

  Optional tissue hint.

- return_validated:

  When `TRUE` (default) request the validated report. Set to `FALSE` to
  retrieve raw LLM output for a single cluster.

- marker_limit:

  Maximum number of markers per cluster.

- offline:

  Logical flag forcing CLI fallback. Defaults to the active
  configuration.

- fallback:

  Attempt CLI fallback automatically when HTTP requests fail.

- config:

  A `GptcaConfig` object.

- cluster_column:

  Name of the column containing cluster IDs in marker data frames or
  Seurat metadata.

- gene_column:

  Name of the gene column in marker data frames.

- markers_column:

  Name of list-column storing marker vectors.

- score_column:

  Column used to rank markers (highest values first).

- top_n:

  Optional number of top markers to keep per cluster before applying
  `marker_limit`.

- object:

  A Seurat object.

- label_col:

  Metadata column used to store primary labels.

- confidence_col:

  Metadata column for confidence values.

- status_col:

  Metadata column storing validation status.

- warnings_col:

  Metadata column storing concatenated warnings.

- annotation:

  A `gptca_annotation` object returned by `gptca_annotate_markers()` or
  `gptca_annotate_seurat()`.

- reduction:

  Name of the reduction containing UMAP embeddings.

- flagged_status:

  Status values highlighted on the plot.

- point_size:

  Point size for scatter plot.

- alpha:

  Point alpha level.

- palette:

  Optional colour palette overriding defaults.

## Value

`gptca_annotate_markers()` and `gptca_annotate_seurat()` return objects
of class `gptca_annotation`. `gptca_add_metadata()` returns the modified
Seurat object. `gptca_plot_umap()` returns a ggplot2 object.

## Examples

``` r
## Not run: 
# markers <- read.csv(system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator"))
# annotations <- gptca_annotate_markers(markers, species = "Homo sapiens")
# if (requireNamespace("Seurat", quietly = TRUE)) {
#   pbmc <- readRDS("path/to/pbmc.rds")
#   pbmc <- gptca_add_metadata(pbmc, annotations)
#   gptca_plot_umap(pbmc)
# }
## End(Not run)
```
