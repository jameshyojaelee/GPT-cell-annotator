# gptcellannotator

`gptcellannotator` connects Seurat workflows to the GPT Cell Annotator
CLI. It is CLI-first (offline by default), with validation guardrails
and plotting helpers.

## Installation

Ensure you have an R 4.3 environment with development headers for
`xml2`, plus Pandoc \>= 2.11 for vignettes. On macOS install the Xcode
command line tools; on Windows install
[Rtools42](https://cran.r-project.org/bin/windows/Rtools/).

### GitHub source (pak or remotes)

``` r
# pak automatically resolves system + R dependencies
if (!require("pak")) install.packages("pak")
pak::pkg_install("github::jameshyojaelee/CellAnnot-GPT@main?subdir=clients/r/gptcellannotator")

# or use remotes / devtools
remotes::install_github("jameshyojaelee/CellAnnot-GPT", subdir = "clients/r/gptcellannotator", build = TRUE)
```

Run the helper script if you prefer a single command:

``` bash
Rscript clients/r/scripts/install_github.R
```

### Local source tarball

``` bash
R CMD build clients/r/gptcellannotator
R CMD INSTALL gptcellannotator_0.1.0.tar.gz
```

The script `clients/r/scripts/install_local_tarball.R` automates those
steps with
[`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html)
while caching dependencies.

### Python CLI

Offline annotations require the Python CLI:

``` bash
pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"
```

## Quick start

``` r
library(gptcellannotator)

cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
gptca_config_set(cfg)

pbmc_markers <- system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator")
markers <- read.csv(pbmc_markers)

annotations <- gptca_annotate_markers(
  markers,
  species = "Homo sapiens",
  return_validated = TRUE
)

annotations$clusters
```

### Seurat pipeline

``` r
library(Seurat)

pbmc <- readRDS("path/to/pbmc.rds")            # replace with your Seurat object
markers <- FindAllMarkers(pbmc)                # positive markers only

annotations <- gptca_annotate_seurat(
  pbmc,
  markers = markers,
  species = "Homo sapiens"
)

pbmc <- gptca_add_metadata(pbmc, annotations)
gptca_plot_umap(pbmc)
```

### Offline fallback

Set `offline = TRUE` (default) or export `GPTCA_CLI_PATH` pointing to
the `gca` binary. The package writes markers to a temp CSV, calls
`gca annotate`, and parses the structured JSON output.

``` r
cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
annotations <- gptca_annotate_markers(markers, config = cfg)
```

## Configuration

- `GPTCA_CLI_PATH` to point at the `gca` binary (auto-detected
  otherwise).
- All exported functions accept a `config` override for ad-hoc sessions.

## Validation and guardrails

Validated responses include:

- Per-cluster `status` (`supported`, `flagged`, `unknown`)
- Canonical markers and ontology identifiers
- Validation notes (`missing_markers`, `contradictory_markers`)
- Warning strings surfaced in Seurat metadata

Use `annotations$raw` to inspect the full JSON payload.

## Getting help

- Open an issue:
  <https://github.com/jameshyojaelee/CellAnnot-GPT/issues>
- Internal docs:
  [`vignette("annotate-seurat", package = "gptcellannotator")`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/articles/annotate-seurat.md)
