# gptcellannotator: Connect Seurat workflows to GPT Cell Annotator

An R interface to the GPT Cell Annotator backend. Provides helpers for
sending marker tables, working with Seurat objects, parsing validation
reports, and plotting annotations.

## Details

The package prefers HTTPS communication via the FastAPI service but can
fall back to the \`gca\` CLI when offline. Set configuration defaults
with
[`gptca_config()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_config.md)
and see
[`vignette("annotate-seurat")`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/articles/annotate-seurat.md)
for a complete workflow.

## Author

GPT Cell Annotator Team \<support@gpt-cell-annotator.org\>

## See also

[`gptca_config()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_config.md),
[`gptca_annotate_markers()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_annotate.md),
[`gptca_plot_umap()`](https://jameshyojaelee.github.io/CellAnnot-GPT/r/reference/gptca_annotate.md)
