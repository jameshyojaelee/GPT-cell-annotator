# Getting Started

This guide gets you from install to your first annotations with the CLI, the Python helper, and the R wrapper.

## 1) Environment setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install "git+https://github.com/jameshyojaelee/CellAnnot-GPT.git#egg=gpt-cell-annotator"
```

If you prefer Poetry, run `poetry install --extras "dev"`.

## 2) Offline demo (CLI)

```bash
# Annotate demo PBMC markers and save the JSON report
gca annotate data/demo/pbmc_markers.csv \
  --species "Homo sapiens" \
  --offline \
  --out-json demo_annotations.json
```

The JSON report matches the schema used by the Python helper and the R wrapper. Add `--out-csv` to emit a flat summary.
Bundled assets are packaged with the library and materialised to `~/.cache/gpt-cell-annotator`
on first use, so the `data/demo/...` paths resolve even if the top-level `data/` directory is
empty. Use `data/` as scratch space for locally generated artifacts.

## 3) Programmatic usage (Python)

```python
from gpt_cell_annotator import annotate_markers

markers = {
    "0": ["MS4A1", "CD79A"],
    "1": ["CD3D", "CD3E"],
}

result = annotate_markers(markers, species="Homo sapiens")
print(result.report_dict())
```

The helper returns annotations plus a structured validation report.

## 4) Seurat / R usage

```r
library(gptcellannotator)
cfg <- gptca_config(offline = TRUE, cli_path = Sys.which("gca"))
gptca_config_set(cfg)

markers <- read.csv(system.file("extdata", "pbmc_markers.csv", package = "gptcellannotator"))
annotations <- gptca_annotate_markers(markers, species = "Homo sapiens")
annotations$clusters
```

Set `offline = FALSE` and supply `OPENAI_API_KEY` only if you need live LLM calls.

## 5) Validation guardrails

Key knobs (env vars or `config.settings.Settings` fields):

- `VALIDATION_MIN_MARKER_OVERLAP` (default 2) - minimum shared markers required to treat a label as supported.
- `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` - bump confidence tiers based on overlap.
- `SYNONYM_ENABLE_ORTHOLOGS` toggles ortholog mapping; `ORTHOLOG_MAPPING_PATH` overrides the mapping file.

Tune these before running the CLI or Python helper to adjust conservatism. Offline mode still applies guardrails against the bundled marker DB.

## 6) Troubleshooting

- Missing `marker_db.parquet` -> re-run the CLI once to re-materialise bundled assets, or set `GCA_MARKER_DB_PATH`.
- Empty annotations in offline mode -> ensure markers are parsed correctly (comma/JSON lists).

See `docs/faq.md` for more scenarios.
