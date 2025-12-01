# Frequently Asked Questions

## Do I need an OpenAI API key?

Only for live LLM annotations. Without `OPENAI_API_KEY`, the engine automatically uses the heuristic mock annotator. You can also force offline mode with `--offline` or `Annotator(force_mock=True)`.

## How do I run completely offline?

```bash
gca build-db --offline
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json
```

Assets are cached under `~/.cache/gpt-cell-annotator` by default. Override with `GPT_CELL_ANNOTATOR_HOME` or pass `--output-dir` to `gca build-db`.

## Where is the marker database?

The CLI copies `marker_db.parquet` to `~/.cache/gpt-cell-annotator/data/processed/`. Point at a different DB with:

- `GCA_MARKER_DB_PATH=/path/to/marker_db.parquet`
- `gca build-db --output-dir /my/db`
- `annotate_anndata(..., marker_db_path="...")`

## How do guardrails work?

Validation checks overlap between suggested markers and the knowledge base. Key knobs:

- `VALIDATION_MIN_MARKER_OVERLAP` (default 2) downgrades labels with low support to `Unknown or Novel`.
- `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` raise confidence tiers with more overlap.
- `RAG_ENABLED`, `RAG_TOP_K`, `RAG_MIN_OVERLAP` gate retrieval augmentation.

Tune via environment variables before running the CLI or Scanpy helpers.

## What about cross-species datasets?

Provide `--species` (CLI) or `dataset_context["species"]` (Python). Ortholog mapping uses the packaged file at `config/orthologs/human_mouse.tsv`; override with `ORTHOLOG_MAPPING_PATH` or disable via `SYNONYM_ENABLE_ORTHOLOGS=false` for debugging.

## How does the R package work offline?

`gptcellannotator` shells out to `gca annotate` when `offline = TRUE` (default). Ensure the Python CLI is installed and `Sys.which("gca")` returns a path. Set `offline = FALSE` with `base_url`/`api_key` only when targeting an external REST service.

## CLI says “marker_db.parquet not found”

Re-run `gca build-db --offline` or point to an existing DB with `GCA_MARKER_DB_PATH` / `--marker-db`.

## Scanpy import errors

Install the Scanpy extra: `pip install "gpt-cell-annotator[scanpy]"`. Keep `PYTHONPATH` empty so the virtualenv takes precedence.

## Cached assets look stale

Delete the cache directory (`rm -rf ~/.cache/gpt-cell-annotator`) or set `GPT_CELL_ANNOTATOR_HOME` to a fresh location before rerunning `gca build-db`.
