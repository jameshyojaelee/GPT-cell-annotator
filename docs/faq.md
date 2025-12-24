# Frequently Asked Questions

## Do I need an OpenAI API key?

Only for live LLM annotations. Without `OPENAI_API_KEY`, the engine automatically uses the heuristic mock annotator. You can also force offline mode with `--offline` or `Annotator(force_mock=True)`.

## How do I run completely offline?

```bash
gca annotate data/demo/pbmc_markers.csv --offline --out-json annotations.json
```

Assets are cached under `~/.cache/gpt-cell-annotator` by default. Override with `GPT_CELL_ANNOTATOR_HOME`.

## Where is the marker database?

The CLI copies `marker_db.parquet` to `~/.cache/gpt-cell-annotator/data/processed/`. Point at a different DB with:

- `GCA_MARKER_DB_PATH=/path/to/marker_db.parquet`
- `GPT_CELL_ANNOTATOR_DATA_DIR=/path/to/db`

## How do guardrails work?

Validation checks overlap between suggested markers and the knowledge base. Key knobs:

- `VALIDATION_MIN_MARKER_OVERLAP` (default 2) downgrades labels with low support to `Unknown or Novel`.
- `CONFIDENCE_OVERLAP_MEDIUM` / `CONFIDENCE_OVERLAP_HIGH` raise confidence tiers with more overlap.

Tune via environment variables before running the CLI or Python helper.

## What about cross-species datasets?

Provide `--species` (CLI) or `species=` (Python). Ortholog mapping uses the packaged file at `config/orthologs/human_mouse.tsv`; override with `ORTHOLOG_MAPPING_PATH` or disable via `SYNONYM_ENABLE_ORTHOLOGS=false` for debugging.

## CLI says "marker_db.parquet not found"

Re-run any CLI command to re-materialise the bundled assets, or point to an existing DB with `GCA_MARKER_DB_PATH`.

## Cached assets look stale

Delete the cache directory (`rm -rf ~/.cache/gpt-cell-annotator`) or set `GPT_CELL_ANNOTATOR_HOME` to a fresh location before rerunning the CLI.
