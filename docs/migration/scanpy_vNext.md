# Migrating to the Scanpy vNext Workflow

This guide highlights the breaking/API changes introduced while delivering the Phase 2 (core API) and Phase 3 (CLI/docs) work for the Scanpy integration.

## TL;DR

- `annotate_anndata` is now purely synchronous and returns a `ScanpyAnnotationResult` whose `.report` is a `ScanpyDatasetReport` (with cache hits, guardrail overrides, offline status, and chunk size metadata).
- The asynchronous helper (`annotate_anndata_async`) and concurrency-oriented `BatchOptions(size=…, concurrency=…)` are gone. Chunk behaviour is controlled by `BatchOptions(chunk_size=…)`.
- `DiskAnnotationCache` uses synchronous `get/set` methods. Marker DB lookups go through the new `MarkerDatabaseCache`, which persists copies under `~/.cache/gpt-cell-annotator` (or `GCA_CACHE_DIR`).
- CLI users now rely on `gca scanpy annotate ... --validate-only` instead of a separate `validate` subcommand. New flags: `--chunk-size`, `--json-report`, `--use-cache`, `--progress`.

## API Changes

| Old | New | Notes |
| --- | --- | --- |
| `annotate_anndata_async` | Removed | Call `annotate_anndata` directly; it no longer enforces “no event loop” restrictions. |
| `BatchOptions(size=32, concurrency=1)` | `BatchOptions(chunk_size=32)` | Chunk size limits how many clusters are sent to the LLM per request. |
| `ScanpyAnnotationResult.report` (`DatasetReport`) + `.stats` | `ScanpyAnnotationResult.report` (`ScanpyDatasetReport`) + `.stats` proxy | Access guardrail/cache information via `.report` or `result.stats`. |
| `DiskAnnotationCache.get/set` (async) | `DiskAnnotationCache.get/set` (sync) | Update custom caches accordingly. |
| `Marker DB` memoisation in-module | `MarkerDatabaseCache` in `gpt_cell_annotator.cache` | Reuse this helper to respect `GCA_CACHE_DIR`. |

### Guardrails & Offline Mode

- `ScanpyDatasetReport.offline_mode` indicates whether the annotator fell back to the heuristic mock (e.g., no `OPENAI_API_KEY`).
- Guardrail overrides (via `GuardrailConfig`) are recorded in the report and echoed through logging (`scanpy.annotate.guardrails.override`).

## CLI Changes

| Old command | Replacement |
| --- | --- |
| `gca scanpy annotate ... --summary-json out.json --stats-json out_stats.json` | `gca scanpy annotate ... --json-report out.json --summary-csv out.csv` |
| `gca scanpy validate ...` | `gca scanpy annotate ... --validate-only --label-column existing_label` |
| `--batch-size`, `--concurrency` | `--chunk-size` |

Additional niceties:
- `--use-cache` scaffolds `~/.cache/gpt-cell-annotator/annotations` automatically.
- `--progress` shows a `tqdm` progress bar when stderr is attached to a TTY (optional dependency).

## Recommended Migration Steps

1. **Update imports**: replace `BatchOptions(size=..., concurrency=...)` with `BatchOptions(chunk_size=...)`. Drop references to `annotate_anndata_async`.
2. **Adjust tests**: expect `result.report.dataset` (instead of raw `DatasetReport`) and leverage fixture helpers in `tests/fixtures/scanpy.py`.
3. **CLI scripts**: swap `gca scanpy validate` invocations for `gca scanpy annotate ... --validate-only`. Update automation to look for `--json-report`.
4. **Caching**: point custom tooling at `gpt_cell_annotator.cache.MarkerDatabaseCache` for marker DB reuse. Honour `GCA_CACHE_DIR` when writing annotation caches.
5. **Docs/Notebooks**: refresh code snippets to use the new API and mention `ScanpyDatasetReport` in narrative text.

## Compatibility Notes

- Python 3.9–3.12, AnnData ≥0.9, Scanpy ≥1.9 are covered in the new tox matrix (`tox.ini`).
- Optional telemetry (structlog, Prometheus) now lives behind the `observability` extra; the default flow sticks to stdlib logging.
- Offline runs (no `OPENAI_API_KEY` or explicit `--offline`) are now surfaced via `.report.offline_mode` and log records.

Need help adapting downstream automation? Open an issue with sample code and we’ll include it in the migration examples.
