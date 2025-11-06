# GPT Cell Annotator × Scanpy Architecture & Roadmap

This document inventories the current Scanpy/AnnData workflow, captures constraints, and stages the upcoming refactor that narrows scope while adding the reliability features requested in the revised mega-plan.

## Current Implementation Snapshot

### Annotate Pipeline (Python API)
- Entry point is `annotate_anndata` (`gpt_cell_annotator/scanpy.py`), which defers to the async variant after asserting we are not already inside an event loop.
- `annotate_anndata_async` optionally recomputes `rank_genes_groups` (via Scanpy) and builds a `{cluster_id -> markers}` payload from `adata.uns["rank_genes_groups"]["names"]`.
- Payloads flow into `_annotate_batches_async`, which:
  - Normalises batching (`BatchOptions`) and uses an `asyncio` semaphore for concurrency.
  - Checks an optional `AnnotationCacheProtocol` (default implementation: `DiskAnnotationCache`) before calling the `backend.llm.annotator.Annotator`.
  - Writes results back to the cache asynchronously and tallies batch/counter telemetry (Prometheus if available).
- `_run_annotation_workflow` feeds annotations through guardrail cross-checking (`backend.validation.crosscheck`) and builds a `DatasetReport` (`backend.validation.report`), returning stats (`cache_hits`, `llm_batches`, overrides, request_id).
- `_apply_annotations_to_obs` writes results into `adata.obs` (`gptca_label`, `gptca_proposed_label`, `gptca_status`, etc.) and mirrors structured metadata for downstream consumption.

### CLI Surface
- `python -m gpt_cell_annotator.scanpy` exposes two subcommands:
  - `annotate`: wraps the API flow above, handles I/O (h5ad/loom), configures presets, optional offline mode (`Annotator(force_mock=True)`), and emits JSON/CSV/statistics outputs.
  - `validate`: re-runs guardrail checks against existing annotations (`validate_anndata`).
- `gpt_cell_annotator/cli.py` registers these via the top-level `gca` CLI (`gca scanpy ...`), materialising bundled assets before delegating.
- Request scoping uses `GCA_REQUEST_ID` when present and enables structlog + Prometheus metrics when optional extras are installed.

### Test Coverage
- `tests/test_scanpy_integration.py` covers:
  - Sync/async parity (`annotate_anndata` vs. `annotate_anndata_async`).
  - CLI annotate/validate round-trips (ensuring outputs, offline flag behaviour).
  - Disk cache hit reduction, guardrail configuration overrides, helper wrappers (`annotate_from_markers`, `annotate_rank_genes`).
- Fixtures construct synthetic AnnData objects and stub the marker DB loader to avoid network/disk dependencies.

## System Diagram

```mermaid
flowchart LR
    A[AnnData Input] --> B{annotate_anndata}
    B -->|ensure rankings| C[scanpy.tl.rank_genes_groups]
    C --> D[Cluster Payload Builder]
    D --> E{DiskAnnotationCache?}
    E -->|hit| F1[Cached Annotation]
    E -->|miss| F2[Annotator.annotate_batch]
    F2 --> G[Guardrail crosscheck_batch]
    G --> H[DatasetReport]
    H --> I[Apply columns to adata.obs]
    H --> J[Summary JSON / CSV stats]
```

## Constraints & Compatibility

| Component | Current baseline | Proposed target | Notes |
| --- | --- | --- | --- |
| Python | `^3.11` (pyproject) | 3.9 – 3.12 matrix | Requires tox/nox orchestration and dependency pin audits. |
| AnnData | `^0.10.7` (direct dependency) | ≥0.9 (LTS) & latest minor | Back-compat verification for `.uns["rank_genes_groups"]` structure. |
| Scanpy | Optional `^1.10.1` extra | ≥1.9 track + latest | Loom IO + rank genes API need smoke tests across versions. |
| Backend services | Existing schema (`schemas/annotation.py`) | Unchanged | Telemetry hooks to remain optional extras. |

### Marker DB Asset Requirements

| Asset | Source / materialisation | Default location | Overrides |
| --- | --- | --- | --- |
| `marker_db.parquet` | `assets.ensure_marker_database` (downloads or bundled) | `settings.data_dir / marker_db.parquet` | `GCA_MARKER_DB_PATH` env var or CLI `--marker-db` |
| Marker cache | `_MARKER_CACHE` in-process map keyed by file mtime | In-memory per interpreter | Flushes only on restart; proposed disk cache will supersede. |
| Annotation cache | `DiskAnnotationCache` (async JSON files) | `cache_dir` supplied by caller | `GCA_CACHE_DIR` planned to pick defaults (`~/.cache/gpt-cell-annotator`). |
| Assets home | `assets.ensure_all_assets` | `~/.cache/gpt-cell-annotator` | `GPT_CELL_ANNOTATOR_ASSETS_HOME`, CLI `--assets-home`. |

## Adjusted Scope Summary

- **Removed (de-scoped)**: async coroutine entrypoint as default, helper wrappers (`annotate_rank_genes`, `annotate_from_markers`), structlog/Prometheus deliverable, standalone `gca scanpy validate` command.
- **Added (new focus)**: compatibility matrix & tox automation, migration guide, marker asset governance checklist, offline/caching first-class features, performance checks (cached vs. uncached).

## Phased Roadmap

### Phase 1 — Packaging & Compatibility
- Establish tox/nox matrix (Py3.9–3.12 × AnnData/Scanpy LTS/latest).
- Restructure extras (`scanpy`, `ui`, `dev`, optional `observability`) with lean defaults.
- Document compatibility guarantees, release checklist, and artifact smoke-tests.

### Phase 2 — Core API & Offline
- Refactor `annotate_anndata` to stay synchronous with internal batching (chunked threads) and relocate caching utilities to `gpt_cell_annotator/cache.py`.
- Enforce guardrail threshold configuration surface (no hidden structlog dependency).
- Introduce offline heuristic fallback in the returned `DatasetReport` and cache hit telemetry.

### Phase 3 — CLI, Docs, Benchmarking
- Collapse `validate` behaviour into `annotate --validate-only`.
- Expand CLI flags (`--chunk-size`, `--offline`, `--json-report`, `--species` presets) and integrate JSON schema output.
- Deliver updated docs (`docs/scanpy_integration.md`, migration guide) and performance benchmark notebook (50k synthetic cells, cached vs uncached).

### Phase 4 — Optional Telemetry & Dashboards
- Reintroduce observability behind an `observability` extra (structlog/Prometheus opt-in).
- Explore Streamlit dashboards / metrics export only if adopted.

## Issue Backlog (Draft Tickets)

### Phase 1
- `[P1-001]` Write Scanpy integration architecture doc & compatibility matrix (this prompt).
- `[P1-002]` Add tox matrix + `make packaging-check`.
- `[P1-003]` Update `RELEASING.md` with Python package checklist.

### Phase 2
- `[P2-001]` Refactor `annotate_anndata` batching to synchronous chunks + cache module.
- `[P2-002]` Implement offline heuristic fallback & typed `DatasetReport` wrapper.
- `[P2-003]` Expand integration tests for caching/offline/guardrail thresholds.

### Phase 3
- `[P3-001]` Merge `validate` into `annotate --validate-only` and update CLI help/tests.
- `[P3-002]` Produce benchmark notebook + doc updates + migration guide.
- `[P3-003]` Update CI workflow to run offline CLI smoke + nbmake benchmark.

### Phase 4
- `[P4-001]` Package observability extras and optional telemetry wiring.
- `[P4-002]` Prototype dashboard / metrics publishing (stretch).
