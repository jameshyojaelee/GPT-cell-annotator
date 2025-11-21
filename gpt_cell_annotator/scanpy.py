"""Scanpy integration helpers, annotation APIs, and CLI entrypoints."""

from __future__ import annotations

import argparse
import json
import logging
import os
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData
from pydantic import BaseModel, ConfigDict, Field

from backend.llm.annotator import Annotator
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import DatasetReport, build_structured_report
from config.settings import Settings, get_settings
from gpt_cell_annotator import assets
from gpt_cell_annotator.cache import MarkerDatabaseCache, default_cache_dir

try:  # Optional dependency used for Loom IO.
    import scanpy as sc
except ImportError:  # pragma: no cover - handled via guards in runtime paths.
    sc = None

try:  # Progress bars are optional; fallback to logging when missing.
    from tqdm import tqdm
except ImportError:  # pragma: no cover - tqdm is an optional extra.
    tqdm = None  # type: ignore[assignment]

logger = logging.getLogger("gpt_cell_annotator.scanpy")

MARKER_DB_COLUMNS = [
    "source",
    "cell_type",
    "ontology_id",
    "gene_symbol",
    "species",
    "tissue",
    "evidence",
    "reference",
    "evidence_score",
]

MARKER_DB_ENV_VAR = "GCA_MARKER_DB_PATH"
REQUEST_ID_ENV_VAR = "GCA_REQUEST_ID"
CACHE_VERSION = "scanpy-v3"

SPECIES_PRESETS: dict[str, dict[str, str | None]] = {
    "human_pbmc": {"species": "Homo sapiens", "tissue": "Peripheral blood"},
    "mouse_pbmc": {"species": "Mus musculus", "tissue": "Peripheral blood"},
    "mouse_brain": {"species": "Mus musculus", "tissue": "Brain"},
}

_MARKER_CACHE_SINGLETON: MarkerDatabaseCache | None = None


def _marker_cache() -> MarkerDatabaseCache:
    global _MARKER_CACHE_SINGLETON
    if _MARKER_CACHE_SINGLETON is None:
        _MARKER_CACHE_SINGLETON = MarkerDatabaseCache()
    return _MARKER_CACHE_SINGLETON


class ScanpyAnnotationError(RuntimeError):
    """Base exception for Scanpy integration failures."""


class MarkerDatabaseError(ScanpyAnnotationError):
    """Raised when the marker database cannot be resolved or loaded."""


def _default_request_id(user_supplied: str | None) -> str:
    env_override = os.environ.get(REQUEST_ID_ENV_VAR)
    if user_supplied:
        return user_supplied
    if env_override:
        return env_override
    return f"scanpy-{os.getpid()}"


@dataclass(slots=True)
class DiskAnnotationCache:
    """Disk-backed cache for Scanpy annotations (synchronous)."""

    root: Path

    def __post_init__(self) -> None:
        self.root = self.root.expanduser()
        self.root.mkdir(parents=True, exist_ok=True)

    def _path_for(self, payload: dict[str, Any]) -> Path:
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        digest = sha256(canonical.encode("utf-8")).hexdigest()
        return self.root / f"{digest}.json"

    def get(self, payload: dict[str, Any]) -> dict[str, Any] | None:
        path = self._path_for(payload)
        if not path.exists():
            return None
        try:
            loaded: Any = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if isinstance(loaded, dict):
            return loaded
        return None

    def set(self, payload: dict[str, Any], value: dict[str, Any]) -> None:
        path = self._path_for(payload)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value), encoding="utf-8")


@dataclass(slots=True)
class GuardrailConfig:
    """Overrides for guardrail behaviour and confidence calibration."""

    min_marker_overlap: int | None = None
    force_unknown_on_fail: bool | None = None
    confidence_overlap_medium: int | None = None
    confidence_overlap_high: int | None = None

    def as_update(self) -> dict[str, Any]:
        payload: dict[str, Any] = {}
        if self.min_marker_overlap is not None:
            payload["validation_min_marker_overlap"] = max(0, self.min_marker_overlap)
        if self.force_unknown_on_fail is not None:
            payload["validation_force_unknown_on_fail"] = bool(self.force_unknown_on_fail)
        if self.confidence_overlap_medium is not None:
            payload["confidence_overlap_medium"] = max(0, self.confidence_overlap_medium)
        if self.confidence_overlap_high is not None:
            payload["confidence_overlap_high"] = max(
                self.confidence_overlap_medium or 0,
                self.confidence_overlap_high,
            )
        return payload


@dataclass(slots=True)
class BatchOptions:
    """Controls chunk size for annotation batches."""

    chunk_size: int = 32

    def normalized(self, total: int) -> int:
        if total <= 0:
            return 1
        return max(1, min(self.chunk_size, total))


class ScanpyDatasetReport(BaseModel):
    """Structured report summarising annotations plus runtime metadata."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    dataset: DatasetReport
    cache_hits: int = Field(default=0, ge=0)
    llm_batches: int = Field(default=0, ge=0)
    offline_mode: bool = False
    guardrail_min_overlap: int = Field(default=1, ge=0)
    guardrail_overrides: dict[str, Any] = Field(default_factory=dict)
    request_id: str = ""
    chunk_size: int = Field(default=1, ge=1)

    @property
    def warnings(self) -> list[str]:
        collected: list[str] = []
        for cluster in self.dataset.clusters:
            collected.extend(cluster.warnings or [])
        return collected

    def stats_dict(self) -> dict[str, Any]:
        return {
            "total_clusters": self.dataset.summary.total_clusters,
            "cache_hits": self.cache_hits,
            "llm_batches": self.llm_batches,
            "offline_mode": self.offline_mode,
            "guardrail_min_marker_overlap": self.guardrail_min_overlap,
            "chunk_size": self.chunk_size,
            "mode": "offline" if self.offline_mode else "live",
        }


@dataclass(slots=True)
class ScanpyAnnotationResult:
    """Return object from annotate_anndata."""

    adata: AnnData
    annotations: list[dict[str, Any]]
    report: ScanpyDatasetReport

    @property
    def stats(self) -> dict[str, Any]:
        return self.report.stats_dict()

    def report_dict(self) -> dict[str, Any]:
        payload = self.report.model_dump()
        payload["dataset"] = self.report.dataset.model_dump()
        return payload


@dataclass(slots=True)
class MarkerAnnotationResult:
    """Result returned when annotating pre-computed marker collections."""

    annotations: list[dict[str, Any]]
    report: ScanpyDatasetReport

    def report_dict(self) -> dict[str, Any]:
        payload = self.report.model_dump()
        payload["dataset"] = self.report.dataset.model_dump()
        return payload


def _unique_ordered(markers: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for marker in markers:
        if not isinstance(marker, str):
            continue
        normalised = marker.strip()
        if not normalised:
            continue
        upper = normalised.upper()
        if upper in seen:
            continue
        seen.add(upper)
        ordered.append(upper)
    return ordered


def _ensure_rankings(
    adata: AnnData,
    cluster_key: str,
    *,
    top_n_markers: int,
    method: str = "wilcoxon",
) -> None:
    rankings = adata.uns.get("rank_genes_groups")
    params = (rankings or {}).get("params", {})
    if rankings is not None and params.get("groupby") == cluster_key:
        return

    if sc is None:
        raise ImportError(
            "scanpy is required to compute marker rankings. Install the "
            'extra with `pip install "gpt-cell-annotator[scanpy]"` or '
            "precompute `rank_genes_groups` before calling annotate_anndata."
        )

    sc.tl.rank_genes_groups(adata, groupby=cluster_key, n_genes=top_n_markers, method=method)


def _markers_from_rankings(
    rankings: Mapping[str, Sequence[str]] | np.ndarray | Mapping[str, Any],
    top_n_markers: int,
) -> dict[str, list[str]]:
    if isinstance(rankings, np.ndarray) and rankings.dtype.names:
        groups = rankings.dtype.names
        structured: dict[str, list[str]] = {}
        for group in groups:
            structured[str(group)] = _unique_ordered(rankings[group][:top_n_markers])
        return structured

    if isinstance(rankings, Mapping):
        derived: dict[str, list[str]] = {}
        for cluster, values in rankings.items():
            derived[str(cluster)] = _unique_ordered(values[:top_n_markers])
        return derived
    raise ValueError("Unsupported structure for rank_genes_groups['names'].")


def _extract_markers(adata: AnnData, cluster_key: str, top_n_markers: int) -> dict[str, list[str]]:
    rankings = adata.uns.get("rank_genes_groups")
    if rankings is None or "names" not in rankings:
        raise ValueError(
            "AnnData object is missing `rank_genes_groups` results for cluster key "
            f"'{cluster_key}'. Run scanpy.tl.rank_genes_groups beforehand or allow "
            "annotate_anndata to compute them."
        )
    return _markers_from_rankings(rankings["names"], top_n_markers)


def _build_cluster_index(cluster_series: pd.Series) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for obs_idx, cluster_id in cluster_series.items():
        mapping.setdefault(str(cluster_id), []).append(str(obs_idx))
    return mapping


def _apply_annotations_to_obs(
    adata: AnnData,
    cluster_key: str,
    result_prefix: str,
    report: DatasetReport,
) -> None:
    label_col = f"{result_prefix}_label"
    confidence_col = f"{result_prefix}_confidence"
    status_col = f"{result_prefix}_status"
    rationale_col = f"{result_prefix}_rationale"
    ontology_col = f"{result_prefix}_ontology_id"
    proposed_col = f"{result_prefix}_proposed_label"
    canonical_col = f"{result_prefix}_canonical_markers"
    mapping_col = f"{result_prefix}_mapping_notes"

    clusters_index = _build_cluster_index(adata.obs[cluster_key].astype(str))
    for column in [
        label_col,
        proposed_col,
        confidence_col,
        status_col,
        rationale_col,
        ontology_col,
        canonical_col,
        mapping_col,
    ]:
        adata.obs[column] = pd.Series(
            ["" for _ in range(adata.n_obs)],
            index=adata.obs.index,
            dtype="object",
        )

    for cluster_report in report.clusters:
        cluster_id = cluster_report.cluster_id
        annotation = cluster_report.annotation
        metadata = annotation.get("metadata") or {}
        mapping_notes = metadata.get("mapping_notes") or annotation.get("mapping_notes") or []
        canonical_markers = metadata.get("canonical_markers") or annotation.get("markers") or []
        mapping_display = " | ".join(
            f"{note.get('source', '?')}→{note.get('target') or 'unmapped'}"
            for note in mapping_notes
        )
        obs_indices = clusters_index.get(cluster_id, [])
        for obs_idx in obs_indices:
            adata.obs.at[obs_idx, label_col] = annotation.get("primary_label") or ""
            adata.obs.at[obs_idx, proposed_col] = annotation.get("proposed_label") or ""
            adata.obs.at[obs_idx, confidence_col] = cluster_report.confidence or ""
            adata.obs.at[obs_idx, status_col] = cluster_report.status or ""
            adata.obs.at[obs_idx, rationale_col] = annotation.get("rationale") or ""
            adata.obs.at[obs_idx, ontology_col] = annotation.get("ontology_id") or ""
            adata.obs.at[obs_idx, canonical_col] = ", ".join(canonical_markers)
            adata.obs.at[obs_idx, mapping_col] = mapping_display


def report_to_dataframe(
    report: DatasetReport | ScanpyDatasetReport | dict[str, Any],
) -> pd.DataFrame:
    """Return a dataframe summarising cluster annotations."""

    dataset: DatasetReport
    if isinstance(report, ScanpyDatasetReport):
        dataset = report.dataset
    elif isinstance(report, DatasetReport):
        dataset = report
    else:
        dataset = DatasetReport.model_validate(report)

    records: list[dict[str, Any]] = []
    for cluster in dataset.clusters:
        annotation = cluster.annotation
        records.append(
            {
                "cluster_id": cluster.cluster_id,
                "primary_label": annotation.get("primary_label"),
                "confidence": cluster.confidence,
                "status": cluster.status,
                "ontology_id": annotation.get("ontology_id"),
                "rationale": annotation.get("rationale"),
                "warnings": "; ".join(cluster.warnings or []),
            }
        )
    return pd.DataFrame.from_records(records)


def _resolve_marker_db_path(marker_db_path: Path | None) -> Path:
    if marker_db_path is not None:
        candidate = marker_db_path.expanduser()
        if candidate.exists():
            return candidate

    env_path = os.environ.get(MARKER_DB_ENV_VAR)
    if env_path:
        candidate = Path(env_path).expanduser()
        if candidate.exists():
            return candidate

    settings = get_settings()
    data_dir = Path(settings.data_dir)
    default_path = data_dir / "marker_db.parquet"
    if default_path.exists():
        return default_path

    ensured_dir = assets.ensure_marker_database(target_dir=data_dir)
    candidate = ensured_dir / "marker_db.parquet"
    if candidate.exists():
        return candidate
    fallback = assets.ensure_marker_database() / "marker_db.parquet"
    if fallback.exists():
        return fallback
    raise MarkerDatabaseError("marker_db.parquet could not be located or materialised.")


def _guardrail_settings(guardrails: GuardrailConfig | None) -> tuple[Settings | None, int]:
    if guardrails is None:
        settings = get_settings()
        return None, max(1, settings.validation_min_marker_overlap)

    update = guardrails.as_update()
    if not update:
        settings = get_settings()
        return None, max(1, settings.validation_min_marker_overlap)

    base = get_settings()
    override = base.model_copy(update=update)
    min_support = max(1, override.validation_min_marker_overlap)
    return override, min_support


def _build_cache_payload(
    cluster: dict[str, Any],
    dataset_context: dict[str, Any],
    *,
    mode: str,
) -> dict[str, Any]:
    return {
        "version": CACHE_VERSION,
        "cluster_id": str(cluster.get("cluster_id")),
        "markers": cluster.get("markers"),
        "context": dataset_context,
        "mode": mode,
    }


def _chunk(payload: Sequence[dict[str, Any]], size: int) -> Iterable[list[dict[str, Any]]]:
    for index in range(0, len(payload), size):
        yield list(payload[index : index + size])


def _build_annotations(
    clusters_payload: Sequence[dict[str, Any]],
    results: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    for cluster in clusters_payload:
        cluster_id = str(cluster["cluster_id"])
        cluster_result = dict(results.get(cluster_id) or {})
        cluster_result.setdefault("primary_label", "Unknown or Novel")
        cluster_result.setdefault("confidence", "Unknown")
        cluster_result.setdefault("rationale", "")
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cluster_id,
                "markers": cluster.get("markers", []),
            }
        )
    return annotations


def _annotate_clusters_payload(
    clusters_payload: Sequence[dict[str, Any]],
    *,
    species: str,
    tissue: str | None = None,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: DiskAnnotationCache | None = None,
    request_id: str | None = None,
    progress_callback: Callable[[int], None] | None = None,
) -> tuple[list[dict[str, Any]], ScanpyDatasetReport]:
    """Shared annotation pipeline for pre-built cluster payloads."""

    annotator_instance = annotator or Annotator()
    offline_mode = annotator_instance.llm_mode != "live"
    if offline_mode:
        logger.info("scanpy.annotate.offline", extra={"reason": "llm_mode=mock"})

    dataset_context: dict[str, str] = {"species": species}
    if tissue:
        dataset_context["tissue"] = tissue

    marker_cache = _marker_cache()
    resolved_path = Path(marker_db_path).expanduser() if marker_db_path is not None else None
    try:
        marker_df = marker_cache.load(
            frame=marker_db,
            source_path=resolved_path,
            resolver=lambda: _resolve_marker_db_path(resolved_path),
            columns=MARKER_DB_COLUMNS,
        )
    except FileNotFoundError as exc:  # pragma: no cover - defensive guard
        raise MarkerDatabaseError(str(exc)) from exc

    request = _default_request_id(request_id)
    options = batch_options or BatchOptions()
    chunk_size = options.normalized(len(clusters_payload))
    cache_hits = 0
    llm_batches = 0
    results: dict[str, dict[str, Any]] = {}
    cache = annotation_cache

    logger.info(
        "scanpy.annotate.start",
        extra={
            "request_id": request,
            "clusters": len(clusters_payload),
            "chunk_size": chunk_size,
            "offline": offline_mode,
        },
    )

    for chunk in _chunk(clusters_payload, chunk_size):
        pending: list[dict[str, Any]] = []
        payloads: list[dict[str, Any]] = []
        for cluster in chunk:
            if cache is None:
                pending.append(cluster)
                continue
            cache_payload = _build_cache_payload(
                cluster,
                dataset_context,
                mode=annotator_instance.llm_mode,
            )
            cached = cache.get(cache_payload)
            if cached:
                cache_hits += 1
                results[str(cluster["cluster_id"])] = dict(cached)
            else:
                pending.append(cluster)
                payloads.append(cache_payload)

        if pending:
            llm_batches += 1
            try:
                batch_result = annotator_instance.annotate_batch(pending, dataset_context)
            except Exception as exc:  # pragma: no cover - delegated to Annotator tests
                logger.error(
                    "scanpy.annotate.batch_failed",
                    extra={"request_id": request, "error": str(exc)},
                )
                raise ScanpyAnnotationError("Annotation batch failed") from exc
            for index, cluster in enumerate(pending):
                cluster_id = str(cluster["cluster_id"])
                annotation = dict(batch_result.get(cluster_id) or {})
                results[cluster_id] = annotation
                if cache is not None:
                    cache_payload = payloads[index]
                    cache.set(cache_payload, annotation)

        if progress_callback:
            progress_callback(len(chunk))
        else:
            logger.debug(
                "scanpy.annotate.progress",
                extra={
                    "request_id": request,
                    "processed": len(results),
                    "cache_hits": cache_hits,
                    "llm_batches": llm_batches,
                },
            )

    annotations = _build_annotations(clusters_payload, results)
    settings_override, min_support = _guardrail_settings(guardrails)
    crosschecked = crosscheck_batch(
        annotations,
        marker_df,
        species=dataset_context.get("species"),
        tissue=dataset_context.get("tissue"),
        min_support=min_support,
    )
    report_model: DatasetReport = build_structured_report(
        annotations,
        crosschecked,
        settings_override=settings_override,
    )
    if guardrails and guardrails.as_update():
        logger.info(
            "scanpy.annotate.guardrails.override",
            extra={"request_id": request, "overrides": guardrails.as_update()},
        )

    report = ScanpyDatasetReport(
        dataset=report_model,
        cache_hits=cache_hits,
        llm_batches=llm_batches,
        offline_mode=offline_mode,
        guardrail_min_overlap=min_support,
        guardrail_overrides=guardrails.as_update() if guardrails else {},
        request_id=request,
        chunk_size=chunk_size,
    )

    logger.info(
        "scanpy.annotate.complete",
        extra={
            "request_id": request,
            "summary": report_model.summary.model_dump(),
            "cache_hits": cache_hits,
            "llm_batches": llm_batches,
        },
    )
    return annotations, report


def annotate_anndata(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    result_prefix: str = "gptca",
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    compute_rankings: bool = True,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: DiskAnnotationCache | None = None,
    request_id: str | None = None,
    progress_callback: Callable[[int], None] | None = None,
) -> ScanpyAnnotationResult:
    """Annotate clusters in an AnnData object using GPT Cell Annotator."""

    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs.")

    if compute_rankings:
        _ensure_rankings(
            adata,
            cluster_key,
            top_n_markers=top_n_markers,
        )

    cluster_markers = _extract_markers(adata, cluster_key, top_n_markers)
    clusters_payload = [
        {"cluster_id": cluster_id, "markers": markers}
        for cluster_id, markers in cluster_markers.items()
    ]

    annotations, report = _annotate_clusters_payload(
        clusters_payload,
        species=species,
        tissue=tissue,
        marker_db=marker_db,
        marker_db_path=marker_db_path,
        annotator=annotator,
        batch_options=batch_options,
        guardrails=guardrails,
        annotation_cache=annotation_cache,
        request_id=request_id,
        progress_callback=progress_callback,
    )

    _apply_annotations_to_obs(adata, cluster_key, result_prefix, report.dataset)
    return ScanpyAnnotationResult(adata=adata, annotations=annotations, report=report)


def annotate_from_markers(
    cluster_markers: Mapping[str, Sequence[str]],
    *,
    species: str,
    tissue: str | None = None,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: DiskAnnotationCache | None = None,
    request_id: str | None = None,
) -> MarkerAnnotationResult:
    clusters_payload = [
        {"cluster_id": str(cluster_id), "markers": _unique_ordered(markers)}
        for cluster_id, markers in cluster_markers.items()
    ]
    annotations, report = _annotate_clusters_payload(
        clusters_payload,
        species=species,
        tissue=tissue,
        marker_db=marker_db,
        marker_db_path=marker_db_path,
        annotator=annotator,
        batch_options=batch_options,
        guardrails=guardrails,
        annotation_cache=annotation_cache,
        request_id=request_id,
    )
    return MarkerAnnotationResult(
        annotations=annotations,
        report=report,
    )


def annotate_rank_genes(
    rank_genes_result: Mapping[str, Sequence[str]] | Mapping[str, Any],
    *,
    species: str,
    tissue: str | None = None,
    top_n_markers: int = 5,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    batch_options: BatchOptions | None = None,
    guardrails: GuardrailConfig | None = None,
    annotation_cache: DiskAnnotationCache | None = None,
    request_id: str | None = None,
) -> MarkerAnnotationResult:
    cluster_markers = _markers_from_rankings(rank_genes_result, top_n_markers)
    return annotate_from_markers(
        cluster_markers,
        species=species,
        tissue=tissue,
        marker_db=marker_db,
        marker_db_path=marker_db_path,
        annotator=annotator,
        batch_options=batch_options,
        guardrails=guardrails,
        annotation_cache=annotation_cache,
        request_id=request_id,
    )


def _cluster_labels_from_obs(
    adata: AnnData,
    cluster_key: str,
    label_column: str,
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for obs_idx, cluster_id in adata.obs[cluster_key].astype(str).items():
        label = adata.obs.at[obs_idx, label_column]
        if isinstance(label, str) and label.strip():
            mapping.setdefault(cluster_id, label.strip())
    return mapping


def validate_anndata(
    adata: AnnData,
    cluster_key: str,
    *,
    species: str,
    label_column: str,
    ontology_column: str | None = None,
    tissue: str | None = None,
    top_n_markers: int = 5,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    compute_rankings: bool = True,
    guardrails: GuardrailConfig | None = None,
) -> ScanpyDatasetReport:
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs.")
    if label_column not in adata.obs:
        raise KeyError(f"Label column '{label_column}' not found in adata.obs.")
    if ontology_column and ontology_column not in adata.obs:
        raise KeyError(f"Ontology column '{ontology_column}' not found in adata.obs.")

    if compute_rankings:
        _ensure_rankings(adata, cluster_key, top_n_markers=top_n_markers)

    cluster_markers = _extract_markers(adata, cluster_key, top_n_markers)
    label_mapping = _cluster_labels_from_obs(adata, cluster_key, label_column)
    ontology_mapping: dict[str, str] = {}
    if ontology_column is not None:
        ontology_mapping = _cluster_labels_from_obs(adata, cluster_key, ontology_column)

    annotations: list[dict[str, Any]] = []
    for cluster_id, markers in cluster_markers.items():
        annotation: dict[str, Any] = {
            "cluster_id": cluster_id,
            "markers": markers,
            "primary_label": label_mapping.get(cluster_id, "Unknown or Novel"),
        }
        if ontology_column is not None:
            annotation["ontology_id"] = ontology_mapping.get(cluster_id)
        annotations.append(annotation)

    marker_cache = _marker_cache()
    resolved_path = Path(marker_db_path).expanduser() if marker_db_path is not None else None
    try:
        marker_df = marker_cache.load(
            frame=marker_db,
            source_path=resolved_path,
            resolver=lambda: _resolve_marker_db_path(resolved_path),
            columns=MARKER_DB_COLUMNS,
        )
    except FileNotFoundError as exc:  # pragma: no cover - defensive guard
        raise MarkerDatabaseError(str(exc)) from exc
    settings_override, min_support = _guardrail_settings(guardrails)
    crosschecked = crosscheck_batch(
        annotations,
        marker_df,
        species=species,
        tissue=tissue,
        min_support=min_support,
    )
    report = build_structured_report(
        annotations,
        crosschecked,
        settings_override=settings_override,
    )
    logger.info(
        "scanpy.validate.complete",
        extra={
            "clusters": len(annotations),
            "flagged": report.summary.flagged_clusters,
            "species": species,
            "tissue": tissue,
        },
    )
    return ScanpyDatasetReport(
        dataset=report,
        cache_hits=0,
        llm_batches=0,
        offline_mode=True,
        guardrail_min_overlap=min_support,
        guardrail_overrides=guardrails.as_update() if guardrails else {},
        request_id=_default_request_id(None),
        chunk_size=len(annotations) or 1,
    )


def _read_adata(path: Path) -> AnnData:
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return ad.read_h5ad(path)
    if suffix == ".loom":
        if sc is None:
            raise ImportError("scanpy is required to read Loom files.")
        return sc.read_loom(path, obsm_keys=None)
    raise ValueError(f"Unsupported AnnData format: {path.suffix}")


def _write_report(report: ScanpyDatasetReport, path: Path) -> None:
    payload = report.model_dump()
    payload["dataset"] = report.dataset.model_dump()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _apply_preset(args: argparse.Namespace) -> None:
    if not args.preset:
        return
    preset = SPECIES_PRESETS.get(args.preset)
    if not preset:
        available = ", ".join(SPECIES_PRESETS)
        raise SystemExit(f"Unknown preset '{args.preset}'. Available presets: {available}")
    if not args.species:
        args.species = preset.get("species")
    if not args.tissue and preset.get("tissue"):
        args.tissue = preset["tissue"]


def _annotation_cache_from_args(args: argparse.Namespace) -> DiskAnnotationCache | None:
    if args.cache_dir:
        return DiskAnnotationCache(Path(args.cache_dir).expanduser())
    if args.use_cache:
        cache_root = default_cache_dir() / "annotations"
        return DiskAnnotationCache(cache_root)
    return None


def cmd_annotate(args: argparse.Namespace) -> int:
    _apply_preset(args)
    adata = _read_adata(args.input)

    if args.validate_only:
        report = validate_anndata(
            adata,
            args.cluster_key,
            species=args.species,
            tissue=args.tissue,
            label_column=args.label_column,
            ontology_column=args.ontology_column,
            top_n_markers=args.top_n_markers,
            marker_db_path=args.marker_db,
            compute_rankings=not args.skip_recompute_markers,
            guardrails=_guardrail_config_from_args(
                args.guardrail_min_overlap,
                args.guardrail_force_unknown,
            ),
        )
        if args.json_report:
            _write_report(report, args.json_report)
        if args.summary_csv:
            df = report_to_dataframe(report)
            args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.summary_csv, index=False)
        return 0

    cache = _annotation_cache_from_args(args)
    annotator = Annotator(force_mock=args.offline) if args.offline else Annotator()
    progress_bar = None

    def _progress_update(delta: int) -> None:
        if progress_bar is not None:
            progress_bar.update(delta)

    if tqdm is not None and args.progress and os.isatty(2):
        progress_bar = tqdm(
            total=adata.obs[args.cluster_key].nunique(),
            desc="Annotating",
            unit="clusters",
        )

    result = annotate_anndata(
        adata,
        args.cluster_key,
        species=args.species,
        tissue=args.tissue,
        top_n_markers=args.top_n_markers,
        result_prefix=args.prefix,
        marker_db_path=args.marker_db,
        annotator=annotator,
        compute_rankings=not args.skip_recompute_markers,
        batch_options=BatchOptions(chunk_size=args.chunk_size),
        guardrails=_guardrail_config_from_args(
            args.guardrail_min_overlap,
            args.guardrail_force_unknown,
        ),
        annotation_cache=cache,
        request_id=args.request_id,
        progress_callback=_progress_update if progress_bar else None,
    )

    if progress_bar is not None:
        progress_bar.close()

    output_path = args.output or args.input
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.adata.write(output_path)

    if args.json_report:
        _write_report(result.report, args.json_report)
    if args.summary_csv:
        df = report_to_dataframe(result.report)
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.summary_csv, index=False)

    return 0


def _guardrail_config_from_args(
    min_overlap: int | None,
    force_unknown: bool | None,
) -> GuardrailConfig | None:
    if min_overlap is None and force_unknown is None:
        return None
    return GuardrailConfig(
        min_marker_overlap=min_overlap,
        force_unknown_on_fail=force_unknown,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m gpt_cell_annotator.scanpy",
        description="Scanpy integration tools for GPT Cell Annotator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate clusters in an AnnData file (h5ad or loom).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    annotate_parser.add_argument("input", type=Path, help="Path to the input AnnData file.")
    annotate_parser.add_argument(
        "--cluster-key",
        required=True,
        help="Column in adata.obs that identifies cluster assignments.",
    )
    annotate_parser.add_argument("--species", help="Species context (overrides preset).")
    annotate_parser.add_argument(
        "--preset",
        choices=sorted(SPECIES_PRESETS),
        help="Use preset species/tissue.",
    )
    annotate_parser.add_argument("--tissue", help="Optional tissue context.")
    annotate_parser.add_argument(
        "--top-n-markers",
        type=int,
        default=5,
        help="Number of top markers per cluster to send for annotation.",
    )
    annotate_parser.add_argument(
        "--chunk-size",
        type=int,
        default=32,
        help="Number of clusters to annotate per LLM batch.",
    )
    annotate_parser.add_argument(
        "--prefix",
        default="gptca",
        help="Prefix for columns written to adata.obs.",
    )
    annotate_parser.add_argument(
        "--output",
        type=Path,
        help="Destination AnnData path. Defaults to overwriting input.",
    )
    annotate_parser.add_argument(
        "--json-report",
        type=Path,
        help="Optional path to write the full annotation report JSON.",
    )
    annotate_parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write a per-cluster summary CSV.",
    )
    annotate_parser.add_argument(
        "--marker-db",
        type=Path,
        help="Override path to marker_db.parquet.",
    )
    annotate_parser.add_argument(
        "--skip-recompute-markers",
        action="store_true",
        help="Assume rank_genes_groups already computed; do not recompute.",
    )
    annotate_parser.add_argument(
        "--offline",
        action="store_true",
        help="Force heuristic mock mode and avoid network access.",
    )
    annotate_parser.add_argument(
        "--guardrail-min-overlap",
        type=int,
        help="Override validation_min_marker_overlap for this run.",
    )
    annotate_parser.add_argument(
        "--guardrail-force-unknown",
        dest="guardrail_force_unknown",
        action="store_true",
        default=None,
        help="Force downgrade to Unknown when guardrails fail.",
    )
    annotate_parser.add_argument(
        "--no-guardrail-force-unknown",
        dest="guardrail_force_unknown",
        action="store_false",
        help="Disable downgrade to Unknown when guardrails fail.",
    )
    annotate_parser.add_argument(
        "--request-id",
        help="Attach a custom request identifier for telemetry.",
    )
    annotate_parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Directory for annotation cache artifacts.",
    )
    annotate_parser.add_argument(
        "--use-cache",
        action="store_true",
        help="Enable persistent disk cache for annotation responses.",
    )
    annotate_parser.add_argument(
        "--progress",
        action="store_true",
        help="Display tqdm progress bar when running in an interactive terminal.",
    )
    annotate_parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Skip LLM annotations and only run guardrail validation.",
    )
    annotate_parser.add_argument(
        "--label-column",
        default="gptca_label",
        help="Existing annotation column to validate when using --validate-only.",
    )
    annotate_parser.add_argument(
        "--ontology-column",
        help="Optional ontology column to validate when using --validate-only.",
    )
    annotate_parser.set_defaults(func=cmd_annotate)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    command = cast(Callable[[argparse.Namespace], int], args.func)
    return command(args)


__all__ = [
    "BatchOptions",
    "DiskAnnotationCache",
    "GuardrailConfig",
    "MarkerAnnotationResult",
    "ScanpyAnnotationResult",
    "ScanpyDatasetReport",
    "annotate_anndata",
    "annotate_from_markers",
    "annotate_rank_genes",
    "main",
    "report_to_dataframe",
    "validate_anndata",
]
