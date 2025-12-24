"""High-level marker annotation helpers for GPT Cell Annotator."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

from backend.llm.annotator import Annotator
from backend.validation.crosscheck import crosscheck_batch
from backend.validation.report import DatasetReport, build_structured_report
from config.settings import Settings, get_settings
from gpt_cell_annotator import assets

MARKER_DB_ENV_VAR = "GCA_MARKER_DB_PATH"


@dataclass(slots=True)
class MarkerAnnotationResult:
    """Result returned when annotating marker lists."""

    annotations: list[dict[str, Any]]
    report: DatasetReport

    def report_dict(self) -> dict[str, Any]:
        return self.report.model_dump()


def report_to_dataframe(report: DatasetReport | dict[str, Any]) -> pd.DataFrame:
    """Return a dataframe summarising cluster annotations."""

    dataset = report if isinstance(report, DatasetReport) else DatasetReport.model_validate(report)
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


def _unique_ordered(markers: Sequence[Any]) -> list[str]:
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


def _resolve_marker_db_path(marker_db_path: str | Path | None) -> Path:
    if marker_db_path is not None:
        candidate = Path(marker_db_path).expanduser()
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

    raise FileNotFoundError("marker_db.parquet could not be located or materialised.")


def _load_marker_db(
    marker_db: pd.DataFrame | None,
    marker_db_path: str | Path | None,
) -> pd.DataFrame:
    if marker_db is not None:
        return marker_db.copy()
    path = _resolve_marker_db_path(marker_db_path)
    return pd.read_parquet(path)


def _guardrail_settings(
    settings: Settings,
    *,
    min_marker_overlap: int | None,
    force_unknown_on_fail: bool | None,
) -> tuple[Settings | None, int]:
    if min_marker_overlap is None and force_unknown_on_fail is None:
        return None, max(1, settings.validation_min_marker_overlap)

    update: dict[str, Any] = {}
    if min_marker_overlap is not None:
        update["validation_min_marker_overlap"] = max(0, min_marker_overlap)
    if force_unknown_on_fail is not None:
        update["validation_force_unknown_on_fail"] = bool(force_unknown_on_fail)

    override = settings.model_copy(update=update)
    min_support = max(1, override.validation_min_marker_overlap)
    return override, min_support


def annotate_markers(
    cluster_markers: Mapping[str, Sequence[str]],
    *,
    species: str,
    tissue: str | None = None,
    marker_db: pd.DataFrame | None = None,
    marker_db_path: str | Path | None = None,
    annotator: Annotator | None = None,
    settings: Settings | None = None,
    min_marker_overlap: int | None = None,
    force_unknown_on_fail: bool | None = None,
) -> MarkerAnnotationResult:
    """Annotate marker lists for clusters and return annotations + validation report."""

    if not cluster_markers:
        raise ValueError("cluster_markers is empty")

    settings = settings or get_settings()
    settings_override, min_support = _guardrail_settings(
        settings,
        min_marker_overlap=min_marker_overlap,
        force_unknown_on_fail=force_unknown_on_fail,
    )

    dataset_context: dict[str, str] = {"species": species}
    if tissue:
        dataset_context["tissue"] = tissue

    clusters_payload = [
        {"cluster_id": str(cluster_id), "markers": _unique_ordered(markers)}
        for cluster_id, markers in cluster_markers.items()
    ]

    annotator_instance = annotator or Annotator(settings=settings)
    raw_result = annotator_instance.annotate_batch(clusters_payload, dataset_context)

    annotations: list[dict[str, Any]] = []
    for cluster in clusters_payload:
        cid = str(cluster["cluster_id"])
        cluster_result = dict(raw_result.get(cid) or {})
        cluster_result.setdefault("primary_label", "Unknown or Novel")
        cluster_result.setdefault("confidence", "Unknown")
        cluster_result.setdefault("rationale", "")
        annotations.append(
            {
                **cluster_result,
                "cluster_id": cid,
                "markers": cluster.get("markers", []),
            }
        )

    marker_df = _load_marker_db(marker_db, marker_db_path)
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

    return MarkerAnnotationResult(annotations=annotations, report=report)


__all__ = ["MarkerAnnotationResult", "annotate_markers", "report_to_dataframe"]
