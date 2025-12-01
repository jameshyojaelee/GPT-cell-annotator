"""Caching utilities for GPT Cell Annotator."""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path

import pandas as pd

CACHE_DIR_ENV_VAR = "GCA_CACHE_DIR"
_DEFAULT_CACHE_SUBDIR = "gpt-cell-annotator"


def default_cache_dir() -> Path:
    """Return the default on-disk cache directory, respecting environment overrides."""

    override = os.environ.get(CACHE_DIR_ENV_VAR)
    if override:
        return Path(override).expanduser()
    return Path.home() / ".cache" / _DEFAULT_CACHE_SUBDIR


def _ensure_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        for column in missing:
            frame[column] = None
    return frame[list(columns)].copy()


@dataclass(slots=True)
class MarkerDatabaseCache:
    """Materialise marker databases in memory/disk for fast reuse."""

    cache_dir: Path = field(default_factory=default_cache_dir)
    persist: bool = True
    _memory: dict[str, tuple[float, pd.DataFrame]] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        self.cache_dir = self.cache_dir.expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        (self.cache_dir / "marker-db").mkdir(parents=True, exist_ok=True)

    def _key(self, source: Path) -> str:
        canonical = str(source.expanduser().resolve())
        return sha256(canonical.encode("utf-8")).hexdigest()

    def _disk_path(self, key: str) -> Path:
        return self.cache_dir / "marker-db" / f"{key}.parquet"

    def _meta_path(self, key: str) -> Path:
        return self.cache_dir / "marker-db" / f"{key}.json"

    def load_frame(self, source: Path, *, columns: Sequence[str]) -> pd.DataFrame:
        """Return a dataframe sourced from `source`, cached in memory/disk when possible."""

        source = source.expanduser()
        key = self._key(source)
        try:
            mtime = source.stat().st_mtime
        except FileNotFoundError:
            mtime = None

        cached = self._memory.get(key)
        if cached and mtime is not None and cached[0] == mtime:
            return cached[1].copy()

        try:
            frame = pd.read_parquet(source)
            trimmed = _ensure_columns(frame, columns)
            if mtime is not None:
                self._memory[key] = (mtime, trimmed)
                if self.persist:
                    self._write_disk_cache(key, trimmed, mtime)
            else:
                self._memory[key] = (-1.0, trimmed)
            return trimmed.copy()
        except FileNotFoundError:
            if not self.persist:
                raise
            disk_frame = self._read_disk_cache(key, columns=columns)
            if disk_frame is not None:
                self._memory[key] = (-1.0, disk_frame)
                return disk_frame.copy()
            raise

    def load(
        self,
        *,
        frame: pd.DataFrame | None = None,
        source_path: Path | None = None,
        resolver: Callable[[], Path] | None = None,
        columns: Sequence[str],
    ) -> pd.DataFrame:
        """Return a dataframe from an in-memory frame or resolved path.

        Args:
            frame: Optional pre-loaded dataframe. When provided, a copy is returned.
            source_path: Candidate on-disk path for the marker DB.
            resolver: Callable returning a fallback path when `source_path` is None.
            columns: Expected columns, missing ones will be added as null.
        """

        if frame is not None:
            return _ensure_columns(frame.copy(), columns)
        if source_path is not None:
            return self.load_frame(source_path, columns=columns)
        if resolver is None:
            raise ValueError("resolver must be provided when frame and source_path are None.")
        resolved = resolver()
        return self.load_frame(resolved, columns=columns)

    def _write_disk_cache(self, key: str, frame: pd.DataFrame, mtime: float) -> None:
        cache_path = self._disk_path(key)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(cache_path, index=False)
        meta_path = self._meta_path(key)
        meta_path.write_text(json.dumps({"mtime": mtime}), encoding="utf-8")

    def _read_disk_cache(self, key: str, *, columns: Sequence[str]) -> pd.DataFrame | None:
        cache_path = self._disk_path(key)
        if not cache_path.exists():
            return None
        frame = pd.read_parquet(cache_path)
        return _ensure_columns(frame, columns)
