from __future__ import annotations

from config import settings


def test_data_dir_honours_gca_env(monkeypatch):
    monkeypatch.setenv("GPT_CELL_ANNOTATOR_DATA_DIR", "/tmp/custom-db")
    settings.get_settings.cache_clear()
    cfg = settings.get_settings()
    assert cfg.data_dir == "/tmp/custom-db"
    settings.get_settings.cache_clear()
