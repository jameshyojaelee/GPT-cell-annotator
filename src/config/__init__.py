"""Config package for GPT Cell Annotator."""

from .settings import (
    CONFIG_DIR,
    DEFAULT_ORTHOLOG_MAPPING,
    DEFAULT_SYNONYM_CONFIG,
    Settings,
    get_settings,
)

__all__ = [
    "CONFIG_DIR",
    "DEFAULT_ORTHOLOG_MAPPING",
    "DEFAULT_SYNONYM_CONFIG",
    "Settings",
    "get_settings",
]
