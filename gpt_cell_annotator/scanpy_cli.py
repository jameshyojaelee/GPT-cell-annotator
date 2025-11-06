"""Entry point alias for the Scanpy-focused CLI."""

from __future__ import annotations

from typing import Sequence

from gpt_cell_annotator.scanpy import main as _scanpy_main


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate to the Scanpy module's main entry point."""

    return _scanpy_main(argv)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
