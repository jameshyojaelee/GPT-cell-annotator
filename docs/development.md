# Development Workflow

This guide captures the standard local workflow for GPT Cell Annotator contributors.

## Environment Setup
- Install Poetry (≥1.8) and ensure Python 3.11 is available on your PATH.
- Bootstrap project dependencies (includes API/UI/Scanpy extras and dev tooling):
  ```bash
  poetry install --no-interaction --extras "dev" --extras "scanpy" --extras "api" --extras "ui"
  ```
- Alternatively, run `make install` which wraps the command above.

## Pre-commit Hooks
- Install hooks into your local Git clone:
  ```bash
  poetry run pre-commit install
  ```
- To run all hooks on demand (same command CI executes), use either:
  ```bash
  poetry run pre-commit run --all-files
  ```
  or `make pre-commit`.

## Test Workflow
- `make test` now runs Ruff + mypy (via `lint`) before executing pytest with `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1`.
- Use `make test` prior to pushing PRs to mirror the default CI gate.

## Focused Checks
- Ruff only: `poetry run ruff check backend frontend evaluation scripts config tests`
- Mypy only: `poetry run mypy backend`
- Pytest only: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest`

## Verification
- `poetry run pre-commit run --all-files` — passed (Python 3.11 via Poetry, 2025-11-04).
- `make test` — lint + pytest pipeline passed (65 tests).
- `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest` — 65 passed, 0 failed (32.83s).
- `poetry run ruff check …` — passed.
- `poetry run mypy backend` — passed (17 source files).
