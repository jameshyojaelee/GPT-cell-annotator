# Development Workflow

## Environment

```bash
poetry install --extras "scanpy" --extras "dev"
# or via venv:
# python -m venv .venv && source .venv/bin/activate
# pip install -e ".[scanpy,dev]"
```

Ensure `PYTHONPATH` includes `src` when running without an editable install; the pytest config already adds it during test runs.

## Common tasks

- Format: `poetry run black src scripts tests`
- Lint: `poetry run ruff check src scripts tests`
- Type check: `poetry run mypy src/backend src/gpt_cell_annotator src/config`
- Tests: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest`

The `make` targets wrap these commands (`make format`, `make ruff`, `make test`).

## Smoke checks

```bash
gca annotate data/demo/pbmc_markers.csv --offline --out-json /tmp/demo.json
gca build-db --offline --output-dir /tmp/gca-db
```

These commands exercise the bundled assets and mock annotator without network access.
