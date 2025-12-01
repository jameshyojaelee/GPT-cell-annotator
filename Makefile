.PHONY: install format lint test typecheck ruff ruff-fix build-marker-db clean pre-commit packaging-check

install:
	poetry install --extras "dev" --extras "scanpy"

format:
	poetry run black src scripts tests

ruff:
	poetry run ruff check src scripts tests

ruff-fix:
	poetry run ruff check --fix src scripts tests

lint: ruff typecheck

typecheck:
	poetry run mypy src/backend src/gpt_cell_annotator src/config

pre-commit:
	poetry run pre-commit run --all-files

test: lint
	PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest

release:
	poetry build
	poetry run twine check dist/*

packaging-check:
	tox -e packaging

build-marker-db:
	poetry run python scripts/build_marker_db.py

clean:
	find src tests -name "__pycache__" -type d -prune -exec rm -rf {} +
	rm -rf __pycache__ .pytest_cache .ruff_cache .mypy_cache dist gptcellannotator.Rcheck
	rm -f marker_db.sqlite
