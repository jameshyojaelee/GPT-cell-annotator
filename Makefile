.PHONY: install format lint test typecheck ruff ruff-fix build-marker-db clean pre-commit packaging-check

install:
	poetry install --extras "dev" --extras "scanpy"

format:
	poetry run black backend scripts config tests gpt_cell_annotator

ruff:
	poetry run ruff check backend scripts config tests gpt_cell_annotator

ruff-fix:
	poetry run ruff check --fix backend scripts config tests gpt_cell_annotator

lint: ruff typecheck

typecheck:
	poetry run mypy backend

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
	rm -rf __pycache__ .pytest_cache .ruff_cache .mypy_cache
	rm -f marker_db.sqlite
