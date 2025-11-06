# Releasing GPT Cell Annotator

This checklist streamlines publishing wheels, source distributions, and container images. It assumes you have push access to PyPI, GitHub Container Registry (GHCR), and the project repository.

## Python Package Release Checklist
Follow this mini checklist whenever you cut a new Python release:
- [ ] Build fresh artifacts with `uv build` (mirrors the tox packaging environment).
- [ ] Run `tox -e packaging` or `make packaging-check` to rebuild, run `twine check`, and smoke-test the `gca scanpy` CLI from the generated wheel.
- [ ] Capture SHA256 checksums: `shasum -a 256 dist/*`.
- [ ] In a throwaway virtual environment, install the wheel with Scanpy extras and run the CLI:
  ```bash
  python -m venv /tmp/gca-scanpy-smoke
  source /tmp/gca-scanpy-smoke/bin/activate
  pip install dist/gpt_cell_annotator-*.whl[scanpy]
  gca scanpy annotate --help
  deactivate && rm -rf /tmp/gca-scanpy-smoke
  ```
- [ ] Log the compatibility tier you validated (AnnData/Scanpy versions) in the release notes.

## 1. Pre-flight
- [ ] Ensure your working tree is clean (`git status`).
- [ ] Update `pyproject.toml` version and run `poetry lock` if dependencies changed.
- [ ] Regenerate assets or marker databases if needed.
- [x] Run the full validation suite:
  ```bash
  make install          # installs api/ui/scanpy/dev extras
  make lint
  make test
  PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest tests/test_scanpy_integration.py -k cli
  ```

## 2. Build & Verify Artifacts
- [x] Create fresh build outputs:
  ```bash
  make release          # poetry build + twine check
  ```
- [x] Inspect the contents of `dist/` to confirm both `*.tar.gz` and `*.whl` exist.
- [x] Smoke-test the wheel in a clean virtual environment:
  ```bash
  python -m venv /tmp/gca-release
  source /tmp/gca-release/bin/activate
  pip install dist/gpt_cell_annotator-*.whl
  gca scanpy annotate --help
  deactivate && rm -rf /tmp/gca-release
  ```

## 3. Publish to PyPI (or TestPyPI)
- [ ] Authenticate (`twine` reads `~/.pypirc`):
  ```bash
  twine upload dist/*               # or use --repository testpypi for staging
  ```
- [ ] Verify installation on a clean machine:
  ```bash
  pip install "gpt-cell-annotator[api,scanpy]"
  gca --version
  ```
- [x] Archive the installation matrix and smoke-test notes under `docs/reports/releases/${VERSION}/install_matrix.md`.

## 4. Build & Push Docker Images
- [ ] Build the release image:
  ```bash
  docker build -t gpt-cell-annotator:${VERSION} .
  ```
- [ ] Tag & push to GHCR:
  ```bash
  docker tag gpt-cell-annotator:${VERSION} ghcr.io/jameshyojaelee/gpt-cell-annotator:${VERSION}
  docker push ghcr.io/jameshyojaelee/gpt-cell-annotator:${VERSION}
  ```
- [ ] Record image digest (`docker inspect --format='{{index .RepoDigests 0}}' ...`) and add to the release notes.

## 5. Git Tag & Release Notes
- [ ] Commit the version bump and changelog updates.
- [ ] Tag the release: `git tag -a v${VERSION} -m "Release v${VERSION}"`
- [ ] Push commits and tags: `git push && git push origin v${VERSION}`
- [ ] Draft GitHub Release notes including:
  - PyPI package link & SHA256 checksums (`shasum -a 256 dist/*`)
  - Docker image digest
  - Highlights (new CLI flags, API changes, docs)

## 6. Post-release
- [ ] Notify downstream integrators (Slack/email) with upgrade instructions.
- [ ] Close the milestone and open a new one for the next cycle.
- [ ] Reset development version (e.g., bump to `0.x.y-dev`).

Happy releasing!

_Release 0.2.0 provenance: see `docs/reports/releases/0.2.0/` for smoke tests, install matrix, hashes, and release notes._
