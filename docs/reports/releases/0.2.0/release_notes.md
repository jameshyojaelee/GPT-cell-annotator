# GPT Cell Annotator 0.2.0 Release Notes

## Highlights
- Scanpy tooling now supports asynchronous batch annotation, offline telemetry, and richer summaries (`gca scanpy`). Full walkthroughs live in the [Scanpy integration guide](../../../scanpy_integration.md).
- Disk-backed caching, guardrail telemetry, and CLI/API parity make offline annotation repeatable; cache layout is demonstrated in the `[full]` smoke tests.
- Seurat users gain a companion R package (`clients/r/gptcellannotator`) with pkgdown docs—see the [Seurat integration playbook](../../../seurat_integration.md).
- Documentation refresh covers install extras, offline workflows, and troubleshooting matrices; start with the [Install guide](../../../install.md).

## Installation
- PyPI wheel/sdist rebuilt via `poetry build`; `poetry run twine check dist/*` passes.
- Tested extras: `[]`, `[api]`, `[ui]`, `[scanpy]`, `[full]` across isolated Python 3.11 virtual environments. Command matrix and observations live in [`install_matrix.md`](install_matrix.md).
- Offline demos rely on bundled assets under `gpt_cell_annotator/_assets`. For Scanpy smoke tests, generate a small `.h5ad` fixture (script snippet in [`smoke_tests.md`](smoke_tests.md)).

## QA Matrix
- Pre-flight (`make install`, `make lint`, `make test`) logs captured in [`smoke_tests.md`](smoke_tests.md).
- API health check returns `{"status":"ok","llm_mode":"mock","cache_enabled":false}`.
- `[full]` install verifies cache bootstrap under a disposable `GPT_CELL_ANNOTATOR_HOME`.
- Additional recommendations: ship a tiny `.h5ad` demo or document fixture creation to streamline Scanpy validation.

## Hashes
- SHA-256 checksums recorded in [`hashes.txt`](hashes.txt):
  - `dist/gpt_cell_annotator-0.2.0-py3-none-any.whl` → `5e2ff7eb619a824e6d07d64e2a8886497b6fb85a8de54f03b10b0e1e01aef34f`
  - `dist/gpt_cell_annotator-0.2.0.tar.gz` → `0f2f3e5866f6a23f44db637117124c8287ac0456370a1e8c491932d44843f12d`

## Next Steps
1. **PyPI** – Obtain credentials and run `twine upload dist/*` (TODO: not executed).
2. **GitHub** – Draft/tag `v0.2.0`, publish release notes, and attach `dist/` artifacts.
3. **Docker** – Build/push GHCR images per `RELEASING.md` (record digest post-push).
4. **Comms** – Notify downstream users, update milestones, and cue 0.3.0 planning once publication completes.
