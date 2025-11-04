# GPT Cell Annotator 0.2.0 – Installation Matrix QA

Date: 2025-11-04  
Engineer: GPT-5 Codex (release engineer)

## Matrix

| Scenario | Command matrix | Exit status | Notes |
| --- | --- | --- | --- |
| Base (`[]`) | `${POETRY_PY311} -m venv /tmp/gca-base-0.2.0`<br>`pip install dist/gpt_cell_annotator-0.2.0-py3-none-any.whl`<br>`gca --version`<br>`gca annotate data/demo/pbmc_markers.csv --offline --out-json /tmp/base.json` | ✅ 0 | CLI reports `0.2.0`; offline annotate writes `/tmp/base.json` with flagged demo markers. |
| API (`[api]`) | `${POETRY_PY311} -m venv /tmp/gca-api-0.2.0`<br>`pip install 'dist/gpt_cell_annotator-0.2.0-py3-none-any.whl[api]'`<br>`gca api --offline --port 8081 &` → `curl -sf http://127.0.0.1:8081/health` → SIGINT | ✅ 0 | Health endpoint returns `{"status":"ok","llm_mode":"mock","cache_enabled":false}`; server shuts down cleanly via `SIGINT`. |
| UI (`[ui]`) | `${POETRY_PY311} -m venv /tmp/gca-ui-0.2.0`<br>`pip install 'dist/gpt_cell_annotator-0.2.0-py3-none-any.whl[ui]'`<br>`python -c "import streamlit"`<br>`gca annotate data/demo/pbmc_markers.csv --offline --out-json /tmp/ui.json` | ✅ 0 | Streamlit imports without warnings; CLI reuse works with UI dependencies present. |
| Scanpy (`[scanpy]`) | `${POETRY_PY311} -m venv /tmp/gca-scanpy-0.2.0`<br>`pip install 'dist/gpt_cell_annotator-0.2.0-py3-none-any.whl[scanpy]'`<br>`gca scanpy annotate /tmp/pbmc_demo.h5ad --cluster-key leiden --species \"Homo sapiens\" --offline --summary-json /tmp/scanpy.json` | ✅ 0 | Repo lacks `.h5ad` demo; generated synthetic AnnData via helper script (categorical `leiden`, ≥2 cells/cluster) before invoking CLI. Summary JSON captured at `/tmp/scanpy.json`. |
| Full (`[full]`) | `${POETRY_PY311} -m venv /tmp/gca-full-0.2.0`<br>`pip install 'dist/gpt_cell_annotator-0.2.0-py3-none-any.whl[full]'`<br>`export GPT_CELL_ANNOTATOR_HOME=/tmp/gca-full-home`<br>`gca annotate data/demo/pbmc_markers.csv --offline --out-json $GPT_CELL_ANNOTATOR_HOME/base.json`<br>`gca scanpy annotate /tmp/pbmc_demo.h5ad --cluster-key leiden --species "Homo sapiens" --offline --summary-json $GPT_CELL_ANNOTATOR_HOME/scanpy.json`<br>`ls -R $GPT_CELL_ANNOTATOR_HOME` | ✅ 0 | New cache home seeds `config/` and `data/` folders plus CLI outputs; both annotate commands succeed with full extra stack. |

## Follow-ups

- Consider adding a lightweight `.h5ad` fixture (or documented generator script) under `data/demo/` to simplify Scanpy smoke tests.
- Pip warns about upgradeable tooling (`pip 25.3`); no action taken during release.
