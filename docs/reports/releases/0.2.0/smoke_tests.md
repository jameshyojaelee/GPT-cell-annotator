# Smoke Tests for v0.2.0 Release

Run log captured on 2025-11-04T15:15:47-05:00.

## 1. Pre-flight Validation

### make install
```console
$ make install
poetry install --extras "dev" --extras "scanpy" --extras "api" --extras "ui"
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Installing dependencies from lock file

No dependencies to install or update

Installing the current project: gpt-cell-annotator (0.2.0)
```

### make lint
```console
$ make lint
poetry run ruff check backend frontend evaluation scripts config tests
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
All checks passed!
poetry run mypy backend
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Success: no issues found in 17 source files
```

### make test
```console
$ make test
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
============================= test session starts ==============================
platform linux -- Python 3.11.14, pytest-8.4.2, pluggy-1.6.0
rootdir: /gpfs/commons/home/jameslee/GPT-cell-annotator
configfile: pyproject.toml
testpaths: tests
asyncio: mode=Mode.STRICT
collected 65 items

tests/test_annotator.py .......                                          [ 10%]
tests/test_api.py .......                                                [ 21%]
tests/test_benchmark_runner.py .....                                     [ 29%]
tests/test_cache.py ...                                                  [ 33%]
tests/test_cli.py ..                                                     [ 36%]
tests/test_frontend.py .......                                           [ 47%]
tests/test_gene_normalization.py ..                                      [ 50%]
tests/test_marker_loader.py ........                                     [ 63%]
tests/test_plots.py .                                                    [ 64%]
tests/test_prompts.py ..                                                 [ 67%]
tests/test_retrieval.py ..                                               [ 70%]
tests/test_scanpy_integration.py ........                                [ 83%]
tests/test_single_cluster.py ....                                        [ 89%]
tests/test_validation.py .......                                         [100%]

=============================== warnings summary ===============================
config/settings.py:11
  /gpfs/commons/home/jameslee/GPT-cell-annotator/config/settings.py:11: PydanticDeprecatedSince20: Support for class-based `config` is deprecated, use ConfigDict instead. Deprecated in Pydantic V2.0 to be removed in V3.0. See Pydantic V2 Migration Guide at https://errors.pydantic.dev/2.12/migration/
    class Settings(BaseSettings):

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================== 65 passed, 1 warning in 35.34s ========================
```

### PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest tests/test_scanpy_integration.py -k cli
```console
$ PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 poetry run pytest tests/test_scanpy_integration.py -k cli
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
============================= test session starts ==============================
platform linux -- Python 3.11.14, pytest-8.4.2, pluggy-1.6.0
rootdir: /gpfs/commons/home/jameslee/GPT-cell-annotator
configfile: pyproject.toml
asyncio: mode=Mode.STRICT
collected 8 items / 6 deselected / 2 selected

tests/test_scanpy_integration.py ..                                      [100%]

=============================== warnings summary ===============================
config/settings.py:11
  /gpfs/commons/home/jameslee/GPT-cell-annotator/config/settings.py:11: PydanticDeprecatedSince20: Support for class-based `config` is deprecated, use ConfigDict instead. Deprecated in Pydantic V2.0 to be removed in V3.0. See Pydantic V2 Migration Guide at https://errors.pydantic.dev/2.12/migration/
    class Settings(BaseSettings):

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
================= 2 passed, 6 deselected, 1 warning in 22.74s ==================
```

## 2. Build & Distribution Checks

### poetry build
```console
$ poetry build
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Building gpt-cell-annotator (0.2.0)
Building sdist
  - Building sdist
  - Built gpt_cell_annotator-0.2.0.tar.gz
Building wheel
  - Building wheel
  - Built gpt_cell_annotator-0.2.0-py3-none-any.whl
```

### poetry run twine check dist/*
```console
$ poetry run twine check dist/*
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Checking dist/gpt_cell_annotator-0.2.0-py3-none-any.whl: PASSED
Checking dist/gpt_cell_annotator-0.2.0.tar.gz: PASSED
```

## 3. Extras Smoke Matrix

### Base extra (`[]`)
```console
$ source /tmp/gca-base-0.2.0/bin/activate && gca --version
0.2.0
$ gca annotate data/demo/pbmc_markers.csv --offline --out-json /tmp/base.json
...rich table output (Unknown or Novel for demo markers)...
$ deactivate
```

### API extra (`[api]`)
```console
$ source /tmp/gca-api-0.2.0/bin/activate
$ gca api --offline --port 8081 &
[server] INFO: Uvicorn running on http://0.0.0.0:8081
$ curl -sf http://127.0.0.1:8081/health
{"status":"ok","llm_mode":"mock","cache_enabled":false}
$ kill -INT <PID> && wait <PID>
$ deactivate
```

### UI extra (`[ui]`)
```console
$ source /tmp/gca-ui-0.2.0/bin/activate
$ python -c "import streamlit; print(streamlit.__version__)"
1.51.0
$ gca annotate data/demo/pbmc_markers.csv --offline --out-json /tmp/ui.json
...rich table output (Unknown or Novel for demo markers)...
$ deactivate
```

### Scanpy extra (`[scanpy]`)
```console
$ source /tmp/gca-scanpy-0.2.0/bin/activate
$ python - <<'PY'
import numpy as np
import pandas as pd
import anndata as ad
obs = pd.DataFrame({"leiden": pd.Categorical(["0","0","1","1","2","2"]),
                    "cell_id": [f"cell{i}" for i in range(6)]})
var = pd.DataFrame(index=["MS4A1","LYZ","GNLY"])
adata = ad.AnnData(X=np.random.rand(6, 3), obs=obs, var=var)
adata.write_h5ad("/tmp/pbmc_demo.h5ad")
PY
$ gca scanpy annotate /tmp/pbmc_demo.h5ad --cluster-key leiden --species "Homo sapiens" --offline --summary-json /tmp/scanpy.json
2025-11-04 15:22:00 [info] scanpy.annotate.complete summary={'total_clusters': 3, 'supported_clusters': 0, 'flagged_clusters': 3, 'unknown_clusters': ['0', '1', '2']}
$ deactivate
```

### Full extra (`[full]`)
```console
$ source /tmp/gca-full-0.2.0/bin/activate
$ export GPT_CELL_ANNOTATOR_HOME=/tmp/gca-full-home
$ gca annotate data/demo/pbmc_markers.csv --offline --out-json $GPT_CELL_ANNOTATOR_HOME/base.json
...rich table output (Unknown or Novel for demo markers)...
$ gca scanpy annotate /tmp/pbmc_demo.h5ad --cluster-key leiden --species "Homo sapiens" --offline --summary-json $GPT_CELL_ANNOTATOR_HOME/scanpy.json
2025-11-04 15:23:25 [info] scanpy.annotate.complete summary={'total_clusters': 3, 'supported_clusters': 0, 'flagged_clusters': 3, 'unknown_clusters': ['0', '1', '2']}
$ ls -R $GPT_CELL_ANNOTATOR_HOME
base.json
config/marker_sources.yaml
data/demo/*.csv
data/processed/marker_db.{parquet,sqlite}
$ deactivate
```

## 4. Hashing Transcript

```console
$ shasum -a 256 dist/* > docs/reports/releases/0.2.0/hashes.txt
$ cat docs/reports/releases/0.2.0/hashes.txt
5e2ff7eb619a824e6d07d64e2a8886497b6fb85a8de54f03b10b0e1e01aef34f  dist/gpt_cell_annotator-0.2.0-py3-none-any.whl
0f2f3e5866f6a23f44db637117124c8287ac0456370a1e8c491932d44843f12d  dist/gpt_cell_annotator-0.2.0.tar.gz
```
