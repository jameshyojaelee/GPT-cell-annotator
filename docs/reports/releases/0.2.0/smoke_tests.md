# Smoke Tests for v0.2.0 Release

This log captures pre-flight commands executed on 2025-11-03T11:53:36-05:00.

## Commands

### make install
poetry install --extras "dev" --extras "scanpy" --extras "api" --extras "ui"
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Installing dependencies from lock file

Package operations: 0 installs, 1 update, 0 removals

  - Downgrading twine (6.2.0 -> 5.1.1)

Installing the current project: gpt-cell-annotator (0.2.0)

### make lint
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

### make test
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
======================== 65 passed, 1 warning in 36.11s ========================

### poetry build
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

### poetry run twine check dist/*
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Checking dist/gpt_cell_annotator-0.2.0-py3-none-any.whl: [31mERROR   [0m InvalidDistribution: Metadata is missing required fields: Name,        
         Version.                                                               
         Make sure the distribution includes the files where those fields are   
         specified, and is using a supported Metadata-Version: 1.0, 1.1, 1.2,   
         2.0, 2.1, 2.2, 2.3.                                                    

### poetry build (after twine upgrade)
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

### poetry run twine check dist/* (after upgrade)
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Traceback (most recent call last):
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 463, in from_name
    return next(iter(cls._prefer_valid(cls.discover(name=name))))
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StopIteration

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/bin/twine", line 8, in <module>
    sys.exit(main())
             ^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/__main__.py", line 33, in main
    error = cli.dispatch(sys.argv[1:])
            ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 99, in dispatch
    version=f"%(prog)s version {twine.__version__} ({dep_versions()})",
                                                     ^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 86, in dep_versions
    "{}: {}".format(*dependency) for dependency in list_dependencies_and_versions()
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 81, in list_dependencies_and_versions
    return [(dep, importlib_metadata.version(dep)) for dep in deps]
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 81, in <listcomp>
    return [(dep, importlib_metadata.version(dep)) for dep in deps]
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 1070, in version
    return distribution(distribution_name).version
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 1043, in distribution
    return Distribution.from_name(distribution_name)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 465, in from_name
    raise PackageNotFoundError(name)
importlib_metadata.PackageNotFoundError: No package metadata was found for pkginfo

### poetry build (after pkginfo restore)
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

### poetry run twine check dist/* (final)
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Traceback (most recent call last):
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 463, in from_name
    return next(iter(cls._prefer_valid(cls.discover(name=name))))
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StopIteration

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/bin/twine", line 8, in <module>
    sys.exit(main())
             ^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/__main__.py", line 33, in main
    error = cli.dispatch(sys.argv[1:])
            ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 99, in dispatch
    version=f"%(prog)s version {twine.__version__} ({dep_versions()})",
                                                     ^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 86, in dep_versions
    "{}: {}".format(*dependency) for dependency in list_dependencies_and_versions()
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 81, in list_dependencies_and_versions
    return [(dep, importlib_metadata.version(dep)) for dep in deps]
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/twine/cli.py", line 81, in <listcomp>
    return [(dep, importlib_metadata.version(dep)) for dep in deps]
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 1070, in version
    return distribution(distribution_name).version
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 1043, in distribution
    return Distribution.from_name(distribution_name)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/gpfs/commons/home/jameslee/.cache/pypoetry/virtualenvs/gpt-cell-annotator-EHIaUXqk-py3.11/lib/python3.11/site-packages/importlib_metadata/__init__.py", line 465, in from_name
    raise PackageNotFoundError(name)
importlib_metadata.PackageNotFoundError: No package metadata was found for pkginfo

### make install (post-dependency updates)
poetry install --extras "dev" --extras "scanpy" --extras "api" --extras "ui"
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Installing dependencies from lock file

Package operations: 1 install, 1 update, 0 removals

  - Installing pkginfo (1.12.1.2)
  - Updating twine (5.1.1 -> 6.2.0)

Installing the current project: gpt-cell-annotator (0.2.0)

### poetry run twine check dist/* (after make install)
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Checking dist/gpt_cell_annotator-0.2.0-py3-none-any.whl: [32mPASSED[0m
Checking dist/gpt_cell_annotator-0.2.0.tar.gz: [32mPASSED[0m

### Base extra validation
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Requirement already satisfied: pip in /scratch/tmp.4noEyhSQqf/base/lib/python3.11/site-packages (24.0)
Collecting pip
  Downloading pip-25.3-py3-none-any.whl.metadata (4.7 kB)
Collecting wheel
  Using cached wheel-0.45.1-py3-none-any.whl.metadata (2.3 kB)
Downloading pip-25.3-py3-none-any.whl (1.8 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.8/1.8 MB 24.3 MB/s eta 0:00:00
Using cached wheel-0.45.1-py3-none-any.whl (72 kB)
Installing collected packages: wheel, pip
  Attempting uninstall: pip
    Found existing installation: pip 24.0
    Uninstalling pip-24.0:
      Successfully uninstalled pip-24.0
Successfully installed pip-25.3 wheel-0.45.1
Processing ./dist/gpt_cell_annotator-0.2.0-py3-none-any.whl
Collecting anndata<0.11.0,>=0.10.7 (from gpt-cell-annotator==0.2.0)
  Using cached anndata-0.10.9-py3-none-any.whl.metadata (6.9 kB)
Collecting httpx<0.28.0,>=0.27.0 (from gpt-cell-annotator==0.2.0)
  Using cached httpx-0.27.2-py3-none-any.whl.metadata (7.1 kB)
Collecting jsonschema<5.0.0,>=4.25.0 (from gpt-cell-annotator==0.2.0)
  Using cached jsonschema-4.25.1-py3-none-any.whl.metadata (7.6 kB)
Collecting numpy<2.0.0,>=1.26.4 (from gpt-cell-annotator==0.2.0)
  Using cached numpy-1.26.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
Collecting openai<2.0.0,>=1.30.0 (from gpt-cell-annotator==0.2.0)
  Using cached openai-1.109.1-py3-none-any.whl.metadata (29 kB)
Collecting pandas<3.0.0,>=2.2.2 (from gpt-cell-annotator==0.2.0)
  Using cached pandas-2.3.3-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl.metadata (91 kB)
Collecting pyarrow<17.0.0,>=16.1.0 (from gpt-cell-annotator==0.2.0)
  Using cached pyarrow-16.1.0-cp311-cp311-manylinux_2_28_x86_64.whl.metadata (3.0 kB)
Collecting pydantic<3.0.0,>=2.7.1 (from gpt-cell-annotator==0.2.0)
  Using cached pydantic-2.12.3-py3-none-any.whl.metadata (87 kB)
Collecting pydantic-settings<3.0.0,>=2.3.4 (from gpt-cell-annotator==0.2.0)
  Using cached pydantic_settings-2.11.0-py3-none-any.whl.metadata (3.4 kB)
Collecting python-dotenv<2.0.0,>=1.0.1 (from gpt-cell-annotator==0.2.0)
  Downloading python_dotenv-1.2.1-py3-none-any.whl.metadata (25 kB)
Collecting pyyaml<7.0.0,>=6.0.2 (from gpt-cell-annotator==0.2.0)
  Using cached pyyaml-6.0.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl.metadata (2.4 kB)
Collecting rich<14.0.0,>=13.7.1 (from gpt-cell-annotator==0.2.0)
  Using cached rich-13.9.4-py3-none-any.whl.metadata (18 kB)
Collecting scikit-learn<2.0.0,>=1.5.0 (from gpt-cell-annotator==0.2.0)
  Using cached scikit_learn-1.7.2-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl.metadata (11 kB)
Collecting sqlalchemy<3.0.0,>=2.0.30 (from gpt-cell-annotator==0.2.0)
  Using cached sqlalchemy-2.0.44-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (9.5 kB)
Collecting structlog<25.0.0,>=24.1.0 (from gpt-cell-annotator==0.2.0)
  Using cached structlog-24.4.0-py3-none-any.whl.metadata (7.3 kB)
Collecting array-api-compat!=1.5,>1.4 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0)
  Using cached array_api_compat-1.12.0-py3-none-any.whl.metadata (2.5 kB)
Collecting h5py>=3.1 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0)
  Using cached h5py-3.15.1-cp311-cp311-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl.metadata (3.0 kB)
Collecting natsort (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0)
  Using cached natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
Collecting packaging>=20.0 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0)
  Using cached packaging-25.0-py3-none-any.whl.metadata (3.3 kB)
Collecting scipy>1.8 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0)
  Downloading scipy-1.16.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl.metadata (62 kB)
Collecting anyio (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached anyio-4.11.0-py3-none-any.whl.metadata (4.1 kB)
Collecting certifi (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached certifi-2025.10.5-py3-none-any.whl.metadata (2.5 kB)
Collecting httpcore==1.* (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached httpcore-1.0.9-py3-none-any.whl.metadata (21 kB)
Collecting idna (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached idna-3.11-py3-none-any.whl.metadata (8.4 kB)
Collecting sniffio (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached sniffio-1.3.1-py3-none-any.whl.metadata (3.9 kB)
Collecting h11>=0.16 (from httpcore==1.*->httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0)
  Using cached h11-0.16.0-py3-none-any.whl.metadata (8.3 kB)
Collecting attrs>=22.2.0 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0)
  Using cached attrs-25.4.0-py3-none-any.whl.metadata (10 kB)
Collecting jsonschema-specifications>=2023.03.6 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0)
  Using cached jsonschema_specifications-2025.9.1-py3-none-any.whl.metadata (2.9 kB)
Collecting referencing>=0.28.4 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0)
  Using cached referencing-0.37.0-py3-none-any.whl.metadata (2.8 kB)
Collecting rpds-py>=0.7.1 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0)
  Downloading rpds_py-0.28.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.1 kB)
Collecting distro<2,>=1.7.0 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0)
  Using cached distro-1.9.0-py3-none-any.whl.metadata (6.8 kB)
Collecting jiter<1,>=0.4.0 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0)
  Using cached jiter-0.11.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (5.2 kB)
Collecting tqdm>4 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0)
  Using cached tqdm-4.67.1-py3-none-any.whl.metadata (57 kB)
Collecting typing-extensions<5,>=4.11 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0)
  Using cached typing_extensions-4.15.0-py3-none-any.whl.metadata (3.3 kB)
Collecting python-dateutil>=2.8.2 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0)
  Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl.metadata (8.4 kB)
Collecting pytz>=2020.1 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0)
  Using cached pytz-2025.2-py2.py3-none-any.whl.metadata (22 kB)
Collecting tzdata>=2022.7 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0)
  Using cached tzdata-2025.2-py2.py3-none-any.whl.metadata (1.4 kB)
Collecting annotated-types>=0.6.0 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0)
  Using cached annotated_types-0.7.0-py3-none-any.whl.metadata (15 kB)
Collecting pydantic-core==2.41.4 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0)
  Using cached pydantic_core-2.41.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (7.3 kB)
Collecting typing-inspection>=0.4.2 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0)
  Using cached typing_inspection-0.4.2-py3-none-any.whl.metadata (2.6 kB)
Collecting markdown-it-py>=2.2.0 (from rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0)
  Using cached markdown_it_py-4.0.0-py3-none-any.whl.metadata (7.3 kB)
Collecting pygments<3.0.0,>=2.13.0 (from rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0)
  Using cached pygments-2.19.2-py3-none-any.whl.metadata (2.5 kB)
Collecting joblib>=1.2.0 (from scikit-learn<2.0.0,>=1.5.0->gpt-cell-annotator==0.2.0)
  Using cached joblib-1.5.2-py3-none-any.whl.metadata (5.6 kB)
Collecting threadpoolctl>=3.1.0 (from scikit-learn<2.0.0,>=1.5.0->gpt-cell-annotator==0.2.0)
  Using cached threadpoolctl-3.6.0-py3-none-any.whl.metadata (13 kB)
Collecting greenlet>=1 (from sqlalchemy<3.0.0,>=2.0.30->gpt-cell-annotator==0.2.0)
  Using cached greenlet-3.2.4-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl.metadata (4.1 kB)
Collecting mdurl~=0.1 (from markdown-it-py>=2.2.0->rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0)
  Using cached mdurl-0.1.2-py3-none-any.whl.metadata (1.6 kB)
Collecting six>=1.5 (from python-dateutil>=2.8.2->pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0)
  Using cached six-1.17.0-py2.py3-none-any.whl.metadata (1.7 kB)
Using cached anndata-0.10.9-py3-none-any.whl (128 kB)
Using cached httpx-0.27.2-py3-none-any.whl (76 kB)
Using cached httpcore-1.0.9-py3-none-any.whl (78 kB)
Using cached jsonschema-4.25.1-py3-none-any.whl (90 kB)
Using cached numpy-1.26.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (18.3 MB)
Using cached openai-1.109.1-py3-none-any.whl (948 kB)
Using cached anyio-4.11.0-py3-none-any.whl (109 kB)
Using cached distro-1.9.0-py3-none-any.whl (20 kB)
Using cached jiter-0.11.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (359 kB)
Using cached pandas-2.3.3-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl (12.8 MB)
Using cached pyarrow-16.1.0-cp311-cp311-manylinux_2_28_x86_64.whl (40.8 MB)
Using cached pydantic-2.12.3-py3-none-any.whl (462 kB)
Using cached pydantic_core-2.41.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (2.1 MB)
Using cached pydantic_settings-2.11.0-py3-none-any.whl (48 kB)
Downloading python_dotenv-1.2.1-py3-none-any.whl (21 kB)
Using cached pyyaml-6.0.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl (806 kB)
Using cached rich-13.9.4-py3-none-any.whl (242 kB)
Using cached pygments-2.19.2-py3-none-any.whl (1.2 MB)
Using cached scikit_learn-1.7.2-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl (9.7 MB)
Using cached sqlalchemy-2.0.44-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (3.3 MB)
Using cached structlog-24.4.0-py3-none-any.whl (67 kB)
Using cached typing_extensions-4.15.0-py3-none-any.whl (44 kB)
Using cached annotated_types-0.7.0-py3-none-any.whl (13 kB)
Using cached array_api_compat-1.12.0-py3-none-any.whl (58 kB)
Using cached attrs-25.4.0-py3-none-any.whl (67 kB)
Using cached greenlet-3.2.4-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl (587 kB)
Using cached h11-0.16.0-py3-none-any.whl (37 kB)
Using cached h5py-3.15.1-cp311-cp311-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl (4.7 MB)
Using cached idna-3.11-py3-none-any.whl (71 kB)
Using cached joblib-1.5.2-py3-none-any.whl (308 kB)
Using cached jsonschema_specifications-2025.9.1-py3-none-any.whl (18 kB)
Using cached markdown_it_py-4.0.0-py3-none-any.whl (87 kB)
Using cached mdurl-0.1.2-py3-none-any.whl (10.0 kB)
Using cached packaging-25.0-py3-none-any.whl (66 kB)
Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl (229 kB)
Using cached pytz-2025.2-py2.py3-none-any.whl (509 kB)
Using cached referencing-0.37.0-py3-none-any.whl (26 kB)
Downloading rpds_py-0.28.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (382 kB)
Downloading scipy-1.16.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl (35.9 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 35.9/35.9 MB 181.3 MB/s  0:00:00
Using cached six-1.17.0-py2.py3-none-any.whl (11 kB)
Using cached sniffio-1.3.1-py3-none-any.whl (10 kB)
Using cached threadpoolctl-3.6.0-py3-none-any.whl (18 kB)
Using cached tqdm-4.67.1-py3-none-any.whl (78 kB)
Using cached typing_inspection-0.4.2-py3-none-any.whl (14 kB)
Using cached tzdata-2025.2-py2.py3-none-any.whl (347 kB)
Using cached certifi-2025.10.5-py3-none-any.whl (163 kB)
Using cached natsort-8.4.0-py3-none-any.whl (38 kB)
Installing collected packages: pytz, tzdata, typing-extensions, tqdm, threadpoolctl, structlog, sniffio, six, rpds-py, pyyaml, python-dotenv, pygments, packaging, numpy, natsort, mdurl, joblib, jiter, idna, h11, greenlet, distro, certifi, attrs, array-api-compat, annotated-types, typing-inspection, sqlalchemy, scipy, referencing, python-dateutil, pydantic-core, pyarrow, markdown-it-py, httpcore, h5py, anyio, scikit-learn, rich, pydantic, pandas, jsonschema-specifications, httpx, pydantic-settings, openai, jsonschema, anndata, gpt-cell-annotator

Successfully installed anndata-0.10.9 annotated-types-0.7.0 anyio-4.11.0 array-api-compat-1.12.0 attrs-25.4.0 certifi-2025.10.5 distro-1.9.0 gpt-cell-annotator-0.2.0 greenlet-3.2.4 h11-0.16.0 h5py-3.15.1 httpcore-1.0.9 httpx-0.27.2 idna-3.11 jiter-0.11.1 joblib-1.5.2 jsonschema-4.25.1 jsonschema-specifications-2025.9.1 markdown-it-py-4.0.0 mdurl-0.1.2 natsort-8.4.0 numpy-1.26.4 openai-1.109.1 packaging-25.0 pandas-2.3.3 pyarrow-16.1.0 pydantic-2.12.3 pydantic-core-2.41.4 pydantic-settings-2.11.0 pygments-2.19.2 python-dateutil-2.9.0.post0 python-dotenv-1.2.1 pytz-2025.2 pyyaml-6.0.3 referencing-0.37.0 rich-13.9.4 rpds-py-0.28.0 scikit-learn-1.7.2 scipy-1.16.3 six-1.17.0 sniffio-1.3.1 sqlalchemy-2.0.44 structlog-24.4.0 threadpoolctl-3.6.0 tqdm-4.67.1 typing-extensions-4.15.0 typing-inspection-0.4.2 tzdata-2025.2
0.2.0
                               GPT Cell Annotator                               
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Cluster     ┃ Primary Label    ┃ Confidence ┃ Status  ┃ Warnings             ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│  \"CD79A\"  │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"CD2\"    │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"GNLY\"   │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"S100A9\" │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"FCGR3A\" │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"TOP2A\"  │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"PF4\"    │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"HBA1\"   │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
└─────────────┴──────────────────┴────────────┴─────────┴──────────────────────┘
-rw-rw-r-- 1 jameslee nslab 14K Nov  3 12:02 /scratch/tmp.4noEyhSQqf/base.json

### API extra validation
The currently activated Python version 3.10.12 is not supported by the project (^3.11).
Trying to find and use a compatible version. 
Using python3 (3.11.14)
Requirement already satisfied: pip in /scratch/tmp.Yo9Yhi3rWv/api/lib/python3.11/site-packages (24.0)
Collecting pip
  Using cached pip-25.3-py3-none-any.whl.metadata (4.7 kB)
Collecting wheel
  Using cached wheel-0.45.1-py3-none-any.whl.metadata (2.3 kB)
Using cached pip-25.3-py3-none-any.whl (1.8 MB)
Using cached wheel-0.45.1-py3-none-any.whl (72 kB)
Installing collected packages: wheel, pip
  Attempting uninstall: pip
    Found existing installation: pip 24.0
    Uninstalling pip-24.0:
      Successfully uninstalled pip-24.0
Successfully installed pip-25.3 wheel-0.45.1
Processing ./dist/gpt_cell_annotator-0.2.0-py3-none-any.whl (from gpt-cell-annotator==0.2.0)
Collecting anndata<0.11.0,>=0.10.7 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached anndata-0.10.9-py3-none-any.whl.metadata (6.9 kB)
Collecting httpx<0.28.0,>=0.27.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached httpx-0.27.2-py3-none-any.whl.metadata (7.1 kB)
Collecting jsonschema<5.0.0,>=4.25.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached jsonschema-4.25.1-py3-none-any.whl.metadata (7.6 kB)
Collecting numpy<2.0.0,>=1.26.4 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached numpy-1.26.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (61 kB)
Collecting openai<2.0.0,>=1.30.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached openai-1.109.1-py3-none-any.whl.metadata (29 kB)
Collecting pandas<3.0.0,>=2.2.2 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pandas-2.3.3-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl.metadata (91 kB)
Collecting pyarrow<17.0.0,>=16.1.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pyarrow-16.1.0-cp311-cp311-manylinux_2_28_x86_64.whl.metadata (3.0 kB)
Collecting pydantic<3.0.0,>=2.7.1 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pydantic-2.12.3-py3-none-any.whl.metadata (87 kB)
Collecting pydantic-settings<3.0.0,>=2.3.4 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pydantic_settings-2.11.0-py3-none-any.whl.metadata (3.4 kB)
Collecting python-dotenv<2.0.0,>=1.0.1 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached python_dotenv-1.2.1-py3-none-any.whl.metadata (25 kB)
Collecting pyyaml<7.0.0,>=6.0.2 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pyyaml-6.0.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl.metadata (2.4 kB)
Collecting rich<14.0.0,>=13.7.1 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached rich-13.9.4-py3-none-any.whl.metadata (18 kB)
Collecting scikit-learn<2.0.0,>=1.5.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached scikit_learn-1.7.2-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl.metadata (11 kB)
Collecting sqlalchemy<3.0.0,>=2.0.30 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached sqlalchemy-2.0.44-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (9.5 kB)
Collecting structlog<25.0.0,>=24.1.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached structlog-24.4.0-py3-none-any.whl.metadata (7.3 kB)
Collecting fastapi<0.117.0,>=0.116.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached fastapi-0.116.2-py3-none-any.whl.metadata (28 kB)
Collecting prometheus-client<0.21.0,>=0.20.0 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached prometheus_client-0.20.0-py3-none-any.whl.metadata (1.8 kB)
Collecting redis<6.0.0,>=5.0.4 (from gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached redis-5.3.1-py3-none-any.whl.metadata (9.2 kB)
Collecting uvicorn<0.31.0,>=0.30.1 (from uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached uvicorn-0.30.6-py3-none-any.whl.metadata (6.6 kB)
Collecting array-api-compat!=1.5,>1.4 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached array_api_compat-1.12.0-py3-none-any.whl.metadata (2.5 kB)
Collecting h5py>=3.1 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached h5py-3.15.1-cp311-cp311-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl.metadata (3.0 kB)
Collecting natsort (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached natsort-8.4.0-py3-none-any.whl.metadata (21 kB)
Collecting packaging>=20.0 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached packaging-25.0-py3-none-any.whl.metadata (3.3 kB)
Collecting scipy>1.8 (from anndata<0.11.0,>=0.10.7->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached scipy-1.16.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl.metadata (62 kB)
Collecting starlette<0.49.0,>=0.40.0 (from fastapi<0.117.0,>=0.116.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached starlette-0.48.0-py3-none-any.whl.metadata (6.3 kB)
Collecting typing-extensions>=4.8.0 (from fastapi<0.117.0,>=0.116.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached typing_extensions-4.15.0-py3-none-any.whl.metadata (3.3 kB)
Collecting anyio (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached anyio-4.11.0-py3-none-any.whl.metadata (4.1 kB)
Collecting certifi (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached certifi-2025.10.5-py3-none-any.whl.metadata (2.5 kB)
Collecting httpcore==1.* (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached httpcore-1.0.9-py3-none-any.whl.metadata (21 kB)
Collecting idna (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached idna-3.11-py3-none-any.whl.metadata (8.4 kB)
Collecting sniffio (from httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached sniffio-1.3.1-py3-none-any.whl.metadata (3.9 kB)
Collecting h11>=0.16 (from httpcore==1.*->httpx<0.28.0,>=0.27.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached h11-0.16.0-py3-none-any.whl.metadata (8.3 kB)
Collecting attrs>=22.2.0 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached attrs-25.4.0-py3-none-any.whl.metadata (10 kB)
Collecting jsonschema-specifications>=2023.03.6 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached jsonschema_specifications-2025.9.1-py3-none-any.whl.metadata (2.9 kB)
Collecting referencing>=0.28.4 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached referencing-0.37.0-py3-none-any.whl.metadata (2.8 kB)
Collecting rpds-py>=0.7.1 (from jsonschema<5.0.0,>=4.25.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached rpds_py-0.28.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.1 kB)
Collecting distro<2,>=1.7.0 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached distro-1.9.0-py3-none-any.whl.metadata (6.8 kB)
Collecting jiter<1,>=0.4.0 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached jiter-0.11.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (5.2 kB)
Collecting tqdm>4 (from openai<2.0.0,>=1.30.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached tqdm-4.67.1-py3-none-any.whl.metadata (57 kB)
Collecting python-dateutil>=2.8.2 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl.metadata (8.4 kB)
Collecting pytz>=2020.1 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pytz-2025.2-py2.py3-none-any.whl.metadata (22 kB)
Collecting tzdata>=2022.7 (from pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached tzdata-2025.2-py2.py3-none-any.whl.metadata (1.4 kB)
Collecting annotated-types>=0.6.0 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached annotated_types-0.7.0-py3-none-any.whl.metadata (15 kB)
Collecting pydantic-core==2.41.4 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pydantic_core-2.41.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (7.3 kB)
Collecting typing-inspection>=0.4.2 (from pydantic<3.0.0,>=2.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached typing_inspection-0.4.2-py3-none-any.whl.metadata (2.6 kB)
Collecting PyJWT>=2.9.0 (from redis<6.0.0,>=5.0.4->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached PyJWT-2.10.1-py3-none-any.whl.metadata (4.0 kB)
Collecting markdown-it-py>=2.2.0 (from rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached markdown_it_py-4.0.0-py3-none-any.whl.metadata (7.3 kB)
Collecting pygments<3.0.0,>=2.13.0 (from rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached pygments-2.19.2-py3-none-any.whl.metadata (2.5 kB)
Collecting joblib>=1.2.0 (from scikit-learn<2.0.0,>=1.5.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached joblib-1.5.2-py3-none-any.whl.metadata (5.6 kB)
Collecting threadpoolctl>=3.1.0 (from scikit-learn<2.0.0,>=1.5.0->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached threadpoolctl-3.6.0-py3-none-any.whl.metadata (13 kB)
Collecting greenlet>=1 (from sqlalchemy<3.0.0,>=2.0.30->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached greenlet-3.2.4-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl.metadata (4.1 kB)
Collecting click>=7.0 (from uvicorn<0.31.0,>=0.30.1->uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached click-8.3.0-py3-none-any.whl.metadata (2.6 kB)
Collecting httptools>=0.5.0 (from uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached httptools-0.7.1-cp311-cp311-manylinux1_x86_64.manylinux_2_28_x86_64.manylinux_2_5_x86_64.whl.metadata (3.5 kB)
Collecting uvloop!=0.15.0,!=0.15.1,>=0.14.0 (from uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached uvloop-0.22.1-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl.metadata (4.9 kB)
Collecting watchfiles>=0.13 (from uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached watchfiles-1.1.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (4.9 kB)
Collecting websockets>=10.4 (from uvicorn[standard]<0.31.0,>=0.30.1; extra == "api" or extra == "full"->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached websockets-15.0.1-cp311-cp311-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (6.8 kB)
Collecting mdurl~=0.1 (from markdown-it-py>=2.2.0->rich<14.0.0,>=13.7.1->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached mdurl-0.1.2-py3-none-any.whl.metadata (1.6 kB)
Collecting six>=1.5 (from python-dateutil>=2.8.2->pandas<3.0.0,>=2.2.2->gpt-cell-annotator==0.2.0->gpt-cell-annotator==0.2.0)
  Using cached six-1.17.0-py2.py3-none-any.whl.metadata (1.7 kB)
Using cached anndata-0.10.9-py3-none-any.whl (128 kB)
Using cached fastapi-0.116.2-py3-none-any.whl (95 kB)
Using cached httpx-0.27.2-py3-none-any.whl (76 kB)
Using cached httpcore-1.0.9-py3-none-any.whl (78 kB)
Using cached jsonschema-4.25.1-py3-none-any.whl (90 kB)
Using cached numpy-1.26.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (18.3 MB)
Using cached openai-1.109.1-py3-none-any.whl (948 kB)
Using cached anyio-4.11.0-py3-none-any.whl (109 kB)
Using cached distro-1.9.0-py3-none-any.whl (20 kB)
Using cached jiter-0.11.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (359 kB)
Using cached pandas-2.3.3-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl (12.8 MB)
Using cached prometheus_client-0.20.0-py3-none-any.whl (54 kB)
Using cached pyarrow-16.1.0-cp311-cp311-manylinux_2_28_x86_64.whl (40.8 MB)
Using cached pydantic-2.12.3-py3-none-any.whl (462 kB)
Using cached pydantic_core-2.41.4-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (2.1 MB)
Using cached pydantic_settings-2.11.0-py3-none-any.whl (48 kB)
Using cached python_dotenv-1.2.1-py3-none-any.whl (21 kB)
Using cached pyyaml-6.0.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl (806 kB)
Using cached redis-5.3.1-py3-none-any.whl (272 kB)
Using cached rich-13.9.4-py3-none-any.whl (242 kB)
Using cached pygments-2.19.2-py3-none-any.whl (1.2 MB)
Using cached scikit_learn-1.7.2-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl (9.7 MB)
Using cached sqlalchemy-2.0.44-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (3.3 MB)
Using cached starlette-0.48.0-py3-none-any.whl (73 kB)
Using cached structlog-24.4.0-py3-none-any.whl (67 kB)
Using cached typing_extensions-4.15.0-py3-none-any.whl (44 kB)
Using cached uvicorn-0.30.6-py3-none-any.whl (62 kB)
Using cached annotated_types-0.7.0-py3-none-any.whl (13 kB)
Using cached array_api_compat-1.12.0-py3-none-any.whl (58 kB)
Using cached attrs-25.4.0-py3-none-any.whl (67 kB)
Using cached click-8.3.0-py3-none-any.whl (107 kB)
Using cached greenlet-3.2.4-cp311-cp311-manylinux_2_24_x86_64.manylinux_2_28_x86_64.whl (587 kB)
Using cached h11-0.16.0-py3-none-any.whl (37 kB)
Using cached h5py-3.15.1-cp311-cp311-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl (4.7 MB)
Using cached httptools-0.7.1-cp311-cp311-manylinux1_x86_64.manylinux_2_28_x86_64.manylinux_2_5_x86_64.whl (456 kB)
Using cached idna-3.11-py3-none-any.whl (71 kB)
Using cached joblib-1.5.2-py3-none-any.whl (308 kB)
Using cached jsonschema_specifications-2025.9.1-py3-none-any.whl (18 kB)
Using cached markdown_it_py-4.0.0-py3-none-any.whl (87 kB)
Using cached mdurl-0.1.2-py3-none-any.whl (10.0 kB)
Using cached packaging-25.0-py3-none-any.whl (66 kB)
Using cached PyJWT-2.10.1-py3-none-any.whl (22 kB)
Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl (229 kB)
Using cached pytz-2025.2-py2.py3-none-any.whl (509 kB)
Using cached referencing-0.37.0-py3-none-any.whl (26 kB)
Using cached rpds_py-0.28.0-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (382 kB)
Using cached scipy-1.16.3-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.whl (35.9 MB)
Using cached six-1.17.0-py2.py3-none-any.whl (11 kB)
Using cached sniffio-1.3.1-py3-none-any.whl (10 kB)
Using cached threadpoolctl-3.6.0-py3-none-any.whl (18 kB)
Using cached tqdm-4.67.1-py3-none-any.whl (78 kB)
Using cached typing_inspection-0.4.2-py3-none-any.whl (14 kB)
Using cached tzdata-2025.2-py2.py3-none-any.whl (347 kB)
Using cached uvloop-0.22.1-cp311-cp311-manylinux2014_x86_64.manylinux_2_17_x86_64.manylinux_2_28_x86_64.whl (3.8 MB)
Using cached watchfiles-1.1.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (456 kB)
Using cached websockets-15.0.1-cp311-cp311-manylinux_2_5_x86_64.manylinux1_x86_64.manylinux_2_17_x86_64.manylinux2014_x86_64.whl (182 kB)
Using cached certifi-2025.10.5-py3-none-any.whl (163 kB)
Using cached natsort-8.4.0-py3-none-any.whl (38 kB)
Installing collected packages: pytz, websockets, uvloop, tzdata, typing-extensions, tqdm, threadpoolctl, structlog, sniffio, six, rpds-py, pyyaml, python-dotenv, PyJWT, pygments, prometheus-client, packaging, numpy, natsort, mdurl, joblib, jiter, idna, httptools, h11, greenlet, distro, click, certifi, attrs, array-api-compat, annotated-types, uvicorn, typing-inspection, sqlalchemy, scipy, referencing, redis, python-dateutil, pydantic-core, pyarrow, markdown-it-py, httpcore, h5py, anyio, watchfiles, starlette, scikit-learn, rich, pydantic, pandas, jsonschema-specifications, httpx, pydantic-settings, openai, jsonschema, fastapi, anndata, gpt-cell-annotator

Successfully installed PyJWT-2.10.1 anndata-0.10.9 annotated-types-0.7.0 anyio-4.11.0 array-api-compat-1.12.0 attrs-25.4.0 certifi-2025.10.5 click-8.3.0 distro-1.9.0 fastapi-0.116.2 gpt-cell-annotator-0.2.0 greenlet-3.2.4 h11-0.16.0 h5py-3.15.1 httpcore-1.0.9 httptools-0.7.1 httpx-0.27.2 idna-3.11 jiter-0.11.1 joblib-1.5.2 jsonschema-4.25.1 jsonschema-specifications-2025.9.1 markdown-it-py-4.0.0 mdurl-0.1.2 natsort-8.4.0 numpy-1.26.4 openai-1.109.1 packaging-25.0 pandas-2.3.3 prometheus-client-0.20.0 pyarrow-16.1.0 pydantic-2.12.3 pydantic-core-2.41.4 pydantic-settings-2.11.0 pygments-2.19.2 python-dateutil-2.9.0.post0 python-dotenv-1.2.1 pytz-2025.2 pyyaml-6.0.3 redis-5.3.1 referencing-0.37.0 rich-13.9.4 rpds-py-0.28.0 scikit-learn-1.7.2 scipy-1.16.3 six-1.17.0 sniffio-1.3.1 sqlalchemy-2.0.44 starlette-0.48.0 structlog-24.4.0 threadpoolctl-3.6.0 tqdm-4.67.1 typing-extensions-4.15.0 typing-inspection-0.4.2 tzdata-2025.2 uvicorn-0.30.6 uvloop-0.22.1 watchfiles-1.1.1 websockets-15.0.1
0.2.0
                               GPT Cell Annotator                               
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
┃ Cluster     ┃ Primary Label    ┃ Confidence ┃ Status  ┃ Warnings             ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
│  \"CD79A\"  │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"CD2\"    │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"GNLY\"   │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"S100A9\" │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"FCGR3A\" │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"TOP2A\"  │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"PF4\"    │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│  \"HBA1\"   │ Unknown or Novel │    Low     │ flagged │ Proposed label not   │
│             │                  │            │         │ found in marker      │
│             │                  │            │         │ knowledge base.; Low │
│             │                  │            │         │ marker overlap with  │
│             │                  │            │         │ knowledge base       │
│             │                  │            │         │ markers.; Only 0     │
│             │                  │            │         │ marker(s) overlap    │
│             │                  │            │         │ with knowledge base  │
│             │                  │            │         │ (minimum 2           │
│             │                  │            │         │ required).; Ontology │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
│             │                  │            │         │ payload.; Label      │
│             │                  │            │         │ absent from marker   │
│             │                  │            │         │ database; Ontology   │
│             │                  │            │         │ identifier missing   │
│             │                  │            │         │ from annotation      │
└─────────────┴──────────────────┴────────────┴─────────┴──────────────────────┘
INFO:     Started server process [3944126]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://127.0.0.1:8081 (Press CTRL+C to quit)
OPENAI_API_KEY not configured; using heuristic MockAnnotator. Outputs are suitable for demos only.
{"timestamp": "2025-11-03T17:03:38.405170Z", "level": "info", "trace_id": "ce5d1333-410b-45b8-b780-2888feeb25cd", "status_code": 200, "path": "/health", "method": "GET", "duration_ms": 4.180473741143942, "message": "request.completed"}
INFO:     127.0.0.1:51110 - "GET /health HTTP/1.1" 200 OK
{"status":"ok","llm_mode":"mock","cache_enabled":false}