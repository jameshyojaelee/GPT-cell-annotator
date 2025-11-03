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
