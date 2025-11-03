# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- Placeholder for future changes. Add new entries below this heading.

## [0.2.0] - 2025-11-03

- Expanded Scanpy integration with asynchronous batching, caching, offline telemetry, and a richer CLI (`gca scanpy`) companion. See `docs/scanpy_integration.md` for the full notebook and CLI walkthrough.
- Hardened annotation validation through structured reports, guardrail overrides, and disk-backed caches that power repeatable offline runs.
- Published an R/Seurat companion package (`clients/r/gptcellannotator`) with REST + CLI workflows, pkgdown documentation, and CI coverage—refer to `docs/seurat_integration.md` for usage guidance.
- Refreshed distribution docs, installation guides, and offline checklists to highlight extras (`[api]`, `[ui]`, `[scanpy]`, `[full]`) and cross-platform setup (`docs/install.md`).
- Automated release assets and QA pipelines, including PyPI-ready wheels/sdists, Docker builds, and benchmark artefact capture under `docs/reports/`.
