# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.46.0] - 2026-09-24

Homopolymer runs without a stutter-table entry (e.g. the dupC C8 run) are no
longer simulated error-free by the ONT read profiles (#132). **Breaking change
for custom read profiles with an `errors` section:** they now need
`molecules.stutter_fallback`. See `docs/about/changelog.md`.

## [0.45.0] - 2026-09-24

Realistic, truth-tracked long-read simulation (read profiles, molecule model,
empirical error channel, per-read truth, fragment-based genomic ONT) and fixes
#97, #102, #104-#124. The maintained changelog is `docs/about/changelog.md`.

## [0.28.1] - 2026-03-16

### Fixed

- **Add biopython to core dependencies** — `biopython` was only listed in optional `dev` and `docs` dependency groups, but is imported at module level in `read_simulator/utils/reference_utils.py` and used in CLI commands, causing `ModuleNotFoundError: No module named 'Bio'` on fresh installs

## [0.20.0] - 2025-10-23

### Added
- Tool version tracking and metadata improvements

### Changed
- Update gitignore configuration

### Fixed
- Configure Trivy to ignore unfixable base image vulnerabilities
- Remove bulk dismissal script (no longer needed)

---

*For complete version history, see git log.*
