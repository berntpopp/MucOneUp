# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.40.0] - 2026-04-04

### Added
- PacBio amplicon read simulation using PBSIM3 template mode (`reads amplicon` command)
- PCR length bias model with exponential decay, calibrated to Madritsch et al. 2026 empirical data
- Deterministic and stochastic (Galton-Watson) PCR bias modes
- Preset profiles (`default`, `no_bias`) for PCR bias configuration
- Primer-based amplicon extraction from diploid VNTR references
- Shared primer binding site utility (refactored from snapshot validator)
- `amplicon_params` configuration section with primer sequences and PCR bias settings
- 69 new tests covering all amplicon simulation components

### Changed
- Updated GitHub Actions to Node.js 24 compatible versions (checkout v6, setup-python v6, setup-uv v7)
- Extended config schema to accept `"amplicon"` and `"pacbio"` as simulator types

---

## [0.19.0] - 2025-10-20

### Added
- MkDocs Material documentation system with GitHub Actions deployment
- Comprehensive user guides (simulation, toxic protein detection, SNaPshot validation)
- Auto-generated CLI documentation via mkdocs-click
- Dark/light mode toggle, instant search, mobile responsive design
- Professional documentation at https://berntpopp.github.io/MucOneUp/

---

## [0.15.0] - 2025-10-18

### Added
- In silico SNaPshot assay validation for MUC1 dupC mutation

---

## [0.14.0] - 2025-10-15

### Added
- Diploid split-simulation for ONT reads

### Fixed
- ONT read simulation bias in diploid references

---

## [0.13.0] - 2025-10-10

### Added
- Toxic protein detection algorithm
- ORF prediction with orfipy integration

---

## Earlier Versions

See [GitHub Releases](https://github.com/berntpopp/MucOneUp/releases) for earlier version history.

---

[Unreleased]: https://github.com/berntpopp/MucOneUp/compare/v0.19.0...HEAD
[0.19.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.19.0
[0.15.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.15.0
[0.14.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.14.0
[0.13.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.13.0
