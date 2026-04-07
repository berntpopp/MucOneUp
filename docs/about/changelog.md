# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Fixed
- ONT amplicon pipeline now reads from dedicated `ont_amplicon_params` config section instead of `pacbio_params` -- was using PacBio Sequel model (ERRHMM-SEQUEL) for ONT reads (#79)
- CLI routes `--model-file` overrides to `ont_amplicon_params` for ONT and `pacbio_params` for PacBio

### Added
- `ont_amplicon_params` config section with ONT-specific defaults (QSHMM-ONT-HQ.model)
- `OntAmpliconConfig` TypedDict for type-safe ONT parameter access
- Config schema validation and path resolution for `ont_amplicon_params`
- Test for ONT platform using `ont_amplicon_params` (not `pacbio_params`)
- Documentation: ONT configuration section in amplicon simulation guide

---

## [0.43.4] - 2026-04-06

### Fixed
- Race condition in VNTR efficiency temp directory — replaced deterministic path with `tempfile.TemporaryDirectory` for process-safe isolation in parallel runs (#82)

---

## [0.43.3] - 2026-04-06

### Fixed
- ONT amplicon: warn when model file does not appear to be an ONT model (#79)
- pbsim3 wrapper: detect and convert numbered SAM files (`prefix_0001.sam`) produced with template mode (#79)
- Config paths (model_file, reseq_model, human_reference) now resolved relative to config file directory, not working directory (#80)
- Default `expected_product_range` widened from `[1500, 6000]` to `[1500, 15000]` to support VNTRs with >99 repeats (#81)

### Added
- 5 new tests for config path resolution, ONT model warning, and numbered SAM detection

---

## [0.43.2] - 2026-04-06

### Added
- `--read-number` CLI option for `reads illumina` to control fragment generation count
- 3 new CLI tests for `--read-number` option

### Changed
- Default `read_number` raised from 10,000 to 100,000 in config.json (was limited by original ReSeq deadlock bug, now resolved by ReSeq2)

---

## [0.43.1] - 2026-04-05

### Changed
- Replaced ReSeq with [ReSeq2](https://github.com/berntpopp/ReSeq2) as the recommended Illumina error modeling tool
- Removed `seqToIllumina` timeout workaround — ReSeq2 v2.0.3 fixes the upstream deadlock bug that limited the original ReSeq to ~10,000 reads
- Default `seqToIllumina` timeout raised from 120s to 600s (safety net only)

### Removed
- `seqtoillumina_timeout` config parameter (no longer needed)

---

## [0.43.0] - 2026-04-05

### Added
- ONT amplicon simulation via `reads amplicon --platform ont`
- `--platform` option for `reads amplicon` command (choices: `pacbio`, `ont`)
- Shared amplicon preparation module (`amplicon_common.py`) used by both platforms
- `Assay_type` field in simulation metadata TSV
- `assay_type` and `ont-amplicon` to config schema
- `keep_intermediate_files` support in ONT amplicon pipeline
- 13 new tests for ONT amplicon pipeline, shared helpers, CLI routing, and metadata

### Changed
- Relaxed pbsim3 `pass_num` validation from `>= 2` to `>= 1` (ONT single-pass)
- Refactored PacBio amplicon pipeline to use shared extraction stages
- Metadata platform field changed from `PacBio-Amplicon` to `PacBio` with separate `assay_type`

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

[Unreleased]: https://github.com/berntpopp/MucOneUp/compare/v0.40.0...HEAD
[0.40.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.40.0
[0.19.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.19.0
[0.15.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.15.0
[0.14.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.14.0
[0.13.0]: https://github.com/berntpopp/MucOneUp/releases/tag/v0.13.0
