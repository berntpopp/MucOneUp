# Changelog

All notable changes to MucOneUp will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.12.0] - 2025-10-18

### Added
- **Progress indicators** for `--simulate-series` mode - shows ETA and completion percentage for long-running series simulations
- **Verbose flag** (`--verbose` / `-v`) as intuitive alias for `--log-level DEBUG`
- **Comprehensive batch processing tests** - DRY parametrized tests for reads/analyze commands
- **ONT test parity** - `test_reads_ont_with_file` mirrors illumina test coverage
- **Shell completion documentation** in README for Bash, Zsh, and Fish shells

### Improved
- Test coverage increased from 148 to 158 CLI-specific tests (+10 tests)
- Code quality: Eliminated code duplication using `nullcontext` pattern (DRY principle)
- Code simplicity: Verbose flag uses simple if-statement instead of callback (KISS principle)
- Better user experience for long-running simulations with visual progress tracking
- Tab completion support documented for all shells

### Changed
- Progress bar appears automatically for series mode (multiple iterations)
- Single simulations continue to run without progress bar (clean output)

### Technical Details
- **DRY**: Used `@pytest.mark.parametrize` for batch tests (5 tests → 1 parametrized test)
- **DRY**: Used `contextlib.nullcontext` to avoid loop duplication in progress bar code
- **KISS**: Verbose flag uses simple precedence logic (3 lines vs callback pattern)
- **YAGNI**: Shell completion uses Click's built-in feature (no custom command needed)

### Developer Notes
- All changes follow Unix philosophy (single responsibility principle)
- SOLID principles maintained throughout
- Zero breaking changes - fully backward compatible
- All 558 existing tests continue to pass

## [0.11.0] - 2025-10-XX

### Note
Version 0.11.0 was the pre-implementation state.

## [0.9.0] - 2025-10-XX

### Changed
- Migrated from argparse to Click for CLI framework
- Restructured CLI following Unix philosophy (single responsibility)
- Improved architecture with SOLID principles

### Added
- Clean separation: `simulate`, `reads`, `analyze` commands
- Batch processing support for reads/analyze commands
- `vntr-stats` analysis command
- Comprehensive CLI test suite (148 tests)

## [Previous versions]

See git history for earlier changes before Click migration.

---

## Versioning Notes

- **0.12.0**: UX improvements (progress, verbose, completion docs, batch testing)
- **0.11.0**: Pre-implementation baseline
- **0.9.0**: Architecture refactor (argparse → Click)
- Earlier versions used different CLI structure
