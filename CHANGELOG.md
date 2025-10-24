# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

#### Illumina Coverage Downsampling Now Works

- **Fixed config key mismatch** that prevented Illumina read downsampling from executing
- Pipeline now correctly reads `coverage` configuration key set by CLI
- Downsampling logic (previously non-functional) now activates properly

**Technical Details:**
- Standardized on `coverage` key across config, CLI, and pipeline
- Fixed log format string bug (changed `%d` to `%.2f` for float values)
- No changes to CLI interface (already used correct key)
- No changes to test fixtures (already used correct key)

### Changed

#### Breaking Change: Config Key Renamed

- **Renamed:** `downsample_coverage` → `coverage` in `read_simulation` section
- **Reason:** Standardize with ONT/PacBio pipelines, fix broken downsampling feature
- **Impact:** User configuration files must be updated

**Migration Required:**

If your `config.json` contains:
```json
{
  "read_simulation": {
    "downsample_coverage": 150
  }
}
```

Change to:
```json
{
  "read_simulation": {
    "coverage": 150
  }
}
```

**Note:** The old key was non-functional due to a bug, so this change enables the feature rather than breaking working functionality. Users who never used downsampling are unaffected.

**Files Changed:**
- `config.json`: Updated key name (line 545)
- `muc_one_up/config.py`: Updated schema validation (line 296)
- `muc_one_up/read_simulator/pipeline.py`: Read correct key (lines 279-331)

**Testing:**
- Added comprehensive test suite: `tests/read_simulator/test_coverage_config.py`
- 20+ test cases covering key standardization, validation, and edge cases
- All existing tests pass (already used correct key)

**Documentation:**
- Implementation plan: `plan/fix-illumina-coverage-downsampling.md`

---

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
