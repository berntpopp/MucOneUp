# Issue #28: Reference Assembly Management

**Status**: ✅ COMPLETED
**Commit**: 7de6267 (feat: Add reference genome assembly management system)
**Date**: 2025-10-19
**Branch**: feat/issue-28-reference-assemblies

---

## Overview

Implemented centralized reference assembly management system following VNtyper's pragmatic config-driven approach. Enables support for multiple reference assemblies (hg19, hg38, GRCh37, GRCh38) with proper FASTA path and VNTR region mapping.

## Implementation Approach

**Chosen Approach**: Simple 3-function extension to existing validation module
**Implementation Plan**: `implementation_plan.md` (FINAL_REVISED_SIMPLE)
**Rejected Approach**: `comprehensive_solution_rejected.md` (over-engineered)

### Why the Simple Approach Won

The FINAL_REVISED_SIMPLE approach (~100 lines) was chosen over the comprehensive solution (500+ lines) because:

- ✅ **KISS**: Simple config-based approach, no unnecessary classes
- ✅ **DRY**: Reuses existing `validate_reference_genome()` function
- ✅ **SOLID**: Extends validation module (Single Responsibility maintained)
- ✅ **Backward Compatible**: Auto-migration from old config format
- ✅ **Maintainable**: Clear, focused code without feature creep

The comprehensive solution violated these principles by:
- ❌ Creating unnecessary `ReferenceManager` class (forbidden in issue)
- ❌ Duplicating existing validation functionality
- ❌ Adding feature creep (mouse genomes, provenance tracking, etc.)
- ❌ Over-engineering with unnecessary abstractions

## Key Changes

### 1. Config Schema (`config.py`)
- Added `reference_genomes` section with assembly-specific configuration
- Auto-migration from old `constants` format for backward compatibility
- Flexible assembly names via JSON schema pattern properties

### 2. Reference Validation (`reference_validation.py`)
Added 3 helper functions (~175 lines):
- `get_reference_path_for_assembly()`: Resolves FASTA path for assembly
- `validate_reference_for_assembly()`: Validates reference and indices
- `get_muc1_region_for_assembly()`: Returns VNTR region coordinates

### 3. Pipeline Integration
- **ONT pipeline**: Uses new helpers with minimap2 validation
- **Illumina pipeline**: Uses new helpers with BWA validation
- Graceful fallback to old `human_reference` config if needed

### 4. Testing
- 23 comprehensive unit tests (all passing)
- 89% coverage for `reference_validation.py`
- No regressions in existing functionality

### 5. Example Configuration
Added `reference_genomes` section to `config.json`:
```json
{
  "reference_genomes": {
    "hg38": {
      "fasta_path": "reference/hg38/hg38.fa",
      "vntr_region": "chr1:155188487-155192239",
      "display_name": "Human Reference Genome (GRCh38/hg38)",
      "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    },
    "hg19": {
      "fasta_path": "reference/hg19/hg19.fa",
      "vntr_region": "chr1:155160963-155162030",
      "display_name": "Human Reference Genome (GRCh37/hg19)",
      "source_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    }
  }
}
```

## CI Validation

All checks passing:
- ✅ Ruff linting: All checks passed
- ✅ Ruff formatting: 85 files formatted
- ✅ Mypy type checking: No issues in 47 source files
- ✅ Full test suite: 676/676 tests passed
- ✅ Coverage: 89% for new module (target: >80%)

## Files Modified

- `config.json` - Added reference_genomes section
- `muc_one_up/bioinformatics/reference_validation.py` - Added 3 helper functions
- `muc_one_up/config.py` - Updated schema + auto-migration
- `muc_one_up/read_simulator/ont_pipeline.py` - Integrated helpers
- `muc_one_up/read_simulator/pipeline.py` - Integrated helpers
- `tests/bioinformatics/test_reference_validation.py` - New comprehensive tests

**Total**: 6 files, 561 insertions, 7 deletions

## Documentation

- `implementation_plan.md` - Final implementation plan (used)
- `comprehensive_solution_rejected.md` - Over-engineered approach (rejected)
- `original_issue.md` - Original issue description from GitHub

## Benefits

1. **Centralized Management**: Single source of truth for assembly configuration
2. **Backward Compatible**: Existing configs work via auto-migration
3. **Extensible**: Easy to add new assemblies (mm10, etc.)
4. **Validated**: Comprehensive test coverage with no regressions
5. **Simple**: Follows KISS/DRY/SOLID principles
6. **Maintainable**: Clear, focused code without unnecessary complexity

## Lessons Learned

- **Over-engineering is real**: First attempt created 500+ lines when 100 sufficed
- **KISS matters**: Simple solutions are easier to test, maintain, and understand
- **DRY prevents bugs**: Reusing existing validation prevented code duplication
- **Planning review crucial**: Ultra-thinking review caught anti-patterns early
- **Tests enable confidence**: 89% coverage ensures robust implementation

## Related Issues

- Implements: #28 (Reference assembly management)
- Related to: #30 (Deterministic simulation - uses reference_genomes config)

## Next Steps

- Merge feature branch to main
- Update user documentation
- Add support for additional assemblies as needed (mm10, etc.)
