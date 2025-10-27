# Feature: VNTR-Biased FASTQ Generation

## Summary
Automatically generates paired-end FASTQ files from VNTR-corrected BAM files, enabling realistic variant caller benchmarking with empirically validated coverage distributions.

## Implementation (Issue #46)

### Core Changes
- **samtools_wrapper.py** (+207 lines):
  - `FastqConversionOptions` dataclass (5 parameters)
  - `convert_bam_to_paired_fastq()` function (160 lines)
  - `_count_fastq_reads()` helper (27 lines)
  - **Key feature**: Read pair integrity validation

- **pipeline.py** (+46 lines):
  - Integration after VNTR bias application
  - Configuration-driven (default: enabled)
  - Non-fatal error handling

### Testing
- **16 unit tests** (all passing)
- **987 total tests** (0 failures)
- Coverage maintained at 7%

### Quality Gates ✅
- Linting: All checks passed
- Type checking: No issues (58 source files)
- CI: All tests passing
- Design: SOLID, DRY, KISS compliant

## Usage

### Default Behavior (Enabled)
```bash
muconeup --config config.json reads illumina sample.fa --seed 42
```

Generates:
- `sample_R1.fastq.gz`, `sample_R2.fastq.gz` (unrealistic coverage)
- `sample_vntr_biased.bam` (realistic coverage)
- `sample_vntr_biased_R1.fastq.gz` ✨ **NEW**
- `sample_vntr_biased_R2.fastq.gz` ✨ **NEW**

### Configuration (Optional)
```json
"vntr_capture_efficiency": {
  "enabled": true,
  "penalty_factor": 0.375,
  "output_fastq": {
    "enabled": true,
    "preserve_read_names": true,
    "singleton_file": null
  }
}
```

## Technical Details

### Design Improvements (v2.0)
- **40% code reduction** (removed over-engineered compression)
- **50% parameter reduction** (10 → 5 via dataclass)
- **Added read pair validation** (prevents data corruption)
- **Simplified compression** (samtools auto-compresses .gz files)

### Key Features
- ✅ Read pair integrity validation (R1/R2 read count matching)
- ✅ Auto-compression based on file extension
- ✅ Parent directory auto-creation
- ✅ Backward compatible (enabled by default)
- ✅ Non-fatal errors (continues with BAM if FASTQ fails)

## Impact
Enables benchmarking variant callers with realistic VNTR coverage (~3.5x VNTR:flanking ratio matching real Twist v2 data).
