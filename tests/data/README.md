# Test Data Files

This directory contains test data files used by the MucOneUp test suite.

## BAM Files

### `sample.001.normal.simulated_vntr_biased.bam`
- **Size**: 297KB (304,128 bytes)
- **Index**: `sample.001.normal.simulated_vntr_biased.bam.bai` (76KB)
- **Source**: Copied from `output/experminents/experiment_1/pair_001/`
- **Purpose**: Test BAM file with unpaired reads from VNTR downsampling
- **Used by**: `tests/test_bam_to_fastq_conversion.py`

#### File Characteristics:
- Total reads: 4,906
- Read1: 2,450
- Read2: 2,456
- Unpaired: 84 reads
- After proper pairing: 2,373 reads per file

This BAM file simulates the exact scenario encountered during VNTR-biased downsampling where some reads lose their mates. It's used to test the BAM-to-FASTQ conversion pipeline's ability to handle unpaired reads correctly using samtools collate and singleton handling.

## Why Commit Test Data?

Test data files are committed to the repository for:
1. **Reproducibility**: Ensures all developers run tests against identical data
2. **CI Compatibility**: GitHub Actions can run tests without generating data
3. **Speed**: Tests run faster without data generation overhead
4. **Reliability**: Known-good test data with verified characteristics

## Size Considerations

The 300KB BAM file is well within acceptable limits for test data:
- Small enough for git (< 1MB)
- Large enough to test real-world scenarios
- Contains sufficient complexity (84 unpaired reads)
- Represents actual VNTR downsampling output

## Maintenance

If test data needs updating:
1. Copy new files from experiment outputs
2. Update this README with new statistics
3. Update test fixtures if file characteristics change
4. Run full test suite to verify compatibility

---

**Last Updated**: 2025-10-28
**MucOneUp Version**: v0.26.0
