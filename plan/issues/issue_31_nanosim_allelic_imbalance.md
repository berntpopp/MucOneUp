# Issue #31: NanoSim Allelic Imbalance Bug

## Problem
NanoSim read simulation produces allelic imbalance when simulating diploid haplotypes with very divergent sizes (e.g., 20 vs 100 repeats). Coverage is not uniformly distributed between haplotypes.

## Root Cause
When NanoSim simulates from a multi-sequence FASTA file, it samples reads based on total sequence length. Longer haplotypes receive proportionally more reads, violating the expected 50:50 diploid coverage.

## Research Background
- **Long-read sequencing bias**: Studies show length-dependent sampling in nanopore sequencing simulators (PEPPER-Margin-DeepVariant, Nature Methods 2021)
- **Diploid simulation**: True diploid samples should have equal allelic coverage regardless of haplotype length differences
- **MUC1 VNTR context**: From bioRxiv 2025.09.06.673538v2, MUC1 VNTR lengths vary 20-125 repeats, making this issue critical for ADTKD-MUC1 studies

## Implementation

### 1. Diagnosis Module (`muc_one_up/read_simulator/allelic_balance_check.py`)
```python
def check_allelic_balance(bam_file: str, haplotype_names: list[str]) -> dict[str, float]:
    """Calculate coverage per haplotype from aligned BAM."""
    # Use pysam to count reads mapping to each haplotype
    # Return {haplotype_name: coverage}
```

### 2. Fix in `ont_pipeline.py`
**Option A: Separate simulation per haplotype (RECOMMENDED)**
```python
# Simulate each haplotype independently with target_coverage / num_haplotypes
for haplotype in haplotypes:
    run_nanosim_simulation(
        reference_fasta=single_haplotype_fa,
        coverage=target_coverage / len(haplotypes)
    )
# Merge FASTQ files
```

**Option B: Weighted sampling** (requires NanoSim modification - NOT recommended)

### 3. Config Parameter
Add to `nanosim_params`:
```json
"enforce_allelic_balance": true,  // Default: true
"allelic_balance_tolerance": 0.1  // Allow 10% deviation
```

### 4. Validation
```python
def validate_allelic_balance(bam: str, tolerance: float = 0.1):
    """Post-simulation check. Warn if coverage ratio exceeds tolerance."""
    # If hap1=45x, hap2=55x with target=50x: OK (within 10%)
    # If hap1=30x, hap2=70x: FAIL (20% deviation)
```

## Testing
- **Unit test**: Mock 2 haplotypes (1kb vs 5kb), verify equal read counts
- **Integration test**: Simulate MUC1 20-repeat vs 100-repeat, check BAM coverage
- **Regression test**: Existing balanced haplotypes should remain unaffected

## Documentation
- Update README with `enforce_allelic_balance` parameter
- Add troubleshooting section for coverage imbalance

## Files Modified
- `muc_one_up/read_simulator/ont_pipeline.py` (split simulation)
- `muc_one_up/read_simulator/allelic_balance_check.py` (NEW)
- `config.json` (new parameters)
- `tests/read_simulator/test_allelic_balance.py` (NEW)
