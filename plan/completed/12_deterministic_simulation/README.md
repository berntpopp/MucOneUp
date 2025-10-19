# Issue #30: Deterministic Read Simulation

**Status**: ✅ COMPLETED
**Commit**: ca769db (feat: Add random seed support for deterministic read simulation)
**PR**: #33
**Date**: 2025-10-19

---

## Overview

Implemented random seed support for deterministic read simulation, enabling reproducible research, fair benchmarking, and consistent debugging. Users can now specify `--seed` parameter to get identical read output across multiple runs.

## Implementation

### Key Changes

1. **CLI Support** (`cli/click_main.py`)
   - Added `--seed` parameter to `reads` command group
   - Seed propagates to both Illumina and ONT pipelines

2. **Illumina Pipeline** (`read_simulator/pipeline.py`)
   - Seed parameter passed to fragment simulation
   - Deterministic fragment generation

3. **ONT Pipeline** (`read_simulator/ont_pipeline.py`)
   - Seed parameter passed to NanoSim wrapper
   - Deterministic ONT read generation

4. **Fragment Simulation** (`read_simulator/fragment_simulation.py`)
   - Uses `random.seed()` for reproducible fragment selection
   - Deterministic insert length sampling

5. **NanoSim Wrapper** (`read_simulator/wrappers/nanosim_wrapper.py`)
   - Passes `--seed` to NanoSim simulator.py
   - Ensures reproducible ONT read characteristics

6. **Configuration Support**
   - Added `seed` parameter to `read_simulation` config section
   - Added `seed` parameter to `nanosim_params` config section

### Testing

- CLI integration tests for `--seed` parameter
- Verified identical output with same seed
- Verified different output with different seeds
- Tested both Illumina and ONT pipelines

### Documentation

- Updated CLAUDE.md with "Deterministic Read Simulation" section
- Added usage examples and caveats
- Documented configuration options

## Usage

```bash
# Illumina reads with seed (reproducible)
muconeup --config config.json reads illumina sample.fa --seed 42

# ONT reads with seed (reproducible)
muconeup --config config.json reads ont sample.fa --seed 42

# Full pipeline with seed
muconeup --config config.json simulate --seed 42 --out-base sample
muconeup --config config.json reads illumina sample.001.simulated.fa --seed 42
```

## Configuration Example

```json
{
  "read_simulation": {
    "simulator": "illumina",
    "seed": 42
  },
  "nanosim_params": {
    "training_data_path": "/path/to/model",
    "coverage": 30,
    "seed": 42
  }
}
```

## Important Caveats

Same seed guarantees identical output ONLY when:
- Using identical input files
- Running on same platform/architecture
- Using same tool versions (NanoSim, Python, reseq)

## Benefits

1. **Reproducible Research**: Identical reads across runs for publication
2. **Fair Benchmarking**: Compare tools with identical input
3. **Debugging**: Consistent behavior for troubleshooting
4. **Testing**: Deterministic unit tests

## Files Modified

- `muc_one_up/cli/click_main.py` - Added --seed parameter
- `muc_one_up/read_simulator/pipeline.py` - Seed propagation
- `muc_one_up/read_simulator/ont_pipeline.py` - Seed propagation
- `muc_one_up/read_simulator/fragment_simulation.py` - Seed usage
- `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py` - Seed to NanoSim
- `tests/cli/test_click_cli.py` - Integration tests
- `CLAUDE.md` - Documentation

## Related Issues

- Implements: #30 (Deterministic simulation)
- Related to: #28 (Reference assembly management)

## References

- Planning: `implementation_plan.md`
- Summary: `implementation_summary.md`
- Quick Reference: `quick_reference.md`
