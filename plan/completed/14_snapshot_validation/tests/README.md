# SNaPshot Validation - Prototype Tests

**Status**: Prototype/Proof-of-Concept Tests
**Purpose**: Validate SNaPshot mechanism on real dupC data

---

## Test Files

### Main Test (Working Prototype)

**`test_complete_snapshot_workflow.py`** (593 lines)
- Complete end-to-end workflow validation
- Tests on real dupC mutant and normal samples
- **Status**: ✓✓ ALL TESTS PASSING
- Implements complete logic with size filtering and mutation prioritization

**Test Coverage**:
1. Mutant sample: Correctly identifies 8C mutation (Black peak)
2. Normal sample: No false positives (all 7C digested)
3. Mechanism validation: 8C survives digest, 7C destroyed
4. Fluorescence detection: Correct Black (dTAMRA) peak for 8C

---

### Supporting Tests

**`test_digest_selection_logic.py`**
- Validates MwoI digest mechanism
- Confirms 7C has restriction site
- Confirms 8C lacks restriction site
- Tests amplicon-specific digest selection

**`test_protocol_validation.py`**
- Protocol compliance tests
- Primer sequence validation
- Amplicon structure verification

**`test_snapshot_validation_real.py`**
- Initial validation tests
- Basic workflow testing

**`test_all_pcr_tools.py`**
- Comparison of PCR simulation tools
- Evaluates different primer binding approaches

---

## Important Notes

### These are PROTOTYPE tests, not production tests

**Limitations**:
1. ⚠️ Contains **hardcoded values** (primers, patterns, sizes)
2. ⚠️ Does NOT use config.json
3. ⚠️ Not following DRY/SOLID principles
4. ⚠️ Violates project architecture standards

### Production Implementation Required

Before production use, the logic from these tests must be:
1. Refactored into modular classes (`PCRSimulator`, `DigestSimulator`, etc.)
2. Config-based (all parameters from config.json)
3. Properly integrated into `muc_one_up/analysis/snapshot_validator.py`
4. Updated tests using config instead of hardcoded values

See **PLAN.md** for complete config-based implementation specification.

---

## Usage

These tests validate the SNaPshot mechanism on real data:

```bash
# Run main test
python -m pytest plan/issue_15_snapshot_validation/tests/test_complete_snapshot_workflow.py -v

# Run specific test
python -m pytest plan/issue_15_snapshot_validation/tests/test_complete_snapshot_workflow.py::test_dupC_mutant_sample -v

# Run all tests in folder
python -m pytest plan/issue_15_snapshot_validation/tests/ -v
```

---

## Test Data

**Required Files**:
- `output/dupC/dupC.001.mut.simulated.fa` - Mutant sample (positive control)
- `output/dupC/dupC.001.normal.simulated.fa` - Normal sample (negative control)

Generate test data:
```bash
muconeup --config config.json simulate --out-base dupC --mutation-name normal,dupC --fixed-lengths 27 --out-dir output/dupC/
```

---

## Expected Results

### Mutant Sample
```
✓ PCR products: 48 (size-filtered 50-65bp)
✓ Digest survivors: 1 (8C amplicon, 56bp)
✓ Mutation detected: YES
✓ Fluorescence: Black (dTAMRA)
✓ Mechanism: 8C lacks MwoI site → survives digest
```

### Normal Sample
```
✓ PCR products: 48 (size-filtered 50-65bp)
✓ Digest survivors: 0 (all 7C digested)
✓ Mutation detected: NO
✓ All products destroyed by MwoI
```

---

## Next Steps

1. Implement config-based production version (see PLAN.md)
2. Create proper unit tests for each modular class
3. Create integration tests using config
4. Deprecate these prototype tests
5. Move to standard tests/ folder once production-ready

---

**Status**: Prototype - Mechanism Validated - Production Implementation Needed
