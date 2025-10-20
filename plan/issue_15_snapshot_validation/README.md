# Issue #15: SNaPshot Validation

**Complete planning and validation for in-silico SNaPshot assay validation**

---

## Quick Links

- **[PLAN.md](PLAN.md)** - Complete implementation plan with config-based architecture
- **[STATUS.md](STATUS.md)** - Current status, achievements, and validation results
- **[tests/](tests/)** - Prototype tests (mechanism validated ✓✓)

---

## Overview

Implement complete in-silico SNaPshot assay validation for MUC1 VNTR mutations, specifically the dupC (8C) mutation. The implementation simulates the complete wet-lab workflow:

```
PCR Amplification → MwoI Digest Selection → SNaPshot Extension → Fluorescence Detection
```

---

## Status Summary

| Component | Status | Notes |
|-----------|--------|-------|
| **Mechanism Validation** | ✅ Complete | 8C disrupts MwoI site - confirmed on real data |
| **Test Implementation** | ✅ Complete | All tests passing (mutant + normal samples) |
| **Architecture Design** | ✅ Complete | Config-based, modular, follows SOLID principles |
| **Production Code** | ⏭️ Pending | Needs config integration (~10 hours) |

---

## Key Findings

### ✅ Validated Mechanism

**MwoI Recognition**: `GCNNNNNNNGC` (exactly 7 bases between GC...GC)

```
7C NORMAL:  GC[CCCCCCC]AGC  → Has MwoI site → DIGESTED ✓
8C MUTANT:  GC[CCCCCCCC]AGC → No MwoI site  → SURVIVES ✓
```

**Result**: The extra C in dupC mutation disrupts the MwoI restriction site, enabling digest-based selection of mutant amplicons!

### ✅ Complete Test Validation

**Test File**: `tests/test_complete_snapshot_workflow.py`

**Results**:
- Mutant sample: ✓ Correctly detects 8C mutation (Black fluorescence)
- Normal sample: ✓ No false positives (all 7C digested)
- Mechanism: ✓ 8C survives, 7C destroyed
- Performance: ✓ < 1 second complete workflow

---

## Implementation Plan

### Architecture (See PLAN.md)

**Modular Classes**:
1. `PCRSimulator` - PCR amplification with multi-product handling
2. `DigestSimulator` - Restriction digest selection
3. `SnapshotExtensionSimulator` - Single-base extension simulation
4. `SnapshotValidator` - Workflow orchestrator

**Design Principles**:
- ✅ Config-based (no hardcoded values)
- ✅ SOLID principles (SRP, OCP, LSP, ISP, DIP)
- ✅ DRY (single source of truth in config)
- ✅ KISS (simple, focused classes)
- ✅ Modular (independently testable components)

### Configuration (Required Before Implementation)

Add to `config.json`:
```json
{
  "snapshot_validation": {
    "dupC": {
      "pcr": {
        "forward_primer": "GGCCGGCCCCGGGCTCCACC",
        "reverse_primer": "TGTCACCTCGGCCCCGGA",
        "reverse_needs_rc": true,
        "max_products": 100,
        "size_range": {"min": 50, "max": 65}
      },
      "digest": {
        "enzyme": "MwoI",
        "recognition_site": "GCNNNNNNNGC"
      },
      "snapshot": {
        "primers": {...},
        "fluorophore_map": {...}
      },
      "validation": {
        "mutant_pattern": "GCCCCCCCCAGC",
        "normal_pattern": "GCCCCCCCAGC"
      }
    }
  }
}
```

---

## Next Steps

### Phase 1: Config Integration (Priority 1, ~2 hours)
- [ ] Add `snapshot_validation` section to config.json
- [ ] Create JSON schema validation
- [ ] Load all parameters from config

### Phase 2: Modularization (Priority 2, ~3 hours)
- [ ] Implement `PCRSimulator` class
- [ ] Implement `DigestSimulator` class
- [ ] Implement `SnapshotExtensionSimulator` class
- [ ] Implement `SnapshotValidator` orchestrator

### Phase 3: Production Code (Priority 3, ~3 hours)
- [ ] Create `muc_one_up/analysis/snapshot_validator.py`
- [ ] Add error handling and logging
- [ ] Write unit tests (config-based)
- [ ] Write integration tests

### Phase 4: CLI Integration (Priority 4, ~2 hours)
- [ ] Add `muconeup analyze snapshot-validate` command
- [ ] JSON output format
- [ ] Update documentation

**Total Estimated Time**: ~10 hours

---

## Usage (After Implementation)

```bash
# Validate dupC mutation
muconeup analyze snapshot-validate \
    output/dupC/dupC.001.mut.simulated.fa \
    --mutation dupC \
    --config config.json \
    --output results.json
```

**Expected Output**:
```json
{
  "haplotype_1": {
    "mutation_detected": true,
    "fluorescence": "Black (dTAMRA)",
    "interpretation": "8C mutation present"
  },
  "haplotype_2": {
    "mutation_detected": false,
    "interpretation": "Normal (7C only)"
  }
}
```

---

## Files Structure

```
plan/issue_15_snapshot_validation/
├── README.md                    # ← This file
├── PLAN.md                      # Complete implementation plan
├── STATUS.md                    # Status, achievements, validation results
└── tests/                       # Prototype tests (mechanism validated)
    ├── README.md
    ├── test_complete_snapshot_workflow.py  # ✓✓ Main test (ALL PASSING)
    ├── test_digest_selection_logic.py
    ├── test_protocol_validation.py
    ├── test_snapshot_validation_real.py
    └── test_all_pcr_tools.py
```

---

## Key References

### Primer Sequences (Validated)
```python
PCR_PRIMER_F = "GGCCGGCCCCGGGCTCCACC"       # 20bp
PCR_PRIMER_R = "TCCGGGGCCGAGGTGACA"         # 18bp (RC)
PCR_FLANKED_7C = "GCCCCCCCAGCCCACGG"        # 17bp (normal)
PCR_FLANKED_8C = "GCCCCCCCCAGCCCACGG"       # 18bp (mutant)
```

### Amplicon Structure
```
[20bp Forward] - [17-18bp Flanked] - [18bp Reverse]
Total: 55bp (7C) or 56bp (8C)
```

### Fluorophore Map
```
A → Green (dR6G)
C → Black (dTAMRA)  ← Expected for 8C mutation
G → Blue (dR110)
T → Red (dROX)
```

---

## Important Notes

1. **User-provided sequence `GCCCCCCCAGCCCACGG` is 7C NORMAL**, not 8C
2. For 8C mutant, add one C: `GCCCCCCCCAGCCCACGG`
3. Protocol's "tagged with 21bp" refers to amplicon structure, not primer length
4. PCR non-specificity is by design (primers bind to all X repeats)
5. Digest provides selection, not PCR specificity
6. Must test PCR amplicon for MwoI sites, not whole genome

---

## Validation Metrics

All metrics achieved in prototype tests:

- ✅ **MwoI site detection**: 100% accuracy (7C has site, 8C doesn't)
- ✅ **Digest selection**: 100% accuracy (7C digested, 8C survives)
- ✅ **Mutation detection**: 100% accuracy (8C correctly identified)
- ✅ **False positive rate**: 0% (normal samples show no 8C)
- ✅ **Performance**: < 1 second for complete workflow

---

**Status**: ✅ Planning Complete - ⏭️ Config Integration Required
**Estimated Time to Production**: ~10 hours
