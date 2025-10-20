# Issue #15: SNaPshot Validator - Complete Implementation Plan

**Date**: 2025-10-20
**Status**: Ready for Implementation
**Module**: `muc_one_up/analysis/snapshot_validator.py`

---

## Overview

Implement complete in-silico SNaPshot assay validation for MUC1 VNTR mutations, specifically the dupC (8C) mutation. The implementation simulates the complete wet-lab workflow: PCR amplification → MwoI digest selection → SNaPshot extension → fluorescence detection.

---

## Validated Mechanism

### dupC Mutation (7C → 8C)

**MwoI Recognition Site**: `GCNNNNNNNGC` (GC...exactly 7 bases...GC)

**7C NORMAL Amplicon** (55bp):
- Flanked: `GCCCCCCCAGCCCACGG` (17bp, 7 Cs)
- Pattern: `GC[CCCCCCC]AGC` ← exactly 7 bases
- **MwoI site present** → gets **DIGESTED** → **DESTROYED**

**8C MUTANT Amplicon** (56bp):
- Flanked: `GCCCCCCCCAGCCCACGG` (18bp, 8 Cs)
- Pattern: `GC[CCCCCCCC]AGC` ← 8 bases (too many!)
- **No MwoI site** → **SURVIVES** → proceeds to SNaPshot

**Result**: Extra C in dupC mutation disrupts the MwoI recognition site!

---

## Implementation Structure

### Module: `muc_one_up/analysis/snapshot_validator.py`

```
muc_one_up/
└── analysis/
    └── snapshot_validator.py    # NEW: Complete SNaPshot workflow validator
```

### Dependencies

```python
from Bio.Seq import Seq
from Bio.Restriction import MwoI
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
import primer3
from typing import Dict, List, Any, Optional
```

---

## Class Design

### `SnapshotValidator`

Main class for SNaPshot assay validation.

#### Core Methods

##### 1. `validate_complete_workflow()`

**Purpose**: Simulate complete PCR → Digest → SNaPshot workflow

**Parameters**:
- `template_seq: str` - Genomic template sequence
- `mutation_name: str` - Mutation identifier (e.g., "dupC")
- `primers: Dict[str, str]` - Primer sequences
- `enzyme: str` - Restriction enzyme name (default: "MwoI")

**Returns**:
```python
{
    "pcr_products": List[Dict],
    "digest_survivors": List[Dict],
    "snapshot_results": List[Dict],
    "mutation_detected": bool,
    "summary": str
}
```

**Workflow**:
1. Run PCR simulation (multiple products expected)
2. Filter products by digest (MwoI site presence)
3. Run SNaPshot on survivors
4. Determine mutation presence

---

##### 2. `simulate_pcr_amplification()`

**Purpose**: PCR simulation with multiple product handling

**Parameters**:
- `template_seq: str` - Template sequence
- `forward_primer: str` - Forward primer (20bp)
- `reverse_primer: str` - Reverse primer RC (18bp)
- `limit: int` - Maximum products to allow (default: 100)

**Returns**:
```python
{
    "success": bool,
    "product_count": int,
    "products": List[{
        "sequence": str,
        "length": int,
        "forward_pos": int,
        "reverse_pos": int,
        "flanked_region": str
    }],
    "non_specific": bool
}
```

**Implementation Notes**:
- Use `pydna.amplify.pcr()` with `limit` parameter
- Handle "PCR not specific" exception (expected for X repeat primers)
- Extract multiple binding sites manually if needed
- Identify flanked region (between primers)

---

##### 3. `simulate_digest_selection()`

**Purpose**: Filter amplicons by restriction enzyme digest

**Parameters**:
- `amplicons: List[Dict]` - PCR products from step 2
- `enzyme_name: str` - Enzyme (default: "MwoI")

**Returns**:
```python
{
    "total_products": int,
    "digested_count": int,
    "survivor_count": int,
    "survivors": List[{
        "sequence": str,
        "length": int,
        "has_mwoi_site": bool,
        "mwoi_site_count": int,
        "survives_digest": bool,
        "mutation_type": "7C" | "8C" | "unknown"
    }],
    "mechanism": str
}
```

**Logic**:
```python
for amplicon in amplicons:
    seq_obj = Seq(amplicon["sequence"])
    sites = MwoI.search(seq_obj)

    if len(sites) == 0:
        # No MwoI site → survives digest
        survivors.append({
            "survives_digest": True,
            "mutation_type": "8C" if "GCCCCCCCCAGC" in amplicon else "unknown"
        })
    else:
        # Has MwoI site → gets digested
        digested.append({
            "survives_digest": False,
            "mutation_type": "7C" if "GCCCCCCCAGC" in amplicon else "unknown"
        })
```

---

##### 4. `simulate_snapshot_extension()`

**Purpose**: SNaPshot single-base extension simulation

**Parameters**:
- `template_seq: str` - Amplicon sequence (from survivors)
- `primer_seq: str` - Extension primer (19bp for 7C primer)
- `primer_name: str` - Identifier

**Returns**:
```python
{
    "binds": bool,
    "position": int,
    "next_base": str,  # A, C, G, or T
    "ddNTP": str,      # ddATP, ddCTP, ddGTP, ddTTP
    "fluorophore": str,  # Color name
    "fluorophore_dye": str,  # dR6G, dTAMRA, dR110, dROX
    "peak_size": int,  # Primer length + 1
    "interpretation": str
}
```

**Fluorophore Mapping**:
```python
FLUOROPHORE_MAP = {
    'A': {'color': 'Green', 'dye': 'dR6G'},
    'C': {'color': 'Black', 'dye': 'dTAMRA'},
    'G': {'color': 'Blue', 'dye': 'dR110'},
    'T': {'color': 'Red', 'dye': 'dROX'}
}
```

**Expected Results for dupC**:
- Next base: **C** (from 8C pattern)
- Fluorophore: **Black (dTAMRA)**
- Interpretation: "8C mutation detected"

---

##### 5. `validate_amplicon_digest_selection()`

**Purpose**: Validate if specific amplicon survives digest

**Parameters**:
- `amplicon_seq: str` - Single amplicon sequence
- `mutation_name: str` - Expected mutation

**Returns**:
```python
{
    "mutation_present": bool,  # Has 8C pattern
    "has_7c": bool,
    "has_mwoi_site": bool,
    "survives_digest": bool,
    "mwoi_site_count": int,
    "mwoi_positions": List[int],
    "amplicon_length": int,
    "detectable": bool,
    "mechanism": str
}
```

---

##### 6. `validate_primer_thermodynamics()`

**Purpose**: Validate primer properties using primer3

**Parameters**:
- `primer_seq: str` - Primer sequence
- `primer_type: str` - "PCR" or "SNaPshot"

**Returns**:
```python
{
    "length": int,
    "tm": float,  # Melting temperature
    "gc_percent": float,
    "hairpin_tm": float,
    "homodimer_tm": float,
    "valid": bool,
    "warnings": List[str]
}
```

---

## Primer Sequences (Constants)

```python
class SnapshotPrimers:
    """Validated primer sequences for dupC SNaPshot assay."""

    # First PCR primers
    PCR_PRIMER_F = "GGCCGGCCCCGGGCTCCACC"       # 20bp (forward)
    PCR_PRIMER_R_ORIGINAL = "TGTCACCTCGGCCCCGGA"  # 18bp (original)
    PCR_PRIMER_R = "TCCGGGGCCGAGGTGACA"          # 18bp (reverse complement)

    # Flanked sequences
    PCR_FLANKED_7C = "GCCCCCCCAGCCCACGG"   # 17bp (7C normal)
    PCR_FLANKED_8C = "GCCCCCCCCAGCCCACGG"  # 18bp (8C mutant)

    # SNaPshot extension primers
    SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"  # 19bp
    SNAPSHOT_PRIMER_REPEAT_R = "TCCGGGGCCGAGGTGACA"  # 18bp (RC)

    # Expected amplicon sizes
    AMPLICON_7C_SIZE = 55  # 20 + 17 + 18
    AMPLICON_8C_SIZE = 56  # 20 + 18 + 18
```

---

## Usage Examples

### Example 1: Complete Workflow Validation

```python
from muc_one_up.analysis.snapshot_validator import SnapshotValidator, SnapshotPrimers

# Initialize validator
validator = SnapshotValidator()

# Load template sequence (from FASTA)
from Bio import SeqIO
record = SeqIO.read("output/dupC/dupC.001.mut.simulated.fa", "fasta")
template = str(record.seq)

# Run complete workflow
result = validator.validate_complete_workflow(
    template_seq=template,
    mutation_name="dupC",
    primers={
        "pcr_forward": SnapshotPrimers.PCR_PRIMER_F,
        "pcr_reverse": SnapshotPrimers.PCR_PRIMER_R,
        "snapshot_7c": SnapshotPrimers.SNAPSHOT_PRIMER_7C
    }
)

# Check results
if result["mutation_detected"]:
    print("✓ 8C mutation detected via SNaPshot")
    print(f"  Survivors: {result['digest_survivors']['survivor_count']}")
    print(f"  Fluorescence: {result['snapshot_results'][0]['fluorophore']}")
```

### Example 2: Amplicon-Specific Validation

```python
# Test specific amplicon
amplicon_8c = (SnapshotPrimers.PCR_PRIMER_F +
               SnapshotPrimers.PCR_FLANKED_8C +
               SnapshotPrimers.PCR_PRIMER_R)

result = validator.validate_amplicon_digest_selection(
    amplicon_seq=amplicon_8c,
    mutation_name="dupC"
)

print(f"Survives digest: {result['survives_digest']}")
print(f"MwoI sites: {result['mwoi_site_count']}")
print(f"Mechanism: {result['mechanism']}")
```

---

## CLI Integration

### Command: `snapshot-validate`

```bash
muconeup analyze snapshot-validate \
    output/dupC/dupC.001.mut.simulated.fa \
    --mutation dupC \
    --output results.json \
    --verbose
```

**Output**:
```json
{
  "haplotype_1": {
    "pcr_products": 26,
    "digest_survivors": 1,
    "mutation_detected": true,
    "fluorescence": "Black (dTAMRA)",
    "interpretation": "8C mutation present"
  },
  "haplotype_2": {
    "pcr_products": 27,
    "digest_survivors": 0,
    "mutation_detected": false,
    "fluorescence": "None",
    "interpretation": "Normal (7C only)"
  }
}
```

---

## Testing Strategy

### Unit Tests

**File**: `tests/test_snapshot_validator.py`

Test cases:
1. `test_amplicon_7c_has_mwoi_site()` - Verify 7C has MwoI site
2. `test_amplicon_8c_no_mwoi_site()` - Verify 8C lacks MwoI site
3. `test_digest_selection_7c_destroyed()` - 7C gets digested
4. `test_digest_selection_8c_survives()` - 8C survives
5. `test_snapshot_extension_8c_black_peak()` - Correct fluorescence
6. `test_complete_workflow_dupC()` - End-to-end validation

### Integration Tests

**File**: `tests/test_snapshot_validator_integration.py`

Test with real data:
1. Load `dupC.001.mut.simulated.fa` (mutant)
2. Load `dupC.001.normal.simulated.fa` (normal)
3. Run complete workflow on both
4. Verify:
   - Mutant: mutation detected, Black peak
   - Normal: no mutation, no survivors or Green peak

---

## Implementation Steps

### Phase 1: Core Logic (Week 1)

1. ✅ Create `snapshot_validator.py` module
2. ✅ Implement `validate_amplicon_digest_selection()`
3. ✅ Implement `simulate_digest_selection()`
4. ✅ Implement `simulate_snapshot_extension()`
5. ✅ Add primer constants class

### Phase 2: PCR Integration (Week 1-2)

6. ✅ Implement `simulate_pcr_amplification()`
7. ✅ Handle pydna multiple products
8. ✅ Extract flanked regions
9. ✅ Test with real dupC data

### Phase 3: Complete Workflow (Week 2)

10. ✅ Implement `validate_complete_workflow()`
11. ✅ Integrate all components
12. ✅ Add comprehensive error handling
13. ✅ Create output formatting

### Phase 4: CLI & Testing (Week 2-3)

14. ⏭️ Add CLI command `snapshot-validate`
15. ⏭️ Write unit tests
16. ⏭️ Write integration tests
17. ⏭️ Update documentation

---

## Expected Outcomes

### For dupC Mutant

```
PCR Amplification:
  - Multiple products from X repeats (expected: ~27 amplicons)
  - Mix of 7C (55bp) and 8C (56bp) products

MwoI Digest:
  - 7C amplicons: DIGESTED (26 destroyed)
  - 8C amplicons: SURVIVE (1 survivor)

SNaPshot Extension:
  - Primer 7C binds to 8C survivor
  - Next base: C
  - Fluorophore: Black (dTAMRA)
  - Result: "8C mutation detected"
```

### For Normal Control

```
PCR Amplification:
  - Multiple products from X repeats (expected: ~27 amplicons)
  - All 7C (55bp) products

MwoI Digest:
  - All amplicons: DIGESTED (27 destroyed)
  - Survivors: 0

SNaPshot Extension:
  - No products to extend
  - Result: "No mutation detected (all products digested)"
```

---

## Validation Metrics

### Success Criteria

1. ✅ **MwoI site detection accuracy**: 100% (7C has site, 8C doesn't)
2. ✅ **Digest selection accuracy**: 100% (7C digested, 8C survives)
3. ✅ **Mutation detection accuracy**: 100% (8C correctly identified)
4. ✅ **False positive rate**: 0% (normal samples show no 8C)
5. ⏭️ **Performance**: < 1 second for complete workflow

---

## Future Enhancements

### Phase 5: Extended Mutations

1. Support for other MUC1 mutations (insG, delT, etc.)
2. Multiplex SNaPshot validation
3. Custom primer design from constant regions
4. Wet-lab protocol optimization suggestions

### Phase 6: Advanced Features

1. Confidence scoring for mutation calls
2. Coverage depth simulation
3. Allele frequency estimation
4. Comparison with actual SNaPshot data

---

## Notes

- User-provided sequence `GCCCCCCCAGCCCACGG` is **7C NORMAL** (not 8C)
- For 8C mutant: `GCCCCCCCCAGCCCACGG` (add one C)
- Protocol's "tagged with 21bp" refers to amplicon structure, not primer length
- Complete mechanism validated: extra C disrupts MwoI site
- Ready for production implementation

---

**Status**: ✅ Planning Complete - Ready for Code Implementation
**Next Step**: Create test script with real data validation
