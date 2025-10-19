# Issue #15: SNaPshot Assay Validation (CORRECTED)

## CORRECTION NOTICE
**Original plan reimplemented `reverse_complement()` which already exists in `translate.py`**

Existing function (USE THIS):
- `muc_one_up/translate.py` contains `reverse_complement(seq: str) -> str`

## Problem Statement
[NO CHANGES - scientifically valid as originally written]

## Implementation

[Most of the implementation is correct, with ONE fix]

### Core Validation Module (`muc_one_up/snapshot_validator.py`)
```python
from typing import Tuple, Optional
from dataclasses import dataclass

# ✅ CORRECTED: Import existing function instead of reimplementing
from .translate import reverse_complement  # ← USE EXISTING FUNCTION

@dataclass
class SnapshotAssayResult:
    """Result of SNaPshot assay validation."""
    accessible: bool
    mutation_site_present: bool
    pcr_amplifiable: bool
    primer_anneals: bool
    pcr_product_size: Optional[int]
    details: str

def validate_snapshot_accessibility(
    haplotype_sequence: str,
    mutation_label: str,
    assay_params: dict,
) -> SnapshotAssayResult:
    """
    Validate if mutation is detectable by SNaPshot assay.

    Returns:
        SnapshotAssayResult with detailed validation status
    """
    # Step 1: Check mutant restriction site present
    mut_site = assay_params["restriction_site_mut"]
    if mut_site not in haplotype_sequence:
        return SnapshotAssayResult(
            accessible=False,
            mutation_site_present=False,
            pcr_amplifiable=False,
            primer_anneals=False,
            pcr_product_size=None,
            details=f"Mutant restriction site {mut_site} not found"
        )

    # Step 2: Verify PCR primers flank mutation site
    fwd_primer = assay_params["pcr_fwd_primer"]
    rev_primer = assay_params["pcr_rev_primer"]

    fwd_pos = haplotype_sequence.find(fwd_primer)
    # ✅ USE EXISTING reverse_complement from translate.py
    rev_pos = haplotype_sequence.find(reverse_complement(rev_primer))
    mut_pos = haplotype_sequence.find(mut_site)

    if fwd_pos == -1 or rev_pos == -1:
        return SnapshotAssayResult(
            accessible=False,
            mutation_site_present=True,
            pcr_amplifiable=False,
            primer_anneals=False,
            pcr_product_size=None,
            details="PCR primers do not anneal to haplotype"
        )

    # Step 3: Validate primer orientation and product size
    if not (fwd_pos < mut_pos < rev_pos):
        return SnapshotAssayResult(
            accessible=False,
            mutation_site_present=True,
            pcr_amplifiable=False,
            primer_anneals=True,
            pcr_product_size=None,
            details="PCR primers do not flank mutation site"
        )

    pcr_product_size = rev_pos - fwd_pos
    min_size = assay_params.get("min_pcr_product_size", 100)
    max_size = assay_params.get("max_pcr_product_size", 5000)

    if not min_size <= pcr_product_size <= max_size:
        return SnapshotAssayResult(
            accessible=False,
            mutation_site_present=True,
            pcr_amplifiable=False,
            primer_anneals=True,
            pcr_product_size=pcr_product_size,
            details=f"PCR product size {pcr_product_size}bp outside valid range"
        )

    # Step 4: Verify SNaPshot extension primer
    snap_primer = assay_params["snapshot_primer"]
    snap_pos = haplotype_sequence.find(snap_primer)

    if snap_pos == -1 or snap_pos >= mut_pos:
        return SnapshotAssayResult(
            accessible=False,
            mutation_site_present=True,
            pcr_amplifiable=True,
            primer_anneals=True,
            pcr_product_size=pcr_product_size,
            details="SNaPshot primer does not anneal upstream of mutation"
        )

    # All checks passed
    return SnapshotAssayResult(
        accessible=True,
        mutation_site_present=True,
        pcr_amplifiable=True,
        primer_anneals=True,
        pcr_product_size=pcr_product_size,
        details=f"Mutation detectable via SNaPshot (PCR product: {pcr_product_size}bp)"
    )

# ❌ REMOVED: Don't reimplement reverse_complement - already exists in translate.py
```

[All other parts of the implementation remain unchanged - they were correct]

## Testing

### Unit Tests
```python
# tests/test_snapshot_validator.py

def test_uses_existing_reverse_complement():
    """Verify we import reverse_complement from translate.py."""
    from muc_one_up.snapshot_validator import reverse_complement
    from muc_one_up.translate import reverse_complement as translate_rc

    # Should be the SAME function object
    assert reverse_complement is translate_rc

def test_dupC_detectable():
    """dupC mutation should disrupt MwoI site and be detectable."""
    haplotype = "GGCCGGCCCCGGGCTCCACC" + \
                "GCCCCCCCCAGC" + \
                "TCCGGGGCCGAGGTGACA"
    assay_params = {
        "restriction_site_mut": "GCCCCCCCCAGC",
        "pcr_fwd_primer": "GGCCGGCCCCGGGCTCCACC",
        "pcr_rev_primer": "TCCGGGGCCGAGGTGACA",  # Forward strand - will be rev-comp'd
        "snapshot_primer": "CGGGCTCCACCGCCCCCCC",
        "min_pcr_product_size": 100,
        "max_pcr_product_size": 5000,
    }

    result = validate_snapshot_accessibility(haplotype, "dupC", assay_params)
    assert result.accessible is True
    assert result.pcr_amplifiable is True

# ... rest of tests unchanged
```

## Files Modified

### MODIFIED (1 function import):
- `muc_one_up/snapshot_validator.py` (import reverse_complement from translate.py)

### DO NOT CREATE:
- ❌ `reverse_complement()` function - Already exists in translate.py

## Compliance with Programming Principles

✅ **DRY**: Reuses existing reverse_complement function

✅ **KISS**: Simple function import

✅ **SOLID**: Single Responsibility maintained

✅ **Modular**: Imports from appropriate module
