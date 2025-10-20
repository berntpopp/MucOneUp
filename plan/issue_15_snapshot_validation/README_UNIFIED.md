# Issue #15: SNaPshot Assay Validation - UNIFIED DOCUMENTATION

**Date**: 2025-10-19
**Status**: 🔴 **BLOCKED - Awaiting Primer Clarification**
**Location**: `/plan/issue_15_snapshot_validation/`

---

## Executive Summary

This issue implements in-silico validation of the SNaPshot assay for detecting MUC1 VNTR mutations, specifically the dupC (8C) mutation.

**Current Status**:
- ✅ Protocol fully analyzed (8-day workflow)
- ✅ Real primer sequences received from user
- ✅ pydna PCR simulation tested and working
- ✅ dupC mutation confirmed in test data (1× 8C in mutant, 0× in normal)
- ✅ MwoI restriction digest simulation working
- ❌ **BLOCKED**: Primer orientation problem (both primers same strand)

**Blocker**: User-provided primers both bind to the forward strand. Need clarification on:
1. Are these complete primer sequences or just binding regions?
2. Where do the tags attach?
3. Which primer should be reverse complement?

---

## Complete Protocol Summary

**Source**: `Protocol MUC1-Englisch_Arif.docx`

### 8-Day Workflow

```
Day 1: Genomic DNA → MwoI digest (3 rounds, 2-3h each @ 60°C)
       ↓
Day 2: PCR amplification with MUC1-Repeat F/R primers (45 cycles)
       PCR conditions: 95°C 5min → (94°C 30s, 67°C 30s, 72°C 30s) × 45 → 72°C 10min
       ↓
Day 3: PCR product → MwoI digest (3 rounds)
       ↓
Day 4: AMPure purification (1.8 vol, 2× wash 70% EtOH, elute 15μl)
       ↓
Day 5: ExoSAP digestion (37°C 20min, 80°C 15min)
       ↓
Day 6: SNaPshot reaction (2 primers: 7C + repeat R)
       Cycling: 94°C 2min → (94°C 10s, 52°C 5s, 60°C 5s) × 55 → 60°C 1min
       ↓
Day 7: SAP digestion (37°C 20min, 80°C 15min)
       ↓
Day 8: Capillary electrophoresis (3500DX sequencer, SnapSHOT program)
       ↓
Day 9: GenMapper analysis
```

### The Mutation: dupC (7C → 8C)

**Wild-Type (7 Cs)**:
```
...GCCCCCC-CAGC...
   ^^^^^^^
   7 C's
```

**Mutant (8 Cs - dupC)**:
```
...GCCCCCCCCAGC...
   ^^^^^^^^
   8 C's (extra C inserted)
```

### MwoI Enzyme

**Recognition Site**: `GCNNNNNNNGC` (GC...7 bases...GC)

**Expected Behavior**:
- WT (7C): **CUTS** - DNA digested
- Mutant (8C): **DOES NOT CUT** - DNA intact

**Reality from Testing**:
- Both mutant and normal have **112 MwoI sites**
- dupC is a frameshift mutation, not site-disrupting
- Global MwoI site count unchanged

---

## Primer Sequences (from User)

### PCR Primers

```python
# User-provided sequences
PCR_PRIMER_1 = "GGCCGGCCCCGGGCTCCACC"  # 20bp
PCR_PRIMER_2 = "TGTCACCTCGGCCCCGGA"    # 18bp
PCR_TAG = "GCCCCCCCAGCCCACGG"          # 17bp (contains 8C: GCCCCCCCCAGC)
```

**Protocol Statement**:
> "The primers MUC1-Repeat F and R are located in 2 contiguous repeats flanking the 7C/8C and are tagged with a 21 bp sequence."

**Problem Identified**:
- Both primers found in X repeat (60bp)
- Primer 1 at position 31: `GGCCGGCCCCGGGCTCCACC`
- Primer 2 at position 8: `TGTCACCTCGGCCCCGGA`
- **Both in FORWARD orientation on same strand** ←←
- Cannot create PCR product (need forward → and reverse ←)

### SNaPshot Extension Primers

```python
SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"     # 19bp
SNAPSHOT_PRIMER_REPEAT_R = "TGTCACCTCGGCCCCGGA"  # 18bp (same as PCR_PRIMER_2)
SNAPSHOT_TAG = "GCCCCACGG"                       # 9bp
```

**Expected Results (from Protocol)**:

**Primer 7C** produces 28 bp peak (19bp + 1pb fluoro):
- **Green** (nucleotide A): Incomplete digestion, no 8C
- **Black** (nucleotide C): **8C mutation present** ✅
- **No peak**: Complete digestion, no 8C

**Primer repeat R** produces 45 bp peak (39pb + 1pb fluoro):
- **Black** (nucleotide C): VNTR amplified (8C or incomplete digestion)

---

## Test Data

### Files Used

```
input/dupC/dupC.001.mut.simulated.fa   (13,601 bp)
output/dupC/dupC.001.normal.simulated.fa (13,600 bp)
```

### Mutation Verification

```
Pattern: GCCCCCCCCAGC (8C)
  Mutant: 1 occurrence  ✓
  Normal: 0 occurrences ✓

Pattern: GCCCCCCCAGC (7C)
  Mutant: 34 occurrences
  Normal: 35 occurrences
```

**dupC mutation confirmed!**

### X Repeat Locations

```
Exact X repeat count:
  Mutant: 26 occurrences
  Normal: 27 occurrences

Primer binding sites (without tags):
  Primer 1 (forward): 49 sites
  Primer 2 (forward): 54 sites
  Primer 2 (reverse complement): 0 sites ← PROBLEM!
```

---

## Testing Results

### Test Files Created

All in `/tests/`:
1. `test_all_pcr_tools.py` - Comprehensive PCR tool comparison
2. `test_snapshot_validation_real.py` - Initial validation tests
3. `test_protocol_validation.py` - Protocol-based primer design tests
4. `test_complete_snapshot_workflow.py` - Full workflow with user primers

### Tools Validated

| Tool | Status | Notes |
|------|--------|-------|
| **pydna** | ✅ Working | v5.5.3, PCR simulation accurate |
| **primer3-py** | ✅ Working | v2.2.0, Tm calculations verified |
| **BioPython** | ✅ Working | v1.85, MwoI digest functional |
| ispcr | ⚠️ Works but awkward | Unmaintained, string output |
| PyPCRtool | ❌ Broken | Import fails |

**Recommendation**: Use pydna (ONLY) - no fallbacks needed

### What Works ✅

1. **pydna PCR simulation**
   - Correctly detects non-specific primers (53 forward, 34 reverse sites)
   - Thermodynamically accurate
   - Fast: 6ms for 13.6kb sequence
   - Error handling: Refuses ambiguous results

2. **MwoI restriction digest**
   - Accurately finds all 112 sites in both sequences
   - Simulates fragmentation correctly
   - Largest fragment: 967 bp

3. **dupC mutation detection**
   - 8C pattern found in mutant only
   - 7C pattern in both (expected - not all are mutation sites)
   - Clear discrimination

4. **SNaPshot extension logic**
   - Primer binding site detection working
   - Next-base determination correct
   - Fluorophore mapping implemented

### What Doesn't Work ❌

1. **Primer orientation**
   - Both primers bind forward strand
   - No reverse complement matches
   - Cannot create PCR product
   - **BLOCKED until resolved**

2. **MwoI site discrimination**
   - Both sequences have 112 sites (same count)
   - dupC frameshift doesn't change global site count
   - Protocol assumes site disruption (may not apply to frameshift)

3. **Multiplex assumption**
   - Initially thought multiplex = multiple mutations
   - Actually: Single mutation (dupC) with 2 validation primers
   - Not a panel for different mutations

---

## Critical Findings

### Finding #1: Repeat-Based Primers Are Not Specific

**Without tags**, primers from X repeat bind to ALL ~27 X repeats:

```
VNTR: X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X-X
      ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑   ↑
      Every X = potential binding site

Result: 49 forward × 54 reverse = 2,646 possible products!
```

**pydna correctly rejects this** as non-specific.

**Solution**: Tags must provide specificity
- Bind to constant regions OR
- Bind to specific repeat junctions OR
- Create unique combined sequence

### Finding #2: dupC Is Frameshift, Not Site-Disrupting

**Two mutation types**:

**Type A (Site-Disrupting)**:
- Mutation creates/destroys restriction site
- Example: Point mutation in GCNNNNNNNGC
- Protocol works: Digest → PCR → SNaPshot

**Type B (Frameshift - like dupC)**:
- Inserts/deletes bases, shifts reading frame
- May not change MwoI site count globally
- Protocol may not work as expected
- Need direct detection without digest pre-selection

**dupC appears to be Type B**

### Finding #3: Primer Orientation Problem

**Current situation**:
```
X repeat (60bp):
Position 8:  TGTCACCTCGGCCCCGGA  ←── Primer 2
Position 31: GGCCGGCCCCGGGCTCCACC ←── Primer 1

Both reading LEFT-TO-RIGHT (forward)
```

**For PCR, need**:
```
Forward:  5'→3' on top strand
Reverse:  5'→3' on bottom strand (binds as RC to top)

Creates: Forward extends →
         Reverse extends ←
         = PCR product
```

**But we have**: Both extend →→ (no product possible)

---

## Current Blocker 🔴

### The Question

**How are the primers supposed to work?**

User provided:
- Primer 1: `GGCCGGCCCCGGGCTCCACC`
- Primer 2: `TGTCACCTCGGCCCCGGA`
- Tag: `GCCCCCCCAGCCCACGG` (17bp, contains 8C)

**But**: Both primers in same orientation!

### Possible Solutions

**Option A**: Primers need tags
```
Full_Forward = [Tag_A] + GGCCGGCCCCGGGCTCCACC
Full_Reverse = [Tag_B] + TGTCACCTCGGCCCCGGA

Where: Tag_A binds upstream constant region
       Tag_B binds downstream constant region
```

**Option B**: One primer should be RC
```
Use: GGCCGGCCCCGGGCTCCACC (forward)
And: TCCGGGGCCGAGGTGACA (RC of primer 2)
```

**Option C**: These are binding regions only
```
Actual primers have additional sequence
Need complete lab primer sequences
```

**Option D**: Different genomic locations
```
Maybe not both from X repeat?
One from constant region?
```

### What We Need

**Please clarify ONE of these**:

1. ✅ Complete primer sequences including tags
2. ✅ Which strand each primer binds to
3. ✅ Genomic positions where primers bind
4. ✅ Expected PCR product size

---

## Implementation Plan (Once Unblocked)

### Module: `muc_one_up/analysis/snapshot_validator.py`

```python
class SnapshotValidator:
    """Complete SNaPshot workflow validation."""

    def validate_complete_workflow(
        self,
        haplotype_seq: str,
        mutation_type: str,
        pcr_primers: Dict[str, str],
        extension_primers: Dict[str, str]
    ) -> SnapshotValidationResult:
        """
        Complete workflow:
        1. Genomic DNA digest (MwoI)
        2. PCR amplification
        3. PCR product digest
        4. SNaPshot extension
        5. Result interpretation
        """
        # Step 1: Digest genomic DNA
        digest_result = self.simulate_mwoi_digest(haplotype_seq)

        # Step 2: PCR on digested DNA
        pcr_result = self.simulate_pcr(
            haplotype_seq,
            pcr_primers["forward"],
            pcr_primers["reverse"]
        )

        # Step 3: Digest PCR product
        if pcr_result.amplified:
            pcr_digest = self.simulate_mwoi_digest(pcr_result.product_sequence)
        else:
            return SnapshotValidationResult(
                accessible=False,
                reason="PCR failed - primers don't amplify"
            )

        # Step 4: SNaPshot extension
        snapshot_results = {}
        for primer_name, primer_seq in extension_primers.items():
            snapshot_results[primer_name] = self.simulate_extension(
                pcr_result.product_sequence,
                primer_seq
            )

        # Step 5: Interpret
        return self.interpret_results(snapshot_results)
```

### CLI Integration

```bash
# Once unblocked
muconeup analyze snapshot-validate \
    dupC.001.mut.simulated.fa \
    --mutation dupC \
    --primers primers.json \
    --output results.json
```

### Configuration Schema

```json
{
  "snapshot_assay": {
    "mutations": {
      "dupC": {
        "type": "frameshift",
        "pcr_primers": {
          "forward": "FULL_SEQUENCE_WITH_TAG",
          "reverse": "FULL_SEQUENCE_WITH_TAG"
        },
        "extension_primers": {
          "primer_7C": "CGGGCTCCACCGCCCCCCC",
          "primer_repeat_R": "TGTCACCTCGGCCCCGGA"
        },
        "restriction_enzyme": "MwoI",
        "expected_results": {
          "mutant_8C": {
            "primer_7C": "Black (C)",
            "primer_repeat_R": "Black (C)"
          },
          "normal_7C": {
            "primer_7C": "No peak (complete digest)",
            "primer_repeat_R": "No peak"
          }
        }
      }
    }
  }
}
```

---

## File Inventory

### Planning Documents (this folder)

| File | Purpose | Status |
|------|---------|--------|
| **README_UNIFIED.md** | This file - unified documentation | ✅ Current |
| Protocol MUC1-Englisch_Arif.docx | Lab protocol source | ✅ Reference |
| issue_15_snapshot_assay_validation.md | Original issue spec | 📦 Archived |
| issue_15_critical_assessment.md | Initial assessment | 📦 Archived |
| issue_15_revised_comprehensive_assessment.md | Corrected assessment | 📦 Archived |
| issue_15_primer3_poc.md | Early POC | 📦 Archived |
| issue_15_pydna_poc_revised.md | Corrected POC | 📦 Archived |
| issue_15_test_results_and_findings.md | Test results | 📦 Archived |
| issue_15_insilico_pcr_tools_final_analysis.md | Tool comparison | ✅ Reference |
| issue_15_pydna_only_implementation.md | Implementation plan | ✅ Reference |
| issue_15_multiplex_pcr_correction.md | Multiplex clarification | ✅ Reference |
| issue_15_protocol_analysis_and_correction.md | Protocol deep-dive | ✅ Reference |
| issue_15_critical_finding_primer_specificity.md | Specificity issue | ✅ Reference |
| issue_15_primer_orientation_problem.md | Current blocker | 🔴 **ACTIVE** |

### Test Scripts (/tests/)

| File | Purpose | Status |
|------|---------|--------|
| test_all_pcr_tools.py | PCR tool comparison | ✅ Complete |
| test_snapshot_validation_real.py | Initial tests | ✅ Complete |
| test_protocol_validation.py | Protocol validation | ✅ Complete |
| test_complete_snapshot_workflow.py | Full workflow | ⚠️ Blocked |

### Test Data

```
output/dupC/dupC.001.mut.simulated.fa    (13,601 bp) - Mutant with 8C
output/dupC/dupC.001.normal.simulated.fa (13,600 bp) - Normal with 7C
```

---

## Timeline

| Date | Event |
|------|-------|
| 2025-10-19 | Issue #15 analysis started |
| 2025-10-19 | Protocol received and analyzed |
| 2025-10-19 | pydna testing completed |
| 2025-10-19 | Real primer sequences received |
| 2025-10-19 | **BLOCKED: Primer orientation problem discovered** |
| **Pending** | **Awaiting primer clarification** |
| **TBD** | Complete implementation |
| **TBD** | Create PR |

---

## Next Steps

### Immediate (User Action Required)

**Clarify primer sequences**:
1. Provide complete primer sequences with tags, OR
2. Specify which primer should be reverse complement, OR
3. Provide genomic binding positions for each primer

### After Unblock (Implementation)

1. Test PCR with corrected primers
2. Validate complete workflow on dupC data
3. Implement `muc_one_up/analysis/snapshot_validator.py`
4. Create CLI command `snapshot-validate`
5. Write comprehensive tests
6. Update documentation
7. Create pull request

---

## Key Contacts

- **User**: bernt (GitHub: berntpopp/MucOneUp)
- **Protocol Author**: Arif (see Protocol MUC1-Englisch_Arif.docx)
- **Issue**: https://github.com/berntpopp/MucOneUp/issues/15

---

## Summary for Quick Reference

**What is this?**: In-silico SNaPshot validation for MUC1 VNTR dupC mutation

**What works?**: pydna PCR simulation, MwoI digest, dupC detection in test data

**What's blocked?**: Primer orientation - both primers same strand, need clarification

**What's needed?**: Complete primer sequences or strand orientation clarification

**When complete?**: Can implement full workflow once primers resolved

---

**Last Updated**: 2025-10-19
**Document Version**: 1.0
**Status**: 🔴 BLOCKED - Awaiting Primer Clarification
