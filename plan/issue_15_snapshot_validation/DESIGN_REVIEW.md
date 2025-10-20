# Issue #15: SNaPshot Validator - Design Review

**Date**: 2025-10-20
**Reviewer**: Architecture Review
**Status**: ⚠️ NEEDS REVISION

---

## Design Principles Assessment

### ✅ SOLID Principles

**Single Responsibility Principle**: ✅ GOOD
- Each method has clear, single purpose
- PCR simulation separate from digest
- Digest separate from SNaPshot extension

**Open/Closed Principle**: ✅ GOOD
- Extensible for other mutations
- Can add new enzymes without modifying core logic

**Liskov Substitution Principle**: ✅ GOOD
- Return types consistent
- Dictionaries allow flexibility

**Interface Segregation Principle**: ✅ GOOD
- Methods have focused interfaces
- No fat interfaces

**Dependency Inversion Principle**: ⚠️ NEEDS WORK
- Currently depends on concrete implementations (hardcoded values)
- Should depend on abstractions (config)

---

## Critical Issues Found

### 🔴 Issue #1: Hardcoded Primer Sequences

**Current Approach** (WRONG):
```python
class SnapshotPrimers:
    PCR_PRIMER_F = "GGCCGGCCCCGGGCTCCACC"  # HARDCODED!
    PCR_PRIMER_R = "TCCGGGGCCGAGGTGACA"
    SNAPSHOT_PRIMER_7C = "CGGGCTCCACCGCCCCCCC"
```

**Problem**:
- Violates DRY (primers repeated in multiple places)
- Not configurable (can't change primers without code change)
- Doesn't integrate with existing config.json structure

**Solution** - Use Config:
```python
# config.json
{
  "snapshot_validation": {
    "dupC": {
      "pcr_primers": {
        "forward": "GGCCGGCCCCGGGCTCCACC",
        "reverse": "TGTCACCTCGGCCCCGGA",
        "reverse_complement": true
      },
      "snapshot_primers": {
        "primer_7c": "CGGGCTCCACCGCCCCCCC",
        "primer_repeat_r": "TGTCACCTCGGCCCCGGA"
      },
      "restriction_enzyme": "MwoI",
      "amplicon_size_range": [50, 65],
      "mutation_patterns": {
        "mutant": "GCCCCCCCCAGC",
        "normal": "GCCCCCCCAGC"
      }
    }
  }
}
```

```python
# Code
class SnapshotValidator:
    def __init__(self, config: Dict[str, Any], mutation_name: str):
        self.config = config
        self.mutation_config = config["snapshot_validation"][mutation_name]

        # Load from config
        self.pcr_forward = self.mutation_config["pcr_primers"]["forward"]
        self.pcr_reverse = self.mutation_config["pcr_primers"]["reverse"]
        self.enzyme = self.mutation_config["restriction_enzyme"]
        self.size_range = self.mutation_config["amplicon_size_range"]
```

---

### 🔴 Issue #2: Hardcoded Mutation Patterns

**Current Approach** (WRONG):
```python
pattern_8c = "GCCCCCCCCAGC"  # HARDCODED!
pattern_7c = "GCCCCCCCAGC"
```

**Problem**:
- Can't validate other mutations without code changes
- Patterns should come from mutation definitions in config

**Solution** - Use Existing Mutation Config:
```python
# config.json already has:
{
  "mutations": {
    "dupC": {
      "allowed_repeats": ["X"],
      "changes": [
        {
          "operation": "insert",
          "position": 0,
          "sequence": "C"
        }
      ]
    }
  }
}

# Derive patterns from mutation definition
def get_mutation_patterns(mutation_def):
    """Derive expected patterns from mutation definition."""
    # For dupC (insert C), we expect 8C vs 7C
    if mutation_def["changes"][0]["operation"] == "insert":
        base = mutation_def["changes"][0]["sequence"]
        # Pattern changes based on mutation
        return {
            "mutant": calculate_mutant_pattern(mutation_def),
            "normal": calculate_normal_pattern(mutation_def)
        }
```

---

### 🔴 Issue #3: Hardcoded Amplicon Size Filter

**Current Approach** (WRONG):
```python
expected_min = 50  # HARDCODED!
expected_max = 65
```

**Problem**:
- Magic numbers
- Different mutations may have different expected sizes
- Should be calculated or configured

**Solution**:
```python
# Calculate from primers + flanked region
primer_f_len = len(self.pcr_forward)
primer_r_len = len(self.pcr_reverse)
flanked_min = 15  # From config
flanked_max = 20  # From config

expected_min = primer_f_len + flanked_min + primer_r_len - 5  # tolerance
expected_max = primer_f_len + flanked_max + primer_r_len + 5
```

---

### 🔴 Issue #4: Hardcoded Enzyme Name

**Current Approach** (WRONG):
```python
def simulate_digest_selection(amplicons, enzyme_name="MwoI"):  # HARDCODED!
```

**Problem**:
- Different mutations may use different enzymes
- Should come from config

**Solution**:
```python
# From config
enzyme_name = self.mutation_config["restriction_enzyme"]

# Or even better, support multiple enzymes
enzymes = self.mutation_config.get("restriction_enzymes", ["MwoI"])
```

---

### 🟡 Issue #5: Modularization Can Be Improved

**Current Structure**:
```
SnapshotValidator
├── validate_complete_workflow()  # Does everything
├── simulate_pcr_amplification()
├── simulate_digest_selection()
└── simulate_snapshot_extension()
```

**Better Structure** (KISS + Modular):
```
# Separate concerns into focused modules

class PCRSimulator:
    """Handles PCR amplification only."""
    def amplify(template, primers, config) -> List[Amplicon]

class DigestSimulator:
    """Handles restriction digest only."""
    def digest(amplicons, enzyme, config) -> List[Amplicon]

class SnapshotExtensionSimulator:
    """Handles SNaPshot extension only."""
    def extend(amplicon, primer, config) -> ExtensionResult

class SnapshotValidator:
    """Orchestrates the complete workflow."""
    def __init__(self, config, mutation_name):
        self.pcr = PCRSimulator(config)
        self.digest = DigestSimulator(config)
        self.extension = SnapshotExtensionSimulator(config)

    def validate(self, template):
        amplicons = self.pcr.amplify(template)
        survivors = self.digest.filter(amplicons)
        results = self.extension.extend_all(survivors)
        return self.interpret(results)
```

**Benefits**:
- Each component testable independently
- Can swap implementations
- Easier to debug
- Follows Single Responsibility Principle

---

## Regression/Bug Prevention

### ✅ Good Practices Already In Place:

1. **Type Hints**: All methods have proper type hints ✅
2. **Error Handling**: Methods return structured results ✅
3. **Testing**: Comprehensive test suite ✅
4. **Documentation**: Well-documented ✅

### ⚠️ Potential Bugs:

1. **Size Filtering Too Aggressive**:
   ```python
   if expected_min <= len(amplicon_seq) <= expected_max:
   ```
   - Could miss valid amplicons if size calculation is off
   - Should log filtered amplicons for debugging

2. **No Validation of Primer RC**:
   ```python
   PCR_PRIMER_R = "TCCGGGGCCGAGGTGACA"  # Assumes already RC'd
   ```
   - Should validate or auto-RC based on config flag

3. **No Handling of Ambiguous Bases**:
   - What if template contains N, R, Y, etc.?
   - Should either skip or handle gracefully

---

## Recommended Revisions

### Priority 1: Config Integration

1. Add `snapshot_validation` section to config.json
2. Load all parameters from config
3. Remove all hardcoded primers, patterns, sizes

### Priority 2: Modularization

1. Split into separate classes:
   - `PCRSimulator`
   - `DigestSimulator`
   - `SnapshotExtensionSimulator`
   - `SnapshotValidator` (orchestrator)

2. Each class independently testable

### Priority 3: DRY Improvements

1. Mutation patterns derived from mutation definitions
2. Expected sizes calculated from primer lengths
3. Enzyme names from config

### Priority 4: Error Handling

1. Add validation for:
   - Config completeness
   - Primer validity (length, composition)
   - Template quality
   - Enzyme availability

2. Graceful degradation:
   - Log warnings for filtered amplicons
   - Return partial results if some steps fail

---

## Revised Config Structure

```json
{
  "snapshot_validation": {
    "dupC": {
      "description": "8C mutation validation",

      "pcr": {
        "forward_primer": "GGCCGGCCCCGGGCTCCACC",
        "reverse_primer": "TGTCACCTCGGCCCCGGA",
        "reverse_needs_rc": true,
        "max_products": 100,
        "size_tolerance": 5
      },

      "digest": {
        "enzyme": "MwoI",
        "recognition_site": "GCNNNNNNNGC",
        "expected_survivors": "mutant_only"
      },

      "snapshot": {
        "primers": {
          "primer_7c": "CGGGCTCCACCGCCCCCCC",
          "primer_repeat_r": "TGTCACCTCGGCCCCGGA"
        },
        "fluorophore_map": {
          "A": {"color": "Green", "dye": "dR6G"},
          "C": {"color": "Black", "dye": "dTAMRA"},
          "G": {"color": "Blue", "dye": "dR110"},
          "T": {"color": "Red", "dye": "dROX"}
        }
      },

      "validation": {
        "mutant_pattern": "GCCCCCCCCAGC",
        "normal_pattern": "GCCCCCCCAGC",
        "expected_mutant_fluorescence": "Black"
      }
    }
  }
}
```

---

## Anti-Pattern Assessment

### ❌ Anti-Patterns Found:

1. **Magic Numbers**: Size filters (50, 65), max products (100)
2. **God Class**: `validate_complete_workflow()` does too much
3. **Hardcoded Configuration**: Primers, patterns, enzymes
4. **Tight Coupling**: Methods depend on hardcoded strings

### ✅ Good Patterns Used:

1. **Dependency Injection**: Methods accept parameters
2. **Return Dictionaries**: Flexible, extensible results
3. **Functional Decomposition**: Separate methods for steps
4. **Type Safety**: Type hints throughout

---

## Final Recommendation

**Status**: ⚠️ **REVISE BEFORE IMPLEMENTATION**

**Required Changes**:
1. ✅ Add config integration (Priority 1)
2. ✅ Remove all hardcoded values
3. ✅ Split into modular classes
4. ✅ Add error handling
5. ✅ Update tests to use config

**Timeline**:
- Config integration: 1-2 hours
- Modularization: 2-3 hours
- Testing updates: 1 hour
- **Total: ~5 hours**

**After Revision**: Ready for implementation ✅

---

## Action Items

- [ ] Create `snapshot_validation` config schema
- [ ] Add config validation tests
- [ ] Refactor into modular classes
- [ ] Update test script to use config
- [ ] Add integration tests
- [ ] Document config options
- [ ] Update IMPLEMENTATION_PLAN.md with revisions

---

**Reviewer**: Architecture Team
**Status**: Needs Revision
**Next Review**: After config integration complete
