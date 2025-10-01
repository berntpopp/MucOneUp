# Type Safety & Validation

**Priority:** 🟠 HIGH
**Estimated Effort:** 3-4 days
**Impact:** MEDIUM-HIGH - Catches bugs early, improves IDE support

## Current State

### Type Hints Coverage

| Status | Count | Files |
|--------|-------|-------|
| **With type hints** | 8 / 20+ | `config.py`, `samtools_wrapper.py`, etc. |
| **Without type hints** | 12+ | `cli.py`, `simulate.py`, most modules |
| **Partial hints** | ~4 | Some signatures typed, others not |

### Problems

1. **No type checking** - No mypy in CI/CD, type errors go unnoticed
2. **Incomplete hints** - Even typed files have gaps
3. **No input validation** - DNA sequences, file paths, etc. not validated
4. **IDE support limited** - Cannot autocomplete or catch errors

## SOLID Principles & Type Safety

### Interface Segregation Principle

Type hints create implicit interfaces:

```python
# ❌ No interface - unclear what's expected
def simulate(config, num, lengths, seed):
    pass

# ✅ Clear interface - types document expectations
def simulate(
    config: Dict[str, Any],
    num_haplotypes: int,
    fixed_lengths: Optional[List[int]] = None,
    seed: Optional[int] = None
) -> List[Tuple[str, List[str]]]:
    pass
```

### Dependency Inversion Principle

Type hints enable proper abstractions:

```python
# ❌ Concrete dependency
def process(samtools_wrapper):
    samtools_wrapper.run_command(...)

# ✅ Abstract interface
from typing import Protocol

class ToolWrapper(Protocol):
    def run_command(self, cmd: List[str]) -> str:
        ...

def process(tool: ToolWrapper):
    tool.run_command(...)
```

## Type Hints Strategy

### Phase 1: Core Types Module

Create `muc_one_up/types.py`:

```python
"""Type definitions for MucOneUp."""

from typing import Dict, List, Tuple, Any, Optional, Union
from pathlib import Path

# Haplotype representation
HaplotypeName = str
RepeatChain = List[str]  # List of repeat symbols like ["1", "2", "7", "8", "9"]
Haplotype = Tuple[HaplotypeName, RepeatChain]
HaplotypeList = List[Haplotype]

# Configuration types
ConfigDict = Dict[str, Any]
RepeatsDict = Dict[str, str]  # Symbol -> DNA sequence
ProbabilitiesDict = Dict[str, Dict[str, float]]  # State -> {Next State -> Probability}

# Mutation types
MutationName = str
MutationTargets = List[Tuple[int, int]]  # [(haplotype_idx, position_idx), ...]
MutatedUnits = Dict[str, str]  # Mutated repeat symbol -> sequence

# SNP types
SNPRecord = Tuple[int, int, str, str]  # (haplotype, position, ref, alt)
SNPList = List[SNPRecord]

# File paths
FilePath = Union[str, Path]

# Sequence types
DNASequence = str
ProteinSequence = str

# Statistics types
SimulationStats = Dict[str, Any]
```

### Phase 2: Add Type Hints to Functions

**Priority order:**
1. Public API functions (called from CLI)
2. Core simulation logic
3. Utility functions
4. Internal helpers

**Example: simulate.py**

```python
# ❌ Before - No type hints
def simulate_diploid(config, num_haplotypes, fixed_lengths=None, seed=None):
    """Simulate diploid haplotypes."""
    pass

# ✅ After - Full type hints
from typing import Dict, List, Tuple, Optional, Any
from .types import HaplotypeList, ConfigDict

def simulate_diploid(
    config: ConfigDict,
    num_haplotypes: int,
    fixed_lengths: Optional[List[int]] = None,
    seed: Optional[int] = None,
) -> HaplotypeList:
    """Simulate diploid haplotypes.

    Args:
        config: Configuration dictionary with repeats, probabilities, etc.
        num_haplotypes: Number of haplotypes to generate (typically 2 for diploid)
        fixed_lengths: Optional fixed repeat counts per haplotype
        seed: Random seed for reproducibility

    Returns:
        List of (haplotype_name, repeat_chain) tuples

    Raises:
        SimulationError: If haplotype generation fails
    """
    pass
```

**Example: mutate.py**

```python
# ✅ With type hints
from typing import Optional
from .types import HaplotypeList, ConfigDict, MutationName, MutationTargets, MutatedUnits

def apply_mutation(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    mutation_targets: Optional[MutationTargets] = None,
    seed: Optional[int] = None,
) -> Tuple[HaplotypeList, MutatedUnits]:
    """Apply mutation to haplotypes.

    Args:
        config: Configuration with mutation definitions
        results: List of haplotypes to mutate
        mutation_name: Name of mutation from config
        mutation_targets: Optional specific targets, otherwise random
        seed: Random seed for target selection

    Returns:
        Tuple of (mutated_haplotypes, mutated_units_dict)

    Raises:
        MutationError: If mutation cannot be applied
    """
    pass
```

### Phase 3: Configure mypy

Create `mypy.ini`:

```ini
[mypy]
# Python version
python_version = 3.8

# Import discovery
mypy_path = muc_one_up
namespace_packages = True
explicit_package_bases = True

# Strictness (gradual typing)
warn_return_any = True
warn_unused_configs = True
warn_unused_ignores = True
warn_redundant_casts = True
warn_no_return = True
warn_unreachable = True

# Initially lenient, tighten over time
disallow_untyped_defs = False  # Start False, migrate to True
disallow_incomplete_defs = False
check_untyped_defs = True
disallow_untyped_calls = False

# Strict settings for new modules
[mypy-muc_one_up.exceptions]
disallow_untyped_defs = True

[mypy-muc_one_up.config]
disallow_untyped_defs = True

[mypy-muc_one_up.types]
disallow_untyped_defs = True

# Ignore third-party libraries without stubs
[mypy-orfipy.*]
ignore_missing_imports = True

[mypy-jsonschema.*]
ignore_missing_imports = True
```

**Run mypy:**
```bash
# Install
pip install mypy

# Check all files
mypy muc_one_up/

# Check specific file
mypy muc_one_up/config.py

# Generate coverage report
mypy muc_one_up/ --html-report mypy-report/
```

### Phase 4: Runtime Type Checking (Optional)

For critical functions, add runtime validation:

```python
# muc_one_up/validation.py
from typing import Any
from .exceptions import ValidationError


def validate_dna_sequence(sequence: str, allow_n: bool = True) -> None:
    """Validate DNA sequence at runtime.

    Args:
        sequence: DNA sequence to validate
        allow_n: Allow 'N' for ambiguous bases

    Raises:
        ValidationError: If sequence contains invalid characters
    """
    if not sequence:
        raise ValidationError("DNA sequence cannot be empty")

    allowed = {"A", "C", "G", "T"}
    if allow_n:
        allowed.add("N")

    invalid = set(sequence.upper()) - allowed
    if invalid:
        raise ValidationError(
            f"Invalid DNA bases: {sorted(invalid)}. "
            f"Allowed: {sorted(allowed)}"
        )


def validate_haplotype_index(index: int, num_haplotypes: int) -> None:
    """Validate haplotype index is within range.

    Args:
        index: 0-based haplotype index
        num_haplotypes: Total number of haplotypes

    Raises:
        ValidationError: If index out of range
    """
    if not 0 <= index < num_haplotypes:
        raise ValidationError(
            f"Haplotype index {index} out of range. "
            f"Must be 0-{num_haplotypes-1}"
        )


def validate_file_exists(filepath: str, description: str = "File") -> None:
    """Validate file exists.

    Args:
        filepath: Path to file
        description: Description for error message

    Raises:
        FileOperationError: If file doesn't exist
    """
    from pathlib import Path
    from .exceptions import FileOperationError

    if not Path(filepath).exists():
        raise FileOperationError(f"{description} not found: {filepath}")
```

**Use in functions:**

```python
def apply_mutation(
    config: ConfigDict,
    results: HaplotypeList,
    mutation_name: MutationName,
    mutation_targets: Optional[MutationTargets] = None,
    seed: Optional[int] = None,
) -> Tuple[HaplotypeList, MutatedUnits]:
    """Apply mutation with validation."""

    # Runtime validation
    if mutation_name not in config.get("mutations", {}):
        raise MutationError(
            f"Mutation '{mutation_name}' not defined. "
            f"Available: {list(config['mutations'].keys())}"
        )

    if mutation_targets:
        for hap_idx, pos_idx in mutation_targets:
            validate_haplotype_index(hap_idx, len(results))

    # Rest of logic...
```

## Bioinformatics-Specific Validation

### Sequence Validation

```python
# muc_one_up/bioinformatics/validation.py
import re
from typing import Set
from ..exceptions import ValidationError

DNA_BASES: Set[str] = {"A", "C", "G", "T", "N"}
DNA_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)


def validate_dna_sequence(sequence: str, allow_ambiguous: bool = True) -> None:
    """Validate DNA sequence contains only valid bases.

    Args:
        sequence: DNA sequence string
        allow_ambiguous: Allow N for ambiguous bases

    Raises:
        ValidationError: If sequence contains invalid characters
    """
    if not sequence:
        raise ValidationError("Sequence cannot be empty")

    bases = set(sequence.upper())
    allowed = DNA_BASES if allow_ambiguous else DNA_BASES - {"N"}
    invalid = bases - allowed

    if invalid:
        raise ValidationError(
            f"Invalid DNA bases: {sorted(invalid)}. "
            f"Allowed: {sorted(allowed)}"
        )


def validate_fasta_format(fasta_path: str) -> None:
    """Validate FASTA file format.

    Args:
        fasta_path: Path to FASTA file

    Raises:
        ValidationError: If FASTA format is invalid
    """
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValidationError(f"FASTA file is empty: {fasta_path}")

    if not lines[0].startswith(">"):
        raise ValidationError(
            f"FASTA must start with header (>): {fasta_path}"
        )

    in_sequence = False
    for i, line in enumerate(lines, 1):
        if line.startswith(">"):
            in_sequence = True
            continue

        if in_sequence and not DNA_PATTERN.match(line):
            raise ValidationError(
                f"Invalid sequence at line {i}: {line[:50]}..."
            )


def validate_repeat_structure(structure: str, valid_symbols: Set[str]) -> None:
    """Validate repeat structure string.

    Args:
        structure: Hyphen-separated repeat structure (e.g., "1-2-7-8-9")
        valid_symbols: Set of valid repeat symbols

    Raises:
        ValidationError: If structure is invalid
    """
    if not structure:
        raise ValidationError("Repeat structure cannot be empty")

    if "-" not in structure:
        raise ValidationError(
            f"Repeat structure must be hyphen-separated: {structure}"
        )

    repeats = structure.split("-")
    for repeat in repeats:
        # Remove mutation marker if present
        clean_repeat = repeat.rstrip("m")

        if clean_repeat not in valid_symbols:
            raise ValidationError(
                f"Invalid repeat symbol '{repeat}' in structure. "
                f"Valid symbols: {sorted(valid_symbols)}"
            )
```

### Reference Genome Validation

```python
# muc_one_up/bioinformatics/reference_validation.py
from pathlib import Path
from typing import List
from ..exceptions import FileOperationError, ValidationError


def validate_reference_genome(
    reference_path: str,
    aligner: str = "bwa",
) -> List[str]:
    """Validate reference genome and indices exist.

    Args:
        reference_path: Path to reference FASTA
        aligner: Aligner to check indices for ('bwa' or 'minimap2')

    Returns:
        List of warnings (empty if all checks pass)

    Raises:
        FileOperationError: If reference not found
        ValidationError: If required indices missing
    """
    ref_path = Path(reference_path)
    warnings = []

    # Check reference exists
    if not ref_path.exists():
        raise FileOperationError(
            f"Reference genome not found: {reference_path}"
        )

    # Check FASTA index
    fai_path = Path(f"{reference_path}.fai")
    if not fai_path.exists():
        warnings.append(
            f"FASTA index missing: {fai_path}. "
            f"Run: samtools faidx {reference_path}"
        )

    # Check aligner index
    if aligner == "bwa":
        required_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        missing = []
        for ext in required_extensions:
            if not Path(f"{reference_path}{ext}").exists():
                missing.append(ext)

        if missing:
            raise ValidationError(
                f"BWA index files missing: {missing}. "
                f"Run: bwa index {reference_path}"
            )

    elif aligner == "minimap2":
        mmi_path = Path(f"{reference_path}.mmi")
        if not mmi_path.exists():
            warnings.append(
                f"Minimap2 index missing: {mmi_path}. "
                f"Run: minimap2 -d {mmi_path} {reference_path}"
            )

    return warnings
```

## KISS Principle: Add Types Incrementally

**Don't try to type everything at once!**

### Week 1: Core modules
1. `types.py` - Type definitions
2. `exceptions.py` - Already typed
3. `config.py` - Expand existing hints
4. `validation.py` - New validation module

### Week 2: Simulation modules
1. `simulate.py`
2. `mutate.py`
3. `probabilities.py`
4. `distribution.py`

### Week 3: Integration modules
1. `snp_integrator.py`
2. `fasta_writer.py`
3. `io.py`
4. `simulation_statistics.py`

### Week 4: Read simulation
1. `read_simulator/pipeline.py`
2. `read_simulator/ont_pipeline.py`
3. `read_simulator/wrappers/*.py`

## Benefits

| Benefit | Before | After |
|---------|--------|-------|
| **IDE Support** | No autocomplete | Full autocomplete & type checking |
| **Bug Detection** | Runtime only | Caught at development time |
| **Documentation** | Docstrings only | Types + docstrings |
| **Refactoring** | Risky | Safe (mypy catches breaks) |
| **Onboarding** | Must read code | Types document interfaces |

## Checklist

### Type Hints Implementation

- [ ] Create `muc_one_up/types.py` with type aliases
- [ ] Add type hints to `config.py`
- [ ] Add type hints to `simulate.py`
- [ ] Add type hints to `mutate.py`
- [ ] Add type hints to `snp_integrator.py`
- [ ] Add type hints to `probabilities.py`
- [ ] Add type hints to `distribution.py`
- [ ] Add type hints to `fasta_writer.py`
- [ ] Add type hints to `io.py`
- [ ] Add type hints to `simulation_statistics.py`
- [ ] Add type hints to `read_simulator/*.py`
- [ ] Add type hints to CLI functions

### Validation Implementation

- [ ] Create `muc_one_up/validation.py`
- [ ] Create `muc_one_up/bioinformatics/validation.py`
- [ ] Add DNA sequence validation
- [ ] Add FASTA format validation
- [ ] Add repeat structure validation
- [ ] Add reference genome validation
- [ ] Add SNP position validation
- [ ] Add mutation target validation
- [ ] Add file path validation

### Configuration & Testing

- [ ] Create `mypy.ini` configuration
- [ ] Run mypy on all modules
- [ ] Fix type errors
- [ ] Add type checking to CI/CD
- [ ] Add tests for validation functions
- [ ] Update docstrings with type info

### Verification

```bash
# Check type coverage
mypy muc_one_up/ --html-report mypy-report/
# Open mypy-report/index.html to see coverage

# Run type checking in CI
mypy muc_one_up/ --strict

# Count functions with type hints
grep -r "def .*) ->" muc_one_up/ --include="*.py" | wc -l
```

## Success Criteria

- ✅ All public functions have type hints
- ✅ `mypy.ini` configured
- ✅ mypy runs without errors in strict mode
- ✅ Type checking in CI/CD
- ✅ DNA sequence validation implemented
- ✅ Reference genome validation implemented
- ✅ All validation functions tested
- ✅ IDE autocomplete working

## Next Steps

After completing type hints and validation:
1. ➡️ Enable **strict mypy** mode module-by-module
2. ➡️ Add **runtime type checking** with pydantic (optional)
3. ➡️ Generate **API documentation** from type hints (Sphinx)
4. ➡️ Add **property-based testing** for validation (Hypothesis)
