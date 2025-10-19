# Testing Strategy

**Priority:** 🟠 HIGH
**Estimated Effort:** 5-7 days
**Impact:** HIGH - Quality assurance, prevents regressions
**Dependencies:** Requires CLI refactoring and exception handling complete

## Current State

### Critical Gaps

| Metric | Current | Target | Gap |
|--------|---------|--------|-----|
| **Test files** | 5 | ~15-20 | -10 to -15 |
| **Test lines** | 470 | ~2,000+ | -1,530 |
| **Test/Code ratio** | 7.7% | 60-80% | -52 to -72% |
| **Coverage measured** | No | Yes | Missing |
| **CLI tests** | 0 | >50 tests | -50 |
| **Integration tests** | 0 | >20 tests | -20 |

### Existing Tests (Good Coverage)

✅ `test_probabilities.py` (~80 lines) - Probability selection logic
✅ `test_mutate.py` (~139 lines) - Mutation application
⚠️ `test_simulate.py` (~120 lines) - Core simulation (partial)
⚠️ `test_config.py` (~47 lines) - Configuration loading (basic)
⚠️ `test_translate.py` (~80 lines) - ORF prediction (basic)

### Missing Tests (Zero Coverage)

❌ CLI interface (`cli.py` - 1,144 lines)
❌ SNP integration (`snp_integrator.py`)
❌ Read simulation pipelines (`read_simulator/*.py`)
❌ External tool wrappers (`read_simulator/wrappers/*.py`)
❌ File I/O (`fasta_writer.py`, `io.py`)
❌ Statistics generation (`simulation_statistics.py`)
❌ Distribution sampling (`distribution.py`)
❌ Toxic protein detection (`toxic_protein_detector.py`)

## KISS Principle: Start Simple

Don't try to achieve 80% coverage in one sprint. Build incrementally:

**Phase 1 (Week 1):** Setup + Critical paths → 30% coverage
**Phase 2 (Week 2):** Core functionality → 50% coverage
**Phase 3 (Week 3):** Edge cases + Integration → 70% coverage
**Phase 4 (Week 4):** Comprehensive coverage + Optimization → 80%+ coverage

## Testing Infrastructure

### Step 1: Create pytest.ini

```ini
# pytest.ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*

# Output settings
addopts =
    --verbose
    --strict-markers
    --tb=short
    --color=yes

    # Coverage settings
    --cov=muc_one_up
    --cov-report=html
    --cov-report=term-missing:skip-covered
    --cov-report=json

    # Warnings
    -W error::UserWarning
    -W ignore::DeprecationWarning

    # Fail if coverage below threshold
    --cov-fail-under=60

# Markers for test organization
markers =
    unit: Unit tests (fast, isolated)
    integration: Integration tests (slower, multiple components)
    slow: Slow tests (> 1 second)
    cli: CLI interface tests
    bioinformatics: Bioinformatics-specific tests (sequence validation, etc.)
    requires_tools: Tests requiring external tools (samtools, bwa, etc.)

# Timeout for slow tests
timeout = 300
```

### Step 2: Create conftest.py with Fixtures

```python
# tests/conftest.py
"""Shared pytest fixtures for MucOneUp tests."""

import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Any

import pytest


# ============================================================================
# Configuration Fixtures
# ============================================================================

@pytest.fixture
def minimal_config() -> Dict[str, Any]:
    """Minimal valid configuration for basic tests."""
    return {
        "repeats": {
            "1": "ATCGATCGATCG",
            "2": "GCTAGCTAGCTA",
            "7": "AAATAAA",
            "8": "TTTCTTT",
            "9": "GGGAGGG",
        },
        "constants": {
            "hg38": {
                "left": "AAAAAAAAAA",
                "right": "TTTTTTTTTT",
                "vntr_start": 10,
                "vntr_end": 50,
            }
        },
        "probabilities": {
            "1": {"2": 1.0},
            "2": {"7": 1.0},
            "7": {"8": 1.0},
            "8": {"9": 1.0},
            "9": {},
        },
        "length_model": {
            "distribution": "normal",
            "min_repeats": 5,
            "max_repeats": 10,
            "mean_repeats": 7,
            "median_repeats": 7,
        },
        "mutations": {
            "dupC": {
                "allowed_repeats": ["1", "2"],
                "strict_mode": False,
                "changes": [
                    {
                        "type": "insert",
                        "position": 1,
                        "sequence": "CCCCCCCC",
                    }
                ],
            }
        },
        "tools": {
            "samtools": "samtools",
            "bwa": "bwa",
            "reseq": "reseq",
        },
        "read_simulation": {
            "simulator": "illumina",
            "human_reference": "/path/to/ref.fa",
            "threads": 1,
            "fragment_size_mean": 300,
            "fragment_size_sd": 50,
            "coverage": 30,
        },
        "reference_assembly": "hg38",
    }


@pytest.fixture
def temp_config_file(tmp_path, minimal_config) -> Path:
    """Create temporary config file."""
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps(minimal_config, indent=2))
    return config_file


# ============================================================================
# File Fixtures
# ============================================================================

@pytest.fixture
def sample_fasta(tmp_path) -> Path:
    """Create sample FASTA file."""
    fasta = tmp_path / "sample.fa"
    fasta.write_text(
        ">haplotype_1\n"
        "ATCGATCGATCGATCGATCGATCG\n"
        ">haplotype_2\n"
        "GCTAGCTAGCTAGCTAGCTAGCTA\n"
    )
    return fasta


@pytest.fixture
def sample_structure_file(tmp_path) -> Path:
    """Create sample structure file."""
    structure = tmp_path / "structure.txt"
    structure.write_text(
        "haplotype_1 1-2-7-8-9\n"
        "haplotype_2 1-2-7-8-9\n"
    )
    return structure


@pytest.fixture
def sample_snp_file(tmp_path) -> Path:
    """Create sample SNP TSV file."""
    snp_file = tmp_path / "snps.tsv"
    snp_file.write_text(
        "1\t10\tA\tG\n"  # Haplotype 1, position 10, A→G
        "2\t20\tC\tT\n"  # Haplotype 2, position 20, C→T
    )
    return snp_file


# ============================================================================
# Output Directory Fixtures
# ============================================================================

@pytest.fixture
def output_dir(tmp_path) -> Path:
    """Create temporary output directory."""
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    return out_dir


# ============================================================================
# Mock Fixtures
# ============================================================================

@pytest.fixture
def mock_external_tools(monkeypatch):
    """Mock external tool execution for faster tests."""

    def mock_subprocess_run(cmd, *args, **kwargs):
        """Mock subprocess.run for external tools."""
        from unittest.mock import Mock

        result = Mock()
        result.returncode = 0
        result.stdout = "Mock output"
        result.stderr = ""
        return result

    import subprocess
    monkeypatch.setattr(subprocess, "run", mock_subprocess_run)


# ============================================================================
# Sequence Fixtures (Bioinformatics)
# ============================================================================

@pytest.fixture
def valid_dna_sequence() -> str:
    """Valid DNA sequence for testing."""
    return "ATCGATCGATCGATCGATCGATCG"


@pytest.fixture
def invalid_dna_sequence() -> str:
    """Invalid DNA sequence (contains non-ACGTN)."""
    return "ATCGXYZATCG"


@pytest.fixture
def sample_haplotypes() -> list:
    """Sample haplotype results for testing."""
    return [
        ("haplotype_1", ["1", "2", "7", "8", "9"]),
        ("haplotype_2", ["1", "2", "7", "8", "9"]),
    ]


# ============================================================================
# Test Data Directory
# ============================================================================

@pytest.fixture
def test_data_dir() -> Path:
    """Path to test data directory."""
    return Path(__file__).parent / "data"


# ============================================================================
# Cleanup
# ============================================================================

@pytest.fixture(autouse=True)
def cleanup_temp_files():
    """Auto-cleanup of temporary files after each test."""
    yield
    # Cleanup happens automatically with tmp_path
```

### Step 3: Create Test Organization

```
tests/
├── conftest.py                      # Shared fixtures
├── data/                            # Test data files
│   ├── sample_config.json
│   ├── invalid_config.json
│   ├── sample.fa
│   └── snps.tsv
├── unit/                            # Unit tests (fast, isolated)
│   ├── test_config.py
│   ├── test_probabilities.py
│   ├── test_distribution.py
│   ├── test_mutate.py
│   ├── test_simulate.py
│   ├── test_snp_integrator.py
│   ├── test_fasta_writer.py
│   ├── test_io.py
│   ├── test_translate.py
│   ├── test_toxic_protein_detector.py
│   └── test_exceptions.py
├── integration/                     # Integration tests
│   ├── test_simulation_workflow.py
│   ├── test_mutation_workflow.py
│   ├── test_read_simulation.py
│   └── test_cli_workflows.py
├── cli/                             # CLI-specific tests
│   ├── test_cli_basic.py
│   ├── test_cli_arguments.py
│   ├── test_cli_simulation.py
│   ├── test_cli_mutation.py
│   └── test_cli_reads.py
└── bioinformatics/                  # Bioinformatics-specific
    ├── test_sequence_validation.py
    ├── test_orf_prediction.py
    └── test_sequence_metrics.py
```

## DRY Principle: Reusable Test Utilities

Create `tests/utils.py`:

```python
"""Test utilities for MucOneUp tests."""

from pathlib import Path
from typing import Dict, Any, List


def assert_valid_fasta(fasta_path: Path) -> None:
    """Assert that FASTA file is valid."""
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]

    assert lines, "FASTA file is empty"
    assert lines[0].startswith(">"), "FASTA must start with header"

    # Check sequence lines contain only valid DNA bases
    for i, line in enumerate(lines):
        if line.startswith(">"):
            continue
        assert set(line.upper()) <= {"A", "C", "G", "T", "N"}, \
            f"Invalid DNA bases at line {i}: {line}"


def assert_structure_file_valid(structure_path: Path, num_haplotypes: int = 2) -> None:
    """Assert that structure file is valid."""
    with open(structure_path) as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    assert len(lines) == num_haplotypes, \
        f"Expected {num_haplotypes} haplotypes, got {len(lines)}"

    for line in lines:
        parts = line.split()
        assert len(parts) == 2, f"Invalid structure line: {line}"
        haplotype_name, structure = parts
        assert haplotype_name.startswith("haplotype"), \
            f"Invalid haplotype name: {haplotype_name}"
        assert "-" in structure, f"Invalid structure format: {structure}"


def count_sequence_bases(sequence: str) -> Dict[str, int]:
    """Count bases in sequence."""
    from collections import Counter
    return dict(Counter(sequence.upper()))


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of sequence."""
    counts = count_sequence_bases(sequence)
    gc = counts.get("G", 0) + counts.get("C", 0)
    total = len(sequence)
    return (gc / total) * 100 if total > 0 else 0.0
```

## Priority Test Coverage Plan

### Phase 1: Critical Paths (Week 1) - Target 30% Coverage

**Files to test:**
1. ✅ `test_exceptions.py` (NEW) - Custom exception hierarchy
2. ✅ `test_config.py` (EXPAND) - Configuration loading
3. ✅ `test_simulate.py` (EXPAND) - Core simulation
4. ✅ `test_mutate.py` (EXISTS) - Keep as-is
5. ✅ `test_probabilities.py` (EXISTS) - Keep as-is

**New tests to write:**

```python
# tests/unit/test_exceptions.py
def test_exception_hierarchy():
    """All custom exceptions inherit from MucOneUpError."""
    pass

def test_external_tool_error_attributes():
    """ExternalToolError preserves tool information."""
    pass


# tests/unit/test_config.py (expand existing)
def test_load_config_success(temp_config_file):
    """Loading valid config succeeds."""
    pass

def test_load_config_missing_file():
    """Loading nonexistent config raises ConfigurationError."""
    pass

def test_load_config_invalid_json(tmp_path):
    """Loading invalid JSON raises ConfigurationError."""
    pass

def test_config_schema_validation(minimal_config):
    """Schema validation catches missing required fields."""
    pass
```

**Success criteria:**
- ✅ pytest runs without errors
- ✅ Coverage report generated
- ✅ Core functionality tested
- ✅ ~30% coverage achieved

### Phase 2: Core Functionality (Week 2) - Target 50% Coverage

**Files to test:**
1. ✅ `test_snp_integrator.py` (NEW)
2. ✅ `test_fasta_writer.py` (NEW)
3. ✅ `test_io.py` (NEW)
4. ✅ `test_distribution.py` (NEW)
5. ✅ `test_simulation_statistics.py` (NEW)

**Example: test_snp_integrator.py**

```python
# tests/unit/test_snp_integrator.py
import pytest
from muc_one_up.snp_integrator import (
    apply_snps_from_file,
    generate_random_snps,
    integrate_snps_to_haplotypes,
)


@pytest.mark.unit
def test_generate_random_snps_density(minimal_config):
    """Random SNPs generated at specified density."""
    snps = generate_random_snps(
        config=minimal_config,
        haplotypes=[("hap1", ["1", "2"])],
        density=1.0,  # 1 SNP per base
        seed=42,
    )
    # Should generate approximately sequence_length SNPs
    assert len(snps) > 0


@pytest.mark.unit
def test_apply_snps_from_file(sample_snp_file, sample_haplotypes, minimal_config):
    """SNPs from file applied correctly."""
    results, snp_info = apply_snps_from_file(
        snp_file=str(sample_snp_file),
        haplotypes=sample_haplotypes,
        config=minimal_config,
    )
    assert len(results) == len(sample_haplotypes)
    assert snp_info["applied"] > 0


@pytest.mark.unit
def test_snp_reference_validation(minimal_config):
    """SNP application validates reference base."""
    # Test that mismatched reference raises error
    pass
```

### Phase 3: Integration & CLI (Week 3) - Target 70% Coverage

**Files to test:**
1. ✅ `test_cli_basic.py` (NEW)
2. ✅ `test_simulation_workflow.py` (NEW)
3. ✅ `test_mutation_workflow.py` (NEW)
4. ✅ `test_read_simulation.py` (NEW)

**Example: CLI testing**

```python
# tests/cli/test_cli_basic.py
from click.testing import CliRunner  # If using Click
from muc_one_up.cli import main


@pytest.mark.cli
def test_cli_help():
    """CLI --help displays usage information."""
    result = main(["--help"])
    assert result == 0
    # Alternatively with Click:
    # runner = CliRunner()
    # result = runner.invoke(cli, ["--help"])
    # assert result.exit_code == 0
    # assert "MucOneUp" in result.output


@pytest.mark.cli
def test_cli_version():
    """CLI --version displays version."""
    pass


@pytest.mark.cli
def test_cli_basic_simulation(temp_config_file, output_dir):
    """Basic simulation completes successfully."""
    pass


@pytest.mark.integration
def test_complete_simulation_workflow(temp_config_file, output_dir):
    """End-to-end simulation workflow."""
    # 1. Run simulation
    # 2. Check FASTA output
    # 3. Check structure file
    # 4. Check stats JSON
    pass
```

### Phase 4: Full Coverage (Week 4) - Target 80%+ Coverage

**Files to test:**
1. ✅ External tool wrappers
2. ✅ Read simulation pipelines
3. ✅ Edge cases
4. ✅ Error scenarios
5. ✅ Bioinformatics validation

## Testing Best Practices (SOLID + KISS)

### Single Responsibility: One Test, One Assertion Concept

```python
# ❌ BAD - Multiple concepts in one test
def test_simulation():
    config = load_config("config.json")
    results = simulate_diploid(config, 2)
    fasta = write_fasta(results)
    stats = generate_stats(results)
    assert len(results) == 2
    assert fasta.exists()
    assert stats["total_length"] > 0
    # What failed if this fails?

# ✅ GOOD - Each test focuses on one concept
def test_simulate_diploid_returns_two_haplotypes(minimal_config):
    results = simulate_diploid(minimal_config, num_haplotypes=2)
    assert len(results) == 2

def test_simulate_diploid_haplotype_structure(minimal_config):
    results = simulate_diploid(minimal_config, num_haplotypes=2)
    name, chain = results[0]
    assert name.startswith("haplotype")
    assert isinstance(chain, list)

def test_write_fasta_creates_file(sample_haplotypes, tmp_path):
    output = tmp_path / "output.fa"
    write_fasta(sample_haplotypes, str(output))
    assert output.exists()
```

### DRY: Use Fixtures, Not Duplication

```python
# ❌ BAD - Setup duplicated in every test
def test_mutation_1():
    config = {"repeats": {...}, "constants": {...}, ...}  # 50 lines
    result = apply_mutation(config, ...)

def test_mutation_2():
    config = {"repeats": {...}, "constants": {...}, ...}  # 50 lines duplicated!
    result = apply_mutation(config, ...)

# ✅ GOOD - Fixture used
@pytest.fixture
def config():
    return {"repeats": {...}, "constants": {...}, ...}

def test_mutation_1(config):
    result = apply_mutation(config, ...)

def test_mutation_2(config):
    result = apply_mutation(config, ...)
```

### KISS: Test Real Behavior, Not Implementation

```python
# ❌ BAD - Testing internal implementation
def test_simulate_uses_random_choice(monkeypatch):
    mock_random = Mock()
    monkeypatch.setattr("random.choice", mock_random)
    simulate_diploid(...)
    assert mock_random.called
    # Fragile - breaks if implementation changes

# ✅ GOOD - Testing observable behavior
def test_simulate_with_seed_is_reproducible(minimal_config):
    results1 = simulate_diploid(minimal_config, seed=42)
    results2 = simulate_diploid(minimal_config, seed=42)
    assert results1 == results2
    # Tests actual requirement: reproducibility
```

## Running Tests

```bash
# All tests
pytest

# Specific category
pytest -m unit
pytest -m integration
pytest -m cli

# Specific file
pytest tests/unit/test_config.py -v

# Specific test
pytest tests/unit/test_config.py::test_load_config_success -v

# With coverage
pytest --cov=muc_one_up --cov-report=html

# Fast tests only (skip slow)
pytest -m "not slow"

# Parallel execution (faster)
pytest -n auto

# Stop on first failure
pytest -x

# Show print statements
pytest -s

# Verbose with tracebacks
pytest -vv --tb=long
```

## Success Criteria

- ✅ `pytest.ini` configured
- ✅ `conftest.py` with reusable fixtures
- ✅ Test coverage >60% (target 80%)
- ✅ All critical paths tested
- ✅ CLI testable (no `sys.exit()`)
- ✅ Fast unit tests (<0.1s each)
- ✅ Integration tests complete (<5s each)
- ✅ Coverage report generated
- ✅ All tests passing in CI/CD

## Next Steps

After achieving 80% test coverage:
1. ➡️ Add **Type Hints** to improve test reliability
2. ➡️ Setup **CI/CD** to run tests automatically
3. ➡️ Add **Mutation Testing** (pytest-mutpy) to verify test quality
4. ➡️ Add **Property-Based Testing** (Hypothesis) for bioinformatics edge cases
