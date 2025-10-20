# Issue #30: Deterministic Simulation - Seed Implementation Plan

**Status**: Ready for Implementation
**Priority**: Phase 1 (Week 1) - Quick Win
**Complexity**: Low-Medium
**Estimated Effort**: 4-6 hours

---

## 📋 Executive Summary

This plan implements deterministic read simulation by adding seed support to all read simulators (Illumina/w-Wessim2, ONT/NanoSim). Following scientific best practices and modern Python patterns, the implementation ensures reproducible results while maintaining backward compatibility.

**Key Principles Applied**:
- ✅ **SOLID**: Optional parameters with sensible defaults (Open/Closed, Liskov Substitution)
- ✅ **DRY**: Centralized seed handling, no duplication
- ✅ **KISS**: Simple `--seed` flag, straightforward propagation
- ✅ **Modular**: Each component handles its own seeding independently

---

## 🎯 Goals

1. **Reproducibility**: Same seed → identical reads across runs
2. **Backward Compatibility**: Default behavior unchanged (no seed = random)
3. **Scientific Rigor**: Proper RNG best practices (avoid global state pollution)
4. **Testing**: Comprehensive validation of deterministic behavior

---

## 📊 Current State Analysis

### ✅ What Already Works
1. **Haplotype simulation** already supports `--seed` (line 142-145 in `cli/click_main.py`)
2. **NanoSim** supports `--seed` parameter natively (confirmed in wrapper)
3. **Python `random` module** used in fragment_simulation.py (line 11)

### ❌ What's Missing
1. **No seed parameter** in `reads illumina` command
2. **No seed parameter** in `reads ont` command
3. **Fragment simulation** doesn't accept seed parameter
4. **NanoSim wrapper** doesn't propagate seed to tool
5. **Config schema** doesn't include seed fields
6. **No tests** for deterministic behavior

---

## 🏗️ Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         CLI Layer                                │
│  muconeup --config X reads illumina file.fa --seed 42          │
│  muconeup --config X reads ont file.fa --seed 42               │
└────────────────────┬────────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Pipeline Layer                                │
│  • simulate_reads_pipeline(config, input_fa)                    │
│  • simulate_ont_reads_pipeline(config, input_fa)                │
│                                                                   │
│  Reads: config["read_simulation"]["seed"]                       │
│         config["nanosim_params"]["seed"]                        │
└────────────────────┬────────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Wrapper Layer                                 │
│  • simulate_fragments(..., seed=None)                           │
│  • run_nanosim_simulation(..., seed=None)                       │
│                                                                   │
│  Implements:                                                     │
│    - random.seed(value) for w-Wessim2                          │
│    - --seed flag for NanoSim                                    │
└─────────────────────────────────────────────────────────────────┘
```

---

## 📝 Detailed Implementation Plan

### Phase 1: Configuration Schema Updates

**File**: `muc_one_up/config.py`
**Lines**: 212-240 (read_simulation section)

#### Changes Required

```python
# Line 212-240: Add seed parameter to read_simulation section
"read_simulation": {
    "type": "object",
    "properties": {
        "simulator": {"type": "string", "enum": ["illumina", "ont"]},
        "reseq_model": {"type": "string"},
        "sample_bam": {"type": "string"},
        # ... existing properties ...
        "seed": {"type": ["number", "null"]},  # ← ADD THIS LINE
    },
    "required": [
        "human_reference",
        "threads",
    ],
    "additionalProperties": False,
},
```

**Line 108-120**: Add seed to nanosim_params section

```python
"nanosim_params": {
    "type": "object",
    "properties": {
        "training_data_path": {"type": "string"},
        "coverage": {"type": "number"},
        "num_threads": {"type": ["number", "null"]},
        "min_read_length": {"type": ["number", "null"]},
        "max_read_length": {"type": ["number", "null"]},
        "other_options": {"type": ["string", "null"]},
        "seed": {"type": ["number", "null"]},  # ← ADD THIS LINE
    },
    "required": ["training_data_path", "coverage"],
    "additionalProperties": False,
},
```

**Rationale**: JSON Schema validation ensures type safety. Using `["number", "null"]` makes seed optional (backward compatibility).

---

### Phase 2: CLI Command Updates

**File**: `muc_one_up/cli/click_main.py`

#### 2.1: Add `--seed` to `reads illumina` command

**Location**: Lines 316-426 (illumina command)

```python
@reads.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
)
@click.option(
    "--coverage",
    type=int,
    default=30,
    show_default=True,
    help="Target sequencing coverage.",
)
@click.option(
    "--threads",
    type=int,
    default=8,
    show_default=True,
    help="Number of threads.",
)
@click.option(  # ← ADD THIS OPTION
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed):  # ← Add seed param
    """Simulate Illumina short reads from one or more FASTA files.

    # ... existing docstring ...

    \b
    Reproducibility:
      Use --seed to generate identical reads across runs:
      muconeup --config X reads illumina file.fa --seed 42
    """
    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure read simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "illumina"
        config["read_simulation"]["coverage"] = coverage
        config["read_simulation"]["threads"] = threads

        # ← ADD SEED HANDLING HERE
        if seed is not None:
            config["read_simulation"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # ... rest of existing code unchanged ...
```

**Lines to modify**:
- Line 347: Add `--seed` option (after `--threads`)
- Line 347: Add `seed` parameter to function signature
- Line 386: Add seed to config (after line 385)

---

#### 2.2: Add `--seed` to `reads ont` command

**Location**: Lines 428-537 (ont command)

```python
@reads.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
# ... existing options ...
@click.option(  # ← ADD THIS OPTION
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed):  # ← Add seed param
    """Simulate Oxford Nanopore long reads from one or more FASTA files.

    # ... existing docstring ...

    \b
    Reproducibility:
      Use --seed to generate identical reads across runs:
      muconeup --config X reads ont file.fa --seed 42
    """
    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure ONT simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "ont"
        config["read_simulation"]["coverage"] = coverage
        if "nanosim_params" not in config:
            config["nanosim_params"] = {}
        config["nanosim_params"]["min_len"] = min_read_length

        # ← ADD SEED HANDLING HERE
        if seed is not None:
            config["nanosim_params"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # ... rest of existing code unchanged ...
```

**Lines to modify**:
- Line 459: Add `--seed` option (after `--min-read-length`)
- Line 459: Add `seed` parameter to function signature
- Line 497: Add seed to config (after line 496)

---

### Phase 3: Fragment Simulation Updates

**File**: `muc_one_up/read_simulator/fragment_simulation.py`

#### 3.1: Update `simulate_fragments` function

**Location**: Lines 297-452

```python
def simulate_fragments(
    ref_fa: str,
    syser_file: str,
    psl_file: str,
    read_number: int,
    fragment_size: int,
    fragment_sd: int,
    min_fragment: int,
    bind: float,
    output_fragments: str,
    seed: int | None = None,  # ← ADD THIS PARAMETER
) -> None:
    """
    Simulate fragments (port of w-Wessim2) and write paired fragment sequences to a FASTA file.

    Args:
        ref_fa: Reference FASTA file.
        syser_file: Systematic errors file.
        psl_file: PSL file.
        read_number: Number of read pairs to generate.
        fragment_size: Mean fragment size.
        fragment_sd: Standard deviation of fragment size.
        min_fragment: Minimum fragment size.
        bind: Minimum fraction (%) for overlap.
        output_fragments: Output FASTA filename for fragments.
        seed: Random seed for reproducibility. If None, uses system randomness.  # ← ADD

    Raises:
        SystemExit: If the simulation fails or produces invalid output.
    """
    # ← ADD SEED INITIALIZATION AT START
    if seed is not None:
        random.seed(seed)
        logging.info(f"Fragment simulation using random seed: {seed}")

    logging.info("Starting fragment simulation (ported w-Wessim2 logic)...")

    # ... rest of function unchanged ...
```

**Lines to modify**:
- Line 306: Add `seed` parameter to function signature
- Line 321: Add `seed` to docstring Args section
- Line 325: Add seed initialization (before logging.info on line 325)

**Rationale**: Using `random.seed()` is acceptable here because:
1. This is a top-level orchestration function (not library code)
2. Setting seed at function start ensures consistent behavior
3. Python's `random` module is used by all helper functions (get_insert_length, pick_on_match, etc.)

---

### Phase 4: Illumina Pipeline Updates

**File**: `muc_one_up/read_simulator/pipeline.py`

#### 4.1: Update `simulate_reads_pipeline` function

**Location**: Lines 203-213 (fragment simulation call)

```python
# Stage 6: Simulate fragments
fragments_fa = str(Path(output_dir) / f"_{output_base}_fragments.fa")
read_number = rs_config.get("read_number", 100000)
fragment_size = rs_config.get("fragment_size", 350)
fragment_sd = rs_config.get("fragment_sd", 50)
min_fragment = rs_config.get("min_fragment", 200)
bind = rs_config.get("binding_min", 0.5)
seed = rs_config.get("seed")  # ← ADD THIS LINE
logging.info("6. Simulating fragments (w-Wessim2)")
simulate_fragments(
    no_ns_fa,
    syser_fq,
    psl_file,
    read_number,
    fragment_size,
    fragment_sd,
    min_fragment,
    bind,
    fragments_fa,
    seed=seed,  # ← ADD THIS PARAMETER
)
```

**Lines to modify**:
- Line 202: Add `seed = rs_config.get("seed")` after `bind` assignment
- Line 213: Add `seed=seed` to `simulate_fragments()` call

---

### Phase 5: NanoSim Wrapper Updates

**File**: `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`

#### 5.1: Update `run_nanosim_simulation` function

**Location**: Lines 20-131

```python
def run_nanosim_simulation(
    nanosim_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    training_model: str,
    coverage: float,
    threads: int = 4,
    min_read_length: int | None = None,
    max_read_length: int | None = None,
    other_options: str = "",
    timeout: int = 3600,
    seed: int | None = None,  # ← ADD THIS PARAMETER
) -> str:
    """
    Run NanoSim simulation to generate Oxford Nanopore reads.

    Args:
        nanosim_cmd: Path to the NanoSim simulator.py command
        reference_fasta: Path to reference FASTA file
        output_prefix: Prefix for output files
        training_model: Path to NanoSim training model
        coverage: Desired coverage
        threads: Number of threads to use (default: 4)
        min_read_length: Minimum read length (optional)
        max_read_length: Maximum read length (optional)
        other_options: Additional NanoSim options (optional)
        timeout: Timeout in seconds (default: 3600)
        seed: Random seed for reproducibility (optional)  # ← ADD

    Returns:
        Path to the generated FASTQ file

    Raises:
        RuntimeError: If NanoSim execution fails
    """
    # ... existing code for output directory setup ...

    # Build the command with required parameters
    # SECURITY: Always use list form, never shell=True
    # Use centralized build_tool_command to safely handle multi-word commands (conda/mamba)
    cmd_list = build_tool_command(
        nanosim_cmd,
        "genome",
        "-rg",
        reference_fasta,
        "-c",
        training_model,
        "-o",
        output_prefix,
        "-t",
        threads,  # build_tool_command handles conversion
        "-x",
        coverage,  # build_tool_command handles conversion
    )

    # ← ADD SEED PARAMETER HANDLING HERE (after line 74)
    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])
        logging.info(f"[NanoSim] Using random seed: {seed}")

    # Add optional parameters
    if min_read_length:
        cmd_list.extend(["--min_len", str(min_read_length)])
    # ... rest of function unchanged ...
```

**Lines to modify**:
- Line 30: Add `seed` parameter to function signature
- Line 45: Add `seed` to docstring Args section
- Line 74: Add seed to command list (after coverage, before min_len check)

---

### Phase 6: ONT Pipeline Updates

**File**: `muc_one_up/read_simulator/ont_pipeline.py`

#### 6.1: Update `simulate_ont_reads_pipeline` function

**Location**: Lines 95-122 (NanoSim call)

```python
# Optional parameters with defaults
threads = ns_params.get("num_threads", rs_config.get("threads", 4))
min_read_length = ns_params.get("min_read_length")
max_read_length = ns_params.get("max_read_length")
other_options = ns_params.get("other_options", "")
seed = ns_params.get("seed")  # ← ADD THIS LINE

# Setup output paths
input_path = Path(input_fa)
# ... existing code ...

# 1. Run NanoSim simulation
logging.info("1. Starting NanoSim simulation")
try:
    fastq_file = run_nanosim_simulation(
        nanosim_cmd=nanosim_cmd,
        reference_fasta=input_fa,
        output_prefix=output_prefix,
        training_model=training_model,
        coverage=float(coverage),
        threads=int(threads),
        min_read_length=min_read_length,
        max_read_length=max_read_length,
        other_options=other_options,
        seed=seed,  # ← ADD THIS PARAMETER
    )
    logging.info("NanoSim simulation completed successfully")
except Exception as e:
    # ... existing error handling ...
```

**Lines to modify**:
- Line 98: Add `seed = ns_params.get("seed")` after `other_options`
- Line 122: Add `seed=seed` to `run_nanosim_simulation()` call

---

## 🧪 Testing Strategy

### Phase 7: Add Comprehensive Tests

#### 7.1: Fragment Simulation Seed Test

**File**: `tests/read_simulator/test_fragment_simulation.py`

**Add new test class** (after line 536):

```python
class TestFragmentSimulationSeeding:
    """Test deterministic behavior with seeds."""

    def test_same_seed_produces_identical_fragments(self, tmp_path):
        """Test that same seed generates identical fragment files."""
        # Arrange: Create input files
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\n" + "ACGT" * 100 + "\n")  # 400 bases

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
            "@chr1 reverse\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
        )

        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t50\t150\t1\t100\t0\t50\n"
        )

        # Act: Run simulation twice with same seed
        output1 = tmp_path / "fragments1.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output1), seed=42
        )

        output2 = tmp_path / "fragments2.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output2), seed=42
        )

        # Assert: Files should be identical
        content1 = output1.read_text()
        content2 = output2.read_text()
        assert content1 == content2, "Same seed should produce identical fragments"

    def test_different_seeds_produce_different_fragments(self, tmp_path):
        """Test that different seeds generate different fragment files."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\n" + "ACGT" * 100 + "\n")

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
            "@chr1 reverse\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
        )

        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t50\t150\t1\t100\t0\t50\n"
        )

        # Act: Run with different seeds
        output1 = tmp_path / "fragments1.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output1), seed=42
        )

        output2 = tmp_path / "fragments2.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output2), seed=123
        )

        # Assert: Files should be different
        content1 = output1.read_text()
        content2 = output2.read_text()
        assert content1 != content2, "Different seeds should produce different fragments"

    def test_no_seed_produces_random_output(self, tmp_path):
        """Test that omitting seed produces non-deterministic output."""
        # Arrange
        ref_fa = tmp_path / "ref.fa"
        ref_fa.write_text(">chr1\n" + "ACGT" * 100 + "\n")

        syser = tmp_path / "test.syser"
        syser.write_text(
            "@chr1 forward\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
            "@chr1 reverse\n" + "N" * 400 + "\n+\n" + "!" * 400 + "\n"
        )

        psl = tmp_path / "test.psl"
        psl.write_text(
            "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery\t150\t0\t100\tchr1\t1000\t50\t150\t1\t100\t0\t50\n"
        )

        # Act: Run twice without seed
        output1 = tmp_path / "fragments1.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output1), seed=None
        )

        output2 = tmp_path / "fragments2.fa"
        simulate_fragments(
            str(ref_fa), str(syser), str(psl), 10, 300, 50, 200, 0.5, str(output2), seed=None
        )

        # Assert: Files are very likely different (probabilistic test)
        content1 = output1.read_text()
        content2 = output2.read_text()
        # Note: There's a tiny chance they could be identical, but astronomically unlikely
        assert content1 != content2, "No seed should produce random (likely different) output"
```

---

#### 7.2: NanoSim Wrapper Seed Test

**File**: `tests/read_simulator/test_nanosim_wrapper.py`

**Add new test** (at end of file):

```python
def test_run_nanosim_with_seed(mocker, tmp_path):
    """Test that seed parameter is correctly passed to NanoSim command."""
    # Arrange
    mock_run_command = mocker.patch("muc_one_up.read_simulator.wrappers.nanosim_wrapper.run_command")
    mock_run_command.return_value = None

    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">chr1\nACGT\n")

    output_prefix = str(tmp_path / "output")
    training_model = "/path/to/model"

    # Create expected output file
    expected_fastq = f"{output_prefix}_aligned_reads.fastq"
    Path(expected_fastq).write_text("@read1\nACGT\n+\nIIII\n")

    # Act
    from muc_one_up.read_simulator.wrappers.nanosim_wrapper import run_nanosim_simulation

    run_nanosim_simulation(
        nanosim_cmd="nanosim",
        reference_fasta=str(ref_fa),
        output_prefix=output_prefix,
        training_model=training_model,
        coverage=30.0,
        threads=4,
        seed=42,
    )

    # Assert
    mock_run_command.assert_called_once()
    called_cmd = mock_run_command.call_args[0][0]

    # Verify --seed is in command
    assert "--seed" in called_cmd
    seed_index = called_cmd.index("--seed")
    assert called_cmd[seed_index + 1] == "42"


def test_run_nanosim_without_seed(mocker, tmp_path):
    """Test that NanoSim runs without seed when not provided."""
    # Arrange
    mock_run_command = mocker.patch("muc_one_up.read_simulator.wrappers.nanosim_wrapper.run_command")
    mock_run_command.return_value = None

    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">chr1\nACGT\n")

    output_prefix = str(tmp_path / "output")
    training_model = "/path/to/model"

    # Create expected output file
    expected_fastq = f"{output_prefix}_aligned_reads.fastq"
    Path(expected_fastq).write_text("@read1\nACGT\n+\nIIII\n")

    # Act
    from muc_one_up.read_simulator.wrappers.nanosim_wrapper import run_nanosim_simulation

    run_nanosim_simulation(
        nanosim_cmd="nanosim",
        reference_fasta=str(ref_fa),
        output_prefix=output_prefix,
        training_model=training_model,
        coverage=30.0,
        threads=4,
        seed=None,  # Explicitly no seed
    )

    # Assert
    mock_run_command.assert_called_once()
    called_cmd = mock_run_command.call_args[0][0]

    # Verify --seed is NOT in command
    assert "--seed" not in called_cmd
```

---

#### 7.3: CLI Integration Test

**File**: `tests/test_click_cli.py` (new file or append to existing)

```python
"""Tests for Click CLI seed parameter handling."""

import json
from unittest.mock import Mock

import pytest
from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


def test_illumina_reads_accepts_seed_parameter(tmp_path, mocker):
    """Test that 'reads illumina' command accepts --seed parameter."""
    # Arrange
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
        },
    }))

    input_fa = tmp_path / "input.fa"
    input_fa.write_text(">seq1\nACGT\n")

    # Mock the simulate_reads_pipeline function
    mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")

    runner = CliRunner()

    # Act
    result = runner.invoke(cli, [
        "--config", str(config_file),
        "reads", "illumina",
        str(input_fa),
        "--seed", "42",
        "--coverage", "10",
    ])

    # Assert
    assert result.exit_code == 0, f"Command failed: {result.output}"

    # Verify seed was passed to config
    mock_simulate.assert_called_once()
    called_config = mock_simulate.call_args[0][0]
    assert called_config["read_simulation"]["seed"] == 42


def test_ont_reads_accepts_seed_parameter(tmp_path, mocker):
    """Test that 'reads ont' command accepts --seed parameter."""
    # Arrange
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
        },
        "nanosim_params": {
            "training_data_path": "/path/to/model",
            "coverage": 30,
        },
    }))

    input_fa = tmp_path / "input.fa"
    input_fa.write_text(">seq1\nACGT\n")

    # Mock the simulate_reads_pipeline function
    mock_simulate = mocker.patch("muc_one_up.read_simulation.simulate_reads")

    runner = CliRunner()

    # Act
    result = runner.invoke(cli, [
        "--config", str(config_file),
        "reads", "ont",
        str(input_fa),
        "--seed", "42",
        "--coverage", "10",
    ])

    # Assert
    assert result.exit_code == 0, f"Command failed: {result.output}"

    # Verify seed was passed to config
    mock_simulate.assert_called_once()
    called_config = mock_simulate.call_args[0][0]
    assert called_config["nanosim_params"]["seed"] == 42
```

---

#### 7.4: Config Schema Validation Test

**File**: `tests/test_config.py`

**Add new tests** (at end of TestLoadConfig class):

```python
def test_read_simulation_accepts_seed(self, tmp_path):
    """Test that read_simulation section accepts seed parameter."""
    config_file = tmp_path / "config.json"
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
            "seed": 42,  # Should be valid
        },
    }
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["read_simulation"]["seed"] == 42


def test_nanosim_params_accepts_seed(self, tmp_path):
    """Test that nanosim_params section accepts seed parameter."""
    config_file = tmp_path / "config.json"
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
        },
        "nanosim_params": {
            "training_data_path": "/path/to/model",
            "coverage": 30,
            "seed": 12345,  # Should be valid
        },
    }
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["nanosim_params"]["seed"] == 12345


def test_seed_can_be_null(self, tmp_path):
    """Test that seed can be explicitly null (for JSON compatibility)."""
    config_file = tmp_path / "config.json"
    config_data = {
        "repeats": {"1": "ACGACT"},
        "constants": {"hg38": {"left": "A", "right": "T", "vntr_region": "chr1:1-100"}},
        "probabilities": {},
        "length_model": {
            "distribution": "normal",
            "min_repeats": 10,
            "max_repeats": 100,
            "mean_repeats": 50,
            "median_repeats": 50,
        },
        "mutations": {},
        "tools": {"samtools": "samtools"},
        "read_simulation": {
            "human_reference": "/ref/hg38.fa",
            "threads": 4,
            "seed": None,  # Explicit null
        },
    }
    config_file.write_text(json.dumps(config_data))

    # Should not raise ValidationError
    config = load_config(str(config_file))
    assert config["read_simulation"]["seed"] is None
```

---

## 📋 Testing Checklist

Before marking this issue as complete, verify:

- [ ] **Unit Tests**
  - [ ] Fragment simulation seed tests pass (3 tests)
  - [ ] NanoSim wrapper seed tests pass (2 tests)
  - [ ] Config schema validation tests pass (3 tests)
  - [ ] CLI integration tests pass (2 tests)

- [ ] **Integration Tests**
  - [ ] Run full Illumina pipeline with seed twice, verify identical output
  - [ ] Run full ONT pipeline with seed twice, verify identical output
  - [ ] Verify different seeds produce different outputs

- [ ] **Regression Tests**
  - [ ] Existing tests still pass (no seed parameter breaks nothing)
  - [ ] Default behavior unchanged (no seed = random output)

- [ ] **Manual Testing**
  ```bash
  # Test 1: Illumina determinism
  muconeup --config config.json reads illumina test.fa --seed 42 --out-base run1
  muconeup --config config.json reads illumina test.fa --seed 42 --out-base run2
  diff run1_R1.fastq.gz run2_R1.fastq.gz  # Should be identical

  # Test 2: ONT determinism
  muconeup --config config.json reads ont test.fa --seed 42 --out-base run1
  muconeup --config config.json reads ont test.fa --seed 42 --out-base run2
  diff run1_ont_aligned_reads.fastq run2_ont_aligned_reads.fastq  # Should be identical

  # Test 3: Different seeds
  muconeup --config config.json reads illumina test.fa --seed 42 --out-base run1
  muconeup --config config.json reads illumina test.fa --seed 123 --out-base run2
  diff run1_R1.fastq.gz run2_R1.fastq.gz  # Should be different
  ```

---

## 🚨 Potential Issues & Mitigation

### Issue 1: Global Random State Pollution

**Problem**: Using `random.seed()` affects global state, which could interfere with other code.

**Mitigation**:
- Fragment simulation is a top-level orchestration function (acceptable use)
- Document seed behavior in function docstrings
- Consider future refactor to use `random.Random()` instance (Phase 4 enhancement)

**Best Practice for Future**:
```python
# Future enhancement (not required for MVP)
if seed is not None:
    rng = random.Random(seed)
    # Pass rng to all helper functions
else:
    rng = random.Random()
```

---

### Issue 2: External Tool Version Differences

**Problem**: NanoSim output may vary across versions even with same seed.

**Mitigation**:
- Document NanoSim version in simulation_stats.json
- Add warning in documentation about version-dependent reproducibility
- Log NanoSim version during execution

**Example Documentation**:
> ⚠️ **Note**: Reproducibility requires identical NanoSim versions. Different versions may produce different results even with the same seed. Use `nanosim --version` to verify version consistency.

---

### Issue 3: Floating Point Precision

**Problem**: Different platforms may have minor floating-point differences.

**Mitigation**:
- Use integer seeds only
- Document platform-specific reproducibility limitations
- This is acceptable for scientific reproducibility (same platform = same results)

---

## 📚 Documentation Updates

### Update 1: CLAUDE.md

**File**: `CLAUDE.md`
**Section**: "Read Simulation" (after line 52)

Add:

```markdown
### Deterministic Read Simulation

Generate reproducible reads using the `--seed` parameter:

```bash
# Illumina reads with seed
muconeup --config config.json reads illumina sample.fa --seed 42

# ONT reads with seed
muconeup --config config.json reads ont sample.fa --seed 42

# Full pipeline with seed (haplotypes + reads)
muconeup --config config.json simulate --seed 42 --out-base sample
muconeup --config config.json reads illumina sample.001.simulated.fa --seed 42
```

**Important**: Same seed guarantees identical output ONLY when:
- Using identical input files
- Running on same platform/architecture
- Using same tool versions (NanoSim, Python)
```

---

### Update 2: config.json Example

**File**: `config.json` (example section)

Add seed examples:

```json
{
  "read_simulation": {
    "simulator": "illumina",
    "human_reference": "/path/to/hg38.fa",
    "threads": 8,
    "coverage": 30,
    "seed": 42
  },
  "nanosim_params": {
    "training_data_path": "/path/to/nanosim_model",
    "coverage": 30,
    "seed": 42
  }
}
```

---

## ⏱️ Implementation Timeline

### Hour 1-2: Configuration & Schema
- [ ] Update `config.py` schema (30 min)
- [ ] Add config validation tests (30 min)
- [ ] Run tests to verify schema changes (30 min)

### Hour 2-3: CLI Updates
- [ ] Add `--seed` to `reads illumina` command (20 min)
- [ ] Add `--seed` to `reads ont` command (20 min)
- [ ] Add CLI integration tests (40 min)

### Hour 3-4: Core Implementation
- [ ] Update `fragment_simulation.py` (30 min)
- [ ] Update `pipeline.py` (15 min)
- [ ] Update `nanosim_wrapper.py` (15 min)
- [ ] Update `ont_pipeline.py` (10 min)

### Hour 4-5: Testing
- [ ] Write fragment simulation seed tests (30 min)
- [ ] Write NanoSim wrapper seed tests (20 min)
- [ ] Run full test suite (20 min)

### Hour 5-6: Documentation & Validation
- [ ] Update CLAUDE.md (20 min)
- [ ] Update config.json examples (10 min)
- [ ] Manual end-to-end testing (40 min)

---

## ✅ Definition of Done

1. **Code Complete**
   - ✅ All 6 files modified as specified
   - ✅ No regressions in existing functionality
   - ✅ Code follows DRY, KISS, SOLID principles

2. **Tests Pass**
   - ✅ All new tests pass (10+ new tests)
   - ✅ All existing tests pass
   - ✅ Manual verification of determinism

3. **Documentation Complete**
   - ✅ CLAUDE.md updated with seed examples
   - ✅ Function docstrings include seed parameter
   - ✅ Example config.json includes seed

4. **Review Complete**
   - ✅ Code reviewed for security issues
   - ✅ No hardcoded seeds in production code
   - ✅ Backward compatibility verified

---

## 🎓 Scientific Best Practices Applied

### ✅ Reproducibility
- Same seed → identical reads (within platform/version constraints)
- Documented version dependencies
- Logged seed values for audit trail

### ✅ Transparency
- Clear documentation of seed behavior
- Explicit logging when seed is used
- Version information in output

### ✅ Flexibility
- Optional parameter (backward compatible)
- Works with existing workflows
- No breaking changes

### ✅ Testability
- Comprehensive unit tests
- Integration tests
- Manual verification procedures

---

## 📊 Success Metrics

After implementation, verify:

1. **Functional Correctness**
   - [ ] `muconeup reads illumina file.fa --seed 42` runs twice produces identical FASTQ
   - [ ] `muconeup reads ont file.fa --seed 42` runs twice produces identical BAM
   - [ ] Different seeds produce different outputs

2. **Code Quality**
   - [ ] Test coverage ≥ 90% for modified functions
   - [ ] No pylint/mypy errors
   - [ ] No security vulnerabilities

3. **Performance**
   - [ ] No measurable performance regression
   - [ ] Seed initialization overhead < 1ms

4. **Usability**
   - [ ] Help text clear and accurate
   - [ ] Error messages informative
   - [ ] Documentation examples work

---

## 🔗 Related Issues/PRs

- Issue #30: Original request for deterministic simulation
- PR #XX: Implementation of seed support (to be created)

---

## 📝 Notes for Implementer

### Priority Order
1. Start with config schema (foundation)
2. CLI updates (user-facing)
3. Core implementation (functionality)
4. Tests (validation)
5. Documentation (completeness)

### Quick Wins
- NanoSim already supports --seed (just need to pass it)
- Python random.seed() is one line
- Most changes are parameter passing (low risk)

### Watch Out For
- Don't forget to propagate seed through all layers
- Test both `seed=42` and `seed=None` cases
- Ensure logging shows seed value when used

---

**Implementation Status**: ⏳ Ready to Start
**Last Updated**: 2025-10-19
**Estimated Completion**: 4-6 hours
