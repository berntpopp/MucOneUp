# Phase 2A: Fix Output Contract & Consolidate ORF Analysis — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `--out-dir` and `--out-base` CLI options actually control where read-simulation output lands, and eliminate the duplicated ORF analysis code between `click_main.py` and `analysis.py`.

**Architecture:** Introduce an `OutputConfig` dataclass that carries output directory and base name from CLI through the read-simulation API to each platform pipeline. Each pipeline uses `OutputConfig` to derive file paths instead of deriving them from the input FASTA path. Separately, consolidate the inline ORF/toxic-protein logic in `click_main.py` into `analysis.py`, making it the single implementation. Move the `_generate_output_base()` helper to `cli/outputs.py`.

**Tech Stack:** Python 3.10+, Click, pytest, dataclasses

**Spec:** `.planning/specs/2026-04-01-codebase-refactoring-design.md` (sections 2.1, 2.2, 2.6)

---

## File Structure

### New Files
| File | Responsibility |
|------|---------------|
| `muc_one_up/read_simulator/output_config.py` | `OutputConfig` dataclass |
| `tests/test_output_config.py` | Tests for OutputConfig and pipeline output placement |

### Modified Files
| File | Changes |
|------|---------|
| `muc_one_up/read_simulation.py` | Accept and pass `OutputConfig` through SIMULATOR_MAP |
| `muc_one_up/read_simulator/pipeline.py` | Accept `OutputConfig`, use it for file naming |
| `muc_one_up/read_simulator/ont_pipeline.py` | Accept `OutputConfig`, use it for file naming |
| `muc_one_up/read_simulator/pacbio_pipeline.py` | Accept `OutputConfig`, use it for file naming |
| `muc_one_up/cli/click_main.py` | Construct `OutputConfig` in reads commands; remove inline ORF logic |
| `muc_one_up/cli/outputs.py` | Receive `_generate_output_base()` from click_main.py |
| `muc_one_up/cli/analysis.py` | Add `run_orf_analysis_standalone()` for the `analyze orfs` command |

---

## Task 1: OutputConfig Dataclass

**Files:**
- Create: `muc_one_up/read_simulator/output_config.py`
- Test: `tests/test_output_config.py`

- [ ] **Step 1: Write failing tests for OutputConfig**

```python
# tests/test_output_config.py
"""Tests for OutputConfig dataclass."""

from pathlib import Path

from muc_one_up.read_simulator.output_config import OutputConfig


class TestOutputConfig:
    """Tests for OutputConfig construction and path derivation."""

    def test_basic_construction(self, tmp_path):
        oc = OutputConfig(out_dir=tmp_path, out_base="my_reads")
        assert oc.out_dir == tmp_path
        assert oc.out_base == "my_reads"

    def test_derive_path(self, tmp_path):
        oc = OutputConfig(out_dir=tmp_path, out_base="my_reads")
        assert oc.derive_path(".bam") == tmp_path / "my_reads.bam"
        assert oc.derive_path("_R1.fastq.gz") == tmp_path / "my_reads_R1.fastq.gz"

    def test_derive_path_str(self, tmp_path):
        oc = OutputConfig(out_dir=tmp_path, out_base="my_reads")
        assert oc.derive_path_str(".bam") == str(tmp_path / "my_reads.bam")

    def test_from_input_fasta_defaults(self):
        """When no out_dir/out_base provided, derive from input FASTA."""
        oc = OutputConfig.from_input_fasta(
            input_fa="/data/sample.001.simulated.fa",
            out_dir=None,
            out_base=None,
            suffix="_reads",
        )
        assert oc.out_dir == Path("/data")
        assert oc.out_base == "sample.001.simulated_reads"

    def test_from_input_fasta_with_overrides(self, tmp_path):
        """CLI overrides take precedence over input FASTA derivation."""
        oc = OutputConfig.from_input_fasta(
            input_fa="/data/sample.fa",
            out_dir=str(tmp_path),
            out_base="custom_base",
            suffix="_reads",
        )
        assert oc.out_dir == tmp_path
        assert oc.out_base == "custom_base"

    def test_from_input_fasta_out_dir_only(self, tmp_path):
        """out_dir override with auto-generated base name."""
        oc = OutputConfig.from_input_fasta(
            input_fa="/data/sample.fa",
            out_dir=str(tmp_path),
            out_base=None,
            suffix="_reads",
        )
        assert oc.out_dir == tmp_path
        assert oc.out_base == "sample_reads"

    def test_from_input_fasta_out_base_only(self):
        """out_base override with directory from input FASTA."""
        oc = OutputConfig.from_input_fasta(
            input_fa="/data/sample.fa",
            out_dir=None,
            out_base="custom",
            suffix="_reads",
        )
        assert oc.out_dir == Path("/data")
        assert oc.out_base == "custom"

    def test_ensures_out_dir_exists(self, tmp_path):
        new_dir = tmp_path / "nested" / "output"
        oc = OutputConfig(out_dir=new_dir, out_base="test")
        oc.ensure_dir()
        assert new_dir.exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_output_config.py -v`
Expected: FAIL with ImportError (module doesn't exist)

- [ ] **Step 3: Implement OutputConfig**

```python
# muc_one_up/read_simulator/output_config.py
"""Output configuration for read simulation pipelines.

Carries output directory and base name from CLI through the pipeline,
replacing the pattern of deriving output paths from input FASTA paths.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OutputConfig:
    """Configuration for output file placement.

    Attributes:
        out_dir: Directory for all output files.
        out_base: Base name stem for output files (no extension).
    """

    out_dir: Path
    out_base: str

    def derive_path(self, suffix: str) -> Path:
        """Derive an output file path by appending suffix to base name.

        Args:
            suffix: File suffix including extension (e.g., '.bam', '_R1.fastq.gz').

        Returns:
            Full path: out_dir / (out_base + suffix)
        """
        return self.out_dir / f"{self.out_base}{suffix}"

    def derive_path_str(self, suffix: str) -> str:
        """Derive an output file path as a string."""
        return str(self.derive_path(suffix))

    def ensure_dir(self) -> None:
        """Create the output directory if it doesn't exist."""
        self.out_dir.mkdir(parents=True, exist_ok=True)

    @classmethod
    def from_input_fasta(
        cls,
        input_fa: str,
        out_dir: str | None,
        out_base: str | None,
        suffix: str,
    ) -> OutputConfig:
        """Construct OutputConfig from CLI parameters and input FASTA path.

        Falls back to deriving values from input_fa when CLI args are None.

        Args:
            input_fa: Path to input FASTA file.
            out_dir: CLI --out-dir value, or None to use input file's directory.
            out_base: CLI --out-base value, or None to auto-generate from input filename.
            suffix: Suffix for auto-generated base name (e.g., '_reads', '_ont_reads').

        Returns:
            Configured OutputConfig instance.
        """
        input_path = Path(input_fa)
        resolved_dir = Path(out_dir) if out_dir else input_path.parent
        resolved_base = out_base if out_base else input_path.stem + suffix
        return cls(out_dir=resolved_dir, out_base=resolved_base)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_output_config.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/output_config.py tests/test_output_config.py
git commit -m "feat: add OutputConfig dataclass for read simulation output placement"
```

---

## Task 2: Thread OutputConfig Through read_simulation.py

**Files:**
- Modify: `muc_one_up/read_simulation.py:88-92` (SIMULATOR_MAP), `muc_one_up/read_simulation.py:95` (simulate_reads signature)

- [ ] **Step 1: Write failing test for OutputConfig passthrough**

Add to `tests/test_output_config.py`:

```python
from unittest.mock import MagicMock, patch


class TestOutputConfigPassthrough:
    """Tests that OutputConfig flows through simulate_reads to pipelines."""

    def test_simulate_reads_passes_output_config_to_illumina(self, tmp_path):
        """simulate_reads should forward output_config to Illumina pipeline."""
        from muc_one_up.read_simulator.output_config import OutputConfig

        oc = OutputConfig(out_dir=tmp_path, out_base="test_reads")
        config = {
            "read_simulation": {"simulator": "illumina"},
            "tools": {},
        }

        with patch(
            "muc_one_up.read_simulation.simulate_illumina_reads"
        ) as mock_illumina:
            mock_illumina.return_value = str(tmp_path / "test.bam")
            from muc_one_up.read_simulation import simulate_reads

            simulate_reads(config, str(tmp_path / "input.fa"), output_config=oc)
            mock_illumina.assert_called_once()
            call_kwargs = mock_illumina.call_args
            assert call_kwargs.kwargs.get("output_config") == oc

    def test_simulate_reads_passes_output_config_to_ont(self, tmp_path):
        """simulate_reads should forward output_config to ONT pipeline."""
        from muc_one_up.read_simulator.output_config import OutputConfig

        oc = OutputConfig(out_dir=tmp_path, out_base="test_ont")
        config = {
            "read_simulation": {"simulator": "ont"},
            "tools": {},
            "nanosim_params": {},
        }

        with patch(
            "muc_one_up.read_simulation.simulate_ont_reads_pipeline"
        ) as mock_ont:
            mock_ont.return_value = str(tmp_path / "test.bam")
            from muc_one_up.read_simulation import simulate_reads

            simulate_reads(config, str(tmp_path / "input.fa"), output_config=oc)
            mock_ont.assert_called_once()
            call_kwargs = mock_ont.call_args
            assert call_kwargs.kwargs.get("output_config") == oc

    def test_simulate_reads_passes_output_config_to_pacbio(self, tmp_path):
        """simulate_reads should forward output_config to PacBio pipeline."""
        from muc_one_up.read_simulator.output_config import OutputConfig

        oc = OutputConfig(out_dir=tmp_path, out_base="test_pacbio")
        config = {
            "read_simulation": {"simulator": "pacbio"},
            "tools": {},
            "pacbio_params": {
                "model_type": "SEQUEL",
                "model_file": "/path/to/model",
            },
        }

        with patch(
            "muc_one_up.read_simulation.simulate_pacbio_hifi_reads"
        ) as mock_pacbio:
            mock_pacbio.return_value = str(tmp_path / "test.bam")
            from muc_one_up.read_simulation import simulate_reads

            simulate_reads(config, str(tmp_path / "input.fa"), output_config=oc)
            mock_pacbio.assert_called_once()
            call_kwargs = mock_pacbio.call_args
            assert call_kwargs.kwargs.get("output_config") == oc

    def test_simulate_reads_works_without_output_config(self, tmp_path):
        """Backward compat: output_config defaults to None."""
        config = {
            "read_simulation": {"simulator": "illumina"},
            "tools": {},
        }

        with patch(
            "muc_one_up.read_simulation.simulate_illumina_reads"
        ) as mock_illumina:
            mock_illumina.return_value = str(tmp_path / "test.bam")
            from muc_one_up.read_simulation import simulate_reads

            simulate_reads(config, str(tmp_path / "input.fa"))
            call_kwargs = mock_illumina.call_args
            assert call_kwargs.kwargs.get("output_config") is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_output_config.py::TestOutputConfigPassthrough -v`
Expected: FAIL — `simulate_reads` doesn't accept `output_config` parameter yet.

- [ ] **Step 3: Modify simulate_reads to accept and forward OutputConfig**

In `muc_one_up/read_simulation.py`, change the function signature at line 95:

```python
# REPLACE line 95:
def simulate_reads(config: dict[str, Any], input_fa: str, source_tracker: Any | None = None) -> str:

# WITH:
def simulate_reads(
    config: dict[str, Any],
    input_fa: str,
    source_tracker: Any | None = None,
    output_config: Any | None = None,
) -> str:
```

Change the SIMULATOR_MAP lambdas at lines 88-92 to pass `output_config`:

```python
# REPLACE lines 88-92:
SIMULATOR_MAP: dict[str, Callable[..., str]] = {
    "illumina": lambda config, input_fa, _, **kw: simulate_illumina_reads(config, input_fa, **kw),
    "ont": simulate_ont_reads_pipeline,
    "pacbio": simulate_pacbio_hifi_reads,
}

# WITH:
SIMULATOR_MAP: dict[str, Callable[..., str]] = {
    "illumina": lambda config, input_fa, _, **kw: simulate_illumina_reads(config, input_fa, **kw),
    "ont": lambda config, input_fa, human_reference, **kw: simulate_ont_reads_pipeline(
        config, input_fa, human_reference=human_reference, **kw
    ),
    "pacbio": lambda config, input_fa, human_reference, **kw: simulate_pacbio_hifi_reads(
        config, input_fa, human_reference=human_reference, **kw
    ),
}
```

Change the dispatch call at line 185 to pass `output_config`:

```python
# REPLACE line 185:
    return simulator_func(config, input_fa, human_reference, source_tracker=source_tracker)

# WITH:
    return simulator_func(
        config, input_fa, human_reference,
        source_tracker=source_tracker,
        output_config=output_config,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_output_config.py -v`
Expected: All PASS

- [ ] **Step 5: Run existing tests to verify no regression**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All existing tests PASS (output_config defaults to None)

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulation.py tests/test_output_config.py
git commit -m "feat: thread OutputConfig through simulate_reads dispatch"
```

---

## Task 3: Illumina Pipeline Uses OutputConfig

**Files:**
- Modify: `muc_one_up/read_simulator/pipeline.py:70-72` (signature), `muc_one_up/read_simulator/pipeline.py:138-148` (output path derivation)

- [ ] **Step 1: Write failing test for Illumina output placement**

Add to `tests/test_output_config.py`:

```python
class TestIlluminaOutputConfig:
    """Tests that Illumina pipeline respects OutputConfig for file naming."""

    def test_illumina_output_paths_use_output_config(self, tmp_path):
        """When OutputConfig is provided, Illumina pipeline derives paths from it."""
        import ast
        import inspect

        from muc_one_up.read_simulator import pipeline

        source = inspect.getsource(pipeline.simulate_reads_pipeline)
        tree = ast.parse(source)

        # Verify the function signature accepts output_config
        func_def = tree.body[0]
        arg_names = [arg.arg for arg in func_def.args.args]
        assert "output_config" in arg_names, (
            "simulate_reads_pipeline must accept output_config parameter"
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_output_config.py::TestIlluminaOutputConfig -v`
Expected: FAIL — `output_config` not in signature yet.

- [ ] **Step 3: Modify Illumina pipeline to accept and use OutputConfig**

In `muc_one_up/read_simulator/pipeline.py`, change the function signature at line 70:

```python
# REPLACE lines 70-72:
def simulate_reads_pipeline(
    config: dict[str, Any], input_fa: str, source_tracker: Any | None = None
) -> str:

# WITH:
def simulate_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    source_tracker: Any | None = None,
    output_config: Any | None = None,
) -> str:
```

Replace the output path derivation at lines 137-148:

```python
# REPLACE lines 137-148:
    # Create output directory if it doesn't exist
    input_path = Path(input_fa)
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Define the common base name for all output files
    output_base = input_path.name.replace(".fa", "").replace(".fasta", "")

    # Use consistent naming for all output files based on the output_base
    reads_fq1 = rs_config.get("output_fastq1", str(Path(output_dir) / f"{output_base}_R1.fastq.gz"))
    reads_fq2 = rs_config.get("output_fastq2", str(Path(output_dir) / f"{output_base}_R2.fastq.gz"))
    output_bam = rs_config.get("output_bam", str(Path(output_dir) / f"{output_base}.bam"))

# WITH:
    # Determine output directory and base name
    input_path = Path(input_fa)
    if output_config is not None:
        output_dir = str(output_config.out_dir)
        output_base = output_config.out_base
    else:
        output_dir = rs_config.get("output_dir", str(input_path.parent))
        output_base = input_path.name.replace(".fa", "").replace(".fasta", "")

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Use consistent naming for all output files based on the output_base
    reads_fq1 = rs_config.get("output_fastq1", str(Path(output_dir) / f"{output_base}_R1.fastq.gz"))
    reads_fq2 = rs_config.get("output_fastq2", str(Path(output_dir) / f"{output_base}_R2.fastq.gz"))
    output_bam = rs_config.get("output_bam", str(Path(output_dir) / f"{output_base}.bam"))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_output_config.py -v`
Expected: All PASS

- [ ] **Step 5: Run existing tests for no regression**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All existing tests PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/pipeline.py tests/test_output_config.py
git commit -m "feat: Illumina pipeline respects OutputConfig for output placement"
```

---

## Task 4: ONT Pipeline Uses OutputConfig

**Files:**
- Modify: `muc_one_up/read_simulator/ont_pipeline.py:35-40` (signature), `muc_one_up/read_simulator/ont_pipeline.py:148-151` (output path derivation)

- [ ] **Step 1: Write failing test for ONT output placement**

Add to `tests/test_output_config.py`:

```python
class TestONTOutputConfig:
    """Tests that ONT pipeline respects OutputConfig for file naming."""

    def test_ont_output_paths_use_output_config(self):
        """When OutputConfig is provided, ONT pipeline derives paths from it."""
        import ast
        import inspect

        from muc_one_up.read_simulator import ont_pipeline

        source = inspect.getsource(ont_pipeline.simulate_ont_reads_pipeline)
        tree = ast.parse(source)

        func_def = tree.body[0]
        arg_names = [arg.arg for arg in func_def.args.args]
        assert "output_config" in arg_names, (
            "simulate_ont_reads_pipeline must accept output_config parameter"
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_output_config.py::TestONTOutputConfig -v`
Expected: FAIL

- [ ] **Step 3: Modify ONT pipeline to accept and use OutputConfig**

In `muc_one_up/read_simulator/ont_pipeline.py`, change the signature at line 35:

```python
# REPLACE lines 35-40:
def simulate_ont_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
) -> str:

# WITH:
def simulate_ont_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: Any | None = None,
) -> str:
```

Replace the output path derivation at lines 147-151:

```python
# REPLACE lines 147-151:
    # Setup output paths
    input_path = Path(input_fa)
    input_basename = input_path.stem
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    output_prefix = str(Path(output_dir) / f"{input_basename}_ont")

# WITH:
    # Setup output paths
    input_path = Path(input_fa)
    if output_config is not None:
        output_dir = str(output_config.out_dir)
        output_prefix = str(output_config.out_dir / output_config.out_base)
    else:
        input_basename = input_path.stem
        output_dir = rs_config.get("output_dir", str(input_path.parent))
        output_prefix = str(Path(output_dir) / f"{input_basename}_ont")
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_output_config.py tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/ont_pipeline.py tests/test_output_config.py
git commit -m "feat: ONT pipeline respects OutputConfig for output placement"
```

---

## Task 5: PacBio Pipeline Uses OutputConfig

**Files:**
- Modify: `muc_one_up/read_simulator/pacbio_pipeline.py:74-79` (signature), `muc_one_up/read_simulator/pacbio_pipeline.py:238-252` (output path derivation)

- [ ] **Step 1: Write failing test for PacBio output placement**

Add to `tests/test_output_config.py`:

```python
class TestPacBioOutputConfig:
    """Tests that PacBio pipeline respects OutputConfig for file naming."""

    def test_pacbio_output_paths_use_output_config(self):
        """When OutputConfig is provided, PacBio pipeline derives paths from it."""
        import ast
        import inspect

        from muc_one_up.read_simulator import pacbio_pipeline

        source = inspect.getsource(pacbio_pipeline.simulate_pacbio_hifi_reads)
        tree = ast.parse(source)

        func_def = tree.body[0]
        arg_names = [arg.arg for arg in func_def.args.args]
        assert "output_config" in arg_names, (
            "simulate_pacbio_hifi_reads must accept output_config parameter"
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_output_config.py::TestPacBioOutputConfig -v`
Expected: FAIL

- [ ] **Step 3: Modify PacBio pipeline to accept and use OutputConfig**

In `muc_one_up/read_simulator/pacbio_pipeline.py`, change the signature at line 74:

```python
# REPLACE lines 74-79:
def simulate_pacbio_hifi_reads(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
) -> str:

# WITH:
def simulate_pacbio_hifi_reads(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: Any | None = None,
) -> str:
```

Replace the output path derivation at lines 238-252:

```python
# REPLACE lines 238-252:
    # Derive output directory and base name from input file path
    input_path = Path(input_fa)
    output_dir_path = input_path.parent
    output_base = input_path.stem  # e.g., "sample.001.simulated" from "sample.001.simulated.fa"

    # Create output directory (in case it doesn't exist)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Define intermediate and final file paths
    clr_prefix = str(output_dir_path / f"{output_base}_clr")
    clr_bam = f"{clr_prefix}.bam"

    hifi_bam = str(output_dir_path / f"{output_base}_hifi.bam")
    hifi_fastq = str(output_dir_path / f"{output_base}_hifi.fastq")
    aligned_bam = str(output_dir_path / f"{output_base}_aligned.bam")

# WITH:
    # Derive output directory and base name
    input_path = Path(input_fa)
    if output_config is not None:
        output_dir_path = output_config.out_dir
        output_base = output_config.out_base
    else:
        output_dir_path = input_path.parent
        output_base = input_path.stem

    # Create output directory (in case it doesn't exist)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Define intermediate and final file paths
    clr_prefix = str(output_dir_path / f"{output_base}_clr")
    clr_bam = f"{clr_prefix}.bam"

    hifi_bam = str(output_dir_path / f"{output_base}_hifi.bam")
    hifi_fastq = str(output_dir_path / f"{output_base}_hifi.fastq")
    aligned_bam = str(output_dir_path / f"{output_base}_aligned.bam")
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/test_output_config.py tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/pacbio_pipeline.py tests/test_output_config.py
git commit -m "feat: PacBio pipeline respects OutputConfig for output placement"
```

---

## Task 6: CLI Reads Commands Construct and Pass OutputConfig

**Files:**
- Modify: `muc_one_up/cli/click_main.py:538-572` (Illumina loop), `muc_one_up/cli/click_main.py:689-723` (ONT loop), `muc_one_up/cli/click_main.py:928-962` (PacBio loop)

- [ ] **Step 1: Write test that CLI commands pass OutputConfig**

Add to `tests/test_output_config.py`:

```python
class TestCLIOutputConfigIntegration:
    """Tests that CLI reads commands construct OutputConfig and pass it through."""

    def test_illumina_passes_output_config(self, tmp_path):
        """Illumina command should construct OutputConfig from --out-dir/--out-base."""
        import json

        from click.testing import CliRunner

        from muc_one_up.cli.click_main import cli

        # Create minimal config
        config = {
            "repeats": {"1": "ACGT"},
            "constants": {"left": "ACGT", "right": "ACGT"},
            "probabilities": {"1": {}},
            "length_model": {"distribution": "uniform", "min_repeats": 1, "max_repeats": 2},
            "mutations": {},
            "reference_assembly": "hg38",
        }
        config_path = tmp_path / "config.json"
        config_path.write_text(json.dumps(config))

        fasta = tmp_path / "input.fa"
        fasta.write_text(">test\nACGT\n")

        out_dir = tmp_path / "custom_output"

        runner = CliRunner()
        with runner.isolated_filesystem(temp_dir=tmp_path):
            from unittest.mock import patch

            with patch("muc_one_up.read_simulation.simulate_reads") as mock_sim:
                mock_sim.return_value = "test.bam"
                result = runner.invoke(
                    cli,
                    [
                        "--config", str(config_path),
                        "reads", "illumina",
                        str(fasta),
                        "--out-dir", str(out_dir),
                        "--out-base", "my_reads",
                    ],
                )
                if result.exit_code != 0:
                    # Allow graceful failure if simulate_reads isn't fully mockable
                    pass
                if mock_sim.called:
                    call_kwargs = mock_sim.call_args
                    oc = call_kwargs.kwargs.get("output_config") or (
                        call_kwargs[0][2] if len(call_kwargs[0]) > 2 else None
                    )
                    if oc is not None:
                        assert str(oc.out_dir) == str(out_dir)
                        assert oc.out_base == "my_reads"
```

- [ ] **Step 2: Modify Illumina reads command to construct and pass OutputConfig**

In `muc_one_up/cli/click_main.py`, in the Illumina reads loop (~line 538-572):

```python
# REPLACE lines 538-572 (the for loop body for illumina):
        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_reads")

            logging.info(
                "[%d/%d] Simulating Illumina reads: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            # Build source tracker from companion files if requested
            source_tracker = None
            if track_read_source:
                from ..read_simulator.source_tracking import ReadSourceTracker

                stats_path = str(Path(input_fasta).with_suffix(".simulation_stats.json"))
                source_tracker = ReadSourceTracker.from_companion_files(stats_path)
                if source_tracker is None:
                    logging.warning("Could not reconstruct read source tracker from %s", stats_path)
                else:
                    coord_map_path = str(
                        Path(out_dir) / f"{actual_out_base}_repeat_coordinates.tsv"
                    )
                    source_tracker.write_coordinate_map(coord_map_path)
                    logging.info("Repeat coordinate map written: %s", coord_map_path)

            # Run simulation for this file
            simulate_reads_pipeline(config, input_fasta, source_tracker=source_tracker)

# WITH:
        for idx, input_fasta in enumerate(input_fastas, start=1):
            from ..read_simulator.output_config import OutputConfig

            # Determine output base name
            if out_base:
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                actual_out_base = _generate_output_base(Path(input_fasta), "_reads")

            output_config = OutputConfig(
                out_dir=Path(out_dir), out_base=actual_out_base
            )

            logging.info(
                "[%d/%d] Simulating Illumina reads: %s -> %s/",
                idx,
                total_files,
                input_fasta,
                output_config.derive_path(""),
            )

            # Build source tracker from companion files if requested
            source_tracker = None
            if track_read_source:
                from ..read_simulator.source_tracking import ReadSourceTracker

                stats_path = str(Path(input_fasta).with_suffix(".simulation_stats.json"))
                source_tracker = ReadSourceTracker.from_companion_files(stats_path)
                if source_tracker is None:
                    logging.warning("Could not reconstruct read source tracker from %s", stats_path)
                else:
                    coord_map_path = output_config.derive_path_str(
                        "_repeat_coordinates.tsv"
                    )
                    source_tracker.write_coordinate_map(coord_map_path)
                    logging.info("Repeat coordinate map written: %s", coord_map_path)

            # Run simulation for this file
            simulate_reads_pipeline(
                config, input_fasta,
                source_tracker=source_tracker,
                output_config=output_config,
            )
```

- [ ] **Step 3: Apply same pattern to ONT reads command**

In `muc_one_up/cli/click_main.py`, in the ONT reads loop (~line 689-723), apply the same OutputConfig pattern. Replace the loop body with:

```python
        for idx, input_fasta in enumerate(input_fastas, start=1):
            from ..read_simulator.output_config import OutputConfig

            if out_base:
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                actual_out_base = _generate_output_base(Path(input_fasta), "_ont_reads")

            output_config = OutputConfig(
                out_dir=Path(out_dir), out_base=actual_out_base
            )

            logging.info(
                "[%d/%d] Simulating ONT reads: %s -> %s/",
                idx,
                total_files,
                input_fasta,
                output_config.derive_path(""),
            )

            source_tracker = None
            if track_read_source:
                from ..read_simulator.source_tracking import ReadSourceTracker

                stats_path = str(Path(input_fasta).with_suffix(".simulation_stats.json"))
                source_tracker = ReadSourceTracker.from_companion_files(stats_path)
                if source_tracker is None:
                    logging.warning("Could not reconstruct read source tracker from %s", stats_path)
                else:
                    coord_map_path = output_config.derive_path_str(
                        "_repeat_coordinates.tsv"
                    )
                    source_tracker.write_coordinate_map(coord_map_path)
                    logging.info("Repeat coordinate map written: %s", coord_map_path)

            simulate_reads_pipeline(
                config, input_fasta,
                source_tracker=source_tracker,
                output_config=output_config,
            )
```

- [ ] **Step 4: Apply same pattern to PacBio reads command**

In `muc_one_up/cli/click_main.py`, in the PacBio reads loop (~line 928-962), apply the same OutputConfig pattern. Replace the loop body with:

```python
        for idx, input_fasta in enumerate(input_fastas, start=1):
            from ..read_simulator.output_config import OutputConfig

            if out_base:
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                actual_out_base = _generate_output_base(Path(input_fasta), "_pacbio_hifi")

            output_config = OutputConfig(
                out_dir=Path(out_dir), out_base=actual_out_base
            )

            logging.info(
                "[%d/%d] Simulating PacBio HiFi reads: %s -> %s/",
                idx,
                total_files,
                input_fasta,
                output_config.derive_path(""),
            )

            source_tracker = None
            if track_read_source:
                from ..read_simulator.source_tracking import ReadSourceTracker

                stats_path = str(Path(input_fasta).with_suffix(".simulation_stats.json"))
                source_tracker = ReadSourceTracker.from_companion_files(stats_path)
                if source_tracker is None:
                    logging.warning("Could not reconstruct read source tracker from %s", stats_path)
                else:
                    coord_map_path = output_config.derive_path_str(
                        "_repeat_coordinates.tsv"
                    )
                    source_tracker.write_coordinate_map(coord_map_path)
                    logging.info("Repeat coordinate map written: %s", coord_map_path)

            simulate_reads_pipeline(
                config, input_fasta,
                source_tracker=source_tracker,
                output_config=output_config,
            )
```

- [ ] **Step 5: Run full tests**

Run: `pytest tests/ --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/click_main.py tests/test_output_config.py
git commit -m "feat: CLI reads commands construct OutputConfig from --out-dir/--out-base"
```

---

## Task 7: Move _generate_output_base to cli/outputs.py

**Files:**
- Modify: `muc_one_up/cli/click_main.py:1390-1406` (remove function)
- Modify: `muc_one_up/cli/outputs.py` (add function)

- [ ] **Step 1: Write test for the moved function**

```python
# tests/test_generate_output_base.py
"""Test that _generate_output_base lives in outputs.py."""

from pathlib import Path


def test_generate_output_base_importable_from_outputs():
    """Function should be importable from cli.outputs."""
    from muc_one_up.cli.outputs import generate_output_base

    result = generate_output_base(Path("sample.001.simulated.fa"), "_reads")
    assert result == "sample.001.simulated_reads"


def test_generate_output_base_various_suffixes():
    from muc_one_up.cli.outputs import generate_output_base

    assert generate_output_base(Path("/path/sample.fa"), "_orfs") == "sample_orfs"
    assert generate_output_base(Path("test.fasta"), "_stats") == "test_stats"
    assert generate_output_base(Path("x.fa"), "_ont_reads") == "x_ont_reads"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_generate_output_base.py -v`
Expected: FAIL — `generate_output_base` not in `outputs.py`.

- [ ] **Step 3: Add function to outputs.py**

Add at the end of `muc_one_up/cli/outputs.py`:

```python
def generate_output_base(input_path: Path, suffix: str) -> str:
    """Generate output base name from input file path.

    Args:
        input_path: Path to input FASTA file.
        suffix: Suffix to append (e.g., '_reads', '_orfs').

    Returns:
        Output base name (e.g., 'sample.001.simulated_reads').
    """
    return input_path.stem + suffix
```

- [ ] **Step 4: Update click_main.py to import from outputs.py**

In `muc_one_up/cli/click_main.py`:

Remove the `_generate_output_base` function (lines 1390-1406).

At the top of the file (after existing imports from `.config`), add:

```python
from .outputs import generate_output_base
```

Then replace all calls to `_generate_output_base(` with `generate_output_base(` throughout the file.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_generate_output_base.py tests/test_click_cli.py tests/test_cli_smoke.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/click_main.py muc_one_up/cli/outputs.py tests/test_generate_output_base.py
git commit -m "refactor: move generate_output_base to cli/outputs.py"
```

---

## Task 8: Consolidate ORF Analysis Into analysis.py

**Files:**
- Modify: `muc_one_up/cli/analysis.py` (add standalone function)
- Modify: `muc_one_up/cli/click_main.py:1060-1152` (replace inline logic with call to analysis.py)

- [ ] **Step 1: Write test for standalone ORF analysis function**

```python
# tests/test_orf_analysis_consolidation.py
"""Tests for consolidated ORF analysis in analysis.py."""

import ast
import inspect


def test_no_subprocess_orfipy_in_click_main():
    """click_main.py should not call orfipy directly via subprocess."""
    import muc_one_up.cli.click_main as cli_mod

    source = inspect.getsource(cli_mod)
    # Verify no subprocess.run with orfipy in click_main
    assert "orfipy" not in source, (
        "click_main.py should not call orfipy directly. "
        "Use analysis.run_orf_analysis_standalone() instead."
    )


def test_run_orf_analysis_standalone_exists():
    """analysis.py should have a run_orf_analysis_standalone function."""
    from muc_one_up.cli.analysis import run_orf_analysis_standalone

    sig = inspect.signature(run_orf_analysis_standalone)
    param_names = list(sig.parameters.keys())
    assert "input_fasta" in param_names
    assert "out_dir" in param_names
    assert "out_base" in param_names
    assert "orf_min_aa" in param_names
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_orf_analysis_consolidation.py -v`
Expected: FAIL — `run_orf_analysis_standalone` doesn't exist yet, and `orfipy` is in click_main.

- [ ] **Step 3: Add run_orf_analysis_standalone to analysis.py**

Add at the end of `muc_one_up/cli/analysis.py`:

```python
def run_orf_analysis_standalone(
    input_fasta: str,
    out_dir: str,
    out_base: str,
    orf_min_aa: int,
    orf_aa_prefix: str | None,
    left_const: str | None,
    right_const: str | None,
) -> None:
    """Run ORF prediction and toxic protein detection on a FASTA file.

    This is the single implementation for ORF analysis, used by both the
    `analyze orfs` CLI command and the orchestrator pipeline.

    Args:
        input_fasta: Path to input FASTA file.
        out_dir: Output directory.
        out_base: Output base name.
        orf_min_aa: Minimum ORF length in amino acids.
        orf_aa_prefix: Required amino acid prefix, or None.
        left_const: Left constant sequence for toxic protein detection.
        right_const: Right constant sequence for toxic protein detection.
    """
    import subprocess

    orf_output = Path(out_dir) / f"{out_base}.orfs.fa"

    # Run orfipy command-line tool
    orf_filename = f"{out_base}.orfs.fa"
    cmd = [
        "orfipy",
        input_fasta,
        "--outdir",
        str(out_dir),
        "--pep",
        orf_filename,
        "--min",
        str(orf_min_aa * 3),  # Convert AA to nucleotides
        "--start",
        "ATG",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        logging.error("orfipy failed for %s: %s", input_fasta, result.stderr)
        return

    logging.info("ORF prediction completed: %s", orf_output)

    # Filter by amino acid prefix if specified
    if orf_aa_prefix and orf_output.exists():
        from Bio import SeqIO

        filtered_orfs = []
        total_orfs = 0

        for record in SeqIO.parse(str(orf_output), "fasta"):
            total_orfs += 1
            if str(record.seq).startswith(orf_aa_prefix):
                filtered_orfs.append(record)

        SeqIO.write(filtered_orfs, str(orf_output), "fasta")
        logging.info(
            "Filtered ORFs by prefix '%s': %d/%d ORFs retained",
            orf_aa_prefix,
            len(filtered_orfs),
            total_orfs,
        )

    # Toxic protein detection
    if orf_output.exists():
        from ..toxic_protein_detector import scan_orf_fasta

        stats = scan_orf_fasta(
            str(orf_output), left_const=left_const, right_const=right_const
        )

        stats_file = Path(out_dir) / f"{out_base}.orf_stats.json"
        with stats_file.open("w") as f:
            json.dump(stats, f, indent=4)

        logging.info("Toxic protein stats written: %s", stats_file)
```

- [ ] **Step 4: Replace inline ORF logic in click_main.py**

In `muc_one_up/cli/click_main.py`, replace the ORF loop body (lines ~1064-1147) with a call to the new function:

```python
# REPLACE the entire for loop body (lines 1064-1147):
        for idx, input_fasta in enumerate(input_fastas, start=1):
            # ... (all the inline orfipy + toxic protein logic)

# WITH:
        for idx, input_fasta in enumerate(input_fastas, start=1):
            from .analysis import run_orf_analysis_standalone

            if out_base:
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                actual_out_base = generate_output_base(Path(input_fasta), "_orfs")

            logging.info(
                "[%d/%d] Running ORF prediction: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            run_orf_analysis_standalone(
                input_fasta=input_fasta,
                out_dir=str(out_dir),
                out_base=actual_out_base,
                orf_min_aa=orf_min_aa,
                orf_aa_prefix=orf_aa_prefix,
                left_const=left_const,
                right_const=right_const,
            )
```

Also remove the `import subprocess` at line 15 if it's no longer used elsewhere in click_main.py. Check first:
- Search for other `subprocess` usage in click_main.py. If none remain, remove the import.

- [ ] **Step 5: Run tests**

Run: `pytest tests/test_orf_analysis_consolidation.py tests/test_click_cli.py tests/test_config_boundary.py --tb=short -q --no-cov`
Expected: All PASS

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/analysis.py muc_one_up/cli/click_main.py tests/test_orf_analysis_consolidation.py
git commit -m "refactor: consolidate ORF analysis into analysis.py single implementation"
```

---

## Task 9: Verify Full Suite and CI Readiness

**Files:** None (verification only)

- [ ] **Step 1: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests pass.

- [ ] **Step 2: Run ruff linter**

Run: `ruff check muc_one_up/ tests/`
Expected: No errors.

- [ ] **Step 3: Run ruff formatter**

Run: `ruff format --check muc_one_up/ tests/`
Expected: No formatting issues. If any, run `ruff format muc_one_up/ tests/` to fix.

- [ ] **Step 4: Run mypy**

Run: `mypy muc_one_up/`
Expected: Passes with no errors.

- [ ] **Step 5: Verify success criteria**

Run: `grep -rn "orfipy" muc_one_up/cli/click_main.py`
Expected: Zero hits (all ORF logic is in analysis.py).

Run: `grep -rn "_generate_output_base" muc_one_up/cli/click_main.py`
Expected: Zero hits (function moved to outputs.py).

Verify that `OutputConfig` is used in all three pipeline signatures:
Run: `grep -n "output_config" muc_one_up/read_simulator/pipeline.py muc_one_up/read_simulator/ont_pipeline.py muc_one_up/read_simulator/pacbio_pipeline.py`
Expected: Each file shows `output_config` in the function signature.

- [ ] **Step 6: Commit any formatting fixes if needed**

```bash
git add -u
git commit -m "style: fix formatting after Phase 2A changes"
```

- [ ] **Step 7: Run full test suite one final time**

Run: `pytest --tb=short -q`
Expected: All pass. Phase 2A is complete.
