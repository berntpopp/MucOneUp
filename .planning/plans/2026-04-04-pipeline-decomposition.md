# Pipeline Decomposition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decompose 4 read simulation pipeline monoliths into shared utilities + well-defined modules, consolidate CLI dispatch, and add E2E tests.

**Architecture:** Composition with shared utility functions — no base class or protocol. Illumina pipeline splits into `stages/` subpackage; ONT/PacBio/Amplicon get private helper extraction. All pipelines use shared `pipeline_utils.py` for output resolution, metadata, and cleanup. Contract tests validate output naming; E2E tests exercise full workflows with real tools.

**Tech Stack:** Python 3.10+, pytest, pytest-mock (mocker), pysam (E2E), Click (CLI)

**Spec:** `.planning/specs/2026-04-04-pipeline-decomposition-design.md`

---

### Task 1: Add E2E marker to pyproject.toml

**Files:**
- Modify: `pyproject.toml:184-191`

- [ ] **Step 1: Add the e2e marker**

In `pyproject.toml`, add the `e2e` marker to the `markers` list and exclude it from default runs:

```toml
markers = [
    "unit: marks tests as unit tests (fast, isolated)",
    "integration: marks tests as integration tests (slower, multiple components)",
    "slow: marks tests as slow-running (> 1 second)",
    "cli: marks tests as CLI interface tests",
    "bioinformatics: marks tests as bioinformatics-specific (sequence validation, etc.)",
    "requires_tools: marks tests that require specific external tools to be available",
    "e2e: marks end-to-end tests requiring real bioinformatics tools (not run in CI)",
]
```

- [ ] **Step 2: Verify pytest still runs**

Run: `pytest --co -q 2>&1 | tail -5`
Expected: Tests are collected, no errors about unknown markers.

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "ci: add e2e pytest marker for real-tool integration tests"
```

---

### Task 2: Create shared pipeline utilities with tests

**Files:**
- Create: `muc_one_up/read_simulator/pipeline_utils.py`
- Create: `tests/read_simulator/test_pipeline_utils.py`

- [ ] **Step 1: Write tests for `resolve_pipeline_outputs`**

```python
"""Tests for shared pipeline utility functions."""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.output_config import OutputConfig
from muc_one_up.read_simulator.pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_human_reference,
    resolve_pipeline_outputs,
)


class TestResolvePipelineOutputs:
    """Test output path resolution logic."""

    def test_with_output_config(self, tmp_path):
        """Output config takes precedence over rs_config."""
        oc = OutputConfig(out_dir=tmp_path / "out", out_base="mybase")
        output_dir, output_base = resolve_pipeline_outputs(
            input_fa="/some/input.fa",
            rs_config={},
            output_config=oc,
        )
        assert output_dir == tmp_path / "out"
        assert output_base == "mybase"

    def test_without_output_config_uses_rs_config(self, tmp_path):
        """Falls back to rs_config output_dir."""
        output_dir, output_base = resolve_pipeline_outputs(
            input_fa=str(tmp_path / "sample.fa"),
            rs_config={"output_dir": str(tmp_path / "custom")},
            output_config=None,
        )
        assert output_dir == tmp_path / "custom"
        assert output_base == "sample"

    def test_without_output_config_or_rs_config(self, tmp_path):
        """Falls back to input file parent directory."""
        output_dir, output_base = resolve_pipeline_outputs(
            input_fa=str(tmp_path / "data" / "sample.fa"),
            rs_config={},
            output_config=None,
        )
        assert output_dir == tmp_path / "data"
        assert output_base == "sample"

    def test_creates_output_directory(self, tmp_path):
        """Output directory is created if it doesn't exist."""
        oc = OutputConfig(out_dir=tmp_path / "new" / "dir", out_base="base")
        output_dir, _ = resolve_pipeline_outputs(
            input_fa="/input.fa",
            rs_config={},
            output_config=oc,
        )
        assert output_dir.exists()
```

- [ ] **Step 2: Write tests for `resolve_human_reference`**

```python
class TestResolveHumanReference:
    """Test human reference resolution with fallback chain."""

    def test_from_reference_genomes(self, tmp_path, mocker):
        """Uses reference_genomes section when available."""
        ref = tmp_path / "hg38.fa"
        ref.write_text(">chr1\nATCG\n")
        for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"hg38.fa{ext}").touch()

        config = {
            "reference_assembly": "hg38",
            "reference_genomes": {"hg38": {"fasta_path": str(ref)}},
        }
        from muc_one_up.read_simulator.assembly_context import AssemblyContext

        ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )

        result = resolve_human_reference(config, ctx, aligner="bwa")
        assert result == str(ref)

    def test_fallback_to_assembly_context(self, tmp_path, mocker):
        """Falls back to assembly_ctx.human_reference."""
        mocker.patch(
            "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
            side_effect=Exception("not found"),
        )
        from muc_one_up.read_simulator.assembly_context import AssemblyContext

        ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="",
            right_constant="",
            human_reference="/fallback/ref.fa",
        )
        result = resolve_human_reference(config={}, assembly_ctx=ctx, aligner="bwa")
        assert result == "/fallback/ref.fa"

    def test_fallback_to_rs_config(self, mocker):
        """Falls back to rs_config human_reference."""
        mocker.patch(
            "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
            side_effect=Exception("not found"),
        )
        from muc_one_up.read_simulator.assembly_context import AssemblyContext

        ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )
        result = resolve_human_reference(
            config={"read_simulation": {"human_reference": "/legacy/ref.fa"}},
            assembly_ctx=ctx,
            aligner="bwa",
        )
        assert result == "/legacy/ref.fa"

    def test_raises_when_no_reference(self, mocker):
        """Raises ConfigurationError when no reference found."""
        mocker.patch(
            "muc_one_up.read_simulator.pipeline_utils.get_reference_path_for_assembly",
            side_effect=Exception("not found"),
        )
        from muc_one_up.exceptions import ConfigurationError
        from muc_one_up.read_simulator.assembly_context import AssemblyContext

        ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )
        with pytest.raises(ConfigurationError):
            resolve_human_reference(config={}, assembly_ctx=ctx, aligner="bwa")
```

- [ ] **Step 3: Write tests for `create_pipeline_metadata` and `cleanup_intermediates`**

```python
class TestCreatePipelineMetadata:
    """Test metadata file creation."""

    def test_writes_metadata_file(self, tmp_path, mocker):
        """Delegates to write_metadata_file with correct args."""
        from datetime import datetime

        mock_write = mocker.patch(
            "muc_one_up.read_simulator.pipeline_utils.write_metadata_file",
            return_value=str(tmp_path / "metadata.tsv"),
        )
        start = datetime(2026, 1, 1, 10, 0, 0)
        end = datetime(2026, 1, 1, 10, 5, 0)

        result = create_pipeline_metadata(
            output_dir=tmp_path,
            output_base="sample",
            config={"tools": {}},
            start_time=start,
            end_time=end,
            platform="Illumina",
            tools_used=["bwa", "samtools"],
        )

        mock_write.assert_called_once_with(
            output_dir=str(tmp_path),
            output_base="sample",
            config={"tools": {}},
            start_time=start,
            end_time=end,
            platform="Illumina",
            tools_used=["bwa", "samtools"],
        )
        assert result == str(tmp_path / "metadata.tsv")


class TestCleanupIntermediates:
    """Test safe intermediate file cleanup."""

    def test_removes_existing_files(self, tmp_path):
        """Removes files that exist."""
        f1 = tmp_path / "a.tmp"
        f2 = tmp_path / "b.tmp"
        f1.touch()
        f2.touch()
        cleanup_intermediates([str(f1), str(f2)])
        assert not f1.exists()
        assert not f2.exists()

    def test_ignores_missing_files(self, tmp_path):
        """Does not raise for files that don't exist."""
        cleanup_intermediates([str(tmp_path / "nonexistent.tmp")])

    def test_ignores_none_and_empty(self):
        """Handles None and empty strings in list."""
        cleanup_intermediates([None, "", None])

    def test_logs_warning_on_failure(self, tmp_path, mocker):
        """Logs warning when file removal fails."""
        mock_log = mocker.patch("muc_one_up.read_simulator.pipeline_utils.logger")
        mocker.patch("pathlib.Path.unlink", side_effect=PermissionError("nope"))
        f = tmp_path / "locked.tmp"
        f.touch()
        cleanup_intermediates([str(f)])
        mock_log.warning.assert_called()
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_pipeline_utils.py -v 2>&1 | tail -5`
Expected: FAIL — `ModuleNotFoundError: No module named 'muc_one_up.read_simulator.pipeline_utils'`

- [ ] **Step 5: Implement `pipeline_utils.py`**

```python
"""Shared utility functions for read simulation pipelines.

Provides cross-cutting helpers used by all four pipeline variants
(Illumina, ONT, PacBio, Amplicon). Each pipeline retains its own
cleanup orchestration and source tracking; these utilities handle
only the leaf operations that are identical across pipelines.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .assembly_context import AssemblyContext
    from .output_config import OutputConfig

from ..bioinformatics.reference_validation import (
    get_reference_path_for_assembly,
    validate_reference_for_assembly,
)
from ..exceptions import ConfigurationError
from .utils import write_metadata_file

logger = logging.getLogger(__name__)


def resolve_pipeline_outputs(
    input_fa: str,
    rs_config: dict[str, Any],
    output_config: OutputConfig | None,
) -> tuple[Path, str]:
    """Resolve output directory and base name for a pipeline run.

    When output_config is provided it takes full precedence.
    Otherwise falls back to rs_config["output_dir"] then
    input file's parent directory.

    Returns:
        (output_dir, output_base) — directory is created if needed.
    """
    input_path = Path(input_fa)

    if output_config is not None:
        output_dir = output_config.out_dir
        output_base = output_config.out_base
    else:
        output_dir = Path(rs_config.get("output_dir", str(input_path.parent)))
        output_base = input_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir, output_base


def resolve_human_reference(
    config: dict[str, Any],
    assembly_ctx: AssemblyContext,
    aligner: str = "bwa",
) -> str:
    """Resolve human reference with 3-path fallback.

    1. reference_genomes config section (preferred)
    2. assembly_ctx.human_reference
    3. config["read_simulation"]["human_reference"] (legacy)

    Raises:
        ConfigurationError: If no reference can be resolved.
    """
    rs_config = config.get("read_simulation", {})

    # Try reference_genomes section first
    try:
        ref_path = get_reference_path_for_assembly(config, assembly_ctx.assembly_name)
        human_reference = str(ref_path)

        warnings = validate_reference_for_assembly(
            config, assembly_ctx.assembly_name, aligner=aligner
        )
        if warnings:
            for warning in warnings:
                logger.warning("%s", warning)

        logger.info(
            "Using reference from config: %s (%s)",
            human_reference,
            assembly_ctx.assembly_name,
        )
        return human_reference
    except Exception as e:
        logger.debug("Could not load reference from reference_genomes: %s", e)

    # Fallback to assembly context or legacy config
    human_reference = assembly_ctx.human_reference or rs_config.get("human_reference")
    if human_reference:
        return human_reference

    raise ConfigurationError(
        "human_reference not specified. Add 'reference_genomes' section or "
        "'human_reference' to config 'read_simulation' section"
    )


def create_pipeline_metadata(
    output_dir: Path,
    output_base: str,
    config: dict[str, Any],
    start_time: datetime,
    end_time: datetime,
    platform: str,
    tools_used: list[str],
) -> str:
    """Write pipeline metadata file.

    Thin wrapper around write_metadata_file that accepts Path for
    output_dir (the pipelines already have it as Path).
    """
    return write_metadata_file(
        output_dir=str(output_dir),
        output_base=output_base,
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform=platform,
        tools_used=tools_used,
    )


def cleanup_intermediates(file_list: list[str | None]) -> None:
    """Remove intermediate files, logging warnings on failure.

    Safe leaf operation: silently skips None, empty strings,
    and non-existent paths. Never raises.
    """
    for file in file_list:
        if not file:
            continue
        file_path = Path(file)
        if file_path.exists():
            try:
                file_path.unlink()
                logger.info("Removed intermediate file: %s", file)
            except Exception as e:
                logger.warning("Failed to remove file %s: %s", file, e)
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_pipeline_utils.py -v`
Expected: All tests PASS.

- [ ] **Step 7: Run linter and type checker**

Run: `ruff check muc_one_up/read_simulator/pipeline_utils.py tests/read_simulator/test_pipeline_utils.py && mypy muc_one_up/read_simulator/pipeline_utils.py`

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/read_simulator/pipeline_utils.py tests/read_simulator/test_pipeline_utils.py
git commit -m "feat: add shared pipeline_utils module with output resolution, reference fallback, metadata, and cleanup helpers"
```

---

### Task 3: Create Illumina stages dataclasses

**Files:**
- Create: `muc_one_up/read_simulator/stages/__init__.py`

- [ ] **Step 1: Write the dataclasses module**

```python
"""Result dataclasses for Illumina pipeline stages.

These frozen dataclasses define the contracts between the orchestrator
(pipeline.py) and the extracted stage modules.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FragmentResult:
    """Result of the fragment preparation stage (Illumina stages 1-8).

    Attributes:
        r1_fastq: Path to the R1 paired-end FASTQ (gzipped).
        r2_fastq: Path to the R2 paired-end FASTQ (gzipped).
        intermediate_files: Paths to intermediate files created during
            fragment preparation, suitable for cleanup.
    """

    r1_fastq: str
    r2_fastq: str
    intermediate_files: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class AlignmentResult:
    """Result of the alignment and refinement stage (Illumina stages 9-11).

    Attributes:
        final_bam: Path to the final BAM file (may be aligned, VNTR-biased,
            or downsampled depending on config).
        intermediate_bams: Paths to BAM files superseded during processing
            (e.g., pre-VNTR-bias BAM, pre-downsampled BAM).
        intermediate_files: Paths to non-BAM intermediate files (depth files).
    """

    final_bam: str
    intermediate_bams: list[str] = field(default_factory=list)
    intermediate_files: list[str] = field(default_factory=list)
```

- [ ] **Step 2: Verify import works**

Run: `python -c "from muc_one_up.read_simulator.stages import FragmentResult, AlignmentResult; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add muc_one_up/read_simulator/stages/__init__.py
git commit -m "feat: add FragmentResult and AlignmentResult dataclasses for Illumina stage contracts"
```

---

### Task 4: Extract fragment preparation (Illumina stages 1-8)

**Files:**
- Create: `muc_one_up/read_simulator/stages/fragment_preparation.py`
- Create: `tests/read_simulator/test_fragment_preparation.py`

- [ ] **Step 1: Write test for prepare_fragments**

```python
"""Tests for Illumina fragment preparation stage."""

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from muc_one_up.read_simulator.assembly_context import AssemblyContext
from muc_one_up.read_simulator.stages import FragmentResult


class TestPrepareFragments:
    """Test the fragment preparation pipeline (stages 1-8)."""

    def test_returns_fragment_result(self, tmp_path, mocker):
        """prepare_fragments returns a FragmentResult with R1, R2, and intermediates."""
        # Mock all wrapper calls
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.replace_Ns")
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.generate_systematic_errors"
        )
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.fa_to_twobit")
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.extract_subset_reference",
            return_value=str(tmp_path / "_collated.bam"),
        )
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.run_pblat")
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.simulate_fragments")
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.create_reads")
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.split_reads")

        from muc_one_up.read_simulator.stages.fragment_preparation import prepare_fragments

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq\nATCG\n")

        tools = {
            "reseq": "reseq",
            "faToTwoBit": "faToTwoBit",
            "samtools": "samtools",
            "pblat": "pblat",
        }
        rs_config = {
            "reseq_model": "/model",
            "sample_bam": "/sample.bam",
            "read_number": 100,
            "fragment_size": 350,
            "fragment_sd": 50,
            "min_fragment": 200,
            "threads": 2,
        }
        assembly_ctx = AssemblyContext(
            assembly_name="hg38",
            left_constant="AAA",
            right_constant="TTT",
            sample_bam="/sample.bam",
        )

        result = prepare_fragments(
            tools=tools,
            rs_config=rs_config,
            input_fa=str(input_fa),
            output_dir=tmp_path,
            output_base="test",
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, FragmentResult)
        assert result.r1_fastq.endswith("_R1.fastq.gz")
        assert result.r2_fastq.endswith("_R2.fastq.gz")
        assert len(result.intermediate_files) > 0

    def test_raises_on_missing_reseq_model(self, tmp_path):
        """Raises ConfigurationError when reseq_model is missing."""
        from muc_one_up.exceptions import ConfigurationError
        from muc_one_up.read_simulator.stages.fragment_preparation import prepare_fragments

        assembly_ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )

        with pytest.raises(ConfigurationError, match="reseq_model"):
            prepare_fragments(
                tools={},
                rs_config={},
                input_fa=str(tmp_path / "input.fa"),
                output_dir=tmp_path,
                output_base="test",
                assembly_ctx=assembly_ctx,
            )

    def test_raises_on_missing_sample_bam(self, tmp_path, mocker):
        """Raises ConfigurationError when no sample BAM is available."""
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.replace_Ns")
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.generate_systematic_errors"
        )
        mocker.patch("muc_one_up.read_simulator.stages.fragment_preparation.fa_to_twobit")

        from muc_one_up.exceptions import ConfigurationError
        from muc_one_up.read_simulator.stages.fragment_preparation import prepare_fragments

        assembly_ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )

        with pytest.raises(ConfigurationError, match="sample BAM"):
            prepare_fragments(
                tools={"reseq": "reseq", "faToTwoBit": "faToTwoBit", "samtools": "samtools", "pblat": "pblat"},
                rs_config={"reseq_model": "/model"},
                input_fa=str(tmp_path / "input.fa"),
                output_dir=tmp_path,
                output_base="test",
                assembly_ctx=assembly_ctx,
            )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_fragment_preparation.py -v 2>&1 | tail -5`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `fragment_preparation.py`**

Extract lines 193-313 from `pipeline.py` into a function. The function takes explicit parameters (not the full config dict) and returns `FragmentResult`:

```python
"""Illumina fragment preparation stage (pipeline stages 1-8).

Transforms an input FASTA into paired-end FASTQ reads through:
Replace Ns -> systematic errors -> 2bit -> subset reference ->
pblat alignment -> fragment simulation -> read creation -> read splitting.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..assembly_context import AssemblyContext

from ...exceptions import ConfigurationError
from ..fragment_simulation import simulate_fragments
from ..wrappers.reseq_wrapper import (
    create_reads,
    generate_systematic_errors,
    replace_Ns,
    split_reads,
)
from ..wrappers.samtools_wrapper import extract_subset_reference
from ..wrappers.ucsc_tools_wrapper import fa_to_twobit, run_pblat
from . import FragmentResult

logger = logging.getLogger(__name__)


def prepare_fragments(
    tools: dict[str, str],
    rs_config: dict[str, Any],
    input_fa: str,
    output_dir: Path,
    output_base: str,
    assembly_ctx: AssemblyContext,
    source_tracker: Any | None = None,
) -> FragmentResult:
    """Run Illumina fragment preparation (stages 1-8).

    Args:
        tools: Dict of tool command paths (reseq, faToTwoBit, samtools, pblat).
        rs_config: The read_simulation config subsection.
        input_fa: Path to the input simulated FASTA.
        output_dir: Output directory (already exists).
        output_base: Base name for output files.
        assembly_ctx: Resolved assembly context.
        source_tracker: Optional source tracker (triggers fragment origins sidecar).

    Returns:
        FragmentResult with R1/R2 FASTQ paths and intermediate file list.
    """
    intermediate_files: list[str] = []

    # Stage 1: Replace Ns
    no_ns_fa = str(output_dir / f"_{output_base}_noNs.fa")
    logger.info("1. Replacing Ns in FASTA")
    replace_Ns(input_fa, no_ns_fa, tools)
    intermediate_files.append(no_ns_fa)

    # Stage 2: Generate systematic errors
    reseq_model = rs_config.get("reseq_model")
    if not reseq_model:
        raise ConfigurationError("reseq_model not specified in config 'read_simulation' section")
    syser_fq = str(output_dir / f"_{output_base}_syser.fq")
    logger.info("2. Generating systematic errors")
    generate_systematic_errors(no_ns_fa, reseq_model, syser_fq, tools)
    intermediate_files.append(syser_fq)

    # Stage 3: Convert FASTA to 2bit
    twobit_file = str(output_dir / f"_{output_base}_noNs.2bit")
    logger.info("3. Converting FASTA to 2bit format")
    fa_to_twobit(no_ns_fa, twobit_file, tools)
    intermediate_files.append(twobit_file)

    # Stage 4: Extract subset reference from BAM
    sample_bam = assembly_ctx.sample_bam
    if not sample_bam:
        sample_bam = rs_config.get("sample_bam")
        if sample_bam:
            logger.warning(
                "Assembly-specific sample BAM not found for %s. Using generic sample_bam: %s",
                assembly_ctx.assembly_name,
                sample_bam,
            )
    if not sample_bam:
        raise ConfigurationError(
            f"No sample BAM specified for {assembly_ctx.assembly_name}. "
            f"Add 'sample_bam_{assembly_ctx.assembly_name}' or 'sample_bam' to config"
        )
    logger.info("Using sample BAM for %s: %s", assembly_ctx.assembly_name, sample_bam)
    subset_ref = str(output_dir / f"_{output_base}_subset_ref.fa")
    logger.info("4. Extracting subset reference from BAM")
    collated_bam = extract_subset_reference(sample_bam, subset_ref, tools)
    intermediate_files.append(subset_ref)
    intermediate_files.append(collated_bam)

    # Stage 5: Run pblat alignment
    psl_file = str(output_dir / f"_{output_base}_alignment.psl")
    threads = rs_config.get("threads", 4)
    pblat_threads = min(rs_config.get("pblat_threads", 24), threads)
    logger.info("5. Running pblat alignment")
    run_pblat(
        twobit_file,
        subset_ref,
        psl_file,
        tools,
        threads=pblat_threads,
        min_score=rs_config.get("pblat_min_score", 95),
        min_identity=rs_config.get("pblat_min_identity", 95),
    )

    # Stage 6: Simulate fragments
    fragments_fa = str(output_dir / f"_{output_base}_fragments.fa")
    read_number = rs_config.get("read_number", 100000)
    fragment_size = rs_config.get("fragment_size", 350)
    fragment_sd = rs_config.get("fragment_sd", 50)
    min_fragment = rs_config.get("min_fragment", 200)
    bind = rs_config.get("binding_min", 0.5)
    seed = rs_config.get("seed")
    logger.info("6. Simulating fragments (Wessim2-style, ported logic)")

    fragment_origins_path = (
        str(output_dir / f"{output_base}_fragment_origins.tsv")
        if source_tracker is not None
        else None
    )
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
        seed=seed,
        fragment_origins_path=fragment_origins_path,
    )
    intermediate_files.append(psl_file)
    intermediate_files.append(fragments_fa)

    # Stage 7: Create reads from fragments
    reads_fq = str(output_dir / f"_{output_base}_reads.fq")
    logger.info("7. Creating reads from fragments")
    seqtoillumina_timeout = rs_config.get("seqtoillumina_timeout", 120)
    create_reads(fragments_fa, reseq_model, reads_fq, threads, tools, timeout=seqtoillumina_timeout)

    # Stage 8: Split reads into paired FASTQ files
    reads_fq1 = str(output_dir / f"{output_base}_R1.fastq.gz")
    reads_fq2 = str(output_dir / f"{output_base}_R2.fastq.gz")
    logger.info("8. Splitting reads into paired FASTQ files")
    split_reads(reads_fq, reads_fq1, reads_fq2)
    intermediate_files.append(reads_fq)

    return FragmentResult(
        r1_fastq=reads_fq1,
        r2_fastq=reads_fq2,
        intermediate_files=intermediate_files,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_fragment_preparation.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/stages/ && mypy muc_one_up/read_simulator/stages/`

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/stages/fragment_preparation.py tests/read_simulator/test_fragment_preparation.py
git commit -m "feat: extract Illumina fragment preparation into stages/fragment_preparation.py (stages 1-8)"
```

---

### Task 5: Extract alignment and refinement (Illumina stages 9-11)

**Files:**
- Create: `muc_one_up/read_simulator/stages/alignment.py`
- Create: `tests/read_simulator/test_alignment_stage.py`

- [ ] **Step 1: Write test for align_and_refine**

```python
"""Tests for Illumina alignment and refinement stage."""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.assembly_context import AssemblyContext
from muc_one_up.read_simulator.stages import AlignmentResult


class TestAlignAndRefine:
    """Test the alignment and refinement pipeline (stages 9-11)."""

    def test_returns_alignment_result_basic(self, tmp_path, mocker):
        """Basic alignment without VNTR bias or downsampling."""
        mocker.patch("muc_one_up.read_simulator.stages.alignment.align_reads")

        from muc_one_up.read_simulator.stages.alignment import align_and_refine

        assembly_ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )

        result = align_and_refine(
            tools={"bwa": "bwa", "samtools": "samtools"},
            rs_config={"vntr_capture_efficiency": {"enabled": False}},
            r1=str(tmp_path / "R1.fastq.gz"),
            r2=str(tmp_path / "R2.fastq.gz"),
            human_ref="/ref.fa",
            output_dir=tmp_path,
            output_base="test",
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, AlignmentResult)
        assert result.final_bam.endswith("test.bam")
        assert result.intermediate_bams == []

    def test_vntr_bias_produces_intermediate_bam(self, tmp_path, mocker):
        """VNTR efficiency bias creates an intermediate BAM."""
        mocker.patch("muc_one_up.read_simulator.stages.alignment.align_reads")

        mock_vntr = mocker.MagicMock()
        mock_vntr.apply_efficiency_bias.return_value = {"reads_kept": 100}
        mocker.patch(
            "muc_one_up.read_simulator.stages.alignment.VNTREfficiencyModel",
            return_value=mock_vntr,
        )
        mocker.patch("shutil.rmtree")

        from muc_one_up.read_simulator.stages.alignment import align_and_refine

        assembly_ctx = AssemblyContext(
            assembly_name="hg38", left_constant="", right_constant=""
        )

        result = align_and_refine(
            tools={"bwa": "bwa", "samtools": "samtools"},
            rs_config={"vntr_capture_efficiency": {"enabled": True}},
            r1=str(tmp_path / "R1.fastq.gz"),
            r2=str(tmp_path / "R2.fastq.gz"),
            human_ref="/ref.fa",
            output_dir=tmp_path,
            output_base="test",
            assembly_ctx=assembly_ctx,
        )

        assert isinstance(result, AlignmentResult)
        assert "vntr_biased" in result.final_bam
        assert len(result.intermediate_bams) == 1  # original aligned BAM
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_alignment_stage.py -v 2>&1 | tail -5`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `alignment.py`**

Extract lines 315-561 from `pipeline.py`. Key change: use explicit variable names instead of mutating `output_bam`:

```python
"""Illumina alignment and refinement stage (pipeline stages 9-11).

Aligns paired-end reads to human reference, optionally applies VNTR
capture efficiency bias, and optionally downsamples to target coverage.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..assembly_context import AssemblyContext

from ...exceptions import ConfigurationError
from ..wrappers.bwa_wrapper import align_reads
from ..wrappers.samtools_wrapper import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    downsample_bam,
    downsample_entire_bam,
)
from . import AlignmentResult

logger = logging.getLogger(__name__)


def align_and_refine(
    tools: dict[str, str],
    rs_config: dict[str, Any],
    r1: str,
    r2: str,
    human_ref: str,
    output_dir: Path,
    output_base: str,
    assembly_ctx: AssemblyContext,
) -> AlignmentResult:
    """Run alignment, optional VNTR bias, and optional downsampling.

    Args:
        tools: Dict of tool command paths (bwa, samtools).
        rs_config: The read_simulation config subsection.
        r1: Path to R1 FASTQ.
        r2: Path to R2 FASTQ.
        human_ref: Path to human reference FASTA.
        output_dir: Output directory.
        output_base: Base name for output files.
        assembly_ctx: Resolved assembly context.

    Returns:
        AlignmentResult with final BAM and intermediate file lists.
    """
    intermediate_bams: list[str] = []
    intermediate_files: list[str] = []
    threads = rs_config.get("threads", 4)

    # Stage 9: Align reads
    aligned_bam = str(output_dir / f"{output_base}.bam")
    logger.info("9. Aligning reads to reference: %s", human_ref)
    align_reads(r1, r2, human_ref, aligned_bam, tools, threads)

    # Stage 10: Apply VNTR capture efficiency bias (optional)
    vntr_config = rs_config.get("vntr_capture_efficiency", {})
    vntr_enabled = vntr_config.get("enabled", True)
    current_bam = aligned_bam

    if vntr_enabled:
        logger.info("10. Applying VNTR capture efficiency bias")
        try:
            from ..vntr_efficiency import VNTREfficiencyModel

            penalty_factor = vntr_config.get("penalty_factor", 0.375)
            vntr_seed = vntr_config.get("seed", rs_config.get("seed", 42))
            vntr_region = vntr_config.get("vntr_region")
            capture_bed = vntr_config.get("capture_bed")
            flanking_size = vntr_config.get("flanking_size", 10000)

            vntr_model = VNTREfficiencyModel(
                penalty_factor=penalty_factor,
                seed=vntr_seed,
                threads=threads,
                vntr_region=vntr_region,
                capture_bed=Path(capture_bed) if capture_bed else None,
                flanking_size=flanking_size,
            )

            vntr_biased_bam = str(
                output_dir / f"{output_base}_vntr_biased.bam"
            )
            temp_dir = output_dir / "_vntr_efficiency_temp"

            vntr_stats = vntr_model.apply_efficiency_bias(
                input_bam=Path(current_bam),
                output_bam=Path(vntr_biased_bam),
                temp_dir=temp_dir,
            )

            if vntr_config.get("validation", {}).get("report_statistics", True):
                stats_file = str(output_dir / f"{output_base}_vntr_efficiency_stats.json")
                with open(stats_file, "w") as f:
                    json.dump(vntr_stats, f, indent=2)
                logger.info("  VNTR efficiency statistics saved to %s", stats_file)

            # Generate FASTQ from VNTR-biased BAM (optional)
            vntr_fastq_config = vntr_config.get("output_fastq", {})
            if vntr_fastq_config.get("enabled", True):
                logger.info("  Generating FASTQ files from VNTR-biased BAM")
                try:
                    from ..wrappers.samtools_wrapper import (
                        FastqConversionOptions,
                        convert_bam_to_paired_fastq,
                    )

                    vntr_base = Path(vntr_biased_bam).name.replace(".bam", "")
                    vntr_fq1 = str(output_dir / f"{vntr_base}_R1.fastq.gz")
                    vntr_fq2 = str(output_dir / f"{vntr_base}_R2.fastq.gz")

                    opts = FastqConversionOptions(
                        output_singleton=vntr_fastq_config.get("singleton_file"),
                        preserve_read_names=vntr_fastq_config.get("preserve_read_names", True),
                        validate_pairs=True,
                        threads=threads,
                        timeout=1800,
                    )
                    fq1_out, fq2_out = convert_bam_to_paired_fastq(
                        samtools_cmd=tools["samtools"],
                        input_bam=vntr_biased_bam,
                        output_fq1=vntr_fq1,
                        output_fq2=vntr_fq2,
                        options=opts,
                    )
                    logger.info("  VNTR-biased FASTQ R1: %s", Path(fq1_out).name)
                    logger.info("  VNTR-biased FASTQ R2: %s", Path(fq2_out).name)
                except Exception as e:
                    logger.warning("  Failed to generate VNTR-biased FASTQ files: %s", e)
                    logger.debug("  Exception details:", exc_info=True)
            else:
                logger.info("  Skipping VNTR-biased FASTQ generation (disabled in config)")

            if temp_dir.exists():
                shutil.rmtree(temp_dir)

            intermediate_bams.append(current_bam)
            current_bam = vntr_biased_bam
            logger.info("  VNTR efficiency bias applied successfully")

        except Exception as e:
            logger.error("VNTR efficiency modeling failed: %s", e)
            logger.warning("Continuing with original BAM file (no bias applied)")
    else:
        logger.info("10. Skipping VNTR capture efficiency (disabled in config)")

    # Stage 11: Optionally downsample
    target_coverage = rs_config.get("coverage")
    if target_coverage:
        logger.info("11. Downsampling to target coverage")
        mode = rs_config.get("downsample_mode", "vntr").strip().lower()

        if mode == "vntr":
            vntr_region = assembly_ctx.vntr_region
            if not vntr_region:
                raise ConfigurationError(
                    f"VNTR region not specified in config for {assembly_ctx.assembly_name}. "
                    f"Add 'vntr_region_{assembly_ctx.assembly_name}' to config"
                )
            current_cov, depth_file = calculate_vntr_coverage(
                tools["samtools"], current_bam, vntr_region, threads,
                str(output_dir), output_base,
            )
            intermediate_files.append(depth_file)
            region_info = vntr_region
        elif mode == "non_vntr":
            bed_file = rs_config.get("sample_target_bed")
            if not bed_file:
                raise ConfigurationError(
                    "For non-VNTR downsampling, 'sample_target_bed' must be provided in config"
                )
            current_cov, depth_file = calculate_target_coverage(
                tools["samtools"], current_bam, bed_file, threads,
                str(output_dir), output_base,
            )
            intermediate_files.append(depth_file)
            region_info = f"BED file: {bed_file}"
        else:
            raise ConfigurationError(
                f"Invalid downsample_mode '{mode}' in config; use 'vntr' or 'non_vntr'"
            )

        if current_cov > target_coverage:
            fraction = min(max(target_coverage / current_cov, 0.0), 1.0)
            logger.info(
                "Downsampling BAM from %.2fx to %.2fx (fraction: %.4f) based on %s",
                current_cov, target_coverage, fraction, region_info,
            )
            intermediate_bams.append(current_bam)
            downsampled_bam = str(output_dir / f"{output_base}_downsampled.bam")

            if mode == "vntr":
                downsample_bam(
                    tools["samtools"], current_bam, downsampled_bam,
                    vntr_region, fraction,
                    rs_config.get("downsample_seed", 42), threads,
                )
            else:
                downsample_entire_bam(
                    tools["samtools"], current_bam, downsampled_bam,
                    fraction, rs_config.get("downsample_seed", 42), threads,
                )
            current_bam = downsampled_bam
        else:
            logger.info(
                "Current coverage (%.2fx) is below the target; no downsampling performed.",
                current_cov,
            )
    else:
        logger.info("11. Skipping downsampling (no target coverage specified)")

    return AlignmentResult(
        final_bam=current_bam,
        intermediate_bams=intermediate_bams,
        intermediate_files=intermediate_files,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_alignment_stage.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/stages/alignment.py && mypy muc_one_up/read_simulator/stages/alignment.py`

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/stages/alignment.py tests/read_simulator/test_alignment_stage.py
git commit -m "feat: extract Illumina alignment and refinement into stages/alignment.py (stages 9-11)"
```

---

### Task 6: Extract source manifest generation (Illumina stage 12a)

**Files:**
- Create: `muc_one_up/read_simulator/stages/source_manifest.py`
- Create: `tests/read_simulator/test_source_manifest.py`

- [ ] **Step 1: Write test for generate_read_manifest**

```python
"""Tests for Illumina source manifest generation."""

from pathlib import Path
from unittest.mock import MagicMock

import pytest


class TestGenerateReadManifest:
    """Test source tracking manifest generation."""

    def test_returns_none_when_no_tracker(self, tmp_path):
        """Returns None when source_tracker is None."""
        from muc_one_up.read_simulator.stages.source_manifest import (
            generate_read_manifest,
        )

        result = generate_read_manifest(
            sidecar_path=str(tmp_path / "origins.tsv"),
            final_bam=str(tmp_path / "out.bam"),
            input_fa=str(tmp_path / "input.fa"),
            source_tracker=None,
            output_dir=tmp_path,
            output_base="test",
        )
        assert result is None

    def test_returns_none_when_no_sidecar(self, tmp_path):
        """Returns None when sidecar_path is None."""
        from muc_one_up.read_simulator.stages.source_manifest import (
            generate_read_manifest,
        )

        tracker = MagicMock()
        result = generate_read_manifest(
            sidecar_path=None,
            final_bam=str(tmp_path / "out.bam"),
            input_fa=str(tmp_path / "input.fa"),
            source_tracker=tracker,
            output_dir=tmp_path,
            output_base="test",
        )
        assert result is None

    def test_generates_manifest(self, tmp_path, mocker):
        """Generates manifest TSV from sidecar + BAM."""
        mock_parse = mocker.patch(
            "muc_one_up.read_simulator.stages.source_manifest.parse_illumina_reads",
            return_value=[MagicMock(read_id="read1")],
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.source_manifest.pysam"
        )

        tracker = MagicMock()
        tracker.annotate_reads.return_value = [MagicMock()]
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">hap1\nATCG\n>hap2\nGCTA\n")

        from muc_one_up.read_simulator.stages.source_manifest import (
            generate_read_manifest,
        )

        result = generate_read_manifest(
            sidecar_path=str(tmp_path / "origins.tsv"),
            final_bam=str(tmp_path / "out.bam"),
            input_fa=str(input_fa),
            source_tracker=tracker,
            output_dir=tmp_path,
            output_base="test",
        )

        assert result is not None
        assert "read_manifest" in result
        tracker.write_manifest.assert_called_once()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/read_simulator/test_source_manifest.py -v 2>&1 | tail -5`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `source_manifest.py`**

Extract lines 575-642 from `pipeline.py`:

```python
"""Illumina source manifest generation (pipeline stage 12a).

Parses fragment origins sidecar, matches to surviving BAM reads,
annotates with haplotype info, and writes the read manifest TSV.
"""

from __future__ import annotations

import contextlib
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def _normalize_read_id(read_id: str | None) -> str:
    """Strip trailing /1 or /2 mate suffix for consistent matching."""
    if not read_id:
        return ""
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def generate_read_manifest(
    sidecar_path: str | None,
    final_bam: str,
    input_fa: str,
    source_tracker: Any | None,
    output_dir: Path,
    output_base: str,
) -> str | None:
    """Generate read source tracking manifest.

    Args:
        sidecar_path: Path to fragment_origins.tsv (None if tracking disabled).
        final_bam: Path to the final BAM file.
        input_fa: Path to the input FASTA (for haplotype mapping).
        source_tracker: ReadSourceTracker instance (None if tracking disabled).
        output_dir: Output directory for manifest.
        output_base: Base name for manifest file.

    Returns:
        Path to manifest TSV, or None if tracking is disabled.
    """
    if source_tracker is None or sidecar_path is None:
        return None

    from ..parsers.illumina_parser import parse_illumina_reads

    logger.info("Generating read source tracking manifest...")

    # Build sequence name to haplotype mapping from input FASTA
    seq_name_to_haplotype: dict[str, int] = {}
    try:
        with open(input_fa) as fa:
            hap_idx = 0
            for line in fa:
                if line.startswith(">"):
                    hap_idx += 1
                    seq_name = line.strip().lstrip(">").split()[0]
                    seq_name_to_haplotype[seq_name] = hap_idx
    except OSError:
        logger.warning("Could not read input FASTA for haplotype mapping")

    # Determine R1 FASTQ path for parser (same naming convention)
    reads_fq1 = str(output_dir / f"{output_base}_R1.fastq.gz")

    origins = parse_illumina_reads(
        sidecar_path=sidecar_path,
        fastq_r1_path=reads_fq1,
        seq_name_to_haplotype=seq_name_to_haplotype,
    )

    # Filter origins to surviving reads in final BAM
    surviving_read_ids: set[str] | None = None
    try:
        import pysam

        with pysam.AlignmentFile(final_bam, "rb") as bam:
            surviving_read_ids = {
                _normalize_read_id(read.query_name)
                for read in bam
                if read.query_name
            }
    except (ImportError, OSError):
        logger.debug("Could not read final BAM for read filtering; including all origins")

    if surviving_read_ids is not None:
        pre_filter_count = len(origins)
        origins = [o for o in origins if _normalize_read_id(o.read_id) in surviving_read_ids]
        if pre_filter_count != len(origins):
            logger.info(
                "Filtered manifest origins: %d -> %d (matched to final BAM)",
                pre_filter_count,
                len(origins),
            )

    annotated = list(source_tracker.annotate_reads(origins))
    manifest_path = str(output_dir / f"{output_base}_read_manifest.tsv.gz")
    source_tracker.write_manifest(annotated, manifest_path)
    logger.info("Read source manifest written: %s (%d reads)", manifest_path, len(annotated))

    # Clean up sidecar after successful manifest generation
    with contextlib.suppress(OSError):
        Path(sidecar_path).unlink(missing_ok=True)

    return manifest_path
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/read_simulator/test_source_manifest.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/stages/source_manifest.py tests/read_simulator/test_source_manifest.py
git commit -m "feat: extract Illumina source manifest generation into stages/source_manifest.py (stage 12a)"
```

---

### Task 7: Rewrite Illumina pipeline.py as thin orchestrator

**Files:**
- Modify: `muc_one_up/read_simulator/pipeline.py`

- [ ] **Step 1: Run existing tests to establish baseline**

Run: `pytest tests/read_simulator/test_pipeline_reference_integration.py -v`
Expected: All 3 tests PASS (these will be our regression gate).

- [ ] **Step 2: Rewrite pipeline.py**

Replace the 615-line function body with calls to the extracted stages. Keep the function signature identical. The new body should be ~80-100 lines:

```python
#!/usr/bin/env python3
"""Read simulation pipeline orchestrator.

Thin orchestrator that delegates to stage modules:
- stages.fragment_preparation: Stages 1-8 (FASTA to paired FASTQ)
- stages.alignment: Stages 9-11 (alignment, VNTR bias, downsampling)
- stages.source_manifest: Stage 12a (read provenance tracking)
- pipeline_utils: Shared output resolution, metadata, cleanup
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from .assembly_context import AssemblyContext
from .pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_human_reference,
    resolve_pipeline_outputs,
)
from .stages.alignment import align_and_refine
from .stages.fragment_preparation import prepare_fragments
from .stages.source_manifest import generate_read_manifest
from .utils import capture_tool_versions, check_external_tools, log_tool_versions

logger = logging.getLogger(__name__)


def simulate_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run the complete Illumina read simulation pipeline.

    Pipeline stages:
    1-8. Fragment preparation (Replace Ns, errors, 2bit, pblat, fragments, reads, split)
    9. Align reads to human reference (BWA MEM)
    10. Apply VNTR capture efficiency bias (optional)
    11. Downsample to target coverage (optional)
    12. Source tracking manifest + metadata + cleanup

    Args:
        config: Full configuration dictionary.
        input_fa: Input simulated FASTA file.
        source_tracker: Optional read source tracker.
        output_config: Optional output configuration.

    Returns:
        Path to the final output BAM file.
    """
    start_time = datetime.now()
    logger.info(
        "Starting read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    # Extract config sections
    tools = config.get("tools", {})
    rs_config = config.get("read_simulation", {})

    # Resolve assembly context
    assembly_ctx = AssemblyContext.from_configs(config, rs_config)

    # Validate external tools
    check_external_tools(tools)

    # Log tool versions
    logger.info("=" * 60)
    logger.info("Tool Version Information")
    logger.info("=" * 60)
    tools_to_check = ["reseq", "faToTwoBit", "samtools", "pblat", "bwa"]
    tools_subset = {k: v for k, v in tools.items() if k in tools_to_check}
    tool_versions = capture_tool_versions(tools_subset)
    log_tool_versions(tool_versions)
    logger.info("=" * 60)

    # Resolve output paths
    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_config, output_config)

    # Override FASTQ/BAM paths from config (legacy support, only when no output_config)
    if output_config is not None:
        reads_fq1 = str(output_dir / f"{output_base}_R1.fastq.gz")
        reads_fq2 = str(output_dir / f"{output_base}_R2.fastq.gz")
        output_bam_path = str(output_dir / f"{output_base}.bam")
    else:
        reads_fq1 = rs_config.get(
            "output_fastq1", str(output_dir / f"{output_base}_R1.fastq.gz")
        )
        reads_fq2 = rs_config.get(
            "output_fastq2", str(output_dir / f"{output_base}_R2.fastq.gz")
        )
        output_bam_path = rs_config.get(
            "output_bam", str(output_dir / f"{output_base}.bam")
        )

    logger.info("Output filenames:")
    logger.info("  BAM: %s", output_bam_path)
    logger.info("  FASTQ pair: %s, %s", reads_fq1, reads_fq2)

    # Stages 1-8: Fragment preparation
    frag_result = prepare_fragments(
        tools=tools,
        rs_config=rs_config,
        input_fa=input_fa,
        output_dir=output_dir,
        output_base=output_base,
        assembly_ctx=assembly_ctx,
        source_tracker=source_tracker,
    )

    # Resolve human reference (3-path fallback)
    human_ref = resolve_human_reference(config, assembly_ctx, aligner="bwa")

    # Stages 9-11: Alignment and refinement
    align_result = align_and_refine(
        tools=tools,
        rs_config=rs_config,
        r1=frag_result.r1_fastq,
        r2=frag_result.r2_fastq,
        human_ref=human_ref,
        output_dir=output_dir,
        output_base=output_base,
        assembly_ctx=assembly_ctx,
    )

    end_time = datetime.now()
    duration = end_time - start_time
    logger.info(
        "Read simulation pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],
    )
    logger.info("Final outputs:")
    logger.info("  Aligned and indexed BAM: %s", align_result.final_bam)
    logger.info("  Paired FASTQ files (gzipped): %s and %s", reads_fq1, reads_fq2)

    # Stage 12a: Source tracking manifest
    fragment_origins_path = (
        str(output_dir / f"{output_base}_fragment_origins.tsv")
        if source_tracker is not None
        else None
    )
    generate_read_manifest(
        sidecar_path=fragment_origins_path,
        final_bam=align_result.final_bam,
        input_fa=input_fa,
        source_tracker=source_tracker,
        output_dir=output_dir,
        output_base=output_base,
    )

    # Stage 12b: Metadata
    create_pipeline_metadata(
        output_dir=output_dir,
        output_base=output_base,
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform="Illumina",
        tools_used=["reseq", "faToTwoBit", "pblat", "bwa", "samtools"],
    )

    # Stage 12c: Cleanup
    keep_intermediates = rs_config.get("keep_intermediate_files", False)
    if not keep_intermediates:
        all_intermediates = frag_result.intermediate_files.copy()
        for bam in align_result.intermediate_bams:
            all_intermediates.append(bam)
            bam_index = f"{bam}.bai"
            if Path(bam_index).exists():
                all_intermediates.append(bam_index)
        all_intermediates.extend(align_result.intermediate_files)

        # Safety: never delete the final BAM or its index
        final_bam_index = f"{align_result.final_bam}.bai"
        all_intermediates = [
            f for f in all_intermediates
            if f not in (align_result.final_bam, final_bam_index)
        ]

        logger.info("=" * 60)
        logger.info("Intermediate File Cleanup")
        logger.info("=" * 60)
        logger.info("Total files to remove: %d", len(all_intermediates))
        logger.info("Preserving final output: %s", Path(align_result.final_bam).name)
        cleanup_intermediates(all_intermediates)
        logger.info("Cleanup completed successfully")
        logger.info("=" * 60)
    else:
        logger.info("=" * 60)
        logger.info("Keeping intermediate files (keep_intermediate_files=true)")
        logger.info("=" * 60)

    return align_result.final_bam


# For direct command-line use (backward compatibility)
if __name__ == "__main__":
    import json
    import sys

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

    if len(sys.argv) != 3:
        print("Usage: python pipeline.py <config.json> <input_fasta>")
        sys.exit(1)

    config_file = sys.argv[1]
    input_fa_arg = sys.argv[2]

    with Path(config_file).open() as fh:
        config = json.load(fh)

    simulate_reads_pipeline(config, input_fa_arg)
```

- [ ] **Step 3: Run existing integration tests**

Run: `pytest tests/read_simulator/test_pipeline_reference_integration.py -v`
Expected: All 3 tests PASS. If any fail, the mock paths need updating to match the new import locations.

- [ ] **Step 4: Update mock paths in existing tests if needed**

The existing tests mock at `muc_one_up.read_simulator.pipeline.replace_Ns` etc. After refactoring, these functions are no longer imported directly in `pipeline.py`. Update mock targets to point to the stage modules, or add pass-through imports. The simplest fix: mock at the wrapper level (e.g., `muc_one_up.read_simulator.wrappers.reseq_wrapper.replace_Ns`).

- [ ] **Step 5: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All existing tests PASS.

- [ ] **Step 6: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/pipeline.py && mypy muc_one_up/read_simulator/pipeline.py`

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/read_simulator/pipeline.py
git commit -m "refactor: rewrite Illumina pipeline.py as thin orchestrator delegating to stages/"
```

---

### Task 8: Refactor ONT pipeline + fix input_basename bug

**Files:**
- Modify: `muc_one_up/read_simulator/ont_pipeline.py`

- [ ] **Step 1: Run existing ONT tests as baseline**

Run: `pytest tests/read_simulator/test_pipeline_reference_integration.py::TestONTPipelineReferenceIntegration -v`
Expected: All 3 tests PASS.

- [ ] **Step 2: Refactor ont_pipeline.py**

Key changes:
1. Use `resolve_pipeline_outputs()` for output path resolution
2. Use `resolve_human_reference()` for reference resolution
3. Use `create_pipeline_metadata()` for metadata
4. Fix `input_basename` bug: derive `output_base` consistently at function top
5. Extract `_resolve_simulation_mode()` helper
6. Move late imports to function top

```python
# At the top of simulate_ont_reads_pipeline(), replace lines 157-171 with:

    # Resolve output paths (fixes input_basename bug — was undefined when output_config set)
    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_config, output_config)
    output_prefix = str(output_dir / output_base)
```

Replace lines 284-312 (reference resolution) with:

```python
    # Resolve reference for alignment
    if human_reference is None:
        try:
            reference_for_alignment = resolve_human_reference(config, assembly_ctx, aligner="minimap2")
        except ConfigurationError:
            logger.warning("Falling back to aligning against simulated reference")
            reference_for_alignment = input_fa
    else:
        reference_for_alignment = human_reference
        logger.info("Using user-provided reference: %s", reference_for_alignment)
```

Replace lines 370-379 (metadata) with:

```python
    create_pipeline_metadata(
        output_dir=output_dir,
        output_base=output_base,
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform="ONT",
        tools_used=["nanosim", "minimap2", "samtools"],
    )
```

Add the helper before the main function:

```python
def _resolve_simulation_mode(
    input_fa: str, ns_params: dict[str, Any]
) -> tuple[bool, float]:
    """Decide simulation mode and effective correction factor.

    Returns:
        (use_split_simulation, correction_factor)
    """
    is_diploid = is_diploid_reference(input_fa)
    enable_split = ns_params.get("enable_split_simulation", True)
    correction_factor = ns_params.get("correction_factor", 0.325)
    return (is_diploid and enable_split, correction_factor)
```

- [ ] **Step 3: Run ONT tests**

Run: `pytest tests/read_simulator/test_pipeline_reference_integration.py::TestONTPipelineReferenceIntegration -v`
Expected: All 3 tests PASS.

- [ ] **Step 4: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/ont_pipeline.py && mypy muc_one_up/read_simulator/ont_pipeline.py`

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/ont_pipeline.py
git commit -m "refactor: ONT pipeline uses shared utils, fix input_basename bug, extract _resolve_simulation_mode"
```

---

### Task 9: Refactor PacBio pipeline

**Files:**
- Modify: `muc_one_up/read_simulator/pacbio_pipeline.py`

- [ ] **Step 1: Refactor pacbio_pipeline.py**

Key changes:
1. Use `resolve_pipeline_outputs()` for output paths
2. Use `create_pipeline_metadata()` for metadata
3. Use `cleanup_intermediates()` as leaf cleanup (keep `finally` block and MAF preservation logic)
4. Extract `_simulate_haplotype_hifi()` helper

Add before the main function:

```python
def _simulate_haplotype_hifi(
    clr_bam: str,
    ccs_cmd: str,
    output_path: str,
    min_passes: int,
    min_rq: float,
    threads: int,
    seed: int | None,
    hap_idx: int,
) -> str:
    """Run CCS consensus on a single haplotype's CLR BAM.

    Args:
        clr_bam: Path to CLR BAM for this haplotype.
        ccs_cmd: CCS command path.
        output_path: Output HiFi BAM path.
        min_passes: Minimum CCS passes.
        min_rq: Minimum read quality.
        threads: Thread count.
        seed: Base seed (hap_idx added for uniqueness).
        hap_idx: 1-based haplotype index.

    Returns:
        Path to HiFi BAM.
    """
    hifi_seed = seed + hap_idx if seed is not None else None
    hifi_bam = run_ccs_consensus(
        ccs_cmd=ccs_cmd,
        input_bam=clr_bam,
        output_bam=output_path,
        min_passes=min_passes,
        min_rq=min_rq,
        threads=threads,
        seed=hifi_seed,
    )
    logging.info("Haplotype %d CCS complete: %s", hap_idx, hifi_bam)
    return hifi_bam
```

Replace lines 244-254 (output resolution) with:

```python
    output_dir_path, output_base = resolve_pipeline_outputs(input_fa, {}, output_config)
```

Replace lines 456-466 (metadata) with:

```python
        create_pipeline_metadata(
            output_dir=output_dir_path,
            output_base=output_base,
            config=config,
            start_time=start_time,
            end_time=end_time,
            platform="PacBio",
            tools_used=["pbsim3", "ccs", "minimap2", "samtools"],
        )
```

Replace the CCS loop (lines 324-345) with calls to `_simulate_haplotype_hifi()`.

Replace `cleanup_files` calls in `finally` block with `cleanup_intermediates`.

- [ ] **Step 2: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests PASS.

- [ ] **Step 3: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/pacbio_pipeline.py && mypy muc_one_up/read_simulator/pacbio_pipeline.py`

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/read_simulator/pacbio_pipeline.py
git commit -m "refactor: PacBio pipeline uses shared utils, extract _simulate_haplotype_hifi"
```

---

### Task 10: Refactor Amplicon pipeline

**Files:**
- Modify: `muc_one_up/read_simulator/amplicon_pipeline.py`

- [ ] **Step 1: Refactor amplicon_pipeline.py**

Key changes:
1. Use `resolve_pipeline_outputs()` for output paths
2. Use `create_pipeline_metadata()` for metadata
3. Use `cleanup_intermediates()` as leaf cleanup (keep `finally` block and tempdir)
4. Extract `_simulate_haplotype_amplicon()` helper

Add before the main function:

```python
def _simulate_haplotype_amplicon(
    clr_bam: str,
    ccs_cmd: str,
    output_path: str,
    min_passes: int,
    min_rq: float,
    threads: int,
    seed: int | None,
    hap_idx: int,
    bam_idx: int,
) -> str:
    """Run CCS consensus on a single amplicon CLR BAM.

    Args:
        clr_bam: Path to CLR BAM.
        ccs_cmd: CCS command path.
        output_path: Output HiFi BAM path.
        min_passes: Minimum CCS passes.
        min_rq: Minimum read quality.
        threads: Thread count.
        seed: Base seed (hap_idx*100 + bam_idx added for uniqueness).
        hap_idx: 1-based haplotype index.
        bam_idx: 1-based BAM index within haplotype.

    Returns:
        Path to HiFi BAM.
    """
    hifi_seed = (seed + hap_idx * 100 + bam_idx) if seed is not None else None
    return run_ccs_consensus(
        ccs_cmd=ccs_cmd,
        input_bam=clr_bam,
        output_bam=output_path,
        min_passes=min_passes,
        min_rq=min_rq,
        threads=threads,
        seed=hifi_seed,
    )
```

Replace lines 109-117 (output resolution) with:

```python
    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_config, output_config)
```

Replace CCS loop (lines 229-243) with calls to `_simulate_haplotype_amplicon()`.

Replace `cleanup_files` calls with `cleanup_intermediates`.

Replace `write_metadata_file` call with `create_pipeline_metadata`.

- [ ] **Step 2: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests PASS.

- [ ] **Step 3: Lint and type check**

Run: `ruff check muc_one_up/read_simulator/amplicon_pipeline.py && mypy muc_one_up/read_simulator/amplicon_pipeline.py`

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/read_simulator/amplicon_pipeline.py
git commit -m "refactor: Amplicon pipeline uses shared utils, extract _simulate_haplotype_amplicon"
```

---

### Task 11: CLI dispatch consolidation

**Files:**
- Modify: `muc_one_up/cli/commands/reads.py`
- Create: `tests/cli/test_apply_pacbio_params.py`

- [ ] **Step 1: Write test for _apply_pacbio_params**

```python
"""Tests for PacBio parameter consolidation in CLI."""

import click
import pytest


class TestApplyPacbioParams:
    """Test the shared PacBio parameter helper."""

    def test_initializes_pacbio_params(self):
        """Creates pacbio_params section if missing."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        _apply_pacbio_params(config, model_type="qshmm", model_file="/model", seed=None)
        assert "pacbio_params" in config
        assert config["pacbio_params"]["model_type"] == "qshmm"

    def test_only_overrides_non_none(self):
        """Only sets params that are not None."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config = {"pacbio_params": {"model_type": "errhmm", "threads": 8}}
        _apply_pacbio_params(
            config, model_type=None, model_file="/model", seed=None, threads=None
        )
        assert config["pacbio_params"]["model_type"] == "errhmm"  # not overwritten
        assert config["pacbio_params"]["model_file"] == "/model"
        assert config["pacbio_params"]["threads"] == 8  # not overwritten

    def test_raises_on_missing_required(self):
        """Raises ClickException when required params missing."""
        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        with pytest.raises(click.ClickException, match="model_type"):
            _apply_pacbio_params(config, model_type=None, model_file=None, seed=None)

    def test_logs_seed(self, caplog):
        """Logs seed when provided."""
        import logging

        from muc_one_up.cli.commands.reads import _apply_pacbio_params

        config: dict = {}
        with caplog.at_level(logging.INFO):
            _apply_pacbio_params(
                config, model_type="qshmm", model_file="/m", seed=42
            )
        assert "42" in caplog.text
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/cli/test_apply_pacbio_params.py -v 2>&1 | tail -5`
Expected: FAIL — `_apply_pacbio_params` doesn't exist yet.

- [ ] **Step 3: Implement `_apply_pacbio_params` and update commands**

Add after `_setup_read_config` in `reads.py`:

```python
def _apply_pacbio_params(
    config: dict[str, Any],
    model_type: str | None,
    model_file: str | None,
    seed: int | None,
    threads: int | None = None,
    pass_num: int | None = None,
    min_passes: int | None = None,
    min_rq: float | None = None,
) -> None:
    """Apply PacBio CLI overrides and validate required params.

    Used by both pacbio() and amplicon() commands to avoid duplication.
    Mutates config in place.

    Raises:
        click.ClickException: If required params are missing after merge.
    """
    if "pacbio_params" not in config:
        config["pacbio_params"] = {}

    params = config["pacbio_params"]
    overrides = {
        "model_type": model_type,
        "model_file": model_file,
        "pass_num": pass_num,
        "min_passes": min_passes,
        "min_rq": min_rq,
    }
    if threads is not None:
        overrides["threads"] = threads

    for key, value in overrides.items():
        if value is not None:
            params[key] = value

    if seed is not None:
        params["seed"] = seed
        logging.info("Using random seed: %d (results will be reproducible)", seed)

    # Validate required params
    required = ["model_type", "model_file"]
    missing = [p for p in required if p not in params]
    if missing:
        raise click.ClickException(
            f"Missing required PacBio parameters: {', '.join(missing)}. "
            f"Provide via CLI options or config.json pacbio_params section."
        )
```

Then update the `pacbio()` command (lines 362-396) to:

```python
    _setup_read_config(config, "pacbio", coverage, seed, seed_config_key="pacbio_params")
    _apply_pacbio_params(config, model_type, model_file, seed, threads=threads if threads != 4 else None, pass_num=pass_num, min_passes=min_passes, min_rq=min_rq)
```

And update the `amplicon()` command (lines 501-520) to:

```python
    _apply_pacbio_params(config, model_type, model_file, seed)
```

- [ ] **Step 4: Run tests**

Run: `pytest tests/cli/test_apply_pacbio_params.py -v`
Expected: All tests PASS.

- [ ] **Step 5: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests PASS.

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/cli/commands/reads.py tests/cli/test_apply_pacbio_params.py
git commit -m "refactor: consolidate PacBio CLI param setup into _apply_pacbio_params"
```

---

### Task 12: Contract tests for output naming and return values

**Files:**
- Create: `tests/read_simulator/test_pipeline_contracts.py`

- [ ] **Step 1: Write contract tests**

```python
"""Contract tests for pipeline output naming, metadata, and return values.

These tests mock all external tool calls and validate only the Python-level
orchestration logic: output file paths, metadata basename, and return values.
They catch regressions in output naming when refactoring pipeline internals.
"""

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from muc_one_up.read_simulator.output_config import OutputConfig


class TestIlluminaOutputNaming:
    """Illumina pipeline output naming contracts."""

    def _make_config(self, tmp_path):
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\nATCG\n")
        for ext in [".fai", ".amb", ".ann", ".bwt", ".pac", ".sa"]:
            (tmp_path / f"ref.fa{ext}").touch()
        return {
            "tools": {
                "reseq": "reseq", "faToTwoBit": "faToTwoBit",
                "samtools": "samtools", "pblat": "pblat", "bwa": "bwa",
            },
            "read_simulation": {
                "reseq_model": str(tmp_path / "model"),
                "sample_bam": str(tmp_path / "sample.bam"),
                "output_dir": str(tmp_path),
                "vntr_capture_efficiency": {"enabled": False},
            },
            "reference_assembly": "hg38",
            "reference_genomes": {"hg38": {"fasta_path": str(ref)}},
        }

    def _mock_all_stages(self, mocker, tmp_path):
        mocker.patch("muc_one_up.read_simulator.pipeline.check_external_tools")
        mocker.patch("muc_one_up.read_simulator.pipeline.capture_tool_versions", return_value={})
        mocker.patch("muc_one_up.read_simulator.pipeline.log_tool_versions")
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.replace_Ns"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.generate_systematic_errors"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.fa_to_twobit"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.extract_subset_reference",
            return_value=str(tmp_path / "_collated.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.run_pblat"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.simulate_fragments"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.create_reads"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.fragment_preparation.split_reads"
        )
        mocker.patch(
            "muc_one_up.read_simulator.stages.alignment.align_reads"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pipeline.cleanup_intermediates"
        )

    def test_output_bam_with_output_config(self, tmp_path, mocker):
        """BAM path respects output_config base name."""
        config = self._make_config(tmp_path)
        self._mock_all_stages(mocker, tmp_path)

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">seq\nATCG\n")
        oc = OutputConfig(out_dir=tmp_path / "custom", out_base="myreads")

        result = simulate_reads_pipeline(config, str(input_fa), output_config=oc)
        assert Path(result).name == "myreads.bam"
        assert Path(result).parent == tmp_path / "custom"

    def test_output_bam_without_output_config(self, tmp_path, mocker):
        """BAM path uses input stem when no output_config."""
        config = self._make_config(tmp_path)
        self._mock_all_stages(mocker, tmp_path)

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "sample_001.fa"
        input_fa.write_text(">seq\nATCG\n")

        result = simulate_reads_pipeline(config, str(input_fa))
        assert Path(result).name == "sample_001.bam"

    def test_returns_string_path(self, tmp_path, mocker):
        """Return value is a string (not Path)."""
        config = self._make_config(tmp_path)
        self._mock_all_stages(mocker, tmp_path)

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">seq\nATCG\n")

        result = simulate_reads_pipeline(config, str(input_fa))
        assert isinstance(result, str)


class TestONTOutputNaming:
    """ONT pipeline output naming contracts — regression for input_basename bug."""

    def test_output_with_output_config(self, tmp_path, mocker):
        """Metadata uses output_config base, not input_basename."""
        config = {
            "tools": {"nanosim": "sim", "minimap2": "mm2", "samtools": "st"},
            "nanosim_params": {"training_data_path": "/model", "coverage": 30},
            "read_simulation": {},
        }
        mock_sim = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.run_nanosim_simulation",
            return_value=str(tmp_path / "reads.fq"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.align_ont_reads_with_minimap2"
        )
        mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.is_diploid_reference",
            return_value=False,
        )
        mock_meta = mocker.patch(
            "muc_one_up.read_simulator.ont_pipeline.create_pipeline_metadata"
        )

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">seq\nATCG\n")
        oc = OutputConfig(out_dir=tmp_path / "out", out_base="custom_ont")

        simulate_ont_reads_pipeline(config, str(input_fa), output_config=oc)

        # Verify metadata uses output_config base, not "test_ont"
        mock_meta.assert_called_once()
        call_kwargs = mock_meta.call_args
        assert call_kwargs[1]["output_base"] == "custom_ont" or call_kwargs[0][1] == "custom_ont"


class TestPacBioOutputNaming:
    """PacBio pipeline output naming contracts."""

    def test_returns_bam_with_reference(self, tmp_path, mocker):
        """Returns BAM path when human_reference is provided."""
        config = {
            "tools": {"pbsim3": "pb", "ccs": "ccs", "samtools": "st", "minimap2": "mm2"},
            "pacbio_params": {
                "model_type": "qshmm", "model_file": "/m",
                "coverage": 5, "pass_num": 3, "min_passes": 3, "min_rq": 0.99,
            },
        }
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.validate_pbsim3_parameters"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.validate_ccs_parameters"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(tmp_path / "clr.bam")],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(tmp_path / "hifi.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(tmp_path / "hifi.fq"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.align_reads_with_minimap2",
            return_value=str(tmp_path / "aligned.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.cleanup_intermediates"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.create_pipeline_metadata"
        )

        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">seq\nATCG\n")

        result = simulate_pacbio_hifi_reads(
            config, str(input_fa), human_reference="/ref.fa"
        )
        assert result.endswith(".bam")

    def test_returns_fastq_without_reference(self, tmp_path, mocker):
        """Returns FASTQ path when human_reference is None."""
        config = {
            "tools": {"pbsim3": "pb", "ccs": "ccs", "samtools": "st", "minimap2": "mm2"},
            "pacbio_params": {
                "model_type": "qshmm", "model_file": "/m",
                "coverage": 5, "pass_num": 3, "min_passes": 3, "min_rq": 0.99,
            },
        }
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.validate_pbsim3_parameters"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.validate_ccs_parameters"
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation",
            return_value=[str(tmp_path / "clr.bam")],
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus",
            return_value=str(tmp_path / "hifi.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq",
            return_value=str(tmp_path / "hifi.fq"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.pacbio_pipeline.cleanup_intermediates"
        )

        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">seq\nATCG\n")

        result = simulate_pacbio_hifi_reads(
            config, str(input_fa), human_reference=None
        )
        assert result.endswith(".fastq")


class TestAmpliconOutputNaming:
    """Amplicon pipeline output naming contracts."""

    def test_returns_bam_with_reference(self, tmp_path, mocker):
        """Returns BAM path when human_reference provided."""
        config = {
            "tools": {"pbsim3": "pb", "ccs": "ccs", "samtools": "st", "minimap2": "mm2"},
            "pacbio_params": {"model_type": "qshmm", "model_file": "/m"},
            "amplicon_params": {
                "forward_primer": "ATCG",
                "reverse_primer": "GCTA",
            },
            "read_simulation": {"coverage": 5},
        }
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.is_diploid_reference",
            return_value=False,
        )
        mock_extract = MagicMock()
        mock_extract.extract.return_value = MagicMock(length=500, fasta_path=str(tmp_path / "amp.fa"))
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.AmpliconExtractor",
            return_value=mock_extract,
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.generate_template_fasta"
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.run_pbsim3_template_simulation",
            return_value=[str(tmp_path / "clr.bam")],
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.run_ccs_consensus",
            return_value=str(tmp_path / "hifi.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.merge_bam_files",
            return_value=str(tmp_path / "merged.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.convert_bam_to_fastq",
            return_value=str(tmp_path / "hifi.fq"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.align_reads_with_minimap2",
            return_value=str(tmp_path / "aligned.bam"),
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.cleanup_intermediates"
        )
        mocker.patch(
            "muc_one_up.read_simulator.amplicon_pipeline.create_pipeline_metadata"
        )

        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        input_fa = tmp_path / "test.fa"
        input_fa.write_text(">seq\nATCG\n")

        result = simulate_amplicon_reads_pipeline(
            config, str(input_fa), human_reference="/ref.fa"
        )
        assert result.endswith(".bam")
```

- [ ] **Step 2: Run tests**

Run: `pytest tests/read_simulator/test_pipeline_contracts.py -v`
Expected: All tests PASS (if run after Tasks 7-10; if run before, update mock paths).

- [ ] **Step 3: Commit**

```bash
git add tests/read_simulator/test_pipeline_contracts.py
git commit -m "test: add contract tests for pipeline output naming, metadata, and return values"
```

---

### Task 13: E2E test fixtures and tests

**Files:**
- Create: `tests/e2e/__init__.py`
- Create: `tests/e2e/conftest.py`
- Create: `tests/e2e/test_pipeline_e2e.py`

- [ ] **Step 1: Create E2E package and conftest**

`tests/e2e/__init__.py`: empty file.

`tests/e2e/conftest.py`:

```python
"""Shared fixtures for E2E pipeline tests.

These tests require real bioinformatics tools installed locally.
They are gated by @pytest.mark.e2e and @pytest.mark.requires_tools.
"""

import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from muc_one_up.simulate import simulate_from_chains
from muc_one_up.type_defs import RepeatUnit


@pytest.fixture
def e2e_config_base() -> dict[str, Any]:
    """Base config with tool paths discovered from PATH."""
    tools = {}
    for tool in ["samtools", "bwa", "minimap2", "reseq", "faToTwoBit", "pblat",
                  "pbsim3", "ccs", "nanosim"]:
        path = shutil.which(tool)
        if path:
            tools[tool] = path
    # pbsim3 may be installed as "pbsim"
    if "pbsim3" not in tools:
        pbsim = shutil.which("pbsim")
        if pbsim:
            tools["pbsim3"] = pbsim

    return {"tools": tools}


@pytest.fixture
def e2e_diploid_fasta(tmp_path, minimal_config) -> Path:
    """Generate a small diploid FASTA from real MUC1 repeat chains.

    Produces ~300-500bp per haplotype using chains 1-2-X-B-6-7-8-9.
    """
    chains = [
        [RepeatUnit.from_str(s) for s in ["1", "2", "X", "B", "6", "7", "8", "9"]],
        [RepeatUnit.from_str(s) for s in ["1", "2", "A", "B", "6p", "7", "8", "9"]],
    ]
    results = simulate_from_chains(chains, minimal_config)

    fasta_path = tmp_path / "e2e_diploid.fa"
    with open(fasta_path, "w") as f:
        for hr in results:
            f.write(f">{hr.name}\n{hr.sequence}\n")

    return fasta_path


@pytest.fixture
def e2e_simulation_stats(tmp_path, e2e_diploid_fasta) -> Path:
    """Write minimal simulation_stats.json companion file."""
    stats = {
        "haplotypes": [
            {"name": "haplotype_1", "chain": "1-2-X-B-6-7-8-9"},
            {"name": "haplotype_2", "chain": "1-2-A-B-6p-7-8-9"},
        ],
        "assembly": "hg38",
    }
    stats_path = tmp_path / "e2e_diploid.simulation_stats.json"
    stats_path.write_text(json.dumps(stats))
    return stats_path
```

- [ ] **Step 2: Write E2E tests**

```python
"""End-to-end pipeline tests with real bioinformatics tools.

These tests simulate a VNTR repeat and produce BAM/FASTQ files from each
pipeline variant. They require real tools installed locally and are gated
by @pytest.mark.e2e and @pytest.mark.requires_tools.

Run with: pytest tests/e2e/ -m e2e -v
"""

import json
from pathlib import Path

import pytest


@pytest.mark.e2e
@pytest.mark.requires_tools("reseq", "bwa", "samtools", "pblat", "faToTwoBit")
class TestIlluminaE2E:
    """Full Illumina pipeline: FASTA -> paired FASTQ -> sorted BAM."""

    def test_illumina_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        """Illumina pipeline produces a valid BAM with reads."""
        import pysam

        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline

        # Build config
        config = {**e2e_config_base, **minimal_config}
        config["read_simulation"] = {
            "simulator": "illumina",
            "reseq_model": _find_reseq_model(),
            "sample_bam": _find_sample_bam(),
            "human_reference": str(e2e_diploid_fasta),  # Align against self
            "output_dir": str(tmp_path / "illumina_out"),
            "read_number": 100,
            "fragment_size": 250,
            "fragment_sd": 30,
            "min_fragment": 150,
            "threads": 2,
            "vntr_capture_efficiency": {"enabled": False},
        }

        result = simulate_reads_pipeline(config, str(e2e_diploid_fasta))

        # Assertions
        assert Path(result).exists(), "BAM file should exist"
        assert Path(f"{result}.bai").exists(), "BAM index should exist"

        with pysam.AlignmentFile(result, "rb") as bam:
            read_count = sum(1 for _ in bam)
        assert read_count > 0, "BAM should contain reads"

        # Metadata file should exist
        metadata_files = list(Path(result).parent.glob("*_metadata.tsv"))
        assert len(metadata_files) > 0, "Metadata TSV should be written"


@pytest.mark.e2e
@pytest.mark.requires_tools("nanosim", "minimap2", "samtools")
class TestONTE2E:
    """Full ONT pipeline: FASTA -> FASTQ -> sorted BAM."""

    def test_ont_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        """ONT pipeline produces a valid BAM with reads."""
        import pysam

        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline

        config = {**e2e_config_base, **minimal_config}
        config["nanosim_params"] = {
            "training_data_path": _find_nanosim_model(),
            "coverage": 5,
            "num_threads": 2,
            "seed": 42,
            "enable_split_simulation": True,
        }
        config["read_simulation"] = {
            "simulator": "ont",
            "output_dir": str(tmp_path / "ont_out"),
        }

        result = simulate_ont_reads_pipeline(
            config, str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "BAM file should exist"
        with pysam.AlignmentFile(result, "rb") as bam:
            read_count = sum(1 for _ in bam)
        assert read_count > 0, "BAM should contain reads"


@pytest.mark.e2e
@pytest.mark.requires_tools("pbsim3", "ccs", "minimap2", "samtools")
class TestPacBioE2E:
    """Full PacBio HiFi pipeline: FASTA -> CLR -> CCS -> BAM."""

    def test_pacbio_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta
    ):
        """PacBio pipeline produces a valid BAM with HiFi reads."""
        import pysam

        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_hifi_reads

        config = {**e2e_config_base}
        config["pacbio_params"] = {
            "model_type": "qshmm",
            "model_file": _find_pbsim3_model(),
            "coverage": 5,
            "pass_num": 3,
            "min_passes": 3,
            "min_rq": 0.99,
            "threads": 2,
            "seed": 42,
        }

        result = simulate_pacbio_hifi_reads(
            config, str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "BAM file should exist"
        with pysam.AlignmentFile(result, "rb") as bam:
            read_count = sum(1 for _ in bam)
        assert read_count > 0, "BAM should contain reads"


@pytest.mark.e2e
@pytest.mark.requires_tools("pbsim3", "ccs", "minimap2", "samtools")
class TestAmpliconE2E:
    """Full Amplicon pipeline: FASTA -> template -> CLR -> CCS -> BAM."""

    def test_amplicon_produces_valid_bam(
        self, tmp_path, e2e_config_base, e2e_diploid_fasta, minimal_config
    ):
        """Amplicon pipeline produces a valid BAM."""
        import pysam

        from muc_one_up.read_simulator.amplicon_pipeline import (
            simulate_amplicon_reads_pipeline,
        )

        config = {**e2e_config_base, **minimal_config}
        config["pacbio_params"] = {
            "model_type": "qshmm",
            "model_file": _find_pbsim3_model(),
            "threads": 2,
            "seed": 42,
        }
        config["amplicon_params"] = {
            "forward_primer": minimal_config["constants"]["hg38"]["left"][:20],
            "reverse_primer": _reverse_complement(
                minimal_config["constants"]["hg38"]["right"][-20:]
            ),
        }
        config["read_simulation"] = {
            "simulator": "amplicon",
            "coverage": 10,
        }

        result = simulate_amplicon_reads_pipeline(
            config, str(e2e_diploid_fasta),
            human_reference=str(e2e_diploid_fasta),
        )

        assert Path(result).exists(), "BAM file should exist"
        with pysam.AlignmentFile(result, "rb") as bam:
            read_count = sum(1 for _ in bam)
        assert read_count > 0, "BAM should contain reads"


# ============================================================================
# Helper functions for locating test resources
# ============================================================================

def _find_reseq_model() -> str:
    """Find a ReSeq model file for testing. Raises pytest.skip if not found."""
    import os
    model_path = os.environ.get("MUCONEUP_RESEQ_MODEL")
    if model_path and Path(model_path).exists():
        return model_path
    pytest.skip("MUCONEUP_RESEQ_MODEL environment variable not set")


def _find_sample_bam() -> str:
    """Find a sample BAM for testing. Raises pytest.skip if not found."""
    import os
    bam_path = os.environ.get("MUCONEUP_SAMPLE_BAM")
    if bam_path and Path(bam_path).exists():
        return bam_path
    pytest.skip("MUCONEUP_SAMPLE_BAM environment variable not set")


def _find_nanosim_model() -> str:
    """Find NanoSim training model. Raises pytest.skip if not found."""
    import os
    model_path = os.environ.get("MUCONEUP_NANOSIM_MODEL")
    if model_path and Path(model_path).exists():
        return model_path
    pytest.skip("MUCONEUP_NANOSIM_MODEL environment variable not set")


def _find_pbsim3_model() -> str:
    """Find pbsim3 model file. Raises pytest.skip if not found."""
    import os
    model_path = os.environ.get("MUCONEUP_PBSIM3_MODEL")
    if model_path and Path(model_path).exists():
        return model_path
    pytest.skip("MUCONEUP_PBSIM3_MODEL environment variable not set")


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement.get(b, b) for b in reversed(seq))
```

- [ ] **Step 3: Verify tests are collected but skipped without tools**

Run: `pytest tests/e2e/ --co -q`
Expected: 4 tests collected.

Run: `pytest tests/e2e/ -v -m e2e`
Expected: Tests are either SKIPPED (tools not found) or PASS (tools available).

- [ ] **Step 4: Commit**

```bash
git add tests/e2e/
git commit -m "test: add E2E tests for all 4 pipeline variants with real tool execution"
```

---

### Task 14: Final verification

- [ ] **Step 1: Run full test suite**

Run: `pytest --tb=short -q`
Expected: All tests PASS (E2E tests skipped if tools not available).

- [ ] **Step 2: Run linter on all changed files**

Run: `ruff check muc_one_up/read_simulator/ muc_one_up/cli/commands/reads.py tests/`

- [ ] **Step 3: Run type checker**

Run: `mypy muc_one_up/`

- [ ] **Step 4: Verify line count reduction**

Run: `wc -l muc_one_up/read_simulator/pipeline.py muc_one_up/read_simulator/ont_pipeline.py muc_one_up/read_simulator/pacbio_pipeline.py muc_one_up/read_simulator/amplicon_pipeline.py`
Expected: pipeline.py ~100 LOC (down from 711), others reduced proportionally.

- [ ] **Step 5: Commit any final fixes**

```bash
git add -A
git commit -m "chore: final cleanup after pipeline decomposition"
```
