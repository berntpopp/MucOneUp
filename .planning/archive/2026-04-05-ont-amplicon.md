# ONT Amplicon Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add ONT amplicon simulation via `reads amplicon --platform ont`, sharing extraction/PCR-bias stages with the existing PacBio amplicon pipeline.

**Architecture:** Extract shared amplicon stages into `amplicon_common.py`, refactor existing PacBio pipeline to use it, create new `ont_amplicon_pipeline.py` that composes shared stages with single-pass pbsim3 + minimap2 map-ont. Add `--platform` CLI flag, `"ont-amplicon"` to dispatch and config schema, and `assay_type` metadata field.

**Tech Stack:** Python 3.10+, Click CLI, pbsim3, samtools, minimap2, pytest

---

### Task 1: Relax pbsim3 wrapper pass_num validation

The wrapper currently rejects `pass_num < 2`. ONT needs `pass_num = 1`.

**Files:**
- Modify: `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py:382-386`
- Test: `tests/read_simulator/test_pbsim3_wrapper.py` (extend if exists, or inline)

- [ ] **Step 1: Write the failing test**

Find or create a test that verifies `pass_num=1` is accepted. Add to the end of the relevant test file. If no dedicated pbsim3 wrapper test file exists, add inline:

```python
# In tests/read_simulator/test_pbsim3_wrapper.py (create if needed)
"""Tests for pbsim3 wrapper validation."""

import pytest

from muc_one_up.exceptions import FileOperationError
from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import run_pbsim3_template_simulation


def test_pass_num_1_accepted(tmp_path):
    """pass_num=1 must be accepted for ONT single-pass simulation."""
    model = tmp_path / "test.model"
    model.write_text("model")
    template = tmp_path / "template.fa"
    template.write_text(">seq\nACGT\n")

    # Should not raise FileOperationError for pass_num validation.
    # It WILL raise for missing pbsim3 executable — that's expected.
    with pytest.raises(Exception) as exc_info:
        run_pbsim3_template_simulation(
            pbsim3_cmd="nonexistent_pbsim3",
            samtools_cmd="samtools",
            template_fasta=str(template),
            model_type="errhmm",
            model_file=str(model),
            output_prefix=str(tmp_path / "out"),
            pass_num=1,
        )
    # The error should NOT be about pass_num validation
    assert "pass_num" not in str(exc_info.value).lower()


def test_pass_num_0_rejected(tmp_path):
    """pass_num=0 must still be rejected."""
    model = tmp_path / "test.model"
    model.write_text("model")
    template = tmp_path / "template.fa"
    template.write_text(">seq\nACGT\n")

    with pytest.raises(FileOperationError, match="pass_num"):
        run_pbsim3_template_simulation(
            pbsim3_cmd="pbsim",
            samtools_cmd="samtools",
            template_fasta=str(template),
            model_type="errhmm",
            model_file=str(model),
            output_prefix=str(tmp_path / "out"),
            pass_num=0,
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/read_simulator/test_pbsim3_wrapper.py::test_pass_num_1_accepted -v`
Expected: FAIL — `FileOperationError: Invalid pass_num: 1. Multi-pass simulation requires pass_num >= 2`

- [ ] **Step 3: Fix the validation**

In `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py`, change lines 382-386 from:

```python
    if pass_num < 2:
        raise FileOperationError(
            f"Invalid pass_num: {pass_num}. Multi-pass simulation requires pass_num >= 2"
        )
```

to:

```python
    if pass_num < 1:
        raise FileOperationError(
            f"Invalid pass_num: {pass_num}. Template simulation requires pass_num >= 1"
        )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/read_simulator/test_pbsim3_wrapper.py -v`
Expected: Both tests pass.

- [ ] **Step 5: Run full suite**

Run: `python -m pytest --tb=short -q`

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py tests/read_simulator/test_pbsim3_wrapper.py
git commit -m "fix: allow pass_num=1 in pbsim3 template mode for ONT

The validation rejected pass_num < 2, but ONT amplicon simulation
needs single-pass (pass_num=1). Relaxed to pass_num >= 1."
```

---

### Task 2: Add assay_type support to metadata writer

**Files:**
- Modify: `muc_one_up/read_simulator/utils/metadata_writer.py:164-165`
- Modify: `tests/read_simulator/test_metadata_writer.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/read_simulator/test_metadata_writer.py`:

```python
    def test_assay_type_written_when_present(
        self, mock_log_versions, mock_capture_versions, tmp_path
    ):
        """Assay_type field from config is written to metadata."""
        config = {
            "tools": {"samtools": "samtools"},
            "read_simulation": {"coverage": 30, "assay_type": "amplicon"},
        }
        path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=datetime(2025, 1, 1, 12, 0, 0),
            end_time=datetime(2025, 1, 1, 13, 0, 0),
            platform="PacBio",
            tools_used=["samtools"],
        )
        content = Path(path).read_text()
        assert "Assay_type\tamplicon" in content

    def test_assay_type_omitted_when_absent(
        self, mock_log_versions, mock_capture_versions, tmp_path
    ):
        """Assay_type field is not written when not in config."""
        config = {
            "tools": {"samtools": "samtools"},
            "read_simulation": {"coverage": 30},
        }
        path = write_metadata_file(
            output_dir=str(tmp_path),
            output_base="test",
            config=config,
            start_time=datetime(2025, 1, 1, 12, 0, 0),
            end_time=datetime(2025, 1, 1, 13, 0, 0),
            platform="Illumina",
            tools_used=["samtools"],
        )
        content = Path(path).read_text()
        assert "Assay_type" not in content
```

Note: Check the existing test file for the correct imports (`datetime`, `Path`) and fixtures (`mock_log_versions`, `mock_capture_versions`). Add any missing imports.

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/read_simulator/test_metadata_writer.py::TestMetadataWriter::test_assay_type_written_when_present -v`
Expected: FAIL — `Assay_type` not in output.

- [ ] **Step 3: Add assay_type writing to metadata_writer.py**

In `muc_one_up/read_simulator/utils/metadata_writer.py`, after the platform-specific block (after line 164, before the tool versions loop), add:

```python
        # Assay type (e.g., "amplicon") — written when present in config
        assay = config.get("read_simulation", {}).get("assay_type")
        if assay:
            f.write(f"Assay_type\t{assay}\n")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/read_simulator/test_metadata_writer.py -v`
Expected: All tests pass including new ones.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/utils/metadata_writer.py tests/read_simulator/test_metadata_writer.py
git commit -m "feat: add assay_type field to simulation metadata

Writes Assay_type (e.g., 'amplicon') to metadata TSV when present
in config['read_simulation']['assay_type']. Used by both PacBio
and ONT amplicon pipelines."
```

---

### Task 3: Extract shared amplicon helpers

Extract stages 1-3 from `amplicon_pipeline.py` into reusable `amplicon_common.py`.

**Files:**
- Create: `muc_one_up/read_simulator/amplicon_common.py`
- Modify: `muc_one_up/read_simulator/amplicon_pipeline.py`
- Create: `tests/read_simulator/test_amplicon_common.py`

- [ ] **Step 1: Write the test**

Create `tests/read_simulator/test_amplicon_common.py`:

```python
"""Tests for shared amplicon extraction and preparation."""

from pathlib import Path

import pytest

from muc_one_up.read_simulator.amplicon_common import AmpliconPrep


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def diploid_fasta_with_primers(tmp_path, muc1_primers):
    """Create a diploid FASTA with primer sites in both haplotypes."""
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())

    vntr1 = "ACGT" * 400
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50

    vntr2 = "ACGT" * 800
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50

    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


class TestAmpliconPrep:
    def test_dataclass_fields(self):
        """AmpliconPrep has required fields."""
        prep = AmpliconPrep(
            allele_templates=[Path("/a.fa")],
            allele_coverages=[100],
            output_dir=Path("/out"),
            output_base="test",
            intermediate_files=[],
            is_diploid=False,
        )
        assert prep.is_diploid is False
        assert len(prep.allele_templates) == 1


class TestExtractAndPrepare:
    def test_diploid_returns_two_templates(
        self, diploid_fasta_with_primers, tmp_path, muc1_primers
    ):
        """Diploid input produces 2 template FASTAs with PCR bias split."""
        from muc_one_up.read_simulator.amplicon_common import extract_and_prepare_amplicons

        prep = extract_and_prepare_amplicons(
            input_fa=str(diploid_fasta_with_primers),
            forward_primer=muc1_primers["forward"],
            reverse_primer=muc1_primers["reverse"],
            total_coverage=100,
            work_dir=tmp_path,
            expected_product_range=None,
            pcr_bias_config={},
            seed=42,
        )

        assert prep.is_diploid is True
        assert len(prep.allele_templates) == 2
        assert len(prep.allele_coverages) == 2
        assert sum(prep.allele_coverages) >= 2  # each allele gets at least 1
        for t in prep.allele_templates:
            assert t.exists()
```

- [ ] **Step 2: Create `amplicon_common.py`**

Create `muc_one_up/read_simulator/amplicon_common.py`:

```python
"""Shared amplicon extraction and preparation stages.

Used by both PacBio and ONT amplicon pipelines. Handles:
- Diploid detection and haplotype extraction
- Primer-based amplicon extraction
- PCR bias coverage split
- Template FASTA generation
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .pcr_bias import PCRBiasModel
from .utils.amplicon_extractor import AmpliconExtractor
from .utils.reference_utils import extract_haplotypes, is_diploid_reference
from .utils.template_generator import generate_template_fasta


@dataclass
class AmpliconPrep:
    """Result of amplicon extraction and preparation stages."""

    allele_templates: list[Path]
    allele_coverages: list[int]
    output_dir: Path
    output_base: str
    intermediate_files: list[str] = field(default_factory=list)
    is_diploid: bool = False


def extract_and_prepare_amplicons(
    *,
    input_fa: str,
    forward_primer: str,
    reverse_primer: str,
    total_coverage: int,
    work_dir: Path,
    expected_product_range: tuple[int, int] | None = None,
    pcr_bias_config: dict[str, Any] | None = None,
    seed: int | None = None,
) -> AmpliconPrep:
    """Extract amplicons and prepare template FASTAs for simulation.

    Stages:
    1. Detect diploid, split haplotypes if needed
    2. Extract amplicon via primer binding sites
    3. Compute PCR bias coverage split (diploid only)
    4. Generate template FASTAs (N copies per allele)

    Args:
        input_fa: Path to input FASTA (haploid or diploid).
        forward_primer: Forward primer sequence.
        reverse_primer: Reverse primer sequence.
        total_coverage: Total number of template molecules.
        work_dir: Working directory for intermediate files.
        expected_product_range: Optional (min, max) expected amplicon size.
        pcr_bias_config: PCR bias config dict (preset, stochastic, etc.).
        seed: Random seed for reproducibility.

    Returns:
        AmpliconPrep with template paths, coverage per allele, etc.
    """
    intermediate_files: list[str] = []

    # STAGE 1: Haplotype detection and extraction
    diploid = is_diploid_reference(input_fa)

    if diploid:
        logging.info("STAGE 1: Extracting haplotypes from diploid reference")
        hap1_fa, hap2_fa = extract_haplotypes(
            input_fa, work_dir, base_name="amplicon_sim"
        )
        haplotype_fastas = [hap1_fa, hap2_fa]
    else:
        logging.info("STAGE 1: Haploid reference — skipping extraction")
        haplotype_fastas = [input_fa]

    # STAGE 2: Amplicon extraction per haplotype
    logging.info("STAGE 2: Extracting amplicons via primer binding sites")

    extractor = AmpliconExtractor(
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        expected_product_range=expected_product_range,
    )

    amplicon_results = []
    for i, hap_fa in enumerate(haplotype_fastas, 1):
        amp_out = str(work_dir / f"amplicon_hap{i}.fa")
        result = extractor.extract(hap_fa, amp_out)
        amplicon_results.append(result)
        logging.info("  Haplotype %d amplicon: %d bp", i, result.length)

    # STAGE 3: PCR bias coverage split
    pcr_model = PCRBiasModel.from_config(pcr_bias_config or {})

    if diploid:
        logging.info("STAGE 3: Computing PCR bias coverage split")
        n1, n2 = pcr_model.compute_coverage_split(
            total_coverage,
            amplicon_results[0].length,
            amplicon_results[1].length,
            seed=seed,
        )
        allele_counts = [max(1, n1), max(1, n2)]
        if n1 == 0 or n2 == 0:
            logging.warning(
                "PCR bias produced 0 reads for one allele at coverage=%d. "
                "Each allele will receive at least 1 read.",
                total_coverage,
            )
        logging.info(
            "  Coverage split: allele1=%d, allele2=%d (total=%d)",
            allele_counts[0],
            allele_counts[1],
            total_coverage,
        )
    else:
        logging.info("STAGE 3: Haploid — full coverage to single allele")
        allele_counts = [total_coverage]

    # STAGE 4: Template FASTA generation
    logging.info("STAGE 4: Generating template FASTAs")

    template_fastas: list[Path] = []
    for i, (amp_result, count) in enumerate(
        zip(amplicon_results, allele_counts, strict=False), 1
    ):
        template_out = work_dir / f"template_hap{i}.fa"
        generate_template_fasta(amp_result.fasta_path, count, str(template_out))
        template_fastas.append(template_out)
        logging.info(
            "  Haplotype %d: %d copies of %d bp amplicon",
            i,
            count,
            amp_result.length,
        )

    return AmpliconPrep(
        allele_templates=template_fastas,
        allele_coverages=allele_counts,
        output_dir=work_dir,
        output_base="amplicon",
        intermediate_files=intermediate_files,
        is_diploid=diploid,
    )
```

- [ ] **Step 3: Run tests**

Run: `python -m pytest tests/read_simulator/test_amplicon_common.py -v`
Expected: All tests pass.

- [ ] **Step 4: Refactor amplicon_pipeline.py to use the shared helper**

In `muc_one_up/read_simulator/amplicon_pipeline.py`, replace stages 1-4 (lines 177-255, from `with tempfile.TemporaryDirectory` through template FASTA generation) with a call to `extract_and_prepare_amplicons()`. The refactored code should:

1. Add import: `from .amplicon_common import extract_and_prepare_amplicons`
2. Inside the `try: with tempfile.TemporaryDirectory(...)` block, replace stages 1-4 with:

```python
            from .amplicon_common import extract_and_prepare_amplicons

            prep = extract_and_prepare_amplicons(
                input_fa=input_fa,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                total_coverage=total_coverage,
                work_dir=temp_path,
                expected_product_range=expected_range_tuple,
                pcr_bias_config=pcr_bias_config,
                seed=seed,
            )

            intermediate_files.extend(prep.intermediate_files)
            template_fastas = [str(t) for t in prep.allele_templates]
            allele_counts = prep.allele_coverages
            diploid = prep.is_diploid
```

3. Remove the now-unused imports: `AmpliconExtractor`, `extract_haplotypes`, `is_diploid_reference`, `generate_template_fasta`, `PCRBiasModel`
4. Keep stages 5-8 unchanged

- [ ] **Step 5: Change PacBio amplicon metadata platform tag**

In `amplicon_pipeline.py`, change the `create_pipeline_metadata` call's `platform` from `"PacBio-Amplicon"` to `"PacBio"`. Also set `assay_type` in config before the metadata call:

```python
            config.setdefault("read_simulation", {})["assay_type"] = "amplicon"
            create_pipeline_metadata(
                output_dir=output_dir,
                output_base=output_base,
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="PacBio",
                tools_used=["pbsim3", "ccs", "minimap2", "samtools"],
            )
```

- [ ] **Step 6: Run existing amplicon tests to verify no regression**

Run: `python -m pytest tests/read_simulator/test_amplicon_pipeline.py -v`
Expected: All existing tests pass (behavior unchanged).

- [ ] **Step 7: Run full suite**

Run: `python -m pytest --tb=short -q`

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/read_simulator/amplicon_common.py muc_one_up/read_simulator/amplicon_pipeline.py tests/read_simulator/test_amplicon_common.py
git commit -m "refactor: extract shared amplicon stages into amplicon_common.py

Stages 1-3 (haplotype extraction, primer-based amplicon extraction,
PCR bias, template generation) are now in extract_and_prepare_amplicons().
PacBio amplicon pipeline refactored to use it. No behavior change.
Metadata platform changed from 'PacBio-Amplicon' to 'PacBio' with
assay_type='amplicon' field."
```

---

### Task 4: Create ONT amplicon pipeline

**Files:**
- Create: `muc_one_up/read_simulator/ont_amplicon_pipeline.py`
- Create: `tests/read_simulator/test_ont_amplicon_pipeline.py`

- [ ] **Step 1: Write the test**

Create `tests/read_simulator/test_ont_amplicon_pipeline.py`:

```python
"""Tests for ONT amplicon simulation pipeline."""

from unittest.mock import patch

import pytest

from muc_one_up.exceptions import ReadSimulationError


@pytest.fixture
def muc1_primers():
    return {
        "forward": "GGAGAAAAGGAGACTTCGGCTACCCAG",
        "reverse": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    }


@pytest.fixture
def diploid_fasta_with_primers(tmp_path, muc1_primers):
    from Bio.Seq import Seq

    fwd = muc1_primers["forward"]
    rev_rc = str(Seq(muc1_primers["reverse"]).reverse_complement())
    vntr1 = "ACGT" * 400
    seq1 = "N" * 50 + fwd + vntr1 + rev_rc + "N" * 50
    vntr2 = "ACGT" * 800
    seq2 = "N" * 50 + fwd + vntr2 + rev_rc + "N" * 50
    fasta = tmp_path / "diploid.fa"
    fasta.write_text(f">hap1\n{seq1}\n>hap2\n{seq2}\n")
    return fasta


@pytest.fixture
def ont_amplicon_config(tmp_path, muc1_primers):
    model_file = tmp_path / "ERRHMM-ONT-HQ.model"
    model_file.write_text("model")
    return {
        "tools": {
            "pbsim3": "pbsim",
            "samtools": "samtools",
            "minimap2": "minimap2",
        },
        "amplicon_params": {
            "forward_primer": muc1_primers["forward"],
            "reverse_primer": muc1_primers["reverse"],
            "pcr_bias": {"preset": "default"},
        },
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": str(model_file),
            "threads": 4,
        },
        "read_simulation": {
            "simulator": "ont-amplicon",
            "coverage": 50,
        },
    }


class TestOntAmpliconPipeline:
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.align_reads_with_minimap2")
    def test_calls_pbsim3_with_pass_num_1(
        self,
        mock_align,
        mock_convert,
        mock_pbsim3,
        diploid_fasta_with_primers,
        ont_amplicon_config,
        tmp_path,
    ):
        """ONT amplicon must call pbsim3 with pass_num=1."""
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_convert.return_value = str(tmp_path / "reads.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["reads.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        simulate_ont_amplicon_pipeline(
            ont_amplicon_config,
            str(diploid_fasta_with_primers),
            human_reference=str(diploid_fasta_with_primers),
        )

        # Verify pass_num=1 in all pbsim3 calls
        for call in mock_pbsim3.call_args_list:
            assert call.kwargs["pass_num"] == 1

    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.run_pbsim3_template_simulation")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.ont_amplicon_pipeline.align_reads_with_minimap2")
    def test_uses_map_ont_preset(
        self,
        mock_align,
        mock_convert,
        mock_pbsim3,
        diploid_fasta_with_primers,
        ont_amplicon_config,
        tmp_path,
    ):
        """ONT amplicon must align with map-ont preset."""
        mock_pbsim3.side_effect = lambda **kw: [kw["output_prefix"] + ".bam"]
        mock_convert.return_value = str(tmp_path / "reads.fastq")
        mock_align.return_value = str(tmp_path / "aligned.bam")

        for name in ["reads.fastq", "aligned.bam"]:
            (tmp_path / name).write_bytes(b"FAKE")

        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        simulate_ont_amplicon_pipeline(
            ont_amplicon_config,
            str(diploid_fasta_with_primers),
            human_reference=str(diploid_fasta_with_primers),
        )

        mock_align.assert_called_once()
        assert mock_align.call_args.kwargs["preset"] == "map-ont"

    def test_rejects_source_tracker(self, ont_amplicon_config, tmp_path):
        """Source tracking not supported — must raise."""
        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        fasta = tmp_path / "test.fa"
        fasta.write_text(">seq\nACGT\n")

        with pytest.raises(ReadSimulationError, match="source tracking"):
            simulate_ont_amplicon_pipeline(
                ont_amplicon_config,
                str(fasta),
                source_tracker="not_none",
            )
```

- [ ] **Step 2: Create `ont_amplicon_pipeline.py`**

Create `muc_one_up/read_simulator/ont_amplicon_pipeline.py`:

```python
"""ONT amplicon read simulation pipeline.

Single-pass amplicon simulation using pbsim3 template mode with ONT
error models. Shares extraction/PCR-bias stages with the PacBio
amplicon pipeline via amplicon_common.py.

Pipeline stages:
1-3. Amplicon extraction and preparation (shared)
4.   PBSIM3 template mode, ONT model, pass_num=1
5.   BAM → FASTQ conversion (no CCS — ONT is single-pass)
6.   Merge per-allele FASTQs (diploid)
7.   Align to reference (minimap2 map-ont)
8.   Write metadata, cleanup
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .output_config import OutputConfig

from ..exceptions import ExternalToolError, FileOperationError, ReadSimulationError
from .amplicon_common import extract_and_prepare_amplicons
from .constants import MINIMAP2_PRESET_ONT
from .pipeline_utils import (
    cleanup_intermediates,
    create_pipeline_metadata,
    resolve_pipeline_outputs,
)
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_template_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq


def simulate_ont_amplicon_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None,
    source_tracker: Any | None = None,
    output_config: OutputConfig | None = None,
) -> str:
    """Run ONT amplicon simulation pipeline.

    Uses pbsim3 template mode with ONT error models (single-pass,
    no CCS consensus). See module docstring for pipeline stages.

    Args:
        config: Pipeline config with tools, amplicon_params, pacbio_params,
            read_simulation sections.
        input_fa: Input FASTA (haploid or diploid).
        human_reference: Optional reference for alignment.
        source_tracker: Not supported — raises if non-None.
        output_config: Optional output path control.

    Returns:
        Path to final output (aligned BAM or FASTQ).
    """
    if source_tracker is not None:
        raise ReadSimulationError(
            "Read source tracking is not yet supported for amplicon simulation. "
            "Remove --track-read-source to proceed."
        )

    from typing import cast

    from ..type_defs import AmpliconConfig, PacbioConfig, ReadSimulationConfig

    amplicon_params = cast(AmpliconConfig, config.get("amplicon_params", {}))
    pacbio_params = cast(PacbioConfig, config.get("pacbio_params", {}))
    tools = config.get("tools", {})
    rs_raw = config.get("read_simulation", {})
    rs_config = cast(ReadSimulationConfig, rs_raw)

    # Tool paths (no CCS needed for ONT)
    pbsim3_cmd = tools.get("pbsim3", "pbsim")
    samtools_cmd = tools.get("samtools", "samtools")
    minimap2_cmd = tools.get("minimap2", "minimap2")

    # Params from config
    model_type = pacbio_params["model_type"]
    model_file = pacbio_params["model_file"]
    threads = pacbio_params.get("threads", 4)
    seed = pacbio_params.get("seed")
    accuracy_mean = pacbio_params.get("accuracy_mean", 0.85)

    forward_primer = amplicon_params["forward_primer"]
    reverse_primer = amplicon_params["reverse_primer"]
    expected_range = amplicon_params.get("expected_product_range")
    expected_range_tuple: tuple[int, int] | None = None
    if expected_range is not None:
        if len(expected_range) != 2:
            raise ValueError(
                "amplicon_params.expected_product_range must contain exactly "
                "2 values: [min_product_size, max_product_size]"
            )
        expected_range_tuple = (int(expected_range[0]), int(expected_range[1]))

    _rs_cov = rs_config.get("coverage")
    _pb_cov = pacbio_params.get("coverage")
    total_coverage: int = int(
        _rs_cov if _rs_cov is not None else (_pb_cov if _pb_cov is not None else 30)
    )

    pcr_bias_config = amplicon_params.get("pcr_bias", {})

    output_dir, output_base = resolve_pipeline_outputs(input_fa, rs_raw, output_config)

    start_time = datetime.now()
    intermediate_files: list[str] = []

    try:
        with tempfile.TemporaryDirectory(prefix="ont_amplicon_sim_") as temp_dir:
            temp_path = Path(temp_dir)

            logging.info("=" * 80)
            logging.info("ONT AMPLICON SIMULATION PIPELINE")
            logging.info("=" * 80)

            # STAGES 1-3: Shared amplicon extraction and preparation
            prep = extract_and_prepare_amplicons(
                input_fa=input_fa,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                total_coverage=total_coverage,
                work_dir=temp_path,
                expected_product_range=expected_range_tuple,
                pcr_bias_config=pcr_bias_config,
                seed=seed,
            )
            intermediate_files.extend(prep.intermediate_files)

            # STAGE 4: PBSIM3 template mode — ONT, single-pass
            logging.info("STAGE 4: Running PBSIM3 template mode (ONT, pass_num=1)")

            allele_bams: list[list[str]] = []
            for i, template_fa in enumerate(prep.allele_templates, 1):
                prefix = str(temp_path / f"ont_hap{i}")
                hap_seed = (seed + i) if seed is not None else None

                bams = run_pbsim3_template_simulation(
                    pbsim3_cmd=pbsim3_cmd,
                    samtools_cmd=samtools_cmd,
                    template_fasta=str(template_fa),
                    model_type=model_type,
                    model_file=model_file,
                    output_prefix=prefix,
                    pass_num=1,
                    accuracy_mean=accuracy_mean,
                    seed=hap_seed,
                )
                allele_bams.append(bams)
                intermediate_files.extend(bams)
                logging.info("  Haplotype %d: %d BAMs", i, len(bams))

            # STAGE 5: BAM → FASTQ conversion (no CCS for ONT)
            logging.info("STAGE 5: Converting BAMs to FASTQ (no CCS)")

            allele_fastqs: list[str] = []
            for i, bam_group in enumerate(allele_bams, 1):
                for j, bam in enumerate(bam_group, 1):
                    fq_out = str(temp_path / f"ont_hap{i}_{j:04d}.fastq")
                    convert_bam_to_fastq(
                        samtools_cmd=samtools_cmd,
                        input_bam=bam,
                        output_fastq=fq_out,
                        threads=threads,
                    )
                    allele_fastqs.append(fq_out)

            # STAGE 6: Merge FASTQs
            merged_fastq = str(output_dir / f"{output_base}_amplicon_ont.fastq")
            logging.info("STAGE 6: Merging %d FASTQs", len(allele_fastqs))

            with open(merged_fastq, "w") as outf:
                for fq in allele_fastqs:
                    with open(fq) as inf:
                        shutil.copyfileobj(inf, outf)

            # STAGE 7: Alignment (optional)
            if human_reference is None:
                logging.info("ONT amplicon simulation complete (no alignment)")
                logging.info("Final output: %s", merged_fastq)
                final_output = merged_fastq
            else:
                logging.info("STAGE 7: Aligning with minimap2 (map-ont)")
                aligned_bam = str(output_dir / f"{output_base}_amplicon_ont.bam")
                aligned_bam = align_reads_with_minimap2(
                    minimap2_cmd=minimap2_cmd,
                    samtools_cmd=samtools_cmd,
                    reference=human_reference,
                    reads_fastq=merged_fastq,
                    output_bam=aligned_bam,
                    preset=MINIMAP2_PRESET_ONT,
                    threads=threads,
                )
                final_output = aligned_bam

            # STAGE 8: Metadata
            end_time = datetime.now()
            duration = end_time - start_time
            logging.info("=" * 80)
            logging.info(
                "ONT amplicon simulation complete (duration: %s)",
                str(duration).split(".")[0],
            )
            logging.info("Final output: %s", final_output)
            logging.info("=" * 80)

            config.setdefault("read_simulation", {})["assay_type"] = "amplicon"
            create_pipeline_metadata(
                output_dir=output_dir,
                output_base=f"{output_base}_amplicon_ont",
                config=config,
                start_time=start_time,
                end_time=end_time,
                platform="ONT",
                tools_used=["pbsim3", "minimap2", "samtools"],
            )

            return final_output

    except (ExternalToolError, FileOperationError) as e:
        logging.error("ONT amplicon pipeline failed: %s", e)
        raise RuntimeError(
            f"ONT amplicon simulation pipeline failed.\nError: {e}\n\n"
            f"Troubleshooting:\n"
            f"  1. Verify tools installed: pbsim3, samtools, minimap2\n"
            f"  2. Check model file: {model_file}\n"
            f"  3. Ensure primers bind in the reference"
        ) from e
    except Exception as e:
        logging.error("Unexpected error in ONT amplicon pipeline: %s", e)
        raise RuntimeError(f"ONT amplicon pipeline failed: {e}") from e
    finally:
        try:
            if intermediate_files:
                cleanup_intermediates(intermediate_files)
        except Exception as cleanup_err:
            logging.warning("Cleanup failed (non-fatal): %s", cleanup_err)
```

- [ ] **Step 3: Run tests**

Run: `python -m pytest tests/read_simulator/test_ont_amplicon_pipeline.py -v`
Expected: All 3 tests pass.

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/read_simulator/ont_amplicon_pipeline.py tests/read_simulator/test_ont_amplicon_pipeline.py
git commit -m "feat: ONT amplicon pipeline — single-pass pbsim3 + map-ont

Uses shared amplicon extraction stages, pbsim3 template mode with
pass_num=1 and ONT error models, no CCS consensus, minimap2 map-ont
preset. Output: *_amplicon_ont.bam / .fastq"
```

---

### Task 5: CLI --platform flag and dispatch integration

**Files:**
- Modify: `muc_one_up/cli/commands/reads.py:439-555`
- Modify: `muc_one_up/read_simulation.py:89-99`
- Modify: `muc_one_up/config.py:285`
- Modify: `tests/test_lazy_imports.py:65`
- Create: `tests/cli/test_amplicon_platform.py`

- [ ] **Step 1: Add `"ont-amplicon"` to config schema**

In `muc_one_up/config.py:285`, change:

```python
"simulator": {"type": "string", "enum": ["illumina", "ont", "pacbio", "amplicon"]},
```

to:

```python
"simulator": {"type": "string", "enum": ["illumina", "ont", "pacbio", "amplicon", "ont-amplicon"]},
```

- [ ] **Step 2: Add `"ont-amplicon"` to dispatch**

In `muc_one_up/read_simulation.py`, add a new branch before the `else` clause (before line 97):

```python
    elif simulator_type == "ont-amplicon":
        from muc_one_up.read_simulator.ont_amplicon_pipeline import (
            simulate_ont_amplicon_pipeline,
        )

        return lambda config, input_fa, human_reference, **kw: simulate_ont_amplicon_pipeline(
            config, input_fa, human_reference=human_reference, **kw
        )
```

Update the `else` branch's valid list:

```python
    else:
        valid = "amplicon, illumina, ont, ont-amplicon, pacbio"
        raise ValueError(f"Unknown simulator: '{simulator_type}'. Valid options: {valid}. ")
```

- [ ] **Step 3: Add `--platform` to CLI amplicon command**

In `muc_one_up/cli/commands/reads.py`, add the `--platform` option to the `amplicon` command. Add it after the existing `--stochastic-pcr` option (around line 463):

```python
@click.option(
    "--platform",
    type=click.Choice(["pacbio", "ont"]),
    default="pacbio",
    show_default=True,
    help="Sequencing platform for amplicon simulation.",
)
```

Add `platform` to the function signature:

```python
def amplicon(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    model_type,
    model_file,
    pcr_preset,
    stochastic_pcr,
    platform,
    seed,
    track_read_source,
):
```

Update the docstring to mention ONT support.

In the function body, after `_setup_read_config`, set the simulator type based on platform:

```python
    if platform == "ont":
        config["read_simulation"]["simulator"] = "ont-amplicon"
    config["read_simulation"]["assay_type"] = "amplicon"
```

- [ ] **Step 4: Update dispatch tests**

In `tests/test_lazy_imports.py:65`, change:

```python
    @pytest.mark.parametrize("sim_type", ["illumina", "ont", "pacbio", "amplicon"])
```

to:

```python
    @pytest.mark.parametrize("sim_type", ["illumina", "ont", "pacbio", "amplicon", "ont-amplicon"])
```

- [ ] **Step 5: Write CLI tests**

Create `tests/cli/test_amplicon_platform.py`:

```python
"""Tests for reads amplicon --platform routing."""

import json
from unittest.mock import patch

from click.testing import CliRunner

from muc_one_up.cli.click_main import cli


def _invoke_amplicon(tmp_path, platform=None, extra_args=None):
    """Invoke reads amplicon and return the config passed to the simulator."""
    model_file = tmp_path / "test.model"
    model_file.write_text("model")

    config_dict = {
        "amplicon_params": {
            "forward_primer": "ACGT",
            "reverse_primer": "TGCA",
        },
        "pacbio_params": {
            "model_type": "errhmm",
            "model_file": str(model_file),
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config_dict))

    fasta = tmp_path / "input.fa"
    fasta.write_text(">seq\nACGT\n")

    captured = {}

    def fake_load(path):
        return json.loads(config_path.read_text())

    def fake_batch(config, *args, **kwargs):
        captured["config"] = config

    runner = CliRunner()
    with (
        patch("muc_one_up.config.load_config_raw", side_effect=fake_load),
        patch(
            "muc_one_up.cli.commands.reads._run_batch_simulation",
            side_effect=fake_batch,
        ),
    ):
        args = ["--config", str(config_path), "reads", "amplicon"]
        if platform:
            args.extend(["--platform", platform])
        if extra_args:
            args.extend(extra_args)
        args.append(str(fasta))
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    return captured["config"]


def test_default_platform_is_pacbio(tmp_path):
    """Without --platform, simulator should be 'amplicon' (PacBio)."""
    config = _invoke_amplicon(tmp_path)
    assert config["read_simulation"]["simulator"] == "amplicon"


def test_platform_ont_sets_ont_amplicon(tmp_path):
    """--platform ont should set simulator to 'ont-amplicon'."""
    config = _invoke_amplicon(tmp_path, platform="ont")
    assert config["read_simulation"]["simulator"] == "ont-amplicon"


def test_assay_type_set_for_both_platforms(tmp_path):
    """Both platforms should set assay_type='amplicon'."""
    for platform in [None, "ont"]:
        config = _invoke_amplicon(tmp_path, platform=platform)
        assert config["read_simulation"]["assay_type"] == "amplicon"
```

- [ ] **Step 6: Run all tests**

Run: `python -m pytest tests/cli/test_amplicon_platform.py tests/test_lazy_imports.py -v`
Expected: All pass.

- [ ] **Step 7: Run full quality gate**

```bash
ruff check muc_one_up/ tests/
ruff format --check muc_one_up/ tests/
python -m mypy muc_one_up/
python -m pytest --tb=short -q
```

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/cli/commands/reads.py muc_one_up/read_simulation.py muc_one_up/config.py tests/test_lazy_imports.py tests/cli/test_amplicon_platform.py
git commit -m "feat: add --platform ont to reads amplicon command

Routes to ONT amplicon pipeline via simulator='ont-amplicon'.
Default --platform pacbio preserves backward compatibility.
Config schema and dispatch updated. assay_type='amplicon' set
for metadata."
```

---

## Final Verification

- [ ] **Run full quality gate**

```bash
ruff check muc_one_up/ tests/
ruff format --check muc_one_up/ tests/
python -m mypy muc_one_up/
python -m pytest --tb=short -q
```

All four commands must pass cleanly.
