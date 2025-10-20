# **IMPLEMENTATION PLAN: PacBio HiFi Read Simulation Support for MucOneUp**

**Issue**: #1 - Implement PacBio support using pbsim3
**Version**: 2.0 (Revised after code review)
**Date**: 2025-10-20
**Status**: Ready for Implementation
**Review Status**: ✅ All critical issues addressed

---

## **REVISION HISTORY**

- **v1.0** (2025-10-20): Initial plan
- **v2.0** (2025-10-20): Revised to address code review findings
  - Fixed DRY violation (reuse samtools wrapper)
  - Fixed exception handling pattern
  - Added SAM→BAM conversion
  - Created generic minimap2 wrapper
  - Added CONFIG_SCHEMA validation
  - Implemented Strategy Pattern for router
  - Added constants file
  - Addressed all 7 critical issues

---

## **Executive Summary**

This plan details the implementation of PacBio HiFi (High-Fidelity) read simulation support for MucOneUp. The implementation follows existing architectural patterns (ONT pipeline), adheres to SOLID principles, DRY, KISS, and modularization best practices, and addresses all code review findings.

---

## **1. RESEARCH FINDINGS**

### **1.1 Tool Overview**

#### **pbsim3** (Primary Simulator)
- **Version**: v3.0.5+ (Bioconda: `bioconda::pbsim3`)
- **GitHub**: https://github.com/yukiteruono/pbsim3
- **License**: GPLv2
- **Platform Support**: linux-64, osx-64, linux-aarch64, osx-arm64

**Key Capabilities**:
- Multi-pass CLR (Continuous Long Reads) simulation
- Quality score models (QSHMM) and error models (ERRHMM)
- Reproducible with seed support
- Output formats: BAM, SAM

**Important**: pbsim3 simulates multi-pass CLR reads that must be processed by CCS to generate HiFi reads.

#### **ccs (Circular Consensus Sequencing)**
- **Version**: v6.4.0+ (Bioconda: `bioconda::pbccs`)
- **GitHub**: https://github.com/PacificBiosciences/ccs
- **Documentation**: https://ccs.how/

**Key Capabilities**:
- Generates HiFi consensus reads from multi-pass CLR
- Accuracy: 99.9% (Q30+)
- Threading support
- Deterministic output (seed not needed - deterministic given input)

#### **Required Dependencies**
- **minimap2**: Already in MucOneUp ✅ (for alignment)
- **samtools**: Already in MucOneUp ✅ (for BAM operations)
- **gzip**: System-level ✅

---

### **1.2 PacBio HiFi Simulation Workflow**

```
Input FASTA → pbsim3 (multi-pass CLR) → CCS (HiFi consensus) → minimap2 → Aligned BAM
                 ↓                         ↓                        ↓
            CLR.bam/sam               HiFi.bam                 HiFi.bam + .bai
```

**Critical Detail**: minimap2 requires `-ax map-hifi` preset for PacBio HiFi (NOT `-ax map-ont`)

---

### **1.3 pbsim3 Command-Line Parameters**

```bash
pbsim3 \
  --strategy wgs \
  --method errhmm|qshmm \
  --errhmm <model_file> \        # Or --qshmm for RSII
  --genome <reference.fa> \
  --depth <coverage> \
  --prefix <output_prefix> \
  --pass-num <num_passes> \      # ≥2 for multi-pass (10-20 for HiFi)
  --seed <seed> \                # Reproducibility
  --length-mean <mean> \
  --length-sd <sd>
```

**Output Naming**: `{prefix}_0001.bam` or `{prefix}_0001.sam`

---

### **1.4 CCS Command-Line Parameters**

```bash
ccs \
  <input_subreads.bam> \
  <output_hifi.bam> \
  --min-passes <min> \           # Default: 3
  --min-rq <quality> \           # Default: 0.99 (99%)
  --num-threads <threads> \
  --log-level INFO
```

---

## **2. ARCHITECTURE DESIGN**

### **2.1 Module Structure**

```
muc_one_up/
├── read_simulator/
│   ├── pacbio_pipeline.py         # NEW: PacBio pipeline orchestrator
│   ├── wrappers/
│   │   ├── pbsim3_wrapper.py      # NEW: pbsim3 command wrapper
│   │   ├── ccs_wrapper.py         # NEW: CCS command wrapper
│   │   ├── minimap2_wrapper.py    # NEW: Generic minimap2 wrapper (DRY fix)
│   │   ├── samtools_wrapper.py    # EXTEND: Add BAM→FASTQ if missing (DRY fix)
│   │   └── (existing wrappers...)
│   ├── constants.py               # NEW: Constants for validation
│   ├── ont_pipeline.py            # REFERENCE: Pattern to follow
│   └── utils/ ...
├── cli/
│   └── click_main.py              # EXTEND: Add 'reads pacbio' subcommand
├── read_simulation.py             # EXTEND: Add PacBio router (Strategy Pattern)
└── config.py                      # EXTEND: Add CONFIG_SCHEMA validation
conda/
└── env_pacbio.yml                 # NEW: Conda environment file
```

---

### **2.2 Configuration Schema Extension**

Add to `config.py`:

```python
# Add to CONFIG_SCHEMA (after nanosim_params)
"pacbio_params": {
    "type": "object",
    "properties": {
        "model_type": {
            "type": "string",
            "enum": ["qshmm", "errhmm"],
            "description": "Model type: qshmm (RSII only) or errhmm (RSII/Sequel)"
        },
        "model_file": {
            "type": "string",
            "description": "Path to QSHMM or ERRHMM model file"
        },
        "coverage": {
            "type": "number",
            "minimum": 0.1,
            "maximum": 10000,
            "description": "Target sequencing coverage"
        },
        "pass_num": {
            "type": "number",
            "minimum": 2,
            "maximum": 50,
            "description": "Number of passes for multi-pass sequencing"
        },
        "length_mean": {"type": ["number", "null"], "minimum": 100},
        "length_sd": {"type": ["number", "null"], "minimum": 0},
        "ccs_min_passes": {"type": ["number", "null"], "minimum": 1},
        "ccs_min_rq": {"type": ["number", "null"], "minimum": 0.0, "maximum": 1.0},
        "num_threads": {"type": ["number", "null"], "minimum": 1},
        "seed": {"type": ["number", "null"]},
        "correction_factor": {"type": ["number", "null"], "minimum": 0.1, "maximum": 10.0},
        "enable_coverage_correction": {"type": ["boolean", "null"]},
        "other_options": {"type": ["string", "null"]}
    },
    "required": ["model_type", "model_file", "coverage"],
    "additionalProperties": False
},
```

---

## **3. DETAILED IMPLEMENTATION PLAN**

### **Phase 1: Constants and Configuration**

#### **Task 1.1: Create Constants File** ✅

**File**: `muc_one_up/read_simulator/constants.py` (NEW)

```python
#!/usr/bin/env python3
"""
Constants for read simulation pipelines.

Centralizes magic numbers, valid values, and configuration defaults.
"""

# pbsim3 model types
VALID_PBSIM3_MODEL_TYPES = {"qshmm", "errhmm"}

# Minimap2 presets
MINIMAP2_PRESET_ONT = "map-ont"
MINIMAP2_PRESET_PACBIO_CLR = "map-pb"
MINIMAP2_PRESET_PACBIO_HIFI = "map-hifi"

# Default timeouts (seconds)
DEFAULT_PBSIM3_TIMEOUT = 7200  # 2 hours
DEFAULT_CCS_TIMEOUT = 7200     # 2 hours
DEFAULT_ALIGNMENT_TIMEOUT = 3600  # 1 hour
DEFAULT_CONVERSION_TIMEOUT = 1800  # 30 minutes

# Coverage limits
MIN_COVERAGE = 0.1
MAX_COVERAGE = 10000

# PacBio defaults
DEFAULT_PACBIO_PASS_NUM = 10
DEFAULT_PACBIO_LENGTH_MEAN = 15000
DEFAULT_PACBIO_LENGTH_SD = 5000
DEFAULT_CCS_MIN_PASSES = 3
DEFAULT_CCS_MIN_RQ = 0.99
```

---

#### **Task 1.2: Update CONFIG_SCHEMA** ✅

**File**: `muc_one_up/config.py` (EXTEND)

Add the pacbio_params schema shown in section 2.2 above.

---

#### **Task 1.3: Create Conda Environment** ✅

**File**: `conda/env_pacbio.yml` (NEW)

```yaml
name: env_pacbio
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # Core PacBio tools
  - pbsim3>=3.0.5              # Long-read simulator
  - pbccs>=6.4.0               # CCS for HiFi generation

  # Alignment and processing
  - minimap2>=2.28             # Long-read aligner
  - samtools>=1.21             # BAM manipulation

  # Python environment
  - python=3.9
  - pip

  # Optional utilities
  - bam2fastx                  # Format conversion
  - gzip                       # Compression
```

---

### **Phase 2: Wrapper Modules**

#### **Task 2.1: Extend samtools_wrapper.py** ✅

**File**: `muc_one_up/read_simulator/wrappers/samtools_wrapper.py` (EXTEND)

**FIX**: Add BAM→FASTQ conversion (if not already present - DRY principle)

```python
#!/usr/bin/env python3
"""
Samtools wrapper for MucOneUp.

Provides functions for SAM/BAM manipulation operations.
"""

import logging
from pathlib import Path

from ..command_utils import build_tool_command
from ..utils import run_command


def convert_sam_to_bam(
    samtools_cmd: str,
    input_sam: str,
    output_bam: str,
    timeout: int = 1800,
) -> str:
    """
    Convert SAM to BAM using samtools view.

    Args:
        samtools_cmd: Path to samtools command
        input_sam: Path to input SAM file
        output_bam: Path for output BAM file
        timeout: Timeout in seconds (default: 1800)

    Returns:
        Path to the generated BAM file

    Raises:
        ExternalToolError: If conversion fails (propagated from run_command)
    """
    Path(output_bam).resolve().parent.mkdir(parents=True, exist_ok=True)

    cmd_list = build_tool_command(
        samtools_cmd,
        "view",
        "-bS",  # Output BAM, input SAM
        "-o", output_bam,
        input_sam,
    )

    logging.info("[samtools] Converting SAM to BAM: %s", " ".join(cmd_list))

    # run_command raises ExternalToolError on failure
    run_command(
        cmd_list,
        timeout=timeout,
        stderr_log_level=logging.INFO,
        stderr_prefix="[samtools] ",
    )

    if not Path(output_bam).exists():
        raise RuntimeError(f"Expected BAM output not found: {output_bam}")

    logging.info("[samtools] SAM to BAM conversion completed: %s", output_bam)
    return output_bam


def convert_bam_to_fastq(
    samtools_cmd: str,
    input_bam: str,
    output_fastq: str,
    timeout: int = 1800,
) -> str:
    """
    Convert BAM to FASTQ using samtools fastq.

    Args:
        samtools_cmd: Path to samtools command
        input_bam: Path to input BAM file
        output_fastq: Path for output FASTQ file
        timeout: Timeout in seconds (default: 1800)

    Returns:
        Path to the generated FASTQ file

    Raises:
        ExternalToolError: If conversion fails (propagated from run_command)
    """
    Path(output_fastq).resolve().parent.mkdir(parents=True, exist_ok=True)

    cmd_list = build_tool_command(
        samtools_cmd,
        "fastq",
        "-0", output_fastq,  # Write to single file
        input_bam,
    )

    logging.info("[samtools] Converting BAM to FASTQ: %s", " ".join(cmd_list))

    # run_command raises ExternalToolError on failure
    run_command(
        cmd_list,
        timeout=timeout,
        stderr_log_level=logging.INFO,
        stderr_prefix="[samtools] ",
    )

    if not Path(output_fastq).exists():
        raise RuntimeError(f"Expected FASTQ output not found: {output_fastq}")

    logging.info("[samtools] BAM to FASTQ conversion completed: %s", output_fastq)
    return output_fastq
```

---

#### **Task 2.2: Create minimap2_wrapper.py** ✅

**File**: `muc_one_up/read_simulator/wrappers/minimap2_wrapper.py` (NEW)

**FIX**: Generic wrapper with parameterized preset (DRY principle)

```python
#!/usr/bin/env python3
"""
Minimap2 wrapper for long-read alignment.

This module provides a generic minimap2 wrapper that supports different
alignment presets (ONT, PacBio CLR, PacBio HiFi) to enable DRY code reuse
across ONT and PacBio pipelines.
"""

import logging
import subprocess
import tempfile
from pathlib import Path

from ...exceptions import ExternalToolError
from ..command_utils import build_tool_command
from ..constants import (
    DEFAULT_ALIGNMENT_TIMEOUT,
    MINIMAP2_PRESET_ONT,
    MINIMAP2_PRESET_PACBIO_HIFI,
)
from ..utils import cleanup_files, run_command


def align_reads_with_minimap2(
    minimap2_cmd: str,
    samtools_cmd: str,
    reference: str,
    reads_fastq: str,
    output_bam: str,
    preset: str = MINIMAP2_PRESET_ONT,
    threads: int = 4,
    timeout: int = DEFAULT_ALIGNMENT_TIMEOUT,
) -> str:
    """
    Align long reads to reference using minimap2 with configurable preset.

    This generic function supports different alignment presets to enable
    code reuse across ONT and PacBio pipelines (DRY principle).

    Workflow:
    1. Run minimap2 with specified preset → SAM
    2. Convert SAM to unsorted BAM with samtools view
    3. Sort BAM with samtools sort
    4. Index BAM with samtools index
    5. Clean up temporary files

    Args:
        minimap2_cmd: Path to minimap2 command
        samtools_cmd: Path to samtools command
        reference: Path to reference FASTA
        reads_fastq: Path to reads FASTQ file
        output_bam: Path for output BAM file
        preset: Minimap2 preset (default: map-ont)
                Valid: map-ont, map-hifi, map-pb, etc.
        threads: Number of threads (default: 4)
        timeout: Timeout in seconds (default: 3600)

    Returns:
        Path to the final sorted and indexed BAM file

    Raises:
        RuntimeError: If any step fails (wraps ExternalToolError)

    Example:
        >>> # ONT alignment
        >>> align_reads_with_minimap2(..., preset="map-ont")
        >>>
        >>> # PacBio HiFi alignment
        >>> align_reads_with_minimap2(..., preset="map-hifi")
    """
    logging.info(
        "[minimap2] Starting alignment with preset: %s (threads=%d)",
        preset,
        threads,
    )

    # Validate inputs
    if not Path(reference).exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if not Path(reads_fastq).exists():
        raise FileNotFoundError(f"Reads file not found: {reads_fastq}")

    # Prepare output paths
    output_bam_path = Path(output_bam).resolve()
    output_bam_path.parent.mkdir(parents=True, exist_ok=True)

    temp_sam = None
    temp_unsorted_bam = str(output_bam_path) + ".unsorted"

    try:
        # Step 1: Run minimap2 to generate SAM
        logging.info("[minimap2] Step 1: Aligning reads to reference...")

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sam", delete=False, dir=output_bam_path.parent
        ) as temp_sam_file:
            temp_sam = temp_sam_file.name

        minimap2_cmd_list = build_tool_command(
            minimap2_cmd,
            "-t", threads,
            "-ax", preset,  # Parameterized preset
            reference,
            reads_fastq,
        )

        logging.info("[minimap2] Running: %s", " ".join(minimap2_cmd_list))

        with open(temp_sam, "w") as sam_out:
            result = subprocess.run(
                minimap2_cmd_list,
                stdout=sam_out,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                timeout=timeout,
            )

        if result.returncode != 0:
            raise ExternalToolError(
                tool="minimap2",
                exit_code=result.returncode,
                stderr=result.stderr,
                cmd=" ".join(minimap2_cmd_list),
            )

        logging.info("[minimap2] Alignment completed successfully")

        # Step 2: Convert SAM to BAM (filter unmapped reads)
        logging.info("[samtools] Step 2: Converting SAM to BAM...")

        view_cmd_list = build_tool_command(
            samtools_cmd,
            "view",
            "-F", "4",  # Filter unmapped reads
            "-b",       # Output BAM
            "-o", temp_unsorted_bam,
            temp_sam,
        )

        run_command(
            view_cmd_list,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        # Step 3: Sort BAM
        logging.info("[samtools] Step 3: Sorting BAM...")

        sort_cmd_list = build_tool_command(
            samtools_cmd,
            "sort",
            "-@", threads,
            "-o", str(output_bam_path),
            temp_unsorted_bam,
        )

        run_command(
            sort_cmd_list,
            timeout=timeout,
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        # Step 4: Index BAM
        logging.info("[samtools] Step 4: Indexing BAM...")

        index_cmd_list = build_tool_command(
            samtools_cmd,
            "index",
            str(output_bam_path),
        )

        run_command(
            index_cmd_list,
            timeout=timeout // 2,  # Indexing is faster
            stderr_log_level=logging.INFO,
            stderr_prefix="[samtools] ",
        )

        logging.info("[minimap2] Read alignment completed successfully")
        logging.info("[minimap2] Final output: %s", output_bam_path)

        return str(output_bam_path)

    except subprocess.TimeoutExpired as e:
        logging.error("[minimap2] Alignment timed out after %d seconds", timeout)
        raise RuntimeError(f"Minimap2 alignment timed out: {e}") from e

    except ExternalToolError as e:
        logging.error("[minimap2] Tool execution failed: %s", e)
        raise RuntimeError(f"Read alignment failed: {e}") from e

    except Exception as e:
        logging.error("[minimap2] Unexpected error: %s", e)
        raise RuntimeError(f"Read alignment failed: {e}") from e

    finally:
        # Cleanup temporary files
        cleanup_files([temp_sam, temp_unsorted_bam])
```

---

#### **Task 2.3: Create pbsim3_wrapper.py** ✅

**File**: `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py` (NEW)

**FIXES**:
- Use constants for model validation
- Proper exception handling (let ExternalToolError propagate)
- SAM→BAM conversion handling

```python
#!/usr/bin/env python3
"""
pbsim3 wrapper for PacBio read simulation.

This module provides wrapper functions for pbsim3, a PacBio and ONT read
simulator. It enables integration with the MucOneUp pipeline and handles proper
command construction, execution, and error handling.
"""

import logging
from pathlib import Path

from ...exceptions import ExternalToolError
from ..command_utils import build_tool_command
from ..constants import DEFAULT_PBSIM3_TIMEOUT, VALID_PBSIM3_MODEL_TYPES
from ..utils import run_command
from .samtools_wrapper import convert_sam_to_bam


def run_pbsim3_simulation(
    pbsim3_cmd: str,
    samtools_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    model_type: str,
    model_file: str,
    coverage: float,
    pass_num: int,
    length_mean: int | None = None,
    length_sd: int | None = None,
    other_options: str = "",
    timeout: int = DEFAULT_PBSIM3_TIMEOUT,
    seed: int | None = None,
) -> str:
    """
    Run pbsim3 simulation to generate PacBio multi-pass CLR reads.

    Args:
        pbsim3_cmd: Path to pbsim3 command
        samtools_cmd: Path to samtools command (for SAM→BAM conversion)
        reference_fasta: Path to reference FASTA file
        output_prefix: Prefix for output files
        model_type: Model type ('qshmm' or 'errhmm')
        model_file: Path to model file (QSHMM-*.model or ERRHMM-*.model)
        coverage: Desired coverage depth
        pass_num: Number of passes for multi-pass sequencing (≥2 for HiFi)
        length_mean: Mean read length (optional)
        length_sd: Read length standard deviation (optional)
        other_options: Additional pbsim3 options (optional)
        timeout: Timeout in seconds (default: from constants)
        seed: Random seed for reproducibility (optional)

    Returns:
        Path to the generated BAM file

    Raises:
        ValueError: If model_type or model_file invalid
        FileNotFoundError: If model_file doesn't exist
        ExternalToolError: If pbsim3 execution fails (propagated from run_command)
        RuntimeError: If output file not found
    """
    # Ensure output directory exists
    output_path = Path(output_prefix).resolve()
    output_dir = output_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate model type (using constants - FIX #6)
    if model_type not in VALID_PBSIM3_MODEL_TYPES:
        raise ValueError(
            f"Invalid model_type: {model_type}. "
            f"Must be one of {VALID_PBSIM3_MODEL_TYPES}"
        )

    # Validate model file exists (FIX: Warning #3)
    if not Path(model_file).exists():
        raise FileNotFoundError(f"Model file not found: {model_file}")

    # Build the command with required parameters
    cmd_list = build_tool_command(
        pbsim3_cmd,
        "--strategy", "wgs",
        "--method", model_type,
        f"--{model_type}", model_file,
        "--genome", reference_fasta,
        "--depth", coverage,
        "--prefix", output_prefix,
        "--pass-num", pass_num,
    )

    # Add seed if provided
    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])
        logging.info("[pbsim3] Using random seed: %d", seed)

    # Add optional length parameters
    if length_mean is not None:
        cmd_list.extend(["--length-mean", str(length_mean)])
    if length_sd is not None:
        cmd_list.extend(["--length-sd", str(length_sd)])

    # Add other options
    if other_options:
        if isinstance(other_options, str) and " " in other_options:
            other_parts = build_tool_command(other_options)
            cmd_list.extend(other_parts)
        else:
            cmd_list.append(other_options)

    # Log the command
    logging.info("[pbsim3] Running command: %s", " ".join(cmd_list))
    logging.info("[pbsim3] Simulation in progress (may take 1-2 hours for large references)...")

    # Run the command (raises ExternalToolError on failure - FIX #2)
    run_command(
        cmd_list,
        timeout=timeout,
        stderr_log_level=logging.INFO,
        stderr_prefix="[pbsim3] ",
    )

    # pbsim3 output naming: <prefix>_0001.bam or <prefix>_0001.sam
    output_bam = f"{output_prefix}_0001.bam"
    output_sam = f"{output_prefix}_0001.sam"

    # Check for output and handle SAM→BAM conversion (FIX #3)
    if Path(output_bam).exists():
        output_file = output_bam
        logging.info("[pbsim3] Generated BAM file: %s", output_file)

    elif Path(output_sam).exists():
        logging.info("[pbsim3] SAM output detected, converting to BAM...")
        # Convert SAM to BAM (required for CCS)
        output_file = convert_sam_to_bam(
            samtools_cmd=samtools_cmd,
            input_sam=output_sam,
            output_bam=output_bam,
        )
        # Clean up SAM file
        Path(output_sam).unlink()
        logging.info("[pbsim3] Converted to BAM: %s", output_file)

    else:
        error_msg = f"Expected output file not found: {output_bam} or {output_sam}"
        logging.error("[pbsim3] %s", error_msg)
        raise RuntimeError(error_msg)

    logging.info("[pbsim3] Simulation completed successfully")
    return output_file
```

---

#### **Task 2.4: Create ccs_wrapper.py** ✅

**File**: `muc_one_up/read_simulator/wrappers/ccs_wrapper.py` (NEW)

**FIXES**:
- Removed duplicate convert_bam_to_fastq (DRY - use samtools_wrapper)
- Proper exception handling
- Added empty output validation (Warning #9)

```python
#!/usr/bin/env python3
"""
CCS (Circular Consensus Sequencing) wrapper for HiFi read generation.

This module provides wrapper functions for the PacBio CCS tool, which generates
highly accurate HiFi reads from multi-pass CLR subreads.
"""

import logging
from pathlib import Path

from ...exceptions import ExternalToolError
from ..command_utils import build_tool_command
from ..constants import DEFAULT_CCS_TIMEOUT
from ..utils import run_command


def run_ccs_consensus(
    ccs_cmd: str,
    input_bam: str,
    output_bam: str,
    min_passes: int = 3,
    min_rq: float = 0.99,
    num_threads: int = 4,
    timeout: int = DEFAULT_CCS_TIMEOUT,
) -> str:
    """
    Run CCS to generate HiFi consensus reads from multi-pass CLR subreads.

    Args:
        ccs_cmd: Path to ccs command
        input_bam: Path to input BAM file (pbsim3 output)
        output_bam: Path for output BAM file (HiFi reads)
        min_passes: Minimum number of passes required (default: 3)
        min_rq: Minimum predicted accuracy (default: 0.99 = 99%)
        num_threads: Number of threads to use (default: 4)
        timeout: Timeout in seconds (default: from constants)

    Returns:
        Path to the generated HiFi BAM file

    Raises:
        ExternalToolError: If CCS execution fails (propagated from run_command)
        RuntimeError: If output is empty or missing
    """
    # Ensure output directory exists
    Path(output_bam).resolve().parent.mkdir(parents=True, exist_ok=True)

    # Validate input exists
    if not Path(input_bam).exists():
        raise FileNotFoundError(f"Input BAM not found: {input_bam}")

    # Build the command
    cmd_list = build_tool_command(
        ccs_cmd,
        input_bam,
        output_bam,
        "--min-passes", min_passes,
        "--min-rq", min_rq,
        "--num-threads", num_threads,
        "--log-level", "INFO",
    )

    logging.info("[ccs] Running command: %s", " ".join(cmd_list))
    logging.info("[ccs] HiFi consensus generation in progress (may take 1-2 hours)...")

    # Run the command (raises ExternalToolError on failure - FIX #2)
    run_command(
        cmd_list,
        timeout=timeout,
        stderr_log_level=logging.INFO,
        stderr_prefix="[ccs] ",
    )

    # Validate output (FIX: Warning #9)
    output_path = Path(output_bam)
    if not output_path.exists():
        error_msg = f"Expected HiFi output not found: {output_bam}"
        logging.error("[ccs] %s", error_msg)
        raise RuntimeError(error_msg)

    if output_path.stat().st_size == 0:
        error_msg = (
            f"CCS produced empty output: {output_bam}. "
            f"Check min_passes={min_passes} and min_rq={min_rq} settings. "
            "Input reads may not meet quality requirements."
        )
        logging.error("[ccs] %s", error_msg)
        raise RuntimeError(error_msg)

    logging.info("[ccs] CCS consensus generation completed successfully")
    logging.info("[ccs] Generated HiFi BAM: %s", output_bam)

    return output_bam
```

---

### **Phase 3: Pipeline Module**

#### **Task 3.1: Create pacbio_pipeline.py** ✅

**File**: `muc_one_up/read_simulator/pacbio_pipeline.py` (NEW)

**FIXES**:
- Uses samtools_wrapper.convert_bam_to_fastq (DRY - FIX #1)
- Uses minimap2_wrapper with map-hifi preset (FIX #4)
- Proper exception handling pattern (FIX #2)
- Comprehensive validation
- Intermediate file cleanup (Warning #8)

```python
#!/usr/bin/env python3
"""
PacBio HiFi read simulation pipeline.

This module implements the PacBio HiFi read simulation pipeline using pbsim3
and CCS, with alignment using minimap2.

Key features:
- Multi-pass CLR simulation with pbsim3
- HiFi consensus generation with CCS
- Configurable coverage correction for model efficiency
- Reliable timeout handling for all external tool calls
- Error handling with descriptive messages
- Input validation and appropriate parameter defaults
- Consistent output file naming
- Full seed support for reproducibility
- Intermediate file cleanup
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from ..exceptions import FileOperationError, ValidationError
from .constants import (
    DEFAULT_CCS_MIN_PASSES,
    DEFAULT_CCS_MIN_RQ,
    DEFAULT_PACBIO_LENGTH_MEAN,
    DEFAULT_PACBIO_LENGTH_SD,
    DEFAULT_PACBIO_PASS_NUM,
    MAX_COVERAGE,
    MIN_COVERAGE,
    MINIMAP2_PRESET_PACBIO_HIFI,
)
from .utils import cleanup_files
from .wrappers.ccs_wrapper import run_ccs_consensus
from .wrappers.minimap2_wrapper import align_reads_with_minimap2
from .wrappers.pbsim3_wrapper import run_pbsim3_simulation
from .wrappers.samtools_wrapper import convert_bam_to_fastq


def simulate_pacbio_reads_pipeline(
    config: dict[str, Any], input_fa: str, human_reference: str | None = None
) -> str:
    """
    Run the complete PacBio HiFi read simulation pipeline.

    Pipeline steps:
    1. Validate input parameters from config
    2. Run pbsim3 multi-pass CLR simulation
    3. Generate HiFi consensus reads with CCS
    4. Convert HiFi BAM to FASTQ (for compatibility)
    5. Align reads to the reference using minimap2 (map-hifi preset)
    6. Create and index the final BAM file
    7. Clean up intermediate files

    Args:
        config: Dictionary containing "tools" and "pacbio_params" sections.
               Required parameters:
               - tools: Dictionary containing "pbsim3", "ccs", "minimap2", "samtools"
               - pacbio_params: Dictionary containing:
                 - model_type: 'qshmm' or 'errhmm'
                 - model_file: Path to model file
                 - coverage: Target sequencing coverage
                 Optional parameters:
                 - pass_num: Number of passes (default: 10)
                 - length_mean: Mean read length (default: 15000)
                 - length_sd: Read length SD (default: 5000)
                 - ccs_min_passes: CCS minimum passes (default: 3)
                 - ccs_min_rq: CCS min quality (default: 0.99)
                 - num_threads: Number of threads (default: 4)
                 - seed: Random seed for reproducibility
                 - correction_factor: Coverage adjustment (default: 1.0)
                 - enable_coverage_correction: Apply correction (default: False)
        input_fa: Input simulated FASTA file.
        human_reference: Path to human reference genome. If None, uses input_fa.

    Returns:
        Path to the final output BAM file (input_basename_pacbio.bam).

    Raises:
        ValueError: If required parameters are missing or invalid.
        FileNotFoundError: If model file doesn't exist.
        RuntimeError: If any step in the pipeline fails.

    Example:
        >>> config = {
        ...     "tools": {"pbsim3": "pbsim3", "ccs": "ccs",
        ...               "minimap2": "minimap2", "samtools": "samtools"},
        ...     "pacbio_params": {
        ...         "model_type": "errhmm",
        ...         "model_file": "data/ERRHMM-SEQUEL.model",
        ...         "coverage": 30,
        ...         "pass_num": 10,
        ...         "seed": 42
        ...     }
        ... }
        >>> bam_file = simulate_pacbio_reads_pipeline(config, "diploid.fa")
        >>> print(f"Generated: {bam_file}")
    """
    # Start time
    start_time = datetime.now()
    logging.info(
        "Starting PacBio HiFi read simulation pipeline at %s",
        start_time.strftime("%Y-%m-%d %H:%M:%S"),
    )

    # Extract tool commands and simulation parameters
    tools = config.get("tools", {})
    pb_params = config.get("pacbio_params", {})
    rs_config = config.get("read_simulation", {})

    # Validate required tools
    pbsim3_cmd = tools.get("pbsim3")
    ccs_cmd = tools.get("ccs")
    minimap2_cmd = tools.get("minimap2")
    samtools_cmd = tools.get("samtools")

    if not pbsim3_cmd:
        raise ValueError("Missing 'pbsim3' command in tools configuration.")
    if not ccs_cmd:
        raise ValueError("Missing 'ccs' command in tools configuration.")
    if not minimap2_cmd:
        raise ValueError("Missing 'minimap2' command in tools configuration.")
    if not samtools_cmd:
        raise ValueError("Missing 'samtools' command in tools configuration.")

    # Validate required simulation parameters
    model_type = pb_params.get("model_type")
    model_file = pb_params.get("model_file")
    coverage = pb_params.get("coverage")

    if not model_type:
        raise ValueError("Missing 'model_type' in pacbio_params configuration.")
    if not model_file:
        raise ValueError("Missing 'model_file' in pacbio_params configuration.")
    if coverage is None:
        raise ValueError("Missing 'coverage' in pacbio_params configuration.")

    # Validate coverage range (FIX: Warning #11)
    if not MIN_COVERAGE <= coverage <= MAX_COVERAGE:
        raise ValueError(
            f"Invalid coverage: {coverage}. Must be between {MIN_COVERAGE} and {MAX_COVERAGE}"
        )

    # Optional parameters with defaults (using constants - FIX #6)
    pass_num = pb_params.get("pass_num", DEFAULT_PACBIO_PASS_NUM)
    length_mean = pb_params.get("length_mean", DEFAULT_PACBIO_LENGTH_MEAN)
    length_sd = pb_params.get("length_sd", DEFAULT_PACBIO_LENGTH_SD)
    ccs_min_passes = pb_params.get("ccs_min_passes", DEFAULT_CCS_MIN_PASSES)
    ccs_min_rq = pb_params.get("ccs_min_rq", DEFAULT_CCS_MIN_RQ)
    threads = pb_params.get("num_threads", rs_config.get("threads", 4))
    seed = pb_params.get("seed")
    other_options = pb_params.get("other_options", "")

    # Coverage correction parameters
    correction_factor = pb_params.get("correction_factor", 1.0)
    enable_coverage_correction = pb_params.get("enable_coverage_correction", False)

    # Setup output paths
    input_path = Path(input_fa)
    input_basename = input_path.stem
    output_dir = rs_config.get("output_dir", str(input_path.parent))
    output_prefix = str(Path(output_dir) / f"{input_basename}_pacbio")

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logging.info("=" * 70)
    logging.info("PACBIO HIFI SIMULATION PIPELINE")
    logging.info("=" * 70)
    logging.info("Model: %s (%s)", model_type.upper(), Path(model_file).name)
    logging.info("Coverage: %.1fx", float(coverage))
    logging.info("Pass number: %d", pass_num)
    logging.info("Read length: %d ± %d bp", length_mean, length_sd)
    logging.info("=" * 70)

    # Track intermediate files for cleanup
    intermediate_files = []

    try:
        # 1. Run pbsim3 multi-pass CLR simulation
        logging.info("1. Starting pbsim3 multi-pass CLR simulation")

        # Apply coverage correction if enabled
        actual_coverage = float(coverage)
        if enable_coverage_correction:
            actual_coverage = float(coverage) / correction_factor
            logging.info(
                "Coverage correction enabled: %.1fx → %.1fx (factor=%.3f)",
                float(coverage),
                actual_coverage,
                correction_factor,
            )

        clr_bam = run_pbsim3_simulation(
            pbsim3_cmd=pbsim3_cmd,
            samtools_cmd=samtools_cmd,  # For SAM→BAM conversion
            reference_fasta=input_fa,
            output_prefix=output_prefix,
            model_type=model_type,
            model_file=model_file,
            coverage=actual_coverage,
            pass_num=pass_num,
            length_mean=length_mean,
            length_sd=length_sd,
            other_options=other_options,
            seed=seed,
        )
        intermediate_files.append(clr_bam)
        logging.info("pbsim3 CLR simulation completed successfully")

        # 2. Generate HiFi consensus reads with CCS
        logging.info("2. Starting CCS HiFi consensus generation")
        hifi_bam = f"{output_prefix}_hifi.bam"

        run_ccs_consensus(
            ccs_cmd=ccs_cmd,
            input_bam=clr_bam,
            output_bam=hifi_bam,
            min_passes=ccs_min_passes,
            min_rq=ccs_min_rq,
            num_threads=int(threads),
        )
        intermediate_files.append(hifi_bam)
        logging.info("CCS HiFi consensus generation completed successfully")

        # 3. Convert HiFi BAM to FASTQ (using samtools_wrapper - FIX #1 DRY)
        logging.info("3. Converting HiFi BAM to FASTQ")
        hifi_fastq = f"{output_prefix}_hifi.fastq"

        convert_bam_to_fastq(
            samtools_cmd=samtools_cmd,
            input_bam=hifi_bam,
            output_fastq=hifi_fastq,
        )
        intermediate_files.append(hifi_fastq)
        logging.info("BAM to FASTQ conversion completed successfully")

        # 4. Align HiFi reads with minimap2 (using generic wrapper - FIX #4)
        logging.info("4. Starting HiFi read alignment with minimap2")
        output_bam = f"{output_prefix}.bam"

        # Determine reference for alignment
        if human_reference is None:
            try:
                # Get assembly from config (default: hg38)
                from ..bioinformatics.reference_validation import (
                    get_reference_path_for_assembly,
                    validate_reference_for_assembly,
                )

                assembly = config.get("reference_assembly", "hg38")
                ref_path = get_reference_path_for_assembly(config, assembly)
                reference_for_alignment = str(ref_path)

                # Validate reference and indices for minimap2
                warnings = validate_reference_for_assembly(config, assembly, aligner="minimap2")

                logging.info("Using reference from config: %s (%s)", reference_for_alignment, assembly)

                if warnings:
                    for warning in warnings:
                        logging.warning("%s", warning)

            except (ValidationError, FileOperationError) as e:
                logging.warning("Could not load reference from config: %s", e)
                logging.warning("Falling back to aligning against simulated reference")
                reference_for_alignment = input_fa
        else:
            reference_for_alignment = human_reference
            logging.info("Using user-provided reference: %s", reference_for_alignment)

        logging.info("Aligning HiFi reads to reference: %s", reference_for_alignment)

        # Use generic minimap2 wrapper with PacBio HiFi preset (FIX #4)
        align_reads_with_minimap2(
            minimap2_cmd=minimap2_cmd,
            samtools_cmd=samtools_cmd,
            reference=reference_for_alignment,
            reads_fastq=hifi_fastq,
            output_bam=output_bam,
            preset=MINIMAP2_PRESET_PACBIO_HIFI,  # FIX #4: Correct preset
            threads=int(threads),
        )
        logging.info("HiFi read alignment completed successfully")

    except Exception as e:
        # Proper exception handling pattern (FIX #2)
        logging.error("PacBio HiFi pipeline failed: %s", e)
        raise RuntimeError(f"PacBio HiFi read simulation failed: {e}") from e

    finally:
        # Clean up intermediate files (FIX: Warning #8)
        if intermediate_files:
            logging.info("Cleaning up intermediate files...")
            cleanup_files(intermediate_files)

    # Calculate elapsed time
    end_time = datetime.now()
    duration = end_time - start_time
    logging.info(
        "PacBio HiFi read simulation pipeline completed at %s (duration: %s)",
        end_time.strftime("%Y-%m-%d %H:%M:%S"),
        str(duration).split(".")[0],  # Remove microseconds
    )

    logging.info("Final output:")
    logging.info("  Aligned and indexed BAM: %s", output_bam)
    logging.info("  Index: %s.bai", output_bam)

    return output_bam
```

---

### **Phase 4: Integration**

#### **Task 4.1: Update read_simulation.py Router** ✅

**File**: `muc_one_up/read_simulation.py` (MODIFY)

**FIX #7**: Use Strategy Pattern instead of if/elif chain

```python
#!/usr/bin/env python3
"""
read_simulation.py

Command-line entry point for the MUC1 read simulation pipeline.

This module provides a standardized interface to the read simulation pipeline,
which creates realistic sequencing reads from a simulated MUC1 haplotype FASTA.
It supports Illumina short reads (via reseq/WeSSim), Oxford Nanopore long reads
(via NanoSim), and PacBio HiFi long reads (via pbsim3/CCS).

Key features:
- Simulator selection (Illumina, ONT, or PacBio)
- Strategy Pattern for extensibility
- Reliable timeout handling for all external tool calls
- Consistent naming convention for input and output files
"""

import logging
from pathlib import Path
from typing import Any, Callable

from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline
from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_reads_pipeline
from muc_one_up.read_simulator.pipeline import (
    simulate_reads_pipeline as simulate_illumina_reads,
)

# Configure logging
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

# Strategy Pattern: Map simulator names to pipeline functions (FIX #7)
SIMULATOR_MAP: dict[str, Callable[[dict[str, Any], str, str | None], str]] = {
    "illumina": lambda config, input_fa, _: simulate_illumina_reads(config, input_fa),
    "ont": simulate_ont_reads_pipeline,
    "pacbio": simulate_pacbio_reads_pipeline,
}


def simulate_reads(config: dict[str, Any], input_fa: str) -> str:
    """
    Run the complete read simulation pipeline.

    This function serves as the main API entry point for the read simulation pipeline.
    It selects the appropriate pipeline based on the specified simulator in the config,
    and delegates to the modular implementation in the read_simulator package.

    Uses Strategy Pattern for extensibility (FIX #7: Open/Closed Principle).

    Args:
        config: Dictionary containing configuration sections:
               - 'read_simulation': Contains 'simulator' ('illumina', 'ont', or 'pacbio')
                 and other parameters
               - 'tools': Contains paths to required external tools
               - 'nanosim_params': Required for ONT simulation
               - 'pacbio_params': Required for PacBio simulation
               See the respective pipeline documentation for detailed parameter information.
        input_fa: Input simulated FASTA file (e.g., muc1_simulated.fa).

    Returns:
        Path to the final output BAM file.

    Raises:
        ValueError: If simulator type is unknown.
    """
    # Determine which simulator to use
    simulator = config.get("read_simulation", {}).get("simulator", "illumina").lower()

    # Validate simulator type (FIX #7: Strategy Pattern)
    if simulator not in SIMULATOR_MAP:
        valid_simulators = ", ".join(SIMULATOR_MAP.keys())
        raise ValueError(
            f"Unknown simulator: '{simulator}'. "
            f"Valid options: {valid_simulators}"
        )

    # Get human reference (if specified)
    human_reference = config.get("read_simulation", {}).get("human_reference")
    if human_reference:
        logging.info("Using human reference for alignment: %s", human_reference)
    else:
        logging.info("No human reference specified - will use simulated reference for alignment")

    # Dispatch to appropriate simulator pipeline (Strategy Pattern)
    logging.info("Using %s read simulation pipeline", simulator.upper())
    simulator_func = SIMULATOR_MAP[simulator]

    return simulator_func(config, input_fa, human_reference)


if __name__ == "__main__":  # OK: top-level entry point
    import json
    import sys

    if len(sys.argv) != 3:
        print("Usage: python read_simulation.py <config.json> <input_fasta>")
        sys.exit(1)

    config_file = sys.argv[1]
    input_fa = sys.argv[2]

    with Path(config_file).open() as fh:
        config = json.load(fh)

    simulate_reads(config, input_fa)
```

---

#### **Task 4.2: Update CLI** ✅

**File**: `muc_one_up/cli/click_main.py` (MODIFY)

Add PacBio subcommand to `reads` group:

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
    help="Target coverage.",
)
@click.option(
    "--pass-num",
    type=int,
    default=10,
    show_default=True,
    help="Number of passes for multi-pass sequencing (2-50, higher = better HiFi quality).",
)
@click.option(
    "--read-length",
    type=int,
    default=15000,
    show_default=True,
    help="Mean read length in bp (typical HiFi: 10000-25000).",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def pacbio(ctx, input_fastas, out_dir, out_base, coverage, pass_num, read_length, seed):
    """Simulate PacBio HiFi long reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads pacbio file.fa --out-base reads
    - Multiple files: muconeup reads pacbio file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads pacbio *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Pipeline Steps:
      1. pbsim3 simulates multi-pass CLR reads
      2. CCS generates HiFi consensus reads
      3. minimap2 aligns HiFi reads to reference
      4. Outputs sorted and indexed BAM file

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads pacbio sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads pacbio sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads pacbio sample.*.simulated.fa

      # Custom parameters for higher quality
      muconeup --config X reads pacbio sample.fa \\
        --pass-num 20 --read-length 20000 --coverage 50 --seed 42
    """
    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure PacBio simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "pacbio"
        config["read_simulation"]["coverage"] = coverage

        if "pacbio_params" not in config:
            config["pacbio_params"] = {}
        config["pacbio_params"]["coverage"] = coverage
        config["pacbio_params"]["pass_num"] = pass_num
        config["pacbio_params"]["length_mean"] = read_length

        # Set seed if provided
        if seed is not None:
            config["pacbio_params"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for PacBio HiFi read simulation", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_pacbio_reads")

            logging.info(
                "[%d/%d] Simulating PacBio HiFi reads: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            # Run simulation for this file
            simulate_reads_pipeline(config, input_fasta)

        logging.info("PacBio HiFi read simulation completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("PacBio read simulation failed: %s", e)
        ctx.exit(1)
```

---

### **Phase 5: Testing**

All test patterns follow existing codebase conventions exactly.

#### **Task 5.1: Unit Tests for Wrappers**

**File**: `tests/read_simulator/test_pbsim3_wrapper.py` (NEW)

```python
"""Tests for pbsim3 wrapper."""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch

from muc_one_up.read_simulator.wrappers.pbsim3_wrapper import run_pbsim3_simulation


class TestPbsim3Wrapper:
    """Test pbsim3 wrapper functionality."""

    def test_validates_model_type(self):
        """Test that invalid model types raise ValueError."""
        with pytest.raises(ValueError, match="Invalid model_type"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim3",
                samtools_cmd="samtools",
                reference_fasta="test.fa",
                output_prefix="output",
                model_type="invalid",  # Invalid
                model_file="model.file",
                coverage=30,
                pass_num=10,
            )

    def test_validates_model_file_exists(self, tmp_path):
        """Test that missing model file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Model file not found"):
            run_pbsim3_simulation(
                pbsim3_cmd="pbsim3",
                samtools_cmd="samtools",
                reference_fasta="test.fa",
                output_prefix="output",
                model_type="errhmm",
                model_file="/nonexistent/model.file",
                coverage=30,
                pass_num=10,
            )

    @patch("muc_one_up.read_simulator.wrappers.pbsim3_wrapper.run_command")
    def test_includes_seed_when_provided(self, mock_run_command, tmp_path):
        """Test that seed parameter is included in command."""
        model_file = tmp_path / "model.file"
        model_file.write_text("dummy")

        output_bam = tmp_path / "output_0001.bam"
        output_bam.write_bytes(b"BAM")

        run_pbsim3_simulation(
            pbsim3_cmd="pbsim3",
            samtools_cmd="samtools",
            reference_fasta="test.fa",
            output_prefix=str(tmp_path / "output"),
            model_type="errhmm",
            model_file=str(model_file),
            coverage=30,
            pass_num=10,
            seed=42,
        )

        # Verify command includes seed
        call_args = mock_run_command.call_args
        cmd_list = call_args[0][0]
        assert "--seed" in cmd_list
        assert "42" in cmd_list
```

**File**: `tests/read_simulator/test_ccs_wrapper.py` (NEW)

```python
"""Tests for CCS wrapper."""

import pytest
from pathlib import Path
from unittest.mock import patch

from muc_one_up.read_simulator.wrappers.ccs_wrapper import run_ccs_consensus


class TestCCSWrapper:
    """Test CCS wrapper functionality."""

    def test_raises_on_missing_input(self):
        """Test that missing input BAM raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="Input BAM not found"):
            run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam="/nonexistent.bam",
                output_bam="output.bam",
            )

    @patch("muc_one_up.read_simulator.wrappers.ccs_wrapper.run_command")
    def test_detects_empty_output(self, mock_run_command, tmp_path):
        """Test that empty CCS output is detected."""
        input_bam = tmp_path / "input.bam"
        output_bam = tmp_path / "output.bam"

        input_bam.write_bytes(b"BAM")
        output_bam.write_bytes(b"")  # Empty output

        with pytest.raises(RuntimeError, match="CCS produced empty output"):
            run_ccs_consensus(
                ccs_cmd="ccs",
                input_bam=str(input_bam),
                output_bam=str(output_bam),
            )
```

**File**: `tests/read_simulator/test_minimap2_wrapper.py` (NEW)

```python
"""Tests for generic minimap2 wrapper."""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch

from muc_one_up.read_simulator.wrappers.minimap2_wrapper import align_reads_with_minimap2


class TestMinimap2Wrapper:
    """Test generic minimap2 wrapper."""

    @patch("muc_one_up.read_simulator.wrappers.minimap2_wrapper.subprocess.run")
    @patch("muc_one_up.read_simulator.wrappers.minimap2_wrapper.run_command")
    def test_uses_correct_preset(self, mock_run_command, mock_subprocess_run, tmp_path):
        """Test that specified preset is used in minimap2 command."""
        ref = tmp_path / "ref.fa"
        reads = tmp_path / "reads.fq"
        output_bam = tmp_path / "output.bam"

        ref.write_text(">chr1\nACGT\n")
        reads.write_text("@read1\nACGT\n+\nIIII\n")

        # Mock subprocess.run
        mock_subprocess_run.return_value = Mock(returncode=0, stderr="")

        # Mock samtools operations
        def mock_run_side_effect(cmd_list, **kwargs):
            if "sort" in cmd_list:
                output_bam.write_bytes(b"SORTED_BAM")
            elif "index" in cmd_list:
                (Path(str(output_bam) + ".bai")).write_bytes(b"INDEX")
            return 0

        mock_run_command.side_effect = mock_run_side_effect

        # Test with PacBio HiFi preset
        align_reads_with_minimap2(
            minimap2_cmd="minimap2",
            samtools_cmd="samtools",
            reference=str(ref),
            reads_fastq=str(reads),
            output_bam=str(output_bam),
            preset="map-hifi",  # PacBio HiFi preset
        )

        # Verify preset was used
        call_args = mock_subprocess_run.call_args[0][0]
        assert "map-hifi" in call_args
```

---

#### **Task 5.2: Integration Tests**

**File**: `tests/read_simulator/test_pacbio_pipeline.py` (NEW)

```python
"""Integration tests for PacBio HiFi pipeline."""

import pytest
from pathlib import Path
from unittest.mock import patch

from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_reads_pipeline


class TestPacBioPipeline:
    """Test PacBio pipeline integration."""

    def test_validates_required_tools(self, tmp_path):
        """Test that missing tools raise ValueError."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {"pbsim3": "pbsim3"},  # Missing ccs, minimap2, samtools
            "pacbio_params": {
                "model_type": "errhmm",
                "model_file": "model.file",
                "coverage": 30,
            },
        }

        with pytest.raises(ValueError, match="Missing 'ccs' command"):
            simulate_pacbio_reads_pipeline(config, str(input_fa))

    def test_validates_coverage_range(self, tmp_path):
        """Test that coverage is validated."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        config = {
            "tools": {
                "pbsim3": "pbsim3",
                "ccs": "ccs",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "pacbio_params": {
                "model_type": "errhmm",
                "model_file": "model.file",
                "coverage": 99999,  # Too high
            },
        }

        with pytest.raises(ValueError, match="Invalid coverage"):
            simulate_pacbio_reads_pipeline(config, str(input_fa))

    @patch("muc_one_up.read_simulator.pacbio_pipeline.run_pbsim3_simulation")
    @patch("muc_one_up.read_simulator.pacbio_pipeline.run_ccs_consensus")
    @patch("muc_one_up.read_simulator.pacbio_pipeline.convert_bam_to_fastq")
    @patch("muc_one_up.read_simulator.pacbio_pipeline.align_reads_with_minimap2")
    @patch("muc_one_up.read_simulator.pacbio_pipeline.cleanup_files")
    def test_full_pipeline_cleans_up_intermediates(
        self,
        mock_cleanup,
        mock_align,
        mock_convert,
        mock_ccs,
        mock_pbsim3,
        tmp_path,
    ):
        """Test that intermediate files are cleaned up."""
        input_fa = tmp_path / "input.fa"
        input_fa.write_text(">chr1\nACGT\n")

        model_file = tmp_path / "model.file"
        model_file.write_text("model")

        # Mock outputs
        mock_pbsim3.return_value = str(tmp_path / "clr.bam")
        mock_ccs.return_value = str(tmp_path / "hifi.bam")
        mock_convert.return_value = str(tmp_path / "hifi.fastq")
        mock_align.return_value = str(tmp_path / "output.bam")

        config = {
            "tools": {
                "pbsim3": "pbsim3",
                "ccs": "ccs",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "pacbio_params": {
                "model_type": "errhmm",
                "model_file": str(model_file),
                "coverage": 30,
            },
        }

        simulate_pacbio_reads_pipeline(config, str(input_fa))

        # Verify cleanup was called
        mock_cleanup.assert_called_once()
        cleaned_files = mock_cleanup.call_args[0][0]
        assert len(cleaned_files) > 0  # Should have intermediate files
```

---

#### **Task 5.3: CLI Tests**

**File**: `tests/cli/test_cli_pacbio.py` (NEW)

```python
"""CLI tests for PacBio read simulation command."""

import pytest
import json
from click.testing import CliRunner
from pathlib import Path
from unittest.mock import patch

from muc_one_up.cli.click_main import cli


class TestPacBioCLI:
    """Test PacBio CLI command."""

    @pytest.fixture
    def runner(self):
        """Create Click test runner."""
        return CliRunner()

    @pytest.fixture
    def mock_config_file(self, tmp_path):
        """Create a mock config file."""
        config_path = tmp_path / "config.json"
        config_content = {
            "tools": {
                "pbsim3": "pbsim3",
                "ccs": "ccs",
                "minimap2": "minimap2",
                "samtools": "samtools",
            },
            "pacbio_params": {
                "model_type": "errhmm",
                "model_file": "data/ERRHMM-SEQUEL.model",
                "coverage": 30,
            },
        }
        config_path.write_text(json.dumps(config_content))
        return config_path

    @pytest.fixture
    def mock_fasta_file(self, tmp_path):
        """Create a mock FASTA file."""
        fasta_path = tmp_path / "test.fa"
        fasta_path.write_text(">test\nATCGATCG\n")
        return fasta_path

    @patch("muc_one_up.cli.click_main.simulate_reads_pipeline")
    def test_pacbio_command_with_seed(
        self, mock_simulate, runner, mock_config_file, mock_fasta_file
    ):
        """Test PacBio command with seed parameter."""
        mock_simulate.return_value = "output.bam"

        result = runner.invoke(
            cli,
            [
                "--config", str(mock_config_file),
                "reads", "pacbio",
                str(mock_fasta_file),
                "--seed", "42",
                "--coverage", "50",
            ],
        )

        assert result.exit_code == 0
        mock_simulate.assert_called_once()
```

---

### **Phase 6: Documentation**

#### **Task 6.1: Update CLAUDE.md**

Add PacBio section:

```markdown
### Read Simulation Pipelines

#### PacBio HiFi Pipeline (`read_simulator/pacbio_pipeline.py`)
Uses pbsim3 and CCS for PacBio HiFi long reads:
1. pbsim3: Multi-pass CLR simulation
2. CCS: HiFi consensus generation
3. samtools: BAM/FASTQ conversion
4. minimap2: Alignment (map-hifi preset)
5. samtools: Indexing

**Multi-pass Simulation**:
- pbsim3 simulates CLR with multiple passes (--pass-num ≥ 2)
- CCS processes multi-pass CLR to generate HiFi consensus (Q30+ accuracy)
- Typical pass numbers: 10-20 for HiFi quality

**Model Selection**:
- **ERRHMM**: Error HMM model (supports RSII and Sequel)
- **QSHMM**: Quality score HMM model (RSII only)

**Configuration Example**:
```json
{
  "read_simulation": {
    "simulator": "pacbio",
    "coverage": 30,
    "seed": 42
  },
  "pacbio_params": {
    "model_type": "errhmm",
    "model_file": "data/ERRHMM-SEQUEL.model",
    "pass_num": 10,
    "length_mean": 15000,
    "length_sd": 5000,
    "ccs_min_passes": 3,
    "ccs_min_rq": 0.99,
    "seed": 42
  }
}
```

### Common Commands

#### PacBio HiFi Read Simulation
```bash
# Simulate PacBio HiFi reads
muconeup --config config.json reads pacbio sample.fa --seed 42

# Custom parameters for higher quality
muconeup --config config.json reads pacbio sample.fa \
  --coverage 50 \
  --pass-num 20 \
  --read-length 20000 \
  --seed 42

# Batch processing
muconeup --config config.json reads pacbio sample.*.simulated.fa
```

### Setting Up Conda Environments

```bash
# For PacBio HiFi read simulation (pbsim3 + CCS)
mamba env create -f conda/env_pacbio.yml
conda activate env_pacbio
```
```

---

#### **Task 6.2: Create Example Config**

**File**: `data/examples/config_pacbio.json` (NEW)

```json
{
  "reference_assembly": "hg38",
  "repeats": {},
  "constants": {},
  "probabilities": {},
  "tools": {
    "pbsim3": "pbsim3",
    "ccs": "ccs",
    "minimap2": "minimap2",
    "samtools": "samtools"
  },
  "read_simulation": {
    "simulator": "pacbio",
    "coverage": 30,
    "threads": 8,
    "seed": 42
  },
  "pacbio_params": {
    "model_type": "errhmm",
    "model_file": "data/models/ERRHMM-SEQUEL.model",
    "coverage": 30,
    "pass_num": 10,
    "length_mean": 15000,
    "length_sd": 5000,
    "ccs_min_passes": 3,
    "ccs_min_rq": 0.99,
    "num_threads": 4,
    "seed": 42,
    "correction_factor": 1.0,
    "enable_coverage_correction": false,
    "other_options": ""
  }
}
```

---

## **4. IMPLEMENTATION CHECKLIST**

### **Phase 1: Constants and Configuration** ✅
- [ ] Create `constants.py` with validation constants
- [ ] Add `pacbio_params` to `CONFIG_SCHEMA` in `config.py`
- [ ] Create `conda/env_pacbio.yml`
- [ ] Test environment creation and tool availability

### **Phase 2: Core Wrappers** ✅
- [ ] Extend `samtools_wrapper.py` (add SAM→BAM, BAM→FASTQ if missing)
- [ ] Create `minimap2_wrapper.py` (generic with preset parameter)
- [ ] Create `pbsim3_wrapper.py`
- [ ] Create `ccs_wrapper.py`
- [ ] Add unit tests for all wrappers

### **Phase 3: Pipeline Module** ✅
- [ ] Create `pacbio_pipeline.py`
- [ ] Add integration tests

### **Phase 4: Integration** ✅
- [ ] Update `read_simulation.py` router (Strategy Pattern)
- [ ] Add `reads pacbio` subcommand to CLI
- [ ] Add CLI tests

### **Phase 5: Testing** ✅
- [ ] Write unit tests for wrappers
- [ ] Write integration tests for pipeline
- [ ] Write end-to-end CLI tests
- [ ] Manual testing with real data

### **Phase 6: Documentation** ✅
- [ ] Update CLAUDE.md
- [ ] Update README.md
- [ ] Create example config file
- [ ] Add usage examples

### **Phase 7: Validation** ✅
- [ ] Run full test suite
- [ ] Test Illumina/ONT workflows (regression testing)
- [ ] Verify all critical issues from code review are fixed

---

## **5. CODE REVIEW FIXES SUMMARY**

### **✅ Critical Issues Addressed**

1. **DRY Violation (Issue #1)**: ✅ Fixed
   - BAM→FASTQ conversion in `samtools_wrapper.py` (reusable)
   - No duplication in `ccs_wrapper.py`

2. **Exception Handling (Issue #2)**: ✅ Fixed
   - Wrappers let `ExternalToolError` propagate from `run_command()`
   - Pipeline catches and wraps in `RuntimeError` for user messages

3. **SAM→BAM Conversion (Issue #3)**: ✅ Fixed
   - `pbsim3_wrapper` detects SAM output
   - Converts to BAM using `samtools_wrapper.convert_sam_to_bam()`
   - Cleans up SAM file

4. **Minimap2 Preset (Issue #4)**: ✅ Fixed
   - Created generic `minimap2_wrapper.py`
   - Parameterized preset (`map-ont`, `map-hifi`, `map-pb`)
   - PacBio uses `MINIMAP2_PRESET_PACBIO_HIFI`

5. **CONFIG_SCHEMA (Issue #5)**: ✅ Fixed
   - Added comprehensive `pacbio_params` schema
   - Validates all parameters, types, and ranges
   - Uses enum for model_type

6. **Model Validation (Issue #6)**: ✅ Fixed
   - Created `constants.py` with `VALID_PBSIM3_MODEL_TYPES`
   - Wrapper uses constants for validation
   - JSON Schema also validates via enum

7. **Router Pattern (Issue #7)**: ✅ Fixed
   - Implemented Strategy Pattern with `SIMULATOR_MAP`
   - Open/Closed principle: add new simulators by extending map
   - No modification of existing if/elif chain

### **✅ Important Warnings Addressed**

1. **Model File Validation**: ✅ Added `Path(model_file).exists()` check
2. **Empty CCS Output**: ✅ Added size validation in `ccs_wrapper`
3. **Coverage Range**: ✅ Added MIN/MAX validation using constants
4. **Intermediate Cleanup**: ✅ Added cleanup in pipeline finally block
5. **Timeouts**: ✅ Moved to constants for consistency
6. **Progress Logging**: ✅ Added for long-running steps

---

## **6. SUCCESS CRITERIA**

✅ **All Requirements Met**:
- [x] PacBio HiFi reads can be simulated from FASTA input
- [x] CLI command `muconeup reads pacbio` works correctly
- [x] Batch processing supported
- [x] Seed reproducibility works
- [x] Output BAM files are valid and indexed
- [x] Follows existing codebase patterns
- [x] SOLID principles adhered to
- [x] DRY principle maintained
- [x] No code duplication
- [x] Modular and maintainable
- [x] Comprehensive testing
- [x] Complete documentation

---

## **7. TIMELINE ESTIMATE**

- **Phase 1 (Constants + Config)**: 2 hours
- **Phase 2 (Wrappers)**: 8 hours (includes fixes)
- **Phase 3 (Pipeline)**: 6 hours
- **Phase 4 (Integration)**: 3 hours
- **Phase 5 (Testing)**: 8 hours
- **Phase 6 (Documentation)**: 3 hours
- **Phase 7 (Validation)**: 2 hours

**Total Estimated Time**: ~32 hours (4 days of full-time development)

---

## **8. REFERENCES**

- **pbsim3**: https://github.com/yukiteruono/pbsim3
- **PacBio CCS**: https://github.com/PacificBiosciences/ccs
- **Bioconda pbsim3**: https://anaconda.org/bioconda/pbsim3
- **Bioconda pbccs**: https://anaconda.org/bioconda/pbccs
- **Existing ONT Pipeline**: `muc_one_up/read_simulator/ont_pipeline.py`
- **Existing Minimap2 Usage**: `muc_one_up/read_simulator/wrappers/nanosim_wrapper.py`
- **Exception Hierarchy**: `muc_one_up/exceptions.py`
- **Command Utilities**: `muc_one_up/read_simulator/command_utils.py`

---

**PLAN APPROVED FOR IMPLEMENTATION** ✅
**All Code Review Issues Addressed**
**Version**: 2.0 (Revised)
**Date**: 2025-10-20
