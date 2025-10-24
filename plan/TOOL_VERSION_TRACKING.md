# Tool Version Tracking Implementation Guide

**Project:** MucOneUp v0.15.0
**Priority:** HIGH (Target: v0.16.0)
**Effort:** 6-10 hours
**Status:** ❌ Not Implemented

---

## 1. Executive Summary

### Current State
MucOneUp uses **9 external bioinformatics tools** without version tracking:

| Tool | Purpose | Platform | Flag | Exit | Binary | Status |
|------|---------|----------|------|------|--------|--------|
| samtools | BAM operations | All | `--version` | 0 | samtools | ✅ |
| minimap2 | Long-read align | ONT/PacBio | `--version` | 0 | minimap2 | ✅ |
| NanoSim | ONT simulation | ONT | `--version` | 0 | **simulator.py** | ✅ |
| pbsim3 | PacBio CLR | PacBio | *(none)* | 255 | **pbsim** ⚠️ | ✅ |
| ccs | PacBio HiFi | PacBio | `--version` | 0 | ccs | ✅ |
| reseq | Illumina sim | Illumina | `--version` | 0 | reseq | ✅ |
| bwa | Short-read align | Illumina | *(none)* | 1 | bwa | ✅ |
| faToTwoBit | FASTA→2bit | Illumina | *(none)* | 255 | faToTwoBit ⚠️ | ✅ |
| pblat | Parallel BLAT | Illumina | *(none)* | 0 | pblat | ✅ |

**⚠️ Critical:** pbsim3 binary is `pbsim` with NO version output; faToTwoBit has NO version info

**Impact:** No reproducibility, debugging challenges, no audit trail.

### Recommended Solution

**Adopt VariantCentrifuge pattern** (battle-tested in clinical pipelines):
- Simple string return: `"samtools 1.17"` or `"N/A"`
- Tool-specific parsers
- Graceful degradation (never crashes)
- Separate TSV metadata file
- KISS principle throughout

---

## 2. Research Summary

### Best Practices from Literature

**Source: Nature Methods 2021 - Reproducible Workflows**
- Modern workflow managers (Nextflow, Snakemake) track tool versions automatically
- Containerization (Docker/Singularity) ensures version consistency
- Provenance metadata is critical for reproducibility

**Source: BMC Bioinformatics 2017 - Provenance Tracking**
- W3C PROV standard for provenance representation
- Retrospective vs. prospective provenance
- Tool versioning as minimum requirement

**Source: PHA4GE Best Practices 2024**
- Code repositories (GitHub/GitLab) for version control
- Automated version capture in CI/CD pipelines
- Quality management standards for clinical applications

### Python Subprocess Best Practices

**Key Findings:**
1. **Try multiple flags:** Different tools use `-v`, `--version`, `-version`, or no flag
2. **Use `check=False`:** Many tools return non-zero for version commands
3. **Parse both stdout AND stderr:** Version info can appear in either stream
4. **Set timeout:** Prevent hanging on tools with interactive prompts
5. **Security:** Never use `shell=True` - always pass command as list

### Real-World Reference Implementation

**VariantCentrifuge** (clinical variant analysis pipeline):
- **Simple design:** One function, string return, tool-specific parsers
- **Graceful failure:** Returns `"N/A"`, logs warning, never crashes
- **TSV output:** Separate metadata file, not mixed with results
- **Production-proven:** Running in clinical diagnostics since 2020+

---

## 3. Tool Version Detection Protocol

### 3.1 Tested Tools (WSL Ubuntu + Conda Environments)

**✅ samtools 1.17** (system install)
```bash
$ samtools --version
samtools 1.17
Using htslib 1.20
Copyright (C) 2023 Genome Research Ltd.
# Exit code: 0
# Stream: stdout
# Parse: First line
```

**Alternative:** `samtools version` (no dashes)
**❌ NOT supported:** `samtools -v` (error: unrecognized command)

---

**✅ minimap2 2.28-r1209** (system install)
```bash
$ minimap2 --version
2.28-r1209
# Exit code: 0
# Stream: stdout
# Parse: Single line, prepend "minimap2 "
```

**❌ NOT supported:** `minimap2 -v` (error: missing option argument - verbosity flag)

---

**⚠️ bwa 0.7.18-r1243** (system install)
```bash
$ bwa
# Outputs to STDERR:
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.18-r1243-dirty
# Exit code: 1 (non-zero!)
# Stream: stderr
# Parse: Regex r"version:\s*([\w\.-]+)" (case-insensitive)
```

**❌ NO version flag:** `bwa --version` and `bwa -v` both fail

---

**✅ NanoSim 3.2.2** (conda: env_nanosim)
```bash
$ source ~/miniforge3/etc/profile.d/conda.sh
$ conda activate env_nanosim
$ simulator.py --version
NanoSim 3.2.2
# Exit code: 0
# Stream: stdout
# Parse: Full line
```

**Alternative:** `simulator.py -v` (also works)
**Note:** Tool name is `simulator.py`, not `nanosim`

---

**⚠️ pbsim3 3.0.5** (conda: env_pacbio, binary name is `pbsim`)
```bash
$ conda activate env_pacbio
$ pbsim
USAGE: pbsim [options]
# Exit code: 255 (non-zero!)
# Stream: stdout
# Parse: No version in usage! Tool is pbsim3 (v3.0.5) but binary is "pbsim"
```

**❌ NO version flag:** `pbsim --version` fails (unrecognized option)
**Critical:** pbsim3 package provides `pbsim` binary, version NOT in output!
**Solution:** Return "pbsim3 (version from conda metadata)" or detect via `conda list`

---

**✅ ccs 6.4.0** (conda: env_pacbio, PacBio HiFi tool)
```bash
$ conda activate env_pacbio
$ ccs --version
ccs 6.4.0 (commit v6.4.0)

Using:
  unanimity : 6.4.0 (commit v6.4.0)
  pbbam     : 2.1.0 (commit v2.0.0-26-g05a8314)
# Exit code: 0
# Stream: stdout
# Parse: First line
```

**❌ NOT supported:** `ccs -v` (crashes with CommandLineParserException)

---

**✅ reseq 1.1** (conda: env_wessim)
```bash
$ conda activate env_wessim
$ reseq --version
ReSeq version 1.1
# Exit code: 0
# Stream: stdout
# Parse: Full line
```

**Alternative:** `reseq --help` (contains "Version: 1.1")

---

**⚠️ faToTwoBit** (conda: env_wessim, UCSC tool)
```bash
$ conda activate env_wessim
$ faToTwoBit
faToTwoBit - Convert DNA from fasta to 2bit format
usage:
   faToTwoBit in.fa [in2.fa in3.fa ...] out.2bit
# Exit code: 255 (non-zero!)
# Stream: stdout
# Parse: No version info in output!
```

**❌ NO version flag or info**
**Solution:** Return "faToTwoBit (version unknown)"

---

**⚠️ pblat** (conda: env_wessim, parallel BLAT)
```bash
$ conda activate env_wessim
$ pblat
pblat - BLAT with parallel supports v. 36x2 fast sequence search command line tool
# Exit code: 0
# Stream: stdout
# Parse: Extract "v. 36x2" from first line
```

**❌ NOT supported:** `pblat --version` (error: not a valid option)
**Solution:** Run with no args, regex parse version from usage line

### 3.2 Version Detection Pattern Matrix

| Pattern | Tools | Command | Exit Code | Stream | Parse Strategy | Status |
|---------|-------|---------|-----------|--------|----------------|--------|
| **Standard --version** | samtools, minimap2, NanoSim, ccs, reseq | `--version` | 0 | stdout | First line or full | ✅ Tested |
| **Usage parsing** | pblat | *(no args)* | 0 | stdout | Regex `v\. ([\w]+)` from usage | ✅ Tested |
| **Stderr parsing** | bwa | *(no args)* | 1 ⚠️ | stderr | Regex `version:\s*([\w\.-]+)` | ✅ Tested |
| **No version info** | faToTwoBit | *(no args)* | 255 | stdout | Return "(version unknown)" | ✅ Tested |
| **Special: pbsim3** | pbsim3 | *(binary: pbsim)* | 255 | stdout | **NO VERSION IN OUTPUT** | ⚠️ Problem |

**Critical Findings:**
1. **NanoSim:** Binary is `simulator.py`, not `nanosim`
2. **pbsim3:** Conda package is `pbsim3` but binary is `pbsim` - NO version in output!
3. **bwa:** Non-zero exit code (1) - must use `check=False`
4. **faToTwoBit:** Exit 255, NO version info - return generic string
5. **pblat:** Must parse version from usage line ("v. 36x2")

### 3.3 Error Handling Protocol

**Principle:** Version checking is **informational**, not **critical**.

```python
# CORRECT: Never crash pipeline
try:
    result = subprocess.run(cmd, capture_output=True, timeout=5, check=False)
    version = parse_func(result.stdout, result.stderr)
    return version if version != "N/A" else "N/A"
except Exception as e:
    logger.warning(f"Failed to get version for {tool}: {e}")
    return "N/A"

# WRONG: Don't do this
try:
    result = subprocess.run(cmd, check=True)  # ❌ Will crash on bwa!
    return result.stdout
except:
    raise  # ❌ Crashes pipeline
```

---

## 4. Implementation Design

### 4.1 Architecture (SOLID Principles)

**Single Responsibility Principle:**
- `get_tool_version()`: One function, one job (get version)
- `capture_tool_versions()`: Batch operations
- `write_metadata_file()`: Output formatting

**Open/Closed Principle:**
- Tool parsers easily added via `tool_map` dictionary
- No modification of core logic needed for new tools

**Dependency Inversion:**
- Functions accept tool commands from config (not hardcoded)
- Uses existing `build_tool_command()` utility

### 4.2 Module Structure (DRY)

```
muc_one_up/read_simulator/utils/
├── tool_version.py          # Version detection (NEW)
│   ├── get_tool_version()
│   ├── capture_tool_versions()
│   └── log_tool_versions()
└── metadata_writer.py       # Output formatting (NEW)
    └── write_metadata_file()
```

**Reuses existing:**
- `command_utils.py::build_tool_command()` (conda/mamba support)
- `command_utils.py::get_tool_executable()` (extract tool name)

### 4.3 Core Implementation

**File: `muc_one_up/read_simulator/utils/tool_version.py`**

```python
#!/usr/bin/env python3
"""Tool version detection for read simulation pipelines."""

import logging
import re
import subprocess
from typing import Callable

from ..command_utils import build_tool_command

logger = logging.getLogger(__name__)


def get_tool_version(tool_cmd: str, tool_name: str) -> str:
    """
    Get version string for a bioinformatics tool.

    Args:
        tool_cmd: Command from config (e.g., "samtools" or "mamba run -n env samtools")
        tool_name: Canonical name for lookup ("samtools", "bwa", etc.)

    Returns:
        Version string like "samtools 1.17" or "N/A" on failure

    Design:
        - KISS: Simple string return, not complex dict
        - DRY: Reuses build_tool_command() for conda/mamba
        - SOLID: Tool-specific parsers via strategy pattern (tool_map)
        - Never crashes (graceful degradation to "N/A")
    """

    # Tool-specific parser functions
    def parse_samtools(stdout: str, stderr: str) -> str:
        """Parse: samtools 1.17\\nUsing htslib..."""
        for line in stdout.splitlines():
            if line.startswith("samtools"):
                return line.strip()
        return "N/A"

    def parse_minimap2(stdout: str, stderr: str) -> str:
        """Parse: 2.28-r1209"""
        line = stdout.strip()
        if line:
            return f"minimap2 {line}"
        return "N/A"

    def parse_bwa(stdout: str, stderr: str) -> str:
        """Parse: Version: 0.7.18-r1243 from stderr"""
        for line in stderr.splitlines():
            if "version:" in line.lower():
                match = re.search(r"version:\s*([\w\.-]+)", line, re.IGNORECASE)
                if match:
                    return f"bwa {match.group(1)}"
        return "N/A"

    def parse_nanosim(stdout: str, stderr: str) -> str:
        """Parse: NanoSim 3.2.2 from simulator.py --version"""
        # Tool binary is "simulator.py", not "nanosim"
        for line in stdout.splitlines():
            line_clean = line.strip()
            if line_clean and "nanosim" in line_clean.lower():
                return line_clean  # Returns "NanoSim 3.2.2"
        return "N/A"

    def parse_pbsim3(stdout: str, stderr: str) -> str:
        """
        Parse pbsim3 version.

        CRITICAL: pbsim3 conda package provides "pbsim" binary with NO version in output!
        Fallback: Return generic string since version not retrievable from tool itself.
        """
        # Check if this is pbsim output (usage message)
        if "pbsim" in stdout.lower() or "USAGE: pbsim" in stdout:
            # pbsim3 doesn't output version - return generic identifier
            return "pbsim3 (installed)"
        return "N/A"

    def parse_ccs(stdout: str, stderr: str) -> str:
        """Parse: ccs 6.4.0 (commit v6.4.0)"""
        for line in stdout.splitlines():
            if line.strip().startswith("ccs"):
                return line.strip()  # Returns "ccs 6.4.0 (commit v6.4.0)"
        return "N/A"

    def parse_reseq(stdout: str, stderr: str) -> str:
        """Parse: ReSeq version 1.1"""
        for line in stdout.splitlines():
            if "reseq" in line.lower() and "version" in line.lower():
                return line.strip()  # Returns "ReSeq version 1.1"
        return "N/A"

    def parse_fatotwobit(stdout: str, stderr: str) -> str:
        """
        Parse faToTwoBit version.

        UCSC tool with NO version flag or version in usage output.
        Return generic string to indicate tool is present.
        """
        if "fatotwobit" in stdout.lower() or "fatotwobit" in stderr.lower():
            return "faToTwoBit (version unknown)"
        return "N/A"

    def parse_pblat(stdout: str, stderr: str) -> str:
        """Parse: pblat - BLAT with parallel supports v. 36x2"""
        for line in stdout.splitlines():
            if "pblat" in line.lower():
                # Extract version using regex: "v. 36x2"
                match = re.search(r"v\.\s*([\w]+)", line)
                if match:
                    return f"pblat v.{match.group(1)}"  # Returns "pblat v.36x2"
                return line.strip()  # Fallback: return full line
        return "N/A"

    # Tool command configuration (Strategy Pattern)
    # NOTE: Tool names must match config.json "tools" keys
    tool_map: dict[str, dict] = {
        "samtools": {"args": ["--version"], "parse": parse_samtools},
        "minimap2": {"args": ["--version"], "parse": parse_minimap2},
        "nanosim": {"args": ["--version"], "parse": parse_nanosim},  # Binary: simulator.py
        "pbsim3": {"args": [], "parse": parse_pbsim3},  # Binary: pbsim, NO version output!
        "ccs": {"args": ["--version"], "parse": parse_ccs},
        "reseq": {"args": ["--version"], "parse": parse_reseq},  # --version works, not just --help
        "bwa": {"args": [], "parse": parse_bwa},  # No flag, stderr output
        "faToTwoBit": {"args": [], "parse": parse_fatotwobit},  # No version info available
        "pblat": {"args": [], "parse": parse_pblat},  # No --version flag, parse usage
    }

    # Lookup tool configuration
    if tool_name not in tool_map:
        logger.warning(f"No version logic for {tool_name}. Returning 'N/A'.")
        return "N/A"

    config = tool_map[tool_name]
    parse_func: Callable[[str, str], str] = config["parse"]
    args: list[str] = config["args"]

    try:
        # Build command (handles conda/mamba from config)
        cmd = build_tool_command(tool_cmd, *args)
        logger.debug(f"Getting version: {' '.join(cmd)}")

        # Run with check=False (tools may return non-zero)
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,  # CRITICAL: Don't raise on non-zero exit
            timeout=5,
        )

        # Parse version from stdout/stderr
        version = parse_func(result.stdout, result.stderr)

        if version == "N/A":
            logger.warning(
                f"Could not parse version for {tool_name}. "
                f"stdout: {result.stdout[:100]}, stderr: {result.stderr[:100]}"
            )
        else:
            logger.debug(f"Got version for {tool_name}: {version}")

        return version

    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout getting version for {tool_name} (5s)")
        return "N/A"
    except Exception as e:
        logger.warning(f"Failed to get version for {tool_name}: {e}")
        return "N/A"


def capture_tool_versions(tools_config: dict[str, str]) -> dict[str, str]:
    """
    Capture versions for all configured tools.

    Args:
        tools_config: {"samtools": "samtools", "reseq": "mamba run -n env reseq"}

    Returns:
        {"samtools": "samtools 1.17", "reseq": "reseq version 1.1"}
    """
    return {
        name: get_tool_version(cmd, name)
        for name, cmd in tools_config.items()
    }


def log_tool_versions(versions: dict[str, str]) -> None:
    """Log versions in formatted table."""
    if not versions:
        return

    # Calculate column widths
    max_tool = max(len(t) for t in versions.keys())
    max_ver = max(len(v) for v in versions.values())
    tw, vw = max(max_tool, 12), max(max_ver, 15)

    # Print table
    logger.info("Tool Versions:")
    logger.info(f"┌{'─' * (tw + 2)}┬{'─' * (vw + 2)}┐")
    logger.info(f"│ {'Tool':<{tw}} │ {'Version':<{vw}} │")
    logger.info(f"├{'─' * (tw + 2)}┼{'─' * (vw + 2)}┤")
    for tool, ver in sorted(versions.items()):
        logger.info(f"│ {tool:<{tw}} │ {ver:<{vw}} │")
    logger.info(f"└{'─' * (tw + 2)}┴{'─' * (vw + 2)}┘")
```

**File: `muc_one_up/read_simulator/utils/metadata_writer.py`**

```python
#!/usr/bin/env python3
"""Metadata file generation for simulation provenance."""

import logging
from datetime import datetime
from pathlib import Path

from ...version import __version__
from .tool_version import capture_tool_versions, log_tool_versions

logger = logging.getLogger(__name__)


def write_metadata_file(
    output_dir: str,
    output_base: str,
    config: dict,
    start_time: datetime,
    end_time: datetime,
    platform: str,
) -> str:
    """
    Write TSV metadata file with provenance information.

    Args:
        output_dir: Output directory
        output_base: Base name for files
        config: Pipeline configuration
        start_time: Pipeline start timestamp
        end_time: Pipeline end timestamp
        platform: "Illumina", "ONT", or "PacBio"

    Returns:
        Path to created metadata file

    Output Format (TSV):
        Parameter	Value
        Tool	MucOneUp
        Version	0.15.0
        Platform	Illumina
        Run_start_time	2025-10-24T12:34:56
        Run_duration_seconds	4216.3
        Seed	42
        Coverage	30
        tool.samtools_version	samtools 1.17
        tool.bwa_version	bwa 0.7.18-r1243

    Design:
        - TSV format (human + machine readable)
        - Separate from simulation_stats.json (separation of concerns)
        - Tool versions namespaced with "tool." prefix
        - Follows VariantCentrifuge proven pattern
    """
    metadata_path = Path(output_dir) / f"{output_base}_metadata.tsv"
    logger.info(f"Writing metadata: {metadata_path}")

    # Capture tool versions
    tool_versions = capture_tool_versions(config.get("tools", {}))
    log_tool_versions(tool_versions)

    # Calculate duration
    duration = (end_time - start_time).total_seconds()

    # Write TSV
    with open(metadata_path, "w", encoding="utf-8") as f:
        f.write("Parameter\tValue\n")

        # Core metadata
        f.write(f"Tool\tMucOneUp\n")
        f.write(f"Version\t{__version__}\n")
        f.write(f"Platform\t{platform}\n")
        f.write(f"Run_start_time\t{start_time.isoformat()}\n")
        f.write(f"Run_end_time\t{end_time.isoformat()}\n")
        f.write(f"Run_duration_seconds\t{duration:.1f}\n")

        # Simulation parameters
        f.write(f"Output_base\t{output_base}\n")
        f.write(f"Output_dir\t{output_dir}\n")

        if "seed" in config:
            f.write(f"Seed\t{config['seed']}\n")

        # Platform-specific parameters
        if platform == "Illumina":
            read_cfg = config.get("read_simulation", {})
            f.write(f"Coverage\t{read_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Fragment_size\t{read_cfg.get('fragment_size_mean', 'N/A')}\n")
        elif platform == "ONT":
            ont_cfg = config.get("nanosim_params", {})
            f.write(f"Coverage\t{ont_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Min_read_length\t{ont_cfg.get('min_read_length', 'N/A')}\n")
            f.write(f"Max_read_length\t{ont_cfg.get('max_read_length', 'N/A')}\n")
        elif platform == "PacBio":
            pb_cfg = config.get("pacbio_params", {})
            f.write(f"Coverage\t{pb_cfg.get('coverage', 'N/A')}\n")
            f.write(f"Pass_num\t{pb_cfg.get('pass_num', 'N/A')}\n")

        # Tool versions (namespaced with "tool." prefix)
        for tool_name, version in sorted(tool_versions.items()):
            f.write(f"tool.{tool_name}_version\t{version}\n")

    logger.info(f"Metadata written: {metadata_path}")
    return str(metadata_path)
```

### 4.4 Pipeline Integration

**Illumina Pipeline: `muc_one_up/read_simulator/pipeline.py`**

```python
from datetime import datetime
from .utils.metadata_writer import write_metadata_file

def run_illumina_pipeline(config, ...):
    start_time = datetime.now()

    # ... existing pipeline code ...

    end_time = datetime.now()

    # Write metadata file (LAST STEP)
    write_metadata_file(
        output_dir=output_dir,
        output_base=output_base,
        config=config,
        start_time=start_time,
        end_time=end_time,
        platform="Illumina",
    )
```

**Similarly for ONT and PacBio pipelines.**

---

## 5. Testing Strategy

### 5.1 Unit Tests

**File: `tests/read_simulator/test_tool_version.py`**

```python
import pytest
from unittest.mock import patch, MagicMock
from muc_one_up.read_simulator.utils.tool_version import get_tool_version


class TestGetToolVersion:
    """Test tool version detection."""

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    def test_samtools_standard_flag(self, mock_run):
        """Test samtools with --version flag."""
        mock_run.return_value = MagicMock(
            stdout="samtools 1.17\nUsing htslib 1.20\n",
            stderr="",
            returncode=0,
        )

        version = get_tool_version("samtools", "samtools")
        assert version == "samtools 1.17"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    def test_bwa_no_flag_stderr(self, mock_run):
        """Test bwa parsing from stderr (no --version flag)."""
        mock_run.return_value = MagicMock(
            stdout="",
            stderr="Program: bwa\nVersion: 0.7.18-r1243\n",
            returncode=1,  # Non-zero!
        )

        version = get_tool_version("bwa", "bwa")
        assert version == "bwa 0.7.18-r1243"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    def test_graceful_failure(self, mock_run):
        """Test graceful degradation on failure."""
        mock_run.side_effect = FileNotFoundError("Tool not found")

        version = get_tool_version("nonexistent", "nonexistent")
        assert version == "N/A"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    def test_timeout_handling(self, mock_run):
        """Test timeout returns N/A."""
        import subprocess
        mock_run.side_effect = subprocess.TimeoutExpired(cmd=["tool"], timeout=5)

        version = get_tool_version("samtools", "samtools")
        assert version == "N/A"

    @patch("muc_one_up.read_simulator.utils.tool_version.subprocess.run")
    def test_unparseable_output(self, mock_run):
        """Test unparseable output returns N/A."""
        mock_run.return_value = MagicMock(
            stdout="random garbage\n",
            stderr="",
            returncode=0,
        )

        version = get_tool_version("samtools", "samtools")
        assert version == "N/A"
```

### 5.2 Integration Tests

**Run with real tools when available:**

```python
@pytest.mark.integration
def test_real_samtools_version():
    """Test with actual samtools installation."""
    version = get_tool_version("samtools", "samtools")
    assert version.startswith("samtools")
    assert version != "N/A"
```

---

## 6. Implementation Checklist

### Phase 1: Core Infrastructure (3-4 hours)
- [ ] Create `muc_one_up/read_simulator/utils/tool_version.py`
  - [ ] Implement `get_tool_version()` with 9 tool parsers
  - [ ] Implement `capture_tool_versions()`
  - [ ] Implement `log_tool_versions()`
- [ ] Create `muc_one_up/read_simulator/utils/metadata_writer.py`
  - [ ] Implement `write_metadata_file()`
- [ ] Write unit tests
  - [ ] Test all parsers with mocked subprocess
  - [ ] Test graceful failure scenarios
  - [ ] Test timeout handling

### Phase 2: Pipeline Integration (2-3 hours)
- [ ] Integrate into `pipeline.py` (Illumina)
  - [ ] Add datetime tracking
  - [ ] Call `write_metadata_file()` at end
- [ ] Integrate into `ont_pipeline.py` (ONT)
- [ ] Integrate into `pacbio_pipeline.py` (PacBio)
- [ ] Test end-to-end with each pipeline

### Phase 3: Documentation (1-2 hours)
- [ ] Update `CLAUDE.md` with version tracking info
- [ ] Add example metadata TSV to `data/examples/`
- [ ] Update README with reproducibility section

### Phase 4: CI/CD (1 hour)
- [ ] Add CI check for metadata file generation
- [ ] Update release notes for v0.16.0

**Total Effort:** 7-10 hours

---

## 7. Design Principles Compliance

### KISS (Keep It Simple, Stupid) ✓
- Simple string return, not complex dict
- One function per responsibility
- Minimal abstraction

### DRY (Don't Repeat Yourself) ✓
- Reuses `build_tool_command()` (no code duplication)
- Tool parsers via strategy pattern (no if/else chains)
- Single metadata writer used by all platforms

### SOLID ✓
- **S**ingle Responsibility: Each function has one job
- **O**pen/Closed: Add tools via `tool_map`, don't modify core
- **L**iskov Substitution: N/A (no inheritance)
- **I**nterface Segregation: Functions accept minimal parameters
- **D**ependency Inversion: Depends on abstractions (config), not concrete tools

### Modularity ✓
- Separate modules: `tool_version.py`, `metadata_writer.py`
- Clear interfaces between modules
- Easy to test in isolation

---

## 8. Example Output

### Console Log
```
2025-10-24 12:34:56 - INFO - Starting Illumina pipeline
2025-10-24 12:34:57 - DEBUG - Getting version: samtools --version
2025-10-24 12:34:57 - DEBUG - Got version for samtools: samtools 1.17
2025-10-24 12:34:57 - DEBUG - Getting version: bwa
2025-10-24 12:34:57 - DEBUG - Got version for bwa: bwa 0.7.18-r1243
2025-10-24 12:34:58 - INFO - Tool Versions:
2025-10-24 12:34:58 - INFO - ┌──────────────┬───────────────────┐
2025-10-24 12:34:58 - INFO - │ Tool         │ Version           │
2025-10-24 12:34:58 - INFO - ├──────────────┼───────────────────┤
2025-10-24 12:34:58 - INFO - │ bwa          │ bwa 0.7.18-r1243  │
2025-10-24 12:34:58 - INFO - │ faToTwoBit   │ N/A               │
2025-10-24 12:34:58 - INFO - │ pblat        │ N/A               │
2025-10-24 12:34:58 - INFO - │ reseq        │ N/A               │
2025-10-24 12:34:58 - INFO - │ samtools     │ samtools 1.17     │
2025-10-24 12:34:58 - INFO - └──────────────┴───────────────────┘
... pipeline runs ...
2025-10-24 13:45:12 - INFO - Writing metadata: output/sample_001_metadata.tsv
```

### Metadata File: `sample_001_metadata.tsv`
```tsv
Parameter	Value
Tool	MucOneUp
Version	0.15.0
Platform	Illumina
Run_start_time	2025-10-24T12:34:56.123456
Run_end_time	2025-10-24T13:45:12.789012
Run_duration_seconds	4216.7
Output_base	sample_001
Output_dir	output
Seed	42
Coverage	30
Fragment_size	350
tool.bwa_version	bwa 0.7.18-r1243
tool.faToTwoBit_version	N/A
tool.pblat_version	N/A
tool.reseq_version	N/A
tool.samtools_version	samtools 1.17
```

---

## 9. Migration Notes

### Backward Compatibility
- **Zero breaking changes:** Existing pipelines work unchanged
- **Additive only:** New metadata file created, no modifications to existing outputs
- **Optional:** If version check fails, pipeline continues (logs warning)

### Deprecation Path
None needed - this is a new feature.

---

## 10. References

### Academic Literature
1. **Wratten et al. 2021** - Reproducible workflows with bioinformatics managers. *Nature Methods* 18(10):1161-1173.
2. **Kanwal et al. 2017** - Investigating reproducibility and provenance. *BMC Bioinformatics* 18:337.
3. **PHA4GE 2024** - Public Health Pipeline Best Practices. GitHub: pha4ge/public-health-pipeline-best-practices

### Production Reference
- **VariantCentrifuge** - Clinical variant analysis pipeline
  - Source: `/mnt/c/development/scholl-lab/variantcentrifuge/variantcentrifuge/utils.py:169-248`
  - Pattern: Simple string return, tool-specific parsers, TSV metadata

### Tools Documentation
- samtools: http://www.htslib.org/doc/samtools.html
- minimap2: https://github.com/lh3/minimap2
- bwa: http://bio-bwa.sourceforge.net/

---

## 11. Decision Record

**Date:** 2025-10-24
**Decision:** Implement VariantCentrifuge pattern for v0.16.0
**Rationale:**
- Battle-tested in clinical production
- KISS/DRY/SOLID compliant
- Zero breaking changes
- Low implementation effort (7-10 hours)

**Approval:** ⏳ Pending team review

---

**End of Document**

*Last Updated: 2025-10-24*
*Author: Bioinformatics Analysis*
*Status: Ready for Implementation*
