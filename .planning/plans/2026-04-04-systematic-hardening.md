# Systematic Codebase Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Harden MucOneUp codebase across CI config, module structure, type safety, and test coverage — shipping as a single v0.39.0 release.

**Architecture:** Four sequential waves: (1) CI/dev-loop quick wins, (2) reads.py DRY + samtools_wrapper split, (3) domain function signatures narrowed to section TypedDicts, (4) error path tests + global randomness fix + pytest markers. All on one feature branch.

**Tech Stack:** Python 3.10+, Click, pytest, mypy, ruff, uv

---

## File Map

**Wave 1 (modify only):**
- `pyproject.toml` — coverage branch, mypy override, marker cleanup
- `.github/workflows/test.yml` — raise threshold 70→75
- `muc_one_up/cli/commands/simulate.py` — fix comment

**Wave 2 (create + modify):**
- Create: `muc_one_up/read_simulator/wrappers/samtools_core.py`
- Create: `muc_one_up/read_simulator/wrappers/samtools_coverage.py`
- Create: `muc_one_up/read_simulator/wrappers/samtools_convert.py`
- Modify: `muc_one_up/read_simulator/wrappers/samtools_wrapper.py` → re-export barrel
- Modify: `muc_one_up/cli/options.py` — shared read options decorator
- Modify: `muc_one_up/cli/commands/reads.py` — use shared decorator + helper

**Wave 3 (modify only):**
- `muc_one_up/assembly.py` — narrow signature
- `muc_one_up/simulate.py` — narrow signatures
- `muc_one_up/mutate.py` — narrow signatures
- `muc_one_up/type_defs.py` — add composite TypedDict if needed
- `muc_one_up/cli/orchestration.py` — update call sites

**Wave 4 (create + modify):**
- Create: `tests/read_simulator/test_common_utils_errors.py`
- Modify: `muc_one_up/read_simulator/fragment_simulation.py` — thread rng
- Modify: `tests/conftest.py` — requires_tools hook
- Modify: `pyproject.toml` — re-add markers
- Modify: `.github/workflows/test.yml` — raise threshold 75→80

---

### Task 1: Create feature branch and capture help-text snapshots

**Files:**
- No source files modified

- [ ] **Step 1: Create feature branch**

```bash
git checkout -b refactor/v0.39-systematic-hardening
```

- [ ] **Step 2: Capture pre-refactor CLI help snapshots for reads subcommands**

```bash
uv run muconeup reads illumina --help > /tmp/reads_illumina_help_before.txt
uv run muconeup reads ont --help > /tmp/reads_ont_help_before.txt
uv run muconeup reads pacbio --help > /tmp/reads_pacbio_help_before.txt
```

- [ ] **Step 3: Verify baseline is green**

Run: `uv run pytest --tb=short -q && uv run mypy muc_one_up/ && uv run ruff check muc_one_up/ tests/`
Expected: 1240+ passed, mypy 0 errors, ruff clean

---

### Task 2: Wave 1a — Raise CI coverage threshold 70→75

**Files:**
- Modify: `.github/workflows/test.yml:106`

- [ ] **Step 1: Edit CI workflow**

In `.github/workflows/test.yml`, change line 106 from:
```yaml
          pytest --cov=muc_one_up --cov-fail-under=70
```
to:
```yaml
          pytest --cov=muc_one_up --cov-fail-under=75
```

- [ ] **Step 2: Verify tests still pass locally with the new threshold**

Run: `uv run pytest --cov=muc_one_up --cov-fail-under=75 --tb=short -q`
Expected: PASS (current coverage is 87%)

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/test.yml
git commit -m "ci: raise coverage threshold from 70% to 75%"
```

---

### Task 3: Wave 1b — Enable branch coverage

**Files:**
- Modify: `pyproject.toml` (coverage.run section)

- [ ] **Step 1: Add branch coverage to pyproject.toml**

In `pyproject.toml`, change the `[tool.coverage.run]` section from:
```toml
[tool.coverage.run]
source = ["muc_one_up"]
omit = [
    "*/tests/*",
    "*/__pycache__/*",
    "*/site-packages/*",
]
```
to:
```toml
[tool.coverage.run]
source = ["muc_one_up"]
branch = true
omit = [
    "*/tests/*",
    "*/__pycache__/*",
    "*/site-packages/*",
]
```

- [ ] **Step 2: Run tests to verify branch coverage works and check the percentage**

Run: `uv run pytest --tb=short -q 2>&1 | tail -5`
Expected: PASS. Note the new coverage percentage — branch coverage is typically lower than line coverage. If it drops below 75%, adjust the CI threshold in `.github/workflows/test.yml` to match (e.g., 70%) and note in the commit message.

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "ci: enable branch coverage measurement"
```

---

### Task 4: Wave 1c — Add mypy override for rfc8785

**Files:**
- Modify: `pyproject.toml` (mypy overrides section)

- [ ] **Step 1: Add rfc8785 mypy override**

In `pyproject.toml`, add after the existing `[[tool.mypy.overrides]]` blocks (after the `tests.*` block around line 160):

```toml
[[tool.mypy.overrides]]
module = "rfc8785.*"
ignore_missing_imports = true
```

- [ ] **Step 2: Verify mypy passes**

Run: `uv run mypy muc_one_up/`
Expected: `Success: no issues found in 75 source files`

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml
git commit -m "chore: add mypy override for rfc8785 stub"
```

---

### Task 5: Wave 1d+1e — Clean up markers and fix stale comment

**Files:**
- Modify: `pyproject.toml` (pytest markers)
- Modify: `muc_one_up/cli/commands/simulate.py:150`

- [ ] **Step 1: Remove unused pytest markers**

In `pyproject.toml`, change the markers list from:
```toml
markers = [
    "unit: marks tests as unit tests (fast, isolated)",
    "integration: marks tests as integration tests (slower, multiple components)",
    "slow: marks tests as slow (> 1 second)",
    "cli: marks tests as CLI interface tests",
    "bioinformatics: marks tests as bioinformatics-specific (sequence validation, etc.)",
    "requires_tools: marks tests requiring external tools (samtools, bwa, etc.)",
]
```
to:
```toml
markers = [
    "unit: marks tests as unit tests (fast, isolated)",
    "integration: marks tests as integration tests (slower, multiple components)",
    "cli: marks tests as CLI interface tests",
    "bioinformatics: marks tests as bioinformatics-specific (sequence validation, etc.)",
]
```

- [ ] **Step 2: Fix stale comment in simulate.py**

In `muc_one_up/cli/commands/simulate.py`, change line 150 from:
```python
    # IMPORTANT: Disable pipeline options (simulate is PURE)
```
to:
```python
    # simulate command generates haplotypes only — no read simulation or ORF output
```

- [ ] **Step 3: Verify**

Run: `uv run pytest --tb=short -q && uv run ruff check muc_one_up/ tests/`
Expected: PASS, ruff clean

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml muc_one_up/cli/commands/simulate.py
git commit -m "chore: remove unused pytest markers, fix stale simulate comment"
```

---

### Task 6: Wave 2b — Split samtools_wrapper.py into 3 modules

**Files:**
- Create: `muc_one_up/read_simulator/wrappers/samtools_core.py`
- Create: `muc_one_up/read_simulator/wrappers/samtools_coverage.py`
- Create: `muc_one_up/read_simulator/wrappers/samtools_convert.py`
- Modify: `muc_one_up/read_simulator/wrappers/samtools_wrapper.py`

This task is mechanical — move functions, fix imports, keep the barrel re-export. The functions themselves do not change.

- [ ] **Step 1: Create samtools_core.py**

Create `muc_one_up/read_simulator/wrappers/samtools_core.py` containing:
- The module docstring and imports (`logging`, `Path`, `ExternalToolError`, `FileOperationError`, `build_tool_command`, `run_command`, `run_pipeline`)
- `extract_subset_reference()` (original lines 27–68)
- `sort_and_index_bam()` (original lines 356–416)
- `merge_bam_files()` (original lines 932–1018)

Copy each function exactly as-is from `samtools_wrapper.py`. The only change is the import header:

```python
"""Samtools core operations: reference extraction, sort/index, merge."""

import logging
from collections.abc import Sequence
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command, run_pipeline
```

- [ ] **Step 2: Create samtools_coverage.py**

Create `muc_one_up/read_simulator/wrappers/samtools_coverage.py` containing:
- `calculate_vntr_coverage()` (original lines 71–137)
- `calculate_target_coverage()` (original lines 140–210)
- `downsample_bam()` (original lines 213–304)
- `downsample_entire_bam()` (original lines 307–353)

Import header:

```python
"""Samtools coverage and downsampling operations."""

import logging
from pathlib import Path

from ...exceptions import FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command
```

- [ ] **Step 3: Create samtools_convert.py**

Create `muc_one_up/read_simulator/wrappers/samtools_convert.py` containing:
- `convert_sam_to_bam()` (original lines 419–494)
- `convert_bam_to_fastq()` (original lines 497–572)
- `FastqConversionOptions` dataclass (original lines 575–600)
- `convert_bam_to_paired_fastq()` (original lines 603–898)
- `_count_fastq_reads()` helper (original lines 901–929)

Import header:

```python
"""Samtools format conversion operations: SAM/BAM/FASTQ."""

import logging
from dataclasses import dataclass
from pathlib import Path

from ...exceptions import ExternalToolError, FileOperationError
from ..command_utils import build_tool_command
from ..utils import run_command, run_pipeline
```

- [ ] **Step 4: Replace samtools_wrapper.py with re-export barrel**

Replace the entire content of `muc_one_up/read_simulator/wrappers/samtools_wrapper.py` with:

```python
"""Samtools wrapper — re-exports from split modules.

All functions are available from this module for backward compatibility.
New code should import from the specific submodules:
- samtools_core: extract_subset_reference, sort_and_index_bam, merge_bam_files
- samtools_coverage: calculate_vntr_coverage, calculate_target_coverage, downsample_bam, downsample_entire_bam
- samtools_convert: convert_sam_to_bam, convert_bam_to_fastq, convert_bam_to_paired_fastq, FastqConversionOptions
"""

from .samtools_convert import (
    FastqConversionOptions,
    convert_bam_to_fastq,
    convert_bam_to_paired_fastq,
    convert_sam_to_bam,
)
from .samtools_core import (
    extract_subset_reference,
    merge_bam_files,
    sort_and_index_bam,
)
from .samtools_coverage import (
    calculate_target_coverage,
    calculate_vntr_coverage,
    downsample_bam,
    downsample_entire_bam,
)

__all__ = [
    "extract_subset_reference",
    "sort_and_index_bam",
    "merge_bam_files",
    "calculate_vntr_coverage",
    "calculate_target_coverage",
    "downsample_bam",
    "downsample_entire_bam",
    "convert_sam_to_bam",
    "convert_bam_to_fastq",
    "convert_bam_to_paired_fastq",
    "FastqConversionOptions",
]
```

- [ ] **Step 5: Verify all tests pass and no module exceeds 500 LOC**

Run: `uv run pytest --tb=short -q && uv run mypy muc_one_up/ && uv run ruff check muc_one_up/ tests/`
Expected: All tests pass, mypy clean, ruff clean.

Run: `wc -l muc_one_up/read_simulator/wrappers/samtools_*.py`
Expected: `samtools_wrapper.py` ~40 lines, all new modules under 500 LOC.

- [ ] **Step 6: Commit**

```bash
git add muc_one_up/read_simulator/wrappers/samtools_core.py \
       muc_one_up/read_simulator/wrappers/samtools_coverage.py \
       muc_one_up/read_simulator/wrappers/samtools_convert.py \
       muc_one_up/read_simulator/wrappers/samtools_wrapper.py
git commit -m "refactor: split samtools_wrapper.py into 3 focused modules

Extract samtools_core (reference, sort, merge), samtools_coverage
(depth, downsampling), samtools_convert (SAM/BAM/FASTQ conversion).
Original module becomes a re-export barrel for backward compatibility."
```

---

### Task 7: Wave 2a — Deduplicate reads.py with shared options

**Files:**
- Modify: `muc_one_up/cli/options.py`
- Modify: `muc_one_up/cli/commands/reads.py`

- [ ] **Step 1: Add shared_read_options decorator to options.py**

Add at the end of `muc_one_up/cli/options.py`:

```python
import functools
import click


def shared_read_options(func):
    """Apply common Click options shared across all read simulation subcommands.

    Adds: input_fastas argument, --out-dir, --out-base, --coverage, --seed,
    --track-read-source. Decorator stacking order is reversed (bottom-up)
    to preserve help text order.
    """
    @click.option(
        "--track-read-source",
        is_flag=True,
        default=False,
        help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
    )
    @click.option(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (same seed = identical reads).",
    )
    @click.option(
        "--coverage",
        type=int,
        default=None,
        help="Target sequencing coverage (overrides config if provided, defaults to config value or 30x).",
    )
    @click.option(
        "--out-base",
        default=None,
        help="Base name for output files (auto-generated if processing multiple files).",
    )
    @click.option(
        "--out-dir",
        default=".",
        show_default=True,
        type=click.Path(file_okay=False),
        help="Output folder.",
    )
    @click.argument(
        "input_fastas",
        nargs=-1,
        required=True,
        type=click.Path(exists=True, dir_okay=False),
    )
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper
```

- [ ] **Step 2: Add _setup_read_config helper to reads.py**

Add this helper function after the `_run_batch_simulation` function in `muc_one_up/cli/commands/reads.py` (around line 104):

```python
def _setup_read_config(
    config: dict[str, Any],
    simulator_type: str,
    coverage: int | None,
    seed: int | None,
    seed_config_key: str = "read_simulation",
) -> None:
    """Initialize read_simulation config section with shared defaults.

    Mutates config in place. Handles:
    - Creating read_simulation section if absent
    - Setting simulator type
    - Applying coverage with fallback to config value or 30x
    - Setting seed and logging it

    Args:
        config: Full configuration dict (mutated in place)
        simulator_type: One of "illumina", "ont", "pacbio"
        coverage: CLI --coverage value (None if not provided)
        seed: CLI --seed value (None if not provided)
        seed_config_key: Config section to store seed in (default: "read_simulation",
                        use "nanosim_params" for ONT, "pacbio_params" for PacBio)
    """
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = simulator_type

    if coverage is not None:
        config["read_simulation"]["coverage"] = coverage
    elif "coverage" not in config.get("read_simulation", {}):
        config["read_simulation"]["coverage"] = 30

    if seed is not None:
        if seed_config_key not in config:
            config[seed_config_key] = {}
        config[seed_config_key]["seed"] = seed
        logging.info(f"Using random seed: {seed} (results will be reproducible)")
```

- [ ] **Step 3: Refactor illumina command**

Replace the illumina command (lines 122–214) with:

```python
@reads.command()
@click.option(
    "--threads",
    type=int,
    default=8,
    show_default=True,
    help="Number of threads.",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed, track_read_source):
    """Simulate Illumina short reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads illumina file.fa --out-base reads
    - Multiple files: muconeup reads illumina file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads illumina *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads illumina sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads illumina sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads illumina sample.*.simulated.fa

      # Compose with shell (Unix philosophy)
      for f in *.fa; do muconeup --config X reads illumina "$f"; done
    """
    require_config(ctx)
    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "illumina", coverage, seed)
    config["read_simulation"]["threads"] = threads

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_reads", "Illumina", track_read_source
    )
```

Add the import at the top of `reads.py`:
```python
from ..options import shared_read_options
```

- [ ] **Step 4: Refactor ont command**

Replace the ont command (lines 217–308) with:

```python
@reads.command()
@click.option(
    "--min-read-length",
    type=int,
    default=100,
    show_default=True,
    help="Minimum read length.",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed, track_read_source):
    """Simulate Oxford Nanopore long reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads ont file.fa --out-base reads
    - Multiple files: muconeup reads ont file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads ont *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads ont sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads ont sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads ont sample.*.simulated.fa
    """
    require_config(ctx)
    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))
    _setup_read_config(config, "ont", coverage, seed, seed_config_key="nanosim_params")

    if "nanosim_params" not in config:
        config["nanosim_params"] = {}
    config["nanosim_params"]["min_len"] = min_read_length

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_ont_reads", "ONT", track_read_source
    )
```

- [ ] **Step 5: Refactor pacbio command**

Replace the pacbio command (lines 311–490). The pacbio command keeps its unique options inline (--pass-num, --min-passes, --min-rq, --model-type, --model-file, --threads) and uses `@shared_read_options` for the 6 common ones:

```python
@reads.command()
@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--model-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to pbsim3 model file (overrides config if provided).",
)
@click.option(
    "--model-type",
    type=click.Choice(["qshmm", "errhmm"]),
    default=None,
    help="pbsim3 model type (overrides config if provided).",
)
@click.option(
    "--min-rq",
    type=float,
    default=None,
    help="Minimum predicted read quality (RQ) score (0.0-1.0). 0.99=Q20 (standard HiFi), 0.999=Q30.",
)
@click.option(
    "--min-passes",
    type=int,
    default=None,
    help="Minimum passes required for CCS HiFi consensus (>=1, overrides config if provided).",
)
@click.option(
    "--pass-num",
    type=int,
    default=None,
    help="Number of passes per molecule for multi-pass CLR simulation (>=2, overrides config if provided).",
)
@shared_read_options
@click.pass_context
@cli_error_handler
def pacbio(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    pass_num,
    min_passes,
    min_rq,
    model_type,
    model_file,
    threads,
    seed,
    track_read_source,
):
    """Simulate PacBio HiFi reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads pacbio file.fa --model-file X.model --out-base reads
    - Multiple files: muconeup reads pacbio file1.fa file2.fa --model-file X.model
    - Glob pattern: muconeup reads pacbio *.simulated.fa --model-file X.model

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Workflow:
      1. Multi-pass CLR simulation (pbsim3)
      2. HiFi consensus generation (CCS)
      3. Read alignment (minimap2 with map-hifi preset)

    \b
    Examples:
      # Single file with standard HiFi settings (Q20)
      muconeup --config X reads pacbio sample.001.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --out-base my_hifi

      # Multiple files with high-accuracy HiFi (Q30)
      muconeup --config X reads pacbio sample.*.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --min-rq 0.999 --min-passes 5

      # Ultra-deep coverage simulation
      muconeup --config X reads pacbio sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --coverage 100 --pass-num 5

    \b
    Model Files:
      Download from: https://github.com/yukiteruono/pbsim3/tree/master/data
      - QSHMM-SEQUEL.model: Sequel II chemistry
      - QSHMM-RSII.model: RS II chemistry
      - ERRHMM-SEQUEL.model: Alternative error model

    \b
    Quality Control:
      - pass_num >=2 required for multi-pass (>=3 recommended)
      - min_passes controls CCS stringency (higher = better quality, lower yield)
      - min_rq=0.99 is Q20 (standard HiFi threshold)
      - min_rq=0.999 is Q30 (ultra-high accuracy)
    """
    require_config(ctx)
    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    # PacBio uses its own config section for most params
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = "pacbio"

    if "pacbio_params" not in config:
        config["pacbio_params"] = {}

    # Override config with CLI params only if provided
    if model_type is not None:
        config["pacbio_params"]["model_type"] = model_type
    if model_file is not None:
        config["pacbio_params"]["model_file"] = model_file
    if coverage is not None:
        config["pacbio_params"]["coverage"] = coverage
    if pass_num is not None:
        config["pacbio_params"]["pass_num"] = pass_num
    if min_passes is not None:
        config["pacbio_params"]["min_passes"] = min_passes
    if min_rq is not None:
        config["pacbio_params"]["min_rq"] = min_rq
    if threads != 4:
        config["pacbio_params"]["threads"] = threads
    if seed is not None:
        config["pacbio_params"]["seed"] = seed
        logging.info(f"Using random seed: {seed} (results will be reproducible)")

    # Validate required config parameters exist
    required_params = ["model_type", "model_file"]
    missing_params = [p for p in required_params if p not in config["pacbio_params"]]
    if missing_params:
        raise click.ClickException(
            f"Missing required PacBio parameters in config: {', '.join(missing_params)}. "
            f"Either add them to config.json pacbio_params section or provide via CLI options."
        )

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_pacbio_hifi", "PacBio HiFi", track_read_source
    )
```

- [ ] **Step 6: Verify help text is identical to pre-refactor snapshots**

```bash
uv run muconeup reads illumina --help > /tmp/reads_illumina_help_after.txt
uv run muconeup reads ont --help > /tmp/reads_ont_help_after.txt
uv run muconeup reads pacbio --help > /tmp/reads_pacbio_help_after.txt
diff /tmp/reads_illumina_help_before.txt /tmp/reads_illumina_help_after.txt
diff /tmp/reads_ont_help_before.txt /tmp/reads_ont_help_after.txt
diff /tmp/reads_pacbio_help_before.txt /tmp/reads_pacbio_help_after.txt
```

Expected: No differences. If decorator stacking order changed the help text order, adjust the decorator order in `shared_read_options` (Click applies decorators bottom-up, so the last decorator in the stack appears first in help).

- [ ] **Step 7: Run full test suite**

Run: `uv run pytest --tb=short -q && uv run mypy muc_one_up/ && uv run ruff check muc_one_up/ tests/`
Expected: All pass.

- [ ] **Step 8: Commit**

```bash
git add muc_one_up/cli/options.py muc_one_up/cli/commands/reads.py
git commit -m "refactor: deduplicate reads.py with shared options decorator

Extract shared_read_options() decorator for 6 common CLI options
and _setup_read_config() helper for config initialization.
Illumina/ONT/PacBio subcommands now use shared decorator.
CLI help text output unchanged."
```

---

### Task 8: Wave 3a — Narrow assembly.py and single-section functions

**Files:**
- Modify: `muc_one_up/assembly.py`
- Modify: `muc_one_up/type_defs.py`
- Modify: `muc_one_up/simulate.py`
- Modify: `muc_one_up/mutate.py`

- [ ] **Step 1: Check what types exist in type_defs.py**

Run: `grep "class\|TypedDict\|= TypeVar\|TypeAlias" muc_one_up/type_defs.py`

Identify existing types: `AssemblyConstants`, `LengthModelDict`, `ProbabilitiesDict`, `MutationDefinition`, `RepeatsDict` (or whatever the repeats type is called). If `RepeatsDict` doesn't exist, add it.

- [ ] **Step 2: Add RepeatsDict if missing**

If `muc_one_up/type_defs.py` does not already have a type alias for the repeats section, add:

```python
RepeatsDict = dict[str, str]
"""Mapping of repeat symbol to DNA sequence."""
```

- [ ] **Step 3: Narrow assembly.py signature**

Change `muc_one_up/assembly.py` from:

```python
from .type_defs import ConfigDict, DNASequence, RepeatUnit
```
```python
def assemble_sequence(chain: list[RepeatUnit], config: ConfigDict) -> DNASequence:
```

to accept explicit sections:

```python
from .type_defs import AssemblyConstants, DNASequence, RepeatUnit
```
```python
def assemble_sequence(
    chain: list[RepeatUnit],
    repeats: dict[str, str],
    constants: AssemblyConstants,
) -> DNASequence:
```

Update the function body to use `repeats` and `constants` directly instead of extracting from config:

```python
    left_const = constants["left"]
    right_const = constants["right"]

    parts: list[str] = [left_const]
    for unit in chain:
        if unit.symbol not in repeats:
            raise KeyError(f"Repeat symbol '{unit.symbol}' not found in repeats")
        parts.append(repeats[unit.symbol])

    if chain and chain[-1].symbol == "9":
        parts.append(right_const)
    elif chain:
        logger.debug("Chain does not end with '9'; right constant omitted.")

    return "".join(parts)
```

- [ ] **Step 4: Update all callers of assemble_sequence**

Search for callers: `grep -rn "assemble_sequence" muc_one_up/`

Update each call site to pass explicit sections. The main callers are:

In `muc_one_up/simulate.py` — `assemble_haplotype_from_chain()`:
```python
def assemble_haplotype_from_chain(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    typed_chain = [RepeatUnit.from_str(s) for s in chain]
    ref_assembly = config.get("reference_assembly", "hg38")
    return assemble_sequence(typed_chain, config["repeats"], config["constants"][ref_assembly])
```

In `muc_one_up/simulate.py` — `simulate_from_chains()`, in the loop body where `assemble_sequence` is called:
```python
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config["constants"][ref_assembly]
    repeats = config["repeats"]
    # ... in the loop:
    seq = assemble_sequence(working_chain, repeats, constants)
```

In `muc_one_up/mutate.py` — `rebuild_haplotype_sequence()`:
```python
def rebuild_haplotype_sequence(chain: RepeatChain, config: ConfigDict) -> DNASequence:
    typed_chain = [RepeatUnit.from_str(s) for s in chain]
    ref_assembly = config.get("reference_assembly", "hg38")
    return assemble_sequence(typed_chain, config["repeats"], config["constants"][ref_assembly])
```

In `muc_one_up/mutate.py` — `apply_mutations()`, where `assemble_sequence` is called after mutation:
```python
    ref_assembly = config.get("reference_assembly", "hg38")
    constants = config["constants"][ref_assembly]
    repeats_dict = config["repeats"]
    # ... then in the mutation loop, pass repeats_dict and constants
```

- [ ] **Step 5: Narrow validate_allowed_repeats in mutate.py**

Change:
```python
def validate_allowed_repeats(mutation_def: MutationDefinition, config: ConfigDict) -> set[str]:
```
to:
```python
def validate_allowed_repeats(mutation_def: MutationDefinition, valid_symbols: set[str]) -> set[str]:
```

Update the body to use `valid_symbols` instead of `set(config["repeats"].keys())`.

Update the caller in `apply_mutations()`:
```python
    valid_symbols = set(config["repeats"].keys())
    allowed_repeats = validate_allowed_repeats(mutation_def, valid_symbols)
```

- [ ] **Step 6: Verify**

Run: `uv run pytest --tb=short -q && uv run mypy muc_one_up/ && uv run ruff check muc_one_up/ tests/`
Expected: All pass. TypedDicts are structurally compatible with plain dicts, so tests pass without modification.

- [ ] **Step 7: Commit**

```bash
git add muc_one_up/assembly.py muc_one_up/simulate.py muc_one_up/mutate.py muc_one_up/type_defs.py
git commit -m "refactor: narrow assembly/mutate signatures to section types

assemble_sequence() now takes repeats + constants instead of full config.
validate_allowed_repeats() takes valid_symbols set instead of config.
rebuild_haplotype_sequence() and assemble_haplotype_from_chain() unpack
config at call boundary."
```

---

### Task 9: Wave 3b — Narrow simulate.py multi-section functions

**Files:**
- Modify: `muc_one_up/simulate.py`
- Modify: `muc_one_up/cli/orchestration.py` (caller updates)

- [ ] **Step 1: Narrow simulate_single_haplotype**

`simulate_single_haplotype()` accesses `config["probabilities"]`, `config["repeats"]`, `config["constants"]`, and `config.get("reference_assembly")`. Since it needs 3+ sections, pass them explicitly:

```python
def simulate_single_haplotype(
    config: ConfigDict,
    target_length: int,
    min_length: int = 10,
    rng: _random_module.Random | None = None,
) -> HaplotypeResult:
```

This function is deeply intertwined with config (it builds sequences inline, not just via `assemble_sequence`). Keep `config: ConfigDict` here but document which sections it reads. The inline assembly in lines 244–282 directly accesses `config["constants"]`, `config["probabilities"]`, and `config["repeats"]` — narrowing would require passing 4 separate params which is worse than one config.

**Decision:** Leave `simulate_single_haplotype` and `simulate_diploid` accepting `config: ConfigDict`. These are the top-level orchestration functions in simulate.py — they legitimately cross 4 config sections. Document the accessed sections in the docstring instead.

- [ ] **Step 2: Update docstrings to document accessed sections**

In `simulate_diploid()`:
```python
    """Simulate multiple haplotypes with configurable parameters.

    Config sections accessed:
        - config["length_model"]: repeat count distribution
        - config["probabilities"]: state transition weights
        - config["repeats"]: repeat symbol to sequence mapping
        - config["constants"][assembly]: flanking constant sequences
```

In `simulate_single_haplotype()`:
```python
    """Build single haplotype by chaining repeats probabilistically.

    Config sections accessed:
        - config["probabilities"]: state transition weights
        - config["repeats"]: repeat symbol to sequence mapping
        - config["constants"][assembly]: flanking constant sequences
        - config["reference_assembly"]: assembly name (default: "hg38")
```

- [ ] **Step 3: Verify**

Run: `uv run pytest --tb=short -q && uv run mypy muc_one_up/ && uv run ruff check muc_one_up/ tests/`
Expected: All pass.

- [ ] **Step 4: Commit**

```bash
git add muc_one_up/simulate.py
git commit -m "docs: document config section access in simulate.py functions

simulate_diploid and simulate_single_haplotype legitimately cross
4 config sections — document rather than artificially narrow."
```

---

### Task 10: Wave 4a — Error path tests for common_utils.py

**Files:**
- Create: `tests/read_simulator/test_common_utils_errors.py`

- [ ] **Step 1: Write parametrized error path tests**

Create `tests/read_simulator/test_common_utils_errors.py`:

```python
"""Error path tests for run_command and run_pipeline."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from muc_one_up.exceptions import ExternalToolError
from muc_one_up.read_simulator.utils.common_utils import RunResult, run_command


class TestRunCommandErrors:
    """Test run_command error handling."""

    @patch("muc_one_up.read_simulator.utils.common_utils.subprocess.Popen")
    def test_nonzero_exit_raises_external_tool_error(self, mock_popen):
        """run_command raises ExternalToolError on non-zero exit."""
        proc = MagicMock()
        proc.returncode = 1
        proc.stdout = iter([])
        proc.stderr = MagicMock()
        proc.stderr.__iter__ = MagicMock(return_value=iter([]))
        proc.wait.return_value = 1
        mock_popen.return_value.__enter__ = MagicMock(return_value=proc)
        mock_popen.return_value.__exit__ = MagicMock(return_value=False)

        with pytest.raises(ExternalToolError, match="failed with exit code 1"):
            run_command(["false"])

    @patch("muc_one_up.read_simulator.utils.common_utils.subprocess.Popen")
    def test_timeout_raises_external_tool_error(self, mock_popen):
        """run_command raises ExternalToolError on timeout."""
        proc = MagicMock()
        proc.wait.side_effect = subprocess.TimeoutExpired(cmd="test", timeout=1)
        proc.kill.return_value = None
        proc.stdout = iter([])
        proc.stderr = MagicMock()
        proc.stderr.__iter__ = MagicMock(return_value=iter([]))
        mock_popen.return_value.__enter__ = MagicMock(return_value=proc)
        mock_popen.return_value.__exit__ = MagicMock(return_value=False)

        with pytest.raises(ExternalToolError, match="timed out"):
            run_command(["sleep", "100"], timeout=1)
```

Note: The exact mock structure depends on how `run_command` uses subprocess internally. Read the actual implementation during execution and adjust mocks to match the real code paths.

- [ ] **Step 2: Run tests to verify they pass**

Run: `uv run pytest tests/read_simulator/test_common_utils_errors.py -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tests/read_simulator/test_common_utils_errors.py
git commit -m "test: add error path tests for run_command/run_pipeline"
```

---

### Task 11: Wave 4b — Fix global randomness in fragment_simulation.py

**Files:**
- Modify: `muc_one_up/read_simulator/fragment_simulation.py`

- [ ] **Step 1: Thread rng through fragment_simulation functions**

In `muc_one_up/read_simulator/fragment_simulation.py`, the functions `pick_on_match()`, `get_insert_length()`, `pick_fragment()` use global `random`. Thread an `rng` parameter through them:

Change `pick_on_match` (around line 210):
```python
def pick_on_match(
    matches: list[tuple[str, int, int, str]],
    rng: random.Random | None = None,
) -> tuple[str, int, int, str]:
```
Body: `return (rng or random).choice(matches)`

Change `get_insert_length` (around line 226):
```python
def get_insert_length(mu: float, sigma: float, lower: int, rng: random.Random | None = None) -> int:
```
Body: replace `random.gauss` with `(rng or random).gauss`

Change `pick_fragment` (around line 244):
```python
def pick_fragment(
    match: tuple[str, int, int, str], ins: int, bind: float, rng: random.Random | None = None,
) -> tuple[str, int, int, str]:
```
Body: replace `random.randint` with `(rng or random).randint`

Change `simulate_fragments` (around line 297):
```python
def simulate_fragments(
    ...,
    seed: int | None = None,
    fragment_origins_path: str | None = None,
) -> None:
```
Replace the global seed block:
```python
    # Initialize RNG for reproducibility
    rng = random.Random(seed) if seed is not None else None
    if seed is not None:
        logging.info(f"Fragment simulation using random seed: {seed}")
```

Then pass `rng=rng` to all calls to `pick_on_match`, `get_insert_length`, `pick_fragment` within `simulate_fragments`.

- [ ] **Step 2: Update existing tests if needed**

Run: `uv run pytest tests/read_simulator/test_fragment_simulation.py -v`
Expected: All existing tests pass (the `rng=None` default preserves backward compat).

- [ ] **Step 3: Verify no global random usage remains**

Run: `grep -n "random\.\(choice\|choices\|seed\|randint\|gauss\|sample\|shuffle\)" muc_one_up/read_simulator/fragment_simulation.py`
Expected: No matches (all usage now goes through `rng` or the `(rng or random)` pattern).

- [ ] **Step 4: Run full suite**

Run: `uv run pytest --tb=short -q && uv run ruff check muc_one_up/ tests/`
Expected: All pass.

- [ ] **Step 5: Commit**

```bash
git add muc_one_up/read_simulator/fragment_simulation.py
git commit -m "refactor: thread rng through fragment_simulation functions

Replace global random.choice/gauss/randint with explicit rng parameter
in pick_on_match, get_insert_length, pick_fragment, simulate_fragments.
Defaults to global random for backward compatibility."
```

---

### Task 12: Wave 4c — Operationalize pytest markers

**Files:**
- Modify: `tests/conftest.py`
- Modify: `pyproject.toml`

- [ ] **Step 1: Add pytest_collection_modifyitems hook to conftest.py**

Add at the top of `tests/conftest.py` (after existing imports):

```python
import shutil
```

Add at the end of `tests/conftest.py`:

```python
# ============================================================================
# Marker-based tool skipping
# ============================================================================


def pytest_collection_modifyitems(config, items):
    """Skip tests marked with @pytest.mark.requires_tools when tools are absent.

    Usage: @pytest.mark.requires_tools("samtools", "bwa")
    Skips the test if any named binary is not found on PATH.
    """
    for item in items:
        for marker in item.iter_markers("requires_tools"):
            for tool in marker.args:
                if not shutil.which(tool):
                    item.add_marker(
                        pytest.mark.skip(reason=f"requires external tool: {tool}")
                    )
                    break
```

- [ ] **Step 2: Re-add markers to pyproject.toml**

In `pyproject.toml`, update the markers list to:

```toml
markers = [
    "unit: marks tests as unit tests (fast, isolated)",
    "integration: marks tests as integration tests (slower, multiple components)",
    "slow: marks tests as slow (> 1 second)",
    "cli: marks tests as CLI interface tests",
    "bioinformatics: marks tests as bioinformatics-specific (sequence validation, etc.)",
    "requires_tools: marks tests requiring external tools (args: tool names, e.g. 'samtools', 'bwa')",
]
```

- [ ] **Step 3: Identify slow tests**

Run: `uv run pytest --durations=20 --no-header -q 2>&1 | head -25`

Mark any tests consistently > 1s with `@pytest.mark.slow`.

- [ ] **Step 4: Verify**

Run: `uv run pytest --tb=short -q`
Expected: All pass. The hook only skips tests that have the marker; no existing tests are affected yet.

- [ ] **Step 5: Commit**

```bash
git add tests/conftest.py pyproject.toml
git commit -m "test: add requires_tools marker hook and re-add pytest markers

pytest_collection_modifyitems skips tests marked @requires_tools('tool')
when the named binary is absent. Re-adds slow and requires_tools markers."
```

---

### Task 13: Wave 4d — Raise CI threshold to 80%

**Files:**
- Modify: `.github/workflows/test.yml:106`

- [ ] **Step 1: Verify current coverage with branch measurement**

Run: `uv run pytest --cov=muc_one_up --cov-fail-under=80 --tb=short -q 2>&1 | tail -3`
Expected: PASS (coverage should still be above 80% with branch coverage).

If coverage is below 80%, set threshold to the nearest 5% below current coverage and note in commit.

- [ ] **Step 2: Update CI threshold**

In `.github/workflows/test.yml`, change:
```yaml
          pytest --cov=muc_one_up --cov-fail-under=75
```
to:
```yaml
          pytest --cov=muc_one_up --cov-fail-under=80
```

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/test.yml
git commit -m "ci: raise coverage threshold from 75% to 80%"
```

---

### Task 14: Version bump and final verification

**Files:**
- Modify: `pyproject.toml`
- Modify: `muc_one_up/version.py`

- [ ] **Step 1: Run full quality gate**

```bash
uv run ruff check muc_one_up/ tests/
uv run ruff format --check muc_one_up/ tests/
uv run mypy muc_one_up/
uv run pytest --tb=short -q
```

Expected: All pass, all clean.

- [ ] **Step 2: Verify success criteria**

```bash
# No module > 500 LOC
wc -l muc_one_up/read_simulator/wrappers/samtools_*.py

# No function takes ConfigDict when it only reads one section
grep -rn "config: ConfigDict" muc_one_up/assembly.py muc_one_up/mutate.py

# No global random in fragment_simulation
grep -n "random\.\(choice\|gauss\|randint\|seed\)" muc_one_up/read_simulator/fragment_simulation.py
```

Expected: samtools modules all < 500 LOC. assembly.py has 0 ConfigDict hits. mutate.py: only `apply_mutations` and `apply_changes_to_repeat` still take ConfigDict (multi-section). fragment_simulation.py: 0 global random hits.

- [ ] **Step 3: Bump version**

```bash
make bump-minor
```

This should bump from 0.38.0 to 0.39.0 in both `pyproject.toml` and `muc_one_up/version.py`.

- [ ] **Step 4: Commit version bump**

```bash
git add pyproject.toml muc_one_up/version.py
git commit -m "chore: bump version to 0.39.0"
```
