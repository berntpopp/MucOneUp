# Issue #1: Add PacBio Read Simulation Support

## CORRECTION NOTICE
**Original plan incorrectly stated ONT was "partially implemented". ONT is FULLY functional.**

Current CLI architecture:
```bash
muconeup reads illumina ...  # ✅ Works
muconeup reads ont ...       # ✅ Works
muconeup reads pacbio ...    # ❌ Missing - THIS IS THE ONLY GAP
```

## Problem
MucOneUp currently supports Illumina (w-Wessim2) and ONT (NanoSim) read simulation. PacBio simulation via pbsim3 is not implemented. Users need PacBio HiFi/CCS support for ultra-accurate long-read sequencing.

## Research Background
- **pbsim3** (NAR Genomics 2022): Simulates PacBio CLR, CCS, HiFi, and ONT reads
- **PacBio CCS**: 10-25kb reads with 99.9% accuracy (comparable to Illumina)
- **Use case**: MUC1 VNTR phasing with long, accurate reads

## Current CLI Architecture (DO NOT CHANGE)
```
muc_one_up/cli/click_main.py follows Unix Philosophy:
  @cli.group() reads                    # Read simulation group
    @reads.command() illumina           # Illumina-specific command
    @reads.command() ont                # ONT-specific command
    @reads.command() pacbio  ← ADD THIS # PacBio-specific command (NEW)
```

**This architecture is CORRECT. Do NOT add --read-type flag to simulate command.**

## Implementation

### 1. PacBio Wrapper (NEW FILE)
```python
# muc_one_up/read_simulator/wrappers/pbsim_wrapper.py

def run_pbsim_simulation(
    pbsim_cmd: str,
    reference_fasta: str,
    output_prefix: str,
    coverage: float,
    read_type: str = "SEQUEL",  # CLR, CCS, SEQUEL, SEQUEL2, RSII
    threads: int = 4,
    seed: int | None = None,  # For reproducibility
    min_read_length: int | None = None,
    max_read_length: int | None = None,
    timeout: int = 3600,
) -> str:
    """Run pbsim3 simulation."""
    cmd_list = build_tool_command(
        pbsim_cmd,
        "--strategy", "wgs",
        "--method", "qshmm",
        "--qshmm", read_type,
        "--genome", reference_fasta,
        "--prefix", output_prefix,
        "--depth", coverage,
        "--threads", threads,
    )

    if seed is not None:
        cmd_list.extend(["--seed", str(seed)])

    # Run command using existing run_command() utility
    # Return path to generated FASTQ
```

### 2. PacBio Pipeline (NEW FILE)
```python
# muc_one_up/read_simulator/pacbio_pipeline.py

from .wrappers.pbsim_wrapper import run_pbsim_simulation
from .wrappers.nanosim_wrapper import align_ont_reads_with_minimap2

def simulate_pacbio_reads_pipeline(
    config: dict[str, Any],
    input_fa: str,
    human_reference: str | None = None
) -> str:
    """
    Run PacBio simulation pipeline.

    Steps:
    1. Run pbsim3 to generate PacBio reads
    2. Align with minimap2 (same as ONT)
    3. Create indexed BAM

    Returns:
        Path to BAM file
    """
    # Extract config
    pbsim_params = config.get("pbsim_params", {})
    # ... similar structure to ont_pipeline.py

    # 1. Run pbsim3
    fastq_file = run_pbsim_simulation(...)

    # 2. Align with minimap2 (reuse ONT alignment)
    align_ont_reads_with_minimap2(...)  # Works for PacBio too

    return output_bam
```

### 3. CLI Integration (MODIFY EXISTING FILE)
```python
# muc_one_up/cli/click_main.py

@reads.command()
@click.argument("input_fastas", nargs=-1, required=True, ...)
@click.option("--coverage", type=int, default=30, ...)
@click.option("--read-type",
              type=click.Choice(["CCS", "CLR", "SEQUEL"]),
              default="CCS",
              help="PacBio read type")
@click.pass_context
def pacbio(ctx, input_fastas, out_dir, coverage, read_type):
    """Simulate PacBio long reads from one or more FASTA files.

    Examples:
      muconeup --config X reads pacbio sample.fa --read-type CCS
      muconeup --config X reads pacbio *.simulated.fa
    """
    from ..read_simulation import simulate_reads as simulate_reads_pipeline

    config_path = Path(ctx.obj["config_path"])
    with config_path.open() as f:
        config = json.load(f)

    # Configure PacBio simulation
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = "pacbio"  # NEW
    config["read_simulation"]["coverage"] = coverage

    if "pbsim_params" not in config:
        config["pbsim_params"] = {}
    config["pbsim_params"]["read_type"] = read_type

    # Process files (same pattern as illumina/ont commands)
    for input_fasta in input_fastas:
        simulate_reads_pipeline(config, input_fasta)
```

### 4. Update Dispatcher (MODIFY EXISTING FILE)
```python
# muc_one_up/read_simulation.py

def simulate_reads(config: dict[str, Any], input_fa: str) -> str:
    """Dispatch to appropriate simulator based on config."""
    simulator = config.get("read_simulation", {}).get("simulator", "illumina")

    if simulator.lower() == "ont":
        from muc_one_up.read_simulator.ont_pipeline import simulate_ont_reads_pipeline
        human_reference = config.get("read_simulation", {}).get("human_reference")
        return simulate_ont_reads_pipeline(config, input_fa, human_reference)

    elif simulator.lower() == "pacbio":  # ← ADD THIS BRANCH
        from muc_one_up.read_simulator.pacbio_pipeline import simulate_pacbio_reads_pipeline
        human_reference = config.get("read_simulation", {}).get("human_reference")
        return simulate_pacbio_reads_pipeline(config, input_fa, human_reference)

    else:  # Default: Illumina
        from muc_one_up.read_simulator.pipeline import simulate_reads_pipeline
        return simulate_reads_pipeline(config, input_fa)
```

### 5. Config Schema Extension
```json
{
  "pbsim_params": {
    "coverage": 30,
    "read_type": "CCS",       // or CLR, SEQUEL, SEQUEL2, RSII
    "num_threads": 4,
    "min_read_length": 1000,
    "max_read_length": 50000
  }
}
```

### 6. Conda Environment (NEW FILE)
```yaml
# conda/env_pbsim.yml
name: env_pbsim
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.12
  - pbsim3>=3.0.2
  - minimap2
  - samtools
```

## Testing

### Unit Tests
```python
# tests/read_simulator/test_pbsim_wrapper.py
def test_pbsim_command_with_seed(mocker):
    """Verify pbsim3 command includes --seed."""
    mock_run = mocker.patch("muc_one_up.read_simulator.utils.run_command")
    mock_run.return_value = ("", "", 0)

    run_pbsim_simulation(
        pbsim_cmd="pbsim",
        reference_fasta="test.fa",
        output_prefix="out",
        coverage=30,
        seed=42
    )

    cmd = mock_run.call_args[0][0]
    assert "--seed" in cmd
    assert "42" in cmd
```

### Integration Test
```python
def test_pacbio_pipeline_integration(tmp_path, config):
    """Test full PacBio pipeline."""
    if not shutil.which("pbsim"):
        pytest.skip("pbsim3 not installed")

    config["pbsim_params"] = {"coverage": 5, "read_type": "CCS"}
    bam = simulate_pacbio_reads_pipeline(config, "test.fa")

    assert Path(bam).exists()
    assert Path(bam + ".bai").exists()
```

## Documentation Updates

### README.md
```markdown
### Read Simulation Platforms

MucOneUp supports three sequencing platforms with separate commands:

**Illumina (w-Wessim2):**
```bash
muconeup --config config.json reads illumina sample.fa
```

**Oxford Nanopore (NanoSim):**
```bash
muconeup --config config.json reads ont sample.fa
```

**PacBio (pbsim3):**
```bash
muconeup --config config.json reads pacbio sample.fa --read-type CCS
```

**Platform Comparison:**
| Platform | Command | Read Length | Accuracy | Best For |
|----------|---------|------------|----------|----------|
| Illumina | `reads illumina` | 150-300bp | 99.9% | High coverage, SNPs |
| ONT | `reads ont` | 10-100kb | 95-99% | Structural variants |
| PacBio CCS | `reads pacbio` | 10-25kb | 99.9% | Long + accurate |
```

## Files Modified/Created

### NEW FILES (3):
- `muc_one_up/read_simulator/wrappers/pbsim_wrapper.py`
- `muc_one_up/read_simulator/pacbio_pipeline.py`
- `conda/env_pbsim.yml`
- `tests/read_simulator/test_pbsim_wrapper.py`
- `tests/read_simulator/test_pacbio_pipeline.py`

### MODIFIED FILES (2):
- `muc_one_up/cli/click_main.py` (add `@reads.command() pacbio`)
- `muc_one_up/read_simulation.py` (add pacbio branch to dispatcher)
- `README.md` (document new command)

### DO NOT MODIFY:
- ❌ `ont_pipeline.py` - already complete
- ❌ CLI architecture - already correct
- ❌ Dispatcher pattern - already sound

## Compliance with Project Style

✅ **Unix Philosophy**: Each command does ONE thing
- `reads illumina` - Illumina reads
- `reads ont` - ONT reads
- `reads pacbio` - PacBio reads

✅ **DRY**: Reuses existing alignment code (`align_ont_reads_with_minimap2`)

✅ **KISS**: Simple command structure, no complex flag combinations

✅ **SOLID**: Single Responsibility - each pipeline handles one platform

✅ **Modular**: New files in appropriate directories, follows existing patterns
