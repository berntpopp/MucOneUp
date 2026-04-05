# Amplicon Simulation Design Spec

**Issue:** [#72 — Add PBSIM3 template mode for amplicon-like long-read simulation](https://github.com/berntpopp/MucOneUp/issues/72)
**Date:** 2026-04-04
**Status:** Draft

## Overview

Add a PacBio amplicon read simulation mode to MucOneUp that produces realistic PCR amplicon-like reads from simulated diploid VNTR references. This enables benchmarking of MUC1 VNTR analysis pipelines (VNTRPipeline, PacMUCI) that work on targeted long-range PCR products rather than WGS data.

The pipeline extracts amplicon regions using primer binding site detection, applies a length-dependent PCR bias model to determine per-allele coverage, and uses PBSIM3's template mode (`--strategy templ`) to generate full-length reads with realistic error profiles.

ONT support is deferred to a follow-up issue.

## Config Contract

The amplicon pipeline reads from two config sections with clear ownership:

- **`amplicon_params`** — Owns amplicon-specific settings only: primer sequences, expected product range, PCR bias configuration. Does not contain any PBSIM3/CCS/minimap2 settings.
- **`pacbio_params`** — Continues to own all PBSIM3, CCS, and minimap2 settings (model_type, model_file, pass_num, min_passes, min_rq, threads, etc.). The amplicon pipeline reads these exactly as the existing PacBio pipeline does.
- **`read_simulation.simulator`** — Must be extended to accept `"amplicon"` in addition to `"illumina"`, `"ont"`, and `"pacbio"`. The config schema in `muc_one_up/config.py` must be updated accordingly.

This means the amplicon CLI command writes to `amplicon_params` for amplicon-specific overrides and to `pacbio_params` for PBSIM3/CCS overrides, following the same pattern as the existing `pacbio` command.

## Pipeline Stages

1. **Haplotype extraction** — Split diploid FASTA into hap1/hap2 using existing `reference_utils.extract_haplotypes()`
2. **Amplicon extraction** — Find primer binding sites in each haplotype, extract the region between primers (inclusive). Hard error if primers not found.
3. **PCR bias computation** — Given amplicon lengths, compute per-allele read counts from total coverage using the exponential decay model.
4. **Template FASTA generation** — Replicate each amplicon N1/N2 times into template FASTAs (one FASTA record per desired read).
5. **PBSIM3 template simulation** — Run `pbsim3 --strategy templ` per allele. Produces one multi-pass CLR read per template record.
6. **CCS consensus** — Run CCS on multi-pass CLR output to produce HiFi reads (existing `ccs_wrapper`).
7. **BAM merge** — Merge per-allele HiFi BAMs (existing `samtools_wrapper.merge_bam_files()`).
8. **Alignment** — Align merged reads to human reference with minimap2 map-hifi preset (optional, skipped if no human reference).

**Haploid/diploid branching:** Before any extraction, the pipeline calls `is_diploid_reference()` (from `reference_utils`) to inspect sequence count. This avoids calling `extract_haplotypes()` on haploid input, which would fail by design. For diploid input (2 sequences), the full pipeline runs. For haploid input (1 sequence), stages 1 and 3 are skipped — full coverage goes to the single allele, amplicon extraction runs directly on the input FASTA.

**Coverage semantics:** `--coverage` targets template molecules before CCS filtering. The pipeline generates N template copies and PBSIM3 produces N multi-pass CLR reads, but CCS filtering (min_rq, min_passes) may reduce the final HiFi read count below N. This is the honest approach — CCS yield is model-dependent and cannot be reliably predicted for oversampling. Users who need a specific final HiFi count should oversample by adjusting `--coverage` upward.

**Read-source tracking:** `--track-read-source` is not supported in amplicon mode (v1). Amplicon reads are simulated from extracted subregions, so template-space coordinates do not map directly to original haplotype coordinates without remapping via primer positions. If `--track-read-source` is passed with `reads amplicon`, the CLI raises a clear error: "Read source tracking is not yet supported for amplicon simulation."

## PCR Bias Model (`pcr_bias.py`)

### Mathematical Foundation

Per-cycle amplification efficiency decays exponentially with amplicon length:

```
E(L) = E_max * exp(-alpha * L)
```

For two competing alleles in the same reaction, the yield ratio after `c` cycles (Suzuki-Giovannoni 1996):

```
ratio = ((1 + E1) / (1 + E2)) ^ c
```

Where:
- `E_max` — Maximum per-cycle efficiency for short amplicons
- `alpha` — Length decay rate (bp^-1), physically motivated by depurination damage and elongation time limits
- `L` — Amplicon length in bp
- `c` — Number of PCR cycles

### Parameters

All parameters are exposed in config and overridable:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `e_max` | float | 0.95 | Max per-cycle efficiency (short amplicons) |
| `alpha` | float | calibrated | Length decay rate (bp^-1) |
| `cycles` | int | 25 | Number of PCR cycles |
| `denaturation_time` | float | 10.0 | Denaturation step duration (seconds) |
| `stochastic` | bool | false | Enable Galton-Watson branching process |

### Presets

- **`"default"`** — KOD HS polymerase, 25 cycles, 98C/10s denaturation. Calibrated so heterozygous alleles differing by >1.5kb produce ~65-70% reads for the shorter allele, matching ~30% allelic dropout observed by Madritsch et al. 2026.
- **`"no_bias"`** — Equal 50/50 split (E_max=1.0, alpha=0). For controlled benchmarking where PCR bias is not the variable under test.

Preset values can be overridden by specifying individual parameters alongside the preset.

### Interface

```python
class PCRBiasModel:
    @classmethod
    def from_preset(cls, name: str, **overrides) -> PCRBiasModel:
        """Load a named preset, with optional parameter overrides."""

    @classmethod
    def from_params(cls, e_max: float, alpha: float, cycles: int,
                    denaturation_time: float = 10.0,
                    stochastic: bool = False) -> PCRBiasModel:
        """Construct from explicit parameters."""

    def compute_coverage_split(
        self, total_coverage: int,
        allele1_length: int, allele2_length: int,
        seed: int | None = None,
    ) -> tuple[int, int]:
        """
        Compute per-allele read counts from total coverage.
        
        In deterministic mode: ratio is computed directly from the model.
        In stochastic mode: Galton-Watson branching process with binomial
        draws per cycle. Seeded for reproducibility.
        
        Returns (reads_allele1, reads_allele2) where sum = total_coverage.
        """
```

### Stochastic Mode

When `stochastic=True`, each PCR cycle is modeled as independent binomial trials:

```python
n1, n2 = N0, N0  # equal starting copies
for cycle in range(cycles):
    n1 += rng.binomial(n1, E1)
    n2 += rng.binomial(n2, E2)
ratio = n1 / (n1 + n2)
```

This produces realistic run-to-run variance. Controlled via `--seed`.

### Calibration

The `alpha` parameter for the `"default"` preset is calibrated against empirical data from Madritsch et al. 2026:
- 30% of heterozygous samples showed complete allelic dropout of the longer allele
- MUC1 VNTR alleles range from ~1.5kb (25 repeats) to ~6kb (100 repeats)
- Target: alleles differing by >1.5kb should produce ~65-70% shorter allele reads

References:
- Suzuki & Giovannoni 1996 — Competitive PCR model
- Booth/Louw et al. 2011 (PMC3148806) — Elongation efficiency framework
- Kebschull & Zador 2015 (PMC4666380) — Stochastic branching process validation
- Madritsch et al. 2026 (Sci Rep 16:762) — MUC1 empirical dropout data
- Lindahl & Nyberg 1972 — Depurination rate constants

## PBSIM3 Template Mode Integration

### New Wrapper Function

Added to `pbsim3_wrapper.py` alongside the existing `run_pbsim3_simulation()`:

```python
def run_pbsim3_template_simulation(
    pbsim3_cmd: str,
    samtools_cmd: str,
    template_fasta: str,
    model_type: str,
    model_file: str,
    output_prefix: str,
    pass_num: int = DEFAULT_PBSIM3_PASS_NUM,
    accuracy_mean: float = DEFAULT_PBSIM3_ACCURACY_MEAN,
    seed: int | None = None,
    timeout: int = DEFAULT_PBSIM3_TIMEOUT,
) -> list[str]:
```

Key differences from WGS mode:
- `--strategy templ` instead of `--strategy wgs`
- `--template` instead of `--genome`
- No `--depth` (one read per FASTA record)
- No `--length-*` parameters (reads span full template length)
- Same `--method`, model, `--pass-num`, accuracy, and seed parameters

### Command Construction

```
pbsim --strategy templ \
  --method qshmm --qshmm /path/to/model \
  --template amplicon_copies.fa \
  --pass-num 3 \
  --accuracy-mean 0.85 \
  --prefix output_prefix \
  --seed 42
```

## Amplicon Extraction

### Shared Primer Utility (`sequence_utils.py`)

Extracted from `PCRSimulator._find_binding_sites()` in `snapshot_validator.py`:

```python
def find_primer_binding_sites(
    template: str,
    primer: str,
    reverse_complement: bool = False,
) -> list[int]:
    """
    Find all binding sites for a primer in template sequence.
    Both template and primer are normalized to uppercase before matching.
    Exact string matching after normalization.
    Returns list of 0-based positions.
    If reverse_complement=True, searches for RC of primer.
    """
```

`PCRSimulator` is refactored to use this shared function. Existing behavior and tests unchanged.

### Amplicon Extractor (`amplicon_extractor.py`)

```python
class AmpliconResult(NamedTuple):
    sequence: str
    length: int
    forward_pos: int
    reverse_pos: int
    fasta_path: str

class AmpliconExtractor:
    def __init__(self, forward_primer: str, reverse_primer: str):
        """
        forward_primer: 5'->3' forward primer sequence
        reverse_primer: 5'->3' reverse primer sequence (will be RC'd for search)
        """

    def extract(self, haplotype_fasta: str, output_path: str) -> AmpliconResult:
        """
        Extract amplicon region from haplotype FASTA.
        
        1. Read sequence from FASTA
        2. Find forward primer binding site (forward strand)
        3. Find reverse primer binding site (search for RC on forward strand)
        4. Extract region from forward primer start to reverse primer end (inclusive)
        5. Write extracted sequence to output FASTA
        
        Raises AmpliconExtractionError if:
        - Forward primer not found
        - Reverse primer not found
        - Multiple binding sites for either primer (ambiguous)
        - Extracted region outside expected_product_range (if configured;
          bounds are inclusive, e.g. [1500, 6000] means 1500 <= length <= 6000)
        """
```

### Error Messages

Clear, actionable error messages:
- "Forward primer 'GGAG...' not found in haplotype_1. Ensure the simulated reference includes sufficient flanking sequence to contain the primer binding site."
- "Multiple forward primer binding sites found at positions [100, 5200]. Amplicon extraction requires unambiguous primer binding."

## Template FASTA Generator (`template_generator.py`)

```python
def generate_template_fasta(
    amplicon_fasta: str,
    num_copies: int,
    output_path: str,
) -> str:
    """
    Replicate amplicon sequence N times into a multi-record FASTA.
    
    Each record named 'amplicon_copy_001', 'amplicon_copy_002', etc.
    PBSIM3 --strategy templ produces one read per record.
    
    Returns path to output FASTA.
    """
```

## Configuration

### Config.json Addition

The amplicon pipeline requires both `amplicon_params` (amplicon-specific) and `pacbio_params` (PBSIM3/CCS settings, shared with the existing PacBio pipeline):

```json
"amplicon_params": {
    "forward_primer": "GGAGAAAAGGAGACTTCGGCTACCCAG",
    "reverse_primer": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    "primer_source": "Wenzel et al. 2018 (PMID: 29520014)",
    "expected_product_range": [1500, 6000],
    "pcr_bias": {
        "preset": "default"
    }
}
```

`expected_product_range` is optional. When omitted, no size validation is performed on extracted amplicons. When present, bounds are inclusive: `[min, max]` means `min <= length <= max`.

PCR bias parameters can be overridden individually:

```json
"pcr_bias": {
    "preset": "default",
    "cycles": 30,
    "stochastic": true
}
```

Or specified without a preset:

```json
"pcr_bias": {
    "e_max": 0.90,
    "alpha": 0.0002,
    "cycles": 20,
    "denaturation_time": 15.0,
    "stochastic": false
}
```

## CLI Command

```bash
muconeup --config config.json reads amplicon sample.simulated.fa \
  --coverage 1000 \
  --model-file /models/QSHMM-SEQUEL.model \
  --stochastic-pcr \
  --pcr-preset default \
  --seed 42
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `INPUT_FASTAS` | argument | required | One or more FASTA files |
| `--coverage` | int | 30 | Total template molecules before CCS (split by PCR bias) |
| `--model-file` | path | from config | PBSIM3 model file |
| `--model-type` | choice | from config | qshmm or errhmm |
| `--pcr-preset` | choice | default | PCR bias preset name |
| `--stochastic-pcr` | flag | off | Enable stochastic PCR bias |
| `--seed` | int | None | Random seed for reproducibility |
| `--out-dir` | path | input dir | Output directory |
| `--out-base` | str | auto | Output base name |
| `--track-read-source` | flag | off | **Not supported in v1** — raises error if used |

Options inherited via `@shared_read_options` decorator: `INPUT_FASTAS`, `--out-dir`, `--out-base`, `--coverage`, `--seed`, `--track-read-source`. Note: `--track-read-source` is accepted by the decorator but rejected at runtime for amplicon mode with a clear error message.

### Haploid Handling

Detected via `is_diploid_reference()` before any extraction attempt:
- **Haploid (1 sequence):** Skip haplotype extraction and PCR bias computation. Full coverage allocated to the single amplicon. Pipeline proceeds directly to amplicon extraction.
- **Diploid (2 sequences):** Full pipeline with haplotype extraction, per-allele amplicon extraction, PCR bias split, and per-allele simulation.

## File Organization

### New Files

```
muc_one_up/read_simulator/
    amplicon_pipeline.py              # Pipeline orchestrator
    pcr_bias.py                       # PCR length bias model + presets
    utils/
        sequence_utils.py             # Shared primer binding site detection
        amplicon_extractor.py         # Primer-based amplicon extraction
        template_generator.py         # Replicate amplicon into template FASTA

tests/read_simulator/
    test_pcr_bias.py
    test_amplicon_extractor.py
    test_template_generator.py
    test_amplicon_pipeline.py
    test_sequence_utils.py
```

### Modified Files

```
muc_one_up/cli/commands/reads.py              # Add amplicon subcommand
muc_one_up/read_simulation.py                 # Add "amplicon" to simulator map
muc_one_up/read_simulator/constants.py        # Add amplicon constants
muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py  # Add template mode function
muc_one_up/analysis/snapshot_validator.py      # Refactor to use sequence_utils
muc_one_up/config.py                          # Add "amplicon" to config schema
config.json                                   # Add amplicon_params section
```

## Constants

Added to `constants.py`:

```python
DEFAULT_AMPLICON_TIMEOUT: Final[int] = 3600
VALID_PCR_PRESETS: Final[set[str]] = {"default", "no_bias"}
DEFAULT_PCR_PRESET: Final[str] = "default"
```

## Testing Strategy

- **Unit tests:** PCR bias model (deterministic mode, stochastic mode, preset loading, equal-length alleles produce 50/50, extreme length differences produce near-dropout), amplicon extractor (primer found, not found, multiple sites, RC handling, size validation), template generator (correct record count, naming convention), sequence utils (binding site detection, RC search)
- **Integration tests:** Pipeline with mocked external tools (PBSIM3, CCS, minimap2, samtools). Verify correct command construction, stage sequencing, diploid vs haploid branching.
- **Refactoring tests:** Existing `snapshot_validator` tests pass unchanged after `_find_binding_sites()` extraction.
- **No end-to-end tests** requiring actual PBSIM3 installation — consistent with existing test patterns.

## Out of Scope

- ONT amplicon simulation (follow-up issue)
- Read-source tracking for amplicon mode (requires coordinate remapping from template-space to original-haplotype-space; follow-up)
- Fuzzy/mismatch-tolerant primer matching
- Multiple amplicon regions per run
- GC-content bias modeling
- PCR chimera simulation
- Primer dimer modeling
