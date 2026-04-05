# ONT Amplicon Simulation Design

## Goal

Extend the existing PacBio amplicon simulation pipeline to support ONT reads using pbsim3 with ONT error models. Users select the platform via `--platform ont` on the existing `reads amplicon` command.

## Background

NanoSim cannot generate full-length amplicon reads (it samples random-length substrings from the reference). pbsim3's template mode (`--strategy templ`) generates one full-length read per template copy and ships with ONT-specific models (`ERRHMM-ONT.model`, `ERRHMM-ONT-HQ.model`, `QSHMM-ONT.model`, `QSHMM-ONT-HQ.model`). This makes pbsim3 the right tool for both PacBio and ONT amplicon simulation.

## Architecture

Extract shared amplicon stages (haplotype extraction, primer-based amplicon extraction, PCR bias, template generation) into `amplicon_common.py`. The existing PacBio pipeline and new ONT pipeline each compose these shared stages with platform-specific read generation and alignment.

```
amplicon_common.py        → shared stages 1-3 (extract, bias, templates)
amplicon_pipeline.py      → PacBio: multi-pass pbsim3 + CCS + map-hifi  (existing, refactored)
ont_amplicon_pipeline.py  → ONT: single-pass pbsim3 + no CCS + map-ont  (new)
```

## 1. Shared Amplicon Helpers

New module: `muc_one_up/read_simulator/amplicon_common.py`

### `AmpliconPrep` dataclass

```python
@dataclass
class AmpliconPrep:
    allele_templates: list[Path]   # template FASTA paths (1 or 2)
    allele_coverages: list[int]    # per-allele read counts from PCR bias
    output_dir: Path
    output_base: str
    intermediate_files: list[str]  # for cleanup
    is_diploid: bool
```

### `extract_and_prepare_amplicons()` function

Extracted from `amplicon_pipeline.py` stages 1-3:
1. Detect diploid reference, split haplotypes if needed
2. Extract amplicon regions via forward/reverse primer binding sites (`AmpliconExtractor`)
3. Compute PCR bias coverage split per allele (`PCRBiasModel`)
4. Generate template FASTAs with N copies per allele (`TemplateGenerator`)
5. Return `AmpliconPrep`

Parameters: `config`, `input_fa`, `output_config`, plus the resolved tool paths and amplicon/pacbio params needed for extraction.

The existing `amplicon_pipeline.py` is refactored to call this helper, then perform PacBio-specific stages. No behavior change.

## 2. pbsim3 Wrapper Change

`muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py:383` currently rejects `pass_num < 2`:

```python
if pass_num < 2:
    raise FileOperationError(
        f"Invalid pass_num: {pass_num}. Multi-pass simulation requires pass_num >= 2"
    )
```

This must be relaxed to allow `pass_num = 1` for ONT single-pass simulation. Change the validation to `pass_num < 1` and update the error message:

```python
if pass_num < 1:
    raise FileOperationError(
        f"Invalid pass_num: {pass_num}. Template simulation requires pass_num >= 1"
    )
```

This is backward-compatible: PacBio callers already pass `pass_num >= 3`. The ONT pipeline passes `pass_num = 1`.

## 3. ONT Amplicon Pipeline

New module: `muc_one_up/read_simulator/ont_amplicon_pipeline.py`

Function: `simulate_ont_amplicon_pipeline(config, input_fa, human_reference=None, source_tracker=None, output_config=None) -> str`

### Pipeline stages

| Stage | What | Tool |
|-------|------|------|
| 1-3 | `extract_and_prepare_amplicons()` | AmpliconExtractor, PCRBiasModel, TemplateGenerator |
| 4 | pbsim3 template mode, ONT model, `--pass-num 1` (per allele) | pbsim3 |
| 5 | Per-allele BAM → FASTQ conversion | samtools |
| 6 | Concatenate per-allele FASTQs (if diploid) | cat (Python-level) |
| 7 | Align merged FASTQ to reference with `map-ont` preset | minimap2 |
| 8 | Write metadata, cleanup intermediates | — |

Note on stage 4-5: `run_pbsim3_template_simulation()` already returns BAM paths (it handles SAM→BAM internally). The ONT pipeline converts each per-allele BAM to FASTQ via `convert_bam_to_fastq()`, then concatenates the FASTQs. This is simpler than the PacBio path which runs CCS on each BAM before merging.

### Key differences from PacBio amplicon

| Parameter | PacBio | ONT |
|-----------|--------|-----|
| `--pass-num` | 3+ (multi-pass CLR) | 1 (single-pass) |
| CCS consensus | Yes (multi-pass → HiFi) | No (skip entirely) |
| minimap2 preset | `map-hifi` | `map-ont` |
| Model files | `QSHMM-SEQUEL.model` etc. | `ERRHMM-ONT-HQ.model` etc. |
| Output names | `*_amplicon_hifi.bam` | `*_amplicon_ont.bam` |

### Output naming

| File | Pattern |
|------|---------|
| Aligned BAM | `{output_base}_amplicon_ont.bam` |
| Merged FASTQ | `{output_base}_amplicon_ont.fastq` |
| Metadata | `{output_base}_amplicon_ont_metadata.tsv` |

This avoids collision with PacBio amplicon names (`*_amplicon_hifi.*`).

### pbsim3 invocation

Uses `run_pbsim3_template_simulation()` from `pbsim3_wrapper.py` with:
- `model_type`: from config (typically `"errhmm"` for ONT)
- `model_file`: user-provided path to ONT model file
- `pass_num`: `1` (single-pass ONT)
- `accuracy_mean`: from config or pbsim3 default
- All other params same as PacBio amplicon

## 4. Metadata Writer

The metadata writer (`muc_one_up/read_simulator/utils/metadata_writer.py:141`) currently only handles `"Illumina"`, `"ONT"`, and `"PacBio"` platform strings. The existing PacBio amplicon pipeline passes `platform="Amplicon"` which already falls through without writing platform-specific parameters.

**Decision:** Use `platform="PacBio"` for PacBio amplicon and `platform="ONT"` for ONT amplicon. Add separate `Assay_type\tamplicon` metadata field for both. This reuses existing platform-specific parameter writing (coverage, pass_num for PacBio; coverage for ONT) without needing new platform branches.

Changes to metadata writer:
- After the platform-specific block, write `Assay_type` if present in config:
  ```python
  assay = config.get("read_simulation", {}).get("assay_type")
  if assay:
      f.write(f"Assay_type\t{assay}\n")
  ```
- The CLI sets `config["read_simulation"]["assay_type"] = "amplicon"` for both PacBio and ONT amplicon commands.
- The existing PacBio amplicon pipeline changes from `platform="Amplicon"` to `platform="PacBio"` with `assay_type="amplicon"` in config.

## 5. CLI and Dispatch

### CLI option

Add `--platform` to `reads amplicon` command in `muc_one_up/cli/commands/reads.py`:

```python
@click.option(
    "--platform",
    type=click.Choice(["pacbio", "ont"]),
    default="pacbio",
    show_default=True,
    help="Sequencing platform for amplicon simulation.",
)
```

- `--platform pacbio` (default): sets `simulator = "amplicon"`, routes to existing `simulate_amplicon_reads_pipeline()`
- `--platform ont`: sets `simulator = "ont-amplicon"`, routes to `simulate_ont_amplicon_pipeline()`
- Both set `config["read_simulation"]["assay_type"] = "amplicon"` for metadata
- `--model-file` required for both platforms (consistent UX)
- `--track-read-source` rejected for both platforms (existing behavior)

### Dispatch

Add `"ont-amplicon"` branch to `_get_simulator()` in `read_simulation.py`:

```python
elif simulator_type == "ont-amplicon":
    from muc_one_up.read_simulator.ont_amplicon_pipeline import (
        simulate_ont_amplicon_pipeline,
    )
    return lambda config, input_fa, human_reference, **kw: simulate_ont_amplicon_pipeline(
        config, input_fa, human_reference=human_reference, **kw
    )
```

### Config schema

Add `"ont-amplicon"` to the simulator enum in `muc_one_up/config.py:285`:

```python
"simulator": {"type": "string", "enum": ["illumina", "ont", "pacbio", "amplicon", "ont-amplicon"]},
```

### Config structure

ONT amplicon uses the same config sections as PacBio amplicon:
- `amplicon_params`: primers, expected_product_range, pcr_bias (shared)
- `pacbio_params`: model_type, model_file (pbsim3 config — name is legacy but applies to both)
- `read_simulation`: simulator, coverage, assay_type
- `tools`: pbsim3, samtools, minimap2 (no `ccs` required for ONT)

## 6. Files Changed

| File | Action | Responsibility |
|------|--------|----------------|
| `muc_one_up/read_simulator/amplicon_common.py` | Create | `AmpliconPrep` dataclass, `extract_and_prepare_amplicons()` |
| `muc_one_up/read_simulator/ont_amplicon_pipeline.py` | Create | ONT amplicon pipeline function |
| `muc_one_up/read_simulator/amplicon_pipeline.py` | Modify | Refactor to call shared helper; change metadata platform to `"PacBio"` |
| `muc_one_up/read_simulator/wrappers/pbsim3_wrapper.py` | Modify | Relax `pass_num` validation from `< 2` to `< 1` |
| `muc_one_up/read_simulator/utils/metadata_writer.py` | Modify | Add `assay_type` field support |
| `muc_one_up/read_simulation.py` | Modify | Add `"ont-amplicon"` to `_get_simulator()` |
| `muc_one_up/config.py` | Modify | Add `"ont-amplicon"` to simulator enum in schema |
| `muc_one_up/cli/commands/reads.py` | Modify | Add `--platform` option, set assay_type in config |
| `tests/read_simulator/test_amplicon_common.py` | Create | Unit tests for shared extraction helper |
| `tests/read_simulator/test_ont_amplicon_pipeline.py` | Create | Unit tests for ONT amplicon pipeline (mocked tools) |
| `tests/read_simulator/test_metadata_writer.py` | Modify | Add test for assay_type field |
| `tests/test_lazy_imports.py` | Modify | Add `"ont-amplicon"` to dispatch parametrize list |
| `tests/cli/test_amplicon_platform.py` | Create | CLI routing tests for `--platform` flag |
| `tests/e2e/test_pipeline_e2e.py` | Modify | Add `TestOntAmpliconE2E` class |

## 7. Testing

**Unit tests:**
- `tests/read_simulator/test_amplicon_common.py`: Verify `extract_and_prepare_amplicons()` returns correct `AmpliconPrep` for diploid/haploid inputs with mocked extraction.
- `tests/read_simulator/test_ont_amplicon_pipeline.py`: Verify pipeline calls pbsim3 with `pass-num=1`, does NOT invoke CCS, uses `map-ont` alignment preset. All tools mocked.
- `tests/read_simulator/test_metadata_writer.py`: Add test verifying `assay_type` is written when present in config.

**Dispatch tests:**
- `tests/test_lazy_imports.py`: Add `"ont-amplicon"` to the `@pytest.mark.parametrize` list in `TestGetSimulatorDispatch.test_valid_simulator_returns_callable`.

**CLI tests:**
- `tests/cli/test_amplicon_platform.py`: Invoke `reads amplicon --platform ont` via CliRunner, verify config sets `simulator="ont-amplicon"` and `assay_type="amplicon"`. Verify default `--platform` routes to PacBio (backward compat).

**E2e tests:**
- `TestOntAmpliconE2E` in `tests/e2e/test_pipeline_e2e.py`: gated by `@pytest.mark.requires_tools("pbsim3", "minimap2", "samtools")` — note: no `ccs` required. Skips if tools/models not installed.

## 8. Out of Scope

- NanoSim-based amplicon simulation (NanoSim cannot generate full-length reads)
- Renaming `pacbio_params` config section (legacy name applies to pbsim3 regardless of platform)
- Read source tracking for amplicon mode (existing limitation, applies to both platforms)
- Basecaller simulation (Guppy/Dorado quality polishing)

## Success Criteria

- `reads amplicon --platform ont --model-file /path/to/ERRHMM-ONT-HQ.model` produces aligned BAM with full-length ONT-error-profile reads
- `reads amplicon` (no `--platform`) works identically to current PacBio behavior
- Existing PacBio amplicon tests pass unchanged after refactor
- pbsim3 wrapper accepts `pass_num=1` without error
- Metadata includes `Assay_type: amplicon` for both platforms
- Config schema validates `"ont-amplicon"` as a valid simulator
- `_get_simulator("ont-amplicon")` returns the ONT amplicon pipeline callable
- mypy green, ruff clean, all tests pass
- No CCS dependency for ONT amplicon (only pbsim3 + samtools + minimap2)
