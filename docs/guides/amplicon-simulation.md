# Amplicon Read Simulation Guide

Simulate PacBio or ONT amplicon reads with realistic PCR length bias from MucOneUp diploid references.

---

## Overview

The `reads amplicon` command produces full-length amplicon reads that mimic targeted long-range PCR products spanning the MUC1 VNTR region. This enables benchmarking of amplicon-based analysis pipelines such as VNTRPipeline (Madritsch et al. 2026) and PacMUCI (Vrbacka et al. 2025).

**Key features:**

- Primer-based amplicon extraction from simulated haplotypes
- PCR length bias modeling (shorter alleles amplify preferentially)
- PBSIM3 template mode for full-length reads with realistic error profiles
- PacBio (multi-pass CCS) and ONT (single-pass) platform support
- Diploid and haploid input support

**How it differs from `reads pacbio`:**

| Feature | `reads pacbio` (WGS) | `reads amplicon` |
|---------|----------------------|------------------|
| Strategy | Random genome sampling | Full-length templates |
| Read length | Variable (length distribution) | Fixed (amplicon length) |
| Coverage model | Uniform across reference | PCR bias (length-dependent) |
| Use case | WGS benchmarking | Amplicon pipeline benchmarking |

---

## Quick Start

```bash
# 1. Generate diploid reference with different VNTR lengths
muconeup --config config.json simulate \
  --fixed-lengths 40 --fixed-lengths 70 \
  --seed 42 --out-base output/sample

# 2a. PacBio amplicon reads (default)
muconeup --config config.json reads amplicon \
  output/sample.001.simulated.fa \
  --model-file /path/to/QSHMM-RSII.model \
  --coverage 500 --seed 42

# 2b. ONT amplicon reads
muconeup --config config.json reads amplicon --platform ont \
  output/sample.001.simulated.fa \
  --model-file /path/to/ERRHMM-ONT-HQ.model \
  --coverage 500 --seed 42
```

**Output:** Aligned BAM with ~500 reads, biased toward the shorter (40-repeat) allele.

---

## Pipeline Stages

The amplicon pipeline runs 8 stages:

```
Input FASTA (diploid)
  |
  v
1. Haplotype extraction ──> hap1.fa, hap2.fa
  |
  v
2. Amplicon extraction ──> amplicon_hap1.fa, amplicon_hap2.fa
   (primer binding site detection)
  |
  v
3. PCR bias computation ──> e.g., 350 reads for hap1, 150 for hap2
   (shorter allele gets more reads)
  |
  v
4. Template FASTA generation ──> 350 copies of hap1, 150 copies of hap2
  |
  v
5. PBSIM3 template mode ──> CLR reads (one per template copy)
  |
  v
6. CCS consensus ──> HiFi reads
  |
  v
7. BAM merge ──> single BAM with reads from both alleles
  |
  v
8. Alignment (minimap2) ──> final aligned BAM
```

For **haploid** input (single sequence), stages 1 and 3 are skipped.

---

## Platform Selection

Use `--platform` to choose between PacBio and ONT amplicon simulation:

```bash
# PacBio HiFi (default)
muconeup --config config.json reads amplicon sample.fa \
  --model-file /path/to/QSHMM-SEQUEL.model

# Oxford Nanopore
muconeup --config config.json reads amplicon --platform ont sample.fa \
  --model-file /path/to/ERRHMM-ONT-HQ.model
```

Stages 1-4 (extraction, PCR bias, template generation) are shared. The platforms differ in read generation and alignment:

| Stage | PacBio (`--platform pacbio`) | ONT (`--platform ont`) |
|-------|------------------------------|------------------------|
| pbsim3 pass_num | 3+ (multi-pass CLR) | 1 (single-pass) |
| Consensus | CCS (multi-pass -> HiFi) | None (skip) |
| minimap2 preset | `map-hifi` | `map-ont` |
| Error model | `QSHMM-SEQUEL.model` etc. | `ERRHMM-ONT-HQ.model` etc. |
| Output suffix | `*_amplicon_hifi.bam` | `*_amplicon_ont.bam` |
| Tool dependencies | pbsim3, ccs, samtools, minimap2 | pbsim3, samtools, minimap2 |

!!! note "ONT does not require CCS"
    The ONT pipeline skips CCS consensus entirely. Each pbsim3 template produces one single-pass read, which is converted directly to FASTQ and aligned. This means `--coverage` maps directly to the number of output reads (before alignment filtering).

---

## PCR Length Bias Model

In long-range PCR, shorter amplicons amplify more efficiently than longer ones. This is modeled as exponential decay of per-cycle efficiency:

```
E(L) = E_max * exp(-alpha * L)
```

For two competing alleles, the yield ratio after `c` PCR cycles:

```
ratio = ((1 + E1) / (1 + E2)) ^ c
```

### Example

With the `default` preset (calibrated to Madritsch et al. 2026):

| Allele lengths | Shorter allele fraction | Longer allele fraction |
|---------------|------------------------|----------------------|
| 2000 vs 2000 bp | 50% | 50% |
| 2000 vs 3500 bp | 70% | 30% |
| 1800 vs 4800 bp | 84% | 16% |

### Presets

| Preset | Description | Use case |
|--------|-------------|----------|
| `default` | KOD HS, 25 cycles, calibrated to ~70% bias at 1.5kb difference | Realistic amplicon benchmarking |
| `no_bias` | Equal 50/50 split regardless of length | Control experiments |

```bash
# Use default preset (realistic PCR bias)
muconeup --config config.json reads amplicon sample.fa \
  --model-file model.model

# No bias (equal allele coverage)
muconeup --config config.json reads amplicon sample.fa \
  --model-file model.model --pcr-preset no_bias
```

### Stochastic Mode

By default, the bias is deterministic (same input = same split). Enable stochastic mode for realistic run-to-run variance:

```bash
muconeup --config config.json reads amplicon sample.fa \
  --model-file model.model --stochastic-pcr --seed 42
```

This uses a Galton-Watson branching process where each PCR cycle is modeled as independent binomial trials per molecule.

### Custom Parameters

Override individual PCR parameters in `config.json`:

```json
"amplicon_params": {
    "forward_primer": "GGAGAAAAGGAGACTTCGGCTACCCAG",
    "reverse_primer": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    "pcr_bias": {
        "preset": "default",
        "cycles": 30,
        "stochastic": true
    }
}
```

Or specify all parameters explicitly (no preset):

```json
"pcr_bias": {
    "e_max": 0.90,
    "alpha": 0.00005,
    "cycles": 20,
    "denaturation_time": 15.0,
    "stochastic": false
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `e_max` | float | 0.95 | Max per-cycle efficiency (short amplicons) |
| `alpha` | float | 0.00005 | Length decay rate (bp^-1) |
| `cycles` | int | 25 | Number of PCR cycles |
| `denaturation_time` | float | 10.0 | Denaturation step duration (seconds) |
| `stochastic` | bool | false | Enable Galton-Watson branching |

---

## Primer Configuration

Amplicon boundaries are defined by primer binding sites in `config.json`:

```json
"amplicon_params": {
    "forward_primer": "GGAGAAAAGGAGACTTCGGCTACCCAG",
    "reverse_primer": "GCCGTTGTGCACCAGAGTAGAAGCTGA",
    "primer_source": "Wenzel et al. 2018 (PMID: 29520014)",
    "expected_product_range": [1500, 15000]
}
```

The primers are the core sequences from Wenzel et al. 2018 (PS2/PS3), used across all published MUC1 VNTR PCR protocols (Wenzel 2018, Madritsch 2026, Vrbacka 2025).

### How extraction works

1. Forward primer is searched on the forward strand (exact match, case-insensitive)
2. Reverse primer is reverse-complemented and searched on the forward strand
3. The region from forward primer start to reverse primer end (inclusive) is extracted
4. If either primer is not found or found multiple times, the pipeline raises an error

### Size validation

`expected_product_range` is optional. When set, the extracted amplicon length is validated:

```json
"expected_product_range": [1500, 15000]
```

This means 1500 <= amplicon_length <= 6000 bp. Amplicons outside this range raise an error.

---

## Coverage Semantics

`--coverage` specifies the total number of **template molecules**. For PacBio, the actual HiFi read count may be lower due to CCS quality filtering (`min_rq`, `min_passes`). For ONT, each template produces one output read (no CCS filtering).

For diploid inputs, the total is split between alleles by the PCR bias model:

```bash
# 500 total template molecules, split by PCR bias
muconeup --config config.json reads amplicon sample.fa \
  --model-file model.model --coverage 500
```

With a 70%/30% bias, this produces ~350 templates for the shorter allele and ~150 for the longer allele.

---

## CLI Reference

```bash
muconeup --config CONFIG reads amplicon [OPTIONS] INPUT_FASTAS...

Required:
  INPUT_FASTAS              One or more FASTA files

Options:
  --platform [pacbio|ont]    Sequencing platform (default: pacbio)
  --model-file PATH         pbsim3 model file (overrides config)
  --model-type [qshmm|errhmm]  pbsim3 model type (overrides config)
  --pcr-preset [default|no_bias]  PCR bias preset
  --stochastic-pcr          Enable stochastic PCR bias
  --coverage INT             Total template molecules (default: 30)
  --seed INT                 Random seed for reproducibility
  --out-dir PATH             Output directory (default: .)
  --out-base TEXT             Output base name
  --track-read-source         Accepted but raises an error (not supported in amplicon mode)
```

!!! warning "Read Source Tracking"
    The `--track-read-source` option is accepted by the CLI parser but raises an error at runtime. Read source tracking is not supported in amplicon mode.

---

## Limitations

- **Read source tracking** (`--track-read-source`) is not supported in amplicon mode
- Primer matching is exact (no mismatch tolerance)
- Single amplicon region per run

---

## References

- Wenzel et al. 2018. *Sci Rep* 8:4170. DOI: [10.1038/s41598-018-22428-0](https://doi.org/10.1038/s41598-018-22428-0)
- Madritsch et al. 2026. *Sci Rep* 16:762. DOI: [10.1038/s41598-025-30441-3](https://doi.org/10.1038/s41598-025-30441-3)
- Vrbacka et al. 2025. *bioRxiv*. DOI: [10.1101/2025.09.06.673538](https://doi.org/10.1101/2025.09.06.673538)
- Suzuki & Giovannoni 1996 -- Competitive PCR bias model
- PBSIM3 template mode: [github.com/yukiteruono/pbsim3](https://github.com/yukiteruono/pbsim3)
