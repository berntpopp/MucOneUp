# Read Source Tracking

## Overview

The `--track-read-source` flag annotates every simulated read with its ground-truth origin: which haplotype it came from, where on the reference it was sampled, and whether it overlaps a mutation, SNP, or specific VNTR repeat unit. This produces a compressed TSV manifest alongside the standard output files.

**Use cases:** benchmarking variant callers, measuring aligner accuracy, quantifying allelic balance, and evaluating VNTR-aware pipelines against known ground truth.

---

## Quick Start

```bash
# Step 1: simulate haplotypes
muconeup --config config.json simulate \
  --out-base sample --fixed-lengths 30 \
  --mutation-name dupC --mutation-targets 1,25

# Step 2: simulate Illumina reads with source tracking
muconeup --config config.json reads illumina sample.001.simulated.fa \
  --track-read-source
```

Works identically with `reads ont` and `reads pacbio`.

---

## Output Files

When `--track-read-source` is enabled, two additional files are produced:

| File | Description |
|------|-------------|
| `{base}_repeat_coordinates.tsv` | Repeat unit coordinate map per haplotype |
| `{base}_read_manifest.tsv.gz` | Gzip-compressed per-read annotations |

### Read Manifest Schema

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | str | Read name from FASTQ/BAM |
| `haplotype` | int | Source haplotype (1 or 2) |
| `ref_start` | int | 0-based start on simulated reference |
| `ref_end` | int | 0-based end on simulated reference |
| `strand` | str | `+` or `-` |
| `overlaps_vntr` | bool | Read overlaps the VNTR region |
| `repeat_units` | str | Comma-separated repeat indices (e.g. `3,4,5`) or `.` |
| `overlaps_mutation` | bool | Read overlaps a mutated repeat |
| `mutation_name` | str | Mutation name (e.g. `dupC`) or `.` |
| `overlaps_snp` | bool | Read overlaps an applied SNP |

Booleans are serialized as lowercase `true`/`false`.

### Repeat Coordinates Schema

| Column | Type | Description |
|--------|------|-------------|
| `haplotype` | int | Haplotype number |
| `index` | int | 1-based repeat unit index |
| `repeat_type` | str | Repeat symbol (e.g. `X`, `Xm`) |
| `start` | int | 0-based start coordinate |
| `end` | int | 0-based end coordinate (exclusive) |
| `is_mutated` | bool | Whether this repeat carries a mutation |
| `mutation_name` | str | Mutation name or `.` |
| `snp_count` | int | Number of SNPs in this repeat |
| `snp_positions` | str | Comma-separated SNP positions or `.` |

---

## Querying the Manifest

### Which reads carry the mutation?

```python
import csv, gzip

with gzip.open("sample.001_read_manifest.tsv.gz", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    mutation_reads = [r for r in reader if r["overlaps_mutation"] == "true"]

print(f"Reads carrying mutation: {len(mutation_reads)}")
for r in mutation_reads:
    print(f"  {r['read_id']}  hap={r['haplotype']}  "
          f"pos={r['ref_start']}-{r['ref_end']}  repeats={r['repeat_units']}")
```

### Allelic balance

```python
import csv, gzip

with gzip.open("sample.001_read_manifest.tsv.gz", "rt") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

vntr_reads = [r for r in rows if r["overlaps_vntr"] == "true"]
mut_reads = [r for r in vntr_reads if r["overlaps_mutation"] == "true"]
print(f"VAF: {len(mut_reads)}/{len(vntr_reads)} = {len(mut_reads)/len(vntr_reads):.1%}")
```

### Command-line filtering

```bash
# Count mutation-carrying reads
zcat sample.001_read_manifest.tsv.gz | awk -F'\t' '$8=="true"' | wc -l

# Extract haplotype 1 reads only
zcat sample.001_read_manifest.tsv.gz | awk -F'\t' '$2==1'
```

---

## Platform-Specific Details

Read origin extraction varies by platform but produces the same manifest format:

| Platform | Source of read positions |
|----------|------------------------|
| **Illumina** | Fragment coordinates captured during Wessim2 simulation |
| **ONT** | NanoSim encodes haplotype and position in read names |
| **PacBio** | pbsim3 MAF alignment files parsed before cleanup |

### Standalone `reads` command

When calling `reads` directly (not via `simulate`), the tracker reconstructs simulation metadata from the companion `simulation_stats.json` file. If the companion file is missing, a warning is emitted and a partial manifest (positions only, no VNTR/mutation annotations) is produced.

---

## Performance

Tracking adds negligible overhead. The manifest is ~100 bytes per read; a 30x simulation of a 10 kb region (~15,000 reads) produces ~1.5 MB uncompressed. The gzip-compressed file is typically under 200 KB.

When `--track-read-source` is omitted (the default), zero tracking code executes.
