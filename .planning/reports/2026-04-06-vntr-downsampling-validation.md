# VNTR Downsampling Validation Report

**Date:** 2026-04-06
**MucOneUp version:** v0.44.0 + fix/89-vntr-downsample-two-step branch
**downsample_mode:** vntr (two-step)
**coverage target:** 150

---

## Datasets

| Dataset | n | Description |
|---------|---|-------------|
| **Real exomes** | 1,043 | CerKiD Berlin Twist v2, chr1 MUC1 region extracts |
| **Simulated** | 50 | MucOneUp, VNTR lengths 30-90 repeats, downsample_mode=vntr, coverage=150 |

All coverage measurements use `samtools depth -a` (includes zero-coverage positions).

- VNTR region: chr1:155,188,487-155,192,239 (hg38, 3,753 bp)
- Non-VNTR BED: `data/bed/twist_v2_muc1_no_vntr_targets.bed` (765 bp, 7 Twist probe regions)
- Real exome BAMs: CerKiD Berlin cohort, pre-sliced chr1 MUC1 region extracts (vntyper-analyses)
- Simulated BAMs: 50 MucOneUp runs with VNTR lengths 30-90

---

## Real Exome Descriptive Statistics (n=1,043)

| Metric | Mean | Median | SD | Min | Q5 | Q25 | Q75 | Q95 | Max |
|--------|------|--------|-----|-----|-----|------|------|------|-----|
| Total reads | 31,851 | 29,784 | 12,310 | 7,328 | 15,921 | 23,074 | 38,076 | 56,547 | 93,964 |
| VNTR mean (x) | 182.0 | 168.7 | 85.8 | 0.0 | 74.7 | 122.4 | 223.0 | 340.4 | 689.6 |
| VNTR median (x) | 111.8 | 106.0 | 51.2 | 0.0 | 43.0 | 75.0 | 140.0 | 206.0 | 402.0 |
| VNTR % uncovered | 10.3 | 9.9 | 2.9 | 0.0 | 6.6 | 8.4 | 11.7 | 14.7 | 31.6 |
| Non-VNTR BED (x) | 110.2 | 107.0 | 25.8 | 55.5 | 73.0 | 91.4 | 125.3 | 155.9 | 240.5 |
| VNTR:BED ratio | 1.7 | 1.6 | 0.8 | 0.0 | 0.8 | 1.2 | 2.0 | 2.7 | 9.0 |

---

## Simulated Descriptive Statistics (n=50)

| Metric | Mean | Median | SD | Min | Q5 | Q25 | Q75 | Q95 | Max |
|--------|------|--------|-----|-----|-----|------|------|------|-----|
| Total reads | 8,844 | 8,827 | 242 | 8,396 | 8,458 | 8,634 | 9,028 | 9,306 | 9,386 |
| VNTR mean (x) | 204.1 | 203.6 | 7.2 | 189.9 | 192.5 | 197.4 | 209.5 | 215.9 | 218.6 |
| VNTR median (x) | 81.2 | 77.0 | 35.2 | 6.0 | 26.0 | 54.0 | 110.0 | 143.0 | 179.0 |
| VNTR % uncovered | 24.9 | 24.4 | 3.0 | 21.7 | 22.2 | 23.4 | 25.7 | 26.6 | 39.7 |
| Non-VNTR BED (x) | 149.8 | 149.8 | 4.9 | 141.1 | 141.6 | 146.0 | 155.1 | 158.6 | 160.0 |
| VNTR:BED ratio | 1.4 | 1.4 | 0.0 | 1.3 | 1.3 | 1.4 | 1.4 | 1.4 | 1.4 |

---

## Side-by-Side Comparison

| Metric                 |   R.Mean |    R.Med |     R.SD |     R.Q5 |    R.Q95 |   S.Mean |    S.Med |     S.SD |     S.Q5 |    S.Q95 |
|------------------------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| Total reads            |   31,851 |   29,784 |   12,310 |   15,921 |   56,547 |    8,844 |    8,827 |      242 |    8,458 |    9,306 |
| VNTR mean (x)          |    182.0 |    168.7 |     85.8 |     74.7 |    340.4 |    204.1 |    203.6 |      7.2 |    192.5 |    215.9 |
| VNTR median (x)        |    111.8 |    106.0 |     51.2 |     43.0 |    206.0 |     81.2 |     77.0 |     35.2 |     26.0 |    143.0 |
| VNTR % uncovered       |     10.3 |      9.9 |      2.9 |      6.6 |     14.7 |     24.9 |     24.4 |      3.0 |     22.2 |     26.6 |
| Non-VNTR BED (x)       |    110.2 |    107.0 |     25.8 |     73.0 |    155.9 |    149.8 |    149.8 |      4.9 |    141.6 |    158.6 |
| VNTR:BED ratio         |      1.7 |      1.6 |      0.8 |      0.8 |      2.7 |      1.4 |      1.4 |      0.0 |      1.3 |      1.4 |

R = Real exomes (n=1,043), S = Simulated (n=50)

---

## Match Assessment

| Metric | Real Median | Sim Median | Sim/Real | Verdict |
|--------|------------|------------|----------|---------|
| Total reads | 29,784 | 8,827 | 0.30 | POOR |
| VNTR mean | 168.7x | 203.6x | 1.21 | GOOD |
| VNTR median | 106.0x | 77.0x | 0.73 | GOOD |
| VNTR % uncovered | 9.9% | 24.4% | 2.48 | POOR |
| Non-VNTR BED | 107.0x | 149.8x | 1.40 | GOOD |
| VNTR:BED ratio | 1.6x | 1.4x | 0.87 | GOOD |

**GOOD** = sim/real ratio within 0.7-1.5x. **ACCEPTABLE** = 0.5-2.0x. **POOR** = outside 2.0x.

---

## Simulated Coverage by VNTR Length

| VNTR len | Reads | VNTR mean | VNTR med | % uncov | BED mean | Ratio |
|----------|------:|----------:|---------:|--------:|---------:|------:|
| 30 | 8,935 | 207.1 | 78 | 24.6 | 153.0 | 1.4 |
| 31 | 8,954 | 205.5 | 7 | 36.6 | 149.6 | 1.4 |
| 32 | 9,013 | 209.4 | 143 | 22.2 | 151.2 | 1.4 |
| 33 | 8,872 | 204.4 | 92 | 23.2 | 148.8 | 1.4 |
| 34 | 9,014 | 206.6 | 114 | 22.0 | 151.5 | 1.4 |
| 36 | 8,745 | 197.7 | 75 | 24.0 | 144.5 | 1.4 |
| 37 | 9,105 | 213.2 | 77 | 24.4 | 155.1 | 1.4 |
| 38 | 8,944 | 203.1 | 179 | 23.1 | 147.3 | 1.4 |
| 39 | 8,672 | 200.0 | 67 | 24.7 | 145.9 | 1.4 |
| 41 | 9,319 | 215.1 | 34 | 22.4 | 154.2 | 1.4 |
| 42 | 8,624 | 198.2 | 144 | 23.3 | 146.1 | 1.4 |
| 43 | 8,732 | 200.4 | 120 | 24.9 | 148.2 | 1.4 |
| 44 | 9,114 | 213.6 | 110 | 24.8 | 156.5 | 1.4 |
| 45 | 8,396 | 189.9 | 82 | 22.8 | 141.1 | 1.3 |
| 47 | 9,306 | 218.6 | 126 | 22.6 | 160.0 | 1.4 |
| 48 | 9,162 | 214.2 | 115 | 23.6 | 155.9 | 1.4 |
| 49 | 8,684 | 195.2 | 86 | 23.2 | 142.2 | 1.4 |
| 50 | 8,434 | 190.4 | 61 | 24.4 | 141.6 | 1.3 |
| 52 | 8,892 | 206.6 | 83 | 25.7 | 150.8 | 1.4 |
| 53 | 8,624 | 195.4 | 6 | 39.7 | 145.1 | 1.3 |
| 54 | 9,105 | 211.4 | 54 | 25.8 | 155.1 | 1.4 |
| 55 | 9,128 | 215.9 | 42 | 26.5 | 158.6 | 1.4 |
| 56 | 8,896 | 203.4 | 125 | 24.2 | 149.3 | 1.4 |
| 58 | 9,386 | 217.4 | 135 | 23.4 | 158.9 | 1.4 |
| 59 | 8,716 | 201.4 | 111 | 24.0 | 148.4 | 1.4 |
| 60 | 8,651 | 198.9 | 102 | 23.8 | 147.0 | 1.4 |
| 61 | 8,904 | 207.6 | 77 | 24.5 | 152.3 | 1.4 |
| 63 | 8,800 | 206.2 | 69 | 24.2 | 151.0 | 1.4 |
| 64 | 8,956 | 209.5 | 73 | 25.0 | 155.6 | 1.3 |
| 65 | 8,795 | 203.6 | 37 | 24.6 | 150.0 | 1.4 |
| 66 | 8,483 | 192.5 | 36 | 25.4 | 142.0 | 1.4 |
| 67 | 8,956 | 203.7 | 43 | 26.2 | 150.4 | 1.4 |
| 69 | 8,458 | 195.0 | 76 | 24.2 | 143.7 | 1.4 |
| 70 | 8,983 | 214.8 | 75 | 25.2 | 157.2 | 1.4 |
| 71 | 8,778 | 200.2 | 94 | 23.7 | 147.0 | 1.4 |
| 72 | 8,488 | 194.3 | 75 | 24.2 | 141.3 | 1.4 |
| 74 | 8,705 | 203.3 | 60 | 23.6 | 150.1 | 1.4 |
| 75 | 8,521 | 197.0 | 127 | 24.8 | 146.3 | 1.3 |
| 76 | 9,115 | 206.0 | 92 | 24.4 | 152.2 | 1.4 |
| 77 | 8,832 | 201.4 | 53 | 26.2 | 148.3 | 1.4 |
| 78 | 8,634 | 201.9 | 26 | 25.7 | 147.5 | 1.4 |
| 80 | 9,028 | 208.9 | 46 | 26.5 | 152.4 | 1.4 |
| 81 | 9,177 | 213.7 | 75 | 25.2 | 156.4 | 1.4 |
| 82 | 8,627 | 197.4 | 77 | 25.5 | 146.0 | 1.4 |
| 83 | 8,609 | 198.7 | 88 | 26.6 | 146.5 | 1.4 |
| 85 | 9,029 | 208.6 | 87 | 23.4 | 154.1 | 1.4 |
| 86 | 8,729 | 204.1 | 75 | 24.7 | 150.9 | 1.4 |
| 87 | 8,617 | 196.1 | 65 | 24.7 | 143.4 | 1.4 |
| 88 | 8,822 | 204.5 | 86 | 23.6 | 150.8 | 1.4 |
| 90 | 8,722 | 202.1 | 81 | 21.7 | 148.3 | 1.4 |

---

## File Locations

### Simulation outputs

Each simulation output directory (`len_{N}/`) contains:

```
sim_{N}.001.simulated.fa                              # Input haplotype FASTA
sim_{N}.001.simulated_reads_downsampled.bam           # Final downsampled BAM
sim_{N}.001.simulated_reads_downsampled.bam.bai       # BAM index
sim_{N}.001.simulated_reads_R1.fastq.gz               # Raw paired FASTQ R1
sim_{N}.001.simulated_reads_R2.fastq.gz               # Raw paired FASTQ R2
sim_{N}.001.simulated_reads_vntr_biased_R1.fastq.gz   # VNTR-biased FASTQ R1
sim_{N}.001.simulated_reads_vntr_biased_R2.fastq.gz   # VNTR-biased FASTQ R2
sim_{N}.001.simulated_reads_vntr_efficiency_stats.json # VNTR efficiency stats
sim_{N}.001.simulated_reads_metadata.tsv              # Simulation metadata
sim_{N}.001.simulation_stats.json                     # Haplotype simulation stats
```

Total: 50 simulations, 658 MB disk.

### Real exome BAMs

CerKiD Berlin cohort, pre-sliced to chr1:155184000-155194000 (from vntyper-analyses screening data).

### Raw data files

Per-sample and per-simulation metrics were saved as JSON during analysis (not committed to repo).

---

## Pipeline Configuration

```json
{
  "downsample_mode": "vntr",
  "coverage": 150,
  "read_number": 100000,
  "vntr_capture_efficiency": {
    "enabled": true,
    "penalty_factor": 0.39
  },
  "vntr_to_flanking_ratio": 1.4
}
```

### Two-step downsampling (vntr mode)

1. **Step 1:** Measure coverage in non-VNTR BED (`sample_target_bed`). Downsample entire BAM so non-VNTR BED = `coverage` target.
2. **Step 2:** Measure VNTR coverage after step 1. If VNTR > non-VNTR x `vntr_to_flanking_ratio`, downsample VNTR region to match ratio. Uses seed+1 to avoid samtools hash-based subsampling no-op.

---

## Known Limitations

1. **Total reads** (sim 8,827 vs real 29,784, ratio 0.30): simulated reads only cover ~4kb around the VNTR; real exomes capture across the full ~10kb region. The read generation pipeline (pBLAT + ReSeq) only produces fragments from the aligned footprint of the simulated FASTA. See issue #89.
2. **VNTR % uncovered** (sim 24.4% vs real 9.9%, ratio 2.48): simulation produces more zero-coverage positions in the VNTR due to the limited read generation footprint and higher mean-to-median ratio (skewed coverage distribution).
3. **Non-VNTR BED** (sim 150x vs real 107x): simulated non-VNTR hits the configured target (150x), while real exomes reflect natural sequencing depth (~107x). These are comparable when the target is set to match real data.
4. **VNTR:BED ratio** (sim 1.4 vs real 1.6): close match, controlled by `vntr_to_flanking_ratio` config parameter (default 1.4, derived from median of 20 CerKiD exomes). Can be tuned in config.
5. **Low variance in simulated metrics**: SD of simulated VNTR mean is 7.2x vs 85.8x in real data. Simulations use fixed `read_number` and deterministic `seed`, producing very consistent results. Real exomes vary by patient, DNA quality, and sequencing run.

---

## Conclusion

The two-step VNTR downsampling approach produces coverage profiles that match real exome data for the key benchmarking metrics:

| Metric | Status | Detail |
|--------|--------|--------|
| VNTR mean coverage | GOOD | 204x sim vs 169x real (1.2x) |
| VNTR:flanking ratio | GOOD | 1.4x sim vs 1.6x real (0.87x) |
| Non-VNTR BED coverage | GOOD | 150x sim, hits target |
| VNTR median coverage | GOOD | 77x sim vs 106x real (0.73x) |
| Total reads | POOR | 8.8k sim vs 29.8k real (0.30x) -- architectural limitation |
| VNTR % uncovered | POOR | 24% sim vs 10% real (2.5x) -- architectural limitation |

The POOR metrics are architectural limitations of the read generation pipeline (issue #89), not the downsampling logic. For variant caller benchmarking, the GOOD metrics (VNTR coverage, ratio, non-VNTR coverage) are the ones that matter most.
