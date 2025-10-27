# Test Datasets

Standardized test datasets for validating MucOneUp workflows and benchmarking performance.

## Download

**Latest Release**: [testdata_40-70.tar.gz](https://github.com/berntpopp/MucOneUp/raw/main/data/testdata_releases/_latest/testdata_40-70.tar.gz)

The `_latest` directory always contains the most recent test dataset with a constant filename for persistent download links.

Extract with: `tar -xzf testdata_40-70.tar.gz`

Check version: `cat VERSION.json` (included in `_latest/` directory)

## Dataset Characteristics

| Property | Value |
|----------|-------|
| **VNTR Structure** | Asymmetric diploid (40/70 repeats) |
| **Variants** | Normal + dupC mutation (position 20) |
| **Platforms** | Illumina (50×), ONT (30×), PacBio HiFi (30×) |
| **Total Files** | ~55 files (~17 MB) |

## Directory Structure

```
testdata_40-70_{version}/
├── references/
│   ├── testdata_40-70.001.simulated.fa           # Normal diploid reference
│   ├── testdata_40-70.002.simulated.fa           # dupC mutated reference
│   ├── *_orf_analysis.json                       # ORF toxicity analysis
│   └── *_snapshot_validation.json                # SNaPshot validation (mutated only)
├── illumina/
│   ├── normal/
│   │   ├── reads_R1.fastq.gz
│   │   ├── reads_R2.fastq.gz
│   │   ├── aligned.bam
│   │   ├── aligned.bam.bai
│   │   └── metadata.tsv
│   └── dupC/
│       └── ...
├── ont/
│   ├── normal/
│   │   ├── reads.fastq                           # Uncompressed FASTQ
│   │   ├── aligned.bam
│   │   └── metadata.tsv
│   └── dupC/
│       └── ...
├── pacbio/
│   ├── normal/
│   │   ├── reads_hifi_0001.fastq.gz
│   │   ├── reads_hifi_0002.fastq.gz
│   │   ├── aligned.bam
│   │   └── vntr_efficiency_stats.json
│   └── dupC/
│       └── ...
├── metadata/
│   ├── dataset_metadata.json
│   ├── tool_versions.txt
│   └── generation.log
└── README.txt
```

## Generation

Regenerate test data using:

```bash
python helpers/generate_test_data.py \
    --version v{VERSION} \
    --config config.json \
    --platforms all \
    --threads 4 \
    --verbose
```

**Requirements**: All platform-specific tools (BWA, NanoSim, pbsim3) must be installed.

**Output**: Tarball in `data/testdata_releases/{date}_{version}/` with automatic `_latest` symlink for persistent downloads.

The script automatically updates the `_latest` symlink to point to the newest release, ensuring documentation links remain valid.

## Use Cases

- **Workflow validation**: Test complete pipeline from reference generation to read simulation
- **Platform comparison**: Compare alignment characteristics across Illumina, ONT, and PacBio
- **Mutation detection**: Validate SNaPshot analysis on dupC variant
- **ORF analysis**: Test toxic protein detection on VNTR references
- **Benchmarking**: Standardized dataset for performance testing

## Validation

Each platform directory includes `metadata.tsv` with read statistics:

- Read counts and length distributions
- Coverage uniformity across VNTR region
- Alignment metrics (MAPQ, error rates)

For Illumina and ONT, SNaPshot validation confirms dupC mutation detection in position 20.
