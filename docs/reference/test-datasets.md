# Test Datasets

Standardized test datasets for validating MucOneUp workflows and benchmarking performance.

## Download

Test datasets are published as [GitHub Releases](https://github.com/vntyper-org/MucOneUp/releases?q=testdata) tagged `testdata-<version>`, separate from the code releases.

**Current dataset**: [`testdata-v0.25.0`](https://github.com/vntyper-org/MucOneUp/releases/tag/testdata-v0.25.0)

```bash
curl -LO https://github.com/vntyper-org/MucOneUp/releases/download/testdata-v0.25.0/testdata_40-70_v0.25.0.tar.gz
curl -LO https://github.com/vntyper-org/MucOneUp/releases/download/testdata-v0.25.0/testdata_40-70_v0.25.0.tar.gz.sha256
sha256sum -c testdata_40-70_v0.25.0.tar.gz.sha256
tar -xzf testdata_40-70_v0.25.0.tar.gz
```

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
│   ├── normal/
│   │   ├── testdata_40-70.001.normal.simulated.fa
│   │   ├── testdata_40-70.001.normal.simulated.structure.tsv
│   │   └── testdata_40-70.001.normal.simulated_toxic.orf_stats.json
│   └── dupC/
│       ├── testdata_40-70.001.mut.simulated.fa
│       ├── testdata_40-70.001.mut.simulated.structure.tsv
│       └── testdata_40-70.001.mut.simulated_toxic.orf_stats.json
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

**Output**: `data/testdata_releases/{date}_{version}/testdata_40-70_{version}.tar.gz` plus a `.sha256` file. This directory is gitignored; tarballs are never committed.

**Publish** the tarball as a GitHub Release (requires an authenticated `gh` CLI):

```bash
make generate-testdata                  # writes data/testdata_releases/<date>_v<version>/
make publish-testdata VERSION=v0.48.0   # creates release testdata-v0.48.0 with tarball + .sha256
```

Then update the **Current dataset** link above.

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
