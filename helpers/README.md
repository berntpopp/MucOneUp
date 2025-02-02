# Helpers

This directory contains helper scripts for the muconeup pipeline. Each helper script automates a specific task required by the pipeline. The following sections describe each helper script and its usage.

---

## 1. Install References Helper

This helper script automates the download, extraction, and (optional) indexing of external reference files required for the muconeup pipeline. It is located in the `helpers` subdirectory of this repository.

### Overview

The script performs the following tasks:
- **Download**: Retrieves reference files from specified URLs.
- **Verification**: Calculates and logs MD5 checksums for downloaded files.
- **Extraction**: Automatically extracts files compressed with gzip.
- **Indexing**: Optionally indexes FASTA files using BWA (or a user-specified executable).
- **Reporting**: Generates a JSON configuration file with the absolute paths of the installed references.
- **Logging**: Logs all actions to both the console and a log file in the output directory.

### Requirements

- Python 3.6 or higher.
- Internet access to download the reference files.
- [BWA](http://bio-bwa.sourceforge.net/) (if indexing of FASTA files is required).

### Installation

No additional installation is required. All dependencies are part of the Python Standard Library.

### Usage

Run the script from the command line. The basic usage is as follows:

```bash
python helpers/install_references.py --output-dir /path/to/destination
```

#### Options

- `--output-dir`: **(Required)** Destination directory where the reference files will be installed.
- `--config`: *(Optional)* Path to a JSON configuration file. If not provided, the script will use the default configuration file (`install_references_config.json`) located alongside the script.
- `--skip-indexing`: *(Optional)* Skip the indexing of FASTA reference files (even if an index command is provided in the config).

#### Example

Download and install references to the `./references` directory, using the default configuration, and skip the indexing step:

```bash
python helpers/install_references.py --output-dir ./references --skip-indexing
```

### Configuration

The script uses a JSON configuration file to determine which references to download. Below is an example configuration (`install_references_config.json`):

```json
{
  "bwa_path": "bwa",
  "references": {
    "GRCh38_no_alt_analysis_set": {
      "url": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
      "target_path": "GRCh38_no_alt_analysis_set.fna.gz",
      "index_command": "bwa index {path}"
    },
    "Hs-Nova-TruSeq": {
      "url": "https://github.com/schmeing/ReSeq-profiles/raw/master/profiles/Hs-Nova-TruSeq.reseq",
      "target_path": "Hs-Nova-TruSeq.reseq"
    }
  }
}
```

#### Configuration Fields

- **`bwa_path`**: The path or name of the BWA executable. This will replace `"bwa"` in the indexing command.
- **`references`**: A mapping of reference names to their details:
  - **`url`**: The URL from which to download the reference file.
  - **`target_path`**: The relative path (from the output directory) where the file will be saved.
  - **`index_command`**: *(Optional)* The command to index the reference file. Use `{path}` as a placeholder for the file path.

### Logging and Output Files

- **Log File**: A log file named `install_references.log` will be created in the output directory.
- **Installed References Configuration**: After running, an `installed_references.json` file will be generated in the output directory. This file maps each reference name to its absolute installed file path.

---

## 2. BAM Anonymizer Helper

The `bam_anonymizer.py` script is designed to prepare sample data for read simulation testing by processing BAM files. It performs the following steps:

- **Subsetting**: Extracts a specific genomic region (default: the hg38 MUC1 VNTR region, `chr1:155184000-155194000`) from an input BAM file.
- **Replacing Read Groups**: Uses GATK AddOrReplaceReadGroups to replace the original read group information with generic values.
- **Anonymization**: Anonymizes the BAM file using GATK's ReadAnonymizer.
- **Reheadering**: Uses `samtools reheader` to remove unwanted header lines (such as @PG, @RG, and @CO) from the anonymized BAM.
- **Output**: Produces a final output BAM file named exactly according to the provided target design (e.g. `twist_v2.bam`).
- **Pseudonymisation Report**: Computes MD5 checksums for the original and final BAM files and writes a pseudonymisation output CSV (named with the target design as a prefix) with the following columns:
  - `old_name`
  - `old_md5`
  - `new_name`
  - `new_md5`

### Requirements

- Python 3.6 or higher.
- External tools:
  - `samtools`
  - `gatk`
- A reference FASTA file is required by GATK.

### Usage

Run the script from the command line. The basic usage is as follows:

```bash
python helpers/bam_anonymizer.py --input-bam /path/to/input.bam --target-design twist_v2 --ref /path/to/reference.fasta
```

#### Options

- `--input-bam`: **(Required)** Path to the input BAM file.
- `--target-design`: **(Required)** Target design name to be used for naming the output file (e.g. twist_v2).
- `--ref`: **(Required)** Path to the reference FASTA file required by GATK.
- `--region`: *(Optional)* Genomic region to subset from the BAM file. Default is `chr1:155184000-155194000` (hg38 MUC1 VNTR region).
- `--output-dir`: *(Optional)* Directory to write output files. Default is the current directory.
- `--keep-intermediates`: *(Optional)* Retain intermediate files (for debugging purposes). If not set, all intermediate files and their index files are deleted.
- `--log-level`: *(Optional)* Logging level (DEBUG, INFO, WARNING, ERROR). Default is INFO.

#### Example

```bash
python helpers/bam_anonymizer.py --input-bam data/sample.bam --target-design twist_v2 --ref reference/GRCh38.fasta --output-dir data --log-level DEBUG
```

### Output

After processing, the script produces:
- A final anonymized and reheadered BAM file named `twist_v2.bam` in the output directory.
- A pseudonymisation CSV file named `twist_v2_pseudonymisation_output.csv` that contains the following columns:
  - `old_name` (original BAM file name)
  - `old_md5` (MD5 checksum of the original BAM)
  - `new_name` (final BAM file name)
  - `new_md5` (MD5 checksum of the final BAM)

### Logging and Cleanup

- **Logging**: The script logs all processing steps. Use `--log-level` to set the desired verbosity.
- **Cleanup**: If `--keep-intermediates` is not set, all intermediate files (and their index files) are removed after processing.
