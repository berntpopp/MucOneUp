# Install References Helper

This helper script automates the download, extraction, and (optional) indexing of external reference files required for the muconeup pipeline. It is located in the `helpers` subdirectory of this repository.

## Overview

The script performs the following tasks:
- **Download**: Retrieves reference files from specified URLs.
- **Verification**: Calculates and logs MD5 checksums for downloaded files.
- **Extraction**: Automatically extracts files compressed with gzip.
- **Indexing**: Optionally indexes FASTA files using BWA (or a user-specified executable).
- **Reporting**: Generates a JSON configuration file with the absolute paths of the installed references.
- **Logging**: Logs all actions to both the console and a log file in the output directory.

## Requirements

- Python 3.6 or higher.
- Internet access to download the reference files.
- [BWA](http://bio-bwa.sourceforge.net/) (if indexing of FASTA files is required).

## Installation

No additional installation is required. All dependencies are part of the Python Standard Library.

## Usage

Run the script from the command line. The basic usage is as follows:

```bash
python helpers/install_references.py --output-dir /path/to/destination
```

### Options

- `--output-dir`: **(Required)** Destination directory where the reference files will be installed.
- `--config`: *(Optional)* Path to a JSON configuration file. If not provided, the script will use the default configuration file (`install_references_config.json`) located alongside the script.
- `--skip-indexing`: *(Optional)* Skip the indexing of FASTA reference files (even if an index command is provided in the config).

### Example

Download and install references to the `./references` directory, using the default configuration, and skip the indexing step:

```bash
python helpers/install_references.py --output-dir ./references --skip-indexing
```

## Configuration

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

### Configuration Fields

- **`bwa_path`**: The path or name of the BWA executable. This will replace `"bwa"` in the indexing command.
- **`references`**: A mapping of reference names to their details:
  - **`url`**: The URL from which to download the reference file.
  - **`target_path`**: The relative path (from the output directory) where the file will be saved.
  - **`index_command`**: *(Optional)* The command to index the reference file. Use `{path}` as a placeholder for the file path.

## Logging and Output Files

- **Log File**: A log file named `install_references.log` will be created in the output directory.
- **Installed References Configuration**: After running, an `installed_references.json` file will be generated in the output directory. This file maps each reference name to its absolute installed file path.
