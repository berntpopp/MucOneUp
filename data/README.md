# Data Folder

This folder contains anonymized BAM files for the MUC1 region generated for various exome designs used as templates for read simulation with [wessim](https://github.com/bioinformatics-centre/wessim). These files have been processed to include only the MUC1 region and have been anonymized for data privacy.

**Reference Assemblies**: All BAM files in this directory are currently mapped to the GRCh38/hg38 reference genome. Files for the GRCh37/hg19 reference assembly will be added in the future.

## Contents

### hg38 (GRCh38) BAM Files

- **`idt_v1_hg38.bam`**: BAM file for the IDT v1 exome design (MUC1 region, hg38).
- **`idt_v1_hg38.bam.bai`**: Index file for `idt_v1_hg38.bam`.
- **`twist_v2_hg38.bam`**: BAM file for the Twist v2 exome design (MUC1 region, hg38).
- **`twist_v2_hg38.bam.bai`**: Index file for `twist_v2_hg38.bam`.
- **`agilent_v8_hg38.bam`**: BAM file for the Agilent SureSelect v8 exome design (MUC1 region, hg38).
- **`agilent_v8_hg38.bam.bai`**: Index file for `agilent_v8_hg38.bam`.

### hg19 (GRCh37) BAM Files

_Coming soon. hg19 versions will be added in the future._

## Usage

These files are intended to be used as input templates for read simulation with wessim. They provide a realistic representation of the MUC1 region as captured by different exome designs:

- **Twist v2**: Represents a widely used exome capture design capturing the MUC1 VNTR region.
- **IDT v1**: Represents one alternative exome capture design capturing the MUC1 VNTR region.
- **Agilent SureSelect v8**: Represents one alternative exome capture design commonly used in different versions. Agilent SureSelct designs typically do not capture the MUC1-VNTR region to sufficient coverage.

To use these files in your read simulation pipeline, point your simulation tool (e.g., wessim) to the desired BAM file and its index.

## Additional Information

- **Anonymization**: The BAM files have been anonymized to protect sensitive sample information.
- **Region**: Only the MUC1 region is included in these BAM files.
- **Index Files**: The accompanying `.bai` index files are required for most alignment and simulation tools.
