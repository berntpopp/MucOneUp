"""Metadata and checksum generation for test datasets."""

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

from .config import DatasetConfig


class MetadataGenerator:
    """Generates manifests, checksums, and documentation."""

    def __init__(self, config: DatasetConfig):
        """
        Initialize metadata generator.

        Args:
            config: Dataset configuration
        """
        self.config = config
        self.logger = logging.getLogger(__name__)

    def generate_manifest(self, output_dir: Path, all_files: list[Path]) -> Path:
        """
        Generate manifest.json with complete metadata.

        Args:
            output_dir: Dataset root directory
            all_files: List of all files in dataset

        Returns:
            Path to created manifest.json
        """
        self.logger.info("Generating manifest.json")

        manifest: dict[str, Any] = {
            "version": self.config.version,
            "generated_with": "MucOneUp",
            "muconeup_version": self.config.version,
            "generation_date": datetime.utcnow().isoformat() + "Z",
            "sample": {
                "base_name": self.config.sample.base_name,
                "vntr_lengths": self.config.sample.vntr_lengths,
                "vntr_h1": 40,
                "vntr_h2": 70,
                "mutation_mode": self.config.sample.mutation_name,
                "mutation_targets": self.config.sample.mutation_targets,
                "seed_base": self.config.sample.seed_base,
            },
            "platforms": {},
            "structure": self._describe_structure(),
            "checksums": {},
            "file_count": len(all_files),
        }

        # Platform metadata
        for platform in self.config.enabled_platforms:
            manifest["platforms"][platform.name] = {
                "enabled": platform.enabled,
                "coverage": platform.coverage,
                "seed": platform.get_seed(self.config.sample.seed_base),
            }

        # Calculate checksums for all files
        for file_path in all_files:
            rel_path = file_path.relative_to(output_dir)
            checksum = self._calculate_checksum(file_path)
            manifest["checksums"][str(rel_path)] = checksum

        manifest_path = output_dir / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)

        self.logger.info(f"✓ Manifest created: {len(manifest['checksums'])} files")

        return manifest_path

    def generate_checksums_file(self, output_dir: Path, manifest: dict[str, Any]) -> Path:
        """
        Generate checksums.txt in sha256sum format.

        Args:
            output_dir: Dataset root directory
            manifest: Manifest dictionary containing checksums

        Returns:
            Path to created checksums.txt
        """
        self.logger.info("Generating checksums.txt")

        checksums_path = output_dir / "checksums.txt"

        with open(checksums_path, "w") as f:
            for rel_path, checksum in sorted(manifest["checksums"].items()):
                f.write(f"{checksum}  {rel_path}\n")

        self.logger.info(f"✓ Checksums file created: {len(manifest['checksums'])} entries")

        return checksums_path

    def generate_readme(self, output_dir: Path) -> Path:
        """
        Generate README.md for the dataset.

        Args:
            output_dir: Dataset root directory

        Returns:
            Path to created README.md
        """
        self.logger.info("Generating README.md")

        readme_content = self._build_readme()

        readme_path = output_dir / "README.md"
        with open(readme_path, "w") as f:
            f.write(readme_content)

        self.logger.info("✓ README created")

        return readme_path

    @staticmethod
    def _calculate_checksum(file_path: Path) -> str:
        """
        Calculate SHA256 checksum for a file.

        Args:
            file_path: Path to file

        Returns:
            Hexadecimal checksum string
        """
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()

    @staticmethod
    def _describe_structure() -> dict[str, str]:
        """
        Describe dataset directory structure.

        Returns:
            Dict mapping directories to descriptions
        """
        return {
            "references/": "Shared diploid references (normal + dupC)",
            "illumina/": "Illumina paired-end sequencing (50x)",
            "ont/": "Oxford Nanopore sequencing (30x)",
            "pacbio/": "PacBio HiFi sequencing (30x)",
        }

    def _build_readme(self) -> str:
        """
        Build README.md content.

        Returns:
            Complete README markdown string
        """
        version = self.config.version

        return f"""# MucOneUp Test Dataset {version}

## Overview

Comprehensive test dataset with **asymmetric VNTR (40/70 repeats)** across three sequencing platforms.

**Generated with**: MucOneUp {version}
**Platforms**: Illumina (50×), ONT (30×), PacBio (30×)
**Sample**: Asymmetric diploid (h1=40, h2=70 repeats)
**Variants**: Normal baseline + dupC mutation

## Structure

```
testdata_40-70_{version}/
├── references/          # Shared diploid references
│   ├── normal/          # Normal baseline
│   └── dupC/            # dupC mutation (h1:r20)
├── illumina/            # Illumina paired-end (50×)
│   ├── normal/
│   └── dupC/
├── ont/                 # Oxford Nanopore (30×)
│   ├── normal/
│   └── dupC/
└── pacbio/              # PacBio HiFi (30×)
    ├── normal/
    └── dupC/
```

## Quick Start

```bash
# Verify integrity
sha256sum -c checksums.txt

# Explore references
cat references/normal/testdata_40-70.001.normal.simulated.structure.tsv
cat references/dupC/testdata_40-70.001.dupC.simulated.structure.tsv

# View Illumina reads
zcat illumina/normal/reads_R1.fastq.gz | head -4

# View ONT reads
zcat ont/normal/reads.fastq.gz | head -4

# View PacBio reads
zcat pacbio/normal/reads.fastq.gz | head -4
```

## Platform Comparison

| Platform | Coverage | Read Type | VNTR-Biased |
|----------|----------|-----------|-------------|
| Illumina | 50× | Paired-end (2×150bp) | ✅ Yes |
| ONT | 30× | Long reads (~3kb) | ✅ Yes |
| PacBio | 30× | HiFi reads (~15kb) | ✅ Yes |

## Reproducibility

```bash
# Generate references (once)
muconeup --config config.json simulate \\
  --out-base testdata_40-70 \\
  --fixed-lengths 40,70 \\
  --mutation-name normal,dupC \\
  --mutation-targets 1,20 \\
  --seed 42000

# Simulate Illumina (seed 42000)
muconeup --config config.json reads illumina \\
  testdata_40-70.001.normal.simulated.fa --seed 42000
muconeup --config config.json reads illumina \\
  testdata_40-70.001.dupC.simulated.fa --seed 42000

# Simulate ONT (seed 42010)
muconeup --config config.json reads ont \\
  testdata_40-70.001.normal.simulated.fa --seed 42010
muconeup --config config.json reads ont \\
  testdata_40-70.001.dupC.simulated.fa --seed 42010

# Simulate PacBio (seed 42020)
muconeup --config config.json reads pacbio \\
  testdata_40-70.001.normal.simulated.fa --seed 42020
muconeup --config config.json reads pacbio \\
  testdata_40-70.001.dupC.simulated.fa --seed 42020
```

## Citation

```bibtex
@software{{muconeup_{version.replace(".", "_").replace("-", "_")},
  author = {{Popp, Bernt}},
  title = {{MucOneUp: MUC1 VNTR Simulation Toolkit}},
  version = {{{version}}},
  year = {{2025}},
  url = {{https://github.com/berntpopp/MucOneUp}}
}}
```

## License

GPL-3.0 (same as MucOneUp)

## Support

- Documentation: https://berntpopp.github.io/MucOneUp/
- Issues: https://github.com/berntpopp/MucOneUp/issues
"""
