"""Sample generation logic for test data."""

import logging
from pathlib import Path

from .config import DatasetConfig, PlatformConfig
from .utils import ensure_dir, run_command


class ReferenceGenerator:
    """Generates diploid references using dual simulation mode."""

    def __init__(self, config: DatasetConfig, muconeup_config: Path):
        """
        Initialize reference generator.

        Args:
            config: Dataset configuration
            muconeup_config: Path to MucOneUp config.json file
        """
        self.config = config
        self.muconeup_config = muconeup_config
        self.logger = logging.getLogger(__name__)

    def generate(self, output_dir: Path) -> dict[str, Path]:
        """
        Generate diploid references for normal and dupC variants.

        Uses dual simulation mode (--mutation-name normal,dupC) to generate
        both variants in a single run for efficiency.

        Args:
            output_dir: Base output directory

        Returns:
            Dict mapping variant name to FASTA path
        """
        self.logger.info("Generating diploid references (dual simulation mode)")

        refs_dir = output_dir / "references"
        ensure_dir(refs_dir)

        sample = self.config.sample

        # Single command generates both variants
        cmd = [
            "muconeup",
            "--config",
            str(self.muconeup_config),
            "simulate",
            "--out-base",
            str(refs_dir / sample.base_name),
            "--fixed-lengths",
            sample.vntr_lengths,
            "--mutation-name",
            sample.mutation_name,
            "--mutation-targets",
            sample.mutation_targets,
            "--seed",
            str(sample.seed_base),
        ]

        run_command(cmd, "Generate diploid references")

        # Organize files by variant
        references = {}

        for variant in ["normal", "dupC"]:
            variant_dir = refs_dir / variant
            ensure_dir(variant_dir)

            # Move files to variant subdirectory
            pattern = f"{sample.base_name}.001.{variant}.*"
            for file in refs_dir.glob(pattern):
                if file.is_file():
                    dest = variant_dir / file.name
                    file.rename(dest)
                    self.logger.debug(f"  Moved: {file.name} -> {variant}/")

                    # Track FASTA file
                    if file.suffix == ".fa":
                        references[variant] = dest

        if len(references) != 2:
            raise RuntimeError(
                f"Expected 2 reference FASTAs, found {len(references)}: {list(references.keys())}"
            )

        self.logger.info(f"✓ Generated references: {list(references.keys())}")

        return references


class ReadSimulator:
    """Simulates sequencing reads for multiple platforms."""

    def __init__(self, config: DatasetConfig, muconeup_config: Path):
        """
        Initialize read simulator.

        Args:
            config: Dataset configuration
            muconeup_config: Path to MucOneUp config.json file
        """
        self.config = config
        self.muconeup_config = muconeup_config
        self.logger = logging.getLogger(__name__)

    def simulate_platform(
        self, platform: PlatformConfig, references: dict[str, Path], output_dir: Path
    ) -> dict[str, list[Path]]:
        """
        Simulate reads for a platform across all variants.

        Args:
            platform: Platform configuration
            references: Dict mapping variant name to reference FASTA
            output_dir: Base output directory

        Returns:
            Dict mapping variant name to list of generated files
        """
        self.logger.info(
            f"Simulating {platform.name.upper()} reads (coverage={platform.coverage}×)"
        )

        platform_dir = output_dir / platform.output_subdir
        ensure_dir(platform_dir)

        generated_files = {}
        seed = platform.get_seed(self.config.sample.seed_base)

        for variant, ref_fa in references.items():
            variant_dir = platform_dir / variant
            ensure_dir(variant_dir)

            self.logger.info(f"  Variant: {variant} (seed={seed})")

            cmd = [
                "muconeup",
                "--config",
                str(self.muconeup_config),
                "reads",
                platform.name,
                str(ref_fa),
                "--seed",
                str(seed),
                "--threads",
                str(self.config.threads),
            ]

            # Platform-specific parameters
            if platform.name in ["ont", "pacbio"]:
                cmd.extend(["--coverage", str(platform.coverage)])

            run_command(cmd, f"Simulate {platform.name} reads for {variant}")

            # Collect and organize generated files
            files = self._collect_and_organize_files(ref_fa, variant_dir, platform.name)

            generated_files[variant] = files
            self.logger.info(f"    ✓ Generated {len(files)} files")

        return generated_files

    def _collect_and_organize_files(
        self, ref_fa: Path, variant_dir: Path, platform: str
    ) -> list[Path]:
        """
        Collect generated read files and organize into directories.

        Args:
            ref_fa: Reference FASTA path
            variant_dir: Target directory for this variant
            platform: Platform name (illumina, ont, pacbio)

        Returns:
            List of organized file paths
        """
        # Extract base name from reference FASTA
        # e.g., "testdata_40-70.001.normal.simulated.fa" -> "testdata_40-70.001.normal"
        base = ref_fa.stem.replace(".simulated", "")
        parent_dir = ref_fa.parent

        files = []

        # Platform-specific file patterns and rename mappings
        if platform == "illumina":
            patterns = [
                f"{base}.simulated_R1.fastq.gz",
                f"{base}.simulated_R2.fastq.gz",
                f"{base}.simulated.bam",
                f"{base}.simulated.bam.bai",
                f"{base}.simulated_vntr_biased.bam",
                f"{base}.simulated_vntr_biased.bam.bai",
                f"{base}.simulated_vntr_biased_R1.fastq.gz",
                f"{base}.simulated_vntr_biased_R2.fastq.gz",
            ]
            rename_map = {
                f"{base}.simulated_R1.fastq.gz": "reads_R1.fastq.gz",
                f"{base}.simulated_R2.fastq.gz": "reads_R2.fastq.gz",
                f"{base}.simulated.bam": "aligned.bam",
                f"{base}.simulated.bam.bai": "aligned.bam.bai",
                f"{base}.simulated_vntr_biased.bam": "vntr_biased.bam",
                f"{base}.simulated_vntr_biased.bam.bai": "vntr_biased.bam.bai",
                f"{base}.simulated_vntr_biased_R1.fastq.gz": "vntr_biased_R1.fastq.gz",
                f"{base}.simulated_vntr_biased_R2.fastq.gz": "vntr_biased_R2.fastq.gz",
            }
        elif platform in ["ont", "pacbio"]:
            patterns = [
                f"{base}.simulated.fastq.gz",
                f"{base}.simulated.bam",
                f"{base}.simulated.bam.bai",
                f"{base}.simulated_vntr_biased.bam",
                f"{base}.simulated_vntr_biased.bam.bai",
                f"{base}.simulated_vntr_biased.fastq.gz",
            ]
            rename_map = {
                f"{base}.simulated.fastq.gz": "reads.fastq.gz",
                f"{base}.simulated.bam": "aligned.bam",
                f"{base}.simulated.bam.bai": "aligned.bam.bai",
                f"{base}.simulated_vntr_biased.bam": "vntr_biased.bam",
                f"{base}.simulated_vntr_biased.bam.bai": "vntr_biased.bam.bai",
                f"{base}.simulated_vntr_biased.fastq.gz": "vntr_biased.fastq.gz",
            }
        else:
            raise ValueError(f"Unknown platform: {platform}")

        # Find and move files
        for pattern in patterns:
            source = parent_dir / pattern
            if source.exists():
                dest_name = rename_map.get(pattern, pattern)
                dest = variant_dir / dest_name
                source.rename(dest)
                files.append(dest)
                self.logger.debug(f"      Moved: {source.name} -> {dest_name}")
            else:
                self.logger.warning(f"      Expected file not found: {pattern}")

        return files
