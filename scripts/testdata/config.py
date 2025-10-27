"""Configuration management for test data generation."""

from dataclasses import dataclass, field


@dataclass
class SampleConfig:
    """Configuration for a single test sample."""

    base_name: str = "testdata_40-70"
    vntr_lengths: str = "40,70"
    mutation_name: str = "normal,dupC"
    mutation_targets: str = "1,20"
    seed_base: int = 42000  # Base seed, incremented per platform


@dataclass
class PlatformConfig:
    """Configuration for a sequencing platform."""

    name: str
    enabled: bool
    coverage: int
    seed_offset: int  # Added to base seed
    output_subdir: str

    def get_seed(self, base_seed: int) -> int:
        """Calculate platform-specific seed."""
        return base_seed + self.seed_offset


@dataclass
class DatasetConfig:
    """Complete dataset generation configuration."""

    version: str
    sample: SampleConfig = field(default_factory=SampleConfig)

    platforms: dict[str, PlatformConfig] = field(
        default_factory=lambda: {
            "illumina": PlatformConfig(
                name="illumina",
                enabled=True,
                coverage=50,
                seed_offset=0,
                output_subdir="illumina",
            ),
            "ont": PlatformConfig(
                name="ont",
                enabled=True,
                coverage=30,
                seed_offset=10,
                output_subdir="ont",
            ),
            "pacbio": PlatformConfig(
                name="pacbio",
                enabled=True,
                coverage=30,
                seed_offset=20,
                output_subdir="pacbio",
            ),
        }
    )

    threads: int = 4
    verbose: bool = False

    @property
    def enabled_platforms(self) -> list[PlatformConfig]:
        """Get list of enabled platforms."""
        return [p for p in self.platforms.values() if p.enabled]


def load_config(version: str, **overrides) -> DatasetConfig:
    """Load configuration with optional overrides."""
    config = DatasetConfig(version=version)

    # Apply overrides
    for key, value in overrides.items():
        if hasattr(config, key):
            setattr(config, key, value)

    return config
