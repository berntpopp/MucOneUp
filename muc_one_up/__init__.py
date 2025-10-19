"""MucOneUp: MUC1 VNTR simulation and analysis toolkit.

MucOneUp is a Python toolkit for simulating realistic MUC1 Variable Number
Tandem Repeat (VNTR) diploid references with customizable mutations and
sequencing read simulation. Supports both Illumina and Oxford Nanopore
sequencing platforms.

Core Capabilities:
    - VNTR haplotype simulation with probability-based repeat transitions
    - Targeted mutation application (insert, delete, replace, delete_insert)
    - SNP integration for haplotype-specific variants
    - Illumina and ONT read simulation with realistic error profiles
    - ORF prediction and toxic protein detection
    - Comprehensive statistics and coverage analysis

Key Modules:
    simulate: Core haplotype generation engine
    mutate: Mutation application and validation
    config: Configuration loading and schema validation
    read_simulation: Illumina and ONT pipeline orchestration
    analysis: VNTR statistics and ORF analysis
    cli: Command-line interface (click-based)

Example:
    Basic simulation workflow::

        from muc_one_up.config import load_config
        from muc_one_up.simulate import simulate_diploid
        from muc_one_up.mutate import apply_mutations

        # Load configuration
        config = load_config("config.json")

        # Generate diploid haplotypes
        results = simulate_diploid(config, num_haplotypes=2, fixed_lengths=[50, 60])

        # Apply mutation
        mutated_results, units = apply_mutations(
            config, results, "dupC", targets=[(1, 25)]
        )

Command-Line Interface:
    Install and use via CLI::

        pip install .
        muconeup --config config.json simulate --out-base test --fixed-lengths 50

See Also:
    - CLAUDE.md: Comprehensive developer documentation
    - README.md: Installation and usage guide
    - config.json: Example configuration file
"""

from .version import __version__

__all__ = ["__version__"]
