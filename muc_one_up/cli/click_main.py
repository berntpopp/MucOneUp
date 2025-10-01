"""
Click-based CLI for MucOneUp.

Design Philosophy: Hybrid approach
- `simulate`: Main workflow command (integrated pipeline)
- `reads` and `analyze`: Standalone utilities

Follows SOLID, DRY, KISS principles by delegating to existing backend modules.
"""

import logging
from pathlib import Path
from typing import Any

import click

from ..exceptions import MucOneUpError
from ..version import __version__
from .config import (
    determine_simulation_mode,
    process_mutation_config,
    setup_configuration,
)
from .main import configure_logging
from .orchestration import run_single_simulation_iteration

# ============================================================================
# Root CLI Group
# ============================================================================


@click.group()
@click.version_option(version=__version__, prog_name="MucOneUp")
@click.option(
    "--config",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to JSON configuration file.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"]),
    default="INFO",
    help="Set logging level.",
)
@click.pass_context
def cli(ctx, config, log_level):
    """MucOneUp - MUC1 VNTR diploid reference simulator.

    \b
    Commands:
      simulate  Main workflow - generate haplotypes and run optional pipelines
      reads     Standalone read simulation utilities (Illumina, ONT)
      analyze   Standalone analysis utilities (ORFs, statistics)
    """
    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config
    ctx.obj["log_level"] = log_level
    configure_logging(log_level)


# ============================================================================
# SIMULATE Command - Main Workflow
# ============================================================================


@cli.command()
@click.option(
    "--out-base",
    default="muc1_simulated",
    show_default=True,
    help="Base name for all output files.",
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder for all files.",
)
@click.option(
    "--num-haplotypes",
    default=2,
    show_default=True,
    type=int,
    help="Number of haplotypes to simulate (typically 2 for diploid).",
)
@click.option(
    "--seed",
    type=int,
    help="Random seed for reproducible simulations.",
)
@click.option(
    "--reference-assembly",
    type=click.Choice(["hg19", "hg38"]),
    help="Reference genome assembly (overrides config file).",
)
@click.option(
    "--output-structure",
    is_flag=True,
    help="Write VNTR structure file.",
)
# Length Options (Mutually Exclusive)
@click.option(
    "--fixed-lengths",
    multiple=True,
    type=str,
    help="Fixed VNTR lengths or ranges (e.g., '60' or '20-40'). Repeat for multiple haplotypes.",
)
@click.option(
    "--input-structure",
    type=click.Path(exists=True, dir_okay=False),
    help="Path to file with predefined VNTR repeat compositions.",
)
@click.option(
    "--simulate-series",
    type=int,
    default=None,
    help="Generate series across fixed-length ranges with given step size.",
)
# Mutation Options
@click.option(
    "--mutation-name",
    help="Mutation to apply. Use 'normal,mutationName' for dual simulation.",
)
@click.option(
    "--mutation-targets",
    multiple=True,
    help="Mutation targets as 'haplotype_index,repeat_index' pairs (1-based).",
)
# SNP Options
@click.option(
    "--snp-input-file",
    type=click.Path(exists=True, dir_okay=False),
    help="TSV file with predefined SNPs.",
)
@click.option(
    "--random-snps",
    is_flag=True,
    help="Enable random SNP generation.",
)
@click.option(
    "--random-snp-density",
    type=float,
    help="SNP density per 1000 bp (requires --random-snps).",
)
@click.option(
    "--random-snp-output-file",
    type=str,
    help="Output file for generated random SNPs (requires --random-snps).",
)
@click.option(
    "--random-snp-region",
    type=click.Choice(["all", "constants_only", "vntr_only"]),
    default="constants_only",
    show_default=True,
    help="Region for random SNPs (requires --random-snps).",
)
@click.option(
    "--random-snp-haplotypes",
    type=click.Choice(["all", "1", "2"]),
    default="all",
    show_default=True,
    help="Haplotypes for random SNPs (requires --random-snps).",
)
# Pipeline Options (NEW: Integrated workflow)
@click.option(
    "--simulate-reads",
    type=click.Choice(["illumina", "ont"]),
    help="Run read simulation pipeline after haplotype generation.",
)
@click.option(
    "--output-orfs",
    is_flag=True,
    help="Run ORF prediction and toxic protein detection.",
)
@click.option(
    "--orf-min-aa",
    type=int,
    default=100,
    show_default=True,
    help="Minimum peptide length in amino acids (requires --output-orfs).",
)
@click.option(
    "--orf-aa-prefix",
    default=None,
    help="Filter ORFs by prefix, e.g., MTSSV (requires --output-orfs).",
)
@click.pass_context
def simulate(ctx, **kwargs):
    """Generate MUC1 VNTR diploid haplotypes with optional pipelines.

    \b
    This is the main workflow command. It can:
    - Generate haplotype sequences
    - Apply mutations
    - Integrate SNPs
    - Simulate reads (--simulate-reads)
    - Predict ORFs (--output-orfs)

    \b
    Examples:
      # Basic simulation
      muconeup --config config.json simulate --out-base output

      # Full pipeline in one command
      muconeup --config config.json simulate --out-base output \\
          --simulate-reads illumina --output-orfs

      # With mutations
      muconeup --config config.json simulate --out-base output \\
          --mutation-name dupC --mutation-targets 1,25 2,30
    """
    try:
        # Convert Click kwargs to argparse-style namespace for backend compatibility
        # This maintains DRY by reusing existing orchestration logic
        args = _make_args_namespace(ctx.obj["config_path"], kwargs)

        # Delegate to existing configuration and orchestration modules (SOLID, DRY)
        config, out_dir, out_base = setup_configuration(args)
        simulation_configs, predefined_chains, structure_mutation_info = determine_simulation_mode(
            args, config
        )
        dual_mutation_mode, mutation_pair, mutation_name = process_mutation_config(
            args, structure_mutation_info
        )

        # Run simulations using existing orchestration function
        for sim_index, fixed_conf in enumerate(simulation_configs, start=1):
            run_single_simulation_iteration(
                args,
                config,
                out_dir,
                out_base,
                sim_index,
                fixed_conf,
                predefined_chains,
                dual_mutation_mode,
                mutation_pair,
                structure_mutation_info,
            )

        logging.info("All simulations completed successfully.")
        ctx.exit(0)

    except KeyboardInterrupt:
        logging.warning("Interrupted by user (Ctrl+C)")
        ctx.exit(130)

    except MucOneUpError as e:
        logging.error("%s: %s", type(e).__name__, e)
        ctx.exit(1)

    except Exception as e:
        logging.exception("Unexpected error: %s", e)
        ctx.exit(2)


# ============================================================================
# READS Command Group - Standalone Utilities
# ============================================================================


@cli.group()
@click.pass_context
def reads(ctx):
    """Standalone read simulation utilities.

    \b
    These commands work with any FASTA file, not just MucOneUp outputs.
    Useful for processing external sequences or re-running read simulation.
    """
    pass


@reads.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    required=True,
    help="Base name for output files.",
)
@click.option(
    "--coverage",
    type=int,
    default=30,
    show_default=True,
    help="Target sequencing coverage.",
)
@click.option(
    "--threads",
    type=int,
    default=8,
    show_default=True,
    help="Number of threads for alignment.",
)
@click.pass_context
def illumina(ctx, input_fasta, out_dir, out_base, coverage, threads):
    """Simulate Illumina short reads from FASTA.

    \b
    Example:
      muconeup --config config.json reads illumina input.fa --out-base reads_output
    """
    try:
        # Load config from context
        import json

        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure read simulation
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "illumina"
        config["read_simulation"]["coverage"] = coverage
        config["read_simulation"]["threads"] = threads

        # Run pipeline
        logging.info("Starting Illumina read simulation for %s", input_fasta)
        simulate_reads_pipeline(config, input_fasta)
        logging.info("Illumina read simulation completed successfully.")
        ctx.exit(0)

    except Exception as e:
        logging.error("Read simulation failed: %s", e)
        ctx.exit(1)


@reads.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    required=True,
    help="Base name for output files.",
)
@click.option(
    "--coverage",
    type=int,
    default=30,
    show_default=True,
    help="Target sequencing coverage.",
)
@click.option(
    "--min-read-length",
    type=int,
    default=100,
    show_default=True,
    help="Minimum read length.",
)
@click.pass_context
def ont(ctx, input_fasta, out_dir, out_base, coverage, min_read_length):
    """Simulate Oxford Nanopore long reads from FASTA.

    \b
    Example:
      muconeup --config config.json reads ont input.fa --out-base reads_output
    """
    try:
        # Load config from context
        import json

        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure read simulation
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "ont"
        config["read_simulation"]["coverage"] = coverage
        config["nanosim_params"]["min_len"] = min_read_length

        # Run pipeline
        logging.info("Starting Oxford Nanopore read simulation for %s", input_fasta)
        simulate_reads_pipeline(config, input_fasta)
        logging.info("Oxford Nanopore read simulation completed successfully.")
        ctx.exit(0)

    except Exception as e:
        logging.error("Read simulation failed: %s", e)
        ctx.exit(1)


# ============================================================================
# ANALYZE Command Group - Standalone Utilities
# ============================================================================


@cli.group()
@click.pass_context
def analyze(ctx):
    """Standalone analysis utilities.

    \b
    These commands work with any FASTA file, not just MucOneUp outputs.
    Useful for analyzing external sequences or re-running analysis.
    """
    pass


@analyze.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    required=True,
    help="Base name for output files.",
)
@click.option(
    "--orf-min-aa",
    type=int,
    default=100,
    show_default=True,
    help="Minimum ORF length in amino acids.",
)
@click.option(
    "--orf-aa-prefix",
    default=None,
    help="Filter ORFs by prefix (e.g., MTSSV).",
)
@click.pass_context
def orfs(ctx, input_fasta, out_dir, out_base, orf_min_aa, orf_aa_prefix):
    """Predict ORFs and detect toxic protein features.

    \b
    Example:
      muconeup --config config.json analyze orfs input.fa --out-base analysis
    """
    # Standalone ORF analysis - delegate to orfipy command-line tool
    try:
        import subprocess

        # Use orfipy command-line tool directly for standalone FASTA analysis
        orf_output = Path(out_dir) / f"{out_base}.orfs.fa"
        logging.info("Running ORF prediction on %s", input_fasta)

        cmd = [
            "orfipy",
            input_fasta,
            "--pep",
            str(orf_output),
            "--min",
            str(orf_min_aa * 3),  # Convert AA to nucleotides
            "--start",
            "ATG",
        ]

        if orf_aa_prefix:
            cmd.extend(["--start-codon-prefix", orf_aa_prefix])

        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            logging.error("orfipy command failed: %s", result.stderr)
            ctx.exit(1)

        logging.info("ORF prediction completed: %s", orf_output)

        # Toxic protein detection (if orfipy succeeded and output exists)
        if orf_output.exists():
            import json

            from ..toxic_protein_detector import scan_orf_fasta

            config_path = Path(ctx.obj["config_path"])
            with config_path.open() as f:
                config = json.load(f)

            left_const = config.get("constants", {}).get("left")
            right_const = config.get("constants", {}).get("right")
            stats = scan_orf_fasta(str(orf_output), left_const=left_const, right_const=right_const)

            stats_file = Path(out_dir) / f"{out_base}.orf_stats.json"
            with stats_file.open("w") as f:
                json.dump(stats, f, indent=4)

            logging.info("Toxic protein detection stats written: %s", stats_file)

        ctx.exit(0)

    except Exception as e:
        logging.error("ORF analysis failed: %s", e)
        ctx.exit(1)


@analyze.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    required=True,
    help="Base name for output files.",
)
@click.pass_context
def stats(ctx, input_fasta, out_dir, out_base):
    """Generate basic sequence statistics for FASTA.

    \b
    Example:
      muconeup --config config.json analyze stats input.fa --out-base stats_output
    """
    # Simple sequence statistics for arbitrary FASTA files
    try:
        import json

        logging.info("Generating basic statistics for %s", input_fasta)

        # Read FASTA manually (simple parsing)
        sequences: list[tuple[str, str]] = []
        current_header: str | None = None
        current_seq: list[str] = []

        with Path(input_fasta).open() as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header is not None:
                        sequences.append((current_header, "".join(current_seq)))
                    current_header = line[1:]
                    current_seq = []
                elif line:
                    current_seq.append(line)

            if current_header is not None:
                sequences.append((current_header, "".join(current_seq)))

        # Compute stats
        stats_data: dict[str, Any] = {
            "input_file": str(input_fasta),
            "num_sequences": len(sequences),
            "sequences": [],
        }

        for header, sequence in sequences:
            gc_count = sequence.upper().count("G") + sequence.upper().count("C")
            gc_content = 100.0 * gc_count / len(sequence) if sequence else 0.0

            stats_data["sequences"].append(
                {"header": header, "length": len(sequence), "gc_content": round(gc_content, 2)}
            )

        stats_file = Path(out_dir) / f"{out_base}.basic_stats.json"
        with stats_file.open("w") as f:
            json.dump(stats_data, f, indent=4)

        logging.info("Statistics written: %s", stats_file)
        ctx.exit(0)

    except Exception as e:
        logging.error("Statistics generation failed: %s", e)
        ctx.exit(1)


# ============================================================================
# Helper Functions
# ============================================================================


def _make_args_namespace(config_path, kwargs):
    """
    Convert Click kwargs to argparse-style namespace.

    This maintains DRY principle by allowing reuse of existing backend functions
    that expect argparse namespace objects.

    Args:
        config_path: Path to configuration file
        kwargs: Dictionary of Click command options

    Returns:
        SimpleNamespace object compatible with existing backend functions
    """
    from types import SimpleNamespace

    # Map Click options to argparse-style attributes
    args = SimpleNamespace(
        config=config_path,
        out_base=kwargs["out_base"],
        out_dir=kwargs["out_dir"],
        num_haplotypes=kwargs["num_haplotypes"],
        seed=kwargs["seed"],
        reference_assembly=kwargs["reference_assembly"],
        output_structure=kwargs["output_structure"],
        # Length options
        fixed_lengths=list(kwargs["fixed_lengths"]) if kwargs["fixed_lengths"] else None,
        input_structure=kwargs["input_structure"],
        simulate_series=kwargs["simulate_series"],
        # Mutation options
        mutation_name=kwargs["mutation_name"],
        mutation_targets=list(kwargs["mutation_targets"]) if kwargs["mutation_targets"] else None,
        # SNP options
        snp_input_file=kwargs["snp_input_file"],
        random_snps=kwargs["random_snps"],
        random_snp_density=kwargs["random_snp_density"],
        random_snp_output_file=kwargs["random_snp_output_file"],
        random_snp_region=kwargs["random_snp_region"],
        random_snp_haplotypes=kwargs["random_snp_haplotypes"],
        # Pipeline options
        simulate_reads=kwargs["simulate_reads"],
        output_orfs=kwargs["output_orfs"],
        orf_min_aa=kwargs["orf_min_aa"],
        orf_aa_prefix=kwargs["orf_aa_prefix"],
    )

    return args


# ============================================================================
# Entry Point
# ============================================================================


def main():
    """Entry point for Click CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
