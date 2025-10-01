"""
Click-based CLI for MucOneUp.

Design Philosophy: Clean Separation (Unix Philosophy)
Each command does ONE thing well:
- `simulate`: ONLY generates haplotypes
- `reads`: ONLY simulates reads from FASTA
- `analyze`: ONLY analyzes FASTA
- `pipeline`: Orchestrates the above for convenience

Follows SOLID principles: Single Responsibility, Dependency Inversion.
"""

import json
import logging
import subprocess
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
from .orchestration import run_single_simulation_iteration

# ============================================================================
# Logging Configuration
# ============================================================================


def configure_logging(level_str: str) -> None:
    """Configure root logging based on the provided level string.

    If level_str is 'NONE', disable logging.
    """
    if level_str.upper() == "NONE":
        logging.disable(logging.CRITICAL)
    else:
        level = getattr(logging, level_str.upper(), logging.INFO)
        root_logger = logging.getLogger()
        if root_logger.handlers:
            for handler in root_logger.handlers[:]:
                root_logger.removeHandler(handler)
        logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")
        logging.info("Logging configured at level: %s", level_str.upper())


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
    Philosophy: Each command does ONE thing (Unix philosophy)
      simulate  - Generate haplotypes ONLY
      reads     - Simulate reads from FASTA
      analyze   - Analyze FASTA (ORFs, stats)
      pipeline  - Run complete workflow (convenience)

    \b
    Examples:
      # Clean separation - compose commands
      muconeup --config X simulate --out-base Y
      muconeup --config X reads illumina Y.001.simulated.fa
      muconeup --config X analyze orfs Y.001.simulated.fa

      # Or use pipeline for convenience
      muconeup --config X pipeline --out-base Y --with-reads --with-orfs
    """
    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config
    ctx.obj["log_level"] = log_level
    configure_logging(log_level)


# ============================================================================
# SIMULATE Command - ONLY Haplotype Generation
# ============================================================================


@cli.command()
@click.option(
    "--out-base",
    default="muc1_simulated",
    show_default=True,
    help="Base name for output files.",
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--num-haplotypes",
    default=2,
    show_default=True,
    type=int,
    help="Number of haplotypes to simulate.",
)
@click.option(
    "--seed",
    type=int,
    help="Random seed for reproducibility.",
)
@click.option(
    "--reference-assembly",
    type=click.Choice(["hg19", "hg38"]),
    help="Reference assembly (overrides config).",
)
@click.option(
    "--output-structure",
    is_flag=True,
    help="Write VNTR structure file.",
)
# Length Options
@click.option(
    "--fixed-lengths",
    multiple=True,
    type=str,
    help="Fixed VNTR lengths or ranges (e.g., '60' or '20-40').",
)
@click.option(
    "--input-structure",
    type=click.Path(exists=True, dir_okay=False),
    help="Predefined VNTR structure file.",
)
@click.option(
    "--simulate-series",
    type=int,
    help="Series step size for fixed-length ranges.",
)
# Mutation Options
@click.option(
    "--mutation-name",
    help="Mutation name. Use 'normal,mutation' for dual simulation.",
)
@click.option(
    "--mutation-targets",
    multiple=True,
    help="Mutation targets as 'hap_idx,rep_idx' pairs (1-based).",
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
    help="SNP density per 1000 bp.",
)
@click.option(
    "--random-snp-output-file",
    type=str,
    help="Output file for random SNPs.",
)
@click.option(
    "--random-snp-region",
    type=click.Choice(["all", "constants_only", "vntr_only"]),
    default="constants_only",
    show_default=True,
    help="Region for random SNPs.",
)
@click.option(
    "--random-snp-haplotypes",
    type=click.Choice(["all", "1", "2"]),
    default="all",
    show_default=True,
    help="Haplotypes for random SNPs.",
)
@click.pass_context
def simulate(ctx, **kwargs):
    """Generate MUC1 VNTR diploid haplotypes.

    \b
    Single Responsibility: ONLY generates haplotype FASTA files.
    Does NOT run read simulation or ORF prediction.
    Use 'pipeline' command or chain commands manually for full workflow.

    \b
    Output:
      - {out_base}.{iteration}.simulated.fa (haplotype sequences)
      - {out_base}.{iteration}.vntr_structure.txt (if --output-structure)
      - {out_base}.{iteration}.simulation_stats.json (statistics)

    \b
    Example:
      muconeup --config config.json simulate --out-base output
    """
    try:
        # Convert Click kwargs to argparse-style namespace (DRY - reuse backend)
        args = _make_args_namespace(ctx.obj["config_path"], kwargs)

        # IMPORTANT: Disable pipeline options (simulate is PURE)
        args.simulate_reads = None
        args.output_orfs = False

        # Delegate to existing orchestration (SOLID - Dependency Inversion)
        config, out_dir, out_base = setup_configuration(args)
        simulation_configs, predefined_chains, structure_mutation_info = determine_simulation_mode(
            args, config
        )
        dual_mutation_mode, mutation_pair, _ = process_mutation_config(
            args, structure_mutation_info
        )

        # Run haplotype generation ONLY
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

        logging.info("Haplotype generation completed successfully.")
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
# READS Command Group - Pure Read Simulation Utilities
# ============================================================================


@cli.group()
def reads():
    """Read simulation utilities.

    \b
    Single Responsibility: Simulate reads from ANY FASTA file.
    Works with MucOneUp outputs or external sequences.
    """


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
    help="Number of threads.",
)
@click.pass_context
def illumina(ctx, input_fasta, out_dir, out_base, coverage, threads):
    """Simulate Illumina short reads from FASTA.

    \b
    Example:
      muconeup --config X reads illumina output.001.simulated.fa --out-base reads_out
    """
    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config
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
        logging.info("Simulating Illumina reads for %s", input_fasta)
        simulate_reads_pipeline(config, input_fasta)
        logging.info("Illumina read simulation completed.")
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
    help="Target coverage.",
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
      muconeup --config X reads ont output.001.simulated.fa --out-base reads_out
    """
    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure ONT simulation
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "ont"
        config["read_simulation"]["coverage"] = coverage
        if "nanosim_params" not in config:
            config["nanosim_params"] = {}
        config["nanosim_params"]["min_len"] = min_read_length

        # Run pipeline
        logging.info("Simulating ONT reads for %s", input_fasta)
        simulate_reads_pipeline(config, input_fasta)
        logging.info("ONT read simulation completed.")
        ctx.exit(0)

    except Exception as e:
        logging.error("Read simulation failed: %s", e)
        ctx.exit(1)


# ============================================================================
# ANALYZE Command Group - Pure Analysis Utilities
# ============================================================================


@cli.group()
def analyze():
    """Analysis utilities.

    \b
    Single Responsibility: Analyze ANY FASTA file.
    Works with MucOneUp outputs or external sequences.
    """


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
      muconeup --config X analyze orfs output.001.simulated.fa --out-base analysis
    """
    try:
        # Use orfipy command-line tool
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
            logging.error("orfipy failed: %s", result.stderr)
            ctx.exit(1)

        logging.info("ORF prediction completed: %s", orf_output)

        # Toxic protein detection
        if orf_output.exists():
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

            logging.info("Toxic protein stats written: %s", stats_file)

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
    """Generate basic sequence statistics.

    \b
    Example:
      muconeup --config X analyze stats output.001.simulated.fa --out-base stats_out
    """
    try:
        logging.info("Generating statistics for %s", input_fasta)

        # Simple FASTA parsing
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
# PIPELINE Command - Convenience Orchestrator
# ============================================================================


@cli.command()
@click.option(
    "--out-base",
    default="muc1_simulated",
    show_default=True,
    help="Base name for all outputs.",
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--with-reads",
    type=click.Choice(["illumina", "ont"]),
    help="Include read simulation step.",
)
@click.option(
    "--with-orfs",
    is_flag=True,
    help="Include ORF prediction step.",
)
# Inherit simulation options
@click.option("--num-haplotypes", default=2, show_default=True, type=int)
@click.option("--seed", type=int)
@click.option("--reference-assembly", type=click.Choice(["hg19", "hg38"]))
@click.option("--output-structure", is_flag=True)
@click.option("--fixed-lengths", multiple=True, type=str)
@click.option("--input-structure", type=click.Path(exists=True, dir_okay=False))
@click.option("--simulate-series", type=int)
@click.option("--mutation-name")
@click.option("--mutation-targets", multiple=True)
@click.option("--snp-input-file", type=click.Path(exists=True, dir_okay=False))
@click.option("--random-snps", is_flag=True)
@click.option("--random-snp-density", type=float)
@click.option("--random-snp-output-file", type=str)
@click.option(
    "--random-snp-region",
    type=click.Choice(["all", "constants_only", "vntr_only"]),
    default="constants_only",
)
@click.option("--random-snp-haplotypes", type=click.Choice(["all", "1", "2"]), default="all")
# Read/ORF options
@click.option("--coverage", type=int, default=30, show_default=True)
@click.option("--threads", type=int, default=8, show_default=True)
@click.option("--orf-min-aa", type=int, default=100, show_default=True)
@click.pass_context
def pipeline(ctx, with_reads, with_orfs, **kwargs):
    """Run complete workflow (convenience command).

    \b
    Orchestrates: simulate → reads → analyze
    Use this for quick full pipelines, or chain individual commands for flexibility.

    \b
    Example:
      muconeup --config X pipeline --out-base Y --with-reads illumina --with-orfs
    """
    try:
        # Step 1: Run simulate command
        logging.info("Pipeline Step 1/3: Generating haplotypes...")
        ctx.invoke(
            simulate,
            out_base=kwargs["out_base"],
            out_dir=kwargs["out_dir"],
            num_haplotypes=kwargs["num_haplotypes"],
            seed=kwargs["seed"],
            reference_assembly=kwargs["reference_assembly"],
            output_structure=kwargs["output_structure"],
            fixed_lengths=kwargs["fixed_lengths"],
            input_structure=kwargs["input_structure"],
            simulate_series=kwargs["simulate_series"],
            mutation_name=kwargs["mutation_name"],
            mutation_targets=kwargs["mutation_targets"],
            snp_input_file=kwargs["snp_input_file"],
            random_snps=kwargs["random_snps"],
            random_snp_density=kwargs["random_snp_density"],
            random_snp_output_file=kwargs["random_snp_output_file"],
            random_snp_region=kwargs["random_snp_region"],
            random_snp_haplotypes=kwargs["random_snp_haplotypes"],
        )

        # Determine output FASTA path
        fasta_file = Path(kwargs["out_dir"]) / f"{kwargs['out_base']}.001.simulated.fa"

        # Step 2: Read simulation (if requested)
        if with_reads:
            logging.info("Pipeline Step 2/3: Simulating %s reads...", with_reads)
            if with_reads == "illumina":
                ctx.invoke(
                    illumina,
                    input_fasta=str(fasta_file),
                    out_dir=kwargs["out_dir"],
                    out_base=f"{kwargs['out_base']}_reads",
                    coverage=kwargs["coverage"],
                    threads=kwargs["threads"],
                )
            else:  # ont
                ctx.invoke(
                    ont,
                    input_fasta=str(fasta_file),
                    out_dir=kwargs["out_dir"],
                    out_base=f"{kwargs['out_base']}_reads",
                    coverage=kwargs["coverage"],
                    min_read_length=100,
                )

        # Step 3: ORF prediction (if requested)
        if with_orfs:
            logging.info("Pipeline Step 3/3: Predicting ORFs...")
            ctx.invoke(
                orfs,
                input_fasta=str(fasta_file),
                out_dir=kwargs["out_dir"],
                out_base=f"{kwargs['out_base']}_orfs",
                orf_min_aa=kwargs["orf_min_aa"],
                orf_aa_prefix=None,
            )

        logging.info("Pipeline completed successfully!")
        ctx.exit(0)

    except Exception as e:
        logging.error("Pipeline failed: %s", e)
        ctx.exit(1)


# ============================================================================
# Helper Functions
# ============================================================================


def _make_args_namespace(config_path, kwargs):
    """Convert Click kwargs to argparse-style namespace (DRY - reuse backend)."""
    from types import SimpleNamespace

    return SimpleNamespace(
        config=config_path,
        out_base=kwargs["out_base"],
        out_dir=kwargs["out_dir"],
        num_haplotypes=kwargs["num_haplotypes"],
        seed=kwargs["seed"],
        reference_assembly=kwargs["reference_assembly"],
        output_structure=kwargs["output_structure"],
        fixed_lengths=list(kwargs["fixed_lengths"]) if kwargs["fixed_lengths"] else None,
        input_structure=kwargs["input_structure"],
        simulate_series=kwargs["simulate_series"],
        mutation_name=kwargs["mutation_name"],
        mutation_targets=list(kwargs["mutation_targets"]) if kwargs["mutation_targets"] else None,
        snp_input_file=kwargs["snp_input_file"],
        random_snps=kwargs["random_snps"],
        random_snp_density=kwargs["random_snp_density"],
        random_snp_output_file=kwargs["random_snp_output_file"],
        random_snp_region=kwargs["random_snp_region"],
        random_snp_haplotypes=kwargs["random_snp_haplotypes"],
        # Pure simulate - no pipeline options
        simulate_reads=None,
        output_orfs=False,
        orf_min_aa=100,
        orf_aa_prefix=None,
    )


# ============================================================================
# Entry Point
# ============================================================================


def main():
    """Entry point for Click CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
