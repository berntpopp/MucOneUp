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
from contextlib import nullcontext
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
# CLI Context Settings
# ============================================================================

# Context settings for custom help flags (-h in addition to --help)
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


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
# Configuration Validation
# ============================================================================


def validate_config_callback(
    ctx: click.Context, param: click.Parameter, value: str | None
) -> str | None:
    """Smart validation for --config option.

    Allows --help and --version to work without config,
    but ensures config is provided when actually executing commands.

    Args:
        ctx: Click context
        param: Click parameter
        value: Config file path (or None)

    Returns:
        Config file path (or None during resilient parsing)
    """
    # Skip validation during resilient parsing (--help, --version)
    if ctx.resilient_parsing:
        return value

    # Config validation will be enforced per-command using require_config()
    # This allows subcommand help to work: muconeup analyze --help
    return value


def require_config(ctx: click.Context) -> None:
    """Ensure config is provided before executing a command.

    Call this at the start of each command that needs config.
    Provides clear error message if config is missing.

    Args:
        ctx: Click context

    Raises:
        click.UsageError: If config is not provided
    """
    if not ctx.obj.get("config_path"):
        raise click.UsageError(
            "Missing required option '--config'.\n"
            "Usage: muconeup --config FILE COMMAND [OPTIONS]\n"
            "Try 'muconeup --help' for more information."
        )


def print_version_callback(ctx: click.Context, param: click.Parameter, value: bool) -> None:
    """Print version and exit when --version or -V is used.

    This callback is marked as eager to ensure it runs before
    required options are validated, allowing --version to work
    without providing --config.

    Args:
        ctx: Click context
        param: Click parameter
        value: Whether the flag was set
    """
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"MucOneUp, version {__version__}")
    ctx.exit()


# ============================================================================
# Root CLI Group
# ============================================================================


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--version",
    "-V",
    is_flag=True,
    callback=print_version_callback,
    expose_value=False,
    is_eager=True,
    help="Show the version and exit.",
)
@click.option(
    "--config",
    callback=validate_config_callback,
    is_eager=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to JSON configuration file.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL", "NONE"]),
    default="INFO",
    help="Set logging level.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Enable verbose output (sets log level to DEBUG).",
)
@click.pass_context
def cli(ctx, config, log_level, verbose):
    """MucOneUp - MUC1 VNTR diploid reference simulator.

    \b
    Philosophy: Each command does ONE thing (Unix philosophy)
      simulate  - Generate haplotypes ONLY
      reads     - Simulate reads from FASTA (supports multiple files)
      analyze   - Analyze FASTA (ORFs, stats) (supports multiple files)

    \b
    Examples:
      # Generate haplotypes
      muconeup --config X simulate --out-base Y

      # Process single file
      muconeup --config X reads illumina Y.001.simulated.fa --out-base reads

      # Process multiple files (batch)
      muconeup --config X reads illumina Y.*.simulated.fa
      muconeup --config X analyze orfs Y.*.simulated.fa

      # Shell composition (Unix philosophy)
      muconeup --config X simulate --fixed-lengths 20-40 --simulate-series 1
      for f in Y.*.fa; do muconeup --config X reads illumina "$f"; done
    """
    # KISS: Simple precedence - verbose overrides log_level
    if verbose:
        log_level = "DEBUG"

    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config  # May be None during help/version
    ctx.obj["log_level"] = log_level

    # Only configure logging if not in resilient parsing mode
    if not ctx.resilient_parsing:
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
    # Validate config is provided
    require_config(ctx)

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
        total_iterations = len(simulation_configs)

        # Show progress bar only for series mode (DRY - no code duplication)
        # Use nullcontext for single iterations to avoid if/else duplication
        progress_ctx = (
            click.progressbar(
                simulation_configs,
                label=f"Simulating {total_iterations} iterations",
                show_eta=True,
                show_pos=True,
            )
            if total_iterations > 1
            else nullcontext(simulation_configs)
        )

        with progress_ctx as configs:
            for sim_index, fixed_conf in enumerate(configs, start=1):
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
        return  # Click handles exit automatically

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
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
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
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed):
    """Simulate Illumina short reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads illumina file.fa --out-base reads
    - Multiple files: muconeup reads illumina file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads illumina *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads illumina sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads illumina sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads illumina sample.*.simulated.fa

      # Compose with shell (Unix philosophy)
      for f in *.fa; do muconeup --config X reads illumina "$f"; done
    """
    # Validate config is provided
    require_config(ctx)

    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure read simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "illumina"
        config["read_simulation"]["coverage"] = coverage
        config["read_simulation"]["threads"] = threads

        # Set seed if provided
        if seed is not None:
            config["read_simulation"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for Illumina read simulation", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_reads")

            logging.info(
                "[%d/%d] Simulating Illumina reads: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            # Run simulation for this file
            simulate_reads_pipeline(config, input_fasta)

        logging.info("Illumina read simulation completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("Read simulation failed: %s", e)
        ctx.exit(1)


@reads.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
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
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed):
    """Simulate Oxford Nanopore long reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads ont file.fa --out-base reads
    - Multiple files: muconeup reads ont file1.fa file2.fa file3.fa
    - Glob pattern: muconeup reads ont *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X reads ont sample.001.fa --out-base my_reads

      # Multiple files (auto-generated output names)
      muconeup --config X reads ont sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X reads ont sample.*.simulated.fa
    """
    # Validate config is provided
    require_config(ctx)

    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure ONT simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "ont"
        config["read_simulation"]["coverage"] = coverage
        if "nanosim_params" not in config:
            config["nanosim_params"] = {}
        config["nanosim_params"]["min_len"] = min_read_length

        # Set seed if provided
        if seed is not None:
            config["nanosim_params"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for ONT read simulation", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_ont_reads")

            logging.info(
                "[%d/%d] Simulating ONT reads: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            # Run simulation for this file
            simulate_reads_pipeline(config, input_fasta)

        logging.info("ONT read simulation completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("Read simulation failed: %s", e)
        ctx.exit(1)


@reads.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
)
@click.option(
    "--coverage",
    type=int,
    default=None,
    help="Target coverage (overrides config if provided).",
)
@click.option(
    "--pass-num",
    type=int,
    default=None,
    help="Number of passes per molecule for multi-pass CLR simulation (≥2, overrides config if provided).",
)
@click.option(
    "--min-passes",
    type=int,
    default=None,
    help="Minimum passes required for CCS HiFi consensus (≥1, overrides config if provided).",
)
@click.option(
    "--min-rq",
    type=float,
    default=None,
    help="Minimum predicted accuracy for HiFi reads (0.0-1.0, overrides config if provided). 0.99 = Q20 (standard HiFi).",
)
@click.option(
    "--model-type",
    type=click.Choice(["qshmm", "errhmm"]),
    default=None,
    help="pbsim3 model type (overrides config if provided).",
)
@click.option(
    "--model-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to pbsim3 model file (overrides config if provided).",
)
@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Number of threads.",
)
@click.option(
    "--seed",
    type=int,
    default=None,
    help="Random seed for reproducibility (same seed = identical reads).",
)
@click.pass_context
def pacbio(
    ctx,
    input_fastas,
    out_dir,
    out_base,
    coverage,
    pass_num,
    min_passes,
    min_rq,
    model_type,
    model_file,
    threads,
    seed,
):
    """Simulate PacBio HiFi reads from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup reads pacbio file.fa --model-file X.model --out-base reads
    - Multiple files: muconeup reads pacbio file1.fa file2.fa --model-file X.model
    - Glob pattern: muconeup reads pacbio *.simulated.fa --model-file X.model

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Workflow:
      1. Multi-pass CLR simulation (pbsim3)
      2. HiFi consensus generation (CCS)
      3. Read alignment (minimap2 with map-hifi preset)

    \b
    Examples:
      # Single file with standard HiFi settings (Q20)
      muconeup --config X reads pacbio sample.001.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --out-base my_hifi

      # Multiple files with high-accuracy HiFi (Q30)
      muconeup --config X reads pacbio sample.*.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --min-rq 0.999 --min-passes 5

      # Ultra-deep coverage simulation
      muconeup --config X reads pacbio sample.fa \\
        --model-file /models/QSHMM-SEQUEL.model \\
        --coverage 100 --pass-num 5

    \b
    Model Files:
      Download from: https://github.com/yukiteruono/pbsim3/tree/master/data
      - QSHMM-SEQUEL.model: Sequel II chemistry
      - QSHMM-RSII.model: RS II chemistry
      - ERRHMM-SEQUEL.model: Alternative error model

    \b
    Quality Control:
      - pass_num ≥2 required for multi-pass (≥3 recommended)
      - min_passes controls CCS stringency (higher = better quality, lower yield)
      - min_rq=0.99 is Q20 (standard HiFi threshold)
      - min_rq=0.999 is Q30 (ultra-high accuracy)
    """
    # Validate config is provided
    require_config(ctx)

    try:
        from ..read_simulation import simulate_reads as simulate_reads_pipeline

        # Load config once (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        # Configure PacBio simulation (shared for all files)
        if "read_simulation" not in config:
            config["read_simulation"] = {}
        config["read_simulation"]["simulator"] = "pacbio"

        # Initialize pacbio_params if not in config (config-first principle)
        if "pacbio_params" not in config:
            config["pacbio_params"] = {}

        # Override config with CLI params only if provided (CLI takes precedence)
        # This follows the ONT pattern: defaults → config → CLI override hierarchy
        if model_type is not None:
            config["pacbio_params"]["model_type"] = model_type
        if model_file is not None:
            config["pacbio_params"]["model_file"] = model_file
        if coverage is not None:
            config["pacbio_params"]["coverage"] = coverage
        if pass_num is not None:
            config["pacbio_params"]["pass_num"] = pass_num
        if min_passes is not None:
            config["pacbio_params"]["min_passes"] = min_passes
        if min_rq is not None:
            config["pacbio_params"]["min_rq"] = min_rq
        if threads != 4:  # threads has explicit default, so keep this check
            config["pacbio_params"]["threads"] = threads
        if seed is not None:
            config["pacbio_params"]["seed"] = seed
            logging.info(f"Using random seed: {seed} (results will be reproducible)")

        # Validate required config parameters exist
        required_params = ["model_type", "model_file"]
        missing_params = [p for p in required_params if p not in config["pacbio_params"]]
        if missing_params:
            raise click.ClickException(
                f"Missing required PacBio parameters in config: {', '.join(missing_params)}. "
                f"Either add them to config.json pacbio_params section or provide via CLI options."
            )

        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for PacBio HiFi read simulation", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_pacbio_hifi")

            logging.info(
                "[%d/%d] Simulating PacBio HiFi reads: %s -> %s",
                idx,
                total_files,
                input_fasta,
                actual_out_base,
            )

            # Run simulation for this file
            simulate_reads_pipeline(config, input_fasta)

        logging.info("PacBio HiFi read simulation completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("PacBio HiFi read simulation failed: %s", e)
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
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
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
def orfs(ctx, input_fastas, out_dir, out_base, orf_min_aa, orf_aa_prefix):
    """Predict ORFs and detect toxic protein features from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup analyze orfs file.fa --out-base analysis
    - Multiple files: muconeup analyze orfs file1.fa file2.fa file3.fa
    - Glob pattern: muconeup analyze orfs *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X analyze orfs sample.001.fa --out-base my_analysis

      # Multiple files (auto-generated output names)
      muconeup --config X analyze orfs sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X analyze orfs sample.*.simulated.fa
    """
    # Validate config is provided
    require_config(ctx)

    try:
        # Load config once for toxic protein detection (DRY principle)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)
        left_const = config.get("constants", {}).get("left")
        right_const = config.get("constants", {}).get("right")

        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for ORF prediction", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_orfs")

            orf_output = Path(out_dir) / f"{actual_out_base}.orfs.fa"
            logging.info(
                "[%d/%d] Running ORF prediction: %s -> %s",
                idx,
                total_files,
                input_fasta,
                orf_output,
            )

            # Run orfipy command-line tool
            # Note: orfipy creates its own output directory, so we must:
            # 1. Use --outdir to specify the output directory
            # 2. Use just the filename (not full path) for --pep
            orf_filename = f"{actual_out_base}.orfs.fa"
            cmd = [
                "orfipy",
                input_fasta,
                "--outdir",
                str(out_dir),
                "--pep",
                orf_filename,
                "--min",
                str(orf_min_aa * 3),  # Convert AA to nucleotides
                "--start",
                "ATG",
            ]

            # Note: Amino acid prefix filtering is done in translate.py after orfipy runs,
            # not via orfipy flags (which don't support this feature)

            result = subprocess.run(cmd, capture_output=True, text=True, check=False)

            if result.returncode != 0:
                logging.error("orfipy failed for %s: %s", input_fasta, result.stderr)
                continue  # Continue with next file instead of exiting

            logging.info("ORF prediction completed: %s", orf_output)

            # Filter by amino acid prefix if specified
            if orf_aa_prefix and orf_output.exists():
                from Bio import SeqIO

                filtered_orfs = []
                total_orfs = 0

                for record in SeqIO.parse(str(orf_output), "fasta"):
                    total_orfs += 1
                    # Check if protein sequence starts with required prefix
                    if str(record.seq).startswith(orf_aa_prefix):
                        filtered_orfs.append(record)

                # Write filtered ORFs back to file
                SeqIO.write(filtered_orfs, str(orf_output), "fasta")
                logging.info(
                    "Filtered ORFs by prefix '%s': %d/%d ORFs retained",
                    orf_aa_prefix,
                    len(filtered_orfs),
                    total_orfs,
                )

            # Toxic protein detection
            if orf_output.exists():
                from ..toxic_protein_detector import scan_orf_fasta

                stats = scan_orf_fasta(
                    str(orf_output), left_const=left_const, right_const=right_const
                )

                stats_file = Path(out_dir) / f"{actual_out_base}.orf_stats.json"
                with stats_file.open("w") as f:
                    json.dump(stats, f, indent=4)

                logging.info("Toxic protein stats written: %s", stats_file)

        logging.info("ORF prediction completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("ORF analysis failed: %s", e)
        ctx.exit(1)


@analyze.command()
@click.argument(
    "input_fastas", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    "--out-dir",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False),
    help="Output folder.",
)
@click.option(
    "--out-base",
    default=None,
    help="Base name for output files (auto-generated if processing multiple files).",
)
@click.pass_context
def stats(ctx, input_fastas, out_dir, out_base):
    """Generate basic sequence statistics from one or more FASTA files.

    Supports batch processing following Unix philosophy:
    - Single file: muconeup analyze stats file.fa --out-base stats
    - Multiple files: muconeup analyze stats file1.fa file2.fa file3.fa
    - Glob pattern: muconeup analyze stats *.simulated.fa

    When processing multiple files, --out-base is auto-generated from input
    filenames unless explicitly provided (which applies to all files).

    \b
    Examples:
      # Single file with custom output name
      muconeup --config X analyze stats sample.001.fa --out-base my_stats

      # Multiple files (auto-generated output names)
      muconeup --config X analyze stats sample.001.fa sample.002.fa

      # Glob pattern (shell expands)
      muconeup --config X analyze stats sample.*.simulated.fa
    """
    try:
        # Warn if --out-base provided for multiple files
        if len(input_fastas) > 1 and out_base:
            logging.warning(
                "--out-base '%s' will be used for all %d files. "
                "Consider omitting --out-base for auto-generated names.",
                out_base,
                len(input_fastas),
            )

        # Process each file (KISS principle - simple iteration)
        total_files = len(input_fastas)
        logging.info("Processing %d FASTA file(s) for statistics generation", total_files)

        for idx, input_fasta in enumerate(input_fastas, start=1):
            # Determine output base name
            if out_base:
                # User provided: use as-is (or append index for multiple files)
                actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
            else:
                # Auto-generate from input filename
                actual_out_base = _generate_output_base(Path(input_fasta), "_stats")

            logging.info("[%d/%d] Generating statistics: %s", idx, total_files, input_fasta)

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

            stats_file = Path(out_dir) / f"{actual_out_base}.basic_stats.json"
            with stats_file.open("w") as f:
                json.dump(stats_data, f, indent=4)

            logging.info("Statistics written: %s", stats_file)

        logging.info("Statistics generation completed for all %d file(s).", total_files)
        return  # Click handles exit automatically

    except Exception as e:
        logging.error("Statistics generation failed: %s", e)
        ctx.exit(1)


# ============================================================================
# ANALYZE vntr-stats - VNTR Transition Probability Analysis
# ============================================================================


@analyze.command("vntr-stats")
@click.argument("input_file", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--structure-column",
    default="vntr",
    show_default=True,
    help="Column name (if header) or 0-based index containing VNTR structure.",
)
@click.option(
    "--delimiter",
    default="\t",
    show_default=True,
    help="Field delimiter for input file.",
)
@click.option(
    "--header",
    is_flag=True,
    help="Specify if input file has header row.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Output JSON file (default: stdout).",
)
@click.pass_context
def vntr_stats(ctx, input_file, structure_column, delimiter, header, output):
    """Analyze VNTR structures and compute transition probabilities.

    Processes a CSV/TSV file containing VNTR structures, calculates statistics
    (min/max/mean/median repeat units), and builds a transition probability
    matrix showing the likelihood of each repeat unit following another.

    The analysis removes duplicate VNTR structures and includes an "END" state
    representing sequence termination. Unknown repeat tokens (not in config)
    trigger warnings but don't cause failure.

    \b
    Examples:
      # Analyze example VNTR database
      muconeup --config X analyze vntr-stats data/examples/vntr_database.tsv --header

      # Use custom column and save to file
      muconeup --config X analyze vntr-stats data.csv \\
        --delimiter "," --structure-column "sequence" --output stats.json

      # Column index without header
      muconeup --config X analyze vntr-stats data.tsv --structure-column 3

      # Pipe to jq for filtering
      muconeup --config X analyze vntr-stats data/examples/vntr_database.tsv \\
        --header | jq '.mean_repeats'

    \b
    Output JSON contains:
      - Statistics: min/max/mean/median repeat counts
      - Probabilities: State transition matrix (including END state)
      - Repeats: Known repeat dictionary from config
    """
    # Validate config is provided
    require_config(ctx)

    try:
        # Load config for known repeats (DRY principle - reuse existing pattern)
        config_path = Path(ctx.obj["config_path"])
        with config_path.open() as f:
            config = json.load(f)

        known_repeats = config.get("repeats", {})
        if not known_repeats:
            logging.error("Configuration file does not contain 'repeats' key")
            ctx.exit(1)

        logging.info("Analyzing VNTR structures from %s", input_file)

        # Import analysis module
        from ..analysis.vntr_statistics import analyze_vntr_sequences

        # Run analysis
        with Path(input_file).open() as f:
            results = analyze_vntr_sequences(f, structure_column, header, known_repeats, delimiter)

        # Build output (matches original helper behavior)
        output_data = {
            "min_repeats": results["min_repeats"],
            "max_repeats": results["max_repeats"],
            "mean_repeats": results["mean_repeats"],
            "median_repeats": results["median_repeats"],
            "repeats": known_repeats,
            "probabilities": results["probabilities"],
        }

        # Write output
        if output:
            output_path = Path(output)
            with output_path.open("w") as out:
                json.dump(output_data, out, indent=2)
            logging.info("VNTR statistics written to %s", output)
        else:
            # Output to stdout (pipeable)
            click.echo(json.dumps(output_data, indent=2))

        logging.info(
            "Analysis complete: %d unique structures, %d repeat types",
            len(results["probabilities"]),
            len([t for t in results["probabilities"] if t != "END"]),
        )

    except Exception as e:
        logging.error("VNTR analysis failed: %s", e, exc_info=True)
        ctx.exit(1)


# ============================================================================
# Helper Functions
# ============================================================================


def _generate_output_base(input_path: Path, suffix: str) -> str:
    """Generate output base name from input file path.

    Args:
        input_path: Path to input FASTA file
        suffix: Suffix to append (e.g., '_reads', '_orfs')

    Returns:
        Output base name (e.g., 'sample.001.simulated_reads')

    Examples:
        >>> _generate_output_base(Path('sample.001.simulated.fa'), '_reads')
        'sample.001.simulated_reads'
        >>> _generate_output_base(Path('/path/sample.fa'), '_orfs')
        'sample_orfs'
    """
    return input_path.stem + suffix


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


@analyze.command()
@click.argument("input_fasta", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--mutation",
    required=True,
    help="Mutation name to validate (e.g., 'dupC').",
)
@click.option(
    "--output",
    type=click.Path(),
    help="Output JSON file for validation results (prints to stdout if not specified).",
)
@click.pass_context
def snapshot_validate(ctx, input_fasta, mutation, output):
    """Validate SNaPshot assay for MUC1 VNTR mutations.

    Simulates complete SNaPshot workflow: PCR → MwoI digest → extension → detection.

    \b
    Examples:
      # Validate dupC mutation in a sample
      muconeup --config config.json analyze snapshot-validate sample.fa --mutation dupC

      # Save results to JSON
      muconeup --config config.json analyze snapshot-validate sample.fa --mutation dupC --output results.json
    """
    # Validate config is provided
    require_config(ctx)

    from Bio import SeqIO

    from ..analysis.snapshot_validator import SnapshotValidator
    from ..config import load_config

    try:
        # Load config
        config_path = Path(ctx.obj["config_path"])
        config = load_config(str(config_path))

        # Check if mutation is configured
        if "snapshot_validation" not in config:
            logging.error("Config missing 'snapshot_validation' section")
            ctx.exit(1)

        if mutation not in config["snapshot_validation"]:
            available = list(config["snapshot_validation"].keys())
            logging.error(
                "Mutation '%s' not found in config. Available: %s",
                mutation,
                available,
            )
            ctx.exit(1)

        # Initialize validator
        validator = SnapshotValidator(config, mutation)
        logging.info("SNaPshot validator initialized for mutation: %s", mutation)

        # Load input FASTA
        records = list(SeqIO.parse(input_fasta, "fasta"))
        if not records:
            logging.error("No sequences found in %s", input_fasta)
            ctx.exit(1)

        logging.info("Loaded %d haplotype(s) from %s", len(records), input_fasta)

        # Run validation on each haplotype
        all_results = {}
        for i, record in enumerate(records, 1):
            haplotype_name = record.id or f"haplotype_{i}"
            template = str(record.seq)

            logging.info(
                "Validating haplotype %d: %s (%d bp)",
                i,
                haplotype_name,
                len(template),
            )

            result = validator.validate_complete_workflow(template)

            all_results[haplotype_name] = {
                "mutation_detected": result["mutation_detected"],
                "pcr_products": result["pcr_results"]["product_count"],
                "digest_survivors": result["digest_results"]["survivor_count"],
                "expected_fluorescence": result.get("expected_fluorescence", ""),
                "summary": result["summary"],
            }

            # Add SNaPshot details if mutation detected
            if result["mutation_detected"] and result.get("snapshot_results"):
                snapshot_info = []
                for snap in result["snapshot_results"]:
                    if snap["extension"].get("binds"):
                        snapshot_info.append(
                            {
                                "fluorophore": snap["extension"].get("fluorophore"),
                                "dye": snap["extension"].get("fluorophore_dye"),
                                "next_base": snap["extension"].get("next_base"),
                            }
                        )
                all_results[haplotype_name]["snapshot_details"] = snapshot_info

            logging.info("  Result: %s", result["summary"])

        # Output results
        output_data = {
            "mutation": mutation,
            "input_file": str(input_fasta),
            "haplotypes": all_results,
            "overall_detection": any(r["mutation_detected"] for r in all_results.values()),
        }

        if output:
            output_path = Path(output)
            with output_path.open("w") as f:
                json.dump(output_data, f, indent=2)
            logging.info("Validation results written to %s", output)
        else:
            click.echo(json.dumps(output_data, indent=2))

        # Exit code: 0 if mutation detected, 1 if not
        if output_data["overall_detection"]:
            logging.info("✓ Mutation %s DETECTED in sample", mutation.upper())
            ctx.exit(0)
        else:
            logging.info("✗ Mutation %s NOT detected in sample", mutation.upper())
            ctx.exit(1)

    except Exception as e:
        logging.error("SNaPshot validation failed: %s", e, exc_info=True)
        ctx.exit(1)


# ============================================================================
# Entry Point
# ============================================================================


def main():
    """Entry point for Click CLI."""
    cli(obj={})


if __name__ == "__main__":
    main()
