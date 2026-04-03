"""Reads commands -- simulate sequencing reads from FASTA files."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import click

from .._common import require_config
from ..error_handling import cli_error_handler

# ============================================================================
# Shared batch helper
# ============================================================================


def _run_batch_simulation(
    config: dict[str, Any],
    input_fastas: tuple[str, ...],
    out_dir: str,
    out_base: str | None,
    suffix: str,
    simulator_label: str,
    track_read_source: bool,
) -> None:
    """Run read simulation over one or more FASTA files.

    Shared loop for illumina, ont, and pacbio commands. Handles output naming,
    source tracker reconstruction, and pipeline dispatch.
    """
    from ...read_simulation import simulate_reads as simulate_reads_pipeline
    from ...read_simulator.output_config import OutputConfig

    # Warn if --out-base provided for multiple files
    if len(input_fastas) > 1 and out_base:
        logging.warning(
            "--out-base '%s' will be used for all %d files. "
            "Consider omitting --out-base for auto-generated names.",
            out_base,
            len(input_fastas),
        )

    total_files = len(input_fastas)
    logging.info("Processing %d FASTA file(s) for %s read simulation", total_files, simulator_label)

    for idx, input_fasta in enumerate(input_fastas, start=1):
        # Determine output base name
        if out_base:
            actual_out_base = f"{out_base}_{idx:03d}" if total_files > 1 else out_base
        else:
            actual_out_base = None

        output_config = OutputConfig.from_input_fasta(
            input_fa=input_fasta,
            out_dir=out_dir,
            out_base=actual_out_base,
            suffix=suffix,
        )

        logging.info(
            "[%d/%d] Simulating %s reads: %s -> %s",
            idx,
            total_files,
            simulator_label,
            input_fasta,
            output_config.out_base,
        )

        # Build source tracker from companion files if requested
        source_tracker = None
        if track_read_source:
            from ...read_simulator.source_tracking import ReadSourceTracker

            # Derive stats path: input may be "base.001.simulated.fa" but stats
            # file is "base.001.simulation_stats.json" (no ".simulated" component).
            input_p = Path(input_fasta)
            if input_p.stem.endswith(".simulated"):
                stats_path = str(
                    input_p.with_name(
                        input_p.stem.removesuffix(".simulated") + ".simulation_stats.json"
                    )
                )
            else:
                stats_path = str(input_p.with_suffix(".simulation_stats.json"))
            source_tracker = ReadSourceTracker.from_companion_files(stats_path)
            if source_tracker is None:
                logging.warning("Could not reconstruct read source tracker from %s", stats_path)
            else:
                coord_map_path = output_config.derive_path_str("_repeat_coordinates.tsv")
                source_tracker.write_coordinate_map(coord_map_path)
                logging.info("Repeat coordinate map written: %s", coord_map_path)

        # Run simulation for this file
        simulate_reads_pipeline(
            config,
            input_fasta,
            source_tracker=source_tracker,
            output_config=output_config,
        )

    logging.info("%s read simulation completed for all %d file(s).", simulator_label, total_files)


# ============================================================================
# READS Command Group - Pure Read Simulation Utilities
# ============================================================================


@click.group()
def reads():
    """Read simulation utilities.

    \b
    Simulate reads from ANY FASTA file.
    Works with MucOneUp outputs or external sequences.
    Requires --config at the top level.
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
    default=None,
    help="Target sequencing coverage (overrides config if provided, defaults to config value or 30x).",
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
@click.option(
    "--track-read-source",
    is_flag=True,
    default=False,
    help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
)
@click.pass_context
@cli_error_handler
def illumina(ctx, input_fastas, out_dir, out_base, coverage, threads, seed, track_read_source):
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
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    # Configure Illumina simulation
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = "illumina"

    if coverage is not None:
        config["read_simulation"]["coverage"] = coverage
    elif "coverage" not in config.get("read_simulation", {}):
        config["read_simulation"]["coverage"] = 30

    config["read_simulation"]["threads"] = threads

    if seed is not None:
        config["read_simulation"]["seed"] = seed
        logging.info(f"Using random seed: {seed} (results will be reproducible)")

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_reads", "Illumina", track_read_source
    )


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
    help="Target coverage (overrides config if provided, defaults to config value or 30x).",
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
@click.option(
    "--track-read-source",
    is_flag=True,
    default=False,
    help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
)
@click.pass_context
@cli_error_handler
def ont(ctx, input_fastas, out_dir, out_base, coverage, min_read_length, seed, track_read_source):
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
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    # Configure ONT simulation
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = "ont"

    if coverage is not None:
        config["read_simulation"]["coverage"] = coverage
    elif "coverage" not in config.get("read_simulation", {}):
        config["read_simulation"]["coverage"] = 30

    if "nanosim_params" not in config:
        config["nanosim_params"] = {}
    config["nanosim_params"]["min_len"] = min_read_length

    if seed is not None:
        config["nanosim_params"]["seed"] = seed
        logging.info(f"Using random seed: {seed} (results will be reproducible)")

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_ont_reads", "ONT", track_read_source
    )


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
    help="Number of passes per molecule for multi-pass CLR simulation (>=2, overrides config if provided).",
)
@click.option(
    "--min-passes",
    type=int,
    default=None,
    help="Minimum passes required for CCS HiFi consensus (>=1, overrides config if provided).",
)
@click.option(
    "--min-rq",
    type=float,
    default=None,
    help="Minimum predicted read quality (RQ) score (0.0-1.0). 0.99=Q20 (standard HiFi), 0.999=Q30.",
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
@click.option(
    "--track-read-source",
    is_flag=True,
    default=False,
    help="Generate read source tracking manifest and coordinate map alongside simulated reads.",
)
@click.pass_context
@cli_error_handler
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
    track_read_source,
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
      - pass_num >=2 required for multi-pass (>=3 recommended)
      - min_passes controls CCS stringency (higher = better quality, lower yield)
      - min_rq=0.99 is Q20 (standard HiFi threshold)
      - min_rq=0.999 is Q30 (ultra-high accuracy)
    """
    require_config(ctx)

    from ...config import load_config_raw

    config = load_config_raw(str(ctx.obj["config_path"]))

    # Configure PacBio simulation
    if "read_simulation" not in config:
        config["read_simulation"] = {}
    config["read_simulation"]["simulator"] = "pacbio"

    if "pacbio_params" not in config:
        config["pacbio_params"] = {}

    # Override config with CLI params only if provided
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
    if threads != 4:
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

    _run_batch_simulation(
        config, input_fastas, out_dir, out_base, "_pacbio_hifi", "PacBio HiFi", track_read_source
    )
